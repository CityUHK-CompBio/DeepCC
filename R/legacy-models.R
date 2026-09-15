# Legacy Keras 2 model support
#
# Models published before the keras3 migration were written by Keras 2.x.
# Keras 3 refuses two entries that appear in those files:
#
#   * `batch_input_shape` on the first layer, where Keras 3 expects an Input
#     layer or `batch_shape`;
#   * `axis` stored as a length-1 list on BatchNormalization layers, where
#     Keras 3 requires a scalar.
#
# The stored weights are unaffected, so a legacy file is rebuilt from the
# architecture it records, with those two entries converted, and then
# populated with the original weights.

#' Convert a class name to snake_case
#' @noRd
legacySnakeCase <- function(x) {
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- gsub("([A-Z]+)([A-Z][a-z])", "\\1_\\2", x)
  tolower(x)
}

#' Map a serialized Keras 2 object to its Keras 3 alias
#' @noRd
legacyObjectAlias <- function(value) {
  if (is.null(value) || is.character(value)) return(value)
  if (is.list(value) && !is.null(value[["class_name"]])) {
    return(legacySnakeCase(value[["class_name"]]))
  }
  value
}

#' Convert a Keras 2 layer config into one Keras 3 accepts
#'
#' @param class_name layer class recorded in the file
#' @param config layer config recorded in the file
#' @param rank number of input dimensions including the batch dimension.
#'   Keras 2 counted axes from the batch dimension, while Keras 3 counts
#'   feature axes and addresses the last one as -1, so `axis` is remapped.
#' @noRd
convertLegacyLayerConfig <- function(class_name, config, rank = NULL) {
  obsolete <- c("batch_input_shape", "batch_shape", "input_shape", "input_dim", "dtype")
  config <- config[setdiff(names(config), obsolete)]

  if (identical(class_name, "BatchNormalization")) {
    axis <- config[["axis"]]
    if (is.list(axis)) axis <- axis[[1L]]
    axis <- as.integer(if (is.null(axis)) -1L else axis)
    if (!is.null(rank) && axis > 0L && axis == rank - 1L) {
      # The normalised axis is the last feature axis. Keras 3 rejects the
      # Keras 2 index here because it leaves the variable shape undefined.
      axis <- -1L
    }
    config[["axis"]] <- axis
  }

  kept <- list()
  for (key in names(config)) {
    value <- config[[key]]
    if (grepl("_initializer$|_regularizer$|_constraint$", key)) {
      value <- legacyObjectAlias(value)
    }
    # NULL entries are unset regularizers, constraints, or activations, whose
    # Keras 3 defaults carry the same meaning.
    if (!is.null(value)) kept[[key]] <- value
  }
  kept
}

#' Look up the keras3 constructor for a legacy layer class
#' @noRd
legacyLayerConstructor <- function(class_name) {
  fn_name <- paste0("layer_", legacySnakeCase(class_name))
  ns <- asNamespace("keras3")
  if (!exists(fn_name, envir = ns, inherits = FALSE)) {
    stop(sprintf("Unsupported layer class '%s' in the saved model.", class_name),
         call. = FALSE)
  }
  get(fn_name, envir = ns)
}

#' Is h5py available in the Python runtime keras3 uses?
#'
#' `reticulate` performs its own interpreter discovery, which can disagree
#' with the environment keras3 was configured with. Asking keras3 for its
#' interpreter first keeps both libraries on the same runtime.
#' @noRd
kerasH5pyAvailable <- function() {
  # config_backend() activates the Python runtime keras3 is configured with,
  # so h5py is looked up there instead of in reticulate's own interpreter
  # discovery, which can resolve to a different environment.
  try(keras3::config_backend(), silent = TRUE)
  available <- tryCatch(reticulate::py_module_available("h5py"),
                        error = function(e) FALSE)
  isTRUE(available)
}

#' Read the architecture recorded in a Keras 2 HDF5 file
#' @noRd
readLegacyKerasConfig <- function(filepath) {
  if (!kerasH5pyAvailable()) {
    stop(paste("Reading this model requires the Python 'h5py' module.",
               "Install it with reticulate::py_install(\"h5py\")."),
         call. = FALSE)
  }
  h5py <- reticulate::import("h5py", delay_load = TRUE)
  json <- reticulate::import("json", delay_load = TRUE)
  handle <- h5py$File(filepath, "r")
  on.exit(handle$close(), add = TRUE)

  raw <- tryCatch(handle$attrs[["model_config"]], error = function(e) NULL)
  if (is.null(raw)) return(NULL)
  reticulate::py_to_r(json$loads(raw))
}

#' Detect whether an HDF5 model file was written by Keras 2
#' @noRd
isLegacyKerasHdf5 <- function(filepath) {
  if (!kerasH5pyAvailable()) return(FALSE)
  version <- tryCatch({
    h5py <- reticulate::import("h5py", delay_load = TRUE)
    handle <- h5py$File(filepath, "r")
    on.exit(handle$close(), add = TRUE)
    as.character(handle$attrs[["keras_version"]])
  }, error = function(e) NULL)
  !is.null(version) && length(version) == 1L && grepl("^2\\.", version)
}

#' Rebuild a Keras 2 model under Keras 3 and restore its weights
#' @noRd
loadLegacyKerasModel <- function(filepath) {
  config <- readLegacyKerasConfig(filepath)
  if (is.null(config)) {
    stop("The saved model does not record an architecture.", call. = FALSE)
  }

  specs <- config[["config"]][["layers"]]
  if (is.null(specs) || !length(specs)) {
    stop("The saved model records no layers.", call. = FALSE)
  }

  shape <- config[["config"]][["build_input_shape"]]
  if (is.null(shape)) shape <- specs[[1L]][["config"]][["batch_input_shape"]]
  shape <- unlist(shape)
  shape <- shape[!is.na(shape)]
  if (!length(shape)) {
    stop("The saved model does not record an input shape.", call. = FALSE)
  }

  input <- keras3::keras_input(shape = as.integer(shape))
  x <- input
  for (spec in specs) {
    ctor <- legacyLayerConstructor(spec[["class_name"]])
    args <- convertLegacyLayerConfig(spec[["class_name"]], spec[["config"]],
                                     rank = length(x$shape))
    x <- do.call(ctor, c(list(x), args))
  }
  model <- keras3::keras_model(input, x)
  keras3::load_model_weights(model, filepath)
  model
}

#' Load a classifier from either a Keras 3 or a legacy Keras 2 file
#'
#' The Keras 3 loader runs first so that keras3 initialises its Python
#' runtime before the compatibility loader needs h5py to read the recorded
#' architecture. Both failures are reported when neither loader succeeds,
#' so the underlying cause stays visible.
#' @noRd
loadClassifierAnyFormat <- function(filepath) {
  # Loading a legacy file through the Keras 3 loader fails only after Keras
  # has emitted noisy Python warnings, so pick the loader by recorded format
  # whenever h5py is available to inspect it.
  legacy_first <- isLegacyKerasHdf5(filepath)
  attempts <- if (legacy_first) {
    list(list(loader = loadLegacyKerasModel,
              label = "Keras 2 compatibility loader"),
         list(loader = keras3::load_model, label = "Keras 3 loader"))
  } else {
    list(list(loader = keras3::load_model, label = "Keras 3 loader"),
         list(loader = loadLegacyKerasModel,
              label = "Keras 2 compatibility loader"))
  }
  failures <- character()
  for (attempt in attempts) {
    result <- tryCatch(attempt$loader(filepath), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    failures <- c(failures, sprintf("  %s: %s", attempt$label, conditionMessage(result)))
  }
  stop(sprintf("Could not load '%s'.\n%s", filepath, paste(failures, collapse = "\n")),
       call. = FALSE)
}
