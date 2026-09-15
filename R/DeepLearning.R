#' Train DeepCC Model
#'
#' This function trains DeepCC Model on the training data using the modern
#' \code{keras3} interface. The network architecture, training recipe, and
#' legacy model format are preserved from DeepCC 0.1.1.
#'
#' @param trainData a data.frame containing functional spectra of training data (each row presents one sample)
#' @param trainLabels a character vector containing lables of training data
#' @param epochs the number of epochs
#' @param dropout dropout rate
#' @param activation_func activation funtion
#' @param validation_split fraction of training data to use for validation
#' @return a trained DeepCC model with \code{classifier}, \code{levels}, and
#'   \code{feature_names} fields
#' @export
#' @examples
#' \dontrun{
#' set.seed(42)
#' eps <- as.data.frame(matrix(rnorm(20*50), nrow=20, ncol=50))
#' colnames(eps) <- paste0("F", seq_len(50))
#' labels <- sample(c("A", "B", "C"), 20, replace=TRUE)
#' deepcc_model <- train_DeepCC_model(eps, labels, epochs=2)
#' }
train_DeepCC_model <- function(trainData, trainLabels, epochs = 100, dropout = 0.4, activation_func = "selu", validation_split = 0.2){
  trainData <- as.matrix(trainData)
  feature_names <- colnames(trainData)

  ind <- !is.na(trainLabels)
  x_train <- trainData[ind, , drop = FALSE]
  y_train <- factor(trainLabels[ind])
  levels <- levels(y_train)
  class <- length(levels)
  y_train <- keras3::to_categorical(as.integer(y_train) - 1L, class)

  keras3::clear_session()

  init_methods <- "glorot_uniform"
  model <- keras3::keras_model_sequential()
  input_layer <- keras3::keras_input(shape = ncol(x_train))
  output_layer <- input_layer |>
    keras3::layer_dense(units = 1024, activation = activation_func, kernel_initializer = init_methods) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_gaussian_dropout(rate = dropout) |>
    keras3::layer_dense(units = 256, activation = activation_func, kernel_initializer = init_methods) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_gaussian_dropout(rate = dropout) |>
    keras3::layer_dense(units = 64, activation = activation_func, kernel_initializer = init_methods) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_gaussian_dropout(rate = dropout) |>
    keras3::layer_dense(units = 64, activation = activation_func, kernel_initializer = init_methods) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_gaussian_dropout(rate = dropout) |>
    keras3::layer_dense(units = 10, activation = activation_func, kernel_initializer = init_methods) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_gaussian_dropout(rate = dropout) |>
    keras3::layer_dense(units = class, activation = 'softmax')
  model <- keras3::keras_model(inputs = input_layer, outputs = output_layer)

  model <- keras3::compile(model,
    loss = "categorical_crossentropy",
    optimizer = keras3::optimizer_adam(learning_rate = 0.001, beta_1 = 0.9, beta_2 = 0.999),
    metrics = c('accuracy')
  )

  history <- keras3::fit(model,
    x_train, y_train,
    epochs = epochs, batch_size = 1024,
    validation_split = validation_split,
    verbose = 0
  )

  model <- keras3::compile(model,
    loss = "categorical_crossentropy",
    optimizer = keras3::optimizer_sgd(learning_rate = 1e-05, momentum = 0.9),
    metrics = c('accuracy')
  )

  history <- keras3::fit(model,
    x_train, y_train,
    epochs = epochs, batch_size = 1024,
    validation_split = validation_split,
    verbose = 0
  )

  list(classifier = model, levels = levels, feature_names = feature_names)
}

#' Save DeepCC Model
#'
#' @param deepcc_model a trained DeepCC model
#' @param prefix file path prefix; outputs \code{prefix.hdf5} and \code{prefix.RData}
#' @export
save_DeepCC_model <- function(deepcc_model, prefix) {
  keras3::save_model(deepcc_model$classifier, paste0(prefix, ".hdf5"))
  levels <- deepcc_model$levels
  feature_names <- deepcc_model$feature_names
  save(levels, feature_names, file = paste0(prefix, ".RData"))
}

#' Load DeepCC Model
#'
#' Loads a saved DeepCC model. Supports both the new format (with metadata)
#' and the legacy 0.1.1 format (without metadata).
#'
#' @param prefix file path prefix
#' @return a DeepCC model with \code{classifier}, \code{levels}, and
#'   optionally \code{feature_names}
#' @export
load_DeepCC_model <- function(prefix){
  load(file = paste0(prefix, ".RData"))
  classifier <- keras3::load_model(paste0(prefix, ".hdf5"))
  model <- list(classifier = classifier, levels = levels)
  if (exists("feature_names", envir = environment())) {
    fn <- get("feature_names", envir = environment())
    if (!is.null(fn)) model$feature_names <- fn
  }
  model
}

#' Get DeepCC Labels
#'
#' This function classifys new data set using trained DeepCC model.
#'
#' @param DeepCCModel a trained DeepCC model
#' @param newData a data.frame containing functional spectra of new data (each presnets one sample)
#' @param cutoff a numeric indicating cutoff of poster probability
#' @param prob_mode a logical flag; if TRUE, return a data.frame with labels and probabilities
#' @param prob_raw a logical flag; if TRUE and prob_mode is TRUE, return the raw probability matrix
#' @return a character vector containing lables of training data
#' @export
get_DeepCC_label <- function(DeepCCModel, newData, cutoff = 0.5, prob_mode = FALSE, prob_raw = FALSE)
{
  newData <- alignNewData(DeepCCModel, newData)
  res <- stats::predict(DeepCCModel$classifier, as.matrix(newData))
  predicted <- apply(res, 1, function(z){
    if (max(z) >= cutoff){
      which.max(z)
    }
    else {
      NA_integer_
    }
  })
  pred <- factor(predicted, levels = seq_along(DeepCCModel$levels),
                 labels = DeepCCModel$levels)
  if (prob_mode) {
    pred <- data.frame(DeepCC = as.character(pred),
                       Probability = round(apply(res, 1, max), digits = 3))
  }
  if (prob_mode && prob_raw) {
    pred <- res
  }
  pred
}

#' Get DeepCC Prob Matrix
#'
#' @param DeepCCModel a trained DeepCC model
#' @param newData a data.frame containing functional spectra of new data
#' @return a matrix containing class probabilities for each sample
#' @export
get_DeepCC_prob <- function(DeepCCModel, newData){
  newData <- alignNewData(DeepCCModel, newData)
  res <- stats::predict(DeepCCModel$classifier, as.matrix(newData))
  colnames(res) <- DeepCCModel$levels
  res
}

#' Get DeepCC Features
#'
#' This function obtains DeepCC Features from functional spectra using the
#' second-to-last layer of the classifier in inference mode.
#'
#' @param DeepCCModel a trained DeepCC model
#' @param fs a data.frame containing functional spectra (each row presents one sample)
#' @return a data.frame containing DeepCC Features extracted from the second-to-last layer
#' @export
get_DeepCC_features <- function(DeepCCModel, fs) {
  fs <- alignNewData(DeepCCModel, fs)
  model <- DeepCCModel$classifier
  intermediate_layer_model <- keras3::keras_model(inputs = model$input,
                                                  outputs = model$layers[[length(model$layers) - 1]]$output)
  df <- stats::predict(intermediate_layer_model, as.matrix(fs))
  rownames(df) <- rownames(fs)
  df
}

#' Align new data columns to model feature order
#' @noRd
alignNewData <- function(DeepCCModel, newData) {
  if (is.null(DeepCCModel$feature_names)) return(newData)
  if (!is.data.frame(newData) && !is.matrix(newData)) stop("newData must be a data.frame or matrix.")
  if (is.null(colnames(newData))) {
    warning("newData has no column names; using positional order for compatibility.")
    return(newData)
  }
  missing <- setdiff(DeepCCModel$feature_names, colnames(newData))
  if (length(missing)) {
    stop(paste("Missing required features:", paste(utils::head(missing, 10), collapse = ", "),
               if (length(missing) > 10) paste0(" (and ", length(missing) - 10, " more)") else ""))
  }
  newData <- newData[, DeepCCModel$feature_names, drop = FALSE]
  if (is.data.frame(newData)) as.data.frame(newData) else newData
}
