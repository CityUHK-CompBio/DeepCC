# Regression tests for models saved by Keras 2 before the keras3 migration.
#
# The published pre-trained models are ~187 MB each, so these tests instead
# synthesise a Keras 2 HDF5 file with the same structure and check that
# loading it reproduces the predictions of the equivalent Keras 3 model.

h5pyAvailable <- function() DeepCC:::kerasH5pyAvailable()

pyArray <- function(x) reticulate::np_array(x)

legacyModelJson <- function(input_dim, hidden, classes) {
  paste0(
    '{"class_name": "Sequential", "config": {"name": "sequential", "layers": [',
    '{"class_name": "Dense", "config": {"name": "dense", "trainable": true, ',
    '"batch_input_shape": [null, ', input_dim, '], "dtype": "float32", "units": ',
    hidden, ', "activation": "selu", "use_bias": true, ',
    '"kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, ',
    '"bias_initializer": {"class_name": "Zeros", "config": {}}, ',
    '"kernel_regularizer": null, "bias_regularizer": null, ',
    '"activity_regularizer": null, "kernel_constraint": null, ',
    '"bias_constraint": null}}, ',
    '{"class_name": "BatchNormalization", "config": {"name": "batch_normalization", ',
    '"trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, ',
    '"epsilon": 0.001, "center": true, "scale": true, ',
    '"beta_initializer": {"class_name": "Zeros", "config": {}}, ',
    '"gamma_initializer": {"class_name": "Ones", "config": {}}, ',
    '"moving_mean_initializer": {"class_name": "Zeros", "config": {}}, ',
    '"moving_variance_initializer": {"class_name": "Ones", "config": {}}, ',
    '"beta_regularizer": null, "gamma_regularizer": null, ',
    '"beta_constraint": null, "gamma_constraint": null}}, ',
    '{"class_name": "GaussianDropout", "config": {"name": "gaussian_dropout", ',
    '"trainable": true, "dtype": "float32", "rate": 0.4}}, ',
    '{"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, ',
    '"dtype": "float32", "units": ', classes, ', "activation": "softmax", ',
    '"use_bias": true, ',
    '"kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, ',
    '"bias_initializer": {"class_name": "Zeros", "config": {}}, ',
    '"kernel_regularizer": null, "bias_regularizer": null, ',
    '"activity_regularizer": null, "kernel_constraint": null, ',
    '"bias_constraint": null}}], ',
    '"build_input_shape": [null, ', input_dim, ']}}'
  )
}

writeLegacyHdf5 <- function(path, layers, layer_order, input_dim, hidden, classes) {
  h5py <- reticulate::import("h5py", delay_load = TRUE)
  handle <- h5py$File(path, "w")
  on.exit(handle$close(), add = TRUE)

  handle$attrs$create("keras_version", "2.3.0-tf")
  handle$attrs$create("backend", "tensorflow")
  handle$attrs$create("model_config", legacyModelJson(input_dim, hidden, classes))

  weights_group <- handle$create_group("model_weights")
  # Keras 2 recorded every layer here, including weight-free ones such as
  # GaussianDropout, and stored the weight paths on each layer group.
  weights_group$attrs$create("layer_names", layer_order)
  for (name in names(layers)) {
    weights <- layers[[name]]
    layer_group <- weights_group$create_group(name)
    # Weight-free layers get an empty group, matching what Keras 2 wrote.
    if (!length(weights)) next
    weight_names <- paste0(name, "/", names(weights), ":0")
    layer_group$attrs$create("weight_names", weight_names)
    for (i in seq_along(weight_names)) {
      layer_group$create_dataset(weight_names[[i]], data = pyArray(weights[[i]]))
    }
  }
  invisible(path)
}

test_that("a Keras 2 model is rebuilt and reproduces Keras 3 predictions", {
  testthat::skip_if_not(h5pyAvailable(), "h5py is not available")

  input_dim <- 5L
  hidden <- 8L
  classes <- 3L

  keras3::clear_session()
  input <- keras3::keras_input(shape = input_dim, name = "input_layer")
  output <- input |>
    keras3::layer_dense(units = hidden, activation = "selu", name = "dense") |>
    keras3::layer_batch_normalization(name = "batch_normalization") |>
    keras3::layer_gaussian_dropout(rate = 0.4, name = "gaussian_dropout") |>
    keras3::layer_dense(units = classes, activation = "softmax", name = "dense_1")
  reference <- keras3::keras_model(input, output)

  # Non-default batch normalisation state so swapped or missing weights show up.
  reference$layers[[3]]$set_weights(unname(list(
    pyArray(seq(0.5, 1.2, length.out = hidden)),
    pyArray(seq(-0.2, 0.2, length.out = hidden)),
    pyArray(seq(-1, 1, length.out = hidden)),
    pyArray(seq(0.5, 2, length.out = hidden))
  )))

  set.seed(11)
  probe <- matrix(rnorm(4L * input_dim), nrow = 4L)
  expected <- stats::predict(reference, probe, verbose = 0)

  layers <- list(
    dense = list(kernel = reference$layers[[2]]$get_weights()[[1]],
                 bias = reference$layers[[2]]$get_weights()[[2]]),
    batch_normalization = list(gamma = reference$layers[[3]]$get_weights()[[1]],
                               beta = reference$layers[[3]]$get_weights()[[2]],
                               moving_mean = reference$layers[[3]]$get_weights()[[3]],
                               moving_variance = reference$layers[[3]]$get_weights()[[4]]),
    gaussian_dropout = list(),
    dense_1 = list(kernel = reference$layers[[5]]$get_weights()[[1]],
                   bias = reference$layers[[5]]$get_weights()[[2]])
  )

  prefix <- file.path(tempdir(), "legacy_synthetic")
  writeLegacyHdf5(paste0(prefix, ".hdf5"), layers, names(layers),
                  input_dim, hidden, classes)
  levels <- c("A", "B", "C")
  save(levels, file = paste0(prefix, ".RData"))

  expect_true(DeepCC:::isLegacyKerasHdf5(paste0(prefix, ".hdf5")))

  model <- load_DeepCC_model(prefix)
  expect_identical(model$levels, levels)

  got <- get_DeepCC_prob(model, probe)
  expect_equal(unname(got), unname(expected), tolerance = 1e-5)
  expect_identical(colnames(got), levels)
})

test_that("legacy axis and initializer entries are converted", {
  dense <- DeepCC:::convertLegacyLayerConfig("Dense", list(
    name = "dense", batch_input_shape = list(NULL, 5L), dtype = "float32",
    units = 8L, activation = "selu",
    kernel_initializer = list(class_name = "GlorotUniform", config = list(seed = NULL)),
    bias_initializer = list(class_name = "Zeros", config = list()),
    kernel_regularizer = NULL, bias_constraint = NULL
  ))
  expect_false("batch_input_shape" %in% names(dense))
  expect_false("dtype" %in% names(dense))
  expect_false("kernel_regularizer" %in% names(dense))
  expect_identical(dense$kernel_initializer, "glorot_uniform")
  expect_identical(dense$bias_initializer, "zeros")

  # Keras 2 counted the batch dimension, so axis 1 is the last feature axis
  # of a rank-2 input and must become -1 for Keras 3.
  bn <- DeepCC:::convertLegacyLayerConfig("BatchNormalization",
                                          list(name = "bn", axis = list(1L)),
                                          rank = 2L)
  expect_identical(bn$axis, -1L)
})

test_that("unloadable files report both loader failures", {
  testthat::skip_if_not(h5pyAvailable(), "h5py is not available")

  prefix <- file.path(tempdir(), "not_a_model")
  writeLines("not an hdf5 file", paste0(prefix, ".hdf5"))
  levels <- c("A")
  save(levels, file = paste0(prefix, ".RData"))

  expect_error(load_DeepCC_model(prefix), "Keras 3 loader")
})
