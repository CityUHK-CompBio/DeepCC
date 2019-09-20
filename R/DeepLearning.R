#' Train DeepCC Model
#'
#' This function trains DeepCC Model on the training data.
#'
#' @param trainData a data.frame containing functional spectra of training data (each row presents one sample)
#' @param trainLabels a character vector containing lables of training data
#' @param epochs the number of epochs
#' @param dropout dropout rate
#' @param activation_func activation function
#' @return a trained DeepCC model
#' @export
#' @examples
#' train_DeepCC_model(tcga.fs, tcga.labels)
train_DeepCC_model <- function(trainData, trainLabels, epochs = 100, dropout = 0.5, activation_func = "selu") {

  ind <- !is.na(y_train)
  x_train <- trainData[ind, ]
  y_train <- factor(trainLabels)
  levels <- levels(y_train)
  y_train <- to_categorical(as.numeric(y_train[ind]) - 1, 4)

  k_clear_session()

  init_methods <- "glorot_uniform"
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 1024, activation = activation_func, input_shape = ncol(x_train), kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = dropout) %>%
    layer_dense(units = 1024, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = dropout) %>%
    layer_dense(units = 256, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = dropout) %>%
    layer_dense(units = 64, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = dropout) %>%
    layer_dense(units = 10, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_dense(units = 4, activation = 'softmax')


  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(),
    metrics = c('accuracy')
  )

  history <- model %>% fit(
    x_train, y_train,
    epochs = epochs, batch_size = 1024,
    #callbacks = callback_tensorboard("../logs/run_deepcc"),
    view_metrics = F,
    validation_split = validation_split
  )

  list(classifier=model, levels=levels)
}

#' @export
save_DeepCC_model <- function(deepcc_model, prefix) {
  deepcc_model$classifier %>% save_model_hdf5(filepath = paste0(prefix, ".hdf5"))
  levels <- deepcc_model$levels
  save(levels, file = paste0(prefix, ".RData"))
}

#' @export
load_DeepCC_model <- function(prefix) {
  load(file = paste0(prefix, ".RData"))
  classifer <- keras::load_model_hdf5(filepath = paste0(prefix, ".hdf5"))
  list(classifier = classifer, levels = levels)
}


#' Get DeepCC Labels
#'
#' This function classifys new data set using trained DeepCC model.
#'
#' @param DeepCCModel a trained DeepCC model
#' @param newData a data.frame containing functional spectra of new data (each row presents one sample)
#' @param cutoff a numeric indicating cutoff of poster probability
#' @return a character vector containing lables of training data
#' @export
#' @examples
#' get_DeepCC_label(deepcc.model, newdata.fs)
get_DeepCC_label <- function (DeepCCModel, newData, cutoff = 0.5, prob.mode = F, prob.raw = F)
{
  res <- keras::predict_proba(DeepCCModel$classifier, newData)
  predicted <- apply(res, 1, function(z) {
    if (max(z) >= cutoff) {
      which.max(z)
    }
    else {
      NA
    }
  })
  pred <- factor(predicted, levels = seq(length(DeepCCModel$levels)),
                 labels = DeepCCModel$levels)
  if (prob.mode) {
    pred <- data.frame(DeepCC = as.character(pred), Probability = round(apply(res,
                                                                              1, max), digits = 3))
  }

  if (prob.mode & prob.raw) {
    pred <- res
  }

  pred
}

#' Get DeepCC prob mat
#'
#' @export
#' @examples
#' getDeepCCLabels(deepcc.model, newdata.fs)
get_DeepCC_prob <- function(DeepCCModel, newData){

  res <- keras::predict_proba(DeepCCModel$classifier, newData)
  colnames(res) <- DeepCCModel$levels

  res
}

#' Get DeepCC Features
#'
#' This function obtains DeepCC Features from functional spectra.
#'
#' @param DeepCCModel a trained DeepCC model
#' @param fs a data.frame containing functional spectra (each row presents one sample)
#' @return a data.frame containing DeepCC Features extracted from the last hidden layer
#' @export
#' @examples
#' getDeepCCFeatures(deepcc.model, fs)
get_DeepCC_features <- function(DeepCCModel, fs){

  model <- DeepCCModel$classifier
  intermediate_layer_model <- keras::keras_model(model$input,
                                          model$layers[[length(model$layers)-1]]$output)
  predict(intermediate_layer_model, fs)
}
