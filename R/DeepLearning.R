#' Train DeepCC Model
#'
#' This function trains DeepCC Model on the training data.
#'
#' @param trainData a data.frame containing functional spectra of training data (each row presents one sample)
#' @param trainLabels a character vector containing lables of training data
#' @param hiddenLayers a numeric vector containing hidden layers architecture
#' @param cores a integer indicating cpu cores used in parallel computing or GPUs (defaut = all cores - 1)
#' @param device a string indicating using CPU or GPU (defaut = "cpu")
#' @return a trained h2o/mxnet deep learning model
#' @export
#' @examples
#' trainDeepCCModel(tcga.fs, tcga.labels)
trainDeepCCModel <- function(trainData, trainLabels, hiddenLayers=c(2000,500,120,30,10), cores = NULL, device = "cpu") {
  if(device == "gpu") {
    device.gpu <- lapply(0:(cores-1), function(i) {
      mx.gpu(i)
    })

    fs <- data.matrix(trainData[!is.na(labels), ])
    labels <- as.factor(na.omit(trainLabels))

    classifier <- mxnet::mx.mlp(fs, as.numeric(labels)-1, array.layout="rowmajor",
                    hidden_node=hiddenLayers, out_node=length(levels(labels)),
                    out_activation="softmax", num.round=100, eval.metric=mx.metric.accuracy,
                    learning.rate=0.1,
                    momentum=0.9,
                    initializer=mx.init.uniform(0.07),
                    device = device.gpu)
  }
    else {
    if(is.null(cores)) cores <- parallel::detectCores() - 1
    h2o::h2o.init(nthreads = cores)
    train_sv <- h2o::as.h2o(data.frame(trainData, class=trainLabels)[!is.na(trainLabels), ])

    classifier <- h2o::h2o.deeplearning(x = 1:ncol(trainData), y=ncol(trainData)+1,
                                        training_frame = train_sv,
                                        activation="Tanh",
                                        epochs = 100,
                                        hidden = hiddenLayers,
                                        ignore_const_cols=F
    )
  }

  classifier
}

#' Get DeepCC Labels
#'
#' This function classifys new data set using trained DeepCC model.
#'
#' @param DeepCCModel a trained h2o deep learning model
#' @param newData a data.frame containing functional spectra of new data (each row presents one sample)
#' @param cutoff a numeric indicating cutoff of poster probability
#' @return a character vector containing lables of training data
#' @export
#' @examples
#' getDeepCCLabels(deepcc.model, newdata.fs)
getDeepCCLabels <- function(DeepCCModel, newData, cutoff=0.5){
  res <- as.data.frame(h2o::h2o.predict(DeepCCModel, h2o::as.h2o(newData)))
  res <- data.frame(res[, -1])

  predicted <- apply(res,1,function(z) {
    if(max(z) >= cutoff ) {
      colnames(res)[which.max(z)]
    } else { NA }
  })
  predicted
}

#' Get DeepCC Features
#'
#' This function obtains DeepCC Features from functional spectra.
#'
#' @param DeepCCModel a trained h2o deep learning model
#' @param fs a data.frame containing functional spectra (each row presents one sample)
#' @return a data.frame containing DeepCC Features extracted from the last hidden layer
#' @export
#' @examples
#' getDeepCCFeatures(deepcc.model, fs)
getDeepCCFeatures <- function(DeepCCModel, fs){
  nLayers <- length(DeepCCModel@allparameters$hidden)

  features <- h2o::h2o.deepfeatures(DeepCCModel,  h2o::as.h2o(fs), layer = nLayers)
  as.data.frame(features)
}
