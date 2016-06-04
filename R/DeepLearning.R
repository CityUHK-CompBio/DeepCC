#' Train DeepCC Model
#'
#' This function trains DeepCC Model on the training data.
#'
#' @param trainData a data.frame containing functional spectra of training data (each row presents one sample)
#' @param trainLabels a character vector containing lables of training data
#' @param hiddenLayers a numeric vector containing hidden layers architecture
#' @param cores a integer indicating cpu cores used in parallel computing (defaut = all cores - 1)
#' @return a trained h2o deep learning model
#' @export
#' @examples
#' trainDeepCCModel(tcga.fs, tcga.labels)
trainDeepCCModel <- function(trainData, trainLabels, hiddenLayers=c(2000,500,120,30,10), cores = parallel::detectCores() - 1) {
    h2o::h2o.init(nthreads = cores)
    train_sv <- h2o::as.h2o(data.frame(trainData, class=trainLabels)[!is.na(trainLabels), ])
    classifier <- h2o::h2o.deeplearning(x = 1:length(MSigDB), y=length(MSigDB)+1,
                                   training_frame = train_sv,
                                   activation="Tanh",
                                   epochs = 100,
                                   hidden = hiddenLayers,
                                   ignore_const_cols=F
    )
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
