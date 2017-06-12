#' Train DeepCC Model
#'
#' This function trains DeepCC Model on the training data.
#'
#' @param trainData a data.frame containing functional spectra of training data (each row presents one sample)
#' @param trainLabels a character vector containing lables of training data
#' @param hiddenLayers a numeric vector containing hidden layers architecture
#' @param gups a integer indicating how many GPUs used in parallel computing using GPUs (defaut = 0)
#' @return a trained DeepCC model
#' @export
#' @examples
#' trainDeepCCModel(tcga.fs, tcga.labels)
trainDeepCCModel <- function(trainData, trainLabels, hiddenLayers=c(2000,500,120,30,10), gpus = 0) {
  fs <- data.matrix(trainData[!is.na(trainLabels), ])
  labels <- as.factor(na.omit(trainLabels))
  levels <- levels(labels)

  if(gpus > 0) {
    device <- lapply(0:(gpus-1), mx.gpu)
    } else {
    device  <- mx.cpu()
    }

  classifier <- mxnet::mx.mlp(fs, as.numeric(labels)-1, array.layout="rowmajor",
                               hidden_node=hiddenLayers, out_node=length(levels),
                               out_activation="softmax", num.round=30, eval.metric=mx.metric.accuracy,
                               learning.rate=0.01,
                               momentum=0.9,
                               optimizer = "sgd",
                               # optimizer = "adadelta", # adadelta is also good, w/o pre-setting learning rate and momentum
                               activation = "tanh",
                               initializer = mx.init.Xavier(),
                               device = device)

  list(classifier=classifier, levels=levels)
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
#' getDeepCCLabels(deepcc.model, newdata.fs)
getDeepCCLabels <- function(DeepCCModel, newData, cutoff=0.5){


  res <- predict(DeepCCModel$classifier, newData, array.layout="rowmajor")

  predicted <- apply(res, 2, function(z) {
    if(max(z) >= cutoff ) {
      which.max(z)
    } else { NA }
  })
  factor(predicted, levels=seq(length(DeepCCModel$levels)), labels=DeepCCModel$levels)

  #predicted <- factor(apply(predict(DeepCCModel$classifier, newData, array.layout="rowmajor"), 2, which.max), levels=seq(length(DeepCCModel$levels)), labels=DeepCCModel$levels)
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
getDeepCCFeatures <- function(DeepCCModel, fs){

  ## modified from fdavidcl on github
  ## ref: https://github.com/dmlc/mxnet/issues/2785

  predictInternal <- function(model, X, ctx=NULL, layer, array.batch.size=128, array.layout="auto") {
    # initialization stuff I probably don't care about ---------------------------
    if (is.null(ctx)) ctx <- mxnet:::mx.ctx.default()
    if (is.array(X) || is.matrix(X)) {
      if (array.layout == "auto") {
        array.layout <- mxnet:::mx.model.select.layout.predict(X, model)
      }
      if (array.layout == "rowmajor") {
        X <- t(X)
      }
    }
    # end initialization ---------------------------------------------------------
    # iterator creation ----------------------------------------------------------
    ## X iterates through the batches of input data
    X <- mxnet:::mx.model.init.iter(X, NULL, batch.size=array.batch.size, is.train=FALSE)
    X$reset()
    if (!X$iter.next()) stop("Cannot predict on empty iterator")
    dlist = X$value()
    # end iterator creation ------------------------------------------------------
    # executor creation ----------------------------------------------------------
    ## mx.simple.bind defined in https://github.com/dmlc/mxnet/blob/master/R-package/R/executor.R#L5
    ## internally calls mx.symbol.bind, defined in https://github.com/dmlc/mxnet/blob/master/R-package/src/executor.cc#L191
    ### see also: https://github.com/dmlc/mxnet/blob/e7514fe1b3265aaf15870b124bb6ed0edd82fa76/R-package/demo/basic_executor.R
    # pexec <- mxnet:::mx.simple.bind(model$symbol, ctx=ctx, data=dim(dlist$data), grad.req="null")
    internals <- model$symbol$get.internals()
    pexec <- mxnet:::mx.simple.bind(internals[[4 * layer]], ctx=ctx, data=dim(dlist$data), grad.req="null")
    # end executor creation ------------------------------------------------------
    # set up arg arrays ----------------------------------------------------------
    internalArgParams <- model$arg.params[1:(2 * layer)]
    # print(names(internalArgParams))

    mxnet:::mx.exec.update.arg.arrays(pexec, internalArgParams, match.name=T)
    mxnet:::mx.exec.update.aux.arrays(pexec, model$aux.params, match.name=T)
    # end set up arg arrays ------------------------------------------------------
    # the rest is left untouched -------------------------------------------------
    packer <- mxnet:::mx.nd.arraypacker()
    X$reset()
    while (X$iter.next()) {
      dlist = X$value()
      mxnet:::mx.exec.update.arg.arrays(pexec, list(data=dlist$data), match.name=T)
      mxnet:::mx.exec.forward(pexec, is.train=FALSE)
      out.pred <- mxnet:::mx.nd.copyto(pexec$ref.outputs[[1]], mxnet:::mx.cpu())
      padded <- X$num.pad()
      oshape <- dim(out.pred)
      ndim <- length(oshape)
      packer$push(mxnet:::mx.nd.slice(out.pred, 0, oshape[[ndim]] - padded))
    }
    X$reset()
    df <- t(packer$get())
    rownames(df) <- colnames(X)
    return(df)
  }

  model <- DeepCCModel$classifier
  layer <- length(model$symbol$get.internals()$arguments)/2 - 2

  predictInternal(model, fs, layer = layer, array.layout = "rowmajor", array.batch.size = 1)
}
