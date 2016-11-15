#' Cross Validation of DeepCC Model
#'
#' This function performs cross validation of DeepCC Mode on the training data.
#'
#' @param fs a data.frame containing functional spectra of patients (each row presents one sample)
#' @param labels a character vector containing training lables
#' @param fold a integer indicating the fold number of cross validation (default: 5)
#' @return a numeric indicating error rate in a single run
#' @export
#' @examples
#' cross_validataion(tcga.fs, tcga.labels)
cross_validataion <- function(fs, labels, fold = 5) {
  fs <- fs[!is.na(labels), ]
  labels <- na.omit(labels)

  n <- length(labels)
  testidx <- sample(n, n*(1/fold))

  trainData.fs <- fs[-testidx, ]
  trainData.labels <- labels[-testidx]

  testData.fs <- fs[testidx, ]
  testData.labels <- labels[testidx]

  deepcc.model <- trainDeepCCModel(trainData.fs, trainData.labels)
  pred.lables <- getDeepCCLabels(deepcc.model, testData.fs, 0)

  error_rate <- 1-mean(as.character(testData.labels) == as.character(pred.lables))
  error_rate
}

#' Get Gene Sets
#'
#' This function extract a list of gene sets from gmt file.
#'
#' @param file filename of the gmt file
#' @return a list containing gene sets by EntrezID
#' @export
#' @examples
#' msigdbv51 <- get_gene_sets("msigdb.v5.1.entrez.gmt")

get_gene_sets <- function(file) {
  msig <- GSEABase::getGmt(file, geneIdType=EntrezIdentifier())
  GSEABase::geneIds(msig)
}
