#' Preprocess Gene List
#'
#' This function preprocess gene list for futhur process.
#'
#' @param geneList a named vector containing the value of gene expression
#' @return a named vector containing the value of gene expression
#' @examples
#' preprocessGeneList(geneList)
preprocessGeneList <- function(geneList) {
  geneList <- geneList[which((!is.na(geneList)) & (names(geneList)!="") & (!is.na(names(geneList))))]
  geneList <- tapply(t(geneList), names(geneList), max)
  geneList[order(geneList, decreasing=TRUE)]
}

#' Calculate Enrichment Score
#'
#' This function calculates enrichment score of a gene list on a specific gene set.
#'
#' @param geneList a named vector containing the values of gene expression
#' @param geneSet a vector containing genes to represent a gene set
#' @return a numeric indicating enrichment score
#' @export
#' @examples
#' calcEnrichmentScore(geneList, geneSet)
calcEnrichmentScore <- function (geneList, geneSet)
{
  hits <- (names(geneList) %in% geneSet)
  nh <- sum(hits)
  N <- length(geneList)
  ES <- 0
  runningES <- rep(0, N)
  if (nh > N) {
    stop("Gene Set is larger than Gene List")
  }
  else {
    if (nh) {
      tmp <- rep(0, N)
      NR = sum(abs(geneList[hits]))
      tmp[hits] <- abs(geneList[hits])/NR
      tmp[!hits] <- -1/(N - nh)
      runningES <- cumsum(tmp)

      ESmax <- max(runningES)
      ESmin <- min(runningES)
      ES <- ifelse(abs(ESmin) > abs(ESmax), ESmin, ESmax)
    }
  }
  ES
}

#' Generate Functional Spectra
#'
#' This function generates functional spectra for given gene expression profiles.
#'
#' @param eps a data.frame containing gene expression profiles (each row presents one sample)
#' @param geneSets a List containg gene sets (defaut = MSigDB v5.0)
#' @param cores a integer indicating cpu cores used in parallel computing (defaut = all cores - 1)
#' @return a data.frame containing functional spectra
#' @importFrom foreach foreach %dopar%
#' @export
#' @examples
#' getFunctionalSpectra(eps)
getFunctionalSpectra <- function(eps, geneSets = MSigDB, cores = parallel::detectCores() - 1) {
  eps <- scale(eps, scale = FALSE)

  doParallel::registerDoParallel(cores)
  res <- do.call(rbind, foreach(idx=1:nrow(eps)) %dopar% {
    geneList <- preprocessGeneList(eps[idx,])
    sapply(geneSets, function(x) calcEnrichmentScore(geneList, x))
  })
  rownames(res) <- rownames(eps)
  res
}


