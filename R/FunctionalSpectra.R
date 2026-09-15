#' Preprocess Gene List
#'
#' This function preprocess gene list for futhur process
#'
#' @param geneList a named vector containing the value of gene expression
#' @return a named vecter containing the value of gene expression
#' @examples
#' \dontrun{
#' geneList <- setNames(rnorm(10), paste0("G", seq_len(10)))
#' preprocessGeneList(geneList)
#' }
preprocessGeneList <- function(geneList) {
  geneList <- geneList[which((!is.na(geneList)) & (names(geneList)!="") & (!is.na(names(geneList))))]
  geneList <- tapply(t(geneList), names(geneList), max)
  geneList[order(geneList, decreasing = TRUE)]
}

#' Calculate Enrichment Score
#'
#' This function calculates enrichment score of a gene list on a specific gene set.
#'
#' @param geneList a named vector containing the values of gene expression
#' @param geneSet a vector containing genes to represent a gene set
#' @return a numeric indicating enrichment score
#' @useDynLib DeepCC
#' @import Rcpp
#' @export
#' @examples
#' geneList <- setNames(rnorm(10), paste0("G", seq_len(10)))
#' geneSet <- c("G1", "G5")
#' calcEnrichmentScore(geneList, geneSet)
calcEnrichmentScore <- function(geneList, geneSet)
{
  calcEnrichmentScoreCPP((names(geneList) %in% geneSet), geneList, 1)
}

#' Resolve gene sets from legacy string or list
#' @noRd
resolveGeneSets <- function(geneSets) {
  if (is.character(geneSets) && length(geneSets) == 1) {
    stop(paste("String geneSets shortcuts (MSigDBv5/v6/v7) have been removed.",
               "Use get_msigdbr() for current MSigDB data, get_gene_sets() for GMT files,",
               "or pass a named list directly."))
  }
  if (!is.list(geneSets)) stop("geneSets must be a character string or a named list.")
  geneSets
}

#' Map gene sets to flat integer indices for batch C++ kernel
#' @noRd
buildFlatGeneSetIndex <- function(geneSets, geneNames) {
  P <- length(geneSets)
  starts <- integer(P + 1)
  members_list <- lapply(seq_len(P), function(i) {
    genes <- geneSets[[i]]
    if (is.null(genes) || length(genes) == 0) return(integer(0))
    idx <- match(unique(as.character(genes)), geneNames) - 1L
    idx[!is.na(idx)]
  })
  starts[-1] <- cumsum(vapply(members_list, length, integer(1)))
  members <- unlist(members_list, use.names = FALSE)
  list(starts = starts, members = members)
}

#' Generate Functional Spectra
#'
#' This function generates functional spectra for given gene expression profiles.
#' Uses a fast batch native kernel by default, with automatic fallback to the
#' legacy per-row parallel path when the input structure prevents batch indexing.
#'
#' @param eps a data.frame containing gene expression profiles (each row presents one sample)
#' @param geneSets a List containing gene sets (default: MSigDB v7)
#' @param scale logical indicating whether to center each gene column (default: TRUE)
#' @param cores integer or NULL; number of CPU cores. NULL uses the native kernel's
#'   thread pool. Set to a specific integer for legacy parallel fallback.
#' @return a data.frame containing functional spectra
#' @seealso  \code{\link{getFunctionalSpectrum}} for a single expression profile.
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @export
#' @examples
#' \dontrun{
#' set.seed(42)
#' eps <- as.data.frame(matrix(rnorm(10*100), nrow=10, ncol=100))
#' colnames(eps) <- paste0("G", seq_len(100))
#' fs <- getFunctionalSpectra(eps, geneSets=list(setA=c("G1","G5","G20")))
#' }
getFunctionalSpectra <- function(eps, geneSets = 'MSigDBv7', scale = TRUE, cores = NULL) {
  geneSets <- resolveGeneSets(geneSets)
  P <- length(geneSets)

  if (scale) eps <- scale(eps, scale = FALSE)
  eps <- as.matrix(eps)
  if (nrow(eps) == 0 || ncol(eps) == 0) stop("eps has zero rows or columns.")

  geneNames <- colnames(eps)
  if (is.null(geneNames)) geneNames <- as.character(seq_len(ncol(eps)))

  # Handle duplicate gene names by keeping max value per name
  dupNames <- duplicated(geneNames) | duplicated(geneNames, fromLast = TRUE)
  if (any(dupNames)) {
    # Reduce to unique names with max value per name (per sample)
    uniqNames <- unique(geneNames)
    epsRed <- matrix(NA_real_, nrow = nrow(eps), ncol = length(uniqNames),
                     dimnames = list(rownames(eps), uniqNames))
    for (gn in uniqNames) {
      cols <- which(geneNames == gn)
      epsRed[, gn] <- apply(eps[, cols, drop = FALSE], 1, max, na.rm = TRUE)
    }
    eps <- epsRed
    geneNames <- uniqNames
  }

  # Try the batch sparse kernel
  useBatch <- TRUE
  flat <- tryCatch({
    buildFlatGeneSetIndex(geneSets, geneNames)
  }, error = function(e) {
    useBatch <<- FALSE
    list(starts = integer(0), members = integer(0))
  })

  if (useBatch && length(flat$starts) == P + 1) {
    # Use the optimized batch kernel
    if (is.null(cores)) {
      res <- calcEnrichmentScoreBatchCPP(eps, flat$starts, flat$members)
    } else {
      res <- calcEnrichmentScoreBatchCPP(eps, flat$starts, flat$members, nthreads = as.integer(cores))
    }
    dimnames(res) <- list(rownames(eps), names(geneSets))
    return(as.data.frame(res))
  }

  # Legacy per-row parallel fallback
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 2L)
  doParallel::registerDoParallel(cores)
  on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  res <- foreach(idx = seq_len(nrow(eps)), .combine = rbind) %dopar% {
    geneList <- preprocessGeneList(eps[idx, ])
    sapply(geneSets, function(x) calcEnrichmentScore(geneList, x))
  }
  rownames(res) <- rownames(eps)
  colnames(res) <- names(geneSets)
  as.data.frame(res)
}

#' Generate Functional Spectrum
#'
#' This function generates functional spectrum for a single gene expression profile.
#'
#' @param expressionProfile a named numeric vector containing gene expression profile
#' @param geneSets a List containing gene sets (default: MSigDB v7)
#' @param refExp a character indicating cancer typer according to TCGA's indentifier, or a named vector reference expression
#' @param logChange a logical flag indicating whether the input data is already in log change form, e.g., for two color microarray, you should turn it on. (default: FALSE)
#' @param inverseRescale a logical flag indicating whether we rescale the reference to the scale of input data. If your single sample is microarray data and the reference is RNA-Seq, you should turn it on. (default: FALSE)
#' @param filter a numeric indicating the cutoff value of expression. (default: -3)
#' @return a numeric vector containing functional spectrum
#' @note You can generate the reference expression profile from your previous data or public data, which is the same(similiar) cancer type and platform.
#' In DeepCC we also prepared average expression profiles of each cancer types in TCGA project as references. To use them, just use the TCGA identifier (COADREAD, BRCA, OV, etc.) to indicate the cancer type.
#' If your single sample is microarray data, we strongly sugguest turn the parameter \code{inverseRescale} on, since TCGA is RNA-Seq, which has very small expression value for low expressed genes, compared with microarray.
#' @seealso \code{\link{getFunctionalSpectra}} for a batch of gene expression profiles.
#' @export
#' @examples
#' \dontrun{
#' ep <- setNames(rnorm(100), paste0("G", seq_len(100)))
#' ref <- setNames(rnorm(100), paste0("G", seq_len(100)))
#' fs <- getFunctionalSpectrum(ep, geneSets=list(setA=c("G1","G5","G20")), refExp=ref)
#' }
getFunctionalSpectrum <- function(expressionProfile, geneSets = 'MSigDBv7', refExp = NULL, logChange = FALSE, inverseRescale = FALSE, filter = -3) {
  expressionProfile <- unlist(expressionProfile)
  if(!logChange) {
    if(is.null(refExp)) stop("Must have a reference expression profile!")
    if(is.character(refExp)) {
      if(!(refExp %in% rownames(ref_eps))) stop(paste(refExp, "is not a support identifier of cancer types!\n Please use one in :", paste(row.names(ref_eps), collapse = " ")))
      refExp <- ref_eps[refExp, ]
    }
    # filter low expressed genes
    refExp <- refExp[refExp > filter]

    common <- intersect(names(expressionProfile), names(refExp))
    if(!inverseRescale) {
      expressionProfile <- stats::predict(stats::lm(refExp[common] ~ expressionProfile[common])) - expressionProfile[common]
    } else {
      expressionProfile <- expressionProfile[common] - stats::predict(stats::lm(expressionProfile[common] ~ refExp[common]))
    }
  }
  geneList <- preprocessGeneList(expressionProfile)

  geneSets <- resolveGeneSets(geneSets)

  res <- vapply(seq_len(length(geneSets)), function(i) calcEnrichmentScore(geneList, geneSets[[i]]), numeric(1))
  names(res) <- names(geneSets)
  res
}
