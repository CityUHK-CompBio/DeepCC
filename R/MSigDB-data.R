#' Bundled MSigDB gene set cache
#'
#' A snapshot of MSigDB gene sets bundled with DeepCC so that functional
#' spectra can be computed without network access and without the optional
#' `msigdbr` package. The cache stores Entrez gene IDs, matching the
#' identifier convention of DeepCC expression input.
#'
#' The bundled snapshot is refreshed only when a new MSigDB milestone
#' release is adopted. Use [get_msigdbr()] to retrieve the current release
#' from the `msigdbr` package instead.
#'
#' @format A named list of gene sets. Each element is a character vector of
#'   Entrez gene IDs. Attributes record the MSigDB release and provenance:
#'   \describe{
#'     \item{db_version}{MSigDB release, e.g. `"2026.1.Hs"`.}
#'     \item{source}{Provenance of the snapshot, e.g. `"msigdbr"`.}
#'     \item{collections}{Named character vector mapping each gene set to its
#'       MSigDB collection.}
#'   }
#' @seealso [get_msigdbr()]
#' @name MSigDB
#' @docType data
#' @keywords datasets
NULL
