#' Get current MSigDB gene sets
#'
#' Retrieves the current MSigDB gene set collection from the `msigdbr`
#' package. Returns a named list suitable for direct use in
#' `getFunctionalSpectra()`.
#'
#' @param collection character; MSigDB collection to retrieve. Use `NULL`
#'   (default) for all collections, or a specific one such as `"C2"`, `"C5"`,
#'   or `"H"`.
#' @return a named list of gene sets, each element a character vector of
#'   Entrez gene IDs
#' @export
#' @examples
#' \dontrun{
#' MSigDBr <- get_msigdbr()
#' fs <- getFunctionalSpectra(eps, geneSets = MSigDBr)
#' }
get_msigdbr <- function(collection = NULL){
  m_df <- msigdbr::msigdbr(species = "Homo sapiens")
  if (!is.null(collection)) {
    m_df <- m_df[m_df[["gs_collection"]] == collection, ]
  }
  m_df_2 <- as.data.frame(m_df)[, c("gs_name", "ncbi_gene")]
  gene_sets <- split(m_df_2[["ncbi_gene"]], m_df_2[["gs_name"]])
  gene_sets <- lapply(gene_sets, as.character)
  gene_sets
}
