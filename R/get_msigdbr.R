#' get MSigDBr from R package `msigdbr`
#'
#' This function defines MSigDBr from a R package, with 25, 724 gene sets
#'
#' @param cores a integer indicating cpu cores used in parallel computing (default = all cores -2 )
#'
#' @return a list containing 25, 724 gene sets, each sets contains multiple entrez_gene
#' @examples
#' \dontrun{
#' MSigDBr <- get_msigdbr()
#' }

get_msigdbr <- function(cores = NULL){
  m_df <- msigdbr::msigdbr(species = "Homo sapiens")
  set_name <- unique(m_df[["gs_name"]])

  m_df_2 <- m_df[, c("gs_name", "entrez_gene")]

  get_list <- function(g_name){
    tmp <- as.character(m_df_2[m_df_2[["gs_name"]] == g_name, "entrez_gene"])
    tmp
  }

  gene_sets <- lapply(seq_along(set_name), function(idx) {
    get_list(set_name[idx])
  })
  names(gene_sets) <- set_name

  gene_sets
}
