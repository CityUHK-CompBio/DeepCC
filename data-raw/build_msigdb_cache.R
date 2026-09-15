# Build the bundled MSigDB milestone cache.
#
# Run manually when a new MSigDB milestone release is adopted:
#   Rscript data-raw/build_msigdb_cache.R
#
# The cache lets DeepCC compute functional spectra without network access and
# without the msigdbr package installed. It stores Entrez IDs because DeepCC
# expression input is keyed by Entrez identifiers.

library(msigdbr)

build_msigdb_cache <- function(output = "data/MSigDB.rda") {
  m <- msigdbr::msigdbr(species = "Homo sapiens")
  df <- as.data.frame(m)[, c("gs_name", "ncbi_gene", "gs_collection")]
  df$ncbi_gene <- as.character(df$ncbi_gene)

  sets <- split(df$ncbi_gene, df$gs_name)
  collections <- vapply(split(df$gs_collection, df$gs_name), function(x) x[[1L]], character(1))

  MSigDB <- sets
  attr(MSigDB, "db_version") <- unique(as.character(m$db_version))
  attr(MSigDB, "source") <- "msigdbr"
  attr(MSigDB, "collections") <- collections

  save(MSigDB, file = output, compress = "xz", version = 2)
  invisible(output)
}

if (identical(environment(), globalenv()) && !interactive()) {
  out <- build_msigdb_cache()
  cat("wrote", out, "\n")
  cat("sets:", length(get(load(out))), "\n")
}
