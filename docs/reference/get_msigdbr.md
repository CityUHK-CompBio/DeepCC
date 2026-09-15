# Get current MSigDB gene sets

Retrieves the current MSigDB gene set collection from the \`msigdbr\`
package. Returns a named list suitable for direct use in
\`getFunctionalSpectra()\`.

## Usage

``` r
get_msigdbr(collection = NULL)
```

## Arguments

- collection:

  character; MSigDB collection to retrieve. Use \`NULL\` (default) for
  all collections, or a specific one such as \`"C2"\`, \`"C5"\`, or
  \`"H"\`.

## Value

a named list of gene sets, each element a character vector of Entrez
gene IDs

## Examples

``` r
if (FALSE) { # \dontrun{
MSigDBr <- get_msigdbr()
fs <- getFunctionalSpectra(eps, geneSets = MSigDBr)
} # }
```
