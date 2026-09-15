# Get current MSigDB gene sets

Retrieves the current MSigDB gene set collection from the \`msigdbr\`
package. Use this when you want the newest release; the bundled
\[MSigDB\] cache is used by default and needs neither \`msigdbr\` nor
network access.

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

## See also

\[MSigDB\] for the bundled offline snapshot

## Examples

``` r
if (FALSE) { # \dontrun{
MSigDBr <- get_msigdbr()
fs <- getFunctionalSpectra(eps, geneSets = MSigDBr)
} # }
```
