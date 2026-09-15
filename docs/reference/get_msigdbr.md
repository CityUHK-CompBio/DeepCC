# get MSigDBr from R package \`msigdbr\`

This function defines MSigDBr from a R package, with 25, 724 gene sets

## Usage

``` r
get_msigdbr(cores = NULL)
```

## Arguments

- cores:

  a integer indicating cpu cores used in parallel computing (default =
  all cores -2 )

## Value

a list containing 25, 724 gene sets, each sets contains multiple
entrez_gene

## Examples

``` r
if (FALSE) { # \dontrun{
MSigDBr <- get_msigdbr()
} # }
```
