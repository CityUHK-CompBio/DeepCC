# Preprocess Gene List

This function preprocess gene list for futhur process

## Usage

``` r
preprocessGeneList(geneList)
```

## Arguments

- geneList:

  a named vector containing the value of gene expression

## Value

a named vecter containing the value of gene expression

## Examples

``` r
if (FALSE) { # \dontrun{
geneList <- setNames(rnorm(10), paste0("G", seq_len(10)))
preprocessGeneList(geneList)
} # }
```
