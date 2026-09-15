# Calculate Enrichment Score

This function calculates enrichment score of a gene list on a specific
gene set.

## Usage

``` r
calcEnrichmentScore(geneList, geneSet)
```

## Arguments

- geneList:

  a named vector containing the values of gene expression

- geneSet:

  a vector containing genes to represent a gene set

## Value

a numeric indicating enrichment score

## Examples

``` r
geneList <- setNames(rnorm(10), paste0("G", seq_len(10)))
geneSet <- c("G1", "G5")
calcEnrichmentScore(geneList, geneSet)
#> [1] 0.6925436
```
