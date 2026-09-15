# Get Gene Sets

This function extract a list of gene sets from gmt file.

## Usage

``` r
get_gene_sets(file)
```

## Arguments

- file:

  filename of the gmt file

## Value

a list containing gene sets by EntrezID

## Examples

``` r
if (FALSE) { # \dontrun{
msigdbv51 <- get_gene_sets("msigdb.v5.1.entrez.gmt")
} # }
```
