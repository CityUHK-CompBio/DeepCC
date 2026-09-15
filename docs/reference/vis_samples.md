# Visualization of samples

This function visualize samples.

## Usage

``` r
vis_samples(data, labels, color, guide_fill = "legend")
```

## Arguments

- data:

  data.frame

- labels:

  groups

- color:

  color

- guide_fill:

  legend

## Value

a ggolot2 object

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
df <- as.data.frame(matrix(rnorm(30*10), nrow=30, ncol=10))
labels <- sample(c("A", "B", "C"), 30, replace=TRUE)
color <- c(A="red", B="blue", C="green")
sample_plot <- vis_samples(df, labels, color)
} # }
```
