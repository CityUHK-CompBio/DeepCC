# Generate Functional Spectra

This function generates functional spectra for given gene expression
profiles. Uses a fast batch native kernel by default, with automatic
fallback to the legacy per-row parallel path when the input structure
prevents batch indexing.

## Usage

``` r
getFunctionalSpectra(eps, geneSets = "MSigDB", scale = TRUE, cores = NULL)
```

## Arguments

- eps:

  a data.frame containing gene expression profiles (each row presents
  one sample)

- geneSets:

  gene sets to score. Either \`"MSigDB"\` for the bundled cache,
  \`"MSigDB\_\<version\>"\` for an explicit bundled release, or a named
  list such as the result of \[get_msigdbr()\]

- scale:

  logical indicating whether to center each gene column (default: TRUE)

- cores:

  integer or NULL; number of CPU cores. NULL uses the native kernel's
  thread pool. Set to a specific integer for legacy parallel fallback.

## Value

a data.frame containing functional spectra

## See also

[`getFunctionalSpectrum`](https://cityuhk-compbio.github.io/DeepCC/reference/getFunctionalSpectrum.md)
for a single expression profile.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
eps <- as.data.frame(matrix(rnorm(10*100), nrow=10, ncol=100))
colnames(eps) <- paste0("G", seq_len(100))
fs <- getFunctionalSpectra(eps, geneSets=list(setA=c("G1","G5","G20")))
} # }
```
