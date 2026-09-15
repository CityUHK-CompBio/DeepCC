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

## Details

A functional spectrum summarises, for one sample, how strongly each gene
set is enriched among that sample's most highly expressed genes.

With \`scale = TRUE\` (the default) each gene is centred across the
samples in \`eps\` by subtracting its mean, so every value becomes a
deviation from the cohort average for that gene. Each sample is then
ranked independently on those centred values and scored with the
weighted running-sum enrichment statistic.

This is not a log fold change. DeepCC does not compare labelled groups,
does not select differentially expressed genes, and needs no group
labels. A score is relative to the other samples supplied in \`eps\`, so
the same sample receives different scores in a different cohort. Supply
the full cohort you wish to compare against, and keep that cohort fixed
between training and prediction.

With \`scale = FALSE\` no centring is applied and each sample is ranked
on its own values, which makes a score independent of the other rows.

A single row with \`scale = TRUE\` centres to zero and therefore returns
all zeros. Use \[getFunctionalSpectrum()\] for genuine single-sample
scoring.

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
