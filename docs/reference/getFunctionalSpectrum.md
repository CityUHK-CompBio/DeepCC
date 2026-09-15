# Generate Functional Spectrum

This function generates functional spectrum for a single gene expression
profile.

## Usage

``` r
getFunctionalSpectrum(
  expressionProfile,
  geneSets = "MSigDBv7",
  refExp = NULL,
  logChange = FALSE,
  inverseRescale = FALSE,
  filter = -3
)
```

## Arguments

- expressionProfile:

  a named numeric vector containing gene expression profile

- geneSets:

  a List containing gene sets (default: MSigDB v7)

- refExp:

  a character indicating cancer typer according to TCGA's indentifier,
  or a named vector reference expression

- logChange:

  a logical flag indicating whether the input data is already in log
  change form, e.g., for two color microarray, you should turn it on.
  (default: FALSE)

- inverseRescale:

  a logical flag indicating whether we rescale the reference to the
  scale of input data. If your single sample is microarray data and the
  reference is RNA-Seq, you should turn it on. (default: FALSE)

- filter:

  a numeric indicating the cutoff value of expression. (default: -3)

## Value

a numeric vector containing functional spectrum

## Note

You can generate the reference expression profile from your previous
data or public data, which is the same(similiar) cancer type and
platform. In DeepCC we also prepared average expression profiles of each
cancer types in TCGA project as references. To use them, just use the
TCGA identifier (COADREAD, BRCA, OV, etc.) to indicate the cancer type.
If your single sample is microarray data, we strongly sugguest turn the
parameter `inverseRescale` on, since TCGA is RNA-Seq, which has very
small expression value for low expressed genes, compared with
microarray.

## See also

[`getFunctionalSpectra`](https://cityuhk-compbio.github.io/DeepCC/reference/getFunctionalSpectra.md)
for a batch of gene expression profiles.

## Examples

``` r
if (FALSE) { # \dontrun{
ep <- setNames(rnorm(100), paste0("G", seq_len(100)))
ref <- setNames(rnorm(100), paste0("G", seq_len(100)))
fs <- getFunctionalSpectrum(ep, geneSets=list(setA=c("G1","G5","G20")), refExp=ref)
} # }
```
