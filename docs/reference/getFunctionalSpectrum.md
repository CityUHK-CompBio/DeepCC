# Generate Functional Spectrum

This function generates functional spectrum for a single gene expression
profile.

## Usage

``` r
getFunctionalSpectrum(
  expressionProfile,
  geneSets = "MSigDB",
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

  gene sets to score. Either \`"MSigDB"\` for the bundled cache,
  \`"MSigDB\_\<version\>"\` for an explicit bundled release, or a named
  list such as the result of \[get_msigdbr()\]

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

## Details

A single sample cannot be centred against a cohort, so a reference
expression profile is required unless the input is already in log-change
form. The reference is first restricted to genes expressed above
\`filter\`, then the sample and reference are related by a linear fit
and the score is computed on the difference between the fitted reference
and the sample. The resulting spectrum is a comparison of the sample
against that reference, not an enrichment of differentially expressed
genes.

Set \`inverseRescale = TRUE\` when the sample is microarray and the
reference is RNA-seq; this reverses the direction of the rescaling to
account for the different expression scales. Set \`logChange = TRUE\`
when the input is already a log-change vector, in which case no
reference is used.

Scores are not comparable across different references, because the
reference defines the baseline being compared against.

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
