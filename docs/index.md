# DeepCC

[![R CMD
check](https://github.com/CityUHK-CompBio/DeepCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CityUHK-CompBio/DeepCC/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://cityuhk-compbio.github.io/DeepCC/LICENSE)
[![R
4.0+](https://img.shields.io/badge/R-%E2%89%A54.0-blue)](https://cran.r-project.org)
[![keras3](https://img.shields.io/badge/keras3-%E2%89%A51.5-blue)](https://cran.r-project.org/package=keras3)

DeepCC is a deep learning-based framework for cancer molecular subtype
classification. It combines functional gene sets (MSigDB) as prior
knowledge with a deep neural network classifier.

> Gao, F., Li, C., Wang, X. DeepCC: a deep learning-based framework for
> cancer classification. *Oncogenesis* 8, 7 (2019). [DOI:
> 10.1038/s41389-019-0157-8](https://www.nature.com/articles/s41389-019-0157-8)

## 2026 modernization

This release modernizes DeepCC for current R environments while
preserving statistical semantics:

| Area | Improvement |
|----|----|
| Enrichment score | Sparse hit-position algorithm replaces full-scan; 39–63× faster with max difference \< 4×10⁻¹³ |
| Deep learning | Migrated from legacy `keras` to modern `keras3` (Keras 3) |
| R compatibility | Requires R ≥ 4.0; fixed namespace imports, parallel registration, and documentation |
| Model metadata | New models record `feature_names` for safe column reordering |
| Package quality | `R CMD check` passes with no errors or warnings |

## Installation

``` r

install.packages("remotes")
remotes::install_github("CityUHK-CompBio/DeepCC")
```

For deep learning operations, install `keras3` and its Python backend:

``` r

install.packages("keras3")
keras3::install_keras()
```

Only functional spectra computation and plotting require the R package
itself; training and prediction load `keras3` on demand.

## Quick start

### Batch functional spectra

``` r

library(DeepCC)

# eps: data.frame (samples × genes), colnames are Entrez IDs
# Use MSigDB v7 by default, or pass your own named list of gene sets
fs <- getFunctionalSpectra(eps, geneSets = "MSigDBv7")
```

### Train and classify

``` r

deepcc_model <- train_DeepCC_model(fs, labels)

# Batch prediction
pred_labels <- get_DeepCC_label(deepcc_model, new_fs)
probs <- get_DeepCC_prob(deepcc_model, new_fs)

# Single sample with TCGA reference
fs_single <- getFunctionalSpectrum(ep, refExp = "COADREAD")
pred_label <- get_DeepCC_label(deepcc_model, fs_single)
```

### Deep features

``` r

features <- get_DeepCC_features(deepcc_model, fs)
```

## Performance

On Apple Silicon (R 4.6.0, Apple clang 21.0.0):

| Workload | Legacy | Sparse kernel | Speedup | Max difference |
|----|----|----|----|----|
| 200 samples × 2,000 genes × 100 sets | 2.83 s | 0.07 s | 39× | 4×10⁻¹⁴ |
| 500 samples × 20,000 genes × 500 sets | 136 s (projected) | 2.15 s | 63× | 3×10⁻¹³ |
| 1 sample × 20,000 genes × 22,596 sets | 62 s (projected) | 4.0 s | 16× | — |

Numbers are from single-process synthetic benchmarks; multi-thread
scaling is flat because per-sample sorting dominates.

## Pre-trained models

CRC models from DeepCC_online are available at
[zero19970/deepcc_model](https://github.com/zero19970/deepcc_model).
Note: HDF5 files in that repository are Git LFS pointers; clone with
`git lfs pull` to obtain actual weights.

## Gene sets

Built-in legacy datasets: MSigDB v5 (10,348 sets), v6 (17,779 sets), and
v7 (22,596 sets). For the latest MSigDB, use
[`msigdbr`](https://cran.r-project.org/package=msigdbr):

``` r

MSigDBr <- get_msigdbr()
fs <- getFunctionalSpectra(eps, geneSets = MSigDBr)
```

## Documentation

Full function reference and tutorials are at
[cityuhk-compbio.github.io/DeepCC](https://cityuhk-compbio.github.io/DeepCC/).

## Citation

If you use DeepCC, please cite:

> Gao, F., Li, C., Wang, X. DeepCC: a deep learning-based framework for
> cancer classification. *Oncogenesis* 8, 7 (2019). [DOI:
> 10.1038/s41389-019-0157-8](https://www.nature.com/articles/s41389-019-0157-8)

## License

MIT — see [LICENSE](https://cityuhk-compbio.github.io/DeepCC/LICENSE).
