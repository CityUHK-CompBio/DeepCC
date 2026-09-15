# DeepCC

[![R-CMD-check](https://github.com/CityUHK-CompBio/DeepCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CityUHK-CompBio/DeepCC/actions/workflows/R-CMD-check.yaml)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![R 4.0+](https://img.shields.io/badge/R-%E2%89%A54.0-blue)](https://cran.r-project.org)

DeepCC is an R package for cancer molecular subtype classification. It
maps gene expression profiles to functional spectra using MSigDB gene
sets, then classifies subtypes with a deep neural network. Single-sample
prediction is supported through platform-specific reference profiles.

> Gao, F., Li, C., Wang, X. DeepCC: a deep learning-based framework for
> cancer classification. *Oncogenesis* 8, 7 (2019).
> [DOI: 10.1038/s41389-019-0157-8](https://www.nature.com/articles/s41389-019-0157-8)

## Installation

```r
install.packages("remotes")
remotes::install_github("CityUHK-CompBio/DeepCC")
```

Deep learning operations require [keras3](https://cran.r-project.org/package=keras3)
and a Python TensorFlow backend. If you plan to train or classify:

```r
install.packages("keras3")
keras3::install_keras()
```

Functional spectra computation and visualization work without any Python
runtime.

## Usage

### Batch functional spectra

```r
library(DeepCC)

# eps: data.frame or matrix (samples × genes), colnames are Entrez IDs
fs <- getFunctionalSpectra(eps)
```

DeepCC ships a bundled MSigDB snapshot (release 2026.1.Hs, 35,361 gene
sets) so functional spectra work offline and without any optional
dependency. To pull the newest release instead, install
[msigdbr](https://cran.r-project.org/package=msigdbr):

```r
# Bundled snapshot (default, offline)
fs <- getFunctionalSpectra(eps)
fs <- getFunctionalSpectra(eps, geneSets = "MSigDB_2026.1.Hs")

# Newest release from msigdbr, optionally filtered to one collection
MSigDBr <- get_msigdbr()
fs <- getFunctionalSpectra(eps, geneSets = MSigDBr)

hallmark <- get_msigdbr(collection = "H")
fs <- getFunctionalSpectra(eps, geneSets = hallmark)
```

Pinning `geneSets` to an explicit release keeps an analysis reproducible
even after MSigDB is updated.

### What the scores mean

A functional spectrum records how strongly each gene set is enriched among
one sample's most highly expressed genes. With the default `scale = TRUE`,
each gene is centred across the samples you supply, each sample is ranked on
those centred values, and the weighted running-sum enrichment statistic is
computed per sample.

This is not a log fold change and no group labels are involved. A score is
relative to the cohort in `eps`, so the same sample receives different scores
in a different cohort; keep the cohort fixed between training and prediction.
A single row with `scale = TRUE` centres to zero, so use
`getFunctionalSpectrum()` for genuine single-sample scoring.

### Train a model

```r
deepcc_model <- train_DeepCC_model(fs, labels)
```

`labels` is a character vector with one label per sample; `NA` entries
are excluded from training. The model records `feature_names` so that
new data with the same columns can be safely reordered at prediction
time.

### Classify new samples

```r
# Batch prediction
pred_labels <- get_DeepCC_label(deepcc_model, new_fs)
prob_matrix <- get_DeepCC_prob(deepcc_model, new_fs)

# Single sample using a TCGA reference profile
fs_single <- getFunctionalSpectrum(ep, refExp = "COADREAD")
pred_label <- get_DeepCC_label(deepcc_model, fs_single)
```

The `cutoff` argument controls label rejection: samples whose maximum
class probability falls below the cutoff receive `NA`.

### Extract deep features

```r
features <- get_DeepCC_features(deepcc_model, fs)
```

Returns the 10-dimensional penultimate layer output, useful for
downstream visualization or clustering.

## Gene sets

DeepCC reads Entrez gene IDs, and every source below returns that
identifier type.

| Source | Sets | Usage |
| --- | --- | --- |
| Bundled MSigDB 2026.1.Hs (offline) | 35,361 | `geneSets = "MSigDB"` (default) |
| Newest MSigDB release | varies | `get_msigdbr()` |
| MSigDB collection subset, e.g. Hallmark | 50 | `get_msigdbr(collection = "H")` |
| GMT file | varies | `get_gene_sets("path.gmt")` |
| Custom named list | varies | pass directly to `geneSets` |

The bundled snapshot is refreshed only when a new MSigDB milestone release
is adopted, so a given DeepCC version always resolves to the same gene sets.

## Pre-trained models

Colorectal cancer models trained on TCGA-COADREAD and CRCSC datasets are
available from the
[deepcc_model repository](https://github.com/zero19970/deepcc_model).
HDF5 files in that repository are Git LFS pointers; use `git lfs pull`
after cloning to obtain actual weights.

Those models were saved by Keras 2. `load_DeepCC_model()` rebuilds them from
their recorded architecture, so they load under the current keras3 without
conversion. Saving a loaded model writes it in the current format.

## Reference profiles

Single-sample classification requires a reference expression profile
from the same cancer type and platform. Built-in TCGA references
(COADREAD, BRCA, OV, etc.) are included. For cross-platform data, use
`inverseRescale = TRUE` when the input is microarray and the reference
is RNA-seq.

## Documentation

Full function reference is available at
[cityuhk-compbio.github.io/DeepCC](https://cityuhk-compbio.github.io/DeepCC/).

## Citation

```r
citation("DeepCC")
```

## License

Apache License 2.0 — see [LICENSE](LICENSE).
