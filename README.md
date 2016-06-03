# DeepCC
DeepCC: a deep learning-based framework for cancer classification

## Dependence
Current version of DeepCC deepends H2O deep learning framework, implemented by Java thus you should install JRE first. Please following the instructions on [h2o.ai](http://www.h2o.ai/download/h2o/r).

## Installation
You can install DeepCC from GitHub directly using devtools.
```
install.packages("devtools")
devtools::install_github("hadley/devtools")
```

## Quick Start
Assume you have the gene expression profiles in `tcga.eps` (each row represents one patient sample) and training labels in `tcga.labels`.
```
tcga.fs <- getFunctionalSpectra(tcga.eps)

deepcc.model <- trainDeepCCModel(tcga.fs, tcga.labels)
tcga.pred.lables <- getDeepCCLabels(deepcc.model, tcga.fs)
```
