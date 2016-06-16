# DeepCC
DeepCC: a deep learning-based framework for cancer classification

## Dependence
Current version of DeepCC depends on H2O deep learning framework, implemented by Java thus you should install JRE first. Please following the instructions on [h2o.ai](http://www.h2o.ai/download/h2o/r).

## Installation
You can install DeepCC from GitHub directly using devtools.
```
install.packages("devtools")
devtools::install_github("gaofeng21cn/DeepCC")
```

## Quick Start
As a case study, you can obtain well organized colorectal cancer data from CRCSC's repository on [Synapse](https://www.synapse.org/#!Synapse:syn2623706/wiki/).

DeepCC only need two input for start.
- a data.frame containing gene expression profiles (each row represents one patient sample and the column names should be Entrez identifiers of genes)
- a character vector containing training labels (`NA` is allowed since DeepCC can ignore them automatically)

Now assume you have the gene expression profiles in `eps` and training labels in `labels`.
```
library(DeepCC)

# get functional spectra from gene expression profiles
# use parameter "cores" to indicate how many cpu cores you what to use, by defaut DeepCC will use all your cores - 1.
fs <- getFunctionalSpectra(eps)

# train DeepCC model
# use parameter "cores" to indicate how many cpu cores you what to use, by defaut DeepCC will use all your cores - 1.
deepcc.model <- trainDeepCCModel(fs, labels)

# obtain deep features
df <- getDeepCCFeatures(deepcc.model, fs)

# classify new data set used trained DeepCC model
# for a batch of samples
new.fs <- getFunctionalSpectra(new.eps)
pred.lables <- getDeepCCLabels(deepcc.model, new.fs)

# for single sample, you have to provide a reference expression profile.
# In DeepCC we prepared average expression profile of each cancer types in TCGA project as reference. To use them, just use the TCGA identifier to indicate the cancer type. Alternatively, you can provide your own reference of the platform and cancer type.
new.fs <- getFunctionalSpectrum(new.ep, refExp = "COADREAD")
pred.lable <- getDeepCCLabels(deepcc.model, new.fs)
```

## Additional tools
#### Benchmark
- `crossValidataion` performs cross validation.

#### Plotting
- `plotKMCurve` plots Kaplan-Meier curve.
