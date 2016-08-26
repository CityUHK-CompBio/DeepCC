# DeepCC
DeepCC: a deep learning-based framework for cancer classification

## Dependence
Current version of DeepCC depends on H2O deep learning framework, implemented by Java thus you should install JRE first. Please following the instructions on [h2o.ai](http://www.h2o.ai/download/h2o/r).

## Installation
You can install DeepCC from GitHub directly using devtools.
```
install.packages("devtools")
devtools::install_github("CityUHK-CompBio/DeepCC")
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
```

After training, now you can use your DeepCC model to classify new sample(s). DeepCC supports both batch sample and single sample classification. The input data should be in the same format as above gene expression profile(s).

```
# classify new data set used trained DeepCC model
# for a batch of samples
new.fs <- getFunctionalSpectra(new.eps)
pred.lables <- getDeepCCLabels(deepcc.model, new.fs)

# for a single sample, you have to provide a reference expression profile.
new.fs <- getFunctionalSpectrum(new.ep, refExp = "COADREAD")
pred.lable <- getDeepCCLabels(deepcc.model, new.fs)
```
Note: You can generate the reference expression profile from your previous data or public data, which is the same(similiar) cancer type and platform.

In DeepCC we also prepared average expression profiles of each cancer types in TCGA project as references. To use them, just use the TCGA identifier (COADREAD, BRCA, OV, etc.) to indicate the cancer type.

Note: if your single sample is microarray data, we strongly sugguest turn the parameter `inverseRescale` on, since TCGA is RNA-Seq, which has very small expression value for low expressed genes, compared with microarray.


## Additional Tools
- `crossValidataion` performs cross validation.
=======
## Additional tools
#### Benchmark
- `cross_validataion` performs cross validation.

#### Plotting
- `plotKMCurve` plots Kaplan-Meier curve.
