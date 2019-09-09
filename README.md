# DeepCC
DeepCC: a deep learning-based framework for cancer classification

### Dependencies
DeepCC dependes on MXNet framework, supporting both multiple CPUs and GPUs. Please install MXNet first according to their instruction.

### Installation
You can install DeepCC from GitHub directly using devtools.
```
install.packages("devtools")
devtools::install_github("CityUHK-CompBio/DeepCC")
```

### Quick start
As a case study, you can obtain well organized colorectal cancer data from CRCSC's repository on [Synapse](https://www.synapse.org/#!Synapse:syn2623706/wiki/).

DeepCC only need two input for start.
- a data.frame containing gene expression profiles (each row represents one patient sample and the column names should be Entrez identifiers of genes)
- a character vector containing training labels (`NA` is allowed since DeepCC can ignore them automatically)

Now assume you have the gene expression profiles in `eps` and training labels in `labels`.
```
library(DeepCC)

# get functional spectra from gene expression profiles
fs <- getFunctionalSpectra(eps)

# train DeepCC model
deepcc.model <- trainDeepCCModel(fs, labels)

# obtain deep features
df <- getDeepCCFeatures(deepcc.model, fs)
```

After training, now you can use your DeepCC model to classify new sample(s). DeepCC can classify samples in a data set, as well as individual samples. The input data should be in the same format as above gene expression profile(s).

```
# classify new data set used trained DeepCC model
# for a batch of samples
new.fs <- getFunctionalSpectra(new.eps)
pred.lables <- getDeepCCLabels(deepcc.model, new.fs)

# for a given single sample, you have to provide a reference expression profile.
new.fs <- getFunctionalSpectrum(new.ep, refExp = "COADREAD")
pred.lable <- getDeepCCLabels(deepcc.model, new.fs)
```
Note: You can generate customized reference expression profile from your previous data or public data, which is the same(similiar) cancer type and platform. Alternatively, you can use pre-defined reference in DeepCC by passing the cancer type (in the format of TCGA cancer types).

### Additional tools
- `cross_validataion` performs cross validation.
- `get_gene_sets` get gene sets list from MSigDB.
- `vis_sample` visualize samples


## Pre-defined data in DeepCC

### List of functional gene sets
By default, DeepCC will use MSigDB v5.0 (10, 348 gene sets) to generate functional spectra. You can also use MSigDB v6.0  (17, 779 gene sets).

### Pre-defined reference
In DeepCC we prepared average expression profiles of each cancer types in TCGA project as references. To use them, just use the TCGA identifier (COADREAD, BRCA, OV, etc.) to indicate the cancer type.

Note: if your single sample is microarray data, we strongly suggest you turn the parameter `inverseRescale` on, since TCGA is RNA-Seq data, which has very small expression value for low expressed genes, compared with microarray. `inverseRescale` can overcome this distribution different a little bit, but it's better to generate your customized reference.


