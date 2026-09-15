# Get DeepCC Features

This function obtains DeepCC Features from functional spectra using the
second-to-last layer of the classifier in inference mode.

## Usage

``` r
get_DeepCC_features(DeepCCModel, fs)
```

## Arguments

- DeepCCModel:

  a trained DeepCC model

- fs:

  a data.frame containing functional spectra (each row presents one
  sample)

## Value

a data.frame containing DeepCC Features extracted from the
second-to-last layer
