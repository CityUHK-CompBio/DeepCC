# Get DeepCC Labels

This function classifys new data set using trained DeepCC model.

## Usage

``` r
get_DeepCC_label(
  DeepCCModel,
  newData,
  cutoff = 0.5,
  prob_mode = FALSE,
  prob_raw = FALSE
)
```

## Arguments

- DeepCCModel:

  a trained DeepCC model

- newData:

  a data.frame containing functional spectra of new data (each presnets
  one sample)

- cutoff:

  a numeric indicating cutoff of poster probability

- prob_mode:

  a logical flag; if TRUE, return a data.frame with labels and
  probabilities

- prob_raw:

  a logical flag; if TRUE and prob_mode is TRUE, return the raw
  probability matrix

## Value

a character vector containing lables of training data
