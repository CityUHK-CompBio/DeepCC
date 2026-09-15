# Cross Validation of DeepCC Model

This function performs cross validation of DeepCC Mode on the training
data.

## Usage

``` r
cross_validation(fs, labels, fold = 5)
```

## Arguments

- fs:

  a data.frame containing functional spectra of patients (each row
  presents one sample)

- labels:

  a character vector containing training lables

- fold:

  a integer indicating the fold number of cross validation (default: 5)

## Value

a numeric indicating error rate in a single run

## Examples

``` r
if (FALSE) { # \dontrun{
cross_validation(fs, labels)
} # }
```
