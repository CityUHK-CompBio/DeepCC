# Train DeepCC Model

This function trains DeepCC Model on the training data using the modern
`keras3` interface. The network architecture, training recipe, and
legacy model format are preserved from DeepCC 0.1.1.

## Usage

``` r
train_DeepCC_model(
  trainData,
  trainLabels,
  epochs = 100,
  dropout = 0.4,
  activation_func = "selu",
  validation_split = 0.2
)
```

## Arguments

- trainData:

  a data.frame containing functional spectra of training data (each row
  presents one sample)

- trainLabels:

  a character vector containing lables of training data

- epochs:

  the number of epochs

- dropout:

  dropout rate

- activation_func:

  activation funtion

- validation_split:

  fraction of training data to use for validation

## Value

a trained DeepCC model with `classifier`, `levels`, and `feature_names`
fields

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
eps <- as.data.frame(matrix(rnorm(20*50), nrow=20, ncol=50))
colnames(eps) <- paste0("F", seq_len(50))
labels <- sample(c("A", "B", "C"), 20, replace=TRUE)
deepcc_model <- train_DeepCC_model(eps, labels, epochs=2)
} # }
```
