# Load DeepCC Model

Loads a saved DeepCC model. Supports both the new format (with metadata)
and the legacy 0.1.1 format (without metadata).

## Usage

``` r
load_DeepCC_model(prefix)
```

## Arguments

- prefix:

  file path prefix

## Value

a DeepCC model with `classifier`, `levels`, and optionally
`feature_names`
