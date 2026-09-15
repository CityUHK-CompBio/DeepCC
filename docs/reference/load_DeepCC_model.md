# Load DeepCC Model

Loads a saved DeepCC model. Models trained under Keras 2 before the
keras3 migration are rebuilt from their recorded architecture, so the
published pre-trained models keep working. Both the new format (with
\`feature_names\` metadata) and the legacy 0.1.1 format are supported.

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

## See also

\[save_DeepCC_model()\] to write a model in the current format
