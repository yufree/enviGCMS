# Compute pooled QC linear index according to run order

Compute pooled QC linear index according to run order

## Usage

``` r
getpqsi(data, order, n = 5)
```

## Arguments

- data:

  peaks intensity list with row as peaks and column as samples

- order:

  run order of pooled QC samples

- n:

  samples numbers used for linear regression

## Value

vector for the peaks proportion with significant changes in linear
regression after FDR control.
