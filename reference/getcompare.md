# Align multiple peaks list to one peak list

Align multiple peaks list to one peak list

## Usage

``` r
getcompare(..., index = 1, ppm = 5, deltart = 5)
```

## Arguments

- ...:

  peaks list, mzrt objects

- index:

  numeric, the index of reference peaks.

- ppm:

  pmd mass accuracy, default 5

- deltart:

  retention time shift table, default 10 seconds

## Value

list object with aligned mzrt objects
