# Merge positive and negative mode data

Merge positive and negative mode data

## Usage

``` r
getpn(pos, neg, ppm = 5, pmd = 2.02, digits = 2, cutoff = 0.9)
```

## Arguments

- pos:

  a list with mzrt profile collected from positive mode. The sample
  order should match the negative mode.

- neg:

  a list with mzrt profile collected from negative mode.The sample order
  should match the positive mode.

- ppm:

  pmd mass accuracy, default 5

- pmd:

  numeric or numeric vector

- digits:

  mass or mass to charge ratio accuracy for pmd, default 2

- cutoff:

  correlation coefficients, default 0.9

## Value

mzrt object with group information from pos mode
