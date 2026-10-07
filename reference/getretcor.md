# Perform peaks list alignment and return features table

Perform peaks list alignment and return features table

## Usage

``` r
getretcor(list, ts = 1, ppm = 10, deltart = 5, FUN)
```

## Arguments

- list:

  each element should be a data.frame with mz, rt and intensity as m/z,
  retention time in seconds and intensity of certain peaks.

- ts:

  template sample index in the list, default 1

- ppm:

  mass accuracy, default 10

- deltart:

  retention time shift table, default 5 seconds

- FUN:

  function to deal with multiple aligned peaks from one sample

## Value

mzrt object without group information
