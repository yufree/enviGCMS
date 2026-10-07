# Density weighted intensity for one sample

Density weighted intensity for one sample

## Usage

``` r
getdwtus(peak, n = 512, log = FALSE)
```

## Arguments

- peak:

  peaks intensity one sample

- n:

  the number of equally spaced points at which the density is to be
  estimated, default 512

- log:

  log transformation

## Value

Density weighted intensity for one sample

## Examples

``` r
data(list)
getdwtus(list$data[,1])
#> [1] 27.6034
```
