# plot density weighted intensity for multiple samples

plot density weighted intensity for multiple samples

## Usage

``` r
plotdwtus(list, n = 512, ...)
```

## Arguments

- list:

  list with data as peaks list, mz, rt and group information

- n:

  the number of equally spaced points at which the density is to be
  estimated, default 512

- ...:

  parameters for \`plot\` function

## Value

Density weighted intensity for multiple samples

## Examples

``` r
data(list)
plotdwtus(list)
```
