# Align mass to charge ratio and/or retention time to remove redundancy

Align mass to charge ratio and/or retention time to remove redundancy

## Usage

``` r
getalign2(mz, rt, ppm = 5, deltart = 5)
```

## Arguments

- mz:

  the mass to charge of reference peaks

- rt:

  retention time of reference peaks

- ppm:

  mass accuracy, default 10

- deltart:

  retention time shift table, default 10 seconds

## Value

index for

## Examples

``` r
mz <- c(221.1171, 221.1170, 229.1546, 233.1497, 271.0790 )
rt <- c(590.8710, 587.3820, 102.9230, 85.8850, 313.8240)
getalign2(mz,rt)
#> [1] 2 3 4 5
```
