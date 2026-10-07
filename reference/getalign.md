# Align two peaks vectors by mass to charge ratio and/or retention time

Align two peaks vectors by mass to charge ratio and/or retention time

## Usage

``` r
getalign(mz1, mz2, rt1 = NULL, rt2 = NULL, ppm = 10, deltart = 10)
```

## Arguments

- mz1:

  the mass to charge of reference peaks

- mz2:

  the mass to charge of peaks to be aligned

- rt1:

  retention time of reference peaks

- rt2:

  retention time of peaks to be aligned

- ppm:

  mass accuracy, default 10

- deltart:

  retention time shift table, default 10 seconds

## Value

data frame with aligned peaks table

## Examples

``` r
mz1 <- c(221.1171, 227.1390, 229.1546, 233.1497, 271.0790 )
mz2 <- c(282.279, 281.113, 227.139, 227.139, 302.207)
rt1 <- c(590.8710, 251.3820, 102.9230, 85.8850, 313.8240)
rt2 <- c(787.08, 160.02, 251.76, 251.76, 220.26)
getalign(mz1,mz2,rt1,rt2)
#>   xid     mz1     mz2     rt1    rt2
#> 1   2 227.139 227.139 251.382 251.76
```
