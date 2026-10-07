# Get mass defect with certain scaled factor

Get mass defect with certain scaled factor

## Usage

``` r
getmassdefect(mass, sf)
```

## Arguments

- mass:

  vector of mass

- sf:

  scaled factors

## Value

dataframe with mass, scaled mass and scaled mass defect

## See also

[`plotkms`](https://yufree.github.io/enviGCMS/reference/plotkms.md)

## Examples

``` r
mass <- c(100.1022,245.2122,267.3144,400.1222,707.2294)
sf <- 0.9988
mf <- getmassdefect(mass,sf)
```
