# plot the kendrick mass defect diagram

plot the kendrick mass defect diagram

## Usage

``` r
plotkms(data, cutoff = 1000)
```

## Arguments

- data:

  vector with the name m/z

- cutoff:

  remove the low intensity

## See also

[`getmassdefect`](https://yufree.github.io/enviGCMS/reference/getmassdefect.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mz <- c(10000,5000,20000,100,40000)
names(mz) <- c(100.1022,245.2122,267.3144,400.1222,707.2294)
plotkms(mz)
} # }
```
