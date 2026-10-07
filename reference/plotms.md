# plot GC/LC-MS data as a heatmap with TIC

plot GC/LC-MS data as a heatmap with TIC

## Usage

``` r
plotms(data, log = FALSE)
```

## Arguments

- data:

  imported data matrix of GC-MS

- log:

  transform the intensity into log based 10

## Value

heatmap

## Examples

``` r
if (FALSE) { # \dontrun{
png('test.png')
plotms(matrix)
dev.off()
} # }
```
