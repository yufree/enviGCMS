# plot GC/LC-MS data as scatter plot

plot GC/LC-MS data as scatter plot

## Usage

``` r
plotmz(data, inscf = 3.5, ...)
```

## Arguments

- data:

  imported data matrix of GC-MS

- inscf:

  Log intensity cutoff for peaks, default 3.5

- ...:

  parameters for \`plot\` function

## Value

scatter plot

## Examples

``` r
if (FALSE) { # \dontrun{
data(matrix)
png('test.png')
plotmz(matrix)
dev.off()
} # }
```
