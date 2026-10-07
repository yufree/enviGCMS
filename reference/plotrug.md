# plot 1-d density for multiple samples

plot 1-d density for multiple samples

## Usage

``` r
plotrug(data, lv = NULL, indexx = NULL, indexy = NULL, ...)
```

## Arguments

- data:

  matrix

- lv:

  factor vector for the column

- indexx:

  index for matrix row

- indexy:

  index for matrix column

- ...:

  parameters for \`title\` function

## Examples

``` r
data(list)
plotrug(list$data)

plotrug(log(list$data), lv = as.factor(list$group$sample_group))
```
