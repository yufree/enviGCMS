# plot intensity of peaks across samples or samples across peaks

plot intensity of peaks across samples or samples across peaks

## Usage

``` r
plotpeak(data, lv = NULL, indexx = NULL, indexy = NULL, ...)
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

## Value

parallel coordinates plot

## Examples

``` r
data(list)
# selected peaks across samples
plotpeak(t(list$data), lv = as.factor(c(rep(1,5),rep(2,nrow(list$data)-5))),1:10,1:10)

# selected samples across peaks
plotpeak(list$data, lv = as.factor(list$group$sample_group),1:10,1:10)
```
