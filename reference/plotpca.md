# plot the PCA for multiple samples

plot the PCA for multiple samples

## Usage

``` r
plotpca(
  data,
  lv = NULL,
  index = NULL,
  center = TRUE,
  scale = TRUE,
  xrange = NULL,
  yrange = NULL,
  pch = NULL,
  ...
)
```

## Arguments

- data:

  data row as peaks and column as samples

- lv:

  group information

- index:

  index for selected peaks

- center:

  parameters for PCA

- scale:

  parameters for scale

- xrange:

  x axis range for return samples, default NULL

- yrange:

  y axis range for return samples, default NULL

- pch:

  default pch would be the first character of group information or
  samples name

- ...:

  other parameters for \`plot\` function

## Value

if xrange and yrange are not NULL, return file name of all selected
samples on 2D score plot

## Examples

``` r
data(list)
plotpca(list$data, lv = as.character(list$group$sample_group))
```
