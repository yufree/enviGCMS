# plot the density for multiple samples

plot the density for multiple samples

## Usage

``` r
plotden(data, lv, index = NULL, name = NULL, lwd = 1, ...)
```

## Arguments

- data:

  data row as peaks and column as samples

- lv:

  group information

- index:

  index for selected peaks

- name:

  name on the figure for samples

- lwd:

  the line width for density plot, default 1

- ...:

  parameters for \`plot\` function

## Examples

``` r
data(list)
plotden(list$data, lv = as.character(list$group$sample_group),ylim = c(0,1))
```
