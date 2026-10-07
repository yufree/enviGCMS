# Plot the heatmap of mzrt profiles

Plot the heatmap of mzrt profiles

## Usage

``` r
plothm(data, lv, index = NULL)
```

## Arguments

- data:

  data row as peaks and column as samples

- lv:

  group information

- index:

  index for selected peaks

## Examples

``` r
data(list)
plothm(list$data, lv = as.factor(list$group$sample_group))
```
