# Relative Log Abundance (RLA) plots

Relative Log Abundance (RLA) plots

## Usage

``` r
plotrla(data, lv, type = "g", ...)
```

## Arguments

- data:

  data row as peaks and column as samples

- lv:

  factor vector for the group information

- type:

  'g' means group median based, other means all samples median based.

- ...:

  parameters for boxplot

## Value

Relative Log Abundance (RLA) plots

## Examples

``` r
data(list)
plotrla(list$data, as.factor(list$group$sample_group))
```
