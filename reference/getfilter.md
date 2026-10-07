# Filter the data based on row and column index

Filter the data based on row and column index

## Usage

``` r
getfilter(list, rowindex = TRUE, colindex = TRUE, name = NULL, type = "o", ...)
```

## Arguments

- list:

  list with data as peaks list, mz, rt and group information

- rowindex:

  logical, row index to keep

- colindex:

  logical, column index to keep

- name:

  file name for csv and/or eic file, default NULL

- type:

  csv format for further analysis, m means Metaboanalyst, a means
  xMSannotator, p means Mummichog(NA values are imputed by
  \`getimputation\`, and F test is used here to generate stats and p
  value), o means full information csv (for \`pmd\` package), default o.
  mapo could output all those format files.

- ...:

  other parameters for \`getcsv\`

## Value

list with remain peaks, and filtered peaks index

## See also

[`getimputation`](https://yufree.github.io/enviGCMS/reference/getimputation.md),
[`getcsv`](https://yufree.github.io/enviGCMS/reference/getcsv.md)

## Examples

``` r
data(list)
li <- getdoe(list)
lif <- getfilter(li,rowindex = li$rsdindex)
```
