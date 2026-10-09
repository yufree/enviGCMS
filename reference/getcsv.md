# Convert an list object to csv file

Convert an list object to csv file

## Usage

``` r
getcsv(list, name, mzdigit = 4, rtdigit = 1, type = "o", target = FALSE, ...)
```

## Arguments

- list:

  list with data as peaks list, mz, rt and group information

- name:

  result name for csv and/or eic file, default NULL

- mzdigit:

  m/z digits of row names of data frame, default 4

- rtdigit:

  retention time digits of row names of data frame, default 1

- type:

  csv format for further analysis, m means Metaboanalyst, a means
  xMSannotator, p means Mummichog(NA values are imputed by
  \`getimputation\`, and F test is used here to generate stats and p
  value), o means full information csv (for \`pmd\` package), default o.
  mapo could output all those format files.

- target:

  logical, preserve original rowname of data or not for target data,
  default FALSE.

- ...:

  other parameters for \`write.table\`

## Value

NULL, csv file

## References

Li, S.; Park, Y.; Duraisingham, S.; Strobel, F. H.; Khan, N.; Soltow, Q.
A.; Jones, D. P.; Pulendran, B. PLOS Computational Biology 2013, 9 (7),
e1003123. Xia, J., Sinelnikov, I.V., Han, B., Wishart, D.S., 2015.
MetaboAnalyst 3.0 - making metabolomics more meaningful. Nucl. Acids
Res. 43, W251-W257.

## Examples

``` r
if (FALSE) { # \dontrun{
data(list)
getcsv(list,name='demo')
} # }
```
