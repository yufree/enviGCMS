# plot scatter plot for rt-mz profile and output gif file for multiple groups

plot scatter plot for rt-mz profile and output gif file for multiple
groups

## Usage

``` r
gifmr(
  list,
  ms = c(100, 500),
  rsdcf = 30,
  inscf = 5,
  imputation = "l",
  name = "test",
  ...
)
```

## Arguments

- list:

  list with data as peaks list, mz, rt and group information

- ms:

  the mass range to plot the data

- rsdcf:

  the rsd cutoff of all peaks in all group

- inscf:

  Log intensity cutoff for peaks across samples. If any peaks show a
  intensity higher than the cutoff in any samples, this peaks would not
  be filtered. default 5

- imputation:

  parameters for \`getimputation\` function method

- name:

  file name for gif file, default test

- ...:

  parameters for \`plot\` function

## Value

gif file

## Examples

``` r
if (FALSE) { # \dontrun{
data(list)
gifmr(list)
} # }
```
