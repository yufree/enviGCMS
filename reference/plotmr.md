# plot the scatter plot for peaks list with threshold

plot the scatter plot for peaks list with threshold

## Usage

``` r
plotmr(
  list,
  rt = NULL,
  ms = NULL,
  inscf = 5,
  rsdcf = 30,
  imputation = "l",
  ...
)
```

## Arguments

- list:

  list with data as peaks list, mz, rt and group information

- rt:

  vector range of the retention time

- ms:

  vector vector range of the m/z

- inscf:

  Log intensity cutoff for peaks across samples. If any peaks show a
  intensity higher than the cutoff in any samples, this peaks would not
  be filtered. default 5

- rsdcf:

  the rsd cutoff of all peaks in all group, default 30

- imputation:

  parameters for \`getimputation\` function method

- ...:

  parameters for \`plot\` function

## Value

data fit the cutoff

## Examples

``` r
data(list)
plotmr(list)
```
