# Write MSP file for NIST search

Write MSP file for NIST search

## Usage

``` r
writeMSP(list, name = "unknown", sep = FALSE)
```

## Arguments

- list:

  a list with spectra information

- name:

  name of the compounds

- sep:

  numeric or logical the numbers of spectra in each file and FALSE to
  include all of the spectra in one msp file

## Value

none a MSP file will be created.

## Examples

``` r
if (FALSE) { # \dontrun{
intensity <- c(10000,20000,10000,30000,5000)
mz <- c(101,143,189,221,234)
writeMSP(list(list(spectra = cbind.data.frame(mz,intensity))), name = 'test')
} # }
```
