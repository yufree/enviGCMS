# Get chemical formula for mass to charge ratio.

Get chemical formula for mass to charge ratio.

## Usage

``` r
getformula(
  mz,
  charge = 0,
  window = 0.001,
  elements = list(C = c(1, 50), H = c(1, 50), N = c(0, 50), O = c(0, 50), P = c(0, 1), S
    = c(0, 1))
)
```

## Arguments

- mz:

  a vector with mass to charge ratio

- charge:

  The charge value of the formula, default 0 for autodetect

- window:

  The window accuracy in the same units as mass

- elements:

  Elements list to take into account.

## Value

list with chemical formula
