# Get the selected isotopologues at certain MS data

Get the selected isotopologues at certain MS data

## Usage

``` r
getisotopologues(formula = "C6H11O6", charge = 1, width = 0.3)
```

## Arguments

- formula:

  the molecular formula.

- charge:

  the charge of that molecular. 1 in EI mode as default

- width:

  the width of the peak width on mass spectrum. 0.3 as default for low
  resolution mass spectrum.

## Examples

``` r
if (FALSE) { # \dontrun{
# show isotopologues
getisotopologues(formula = 'C6H11O6', charge = 1, width = 0.3)
} # }
```
