# Get the selected isotopologues at certain MS data

Get the selected isotopologues at certain MS data

## Usage

``` r
getisotopologues(formula = "C6H11O6", charge = 1, width = 0.3, cutoff = 0.05)
```

## Arguments

- formula:

  the molecular formula. 'C6H11O6' as default

- charge:

  the charge of that molecular. 1 in EI mode as default

- width:

  the width of the peak width on mass spectrum. 0.3 as default for low
  resolution mass spectrum.

- cutoff:

  numeric, minimum relative abundance (fraction of base peak) to
  consider, default 0.05 (5%).

## Examples

``` r
if (FALSE) { # \dontrun{
# show isotopologues
getisotopologues(formula = 'C6H11O6', charge = 1, width = 0.3)
} # }
```
