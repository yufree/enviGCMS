# Check Chemical Formula against Fiehn Seven Golden Rules

Evaluates whether a candidate chemical formula satisfies heuristic
chemical rules based on Senior's valence rules (Senior 1951),
Hydrogen-to-Carbon (H/C) ratios, heteroatom ratios (N/C, O/C, P/C, S/C,
Halogen/C), and Double Bond Equivalents (DBE) as described by Kind &
Fiehn (2007).

## Usage

``` r
checkGoldenRules(formula, z = 0, min_hc = 0.1, max_hc = 6, max_dbe = 40)
```

## Arguments

- formula:

  Character string representing molecular formula (e.g. `"C6H12O6"`).

- z:

  Integer charge state (default: 0).

- min_hc:

  Numeric minimum H/C ratio (default: 0.1).

- max_hc:

  Numeric maximum H/C ratio (default: 6.0).

- max_dbe:

  Numeric maximum DBE (default: 40.0).

## Value

A list containing booleans indicating whether the formula passed Senior
rules, H/C ratio, heteroatom ratios, DBE, nitrogen rule, and overall
validity.

## References

Kind, T., Fiehn, O. Seven Golden Rules for heuristic filtering of
molecular formulas obtained by accurate mass spectrometry. BMC
Bioinformatics 8, 105 (2007).
[doi:10.1186/1471-2105-8-105](https://doi.org/10.1186/1471-2105-8-105)

Senior, J. K. Partitions and their connection with the problems of
chemical structure. J. Chem. Phys. 19, 865-873 (1951).

## Examples

``` r
checkGoldenRules("C6H12O6")
#> $valid
#> [1] TRUE
#> 
#> $senior_rules
#> [1] TRUE
#> 
#> $hc_ratio
#> [1] 2
#> 
#> $hc_pass
#> [1] TRUE
#> 
#> $heteroatom_pass
#> [1] TRUE
#> 
#> $dbe
#> [1] 1
#> 
#> $dbe_pass
#> [1] TRUE
#> 
#> $nitrogen_rule
#> [1] TRUE
#> 
checkGoldenRules("C3H32N12O3") # chemically impossible formula
#> $valid
#> [1] FALSE
#> 
#> $senior_rules
#> [1] TRUE
#> 
#> $hc_ratio
#> [1] 10.66667
#> 
#> $hc_pass
#> [1] FALSE
#> 
#> $heteroatom_pass
#> [1] TRUE
#> 
#> $dbe
#> [1] -6
#> 
#> $dbe_pass
#> [1] FALSE
#> 
#> $nitrogen_rule
#> [1] TRUE
#> 
```
