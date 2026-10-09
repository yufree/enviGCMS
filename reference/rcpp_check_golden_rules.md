# Check Formula against Fiehn Seven Golden Rules

Tests molecular formula against Senior's valences, H/C ratio, heteroatom
ratios, and DBE.

## Usage

``` r
rcpp_check_golden_rules(
  formula,
  z = 0L,
  min_hc = 0.1,
  max_hc = 6,
  max_dbe = 40
)
```

## Arguments

- formula:

  Character molecular formula string.

- z:

  Charge state (default: 0).

- min_hc:

  Minimum H/C ratio (default: 0.1).

- max_hc:

  Maximum H/C ratio (default: 6.0).

- max_dbe:

  Maximum DBE (default: 40.0).

## Value

A list with rule pass/fail booleans and summary details.
