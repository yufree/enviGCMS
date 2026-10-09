# Find PFCs based on mass defect analysis

Find PFCs based on mass defect analysis

## Usage

``` r
findpfc(list)
```

## Arguments

- list:

  list with data as peaks list, mz, rt and group information, retention
  time should be in seconds

## Value

list list with potential PFCs compounds index

## References

Liu, Y.; D'Agostino, L. A.; Qu, G.; Jiang, G.; Martin, J. W.
High-Resolution Mass Spectrometry (HRMS) Methods for Nontarget Discovery
and Characterization of Poly- and per-Fluoroalkyl Substances (PFASs) in
Environmental and Human Samples. TrAC Trends in Analytical Chemistry
2019, 121, 115420.

## Examples

``` r
data(list)
pfc <- findpfc(list)
```
