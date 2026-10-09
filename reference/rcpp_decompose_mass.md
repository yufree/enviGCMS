# Decompose Mass into Candidate Chemical Formulas using the HORIZON Algorithm

Finds all combinations of elements matching a mass within given error
bounds using the HORIZON (Heavy-first Ordered Recursive Inference with
Zero-loop Optimal Navigation) algorithm.

## Usage

``` r
rcpp_decompose_mass(
  mass,
  ppm = 2,
  mzabs = 1e-04,
  elements = NULL,
  minElements = "",
  maxElements = "",
  z = 0L,
  maxisotopes = 10L,
  golden_rules = TRUE,
  resolution = 0,
  fwhm = 0
)
```

## Arguments

- mass:

  Target mass or m/z value.

- ppm:

  Tolerance in ppm (default: 2.0).

- mzabs:

  Absolute mass deviation in Da (default: 0.0001).

- elements:

  Allowed chemical elements (string or vector).

- minElements:

  Lower bounds on elements (e.g. "C0H0", "C1").

- maxElements:

  Upper bounds on elements (e.g. "C50H50").

- z:

  Charge state (default: 0).

- maxisotopes:

  Maximum number of isotopes to compute (default: 10).

- golden_rules:

  Whether to apply Fiehn Seven Golden Rules (default: TRUE).

- resolution:

  Mass resolution for isotope merging (default: 0).

- fwhm:

  Full width at half maximum in Da for isotope peak merging (default:
  0).

## Value

A list compatible with Rdisop::decomposeMass.
