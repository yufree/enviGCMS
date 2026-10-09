# Decompose a Vector of Masses in Parallel using the HORIZON Algorithm

Vectorized, multi-threaded chemical formula decomposition across
thousands of masses using the HORIZON (Heavy-first Ordered Recursive
Inference with Zero-loop Optimal Navigation) algorithm.

## Usage

``` r
rcpp_decompose_masses(
  masses,
  ppm = 2,
  mzabs = 1e-04,
  elements = NULL,
  minElements = "",
  maxElements = "",
  z = 0L,
  maxisotopes = 10L,
  golden_rules = TRUE,
  resolution = 0,
  fwhm = 0,
  nthreads = 1L
)
```

## Arguments

- masses:

  Numeric vector of target masses.

- ppm:

  Tolerance in ppm (default: 2.0).

- mzabs:

  Absolute mass deviation in Da (default: 0.0001).

- elements:

  Allowed chemical elements.

- minElements:

  Lower bounds on elements.

- maxElements:

  Upper bounds on elements.

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

- nthreads:

  Number of OpenMP threads to use (default: 1).

## Value

A list of decomposition results corresponding to each input mass.
