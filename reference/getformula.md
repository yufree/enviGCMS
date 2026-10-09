# Get chemical formula for mass to charge ratio using the HORIZON algorithm

Rapidly decomposes accurate mass-to-charge ratios into candidate
chemical formulas using the native HORIZON (Heavy-first Ordered
Recursive Inference with Zero-loop Optimal Navigation) engine. Supports
chemical plausibility filtering via Double Bond Equivalents (DBE),
Nitrogen rule, and Kind & Fiehn Seven Golden Rules.

## Usage

``` r
getformula(
  mz,
  charge = 0,
  window = 0.001,
  elements = list(C = c(1, 50), H = c(1, 50), N = c(0, 50), O = c(0, 50), P = c(0, 1), S
    = c(0, 1)),
  golden_rules = TRUE,
  detailed = FALSE,
  resolution = 0,
  nthreads = 1
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

- golden_rules:

  logical, whether to filter by Fiehn Seven Golden Rules (default: TRUE)

- detailed:

  logical, if TRUE return a list of data.frames with formula, exact
  mass, error, DBE, H/C ratio, Senior rules, and validation details; if
  FALSE (default), return a list of character vectors of valid formulas

- resolution:

  numeric, mass spectrometer resolving power for isotope peak merging
  (default: 0)

- nthreads:

  integer, number of OpenMP parallel threads for batch processing
  (default: 1)

## Value

list of chemical formulas (or data.frames if detailed = TRUE)
