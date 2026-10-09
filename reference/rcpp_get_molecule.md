# Calculate Exact Mass and Isotope Distribution using the HORIZON Engine

Computes exact monoisotopic mass, DBE, and theoretical isotopic pattern
using the HORIZON engine.

## Usage

``` r
rcpp_get_molecule(formula, z = 0L, maxisotopes = 10L, resolution = 0, fwhm = 0)
```

## Arguments

- formula:

  Character string of the molecular formula (e.g. "C6H12O6", "H2O").

- z:

  Integer charge state (default: 0).

- maxisotopes:

  Maximum number of isotopic peaks to return (default: 10).

- resolution:

  Mass resolution (m / FWHM) for merging fine isotope peaks (default:
  0).

- fwhm:

  Full width at half maximum in Da for isotope peak merging (default:
  0).

## Value

A list compatible with Rdisop::getMolecule.
