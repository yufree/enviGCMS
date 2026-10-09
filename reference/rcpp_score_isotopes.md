# Calculate Isotopic Pattern Similarity Scores

Computes cosine similarity, weighted dot product, and log-likelihood
between experimental peaks and theoretical isotopic distribution.

## Usage

``` r
rcpp_score_isotopes(
  obs_mz,
  obs_int,
  theo_mz,
  theo_int,
  tolerance = 0.005,
  ppm = FALSE
)
```

## Arguments

- obs_mz:

  Numeric vector of observed m/z values.

- obs_int:

  Numeric vector of observed peak intensities.

- theo_mz:

  Numeric vector of theoretical isotopic m/z values.

- theo_int:

  Numeric vector of theoretical isotopic abundances.

- tolerance:

  Mass matching tolerance (default: 0.005 Da or 10 ppm).

- ppm:

  Logical indicating if tolerance is in ppm (default: FALSE).

## Value

A list containing cosine score, weighted cosine score, log-likelihood,
and matching details.
