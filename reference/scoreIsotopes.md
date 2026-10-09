# Calculate Isotopic Pattern Similarity Scores

Computes modern spectral similarity metrics including Cosine similarity,
intensity-weighted dot product, and log-likelihood between observed mass
spectral peaks and theoretical isotopic distributions.

## Usage

``` r
scoreIsotopes(
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

  Numeric vector of observed intensities.

- theo_mz:

  Numeric vector of theoretical isotopic m/z values.

- theo_int:

  Numeric vector of theoretical isotopic abundances.

- tolerance:

  Numeric mass tolerance for peak alignment (default: 0.005).

- ppm:

  Logical indicating whether `tolerance` is in ppm (default: FALSE).

## Value

A list containing:

- cosine:

  Cosine similarity (unweighted dot product of normalized spectra).

- weighted_cosine:

  Mass- and intensity-weighted dot product (MassBank/NIST style).

- log_likelihood:

  Multinomial log-likelihood of observed peak distribution.

- matched_peaks:

  Integer count of successfully matched peaks.

## Examples

``` r
obs_m <- c(180.0634, 181.0668)
obs_i <- c(100, 6.5)
theo_m <- c(180.0634, 181.0668)
theo_i <- c(1.0, 0.066)
scoreIsotopes(obs_m, obs_i, theo_m, theo_i)
#> $cosine
#> [1] 0.9999995
#> 
#> $weighted_cosine
#> [1] 0.9999983
#> 
#> $log_likelihood
#> [1] -0.2298068
#> 
#> $matched_peaks
#> [1] 2
#> 
```
