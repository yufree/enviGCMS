# Batch High Resolution Mass Filtering for all compounds in an MSP file

Applies [`HRMF`](https://yufree.github.io/enviGCMS/reference/HRMF.md) to
every compound in an MSP file that has a chemical formula annotation.

## Usage

``` r
getHRMF(
  file,
  charge = 1,
  adduct = NULL,
  mass_accuracy = 5,
  intensity_cutoff = 1,
  IR_RelAb_cutoff = 1
)
```

## Arguments

- file:

  character. Path to an MSP file.

- charge:

  integer. Charge state (see
  [`HRMF`](https://yufree.github.io/enviGCMS/reference/HRMF.md)).
  Default 1. Ignored when `adduct` is supplied.

- adduct:

  character vector of adduct names or named numeric vector of custom
  mass deltas (see
  [`HRMF`](https://yufree.github.io/enviGCMS/reference/HRMF.md)).
  Default NULL.

- mass_accuracy:

  numeric. Mass accuracy in ppm. Default 5.

- intensity_cutoff:

  numeric. Minimum absolute intensity. Default 1.

- IR_RelAb_cutoff:

  numeric. Relative abundance cutoff (%) for theoretical isotopologues.
  Default 1.

## Value

A named list where each element contains the HRMF results for one
compound (with a valid formula) from the MSP file.

## See also

[`HRMF`](https://yufree.github.io/enviGCMS/reference/HRMF.md),
[`getMSP`](https://yufree.github.io/enviGCMS/reference/getMSP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
results <- getHRMF("library.msp")

# LC-MS/MS library with protonated ions
results <- getHRMF("library.msp", adduct = "[M+H]+")
} # }
```
