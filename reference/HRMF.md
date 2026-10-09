# High Resolution Mass Filtering (HRMF) for GC/LC HRMS data

Performs high-resolution mass filtering by matching experimental mass
spectral peaks against theoretical isotope patterns derived from one or
more candidate chemical formulae. Calculates forward HRMF score, reverse
(MSP) score, and Figure of Merit (FoM) for each candidate.

## Usage

``` r
HRMF(
  msp,
  formula,
  charge = 1,
  adduct = NULL,
  mass_accuracy = 5,
  intensity_cutoff = 1,
  IR_RelAb_cutoff = 1,
  detailed = FALSE
)
```

## Arguments

- msp:

  list. A single compound entry as returned by
  [`getMSP`](https://yufree.github.io/enviGCMS/reference/getMSP.md),
  containing at least a `spectra` element with columns `mz` and
  `intensity`.

- formula:

  character vector. One or more candidate chemical formulae to evaluate.

- charge:

  integer. Charge state: 1 for positive, -1 for negative, 0 for neutral.
  Default 1 (radical cation in EI). Ignored when `adduct` is supplied.

- adduct:

  character vector of adduct names for singly-charged ESI ions (e.g.
  `"[M+H]+"`, `c("[M+H]+", "[M+Na]+")`), or a named numeric vector of
  custom mass deltas in Da with names ending in "+" or "-" (e.g.
  `c("[M+MeOH]+" = 33.033491)`). Default NULL keeps the legacy
  radical-ion/neutral behaviour governed by `charge`. Built-in adducts:
  "\[M+H\]+", "\[M+Na\]+", "\[M+NH4\]+", "\[M+K\]+", "\[M+H-H2O\]+",
  "\[M-H\]-", "\[M+Cl\]-", "\[M+HCOO\]-", "\[M+CH3COO\]-".

- mass_accuracy:

  numeric. Mass accuracy in ppm. Default 5.

- intensity_cutoff:

  numeric. Minimum absolute intensity to retain a peak. Default 1.

- IR_RelAb_cutoff:

  numeric. Relative abundance cutoff (%) for theoretical isotopologues.
  Default 1.

- detailed:

  logical. If TRUE, return detailed list with all_ions, compound, and
  HRMF_scores for each formula. If FALSE (default), return a summary
  data.frame.

## Value

If `detailed = FALSE`, a data.frame with one row per candidate formula
(and per adduct, if supplied) and columns: Candidate, Adduct,
peak_count_forw, df_theortomsp, HRMF_theor_score, peak_count_rev,
df_msptotheor, HRMF_msp_score, FoM. If `detailed = TRUE`, a named list
of detailed results per formula/adduct.

## Details

The method is based on Kwiecien et al. (2015)
[doi:10.1021/acs.analchem.5b01503](https://doi.org/10.1021/acs.analchem.5b01503)
.

Unlike the original MSxplorer implementation, this version uses the
native HORIZON (Heavy-first Ordered Recursive Inference with Zero-loop
Optimal Navigation) engine via Rcpp for both high-performance formula
decomposition and isotope pattern calculation, requiring no external
dependencies.

The input `msp` should be a single entry from the list returned by
[`getMSP`](https://yufree.github.io/enviGCMS/reference/getMSP.md). For
batch processing of entire MSP files, see
[`getHRMF`](https://yufree.github.io/enviGCMS/reference/getHRMF.md).

For EI (GC-MS) radical ions use the default `charge` argument. For ESI
(LC-MS/MS) data with even-electron adduct ions such as \[M+H\]+ or
\[M-H\]-, supply `adduct` so that fragment m/z values are converted
internally; the adduct should describe how fragment ions are formed
(usually the same as the precursor adduct for singly-charged
fragmentation).

## See also

[`getHRMF`](https://yufree.github.io/enviGCMS/reference/getHRMF.md) for
batch processing,
[`getMSP`](https://yufree.github.io/enviGCMS/reference/getMSP.md) for
reading MSP files.

## Examples

``` r
if (FALSE) { # \dontrun{
# Read MSP file and run HRMF on the first compound (EI radical cation)
msp_data <- getMSP("spectrum.msp")
result <- HRMF(msp_data[[1]], formula = "C8H11NO")

# Compare multiple candidates
result <- HRMF(msp_data[[1]], formula = c("C8H11NO", "C7H9NO2"))

# LC-MS/MS with protonated fragments
result <- HRMF(msp_data[[1]], formula = "C8H10N4O2", adduct = "[M+H]+")

# Multiple candidate adducts are scored separately
result <- HRMF(msp_data[[1]], formula = "C8H10N4O2",
               adduct = c("[M+H]+", "[M+Na]+"))
} # }
```
