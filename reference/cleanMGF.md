# Clean an MGF file by keeping only formula-explainable peaks

For each MS2 spectrum, the precursor neutral mass is derived from
`PEPMASS`, `CHARGE` and the supplied `adduct`, candidate parent formulae
are enumerated with Rdisop, and fragment peaks that can be explained as
sub-formulae of the best candidate (the one explaining the most peaks)
are kept. Spectra without `PEPMASS` or without any candidate formula are
written back unchanged.

## Usage

``` r
cleanMGF(
  file,
  out_file = NULL,
  adduct = "[M+H]+",
  charge = NULL,
  mass_accuracy = 5,
  elements = "CHNOPS",
  max_candidates = 3,
  intensity_cutoff = 0,
  min_peaks = 2
)
```

## Arguments

- file:

  character. Path to an MGF file.

- out_file:

  character. Output MGF path. Default NULL writes
  `<file stem>_clean.mgf` next to the input file.

- adduct:

  character. Single adduct name (see
  [`HRMF`](https://yufree.github.io/enviGCMS/reference/HRMF.md)), e.g.
  "\[M+H\]+". Fragments are assumed singly-charged with this adduct.

- charge:

  integer. Precursor charge state used when the `CHARGE` header is
  missing. Default NULL uses 1 (from the adduct sign).

- mass_accuracy:

  numeric. Mass accuracy in ppm for both precursor decomposition and
  fragment matching. Default 5.

- elements:

  character. Elements allowed in precursor formula enumeration, e.g.
  "CHNOPS". Default "CHNOPS".

- max_candidates:

  integer. Maximum number of precursor formula candidates to evaluate
  per spectrum. Default 3.

- intensity_cutoff:

  numeric. Discard peaks below this absolute intensity before
  annotation. Default 0 (keep all).

- min_peaks:

  integer. Drop spectra left with fewer than this many peaks after
  cleaning. Default 2.

## Value

Invisibly, a data.frame with one row per spectrum: title, pepmass,
charge, chosen formula, peak counts before/after, kept fraction and
status ("cleaned", "kept" = unchanged, or "dropped").

## Details

The cleaned file preserves the original header lines and peak line
formatting; the chosen parent formula is appended as a `FORMULA=` header
when a candidate is found.

## See also

[`HRMF`](https://yufree.github.io/enviGCMS/reference/HRMF.md),
[`getMSP`](https://yufree.github.io/enviGCMS/reference/getMSP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
summary <- cleanMGF("ms2.mgf", adduct = "[M+H]+")
summary[, c("title", "formula", "n_before", "n_kept", "status")]
} # }
```
