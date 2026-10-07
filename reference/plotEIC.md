# Plot extracted ion chromatograms from raw mzML/mzXML data

Reads mzML/mzXML raw data and plots extracted ion chromatograms (EIC)
for a list of target m/z values at MS1 or MS2 level. Uses RaMS for
lightweight data access and plotly for interactive visualization.

## Usage

``` r
plotEIC(filepath, featlist, diff = 0.005)
```

## Arguments

- filepath:

  character. Path to an mzML or mzXML file.

- featlist:

  data.frame with columns: `name` (feature label), `mz` (target m/z),
  `ms_level` (either `"ms1"` or `"ms2"`).

- diff:

  numeric. Half-width of the m/z extraction window in Da. Default 0.005.

## Value

A plotly object with overlaid extracted ion chromatograms.

## Examples

``` r
if (FALSE) { # \dontrun{
featlist <- data.frame(
  name = c("target1", "target2"),
  mz = c(300.1234, 350.5678),
  ms_level = c("ms1", "ms2")
)
plotEIC("sample.mzML", featlist)
} # }
```
