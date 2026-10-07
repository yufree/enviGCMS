# Extract top MS2 ions from MS1 EIC interactively

A Shiny app that displays extracted ion chromatograms from MS1 data.
Click on a peak to select a retention time window, then extract the most
intense MS2 ions from that window and overlay their chromatograms.

## Usage

``` r
plotTopMS2Peaks(
  filepath,
  featlist,
  numTopIons = 10,
  diff = 0.01,
  rtWindow = 0.3
)
```

## Arguments

- filepath:

  character. Path to an mzML or mzXML file.

- featlist:

  data.frame with columns: `name`, `mz`, `ms_level`.

- numTopIons:

  integer. Number of most intense MS1 ions to extract. Default 10.

- diff:

  numeric. Half-width of m/z extraction window (Da). Default 0.01.

- rtWindow:

  numeric. Half-width of RT window (seconds) around clicked peak.
  Default 0.3.

## Value

Opens an interactive Shiny app in the browser.

## Examples

``` r
if (FALSE) { # \dontrun{
targets <- data.frame(
  name = c("precursor1"),
  mz = c(500.1234),
  ms_level = c("ms1")
)
plotTopMS2Peaks("sample.mzML", targets, numTopIons = 5)
} # }
```
