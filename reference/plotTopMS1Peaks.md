# Extract top MS1 ions from MS2 EIC interactively

A Shiny app that displays extracted ion chromatograms from MS2 data.
Click on a peak to select a retention time window, then extract the most
intense MS1 ions from that window and overlay their chromatograms.

## Usage

``` r
plotTopMS1Peaks(
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
frags <- data.frame(
  name = c("Br79", "Br81"),
  mz = c(78.9183, 80.9163),
  ms_level = c("ms2", "ms2")
)
plotTopMS1Peaks("sample.mzML", frags, numTopIons = 3)
} # }
```
