# getintegration was mainly used for get the integration of certain ion's chromatogram data and plot the data

getintegration was mainly used for get the integration of certain ion's
chromatogram data and plot the data

## Usage

``` r
getintegration(
  data,
  rt = c(8.3, 9),
  n = 5,
  m = 5,
  slope = c(2, 2),
  baseline = 10,
  noslope = TRUE,
  smoothit = TRUE,
  half = FALSE
)
```

## Arguments

- data:

  file should be a dataframe with the first column RT and second column
  intensity of the SIM ions.

- rt:

  a rough RT range contained only one peak to get the area

- n:

  points in the moving average smooth box, default value is 5

- m:

  numbers of points for regression to get the slope

- slope:

  the threshold value for start/stop peak as percentage of max slope

- baseline:

  numbers of the points for the baseline of the signal

- noslope:

  logical, if using a horizon line to get area or not

- smoothit:

  logical, if using an average smooth box or not. If using, n will be
  used

- half:

  logical, if using the left half peak to calculate the area

## Value

integration data such as peak area, peak height, signal and the slope
data.

## Examples

``` r
if (FALSE) { # \dontrun{
list <- getintegration(data)
} # }
```
