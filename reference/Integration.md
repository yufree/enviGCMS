# Just integrate data according to fixed rt and fixed noise area

Just integrate data according to fixed rt and fixed noise area

## Usage

``` r
integration(data, rt = c(8.3, 9), brt = c(8.3, 8.4), smoothit = TRUE)
```

## Arguments

- data:

  file should be a dataframe with the first column RT and second column
  intensity of the SIM ions.

- rt:

  a rough RT range contained only one peak to get the area

- brt:

  a rough RT range contained only one peak and enough noises to get the
  area

- smoothit:

  logical, if using an average smooth box or not. If using, n will be
  used

## Value

area integration data

## Examples

``` r
if (FALSE) { # \dontrun{
area <- integration(data)
} # }
```
