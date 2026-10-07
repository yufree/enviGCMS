# Plot the response group of GC-MS

Plot the response group of GC-MS

## Usage

``` r
plotgroup(data, threshold = 2)
```

## Arguments

- data:

  imported data matrix of GC-MS

- threshold:

  the threshold of the response (log based 10) to separate the group

## Value

list linear regression model for the data matrix

## Examples

``` r
if (FALSE) { # \dontrun{
data(matrix)
plotgroup(matrix)
} # }
```
