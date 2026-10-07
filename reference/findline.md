# find line of the regression model for GC-MS

find line of the regression model for GC-MS

## Usage

``` r
findline(data, threshold = 2, temp = c(100, 320))
```

## Arguments

- data:

  imported data matrix of GC-MS

- threshold:

  the threshold of the response (log based 10)

- temp:

  the scale of the oven temperature (constant rate)

## Value

list linear regression model for the matrix

## Examples

``` r
if (FALSE) { # \dontrun{
data(matrix)
findline(matrix)
} # }
```
