# Get the high order unit based Mass Defect

Get the high order unit based Mass Defect

## Usage

``` r
getmdh(mz, cus = c("CH2,H2"), method = "round")
```

## Arguments

- mz:

  numeric vector for exact mass

- cus:

  chemical formula or reaction

- method:

  you could use \`round\`, \`floor\` or \`ceiling\`

## Value

high order Mass Defect with details

## Examples

``` r
if (FALSE) { # \dontrun{
getmdh(getmass('C2H4'))
} # }
```
