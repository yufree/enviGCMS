# Find lipid class of metabolites base on referenced Kendrick mass defect

Find lipid class of metabolites base on referenced Kendrick mass defect

## Usage

``` r
findlipid(list, mode = "pos")
```

## Arguments

- list:

  list with data as peaks list, mz, rt and group information, retention
  time should be in seconds

- mode:

  'pos' for positive mode, 'neg' for negative mode and 'none' for
  neutral mass, only support \[M+H\] and \[M-H\] for each mode

## Value

list list with dataframe with the lipid referenced Kendrick mass
defect(RKMD) and logical for class

## References

Method for the Identification of Lipid Classes Based on Referenced
Kendrick Mass Analysis. Lerno LA, German JB, Lebrilla CB. Anal Chem.
2010 May 15;82(10):4236–45.

## Examples

``` r
data(list)
RKMD <- findlipid(list)
```
