# Annotation of MS1 data by compounds database by predefined paired mass distance

Annotation of MS1 data by compounds database by predefined paired mass
distance

## Usage

``` r
getms1anno(pmd, mz, ppm = 10, db = NULL)
```

## Arguments

- pmd:

  adducts formula or paired mass distance for ions

- mz:

  unknown mass to charge ratios vector

- ppm:

  mass accuracy

- db:

  compounds database as dataframe. Two required columns are name and
  monoisotopic molecular weight with column names of name and mass

## Value

list or data frame
