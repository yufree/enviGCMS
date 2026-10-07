# Get the overlap peaks by mass and retention time range

Get the overlap peaks by mass and retention time range

## Usage

``` r
getoverlappeak(list1, list2)
```

## Arguments

- list1:

  list with data as peaks list, mz, rt, mzrange, rtrange and group
  information to be overlapped

- list2:

  list with data as peaks list, mz, rt, mzrange, rtrange and group
  information to overlap

## Value

logical index for list 1's peaks

## See also

[`getimputation`](https://yufree.github.io/enviGCMS/reference/getimputation.md),[`getdoe`](https://yufree.github.io/enviGCMS/reference/getdoe.md)
