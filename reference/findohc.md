# Screen organohalogen compounds

Screen organohalogen compounds by retention time, mass defect analysis
and isotope relationship modified by literature report. Also support
compounds with \[M\] and \[M+2\] ratio cutoff.

## Usage

``` r
findohc(
  list,
  sf = 78/77.91051,
  step = 0.001,
  stepsd1 = 0.003,
  stepsd2 = 0.005,
  mzc = 700,
  cutoffint = 1000,
  cutoffr = 0.4,
  clustercf = 10
)
```

## Arguments

- list:

  list with data as peaks list, mz, rt and group information, retention
  time should be in seconds

- sf:

  scale factor, default 78/77.91051(Br)

- step:

  mass defect step, default 0.001

- stepsd1:

  mass defect uncertainty for lower mass, default 0.003

- stepsd2:

  mass defect uncertainty for higher mass, default 0.005

- mzc:

  threshold of lower mass and higher mass, default 700

- cutoffint:

  the cutoff of intensity, default 1000

- cutoffr:

  the cutoff of \[M\] and \[M+2\] ratio, default 0.4

- clustercf:

  the cutoff of cluster analysis to separate two different ions groups
  for retention time, default 10

## Value

list with filtered organohalogen compounds

## References

Identification of Novel Brominated Compounds in Flame Retarded Plastics
Containing TBBPA by Combining Isotope Pattern and Mass Defect Cluster
Analysis Ana Ballesteros-Gomez, Joaquin Ballesteros, Xavier Ortiz,
Willem Jonker, Rick Helmus, Karl J. Jobst, John R. Parsons, and Eric J.
Reiner Environmental Science & Technology 2017 51 (3), 1518-1526 DOI:
10.1021/acs.est.6b03294
