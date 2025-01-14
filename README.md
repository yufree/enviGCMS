enviGCMS: GC-MS Data Analysis for Environmental Science
================

[![CRAN status](http://www.r-pkg.org/badges/version/enviGCMS)](https://cran.r-project.org/package=enviGCMS) [![Download counter](http://cranlogs.r-pkg.org/badges/enviGCMS)](https://cran.r-project.org/package=enviGCMS) [![](https://cranlogs.r-pkg.org/badges/grand-total/enviGCMS)](https://cran.r-project.org/package=enviGCMS) [![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

`enviGCMS` provides functions for GC/LC-MS data analysis for environmental sciences.

Installation
------------

You can either use the stable version of `enviGCMS` from CRAN,

``` {r}
install.packages("enviGCMS")
```

or the current development snapshot from this GitHub repository:

``` {r}
remotes::install_github("yufree/enviGCMS")
```

This package has removed the depends/suggests on xcms package as it is a package with too many unstable depends. If you need to use the functions in xcms, please find the release version of enviGCMS 0.7.4.

Usage
-----

Check this [vignette](http://yufree.github.io/enviGCMS/articles/GCMSDA.html) for Data analysis of GC-MS and LC-MS in Environmental Science.

Check this [vignette](http://yufree.github.io/enviGCMS/articles/PooledQC.html) for Pooled QC analysis in Environmental Science.

Detailed usage of functions in Environmental analysis could be found in this [paper](https://doi.org/10.1016/j.talanta.2016.11.046) and the vignettes in this package.

