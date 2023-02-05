enviGCMS: GC-MS Data Analysis for Environmental Science
================

[![CRAN status](http://www.r-pkg.org/badges/version/enviGCMS)](https://cran.r-project.org/package=enviGCMS) [![Download counter](http://cranlogs.r-pkg.org/badges/enviGCMS)](https://cran.r-project.org/package=enviGCMS) [![](https://cranlogs.r-pkg.org/badges/grand-total/enviGCMS)](https://cran.r-project.org/package=enviGCMS) [![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![Build status](https://api.travis-ci.org/yufree/enviGCMS.svg?branch=master)](https://travis-ci.org/yufree/enviGCMS)

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

Usage
-----

Check this [vignette](http://yufree.github.io/enviGCMS/articles/GCMSDA.html) for Data analysis of GC-MS and LC-MS in Environmental Science.

Check this [vignette](http://yufree.github.io/enviGCMS/articles/PooledQC.html) for Pooled QC analysis in Environmental Science.

- get the mean and RSD of one sample for 5 technique replicate

~~~
# enviGCMS use functions in xcms to import the data, just type the path to your single sample
data1 <- enviGCMS:::getmd(‘sample1-1’)
data2 <- enviGCMS:::getmd(‘sample1-2’)
data3 <- enviGCMS:::getmd(‘sample1-3’)
data4 <- enviGCMS:::getmd(‘sample1-4’)
data5 <- enviGCMS:::getmd(‘sample1-5’)
~~~

- get the mean

~~~
data <- (data1+data2+data3+data4+data5)/5
~~~

- get the standard deviation

~~~
datasd <- sqrt(((data1-data)^2+(data2-data)^2+(data3-data)^2+(data4-data)^2+(data5-data)^2)/4)
~~~

- get the RSD

~~~
databrsd <- datasd/data
~~~

- plot the smooth scatter

~~~
plotsms(datarsd)
~~~

- plot the heatmap

~~~
plotms(data)
~~~

- plot the mz-rt scatter plot

~~~
plotmz(data)
~~~

- plot the boundary model

~~~
findline(data)
~~~

Detailed usage of those functions in Environmental analysis could be found in this [paper](https://www.sciencedirect.com/science/article/pii/S0039914016309298) and the vignettes in this package.

