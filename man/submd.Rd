% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{submd}
\alias{submd}
\title{Get the differences of two GC/LC-MS data}
\usage{
submd(data1, data2, mzstep = 0.1, rtstep = 0.01)
}
\arguments{
\item{data1}{data file path of first data}

\item{data2}{data file path of second data}

\item{mzstep}{the m/z step for generating matrix data from raw mass spectral data}

\item{rtstep}{the alignment accuracy of retention time, e.g. 0.01 means the retention times of combined data should be the same at the accuracy 0.01s. Higher rtstep would return less scans for combined data}
}
\value{
list four matrix with the row as scantime in second and column as m/z, the first matrix refer to data 1, the second matrix refer to data 2, the third matrix refer to data1 - data2 while the fourth refer to data2 - data1, minus values are imputed by 0
}
\description{
Get the differences of two GC/LC-MS data
}
\examples{
\dontrun{
library(faahKO)
cdfpath <- system.file('cdf', package = 'faahKO')
cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
matrix <- submd(cdffiles[1],cdffiles[7])
}
}
