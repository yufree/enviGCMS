% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anno.R
\name{getalign2}
\alias{getalign2}
\title{Align mass to charge ratio and/or retention time to remove redundancy}
\usage{
getalign2(mz, rt, ppm = 5, deltart = 5)
}
\arguments{
\item{mz}{the mass to charge of reference peaks}

\item{rt}{retention time of reference peaks}

\item{ppm}{mass accuracy, default 10}

\item{deltart}{retention time shift table, default 10 seconds}
}
\value{
index for
}
\description{
Align mass to charge ratio and/or retention time to remove redundancy
}
\examples{
mz <- c(221.1171, 221.1170, 229.1546, 233.1497, 271.0790 )
rt <- c(590.8710, 587.3820, 102.9230, 85.8850, 313.8240)
getalign2(mz,rt)
}
