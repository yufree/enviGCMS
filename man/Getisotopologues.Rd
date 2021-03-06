% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peaks.R
\name{Getisotopologues}
\alias{Getisotopologues}
\title{Get the selected isotopologues at certain MS data}
\usage{
Getisotopologues(formula = "C12OH6Br4", charge = 1, width = 0.3)
}
\arguments{
\item{formula}{the molecular formula. C12OH6Br4 means BDE-47 as default}

\item{charge}{the charge of that molecular. 1 in EI mode as default}

\item{width}{the width of the peak width on mass spectrum. 0.3 as default for low resolution mass spectrum.}
}
\description{
Get the selected isotopologues at certain MS data
}
\examples{
# show isotopologues for BDE-47
Getisotopologues(formula = 'C12OH6Br4')
}
