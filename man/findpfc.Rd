% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mda.R
\name{findpfc}
\alias{findpfc}
\title{Find PFCs based on mass defect analysis}
\usage{
findpfc(list)
}
\arguments{
\item{list}{list with data as peaks list, mz, rt and group information, retention time should be in seconds}
}
\value{
list list with potential PFCs compounds index
}
\description{
Find PFCs based on mass defect analysis
}
\examples{
data(list)
pfc <- findpfc(list)
}
\references{
Liu, Y.; D’Agostino, L. A.; Qu, G.; Jiang, G.; Martin, J. W. High-Resolution Mass Spectrometry (HRMS) Methods for Nontarget Discovery and Characterization of Poly- and per-Fluoroalkyl Substances (PFASs) in Environmental and Human Samples. TrAC Trends in Analytical Chemistry 2019, 121, 115420.
}
