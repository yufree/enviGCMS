% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotmzrt.R
\name{plotridges}
\alias{plotridges}
\title{Relative Log Abundance Ridge (RLAR) plots for samples or peaks}
\usage{
plotridges(data, lv, type = "g")
}
\arguments{
\item{data}{data row as peaks and column as samples}

\item{lv}{factor vector for the group information of samples}

\item{type}{'g' means group median based, other means all samples median based.}
}
\value{
Relative Log Abundance Ridge(RLA) plots
}
\description{
Relative Log Abundance Ridge (RLAR) plots for samples or peaks
}
\examples{
data(list)
plotridges(list$data, as.factor(list$group$sample_group))
}
