% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotmzrt.R
\name{plotpca}
\alias{plotpca}
\title{plot the PCA for multiple samples}
\usage{
plotpca(
  data,
  lv = NULL,
  index = NULL,
  center = TRUE,
  scale = TRUE,
  xrange = NULL,
  yrange = NULL,
  pch = NULL,
  ...
)
}
\arguments{
\item{data}{data row as peaks and column as samples}

\item{lv}{group information}

\item{index}{index for selected peaks}

\item{center}{parameters for PCA}

\item{scale}{parameters for scale}

\item{xrange}{x axis range for return samples, default NULL}

\item{yrange}{y axis range for return samples, default NULL}

\item{pch}{default pch would be the first character of group information or samples name}

\item{...}{other parameters for `plot` function}
}
\value{
if xrange and yrange are not NULL, return file name of all selected samples on 2D score plot
}
\description{
plot the PCA for multiple samples
}
\examples{
data(list)
plotpca(list$data, lv = as.character(list$group$sample_group))
}
