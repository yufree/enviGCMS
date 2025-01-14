#' Short-Chain Chlorinated Paraffins(SCCPs) peaks information for quantitative analysis
#'
#' A dataset containing the ions, formula, Cl%, peak abundances for all SCCPs compounds
#' @docType data
#' @usage data(sccp)
#' @format A data frame with 24 rows and 8 variables:
#' \describe{
#'   \item{Cln}{Chlorine atom numbers}
#'   \item{Cn}{Carbon atom numbers}
#'   \item{formula}{molecular formula}
#'   \item{Hn}{hydrogen atom numbers}
#'   \item{ions}{[M-Cl]- ions}
#'   \item{mz}{m/z for the isotopologues with highest intensity}
#'   \item{intensity}{abundance of the isotopologues with highest intensity}
#'   \item{Clp}{Chlorine contents}
#' }
"sccp"
