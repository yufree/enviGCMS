#' @useDynLib enviGCMS, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Internal wrapper for mass decomposition using the HORIZON algorithm
#'
#' Powered by the HORIZON (Heavy-first Ordered Recursive Inference with
#' Zero-loop Optimal Navigation) engine. Supports both single mass decomposition
#' and vectorized multi-threaded batch decomposition with Fiehn Seven Golden Rules
#' filtering and high-resolution instrument peak merging.
#'
#' @param mass Numeric mass or vector of masses to decompose.
#' @param ppm Allowed deviation in ppm (default: 2.0).
#' @param mzabs Absolute deviation in Da (default: 0.0001).
#' @param elements Allowed chemical elements (character string or list).
#' @param minElements Lower bounds formula string.
#' @param maxElements Upper bounds formula string.
#' @param z Charge state (default: 0).
#' @param maxisotopes Maximum number of isotopes (default: 10).
#' @param golden_rules Logical, whether to apply Fiehn Seven Golden Rules (default: TRUE).
#' @param resolution Instrument resolving power (m / FWHM) for merging fine isotope peaks (default: 0).
#' @param fwhm Peak width at half maximum in Da for merging isotope peaks (default: 0).
#' @param nthreads Number of parallel threads for multi-mass decomposition (default: 1).
#' @return List of decomposition results compatible with Rdisop.
#' @noRd
.decomposeMass <- function(mass, ppm = 2.0, mzabs = 0.0001, elements = NULL,
                           minElements = "", maxElements = "", z = 0,
                           maxisotopes = 10, golden_rules = TRUE,
                           resolution = 0.0, fwhm = 0.0, nthreads = 1) {
    if (is.null(minElements)) minElements <- ""
    if (is.null(maxElements)) maxElements <- ""

    if (length(mass) <= 1) {
        rcpp_decompose_mass(
            mass = if (length(mass) == 0) 0.0 else as.numeric(mass),
            ppm = as.numeric(ppm),
            mzabs = as.numeric(mzabs),
            elements = elements,
            minElements = as.character(minElements),
            maxElements = as.character(maxElements),
            z = as.integer(z),
            maxisotopes = as.integer(maxisotopes),
            golden_rules = as.logical(golden_rules),
            resolution = as.numeric(resolution),
            fwhm = as.numeric(fwhm)
        )
    } else {
        rcpp_decompose_masses(
            masses = as.numeric(mass),
            ppm = as.numeric(ppm),
            mzabs = as.numeric(mzabs),
            elements = elements,
            minElements = as.character(minElements),
            maxElements = as.character(maxElements),
            z = as.integer(z),
            maxisotopes = as.integer(maxisotopes),
            golden_rules = as.logical(golden_rules),
            resolution = as.numeric(resolution),
            fwhm = as.numeric(fwhm),
            nthreads = as.integer(nthreads)
        )
    }
}

#' Internal wrapper for molecule mass and isotope calculation using the HORIZON engine
#'
#' @param formula Character formula string.
#' @param elements Allowed chemical elements (optional).
#' @param z Charge state (default: 0).
#' @param maxisotopes Maximum number of isotopes (default: 10).
#' @param resolution Instrument resolving power (default: 0).
#' @param fwhm Peak width at half maximum in Da (default: 0).
#' @return List of molecule properties compatible with Rdisop.
#' @noRd
.getMolecule <- function(formula, elements = NULL, z = 0, maxisotopes = 10,
                         resolution = 0.0, fwhm = 0.0) {
    rcpp_get_molecule(
        formula = as.character(formula),
        z = as.integer(z),
        maxisotopes = as.integer(maxisotopes),
        resolution = as.numeric(resolution),
        fwhm = as.numeric(fwhm)
    )
}

#' Check Chemical Formula against Fiehn Seven Golden Rules
#'
#' Evaluates whether a candidate chemical formula satisfies heuristic chemical
#' rules based on Senior's valence rules (Senior 1951), Hydrogen-to-Carbon (H/C) ratios,
#' heteroatom ratios (N/C, O/C, P/C, S/C, Halogen/C), and Double Bond Equivalents (DBE)
#' as described by Kind & Fiehn (2007).
#'
#' @param formula Character string representing molecular formula (e.g. \code{"C6H12O6"}).
#' @param z Integer charge state (default: 0).
#' @param min_hc Numeric minimum H/C ratio (default: 0.1).
#' @param max_hc Numeric maximum H/C ratio (default: 6.0).
#' @param max_dbe Numeric maximum DBE (default: 40.0).
#' @return A list containing booleans indicating whether the formula passed Senior rules,
#'   H/C ratio, heteroatom ratios, DBE, nitrogen rule, and overall validity.
#' @references
#' Kind, T., Fiehn, O. Seven Golden Rules for heuristic filtering of molecular formulas
#' obtained by accurate mass spectrometry. BMC Bioinformatics 8, 105 (2007).
#' \doi{10.1186/1471-2105-8-105}
#'
#' Senior, J. K. Partitions and their connection with the problems of chemical
#' structure. J. Chem. Phys. 19, 865-873 (1951).
#' @examples
#' checkGoldenRules("C6H12O6")
#' checkGoldenRules("C3H32N12O3") # chemically impossible formula
#' @export
checkGoldenRules <- function(formula, z = 0, min_hc = 0.1, max_hc = 6.0, max_dbe = 40.0) {
    rcpp_check_golden_rules(
        formula = as.character(formula),
        z = as.integer(z),
        min_hc = as.numeric(min_hc),
        max_hc = as.numeric(max_hc),
        max_dbe = as.numeric(max_dbe)
    )
}

#' Calculate Isotopic Pattern Similarity Scores
#'
#' Computes modern spectral similarity metrics including Cosine similarity,
#' intensity-weighted dot product, and log-likelihood between observed mass
#' spectral peaks and theoretical isotopic distributions.
#'
#' @param obs_mz Numeric vector of observed m/z values.
#' @param obs_int Numeric vector of observed intensities.
#' @param theo_mz Numeric vector of theoretical isotopic m/z values.
#' @param theo_int Numeric vector of theoretical isotopic abundances.
#' @param tolerance Numeric mass tolerance for peak alignment (default: 0.005).
#' @param ppm Logical indicating whether \code{tolerance} is in ppm (default: FALSE).
#' @return A list containing:
#'   \item{cosine}{Cosine similarity (unweighted dot product of normalized spectra).}
#'   \item{weighted_cosine}{Mass- and intensity-weighted dot product (MassBank/NIST style).}
#'   \item{log_likelihood}{Multinomial log-likelihood of observed peak distribution.}
#'   \item{matched_peaks}{Integer count of successfully matched peaks.}
#' @examples
#' obs_m <- c(180.0634, 181.0668)
#' obs_i <- c(100, 6.5)
#' theo_m <- c(180.0634, 181.0668)
#' theo_i <- c(1.0, 0.066)
#' scoreIsotopes(obs_m, obs_i, theo_m, theo_i)
#' @export
scoreIsotopes <- function(obs_mz, obs_int, theo_mz, theo_int, tolerance = 0.005, ppm = FALSE) {
    rcpp_score_isotopes(
        obs_mz = as.numeric(obs_mz),
        obs_int = as.numeric(obs_int),
        theo_mz = as.numeric(theo_mz),
        theo_int = as.numeric(theo_int),
        tolerance = as.numeric(tolerance),
        ppm = as.logical(ppm)
    )
}
