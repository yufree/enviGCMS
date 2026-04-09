#' Parse chemical formula into element-count pairs
#'
#' @param formula character string of a chemical formula (e.g. "C8H11NO")
#' @return data.frame with columns \code{element} and \code{count}
#' @noRd
.parse_formula <- function(formula) {
        # Match element symbols (uppercase + optional lowercase) followed by optional count
        elements <- gregexpr("[A-Z][a-z]?\\d*", formula)
        tokens <- regmatches(formula, elements)[[1]]
        elem <- sub("\\d+$", "", tokens)
        count <- as.integer(sub("^[A-Za-z]+", "", tokens))
        count[is.na(count)] <- 1L
        data.frame(element = elem, count = count, stringsAsFactors = FALSE)
}

#' Electron mass constant (Da)
#' @noRd
.ELECTRON_MASS <- 0.00054857990924

#' Match query m/z values against target m/z ranges
#'
#' For each query m/z, find the target row whose (min, max) range contains it.
#' Returns columns from \code{target_values} for the last match (mimicking
#' the original nested-loop semantics where later matches overwrite earlier ones).
#'
#' @param query_mz numeric vector of m/z values to look up.
#' @param target_min numeric vector of lower bounds for target ranges.
#' @param target_max numeric vector of upper bounds for target ranges.
#' @param target_values data.frame with one row per target entry; columns are
#'   extracted on match.
#' @return A data.frame with \code{length(query_mz)} rows and the same columns
#'   as \code{target_values}, filled with 0 / "" defaults where no match occurred.
#' @noRd
.match_mz_peaks <- function(query_mz, target_min, target_max, target_values) {
        n_query  <- length(query_mz)
        n_target <- length(target_min)

        # Pre-fill with type-appropriate defaults
        out <- as.data.frame(
                lapply(target_values, function(col) {
                        if (is.character(col)) rep("", n_query)
                        else rep(0, n_query)
                }),
                stringsAsFactors = FALSE
        )

        for (i in seq_len(n_query)) {
                for (j in seq_len(n_target)) {
                        if (query_mz[i] > target_min[j] & query_mz[i] < target_max[j]) {
                                out[i, ] <- target_values[j, , drop = FALSE]
                        }
                }
        }
        out
}

#' Calculate HRMF forward, reverse and FoM scores
#'
#' @param all_ions data.frame of theoretical ions with matched experimental data.
#' @param comp_scored data.frame of experimental peaks with matched theoretical data.
#' @param formula_name character scalar, the candidate formula label.
#' @return A one-row data.frame with HRMF score columns.
#' @noRd
.calc_hrmf_scores <- function(all_ions, comp_scored, formula_name) {
        # Forward HRMF score (theoretical -> experimental)
        all_iso <- nrow(all_ions)
        n_matched_forw <- sum(!is.na(all_ions$MassError_ppm))
        HRMF_theor_score <- round(
                sum(all_ions$Detected_mz * all_ions$Detected_RelAb, na.rm = TRUE) /
                        sum(all_ions$Iso_mz * all_ions$Theor_RelAb, na.rm = TRUE), 2
        )

        # Reverse HRMF score (experimental -> theoretical)
        all_iso_cmp <- nrow(comp_scored)
        n_matched_rev <- sum(!is.na(comp_scored$MassError_ppm_rev))
        HRMF_msp_score <- round(
                sum(comp_scored$Theor_mz * comp_scored$Theor_RelAb, na.rm = TRUE) /
                        sum(comp_scored$mz * comp_scored$Detected_RelAb, na.rm = TRUE), 2
        )

        # Figure of Merit (FoM)
        det_relab <- all_ions$Detected_RelAb
        det_relab[is.na(det_relab)] <- 0
        theor_relab <- all_ions$Theor_RelAb
        denom <- pmax(theor_relab, det_relab)
        denom[denom == 0] <- 1
        fom_terms <- abs(theor_relab - det_relab) / denom
        FoM <- round(1 - mean(fom_terms), 2)

        data.frame(
                Candidate = formula_name,
                peak_count_forw = all_iso,
                df_theortomsp = round(n_matched_forw / all_iso * 100, 0),
                HRMF_theor_score = HRMF_theor_score,
                peak_count_rev = all_iso_cmp,
                df_msptotheor = round(n_matched_rev / all_iso_cmp * 100, 0),
                HRMF_msp_score = HRMF_msp_score,
                FoM = FoM,
                stringsAsFactors = FALSE
        )
}

#' Annotate experimental peaks with sub-formula matches and isotope patterns
#'
#' For each experimental peak, uses Rdisop to decompose the mass into
#' sub-formulae consistent with the candidate formula constraints. Returns
#' an \code{all_ions} data.frame of theoretical isotopologues, or NULL if
#' no matches are found.
#'
#' @param compound data.frame with columns mz, intensity, mz_min, mz_max.
#' @param element_str character, concatenated element symbols.
#' @param min_str character, Rdisop minElements string.
#' @param max_str character, Rdisop maxElements string.
#' @param me numeric, electron mass adjustment.
#' @param charge integer, charge state.
#' @param mass_accuracy numeric, mass accuracy in ppm.
#' @param IR_RelAb_cutoff numeric, relative abundance cutoff for isotopologues.
#' @param formula_label character, formula name for messages.
#' @return A data.frame of theoretical isotopologues or NULL.
#' @noRd
.annotate_peaks <- function(compound, element_str, min_str, max_str,
                            me, charge, mass_accuracy, IR_RelAb_cutoff,
                            formula_label) {
        windows <- mass_accuracy / 1e6 * compound$mz
        match_comp <- vector("list", nrow(compound))

        for (i in seq_len(nrow(compound))) {
                match_list <- list(annotated = FALSE, annodf = NULL, isopat = NULL)

                mfSet <- tryCatch(
                        Rdisop::decomposeMass(
                                compound$mz[i] + me,
                                mzabs = windows[i],
                                z = charge,
                                elements = element_str,
                                minElements = min_str,
                                maxElements = max_str
                        ),
                        error = function(e) list(formula = character(0), valid = character(0))
                )

                valid_idx <- which(mfSet$valid == "Valid")

                if (length(valid_idx) > 0) {
                        matched_formula <- mfSet$formula[valid_idx[1]]
                        matched_mass <- mfSet$exactmass[valid_idx[1]]

                        match_list$annotated <- TRUE
                        match_list$annodf <- data.frame(
                                MonoIso_mz = matched_mass,
                                MonoIonFormula = matched_formula,
                                MassError_ppm = round(
                                        (compound$mz[i] - matched_mass) / matched_mass * 1e6, 2
                                ),
                                stringsAsFactors = FALSE
                        )

                        mol <- tryCatch(
                                Rdisop::getMolecule(matched_formula,
                                                    z = charge,
                                                    maxisotopes = 20),
                                error = function(e) NULL
                        )

                        if (!is.null(mol) && !is.null(mol$isotopes[[1]])) {
                                iso_mat <- t(mol$isotopes[[1]])
                                iso_df <- data.frame(
                                        Iso_mz = iso_mat[, 1],
                                        Abundance = round(iso_mat[, 2] / max(iso_mat[, 2]) * 100, 1),
                                        stringsAsFactors = FALSE
                                )
                                iso_df <- iso_df[iso_df$Abundance >= IR_RelAb_cutoff, , drop = FALSE]
                                iso_df$MonoIsoFormula <- matched_formula
                                iso_df$Isoformula <- matched_formula
                                match_list$isopat <- iso_df[, c("MonoIsoFormula", "Iso_mz",
                                                                "Abundance", "Isoformula")]
                        }
                }
                match_comp[[i]] <- match_list
        }

        annotated_mask <- vapply(match_comp, function(x) x$annotated, logical(1))

        if (!any(annotated_mask)) {
                message("No match for any peaks for formula: ", formula_label,
                        " within ", mass_accuracy, " ppm. Skipping.")
                return(NULL)
        }

        iso_list <- lapply(match_comp[annotated_mask], function(x) x$isopat)
        iso_list <- iso_list[!vapply(iso_list, is.null, logical(1))]

        if (length(iso_list) == 0) {
                message("Isotope pattern calculation failed for formula: ",
                        formula_label, ". Skipping.")
                return(NULL)
        }

        do.call(rbind, iso_list)
}

#' Compute relative abundances per mono-isotopic group
#'
#' @param df data.frame with columns MonoIsoFormula, intensity (or
#'   Detected_int), and Abundance (theoretical relative abundance).
#' @param int_col character, name of the intensity column.
#' @param has_abundance logical, if TRUE also compute Expected_int from Abundance.
#' @return The input data.frame with added Expected_int, Detected_RelAb columns.
#' @noRd
.calc_group_relab <- function(df, int_col, has_abundance = FALSE) {
        mono_groups <- unique(df$MonoIsoFormula)
        if (has_abundance) df$Expected_int <- 0
        df$Detected_RelAb <- 0

        for (mg in mono_groups) {
                if (mg == "") next
                idx <- df$MonoIsoFormula == mg
                max_int <- max(df[[int_col]][idx], na.rm = TRUE)
                if (max_int > 0) {
                        if (has_abundance) {
                                df$Expected_int[idx] <- round(
                                        df$Abundance[idx] / 100 * max_int, 0
                                )
                        }
                        df$Detected_RelAb[idx] <- round(
                                df[[int_col]][idx] / max_int * 100, 1
                        )
                }
        }
        df
}

#' High Resolution Mass Filtering (HRMF) for GC/LC HRMS data
#'
#' Performs high-resolution mass filtering by matching experimental mass spectral
#' peaks against theoretical isotope patterns derived from one or more candidate
#' chemical formulae. Calculates forward HRMF score, reverse (MSP) score, and
#' Figure of Merit (FoM) for each candidate.
#'
#' The method is based on Kwiecien et al. (2015) \doi{10.1021/acs.analchem.5b01503}.
#'
#' @param msp list. A single compound entry as returned by \code{\link{getMSP}},
#'   containing at least a \code{spectra} element with columns \code{mz} and \code{intensity}.
#' @param formula character vector. One or more candidate chemical formulae to evaluate.
#' @param charge integer. Charge state: 1 for positive, -1 for negative, 0 for neutral.
#'   Default 1 (radical cation in EI).
#' @param mass_accuracy numeric. Mass accuracy in ppm. Default 5.
#' @param intensity_cutoff numeric. Minimum absolute intensity to retain a peak. Default 1.
#' @param IR_RelAb_cutoff numeric. Relative abundance cutoff (\%) for theoretical
#'   isotopologues. Default 1.
#' @param detailed logical. If TRUE, return detailed list with all_ions, compound,
#'   and HRMF_scores for each formula. If FALSE (default), return a summary data.frame.
#'
#' @return If \code{detailed = FALSE}, a data.frame with one row per candidate formula
#'   and columns: Candidate, peak_count_forw, df_theortomsp, HRMF_theor_score,
#'   peak_count_rev, df_msptotheor, HRMF_msp_score, FoM.
#'   If \code{detailed = TRUE}, a named list of detailed results per formula.
#'
#' @details
#' Unlike the original MSxplorer implementation, this version uses \pkg{Rdisop}
#' (already a dependency of enviGCMS) for both formula decomposition and isotope
#' pattern calculation, instead of \pkg{rcdk}/\pkg{rJava}/\pkg{enviPat},
#' requiring no additional dependencies.
#'
#' The input \code{msp} should be a single entry from the list returned by
#' \code{\link{getMSP}}. For batch processing of entire MSP files, see
#' \code{\link{getHRMF}}.
#'
#' @seealso \code{\link{getHRMF}} for batch processing, \code{\link{getMSP}} for
#'   reading MSP files.
#'
#' @examples
#' \dontrun{
#' # Read MSP file and run HRMF on the first compound
#' msp_data <- getMSP("spectrum.msp")
#' result <- HRMF(msp_data[[1]], formula = "C8H11NO")
#'
#' # Compare multiple candidates
#' result <- HRMF(msp_data[[1]], formula = c("C8H11NO", "C7H9NO2"))
#' }
#' @export
HRMF <- function(msp, formula, charge = 1, mass_accuracy = 5,
                 intensity_cutoff = 1, IR_RelAb_cutoff = 1,
                 detailed = FALSE) {

        # Extract spectra from getMSP format
        if (is.null(msp$spectra)) {
                stop("Input 'msp' must contain a 'spectra' element with columns 'mz' and 'intensity'.",
                     call. = FALSE)
        }
        compound <- data.frame(mz = msp$spectra$mz,
                               intensity = msp$spectra$intensity,
                               stringsAsFactors = FALSE)

        # Apply mass accuracy windows and intensity filter
        compound$mz_min <- compound$mz - mass_accuracy / 1e6 * compound$mz
        compound$mz_max <- compound$mz + mass_accuracy / 1e6 * compound$mz
        compound <- compound[compound$intensity > intensity_cutoff, , drop = FALSE]

        if (nrow(compound) == 0) {
                warning("No peaks remain after intensity filtering.", call. = FALSE)
                return(NULL)
        }

        # Electron mass adjustment for charged species
        if (charge > 0) {
                me <- .ELECTRON_MASS
        } else if (charge == 0) {
                me <- 0
        } else {
                me <- -.ELECTRON_MASS
        }

        all_hrmf <- list()

        # Loop through candidate formulae
        for (a in seq_along(formula)) {

                # Parse the candidate formula into elements and counts
                atoms <- .parse_formula(formula[a])
                element_str <- paste(atoms$element, collapse = "")
                min_str <- paste0(paste0(atoms$element, "0"), collapse = "")
                max_str <- paste0(paste0(atoms$element, atoms$count), collapse = "")

                # Annotate experimental peaks with sub-formulae and isotope patterns
                all_ions <- .annotate_peaks(compound, element_str, min_str, max_str,
                                            me, charge, mass_accuracy, IR_RelAb_cutoff,
                                            formula[a])
                if (is.null(all_ions)) next

                # Match theoretical isotope m/z against experimental peaks
                matched_exp <- .match_mz_peaks(
                        query_mz      = all_ions$Iso_mz,
                        target_min    = compound$mz_min,
                        target_max    = compound$mz_max,
                        target_values = compound[, c("mz", "intensity"), drop = FALSE]
                )
                all_ions$Detected_mz  <- matched_exp$mz
                all_ions$Detected_int <- matched_exp$intensity

                # Calculate mass errors and relative abundances for theoretical ions
                windows2 <- mass_accuracy / 1e6 * all_ions$Iso_mz
                all_ions$Iso_mz_min <- all_ions$Iso_mz - windows2
                all_ions$Iso_mz_max <- all_ions$Iso_mz + windows2

                all_ions$MassError_ppm <- ifelse(
                        all_ions$Detected_mz > 1,
                        round((all_ions$Detected_mz - all_ions$Iso_mz) /
                                      all_ions$Iso_mz * 1e6, 1),
                        NA
                )

                # Expected intensity and relative abundances per mono-isotopic group
                all_ions <- .calc_group_relab(all_ions, "Detected_int",
                                             has_abundance = TRUE)

                # Rename for clarity
                all_ions$Theor_RelAb <- all_ions$Abundance
                all_ions$Abundance <- NULL

                # Filter by expected intensity
                all_ions <- all_ions[all_ions$Expected_int > intensity_cutoff, , drop = FALSE]

                if (nrow(all_ions) == 0) {
                        message("No ions passed intensity filter for formula: ",
                                formula[a], ". Skipping.")
                        next
                }

                # Match experimental peaks back to theoretical (reverse direction)
                comp_scored <- compound
                matched_theor <- .match_mz_peaks(
                        query_mz      = comp_scored$mz,
                        target_min    = all_ions$Iso_mz_min,
                        target_max    = all_ions$Iso_mz_max,
                        target_values = all_ions[, c("Iso_mz", "Isoformula",
                                                     "MonoIsoFormula", "Theor_RelAb"),
                                                 drop = FALSE]
                )
                comp_scored$Theor_mz       <- matched_theor$Iso_mz
                comp_scored$IsoFormula     <- matched_theor$Isoformula
                comp_scored$MonoIsoFormula <- matched_theor$MonoIsoFormula
                comp_scored$Theor_RelAb    <- matched_theor$Theor_RelAb

                comp_scored$MassError_ppm_rev <- ifelse(
                        comp_scored$Theor_mz > 1,
                        round((comp_scored$mz - comp_scored$Theor_mz) /
                                      comp_scored$Theor_mz * 1e6, 1),
                        NA
                )

                # Relative abundance for compound peaks per mono-isotopic group
                comp_scored <- .calc_group_relab(comp_scored, "intensity")

                HRMF_scores <- .calc_hrmf_scores(all_ions, comp_scored, formula[a])

                all_hrmf[[formula[a]]] <- list(
                        all_ions = all_ions,
                        compound = comp_scored,
                        HRMF_scores = HRMF_scores
                )
        }

        if (length(all_hrmf) == 0) {
                message("No formulae matched any peaks.")
                return(NULL)
        }

        if (detailed) {
                return(all_hrmf)
        } else {
                scores <- do.call(rbind, lapply(all_hrmf, function(x) x$HRMF_scores))
                rownames(scores) <- NULL
                return(scores)
        }
}


#' Batch High Resolution Mass Filtering for all compounds in an MSP file
#'
#' Applies \code{\link{HRMF}} to every compound in an MSP file that has a
#' chemical formula annotation.
#'
#' @param file character. Path to an MSP file.
#' @param charge integer. Charge state (see \code{\link{HRMF}}). Default 1.
#' @param mass_accuracy numeric. Mass accuracy in ppm. Default 5.
#' @param intensity_cutoff numeric. Minimum absolute intensity. Default 1.
#' @param IR_RelAb_cutoff numeric. Relative abundance cutoff (\%) for
#'   theoretical isotopologues. Default 1.
#'
#' @return A named list where each element contains the HRMF results for one
#'   compound (with a valid formula) from the MSP file.
#'
#' @seealso \code{\link{HRMF}}, \code{\link{getMSP}}
#'
#' @examples
#' \dontrun{
#' results <- getHRMF("library.msp")
#' }
#' @export
getHRMF <- function(file, charge = 1, mass_accuracy = 5,
                    intensity_cutoff = 1, IR_RelAb_cutoff = 1) {

        compounds <- getMSP(file)

        all_hrmf <- list()

        for (a in seq_along(compounds)) {
                comp <- compounds[[a]]
                formula <- comp$formula

                # Skip compounds without formula annotation
                if (is.null(formula) || length(formula) == 0 ||
                    is.na(formula) || formula == "") {
                        message("Compound ", a, " (",
                                ifelse(is.null(comp$name), "unknown", comp$name),
                                "): no formula, skipping.")
                        next
                }

                result <- tryCatch(
                        HRMF(msp = comp,
                             formula = formula,
                             charge = charge,
                             mass_accuracy = mass_accuracy,
                             intensity_cutoff = intensity_cutoff,
                             IR_RelAb_cutoff = IR_RelAb_cutoff,
                             detailed = TRUE),
                        error = function(e) {
                                message("Error processing compound ", a,
                                        " (", comp$name, "): ", e$message)
                                NULL
                        }
                )

                if (!is.null(result)) {
                        name_label <- ifelse(is.null(comp$name) || comp$name == "",
                                             paste0("compound_", a), comp$name)
                        all_hrmf[[name_label]] <- result
                }
        }

        return(all_hrmf)
}
