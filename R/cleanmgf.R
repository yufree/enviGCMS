#' Parse an MGF file into a list of spectra blocks
#'
#' Each element contains the raw header lines (KEY=VALUE) and the parsed
#' peak table with the original peak line strings preserved.
#'
#' @param file character. Path to an MGF file.
#' @return A list of lists with elements \code{headers} (character vector)
#'   and \code{peaks} (data.frame with columns mz, intensity, line).
#' @noRd
.parse_mgf <- function(file) {
        lines <- readLines(file, warn = FALSE)
        begin <- grep("^BEGIN IONS", lines, ignore.case = TRUE)
        end <- grep("^END IONS", lines, ignore.case = TRUE)
        if (length(begin) == 0) {
                stop("No BEGIN IONS blocks found; is this an MGF file?",
                     call. = FALSE)
        }
        if (length(begin) != length(end)) {
                stop("Unbalanced BEGIN/END IONS blocks in MGF file.",
                     call. = FALSE)
        }

        lapply(seq_along(begin), function(i) {
                block <- trimws(lines[(begin[i] + 1):(end[i] - 1)])
                block <- block[block != "" & !grepl("^#", block)]
                hdr_idx <- grepl("=", block)
                headers <- block[hdr_idx]
                peak_lines <- block[!hdr_idx]

                peaks <- data.frame(mz = numeric(0), intensity = numeric(0),
                                    line = character(0),
                                    stringsAsFactors = FALSE)
                if (length(peak_lines) > 0) {
                        toks <- strsplit(peak_lines, "[ \t]+")
                        keep <- vapply(toks, function(t) {
                                length(t) >= 2 &&
                                        !is.na(suppressWarnings(as.numeric(t[1]))) &&
                                        !is.na(suppressWarnings(as.numeric(t[2])))
                        }, logical(1))
                        if (any(keep)) {
                                m <- vapply(toks[keep], function(t) {
                                        as.numeric(t[1:2])
                                }, numeric(2))
                                peaks <- data.frame(
                                        mz = m[1, ], intensity = m[2, ],
                                        line = unlist(peak_lines[keep]),
                                        stringsAsFactors = FALSE
                                )
                        }
                }
                list(headers = headers, peaks = peaks)
        })
}

#' Extract a header value from an MGF block
#'
#' @param headers character vector of KEY=VALUE lines.
#' @param key character, header key (case-insensitive).
#' @return character value or NA when absent.
#' @noRd
.mgf_header <- function(headers, key) {
        hit <- grep(paste0("^", key, "="), headers, ignore.case = TRUE)
        if (length(hit) == 0) return(NA_character_)
        sub(paste0("^", key, "="), "", headers[hit[1]], ignore.case = TRUE)
}

#' Parse an MGF CHARGE header value such as "1+", "2+" or "+2"
#'
#' @param x character, raw header value.
#' @param default integer, fallback when missing or unparsable.
#' @return integer charge (negative for negative-mode values).
#' @noRd
.mgf_charge <- function(x, default = 1L) {
        if (is.na(x)) return(default)
        neg <- grepl("-", x)
        num <- suppressWarnings(as.numeric(gsub("[+-]", "", x)))
        if (is.na(num) || num == 0) return(default)
        as.integer(if (neg) -num else num)
}

#' Flag fragment peaks explainable as sub-formulae of a parent formula
#'
#' Uses the same ion-m/z conversion as \code{\link{HRMF}}: query =
#' fragment m/z - (adduct delta + electron mass) at the radical-ion
#' convention.
#'
#' @param mz numeric vector of fragment m/z.
#' @param parent_formula character, candidate parent formula.
#' @param delta numeric, adduct mass delta (ion m/z = neutral + delta).
#' @param charge integer, fragment charge sign (1 or -1).
#' @param mass_accuracy numeric, matching window in ppm.
#' @return logical vector, TRUE where a sub-formula was found.
#' @noRd
.explain_peaks <- function(mz, parent_formula, delta, charge, mass_accuracy) {
        atoms <- .parse_formula(parent_formula)
        element_str <- paste(atoms$element, collapse = "")
        min_str <- paste0(paste0(atoms$element, "0"), collapse = "")
        max_str <- paste0(paste0(atoms$element, atoms$count), collapse = "")
        me <- if (charge > 0) .ELECTRON_MASS else -.ELECTRON_MASS
        offset <- delta + me
        qmz <- mz - offset

        vapply(qmz, function(m) {
                mf <- tryCatch(
                        .decomposeMass(
                                m,
                                mzabs = mass_accuracy / 1e6 * max(abs(m), 1),
                                z = charge,
                                elements = element_str,
                                minElements = min_str,
                                maxElements = max_str
                        ),
                        error = function(e) NULL
                )
                !is.null(mf) && length(mf$formula) > 0
        }, logical(1))
}

#' Clean an MGF file by keeping only formula-explainable peaks
#'
#' For each MS2 spectrum, the precursor neutral mass is derived from
#' \code{PEPMASS}, \code{CHARGE} and the supplied \code{adduct}, candidate
#' parent formulae are enumerated with the native HORIZON engine, and fragment peaks that
#' can be explained as sub-formulae of the best candidate (the one explaining
#' the most peaks) are kept. Spectra without \code{PEPMASS} or without any
#' candidate formula are written back unchanged.
#'
#' @param file character. Path to an MGF file.
#' @param out_file character. Output MGF path. Default NULL writes
#'   \code{<file stem>_clean.mgf} next to the input file.
#' @param adduct character. Single adduct name (see \code{\link{HRMF}}),
#'   e.g. "[M+H]+". Fragments are assumed singly-charged with this adduct.
#' @param charge integer. Precursor charge state used when the \code{CHARGE}
#'   header is missing. Default NULL uses 1 (from the adduct sign).
#' @param mass_accuracy numeric. Mass accuracy in ppm for both precursor
#'   decomposition and fragment matching. Default 5.
#' @param elements character. Elements allowed in precursor formula
#'   enumeration, e.g. "CHNOPS". Default "CHNOPS".
#' @param max_candidates integer. Maximum number of precursor formula
#'   candidates to evaluate per spectrum. Default 3.
#' @param intensity_cutoff numeric. Discard peaks below this absolute
#'   intensity before annotation. Default 0 (keep all).
#' @param min_peaks integer. Drop spectra left with fewer than this many
#'   peaks after cleaning. Default 2.
#'
#' @return Invisibly, a data.frame with one row per spectrum: title, pepmass,
#'   charge, chosen formula, peak counts before/after, kept fraction and
#'   status ("cleaned", "kept" = unchanged, or "dropped").
#'
#' @details The cleaned file preserves the original header lines and peak
#'   line formatting; the chosen parent formula is appended as a
#'   \code{FORMULA=} header when a candidate is found.
#'
#' @seealso \code{\link{HRMF}}, \code{\link{getMSP}}
#'
#' @examples
#' \dontrun{
#' summary <- cleanMGF("ms2.mgf", adduct = "[M+H]+")
#' summary[, c("title", "formula", "n_before", "n_kept", "status")]
#' }
#' @export
cleanMGF <- function(file, out_file = NULL, adduct = "[M+H]+",
                     charge = NULL, mass_accuracy = 5,
                     elements = "CHNOPS", max_candidates = 3,
                     intensity_cutoff = 0, min_peaks = 2) {

        if (length(adduct) != 1) {
                stop("'adduct' must be a single adduct name; use HRMF() to ",
                     "score multiple adducts.", call. = FALSE)
        }
        modes <- .resolve_adducts(adduct)
        mode <- modes[1, ]
        if (is.null(charge)) charge <- mode$charge

        spectra <- .parse_mgf(file)
        out_lines <- character(0)
        report <- vector("list", length(spectra))

        for (i in seq_along(spectra)) {
                sp <- spectra[[i]]
                title <- utils::head(.mgf_header(sp$headers, "TITLE"), 1)
                if (is.na(title)) title <- paste0("spectrum_", i)

                pepmass_raw <- .mgf_header(sp$headers, "PEPMASS")
                prec <- if (is.na(pepmass_raw)) NA_real_ else
                        suppressWarnings(as.numeric(
                                strsplit(pepmass_raw, "[ \t]+")[[1]][1]
                        ))

                rep <- list(title = title, pepmass = prec,
                            charge = NA_integer_, formula = NA_character_,
                            n_before = nrow(sp$peaks), n_kept = nrow(sp$peaks),
                            kept_fraction = 1, status = "kept")

                # Spectra without usable PEPMASS or peaks pass through unchanged
                if (is.na(prec) || nrow(sp$peaks) == 0) {
                        out_lines <- c(out_lines, "BEGIN IONS",
                                       sp$headers, sp$peaks$line, "END IONS", "")
                        rep$status <- if (is.na(prec)) "kept (no pepmass)" else "kept"
                        report[[i]] <- rep
                        next
                }

                z <- .mgf_charge(.mgf_header(sp$headers, "CHARGE"),
                                 default = as.integer(charge))
                rep$charge <- z

                peaks <- sp$peaks
                if (intensity_cutoff > 0) {
                        peaks <- peaks[peaks$intensity >= intensity_cutoff, , drop = FALSE]
                }

                # Neutral precursor mass: [M + n*adduct]^n+ -> M = m/z*n - n*delta
                neutral <- prec * z - mode$delta * z

                cand <- tryCatch(
                        .decomposeMass(
                                neutral,
                                mzabs = mass_accuracy / 1e6 * neutral,
                                z = 0,
                                elements = elements
                        ),
                        error = function(e) NULL
                )

                if (is.null(cand) || length(cand$formula) == 0) {
                        out_lines <- c(out_lines, "BEGIN IONS",
                                       sp$headers, sp$peaks$line, "END IONS", "")
                        rep$status <- "kept (no formula candidates)"
                        report[[i]] <- rep
                        next
                }

                cand_idx <- which(cand$valid == "Valid")
                if (length(cand_idx) == 0) cand_idx <- seq_along(cand$formula)
                # Drop chemically absurd candidates (e.g. strongly negative
                # double-bond equivalents); keep all as fallback
                sane_idx <- cand_idx[!is.na(cand$DBE[cand_idx]) &
                                              cand$DBE[cand_idx] >= -0.5]
                if (length(sane_idx) > 0) cand_idx <- sane_idx
                cand_formulas <- cand$formula[cand_idx]

                best_mask <- rep(FALSE, nrow(peaks))
                best_formula <- NA_character_
                best_frac <- -1
                for (ci in seq_len(min(max_candidates, length(cand_formulas)))) {
                        mask <- .explain_peaks(peaks$mz, cand_formulas[ci],
                                               mode$delta, mode$charge,
                                               mass_accuracy)
                        frac <- if (length(mask)) mean(mask) else 0
                        if (frac > best_frac) {
                                best_frac <- frac
                                best_mask <- mask
                                best_formula <- cand_formulas[ci]
                        }
                }

                kept_peaks <- peaks[best_mask, , drop = FALSE]
                rep$formula <- best_formula
                rep$n_before <- nrow(sp$peaks)
                rep$n_kept <- nrow(kept_peaks)
                rep$kept_fraction <- if (nrow(sp$peaks) > 0) {
                        round(rep$n_kept / nrow(sp$peaks), 3)
                } else 1

                if (nrow(kept_peaks) < min_peaks) {
                        rep$status <- "dropped"
                        report[[i]] <- rep
                        next
                }

                new_headers <- sp$headers
                if (!any(grepl("^FORMULA=", new_headers, ignore.case = TRUE)) &&
                    !is.na(best_formula)) {
                        new_headers <- c(new_headers,
                                         paste0("FORMULA=", best_formula))
                }
                out_lines <- c(out_lines, "BEGIN IONS", new_headers,
                               kept_peaks$line, "END IONS", "")
                rep$status <- "cleaned"
                report[[i]] <- rep
        }

        if (is.null(out_file)) {
                stem <- sub("\\.[Mm][Gg][Ff]$", "", file)
                out_file <- paste0(stem, "_clean.mgf")
        }
        con <- file(out_file, open = "wt", encoding = "UTF-8")
        writeLines(out_lines, con = con)
        close(con)

        res <- do.call(rbind, lapply(report, function(r) {
                data.frame(title = r$title, pepmass = r$pepmass,
                           charge = r$charge, formula = r$formula,
                           n_before = r$n_before, n_kept = r$n_kept,
                           kept_fraction = r$kept_fraction,
                           status = r$status, stringsAsFactors = FALSE)
        }))
        rownames(res) <- NULL
        invisible(res)
}
