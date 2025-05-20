#' Covert the peaks list csv file into list
#' @param path the path to your csv file
#' @return list with rtmz profile and group information as the first row
#' @seealso \code{\link{getmzrt}}
#' @export
getmzrtcsv <- function(path) {
        dataraw <- utils::read.csv(path, skip = 1)
        sample_name <-
                names(utils::read.csv(path, nrows = 1)[-(1:3)])
        mz <- dataraw[, 2]
        rt <- dataraw[, 3]
        data <- dataraw[, -c(1:3)]
        colnames(data) <- sample_name
        sample_group <-
                c(t(utils::read.csv(path, nrows = 1)[-(1:3)]))
        group <-
                cbind.data.frame(sample_name, sample_group, stringsAsFactors = FALSE)
        rownames(data) <- dataraw[, 1]
        re <- list(
                data = data,
                mz = mz,
                group = group,
                rt = rt
        )
        class(re) <- "mzrt"
        return(re)
}
#' Write MSP file for NIST search
#' @param list a list with spectra information
#' @param name name of the compounds
#' @param sep numeric or logical the numbers of spectra in each file and FALSE to include all of the spectra in one msp file
#' @return none a MSP file will be created.
#' @examples
#' \dontrun{
#' ins <- c(10000,20000,10000,30000,5000)
#' mz <- c(101,143,189,221,234)
#' writeMSP(list(list(spectra = cbind.data.frame(mz,ins))), name = 'test')
#' }
#' @export
writeMSP <- function(list, name = 'unknown', sep = FALSE) {
        writemsp <- function(list) {
                mz <- paste(list$spectra$mz, list$spectra$ins)
                nPeaks <- length(mz)
                cat(
                        "BEGIN IONS",
                        paste("Name:", list$name),
                        paste("RetentionIndex:", list$rti),
                        paste("Formula:", list$formula),
                        paste("IonMode:", list$ionmode),
                        paste("CHARGE:", list$charge),
                        paste("PrecursorMz:", list$prec),
                        paste("Collision_energy:", list$ce),
                        paste("Num Peaks:", nPeaks),
                        paste("Instrument_type:", list$instr),
                        paste("Spectrum_type:", list$msm),
                        paste(mz),
                        "END IONS",
                        file = zz,
                        sep = "\n"
                )
        }
        if (sep) {
                for (i in 1:floor(length(list) / sep)) {
                        zz <-
                                file(file.path(paste(
                                        name, i, ".msp",                                            sep = ""
                                )), "w")
                        idx <- c(1:sep) + sep * (i - 1)
                        sapply(list[idx], writemsp)
                        close(zz)
                }
                zz <-
                        file(file.path(paste(
                                name, ceiling(length(list) / sep), ".msp",                                            sep = ""
                        )), "w")
                idx <- sep * floor(length(list) / sep) + 1:length(list)
                sapply(list[idx], writemsp)
                close(zz)
                message(paste0("MSP files have been generated."))
        } else{
                zz <-
                        file(file.path(paste(name, ".msp",                                            sep = "")), "w")
                sapply(list, writemsp)
                close(zz)
                message(paste0("A data file ",
                               name,
                               ".MSP has been generated."))
        }
}

#' read in MSP file as list for ms/ms or ms(EI) annotation
#' @param file the path to your MSP file
#' @return list a list with MSP information for annotation
#' @export
getMSP <- function(file) {
  msp <- readLines(file, warn = FALSE)
  msp <- msp[msp != ""]
  
  ncomp <- grep('^BEGIN IONS', msp, ignore.case = TRUE)
  if (length(ncomp) == 0) {
    ncomp <- grep("^Name", msp, ignore.case = TRUE)
  }
  splitFactorTmp <- rep(seq_along(ncomp), diff(c(ncomp, length(msp) + 1)))
  li <- split(msp, f = splitFactorTmp)
  
  # Precompile regex patterns for speed
  patterns <- list(
    name = '^NAME: |^TITLE=',
    charge = '^CHARGE=',
    ionmode = '^ION MODE:|^MODE:|^IONMODE:|^Ion_mode:',
    prec = '^PRECURSORMZ: |^PRECURSOR M/Z: |^PRECURSOR MZ: |^PEPMASS: |^PrecursorMZ: |^PEPMASS=',
    formula = '^FORMULA: |^Formula: ',
    exactmass = '^ExactMass: ',
    inchikey = '^InChIKey: ',
    np = '^Num Peaks: ',
    ce = 'COLLISIONENERGY: |Collision_energy: ',
    rt = 'RETENTIONINDEX: |RTINSECONDS: |RTINSECONDS=|retention time\\s*=\\s*[^\\\"]+',
    column = 'column\\s*=\\s*[^\\\"]+',
    instr = 'Instrument_type: ',
    msm = 'Spectrum_type: '
  )
  process_single_entry <- function(entry_lines, patterns) {
    extracted_values <- character(length(patterns))
    names(extracted_values) <- names(patterns)
    for (nm in names(patterns)) {
      pattern_to_search <- patterns[[nm]]
      idx_matches <- grep(pattern_to_search,
                          entry_lines,
                          ignore.case = TRUE,
                          perl = TRUE)
      
      if (!length(idx_matches)) {
        next
      }
      
      line_content <- entry_lines[idx_matches[1]]
      val <- NA_character_
      
      if (nm == "rt") {
        rt_match_comment <- regexpr(
          'retention time\\s*=\\s*[^\\\"]+',
          line_content,
          ignore.case = TRUE,
          perl = TRUE
        )
        if (grepl(
          "^(RETENTIONINDEX:|RTINSECONDS:|RTINSECONDS=)",
          line_content,
          ignore.case = TRUE
        )) {
          temp_val <- gsub(
            "^(RETENTIONINDEX:|RTINSECONDS:|RTINSECONDS=)\\s*",
            "",
            line_content,
            ignore.case = TRUE
          )
          val <- temp_val
        } else if (rt_match_comment != -1) {
          matched_text_list <- regmatches(line_content, rt_match_comment)
          actual_match <- matched_text_list[[1]]
          temp_val <- gsub("^retention time\\s*=",
                           "",
                           actual_match,
                           ignore.case = TRUE)
          x <- sapply(strsplit(temp_val, '\\ '), function (x)
            x[2])
          v <- sapply(strsplit(temp_val, '\\ '), function (x)
            as.numeric(x[1]))
          if (!is.na(x)&(x == 'min' | x == 'minute')) {
            val <- v * 60
          } else if (!is.na(x)&(x == 's' | x == 'sec')) {
            val <- v
          } else{
            val <- NA
          }
        }
      } else if (nm == "column") {
        col_match_obj <- regexpr(pattern_to_search,
                                 line_content,
                                 ignore.case = TRUE,
                                 perl = TRUE)
        if (col_match_obj != -1) {
          matched_text_list <- regmatches(line_content, col_match_obj)
          if (length(matched_text_list) > 0 &&
              length(matched_text_list[[1]]) > 0) {
            actual_match <- matched_text_list[[1]]
            val <- sub('^column\\s*=\\s*',
                       '',
                       actual_match,
                       ignore.case = TRUE)
          }
        }
      } else {
        temp_val <- gsub(pattern_to_search, '', line_content, ignore.case = TRUE)
        if (nm == "prec" && !is.na(temp_val)) {
          val <- strsplit(trimws(temp_val), "[ \t]+")[[1]][1]
        } else {
          val <- temp_val
        }
      }
      if (!is.null(val) && length(val) > 0 && !is.na(val)) {
        extracted_values[nm] <- trimws(val)
      }
    }
    final_fields <- as.list(extracted_values)
    
    # Parse specific fields to numeric after initial character extraction
    final_fields$prec       <- suppressWarnings(as.numeric(final_fields$prec))
    final_fields$exactmass  <- suppressWarnings(as.numeric(final_fields$exactmass))
    final_fields$rt  <- suppressWarnings(as.numeric(final_fields$rt))
    np_val                  <- suppressWarnings(as.numeric(final_fields$np))
    
    # Get masses and intensities
    massIntIndx <- which(grepl('^[0-9]', entry_lines) &
                           !grepl(': ', entry_lines))
    
    process_peaks_flag <- (length(massIntIndx) > 0) &&
      ((!is.na(np_val) &&
          np_val > 0) ||
         is.na(np_val)) # Process if NumPeaks > 0 or if NumPeaks is unknown
    
    if (process_peaks_flag) {
      peak_lines <- entry_lines[massIntIndx]
      # Robustly split by space, tab, or semicolon, and handle potential annotations
      massesInts_str <- unlist(strsplit(peak_lines, '[ \t;]+'))
      massesInts <- suppressWarnings(as.numeric(massesInts_str))
      massesInts <- massesInts[!is.na(massesInts)] # Remove NAs from non-numeric parts
      
      if (length(massesInts) > 0 &&
          (length(massesInts) %% 2 == 0)) {
        # Ensure pairs and non-empty
        mz <- massesInts[seq(1, length(massesInts), by = 2)]
        ins <- massesInts[seq(2, length(massesInts), by = 2)]
        
        if (length(ins) > 0 &&
            any(ins > 0, na.rm = TRUE) &&
            (max(ins, na.rm = TRUE) > 0)) {
          # Check max(ins) > 0 before division
          ins <- ins / max(ins, na.rm = TRUE) * 100
        } else if (length(ins) > 0) {
          # All intensities are zero, NA, or max is 0
          ins <- rep(0, length(ins))
        } else {
          mz <- numeric(0) # Ensure consistency if ins is empty
        }
        final_fields$spectra <- data.frame(mz = mz, ins = ins)
      } else {
        # If parsing m/z and intensity fails (e.g. odd number of values)
        final_fields$spectra <- data.frame(mz = numeric(0), ins = numeric(0))
      }
    } else {
      final_fields$spectra <- data.frame(mz = numeric(0), ins = numeric(0)) # Empty spectra if no peaks
    }
    return(final_fields)
  }
  li_processed <- BiocParallel::bplapply(
    li,
    FUN = process_single_entry,
    patterns = patterns,
    BPPARAM = BiocParallel::bpparam()
  )
  return(li_processed)
}

#' Get chemical formula for mass to charge ratio.
#' @param mz a vector with mass to charge ratio
#' @param charge The charge value of the formula, default 0 for autodetect
#' @param window The window accuracy in the same units as mass
#' @param elements Elements list to take into account.
#' @return list with chemical formula
#' @export
getformula <-
        function(mz,
                 charge = 0,
                 window = 0.001,
                 elements = list(
                         C = c(1, 50),
                         H = c(1, 50),
                         N = c(0, 50),
                         O = c(0, 50),
                         P = c(0, 1),
                         S = c(0, 1)
                 )) {
                list <- list()
                for (i in seq_along(mz)) {
                        element <- paste(names(elements),
                                         sep = "",
                                         collapse = "")
                        minelement <-
                                vapply(elements, function(x)
                                        x[1], 1)
                        maxelement <-
                                vapply(elements, function(x)
                                        x[2], 1)
                        minelement2 <-
                                paste(paste0(names(minelement), minelement),
                                      sep = "",
                                      collapse = "")
                        maxelement2 <-
                                paste(paste0(names(maxelement), maxelement),
                                      sep = "",
                                      collapse = "")
                        mfSet <-
                                Rdisop::decomposeMass(
                                        mz[i],
                                        mzabs = window,
                                        z = charge,
                                        elements = element,
                                        minElements = minelement2,
                                        maxElements = maxelement2
                                )
                        formula <- mfSet$formula
                        valid <- mfSet$valid
                        list[[i]] <- formula[valid == 'Valid']

                }
                return(list)
        }
