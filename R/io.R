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
    if(length(ncomp) == 0) {
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
        inchikey = '^InChIKey: ',
        np = '^Num Peaks: ',
        ce = 'COLLISIONENERGY: |Collision_energy: ',
        rt = 'RETENTIONINDEX: |RTINSECONDS: |RTINSECONDS=|retention time=',
        column = 'column=',
        instr = 'Instrument_type: ',
        msm = 'Spectrum_type: '
    )
    
    getmsp <- function(x) {
        # Extract all fields in one pass
        fields <- vapply(names(patterns), function(nm) {
            idx <- grep(patterns[[nm]], x, ignore.case = TRUE)
            if(length(idx)) {
                gsub(patterns[[nm]], '', x[idx[1]], ignore.case = TRUE)
            } else {
                NA_character_
            }
        }, character(1))
        
        # Parse numeric fields
        fields["prec"] <- as.numeric(fields["prec"])
        fields["rt"] <- as.numeric(fields["rt"])
        np_val <- as.numeric(fields["np"])
        
        # Get masses and intensities efficiently
        massIntIndx <- which(grepl('^[0-9]', x) & !grepl(': ', x))
        if(length(massIntIndx) && (!is.na(np_val) && np_val > 0 || all(is.na(fields["np"])))) {
            massesInts <- as.numeric(unlist(strsplit(x[massIntIndx], '[ \t]+')))
            mz <- massesInts[seq(1, length(massesInts), 2)]
            ins <- massesInts[seq(2, length(massesInts), 2)]
            ins <- ins / max(ins) * 100
            spectra <- data.frame(mz = mz, ins = ins)
            fields <- c(fields, list(spectra = spectra))
        }
        return(fields)
    }
    
    li <- lapply(li, getmsp)
    return(li)
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
