#' Get the mzrt profile and group information for batch correction and plot as a list directly from path with default setting
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the references
#' @param minfrac minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group, default 0.67
#' @param ... arguments for xcmsSet function
#' @return list with rtmz profile and group infomation
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' }
#' @seealso \code{\link{getdata}},\code{\link{getupload}}, \code{\link{getmzrt}}, \code{\link{getdoe}}
getmr <-
        function(path,
                 index = F,
                 BPPARAM = BiocParallel::SnowParam(),
                 pmethod = "hplcorbitrap",
                 minfrac = 0.67,
                 ...) {
                xset <- getdata(
                        path = path,
                        index = index,
                        BPPARAM = BPPARAM,
                        pmethod = pmethod,
                        minfrac = minfrac
                )
                list <- getmzrt(xset, ...)
                return(list)
        }

#' Covert the peaks list csv file into list
#' @param path the path to your csv file
#' @return list with rtmz profile and group infomation as the first row
#' @seealso \code{\link{getmzrt}}
#' @export
getmzrtcsv <- function(path) {
        dataraw <- utils::read.csv(path, skip = 1)
        mz <- dataraw[, 2]
        rt <- dataraw[, 3]
        data <- dataraw[,-c(1:3)]
        group <- c(t(utils::read.csv(path, nrows = 1)[-(1:3)]))
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
#' Write MSP files for NIST search
#' @param mz a intensity vector, who name is the mass in m/z
#' @param outfilename the name of the MSP file, default is 'unknown'
#' @return none a MSP file will be created at the subfolder working dictionary with name 'MSP'
#' @examples
#' \dontrun{
#' mz <- c(10000,20000,10000,30000,5000)
#' names(mz) <- c(101,143,189,221,234)
#' writeMSP(mz,'test')
#' }
#' @export
writeMSP <- function(mz, outfilename = "unknown") {
        mz <- paste(names(mz), round(mz))
        dir.create("MSP")
        zz <- file(file.path("MSP", paste(outfilename, ".msp",
                                          sep = "")), "w")
        nPeaks <- length(mz)
        cat(
                "Name: unknown",
                paste("Num Peaks: ", nPeaks),
                file = zz,
                sep = "\n"
        )
        while (length(mz) >= 5) {
                cat(paste(mz[1:5]),
                    "",
                    file = zz,
                    sep = "; ")
                cat(paste("\n"), file = zz)
                mz <- mz[6:length(mz)]
        }
        if (!is.na(mz[1])) {
                cat(paste(mz),
                    "",
                    file = zz,
                    sep = "; ")
                cat(paste("\n"), file = zz)
        }
        close(zz)
        print(
                paste(
                        "A data file",
                        outfilename,
                        ".MSP has been generated in the folder:",
                        "MSP",
                        cat("\n")
                )
        )
}

#' get the data of QC compound for a group of data
#' @param path data path for your QC samples
#' @param mzrange mass of the QC compound
#' @param rtrange retention time of the QC compound
#' @param index index of the files contained QC compounds, default is all of the compounds
#' @return number vector, each number indicate the peak area of that mass and retention time range
getQCraw <- function(path, mzrange, rtrange, index = NULL) {
        cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
        if (index) {
                cdffiles <- cdffiles[index]
        }
        nsamples <- length(cdffiles)
        area <- numeric()
        for (i in 1:nsamples) {
                RAW <- xcms::xcmsRaw(cdffiles[i])
                peak <- xcms::rawEIC(RAW, mzrange, rtrange)
                area[i] <- sum(peak$intensity)
        }
        return(area)
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
                for (i in 1:length(mz)) {
                        element <- paste(names(elements), sep="", collapse="")
                        minelement <- sapply(elements, function(x) x[1])
                        maxelement <- sapply(elements, function(x) x[2])
                        minelement2 <- paste(paste0(names(minelement),minelement), sep="", collapse="")
                        maxelement2 <- paste(paste0(names(maxelement),maxelement), sep="", collapse="")
                        mfSet <-
                                Rdisop::decomposeMass(mz[i],mzabs = window,z = charge, elements = element, minElements = minelement2, maxElements = maxelement2)
                        formula <- mfSet$formula
                        valid <- mfSet$valid
                        list[[i]] <- formula[valid=='Valid']

                }
                return(list)
        }

#' read in MSP file as list for ms/ms annotation
#' @param file the path to your MSP file
#' @return list a list with MSP information for MS/MS annotation
#' @export
getMSP <- function(file){
        # this part is modified from compMS2Miner's code: https://github.com/WMBEdmands/compMS2Miner/blob/ee20d3d632b11729d6bbb5b5b93cd468b097251d/R/metID.matchSpectralDB.R
        msp <- readLines(file)
        # remove empty lines
        msp <- msp[msp != '']
        ncomp <- grep('^NAME:', msp, ignore.case = TRUE)
        splitFactorTmp <- rep(1:length(ncomp), diff(c(ncomp, length(msp) + 1)))
        # matrix of masses and intensities
        massIntIndx <- which(grepl('^[0-9]', msp) & !grepl(': ', msp))
        massesInts <- unlist(strsplit(msp[massIntIndx], '\t| '))
        massesInts <- as.numeric(massesInts[grep('^[0-9].*[0-9]$|^[0-9]$', massesInts)])
        # if any NAs remove from indx
        massesTmp <-  massesInts[seq(1, length(massesInts), 2)]
        relIntTmp <-  massesInts[seq(2, length(massesInts), 2)]
        # label with entry numbers
        names(massesTmp) <- splitFactorTmp[massIntIndx]
        names(relIntTmp) <- splitFactorTmp[massIntIndx]
        # label with entry numbers
        names(massesTmp) <- splitFactorTmp[massIntIndx]
        names(relIntTmp) <- splitFactorTmp[massIntIndx]
        # entry info
        entryInfoIndx <- setdiff(1:length(msp), massIntIndx)
        entryInfo <- msp[entryInfoIndx]
        names(entryInfo) <- paste0(splitFactorTmp[entryInfoIndx], '_', gsub(':.+', '', entryInfo))
        entryInfo <- gsub('.+: ', '', entryInfo, ignore.case = TRUE)
        # remove duplicate entries
        entryInfo <- entryInfo[duplicated(names(entryInfo)) == FALSE]
        # identify pertinent msp file entries
        mandEntryData <- data.frame(matrix('', nrow=length(ncomp), ncol=10), stringsAsFactors = FALSE)
        colnames(mandEntryData) <- c('DBname', 'DBIonMode', 'DBSMILES', 'DBINCHI',
                                     'DBprecursorMasses', 'DBesiType', 'FORMULA','RETENTIONTIME','COLLISIONENERGY','NumPeaks')
        # name of entry
        mandFieldsRegEx <- c('_NAME$', '_ION MODE$|_MODE$|_IONMODE', '_SMILES$', '_INCHI$|_INCHIKEY$',
                             '_PRECURSORMZ$|_PRECURSOR M/Z$|_PRECURSOR MZ$|_PEPMASS$',
                             '_PRECURSORTYPE$|_PRECURSOR TYPE$|_ADDUCT$|_ION TYPE$|_IONTYPE$',
                             '_FORMULA$','_RETENTIONTIME$','_COLLISIONENERGY$','_Num Peaks$')
        names(mandFieldsRegEx) <- colnames(mandEntryData)
        for(mF in 1:length(mandFieldsRegEx)){
                entInfoTmp <- entryInfo[grep(mandFieldsRegEx[mF], names(entryInfo), ignore.case = TRUE)]
                entryNosTmp <- gsub(mandFieldsRegEx[mF], '', names(entInfoTmp), ignore.case=TRUE)
                indxTmp <- row.names(mandEntryData) %in% entryNosTmp
                mandEntryData[indxTmp, names(mandFieldsRegEx[mF])] <- entInfoTmp
        }
        mandEntryData$DBprecursorMasses <- suppressWarnings(as.numeric(mandEntryData$DBprecursorMasses))
        mandEntryData$NumPeaks <- suppressWarnings(as.numeric(mandEntryData$NumPeaks))
        missingVals <- mandEntryData$NumPeaks!=0
        mandEntryData <- mandEntryData[missingVals, , drop=FALSE]
        mandEntryData$ID <- rownames(mandEntryData)
        getlist <- function(x) {
                mz <- massesTmp[names(massesTmp) == x['ID']]
                ins <- relIntTmp[names(relIntTmp) == x['ID']]
                msms <- cbind.data.frame(mz=mz,ins=ins)
                return(list(msms=msms,name=x['DBname'],ionmode = as.character(x['DBIonMode']),smiles = as.character(x['DBSMILES']), inchi = as.character(x['DBINCHI']), prec = suppressWarnings(as.numeric(x['DBprecursorMasses'])), DBesiType = as.character(x['DBesiType']), formula = as.character(x['FORMULA']), rt=as.character(x['RETENTIONTIME']), ce=as.character(x['COLLISIONENERGY']), npeaks=suppressWarnings(as.numeric(x['NumPeaks']))))
        }
        li <- apply(mandEntryData,1,getlist)
        return(li)
}






