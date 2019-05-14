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
#' @return list with rtmz profile and group infomation
#' @seealso \code{\link{getmzrt}}
getmzrtcsv <- function(path) {
        dataraw <- utils::read.csv(path, skip = 1)
        mz <- dataraw[, 2]
        rt <- dataraw[, 3]
        data <- dataraw[,-c(1:3)]
        group <- data.frame(t(utils::read.csv(path, nrows = 1)[-(1:3)]))
        colnames(group) <- c(1:ncol(group))
        colnames(data) <- rownames(group)
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
#' mz <- c(10000,20000,10000,30000,5000)
#' names(mz) <- c(101,143,189,221,234)
#' writeMSP(mz,'test')
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
#' @param charge The charge value of the formula.
#' @param window The window accuracy in the same units as mass
#' @param elements Elements list to take into account.
#' @return list with chemical formula
#' @export
getformula <-
        function(mz,
                 charge = 1,
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
                        a <- mz[i]
                        mfSet <-
                                rcdk::generate.formula.iter(
                                        a,
                                        charge = charge,
                                        window = window,
                                        elements = elements,
                                        validation = T
                                )
                        hit <- itertools::ihasNext(mfSet)
                        re <- NULL
                        while (itertools::hasNext(hit)) {
                                temp <- iterators::nextElem(hit)
                                re <- c(re, temp)
                        }
                        list[[i]] <- re
                }
                aa <- lapply(list, unlist)
                getvalid <- function(x) {
                        a <- NULL
                        for (i in x) {
                                re <- rcdk::isvalid.formula(rcdk::get.formula(i))
                                a <- c(a, re)
                        }
                        return(x[a])
                }
                bb <- lapply(aa, getvalid)

                return(bb)
        }
