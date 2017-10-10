#' Get xcmsset object in one step with optimized methods.
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the references
#' @param minfrac minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group, default 0.67
#' @param ... arguments for xcmsSet function
#' @details the parameters are extracted from the papers. If you use name other than the name above, you will use the default setting of XCMS. Also I suggest IPO packages or apLCMS packages to get reasonable data for your own instrumental. If you want to summit the results to a paper, remember to include those parameters.
#' @return a xcmsset object for that path or selected samples
#' @references Patti, G. J.; Tautenhahn, R.; Siuzdak, G. Nat. Protocols 2012, 7 (3), 508–516.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' }
#' @export
getdata <-
        function(path,
                 index = F,
                 BPPARAM = BiocParallel::SnowParam(),
                 pmethod = "hplcorbitrap",
                 minfrac = 0.67,
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE,
                                       full.names = TRUE)
                if (index) {
                        cdffiles <- cdffiles[index]
                }
                if (pmethod == "hplcorbitrap") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 2.5,
                                peakwidth = c(10,
                                              60),
                                prefilter = c(3, 5000),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 5,
                                        mzwid = 0.015,
                                        minfrac = min
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 5,
                                                mzwid = 0.015,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == "uplcorbitrap") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 2.5,
                                peakwidth = c(5,
                                              20),
                                prefilter = c(3, 5000),
                                ...
                        )
                        xset <-
                                xcms::group(
                                        xset,
                                        bw = 2,
                                        mzwid = 0.015,
                                        minfrac = minfrac
                                )
                        xset2 <- xcms::retcor(xset, 'obiwarp')
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <-
                                xcms::group(
                                        xset2,
                                        bw = 2,
                                        mzwid = 0.015,
                                        minfrac = minfrac
                                )
                        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                } else if (pmethod == "hplcqtof") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 30,
                                peakwidth = c(10,
                                              60),
                                prefilter = c(0, 0),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 5,
                                        mzwid = 0.025,
                                        minfrac = minfrac
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 5,
                                                mzwid = 0.025,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == "hplchqtof") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 15,
                                peakwidth = c(10,
                                              60),
                                prefilter = c(0, 0),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 5,
                                        mzwid = 0.015,
                                        minfrac = minfrac
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 5,
                                                mzwid = 0.015,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == "uplcqtof") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 30,
                                peakwidth = c(5,
                                              20),
                                prefilter = c(0, 0),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 2,
                                        mzwid = 0.025,
                                        minfrac = minfrac
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 2,
                                                mzwid = 0.025,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == "uplchqtof") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 15,
                                peakwidth = c(5,
                                              20),
                                prefilter = c(0, 0),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 2,
                                        mzwid = 0.015,
                                        minfrac = minfrac
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 2,
                                                mzwid = 0.015,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else {
                        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                              ...)
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(xset, minfrac = minfrac)
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(xset2, minfrac = minfrac)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                }
                return(xset3)
        }

#' Get XCMSnExp object in one step from structured folder path for xcms 3.
#' @param path the path to your data
#' @param index the index of the files
#' @param snames sample names. By default the file name without extension is used
#' @param sclass sample classes.
#' @param phenoData data.frame or NAnnotatedDataFrame defining the sample names and classes and other sample related properties. If not provided, the argument sclass or the subdirectories in which the samples are stored will be used to specify sample grouping.
#' @param BPPARAM used for BiocParallel package
#' @param ppp parameters for peaks picking, e.g. xcms::CentWaveParam()
#' @param rtp parameters for retention time correction, e.g. xcms::ObiwarpParam()
#' @param gpp parameters for peaks grouping, e.g. xcms::PeakDensityParam()
#' @param fpp parameters for peaks filling, e.g. xcms::FillChromPeaksParam(), PeakGroupsParam()
#' @details This is a wrap function for metabolomics data process for xcms 3.
#' @return a XCMSnExp object with processed data
#' @export
getdata2 <- function(path,
                     index = F,
                     snames = NULL,
                     sclass = NULL,
                     phenoData = NULL,
                     BPPARAM = BiocParallel::SnowParam(),
                     ppp = xcms::CentWaveParam(
                             ppm = 5,
                             peakwidth = c(5, 25),
                             prefilter = c(3, 5000)
                     ),
                     rtp = xcms::PeakGroupsParam(minFraction = 0.67),
                     gpp = xcms::PeakDensityParam(minFraction = 0.67,
                                                  bw = 2,
                                                  binSize = 0.025),
                     fpp = xcms::FillChromPeaksParam()) {
        files <- list.files(path, recursive = TRUE,
                            full.names = TRUE)
        if (index) {
                files <- files[index]
        }

        fromPaths <- xcms::phenoDataFromPaths(files)
        if (is.null(snames)) {
                snames <- rownames(fromPaths)
        } else {
                rownames(fromPaths) <- snames
        }

        if (is.null(snames)) {
                snames <- rownames(fromPaths)
        } else {
                rownames(fromPaths) <- snames
        }
        pdata <- phenoData
        if (is.null(pdata)) {
                pdata <- sclass
                if (is.null(pdata))
                        pdata <-
                                methods::new("NAnnotatedDataFrame", fromPaths)
        } else{
                if (class(pdata) == "data.frame")
                        pdata <-
                                methods::new("NAnnotatedDataFrame", fromPaths)
                if (class(pdata) != "NAnnotatedDataFrame")
                        stop("phenoData has to be a data.frame or NAnnotatedDataFrame!")
        }
        raw_data <- MSnbase::readMSData2(files, pdata = pdata)
        xod <-
                xcms::findChromPeaks(raw_data, param = ppp, BPPARAM = BPPARAM)
        xod <- xcms::groupChromPeaks(xod, param = gpp)
        xod <- xcms::adjustRtime(xod, param = rtp)
        xod <- xcms::groupChromPeaks(xod, param = gpp)
        xod <-
                xcms::fillChromPeaks(xod, param = fpp, BPPARAM = BPPARAM)
        return(xod)
}

#' Get the csv files to be submitted to Metaboanalyst
#' @param xset the xcmsset object which you want to submitted to Metaboanalyst
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param name file name
#' @return dataframe with data needed for Metaboanalyst if your want to perform local analysis.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getupload(xset)
#' }
#' @export
getupload <- function(xset,
                      method = "medret",
                      intensity = "into",
                      name = "Peaklist") {
        peakIntensities <- xcms::groupval(xset, method,
                                          intensity)
        peakIntensities[is.na(peakIntensities)] = 0

        data <-
                rbind(group = as.character(xcms::phenoData(xset)$class),
                      peakIntensities)
        data <- data[!duplicated(rownames(data)),]
        filename <- paste0(name, ".csv")
        utils::write.csv(data, file = filename)
        return(data)
}

#' Get the mzrt profile and group information for batch correction and plot as a list
#' @param xset xcmsSet objects
#' @return list with rtmz profile and group infomation
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getmzrt(xset)
#' }

getmzrt <- function(xset){
        data <- xcms::groupval(xset, value = "into")
        idx <- stats::complete.cases(data)
        data0 <- data[idx,]
        # group info
        group <- xcms::phenoData(xset)
        # peaks info
        peaks <- as.data.frame(xcms::groups(xset))
        mz <- peaks$mzmed[idx]
        rt <- peaks$rtmed[idx]
        # return as list
        result <- list(data = data0, group = group, mz = mz, rt = rt)
        return(result)
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
    zz <- file(file.path("MSP", paste(outfilename,
        ".msp", sep = "")), "w")
    nPeaks <- length(mz)
    cat("Name: unknown", paste("Num Peaks: ", nPeaks),
        file = zz, sep = "\n")
    while (length(mz) >= 5) {
        cat(paste(mz[1:5]), "", file = zz, sep = "; ")
        cat(paste("\n"), file = zz)
        mz <- mz[6:length(mz)]
    }
    if (!is.na(mz[1])) {
        cat(paste(mz), "", file = zz, sep = "; ")
        cat(paste("\n"), file = zz)
    }
    close(zz)
    print(paste("A data file", outfilename, ".MSP has been generated in the folder:",
        "MSP", cat("\n")))
}
