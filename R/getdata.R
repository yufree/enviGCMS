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
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' }
#' @seealso \code{\link{getdata2}}, \code{\link{getmzrt}}
#' @export
getdata <-
        function(path,
                 index = F,
                 BPPARAM = BiocParallel::SnowParam(),
                 pmethod = "hplcorbitrap",
                 minfrac = 0.67,
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
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
                                xset2 <- xcms::retcor(xset, "obiwarp")
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
                        xset <- xcms::group(
                                xset,
                                bw = 2,
                                mzwid = 0.015,
                                minfrac = minfrac
                        )
                        xset2 <- xcms::retcor(xset, "obiwarp")
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <- xcms::group(
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
                                xset2 <- xcms::retcor(xset, "obiwarp")
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
                                xset2 <- xcms::retcor(xset, "obiwarp")
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
                                xset2 <- xcms::retcor(xset, "obiwarp")
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
                                xset2 <- xcms::retcor(xset, "obiwarp")
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
                                xset2 <- xcms::retcor(xset, "obiwarp")
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
#' @param mode 'inMemory' or 'onDisk' see `?MSnbase::readMSData` for details, default 'onDisk'
#' @param ppp parameters for peaks picking, e.g. xcms::CentWaveParam()
#' @param rtp parameters for retention time correction, e.g. xcms::ObiwarpParam()
#' @param gpp parameters for peaks grouping, e.g. xcms::PeakDensityParam()
#' @param fpp parameters for peaks filling, e.g. xcms::FillChromPeaksParam(), PeakGroupsParam()
#' @details This is a wrap function for metabolomics data process for xcms 3.
#' @return a XCMSnExp object with processed data
#' @seealso \code{\link{getdata}},\code{\link{getmzrt}}
#' @export
getdata2 <- function(path,
                     index = F,
                     snames = NULL,
                     sclass = NULL,
                     phenoData = NULL,
                     BPPARAM = BiocParallel::SnowParam(),
                     mode = "onDisk",
                     ppp = xcms::CentWaveParam(
                             ppm = 5,
                             peakwidth = c(5, 25),
                             prefilter = c(3, 5000)
                     ),
                     rtp = xcms::ObiwarpParam(binSize = 1),
                     gpp = xcms::PeakDensityParam(
                             sampleGroups = 1,
                             minFraction = 0.67,
                             bw = 2,
                             binSize = 0.025
                     ),
                     fpp = xcms::FillChromPeaksParam()) {
        files <- list.files(path, recursive = TRUE, full.names = TRUE)
        if (index) {
                files <- files[index]
        }

        fromPaths <- xcms::phenoDataFromPaths(files)
        n <- dim(fromPaths)[2]
        sample_group = NULL
        if (n > 1) {
                sample_group <- fromPaths[, 1]
                for (i in 2:n) {
                        sample_group <- paste(sample_group, fromPaths[,
                                                                      i], sep = "_")
                }
        } else {
                sample_group = fromPaths[, 1]
        }
        sample_group = data.frame(sample_group)

        if (is.null(snames)) {
                snames <- rownames(fromPaths)
        } else {
                rownames(sample_group) <- snames
        }

        pdata <- phenoData
        if (is.null(pdata)) {
                pdata <- sclass
                if (is.null(pdata))
                        pdata <- methods::new("NAnnotatedDataFrame",
                                              sample_group)
        } else {
                if (class(pdata) == "data.frame")
                        pdata <- methods::new("NAnnotatedDataFrame",
                                              sample_group)
                if (class(pdata) != "NAnnotatedDataFrame")
                        stop("phenoData has to be a data.frame or NAnnotatedDataFrame!")
        }
        raw_data <- MSnbase::readMSData(files, pdata = pdata,
                                        mode = mode)
        gpp@sampleGroups <- pdata$sample_group
        xod <- xcms::findChromPeaks(raw_data, param = ppp,
                                    BPPARAM = BPPARAM)
        xod <- xcms::adjustRtime(xod, param = rtp)
        xod <- xcms::groupChromPeaks(xod, param = gpp)
        xod <-
                xcms::fillChromPeaks(xod, param = fpp, BPPARAM = BPPARAM)
        return(xod)
}
