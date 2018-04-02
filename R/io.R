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
#' @seealso \code{\link{getdata2}},\code{\link{getupload}}, \code{\link{getmzrt}}
#' @export
getdata <- function(path, index = F, BPPARAM = BiocParallel::SnowParam(),
    pmethod = "hplcorbitrap", minfrac = 0.67, ...) {
    cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
    if (index) {
        cdffiles <- cdffiles[index]
    }
    if (pmethod == "hplcorbitrap") {
        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
            method = "centWave", ppm = 2.5, peakwidth = c(10,
                60), prefilter = c(3, 5000), ...)
        if (index & length(index) == 1) {
            xset3 <- xset
        } else {
            xset <- xcms::group(xset, bw = 5, mzwid = 0.015,
                minfrac = min)
            xset2 <- xcms::retcor(xset, "obiwarp")
            # you need group the peaks again for this corrected
            # data
            xset2 <- xcms::group(xset2, bw = 5, mzwid = 0.015,
                minfrac = minfrac)
            xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
        }
    } else if (pmethod == "uplcorbitrap") {
        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
            method = "centWave", ppm = 2.5, peakwidth = c(5,
                20), prefilter = c(3, 5000), ...)
        xset <- xcms::group(xset, bw = 2, mzwid = 0.015,
            minfrac = minfrac)
        xset2 <- xcms::retcor(xset, "obiwarp")
        # you need group the peaks again for this corrected
        # data
        xset2 <- xcms::group(xset2, bw = 2, mzwid = 0.015,
            minfrac = minfrac)
        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
    } else if (pmethod == "hplcqtof") {
        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
            method = "centWave", ppm = 30, peakwidth = c(10,
                60), prefilter = c(0, 0), ...)
        if (index & length(index) == 1) {
            xset3 <- xset
        } else {
            xset <- xcms::group(xset, bw = 5, mzwid = 0.025,
                minfrac = minfrac)
            xset2 <- xcms::retcor(xset, "obiwarp")
            # you need group the peaks again for this corrected
            # data
            xset2 <- xcms::group(xset2, bw = 5, mzwid = 0.025,
                minfrac = minfrac)
            xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
        }
    } else if (pmethod == "hplchqtof") {
        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
            method = "centWave", ppm = 15, peakwidth = c(10,
                60), prefilter = c(0, 0), ...)
        if (index & length(index) == 1) {
            xset3 <- xset
        } else {
            xset <- xcms::group(xset, bw = 5, mzwid = 0.015,
                minfrac = minfrac)
            xset2 <- xcms::retcor(xset, "obiwarp")
            # you need group the peaks again for this corrected
            # data
            xset2 <- xcms::group(xset2, bw = 5, mzwid = 0.015,
                minfrac = minfrac)
            xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
        }
    } else if (pmethod == "uplcqtof") {
        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
            method = "centWave", ppm = 30, peakwidth = c(5,
                20), prefilter = c(0, 0), ...)
        if (index & length(index) == 1) {
            xset3 <- xset
        } else {
            xset <- xcms::group(xset, bw = 2, mzwid = 0.025,
                minfrac = minfrac)
            xset2 <- xcms::retcor(xset, "obiwarp")
            # you need group the peaks again for this corrected
            # data
            xset2 <- xcms::group(xset2, bw = 2, mzwid = 0.025,
                minfrac = minfrac)
            xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
        }
    } else if (pmethod == "uplchqtof") {
        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
            method = "centWave", ppm = 15, peakwidth = c(5,
                20), prefilter = c(0, 0), ...)
        if (index & length(index) == 1) {
            xset3 <- xset
        } else {
            xset <- xcms::group(xset, bw = 2, mzwid = 0.015,
                minfrac = minfrac)
            xset2 <- xcms::retcor(xset, "obiwarp")
            # you need group the peaks again for this corrected
            # data
            xset2 <- xcms::group(xset2, bw = 2, mzwid = 0.015,
                minfrac = minfrac)
            xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
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
            xset2 <- xcms::group(xset2, minfrac = minfrac)
            xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
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
#' @seealso \code{\link{getdata}},\code{\link{getupload2}}, \code{\link{getmzrt2}}
#' @export
getdata2 <- function(path, index = F, snames = NULL, sclass = NULL,
    phenoData = NULL, BPPARAM = BiocParallel::SnowParam(),
    mode = "onDisk", ppp = xcms::CentWaveParam(ppm = 5,
        peakwidth = c(5, 25), prefilter = c(3, 5000)),
    rtp = xcms::PeakGroupsParam(minFraction = 0.67), gpp = xcms::PeakDensityParam(sampleGroups = 1,
        minFraction = 0.67, bw = 2, binSize = 0.025), fpp = xcms::FillChromPeaksParam()) {
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
    xod <- xcms::groupChromPeaks(xod, param = gpp)
    xod <- xcms::adjustRtime(xod, param = rtp)
    xod <- xcms::groupChromPeaks(xod, param = gpp)
    xod <- xcms::fillChromPeaks(xod, param = fpp, BPPARAM = BPPARAM)
    return(xod)
}

#' Get the csv files to be submitted to Metaboanalyst
#' @param xset the xcmsset object which you want to submitted to Metaboanalyst
#' @param method parameter for groupval function
#' @param value parameter for groupval function
#' @param name file name
#' @return dataframe with data needed for Metaboanalyst if your want to perform local analysis.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getupload(xset)
#' }
#' @seealso \code{\link{getdata}},\code{\link{getupload2}}, \code{\link{getmzrt}}
#' @export
getupload <- function(xset, method = "medret", value = "into",
    name = "Peaklist") {
    peakIntensities <- xcms::groupval(xset, method = method,
        value = value)

    data <- rbind(group = as.character(xcms::phenoData(xset)$class),
        peakIntensities)
    # peaks info
    peaks <- as.data.frame(xcms::groups(xset))
    mz <- peaks$mzmed
    rt <- peaks$rtmed
    rownames(data) <- c("group", paste0(round(mz, 4), "/",
        round(rt, 4)))
    filename <- paste0(name, ".csv")
    utils::write.csv(data, file = filename)
    return(data)
}

#' Get the csv files to be submitted to Metaboanalyst
#' @param xset a XCMSnExp object with processed data which you want to submitted to Metaboanalyst
#' @param value value for `xcms::featureValues`
#' @param name file name
#' @return dataframe with data needed for Metaboanalyst if your want to perform local analysis.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata2(cdfpath,
#' ppp = xcms::MatchedFilterParam(),
#' rtp = xcms::ObiwarpParam(),
#' gpp = xcms::PeakDensityParam())
#' getupload2(xset)
#' }
#' @seealso \code{\link{getdata2}},\code{\link{getupload}}, \code{\link{getmzrt2}}
#' @export
getupload2 <- function(xset, value = "into", name = "Peaklist") {
    data <- xcms::featureValues(xset, value = value)
    data <- rbind(group = as.character(xset@phenoData@data),
        data)
    # peaks info
    peaks <- xcms::featureDefinitions(xset)
    mz <- peaks$mzmed
    rt <- peaks$rtmed
    rownames(data) <- c("group", paste0(round(mz, 4), "/",
        round(rt, 4)))
    filename <- paste0(name, ".csv")
    utils::write.csv(data, file = filename)
    return(data)
}

#' Get the mzrt profile and group information for batch correction and plot as a list
#' @param xset xcmsSet objects
#' @param name file name for csv file, default NULL
#' @return list with rtmz profile and group infomation
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getmzrt(xset)
#' }
#' @seealso \code{\link{getdata}},\code{\link{getupload}}, \code{\link{getmzrt2}}, \code{\link{getdoe}},\code{\link{getmzrt}}
#' @export
getmzrt <- function(xset, name = NULL) {
    data <- xcms::groupval(xset, value = "into")
    # group info
    group <- xcms::phenoData(xset)
    # peaks info
    peaks <- as.data.frame(xcms::groups(xset))
    mz <- peaks$mzmed
    rt <- peaks$rtmed
    # return as list
    result <- list(data = data, group = group, mz = mz,
        rt = rt)
    if (!is.null(name)) {
            data <- cbind(mz = list$mz, rt = list$rt, list$data)
            data <- t(cbind(group = t(cbind(mz = 'mz',rt = 'rt',t(list$group))), t(data)))
        utils::write.csv(data, file = paste0(name, ".csv"))
    }
    return(result)
}

#' Get the mzrt profile and group information for batch correction and plot as a list for xcms 3 object
#' @param xset a XCMSnExp object with processed data
#' @param name file name for csv file, default NULL
#' @return list with rtmz profile and group infomation
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata2(cdfpath,
#' ppp = xcms::MatchedFilterParam(),
#' rtp = xcms::ObiwarpParam(),
#' gpp = xcms::PeakDensityParam())
#' getmzrt2(xset)
#' }
#' @seealso \code{\link{getdata2}},\code{\link{getupload2}}, \code{\link{getmzrt}}, \code{\link{getdoe}},\code{\link{getmzrtcsv}}
#' @export
getmzrt2 <- function(xset, name = NULL) {
    data <- xcms::featureValues(xset, value = "into")
    # group info
    group <- xset@phenoData@data
    # peaks info
    peaks <- xcms::featureDefinitions(xset)
    mz <- peaks$mzmed
    rt <- peaks$rtmed
    # return as list
    result <- list(data = data, group = group, mz = mz,
        rt = rt)
    if (!is.null(name)) {
            data <- cbind(mz = list$mz, rt = list$rt, list$data)
            data <- t(cbind(group = t(cbind(mz = 'mz',rt = 'rt',t(list$group))), t(data)))
            utils::write.csv(data, file = paste0(name, ".csv"))
    }
    return(result)
}

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
#' @export
getmr <- function(path, index = F, BPPARAM = BiocParallel::SnowParam(),
    pmethod = "hplcorbitrap", minfrac = 0.67, ...) {
    xset <- getdata(path = path, index = index, BPPARAM = BPPARAM,
        pmethod = pmethod, minfrac = minfrac)
    list <- getmzrt(xset)
    return(list)
}

#' Covert the peaks list csv file into list
#' @param path the path to your csv file
#' @return list with rtmz profile and group infomation
#' @seealso \code{\link{getmzrt}}, \code{\link{getmzrt2}}
#' @export
getmzrtcsv <- function(path) {
    dataraw <- utils::read.csv(path, skip = 1)
    mz <- dataraw[, 2]
    rt <- dataraw[, 3]
    data <- dataraw[, -c(1:3)]
    group <- data.frame(t(utils::read.csv(path, nrows = 1)[-(1:3)]))
    colnames(group) <- c(1:ncol(group))
    colnames(data) <- group
    rownames(data) <- dataraw[, 1]
    return(list(data = data, mz = mz, group = group, rt = rt))
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

#' get the data of QC compound for a group of data
#' @param path data path for your QC samples
#' @param mzrange mass of the QC compound
#' @param rtrange retention time of the QC compound
#' @param index index of the files contained QC compounds, default is all of the compounds
#' @return number vector, each number indicate the peak area of that mass and retention time range
#' @export
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

#' plot EIC and boxplot for all peaks and return diffreport
#' @param xset xcmsset object
#' @param name filebase of the sub dir
#' @param test 't' means two-sample welch t-test, 't.equalvar' means two-sample welch t-test with equal variance, 'wilcoxon' means rank sum wilcoxon test, 'f' means F-test, 'pairt' means paired t test, 'blockf' means Two-way analysis of variance, default 't'
#' @param nonpara 'y' means using nonparametric ranked data, 'n' means original data
#' @param ... other parameters for `diffreport`
#' @return diffreport and pdf figure for EIC and boxplot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plote(xset)
#' }
#' @export
plote <- function(xset, name = "test", test = "t", nonpara = "n",
    ...) {
    gt <- xcms::groups(xset)
    a <- xcms::diffreport(xset, filebase = name, eicmax = nrow(gt),
        nonpara = nonpara, ...)
    return(a)
}
