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

        fromPaths <- xcms:::phenoDataFromPaths(files)
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

#' Get the features from t test, with p value, q value, rsd and power restriction
#' @param xod xcmsSet objects
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param rsdt rsd threshold to filter the data
#' @return dataframe with peaks fit the setting above
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getfeaturest(xset)
#' }
#' @export

getfeaturest <- function(xod,
                         power = 0.8,
                         pt = 0.05,
                         qt = 0.05,
                         n = 3,
                         rsdt = 30) {
        data <- xcms::groupval(xod, value = "into")
        idx <- stats::complete.cases(data)
        data1 <- as.data.frame(xcms::groups(xod))
        lv <- xod@phenoData[, 1]
        data0 <- data[idx,]
        sd <- genefilter::rowSds(data0[, 1:n])
        mean <- stats::aggregate(t(data0), list(lv), mean)
        mean <- t(mean[,-1])
        rsd <- sd / mean[, 1] * 100
        sd <- sd[rsd < rsdt]
        mz <- data1$mzmed[idx & rsd < rsdt]
        rt <- data1$rtmed[idx & rsd < rsdt]
        data0 <- data0[rsd < rsdt,]
        ar <- genefilter::rowttests(data0, fac = lv)
        dm <- ar$dm
        p <- ar$p.value
        q <- stats::p.adjust(p, method = "BH")
        m <- nrow(data0)
        df <- cbind.data.frame(sd, dm, p, q, mz, rt, data0)
        df <- df[order(df$p),]
        df$alpha <- c(1:m) * pt / m
        rp <- vector()
        for (i in c(1:nrow(df))) {
                r <- stats::power.t.test(
                        delta = df$dm[i],
                        sd = df$sd[i],
                        sig.level = df$alpha[i],
                        n = n
                )
                rp[i] <- r$power
        }
        df <- cbind(power = rp, df)
        df <- df[df$power > power,]
        return(df)
}

#' Get the features from anova, with p value, q value, rsd and power restriction
#' @param xod xcmsSet objects
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param ng group numbers
#' @param rsdt rsd threshold to filter the data
#' @return dataframe with peaks fit the setting above
#' @export

getfeaturesanova <- function(xod,
                             power = 0.8,
                             pt = 0.05,
                             qt = 0.05,
                             n = 3,
                             ng = 3,
                             rsdt = 30) {
        data <- xcms::groupval(xod, value = "into")
        idx <- stats::complete.cases(data)
        data1 <- as.data.frame(xcms::groups(xod))
        lv <- xod@phenoData[, 1]
        data0 <- data[idx,]

        sd <- genefilter::rowSds(data0[, 1:n])
        mean <- stats::aggregate(t(data0), list(lv), mean)
        mean <- t(mean[,-1])
        sd2 <- genefilter::rowSds(mean)
        rsd <- sd / mean[, 1] * 100
        sd2 <- sd2[rsd < rsdt]
        sd <- sd[rsd < rsdt]
        mz <- data1$mzmed[idx & rsd < rsdt]
        rt <- data1$rtmed[idx & rsd < rsdt]
        data0 <- data0[rsd < rsdt,]
        ar <- genefilter::rowFtests(data0, lv)
        p <- ar$p.value
        q <- stats::p.adjust(p, method = "BH")
        m <- nrow(data0)
        df <- cbind.data.frame(sd, sd2, p, q, mz, rt, data0)
        df <- df[order(df$p),]
        df$alpha <- c(1:m) * pt / m
        rp <- vector()
        for (i in c(1:nrow(df))) {
                r <- stats::power.anova.test(
                        groups = ng,
                        between.var = df$sd2[i],
                        within.var = df$sd[i],
                        sig.level = df$alpha[i],
                        n = n
                )
                rp[i] <- r$power
        }
        df <- cbind(power = rp, df)
        df <- df[df$power > power,]
        return(df)
}

#' output the similarity of two dataset
#' @param xset1 the first dataset
#' @param xset2 the second dateset
#' @return similarity on retention time and rsd %
#' @export
getsim <- function(xset1, xset2) {
        data1 <- gettechrep(xset1)[, c("mzmed", "rtmed",
                                       "rsd")]
        data2 <- gettechrep(xset2)[, c("mzmed", "rtmed",
                                       "rsd")]

        # data1$weight <- ifelse(data1$rsd > 100, 0, 0.2)
        # data1$weight[data1$rsd < 80] <- 0.4
        # data1$weight[data1$rsd < 60] <- 0.6
        # data1$weight[data1$rsd < 40] <- 0.8
        # data1$weight[data1$rsd < 20] <- 1 data1$rtorder
        # <- order(data1$rtmed)
        data1$mzmedn <- round(data1$mzmed, 0.1)

        # data2$weight <- ifelse(data2$rsd > 100, 0, 0.2)
        # data2$weight[data2$rsd < 80] <- 0.4
        # data2$weight[data2$rsd < 60] <- 0.6
        # data2$weight[data2$rsd < 40] <- 0.8
        # data2$weight[data2$rsd < 20] <- 1 data2$rtorder
        # <- order(data2$rtmed)
        data2$mzmedn <- round(data2$mzmed, 0.1)

        data <- merge(data1, data2, by = "mzmedn")
        data <- data[stats::complete.cases(data),]
        cor1 <- cor(data$rtmed.x, data$rtmed.y)
        cor2 <- cor(data$rsd.x, data$rsd.y)
        cor <- c(cor1, cor2)
        return(cor)
}

#' get the data of QC compound for a group of data
#' @param path data path for your QC samples
#' @param mzrange mass of the QC compound
#' @param rtrange retention time of the QC compound
#' @param index index of the files contained QC compounds, default is all of the compounds
#' @return number vector, each number indicate the peak area of that mass and retention time range
#' @export
getQCraw <- function(path, mzrange, rtrange, index = NULL) {
        cdffiles <- list.files(path, recursive = TRUE,
                               full.names = TRUE)
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

#' Isotope extraction for single group of samples with certain mass diff
#' @param xcmsSet  a xcmsSet object
#' @param massdiff mass defection
#' @param rtwindow retention time range
#' @param mzwindow mass charge ratio window
#' @param ppm resolution of the mass spectrum
#' @return table with mass, retention time, scaled mass and scaled mass defect
getmassdiff <- function(xcmsSet, massdiff, rtwindow,
                        mzwindow, ppm) {
        # get group infomation
        groups = data.frame(xcmsSet@groups)
        peakIntensities = xcms::groupval(xcmsSet, "medret",
                                         "inio")
        # order peaks by rt
        peakIntensities = peakIntensities[order(groups$rtmed),]
        groups <- groups[order(groups$rtmed),]
        groups$peakins <- apply(peakIntensities, 1, mean)
        result <- NULL
        # search:
        for (i in 1:nrow(groups)) {
                bin = groups[groups$rtmed - groups$rtmed[i] >=
                                     0 &
                                     groups$rtmed - groups$rtmed[i] <= rtwindow,]
                if (nrow(bin) > 1) {
                        dis <- stats::dist(bin$mzmed, method = "manhattan") / massdiff
                        df <-
                                data.frame(
                                        ms1 = bin$mzmed[which(lower.tri(dis),
                                                              arr.ind = T)[, 1]],
                                        ms2 = bin$mzmed[which(lower.tri(dis),
                                                              arr.ind = T)[, 2]],
                                        diff = as.numeric(dis)
                                )
                        df$rdiff <- round(df$diff)
                        dfn <-
                                df[df$diff <= df$rdiff * (1 + ppm / 1e+06) +
                                           (df$ms1 * ppm / 1e+06) / (massdiff * (1 -
                                                                                         ppm /
                                                                                         1e+06)) && df$diff >= df$rdiff *
                                           (1 - ppm / 1e+06) - (df$ms1 * ppm /
                                                                        1e+06) / (massdiff *
                                                                                          (1 + ppm /
                                                                                                   1e+06)),]
                        dfn$msdiff <- abs(dfn$ms1 - dfn$ms2)
                        dfn <- dfn[dfn$msdiff < mzwindow,]
                        # candidate number of labeled atoms
                        result <- rbind(result, bin[bin$mzmed %in%
                                                            dfn$ms1 |
                                                            bin$mzmed %in% dfn$ms2,])
                        result <- result[rownames(unique(result[,
                                                                c("mzmed", "rtmed")])),]
                }
        }
        result$sm <- result$mzmed * massdiff
        result$smd <- ceiling(result$sm) - result$sm
        return(result)
}
