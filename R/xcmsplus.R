#' Get xcmsset object in one step with optimized methods.
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the references
#' @param ... arguments for xcmsSet function
#' @details the parameters are extracted from the papers. If you use name other than the name above, you will use the default setting of XCMS. Also I suggest IPO packages or apLCMS packages to get reasonable data for your own instrumental. If you want to summit the results to a paper, remember to include those parameters.
#' @return a xcmsset object for that path or selected samples
#' @references Patti, G. J.; Tautenhahn, R.; Siuzdak, G. Nat. Protocols 2012, 7 (3), 508–516.
#' @export
getdata <- function(path, index = F, BPPARAM = BiocParallel::SnowParam(workers = 12),
                    pmethod = "hplcorbitrap", ...) {
        cdffiles <- list.files(path, recursive = TRUE,
                               full.names = TRUE)
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
                        xset <- xcms::group(xset, bw = 5, mzwid = 0.015)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <- xcms::group(xset2, bw = 5, mzwid = 0.015)
                        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                }
        } else if (pmethod == "uplcorbitrap") {
                xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                      method = "centWave", ppm = 2.5, peakwidth = c(5,
                                                                                    20), prefilter = c(3, 5000), ...)
                xset <- xcms::group(xset, bw = 2, mzwid = 0.015)
                xset2 <- xcms::retcor(xset)
                # you need group the peaks again for this corrected
                # data
                xset2 <- xcms::group(xset2, bw = 2, mzwid = 0.015)
                xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
        } else if (pmethod == "hplcqtof") {
                xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                      method = "centWave", ppm = 30, peakwidth = c(10,
                                                                                   60), prefilter = c(0, 0), ...)
                if (index & length(index) == 1) {
                        xset3 <- xset
                } else {
                        xset <- xcms::group(xset, bw = 5, mzwid = 0.025)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <- xcms::group(xset2, bw = 5, mzwid = 0.025)
                        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                }
        } else if (pmethod == "hplchqtof") {
                xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                      method = "centWave", ppm = 15, peakwidth = c(10,
                                                                                   60), prefilter = c(0, 0), ...)
                if (index & length(index) == 1) {
                        xset3 <- xset
                } else {
                        xset <- xcms::group(xset, bw = 5, mzwid = 0.015)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <- xcms::group(xset2, bw = 5, mzwid = 0.015)
                        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                }
        } else if (pmethod == "uplcqtof") {
                xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                      method = "centWave", ppm = 30, peakwidth = c(5,
                                                                                   20), prefilter = c(0, 0), ...)
                if (index & length(index) == 1) {
                        xset3 <- xset
                } else {
                        xset <- xcms::group(xset, bw = 2, mzwid = 0.025)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <- xcms::group(xset2, bw = 2, mzwid = 0.025)
                        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                }
        } else if (pmethod == "uplchqtof") {
                xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                      method = "centWave", ppm = 15, peakwidth = c(5,
                                                                                   20), prefilter = c(0, 0), ...)
                if (index & length(index) == 1) {
                        xset3 <- xset
                } else {
                        xset <- xcms::group(xset, bw = 2, mzwid = 0.015)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <- xcms::group(xset2, bw = 2, mzwid = 0.015)
                        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                }
        } else {
                xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                      ...)
                if (index & length(index) == 1) {
                        xset3 <- xset
                } else {
                        xset <- xcms::group(xset)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <- xcms::group(xset2)
                        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                }
        }
        return(xset3)
}

#' Get the csv files to be submitted to Metaboanalyst
#' @param xset the xcmsset object which you want to submitted to Metaboanalyst
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param name file name
#' @return dataframe with data needed for Metaboanalyst if your want to perform local analysis.
#' @export
getupload <- function(xset, method = "medret", intensity = "into",
                      name = "Peaklist") {
        peakIntensities <- xcms::groupval(xset, method,
                                          intensity)
        if (intensity == "intb") {
                peakIntensities[is.na(peakIntensities)] = 0
        }
        data <- rbind(group = as.character(xcms::phenoData(xset)$class),
                      peakIntensities)
        data <- data[!duplicated(rownames(data)), ]
        filename <- paste0(name, ".csv")
        utils::write.csv(data, file = filename)
        return(data)
}
#' Get the report for technique replicates.
#' @param xset the xcmsset object which for all of your technique replicates for one sample
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @return dataframe with mean, standard deviation and RSD for those technique replicates combined with raw data
#' @export
gettechrep <- function(xset, method = "medret", intensity = "into") {
        data <- t(xcms::groupval(xset, method, intensity))
        lv <- xset@phenoData[, 1]
        mean <- stats::aggregate(data, list(lv), mean)
        sd <- stats::aggregate(data, list(lv), sd)
        suppressWarnings(rsd <- sd/mean * 100)
        result <- data.frame(cbind(t(mean[, -1]), t(sd[,
                                                       -1]), t(rsd[, -1])))
        colnames(result) <- c("mean", "sd", "rsd")
        datap <- xcms::groups(xset)
        report <- cbind.data.frame(datap, result)
        return(report)
}

#' Get the report for samples with technique replicates
#' @param xset the xcmsset object all of samples with technique replicates
#' @param anno logical if set as True, it will return the table for further annotation, default false
#' @param peaklist logical if set as True, it will return csv files for metaboanalyst, default false
#' @param file file name for the peaklist
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @return dataframe with mean, standard deviation and RSD for those technique replicates combined with raw data for all of the samples if anno and peaklist are defaults false.
#' @export
gettechbiorep <- function(xset, anno = F, peaklist = F,
                          file = NULL, method = "medret", intensity = "into") {
        data <- t(xcms::groupval(xset, method, intensity))
        lv <- xset@phenoData[, 1]
        lv2 <- xset@phenoData[, 2]
        mean <- stats::aggregate(data, list(lv, lv2), mean)
        sd <- stats::aggregate(data, list(lv, lv2), sd)
        suppressWarnings(rsd <- sd/mean * 100)
        result <- data.frame(cbind(t(mean[, -c(1:2)]),
                                   t(sd[, -c(1:2)]), t(rsd[, -c(1:2)])))
        name <- unique(c(paste0(lv, lv2)))
        colnames(result) <- c(paste0(name, "mean"), paste0(name,
                                                           "sd"), paste0(name, "rsd%"))
        datap <- xcms::groups(xset)
        report <- cbind.data.frame(datap, result)
        if (anno) {
                anno <- as.data.frame(cbind(xset@groups[, 1],
                                            xset@groups[, 4], cbind(t(mean[, -c(1:2)]))))
                colnames(anno) <- c("mz", "time", name)
                return(anno)
        } else if (peaklist) {
                result <- data.frame(t(mean[, -c(1:2)]))
                data <- rbind(group = name, result)
                utils::write.csv(data, file = file)
        } else {
                return(report)
        }
}


#' plot the scatter plot for xcmsset (or two) objects with threshold
#' @param data1 the first xcmsset data set
#' @param data2 the second xcmsset data set
#' @param threshold the threshold of the response (log based 10)
#' @param ms the mass range to plot the data
#' @param ... parameters for `plot` function
#' @return NULL
#' @export
plotmr <- function(data1,
                   data2 = NULL,
                   threshold = 5,
                   ms = c(100, 1000),
                   ...) {
        data1mzrt <- xcms::groups(data1)
        mz <- data1mzrt[,1]
        rt <- data1mzrt[,4]
        data1into <- xcms::groupval(data1,'medret','into')
        into <- apply(data1into,1,mean)
        data1 <- cbind.data.frame(mz,rt,into)

        graphics::plot(
                data1$mz[log10(data1$into+1) > threshold] ~ data1$rt[log10(data1$into+1) > threshold],
                xlab = "Retention Time",
                ylab = "m/z",
                ylim = ms,
                cex = log10(data1$into+1) - threshold + 1,
                col = grDevices::rgb(0,0, 1, 0.1),
                pch = 19,
                ...
        )
        if(!is.null(data2)){
                data2mzrt <- xcms::groups(data2)
                mz <- data2mzrt[,1]
                rt <- data2mzrt[,4]
                data2into <- xcms::groupval(data2,'medret','into')
                into <- apply(data2into,1,mean)
                data2 <- cbind.data.frame(mz,rt,into)
                graphics::points(data2$mz[log10(data2$into+1) > threshold] ~ data2$rt[log10(data2$into+1) > threshold],
                                 cex = log10(data2$into+1) - threshold + 1,
                                 col = grDevices::rgb(1,0, 0, 0.1),
                                 pch = 19)
        }
}

#' plot EIC and boxplot for all peaks and return diffreport
#' @param xset xcmsset object
#' @param name filebase of the sub dir
#' @param test 't' means two-sample welch t-test, 't.equalvar' means two-sample welch t-test with equal variance, 'wilcoxon' means rank sum wilcoxon test, 'f' means F-test, 'pairt' means paired t test, 'blockf' means Two-way analysis of variance, default 't'
#' @param nonpara 'y' means using nonparametric ranked data, 'n' means original data
#' @param ... other parameters for `diffreport`
#' @return diffreport and pdf figure for EIC and boxplot
#' @export
plote <- function(xset,
                  name = "test",
                  test = "t",
                  nonpara = "n",
                  ...) {
        gt <- xcms::groups(xset)
        a <-
                xcms::diffreport(
                        xset,
                        filebase = name,
                        eicmax = nrow(gt),
                        nonpara = nonpara,
                        ...
                )
        return(a)
}

#' Get the features from t test, with p value, q value, rsd and power restriction
#' @param xod xcmsSet objects
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param rsdt rsd threshold to filter the data
#' @return dataframe with peaks fit the setting above
#' @export

getfeaturest <- function(xod, power = 0.8, pt = 0.05,
    qt = 0.05, n = 3, rsdt = 50) {
    data <- xcms::groupval(xod, value = "into")
    idx <- stats::complete.cases(data)
    data1 <- as.data.frame(xcms::groups(xod))
    lv <- xod@phenoData[, 1]
    data0 <- data[idx, ]
    sd <- genefilter::rowSds(data0[, 1:n])
    mean <- stats::aggregate(t(data0), list(lv), mean)
    mean <- t(mean[, -1])
    rsd <- sd/mean[, 1] * 100
    sd <- sd[rsd < rsdt]
    mz <- data1$mzmed[idx & rsd < rsdt]
    rt <- data1$rtmed[idx & rsd < rsdt]
    data0 <- data0[rsd < rsdt, ]
    ar <- genefilter::rowttests(data0, fac = lv)
    dm <- ar$dm
    p <- ar$p.value
    q <- stats::p.adjust(p, method = "BH")
    m <- nrow(data0)
    df <- cbind.data.frame(sd, dm, p, q, mz, rt, data0)
    df <- df[order(df$p), ]
    df$alpha <- c(1:m) * pt/m
    rp <- vector()
    for (i in c(1:nrow(df))) {

        r <- stats::power.t.test(delta = df$dm[i],
            sd = df$sd[i], sig.level = df$alpha[i],
            n = n)
        rp[i] <- r$power
    }
    df <- cbind(power = rp, df)
    df <- df[df$power > power, ]
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

getfeaturesanova <- function(xod, power = 0.8, pt = 0.05,
    qt = 0.05, n = 3, ng = 3, rsdt = 50) {
    data <- xcms::groupval(xod, value = "into")
    idx <- stats::complete.cases(data)
    data1 <- as.data.frame(xcms::groups(xod))
    lv <- xod@phenoData[, 1]
    data0 <- data[idx, ]

    sd <- genefilter::rowSds(data0[, 1:n])
    mean <- stats::aggregate(t(data0), list(lv), mean)
    mean <- t(mean[, -1])
    sd2 <- genefilter::rowSds(mean)
    rsd <- sd/mean[, 1] * 100
    sd2 <- sd2[rsd < rsdt]
    sd <- sd[rsd < rsdt]
    mz <- data1$mzmed[idx & rsd < rsdt]
    rt <- data1$rtmed[idx & rsd < rsdt]
    data0 <- data0[rsd < rsdt, ]
    ar <- genefilter::rowFtests(data0, lv)
    p <- ar$p.value
    q <- stats::p.adjust(p, method = "BH")
    m <- nrow(data0)
    df <- cbind.data.frame(sd, sd2, p, q, mz, rt, data0)
    df <- df[order(df$p), ]
    df$alpha <- c(1:m) * pt/m
    rp <- vector()
    for (i in c(1:nrow(df))) {
        r <- stats::power.anova.test(groups = ng, between.var = df$sd2[i],
            within.var = df$sd[i], sig.level = df$alpha[i],
            n = n)
        rp[i] <- r$power
    }
    df <- cbind(power = rp, df)
    df <- df[df$power > power, ]
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
        data <- data[stats::complete.cases(data), ]
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
        peakIntensities = peakIntensities[order(groups$rtmed),
                                          ]
        groups <- groups[order(groups$rtmed), ]
        groups$peakins <- apply(peakIntensities, 1, mean)
        result <- NULL
        # search:
        for (i in 1:nrow(groups)) {
                bin = groups[groups$rtmed - groups$rtmed[i] >=
                                     0 & groups$rtmed - groups$rtmed[i] <= rtwindow,
                             ]
                if (nrow(bin) > 1) {
                        dis <- stats::dist(bin$mzmed, method = "manhattan")/massdiff
                        df <- data.frame(ms1 = bin$mzmed[which(lower.tri(dis),
                                                               arr.ind = T)[, 1]], ms2 = bin$mzmed[which(lower.tri(dis),
                                                                                                         arr.ind = T)[, 2]], diff = as.numeric(dis))
                        df$rdiff <- round(df$diff)
                        dfn <- df[df$diff <= df$rdiff * (1 + ppm/1e+06) +
                                          (df$ms1 * ppm/1e+06)/(massdiff * (1 -
                                                                                    ppm/1e+06)) && df$diff >= df$rdiff *
                                          (1 - ppm/1e+06) - (df$ms1 * ppm/1e+06)/(massdiff *
                                                                                          (1 + ppm/1e+06)), ]
                        dfn$msdiff <- abs(dfn$ms1 - dfn$ms2)
                        dfn <- dfn[dfn$msdiff < mzwindow, ]
                        # candidate number of labeled atoms
                        result <- rbind(result, bin[bin$mzmed %in%
                                                            dfn$ms1 | bin$mzmed %in% dfn$ms2, ])
                        result <- result[rownames(unique(result[,
                                                                c("mzmed", "rtmed")])), ]
                }
        }
        result$sm <- result$mzmed * massdiff
        result$smd <- ceiling(result$sm) - result$sm
        return(result)
}
