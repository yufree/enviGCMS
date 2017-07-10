#' Get xcmsset object in one step with optimized methods.
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the references
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
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getupload(xset)
#' }
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

#' Get the peak list with blank samples' peaks removed
#' @param xset the xcmsset object with blank and certain group samples' data
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 1000
#' @return diff report
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getbgremove(xset)
#' }
#' @export
getbgremove <- function(xset, method = "medret", intensity = "into", file = NULL, rsdcf = 30, inscf = 1000){
        data0 <- xcms::groupval(xset, method, intensity)
        if (intensity == "intb") {
                data0[is.na(data0)] = 0
        }
        data <- t(data0)
        lv <- xset@phenoData[, 1]
        mean <- stats::aggregate(data, list(lv), mean)
        sd <- stats::aggregate(data, list(lv), sd)
        suppressWarnings(rsd <- sd/mean * 100)
        result <- data.frame(cbind(t(mean[,-1]), t(rsd[, -1])))
        index <- t(rsd[, -1]) < rsdcf & t(mean[, -1]) > inscf

        diff <- result[,2] - result[,1]
        datap <- xcms::groups(xset)
        report <- cbind.data.frame(xset@groups[, 1], xset@groups[, 4], diff)
        colnames(report) <- c("mz", "time", "diff")

        N <- sum(index)
        L <- length(index)

        report <- report[index,]

        message(paste(N,'out of',L,'peaks found with rsd cutoff',rsdcf,'and intensity cutoff', inscf))
        if (!is.null(file)) {
                utils::write.csv(report, file = paste0(file,'.csv'), row.names = F)}
        return(report)
}


#' Get the report for technique replicates.
#' @param xset the xcmsset object which for all of your technique replicates for one sample
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 1000
#' @return dataframe with mean, standard deviation and RSD for those technique replicates combined with raw data
#' @export
gettechrep <- function(xset, method = "medret", intensity = "into", file = NULL, rsdcf = 30, inscf = 1000) {
        data0 <- xcms::groupval(xset, method, intensity)
        if (intensity == "intb") {
                data0[is.na(data0)] = 0
        }
        data <- t(data0)
        lv <- xset@phenoData[, 1]
        mean <- stats::aggregate(data, list(lv), mean)
        sd <- stats::aggregate(data, list(lv), sd)
        suppressWarnings(rsd <- sd/mean * 100)
        result <- data.frame(cbind(t(mean[, -1]), t(sd[,
                                                       -1]), t(rsd[, -1])))
        index <- t(rsd[, -1]) < rsdcf & t(mean[, -1]) > inscf
        colnames(result) <- c("mean", "sd", "rsd")
        datap <- xcms::groups(xset)
        report <- cbind.data.frame(datap, result)
        N <- sum(index)
        L <- length(index)

        message(paste(N,'out of',L,'peaks found with rsd cutoff',rsdcf,'and intensity cutoff', inscf))

        report <- report[index,]

        if (!is.null(file)) {
                anno <- cbind.data.frame(xset@groups[, 1], xset@groups[, 4], t(mean[, -1]))
                colnames(anno) <- c("mz", "time", as.character(t(mean[, 1])))
                utils::write.csv(anno, file = paste0(file,'.csv'), row.names = F)
                anno <- anno[index,]
                return(anno)
        } else {
                return(report)
        }
}

#' Get the report for biological replicates.
#' @param xset the xcmsset object which for all of your technique replicates for bio replicated sample in single group
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 0
#' @return dataframe with mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data
#' @export
getbiotechrep <- function(xset, method = "medret", intensity = "into", file = NULL, rsdcf = 30, inscf = 1000) {
        data0 <- xcms::groupval(xset, method, intensity)
        if (intensity == "intb") {
                data0[is.na(data0)] = 0
        }
        data <- t(data0)
        lv <- xset@phenoData[, 1]
        mean <- stats::aggregate(data, list(lv), mean)
        sd <- stats::aggregate(data, list(lv), sd)
        suppressWarnings(rsd <- sd/mean * 100)
        result <- data.frame(cbind(t(mean[, -1]), t(sd[,-1]), t(rsd[, -1])))
        colnames(result) <- c(paste0(t(mean[, 1]), "mean"), paste0(t(sd[,1]),"sd"), paste0(t(mean[, 1]), "rsd%"))

        meanB <- apply(t(mean[, -1]),1,mean)
        sdB <- apply(t(mean[,-1]),1,sd)
        rsdB <- sdB/meanB * 100

        index <- rsdB < rsdcf & meanB > inscf

        datap <- xcms::groups(xset)
        report <- cbind.data.frame(datap, result,meanB,sdB,rsdB)
        report <- report[index,]

        N <- sum(index)
        L <- length(index)

        message(paste(N,'out of',L,'peaks found with rsd cutoff',rsdcf,'and intensity cutoff', inscf))

        report <- report[index,]

        if (!is.null(file)) {
                anno <- cbind.data.frame(xset@groups[, 1], xset@groups[, 4], t(mean[, -1]))
                colnames(anno) <- c("mz", "time", as.character(t(mean[, 1])))
                utils::write.csv(anno, file = paste0(file,'.csv'), row.names = F)
                anno <- anno[index,]
                return(anno)
        } else {
                return(report)
        }
}

#' Get the report for samples with biological and technique replicates in different groups
#' @param xset the xcmsset object all of samples with technique replicates
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for the peaklist to MetaboAnalyst
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 1000
#' @return dataframe with mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data in different groups if file are defaults NULL.
#' @export
getgrouprep <- function(xset, file = NULL, method = "medret", intensity = "into", rsdcf = 30, inscf = 1000) {
        data0 <- xcms::groupval(xset, method, intensity)
        if (intensity == "intb") {
                data0[is.na(data0)] = 0
        }
        data <- t(data0)
        lv <- xset@phenoData[, 1]
        lv2 <- xset@phenoData[, 2]
        mean <- stats::aggregate(data, list(lv, lv2), mean)
        sd <- stats::aggregate(data, list(lv, lv2), sd)
        suppressWarnings(rsd <- sd/mean * 100)
        result <- cbind.data.frame(t(mean[, -c(1:2)]), t(sd[, -c(1:2)]), t(rsd[, -c(1:2)]))
        name <- paste0(t(mean[, 1]),t(mean[, 2]))
        colnames(result) <- c(paste0(name, "mean"), paste0(name, "sd"), paste0(name, "rsd%"))

        meanB <- stats::aggregate(mean[, -c(1:2)], list(as.vector(t(mean[, 1]))), mean)
        sdB <- stats::aggregate(mean[, -c(1:2)], list(as.vector(t(mean[, 1]))), sd)
        suppressWarnings(rsdB <- sdB[,-1]/meanB[,-1] * 100)

        resultB <- cbind.data.frame(t(meanB[,-1]), t(sdB[,-1]), t(rsdB))
        nameB <- as.vector(t(meanB[, 1]))
        colnames(resultB) <- c(paste0(nameB, "mean"), paste0(nameB,"sd"), paste0(nameB, "rsd%"))

        datap <- xcms::groups(xset)
        report <- cbind.data.frame(datap, result, resultB)

        index <- as.vector(apply(rsdB,2,function(x) all(x<rsdcf))) & as.vector(apply(meanB[,-1],2,function(x) all(x>inscf)))

        N <- sum(index)
        L <- length(index)

        message(paste(N,'out of',L,'peaks found with rsd cutoff',rsdcf,'and intensity cutoff', inscf))

        report <- report[index,]

        if (!is.null(file)) {
                result <- data.frame(t(mean[, -c(1:2)]))[index,]
                data <- rbind(group = as.character(mean[, 1]), result)
                utils::write.csv(data, file = paste0(file,'.csv'))
                return(data)
        } else {
                return(report)
        }
}

#' Get the time series or two factor DoE report for samples with biological and technique replicates in different groups
#' @param xset the xcmsset object all of samples with technique replicates in time series or two factor DoE
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for the peaklist to MetaboAnalyst
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 1000
#' @return dataframe with time series or two factor DoE mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data in different groups if file are defaults NULL.
#' @export
gettimegrouprep <- function(xset, file = NULL, method = "medret", intensity = "into", rsdcf = 30, inscf = 1000) {
        data0 <- xcms::groupval(xset, method, intensity)
        if (intensity == "intb") {
                data0[is.na(data0)] = 0
        }
        data <- t(data0)
        lv <- xset@phenoData[, 1]
        lv2 <- xset@phenoData[, 2]
        lv3 <- xset@phenoData[, 3]
        mean <- stats::aggregate(data, list(lv, lv2, lv3), mean)
        sd <- stats::aggregate(data, list(lv, lv2, lv3), sd)
        suppressWarnings(rsd <- sd/mean * 100)
        result <- cbind.data.frame(t(mean[, -c(1:3)]), t(sd[, -c(1:3)]), t(rsd[, -c(1:3)]))
        name <- paste0(t(mean[, 1]),t(mean[, 2]),t(mean[, 3]))
        colnames(result) <- c(paste0(name, "mean"), paste0(name, "sd"), paste0(name, "rsd%"))

        meanB <- stats::aggregate(mean[, -c(1:3)], list(as.vector(t(mean[, 1])),as.vector(t(mean[, 2]))), mean)
        sdB <- stats::aggregate(mean[, -c(1:3)], list(as.vector(t(mean[, 1])),as.vector(t(mean[, 2]))), sd)
        suppressWarnings(rsdB <- sdB[,-c(1:2)]/meanB[,-c(1:2)] * 100)

        resultB <- cbind.data.frame(t(meanB[,-c(1:2)]), t(sdB[,-c(1:2)]), t(rsdB))
        nameB <- paste0(t(meanB[, 1]), t(meanB[, 2]))
        colnames(resultB) <- c(paste0(nameB, "mean"), paste0(nameB,"sd"), paste0(nameB, "rsd%"))

        datap <- xcms::groups(xset)
        report <- cbind.data.frame(datap, result, resultB)

        index <- as.vector(apply(rsdB,2,function(x) all(x<rsdcf))) & as.vector(apply(meanB[,-c(1,2)],2,function(x) all(x>inscf)))

        N <- sum(index)
        L <- length(index)

        message(paste(N,'out of',L,'peaks found with rsd cutoff',rsdcf,'and intensity cutoff', inscf))

        report <- report[index,]

        if (!is.null(file)) {
                result <- data.frame(t(meanB[, -c(1:2)]))[index,]
                data <- rbind(time = as.character(meanB[, 1]), group = as.character(meanB[, 2]), result)
                utils::write.csv(data, file = paste0(file,'.csv'))
                return(data)
        } else {
                return(report)
        }
}

#' plot the scatter plot for xcmsset (or two) objects with threshold
#' @param data1 the first xcmsset object
#' @param data2 the second xcmsset object
#' @param threshold the threshold of the response (log based 10)
#' @param ms the mass range to plot the data
#' @param ... parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plotmr(xset)
#' }
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
                col = grDevices::rgb(0,0, 1, 0.2),
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
                                 col = grDevices::rgb(1,0, 0, 0.2),
                                 pch = 19)
        }
}

#' plot the diff scatter plot for one xcmsset objects with threshold and two groups
#' @param data xcmsset object with two groups
#' @param threshold the threshold of the response (log based 10)
#' @param ms the mass range to plot the data
#' @param ... parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plotmrc(xset)
#' }
#' @export
plotmrc <- function(data,
                   threshold = 5,
                   ms = c(100, 1000),
                   ...) {
        mzrt <- xcms::groups(data)
        mz <- mzrt[,1]
        rt <- mzrt[,4]
        into <- xcms::groupval(data,'medret','into')
        lv <- data@phenoData[, 1]
        mean <- t(stats::aggregate(t(into), list(lv), mean)[,-1])
        lvname <- stats::aggregate(t(into), list(lv), mean)[,1]
        diff1 <- mean[,1]-mean[,2]
        diff2 <- mean[,2]-mean[,1]
        diff1[diff1 < 0] <- 0
        diff2[diff2 < 0] <- 0
        name1 <- paste0(lvname[1],"-",lvname[2])
        name2 <- paste0(lvname[2],"-",lvname[1])
        data <- cbind.data.frame(mz,rt,diff1,diff2)

        graphics::plot(
                data$mz[log10(data$diff1+1) > threshold] ~ data$rt[log10(data$diff1+1) > threshold],
                xlab = "Retention Time",
                ylab = "m/z",
                ylim = ms,
                cex = log10(data$diff1+1) - threshold + 1,
                col = grDevices::rgb(0,0, 1, 0.618),
                pch = 19,
                ...
        )

        graphics::points(
                data$mz[log10(data$diff2+1) > threshold] ~ data$rt[log10(data$diff2+1) > threshold],
                cex = log10(data$diff2+1) - threshold + 1,
                col = grDevices::rgb(1,0, 0, 0.618),
                pch = 19,
                ...
        )

        graphics::legend('topright', legend = c(name1,name2),pch = 19,col = c(grDevices::rgb(0,0, 1, 0.618),grDevices::rgb(1,0, 0, 0.618)),bty = 'n')

}

#' plot the rsd influnces of data
#' @param xset xcmsset data
#' @param ... other parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plotrsd(xset)
#' }
#' @export
plotrsd <- function(xset,...){
        df <- getbiotechrep(xset)
        mz <- df$mzmed
        rt <- df$rtmin
        cex <- df$rsd
        graphics::plot(mz~rt,cex = scale(cex),
             xlab = 'retention time', ylab = 'm/z',
             ...)
}

#' plot the PCA of xcmsset
#' @param xset xcmsset data
#' @param lv group information
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param ... other parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plotpca(xset)
#' }
#' @export

plotpca <- function(xset,lv = NULL,center = T, scale = T,...){
        data <- xcms::groupval(xset,'medret','into')
        if (is.null(lv)) {
                pch = colnames(data)
        } else {
                pch = lv
        }
        pcao <- stats::prcomp(t(data), center = center,
                              scale = scale)
        pcaoVars = signif(((pcao$sdev)^2)/(sum((pcao$sdev)^2)),
                          3) * 100
        graphics::plot(pcao$x[, 1], pcao$x[, 2], xlab = paste("PC1:",
                                                              pcaoVars[1], "% of Variance Explained"), ylab = paste("PC2:",
                                                                                                                    pcaoVars[2], "% of Variance Explained"), pch = pch,
                       cex = 2, ...)
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
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plote(xset)
#' }
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
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getfeaturest(xset)
#' }
#' @export

getfeaturest <- function(xod, power = 0.8, pt = 0.05,
    qt = 0.05, n = 3, rsdt = 30) {
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
    qt = 0.05, n = 3, ng = 3, rsdt = 30) {
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
