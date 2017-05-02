#' import data and return the annotated matrix for GC-MS
#' @param data file type which xcmsRaw could handle
#' @param profstep the m/z step for generating matrix data from raw mass spectral data.
#' @param time round numbers of retention time, default is 1
#' @return matrix with the row as increasing m/z second and column as increasing scantime
#' @export
getmd <- function(data,profstep = 1,time = 1){
        data <- xcms::xcmsRaw(data,profstep)
        z1 <- data@env$profile
        zf <- as.factor(round(data@scantime,time))
        df <- aggregate(t(z1), list(zf), sum)[-1]
        rownames(df) <- unique(round(data@scantime,time))
        colnames(df) <- seq(data@mzrange[1],data@mzrange[2],by = profstep)
        return(t(as.matrix(df)))
}

#' Combine two or more matrix
#'
#' @param data1 matrix lower mass range
#' @param data2 matrix higher mass range
#' @param ... matrix even higher mass range
#' @return matrix with the row as scantime in second and column as m/z
#' @export
combinemd <- function(data1,data2,...){
        if(missing(...)){
                z1 <- getmd(data1)
                z2 <- getmd(data2)
                ind <- intersect(colnames(z1),colnames(z2))
                z <- rbind(z1[as.character(ind)],z2[,as.character(ind)])
                rownames(z) <- c(seq(min(rownames(z1)),max(rownames(z1))),seq(min(colnames(z2)),max(colnames(z2))))
                colnames(z) <- ind
                return(z)
        }
        else{
                combinemd(data1, combinemd(data2, ...))
        }
}

#' Subset the data mass spectrum of certain retention time and plot them
#' @param data imported data matrix of GC-MS
#' @param rt vector range of the retention time min
#' @param ms vector range of the m/z
#' @return data matrix
#' @export
getsubmd <- function(data,rt,ms){
        mzindexstart <- as.numeric(head(rownames(data),1))
        rtindexstart <- as.numeric(head(colnames(data),1))
        rts <- rt*60
        rt1 <- which(round(as.numeric(colnames(data))) == round(rts[1]))[1]
        rt2 <- which(round(as.numeric(colnames(data))) == round(rts[2]))[1]
        mzs <- ms-mzindexstart+1
        mz1 <- min(mzs)
        mz2 <- max(mzs)
        data <- data[mz1:mz2,]
        data <- t(data)[rt1:rt2,]
        return(t(data))
}

#' Write MSP files for NIST search
#' @param mz a intensity vector, who name is the mass in m/z
#' @param outfilename the name of the MSP file, default is 'unknown'
#' @return none a MSP file will be created at the subfolder working dictionary with name 'MSP'
#' @export
writeMSP<-function(mz, outfilename="unknown"){
        mz <- paste(names(mz),round(mz))
        dir.create('MSP')
        zz <- file(file.path('MSP',paste(outfilename,".msp",sep="")), "w")
        nPeaks <- length(mz)
        cat("Name: unknown", paste("Num Peaks: ",nPeaks),  file = zz, sep = "\n")
        while (length(mz) >=5 ){
                cat(paste(mz[1:5]),"", file = zz, sep="; ")
                cat(paste("\n"), file = zz)
                mz<-mz[6:length(mz)]
        }
        if(!is.na(mz[1])){
                cat(paste(mz),"", file = zz, sep="; ")
                cat(paste("\n"), file = zz)
        }
        close(zz)
        print(paste("A data file",outfilename,".MSP has been generated in the folder:", 'MSP', cat("\n")))
}

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
getdata <-
        function(path,
                 index = F,
                 BPPARAM = BiocParallel::SnowParam(workers = 12),
                 pmethod = 'hplcorbitrap',
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
                if (index) {
                        cdffiles <- cdffiles[index]
                }
                if (pmethod == 'hplcorbitrap') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 2.5,
                                        peakwidth = c(10, 60),
                                        prefilter = c(3, 5000),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplcorbitrap') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 2.5,
                                        peakwidth = c(5, 20),
                                        prefilter = c(3, 5000),
                                        ...
                                )
                        xset <-
                                xcms::group(xset, bw = 2, mzwid = 0.015)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected data
                        xset2 <-
                                xcms::group(xset2, bw = 2, mzwid = 0.015)
                        xset3 <-
                                xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                } else if (pmethod == 'hplcqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 30,
                                        peakwidth = c(10, 60),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.025)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.025)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'hplchqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 15,
                                        peakwidth = c(10, 60),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplcqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 30,
                                        peakwidth = c(5, 20),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 2,
                                                    mzwid = 0.025)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 2,
                                                    mzwid = 0.025)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplchqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 15,
                                        peakwidth = c(5, 20),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 2,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 2,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else{
                        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM, ...)
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <- xcms::group(xset2)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
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
getupload <-
        function(xset,
                 method = "medret",
                 intensity = 'into',
                 name = 'Peaklist') {
                peakIntensities <- xcms::groupval(xset, method, intensity)
                if (intensity == "intb") {
                        peakIntensities[is.na(peakIntensities)] = 0
                }
                data <-
                        rbind(group = as.character(xcms::phenoData(xset)$class), peakIntensities)
                data <- data[!duplicated(rownames(data)),]
                filename <- paste0(name, '.csv')
                utils::write.csv(data, file = filename)
                return(data)
        }
#' Get the report for technique replicates.
#' @param xset the xcmsset object which for all of your technique replicates for one sample
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @return dataframe with mean, standard deviation and RSD for those technique replicates combined with raw data
#' @export
gettechrep <- function(xset,
                       method =  'medret',
                       intensity = 'into') {
        data <- t(xcms::groupval(xset, method, intensity))
        lv <- xset@phenoData[, 1]
        mean <- stats::aggregate(data, list(lv), mean)
        sd <- stats::aggregate(data, list(lv), sd)
        suppressWarnings(rsd <- sd / mean * 100)
        result <-
                data.frame(cbind(t(mean[, -1]), t(sd[, -1]), t(rsd[, -1])))
        colnames(result) <-
                c('mean','sd','rsd')
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
gettechbiorep <-
        function(xset,
                 anno = F,
                 peaklist = F,
                 file = NULL,
                 method =  'medret',
                 intensity = 'into') {
                data <- t(xcms::groupval(xset, method, intensity))
                lv <- xset@phenoData[, 1]
                lv2 <- xset@phenoData[, 2]
                mean <- stats::aggregate(data, list(lv, lv2), mean)
                sd <- stats::aggregate(data, list(lv, lv2), sd)
                suppressWarnings(rsd <- sd / mean * 100)
                result <-
                        data.frame(cbind(t(mean[, -c(1:2)]), t(sd[, -c(1:2)]), t(rsd[, -c(1:2)])))
                name <- unique(c(paste0(lv, lv2)))
                colnames(result) <-
                        c(paste0(name, 'mean'),
                          paste0(name, 'sd'),
                          paste0(name, 'rsd%'))
                datap <- xcms::groups(xset)
                report <- cbind.data.frame(datap, result)
                if (anno) {
                        anno <-
                                as.data.frame(cbind(xset@groups[, 1], xset@groups[, 4], cbind(t(
                                        mean[, -c(1:2)]
                                ))))
                        colnames(anno) <- c('mz', 'time', name)
                        return(anno)
                } else if (peaklist) {
                        result <- data.frame(t(mean[, -c(1:2)]))
                        data <- rbind(group = name, result)
                        utils::write.csv(data, file = file)
                } else{
                        return(report)
                }
        }

#' output the similarity of two dataset
#' @param xset1 the first dataset
#' @param xset2 the second dateset
#' @return similarity
#' @export
getsim <- function(xset1,xset2){
        data1 <- gettechrep(xset1)[,c('mzmed','rtmed','rsd')]
        data2 <- gettechrep(xset2)[,c('mzmed','rtmed','rsd')]

        # data1$weight <- ifelse(data1$rsd > 100, 0, 0.2)
        # data1$weight[data1$rsd < 80] <- 0.4
        # data1$weight[data1$rsd < 60] <- 0.6
        # data1$weight[data1$rsd < 40] <- 0.8
        # data1$weight[data1$rsd < 20] <- 1
        # data1$rtorder <- order(data1$rtmed)
        data1$mzmedn <- round(data1$mzmed,0.1)

        # data2$weight <- ifelse(data2$rsd > 100, 0, 0.2)
        # data2$weight[data2$rsd < 80] <- 0.4
        # data2$weight[data2$rsd < 60] <- 0.6
        # data2$weight[data2$rsd < 40] <- 0.8
        # data2$weight[data2$rsd < 20] <- 1
        # data2$rtorder <- order(data2$rtmed)
        data2$mzmedn <- round(data2$mzmed,0.1)

        data <- merge(data1, data2, by = 'mzmedn')
        data <- data[complete.cases(data),]
        cor1 <- cor(data$rtmed.x, data$rtmed.y)
        cor2 <- cor(data$rsd.x, data$rsd.y)
        cor <- c(cor1,cor2)
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
