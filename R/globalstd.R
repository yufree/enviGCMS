#' Filter ions/peaks based on retention time hierarchical clustering, paired mass differences(PMD) and PMD frequency analysis.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 9
#' @param freqcutoff cutoff of the paired mass difference frequency, default 30
#' @param pmdcutoff cutoff of the largest paired mass differences, default 100
#' @return list with tentative isotope, adducts, and neutral loss peaks' index, retention time clusters.
#' @seealso \code{\link{getstd}},\code{\link{getsda}},\code{\link{plotpaired}}
#' @export
getpaired <- function(list, rtcutoff = 9, freqcutoff = 30, pmdcutoff = 100){
        # paired mass diff analysis
        groups <- cbind.data.frame(mz = list$mz, rt = list$rt, list$data)
        resultstd <- resultdiffstd <- resultsolo <- resultiso <- result <- NULL

        dis <- stats::dist(list$rt, method = "manhattan")
        fit <- stats::hclust(dis)
        rtcluster <- stats::cutree(fit, h=rtcutoff)
        n <- length(unique(rtcluster))
        message(paste(n, 'retention time cluster found.'))
        # search:
        for (i in 1:length(unique(rtcluster))) {
                # find the mass within RT
                rtxi <- list$rt[rtcluster == i]
                bin = groups[groups$rt %in% rtxi, ]
                medianrtxi <- stats::median(rtxi)

                if (nrow(bin) > 1) {
                        # get mz diff
                        dis <- stats::dist(bin$mz, method = "manhattan")
                        # get intensity cor
                        cor <- stats::cor(t(bin[,-c(1,2)]))

                        df <- data.frame(ms1 = bin$mz[which(lower.tri(dis),arr.ind = T)[, 1]], ms2 = bin$mz[which(lower.tri(dis),arr.ind = T)[, 2]], diff = as.numeric(dis), rt = medianrtxi, rtg = i, cor = cor[lower.tri(cor)])

                        isoindex <- (df$diff<1.01 & df$diff>0.99)|(df$diff<2.01 & df$diff>1.99)
                        dfiso <- df[isoindex,]
                        if(nrow(dfiso)>0){
                                resultiso <- rbind(resultiso,dfiso)
                        }
                        dfdiff <- df[!isoindex,]

                        result <- rbind(result,dfdiff)
                }else{
                        solo <- cbind(bin,rtg = i, cor = 1)
                        resultsolo <- rbind(solo,resultsolo)
                }
        }

        # # mass defect filter
        # mdst <- round(14.01565)/14.01565
        # msdefect <- round(result$diff*mdst) - result$diff*mdst
        # index <- abs(msdefect)<0.05
        # result <- result[index,]
        # get the high freq ions pairs
        result$diff2 <- round(result$diff,2)
        if(nrow(result)>0){
                freq <- table(result$diff2)[order(table(result$diff2),decreasing = T)]
                resultdiff <- result[result$diff2 %in% as.numeric(names(freq[freq>freqcutoff])),]
                resultdiff <- resultdiff[resultdiff$diff2<pmdcutoff,]
        }

        # filter the list
        # get the rt cluster
        list$rtcluster <- rtcluster
        # get the data index by rt groups with single ions
        if(!is.null(resultsolo)){
                list$soloindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(round(resultsolo$mz,4), resultsolo$rtg)
                list$solo <- resultsolo
        }
        # get the data index by rt groups with isotope ions
        if(!is.null(resultiso)){
        list$isoindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(c(round(resultiso$ms1,4),round(resultiso$ms2,4)), c(resultiso$rtg,resultiso$rtg))
        list$iso <- resultiso
        }
        # get the data index by rt groups with high freqences PMD
        list$diffindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(c(round(result$ms1,4),round(result$ms2,4)), c(result$rtg,result$rtg))
        list$diff <- result

        list$pairedindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(c(round(resultdiff$ms1,4),round(resultdiff$ms2,4)), c(resultdiff$rtg,resultdiff$rtg))
        list$paired <- resultdiff

        # show message about std mass
        n <- sum(list$pairedindex)
        message(paste(n, 'paired mass found.'))
        # return results
        return(list)
}
#' Find the standard mass for each retention time hierarchical clustering based on PMD relationship within each retention time cluster and isotope and return the index of the std data for each retention time cluster.
#' @param list a list from getpaired function
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @return list with std mass index
#' @seealso \code{\link{getpaired}},\code{\link{getsda}},\code{\link{plotstd}}
#' @export
getstd <- function(list, corcutoff = NULL){
        resultstd2A <- resultstd2B1 <- resultstd2B2 <- resultstd2B3 <- NULL
        # filter high freq ions and find std mass
        resultdiff <- list$paired
        resultiso <- list$iso

        if(!is.null(corcutoff)){
                resultdiff <- resultdiff[resultdiff$cor > corcutoff,]
                resultiso <- resultiso[resultiso$cor > corcutoff,]
        }
        # filter the mass from mass pairs within retention time group
        # group 1: RT groups with solo peak
        resultstd1 <- NULL
        if(!is.null(list$solo)){
        resultstd1 <- cbind(list$solo$mz,list$solo$rt,list$solo$rtg)
        }
        # group 2: RT groups with multiple peaks
        # group 2A: RT groups with multiple peaks while no isotope/paired relationship
        index2A <- !(unique(list$rtcluster) %in% unique(resultdiff$rtg)|unique(list$rtcluster) %in% unique(resultiso$rtg))
        rtg2A <- unique(list$rtcluster)[index2A]
        for(i in 1:length(rtg2A)){
                mass <- list$mz[list$rtcluster == rtg2A[i]]
                rt <- list$rt[list$rtcluster == rtg2A[i]]
                mass <- max(mass)
                suppressWarnings(resultstdtemp <- c(mass, stats::median(rt), rtg2A[i]))
                suppressWarnings(resultstd2A <- rbind(resultstd2A,resultstdtemp))
        }
        # group 2B: RT groups with multiple peaks with isotope/paired relationship
        # index2B <- (unique(list$rtcluster) %in% unique(resultdiff$rtg))|(unique(list$rtcluster) %in% unique(resultiso$rtg))
        # group 2B1: RT groups with multiple peaks with isotope without paired relationship
        index2B1 <- !(unique(list$rtcluster) %in% unique(resultdiff$rtg))&unique(list$rtcluster) %in% unique(resultiso$rtg)
        rtg2B1 <- unique(list$rtcluster)[index2B1]
        for(i in 1:length(rtg2B1)){
                # filter the isotope peaks
                dfiso <- resultiso[resultiso$rtg == rtg2B1[i],]
                if(nrow(dfiso)>0){
                        massstd <- apply(dfiso,1,function(x) min(x[1],x[2]))
                        massstdmax <- apply(dfiso,1,function(x) max(x[1],x[2]))
                        mass <- unique(massstd[!(massstd %in% massstdmax)])
                suppressWarnings(resultstdtemp <- cbind(mz = c(mass), rt = dfiso$rt, rtg = dfiso$rtg))
                resultstd2B1 <- rbind(resultstd2B1,resultstdtemp)
                }
        }
        # group 2B2: RT groups with multiple peaks with paired relationship without isotope
        index2B2 <- (unique(list$rtcluster) %in% unique(resultdiff$rtg))&!(unique(list$rtcluster) %in% unique(resultiso$rtg))
        rtg2B2 <- unique(list$rtcluster)[index2B2]
        for(i in 1:length(rtg2B2)){
                # filter the paired peaks
                df <- resultdiff[resultdiff$rtg == rtg2B2[i],]
                if(nrow(df)>0){
                        mass <- apply(df,1,function(x) min(x[1],x[2]))
                        mass <- unique(mass)
                        suppressWarnings(resultstdtemp <- cbind(mz = c(mass), rt = df$rt, rtg = df$rtg))
                        suppressWarnings(resultstd2B2 <- rbind(resultstd2B2,resultstdtemp))
                }
        }
        # group 2B3: RT groups with multiple peaks with paired relationship and isotope
        index2B3 <- (unique(list$rtcluster) %in% unique(resultdiff$rtg))&(unique(list$rtcluster) %in% unique(resultiso$rtg))
        rtg2B3 <- unique(list$rtcluster)[index2B3]

        for(i in 1:length(rtg2B3)){
                # filter the isotope peaks
                dfiso <- resultiso[resultiso$rtg == rtg2B3[i],]
                dfpaired <- resultdiff[resultdiff$rtg == rtg2B3[i],]
                if(nrow(dfiso)>0&nrow(dfpaired)>0){
                        # remove peaks with more than one isotopes
                        massstd <- apply(dfiso,1,function(x) min(x[1],x[2]))
                        massstdmax <- apply(dfiso,1,function(x) max(x[1],x[2]))
                        massstd <- unique(massstd[!(massstd %in% massstdmax)])
                        dis <- stats::dist(massstd, method = "manhattan")
                        df <- data.frame(ms1 = massstd[which(lower.tri(dis),arr.ind = T)[, 1]], ms2 = massstd[which(lower.tri(dis),arr.ind = T)[, 2]], diff = round(as.numeric(dis),2))
                        # remove the adducts
                        if(sum((df$diff %in% dfpaired$diff2))>0){
                                massstd <- unique(apply(df[df$diff %in% dfpaired$diff2,],1,function(x) min(x[1],x[2])))
                                massused <- unique(c(df$ms1,df$ms2))

                                massadd <- unique(c(df$ms1[df$diff %in% dfpaired$diff2],df$ms2[df$diff %in% dfpaired$diff2]))
                                massextra <- massused[!(massused %in% massadd)]
                                mass <- c(massextra,massstd)
                        }else{
                                mass <- massstd
                        }
                        suppressWarnings(resultstdtemp <- cbind(mz = c(mass), rt = dfiso$rt, rtg = dfiso$rtg))
                        suppressWarnings(resultstd2B3 <- rbind(resultstd2B3,resultstdtemp))
                }else if(nrow(dfiso)>0){
                        # remove peaks with more than one peaks
                        massstd <- apply(dfiso,1,function(x) min(x[1],x[2]))
                        massstdmax <- apply(dfiso,1,function(x) max(x[1],x[2]))
                        mass <- unique(massstd[!(massstd %in% massstdmax)])
                        suppressWarnings(resultstdtemp <- cbind(mz = c(mass), rt = dfiso$rt, rtg = dfiso$rtg))
                        suppressWarnings(resultstd2B3 <- rbind(resultstd2B3,resultstdtemp))
                }else{
                        mass <- apply(dfpaired,1,function(x) min(x[1],x[2]))
                        mass <- unique(mass)
                        suppressWarnings(resultstdtemp <- cbind(mz = c(mass), rt = dfpaired$rt, rtg = dfpaired$rtg))
                        suppressWarnings(resultstd2B3 <- rbind(resultstd2B3,resultstdtemp))
                }

        }

        # Combine the peaks from rt groups with single ion
        resultstd <- rbind(resultstd1,resultstd2A,resultstd2B3,resultstd2B2,resultstd2B1)
        resultstd <- unique(resultstd)
        colnames(resultstd) <- c('mz','rt','rtg')
        resultstd <- as.data.frame(resultstd)
        # show message about std mass
        n <- nrow(resultstd)
        message(paste(n, 'std mass found.'))
        # return the data
        list$stdmassindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(round(resultstd$mz,4), resultstd$rtg)
        list$stdmass <- resultstd
        return(list)
}

#' Perform structure directed analysis for peaks list.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the retention time in secounds, default 3
#' @param freqcutoff cutoff of the paired mass difference frequency, default 10
#' @param pmdcutoff cutoff of the largest paired mass differences, default 500
#' @return list with tentative isotope, adducts, and neutral loss peaks' index, retention time clusters.
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{plotpaired}}
getsda <- function(list, rtcutoff = 3, freqcutoff = 10, pmdcutoff = 500){
        if(is.null(list$stdmass)&is.null(list$paired)){
                mz <- list$mz
                rt <- list$rt
                data <- list$data
        }else if(is.null(list$stdmass)){
                mz <- list$mz[list$pairedindex]
                rt <- list$rt[list$pairedindex]
                data <- list$data[list$pairedindex,]
        }else{
                mz <- list$mz[list$stdmassindex]
                rt <- list$rt[list$stdmassindex]
                data <- list$data[list$stdmassindex,]
        }
        # PMD analysis
        dis <- stats::dist(mz, method = "manhattan")
        disrt <- stats::dist(rt, method = "manhattan")
        df <- data.frame(ms1 = mz[which(lower.tri(dis),arr.ind = T)[, 1]], ms2 = mz[which(lower.tri(dis),arr.ind = T)[, 2]], diff = as.numeric(dis), rt1 = rt[which(lower.tri(disrt),arr.ind = T)[, 1]],rt2 = rt[which(lower.tri(disrt),arr.ind = T)[, 2]],diffrt = as.numeric(disrt))
        df$diff2 <- round(df$diff,2)
        df <- df[df$diffrt>rtcutoff&df$diff2<pmdcutoff,]
        freq <- table(df$diff2)[order(table(df$diff2),decreasing = T)]
        list$sda <- df[(df$diff2 %in% as.numeric(names(freq[freq>=freqcutoff]))),]

        # show message about std mass
        sub <- names(table(list$sda$diff2))
        n <- length(sub)
        message(paste(n, 'groups were found as high frequency PMD group.','\n'))
        message(paste(sub, 'were found as high frequency PMD.','\n'))
        return(list)
}

#' GlobalStd algorithm
#' @param list a peaks list with mass to charge
#' @param rtcutoff cutoff of the distances in cluster
#' @param freqcutoff cutoff of the mass differences frequency
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param pmdcutoff cutoff of the largest mass differences, default 100
#' @param rtcutoff2 cutoff of the retention time in secounds, default 3
#' @param freqcutoff2 cutoff of the paired mass difference frequency, default 10
#' @param pmdcutoff2 cutoff of the largest paired mass differences, default 500
#' @return list with GlobalStd algorithm processed data.
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{plotstd}},\code{\link{plotstdsda}},\code{\link{plotstdrt}}
#' @export
globalstd <- function(list, rtcutoff = 9, freqcutoff = 30, corcutoff = NULL,pmdcutoff = 100, rtcutoff2 = 3, freqcutoff2 = 10, pmdcutoff2 = 500){
        list <- getpaired(list, rtcutoff = rtcutoff, freqcutoff = freqcutoff, pmdcutoff = pmdcutoff)
        list2 <- getstd(list,corcutoff = corcutoff)
        list3 <- getsda(list2, rtcutoff = rtcutoff2, freqcutoff = freqcutoff, pmdcutoff = pmdcutoff2)
        return(list3)
}

#' Plot the mass pairs and high frequency mass differences
#' @param list a list from getpaired function
#' @return NULL
#' @seealso \code{\link{getpaired}}, \code{\link{globalstd}}
#' @export
plotpaired <- function(list){
        paired <- list$paired
        diffgroup <- as.numeric(as.factor(paired$diff2))
        col <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu"))))(length(unique(paired$diff2)))
        graphics::par(mfrow = c(2,1),mar = c(4,4,2,1)+0.1)
        graphics::plot(range(paired$rt),range(paired$ms1,paired$ms2),type = 'n', xlab = 'retention time(s)', ylab = 'm/z')
        graphics::segments(paired$rt,paired$ms1,paired$rt,paired$ms2,col = col[diffgroup],lwd = 1.5)
        graphics::barplot(table(list$paired$diff2),col = col,ylab = 'Frequency', las=2, xlab = 'Paired mass difference',cex.names=0.618)
}

#' Plot the std mass from GlobalStd algorithm
#' @param list a list from getstd function
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}}
#' @export
plotstd <- function(list){
        std <- list$stdmass
        graphics::par(mfrow = c(1,2),mar = c(4,4,2,1)+0.1)
        col <- grDevices::rgb(0,0,1, alpha = 0.318)
        graphics::plot(list$rt,list$mz,xlab = 'retention time(s)', ylab = 'm/z', pch = 19, col =col,main = 'all peaks')
        graphics::plot(std$rt,std$mz,xlab = 'retention time(s)', ylab = 'm/z', pch = 19, col = col, main = 'GlobalStd peaks')
}

#' Plot the std mass from GlobalStd algorithm in certain retention time groups
#' @param list a list from getstd function
#' @param rtcluster retention time group index
#' @param ... other parameters for plot function
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}},\code{\link{plotstd}},\code{\link{plotpaired}},\code{\link{plotstdmd}}
#' @export
#'
plotstdrt <- function(list,rtcluster,...){
        data <- list$data[list$rtcluster == rtcluster,]
        if(length(data)>ncol(list$data)){
                msdata <- apply(data,1,mean)
        }else{
                msdata <- mean(data)
        }
        mz <- list$mz[list$rtcluster == rtcluster]
        rt <- median(list$rt[list$rtcluster == rtcluster])
        graphics::plot(mz,msdata,type = 'h',xlab = paste('m/z','@',rt,'s'), ylab = 'Intensity',...)
        stdmz <- list$stdmass$mz[list$stdmass$rtg == rtcluster]
        index <- round(mz,4) %in% round(stdmz,4)

        graphics::points(mz[index],msdata[index],type = 'h',lwd = 2, col = 'red')
}

#' Plot the std mass from GlobalStd algorithm in structure directed analysis(SDA) groups
#' @param list a list from getsda function
#' @param index index for PMD value
#' @param ... other parameters for plot function
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}},\code{\link{plotstd}},\code{\link{plotpaired}},\code{\link{plotstdrt}}
#' @export
plotstdsda <- function(list,index = NULL,...){
        sda <- list$sda
        diffgroup <- as.numeric(as.factor(sda$diff2))
        if(is.null(index)){
                col <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu"))))(length(unique(diffgroup)))
                graphics::par(mfrow = c(2,1),mar = c(4,4,2,1)+0.1)
                graphics::plot(range(sda$rt1,sda$rt2),range(sda$ms1,sda$ms2),type = 'n', xlab = 'retention time(s)', ylab = 'm/z',...)
                graphics::segments(sda$rt1,sda$ms1,sda$rt2,sda$ms2,col = col[diffgroup],lwd = 1.5)
                points(sda$rt1,sda$ms1,pch = 19, col = grDevices::rgb(0,0,1, alpha = 0.318))
                points(sda$rt2,sda$ms2,pch = 19, col = grDevices::rgb(0,0,1, alpha = 0.318))
                graphics::barplot(table(sda$diff2),col = col,ylab = 'Frequency', las=2, xlab = 'Paired mass differences',cex.names=0.618)

        }else{
                sda <- list$sda[index,]
                graphics::plot(range(list$sda$rt1,list$sda$rt2),range(list$sda$ms1,list$sda$ms2),type = 'n', xlab = 'retention time(s)', ylab = 'm/z',main = paste(sda$diff2[1],'group'),...)
                graphics::segments(sda$rt1,sda$ms1,sda$rt2,sda$ms2,lwd = 1.5,pch = 19)
                points(sda$rt1,sda$ms1,pch = 19, col = grDevices::rgb(0,0,1, alpha = 0.318))
                points(sda$rt2,sda$ms2,pch = 19, col = grDevices::rgb(0,0,1, alpha = 0.318))
        }

}
