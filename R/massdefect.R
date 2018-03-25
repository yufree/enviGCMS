#' Get mass defect with certain scaled factor
#' @param mass vector of mass
#' @param sf scaled factors
#' @return dataframe with mass, scaled mass and scaled mass defect
#' @examples
#' mass <- c(100.1022,245.2122,267.3144,400.1222,707.2294)
#' sf <- 0.9988
#' mf <- getmassdefect(mass,sf)
#' @seealso \code{\link{globalstd}},\code{\link{plotkms}}
#' @export

getmassdefect <- function(mass, sf) {
    sm <- mass * sf
    sd <- round(sm) - sm
    df <- as.data.frame(cbind(mass, sm, sd))
    graphics::plot(df$sd ~ df$sm, xlab = "m/z", ylab = "scaled MD")
    return(df)
}
#' define the Mode function
#' @param x vector
#' @return Mode of the vector
#' @export
Mode = function(x){
        ta = table(x)
        tam = max(ta)
        if (all(ta == tam))
                mod = x
        else
                if(is.numeric(x))
                        mod = as.numeric(names(ta)[ta == tam])
        else
                mod = names(ta)[ta == tam]
        return(mod)
}
#' Filter ions/peaks based on retention time hierarchical clustering, paired mass differences(PMD) and PMD frequency analysis.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in cluster, default 9
#' @param freqcutoff cutoff of the mass differences frequency, default 20
#' @return list with tentative isotope, adducts, and neutral loss peaks' index, retention time clusters.
#' @seealso \code{\link{getmassdefect}},\code{\link{getstd}},\code{\link{getstd}},\code{\link{getmdg}},\code{\link{plotpaired}}
#' @export
getpaired <- function(list, rtcutoff = 9, freqcutoff = 20){
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
        }

        # filter the list
        # get the rt cluster
        list$rtcluster <- rtcluster
        # get the data index by rt groups with single ions
        list$soloindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(round(resultsolo$mz,4), resultsolo$rtg)
        list$solo <- resultsolo
        # get the data index by rt groups with isotope ions
        list$isoindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(c(round(resultiso$ms1,4),round(resultiso$ms2,4)), c(resultiso$rtg,resultiso$rtg))
        list$iso <- resultiso
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
#' @seealso \code{\link{getpaired}},\code{\link{getmdg}},\code{\link{plotstd}}
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
        resultstd1 <- cbind(list$solo$mz,list$solo$rt,list$solo$rtg)
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
                df <- resultiso[resultiso$rtg == rtg2B1[i],]
                if(nrow(df)>0){
                mass <- apply(df,1,function(x) min(x[1],x[2]))
                mass <- unique(mass)
                suppressWarnings(resultstdtemp <- cbind(mz = c(mass), rt = df$rt, rtg = df$rtg))
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
                        massstd <- apply(dfiso,1,function(x) min(x[1],x[2]))
                        massstd <- unique(massstd)
                        dis <- stats::dist(massstd, method = "manhattan")
                        df <- data.frame(ms1 = massstd[which(lower.tri(dis),arr.ind = T)[, 1]], ms2 = massstd[which(lower.tri(dis),arr.ind = T)[, 2]], diff = round(as.numeric(dis),2))
                        # remove the adducts
                        if(sum((df$diff %in% dfpaired$diff2))>0){
                                mass <- unique(df$ms1[df$diff %in% dfpaired$diff2])
                                massadd <- unique(c(df$ms1[df$diff %in% dfpaired$diff2],df$ms2[df$diff %in% dfpaired$diff2]))
                                massextra <- massstd[!(massstd %in% massadd)]
                                mass <- c(massextra,mass)
                        }else{
                                mass <- massstd
                        }
                        suppressWarnings(resultstdtemp <- cbind(mz = c(mass), rt = dfiso$rt, rtg = dfiso$rtg))
                        suppressWarnings(resultstd2B3 <- rbind(resultstd2B3,resultstdtemp))
                }else if(nrow(dfiso)>0){
                        mass <- apply(dfiso,1,function(x) min(x[1],x[2]))
                        mass <- unique(mass)
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

#' Group peaks by the mass defect interval group for different substructures
#' @param list a peaks list with mass to charge
#' @param submass mass vector of sub structure of homologous series
#' @param mdgn mass defect groups numbers for interval, 20 means 0.05 inteval on mass defect scale from -0.5 to 0.5
#' @param lv group info for the data
#' @return list with mass defect analysis dataframe.
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{plotstd}},\code{\link{plotstdmd}},\code{\link{plotstdrt}}
#' @export
getmdg <- function(list, submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdgn = 20, lv = NULL){
        if(is.null(list$stdmass)&is.null(list$paired)){
                mz <- list$mz
                rt <- list$rt
                data <- list$data
                colnames(data) <- lv
        }else if(is.null(list$stdmass)){
                mz <- list$mz[list$pairedindex]
                rt <- list$rt[list$pairedindex]
                data <- list$data[list$pairedindex,]
                colnames(data) <- lv
        }else{
                mz <- list$mz[list$stdmassindex]
                rt <- list$rt[list$stdmassindex]
                data <- list$data[list$stdmassindex,]
                colnames(data) <- lv
        }
        # perform mass defect analysis for std mass
        mda <- cbind.data.frame(mz = mz, rt = rt, data)
        for(i in 1:length(submass)){
                mdst <- round(submass[i])/submass[i]
                msdefect <- round(mz*mdst) - mz*mdst

                mdg <- cut(msdefect, seq(from = -.5, to = .5, by = 1/mdgn),include.lowest = T)
                mdg2 <- mdg[!is.na(mdg)]
                index <- mdg %in% Mode(mdg2)

                name <- c(submass[i],paste0(submass[i],'gi'), 'majormdg')
                md <- cbind.data.frame(msdefect,mdg,index)
                colnames(md) <- name
                mda <- cbind.data.frame(mda,md)
        }
        # get the data
        list$mda <- mda
        return(list)
}
#' GlobalStd algorithm
#' @param list a peaks list with mass to charge
#' @param rtcutoff cutoff of the distances in cluster
#' @param freqcutoff cutoff of the mass differences frequency
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param submass mass vector of sub structure of homologous series
#' @param mdgn mass defect groups numbers for interval, 20 means 0.05 inteval on mass defect scale from -0.5 to 0.5
#' @param lv group info for the data
#' @return list with GlobalStd algorithm processed data.
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getmdg}},\code{\link{plotstd}},\code{\link{plotstdmd}},\code{\link{plotstdrt}}
#' @export
globalstd <- function(list, rtcutoff = 9, freqcutoff = 20, corcutoff = NULL,submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834),  mdgn = 20, lv = NULL){
        list <- getpaired(list, rtcutoff = rtcutoff, freqcutoff = freqcutoff)
        list2 <- getstd(list,corcutoff = corcutoff)
        list3 <- getmdg(list2, submass = submass, mdgn = mdgn, lv = NULL)
        return(list3)
}
#' Paired correlationship among peak list based on cluster analysis
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in cluster
#' @param submass mass vector of sub structure of homologous series
#' @param mdcutoff mass defect cluster cutoff
#' @return list with retention time cluster, std mass defect analysis dataframe based on max average correlation
#' @seealso \code{\link{getmassdefect}},\code{\link{plotkms}},\code{\link{getpaired}}
getcorstd <- function(list, rtcutoff = 9, submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdcutoff = 0.02){

        # paired mass diff analysis
        groups <- cbind.data.frame(mz = list$mz, rt = list$rt, list$data)
        resultstd <- resultsolo <- resultiso <- result <- NULL

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
                        cor <- stats::cor(t(bin[,-c(1,2)]))
                        cormean <- apply(cor,1,mean)
                        corindex <- which.max(cormean)
                        df <- cbind(bin[corindex,],rtg = i)
                        result <- rbind(result,df)
                }else{
                        solo <- cbind(bin,rtg = i)
                        resultsolo <- rbind(solo,resultsolo)
                }
        }

        resultstd <- rbind(result,resultsolo)
        resultstd <- unique(resultstd)

        # perform mass defect analysis for std mass
        for(i in 1:length(submass)){
                mdst <- round(submass[i])/submass[i]
                msdefect <- round(resultstd$mz*mdst) - resultstd$mz*mdst
                dis <- stats::dist(msdefect, method = "manhattan")
                fit <- stats::hclust(dis)
                mdcluster <- stats::cutree(fit, h=mdcutoff)
                n <- length(unique(mdcluster))
                message(paste(n, 'mass defect clusters found for mass', submass[i], 'substructures' ))
                name <- c(colnames(resultstd),submass[i],paste0(submass[i],'g'))
                resultstd <- cbind.data.frame(resultstd,msdefect,mdcluster)
                colnames(resultstd) <- name
        }

        # filter the list
        list$rtcluster <- rtcluster
        list$stdmassindex <- (round(list$mz,4) %in% round(resultstd$mz,4))
        list$stdmass <- resultstd
        return(list)
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
        graphics::barplot(table(list$paired$diff2),col = col[unique(diffgroup)],ylab = 'Frequency', las=2, xlab = 'mass differences')
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
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}},\code{\link{plotstd}},\code{\link{plotpaired}},\code{\link{plotstdmd}}
#' @export
#'
plotstdrt <- function(list,rtcluster){
        data <- list$data[list$rtcluster == rtcluster,]
        if(length(data)>ncol(list$data)){
                msdata <- apply(data,1,mean)
        }else{
                msdata <- mean(data)
        }
        mz <- list$mz[list$rtcluster == rtcluster]
        graphics::plot(mz,msdata,type = 'h',xlab = 'm/z', ylab = 'Intensity')
        stdmz <- list$stdmass$mz[list$stdmass$rtg == rtcluster]
        index <- round(mz,4) %in% round(stdmz,4)

        graphics::points(mz[index],msdata[index],type = 'h',lwd = 2, col = 'red')
}

#' Plot the std mass from GlobalStd algorithm in certain mass defect groups
#' @param list a list from getmdg function
#' @param mdindex mass defect index
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}},\code{\link{plotstd}},\code{\link{plotpaired}},\code{\link{plotstdrt}}
#' @export
#'
plotstdmd <- function(list,mdindex){
        mda <- list$mda
        graphics::plot(mda$rt[mdindex],mda$mz[mdindex],
             xlab = 'retention time(s)', ylab = 'm/z',
             pch = 19, col = 'blue',
             xlim = range(mda$rt),ylim = range(mda$mz))
}
#' plot the kendrick mass defect diagram
#' @param data vector with the name m/z
#' @param cutoff remove the low intensity
#' @return NULL
#' @seealso \code{\link{getmassdefect}},\code{\link{getpaired}}
#' @examples
#' \dontrun{
#' mz <- c(10000,5000,20000,100,40000)
#' names(mz) <- c(100.1022,245.2122,267.3144,400.1222,707.2294)
#' plotkms(mz)
#' }
#' @export
plotkms <- function(data, cutoff = 1000) {
    data <- data[data > cutoff]
    mz <- as.numeric(names(data))
    km <- mz * 14/14.01565
    kmd <- round(km) - km
    graphics::smoothScatter(kmd ~ round(km), xlab = "Kendrick nominal mass",
        ylab = "Kendrick mass defect")
}
