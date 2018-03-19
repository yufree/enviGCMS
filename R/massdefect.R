#' Get mass defect with certain scaled factor
#' @param mass vector of mass
#' @param sf scaled factors
#' @return dataframe with mass, scaled mass and scaled mass defect
#' @examples
#' mass <- c(100.1022,245.2122,267.3144,400.1222,707.2294)
#' sf <- 0.9988
#' mf <- getmassdefect(mass,sf)
#' @seealso \code{\link{getpaired}},\code{\link{plotkms}}
#' @export

getmassdefect <- function(mass, sf) {
    sm <- mass * sf
    sd <- ceiling(sm) - sm
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
#' Filter ions/peaks based on retention time hierarchical clustering, paired mass differences(PMD), PMD frequency analysis and correlation coefficient.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in cluster, default 9
#' @param isocutoff cutoff to find the isotope relationship, default 3
#' @param freqcutoff cutoff of the mass differences frequency, default 20
#' @param corcutoff cutoff of the correlation coefficient, default 0.8
#' @param rtpeaks cutoff of the max peaks within each retention time cluster, default 100
#' @return list with tentative isotope, adducts, and neutral loss peaks' index, retention time cluster.
#' @seealso \code{\link{getmassdefect}},\code{\link{plotkms}},\code{\link{getcorstd}}
#' @export
getpaired <- function(list, rtcutoff = 9, isocutoff = 3, freqcutoff = 20, corcutoff = 0.8,rtpeaks = 100){
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
                        # avoid groups with peaks more than rtpeaks and use the max rtpeaks peaks
                        if(nrow(bin)>rtpeaks){
                                bin <- bin[order(bin$mz,decreasing=TRUE)[1:rtpeaks],]
                        }
                        # get mz diff
                        dis <- stats::dist(bin$mz, method = "manhattan")
                        # get intensity cor
                        cor <- stats::cor(t(bin[,-c(1,2)]))

                        df <- data.frame(ms1 = bin$mz[which(lower.tri(dis),arr.ind = T)[, 1]], ms2 = bin$mz[which(lower.tri(dis),arr.ind = T)[, 2]], diff = as.numeric(dis), rt = medianrtxi, rtg = i, cor = cor[lower.tri(cor)])

                        dfiso <- df[df$diff<isocutoff,]
                        if(nrow(dfiso)>0){
                                resultiso <- rbind(resultiso,dfiso)
                        }
                        dfdiff <- df[df$diff>=isocutoff,]

                        result <- rbind(result,dfdiff)
                }else{
                        solo <- cbind(bin,rtg = i, cor = 1)
                        resultsolo <- rbind(solo,resultsolo)
                }
        }

        result$diff2 <- round(result$diff,2)
        if(nrow(result)>0){
                # get the high freq and high correlated ions pairs
                freq <- table(result$diff2)[order(table(result$diff2),decreasing = T)]
                resultdiff <- result[result$diff2 %in% as.numeric(names(freq[freq>freqcutoff])),]
                resultdiff <- resultdiff[resultdiff$cor > corcutoff,]
        }

        resultiso <- resultiso[resultiso$cor>corcutoff,]

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
#' @return list with std mass index
#' @seealso \code{\link{getpaired}},\code{\link{getmdg}}
#' @export
getstd <- function(list){
        resultdiffstd <- NULL
        # filter high freq ions and find std mass
        resultdiff <- list$paired
        resultsolo <- cbind.data.frame(mz = list$solo$mz, rt = list$solo$rt, rtg = list$solo$rtg)
        resultiso <- list$iso
        n <- unique(resultdiff$rtg)

        # filter the Mode mass from mass pairs within retention time group
        for(i in 1:length(n)){
                df <- resultdiff[resultdiff$rtg == n[i],]
                # massstd <- Mode(c(df$ms1,df$ms2))
                massstd <- min(c(df$ms1[which.max(df$cor)],df$ms2[which.max(df$cor)]))
                suppressWarnings(resultdiffstdtemp <- cbind(mz = c(massstd), rt = df$rt, rtg = df$rtg))
                resultdiffstd <- rbind(resultdiffstd,resultdiffstdtemp)
        }
        # filter the isotope mass from mass pairs within retention time group
        m <- unique(resultiso$rtg)
        for(i in 1:length(m)){
                df <- resultiso[resultiso$rtg == m[i],]
                massstd <- apply(df,1,function(x) min(x[1],x[2]))
                suppressWarnings(resultdiffstdtemp2 <- cbind(mz = c(massstd), rt = df$rt, rtg = df$rtg))
                resultdiffstd <- rbind(resultdiffstd,resultdiffstdtemp2)
        }
        # Combine the peaks from rt groups with single ion
        resultstd <- rbind(resultdiffstd,resultsolo)
        resultstd <- unique(resultstd)
        # show message about std mass
        n <- nrow(resultstd)
        message(paste(n, 'std mass found.'))
        # return the data
        list$stdmassindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(round(resultstd$mz,4), resultstd$rtg)
        list$stdmass <- resultstd
        return(list)
}

#' Group peaks by the mass defect hierarchical clustering group and interval group for different substructures
#' @param list a peaks list with mass to charge
#' @param submass mass vector of sub structure of homologous series
#' @param mdcutoff mass defect cluster cutoff
#' @param mdgn mass defect groups numbers for interval, 50 means 0.02 inteval on mass defect scale from -0.5 to 0.5
#' @param lv group info for the data
#' @return list with mass defect analysis dataframe.
#' @seealso \code{\link{getpaired}},\code{\link{getstd}}
#' @export
getmdg <- function(list, submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdcutoff = 0.02, mdgn = 50, lv = NULL){
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
                dis <- stats::dist(msdefect, method = "manhattan")
                fit <- stats::hclust(dis)
                mdcluster <- stats::cutree(fit, h=mdcutoff)
                n <- length(unique(mdcluster))
                message(paste(n, 'mass defect clusters found for mass', submass[i], 'substructures' ))
                index1 <- mdcluster %in% Mode(mdcluster)

                mdg <- cut(msdefect, seq(from = -.5, to = .5, by = 1/mdgn),include.lowest = T)
                mdg2 <- mdg[!is.na(mdg)]
                index2 <- mdg %in% Mode(mdg2)

                name <- c(submass[i],paste0(submass[i],'g'),paste0(submass[i],'gi'),'majormdc', 'majormdg')
                md <- cbind.data.frame(msdefect,mdcluster,mdg,index1,index2)
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
#' @param isocutoff cutoff to find the isotope relationship
#' @param freqcutoff cutoff of the mass differences frequency
#' @param corcutoff cutoff of the correlation coefficient, default 0.8
#' @param submass mass vector of sub structure of homologous series
#' @param mdcutoff mass defect cluster cutoff
#' @param mdgn mass defect groups numbers for interval, 50 means 0.02 inteval on mass defect scale from -0.5 to 0.5
#' @param lv group info for the data
#' @return list with GlobalStd algorithm processed data.
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getmdg}}
#' @export
globalstd <- function(list, rtcutoff = 9, isocutoff = 3, freqcutoff = 20, corcutoff = 0.8,submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdcutoff = 0.02, mdgn = 50, lv = NULL){
        list <- getpaired(list, rtcutoff = rtcutoff, isocutoff = isocutoff, freqcutoff = freqcutoff,corcutoff = corcutoff)
        list2 <- getstd(list)
        list3 <- getmdg(list2, submass = submass, mdcutoff = mdcutoff, mdgn = mdgn, lv = NULL)
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
        par(mfrow = c(1,2),mar = c(4,4,2,1)+0.1)
        plot(range(paired$rt),range(paired$ms1,paired$ms2),type = 'n', xlab = 'retention time(s)', ylab = 'm/z')
        segments(paired$rt,paired$ms1,paired$rt,paired$ms2,col = col[diffgroup],lwd = 1.5)
        barplot(table(list$paired$diff2),col = col[unique(diffgroup)],ylab = 'Frequency', las=2, xlab = 'mass differences')
}

#' Plot the std mass from GlobalStd algorithm
#' @param list a list from getstd function
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}}
#' @export
plotstd <- function(list){
        std <- list$stdmass
        col <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu"))))(length(unique(paired$diff2)))
        par(mfrow = c(1,2),mar = c(4,4,2,1)+0.1)
        col <- grDevices::rgb(0,0,1, alpha = 0.318)
        plot(list$rt,list$mz,xlab = 'retention time(s)', ylab = 'm/z', pch = 19, col =col,main = 'all peaks')
        plot(std$rt,std$mz,xlab = 'retention time(s)', ylab = 'm/z', pch = 19, col =col,main = 'GlobalStd peaks')
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
        plot(mz,msdata,type = 'h',xlab = 'm/z', ylab = 'Intensity')
        stdmz <- list$stdmass$mz[list$stdmass$rtg == rtcluster]
        index <- round(mz,4) %in% round(stdmz,4)

        points(mz[index],msdata[index],type = 'h',lwd = 2, col = 'red')
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
        plot(mda$rt[mdindex],mda$mz[mdindex],xlab = 'retention time(s)', ylab = 'm/z', pch = 19, col = 'blue')
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
