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

#' Paired mass diff relationship among peak list based on cluster analysis
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in cluster
#' @param isocutoff cutoff to find the isotope relationship
#' @param freqcutoff cutoff of the mass differences frequency
#' @param submass mass vector of sub structure of homologous series
#' @param mdcutoff mass defect cluster cutoff
#' @return list with tentative isotope, adducts, and neutral loss peaks' index, retention time cluster, std mass defect analysis dataframe. The isotope index might be less than the table and the mass pairs are the reason.
#' @seealso \code{\link{getmassdefect}},\code{\link{plotkms}},\code{\link{getcorstd}}
#' @export
getpaired <- function(list, rtcutoff = 9, isocutoff = 3, freqcutoff = 20, submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdcutoff = 0.02){

        # paired mass diff analysis
        groups <- cbind.data.frame(mz = list$mz, rt = list$rt)
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
                        df <- data.frame(ms1 = bin$mz[which(lower.tri(dis),arr.ind = T)[, 1]], ms2 = bin$mz[which(lower.tri(dis),arr.ind = T)[, 2]], diff = as.numeric(dis), rt = medianrtxi, rtg = i)

                                dfiso <- df[df$diff<isocutoff,]
                                if(nrow(dfiso)>0){
                                        resultiso <- rbind(resultiso,dfiso)
                                }
                                dfdiff <- df[df$diff>=isocutoff,]
                                result <- rbind(result,dfdiff)
                }else{
                        solo <- cbind(bin,rtg = i)
                        resultsolo <- rbind(solo,resultsolo)
                }
        }


        result$diff2 <- round(result$diff,2)
        if(nrow(result)>0){
                # get the high freq ions pair
                freq <- table(result$diff2)[order(table(result$diff2),decreasing = T)]
                resultdiff <- result[result$diff2 %in% as.numeric(names(freq[freq>freqcutoff])),]
        }
        # filter high freq ions and find std mass
        n <- unique(resultdiff$rtg)
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

        for(i in 1:length(n)){
                df <- resultdiff[resultdiff$rtg == n[i],]
                massstd <- Mode(c(df$ms1,df$ms2))
                suppressWarnings(resultdiffstdtemp <- cbind(mz = massstd, rt = df$rt, rtg = df$rtg))
                resultdiffstd <- rbind(resultdiffstd,resultdiffstdtemp)
        }

        resultstd <- rbind(resultdiffstd,resultsolo)
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

        list$soloindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(round(resultsolo$mz,4), resultsolo$rtg)
        list$solo <- resultsolo

        list$diffindex <- list$mz %in% c(result$ms1,result$ms2)
        list$diff <- result

        list$paired <- resultdiff

        list$isoindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(c(round(resultiso$ms1,4),round(resultiso$ms2,4)), c(resultiso$rtg,resultiso$rtg))
        list$iso <- resultiso

        list$stdmassindex <- paste(round(list$mz,4),list$rtcluster) %in% paste(round(resultstd$mz,4), resultstd$rtg)
        list$stdmass <- resultstd
        return(list)
}

#' Paired correlationship among peak list based on cluster analysis
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in cluster
#' @param submass mass vector of sub structure of homologous series
#' @param mdcutoff mass defect cluster cutoff
#' @return list with retention time cluster, std mass defect analysis dataframe based on max average correlation
#' @seealso \code{\link{getmassdefect}},\code{\link{plotkms}},\code{\link{getpaired}}
#' @export
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
