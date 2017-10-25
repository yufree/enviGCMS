#' Get mass defect with certain scaled factor
#' @param mass vector of mass
#' @param sf scaled factors
#' @return dataframe with mass, scaled mass and scaled mass defect
#' @examples
#' mass <- c(100.1022,245.2122,267.3144,400.1222,707.2294)
#' sf <- 0.9988
#' mf <- getmassdefect(mass,sf)
#' @seealso \code{\link{getmassdiff}},\code{\link{plotkms}}
#' @export

getmassdefect <- function(mass, sf) {
    sm <- mass * sf
    sd <- ceiling(sm) - sm
    df <- as.data.frame(cbind(mass, sm, sd))
    graphics::plot(df$sd ~ df$sm, xlab = "m/z", ylab = "scaled MD")
    return(df)
}

#' Isotope extraction for single group of samples with certain mass diff
#' @param list  a list with mzrt profile, mz and rt
#' @param massdiff mass defect
#' @param rtwindow retention time range
#' @param mzwindow mass charge ratio window
#' @param ppm resolution of the mass spectrum
#' @seealso \code{\link{getmassdefect}},\code{\link{plotkms}}
#' @return dataframe with mass, retention time, scaled mass and scaled mass defect
getmassdiff <- function(list, massdiff, rtwindow, mzwindow, 
    ppm) {
    # get intensity infomation
    peakIntensities = list$data
    groups = cbind.data.frame(mz = list$mz, rt = list$rt)
    # order peaks by rt
    peakIntensities = peakIntensities[order(list$rt), ]
    groups <- groups[order(list$rt), ]
    groups$peakins <- apply(peakIntensities, 1, mean)
    result <- NULL
    # search:
    for (i in 1:nrow(groups)) {
        bin = groups[list$rt - list$rt[i] >= 0 & list$rt - 
            list$rt[i] <= rtwindow, ]
        if (nrow(bin) > 1) {
            dis <- stats::dist(bin$mz, method = "manhattan")/massdiff
            df <- data.frame(ms1 = bin$mz[which(lower.tri(dis), 
                arr.ind = T)[, 1]], ms2 = bin$mz[which(lower.tri(dis), 
                arr.ind = T)[, 2]], diff = as.numeric(dis))
            df$rdiff <- round(df$diff)
            dfn <- df[df$diff <= df$rdiff * (1 + ppm/1e+06) + 
                (df$ms1 * ppm/1e+06)/(massdiff * (1 - ppm/1e+06)) && 
                df$diff >= df$rdiff * (1 - ppm/1e+06) - (df$ms1 * 
                  ppm/1e+06)/(massdiff * (1 + ppm/1e+06)), 
                ]
            dfn$msdiff <- abs(dfn$ms1 - dfn$ms2)
            dfn <- dfn[dfn$msdiff < mzwindow, ]
            # candidate number of labeled atoms
            result <- rbind(result, bin[bin$mz %in% dfn$ms1 | 
                bin$mz %in% dfn$ms2, ])
            result <- result[rownames(unique(result[, c("mz", 
                "rt")])), ]
        }
    }
    result$sm <- result$mz * massdiff
    result$smd <- ceiling(result$sm) - result$sm
    return(result)
}


#' plot the kendrick mass defect diagram
#' @param data vector with the name m/z
#' @param cutoff remove the low intensity
#' @return NULL
#' @seealso \code{\link{getmassdefect}},\code{\link{getmassdiff}}
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
