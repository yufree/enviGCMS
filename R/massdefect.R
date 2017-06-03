#' Get mass defect with certain scaled factor
#' @param mass vector of mass
#' @param sf scaled factors
#' @return dataframe with mass, scaled mass and scaled mass defect
getmassdefect <- function(mass, sf) {
    sm <- mass * sf
    sd <- ceiling(sm) - sm
    df <- as.data.frame(cbind(mass, sm, sd))
    graphics::plot(df$sd ~ df$sm, xlab = "m/z", ylab = "scaled MD")
    return(df)
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

#' plot the kendrick mass defect diagram
#' @param data vector with the name m/z
#' @param cutoff remove the low intensity
#' @return NULL
#' @export
plotkms <- function(data, cutoff = 1000) {
    data <- data[data > cutoff]
    mz <- as.numeric(names(data))
    km <- mz * 14/14.01565
    kmd <- round(km) - km
    graphics::smoothScatter(kmd ~ round(km), xlab = "Kendrick nominal mass", 
        ylab = "Kendrick mass defect")
}
#' to do
#' Van Krevelen diagram
