#' Impute the peaks list data
#' @param list list with data as peaks list, mz, rt and group information
#' @param method 'r' means remove, 'l' means use half the minimum of the values across the peaks list, 'mean' means mean of the values across the samples, 'median' means median of the values across the samples, '0' means 0, '1' means 1. Default 'l'.
#' @return list with imputed peaks
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' getimputation(list)
#' }
#' @export
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}},\code{\link{getmzrt2}}, \code{\link{getdoe}}, \code{\link{getmr}}
getimputation <- function(list, method = "l") {
    data <- list$data
    mz <- list$mz
    rt <- list$rt
    
    if (method == "r") {
        idx <- stats::complete.cases(data)
        data <- data[idx, ]
        mz <- mz[idx]
        rt <- rt[idx]
    } else if (method == "l") {
        impute <- min(data, na.rm = T)/2
        data[is.na(data)] <- impute
    } else if (method == "mean") {
        for (i in 1:ncol(data)) {
            data[is.na(data[, i]), i] <- mean(data[, i], 
                na.rm = TRUE)
        }
    } else if (method == "median") {
        for (i in 1:ncol(data)) {
            data[is.na(data[, i]), i] <- stats::median(data[, 
                i], na.rm = TRUE)
        }
    } else if (method == "1") {
        data[is.na(data)] <- 1
    } else if (method == "0") {
        data[is.na(data)] <- 0
    } else {
        data <- data
    }
    list$data <- data
    list$mz <- mz
    list$rt <- rt
    return(list)
}

#' Filter the data based on DoE, rsd, intensity
#' @param list list with data as peaks list, mz, rt and group information
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param tr logical. TRUE means dataset with technical replicates at the base level folder
#' @param rsdcft the rsd cutoff of all peaks in technical replicates
#' @param index the index of peaks considered, default NULL
#' @return list with group infomation, filtered peaks and index
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' getdoe(list)
#' }
#' @export
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}},\code{\link{getmzrt2}}, \code{\link{getimputation}}, \code{\link{getmr}}
getdoe <- function(list, inscf = 5, rsdcf = 100, rsdcft = 30, 
    imputation = "l", tr = F, index = NULL) {
    list <- getimputation(list, method = imputation)
    # use index
    if (!is.null(index)) {
        list$data <- list$data[index, ]
        list$mz <- list$mz[index]
        list$rt <- list$rt[index]
    }
    # remove the technical replicates and use biological
    # replicates instead
    if (tr) {
        data <- list$data
        lv <- list$group
        # group base on levels
        cols <- colnames(lv)
        mlv <- do.call(paste, c(lv[cols]))
        # get the rsd of the technical replicates
        meant <- stats::aggregate(t(data), list(mlv), mean)
        sdt <- stats::aggregate(t(data), list(mlv), sd)
        suppressWarnings(rsd <- sdt[, -1]/meant[, -1] * 
            100)
        data <- t(meant[, -1])
        colnames(data) <- unique(mlv)
        rsd <- t(rsd)
        # filter the data based on rsd of the technical
        # replicates
        indext <- as.vector(apply(rsd, 1, function(x) all(x < 
            rsdcft)))
        data <- data[indext, ]
        # data with mean of the technical replicates
        list$data <- data
        # get new group infomation
        ng <- NULL
        if (ncol(lv) > 1) {
            for (i in 1:(ncol(lv) - 1)) {
                lvi <- sapply(strsplit(unique(mlv), split = " ", 
                  fixed = TRUE), `[`, i)
                ng <- cbind(ng, lvi)
            }
            list$group <- data.frame(ng)
        } else {
            list$group <- data.frame(unique(mlv))
        }
        # save the index
        list$indext <- indext
    }
    
    # filter the data based on rsd/intensity
    data <- list$data
    lv <- list$group
    cols <- colnames(lv)
    # one peak for metabolomics is hard to happen
    if (nrow(lv) > 1) {
        if (ncol(lv) > 1) {
            mlv <- do.call(paste0, c(lv[cols], sep = ""))
        } else {
            mlv <- unlist(lv)
        }
        mean <- stats::aggregate(t(data), list(mlv), mean)
        sd <- stats::aggregate(t(data), list(mlv), sd)
        suppressWarnings(rsd <- sd[, -1]/mean[, -1] * 100)
        mean <- t(mean[, -1])
        sd <- t(sd[, -1])
        rsd <- t(rsd)
        colnames(rsd) <- colnames(sd) <- colnames(mean) <- unique(mlv)
        index <- as.vector(apply(rsd, 1, function(x) all(x < 
            rsdcf))) & as.vector(apply(mean, 1, function(x) any(x > 
            10^(inscf))))
        list$groupmean <- mean
        list$groupsd <- sd
        list$grouprsd <- rsd
        list$datafiltered <- data[index, ]
        list$mzfiltered <- list$mz[index]
        list$rtfiltered <- list$rt[index]
        list$groupmeanfiltered <- mean[index, ]
        list$groupsdfiltered <- sd[index, ]
        list$grouprsdfiltered <- rsd[index, ]
        list$index <- index
        return(list)
    } else {
        index <- data > 10^(inscf)
        list$datafiltered <- data[index]
        list$mzfiltered <- list$mz[index]
        list$rtfiltered <- list$rt[index]
        list$index <- index
        message("Only technical replicates were shown for ONE sample !!!")
        return(list)
    }
}

#' Get the features from t test, with p value, q value, rsd and power restriction
#' @param list list with data as peaks list, mz, rt and group information (two groups)
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param index the index of peaks considered, default NULL
#' @return dataframe with peaks fit the setting above
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' getfeaturest(list)
#' }
#' @export

getfeaturest <- function(list, power = 0.8, pt = 0.05, 
    qt = 0.05, n = 3, inscf = 5, rsdcf = 30, imputation = "l", 
    index = NULL) {
    list <- getdoe(list, inscf = inscf, rsdcf = rsdcf, 
        imputation = imputation, index = index)
    data <- list$datafiltered
    lv <- list$group$class
    sdn <- list$groupsdfiltered[, 1]
    mz <- list$mzfiltered
    rt <- list$rtfiltered
    
    ar <- genefilter::rowttests(data, fac = lv)
    dm <- ar$dm
    p <- ar$p.value
    q <- stats::p.adjust(p, method = "BH")
    m <- nrow(data)
    df <- cbind.data.frame(sdn, dm, p, q, mz, rt, data)
    df <- df[order(df$p), ]
    df$alpha <- c(1:m) * pt/m
    rp <- vector()
    for (i in c(1:nrow(df))) {
        r <- stats::power.t.test(delta = df$dm[i], sd = df$sd[i], 
            sig.level = df$alpha[i], n = n)
        rp[i] <- r$power
    }
    df <- cbind(power = rp, df)
    df <- df[df$power > power, ]
    return(df)
}

#' Get the features from anova, with p value, q value, rsd and power restriction
#' @param list list with data as peaks list, mz, rt and group information (more than two groups)
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param ng group numbers
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param index the index of peaks considered, default NULL
#' @return dataframe with peaks fit the setting above
#' @export

getfeaturesanova <- function(list, power = 0.8, pt = 0.05, 
    qt = 0.05, n = 3, ng = 3, rsdcf = 100, inscf = 5, imputation = "l", 
    index = NULL) {
    list <- getdoe(list, inscf = inscf, rsdcf = rsdcf, 
        imputation = imputation, index = index)
    data <- list$datafiltered
    lv <- list$group$class
    sdn <- list$groupsdfiltered[, 1]
    mz <- list$mzfiltered
    rt <- list$rtfiltered
    rsd <- list$grouprsdfiltered
    
    sdn <- genefilter::rowSds(data[, 1:n])
    sdg <- genefilter::rowSds(list$groupmeanfiltered)
    
    ar <- genefilter::rowFtests(data, lv)
    p <- ar$p.value
    q <- stats::p.adjust(p, method = "BH")
    m <- nrow(data)
    df <- cbind.data.frame(sdn, sdg, rsd, p, q, mz, rt, 
        data)
    df <- df[order(df$p), ]
    df$alpha <- c(1:m) * pt/m
    rp <- vector()
    for (i in c(1:nrow(df))) {
        r <- stats::power.anova.test(groups = ng, between.var = df$sdg[i], 
            within.var = df$sdn[i], sig.level = df$alpha[i], 
            n = n)
        rp[i] <- r$power
    }
    df <- cbind(power = rp, df)
    df <- df[df$power > power, ]
    return(df)
}

#' plot the scatter plot for peaks list with threshold
#' @param list list with data as peaks list, mz, rt and group information
#' @param rt vector range of the retention time
#' @param ms vector vector range of the m/z
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param index the index of peaks considered, default NULL
#' @param ... parameters for `plot` function
#' @return data fit the cutoff
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' plotmr(list)
#' }
#' @export
plotmr <- function(list, rt = NULL, ms = NULL, inscf = 5, 
    rsdcf = 30, imputation = "l", index = NULL, ...) {
    graphics::par(mar = c(5, 4.2, 6.1, 2.1), xpd = TRUE)
    list <- getdoe(list, rsdcf = rsdcf, inscf = inscf, 
        imputation = imputation, index = index)
    data <- list$groupmeanfiltered
    dataname <- colnames(data)
    mz <- list$mzfiltered
    RT <- list$rtfiltered
    suppressWarnings(if (nrow(data) > 0) {
        n <- dim(data)[2]
        col <- grDevices::rainbow(n, alpha = 0.318)
        
        graphics::plot(mz ~ RT, xlab = "Retention Time(s)", 
            ylab = "m/z", type = "n", pch = 19, ylim = ms, 
            xlim = rt, ...)
        
        for (i in 1:n) {
            cex = as.numeric(cut((log10(data[, i] + 1) - 
                inscf), breaks = c(0, 1, 2, 3, 4, Inf)/2))/2
            cexlab = c(paste0(inscf, "-", inscf + 0.5), 
                paste0(inscf + 0.5, "-", inscf + 1), paste0(inscf + 
                  1, "-", inscf + 1.5), paste0(inscf + 
                  1.5, "-", inscf + 2), paste0(">", inscf + 
                  2))
            if (!is.null(ms) & !is.null(rt)) {
                graphics::points(x = RT[RT > rt[1] & RT < 
                  rt[2] & mz > ms[1] & mz < ms[2]], y = mz[RT > 
                  rt[1] & RT < rt[2] & mz > ms[1] & mz < 
                  ms[2]], cex = cex, col = col[i], pch = 19)
            } else {
                graphics::points(x = RT, y = mz, cex = cex, 
                  col = col[i], pch = 19)
            }
            
        }
        graphics::legend("topright", legend = dataname, 
            col = col, pch = 19, horiz = T, bty = "n", 
            inset = c(0, -0.25))
        graphics::legend("topleft", legend = cexlab, title = "Intensity in Log scale", 
            pt.cex = c(1, 2, 3, 4, 5)/4, pch = 19, bty = "n", 
            horiz = T, cex = 0.7, col = grDevices::rgb(0, 
                0, 0, 0.318), inset = c(0, -0.25))
    } else {
        graphics::plot(1, xlab = "Retention Time(s)", ylab = "m/z", 
            main = "No peaks found", ylim = ms, xlim = rt, 
            type = "n", pch = 19, ...)
    })
}

#' plot the diff scatter plot for one xcmsset objects with threshold between two groups
#' @param list list with data as peaks list, mz, rt and group information
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param index the index of peaks considered, default NULL
#' @param ... parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' plotmrc(list)
#' }
#' @export
plotmrc <- function(list, ms = c(100, 800), inscf = 5, 
    rsdcf = 30, imputation = "l", index = NULL, ...) {
    list <- getdoe(list, rsdcf = rsdcf, inscf = inscf, 
        imputation = imputation, index = NULL)
    data <- list$groupmeanfiltered
    dataname <- colnames(data)
    mz <- list$mzfiltered
    rt <- list$rtfiltered
    suppressWarnings(if (!is.na(data)) {
        
        diff1 <- data[, 1] - data[, 2]
        diff2 <- data[, 2] - data[, 1]
        diff1[diff1 < 0] <- 0
        diff2[diff2 < 0] <- 0
        name1 <- paste0(dataname[1], "-", dataname[2])
        name2 <- paste0(dataname[2], "-", dataname[1])
        
        cex1 = as.numeric(cut((log10(diff1 + 1) - inscf), 
            breaks = c(0, 1, 2, 3, 4, Inf)/2))/2
        cex2 = as.numeric(cut((log10(diff2 + 1) - inscf), 
            breaks = c(0, 1, 2, 3, 4, Inf)/2))/2
        cexlab = c(paste0(inscf, "-", inscf + 0.5), paste0(inscf + 
            0.5, "-", inscf + 1), paste0(inscf + 1, "-", 
            inscf + 1.5), paste0(inscf + 1.5, "-", inscf + 
            2), paste0(">", inscf + 2))
        
        graphics::plot(mz ~ rt, xlab = "Retention Time(s)", 
            ylab = "m/z", ylim = ms, cex = cex1, col = grDevices::rgb(0, 
                0, 1, 0.618), pch = 19, ...)
        
        graphics::points(mz ~ rt, cex = cex2, col = grDevices::rgb(1, 
            0, 0, 0.618), pch = 19)
        
        graphics::legend("topleft", legend = cexlab, title = "Intensity in Log scale", 
            pt.cex = c(1, 2, 3, 4, 5)/2, pch = 19, col = grDevices::rgb(0, 
                0, 0, 0.618), bty = "n", horiz = T, inset = c(0, 
                -0.25))
        graphics::legend("topright", legend = c(name1, 
            name2), pch = 19, col = c(grDevices::rgb(0, 
            0, 1, 0.618), grDevices::rgb(1, 0, 0, 0.618)), 
            bty = "n", horiz = T, inset = c(0, -0.25))
    } else {
        graphics::plot(1, xlab = "Retention Time(s)", ylab = "m/z", 
            main = "No peaks found", ylim = ms, type = "n", 
            pch = 19, ...)
    })
    
}

#' plot the rsd influnces of data in different groups
#' @param list list with data as peaks list, mz, rt and group information
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param index the index of peaks considered, default NULL
#' @param ... other parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' plotrsd(list)
#' }
#' @export
plotrsd <- function(list, ms = c(100, 800), inscf = 5, 
    rsdcf = 100, imputation = "l", index = NULL, ...) {
    cexlab = c("<20%", "20-40%", "40-60%", "60-80%", ">80%")
    list <- getdoe(list, rsdcf = rsdcf, inscf = inscf, 
        imputation = imputation, index = NULL)
    data <- list$groupmeanfiltered
    dataname <- colnames(data)
    mz <- list$mzfiltered
    rt <- list$rtfiltered
    rsd <- list$grouprsdfiltered
    
    n <- dim(rsd)[2]
    col <- grDevices::rainbow(n, alpha = 0.318)
    
    graphics::plot(mz ~ rt, xlab = "Retention Time(s)", 
        ylab = "m/z", main = "RSD(%) distribution", type = "n", 
        pch = 19, ...)
    
    for (i in 1:n) {
        cex = as.numeric(cut(rsd[, i], breaks = c(0, 20, 
            40, 60, 80, Inf)))/2
        graphics::points(x = rt, y = mz, ylim = ms, cex = cex, 
            col = col[i], pch = 19)
    }
    graphics::legend("topright", legend = dataname, col = col, 
        horiz = T, pch = 19, bty = "n", inset = c(0, -0.25))
    graphics::legend("topleft", legend = cexlab, pt.cex = c(1, 
        2, 3, 4, 5)/2, pch = 19, bty = "n", horiz = T, 
        cex = 0.8, col = grDevices::rgb(0, 0, 0, 0.318), 
        inset = c(0, -0.25))
}

#' plot scatter plot for rt-mz profile and output gif file for mutiple groups
#' @param list list with data as peaks list, mz, rt and group information
#' @param file file name for gif file, default NULL
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param index the index of peaks considered, default NULL
#' @param ... parameters for `plot` function
#' @return gif file
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' gifmr(list)
#' }
#' @export
gifmr <- function(list, ms = c(100, 500), rsdcf = 30, inscf = 5, 
    imputation = "i", index = NULL, file = "test") {
    list <- getdoe(list, rsdcf = rsdcf, inscf = inscf, 
        imputation = imputation, index = NULL)
    data <- list$groupmeanfiltered
    mz <- list$mzfiltered
    rt <- list$rtfiltered
    filename = paste0(file, ".gif")
    mean <- apply(data, 1, mean)
    
    graphics::plot(mz ~ rt, xlab = "Retention Time(s)", 
        ylab = "m/z", pch = 19, cex = log10(mean + 1) - 
            4, col = grDevices::rgb(0, 0, 1, 0.2), main = "All peaks")
    
    col <- grDevices::rainbow(dim(data)[2], alpha = 0.318)
    animation::saveGIF({
        for (i in 1:dim(data)[2]) {
            name <- colnames(data)[i]
            value <- data[, i]
            graphics::plot(mz ~ rt, xlab = "Retention Time(s)", 
                ylab = "m/z", pch = 19, cex = log10(value + 
                  1) - 4, col = col[i], ylim = ms, main = name)
        }
    }, movie.name = filename, ani.width = 800, ani.height = 500)
}

#' plot the PCA of list
#' @param data mzrt profile with row peaks and column samples
#' @param lv group information
#' @param index index for selected peaks
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param ... other parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' data <- list$data
#' lv <- as.character(list$group$class)
#' plotpca(data, lv)
#' }
#' @export

plotpca <- function(data, lv = NULL, index = NULL, center = T, 
    scale = T, ...) {
    
    if (!is.null(index)) {
        data <- data[index, ]
    }
    
    if (is.null(lv)) {
        pch = colnames(data)
    } else {
        pch = lv
    }
    pcao <- stats::prcomp(t(data), center = center, scale = scale)
    pcaoVars = signif(((pcao$sdev)^2)/(sum((pcao$sdev)^2)), 
        3) * 100
    graphics::plot(pcao$x[, 1], pcao$x[, 2], xlab = paste("PC1:", 
        pcaoVars[1], "% of Variance Explained"), ylab = paste("PC2:", 
        pcaoVars[2], "% of Variance Explained"), cex = 2, 
        pch = pch, ...)
}

#' Plot the heatmap of mzrt profiles
#' @param data mzrt profile with row peaks and column samples
#' @param lv group information
#' @param index index for selected peaks
#' @return NULL
#' @export
plothm <- function(data, lv, index = NULL) {
    icolors <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 
        "RdYlBu"))))(100)
    zlim <- range(data)
    if (!is.null(index)) {
        data <- data[index, order(lv)]
    } else {
        data <- data[, order(lv)]
    }
    plotchange <- function(zlim) {
        breaks <- seq(zlim[1], zlim[2], round((zlim[2] - 
            zlim[1])/10))
        poly <- vector(mode = "list", length(icolors))
        graphics::plot(1, 1, t = "n", xlim = c(0, 1), ylim = zlim, 
            xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
            ylab = "", xlab = "", frame.plot = F)
        graphics::axis(4, at = breaks, labels = round(breaks), 
            las = 1, pos = 0.4, cex.axis = 0.8)
        p <- graphics::par("usr")
        graphics::text(p[2] + 2, mean(p[3:4]), labels = "intensity", 
            xpd = NA, srt = -90)
        bks <- seq(zlim[1], zlim[2], length.out = (length(icolors) + 
            1))
        for (i in seq(poly)) {
            graphics::polygon(c(0.1, 0.1, 0.3, 0.3), c(bks[i], 
                bks[i + 1], bks[i + 1], bks[i]), col = icolors[i], 
                border = NA)
        }
    }
    pos <- cumsum(as.numeric(table(lv)/sum(table(lv)))) - 
        as.numeric(table(lv)/sum(table(lv)))/2
    posv <- cumsum(as.numeric(table(lv)/sum(table(lv))))[1:(nlevels(lv) - 
        1)]
    
    
    graphics::layout(matrix(rep(c(1, 1, 1, 2), 10), 10, 
        4, byrow = TRUE))
    graphics::par(mar = c(3, 6, 2, 1))
    graphics::image(t(data), col = icolors, xlab = "samples", 
        main = "peaks intensities on log scale", xaxt = "n", 
        yaxt = "n", zlim = zlim)
    graphics::axis(1, at = pos, labels = levels(lv), cex.axis = 0.8)
    graphics::axis(2, at = seq(0, 1, 1/(nrow(data) - 1)), 
        labels = rownames(data), cex.axis = 1, las = 2)
    graphics::abline(v = posv)
    graphics::par(mar = c(3, 1, 2, 6))
    plotchange(zlim)
}
