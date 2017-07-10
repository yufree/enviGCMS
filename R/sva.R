#' Surrogate variable analysis(SVA) to correct the unknown batch effects
#' @param xset xcmsset object
#' @param lv group information
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @details this is used for reviesed version of SVA to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, signal part, random errors part, batch part, p-values, q-values, mass, rt, Posterior Probabilities of Surrogate variables and Posterior Probabilities of Mod. If no surrogate variable found, corresponding part would miss.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' }
#' @seealso \code{\link{svapca}}, \code{\link{svaplot}}, \code{\link{svabatch}}
#' @export
svacor <- function(xset, lv = NULL, method = "medret",
    intensity = "into") {
    data <- xcms::groupval(xset, method, intensity)
    if (intensity == "intb") {
        data[is.na(data)] = 0
    }
    if (is.null(lv)) {
        lv <- xset@phenoData[, 1]
    }
    mz <- xset@groups[, 1]
    rt <- xset@groups[, 4]
    mod <- stats::model.matrix(~lv)
    mod0 <- as.matrix(c(rep(1, ncol(data))))
    svafit <- sva::sva(data, mod)
    if (svafit$n.sv == 0) {
        svaX <- stats::model.matrix(~lv)
        lmfit <- limma::lmFit(data, svaX)
        signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[,
            1:nlevels(lv)])
        error <- data - signal
        rownames(signal) <- rownames(error) <- rownames(data)
        colnames(signal) <- colnames(error) <- colnames(data)
        pValues = sva::f.pvalue(data, mod, mod0)
        qValues = qvalue::qvalue(pValues)
        qValues = qValues$qvalues

        li <- list(data, signal, error, pValues, qValues,
            mz, rt)
        names(li) <- c("data", "signal", "error", "p-values",
            "q-values", "mz", "rt")

    } else {
        message("Data is correcting ...")
        svaX <- stats::model.matrix(~lv + svafit$sv)
        lmfit <- limma::lmFit(data, svaX)
        batch <- lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) +
            svafit$n.sv)] %*% t(svaX[, (nlevels(lv) +
            1):(nlevels(lv) + svafit$n.sv)])
        signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[,
            1:nlevels(lv)])
        error <- data - signal - batch
        datacor <- signal + error
        svaX2 <- stats::model.matrix(~lv)
        lmfit2 <- limma::lmFit(data, svaX2)
        signal2 <- lmfit2$coef[, 1:nlevels(lv)] %*%
            t(svaX2[, 1:nlevels(lv)])
        error2 <- data - signal2
        rownames(signal2) <- rownames(error2) <- rownames(datacor) <- rownames(signal) <- rownames(batch) <- rownames(error) <- rownames(data)
        colnames(signal2) <- colnames(error2) <- colnames(datacor) <- colnames(signal) <- colnames(batch) <- colnames(error) <- colnames(data)

        modSv = cbind(mod, svafit$sv)
        mod0Sv = cbind(mod0, svafit$sv)
        pValuesSv = sva::f.pvalue(data, modSv, mod0Sv)
        qValuesSv = qvalue::qvalue(pValuesSv)
        qValuesSv = qValuesSv$qvalues

        pValues = sva::f.pvalue(data, mod, mod0)
        qValues = qvalue::qvalue(pValues)
        qValues = qValues$qvalues
        li <- list(data, datacor, signal, batch, error,
            signal2, error2, pValues, qValues, pValuesSv,
            qValuesSv, svafit$pprob.gam, svafit$pprob.b,
            mz, rt)
        names(li) <- c("data", "dataCorrected", "signal",
            "batch", "error", "signal2", "error2",
            "p-values", "q-values", "p-valuesCorrected",
            "q-valuesCorrected", "PosteriorProbabilitiesSurrogate",
            "PosteriorProbabilitiesMod", "mz", "rt")
        message("Done!")
    }
    return(li)
}

#' Principal component analysis(PCA) for SVA corrected data and raw data
#' @param list results from svacor function
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param lv group information
#' @return plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' svapca(df)
#' }
#' @seealso \code{\link{svacor}}, \code{\link{svaplot}}, \code{\link{svabatch}}
#' @export
svapca <- function(list, center = T, scale = T, lv = NULL) {
    data <- list$data
    Signal <- list$signal
    Batch <- list$batch
    error <- list$error
    datacor <- list$dataCorrected
    if (is.null(lv)) {
        pch = colnames(data)
    } else {
        pch = lv
    }

    graphics::par(mfrow = c(2, 5), mar = c(4, 4, 2.6,
        1))

    pcao <- stats::prcomp(t(data), center = center,
        scale = scale)
    pcaoVars = signif(((pcao$sdev)^2)/(sum((pcao$sdev)^2)),
        3) * 100
    graphics::plot(pcao, type = "l", main = "PCA")

    pca <- stats::prcomp(t(Signal), center = TRUE,
        scale = TRUE)
    pcaVars = signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),
        3) * 100
    graphics::plot(pca, type = "l", main = "PCA-signal")

    pcab <- stats::prcomp(t(Batch), center = center,
        scale = scale)
    pcabVars = signif(((pcab$sdev)^2)/(sum((pcab$sdev)^2)),
        3) * 100
    graphics::plot(pcab, type = "l", main = "PCA-batch")

    pcae <- stats::prcomp(t(error), center = center,
        scale = scale)
    pcaeVars = signif(((pcae$sdev)^2)/(sum((pcae$sdev)^2)),
        3) * 100
    graphics::plot(pcae, type = "l", main = "PCA-error")

    pcac <- stats::prcomp(t(datacor), center = center,
        scale = scale)
    pcacVars = signif(((pcac$sdev)^2)/(sum((pcac$sdev)^2)),
        3) * 100
    graphics::plot(pcac, type = "l", main = "PCA-corrected")

    graphics::plot(pcao$x[, 1], pcao$x[, 2], xlab = paste("PC1:",
        pcaoVars[1], "% of Variance Explained"), ylab = paste("PC2:",
        pcaoVars[2], "% of Variance Explained"), pch = pch,
        cex = 2, main = "PCA")

    graphics::plot(pca$x[, 1], pca$x[, 2], xlab = paste("PC1:",
        pcaVars[1], "% of Variance Explained"), ylab = paste("PC2:",
        pcaVars[2], "% of Variance Explained"), pch = pch,
        cex = 2, main = "PCA-signal")

    graphics::plot(pcab$x[, 1], pcab$x[, 2], xlab = paste("PC1:",
        pcabVars[1], "% of Variance Explained"), ylab = paste("PC2:",
        pcabVars[2], "% of Variance Explained"), pch = pch,
        cex = 2, main = "PCA-batch")

    graphics::plot(pcae$x[, 1], pcae$x[, 2], xlab = paste("PC1:",
        pcaeVars[1], "% of Variance Explained"), ylab = paste("PC2:",
        pcaeVars[2], "% of Variance Explained"), pch = pch,
        cex = 2, main = "PCA-error")

    graphics::plot(pcac$x[, 1], pcac$x[, 2], xlab = paste("PC1:",
        pcacVars[1], "% of Variance Explained"), ylab = paste("PC2:",
        pcacVars[2], "% of Variance Explained"), pch = pch,
        cex = 2, main = "PCA-corrected")
}
#' Filter the data with p value and q value
#' @param list results from svacor function
#' @param pqvalues method for ANOVA or SVA
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @return data, corrected data, mz and retention for fileted data
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' svadata(df)
#' }
#' @export
svadata <- function(list, pqvalues = "sv", pt = 0.05,
    qt = 0.05) {
    data <- list$data
    signal2 <- list$signal2
    datacor <- list$dataCorrected
    pValues <- list$"p-values"
    qValues <- list$"q-values"
    pValuesSv <- list$"p-valuesCorrected"
    qValuesSv <- list$"q-valuesCorrected"
    mz <- list$mz
    rt <- list$rt
    if (is.null(signal2)) {
        if (pqvalues == "anova" & sum(pValues < pt &
            qValues < qt) != 0) {
            message("No SV while p-values and q-values have results")
            data <- data[pValues < pt & qValues < qt,
                ]
            mz <- mz[pValues < pt & qValues < qt]
            rt <- rt[pValues < pt & qValues < qt]
            li <- list(data, pValues < pt & qValues <
                qt, mz, rt)
            names(li) <- c("data", "pqvalues", "mz",
                "rt")
            return(li)
        } else {
            message("No SV while p-values and q-values have no results")
        }
    } else {
        if (pqvalues == "anova" & sum(pValues < pt &
            qValues < qt) != 0) {
            message("Have SVs while p-values and q-values have results")
            data <- data[pValues < pt & qValues < qt,
                ]
            mz <- mz[pValues < pt & qValues < qt]
            rt <- rt[pValues < pt & qValues < qt]
            li <- list(data, pValues < pt & qValues <
                qt, mz, rt)
            names(li) <- c("data", "pqvalues", mz,
                rt)
            return(li)
        } else if (pqvalues == "anova") {
            message("Have SVs while p-values and q-values have no results")
        } else if (pqvalues == "sv" & sum(pValuesSv <
            pt & qValuesSv < qt) != 0) {
            message("SVs corrected while p-values and q-values have results")
            data <- data[pValuesSv < pt & qValuesSv <
                qt, ]
            datacor <- datacor[pValuesSv < pt & qValuesSv <
                qt, ]
            mz <- mz[pValuesSv < pt & qValuesSv < qt]
            rt <- rt[pValuesSv < pt & qValuesSv < qt]
            li <- list(datacor, data, pValuesSv < pt &
                qValuesSv < qt, mz, rt)
            names(li) <- c("dataCorrected", "data",
                "pqvalues", "mz", "rt")
            return(li)
        } else {
            message("SVs corrected while p-values and q-values have no results")
        }
    }
}
#' Filter the data with p value and q value and show them
#' @param list results from svacor function
#' @param pqvalues method for ANOVA or SVA
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @param lv group information
#' @param index index for selected peaks
#' @return heatmap for the data
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' svaplot(df)
#' }
#' @seealso \code{\link{svacor}}, \code{\link{svapca}}, \code{\link{svabatch}}
#' @export
svaplot <- function(list, pqvalues = "sv", pt = 0.05,
    qt = 0.05, lv = NULL, index = NULL) {
    data <- list$data
    signal <- list$signal
    signal2 <- list$signal2
    batch <- list$batch
    error <- list$error
    error2 <- list$error2
    datacor <- list$dataCorrected
    pValues <- list$"p-values"
    qValues <- list$"q-values"
    pValuesSv <- list$"p-valuesCorrected"
    qValuesSv <- list$"q-valuesCorrected"
    if (!is.null(index)) {
        data <- data[index, ]
        signal <- signal[index, ]
        signal2 <- signal2[index, ]
        batch <- batch[index, ]
        error <- error[index, ]
        error2 <- error2[index, ]
        datacor <- datacor[index, ]
        pValues <- pValues[index]
        qValues <- qValues[index]
        pValuesSv <- pValuesSv[index]
        qValuesSv <- qValuesSv[index]
    }

    if (is.null(lv)) {
        lv <- as.factor(colnames(data))
    }
    pos <- cumsum(as.numeric(table(lv)/sum(table(lv)))) -
        as.numeric(table(lv)/sum(table(lv)))/2
    posv <- cumsum(as.numeric(table(lv)/sum(table(lv))))[1:(nlevels(lv) -
        1)]

    icolors <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,
        "RdYlBu"))))(100)

    plotchange <- function(zlim) {
        breaks <- seq(zlim[1], zlim[2], round((zlim[2] -
            zlim[1])/10))
        poly <- vector(mode = "list", length(icolors))
        graphics::plot(1, 1, t = "n", xlim = c(0, 1),
            ylim = zlim, xaxt = "n", yaxt = "n", xaxs = "i",
            yaxs = "i", ylab = "", xlab = "", frame.plot = F)
        graphics::axis(4, at = breaks, labels = round(breaks),
            las = 1, pos = 0.4, cex.axis = 0.8)
        p <- graphics::par("usr")
        graphics::text(p[2] + 2, mean(p[3:4]), labels = "intensity",
            xpd = NA, srt = -90)
        bks <- seq(zlim[1], zlim[2], length.out = (length(icolors) +
            1))
        for (i in seq(poly)) {
            graphics::polygon(c(0.1, 0.1, 0.3, 0.3),
                c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                col = icolors[i], border = NA)
        }
    }
    plotimage1 <- function(data, signal, error, zlim) {
        graphics::image(t(data), col = icolors, xlab = "samples",
            main = "peaks", xaxt = "n", yaxt = "n",
            zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv),
            cex.axis = 0.8)
        graphics::axis(2, at = seq(0, 1, 1/(nrow(data) -
            1)), labels = rownames(data), cex.axis = 1,
            las = 2)
        graphics::abline(v = posv)

        graphics::image(t(signal), col = icolors, xlab = "samples",
            main = "peaks-signal", xaxt = "n", yaxt = "n",
            zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv),
            cex.axis = 0.8)
        graphics::axis(2, at = seq(0, 1, 1/(nrow(signal) -
            1)), labels = rownames(signal), cex.axis = 1,
            las = 2)
        graphics::abline(v = posv)

        graphics::image(t(error), col = icolors, xlab = "samples",
            main = "peaks-error", xaxt = "n", yaxt = "n",
            zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv),
            cex.axis = 0.8)
        graphics::axis(2, at = seq(0, 1, 1/(nrow(error) -
            1)), labels = rownames(error), cex.axis = 1,
            las = 2)
        graphics::abline(v = posv)
    }
    plotimage2 <- function(data, signal, batch, error,
        datacor, zlim) {

        graphics::image(t(data), col = icolors, xlab = "samples",
            main = "peaks", xaxt = "n", yaxt = "n",
            zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv))
        graphics::axis(2, at = seq(0, 1, 1/(nrow(data) -
            1)), labels = rownames(data), las = 1)
        graphics::abline(v = posv)

        graphics::image(t(signal), col = icolors, xlab = "samples",
            main = "peaks-signal", xaxt = "n", yaxt = "n",
            zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv),
            cex.axis = 1)
        graphics::axis(2, at = seq(0, 1, 1/(nrow(signal) -
            1)), labels = rownames(signal), las = 1)
        graphics::abline(v = posv)

        graphics::image(t(batch), col = icolors, xlab = "samples",
            main = "peaks-batch", xaxt = "n", yaxt = "n",
            zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv))
        graphics::axis(2, at = seq(0, 1, 1/(nrow(batch) -
            1)), labels = rownames(batch), las = 1)
        graphics::abline(v = posv)

        graphics::image(t(error), col = icolors, xlab = "samples",
            main = "peaks-error", xaxt = "n", yaxt = "n",
            zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv))
        graphics::axis(2, at = seq(0, 1, 1/(nrow(error) -
            1)), labels = rownames(error), las = 1)
        graphics::abline(v = posv)

        graphics::image(t(datacor), col = icolors,
            xlab = "samples", main = "peaks-corrected",
            xaxt = "n", yaxt = "n", zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv))
        graphics::axis(2, at = seq(0, 1, 1/(nrow(datacor) -
            1)), labels = rownames(datacor), las = 1)
        graphics::abline(v = posv)

    }

    if (is.null(signal2)) {
        if (pqvalues == "anova" & sum(pValues < pt &
            qValues < qt) != 0) {
            message("No SV while p-values and q-values have results")
            graphics::layout(matrix(rep(c(1, 1, 2,
                2, 3, 3, 4, 4, 5), 9), 9, 9, byrow = TRUE))
            data <- data[pValues < pt & qValues < qt,
                ]
            signal <- signal[pValues < pt & qValues <
                qt, ]
            error <- error[pValues < pt & qValues <
                qt, ]

            zlim <- range(c(data, signal, error))
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage1(data, signal, error, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)

            li <- list(data, pValues < pt & qValues <
                qt)
            names(li) <- c("data", "pqvalues")
            return(li)
        } else {
            message("No SV while p-values and q-values have no results")
            graphics::layout(matrix(rep(c(1, 1, 1,
                2, 2, 3, 3, 4), 8), 8, 8, byrow = TRUE))
            zlim <- range(c(data, signal, error))
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage1(data, signal, error, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
        }
    } else {
        if (pqvalues == "anova" & sum(pValues < pt &
            qValues < qt) != 0) {
            message("Have SVs while p-values and q-values have results")
            graphics::layout(matrix(rep(c(1, 1, 2,
                2, 3, 3, 4), 7), 7, 7, byrow = TRUE))
            data <- data[pValues < pt & qValues < qt,
                ]
            signal <- signal2[pValues < pt & qValues <
                qt, ]
            error <- error2[pValues < pt & qValues <
                qt, ]
            zlim <- range(c(data, signal, error))

            graphics::par(mar = c(3, 6, 2, 1))
            plotimage1(data, signal, error, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
            li <- list(data, pValues < pt & qValues <
                qt)
            names(li) <- c("data", "pqvalues")
            return(li)
        } else if (pqvalues == "anova") {
            message("Have SVs while p-values and q-values have no results")
            graphics::layout(matrix(rep(c(1, 1, 1,
                2, 2, 3, 3, 4), 8), 8, 8, byrow = TRUE))
            zlim <- range(c(data, signal2, error2))

            graphics::par(mar = c(3, 6, 2, 1))
            plotimage1(data, signal2, error2, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
        } else if (pqvalues == "sv" & sum(pValuesSv <
            pt & qValuesSv < qt) != 0) {
            message("SVs corrected while p-values and q-values have results")
            graphics::layout(matrix(rep(c(1, 1, 2,
                2, 3, 3, 4, 4, 5, 5, 6), 11), 11, 11,
                byrow = TRUE))
            data <- data[pValuesSv < pt & qValuesSv <
                qt, ]
            signal <- signal[pValuesSv < pt & qValuesSv <
                qt, ]
            batch <- batch[pValuesSv < pt & qValuesSv <
                qt, ]
            error <- error[pValuesSv < pt & qValuesSv <
                qt, ]
            datacor <- datacor[pValuesSv < pt & qValuesSv <
                qt, ]
            zlim <- range(c(data, signal, batch, error,
                datacor))
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage2(data, signal, batch, error,
                datacor, zlim)
            graphics::par(mar = c(3, 1, 2, 5))
            plotchange(zlim)
            li <- list(datacor, data, pValuesSv < pt &
                qValuesSv < qt)
            names(li) <- c("dataCorrected", "data",
                "pqvalues")
            return(li)
        } else {
            message("SVs corrected while p-values and q-values have no results")
            graphics::layout(matrix(rep(c(1, 1, 1,
                2, 2, 3, 3, 4, 4, 5, 5, 5, 6), 13),
                13, 13, byrow = TRUE))
            zlim <- range(c(signal, data, batch, error,
                datacor))
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage2(data, signal, batch, error,
                datacor, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
        }
    }
}

#' Get the corrected data after SVA for metabolanalyst
#' @param xset xcmsset object
#' @param lv group information
#' @return csv files for both raw and corrected data for metabolanalyst if SVA could be applied
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' svaupload(xset3)
#' }
#' @export
svaupload <- function(xset, lv = NULL) {
    raw <- svacor(xset, lv = lv)
    if (is.null(raw$dataCorrected)) {
        data <- raw$data
        data <- rbind(group = as.character(xset@phenoData[,
            1]), data)
        utils::write.csv(data, file = "Peaklist.csv")
        return(data)
    } else {
        datacor <- raw$dataCorrected
        data <- raw$data
        data <- rbind(group = as.character(xset@phenoData[,
            1]), data)
        utils::write.csv(datacor, file = "Peaklistcor.csv")
        utils::write.csv(data, file = "Peaklist.csv")
        return(list(data, datacor))
    }
}

#' Plot the influnces of DoE and Batch effects on each peaks
#' @param df data output from `svacor` function
#' @param dfsv data output from `svaplot` function for corrected data
#' @param dfanova data output from `svaplot` function for raw data
#' @return influnces plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' dfsv <- svaplot(xset3)
#' dfanova <- svaplot(xset3, pqvalues = "anova")
#' svabatch(df,dfsv,dfanova)
#' }
#' @seealso \code{\link{svacor}}, \code{\link{svaplot}}, \code{\link{svapca}}
#' @export
svabatch <- function(df, dfsv, dfanova) {
    graphics::par(mfrow = c(1, 2))
    p0 <- -log10(df$`p-valuesCorrected`)
    p0[is.infinite(p0)] <- max(p0[!is.infinite(p0)]) +
        1

    p1 <- -log10(df$`p-values`)
    p1[is.infinite(p1)] <- max(p1[!is.infinite(p1)]) +
        1

    big <- p0 - p1

    graphics::plot(df$PosteriorProbabilitiesSurrogate[big <
        0] ~ df$PosteriorProbabilitiesMod[big < 0],
        xlab = "Influences  from primary variable",
        ylab = "Influences  from surrogate variable",
        main = "Peaks with larger p-values", pch = 19,
        cex = 1.2, ylim = c(0, 1), xlim = c(0, 1.02),
        col = grDevices::rgb(0, 0, 0, 0.1))
    graphics::points(df$PosteriorProbabilitiesSurrogate[big <
        0 & dfsv$pqvalues & dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[big <
        0 & dfsv$pqvalues & dfanova$pqvalue], pch = 19,
        col = "red")
    graphics::points(df$PosteriorProbabilitiesSurrogate[big <
        0 & !dfsv$pqvalues & dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[big <
        0 & !dfsv$pqvalues & dfanova$pqvalue], pch = 19,
        col = "green")
    graphics::plot(df$PosteriorProbabilitiesSurrogate[big >
        0] ~ df$PosteriorProbabilitiesMod[big > 0],
        xlab = "Influences  from primary variable",
        ylab = "Influences  from surrogate variable",
        main = "Peaks with smaller p-values", pch = 19,
        cex = 1.2, ylim = c(0, 1), xlim = c(0, 1.02),
        col = grDevices::rgb(0, 0, 0, 0.1))
    graphics::points(df$PosteriorProbabilitiesSurrogate[big >
        0 & dfsv$pqvalues & dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[big >
        0 & dfsv$pqvalues & dfanova$pqvalue], pch = 19,
        col = "red")
    graphics::points(df$PosteriorProbabilitiesSurrogate[big >
        0 & dfsv$pqvalues & !dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[big >
        0 & dfsv$pqvalues & !dfanova$pqvalue], pch = 19,
        col = "blue")

}
