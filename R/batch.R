#' Use Surrogate Variable Analysis(SVA) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details this is used for SVA to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, signal part, random errors part, batch part, p-values, q-values, mass, rt, Posterior Probabilities of Surrogate variables and Posterior Probabilities of Mod. If no surrogate variable found, corresponding part would miss.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- svacor(list$data,list$group$class)
#' }
#' @seealso \code{\link{svapca}}, \code{\link{svaplot}}, \code{\link{svabatch}}
#' @export
svacor <- function(data, lv) {
    mod <- stats::model.matrix(~lv)
    mod0 <- as.matrix(c(rep(1, ncol(data))))
    svafit <- sva::sva(data, mod)
    if (svafit$n.sv == 0) {
        message("No surrogate variable found")
        svaX <- stats::model.matrix(~lv)
        lmfit <- limma::lmFit(data, svaX)
        signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 
            1:nlevels(lv)])
        error <- data - signal
        rownames(signal) <- rownames(error) <- rownames(data)
        colnames(signal) <- colnames(error) <- colnames(data)
        # find the peaks with significant differences by F test
        # with BH correction for fdr control
        pValues = sva::f.pvalue(data, mod, mod0)
        qValues = stats::p.adjust(pValues, method = "BH")
        # get the results as list
        li <- list(data, signal, error, pValues, qValues)
        names(li) <- c("data", "signal", "error", "p-values", 
            "q-values")
        
    } else {
        message("Data is correcting ...")
        svaX <- stats::model.matrix(~lv + svafit$sv)
        lmfit <- limma::lmFit(data, svaX)
        # data decomposition with sv
        batch <- lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) + 
            svafit$n.sv)] %*% t(svaX[, (nlevels(lv) + 1):(nlevels(lv) + 
            svafit$n.sv)])
        signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 
            1:nlevels(lv)])
        error <- data - signal - batch
        datacor <- signal + error
        svaX2 <- stats::model.matrix(~lv)
        lmfit2 <- limma::lmFit(data, svaX2)
        # data decomposition without sv
        signal2 <- lmfit2$coef[, 1:nlevels(lv)] %*% t(svaX2[, 
            1:nlevels(lv)])
        error2 <- data - signal2
        rownames(signal2) <- rownames(error2) <- rownames(datacor) <- rownames(signal) <- rownames(batch) <- rownames(error) <- rownames(data)
        colnames(signal2) <- colnames(error2) <- colnames(datacor) <- colnames(signal) <- colnames(batch) <- colnames(error) <- colnames(data)
        
        # find the peaks with significant differences by F test
        # with BH correction for fdt control with surrogate
        # variables
        modSv = cbind(mod, svafit$sv)
        mod0Sv = cbind(mod0, svafit$sv)
        pValuesSv = sva::f.pvalue(data, modSv, mod0Sv)
        qValuesSv = stats::p.adjust(pValuesSv, method = "BH")
        # find the peaks with significant differences by F test
        # with BH correction for fdt control without surrogate
        # variables
        pValues = sva::f.pvalue(data, mod, mod0)
        qValues = stats::p.adjust(pValues, method = "BH")
        # get the results as list
        li <- list(data, datacor, signal, batch, error, 
            signal2, error2, pValues, qValues, pValuesSv, 
            qValuesSv, svafit$pprob.gam, svafit$pprob.b)
        names(li) <- c("data", "dataCorrected", "signal", 
            "batch", "error", "signal2", "error2", "p-values", 
            "q-values", "p-valuesCorrected", "q-valuesCorrected", 
            "PosteriorProbabilitiesSurrogate", "PosteriorProbabilitiesMod")
        message("Done!")
    }
    return(li)
}

#' Use Independent Surrogate Variable Analysis(ISVA) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details this is used for reviesed version of ISVA to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, signal part, random errors part, batch part, p-values, q-values, mass, rt. If no surrogate variable found, corresponding part would miss.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' list <- isvacor(list$data,list$group$class)
#' }
#' @seealso \code{\link{svacor}},\code{\link{svapca}}, \code{\link{svaplot}}, \code{\link{svabatch}}
#' @export
isvacor <- function(data, lv) {
    mod <- stats::model.matrix(~lv)
    mod0 <- as.matrix(c(rep(1, ncol(data))))
    
    isvafit <- isva::DoISVA(data, lv, factor.log = T)
    if (isvafit$nsv == 0) {
        message("No surrogate variable found")
        svaX <- stats::model.matrix(~lv)
        lmfit <- limma::lmFit(data, svaX)
        signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 
            1:nlevels(lv)])
        error <- data - signal
        rownames(signal) <- rownames(error) <- rownames(data)
        colnames(signal) <- colnames(error) <- colnames(data)
        # find the peaks with significant differences by F test
        # with BH correction for fdt control
        pValues = sva::f.pvalue(data, mod, mod0)
        qValues = stats::p.adjust(pValues, method = "BH")
        # get the results as list
        li <- list(data, signal, error, pValues, qValues)
        names(li) <- c("data", "signal", "error", "p-values", 
            "q-values")
        
    } else {
        message("Data is correcting ...")
        svaX <- stats::model.matrix(~lv + isvafit$isv)
        lmfit <- limma::lmFit(data, svaX)
        # data decomposition with sv
        batch <- lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) + 
            isvafit$nsv)] %*% t(svaX[, (nlevels(lv) + 1):(nlevels(lv) + 
            isvafit$nsv)])
        signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 
            1:nlevels(lv)])
        error <- data - signal - batch
        datacor <- signal + error
        # data decomposition without sv
        svaX2 <- stats::model.matrix(~lv)
        lmfit2 <- limma::lmFit(data, svaX2)
        signal2 <- lmfit2$coef[, 1:nlevels(lv)] %*% t(svaX2[, 
            1:nlevels(lv)])
        error2 <- data - signal2
        rownames(signal2) <- rownames(error2) <- rownames(datacor) <- rownames(signal) <- rownames(batch) <- rownames(error) <- rownames(data)
        colnames(signal2) <- colnames(error2) <- colnames(datacor) <- colnames(signal) <- colnames(batch) <- colnames(error) <- colnames(data)
        
        modSv = cbind(mod, isvafit$isv)
        mod0Sv = cbind(mod0, isvafit$isv)
        pValuesSv = sva::f.pvalue(data, modSv, mod0Sv)
        qValuesSv = stats::p.adjust(pValuesSv, method = "BH")
        
        pValues = sva::f.pvalue(data, mod, mod0)
        qValues = stats::p.adjust(pValues, method = "BH")
        li <- list(data, datacor, signal, batch, error, 
            signal2, error2, pValues, qValues, pValuesSv, 
            qValuesSv)
        names(li) <- c("data", "dataCorrected", "signal", 
            "batch", "error", "signal2", "error2", "p-values", 
            "q-values", "p-valuesCorrected", "q-valuesCorrected")
        message("Done!")
    }
    return(li)
}

#' Principal component analysis(PCA) for SVA/ISVA corrected data and raw data
#' @param list results from `svacor` or `isvacor` function
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param lv factor vector for the group infomation
#' @return plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- svacor(list$data,list$group$class)
#' svapca(li)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{svaplot}}, \code{\link{svabatch}}
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
        pch = as.character(lv)
    }
    
    graphics::par(mfrow = c(2, 5), mar = c(4, 4, 2.6, 1))
    
    pcao <- stats::prcomp(t(data), center = center, scale = scale)
    pcaoVars = signif(((pcao$sdev)^2)/(sum((pcao$sdev)^2)), 
        3) * 100
    graphics::plot(pcao, type = "l", main = "PCA")
    
    pca <- stats::prcomp(t(Signal), center = TRUE, scale = TRUE)
    pcaVars = signif(((pca$sdev)^2)/(sum((pca$sdev)^2)), 
        3) * 100
    graphics::plot(pca, type = "l", main = "PCA-signal")
    
    pcab <- stats::prcomp(t(Batch), center = center, scale = scale)
    pcabVars = signif(((pcab$sdev)^2)/(sum((pcab$sdev)^2)), 
        3) * 100
    graphics::plot(pcab, type = "l", main = "PCA-batch")
    
    pcae <- stats::prcomp(t(error), center = center, scale = scale)
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

#' Filter the data with p value and q value and show them
#' @param list results from `svacor` or `isvacor` function
#' @param lv factor vector for the group infomation
#' @param pqvalues method for ANOVA or SVA
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @param index index for selected peaks
#' @return heatmap for the data
#' @examples
#' \dontrun{
#' sim <- mzrtsim()
#' li <- svacor(log(sim$data), as.factor(sim$con))
#' svaplot(li,as.factor(sim$con))
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{svapca}}, \code{\link{svabatch}}
#' @export
svaplot <- function(list, lv, pqvalues = "sv", pt = 0.05, 
    qt = 0.05, index = NULL) {
    data <- list$data[, order(lv)]
    signal <- list$signal[, order(lv)]
    signal2 <- list$signal2[, order(lv)]
    batch <- list$batch[, order(lv)]
    error <- list$error[, order(lv)]
    error2 <- list$error2[, order(lv)]
    datacor <- list$dataCorrected[, order(lv)]
    pValues <- list$"p-values"
    qValues <- list$"q-values"
    pValuesSv <- list$"p-valuesCorrected"
    qValuesSv <- list$"q-valuesCorrected"
    if (!is.null(index)) {
        data <- data[index, order(lv)]
        signal <- signal[index, order(lv)]
        signal2 <- signal2[index, order(lv)]
        batch <- batch[index, order(lv)]
        error <- error[index, order(lv)]
        error2 <- error2[index, order(lv)]
        datacor <- datacor[index, order(lv)]
        pValues <- pValues[index]
        qValues <- qValues[index]
        pValuesSv <- pValuesSv[index]
        qValuesSv <- qValuesSv[index]
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
    plotimage1 <- function(data, signal, error, zlim) {
        graphics::image(t(data), col = icolors, xlab = "samples", 
            main = "peaks", xaxt = "n", yaxt = "n", zlim = zlim)
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
            main = "peaks", xaxt = "n", yaxt = "n", zlim = zlim)
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
        
        graphics::image(t(datacor), col = icolors, xlab = "samples", 
            main = "peaks-corrected", xaxt = "n", yaxt = "n", 
            zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv))
        graphics::axis(2, at = seq(0, 1, 1/(nrow(datacor) - 
            1)), labels = rownames(datacor), las = 1)
        graphics::abline(v = posv)
        
    }
    
    if (is.null(signal2)) {
        if (pqvalues == "anova" & sum(pValues < pt & qValues < 
            qt) != 0) {
            message("No SV while p-values and q-values have results")
            graphics::layout(matrix(rep(c(1, 1, 2, 2, 3, 
                3, 4, 4, 5), 9), 9, 9, byrow = TRUE))
            data <- data[pValues < pt & qValues < qt, ]
            signal <- signal[pValues < pt & qValues < qt, 
                ]
            error <- error[pValues < pt & qValues < qt, 
                ]
            
            zlim <- range(c(data, signal, error))
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage1(data, signal, error, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
            
            li <- list(data, pValues < pt & qValues < qt)
            names(li) <- c("data", "pqvalues")
            return(li)
        } else {
            message("No SV while p-values and q-values have no results")
            graphics::layout(matrix(rep(c(1, 1, 1, 2, 2, 
                3, 3, 4), 8), 8, 8, byrow = TRUE))
            zlim <- range(c(data, signal, error))
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage1(data, signal, error, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
        }
    } else {
        if (pqvalues == "anova" & sum(pValues < pt & qValues < 
            qt) != 0) {
            message("Have SVs while p-values and q-values have results")
            graphics::layout(matrix(rep(c(1, 1, 2, 2, 3, 
                3, 4), 7), 7, 7, byrow = TRUE))
            data <- data[pValues < pt & qValues < qt, ]
            signal <- signal2[pValues < pt & qValues < 
                qt, ]
            error <- error2[pValues < pt & qValues < qt, 
                ]
            zlim <- range(c(data, signal, error))
            
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage1(data, signal, error, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
            li <- list(data, pValues < pt & qValues < qt)
            names(li) <- c("data", "pqvalues")
            return(li)
        } else if (pqvalues == "anova") {
            message("Have SVs while p-values and q-values have no results")
            graphics::layout(matrix(rep(c(1, 1, 1, 2, 2, 
                3, 3, 4), 8), 8, 8, byrow = TRUE))
            zlim <- range(c(data, signal2, error2))
            
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage1(data, signal2, error2, zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
        } else if (pqvalues == "sv" & sum(pValuesSv < pt & 
            qValuesSv < qt) != 0) {
            message("SVs corrected while p-values and q-values have results")
            graphics::layout(matrix(rep(c(1, 1, 2, 2, 3, 
                3, 4, 4, 5, 5, 6), 11), 11, 11, byrow = TRUE))
            data <- data[pValuesSv < pt & qValuesSv < qt, 
                ]
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
            plotimage2(data, signal, batch, error, datacor, 
                zlim)
            graphics::par(mar = c(3, 1, 2, 5))
            plotchange(zlim)
            li <- list(datacor, data, pValuesSv < pt & 
                qValuesSv < qt)
            names(li) <- c("dataCorrected", "data", "pqvalues")
            return(li)
        } else {
            message("SVs corrected while p-values and q-values have no results")
            graphics::layout(matrix(rep(c(1, 1, 1, 2, 2, 
                3, 3, 4, 4, 5, 5, 5, 6), 13), 13, 13, byrow = TRUE))
            zlim <- range(c(signal, data, batch, error, 
                datacor))
            graphics::par(mar = c(3, 6, 2, 1))
            plotimage2(data, signal, batch, error, datacor, 
                zlim)
            graphics::par(mar = c(3, 1, 2, 6))
            plotchange(zlim)
        }
    }
}


#' Plot the influnces of DoE and Batch effects on each peaks
#' @param df data output from `svacor` or `isvacor` function
#' @param dfsv data output from `svaplot` function for corrected data
#' @param dfanova data output from `svaplot` function for raw data
#' @return influnces plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- svacor(list$data,list$group$class)
#' dfsv <- svaplot(li,list$group$class)
#' dfanova <- svaplot(li, list$group$class, pqvalues = 'anova')
#' svabatch(li,dfsv,dfanova)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{svaplot}}, \code{\link{svapca}}
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
        0] ~ df$PosteriorProbabilitiesMod[big < 0], xlab = "Influences  from primary variable", 
        ylab = "Influences  from surrogate variable", main = "Peaks with larger p-values", 
        pch = 19, cex = 1.2, ylim = c(0, 1), xlim = c(0, 
            1.02), col = grDevices::rgb(0, 0, 0, 0.1))
    graphics::points(df$PosteriorProbabilitiesSurrogate[big < 
        0 & dfsv$pqvalues & dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[big < 
        0 & dfsv$pqvalues & dfanova$pqvalue], pch = 19, 
        col = "red")
    graphics::points(df$PosteriorProbabilitiesSurrogate[big < 
        0 & !dfsv$pqvalues & dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[big < 
        0 & !dfsv$pqvalues & dfanova$pqvalue], pch = 19, 
        col = "green")
    graphics::plot(df$PosteriorProbabilitiesSurrogate[big > 
        0] ~ df$PosteriorProbabilitiesMod[big > 0], xlab = "Influences  from primary variable", 
        ylab = "Influences  from surrogate variable", main = "Peaks with smaller p-values", 
        pch = 19, cex = 1.2, ylim = c(0, 1), xlim = c(0, 
            1.02), col = grDevices::rgb(0, 0, 0, 0.1))
    graphics::points(df$PosteriorProbabilitiesSurrogate[big > 
        0 & dfsv$pqvalues & dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[big > 
        0 & dfsv$pqvalues & dfanova$pqvalue], pch = 19, 
        col = "red")
    graphics::points(df$PosteriorProbabilitiesSurrogate[big > 
        0 & dfsv$pqvalues & !dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[big > 
        0 & dfsv$pqvalues & !dfanova$pqvalue], pch = 19, 
        col = "blue")
    
}

#' Relative Log Abundance (RLA) plots
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @param type 'g' means group median based, other means all samples median based.
#' @return Relative Log Abundance (RLA) plots
#' @examples
#' sim <- mzrtsim()
#' rlaplot(sim$data, as.factor(sim$con))
#' @export
rlaplot <- function(data, lv, type = "g") {
    data <- log(data)
    data[is.nan(data)] <- 0
    outmat = NULL
    
    if (type == "g") {
        for (lvi in levels(lv)) {
            submat <- data[, lv == lvi]
            median <- apply(submat, 1, median)
            tempmat <- sweep(submat, 1, median, "-")
            outmat <- cbind(outmat, tempmat)
        }
    } else {
        median <- apply(data, 1, median)
        outmat <- sweep(data, 1, median, "-")
        
    }
    
    outmat <- outmat[, order(lv)]
    graphics::boxplot(outmat, col = lv[order(lv)])
    graphics::abline(h = 0)
}
#' Relative Log Abundance Ridge(RLAR) plots
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @param type 'g' means group median based, other means all samples median based.
#' @return Relative Log Abundance (RLA) plots
#' @examples
#' sim <- mzrtsim()
#' ridgesplot(sim$data, as.factor(sim$con))
#' @export
ridgesplot <- function(data, lv, type = "g") {
    data <- log(data)
    data[is.nan(data)] <- 0
    outmat = NULL
    
    if (type == "g") {
        for (lvi in levels(lv)) {
            submat <- data[, lv == lvi]
            median <- apply(submat, 1, median)
            tempmat <- sweep(submat, 1, median, "-")
            outmat <- cbind(outmat, tempmat)
        }
    } else {
        median <- apply(data, 1, median)
        outmat <- sweep(data, 1, median, "-")
        
    }
    
    ov <- reshape2::melt(outmat)
    colnames(outmat) <- lv[order(lv)]
    ov2 <- reshape2::melt(outmat)
    ov2$group <- as.factor(ov2$Var2)
    ggplot2::ggplot(ov, ggplot2::aes(x = ov$value, y = ov$Var2, 
        fill = ov2$group)) + ggridges::geom_density_ridges(stat = "binline", 
        bins = 100) + ggplot2::xlim(-0.5, 0.5) + ggplot2::scale_fill_discrete(name = "Group") + 
        ggplot2::labs(x = "Relative Log Abundance", y = "Samples")
}



