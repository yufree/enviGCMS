#' plot the scatter plot for xcmsset objects with threshold
#' @param xset the xcmsset object
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks, default 5
#' @param rsdcf the rsd cutoff of all peaks
#' @param ... parameters for `plot` function
#' @return data fit the cutoff
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plotmr(xset)
#' }
#' @export
plotmr <- function(xset,
                   ms = c(100, 800),
                   inscf = 5,
                   rsdcf = 30,
                   ...) {
        graphics::par(mar=c(5, 4.2, 6.1, 2.1), xpd=TRUE)
        data <- getbiorep(xset, rsdcf = rsdcf, inscf = inscf)
        suppressWarnings(if (!is.na(data)) {
                datamean <- data[, grepl('*mean', colnames(data))]
                dataname <- unique(xcms::sampclass(xset))
                n <- dim(datamean)[2]
                col <- grDevices::rainbow(n, alpha = 0.318)

                graphics::plot(
                        data$mzmed ~ data$rtmed,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        type = 'n',
                        pch = 19,
                        ...
                )

                for (i in 1:n) {
                        cex = as.numeric(cut((log10(datamean[, i] + 1)-inscf), breaks=c(0,1,2,3,4,Inf)/2))/2
                        cexlab = c(paste0(inscf,'-',inscf+0.5),paste0(inscf+0.5,'-',inscf+1),paste0(inscf+1,'-',inscf+1.5),paste0(inscf+1.5,'-',inscf+2),paste0('>',inscf+2))
                        graphics::points(
                                x = data$rtmed,
                                y = data$mzmed,
                                cex = cex,
                                col = col[i],
                                pch = 19,
                                ylim = ms
                        )
                }
                graphics::legend(
                        'topright',
                        legend = dataname,
                        col = col,
                        pch = 19,
                        horiz = T,
                        bty = 'n',
                        inset = c(0,-0.25)
                )
                graphics::legend(
                        'topleft',
                        legend = cexlab,
                        title = 'Intensity in Log scale',
                        pt.cex = c(1,2,3,4,5)/2,
                        pch = 19,
                        bty = 'n',
                        horiz = T,
                        cex = 0.7,
                        col = grDevices::rgb(0,0,0,0.318),inset = c(0,-0.25)
                )
        } else{
                graphics::plot(
                        1,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        main = "No peaks found",
                        ylim = ms,
                        type = 'n',
                        pch = 19,
                        ...
                )
        })
}

#' plot the diff scatter plot for one xcmsset objects with threshold and two groups
#' @param xset xcmsset object with two groups
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks, default 5
#' @param rsdcf the rsd cutoff of all peaks
#' @param ... parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plotmrc(xset)
#' }
#' @export
plotmrc <- function(xset,
                    ms = c(100, 800),
                    inscf = 5,
                    rsdcf = 30,
                    ...)  {
        data <- getbiorep(xset, rsdcf = rsdcf, inscf = inscf)
        suppressWarnings(if (!is.na(data)) {
                datamean <- data[, grepl('*mean', colnames(data))]
                dataname <- unique(xcms::sampclass(xset))

                diff1 <- datamean[, 1] - datamean[, 2]
                diff2 <- datamean[, 2] - datamean[, 1]
                diff1[diff1 < 0] <- 0
                diff2[diff2 < 0] <- 0
                name1 <- paste0(dataname[1], "-", dataname[2])
                name2 <- paste0(dataname[2], "-", dataname[1])

                cex1 = as.numeric(cut((log10(diff1 + 1)-inscf), breaks=c(0,1,2,3,4,Inf)/2))/2
                cex2 = as.numeric(cut((log10(diff2 + 1)-inscf), breaks=c(0,1,2,3,4,Inf)/2))/2
                cexlab = c(paste0(inscf,'-',inscf+0.5),paste0(inscf+0.5,'-',inscf+1),paste0(inscf+1,'-',inscf+1.5),paste0(inscf+1.5,'-',inscf+2),paste0('>',inscf+2))

                graphics::plot(
                        data$mzmed ~ data$rtmed,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        ylim = ms,
                        cex = cex1,
                        col = grDevices::rgb(0, 0, 1, 0.618),
                        pch = 19,
                        ...
                )

                graphics::points(
                        data$mzmed ~ data$rtmed,
                        cex = cex2,
                        col = grDevices::rgb(1, 0, 0, 0.618),
                        pch = 19
                )

                graphics::legend(
                        'top',
                        legend = cexlab,
                        title = 'Intensity in Log scale',
                        pt.cex = c(1,2,3,4,5)/2,
                        pch = 19,
                        col = grDevices::rgb(0, 0, 0, 0.618),
                        bty = 'n',
                        horiz = T
                )
                graphics::legend(
                        'bottom',
                        legend = c(name1, name2),
                        pch = 19,
                        col = c(
                                grDevices::rgb(0, 0, 1, 0.618),
                                grDevices::rgb(1, 0, 0, 0.618)
                        ),
                        bty = 'n',
                        horiz = T
                )
        } else{
                graphics::plot(
                        1,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        main = "No peaks found",
                        ylim = ms,
                        type = 'n',
                        pch = 19,
                        ...
                )
        })

}

#' plot the rsd influnces of data
#' @param xset xcmsset data
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks, default 5
#' @param rsdcf the rsd cutoff of all peaks
#' @param type 't' means biological samples with technique replicates; 'b' means biological samples without technique replicates in one group; 'g' means biological samples without technique replicates in mutiple groups. default 'g'
#' @param ... other parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath[1:6], pmethod = ' ')
#' plotrsd(xset)
#' }
#' @export
plotrsd <- function(xset,
                    ms = c(0, 800),
                    type = 'g',
                    inscf = 5,
                    rsdcf = 30,
                    ...) {
        cexlab = c('<20%','20-40%','40-60%','60-80%','>80%')
        if (type == 't') {
                df <- gettechrep(xset, inscf = inscf, rsdcf = rsdcf)
                suppressWarnings(if (!is.na(df)) {
                        mz <- df$mzmed
                        rt <- df$rtmin
                        cex = as.numeric(cut(df$rsd, breaks=c(0,20,40,60,80,Inf)))/2
                        graphics::plot(
                                mz ~ rt,
                                cex = cex,
                                ylim = ms,
                                xlab = 'retention time',
                                ylab = 'm/z',
                                main = 'RSD(%) distribution',
                                col = grDevices::rgb(0,0,0,0.318),
                                pch = 19,
                                ...
                        )
                        graphics::legend(
                                'top',
                                legend = cexlab,
                                pt.cex = c(1,2,3,4,5)/2,
                                pch = 19,
                                bty = 'n',
                                horiz = T,
                                cex = 0.8,
                                col = grDevices::rgb(0,0,0,0.318)
                        )
                } else{
                        graphics::plot(
                                1,
                                xlab = "Retention Time(s)",
                                ylab = "m/z",
                                main = "No peaks found",
                                ylim = ms,
                                type = 'n',
                                pch = 19,
                                ...
                        )
                })
        } else if (type == 'b') {
                df <- getbiotechrep(xset, inscf = inscf, rsdcf = rsdcf)
                suppressWarnings(if (!is.na(df)) {
                        mz <- df$mzmed
                        rt <- df$rtmin
                        cex = as.numeric(cut(df$rsdB, breaks=c(0,20,40,60,80,Inf)))/2
                        graphics::plot(
                                mz ~ rt,
                                ylim = ms,
                                xlab = 'retention time',
                                ylab = 'm/z',
                                main = 'RSD(%) distribution',
                                cex = cex,
                                col = grDevices::rgb(0,0,0,0.318),
                                pch = 19,
                                ...
                        )
                        graphics::legend(
                                'top',
                                legend = cexlab,
                                pt.cex = c(1,2,3,4,5)/2,
                                cex = 0.8,
                                pch = 19,
                                bty = 'n',
                                horiz = T,
                                col = grDevices::rgb(0,0,0,0.318)

                        )
                } else{
                        graphics::plot(
                                1,
                                xlab = "Retention Time(s)",
                                ylab = "m/z",
                                main = "No peaks found",
                                ylim = ms,
                                type = 'n',
                                pch = 19,
                                ...
                        )
                })
        } else{
                df <- getbiorep(xset, inscf = inscf, rsdcf = rsdcf)
                mz <- df$mzmed
                rt <- df$rtmin
                datarsd <- df[, grepl('*rsd%', colnames(df))]
                dataname <- unique(xcms::sampclass(xset))
                n <- dim(datarsd)[2]
                col <- grDevices::rainbow(n, alpha = 0.318)

                graphics::plot(
                        mz ~ rt,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        main = "RSD(%) distribution",
                        ylim = ms,
                        type = 'n',
                        pch = 19,
                        ...
                )

                for (i in 1:n) {
                        cex = as.numeric(cut(datarsd[, i], breaks=c(0,20,40,60,80,Inf)))/2
                        graphics::points(
                                x = rt,
                                y = mz,
                                ylim = ms,
                                cex = cex,
                                col = col[i],
                                pch = 19
                        )
                }
                graphics::legend(
                        'bottom',
                        legend = dataname,
                        col = col,
                        horiz = T,
                        pch = 19,
                        bty = 'n'
                )
                graphics::legend(
                        'top',
                        legend = cexlab,
                        pt.cex = c(1,2,3,4,5)/2,
                        pch = 19,
                        bty = 'n',
                        horiz = T,
                        cex = 0.8,
                        col = grDevices::rgb(0,0,0,0.318)
                )

        }


}

#' plot the PCA of xcmsset
#' @param xset xcmsset data
#' @param lv group information
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param ... other parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plotpca(xset)
#' }
#' @export

plotpca <- function(xset,
                    lv = NULL,
                    center = T,
                    scale = T,
                    ...) {
        data <- xcms::groupval(xset, 'medret', 'into')
        data <- data[stats::complete.cases(data),]
        if (is.null(lv)) {
                pch = colnames(data)
        } else {
                pch = lv
        }
        pcao <- stats::prcomp(t(data), center = center,
                              scale = scale)
        pcaoVars = signif(((pcao$sdev) ^ 2) / (sum((pcao$sdev) ^ 2)),
                          3) * 100
        graphics::plot(
                pcao$x[, 1],
                pcao$x[, 2],
                xlab = paste("PC1:",
                             pcaoVars[1], "% of Variance Explained"),
                ylab = paste("PC2:",
                             pcaoVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                ...
        )
}
#' plot EIC and boxplot for all peaks and return diffreport
#' @param xset xcmsset object
#' @param name filebase of the sub dir
#' @param test 't' means two-sample welch t-test, 't.equalvar' means two-sample welch t-test with equal variance, 'wilcoxon' means rank sum wilcoxon test, 'f' means F-test, 'pairt' means paired t test, 'blockf' means Two-way analysis of variance, default 't'
#' @param nonpara 'y' means using nonparametric ranked data, 'n' means original data
#' @param ... other parameters for `diffreport`
#' @return diffreport and pdf figure for EIC and boxplot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plote(xset)
#' }
#' @export
plote <- function(xset,
                  name = "test",
                  test = "t",
                  nonpara = "n",
                  ...) {
        gt <- xcms::groups(xset)
        a <-
                xcms::diffreport(
                        xset,
                        filebase = name,
                        eicmax = nrow(gt),
                        nonpara = nonpara,
                        ...
                )
        return(a)
}

#' plot scatter plot for rt-mz profile and output gif file for mutiple groups
#' @param xset xcmsset object with mutiple groups
#' @param file file name for further annotation, default NULL
#' @param inscf log intensity cutoff for peaks, default 5
#' @param rsdcf the rsd cutoff of all peaks
#' @param ... parameters for `plot` function
#' @return gif file and csv file
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' gifmr(xset)
#' }
#' @export
gifmr <- function(xset,
                  rsdcf = 30,
                  inscf = 5,
                  file = 'test') {
        data <- getbiorep(xset,
                          rsdcf = rsdcf,
                          inscf = inscf,
                          file = file)
        dmx <- data[,-c(1, 2)]
        mean <- apply(dmx, 1, mean)
        graphics::plot(
                data$mz ~ data$time,
                xlab = "Retention Time(s)",
                ylab = "m/z",
                pch = 19,
                cex = log10(mean + 1) - 4,
                col = grDevices::rgb(0, 0, 1, 0.2),
                main = 'All peaks'
        )
        filename = paste0(file, '.gif')
        col <- grDevices::rainbow(dim(dmx)[2], alpha = 0.318)
        animation::saveGIF({
                for (i in 1:dim(dmx)[2]) {
                        name <- colnames(dmx)[i]
                        value <- as.numeric(t(dmx)[i,])
                        graphics::plot(
                                data$mz ~ data$time,
                                xlab = "Retention Time(s)",
                                ylab = "m/z",
                                pch = 19,
                                cex = log10(value + 1) - 4,
                                col = col[i],
                                main = name
                        )
                }
        }, movie.name = filename, ani.width = 800, ani.height = 500)
}
