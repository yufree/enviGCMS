#' plot the scatter plot for peaks list with threshold
#' @param list list with data as peaks list, mz, rt and group information
#' @param rt vector range of the retention time
#' @param ms vector vector range of the m/z
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group, default 30
#' @param imputation parameters for `getimputation` function method
#' @param ... parameters for `plot` function
#' @return data fit the cutoff
#' @examples
#' data(list)
#' plotmr(list)
#' @export
plotmr <- function(list,
                   rt = NULL,
                   ms = NULL,
                   inscf = 5,
                   rsdcf = 30,
                   imputation = "l",
                   ...) {
        graphics::par(mar = c(5, 4.2, 6.1, 2.1), xpd = TRUE)
        list <- getdoe(list,
                       rsdcf = rsdcf,
                       inscf = inscf,
                       imputation = imputation)
        lif <-
                getfilter(list, rowindex = list$rsdindex & list$insindex)
        data <- lif$groupmean
        dataname <- colnames(data)
        mz <- lif$mz
        RT <- lif$rt
        suppressWarnings(if (sum(nrow(data) > 0) > 0) {
                n <- dim(data)[2]
                col <- grDevices::rainbow(n, alpha = 0.318)

                graphics::plot(
                        mz ~ RT,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        type = "n",
                        pch = 19,
                        ylim = ms,
                        xlim = rt,
                        ...
                )

                for (i in 1:n) {
                        cex = as.numeric(cut((log10(data[, i] + 1) -
                                                      inscf),
                                             breaks = c(0, 1, 2, 3, 4, Inf) / 2
                        )) / 2
                        cexlab = c(
                                paste0(inscf, "-", inscf + 0.5),
                                paste0(inscf + 0.5, "-", inscf + 1),
                                paste0(inscf +
                                               1, "-", inscf + 1.5),
                                paste0(inscf +
                                               1.5, "-", inscf + 2),
                                paste0(">", inscf +
                                               2)
                        )
                        if (!is.null(ms) & !is.null(rt)) {
                                graphics::points(
                                        x = RT[RT > rt[1] & RT <
                                                       rt[2] &
                                                       mz > ms[1] &
                                                       mz < ms[2]],
                                        y = mz[RT >
                                                       rt[1] &
                                                       RT < rt[2] &
                                                       mz > ms[1] & mz <
                                                       ms[2]],
                                        cex = cex,
                                        col = col[i],
                                        pch = 19
                                )
                        } else {
                                graphics::points(
                                        x = RT,
                                        y = mz,
                                        cex = cex,
                                        col = col[i],
                                        pch = 19
                                )
                        }

                }
                graphics::legend(
                        "topright",
                        legend = dataname,
                        col = col,
                        pch = 19,
                        horiz = T,
                        bty = "n",
                        inset = c(0,-0.25)
                )
                graphics::legend(
                        "topleft",
                        legend = cexlab,
                        title = "Intensity in Log scale",
                        pt.cex = c(1, 2, 3, 4, 5) / 4,
                        pch = 19,
                        bty = "n",
                        horiz = T,
                        cex = 0.7,
                        col = grDevices::rgb(0,
                                             0, 0, 0.318),
                        inset = c(0,-0.25)
                )
        } else if (NROW(data) > 0) {
                graphics::plot(
                        mz ~ RT,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        type = "n",
                        pch = 19,
                        ylim = ms,
                        xlim = rt,
                        ...
                )
                cex = as.numeric(cut((log10(
                        data + 1
                ) -
                        inscf), breaks = c(0, 1, 2, 3, 4, Inf) /
                        2)) / 2
                cexlab = c(
                        paste0(inscf, "-", inscf + 0.5),
                        paste0(inscf + 0.5, "-", inscf + 1),
                        paste0(inscf +
                                       1, "-", inscf + 1.5),
                        paste0(inscf +
                                       1.5, "-", inscf + 2),
                        paste0(">", inscf +
                                       2)
                )
                col <- grDevices::rgb(0, 0, 0, alpha = 0.318)
                if (!is.null(ms) & !is.null(rt)) {
                        graphics::points(
                                x = RT[RT > rt[1] & RT <
                                               rt[2] &
                                               mz > ms[1] &
                                               mz < ms[2]],
                                y = mz[RT >
                                               rt[1] &
                                               RT < rt[2] &
                                               mz > ms[1] & mz <
                                               ms[2]],
                                ,
                                col = col,
                                cex = cex,
                                pch = 19
                        )
                } else {
                        graphics::points(
                                x = RT,
                                y = mz,
                                cex = cex,
                                pch = 19,
                                col = col
                        )
                }
                graphics::legend(
                        "topright",
                        legend = unique(list$group),
                        col = col,
                        pch = 19,
                        horiz = T,
                        bty = "n",
                        inset = c(0,-0.25)
                )
                graphics::legend(
                        "topleft",
                        legend = cexlab,
                        title = "Intensity in Log scale",
                        pt.cex = c(1, 2, 3, 4, 5) / 4,
                        pch = 19,
                        bty = "n",
                        horiz = T,
                        cex = 0.7,
                        col = grDevices::rgb(0,
                                             0, 0, 0.318),
                        inset = c(0,-0.25)
                )

        } else
        {
                graphics::plot(
                        1,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        main = "No peaks found",
                        ylim = ms,
                        xlim = rt,
                        type = "n",
                        pch = 19,
                        ...
                )
        })
}

#' plot the diff scatter plot for one xcmsset objects with threshold between two groups
#' @param list list with data as peaks list, mz, rt and group information
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param ... parameters for `plot` function
#' @return NULL
#' @examples
#' data(list)
#' plotmrc(list)
#' @export
plotmrc <- function(list,
                    ms = c(100, 800),
                    inscf = 5,
                    rsdcf = 30,
                    imputation = "l",
                    ...) {
        list <- getdoe(
                list,
                rsdcf = rsdcf,
                inscf = inscf,
                imputation = imputation,
        )
        lif <-
                getfilter(list, rowindex = list$rsdindex & list$insindex)
        data <- lif$groupmean
        dataname <- colnames(data)
        mz <- lif$mz
        rt <- lif$rt
        suppressWarnings(if (!is.na(data[1,1])) {
                diff1 <- data[, 1] - data[, 2]
                diff2 <- data[, 2] - data[, 1]
                diff1[diff1 < 0] <- 0
                diff2[diff2 < 0] <- 0
                name1 <- paste0(dataname[1], "-", dataname[2])
                name2 <- paste0(dataname[2], "-", dataname[1])

                cex1 = as.numeric(cut((log10(
                        diff1 + 1
                ) - inscf),
                breaks = c(0, 1, 2, 3, 4, Inf) / 2)) / 2
                cex2 = as.numeric(cut((log10(
                        diff2 + 1
                ) - inscf),
                breaks = c(0, 1, 2, 3, 4, Inf) / 2)) / 2
                cexlab = c(
                        paste0(inscf, "-", inscf + 0.5),
                        paste0(inscf +
                                       0.5, "-", inscf + 1),
                        paste0(inscf + 1, "-",
                               inscf + 1.5),
                        paste0(inscf + 1.5, "-", inscf +
                                       2),
                        paste0(">", inscf + 2)
                )

                graphics::plot(
                        mz ~ rt,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        ylim = ms,
                        cex = cex1,
                        col = grDevices::rgb(0,
                                             0, 1, 0.618),
                        pch = 19,
                        ...
                )

                graphics::points(
                        mz ~ rt,
                        cex = cex2,
                        col = grDevices::rgb(1,
                                             0, 0, 0.618),
                        pch = 19
                )

                graphics::legend(
                        "topleft",
                        legend = cexlab,
                        title = "Intensity in Log scale",
                        pt.cex = c(1, 2, 3, 4, 5) / 2,
                        pch = 19,
                        col = grDevices::rgb(0,
                                             0, 0, 0.618),
                        bty = "n",
                        horiz = T,
                        inset = c(0,-0.25)
                )
                graphics::legend(
                        "topright",
                        legend = c(name1,
                                   name2),
                        pch = 19,
                        col = c(
                                grDevices::rgb(0,
                                               0, 1, 0.618),
                                grDevices::rgb(1, 0, 0, 0.618)
                        ),
                        bty = "n",
                        horiz = T,
                        inset = c(0,-0.25)
                )
        } else {
                graphics::plot(
                        1,
                        xlab = "Retention Time(s)",
                        ylab = "m/z",
                        main = "No peaks found",
                        ylim = ms,
                        type = "n",
                        pch = 19,
                        ...
                )
        })

}

#' plot the rsd influnces of data in different groups
#' @param list list with data as peaks list, mz, rt and group information
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param ... other parameters for `plot` function
#' @return NULL
#' @examples
#' data(list)
#' plotrsd(list)
#' @export
plotrsd <- function(list,
                    ms = c(100, 800),
                    inscf = 5,
                    rsdcf = 100,
                    imputation = "l",
                    ...) {
        graphics::par(mar = c(5, 4.2, 6.1, 2.1), xpd = TRUE)
        cexlab = c("<20%", "20-40%", "40-60%", "60-80%", ">80%")
        list <- getdoe(list,
                       rsdcf = rsdcf,
                       inscf = inscf,
                       imputation = imputation)
        lif <-
                getfilter(list, rowindex = list$rsdindex & list$insindex)
        data <- lif$groupmean
        dataname <- colnames(data)
        mz <- lif$mz
        rt <- lif$rt
        rsd <- lif$grouprsd

        n <- dim(rsd)[2]
        col <- grDevices::rainbow(n, alpha = 0.318)

        graphics::plot(
                mz ~ rt,
                xlab = "Retention Time(s)",
                ylab = "m/z",
                main = "RSD(%) distribution",
                type = "n",
                pch = 19,
                ...
        )

        for (i in 1:n) {
                cex = as.numeric(cut(rsd[, i], breaks = c(0, 20,
                                                          40, 60, 80, Inf))) / 2
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
                "topright",
                legend = dataname,
                col = col,
                horiz = T,
                pch = 19,
                bty = "n",
                inset = c(0,-0.25)
        )
        graphics::legend(
                "topleft",
                legend = cexlab,
                pt.cex = c(1,
                           2, 3, 4, 5) / 2,
                pch = 19,
                bty = "n",
                horiz = T,
                cex = 0.8,
                col = grDevices::rgb(0, 0, 0, 0.318),
                inset = c(0,-0.25)
        )
}

#' plot scatter plot for rt-mz profile and output gif file for mutiple groups
#' @param list list with data as peaks list, mz, rt and group information
#' @param name file name for gif file, default test
#' @param ms the mass range to plot the data
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param ... parameters for `plot` function
#' @return gif file
#' @examples
#' \dontrun{
#' data(list)
#' gifmr(list)
#' }
#' @export
gifmr <- function(list,
                  ms = c(100, 500),
                  rsdcf = 30,
                  inscf = 5,
                  imputation = "i",
                  name = "test",
                  ...) {
        list <- getdoe(list,
                       rsdcf = rsdcf,
                       inscf = inscf,
                       imputation = imputation)
        lif <-
                getfilter(list, rowindex = list$rsdindex & list$insindex)
        data <- lif$groupmean
        mz <- lif$mz
        rt <- lif$rt
        filename = paste0(name, ".gif")
        mean <- apply(data, 1, mean)

        graphics::plot(
                mz ~ rt,
                xlab = "Retention Time(s)",
                ylab = "m/z",
                pch = 19,
                cex = log10(mean + 1) -
                        4,
                col = grDevices::rgb(0, 0, 1, 0.2),
                main = "All peaks",
                ...
        )

        col <- grDevices::rainbow(dim(data)[2], alpha = 0.318)
        animation::saveGIF({
                for (i in 1:dim(data)[2]) {
                        name <- colnames(data)[i]
                        value <- data[, i]
                        graphics::plot(
                                mz ~ rt,
                                xlab = "Retention Time(s)",
                                ylab = "m/z",
                                pch = 19,
                                cex = log10(value +
                                                    1) - 4,
                                col = col[i],
                                ylim = ms,
                                main = name,
                                ...
                        )
                }
        }, movie.name = filename, ani.width = 800, ani.height = 500)
}

#' plot the PCA for multiple samples
#' @param data mzrt profile with row peaks and column samples
#' @param lv group information
#' @param index index for selected peaks
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param xrange x axis range for return samples, default NULL
#' @param yrange y axis range for return samples, default NULL
#' @param pch deault pch would be the first charactor of group information or samples name
#' @param ... other parameters for `plot` function
#' @return if xrange and yrange are not NULL, return file name of all selected samples on 2D score plot
#' @examples
#' data(list)
#' plotpca(list$data, lv = as.character(list$group))
#' @export
plotpca <- function(data,
                    lv = NULL,
                    index = NULL,
                    center = T,
                    scale = T,
                    xrange = NULL,
                    yrange = NULL,
                    pch = NULL,
                    ...) {
        if (!is.null(index)) {
                data <- data[index,]
        }

        if (is.null(lv)) {
                pch0 = colnames(data)
        } else {
                pch0 = lv
        }
        pcao <-
                stats::prcomp(t(data), center = center, scale = scale)
        pcaoVars = signif(((pcao$sdev) ^ 2) / (sum((pcao$sdev) ^ 2)),
                          3) * 100
        if(!is.null(pch)){
                graphics::plot(
                        pcao$x[, 1],
                        pcao$x[, 2],
                        xlab = paste("PC1:",
                                     pcaoVars[1], "% of Variance Explained"),
                        ylab = paste("PC2:",
                                     pcaoVars[2], "% of Variance Explained"),
                        cex = 2,
                        pch = pch,
                        ...
                )
        }else{
                graphics::plot(
                        pcao$x[, 1],
                        pcao$x[, 2],
                        xlab = paste("PC1:",
                                     pcaoVars[1], "% of Variance Explained"),
                        ylab = paste("PC2:",
                                     pcaoVars[2], "% of Variance Explained"),
                        cex = 2,
                        pch = pch0,
                        ...
                )
        }
        # pch = ifelse(hasArg(pch),pch,pch0)

        if (!(is.null(xrange) & is.null(yrange))) {
                return(colnames(data)[pcao$x[, 1] > xrange[1] &
                                              pcao$x[, 1] < xrange[2] &
                                              pcao$x[, 2] > yrange[1] &
                                              pcao$x[, 2] < yrange[2]])
        }
}

#' Plot the heatmap of mzrt profiles
#' @param data mzrt profile with row peaks and column samples
#' @param lv group information
#' @param index index for selected peaks
#' @return NULL
#' @examples
#' data(list)
#' plothm(list$data, lv = as.factor(list$group))
#' @export
plothm <- function(data, lv, index = NULL) {
        icolors <-
                (grDevices::colorRampPalette(rev(
                        RColorBrewer::brewer.pal(11,
                                                 "RdYlBu")
                )))(100)
        zlim <- range(data)
        if (!is.null(index)) {
                data <- data[index, order(lv)]
        } else {
                data <- data[, order(lv)]
        }
        plotchange <- function(zlim) {
                breaks <- seq(zlim[1], zlim[2], round((zlim[2] -
                                                               zlim[1]) / 10))
                poly <- vector(mode = "list", length(icolors))
                graphics::plot(
                        1,
                        1,
                        t = "n",
                        xlim = c(0, 1),
                        ylim = zlim,
                        xaxt = "n",
                        yaxt = "n",
                        xaxs = "i",
                        yaxs = "i",
                        ylab = "",
                        xlab = "",
                        frame.plot = F
                )
                graphics::axis(
                        4,
                        at = breaks,
                        labels = round(breaks),
                        las = 1,
                        pos = 0.4,
                        cex.axis = 0.8
                )
                p <- graphics::par("usr")
                graphics::text(
                        p[2] + 2,
                        mean(p[3:4]),
                        labels = "intensity",
                        xpd = NA,
                        srt = -90
                )
                bks <-
                        seq(zlim[1], zlim[2], length.out = (length(icolors) +
                                                                    1))
                for (i in seq(poly)) {
                        graphics::polygon(
                                c(0.1, 0.1, 0.3, 0.3),
                                c(bks[i],
                                  bks[i + 1], bks[i + 1], bks[i]),
                                col = icolors[i],
                                border = NA
                        )
                }
        }
        pos <- cumsum(as.numeric(table(lv) / sum(table(lv)))) -
                as.numeric(table(lv) / sum(table(lv))) / 2
        posv <-
                cumsum(as.numeric(table(lv) / sum(table(lv))))[1:(nlevels(lv) -
                                                                          1)]


        graphics::layout(matrix(rep(c(1, 1, 1, 2), 10), 10,
                                4, byrow = TRUE))
        graphics::par(mar = c(3, 6, 2, 1))
        graphics::image(
                t(data),
                col = icolors,
                xlab = "samples",
                main = "peaks intensities on log scale",
                xaxt = "n",
                yaxt = "n",
                zlim = zlim
        )
        graphics::axis(1,
                       at = pos,
                       labels = levels(lv),
                       cex.axis = 0.8)
        graphics::axis(
                2,
                at = seq(0, 1, 1 / (nrow(data) - 1)),
                labels = rownames(data),
                cex.axis = 1,
                las = 2
        )
        graphics::abline(v = posv)
        graphics::par(mar = c(3, 1, 2, 6))
        plotchange(zlim)
}

#' plot the density for multiple samples
#' @param data mzrt profile with row peaks and column samples
#' @param lv group information
#' @param index index for selected peaks
#' @param name name on the figure for samples
#' @param lwd the line width for density plot, default 1
#' @param ... parameters for `plot` function
#' @return NULL
#' @examples
#' data(list)
#' plotden(list$data, lv = as.character(list$group),ylim = c(0,1))
#' @export
plotden <- function(data,
                    lv,
                    index = NULL,
                    name = NULL,
                    lwd = 1,
                    ...) {
        if (!is.null(index)) {
                data <- data[index,]
        }
        xlim <- range(log10(data + 1), na.rm = T)
        if (is.null(lv)) {
                col <- as.numeric(as.factor(colnames(data)))
                coli <- unique(colnames(data))
        } else {
                col <- as.numeric(as.factor(lv))
                coli <- unique(lv)
        }
        graphics::plot(
                1,
                1,
                typ = 'n',
                main = paste0('Metabolites profiles changes in ',
                              name,
                              ' samples'),
                xlab = 'Intensity(log based 10)',
                ylab = 'Density',
                xlim = c(xlim[1], xlim[2] + 1),
                ...
        )
        for (i in 1:(ncol(data))) {
                graphics::lines(stats::density(log10(data[, i] + 1)),
                                col = col[i],
                                lwd = lwd)
        }
        graphics::legend(
                "topright",
                legend = coli,
                col = unique(col),
                pch = 19,
                bty = "n"
        )
}
#' Relative Log Abundance (RLA) plots
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @param type 'g' means group median based, other means all samples median based.
#' @return Relative Log Abundance (RLA) plots
#' @examples
#' data(list)
#' plotrla(list$data, as.factor(list$group))
#' @export
plotrla <- function(data, lv, type = "g") {
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
#' Relative Log Abundance Ridge (RLAR) plots
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @param type 'g' means group median based, other means all samples median based.
#' @return Relative Log Abundance Ridge(RLA) plots
#' @examples
#' data(list)
#' plotridges(list$data, as.factor(list$group))
#' @export
plotridges <- function(data, lv, type = "g") {
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
        ggplot2::ggplot(ov,
                        ggplot2::aes(
                                x = ov$value,
                                y = ov$Var2,
                                fill = ov2$group
                        )) + ggridges::geom_density_ridges(stat = "binline",
                                                           bins = 100) + ggplot2::xlim(-0.5, 0.5) + ggplot2::scale_fill_discrete(name = "Group") +
                ggplot2::labs(x = "Relative Log Abundance", y = "Samples")
}
#' plot density weighted intensity for multiple samples
#' @param list list with data as peaks list, mz, rt and group information
#' @param n the number of equally spaced points at which the density is to be estimated, default 512
#' @param ... parameters for `plot` function
#' @return Density weighted intensity for multiple samples
#' @examples
#' data(list)
#' plotdwtus(list)
#' @export
plotdwtus <- function(list,n=512,...){
        dwtus <- apply(list$data, 2, function(x) getdwtus(x,n = n))
        if(!is.null(list$order)){
                graphics::plot(dwtus~as.numeric(list$order), xlab='Run order', ylab = 'DWTUS', col = as.numeric(as.factor(list$group)),...)
                graphics::legend('topright',legend = unique(list$group),col = unique(as.numeric(as.factor(list$group))),pch = 19,bty = 'n')
        }else{
                graphics::plot(dwtus~as.numeric(as.factor(list$group)),xlab='Group', ylab = 'DWTUS',col = as.numeric(as.factor(list$group)),...)
                graphics::legend('topright',legend = unique(list$group),col = unique(as.numeric(as.factor(list$group))),pch = 19,bty = 'n')
        }
}
