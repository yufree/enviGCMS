#' Demo data
#' @docType data
#' @usage data(list)
#' @format A list object with data, mass to charge ratio, retention time and group information. The list is generated from faahKO package.
"list"

#' Demo raw data matrix
#' @docType data
#' @usage data(matrix)
#' @format A matrix object from raw mass spectrometry data. The list is generated from faahKO package.
"matrix"

#' define the Mode function
#' @param x vector
#' @return Mode of the vector
#' @export
Mode <- function(x) {
        ta <- table(x)
        tam <- max(ta)
        if (all(ta == tam))
                mod <- x
        else if (is.numeric(x))
                mod <- as.numeric(names(ta)[ta == tam])
        else
                mod <- names(ta)[ta == tam]
        return(mod)
}

#' filter data by average moving box
#'
#' @param x a vector
#' @param n A number to identify the size of the moving box.
#' @return The filtered data
#' @examples
#' ma(rnorm(1000),5)
#' @export
ma <- function(x, n) {
        stats::filter(x, rep(1 / n, n), circular = TRUE)
}

#' plot GC/LC-MS data as a heatmap with TIC
#'
#' @param data imported data matrix of GC-MS
#' @param log transform the intensity into log based 10
#' @return heatmap
#' @examples
#' \dontrun{
#' png('test.png')
#' plotms(matrix)
#' dev.off()
#' }
#' @export
plotms <- function(data, log = FALSE) {
        # get the mz and rt range and rotate the matrix to
        # adapt the image function
        indmz <- as.numeric(rownames(data))
        indrt <- as.numeric(colnames(data))
        col <-
                (grDevices::colorRampPalette(rev(
                        RColorBrewer::brewer.pal(11,
                                                 "RdYlBu")
                )))(100)
        if (log) {
                z <- log10(t(data) + 1)
        } else {
                z <- t(data)
        }
        # show the intensity scale in log 10 based scale
        graphics::par(mar = c(2, 5, 1, 4), fig = c(0, 1, 0.9,
                                                   1))
        zlim <- range(z, na.rm = TRUE)
        breaks <- seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
        poly <- vector(mode = "list", length(col))
        graphics::plot(
                1,
                1,
                t = "n",
                ylim = c(0, 1),
                xlim = zlim,
                xaxt = "n",
                yaxt = "n",
                xaxs = "i",
                yaxs = "i",
                ylab = "",
                xlab = ""
        )
        graphics::mtext(
                "intensity",
                side = 2,
                line = 0.5,
                las = 1,
                cex = 1
        )
        if (log) {
                graphics::axis(
                        1,
                        at = breaks,
                        labels = round(10 ^ (breaks)),
                        las = 1
                )
        } else {
                graphics::axis(
                        1,
                        at = breaks,
                        labels = round(breaks),
                        las = 1
                )
        }

        bks <- seq(zlim[1], zlim[2], length.out = (length(col) +
                                                           1))
        for (i in seq(poly)) {
                graphics::polygon(
                        c(bks[i], bks[i + 1], bks[i +
                                                          1], bks[i]),
                        c(0, 0, 1, 1),
                        col = col[i],
                        border = NA
                )
        }
        # show the heatmap
        graphics::par(
                mar = c(4, 5, 0, 4),
                fig = c(0, 1, 0,
                        0.9),
                new = TRUE
        )
        graphics::image(
                z,
                ylab = "",
                axes = FALSE,
                col = col,
                useRaster = TRUE
        )
        # display the m/z as y
        mzy <- seq(0, 1, length.out = length(indmz))
        graphics::axis(4, at = mzy[indmz %% 100 == 0][-c(1, sum(indmz %% 100 ==
                                                                        0))], labels = c(rep("", length(indmz[indmz %% 100 ==
                                                                                                                      0])))[-c(1, sum(indmz %% 100 == 0))])
        # add the lable for double y axis
        graphics::text(
                graphics::par("usr")[2] * 1.11,
                mean(graphics::par("usr")[3:4]),
                "m/z",
                srt = -90,
                xpd = TRUE,
                pos = 4
        )
        graphics::text(
                graphics::par("usr")[2] * 1.05,
                mzy[indmz %% 100 ==
                            0][-c(1, sum(indmz %% 100 == 0))],
                labels = indmz[indmz %% 100 ==
                                       0][-c(1, sum(indmz %% 100 == 0))],
                srt = 270,
                xpd = TRUE
        )
        graphics::par(
                mar = c(4, 5, 4, 4),
                fig = c(0, 1, 0,
                        0.9),
                new = TRUE,
                cex.lab = 1
        )
        data[is.na(data)] <- 0
        data <- apply(data, 2, sum)
        graphics::plot(
                data ~ indrt,
                type = "l",
                ylab = "intensity",
                xlab = "retention time(s)",
                frame.plot = FALSE,
                xaxs = "i",
                yaxs = "i"
        )
}

#' plot GC/LC-MS data as scatter plot
#'
#' @param data imported data matrix of GC-MS
#' @param inscf Log intensity cutoff for peaks, default 3.5
#' @param ... parameters for `plot` function
#' @return scatter plot
#' @examples
#' \dontrun{
#' data(matrix)
#' png('test.png')
#' plotmz(matrix)
#' dev.off()
#' }
#' @export
plotmz <- function(data, inscf = 3.5, ...) {
        mz <- as.numeric(rownames(data))
        rt <- as.numeric(colnames(data))
        z <- log10(data + 1)
        cex <- as.numeric(cut(z - inscf, breaks = c(0, 1, 2,
                                                    3, 4, Inf) / 2)) / 2
        cexlab <- c(
                paste0(inscf, "-", inscf + 0.5),
                paste0(inscf +
                               0.5, "-", inscf + 1),
                paste0(inscf + 1, "-", inscf +
                               1.5),
                paste0(inscf + 1.5, "-", inscf + 2),
                paste0(">",
                       inscf + 2)
        )

        z[z < inscf] <- NA
        corr <- which(!is.na(z), arr.ind = TRUE)
        mz0 <- mz[corr[, 1]]
        rt0 <- rt[corr[, 2]]
        int <- z[which(!is.na(z))]

        graphics::par(mar = c(5, 4.2, 6.1, 2.1), xpd = TRUE)
        graphics::plot(
                mz0 ~ rt0,
                pch = 19,
                cex = cex,
                col = grDevices::rgb(0,
                                     0, 1, 0.1),
                xlab = "retention time(s)",
                ylab = "m/z",
                ...
        )
        graphics::legend(
                "top",
                legend = cexlab,
                title = "Intensity in Log scale",
                pt.cex = c(1, 2, 3, 4, 5) / 2,
                pch = 19,
                bty = "n",
                horiz = TRUE,
                cex = 0.7,
                col = grDevices::rgb(0, 0,
                                     1, 0.1),
                inset = c(0,-0.25)
        )
}

#' plot GC-MS data as a heatmap for constant speed of temperature rising
#' @param data imported data matrix of GC-MS
#' @param log transform the intensity into log based 10
#' @param temp temperature range for constant speed
#' @return heatmap
#' @examples
#' \dontrun{
#' plott(matrix)
#' }
#' @export
plott <- function(data,
                  log = FALSE,
                  temp = c(100, 320)) {
        indmz <- as.numeric(rownames(data))
        indrt <- as.numeric(colnames(data))
        col <-
                (grDevices::colorRampPalette(rev(
                        RColorBrewer::brewer.pal(11,
                                                 "RdYlBu")
                )))(100)
        if (log) {
                z <- log10(t(data) + 1)
        } else {
                z <- t(data)
        }
        graphics::par(
                mar = c(2, 5, 1, 4),
                fig = c(0, 1, 0.9,
                        1),
                new = FALSE
        )
        # get the mz and rt range and rotate the matrix to
        zlim <- range(z, na.rm = TRUE)
        breaks <- seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
        poly <- vector(mode = "list", length(col))
        graphics::plot(
                1,
                1,
                t = "n",
                ylim = c(0, 1),
                xlim = zlim,
                xaxt = "n",
                yaxt = "n",
                xaxs = "i",
                yaxs = "i",
                ylab = "",
                xlab = ""
        )
        graphics::mtext(
                "intensity",
                side = 2,
                line = 0.5,
                las = 1,
                cex = 1
        )
        if (log) {
                graphics::axis(
                        1,
                        at = breaks,
                        labels = round(10 ^ (breaks)),
                        las = 1
                )
        } else {
                graphics::axis(
                        1,
                        at = breaks,
                        labels = round(breaks),
                        las = 1
                )
        }

        bks <- seq(zlim[1], zlim[2], length.out = (length(col) +
                                                           1))
        for (i in seq(poly)) {
                graphics::polygon(
                        c(bks[i], bks[i + 1], bks[i +
                                                          1], bks[i]),
                        c(0, 0, 1, 1),
                        col = col[i],
                        border = NA
                )
        }
        # show the heatmap
        graphics::par(
                mar = c(4, 5, 0, 4),
                fig = c(0, 1, 0,
                        0.9),
                new = TRUE
        )
        graphics::image(
                z,
                xlab = expression("Temperature (" *
                                          degree * C * ")"),
                ylab = "m/z",
                axes = FALSE,
                col = col,
                useRaster = TRUE
        )
        # display the temperature as x
        rtx <- seq(0, 1, length.out = length(indrt))
        temp <- round(seq(temp[1], temp[2], length.out = length(indrt)),
                      0)
        graphics::axis(1, at = rtx[temp %% 20 == 0], labels = temp[temp %% 20 ==
                                                                           0])
        # display the m/z as y
        mzy <- seq(0, 1, length.out = length(indmz))
        graphics::axis(2,
                       at = mzy[indmz %% 100 == 0],
                       labels = indmz[indmz %% 100 ==
                                              0],
                       las = 2)
}
#' Plot the background of data
#' @param data imported data matrix of GC-MS
#' @return NULL
#' @examples
#' \dontrun{
#' plotsub(matrix)
#' }
#' @export
plotsub <- function(data) {
        datan <- apply(data, 1, diff)
        datan[datan < 0] <- 0
        datan <- t(datan)
        plotms(datan)
}
#' Plot Total Ion Chromatogram (TIC)
#' @param data imported data matrix of GC-MS
#' @param n logical smooth or not
#' @return plot
#' @examples
#' \dontrun{
#' plottic(matrix)
#' }
#' @export
plottic <- function(data, n = FALSE) {
        data <- apply(data, 2, sum)
        if (n) {
                data <- ma(data, n)
        }
        x <- as.numeric(names(data))
        graphics::plot(
                data ~ x,
                type = "l",
                ylab = "intensity",
                xlab = "retention time(s)",
                frame.plot = FALSE
        )
}
#' plot the information of integration
#' @param list list from getinteagtion
#' @param name the title of the plot
#' @return NULL
#' @examples
#' \dontrun{
#' list <- getinteagtion(rawdata)
#' plotint(list)
#' }
#' @export
plotint <- function(list, name = NULL) {
        area <- list$area
        height <- list$height
        peakdata <- list$peakdata
        RTrange <- list$RTrange
        signal <- list$signal
        slopedata <- list$slopedata
        baseline <- peakdata[1]
        rtstart <- peakdata[2]
        rtend <- peakdata[3]
        rtpeak <- peakdata[4]
        scanstart <- peakdata[5]
        scanend <- peakdata[6]
        scanpeak <- peakdata[7]
        sigstart <- peakdata[8]
        sigend <- peakdata[9]
        sigpeak <- peakdata[10]
        sigpeakbase <- peakdata[11]
        lengthsig <- peakdata[12]
        graphics::plot(
                RTrange,
                signal,
                xlab = "time (min)",
                ylab = "intensity",
                "l",
                ylim = c(-0.02 * max(signal),
                         1.02 * max(signal)),
                main = paste(name, "Peak")
        )
        graphics::lines(c(rtstart, rtend), c(sigstart, sigend),
                        "l", col = "red")
        graphics::lines(c(RTrange[scanstart - baseline + 1],
                          rtstart), c(sigstart, sigstart), "l", col = "darkgreen")
        graphics::lines(c(rtend, RTrange[scanend + baseline -
                                                 1]), c(sigend, sigend), "l", col = "darkgreen")
        graphics::lines(c(rtstart, rtstart),
                        c(0.8 * sigstart,
                          sigstart * 1.2),
                        "l",
                        col = "blue")
        graphics::lines(c(rtend, rtend), c(0.8 * sigend, sigend *
                                                   1.2), "l", col = "blue")
        graphics::lines(c(rtpeak, rtpeak), c(sigpeak, sigpeakbase),
                        "l", col = "blue")

        # print RT, heights & areas
        graphics::text(rtstart, sigpeak * 0.2, as.character(round(rtstart,
                                                                  3)), col = "darkgreen")
        graphics::text(rtend, sigpeak * 0.3, as.character(round(rtend,
                                                                3)), col = "darkgreen")
        graphics::text(rtpeak - 0.1,
                       0.9 * sigpeak,
                       paste("RT:",
                             as.character(round(rtpeak, 3))),
                       col = "darkgreen")
        graphics::text(rtpeak + 0.1, sigpeak * 0.7, paste("area:",
                                                          format(area, digits = 2)), col = "red")
        graphics::text(rtpeak + 0.1,
                       sigpeak * 0.9,
                       paste("height:",
                             format(sigpeak, digits = 2)),
                       col = "red")
}
#' plot the slope information of integration
#' @param list list from getintegration
#' @param name the title of the plot
#' @return NULL
#' @examples
#' \dontrun{
#' list <- getinteragtion(rawdata)
#' plotintslope(list)
#' }
#' @export
plotintslope <- function(list, name = NULL) {
        area <- list$area
        height <- list$height
        peakdata <- list$peakdata
        RTrange <- list$RTrange
        signal <- list$signal
        slopedata <- list$slopedata
        baseline <- peakdata[1]
        rtstart <- peakdata[2]
        rtend <- peakdata[3]
        rtpeak <- peakdata[4]
        scanstart <- peakdata[5]
        scanend <- peakdata[6]
        scanpeak <- peakdata[7]
        sigstart <- peakdata[8]
        sigend <- peakdata[9]
        sigpeak <- peakdata[10]
        sigpeakbase <- peakdata[11]
        lengthsig <- peakdata[12]
        graphics::plot(
                RTrange,
                slopedata,
                xlab = "time (min)",
                ylab = "slope",
                "l",
                main = paste(name, "Slope")
        )
        graphics::lines(c(rtstart, rtstart), c(-0.1 * max(slopedata),
                                               0.1 * max(slopedata)), "l", col = "blue")
        graphics::lines(c(rtend, rtend), c(-0.1 * max(slopedata),
                                           0.1 * max(slopedata)), "l", col = "blue")
        graphics::lines(c(rtpeak, rtpeak), c(-0.5 * max(slopedata),
                                             0.5 * max(slopedata)), "l", col = "blue")
}

#' find line of the regression model for GC-MS
#' @param data imported data matrix of GC-MS
#' @param threshold the threshold of the response (log based 10)
#' @param temp the scale of the oven temperature (constant rate)
#' @return list linear regression model for the matrix
#' @examples
#' \dontrun{
#' data(matrix)
#' findline(matrix)
#' }
#' @export
findline <- function(data,
                     threshold = 2,
                     temp = c(100,
                              320)) {
        y0 <- as.numeric(rownames(data))
        x <- as.numeric(colnames(data))
        # get the group
        group <- ifelse(log10(data) > threshold, 1, 0)
        # get the difference matrix
        diffmatrix <- apply(group, 2, diff)
        # get the points with the smallest differences at the
        # smallest m/z
        difftemp <- apply(diffmatrix, 2, which.min)
        y <- y0[difftemp]
        data <- data.frame(y, x)
        # remove the meaningless bottom
        data <- data[data$y > min(y0) + 1,]
        rangemz <- range(y0)
        rangert <- range(x)
        graphics::plot(
                data$y ~ data$x,
                xlab = expression("Temperature (" *
                                          degree * C * ")"),
                ylab = "m/z",
                pch = 19,
                xlim = rangert,
                ylim = rangemz,
                xaxt = "n",
                yaxt = "n",
                main = "",
                frame.plot = FALSE
        )
        # display the temperature as x
        temp <- round(seq(temp[1], temp[2], length.out = length(x)))
        graphics::axis(1, at = x[temp %% 20 == 0], labels = temp[temp %% 20 ==
                                                                         0])
        # display the m/z as y
        mzy <- seq(min(y0), max(y0), length.out = length(y0))
        graphics::axis(2,
                       at = mzy[y0 %% 100 == 0],
                       labels = y0[y0 %% 100 ==
                                           0],
                       las = 2)
        graphics::abline(stats::lm(data$y ~ data$x),
                         col = "red",
                         lwd = 5)
        graphics::lines(stats::lowess(data$y ~ data$x),
                        col = "blue",
                        lwd = 5)
        slope <- (max(temp) - min(temp)) / (max(data$x) - min(data$x))
        intercept <- min(temp)
        data$x0 <- slope * (data$x - min(data$x)) + intercept
        fit <- stats::lm(data$y ~ data$x0)
        rmse <- round(sqrt(mean(stats::resid(fit) ^ 2)), 2)
        coefs <- stats::coef(fit)
        b0 <- round(coefs[1], 2)
        b1 <- round(coefs[2], 2)
        r2 <- round(summary(fit)$r.squared, 2)
        pv <- stats::anova(fit)$"Pr(>F)"[1]
        eqn <- bquote(italic(m / z) == .(b0) + .(b1) * "*" *
                              italic(Temprature))
        eqn2 <- bquote(r ^ 2 == .(r2) * "," ~ ~ p == .(pv))
        # plot the equation
        posx <- min(x) + (rangert[2] - rangert[1]) * 0.05
        posy <- rangemz[2] - (rangemz[2] - rangemz[1]) * 0.05

        graphics::text(posx, posy, adj = c(0, 0), cex = 1,
                       eqn)
        graphics::text(
                posx,
                posy - (rangemz[2] - rangemz[1]) *
                        0.05,
                adj = c(0, 0),
                cex = 1,
                eqn2
        )
        graphics::legend(
                "topright",
                c("OLS", "lowess"),
                box.lty = 0,
                pch = c(-1,-1),
                lty = c(1, 1),
                lwd = c(2, 2),
                col = c("red", "blue")
        )
        return(fit)
}

#' Plot the response group of GC-MS
#' @param data imported data matrix of GC-MS
#' @param threshold the threshold of the response (log based 10) to separate the group
#' @return list linear regression model for the data matrix
#' @examples
#' \dontrun{
#' data(matrix)
#' plotgroup(matrix)
#' }
#' @export
plotgroup <- function(data, threshold = 2) {
        group <- ifelse(log10(data) > threshold, 1, 0)
        ind <- as.numeric(rownames(data))
        m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2)
        graphics::layout(m)
        graphics::hist(
                log10(data),
                breaks = 100,
                main = "",
                xlab = "Intensity"
        )
        graphics::abline(
                v = threshold,
                lwd = 5,
                lty = 2,
                col = "red"
        )
        graphics::legend(
                "topright",
                "threshold",
                box.lty = 0,
                pch = -1,
                lty = 2,
                lwd = 2,
                col = "red"
        )

        graphics::image(
                t(group),
                xlab = "retention time(min)",
                ylab = "m/z",
                axes = FALSE,
                col = grDevices::heat.colors(2),
                useRaster = TRUE
        )
        indmz <- as.numeric(rownames(data))
        indrt <- as.numeric(colnames(data))
        # display the RT as x
        rtx <- seq(0, 1, length.out = length(indrt))
        graphics::axis(1,
                       at = c(0, rtx[indrt %% 300 == 0], 1),
                       labels = c("", indrt[indrt %% 300 == 0], ""))
        # display the m/z as y
        mzy <- seq(0, 1, length.out = length(indmz))
        graphics::axis(2,
                       at = mzy[indmz %% 100 == 0],
                       labels = indmz[indmz %% 100 ==
                                              0],
                       las = 2)
        findline(data)
}

#' Plot the intensity distribution of GC-MS
#' @param meanmatrix mean data matrix of GC-MS(n=5)
#' @param rsdmatrix standard deviation matrix of GC-MS(n=5)
#' @return NULL
#' @examples
#' \dontrun{
#' plotsms(meanmatrix,rsdmatrix)
#' }
#' @export
plotsms <- function(meanmatrix, rsdmatrix) {
        graphics::par(
                mar = c(4.2, 4.2, 0, 1.5),
                fig = c(0,
                        1, 0, 0.8),
                new = FALSE,
                cex.axis = 1.5,
                cex.lab = 1.5
        )
        graphics::smoothScatter(
                y = rsdmatrix * 100,
                x = log10(c(meanmatrix)),
                main = "",
                xlab = "Intensity",
                ylab = "Relative Standard Deviation(%)",
                xaxt = "n",
                frame.plot = FALSE
        )
        graphics::abline(h = 20,
                         lty = 2,
                         col = "red")
        graphics::abline(h = 10, col = "red")
        graphics::axis(
                1,
                at = c(0, 1, 2, 3, 4, 5, 6, 7),
                labels = c(
                        "1",
                        "10",
                        "100",
                        "1000",
                        "10000",
                        "100000",
                        "1000000",
                        "100000000"
                )
        )
        graphics::par(
                mar = c(0, 4.2, 1, 1.5),
                oma = c(0, 0,
                        0, 0),
                fig = c(0, 1, 0.8, 1),
                new = TRUE,
                cex.axis = 1.5,
                cex.lab = 1.5
        )
        graphics::hist(
                log10(meanmatrix),
                breaks = 100,
                xlab = "Intensity",
                main = "",
                xaxt = "n"
        )
        graphics::axis(
                1,
                at = c(0, 1, 2, 3, 4, 5, 6, 7),
                labels = c(
                        "1",
                        "10",
                        "100",
                        "1000",
                        "10000",
                        "100000",
                        "1000000",
                        "100000000"
                )
        )
}

#' plot the density of the GC-MS data with EM algorithm to separate the data into two log normal distribution.
#' @param data imported data matrix of GC-MS
#' @return NULL
#' @examples
#' \dontrun{
#' # generate a matrix from raw data with row as m/z and column as retention time
#' plothist(matrix)
#' }
#' @export
plothist <- function(data) {
        data1 <- sample(data, 1e+05)
        mixmdl <- mixtools::normalmixEM(log10(data1))
        graphics::plot(mixmdl,
                       which = 2,
                       breaks = 100,
                       xlab2 = "Intensity")
        graphics::lines(stats::density(log10(data1)),
                        lty = 2,
                        lwd = 2)
        graphics::legend(
                "topright",
                c("noise", "signal", "density"),
                box.lty = 0,
                pch = c(-1,-1,-1),
                lty = c(1, 1,
                        2),
                lwd = c(2, 2, 2),
                col = c("red", "green",
                        "black")
        )
}
#' plot the calibration curve with error bar, r squared and equation.
#' @param x concentration
#' @param y response
#' @param upper upper error bar
#' @param lower lower error bar
#' @param ... parameters for `plot` function
#' @return NULL
#' @examples
#' \dontrun{
#' plotcc(x,y,upper)
#' }
#' @export
plotcc <- function(x, y, upper, lower = upper, ...) {
        graphics::plot(x, y, ...)
        graphics::arrows(
                x,
                y + upper,
                x,
                y - lower,
                angle = 90,
                code = 3,
                length = 0.1
        )
        fit <- stats::lm(y ~ x)
        graphics::abline(fit)
        m1 <-  summary(fit)
        graphics::mtext(paste0("R squared: ", round(m1$r.squared, 4)), adj = 0)
        graphics::mtext(paste0(
                "y ~ ",
                round(m1$coefficients[2], 2),
                "x",
                " + ",
                round(m1$coefficients[1], 2)
        ),
        adj = 1)
}
#' Just integrate data according to fixed rt and fixed noise area
#' @param data file should be a dataframe with the first column RT and second column intensity of the SIM ions.
#' @param rt a rough RT range contained only one peak to get the area
#' @param brt a rough RT range contained only one peak and enough noises to get the area
#' @param smoothit logical, if using an average smooth box or not. If using, n will be used
#' @return area integration data
#' @examples
#' \dontrun{
#' area <- Integration(data)
#' }
#' @export
Integration <- function(data,
                        rt = c(8.3, 9),
                        brt = c(8.3,
                                8.4),
                        smoothit = TRUE) {
        # subset the data
        subdata <- data[data[, 1] > rt[2] & data[, 1] < rt[1],]
        # get the signal and the RT
        RTrange <- subdata[, 1]
        signal <- subdata[, 2]
        # subset the noise
        subnoise <- data[data[, 1] > brt[2] & data[, 1] < brt[1],]
        # get the noise and the RT
        RTrange2 <- subnoise[, 1]
        noise <- subnoise[, 2]
        # data smooth
        if (smoothit) {
                signal <- ma(signal, n = 5)
                noise <- ma(noise, n = 5)
        }
        baseline <- mean(noise)
        signaln <- signal - baseline
        area <- 0
        # calculate area; using a Riemann integral (dimension:
        # intensity x delta time)
        for (i in 1:(length(RTrange) - 1)) {
                area <- area + signaln[i] * (RTrange[i + 1] - RTrange[i]) *
                        60
        }
        return(area)
}

#' GetIntegration was mainly used for get the integration of certain ion's chromatogram data and plot the data
#' @param data file should be a dataframe with the first column RT and second column intensity of the SIM ions.
#' @param rt a rough RT range contained only one peak to get the area
#' @param n points in the moving average smooth box, default value is 5
#' @param m numbers of points for regression to get the slope
#' @param slope the threshold value for start/stop peak as percentage of max slope
#' @param baseline numbers of the points for the baseline of the signal
#' @param noslope logical, if using a horizon line to get area or not
#' @param smoothit logical, if using an average smooth box or not. If using, n will be used
#' @param half logical, if using the left half peak to calculate the area
#' @return integration data such as peak area, peak height, signal and the slope data.
#' @examples
#' \dontrun{
#' list <- GetIntegration(data)
#' }
#' @export
GetIntegration <- function(data,
                           rt = c(8.3, 9),
                           n = 5,
                           m = 5,
                           slope = c(2, 2),
                           baseline = 10,
                           noslope = TRUE,
                           smoothit = TRUE,
                           half = FALSE) {
        # subset the data
        subdata <- data[data[, 1] > rt[1] & data[, 1] < rt[2],]
        # get the signal and the RT
        RTrange <- subdata[, 1]
        signal <- subdata[, 2]
        # data smooth
        if (smoothit) {
                signal <- ma(signal, n)
        }
        # get the slope data
        RTrangemsec <- RTrange * 60 * 1000
        slopedata <- signal
        back <- round(m / 2)
        forth <- m - back
        if (m > 2) {
                for (i in (back + 1):(length(signal) - forth)) {
                        slopedata[i] <- stats::coef(stats::lm(signal[(i -
                                                                              back + 1):(i + forth)] ~ RTrangemsec[(i -
                                                                                                                            back + 1):(i + forth)]))[2]
                }
                slopedata[1:back] <-
                        slopedata[back + 1]  # first few points
                slopedata[(length(signal) - forth - 1):length(signal)] <-
                        slopedata[(length(signal) -
                                           forth)]  # last few points
        } else {
                # for n = 2 points; without linear regression (much
                # faster) time per scan in millisec
                delta_t <-
                        (t[length(RTrangemsec)] - RTrangemsec[1] / (length(RTrangemsec) -
                                                                            1))
                for (i in 2:length(signal)) {
                        slopedata[i] <- (signal[i] - signal[i - 1]) / delta_t
                }
                slopedata[1] <- 0
        }
        # search for peak start
        i <- baseline
        while ((slopedata[i] <= (slope[1] / 100 * max(slopedata))) &
               (i < length(signal)))
                i <- i + 1
        rtstart <- RTrange[i]  # peak start found
        scanstart <- i  # (slope > threshold)
        sigstart <-
                mean(signal[(i - baseline + 1):i])  # baseline intensity found
        # search for peak top
        i <- which.max(slopedata)  # jump to slope max.
        while ((slopedata[i] >= 0) & (i < (length(signal) -
                                           baseline)))
                i <- i + 1
        rtpeak <- RTrange[i]  # peak top found
        scanpeak <- i  # (slope = 0)
        sigpeak <- signal[i]
        # search for peak end
        i <- which.min(slopedata)  # jump to slope min.
        while ((slopedata[i] <= -(slope[2] / 100 * max(slopedata))) &
               (i < (length(signal) - baseline)))
                i <- i + 1
        rtend <- RTrange[i]  # peak end found
        scanend <- i  # (-slope < threshold)
        sigend <- mean(signal[i:(i + baseline - 1)])
        # if background without slope
        if (!noslope)
                sigend <- sigstart
        # subtract background from signal
        background <- signal
        for (i in scanstart:scanend) {
                # get background
                background[i] <-
                        sigstart + (sigend - sigstart) / (scanend -
                                                                  scanstart) * (i - scanstart)
        }
        subsignal <- signal - background
        # get the length of the signal
        lengthsig <- length(scanstart:scanend)
        # calculate area; using a Riemann integral (dimension:
        # intensity x min)
        area <- 0
        scantime <-
                (RTrange[scanend] - RTrange[scanstart]) / (scanend -
                                                                   scanstart) * 60  # time per scan in second
        # when half peak
        if (half == TRUE) {
                for (i in scanstart:scanpeak)
                        area <- area + subsignal[i] *
                                scantime
        } else {
                for (i in scanstart:scanend)
                        area <- area + subsignal[i] *
                                scantime
        }
        bgstart <- RTrange[scanstart - baseline + 1]
        bgend <- RTrange[scanend + baseline - 1]
        # calculate height
        sigpeakbase <- sigstart + (sigend - sigstart) / (scanend -
                                                                 scanstart) * (scanpeak - scanstart)
        height <- sigpeak - sigpeakbase
        # SNR
        snrnoise <- abs(diff(range(signal[(scanstart - baseline +
                                                   1):scanstart])))
        SNR <- height / snrnoise
        # collect the data for plot peak and slope
        peakdata <- c(
                baseline,
                rtstart,
                rtend,
                rtpeak,
                scanstart,
                scanend,
                scanpeak,
                sigstart,
                sigend,
                sigpeak,
                sigpeakbase,
                lengthsig,
                SNR
        )
        names(peakdata) <- c(
                "baseline",
                "peak start RT",
                "peak end RT",
                "peak RT",
                "baseline start RT ID",
                "baseline end RT ID",
                "baseline peak RT ID",
                "start RT intensity",
                "end RT intensity",
                "peak RT intensity",
                "peak baseline",
                "points",
                "SNR"
        )
        # return the result as a list
        list <- list(
                area = area,
                height = height,
                peakdata = peakdata,
                RTrange = RTrange,
                signal = signal,
                slopedata = slopedata
        )
        return(list)
}

#' Get the selected isotopologues at certain MS data
#' @param formula the molecular formula.
#' @param charge the charge of that molecular. 1 in EI mode as default
#' @param width the width of the peak width on mass spectrum. 0.3 as default for low resolution mass spectrum.
#' @examples
#' \dontrun{
#' # show isotopologues
#' Getisotopologues(formula = 'C6H11O6', charge = 1, width = 0.3)
#' }
#' @export
Getisotopologues <- function(formula = "C6H11O6",
                             charge = 1,
                             width = 0.3) {
        # input the formula and charge for your molecular,
        formula <-  Rdisop::getMolecule(formula, z = charge, maxisotopes = 20)
        # get the isotopes pattern of your molecular with high
        # abundances. Here we suggest more than 10% abundance
        # of your base peak would meet the SNR
        isotopes <- data.frame(t(formula$isotopes[[1]]))
        isotopes <- isotopes[isotopes[, 2] > 0.1,]
        # order the intensity by the abundance
        findpairs <-
                isotopes[order(isotopes[, 2], decreasing = TRUE),]
        # find the most similar pairs with high abundance
        df <- outer(findpairs[, 1], findpairs[, 1], "/")
        rownames(df) <- colnames(df) <- findpairs[, 1]
        diag(df) <- df[upper.tri(df)] <- 0
        t <- which(df == max(df), arr.ind = TRUE)
        isotopologues1 <- as.numeric(rownames(df)[t[1]])
        isotopologues2 <- as.numeric(colnames(df)[t[2]])
        isotopologuesL <- min(isotopologues1, isotopologues2)
        isotopologuesH <- max(isotopologues1, isotopologues2)
        # get the caculated ratio at certain resolution
        isotopes2 <-
                data.frame(t(formula$isotopes[[1]]))
        ratio <- sum(isotopes2[isotopes2[, 1] > isotopologuesL -
                                       width &
                                       isotopes2[, 1] < isotopologuesL + width,
                               2]) / sum(isotopes2[isotopes2[, 1] > isotopologuesH -
                                                           width &
                                                           isotopes2[, 1] < isotopologuesH + width,
                                                   2])
        peak <-
                c(
                        round(isotopologuesL, digits = 5),
                        round(isotopologuesH,
                              digits = 5),
                        round(ratio, digits = 5)
                )
        # peak <- as.character(peak)
        names(peak) <- c("light isotopologue",
                         "high isotopologue",
                         "caculated ratio")
        return(data.frame(peak))
}
