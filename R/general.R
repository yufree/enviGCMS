#' Demo data
#' @docType data
#' @usage data(list)
#' @format A list object with data, mass to charge ratio, retention time and group information. The list is generated from faahKO package by `getmr` function.
"list"

#' define the Mode function
#' @param x vector
#' @return Mode of the vector
#' @export
Mode = function(x) {
        ta = table(x)
        tam = max(ta)
        if (all(ta == tam))
                mod = x
        else if (is.numeric(x))
                mod = as.numeric(names(ta)[ta == tam])
        else
                mod = names(ta)[ta == tam]
        return(mod)
}
#' filter data by average moving box
#'
#' @param x a vector
#' @param n A number to indentify the size of the moving box.
#' @return The filtered data
#' @examples
#' ma(rnorm(1000),5)
#' @export
ma <- function(x, n) {
        stats::filter(x, rep(1 / n, n), circular = T)
}

#' Import data and return the annotated matrix for GC/LC-MS by m/z range and retention time
#' @param data file type which xcmsRaw could handle
#' @param mzstep the m/z step for generating matrix data from raw mass spectral data
#' @param mzrange vector range of the m/z, default all
#' @param rtrange vector range of the retention time, default all
#' @return matrix with the row as increasing m/z second and column as increasing scantime
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' matrix <- getmd(cdffiles[1])
#' }
getmd <- function(data,
                  mzstep = 0.1,
                  mzrange = F,
                  rtrange = F) {
        data <- xcms::xcmsRaw(data, profstep = mzstep)
        pf <- xcms::profMat(data)
        rownames(pf) <- mz <- xcms::profMz(data)
        colnames(pf) <- rt <- data@scantime
        if (mzrange[1]) {
                pf <- pf[mz > mzrange[1] & mz < mzrange[2],]
        }
        if (rtrange[1]) {
                pf <- pf[, rt > rtrange[1] & rt < rtrange[2]]
        }
        return(pf)
}

#' Combine two data with similar retention time while different mass range
#'
#' @param data1 data file path of lower mass range
#' @param data2 data file path of higher mass range
#' @param mzstep the m/z step for generating matrix data from raw mass spectral data
#' @param rtstep the alignment accuracy of retention time, e.g. 0.01 means the retention times of combined data should be the same at the accuracy 0.01s. Higher rtstep would return less scans for combined data
#' @return matrix with the row as scantime in second and column as m/z
#' @examples
#' \dontrun{
#' # mz100_200 and mz201_300 were the path to the raw data
#' matrix <- getmd(mz100_200,mz201_300)
#' }
#' @export
cbmd <- function(data1,
                 data2,
                 mzstep = 0.1,
                 rtstep = 0.01) {
        z1 <- getmd(data1, mzstep = mzstep)
        z2 <- getmd(data2, mzstep = mzstep)
        roundrt <- -log10(rtstep)
        z <- rbind(z1[, which(round(as.numeric(colnames(z1)),
                                    roundrt) %in% round(as.numeric(colnames(z2)), roundrt))],
                   z2[, which(round(as.numeric(colnames(z2)), roundrt) %in%
                                      round(as.numeric(colnames(z1)), roundrt))])
        rownames(z) <- c(rownames(z1), rownames(z2))
        colnames(z) <- intersect(round(as.numeric(colnames(z1)),
                                       roundrt),
                                 round(as.numeric(colnames(z2)), roundrt))
        return(z)
}

#' Get the differences of two GC/LC-MS data
#'
#' @param data1 data file path of first data
#' @param data2 data file path of second data
#' @param mzstep the m/z step for generating matrix data from raw mass spectral data
#' @param rtstep the alignment accuracy of retention time, e.g. 0.01 means the retention times of combined data should be the same at the accuracy 0.01s. Higher rtstep would return less scans for combined data
#' @return list four matrix with the row as scantime in second and column as m/z, the first matrix refer to data 1, the second matrix refer to data 2, the third matrix refer to data1 - data2 while the fourth refer to data2 - data1, minus values are imputed by 0
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' matrix <- submd(cdffiles[1],cdffiles[7])
#' }
#' @export
submd <- function(data1,
                  data2,
                  mzstep = 0.1,
                  rtstep = 0.01) {
        z1 <- getmd(data1, mzstep = mzstep)
        z2 <- getmd(data2, mzstep = mzstep)
        roundmz <- -log10(mzstep)
        roundrt <- -log10(rtstep)

        z10 <- z1[which(round(as.numeric(rownames(z1)), roundmz) %in%
                                round(as.numeric(rownames(z2)), roundmz)), which(round(as.numeric(colnames(z1)),
                                                                                       roundrt) %in% round(as.numeric(colnames(z2)), roundrt))]
        z20 <- z2[which(round(as.numeric(rownames(z2)), roundmz) %in%
                                round(as.numeric(rownames(z1)), roundmz)), which(round(as.numeric(colnames(z2)),
                                                                                       roundrt) %in% round(as.numeric(colnames(z1)), roundrt))]

        z <- z10 - z20
        z[z < 0] <- 0
        z0 <- z20 - z10
        z0[z0 < 0] <- 0

        rownames(z0) <-
                rownames(z) <- intersect(round(as.numeric(rownames(z1)),
                                               roundmz),
                                         round(as.numeric(rownames(z2)), roundmz))
        colnames(z0) <-
                colnames(z) <- intersect(round(as.numeric(colnames(z1)),
                                               roundrt),
                                         round(as.numeric(colnames(z2)), roundrt))
        xlim <- range(as.numeric(colnames(z)))
        ylim <- range(as.numeric(rownames(z)))

        graphics::par(mfrow = c(2, 2))

        plotmz(z10,
               main = "data 1",
               xlim = xlim,
               ylim = ylim)
        plotmz(z20,
               main = "data 2",
               xlim = xlim,
               ylim = ylim)
        plotmz(z,
               main = "data1 - data2",
               xlim = xlim,
               ylim = ylim)
        plotmz(z0,
               main = "data2 - data1",
               xlim = xlim,
               ylim = ylim)

        li <- list(
                data1 = z10,
                data2 = z20,
                data1s2 = z,
                data2s1 = z0
        )
        return(li)
}


#' plot GC/LC-MS data as a heatmap with TIC
#'
#' @param data imported data matrix of GC-MS
#' @param log transform the intensity into log based 10
#' @return heatmap
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' matrix <- getmd(cdffiles[1])
#' png('test.png')
#' plotms(matrix)
#' dev.off()
#' }
#' @export
plotms <- function(data, log = F) {
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
        zlim <- range(z, na.rm = T)
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
                new = T
        )
        graphics::image(
                z,
                ylab = "",
                axes = F,
                col = col,
                useRaster = T
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
                new = T,
                cex.lab = 1
        )
        data[is.na(data)] <- 0
        data <- apply(data, 2, sum)
        graphics::plot(
                data ~ indrt,
                type = "l",
                ylab = "intensity",
                xlab = "retention time(s)",
                frame.plot = F,
                xaxs = "i",
                yaxs = "i"
        )
}

#' plot GC/LC-MS data as scatter plot
#'
#' @param data imported data matrix of GC-MS
#' @param inscf Log intensity cutoff for peaks, default 5
#' @param ... parameters for `plot` function
#' @return scatter plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' matrix <- getmd(cdffiles[1])
#' png('test.png')
#' plotmz(matrix)
#' dev.off()
#' }
#' @export
plotmz <- function(data, inscf = 5, ...) {
        mz <- as.numeric(rownames(data))
        rt <- as.numeric(colnames(data))
        z <- log10(data + 1)
        cex = as.numeric(cut(z - inscf, breaks = c(0, 1, 2,
                                                   3, 4, Inf) / 2)) / 2
        cexlab = c(
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
                horiz = T,
                cex = 0.7,
                col = grDevices::rgb(0, 0,
                                     1, 0.1),
                inset = c(0,-0.25)
        )
}

#' plot GC-MS data as a heatmap for constant speed of temperature rising
#' @param data imported data matrix of GC-MS
#' @param log transform the intensity into log based 10
#' @param temp temprature range for constant speed
#' @return heatmap
#' @examples
#' \dontrun{
#' matrix <- getmd(rawdata)
#' plott(matrix)
#' }
#' @export
plott <- function(data,
                  log = F,
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
                new = F
        )
        # get the mz and rt range and rotate the matrix to
        zlim <- range(z, na.rm = T)
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
                new = T
        )
        graphics::image(
                z,
                xlab = expression("Temperature (" *
                                          degree * C * ")"),
                ylab = "m/z",
                axes = F,
                col = col,
                useRaster = T
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

#' Plot the backgrond of data
#' @param data imported data matrix of GC-MS
#' @return NULL
#' @examples
#' \dontrun{
#' matrix <- getmd(rawdata)
#' plotsub(matrix)
#' }
#' @export
plotsub <- function(data) {
        datan <- apply(data, 1, diff)
        datan[datan < 0] <- 0
        datan <- t(datan)
        plotms(datan)
}

#' Plot mass spectrum of certain retention time and return mass spectrum vector (MSP file) for NIST search
#' @param data imported data matrix of GC-MS
#' @param rt vector range of the retention time
#' @param ms vector range of the m/z
#' @return plot, vector and MSP files for NIST search
#' @examples
#' \dontrun{
#' matrix <- getmd(rawdata)
#' plotrtms(matrix,rt = c(500,1000),ms = (300,500))
#' }
#' @export
plotrtms <- function(data, rt, ms) {
        data <- getmd(data, rt, ms)
        temp <- apply(data, 1, mean)
        graphics::plot(
                temp,
                xaxs = "i",
                yaxs = "i",
                type = "h",
                ylab = "intensity",
                xlab = "m/z",
                xaxt = "n",
                frame.plot = F
        )
        mz <- as.numeric(names(temp))
        index <- seq(1, length(temp))
        graphics::axis(1, at = index[mz %% 100 == 0], labels = mz[mz %% 100 ==
                                                                          0])
        writeMSP(temp)
        return(temp)
}

#' Plot EIC of certain m/z and return dataframe for intergration
#' @param data imported data matrix of GC-MS
#' @param ms m/z to be extracted
#' @param rt vector range of the retention time
#' @param n logical smooth or not
#' @return dataframe with  with the first column RT and second column intensity of the SIM ions.
#' @examples
#' \dontrun{
#' matrix <- getmd(rawdata)
#' plotmsrt(matrix,rt = c(500,1000),ms = 300)
#' }
#' @export
plotmsrt <- function(data, ms, rt, n = F) {
        data <- getmd(data, rt, c(ms, ms + 1))[1,]
        if (n) {
                data <- ma(data, n)
        }
        x <- as.numeric(names(data))
        graphics::plot(
                data ~ x,
                type = "l",
                xaxs = "i",
                yaxs = "i",
                main = bquote("m/z = " ~ .(ms)),
                ylab = "intensity",
                xlab = "retention time(s)",
                frame.plot = F
        )
        data <- data.frame(x, data)
        colnames(data) <- c("RT", "Intensity")
        return(data)
}

#' Plot Total Ion Chromatogram (TIC)
#' @param data imported data matrix of GC-MS
#' @param n logical smooth or not
#' @return plot
#' @examples
#' \dontrun{
#' matrix <- getmd(rawdata)
#' plottic(matrix)
#' }
#' @export
plottic <- function(data, n = F) {
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
                frame.plot = F
        )
}
#' plot the information of intergretion
#' @param list list from getinteragtion
#' @param name the title of the plot
#' @return NULL
#' @examples
#' \dontrun{
#' list <- getinteragtion(rawdata)
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
#' plot the slope information of intergretion
#' @param list list from getinteragtion
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
#' @param temp the scale of the oven temprature(constant rate)
#' @return list linear regression model for the matrix
#' @examples
#' \dontrun{
#' data <- getmd(rawdata)
#' findline(data)
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
                frame.plot = F
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
#' @param threshold the threshold of the response (log based 10) to seperate the group
#' @return list linear regression model for the data matrix
#' @examples
#' \dontrun{
#' data <- getmd(rawdata)
#' plotgroup(data)
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
                axes = F,
                col = grDevices::heat.colors(2),
                useRaster = T
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
#' data1 <- getmd(‘sample1-1’)
#' data2 <- getmd(‘sample1-2’)
#' data3 <- getmd(‘sample1-3’)
#' data4 <- getmd(‘sample1-4’)
#' data5 <- getmd(‘sample1-5’)
#' data <- (data1+data2+data3+data4+data5)/5
#' datasd <- sqrt(((data1-data)^2+(data2-data)^2+(data3-data)^2+(data4-data)^2+(data5-data)^2)/4)
#' databrsd <- datasd/data
#' plotsms(meanmatrix,rsdmatrix)
#' }
#' @export
plotsms <- function(meanmatrix, rsdmatrix) {
        graphics::par(
                mar = c(4.2, 4.2, 0, 1.5),
                fig = c(0,
                        1, 0, 0.8),
                new = F,
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
                frame.plot = F
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
                new = T,
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

#' plot the density of the GC-MS data with EM algorithm to seperate the data into two log normal distribution.
#' @param data imported data matrix of GC-MS
#' @return NULL
#' @examples
#' \dontrun{
#' matrix <- getmd(rawdata)
#' plothist(matrix)
#' }
#' @export
plothist <- function(data) {
        data1 = sample(data, 1e+05)
        mixmdl = mixtools::normalmixEM(log10(data1))
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
