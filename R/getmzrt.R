#' Check class
#'
#' @noRd
is.mzrt <- function(x)
        inherits(x, "mzrt")
#' c method for mzrt object
#'
#' @noRd
c.mzrt <- function(x, ...) {
        li <- list(x, ...)
        re <- list()
        re$mz <- do.call(c, lapply(li, function(x)
                x$mz))
        re$rt <- do.call(c, lapply(li, function(x)
                x$rt))
        re$data <- do.call(rbind, lapply(li, function(x)
                x$data))
        re$mzrange <-
                do.call(rbind, lapply(li, function(x)
                        x$mzrange))
        re$rtrange <-
                do.call(rbind, lapply(li, function(x)
                        x$rtrange))
        re$groupmean <-
                do.call(rbind, lapply(li, function(x)
                        x$groupmean))
        re$groupsd <-
                do.call(rbind, lapply(li, function(x)
                        x$groupsd))
        re$grouprsd <-
                do.call(rbind, lapply(li, function(x)
                        x$grouprsd))
        re$group <- x$group
        class(re) <- "mzrt"
        return(re)
}
#' Convert an list object to csv file.
#' @param list list with data as peaks list, mz, rt and group information
#' @param name result name for csv and/or eic file, default NULL
#' @param mzdigit m/z digits of row names of data frame, default 4
#' @param rtdigit retention time digits of row names of data frame, default 1
#' @param type csv format for further analysis, m means  Metaboanalyst, a means xMSannotator, p means Mummichog(NA values are imputed by `getimputation`, and F test is used here to generate stats and p value), o means full information csv (for `pmd` package), default o. mapo could output all those format files.
#' @param target logical, preserve original rowname of data or not for target data, default FALSE.
#' @param ... other parameters for `write.table`
#' @return NULL, csv file
#' @references Li, S.; Park, Y.; Duraisingham, S.; Strobel, F. H.; Khan, N.; Soltow, Q. A.; Jones, D. P.; Pulendran, B. PLOS Computational Biology 2013, 9 (7), e1003123.
#' Xia, J., Sinelnikov, I.V., Han, B., Wishart, D.S., 2015. MetaboAnalyst 3.0—making metabolomics more meaningful. Nucl. Acids Res. 43, W251–W257.
#' @examples
#' \dontrun{
#' data(list)
#' getcsv(list,name='demo')
#' }
#' @export
getcsv <-
        function(list,
                 name,
                 mzdigit = 4,
                 rtdigit = 1,
                 type = 'o',
                 target = FALSE,
                 ...) {
                if (!is.null(name)) {
                        if (grepl('m', type)) {
                                sample_group <- list$group[,-1]
                                data <-
                                        rbind(sample_group, list$data)
                                rownames(data) <-
                                        c("group",
                                          paste0(
                                                  "M",
                                                  round(list$mz, mzdigit),
                                                  "T",
                                                  round(list$rt, mzdigit)
                                          ))
                                filename <-
                                        paste0(name, "metaboanalyst.csv")
                                utils::write.csv(data, file = filename, ...)
                        }
                        if (grepl('a', type)) {
                                mz <- list$mz
                                time <- list$rt
                                data <-
                                        as.data.frame(cbind(mz, time, list$data))
                                rownames(data) <-
                                        paste0("M",
                                               round(mz, mzdigit),
                                               "T",
                                               round(time, rtdigit))
                                data <- unique(data)
                                filename <-
                                        paste0(name, "xMSannotator.csv")
                                utils::write.csv(data, file = filename, ...)
                        }
                        if (grepl('p', type)) {
                                lv <- list$group[,-1]
                                lif <-
                                        getimputation(list, method = "l")
                                ar <-
                                        apply(lif$data, 1, function(x)
                                                stats::anova(stats::lm(x ~ lv)))

                                df <-
                                        cbind.data.frame(
                                                m.z = list$mz,
                                                rt = list$rt,
                                                p.value = vapply(ar, function(x)
                                                        x$`Pr(>F)`[1], 1),
                                                t.score = vapply(ar, function(x)
                                                        x$`F value`[1], 1)
                                        )
                                filename <-
                                        paste0(name, 'mummichog.txt')
                                utils::write.table(
                                        df,
                                        file = filename,
                                        sep = "\t",
                                        row.names = FALSE,
                                        ...
                                )
                        }
                        if (grepl('o', type)) {
                                data <- cbind(mz = list$mz,
                                              rt = list$rt,
                                              list$data)
                                colname <- colnames(data)
                                groupz <- list$group[,-1]
                                if (ncol(list$group) > 2) {
                                        cols <- colnames(groupz)
                                        groupz <-
                                                do.call(paste, c(groupz[cols], sep = "_"))
                                }
                                groupt <- c('mz', 'rt', groupz)
                                data <- rbind(groupt, data)
                                if (!target) {
                                        rownames(data) <-
                                                c('group',
                                                  paste0(
                                                          "M",
                                                          round(list$mz, mzdigit),
                                                          "T",
                                                          round(list$rt, rtdigit)
                                                  ))
                                }
                                colnames(data) <- colname
                                filename <- paste0(name, "mzrt.csv")
                                utils::write.csv(data, file = filename, ...)
                        }
                }
        }
#' Get a mzrt list and/or save mz and rt range as csv file.
#' @param list list with data as peaks list, mz, rt and group information
#' @param name result name for csv and/or eic file, default NULL
#' @param ... other parameters for `write.table`
#' @return NULL, csv file
#' @export
getrangecsv <- function(list, name, ...) {
        df <- cbind.data.frame(list$mzrange, list$rtrange)
        filename <- paste0(name, "mzrtrange.csv")
        utils::write.csv(df, file = filename, ...)
}

#' Impute the peaks list data
#' @param list list with data as peaks list, mz, rt and group information
#' @param method 'r' means remove, 'l' means use half the minimum of the values across the peaks list, 'mean' means mean of the values across the samples, 'median' means median of the values across the samples, '0' means 0, '1' means 1. Default 'l'.
#' @return list with imputed peaks
#' @examples
#' data(list)
#' getimputation(list)
#' @export
#' @seealso \code{\link{getdoe}}
getimputation <- function(list, method = "l") {
        data <- list$data
        mz <- list$mz
        rt <- list$rt

        if (method == "r") {
                idx <- stats::complete.cases(data)
                data <- data[idx,]
                mz <- mz[idx]
                rt <- rt[idx]
        } else if (method == "l") {
                impute <- min(data, na.rm = TRUE) / 2
                data[is.na(data)] <- impute
        } else if (method == "mean") {
                for (i in seq_len(ncol(data))) {
                        data[is.na(data[, i]), i] <- mean(data[, i],
                                                          na.rm = TRUE)
                }
        } else if (method == "median") {
                for (i in seq_len(ncol(data))) {
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
#' Filter the data based on row and column index
#' @param list list with data as peaks list, mz, rt and group information
#' @param rowindex logical, row index to keep
#' @param colindex logical, column index to keep
#' @param name file name for csv and/or eic file, default NULL
#' @param type csv format for further analysis, m means  Metaboanalyst, a means xMSannotator, p means Mummichog(NA values are imputed by `getimputation`, and F test is used here to generate stats and p value), o means full information csv (for `pmd` package), default o. mapo could output all those format files.
#' @param ... other parameters for `getcsv`
#' @return list with remain peaks, and filtered peaks index
#' @examples
#' data(list)
#' li <- getdoe(list)
#' lif <- getfilter(li,rowindex = li$rsdindex)
#' @export
#' @seealso \code{\link{getimputation}}, \code{\link{getcsv}}
getfilter <-
        function(list,
                 rowindex = TRUE,
                 colindex = TRUE,
                 name = NULL,
                 type = 'o',
                 ...) {
                if (!is.null(rowindex) & !is.null(list$rowindex)) {
                        rowindex <- rowindex & list$rowindex
                }
                else if (is.null(rowindex) &
                         !is.null(list$rowindex)) {
                        rowindex <- list$rowindex
                }
                list$data <- list$data[rowindex, ]
                list$mz <- list$mz[rowindex]
                list$rt <- list$rt[rowindex]
                list$mzrange <- list$mzrange[rowindex, ]
                list$rtrange <- list$rtrange[rowindex, ]
                list$groupmean <- list$groupmean[rowindex, ]
                list$groupsd <- list$groupsd[rowindex, ]
                list$grouprsd <- list$grouprsd[rowindex, ]
                # list$rowindex <- rowindex

                if (!is.null(colindex) & !is.null(list$colindex)) {
                        colindex <- colindex & list$colindex
                }
                else if (is.null(colindex) &
                         !is.null(list$colindex)) {
                        colindex <- list$colindex
                }
                list$data <- list$data[, colindex]
                if (NCOL(list$group) > 1) {
                        list$group <- list$group[colindex,]
                } else{
                        list$group <- list$group[colindex]
                }
                # list$colindex <- colindex
                getcsv(list, name = name, type = type, ...)
                class(list) <- 'mzrt'
                return(list)
        }
#' Generate the group level rsd and average intensity based on DoE,
#' @param list list with data as peaks list, mz, rt and group information
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param tr logical. TRUE means dataset with technical replicates at the base level folder
#' @param rsdcft the rsd cutoff of all peaks in technical replicates
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation.
#' @return list with group mean, standard deviation, and relative standard deviation for all peaks, and filtered peaks index
#' @examples
#' data(list)
#' getdoe(list)
#' @export
#' @seealso \code{\link{getimputation}}, \code{\link{getpower}}
getdoe <- function(list,
                   inscf = 5,
                   rsdcf = 100,
                   rsdcft = 30,
                   imputation = "l",
                   tr = FALSE,
                   BPPARAM = BiocParallel::bpparam()) {
        list <- getimputation(list, method = imputation)
        # remove the technical replicates and use biological
        # replicates instead
        if (tr) {
                data <- list$data
                lv <- list$group[,-1, drop = FALSE]
                # group base on levels
                cols <- colnames(lv)
                mlv <- do.call(paste, c(lv[cols]))
                # get the rsd of the technical replicates
                meant <-
                        BiocParallel::bpaggregate(t(data), list(mlv),  mean, BPPARAM = BPPARAM)
                idx <- match(unique(mlv), meant[, 1])
                meant <- meant[idx,]
                sdt <-
                        BiocParallel::bpaggregate(t(data), list(mlv), stats::sd, BPPARAM = BPPARAM)
                idx <- match(unique(mlv), sdt[, 1])
                sdt <- sdt[idx,]
                suppressWarnings(rsd <- sdt[,-1] / meant[,-1] *
                                         100)
                data <- t(meant[,-1])
                colnames(data) <- unique(mlv)
                rsd <- t(rsd)
                # filter the data based on rsd of the technical
                # replicates
                indext <- as.vector(apply(rsd, 1, function(x)
                        all(x <
                                    rsdcft)))
                indext <- indext & (!is.na(indext))
                data <- data[indext,]
                # data with mean of the technical replicates
                list$data <- data
                list$mz <- list$mz[indext]
                list$rt <- list$rt[indext]
                # get new group information
                ng <- NULL
                if (ncol(lv) > 1) {
                        for (i in 1:(ncol(lv) - 1)) {
                                lvi <- vapply(
                                        strsplit(
                                                unique(mlv),
                                                split = " ",
                                                fixed = TRUE
                                        ),
                                        `[`,
                                        i,
                                        'g'
                                )
                                ng <- cbind(ng, lvi)
                        }
                        list$group <-
                                cbind.data.frame(
                                        sample_name = unique(mlv),
                                        ng,
                                        stringsAsFactors = FALSE
                                )
                } else {
                        list$group <-
                                data.frame(
                                        sample_name = unique(mlv),
                                        sample_group = unique(mlv),
                                        stringsAsFactors = FALSE
                                )
                }
                # save the index
                list$techindex <- indext
        }

        # filter the data based on rsd/intensity
        data <- list$data
        lv <- list$group[,-1, drop = FALSE]
        cols <- colnames(lv)
        # one peak for metabolomics is hard to happen
        if (sum(NROW(lv) > 1) != 0) {
                if (sum(NCOL(lv) > 1)) {
                        mlv <- do.call(paste0, c(lv[cols], sep = ""))
                } else {
                        mlv <- unlist(lv)
                }
                mean <-
                        BiocParallel::bpaggregate(t(data), list(mlv), mean, BPPARAM = BPPARAM)
                idx <- match(unique(mlv), mean[, 1])
                mean <- mean[idx,]
                sd <-
                        BiocParallel::bpaggregate(t(data), list(mlv), stats::sd, BPPARAM = BPPARAM)
                idx <- match(unique(mlv), sd[, 1])
                sd <- sd[idx,]
                suppressWarnings(rsd <- sd[,-1] / mean[,-1] * 100)
                mean <- t(mean[,-1])
                sd <- t(sd[,-1])
                rsd <- t(rsd)
                colnames(rsd) <-
                        colnames(sd) <-
                        colnames(mean) <- unique(mlv)
                indexrsd <- as.vector(apply(rsd, 1, function(x)
                        all(x <
                                    rsdcf)))
                indexins <- as.vector(apply(mean, 1, function(x)
                        any(x >
                                    10 ^ (
                                            inscf
                                    ))))
                list$groupmean <- mean
                list$groupsd <- sd
                list$grouprsd <- rsd
                indexrsd[is.na(indexrsd)] <- FALSE
                list$rsdindex <- indexrsd
                list$insindex <- indexins
                return(list)
        } else {
                indexins <- data > 10 ^ (inscf)
                list$groupmean <- apply(data, 1, mean)
                list$groupsd <- apply(data, 1, sd)
                suppressWarnings(list$grouprsd <-
                                         list$groupsd / list$groupmean * 100)
                list$insindex <- indexins
                message("Only technical replicates were shown for ONE sample !!!")
                return(list)
        }
}

#' Get the index with power restriction for certain study with BH adjusted p-value and certain power.
#' @param list list with data as peaks list, mz, rt and group information
#' @param pt p value threshold, default 0.05
#' @param qt q value threshold, BH adjust, default 0.05
#' @param powert power cutoff, default 0.8
#' @param imputation parameters for `getimputation` function method
#' @return list with current power and sample numbers for each peaks
#' @examples
#' data(list)
#' getpower(list)
#' @export
#' @seealso \code{\link{getimputation}}, \code{\link{getdoe}}
getpower <-
        function(list,
                 pt = 0.05,
                 qt = 0.05,
                 powert = 0.8,
                 imputation = "l") {
                group <- as.factor(list$group[,-1])
                g <- unique(group)
                ng <- length(g)
                n <- min(table(group))
                list <- getdoe(list, imputation = imputation)
                sd <- apply(list$groupsd, 1, mean)
                if (ng == 2) {
                        ar <- apply(list$data, 1, function(x)
                                stats::t.test(x ~ group))
                        dm <-
                                vapply(ar, function(x)
                                        x$estimate[1] - x$estimate[2], 1)
                        m <- nrow(list$data)
                        p <- vapply(ar, function(x)
                                x$p.value, 1)
                        q <- stats::p.adjust(p, method = "BH")
                        qc <- c(1:m) * pt / m
                        cf <- qc[match(order(qc), order(q))]
                        re <- stats::power.t.test(
                                delta = dm,
                                sd = sd,
                                sig.level = cf,
                                n = n
                        )
                        n <- vector()
                        for (i in 1:m) {
                                re2 <- try(stats::power.t.test(
                                        delta = dm[i],
                                        sd = sd[i],
                                        sig.level = cf[i],
                                        power = powert
                                ),
                                silent = TRUE)
                                if (inherits(re2, "try-error"))
                                        n[i] <- NA
                                else
                                        n[i] <- re2$n
                        }

                        list$power <- re$power
                        list$n <- n
                } else{
                        sdg <- apply(list$groupmean, 1, function(x)
                                stats::sd(x))
                        ar <-
                                apply(list$data, 1, function(x)
                                        stats::anova(stats::lm(x ~ group)))
                        p <- vapply(ar, function(x)
                                x$`Pr(>F)`[1], 1)
                        m <- nrow(list$data)
                        q <- stats::p.adjust(p, method = "BH")
                        qc <- c(1:m) * pt / m
                        cf <- qc[match(order(qc), order(q))]
                        re <- stats::power.anova.test(
                                groups = ng,
                                between.var = sdg,
                                within.var = sd,
                                sig.level = cf,
                                n = n
                        )
                        n <- vector()
                        for (i in 1:m) {
                                re2 <- try(stats::power.anova.test(
                                        groups = ng,
                                        between.var = sdg[i],
                                        within.var = sd[i],
                                        sig.level = cf[i],
                                        power = powert
                                ),
                                silent = TRUE)
                                if (inherits(re2, "try-error"))
                                        n[i] <- NA
                                else
                                        n[i] <- re2$n
                        }
                        list$power <- re$power
                        list$n <- n
                }
                return(list)
        }

#' Density weighted intensity for one sample
#' @param peak peaks intensity one sample
#' @param n the number of equally spaced points at which the density is to be estimated, default 512
#' @param log log transformation
#' @return Density weighted intensity for one sample
#' @examples
#' data(list)
#' getdwtus(list$data[,1])
#' @export
#'
getdwtus <- function(peak, n = 512, log = FALSE) {
        if (log) {
                peak <- log(peak + 1)
        }
        sum <-
                sum(stats::density(peak, bw = 'sj', n = n)$x * stats::density(peak, bw =
                                                                                      'sj', n = n)$y)
        return(sum)
}

#' Compute pooled QC linear index according to run order
#' @param data peaks intensity list with row as peaks and column as samples
#' @param order run order of pooled QC samples
#' @param n samples numbers used for linear regression
#' @return vector for the peaks proportion with significant changes in linear regression after FDR control.
#' @export

getpqsi <- function(data, order, n = 5) {
        data <- data[, order(order)]
        porp <- seq_len(ncol(data))
        for (i in n:ncol(data)) {
                p <-
                        apply(data[, c((i - n + 1):i)], 1, function(x)
                                summary(stats::lm(x ~ c((i - n + 1):i
                                )))$coefficients[2, 4])
                # FDR control
                q <- stats::p.adjust(p, method = 'BH')
                porp[i] <- sum(q < 0.1) / nrow(data)
        }
        return(porp[-c(1:(n - 1))])
}

#' Perform peaks list alignment and return features table
#' @param list each element should be a data.frame with mz, rt and ins as m/z, retention time in seconds and intensity of certain peaks.
#' @param ts template sample index in the list, default 1
#' @param ppm mass accuracy, default 10
#' @param deltart retention time shift table, default 5 seconds
#' @param FUN function to deal with multiple aligned peaks from one sample
#' @return mzrt object without group information
#' @export
getretcor <- function(list,
                      ts = 1,
                      ppm = 10,
                      deltart = 5,
                      FUN) {
        nli <- list[-ts]
        csd <- list[[ts]]
        i <- 1
        df1 <- csd
        ins <- df1$ins
        while (i <= length(nli)) {
                df2 <- list[[i]]
                df <-
                        enviGCMS::getalign(df1$mz,
                                           df2$mz,
                                           df1$rt,
                                           df2$rt,
                                           ppm = ppm,
                                           deltart = deltart)
                mr2 <- paste0(df2$mz, '@', df2$rt)
                mrx <- paste0(df$mz2, '@', df$rt2)

                df$ins2 <- df2$ins[match(mrx, mr2)]
                dfx <- df[!duplicated(df$xid),]
                dfx$ins2 <-
                        stats::aggregate(df$ins2, by = list(df$xid), FUN)[, 2]
                df1 <- cbind.data.frame(mz = dfx$mz1, rt = dfx$rt1)
                if (length(dim(ins)) > 1) {
                        insn <- ins[df$xid,]
                        ins <-
                                cbind.data.frame(insn[!duplicated(df$xid),], dfx$ins2)
                } else{
                        insn <- ins[df$xid]
                        ins <-
                                cbind.data.frame(ins1 = insn[!duplicated(df$xid)], dfx$ins2)
                }
                i <- i + 1
        }
        colnames(ins) <- paste0('ins', seq_along(list))
        li <- list(data = ins,
                   mz = df1[, 1],
                   rt = df1[, 2])
        return(li)
}

#' Merge positive and negative mode data
#' @param pos a list with mzrt profile collected from positive mode. The sample order should match the negative mode.
#' @param neg a list with mzrt profile collected from negative mode.The sample order should match the positive mode.
#' @param ppm pmd mass accuracy, default 5
#' @param pmd numeric or numeric vector
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param cutoff correlation coefficients, default 0.9
#' @return mzrt object with group information from pos mode
#' @export
getpn <- function(pos,
                  neg,
                  ppm = 5,
                  pmd = 2.02,
                  digits = 2,
                  cutoff = 0.9) {
        df <- NULL
        x <- rep(NA, length(pos$mz))
        for (i in seq_along(pos$mz)) {
                diff <- round((pos$mz[i] - neg$mz), digits)
                mzr <-
                        data.table::as.data.table(cbind.data.frame(
                                mzmin = pmd * (1 - ppm * 10e-6),
                                mzmax = pmd * (1 + ppm * 10e-6)
                        ))
                colnames(mzr) <- c("min", "max")
                data.table::setkey(mzr, min, max)
                diff <-
                        data.table::data.table(min = diff, max = diff)
                overlapms <-
                        data.table::foverlaps(
                                diff,
                                mzr,
                                type = 'within',
                                nomatch = 0,
                                which = TRUE
                        )
                if (nrow(overlapms) != 0) {
                        if (nrow(overlapms) > 1) {
                                cor <-
                                        apply(neg$data[overlapms$xid,], 1, function(x)
                                                suppressWarnings(
                                                        cor(
                                                                as.numeric(x),
                                                                as.numeric(pos$data[i,])
                                                        )
                                                ))
                        } else{
                                cor <-
                                        suppressWarnings(cor(
                                                as.numeric(pos$data[i,]),
                                                as.numeric(neg$data[overlapms$xid,])
                                        ))
                        }

                        t <-
                                cbind.data.frame(
                                        pos = pos$mz[i],
                                        rt = pos$rt[i],
                                        neg = neg$mz[overlapms$xid],
                                        rt = neg$rt[overlapms$xid],
                                        diffmz = pos$mz[i] - neg$mz[overlapms$xid],
                                        diffrt = pos$rt[i] - neg$rt[overlapms$xid],
                                        cor = cor
                                )
                        df <- rbind.data.frame(df, t)
                }
        }
        merge2 <- df[df$cor > cutoff,]
        idxp <-
                match(paste(merge2$pos, merge2[, 2]), paste(pos$mz, pos$rt))
        idxn <-
                match(paste(merge2$neg, merge2[, 4]), paste(neg$mz, neg$rt))
        pos$anno[idxp]
        neg$anno[idxn]
        colnames(neg$data) <- colnames(pos$data)


        if (is.null(pos$anno) & is.null(neg$anno)) {
                meta <-
                        list(
                                data = rbind(pos$data[-idxp,], neg$data),
                                group = pos$group,
                                mz = c(pos$mz[-idxp], neg$mz),
                                rt = c(pos$rt[-idxp], neg$rt),
                                mname = c(
                                        paste0('P_', rownames(pos$data)[-idxp]),
                                        paste0('N_', rownames(neg$data))
                                )
                        )
        } else{
                meta <-
                        list(
                                data = rbind(pos$data[-idxp,], neg$data),
                                group = pos$group,
                                mz = c(pos$mz[-idxp], neg$mz),
                                rt = c(pos$rt[-idxp], neg$rt),
                                mname = c(
                                        paste0('P_', rownames(pos$data)[-idxp]),
                                        paste0('N_', rownames(neg$data))
                                ),
                                anno = c(pos$anno[-idxp], neg$anno)
                        )
        }
        return(meta)
}

#' Align multiple peaks list to one peak list
#' @param ... peaks list, mzrt objects
#' @param index numeric, the index of reference peaks.
#' @param ppm pmd mass accuracy, default 5
#' @param deltart retention time shift table, default 10 seconds
#' @return list object with aligned mzrt objects
#' @export
getcompare <- function(...,
                       index = 1,
                       ppm = 5,
                       deltart = 5) {
        li <- list(...)
        ref <- li[[index]]
        lire <- li[-index]
        z <- lapply(lire, function(x) {
                over <-
                        enviGCMS::getalign(
                                mz1 = x$mz,
                                mz2 = ref$mz,
                                rt1 = x$rt,
                                rt2 = ref$rt,
                                ppm = ppm,
                                deltart = deltart
                        )
                over2 <-
                        enviGCMS::getalign(
                                mz1 = ref$mz,
                                mz2 = x$mz,
                                rt1 = ref$rt,
                                rt2 = x$rt,
                                ppm = ppm,
                                deltart = deltart
                        )
                refi <-
                        enviGCMS::getfilter(x, rowindex = -unique(over$xid))
                refo <-
                        enviGCMS::getfilter(ref, rowindex = unique(over2$xid))
                ni <- c(refi, refo)
                return(ni)
        })
        re <- append(list(ref), z)
        return(re)
}
