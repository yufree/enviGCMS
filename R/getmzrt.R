#' Convert an XCMSnExp object to an mzrt S3 object.
#'
#' @noRd
.XCMSnExp2mzrt <-
        function(XCMSnExp,
                 method = "medret",
                 value = "into",
                 mzdigit = 4,
                 rtdigit = 1)
        {
                data <- xcms::featureValues(XCMSnExp, value = value, missing = 0)
                group <- data.frame(apply(XCMSnExp@phenoData@data, 2, as.character),stringsAsFactors = F)

                # peaks info
                peaks <- xcms::featureDefinitions(XCMSnExp)
                mz <- peaks$mzmed
                rt <- peaks$rtmed
                mzrange <- peaks[, c("mzmin", "mzmax")]
                rtrange <- peaks[, c("rtmin", "rtmax")]
                rownames(data) <-
                        paste0("M", round(mz, mzdigit), "T", round(rt, rtdigit))
                mzrt <-
                        list(
                                data = data,
                                group = group,
                                mz = mz,
                                rt = rt,
                                mzrange = mzrange,
                                rtrange = rtrange
                        )
                class(mzrt) <- "mzrt"
                return(mzrt)
        }
#' Convert an xcmsSet object to an mzrt S3 object.
#'
#' @noRd
.xcmsSet2mzrt <-
        function(xcmsSet,
                 method = "medret",
                 value = "into",
                 mzdigit = 4,
                 rtdigit = 1)
        {
                data <- xcms::groupval(xcmsSet, method = method,
                                       value = value)
                rownames(data) <- xcms::groupnames(xcmsSet,
                                                   template = paste0(
                                                           'M.',
                                                           10 ^ (mzdigit - 1),
                                                           'T.',
                                                           10 ^ (rtdigit - 1)
                                                   ))
                # peaks info
                peaks <- as.data.frame(xcms::groups(xcmsSet))
                mz <- peaks$mzmed
                rt <- peaks$rtmed
                if(dim(xcms::phenoData(xcmsSet))[2]>1){
                        group <- xcms::phenoData(xcmsSet)
                }else{
                        group <- as.character(xcms::phenoData(xcmsSet)$class)
                }

                mzrange <- peaks[, c("mzmin", "mzmax")]
                rtrange <- peaks[, c("rtmin", "rtmax")]
                mzrt <-
                        list(
                                data = data,
                                group = group,
                                mz = mz,
                                rt = rt,
                                mzrange = mzrange,
                                rtrange = rtrange
                        )
                class(mzrt) <- "mzrt"
                return(mzrt)
        }

#' Convert an list object to csv file.
#' @param list list with data as peaks list, mz, rt and group information
#' @param name result name for csv and/or eic file, default NULL
#' @param mzdigit m/z digits of row names of data frame, default 4
#' @param rtdigit retention time digits of row names of data frame, default 1
#' @param type csv formate for furthor analysis, m means  Metaboanalyst, a means xMSannotator, p means Mummichog(NA values are imputed by `getimputation`, and F test is used here to generate stats and p vlaue), o means full infomation csv (for `pmd` package), default o. mapo could output all those format files.
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
                 ...) {
                if (!is.null(name)) {
                        if (grepl('m', type)) {
                                data <- rbind(list$group, list$data)
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
                                filename <- paste0(name, "xMSannotator.csv")
                                utils::write.csv(data, file = filename, ...)
                        }
                        if (grepl('p', type)) {
                                lv <- list$group
                                lif <- getimputation(list, method = "l")
                                fstats <-
                                        genefilter::rowFtests(as.matrix(lif$data), fac = as.factor(lv))
                                df <-
                                        cbind.data.frame(
                                                m.z = list$mz,
                                                rt = list$rt,
                                                p.value = fstats$p.value,
                                                t.score = fstats$statistic
                                        )
                                filename <- paste0(name, 'mummichog.txt')
                                utils::write.table(
                                        df,
                                        file = filename,
                                        sep = "\t",
                                        row.names = F,
                                        ...
                                )
                        }
                        if (grepl('o', type)) {
                                data <- cbind(mz = list$mz,
                                              rt = list$rt,
                                              list$data)
                                colname <- colnames(data)
                                groupt = c('mz', 'rt', list$group)
                                data <- rbind(groupt, data)
                                rownames(data) <-
                                        c('group',
                                          paste0(
                                                  "M",
                                                  round(list$mz, mzdigit),
                                                  "T",
                                                  round(list$rt, rtdigit)
                                          ))
                                colnames(data) <- colname
                                filename <- paste0(name, "mzrt.csv")
                                utils::write.csv(data, file = filename, ...)
                        }
                }
        }
#' Get the mzrt profile and group information as a mzrt list and/or save them as csv or rds for further analysis.
#' @param xset xcmsSet/XCMSnExp objects
#' @param name file name for csv and/or eic file, default NULL
#' @param mzdigit m/z digits of row names of data frame, default 4
#' @param rtdigit retention time digits of row names of data frame, default 1
#' @param method parameter for groupval or featureDefinitions function, default medret
#' @param value parameter for groupval or featureDefinitions function, default into
#' @param eic logical, save xcmsSet and xcmsEIC objects for further investigation with the same name of files, you will need raw files in the same directory as defined in xcmsSet to extract the EIC based on the binned data. You could use `plot` to plot EIC for specific peaks. For example, `plot(xcmsEIC,xcmsSet,groupidx = 'M123.4567T278.9')` could show the EIC for certain peaks with m/z 206 and retention time 2789. default F
#' @param type csv formate for furthor analysis, m means  Metaboanalyst, a means xMSannotator, p means Mummichog(NA values are imputed by `getimputation`, and F test is used here to generate stats and p vlaue), o means full infomation csv (for `pmd` package), default o. mapo could output all those format files.
#' @return mzrt object, a list with mzrt profile and group infomation
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getmzrt(xset, name = 'demo', type = 'mapo')
#' }
#' @seealso \code{\link{getdata}},\code{\link{getdata2}}, \code{\link{getdoe}}, \code{\link{getcsv}}, \code{\link{getfilter}}
#' @references
#' Smith, C.A., Want, E.J., O’Maille, G., Abagyan, R., Siuzdak, G., 2006. XCMS: Processing Mass Spectrometry Data for Metabolite Profiling Using Nonlinear Peak Alignment, Matching, and Identification. Anal. Chem. 78, 779–787.
#' @export

getmzrt <-
        function(xset,
                 name = NULL,
                 mzdigit = 4,
                 rtdigit = 1,
                 method = "medret",
                 value = "into",
                 eic = F,
                 type = 'o') {
                if (class(xset) == 'xcmsSet') {
                        if (eic) {
                                eic <-
                                        xcms::getEIC(
                                                xset,
                                                rt = "corrected",
                                                groupidx = 1:nrow(xset@groups)
                                        )
                                eic@groupnames <-
                                        xcms::groupnames(xset,
                                                         template = paste0(
                                                                 'M.',
                                                                 10 ^ (mzdigit - 1),
                                                                 'T.',
                                                                 10 ^ (rtdigit - 1)
                                                         ))
                                saveRDS(eic, file = paste0(name, 'eic.rds'))
                                saveRDS(xset, file = paste0(name, 'xset.rds'))
                        }
                        result <-
                                .xcmsSet2mzrt(
                                        xset,
                                        mzdigit = mzdigit,
                                        rtdigit = rtdigit,
                                        method = method,
                                        value = value
                                )
                }
                else if (class(xset) == 'XCMSnExp') {
                        xset2 <- methods::as(xset, 'xcmsSet')
                        if (eic) {
                                eic <-
                                        xcms::getEIC(
                                                xset2,
                                                rt = "corrected",
                                                groupidx = 1:nrow(xset2@groups)
                                        )
                                saveRDS(eic, file = paste0(name, '.rds'))
                                saveRDS(xset2, file = paste0(name, 'xset.rds'))
                        }
                        result <-
                                .XCMSnExp2mzrt(
                                        xset,
                                        mzdigit = mzdigit,
                                        rtdigit = rtdigit,
                                        method = method,
                                        value = value
                                )
                }
                getcsv(
                        list = result,
                        name = name,
                        mzdigit = mzdigit,
                        rtdigit = rtdigit,
                        type = type
                )
                return(result)
        }
#' Impute the peaks list data
#' @param list list with data as peaks list, mz, rt and group information
#' @param method 'r' means remove, 'l' means use half the minimum of the values across the peaks list, 'mean' means mean of the values across the samples, 'median' means median of the values across the samples, '0' means 0, '1' means 1. Default 'l'.
#' @return list with imputed peaks
#' @examples
#' data(list)
#' getimputation(list)
#' @export
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}},\code{\link{getdoe}}, \code{\link{getmr}}
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
                impute <- min(data, na.rm = T) / 2
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
#' Filter the data based on row and column index
#' @param list list with data as peaks list, mz, rt and group information
#' @param rowindex logical, row index to keep
#' @param colindex logical, column index to keep
#' @param name file name for csv and/or eic file, default NULL
#' @param type csv formate for furthor analysis, m means  Metaboanalyst, a means xMSannotator, p means Mummichog(NA values are imputed by `getimputation`, and F test is used here to generate stats and p vlaue), o means full infomation csv (for `pmd` package), default o. mapo could output all those format files.
#' @param ... other parameters for `getcsv`
#' @return list with remain peaks, and filtered peaks index
#' @examples
#' data(list)
#' li <- getdoe(list)
#' lif <- getfilter(li,rowindex = li$rsdindex)
#' @export
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}}, \code{\link{getcsv}}
getfilter <-
        function(list,
                 rowindex = T,
                 colindex = T,
                 name = NULL,
                 type = 'o',
                 ...) {
                if (!is.null(rowindex) & !is.null(list$rowindex)) {
                        rowindex <- rowindex & list$rowindex
                }
                else if (is.null(rowindex) & !is.null(list$rowindex)) {
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
                else if (is.null(colindex) & !is.null(list$colindex)) {
                        colindex <- list$colindex
                }
                list$data <- list$data[, colindex]
                if(NCOL(list$group)>1) {
                        list$group <- list$group[colindex,]
                }else{
                        list$group <- list$group[colindex]
                        }
                # list$colindex <- colindex
                getcsv(list, name = name, type = type, ...)
                return(list)
        }
#' Filter the data based on DoE, rsd, intensity
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
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}},\code{\link{getpower}}
getdoe <- function(list,
                   inscf = 5,
                   rsdcf = 100,
                   rsdcft = 30,
                   imputation = "l",
                   tr = F, BPPARAM = BiocParallel::bpparam()) {
        list <- getimputation(list, method = imputation)
        # remove the technical replicates and use biological
        # replicates instead
        if (tr) {
                data <- list$data
                lv <- list$group
                # group base on levels
                cols <- colnames(lv)
                mlv <- do.call(paste, c(lv[cols]))
                # get the rsd of the technical replicates
                meant <- BiocParallel::bpaggregate(t(data), list(mlv),  mean,BPPARAM = BPPARAM)
                idx <- match(unique(mlv), meant[,1])
                meant <- meant[idx,]
                sdt <- BiocParallel::bpaggregate(t(data), list(mlv), stats::sd,BPPARAM = BPPARAM)
                idx <- match(unique(mlv), sdt[,1])
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
                # get new group infomation
                ng <- NULL
                if (ncol(lv) > 1) {
                        for (i in 1:(ncol(lv) - 1)) {
                                lvi <- sapply(strsplit(
                                        unique(mlv),
                                        split = " ",
                                        fixed = TRUE
                                ),
                                `[`,
                                i)
                                ng <- cbind(ng, lvi)
                        }
                        list$group <- data.frame(ng)
                } else {
                        list$group <- data.frame(unique(mlv))
                }
                # save the index
                list$techindex <- indext
        }

        # filter the data based on rsd/intensity
        data <- list$data
        lv <- list$group
        cols <- colnames(lv)
        # one peak for metabolomics is hard to happen
        if (sum(NROW(lv) > 1) != 0) {
                if (sum(NCOL(lv) > 1)) {
                        mlv <- do.call(paste0, c(lv[cols], sep = ""))
                } else {
                        mlv <- unlist(lv)
                }
                mean <- BiocParallel::bpaggregate(t(data), list(mlv), mean,BPPARAM = BPPARAM)
                idx <- match(unique(mlv), mean[,1])
                mean <- mean[idx,]
                sd <- BiocParallel::bpaggregate(t(data), list(mlv), stats::sd,BPPARAM = BPPARAM)
                idx <- match(unique(mlv), sd[,1])
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
                indexrsd[is.na(indexrsd)] <- F
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
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}},\code{\link{getdoe}}
getpower <-
        function(list,
                 pt = 0.05,
                 qt = 0.05,
                 powert = 0.8,
                 imputation = "l") {
                group <- as.factor(list$group)
                g <- unique(group)
                ng <- length(g)
                n <- min(table(group))
                list <- getdoe(list, imputation = imputation)
                sd <- apply(list$groupmean, 1, mean)
                if (ng == 2) {
                        ar <- genefilter::rowttests(list$data, fac = group)
                        dm <- ar$dm
                        m <- nrow(list$data)
                        p <- ar$p.value
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
                                silent = T)
                                if (inherits(re2, "try-error"))
                                        n[i] <- NA
                                else
                                        n[i] <- re2$n
                        }

                        list$power <- re$power
                        list$n <- n
                } else{
                        sdg <- genefilter::rowSds(list$groupmean)
                        ar <-
                                genefilter::rowFtests(list$data, group)
                        p <- ar$p.value
                        q <- stats::p.adjust(p, method = "BH")
                        m <- nrow(list$data)
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
                                        between.var = sdg,
                                        within.var = sd,
                                        sig.level = cf,
                                        power = powert
                                ),
                                silent = T)
                                if (inherits(re2, "try-error"))
                                        n[i] <- NA
                                else
                                        n[i] <- re2$n
                        }
                        list$power <- re$power
                        list$n <- re2$n
                }
                return(list)
        }

#' Get the overlap peaks by mass and retention time range
#' @param list1 list with data as peaks list, mz, rt, mzrange, rtrange and group information to be overlapped
#' @param list2 list with data as peaks list, mz, rt, mzrange, rtrange and group information to overlap
#' @return logical index for list 1's peaks
#' @export
#' @seealso \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}},\code{\link{getdoe}}, \code{\link{getoverlapmass}},\code{\link{getoverlaprt}}
getoverlappeak <- function(list1, list2) {
        mz1 <- data.table::as.data.table(list1$mzrange)
        rt1 <- data.table::as.data.table(list1$rtrange)
        mz2 <- data.table::as.data.table(list2$mzrange)
        rt2 <- data.table::as.data.table(list2$rtrange)
        colnames(mz1) <-
                colnames(mz2) <-
                colnames(rt1) <- colnames(rt2) <- c('min', 'max')
        data.table::setkey(mz2, min, max)
        data.table::setkey(rt2, min, max)
        overlapms <-
                data.table::foverlaps(mz1, mz2, which = TRUE, mult = 'first')
        overlaprt <-
                data.table::foverlaps(rt1, rt2, which = TRUE, mult = 'first')
        index <- (!is.na(overlapms)) & (!is.na(overlaprt))
        return(index)
}
#' Get the overlap peaks by mass range
#' @param mzrange1 mass range 1 to be overlapped
#' @param mzrange2 mass range 2 to overlap
#' @return logical index for mzrange1's peaks
#' @export
#' @seealso \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}},\code{\link{getdoe}}, \code{\link{getoverlappeak}},\code{\link{getoverlaprt}}
getoverlapmass <- function(mzrange1, mzrange2) {
        mz1 <- data.table::as.data.table(mzrange1)
        mz2 <- data.table::as.data.table(mzrange2)
        colnames(mz1) <- colnames(mz2) <- c('min', 'max')
        data.table::setkey(mz2, min, max)
        overlapms <-
                data.table::foverlaps(mz1, mz2, which = TRUE, mult = 'first')

        index <- (!is.na(overlapms))
        return(index)
}
#' Get the overlap peaks by retention time
#' @param rtrange1 mass range 1 to be overlapped
#' @param rtrange2 mass range 2 to overlap
#' @return logical index for rtrange1's peaks
#' @export
#' @seealso \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}},\code{\link{getdoe}}, \code{\link{getoverlapmass}},\code{\link{getoverlappeak}}
getoverlaprt <- function(rtrange1, rtrange2) {
        rt1 <- data.table::as.data.table(rtrange1)
        rt2 <- data.table::as.data.table(rtrange2)
        colnames(rt1) <- colnames(rt2) <- c('min', 'max')
        data.table::setkey(rt2, min, max)
        overlapms <-
                data.table::foverlaps(rt1, rt2, which = TRUE, mult = 'first')

        index <- (!is.na(overlapms))
        return(index)
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
getdwtus <- function(peak,n=512,log=F){
        if(log){
                peak <- log(peak+1)
        }
        sum <- sum(stats::density(peak,bw='sj',n=n)$x*stats::density(peak,bw='sj',n=n)$y)
        return(sum)
}

#' Compute pooled QC linear index according to run order
#' @param data peaks intensity list with row as peaks and column as samples
#' @param order run order of pooled QC samples
#' @param n samples numbers used for linear regression
#' @return vector for the peaks proportion with significant changes in linear regression after FDR control.
#' @export

getpqsi <- function(data, order, n=5){
        data <- data[,order(order)]
        porp <- 1:ncol(data)
        for(i in n:ncol(data)){
                p <- apply(data[,c((i-n+1):i)],1,function(x) summary(stats::lm(x~c((i-n+1):i)))$coefficients[2,4])
                # FDR control
                q <- stats::p.adjust(p,method='BH')
                porp[i] <- sum(q<0.1)/nrow(data)
        }
        return(porp[-c(1:(n-1))])
}
