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
                group <- data.frame(apply(XCMSnExp@phenoData@data, 2, as.character),stringsAsFactors = FALSE)

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
                        if(sum(grepl('sample_name',colnames(xcms::phenoData(xcmsSet))))>0){
                                group <- xcms::phenoData(xcmsSet)
                                colnames(data) <- group$sample_name
                        }else{
                                sample_name <- rownames(xcms::phenoData(xcmsSet))
                                sample_group <- xcms::phenoData(xcmsSet)
                                group <- cbind.data.frame(sample_name,sample_group)
                        }
                }else{
                        sample_name <- rownames(xcms::phenoData(xcmsSet))
                        sample_group <- as.character(xcms::phenoData(xcmsSet)$class)
                        group <- cbind.data.frame(sample_name,sample_group)
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
                                data <- rbind(sample_group, list$data)
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
                                lv <- list$group[,-1]
                                lif <- getimputation(list, method = "l")
                                ar <-
                                        apply(lif$data,1, function(x) stats::anova(stats::lm(x~lv)))

                                df <-
                                        cbind.data.frame(
                                                m.z = list$mz,
                                                rt = list$rt,
                                                p.value = vapply(ar,function(x) x$`Pr(>F)`[1],1),
                                                t.score = vapply(ar,function(x) x$`F value`[1],1)
                                        )
                                filename <- paste0(name, 'mummichog.txt')
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
                                groupt <- c('mz', 'rt', list$group[,-1])
                                data <- rbind(groupt, data)
                                if(!target){
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
getrangecsv <- function(list,name,...){
        df <- cbind.data.frame(list$mzrange,list$rtrange)
        filename <- paste0(name, "mzrtrange.csv")
        utils::write.csv(df, file = filename, ...)
}
#' Get the mzrt profile and group information as a mzrt list and/or save them as csv or rds for further analysis.
#' @param xset xcmsSet/XCMSnExp objects
#' @param name file name for csv and/or eic file, default NULL
#' @param mzdigit m/z digits of row names of data frame, default 4
#' @param rtdigit retention time digits of row names of data frame, default 1
#' @param method parameter for groupval or featureDefinitions function, default medret
#' @param value parameter for groupval or featureDefinitions function, default into
#' @param eic logical, save xcmsSet and xcmsEIC objects for further investigation with the same name of files, you will need raw files in the same directory as defined in xcmsSet to extract the EIC based on the binned data. You could use `plot` to plot EIC for specific peaks. For example, `plot(xcmsEIC,xcmsSet,groupidx = 'M123.4567T278.9')` could show the EIC for certain peaks with m/z 206 and retention time 2789. default F
#' @param type csv format for further analysis, m means  Metaboanalyst, a means xMSannotator, p means Mummichog(NA values are imputed by `getimputation`, and F test is used here to generate stats and p value), o means full information csv (for `pmd` package), default o. mapo could output all those format files.
#' @return mzrt object, a list with mzrt profile and group information
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
                 eic = FALSE,
                 type = 'o') {
                if (class(xset) == 'xcmsSet') {
                        if (eic) {
                                eic <-
                                        xcms::getEIC(
                                                xset,
                                                rt = "corrected",
                                                groupidx = seq_len(nrow(xset@groups))
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
                                                groupidx = seq_len(nrow(xset2@groups))
                                        )
                                saveRDS(eic, file = paste0(name, 'eic.rds'))
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
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}}, \code{\link{getcsv}}
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
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}},\code{\link{getpower}}
getdoe <- function(list,
                   inscf = 5,
                   rsdcf = 100,
                   rsdcft = 30,
                   imputation = "l",
                   tr = FALSE, BPPARAM = BiocParallel::bpparam()) {
        list <- getimputation(list, method = imputation)
        # remove the technical replicates and use biological
        # replicates instead
        if (tr) {
                data <- list$data
                lv <- list$group[,-1,drop=FALSE]
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
                # get new group information
                ng <- NULL
                if (ncol(lv) > 1) {
                        for (i in 1:(ncol(lv) - 1)) {
                                lvi <- vapply(strsplit(
                                        unique(mlv),
                                        split = " ",
                                        fixed = TRUE
                                ),
                                `[`,
                                i,'g')
                                ng <- cbind(ng, lvi)
                        }
                        list$group <- cbind.data.frame(sample_name = unique(mlv), ng,stringsAsFactors = FALSE)
                } else {
                        list$group <- data.frame(sample_name = unique(mlv),sample_group = unique(mlv),stringsAsFactors = FALSE)
                }
                # save the index
                list$techindex <- indext
        }

        # filter the data based on rsd/intensity
        data <- list$data
        lv <- list$group[,-1,drop=FALSE]
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
#' @seealso \code{\link{getdata2}},\code{\link{getdata}}, \code{\link{getmzrt}}, \code{\link{getimputation}}, \code{\link{getmr}},\code{\link{getdoe}}
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
                        ar <- apply(list$data, 1, function(x) stats::t.test(x~group))
                        dm <- vapply(ar,function(x) x$estimate[1]-x$estimate[2],1)
                        m <- nrow(list$data)
                        p <- vapply(ar,function(x) x$p.value,1)
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
                        sdg <- apply(list$groupmean, 1, function(x) stats::sd(x))
                        ar <-
                                apply(list$data,1, function(x) stats::anova(stats::lm(x~group)))
                        p <- vapply(ar,function(x) x$`Pr(>F)`[1],1)
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
getdwtus <- function(peak,n=512,log=FALSE){
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
        porp <- seq_len(ncol(data))
        for(i in n:ncol(data)){
                p <- apply(data[,c((i-n+1):i)],1,function(x) summary(stats::lm(x~c((i-n+1):i)))$coefficients[2,4])
                # FDR control
                q <- stats::p.adjust(p,method='BH')
                porp[i] <- sum(q<0.1)/nrow(data)
        }
        return(porp[-c(1:(n-1))])
}

#' Perform peaks list alignment and return features table
#' @param list each element should be a data.frame with mz, rt and ins as m/z, retention time in seconds and intensity of certain peaks.
#' @param ts template sample index in the list, default 1
#' @param ppm mass accuracy, default 10
#' @param deltart retention time shift table, default 5 seconds
#' @param FUN function to deal with multiple aligned peaks from one sample
#' @return mzrt object without group information
#' @export
getretcor <- function(list, ts=1, ppm=10, deltart=5, FUN){
        nli <- list[-ts]
        csd <- list[[ts]]
        i <- 1
        df1 <- csd
        ins <- df1$ins
        while(i<=length(nli)){
                df2 <- list[[i]]
                df <- enviGCMS::getalign(df1$mz,df2$mz,df1$rt,df2$rt,ppm=ppm,deltart=deltart)
                mr2 <- paste0(df2$mz,'@',df2$rt)
                mrx <- paste0(df$mz2,'@',df$rt2)

                df$ins2 <- df2$ins[match(mrx,mr2)]
                dfx <- df[!duplicated(df$xid),]
                dfx$ins2 <- stats::aggregate(df$ins2,by=list(df$xid),FUN)[,2]
                df1 <- cbind.data.frame(mz=dfx$mz1,rt=dfx$rt1)
                if(length(dim(ins))>1){
                        insn <- ins[df$xid,]
                        ins <- cbind.data.frame(insn[!duplicated(df$xid),],dfx$ins2)
                }else{
                        insn <- ins[df$xid]
                        ins <- cbind.data.frame(ins1=insn[!duplicated(df$xid)],dfx$ins2)
                }
                i <- i+1
        }
        colnames(ins) <- paste0('ins',seq_along(list))
        li <- list(data = ins, mz=df1[,1], rt=df1[,2])
        return(li)
}
