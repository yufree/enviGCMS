#' Get the peak list with blank samples' peaks removed
#' @param xset the xcmsset object with blank and certain group samples' data
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf Log intensity cutoff for peaks, default 5
#' @return diff report
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getbgremove(xset)
#' }
#' @export
getbgremove <-
        function(xset,
                 method = "medret",
                 intensity = "into",
                 file = NULL,
                 rsdcf = 30,
                 inscf = 5) {
                data0 <- xcms::groupval(xset, method, intensity)
                data0[is.na(data0)] = 0
                data <- t(data0)
                lv <- xset@phenoData[, 1]
                mean <- stats::aggregate(data, list(lv), mean)
                sd <- stats::aggregate(data, list(lv), sd)
                suppressWarnings(rsd <- sd / mean * 100)
                result <-
                        data.frame(cbind(t(mean[,-1]), t(rsd[, -1])))
                rsd0 <- as.matrix(rsd[,-1])
                rsd0[is.nan(rsd0)] <- 0
                index <-
                        t(rsd0) < rsdcf & t(mean[, -1]) > 10 ^ (inscf)

                diff <- result[, 2] - result[, 1]
                report <-
                        cbind.data.frame(xset@groups[, 1], xset@groups[, 4], diff)
                colnames(report) <- c("mz", "time", "diff")

                N <- sum(index)
                L <- length(index)

                report <- report[index,]

                message(
                        paste(
                                N,
                                'out of',
                                L,
                                'peaks found with rsd cutoff',
                                rsdcf,
                                'and Log intensity cutoff',
                                inscf
                        )
                )
                if (!is.null(file)) {
                        utils::write.csv(report,
                                         file = paste0(file, '.csv'),
                                         row.names = F)
                }
                return(report)
        }


#' Get the report for technique replicates.
#' @param xset the xcmsset object which for all of your technique replicates for one sample
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf Log intensity cutoff for peaks, default 5
#' @return dataframe with mean, standard deviation and RSD for those technique replicates combined with raw data
#' @examples
#' \dontrun{
#' data(fish4)
#' fish <- as(fish4,"xcmsSet")
#' gettechrep(fish)
#' }
#' @export
gettechrep <-
        function(xset,
                 method = "medret",
                 intensity = "into",
                 file = NULL,
                 rsdcf = 30,
                 inscf = 5) {
                data0 <- xcms::groupval(xset, method, intensity)
                data0[is.na(data0)] = 0
                data <- t(data0)
                lv <- xset@phenoData[, 1]
                mean <- stats::aggregate(data, list(lv), mean)
                sd <- stats::aggregate(data, list(lv), sd)
                suppressWarnings(rsd <- sd / mean * 100)

                result <-
                        data.frame(cbind(t(mean[, -1]), t(sd[, -1]), t(rsd[, -1])))
                rsd0 <- as.matrix(rsd[,-1])
                rsd0[is.nan(rsd0)] <- 0
                index <-
                        t(rsd0) < rsdcf & t(mean[, -1]) > 10 ^ (inscf)
                colnames(result) <- c("mean", "sd", "rsd")
                datap <- xcms::groups(xset)
                report <- cbind.data.frame(datap, result)
                N <- sum(index)
                L <- length(index)

                if (N == 0) {
                        message(paste('No peaks found'))
                        return(NA)
                } else{
                        message(
                                paste(
                                        N,
                                        'out of',
                                        L,
                                        'peaks found with rsd cutoff',
                                        rsdcf,
                                        'and Log intensity cutoff',
                                        inscf
                                )
                        )

                        report <- report[index,]

                        if (!is.null(file)) {
                                anno <-
                                        cbind.data.frame(xset@groups[, 1], xset@groups[, 4], t(mean[, -1]))
                                colnames(anno) <-
                                        c("mz", "time", as.character(t(mean[, 1])))
                                utils::write.csv(anno,
                                                 file = paste0(file, '.csv'),
                                                 row.names = F)
                                anno <- anno[index,]
                                return(anno)
                        } else {
                                return(report)
                        }
                }
        }

#' Get the report for biological replicates.
#' @param xset the xcmsset object which for all of your technique replicates for bio replicated sample in single group
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf Log intensity cutoff for peaks, default 5
#' @return dataframe with mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data
#' @export
getbiotechrep <-
        function(xset,
                 method = "medret",
                 intensity = "into",
                 file = NULL,
                 rsdcf = 30,
                 inscf = 5) {
                data0 <- xcms::groupval(xset, method, intensity)
                data0[is.na(data0)] = 0
                data <- t(data0)
                lv <- xset@phenoData[, 1]
                mean <- stats::aggregate(data, list(lv), mean)
                sd <- stats::aggregate(data, list(lv), sd)
                suppressWarnings(rsd <- sd / mean * 100)
                result <-
                        data.frame(cbind(t(mean[, -1]), t(sd[,-1]), t(rsd[, -1])))
                colnames(result) <-
                        c(paste0(t(mean[, 1]), "mean"),
                          paste0(t(sd[, 1]), "sd"),
                          paste0(t(mean[, 1]), "rsd%"))

                meanB <- apply(t(mean[, -1]), 1, mean)
                sdB <- apply(t(mean[,-1]), 1, sd)
                rsdB <- sdB / meanB * 100

                rsd0 <- as.matrix(rsdB)
                rsd0[is.nan(rsd0)] <- 0

                index <- t(rsd0) < rsdcf & meanB > 10 ^ (inscf)

                datap <- xcms::groups(xset)
                report <-
                        cbind.data.frame(datap, result, meanB, sdB, rsdB)

                N <- sum(index)
                L <- length(index)

                if (N == 0) {
                        message(paste('No peaks found'))
                        return(NA)
                } else{
                        message(
                                paste(
                                        N,
                                        'out of',
                                        L,
                                        'peaks found with rsd cutoff',
                                        rsdcf,
                                        'and Log intensity cutoff',
                                        inscf
                                )
                        )

                        report <- report[index,]

                        if (!is.null(file)) {
                                anno <-
                                        cbind.data.frame(xset@groups[, 1], xset@groups[, 4], t(mean[, -1]))
                                colnames(anno) <-
                                        c("mz", "time", as.character(t(mean[, 1])))
                                utils::write.csv(anno,
                                                 file = paste0(file, '.csv'),
                                                 row.names = F)
                                anno <- anno[index,]
                                return(anno)
                        } else {
                                return(report)
                        }
                }
        }

#' Get the report for biological replicates without technique replicates.
#' @param xset the xcmsset object which for bio replicated sample in different groups
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf Log intensity cutoff for peaks, default 5
#' @return dataframe with mean, standard deviation and RSD for biological replicates combined with raw data
#' @export
getbiorep <-
        function(xset,
                 method = "medret",
                 intensity = "into",
                 file = NULL,
                 rsdcf = 30,
                 inscf = 5) {
                data0 <- xcms::groupval(xset, method, intensity)
                data0[is.na(data0)] = 0
                data <- t(data0)
                lv <- xset@phenoData[, 1]
                mean <- stats::aggregate(data, list(lv), mean)
                sd <- stats::aggregate(data, list(lv), sd)
                suppressWarnings(rsd <- sd / mean * 100)
                result <-
                        data.frame(cbind(t(mean[, -1]), t(sd[,-1]), t(rsd[, -1])))
                colnames(result) <-
                        c(paste0(t(mean[, 1]), "mean"),
                          paste0(t(sd[, 1]), "sd"),
                          paste0(t(mean[, 1]), "rsd%"))

                rsd0 <- as.matrix(rsd[,-1])
                rsd0[is.nan(rsd0)] <- 0
                indexrsd <-
                        apply(t(rsd0), 1, function(x)
                                any(x < rsdcf))
                indexmean <-
                        apply(t(mean[,-1]), 1, function(x)
                                any(x > 10 ^ (inscf)))

                index <- indexrsd & indexmean

                datap <- xcms::groups(xset)
                report <- cbind.data.frame(datap, result)

                N <- sum(index)
                L <- length(index)

                if (N == 0) {
                        message(paste('No peaks found'))
                        return(NA)
                } else{
                        message(
                                paste(
                                        N,
                                        'out of',
                                        L,
                                        'peaks found with rsd cutoff',
                                        rsdcf,
                                        'and Log intensity cutoff',
                                        inscf
                                )
                        )
                        report <- report[index,]

                        if (!is.null(file)) {
                                anno <-
                                        cbind.data.frame(xset@groups[, 1], xset@groups[, 4], t(mean[, -1]))
                                colnames(anno) <-
                                        c("mz", "time", as.character(t(mean[, 1])))
                                utils::write.csv(anno,
                                                 file = paste0(file, '.csv'),
                                                 row.names = F)
                                anno <- anno[index,]
                                return(anno)
                        } else {
                                return(report)
                        }
                }
        }


#' Get the report for samples with biological and technique replicates in different groups
#' @param xset the xcmsset object all of samples with technique replicates
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for the peaklist to MetaboAnalyst
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf Log intensity cutoff for peaks, default 5
#' @return dataframe with mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data in different groups if file are defaults NULL.
#' @export
getgrouprep <-
        function(xset,
                 file = NULL,
                 method = "medret",
                 intensity = "into",
                 rsdcf = 30,
                 inscf = 5) {
                data0 <- xcms::groupval(xset, method, intensity)
                data0[is.na(data0)] = 0
                data <- t(data0)
                lv <- xset@phenoData[, 1]
                lv2 <- xset@phenoData[, 2]
                mean <- stats::aggregate(data, list(lv, lv2), mean)
                sd <- stats::aggregate(data, list(lv, lv2), sd)
                suppressWarnings(rsd <- sd / mean * 100)
                result <-
                        cbind.data.frame(t(mean[, -c(1:2)]), t(sd[, -c(1:2)]), t(rsd[, -c(1:2)]))
                name <- paste0(t(mean[, 1]), t(mean[, 2]))
                colnames(result) <-
                        c(paste0(name, "mean"),
                          paste0(name, "sd"),
                          paste0(name, "rsd%"))

                meanB <-
                        stats::aggregate(mean[, -c(1:2)], list(as.vector(t(mean[, 1]))), mean)
                sdB <-
                        stats::aggregate(mean[, -c(1:2)], list(as.vector(t(mean[, 1]))), sd)
                suppressWarnings(rsdB <-
                                         sdB[,-1] / meanB[,-1] * 100)

                resultB <-
                        cbind.data.frame(t(meanB[,-1]), t(sdB[,-1]), t(rsdB))
                nameB <- as.vector(t(meanB[, 1]))
                colnames(resultB) <-
                        c(paste0(nameB, "mean"),
                          paste0(nameB, "sd"),
                          paste0(nameB, "rsd%"))

                datap <- xcms::groups(xset)
                report <- cbind.data.frame(datap, result, resultB)

                rsd0 <- as.matrix(rsdB)
                rsd0[is.nan(rsd0)] <- 0
                index <-
                        as.vector(apply(rsd0, 2, function(x)
                                all(x < rsdcf))) &
                        as.vector(apply(meanB[,-1], 2, function(x)
                                all(x > 10 ^ (inscf))))

                N <- sum(index)
                L <- length(index)

                if (N == 0) {
                        message(paste('No peaks found'))
                        return(NA)
                } else{
                        message(
                                paste(
                                        N,
                                        'out of',
                                        L,
                                        'peaks found with rsd cutoff',
                                        rsdcf,
                                        'and Log intensity cutoff',
                                        inscf
                                )
                        )

                        report <- report[index,]

                        if (!is.null(file)) {
                                result <- data.frame(t(mean[, -c(1:2)]))[index,]
                                data <-
                                        rbind(group = as.character(mean[, 1]), result)
                                utils::write.csv(data, file = paste0(file, '.csv'))
                                return(data)
                        } else {
                                return(report)
                        }
                }
        }

#' Get the time series or two factor DoE report for samples with biological and technique replicates in different groups
#' @param xset the xcmsset object all of samples with technique replicates in time series or two factor DoE
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for the peaklist to MetaboAnalyst
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf Log intensity cutoff for peaks, default 5
#' @return dataframe with time series or two factor DoE mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data in different groups if file are defaults NULL.
#' @export
gettimegrouprep <-
        function(xset,
                 file = NULL,
                 method = "medret",
                 intensity = "into",
                 rsdcf = 30,
                 inscf = 5) {
                data0 <- xcms::groupval(xset, method, intensity)
                data0[is.na(data0)] = 0
                data <- t(data0)
                lv <- xset@phenoData[, 1]
                lv2 <- xset@phenoData[, 2]
                lv3 <- xset@phenoData[, 3]
                mean <-
                        stats::aggregate(data, list(lv, lv2, lv3), mean)
                sd <- stats::aggregate(data, list(lv, lv2, lv3), sd)
                suppressWarnings(rsd <- sd / mean * 100)
                result <-
                        cbind.data.frame(t(mean[, -c(1:3)]), t(sd[, -c(1:3)]), t(rsd[, -c(1:3)]))
                name <-
                        paste0(t(mean[, 1]), t(mean[, 2]), t(mean[, 3]))
                colnames(result) <-
                        c(paste0(name, "mean"),
                          paste0(name, "sd"),
                          paste0(name, "rsd%"))

                meanB <-
                        stats::aggregate(mean[, -c(1:3)], list(as.vector(t(mean[, 1])), as.vector(t(mean[, 2]))), mean)
                sdB <-
                        stats::aggregate(mean[, -c(1:3)], list(as.vector(t(mean[, 1])), as.vector(t(mean[, 2]))), sd)
                suppressWarnings(rsdB <-
                                         sdB[,-c(1:2)] / meanB[,-c(1:2)] * 100)

                resultB <-
                        cbind.data.frame(t(meanB[,-c(1:2)]), t(sdB[,-c(1:2)]), t(rsdB))
                nameB <- paste0(t(meanB[, 1]), t(meanB[, 2]))
                colnames(resultB) <-
                        c(paste0(nameB, "mean"),
                          paste0(nameB, "sd"),
                          paste0(nameB, "rsd%"))

                datap <- xcms::groups(xset)
                report <- cbind.data.frame(datap, result, resultB)

                rsd0 <- as.matrix(rsdB)
                rsd0[is.nan(rsd0)] <- 0
                index <-
                        as.vector(apply(rsd0, 2, function(x)
                                all(x < rsdcf))) &
                        as.vector(apply(meanB[,-c(1, 2)], 2, function(x)
                                all(x > 10 ^ (inscf))))

                N <- sum(index)
                L <- length(index)

                if (N == 0) {
                        message(paste('No peaks found'))
                        return(NA)
                } else{
                        message(
                                paste(
                                        N,
                                        'out of',
                                        L,
                                        'peaks found with rsd cutoff',
                                        rsdcf,
                                        'and Log intensity cutoff',
                                        inscf
                                )
                        )

                        report <- report[index,]

                        if (!is.null(file)) {
                                result <- data.frame(t(meanB[, -c(1:2)]))[index,]
                                data <-
                                        rbind(
                                                time = as.character(meanB[, 1]),
                                                group = as.character(meanB[, 2]),
                                                result
                                        )
                                utils::write.csv(data, file = paste0(file, '.csv'))
                                return(data)
                        } else {
                                return(report)
                        }
                }
        }

#' Get the two factor DoE report for samples with biological replicates in different groups
#' @param xset the xcmsset object all of samples with technique replicates in time series or two factor DoE
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for the peaklist to MetaboAnalyst
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf Log intensity cutoff for peaks, default 5
#' @return dataframe with two factor DoE mean, standard deviation and RSD for those biological replicates combined with raw data in different groups if file are defaults NULL.
#' @export
getgroup2rep <-
        function(xset,
                 file = NULL,
                 method = "medret",
                 intensity = "into",
                 rsdcf = 30,
                 inscf = 5) {
                data0 <- xcms::groupval(xset, method, intensity)
                data0[is.na(data0)] = 0
                data <- t(data0)
                lv <- xset@phenoData[, 1]
                lv2 <- xset@phenoData[, 2]
                mean <- stats::aggregate(data, list(lv, lv2), mean)
                sd <- stats::aggregate(data, list(lv, lv2), sd)
                suppressWarnings(rsd <- sd / mean * 100)

                result <-
                        cbind.data.frame(t(mean[, -c(1:2)]), t(sd[, -c(1:2)]), t(rsd[, -c(1:2)]))
                name <- paste0(t(mean[, 1]), t(mean[, 2]))
                colnames(result) <-
                        c(paste0(name, "mean"),
                          paste0(name, "sd"),
                          paste0(name, "rsd%"))
                datap <- xcms::groups(xset)
                report <- cbind.data.frame(datap, result)
                rsd0 <- as.matrix(rsd[,-c(1,2)])
                rsd0[is.nan(rsd0)] <- 0
                index <-
                        as.vector(apply(rsd0, 2, function(x)
                                any(x < rsdcf))) &
                        as.vector(apply(mean[,-c(1,2)], 2, function(x)
                                any(x > 10 ^ (inscf))))

                N <- sum(index)
                L <- length(index)

                message(
                                paste(
                                        N,
                                        'out of',
                                        L,
                                        'peaks found with rsd cutoff',
                                        rsdcf,
                                        'and Log intensity cutoff',
                                        inscf
                                ))

                        report <- report[index,]

                        if (!is.null(file)) {
                                result <- data.frame(t(mean[, -c(1:2)]))[index,]
                                data <-
                                        rbind(group = as.character(mean[, 1]), result)
                                utils::write.csv(data, file = paste0(file, '.csv'))
                                return(data)
                        } else {
                                return(report)
                        }
        }
