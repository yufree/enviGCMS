#' Align two peaks vectors by mass to charge ratio and/or retention time
#' @param mz1 the mass to charge of reference peaks
#' @param mz2 the mass to charge of peaks to be aligned
#' @param ppm mass accuracy, default 10
#' @param deltart retention time shift table, default 10 seconds
#' @param rt1 retention time of reference peaks
#' @param rt2 retention time of peaks to be aligned
#' @return data frame with aligned peaks table
#' @examples
#' mz1 <- c(221.1171, 227.1390, 229.1546, 233.1497, 271.0790 )
#' mz2 <- c(282.279, 281.113, 227.139, 227.139, 302.207)
#' rt1 <- c(590.8710, 251.3820, 102.9230, 85.8850, 313.8240)
#' rt2 <- c(787.08, 160.02, 251.76, 251.76, 220.26)
#' getalign(mz1,mz2,rt1,rt2)
#' @export
getalign <-
        function(mz1,
                 mz2,
                 rt1 = NULL,
                 rt2 = NULL,
                 ppm = 10,
                 deltart = 10) {
                mza <-
                        data.table::as.data.table(cbind.data.frame(
                                mzmin = mz1 * (1 - ppm * 1e-6),
                                mzmax = mz1 * (1 + ppm * 1e-6)
                        ))
                mzb <-
                        data.table::as.data.table(cbind.data.frame(
                                mzmin = mz2 * (1 - ppm * 1e-6),
                                mzmax = mz2 * (1 + ppm * 1e-6)
                        ))
                colnames(mza) <- colnames(mzb) <- c("min", "max")
                data.table::setkey(mzb, min, max)
                overlapms <- data.table::foverlaps(mza, mzb, which = TRUE)

                overlapms$mz1 <- mz1[overlapms$xid]
                overlapms$mz2 <-
                        mz2[order(mz2, decreasing = FALSE)][overlapms$yid]
                if (!is.null(rt1) & !is.null(rt2)) {
                        overlapms$rt1 <- rt1[overlapms$xid]
                        overlapms$rt2 <-
                                rt2[order(mz2, decreasing = FALSE)][overlapms$yid]
                        drt <- abs(overlapms$rt1 - overlapms$rt2)
                        overlapms <- overlapms[drt < deltart, -c(2)]
                        re <-
                                data.frame(overlapms[!duplicated(overlapms) &
                                                             stats::complete.cases(overlapms), ])

                        if (nrow(re) > 0) {
                                return(re)
                        } else{
                                message('No result could be found.')
                        }

                } else{
                        message('No retention time information!')
                        over <- overlapms[, -2]
                        over2 <-
                                over[stats::complete.cases(over) & !duplicated(over), ]
                        re <- data.frame(over2)
                        if (nrow(re) > 0) {
                                return(re)
                        } else{
                                message('No result could be found.')
                        }
                }
        }
#' Align mass to charge ratio and/or retention time to remove redundancy
#' @param mz the mass to charge of reference peaks
#' @param rt retention time of reference peaks
#' @param ppm mass accuracy, default 10
#' @param deltart retention time shift table, default 10 seconds
#' @return index for
#' @examples
#' mz <- c(221.1171, 221.1170, 229.1546, 233.1497, 271.0790 )
#' rt <- c(590.8710, 587.3820, 102.9230, 85.8850, 313.8240)
#' getalign2(mz,rt)
#' @export
getalign2 <-
        function (mz,
                  rt,
                  ppm = 5,
                  deltart = 5) {
                yyy <- paste0(mz, rt)
                mzn <- mz[!duplicated(yyy)]
                rtn <- rt[!duplicated(yyy)]
                mzr <-
                        data.table::as.data.table(cbind.data.frame(
                                mzmin = mzn *         (1 - ppm * 1e-06),
                                mzmax = mzn * (1 + ppm * 1e-06),
                                mz = mzn,
                                rt = rtn
                        ))
                colnames(mzr) <- c("min", "max", "mz", "rt")
                data.table::setkey(mzr, min, max)
                overlapms <- data.table::foverlaps(mzr, mzr, which = TRUE)
                overlapms$mz1 <- mzr$mz[overlapms$xid]
                overlapms$rt1 <- mzr$rt[overlapms$xid]
                overlapms$mz2 <- mzr$mz[overlapms$yid]
                overlapms$rt2 <- mzr$rt[overlapms$yid]
                overlapms$drt <- abs(overlapms$rt1 - overlapms$rt2)
                overlapms <- overlapms[overlapms$drt < deltart,]
                overlap1 <- overlapms[overlapms$xid == overlapms$yid,]
                overlap2 <- overlapms[overlapms$xid != overlapms$yid,]
                overlap3 <- overlap1[!overlap1$xid %in% overlap2$xid,]
                x <- igraph::graph_from_data_frame(overlap2, directed = FALSE)
                y <-
                        names(igraph::components(x)$membership[!duplicated(igraph::components(x)$membership)])
                idx <- c(overlap3$xid, as.numeric(y))
                idxx <- idx[order(idx, decreasing = FALSE)]
                if (length(idxx) > 0) {
                        mzu <- mzr$mz[idxx]
                        rtu <- mzr$rt[idxx]
                        xxx <- paste0(mzu, rtu)
                        idxx <- which(yyy %in% xxx & !duplicated(yyy))
                        return(idxx)
                } else {
                        message("No result could be found.")
                }
        }

#' Get the overlap peaks by mass and retention time range
#' @param list1 list with data as peaks list, mz, rt, mzrange, rtrange and group information to be overlapped
#' @param list2 list with data as peaks list, mz, rt, mzrange, rtrange and group information to overlap
#' @return logical index for list 1's peaks
#' @export
#' @seealso \code{\link{getimputation}},\code{\link{getdoe}}
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
        if (length(index) > 0) {
                return(index)
        } else{
                message('No result could be found.')
        }
}

#' Annotation of MS1 data by compounds database by predefined paired mass distance
#' @param pmd adducts formula or paired mass distance for ions
#' @param mz unknown mass to charge ratios vector
#' @param ppm mass accuracy
#' @param db compounds database as dataframe. Two required columns are name and monoisotopic molecular weight with column names of name and mass
#' @return list or data frame
#' @export
getms1anno <- function(pmd, mz, ppm = 10, db = NULL) {
        # hr <- get(hr)
        if (is.character(pmd)) {
                pmds <- unlist(Map(enviGCMS::getmass, pmd))
        } else{
                pmds <- pmd
        }

        if (length(pmd) > 1) {
                pmdmt <- outer(db$mass, pmds, '+')
                rownames(pmdmt) <- db$name
                li <- list()
                for (i in seq_along(pmd)) {
                        re <- enviGCMS::getalign(pmdmt[, i], mz, ppm = ppm)
                        re2 <- db[re$xid, ]
                        re3 <- cbind.data.frame(re[, -1], re2)
                        colnames(re3)[1] <- pmd[i]
                        li[[i]] <- re3
                }
                return(li)
        } else{
                pmdmt <- outer(db$mass, pmds, '+')
                rownames(pmdmt) <- db$name
                re <- enviGCMS::getalign(pmdmt, mz, ppm = ppm)
                re2 <- db[re$xid, ]
                re3 <- cbind.data.frame(re[, -1], re2)
                colnames(re3)[1] <- pmd
                return(re3)
        }
}
