#' Short-Chain Chlorinated Paraffins(SCCPs) peaks infomation for quantitative analysis
#'
#' A dataset containing the ions, formula, Cl%, peak abundances for all SCCPs compounds
#' @docType data
#' @usage data(sccp)
#' @format A data frame with 24 rows and 8 variables:
#' \describe{
#'   \item{Cln}{Chlorine atom numbers}
#'   \item{Cn}{Carbon atom numbers}
#'   \item{formula}{molecular formula}
#'   \item{Hn}{hydrogen atom numbers}
#'   \item{ions}{[M-Cl]- ions}
#'   \item{mz}{m/z for the isotopologues with highest intensity}
#'   \item{intensity}{abundance of the isotopologues with highest intensity}
#'   \item{Clp}{Chlorine contents}
#' }
"sccp"

#' Get the peak information from SCCPs standards
#' @param data list from `xcmsRaw` function
#' @param ismz internal standards m/z
#' @param ppm resolution of mass spectrum
#' @param con concentration of standards
#' @param rt retention time range of sccps
#' @param rts retention time range of internal standards
#' @return list with peak information
#' @seealso \code{\link{getarea}},\code{\link{getsccp}}
#' @export

getareastd <- function(data = NULL, ismz = 323, ppm = 5, 
    con = 2000, rt = NULL, rts = NULL) {
    mz <- sccp$mz
    mzh <- mz + mz * 1e-06 * ppm
    mzl <- mz - mz * 1e-06 * ppm
    
    mzhis <- ismz + ismz * 1e-06 * ppm
    mzlis <- ismz + ismz * 1e-06 * ppm
    if (is.null(rts)) {
        eicis <- xcms::getEIC(data, mz = c(mzhis, mzlis))
        dfis <- eicis@eic$xcmsRaw[[1]]
        areais <- sum(diff(dfis[, 1]) * dfis[-1, 2])
    } else {
        eicis <- xcms::getEIC(data, mz = c(mzhis, mzlis))
        dfis <- eicis@eic$xcmsRaw[[1]]
        dfis <- dfis[dfis[, 1] > rts[1] & dfis[, 1] < rts[2], 
            ]
        areais <- sum(diff(dfis[, 1]) * dfis[-1, 2])
    }
    
    area <- vector()
    if (is.null(rt)) {
        for (i in 1:length(mz)) {
            eici <- xcms::getEIC(data, mz = c(mzh[i], mzl[i]))
            df <- eici@eic$xcmsRaw[[1]]
            area[i] <- sum(diff(df[, 1]) * df[-1, 2])
        }
    } else {
        for (i in 1:length(mz)) {
            eici <- xcms::getEIC(data, mz = c(mzh[i], mzl[i]))
            df <- eici@eic$xcmsRaw[[1]]
            df <- df[df[, 1] > rt[1] & df[, 1] < rt[2], 
                ]
            area[i] <- sum(diff(df[, 1]) * df[-1, 2])
        }
    }
    
    
    rarea <- area/(areais * sccp$Cln)
    rrares <- rarea/sum(rarea)
    pCl <- rrares * sccp$Clp
    
    sumpCl <- sum(pCl)
    sumrarea <- sum(rarea)/con
    
    ccomp <- stats::aggregate(rrares, by = list(sccp$Cn), 
        sum)
    colnames(ccomp) <- c("nC", "Formula group abundance")
    clcomp <- stats::aggregate(rrares, by = list(sccp$Cln), 
        sum)
    colnames(clcomp) <- c("nCl", "Formula group abundance")
    list <- list(sumpCl = sumpCl, sumrarea = sumrarea, 
        ccomp = ccomp, clcomp = clcomp)
    return(list)
}
#' Get the peak information from sampels for SCCPs detection
#' @param data list from `xcmsRaw` function
#' @param ismz internal standards m/z
#' @param ppm resolution of mass spectrum
#' @param rt retention time range of sccps
#' @param rts retention time range of internal standards
#' @return list with peak information
#' @seealso \code{\link{getareastd}},\code{\link{getsccp}}
#' @export
getarea <- function(data, ismz = 323, ppm = 5, rt = NULL, 
    rts = NULL) {
    mz <- sccp$mz
    mzh <- mz + mz * 1e-06 * ppm
    mzl <- mz - mz * 1e-06 * ppm
    
    mzhis <- ismz + ismz * 1e-06 * ppm
    mzlis <- ismz + ismz * 1e-06 * ppm
    if (is.null(rts)) {
        eicis <- xcms::getEIC(data, mz = c(mzhis, mzlis))
        dfis <- eicis@eic$xcmsRaw[[1]]
        areais <- sum(diff(dfis[, 1]) * dfis[-1, 2])
    } else {
        eicis <- xcms::getEIC(data, mz = c(mzhis, mzlis))
        dfis <- eicis@eic$xcmsRaw[[1]]
        dfis <- dfis[dfis[, 1] > rts[1] & dfis[, 1] < rts[2], 
            ]
        areais <- sum(diff(dfis[, 1]) * dfis[-1, 2])
    }
    
    area <- vector()
    if (is.null(rt)) {
        for (i in 1:length(mz)) {
            eici <- xcms::getEIC(data, mz = c(mzh[i], mzl[i]))
            df <- eici@eic$xcmsRaw[[1]]
            area[i] <- sum(diff(df[, 1]) * df[-1, 2])
        }
    } else {
        for (i in 1:length(mz)) {
            eici <- xcms::getEIC(data, mz = c(mzh[i], mzl[i]))
            df <- eici@eic$xcmsRaw[[1]]
            df <- df[df[, 1] > rt[1] & df[, 1] < rt[2], 
                ]
            area[i] <- sum(diff(df[, 1]) * df[-1, 2])
        }
    }
    
    
    rarea <- area/(areais * sccp$Cln)
    rrares <- rarea/sum(rarea)
    pCl <- rrares * sccp$Clp
    
    sumpCl <- sum(pCl)
    sumrarea <- sum(rarea)
    
    ccomp <- stats::aggregate(rrares, by = list(sccp$Cn), 
        sum)
    colnames(ccomp) <- c("nC", "Formula group abundance")
    clcomp <- stats::aggregate(rrares, by = list(sccp$Cln), 
        sum)
    colnames(clcomp) <- c("nCl", "Formula group abundance")
    list <- list(sumpCl = sumpCl, sumrarea = sumrarea, 
        ccomp = ccomp, clcomp = clcomp)
    return(list)
}
#' Quantitative analysis for short-chain chlorinated paraffins(SCCPs)
#' @param pathstds mzxml file path for SCCPs standards
#' @param pathsample mzxml file path for samples
#' @param ismz internal standards m/z
#' @param ppm resolution of mass spectrum
#' @param con concentration of standards
#' @param rt retention time range of sccps
#' @param rts retention time range of internal standards
#' @param log log transformation for response factor
#' @return list with peak information
#' @seealso \code{\link{getareastd}},\code{\link{getarea}}
#' @export
getsccp <- function(pathstds, pathsample, ismz = 323, ppm = 5, 
    con = 2000, rt = NULL, rts = NULL, log = T) {
    pathstd <- list.files(path = pathstds, full.names = T, 
        recursive = T)
    pathsamp <- list.files(path = pathsample, full.names = T, 
        recursive = T)
    nstd <- length(pathstd)
    nsamp <- length(pathsamp)
    # process SCCPs standards
    liststd <- list()
    for (i in 1:nstd) {
        file <- xcms::xcmsRaw(pathstd[i], profstep = 0.1)
        liststd[[i]] <- getareastd(file, ismz = ismz, ppm = ppm, 
            con = con, rt = rt, rts = rts)
    }
    pCl <- sapply(liststd, function(x) x$sumpCl)
    rarea <- sapply(liststd, function(x) x$sumrarea)
    
    # get the slope and intercept
    if (log) {
        lmfit <- stats::lm(log(rarea) ~ pCl)
    } else {
        lmfit <- stats::lm(rarea ~ pCl)
    }
    intercept <- lmfit$coefficients[1]
    slope <- lmfit$coefficients[2]
    
    # process SCCPs samples
    listsamp <- list()
    for (i in 1:nsamp) {
        file <- xcms::xcmsRaw(pathsamp[i], profstep = 0.1)
        listsamp[[i]] <- getarea(file, ismz = ismz, ppm = ppm, 
            rt = rt, rts = rts)
    }
    
    pCls <- sapply(listsamp, function(x) x$sumpCl)
    rareas <- sapply(listsamp, function(x) x$sumrarea)
    # get the concertration
    if (log) {
        rareasc <- exp(pCls * slope + intercept)
    } else {
        rareasc <- pCls * slope + intercept
    }
    
    con <- rareas/rareasc
    
    # extract the composition
    
    ccomp <- lapply(liststd, function(x) x$ccomp)
    clcomp <- lapply(liststd, function(x) x$clcomp)
    
    ccomps <- lapply(listsamp, function(x) x$ccomp)
    clcomps <- lapply(listsamp, function(x) x$clcomp)
    
    return(list(cons = con, Ccomp = ccomp, Clcomp = clcomp, 
        Ccomps = ccomps, Clcomps = clcomps, lmfit = lmfit))
}
