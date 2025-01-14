options(shiny.maxRequestSize=100*1024^2)
library(shiny)
library(xcms)
library(enviGCMS)
data("sccp")
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

getareastd <- function(data = NULL,
                       ismz = 323,
                       ppm = 5,
                       con = 2000,
                       rt = NULL,
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
                dfis <-
                        dfis[dfis[, 1] > rts[1] & dfis[, 1] < rts[2], ]
                areais <- sum(diff(dfis[, 1]) * dfis[-1, 2])
        }

        area <- vector()
        if (is.null(rt)) {
                for (i in seq_along(mz)) {
                        eici <- xcms::getEIC(data, mz = c(mzh[i], mzl[i]))
                        df <- eici@eic$xcmsRaw[[1]]
                        area[i] <- sum(diff(df[, 1]) * df[-1, 2])
                }
        } else {
                for (i in seq_along(mz)) {
                        eici <- xcms::getEIC(data, mz = c(mzh[i], mzl[i]))
                        df <- eici@eic$xcmsRaw[[1]]
                        df <- df[df[, 1] > rt[1] & df[, 1] < rt[2], ]
                        area[i] <- sum(diff(df[, 1]) * df[-1, 2])
                }
        }


        rarea <- area / (areais * sccp$Cln)
        rrares <- rarea / sum(rarea)
        pCl <- rrares * sccp$Clp

        sumpCl <- sum(pCl)
        sumrarea <- sum(rarea) / con

        ccomp <- stats::aggregate(rrares, by = list(sccp$Cn),
                                  sum)
        colnames(ccomp) <- c("nC", "Formula group abundance")
        clcomp <- stats::aggregate(rrares, by = list(sccp$Cln),
                                   sum)
        colnames(clcomp) <- c("nCl", "Formula group abundance")
        list <- list(
                sumpCl = sumpCl,
                sumrarea = sumrarea,
                ccomp = ccomp,
                clcomp = clcomp
        )
        return(list)
}
#' Get the peak information from samples for SCCPs detection
#' @param data list from `xcmsRaw` function
#' @param ismz internal standards m/z
#' @param ppm resolution of mass spectrum
#' @param rt retention time range of sccps
#' @param rts retention time range of internal standards
#' @return list with peak information
#' @seealso \code{\link{getareastd}},\code{\link{getsccp}}
#' @export
getarea <- function(data,
                    ismz = 323,
                    ppm = 5,
                    rt = NULL,
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
                dfis <-
                        dfis[dfis[, 1] > rts[1] & dfis[, 1] < rts[2], ]
                areais <- sum(diff(dfis[, 1]) * dfis[-1, 2])
        }

        area <- vector()
        if (is.null(rt)) {
                for (i in seq_along(mz)) {
                        eici <- xcms::getEIC(data, mz = c(mzh[i], mzl[i]))
                        df <- eici@eic$xcmsRaw[[1]]
                        area[i] <- sum(diff(df[, 1]) * df[-1, 2])
                }
        } else {
                for (i in seq_along(mz)) {
                        eici <- xcms::getEIC(data, mz = c(mzh[i], mzl[i]))
                        df <- eici@eic$xcmsRaw[[1]]
                        df <- df[df[, 1] > rt[1] & df[, 1] < rt[2], ]
                        area[i] <- sum(diff(df[, 1]) * df[-1, 2])
                }
        }


        rarea <- area / (areais * sccp$Cln)
        rrares <- rarea / sum(rarea)
        pCl <- rrares * sccp$Clp

        sumpCl <- sum(pCl)
        sumrarea <- sum(rarea)

        ccomp <- stats::aggregate(rrares, by = list(sccp$Cn),
                                  sum)
        colnames(ccomp) <- c("nC", "Formula group abundance")
        clcomp <- stats::aggregate(rrares, by = list(sccp$Cln),
                                   sum)
        colnames(clcomp) <- c("nCl", "Formula group abundance")
        list <- list(
                sumpCl = sumpCl,
                sumrarea = sumrarea,
                ccomp = ccomp,
                clcomp = clcomp
        )
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
getsccp <- function(pathstds,
                    pathsample,
                    ismz = 323,
                    ppm = 5,
                    con = 2000,
                    rt = NULL,
                    rts = NULL,
                    log = TRUE) {
        pathstd <- list.files(path = pathstds,
                              full.names = TRUE,
                              recursive = TRUE)
        pathsamp <- list.files(path = pathsample,
                               full.names = TRUE,
                               recursive = TRUE)
        nstd <- length(pathstd)
        nsamp <- length(pathsamp)
        # process SCCPs standards
        liststd <- list()
        for (i in 1:nstd) {
                file <- xcms::xcmsRaw(pathstd[i], profstep = 0.1)
                liststd[[i]] <- getareastd(
                        file,
                        ismz = ismz,
                        ppm = ppm,
                        con = con,
                        rt = rt,
                        rts = rts
                )
        }
        pCl <- vapply(liststd, function(x)
                x$sumpCl, 1)
        rarea <- vapply(liststd, function(x)
                x$sumrarea, 1)

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
                listsamp[[i]] <- getarea(
                        file,
                        ismz = ismz,
                        ppm = ppm,
                        rt = rt,
                        rts = rts
                )
        }

        pCls <- vapply(listsamp, function(x)
                x$sumpCl, 1)
        rareas <- vapply(listsamp, function(x)
                x$sumrarea, 1)
        # get the concentration
        if (log) {
                rareasc <- exp(pCls * slope + intercept)
        } else {
                rareasc <- pCls * slope + intercept
        }

        con <- rareas / rareasc

        # extract the composition

        ccomp <- lapply(liststd, function(x)
                x$ccomp)
        clcomp <- lapply(liststd, function(x)
                x$clcomp)

        ccomps <- lapply(listsamp, function(x)
                x$ccomp)
        clcomps <- lapply(listsamp, function(x)
                x$clcomp)

        return(
                list(
                        cons = con,
                        Ccomp = ccomp,
                        Clcomp = clcomp,
                        Ccomps = ccomps,
                        Clcomps = clcomps,
                        lmfit = lmfit
                )
        )
}

shinyServer(function(input, output) {
        liststd <- eventReactive(input$go, {
                if (!is.null(input$Standards)){
                        n <- length(input$Standards$datapath)
                        list <- list()
                        for(i in 1:n){

                                file <- xcmsRaw(input$Standards$datapath[i],profstep = 0.1)
                                list[[i]] <- getareastd(file,ismz = input$ISmz,ppm = input$ppm, con = input$con, rt = input$SCCPrt, rts = input$ISrt)

                        }}
                list
        })

        output$reg <- renderTable({
                li <- liststd()
                pCl <- vapply(li,function(x) x$sumpCl,1)
                rarea <- vapply(li,function(x) x$sumrarea,1)
                t <- summary(lm(rarea~pCl))
                t$coefficients
        })

        output$reg2 <- renderTable({
                li <- liststd()
                pCl <- vapply(li,function(x) x$sumpCl,1)
                rarea <- vapply(li,function(x) x$sumrarea,,1)
                t <- summary(lm(rarea~pCl))
                t$coefficients
        })

        output$plotstd <- renderPlot({
                li <- liststd()
                pCl <- vapply(li,function(x) x$sumpCl,1)
                rarea <- vapply(li,function(x) x$sumrarea,1)
                plot(rarea~pCl,xlab = 'Chlorine content %', ylab = 'Response Factor', pch = 19)
        })

        output$plotcomp <- renderPlot({
                li <- liststd()
                ccomp <- lapply(li,function(x) x$ccomp)
                clcomp <- lapply(li,function(x) x$clcomp)
                par(mfrow = c(length(ccomp),2),mar=c(1,1,1,1))

                for(i in seq_along(ccomp)){
                        ccompi <- ccomp[[i]]
                        clcompi <- clcomp[[i]]
                        barplot(ccompi[,2],names.arg = ccompi[,1],main = paste0('Standard',i,"'s C Composition"))
                        barplot(clcompi[,2],names.arg = clcompi[,1], main = paste0('Standard',i,"'s Cl Composition"))
                }
        })

        listsample <- eventReactive(input$go2, {
                if (!is.null(input$Samples)){
                        n <- length(input$Samples$datapath)
                        list <- list()
                        for(i in 1:n){

                                file <- xcmsRaw(input$Samples$datapath[i],profstep = 0.1)
                                list[[i]] <- getarea(file,ismz = input$ISmz,ppm = input$ppm, rt = input$SCCPrt, rts = input$ISrt)

                        }}
                list
        })

        output$plotcomps <- renderPlot({
                li <- listsample()
                ccomp <- lapply(li,function(x) x$ccomp)
                clcomp <- lapply(li,function(x) x$clcomp)
                par(mfrow = c(length(ccomp),2),mar=c(1,1,1,1))

                for(i in seq_along(ccomp)){
                        ccompi <- ccomp[[i]]
                        clcompi <- clcomp[[i]]
                        barplot(ccompi[,2],names.arg = ccompi[,1],main = paste0('Sample',i,"'s C Composition"))
                        barplot(clcompi[,2],names.arg = clcompi[,1], main = paste0('Sample',i,"'s Cl Composition"))
                }
        })

        output$results <- renderPrint({
                li <- listsample()
                pCl <- vapply(li,function(x) x$sumpCl,1)
                rarea <- vapply(li,function(x) x$sumrarea,1)
                if(input$log){
                        cons <- rarea/exp((pCl*input$slope+input$inc))
                }else{
                        cons <- rarea/(pCl*input$slope+input$inc)
                }

                round(cons,2)
        })

        output$data <- renderDataTable({
                sccp
        })

})
