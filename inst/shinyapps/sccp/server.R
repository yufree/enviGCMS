options(shiny.maxRequestSize=100*1024^2)
library(shiny)
library(xcms)
library(enviGCMS)
data("sccp")

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
                pCl <- sapply(li,function(x) x$sumpCl)
                rarea <- sapply(li,function(x) x$sumrarea)
                t <- summary(lm(rarea~pCl))
                t$coefficients
        })

        output$reg2 <- renderTable({
                li <- liststd()
                pCl <- sapply(li,function(x) x$sumpCl)
                rarea <- sapply(li,function(x) x$sumrarea)
                t <- summary(lm(rarea~pCl))
                t$coefficients
        })

        output$plotstd <- renderPlot({
                li <- liststd()
                pCl <- sapply(li,function(x) x$sumpCl)
                rarea <- sapply(li,function(x) x$sumrarea)
                plot(rarea~pCl,xlab = 'Chlorine content %', ylab = 'Response Factor', pch = 19)
        })

        output$plotcomp <- renderPlot({
                li <- liststd()
                ccomp <- lapply(li,function(x) x$ccomp)
                clcomp <- lapply(li,function(x) x$clcomp)
                par(mfrow = c(length(ccomp),2),mar=c(1,1,1,1))

                for(i in 1:length(ccomp)){
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

                for(i in 1:length(ccomp)){
                        ccompi <- ccomp[[i]]
                        clcompi <- clcomp[[i]]
                        barplot(ccompi[,2],names.arg = ccompi[,1],main = paste0('Sample',i,"'s C Composition"))
                        barplot(clcompi[,2],names.arg = clcompi[,1], main = paste0('Sample',i,"'s Cl Composition"))
                }
        })

        output$results <- renderPrint({
                li <- listsample()
                pCl <- sapply(li,function(x) x$sumpCl)
                rarea <- sapply(li,function(x) x$sumrarea)
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
