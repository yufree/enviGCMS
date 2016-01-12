#' function for molecular isotope ratio
#' load the rcdk package
library(rcdk)

#' source("http://bioconductor.org/biocLite.R")
#' biocLite("xcms")
library(xcms)
#' Distiguished the distribution by MLE
library(mixtools)

#' filter data by average moving box
#'
#' @param x a vector
#' @param n A number to indentify the size of the moving box.
#' @return The filtered data
#' @examples
#' ma(c(1:10), 2)
ma <- function(x,n){
        filter(x,rep(1/n,n),circular = T)
}

#' plot GC-MS data as a heatmap
#'
#' @param data imported data matrix of GC-MS
#' @param col custumized color
#' @return heatmap
plotms <- function(data,col = heat.colors(108),...){
        .pardefault <- par(no.readonly = T)
        layout(matrix(c(1,2), nrow=2), heights = c(1,5))
        # get the mz and rt range and rotate the matrix to adapt the image function
        indmz <- as.numeric(rownames(data))
        indrt <- as.numeric(colnames(data))
        z <- log10(t(data+1))
        # show the intensity scale in log 10 based scale
        par(mar=c(2,4,1,2))
        zlim <- range(z)
        breaks <- seq(zlim[1], zlim[2], round((zlim[2]-zlim[1])/7))
        poly <- vector(mode="list", length(col))
        plot(1,1,
             t="n",
             ylim=c(0,1), xlim=range(z),
             xaxt='n', yaxt='n',
             xaxs="i", yaxs="i",
             ylab = '',xlab = '',...)
        mtext('Intensity',side = 2,line = 0.5,las = 1)
        axis(1,at=breaks,labels = 10^(breaks),las=1)
        bks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
        for(i in seq(poly)){
                polygon(c(bks[i], bks[i+1], bks[i+1], bks[i]), c(0,0,1,1), col=col[i], border=NA)
        }
        # show the heatmap
        par(mar=c(4,4,1,2))
        image(z,
              xlab = 'retention time(min)',
              ylab = 'm/z',
              axes = F,
              col = col,
              useRaster = T)
        # display the RT as x
        rtx <- seq(0,1,length.out=length(indrt))
        axis(1, at=c(0,rtx[indrt%%300==0],1), labels= c('',indrt[indrt%%300==0]/60,''))
        # display the m/z as y
        mzy <- seq(0,1,length.out=length(indmz))
        axis(2, at=mzy[indmz%%100==0], labels= indmz[indmz%%100==0], las=2)
        par(.pardefault)
}

#' plot GC-MS data as a heatmap for constant speed of temperature rising
#' @param data imported data matrix of GC-MS
#' @param col custumized color
#' @param temp temprature range for constant speed
#' @return heatmap
plott <- function(data,col = heat.colors(108),temp = c(100,320),...){
        .pardefault <- par(no.readonly = T)
        layout(matrix(c(1,2), nrow=2), heights = c(1,5))
        # get the mz and rt range and rotate the matrix to adapt the image function
        indmz <- as.numeric(rownames(data))
        indrt <- as.numeric(colnames(data))
        z <- log10(t(data+1))
        # show the intensity scale in log 10 based scale
        par(mar=c(2,4,1,2))
        zlim <- range(z)
        breaks <- seq(zlim[1], zlim[2], round((zlim[2]-zlim[1])/7))
        poly <- vector(mode="list", length(col))
        plot(1,1,
             t="n",
             ylim=c(0,1), xlim=range(z),
             xaxt='n', yaxt='n',
             xaxs="i", yaxs="i",
             ylab = '',xlab = '',...)
        mtext('Intensity',side = 2,line = 0.5,las = 1)
        axis(1,at=breaks,labels = 10^(breaks),las=1)
        bks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
        for(i in seq(poly)){
                polygon(c(bks[i], bks[i+1], bks[i+1], bks[i]), c(0,0,1,1), col=col[i], border=NA)
        }
        # show the heatmap
        par(mar=c(4,4,1,2))
        image(z,
              xlab = 'Temperature(°C)',
              ylab = 'm/z',
              axes = F,
              col = col,
              useRaster = T)
        # display the temperature as x
        rtx <- seq(0,1,length.out=length(indrt))
        temp <- round(seq(temp[1],temp[2],length.out=length(indrt)),0)
        axis(1, at=rtx[temp%%20==0], labels= temp[temp%%20==0])
        # display the m/z as y
        mzy <- seq(0,1,length.out=length(indmz))
        axis(2, at=mzy[indmz%%100==0], labels= indmz[indmz%%100==0], las=2)
        par(.pardefault)
}

#' Plot mass spectrum of certain retention time and return mass spectrum vector (MSP file) for NIST search
#' @param data imported data matrix of GC-MS
#' @param rt vector range of the retention time
#' @param ms vector range of the m/z
#' @return plot and MSP files for NIST search
plotrtms <- function(data,rt,ms){
        data <- getsubmd(data,rt,ms)
        temp <- apply(data,1,mean)
        plot(temp,
             type = 'h',
             ylab = 'intensity',
             xlab = 'm/z',
             xaxt = "n",
             frame.plot = F
        )
        mz <- as.numeric(names(temp))
        index <- seq(1,length(temp))
        axis(1, at=index[mz%%100==0], labels= mz[mz%%100==0])
        writeMSP(temp)
}

#' Plot EIC of certain m/z and return dataframe for intergration
#' @param data imported data matrix of GC-MS
#' @param ms m/z to be extracted
#' @param rt vector range of the retention time
#' @param n logical smooth or not
#' @return dataframe with  with the first column RT and second column intensity of the SIM ions.
plotmsrt <- function(data,ms,rt=c(3.1,25),n=F){
        data <- getsubmd(data,rt,c(ms,ms+1))[1,]
        if(n){
                data <- ma(data,n)
        }
        x <- as.numeric(names(data))/60
        plot(data~x,
             type = 'l',
             main = bquote('m/z = '~.(ms)),
             ylab = 'intensity',
             xlab = 'retention time(min)',
             frame.plot = F
        )
        data <- data.frame(x,data)
        colnames(data) <- c('RT','Intensity')
        return(data)
}

#' Plot TIC
#' @param data imported data matrix of GC-MS
#' @param ms vector range of the m/z
#' @param rt vector range of the retention time
#' @param n logical smooth or not
#' @return plot
plottic <- function(data,rt=c(3.1,25),ms=c(100,1000),n=F){
        data <- getsubmd(data,rt,ms)
        data <- apply(data,2,sum)
        if(n){
                data <- ma(data,n)
        }
        x <- as.numeric(names(data))/60
        plot(data~x,
             type = 'l',
             ylab = 'intensity',
             xlab = 'retention time(min)',
             frame.plot = F
        )
}
#' plot the information of intergretion
#' @param list list from getinteragtion
#' @param name the title of the plot
#' @return NULL
plotint <- function(list,name=NULL){
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
        plot(RTrange,signal,
             xlab="time (min)",ylab="intensity","l",
             ylim = c(-0.02*max(signal),1.02*max(signal)),
             main= paste(name,'Peak'))
        lines(c(rtstart, rtend),c(sigstart, sigend),"l", col = "red")
        lines(c(RTrange[scanstart-baseline+1], rtstart), c(sigstart, sigstart),"l", col = "darkgreen")
        lines(c(rtend,RTrange[scanend+baseline-1]), c(sigend, sigend),"l", col = "darkgreen")
        lines(c(rtstart,rtstart),c(0.8*sigstart, sigstart*1.2), "l", col="blue")
        lines(c(rtend,rtend), c(0.8*sigend, sigend*1.2), "l", col="blue")
        lines(c(rtpeak,rtpeak), c(sigpeak, sigpeakbase), "l", col="blue")

        # print RT, heights & areas
        text(rtstart,sigpeak*.2,as.character(round(rtstart,3)),col="darkgreen")
        text(rtend,sigpeak*.3,as.character(round(rtend,3)),col="darkgreen")
        text(rtpeak-0.1,0.9*sigpeak,paste('RT:',as.character(round(rtpeak,3))),col="darkgreen")
        text(rtpeak+0.1,sigpeak*.7,paste("area:", format(area,digits=2)),col="red")
        text(rtpeak+0.1,sigpeak*.9,paste("height:",format(sigpeak, digits=2)),col="red")
}
#' plot the information of intergretion
#' @param list list from getinteragtion
#' @param name the title of the plot
#' @return NULL
plotintslope <- function(list,name=NULL){
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
        plot(RTrange,slopedata,xlab="time (min)",ylab="slope","l",main=paste(name,'Slope'))
        lines(c(rtstart,rtstart),c(-.1*max(slopedata), .1*max(slopedata)), "l", col="blue")
        lines(c(rtend,rtend), c(-.1*max(slopedata), .1*max(slopedata)), "l", col="blue")
        lines(c(rtpeak,rtpeak), c(-.5*max(slopedata), .5*max(slopedata)), "l", col="blue")
}
