#' GetIntegration was mainly used for get the intergration of certain ion's chromatogram data and plot the data
#' @param data file should be a dataframe with the first column RT and second column intensity of the SIM ions.
#' @param minrt a rough smallest RT range contained only one peak to get the area
#' @param maxrt a rough largest RT range contained only one peak to get the area
#' @param n points in the moving average smooth box, default value is 5
#' @param m numbers of points for regression to get the slope
#' @param startslope the threshold value for start peak (%max slope)
#' @param stopslope the threshold value for stop peak (%max slope)
#' @param baseline numbers of the points for the baseline of the signal
#' @param noslope logical, if using a horizon line to get area or not
#' @param smoothit logical, if using an average smooth box or not. If using, n will be used
#' @param half logical, if using the left half peak to caculate the area
#' @return list intergration data such as peak area, peak hight, signal and the slope data.
#' @examples
#' NULL
#
GetIntegration <- function(data, minrt = 8.3 , maxrt = 9, n=5, m=5, startslope = 2, stopslope = 2, baseline = 20, noslope = T, smoothit = T,name = NULL,half = F,...){
        # subset the data
        subdata <- data[data[,1] > minrt&data[,1] < maxrt,]
        # get the signal and the RT
        RTrange <- subdata[,1]
        signal <- subdata[,2]
        # data smooth
        if(smoothit) {
                ma <- function(x,n){filter(x,rep(1/n,n),circular = T)}
                signal <- ma(signal,n)
        }
        # get the slope data
        RTrangemsec <- RTrange*60*1e3
        slopedata <- signal
        back <- round(m/2)
        forth <- m-back
        if(m>2){
                for(i in (back+1):(length(signal)-forth)){
                        slopedata[i] <- coef(lm(signal[(i-back+1):(i+forth)] ~ RTrangemsec[(i-back+1):(i+forth)]))[2]
                }
                slopedata[1:back] <- slopedata[back+1] # first few points
                slopedata[(length(signal)-forth-1):length(signal)] <- slopedata[(length(signal)-forth)]# last few points
        }else{
                # for n = 2 points; without linear regression (much faster)
                # time per scan in millisec
                delta_t <- (t[length(RTrangemsec)]-RTrangemsec[1]/(length(RTrangemsec)-1))
                for(i in 2:length(signal)) {
                        slopedata[i] <- (signal[i]-signal[i-1])/delta_t
                }
                slopedata[1] <- 0
        }
        # search for peak start
        i <- baseline
        while((slopedata[i]<=(startslope/100*max(slopedata)))&(i<length(signal))) i <- i+1
        rtstart <- RTrange[i] # peak start found
        scanstart <- i # (slope > threshold)
        sigstart <- mean(signal[(i-baseline+1):i]) # baseline intensity found
        # search for peak top
        i <- which.max(slopedata) # jump to slope max.
        while((slopedata[i]>=0)&(i<(length(signal)-baseline))) i <- i+1
        rtpeak <- RTrange[i] # peak top found
        scanpeak <- i # (slope = 0)
        sigpeak <- signal[i]
        # search for peak end
        i <- which.min(slopedata) # jump to slope min.
        while((slopedata[i]<= -(stopslope/100*max(slopedata))) & (i<(length(signal)-baseline))) i <- i+1
        rtend <- RTrange[i] # peak end found
        scanend <- i # (-slope < threshold)
        sigend <- mean(signal[i:(i+ baseline-1)])
        # if background without slope
        if(!noslope) sigend <- sigstart
        # subtract background from signal
        background <- signal
        for(i in scanstart:scanend) { # get background
                background[i] <- sigstart+(sigend-sigstart)/(scanend-scanstart)*(i-scanstart)
        }
        subsignal <- signal-background
        # get the length of the signal
        lengthsig <- length(scanstart:scanend)
        # calculate area; using a Riemann integral (dimension: intensity x min)
        area <- 0
        scantime <- (RTrange[scanend]-RTrange[scanstart])/(scanend-scanstart)*60 # time per scan in second
        # when half peak
        if(half==T) {
                for(i in scanstart:scanpeak) area <- area + subsignal[i]*scantime
        }else{
                for(i in scanstart:scanend) area <- area + subsignal[i]*scantime
        }
        bgstart <- RTrange[scanstart-baseline+1]
        bgend <- RTrange[scanend+baseline-1]
        # calculate height
        sigpeakbase <- sigstart+(sigend-sigstart)/(scanend-scanstart)*(scanpeak-scanstart)
        height <- sigpeak - sigpeakbase
        # collect the data for plot peak and slope
        peakdata <- c(baseline,rtstart,rtend,rtpeak,scanstart,scanend,scanpeak,sigstart,sigend,sigpeak,sigpeakbase,lengthsig)
        names(peakdata) <- c('baseline','peak start RT','peak end RT','peak RT','baseline start RT ID','baseline end RT ID','baseline peak RT ID','start RT intensity','end RT intensity','peak RT intensity','peak baseline','points')
        # return the result as a list
        list <- list(area=area,height=height,peakdata=peakdata,RTrange=RTrange,signal=signal,slopedata=slopedata)
        return(list)
}
#' Get the MIR and related information from the files
batch <- function(file,mz1,mz2,...){
        data1 <- xcmsRaw(file)
        df <- data1@env$profile
        rt <- data1@scantime/60
        rownames(df) <- seq(data1@mzrange[1],data1@mzrange[2])
        xa <- df[as.character(mz1),]
        xa <- data.frame(rt,xa)
        xb <- df[as.character(mz2),]
        xb <- data.frame(rt,xb)
        name1 <- paste('m/z:',mz1)
        name2 <- paste('m/z:',mz2)
        xl <- GetIntegration(xa,name1,...)
        xh <- GetIntegration(xb,name2,...)
        par(mfrow = c(2,2),mar = c(2,2,2,2),oma=c(0,0,0,0))
        plotint(xl,name1)
        plotintslope(xl,name1)
        plotint(xh,name2)
        plotintslope(xh,name2)
        list <- list(xl,xh)
        area <- sapply(list,function(x) x$area)
        height <- sapply(list,function(x) x$height)
        arearatio <- area[1]/area[2]
        heightratio <- height[1]/height[2]
        points <- round(mean(sapply(list,function(x) x$peakdata[12])))
        ratio <- c(arearatio,heightratio,points)
        names(ratio) <- c('area ratio', 'height ratio','points')
        return(ratio)
}
#' Just intergrate data according to fixed rt and fixed noise area
#' @param data
#' @param minrt
#' @param maxrt
#' @param minbrt
#' @param maxbrt
#' @param smoothit
#' @return area intergration data
Integration <- function(data, minrt = 8.3 , maxrt = 9,minbrt = 8.3, maxbrt = 8.4,smoothit = T){
        # subset the data
        subdata <- data[data[,1] > minrt&data[,1] < maxrt,]
        # get the signal and the RT
        RTrange <- subdata[,1]
        signal <- subdata[,2]
        # subset the noise
        subnoise <- data[data[,1]>minbrt&data[,1]<maxbrt,]
        # get the noise and the RT
        RTrange2 <- subnoise[,1]
        noise <- subnoise[,2]
        # data smooth
        if(smoothit) {
                ma <- function(x,n){filter(x,rep(1/n,n),circular = T)}
                signal <- ma(signal,n=5)
                noise <- ma(noise,n=5)
        }
        baseline <- mean(noise)
        signaln <- signal-baseline
        area <- 0
        # calculate area; using a Riemann integral (dimension: intensity x delta time)
        for(i in 1:(length(RTrange)-1)) {
                area <- area + signaln[i]*(RTrange[i+1]-RTrange[i])*60
        }
        return(area)
}
#' Get the MIR from the file
qbatch <- function(file,mz1,mz2,minrt = 8.65,maxrt = 8.74,minbrt = 8.85, maxbrt = 8.87){
        data1 <- xcmsRaw(file)
        df <- data1@env$profile
        rt <- data1@scantime/60
        rownames(df) <- seq(data1@mzrange[1],data1@mzrange[2])
        xa <- df[as.character(mz1),]
        xa <- data.frame(rt,xa)
        xb <- df[as.character(mz2),]
        xb <- data.frame(rt,xb)
        xl <- Integration(xa,minrt,maxrt,minbrt,maxbrt)
        xh <- Integration(xb,minrt,maxrt,minbrt,maxbrt)
        arearatio <- xl/xh
        return(arearatio)
}

#' Get the selected isotopologues at certain MS data
Getisotopologues <- function(formula = 'C12OH6Br4', charge = '1',width = 0.3){
        # input the forlmula and charge for your molecular, this demo was for BDE-47
        formula <- get.formula(formula, charge)
        # get the isotopes pattern of your molecular with high abandances. Here
        # we suggest more than 10% abundance of your base peak would meet the SNR
        isotopes <- data.frame(get.isotopes.pattern(formula,minAbund=0.1))
        # order the intensity by the abandance
        findpairs <- isotopes[order(isotopes[,2],decreasing = T),]
        # find the most similar pairs with high abandance
        df <- outer(findpairs[,2],findpairs[,2],'/')
        rownames(df) <- colnames(df) <- findpairs$mass
        diag(df) <- df[upper.tri(df)] <- 0
        t <- which(df==max(df),arr.ind=T)
        isotopologues1 <- as.numeric(rownames(df)[t[1]])
        isotopologues2 <- as.numeric(colnames(df)[t[2]])
        isotopologuesL <- min(isotopologues1,isotopologues2)
        isotopologuesH <- max(isotopologues1,isotopologues2)
        # get the caculated ratio at certain resolution
        isotopes2 <- get.isotopes.pattern(formula,minAbund=0.00000001)
        ratio <- sum(isotopes2[isotopes2[,1]>isotopologuesL-width & isotopes2[,1]<isotopologuesL+width,2])/sum(isotopes2[isotopes2[,1]>isotopologuesH-width & isotopes2[,1]<isotopologuesH+width,2])
        peak <- c(round(isotopologuesL,digits = 1),round(isotopologuesH,digits = 1),round(ratio,digits = 5))
        # peak <- as.character(peak)
        names(peak) <- c('light isotopologue','high isotopologue','caculated ratio')
        return(data.frame(peak))
}
