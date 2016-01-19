#' findline of the regression model
#' @export
findline <- function(data,threshold=2,temp = c(100,320)){
        y0 <- as.numeric(rownames(data))
        x <- as.numeric(colnames(data))/60
        group <- ifelse(log10(data)>threshold,1,0)
        diffmatrix <- apply(group,2,diff)
        difftemp <- apply(diffmatrix,2,which.min)
        y <- y0[difftemp]
        data <- data.frame(y,x)
        data <- data[data$y>101,]
        rangemz <- range(y0)
        rangert <- range(x)
        plot(data$y~data$x,
             xlab = 'Temperature(°C)',
             ylab = 'm/z',
             pch = 19,
             xlim = rangert,
             ylim = rangemz,
             xaxt = "n",
             yaxt = 'n',
             main = 'Regression Model for Signals',
             frame.plot = F)
        # display the temperature as x
        rtx <- seq(min(x),max(x),length.out=length(x))
        temp <- round(seq(temp[1],temp[2],length.out=length(x)))
        axis(1, at=x[temp%%20==0], labels= temp[temp%%20==0])
        # display the m/z as y
        mzy <- seq(min(y0),max(y0),length.out=length(y0))
        axis(2, at=mzy[y0%%100==0], labels= y0[y0%%100==0], las=2)
        abline(lm(data$y~data$x), col="red",lwd = 5)
        lines(lowess(data$y~data$x), col="blue",lwd = 5)
        slope <- (max(temp)-min(temp))/(max(data$x)-min(data$x))
        intercept <- min(data$x)
        data$x0 <- slope*data$x+intercept
        fit <- lm(data$y~data$x0)
        rmse <- round(sqrt(mean(resid(fit)^2)), 2)
        coefs <- coef(fit)
        b0 <- round(coefs[1],2)
        b1 <- round(coefs[2],2)
        r2 <- round(summary(fit)$r.squared, 2)
        pv <- round(anova(fit)$'Pr(>F)'[1], 2)
        eqn <- bquote(italic(Mass(m/z)) == .(b0) + .(b1)*"*"*italic(Temprature))
        eqn2 <- bquote(r^2 == .(r2) * "," ~~ p < 0.01)
        text(5, 950, adj=c(0,0), cex = 1,eqn)
        text(5, 900, adj=c(0,0), cex = 1,eqn2)
        legend("topright",c('OLS','LOWESS'),box.lty=0,pch = c(-1,-1),lty=c(1,1),lwd=c(2,2),col=c('red','blue'))
        return(fit)
}

#' substrict two GC-MS data
#' @export
comparems <- function(data1,data2,...){
        z <- data2-data1
        z[z<0] <- 0
        par(mfrow=c(2,3))

        plot(rowSums(data1),
             type = 'l',
             xlab = 'retention time',
             ylab = 'intensity',
             xaxt = 'n',
             frame.plot = F,
             main = '35ev')
        ind = as.numeric(rownames(data1))
        axis( 1, at=seq(0,1500,300 ), labels= round((ind[1]+(max(ind)-min(ind))/5*(0:5))/60,2))
        plot(rowSums(data2),
             type = 'l',
             xlab = 'retention time',
             ylab = 'intensity',
             xaxt = 'n',
             frame.plot = F,
             main = '70ev')
        ind = as.numeric(rownames(data2))
        axis( 1, at=seq(0,1500,300 ), labels= round((ind[1]+(max(ind)-min(ind))/5*(0:5))/60,2))
        plot(rowSums(z),
             type = 'l',
             xlab = 'retention time',
             ylab = 'intensity',
             xaxt = 'n',
             frame.plot = F,
             main = 'New TIC')
        ind = as.numeric(rownames(z))
        axis( 1, at=seq(0,1500,300 ), labels= round((ind[1]+(max(ind)-min(ind))/5*(0:5))/60,2))

        findline(data1)
        findline(data2)
        findline(z)
        return(z)
}

#' Plot the response group of GC-MS
#' @export
plotgroup <- function(data,threshold=2){
        .pardefault <- par(no.readonly = T)
        group <- ifelse(log10(data)>threshold,1,0)
        ind <- as.numeric(rownames(data))
        par(mfrow=c(1,2))
        image(t(group),
              xlab = 'retention time(min)',
              ylab = 'm/z',
              axes = F,
              col = heat.colors(2),
              useRaster = T)
        indmz <- as.numeric(rownames(data))
        indrt <- as.numeric(colnames(data))
        # display the RT as x
        rtx <- seq(0,1,length.out=length(indrt))
        axis(1, at=c(0,rtx[indrt%%300==0],1), labels= c('',indrt[indrt%%300==0]/60,''))
        # display the m/z as y
        mzy <- seq(0,1,length.out=length(indmz))
        axis(2, at=mzy[indmz%%100==0], labels= indmz[indmz%%100==0], las=2)

        hist(log10(data),
             breaks=100,
             main='',
             xlab='Distribution of Intensity(log base 10)')
        abline(v=threshold,lwd=5,lty=2,col='red')
        legend("topright","threshold",box.lty=0,pch = -1,lty=2,lwd=2,col='red')
        par(.pardefault)
}

#' Plot the backgrond of data
#' @export
plotsub <- function(data,...){
        datan <- apply(data,1,diff)
        datan[datan<0] <- 0
        datan <- t(datan)
        plotms(datan,...)
        findline(datan)
}
#' Plot the intensity distribution of GC-MS
#' @export
plotsms <- function(meanmatrix,sdmatrix,...){
        smoothScatter(log10(c(sdmatrix))~log10(c(meanmatrix)),
                      xlab = 'Mean(log scale based 10)',
                      ylab = 'Standard Deviation(log scale based 10)',
                      frame.plot = F,...)
        abline(a = 0, b = 1)
        abline(v=2)
}
#' @export
plothist <- function(data){
        data1=sample(data,100000)
        mixmdl = normalmixEM(log10(data1))
        plot(mixmdl,which=2,breaks=100,
             xlab2 = 'Intensity(log scale based 10)')
        lines(density(log10(data1)), lty=2, lwd=2)
        legend("topright",c('noise','signal','density'),box.lty=0,pch = c(-1,-1,-1),lty=c(1,1,2),lwd=c(2,2,2),col=c('red','green','black'))
}
