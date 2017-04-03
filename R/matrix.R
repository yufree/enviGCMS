#' find line of the regression model
#' @param data imported data matrix of GC-MS
#' @param threshold the threshold of the response (log based 10)
#' @param temp the scale of the oven temprature(constant rate)
#' @return list linear regression model for the matrix
#' @export
findline <- function(data,threshold=2,temp = c(100,320)){
        y0 <- as.numeric(rownames(data))
        x <- as.numeric(colnames(data))/60
        # get the group
        group <- ifelse(log10(data)>threshold,1,0)
        # get the difference matrix
        diffmatrix <- apply(group,2,diff)
        # get the points with the smallest differences at the smallest m/z
        difftemp <- apply(diffmatrix,2,which.min)
        y <- y0[difftemp]
        data <- data.frame(y,x)
        # remove the meaningless bottom
        data <- data[data$y>101,]
        rangemz <- range(y0)
        rangert <- range(x)
        plot(data$y~data$x,
             xlab = 'Temperature(\u00b0C)',
             ylab = 'm/z',
             pch = 19,
             xlim = rangert,
             ylim = rangemz,
             xaxt = "n",
             yaxt = 'n',
             main = '',
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
        intercept <- min(temp)
        data$x0 <- slope*(data$x-min(data$x))+intercept
        fit <- lm(data$y~data$x0)
        rmse <- round(sqrt(mean(resid(fit)^2)), 2)
        coefs <- coef(fit)
        b0 <- round(coefs[1],2)
        b1 <- round(coefs[2],2)
        r2 <- round(summary(fit)$r.squared, 2)
        pv <- anova(fit)$'Pr(>F)'[1]
        eqn <- bquote(italic(m/z) == .(b0) + .(b1)*"*"*italic(Temprature))
        eqn2 <- bquote(r^2 == .(r2) * "," ~~ p == .(pv))
        text(5, 950, adj=c(0,0), cex = 1,eqn)
        text(5, 900, adj=c(0,0), cex = 1,eqn2)
        legend("topright",c('OLS','LOWESS'),box.lty=0,pch = c(-1,-1),lty=c(1,1),lwd=c(2,2),col=c('red','blue'))
        return(fit)
}

#' substrict two GC-MS data
#' @param data1 the first data
#' @param data2 the second data and will substrict the first one
#' @param lab1 character the title of the first data
#' @param lab2 character the title of the second data
#' @return list linear regression model for the diff matrix
#' @export
comparems <- function(data1,data2,lab1='',lab2=''){
        z <- data2-data1
        z[z<0] <- 0
        par(mfrow=c(2,3))
        plot(rowSums(data1),
             type = 'l',
             xlab = 'retention time',
             ylab = 'intensity',
             xaxt = 'n',
             frame.plot = F,
             main = lab1)
        ind = as.numeric(rownames(data1))
        axis( 1, at=seq(0,1500,300 ), labels= round((ind[1]+(max(ind)-min(ind))/5*(0:5))/60,2))
        plot(rowSums(data2),
             type = 'l',
             xlab = 'retention time',
             ylab = 'intensity',
             xaxt = 'n',
             frame.plot = F,
             main = lab2)
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
#' @param data imported data matrix of GC-MS
#' @param threshold the threshold of the response (log based 10) to seperate the group
#' @return list linear regression model for the data matrix
#' @export
plotgroup <- function(data,threshold=2){
        group <- ifelse(log10(data)>threshold,1,0)
        ind <- as.numeric(rownames(data))
        m <- matrix(c(1, 3, 2,  3), nrow = 2, ncol = 2)
        layout(m)
        hist(log10(data),
             breaks=100,
             main='',
             xlab='Intensity')
        abline(v=threshold,lwd=5,lty=2,col='red')
        legend('topright',"threshold",box.lty=0,pch = -1,lty=2,lwd=2,col='red')

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
        findline(data)
}

#' Plot the backgrond of data
#' @param data imported data matrix of GC-MS
#' @return NULL
#' @export
plotsub <- function(data){
        datan <- apply(data,1,diff)
        datan[datan<0] <- 0
        datan <- t(datan)
        plotms(datan)
}
#' Plot the intensity distribution of GC-MS
#' @param meanmatrix mean data matrix of GC-MS(n=5)
#' @param rsdmatrix standard deviation matrix of GC-MS(n=5)
#' @return NULL
#' @export
plotsms <- function(meanmatrix,rsdmatrix){
        par(mar=c(4.2,4.2,0,1.5),fig=c(0,1,0,0.8), new=F,cex.axis=1.5,cex.lab=1.5)
        smoothScatter(y=rsdmatrix*100,x=log10(c(meanmatrix)),
                      main = '',
                      xlab = 'Intensity',
                      ylab = 'Relative Standard Deviation(%)',
                      xaxt="n",
                      frame.plot = F)
        abline(h=20,lty = 2,col='red')
        abline(h=10,col='red')
        axis(1, at=c(0,1,2,3,4,5,6,7), labels= c('1','10','100','1000','10000','100000','1000000','100000000'))
        par(mar=c(0,4.2,1,1.5),oma=c(0,0,0,0),fig=c(0,1,0.8,1), new=T,cex.axis=1.5,cex.lab=1.5)
        hist(log10(meanmatrix),breaks = 100,xlab = 'Intensity',main = '',xaxt='n')
        axis(1, at=c(0,1,2,3,4,5,6,7), labels= c('1','10','100','1000','10000','100000','1000000','100000000'))
}
#' plot the density of the GC-MS data with EM algorithm to seperate the data into two log normal distribution.
#' @param data imported data matrix of GC-MS
#' @return NULL
#' @export
plothist <- function(data){
        data1=sample(data,100000)
        mixmdl = mixtools::normalmixEM(log10(data1))
        plot(mixmdl,which=2,breaks=100,
             xlab2 = 'Intensity')
        lines(density(log10(data1)), lty=2, lwd=2)
        legend("topright",c('noise','signal','density'),box.lty=0,pch = c(-1,-1,-1),lty=c(1,1,2),lwd=c(2,2,2),col=c('red','green','black'))
}
