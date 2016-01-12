#' findline of the regression model
findline <- function(data,threshold=2){
        x <- as.numeric(rownames(data))[-1]
        y0 <- as.numeric(colnames(data))[-1]
        group <- ifelse(log10(data)>threshold,0,1)
        diffmatrix <- apply(group,2,diff)
        difftemp <- apply(diffmatrix,1,which.max)
        y <- y0[difftemp]
        data <- data.frame(y,x)
        data <- data[data$y>101,]
        plot(data$y~data$x,
             xlab = 'Retention time(min)',
             ylab = 'm/z',
             pch = 19,
             xaxt = "n",
             yaxt = 'n',
             main = 'Regression Model for Signals',
             frame.plot = F)
        axis( 2, at=seq(0,1000,100), labels= 100+100*c(0:10), las= 2 )
        axis( 1, at=seq(min(data$x),max(data$x),(max(data$x)-min(data$x))/22 ), labels= (100+(220)/22*(0:22)))
        abline(lm(data$y~data$x), col="red",lwd = 5)
        lines(lowess(data$y~data$x), col="blue",lwd = 5)
        fit <- lm(data$y~data$x)
        rmse <- round(sqrt(mean(resid(fit)^2)), 2)
        coefs <- coef(fit)
        b0 <- round(coefs[1],2)
        b1 <- round(coefs[2],2)
        r2 <- round(summary(fit)$r.squared, 2)
        pv <- round(anova(fit)$'Pr(>F)'[1], 2)
        eqn <- bquote(italic(Mass(m/z)) == .(b0) + .(b1)*"*"*italic(Temperature(~degree~C)))
        eqn2 <- bquote(r^2 == .(r2) * "," ~~ p < 0.01)
        text(200, 700, adj=c(0,0), cex = 1,eqn)
        text(200, 650, adj=c(0,0), cex = 1,eqn2)
        legend("bottomright",c('OLS','LOWESS'),box.lty=0,pch = c(-1,-1),lty=c(1,1),lwd=c(2,2),col=c('red','blue'))
        return(fit)
}

#' substrict two GC-MS data
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
plotgroup <- function(data,threshold=2){
        group <- ifelse(log10(data)>threshold,1,-1)
        ind <- as.numeric(rownames(data))
        par(mfrow=c(1,2))
        image(group,
              xlab = 'retention time(min)',
              ylab = 'm/z',
              axes = F,
              col = heat.colors(2),
              useRaster = T)
        axis( 2, at=seq(0,1,1/9), labels= 100+100*c(0:9), las= 2 )
        axis( 1, at=seq(0,1,0.1 ), labels= round((ind[1]+(max(ind)-min(ind))/10*(0:10))/60,2))
        hist(log10(data),
             breaks=100,
             main='',
             xlab='Distribution of Intensity(log base 10)')
        abline(v=threshold,lwd=5,lty=2,col='red')
        legend("topright","threshold",box.lty=0,pch = -1,lty=2,lwd=2,col='red')
}

#' Plot the backgrond of data
plotsub <- function(data,...){
        datan <- apply(data,1,diff)
        datan[datan<0] <- 0
        datan <- t(datan)
        plotms(datan,...)
        findline(datan)
}
#' Plot the intensity distribution of GC-MS
plotsms <- function(meanmatrix,sdmatrix,...){
        smoothScatter(log10(c(sdmatrix))~log10(c(meanmatrix)),
                      xlab = 'Mean(log scale based 10)',
                      ylab = 'Standard Deviation(log scale based 10)',
                      frame.plot = F,...)
        abline(a = 0, b = 1)
        abline(v=2)
}
plothist <- function(data){
        data1=sample(data,100000)
        mixmdl = normalmixEM(log10(data1))
        plot(mixmdl,which=2,breaks=100,
             xlab2 = 'Intensity(log scale based 10)')
        lines(density(log10(data1)), lty=2, lwd=2)
        legend("topright",c('noise','signal','density'),box.lty=0,pch = c(-1,-1,-1),lty=c(1,1,2),lwd=c(2,2,2),col=c('red','green','black'))
}
