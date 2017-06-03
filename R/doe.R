#' Get the features from t test, with p value, q value, rsd and power restriction
#' @param xod xcmsSet objects
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param rsdt rsd threshold to filter the data
#' @return dataframe with peaks fit the setting above
#' @export

getfeaturest <- function(xod, power = 0.8, pt = 0.05, 
    qt = 0.05, n = 3, rsdt = 50) {
    data <- xcms::groupval(xod, value = "into")
    idx <- stats::complete.cases(data)
    data1 <- as.data.frame(xcms::groups(xod))
    lv <- xod@phenoData[, 1]
    data0 <- data[idx, ]
    sd <- genefilter::rowSds(data0[, 1:n])
    mean <- stats::aggregate(t(data0), list(lv), mean)
    mean <- t(mean[, -1])
    rsd <- sd/mean[, 1] * 100
    sd <- sd[rsd < rsdt]
    mz <- data1$mzmed[idx & rsd < rsdt]
    rt <- data1$rtmed[idx & rsd < rsdt]
    data0 <- data0[rsd < rsdt, ]
    ar <- genefilter::rowttests(data0, fac = lv)
    dm <- ar$dm
    p <- ar$p.value
    q <- stats::p.adjust(p, method = "BH")
    m <- nrow(data0)
    df <- cbind.data.frame(sd, dm, p, q, mz, rt, data0)
    df <- df[order(df$p), ]
    df$alpha <- c(1:m) * pt/m
    rp <- vector()
    for (i in c(1:nrow(df))) {
        
        r <- stats::power.t.test(delta = df$dm[i], 
            sd = df$sd[i], sig.level = df$alpha[i], 
            n = n)
        rp[i] <- r$power
    }
    df <- cbind(power = rp, df)
    df <- df[df$power > power, ]
    return(df)
}

#' Get the features from anova, with p value, q value, rsd and power restriction
#' @param xod xcmsSet objects
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param ng group numbers
#' @param rsdt rsd threshold to filter the data
#' @return dataframe with peaks fit the setting above
#' @export

getfeaturesanova <- function(xod, power = 0.8, pt = 0.05, 
    qt = 0.05, n = 3, ng = 3, rsdt = 50) {
    data <- xcms::groupval(xod, value = "into")
    idx <- stats::complete.cases(data)
    data1 <- as.data.frame(xcms::groups(xod))
    lv <- xod@phenoData[, 1]
    data0 <- data[idx, ]
    
    sd <- genefilter::rowSds(data0[, 1:n])
    mean <- stats::aggregate(t(data0), list(lv), mean)
    mean <- t(mean[, -1])
    sd2 <- genefilter::rowSds(mean)
    rsd <- sd/mean[, 1] * 100
    sd2 <- sd2[rsd < rsdt]
    sd <- sd[rsd < rsdt]
    mz <- data1$mzmed[idx & rsd < rsdt]
    rt <- data1$rtmed[idx & rsd < rsdt]
    data0 <- data0[rsd < rsdt, ]
    ar <- genefilter::rowFtests(data0, lv)
    p <- ar$p.value
    q <- stats::p.adjust(p, method = "BH")
    m <- nrow(data0)
    df <- cbind.data.frame(sd, sd2, p, q, mz, rt, data0)
    df <- df[order(df$p), ]
    df$alpha <- c(1:m) * pt/m
    rp <- vector()
    for (i in c(1:nrow(df))) {
        r <- stats::power.anova.test(groups = ng, between.var = df$sd2[i], 
            within.var = df$sd[i], sig.level = df$alpha[i], 
            n = n)
        rp[i] <- r$power
    }
    df <- cbind(power = rp, df)
    df <- df[df$power > power, ]
    return(df)
}
