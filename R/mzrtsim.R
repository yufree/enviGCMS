#' Generate simulated count data with batch effects for npeaks
#'
#' @param npeaks Number of peaks to simulate
#' @param ncomp percentage of compounds
#' @param ncond Number of conditions to simulate
#' @param ncpeaks percentage of peaks influenced by conditions
#' @param nbatch Number of batches to simulate
#' @param nbpeaks percentage of peaks influenced by batchs
#' @param npercond Number of samples per condition to simulate
#' @param nperbatch Number of samples per batch to simulate
#' @param shape shape for Weibull distribution of sample mean
#' @param scale scale for Weibull distribution of sample mean
#' @param shapersd shape for Weibull distribution of sample rsd
#' @param scalersd scale for Weibull distribution of sample rsd
#' @param seed Random seed for reproducibility
#' @details the numbers of batch columns should be the same with the condition columns.
#' @return list with rtmz data matrix, row index of peaks influenced by conditions, row index of peaks influenced by batchs, column index of conditions, column of batchs, raw condition matrix, raw batch matrix
#' @export
#' @examples
#' sim <- mzrtsim()
mzrtsim <- function(npeaks = 1000,
                    ncomp = 0.8,
                    ncond = 2,
                    ncpeaks = 0.05,
                    nbatch = 3,
                    nbpeaks = 0.1,
                    npercond = 10,
                    nperbatch = c(8,5,7),
                    shape = 2,
                    scale = 3,
                    shapersd = 1,
                    scalersd = 0.18,
                    seed = 42) {
        set.seed(seed)
        batch <- rep(1:nbatch, nperbatch)
        condition <- rep(1:ncond, npercond)
        # check the col numbers
        if(length(batch) != length(condition)){
                stop('Try to use the same numbers for both batch and condition columns')
        }
        ncol <- length(batch)
        # change the column order
        batch <- batch[sample(ncol)]
        condition <- condition[sample(ncol)]
        # generate the colname
        bc <- paste0('C', condition, 'B', batch)
        # get the compounds numbers
        ncomp <- npeaks*ncomp
        # get the matrix
        matrix <- matrix(0, nrow = ncomp, ncol = ncol)
        colnames(matrix) <- bc

        # generate the base peaks
        samplem <- 10^(stats::rweibull(ncomp, shape = shape, scale = scale))
        # generate the rsd for base peaks
        samplersd <- stats::rweibull(ncomp,shape = shapersd,scale = scalersd)

        for (i in 1:ncomp) {
                samplei <- abs(stats::rnorm(ncol,mean = samplem[i], sd = samplem[i]*samplersd[i]))
                matrix[i,] <- samplei
        }

        # get the multi-peaks
        nmpeaks <- npeaks - ncomp
        indexm <- sample(1:ncomp,nmpeaks,replace = T)
        change <- exp(stats::rnorm(nmpeaks))
        mpeaks <- matrix[indexm,]*change

        # get the matrix
        matrix0 <- matrix <- rbind(matrix,mpeaks)

        # get the numbers of signal and batch peaks
        ncpeaks <- npeaks * ncpeaks
        nbpeaks <- npeaks * nbpeaks
        # simulation of condition
        index <- sample(1:npeaks, ncpeaks)
        matrixc <- matrix[index, ]
        changec <- NULL
        for (i in 1:ncond) {
                colindex <- condition == i
                change <- exp(stats::rnorm(ncpeaks))
                matrixc[, colindex] <- matrixc[, colindex] * change
                changec <- cbind(changec,change)
        }
        matrix[index, ] <- matrixc
        # simulation of batch
        indexb <- sample(1:npeaks, nbpeaks)
        matrixb <- matrix[indexb, ]
        matrixb0 <- matrix0[indexb,]
        changeb <- NULL
        for (i in 1:nbatch) {
                colindex <- batch == i
                change <- exp(stats::rnorm(nbpeaks))
                matrixb[, colindex] <-
                        matrixb[, colindex] * change
                matrixb0[, colindex] <-
                        matrixb0[, colindex] * change
                changeb <- cbind(changeb,change)
        }
        matrix[indexb, ] <- matrixb
        # add row names
        rownames(matrix) <- paste0('P', c(1:npeaks))
        return(list(
                data = matrix,
                conp = index,
                batchp = indexb,
                con = condition,
                batch = batch,
                cmatrix = matrixc,
                changec = changec,
                bmatrix = matrixb0,
                changeb = changeb,
                matrix = matrix0,
                compmatrix = matrix[1:ncomp,]
        ))
}

#' Simulation from data by sample mean and rsd from Empirical Cumulative Distribution or Bootstrap sampling
#'
#' @param data matrix with row peaks and column samples
#' @param type 'e' means simulation from data by sample mean and rsd from Empirical Cumulative Distribution; 'f' means simulation from data by sample mean and rsd from Bootstrap sampling
#' @param npeaks Number of peaks to simulate
#' @param ncomp percentage of compounds
#' @param ncond Number of conditions to simulate
#' @param ncpeaks percentage of peaks influenced by conditions
#' @param nbatch Number of batches to simulate
#' @param nbpeaks percentage of peaks influenced by batchs
#' @param npercond Number of samples per condition to simulate
#' @param nperbatch Number of samples per batch to simulate
#' @param seed Random seed for reproducibility
#' @details the numbers of batch columns should be the same with the condition columns.
#' @return list with rtmz data matrix, row index of peaks influenced by conditions, row index of peaks influenced by batchs, column index of conditions, column of batchs, raw condition matrix, raw batch matrix
#' @export
#' @examples
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' list <- getmr(cdfpath, pmethod = ' ')
#' sim <- simmzrt(list$data)
simmzrt <- function(data,
                    type = 'e',
                    npeaks = 1000,
                    ncomp = 0.8,
                    ncond = 2,
                    ncpeaks = 0.05,
                    nbatch = 3,
                    nbpeaks = 0.1,
                    npercond = 10,
                    nperbatch = c(8,5,7),
                    seed = 42) {
        set.seed(seed)
        batch <- rep(1:nbatch, nperbatch)
        condition <- rep(1:ncond, npercond)
        # check the col numbers
        if(length(batch) != length(condition)){
                stop('Try to use the same numbers for both batch and condition columns')
        }

        ncol <- length(batch)
        # change the column order
        batch <- batch[sample(ncol)]
        condition <- condition[sample(ncol)]
        # generate the colname
        bc <- paste0('C', condition, 'B', batch)
        # get the compounds numbers
        ncomp <- npeaks*ncomp
        # get the matrix
        matrix <- matrix(0, nrow = ncomp, ncol = ncol)
        colnames(matrix) <- bc
        # get the mean, sd, and rsd from the data
        datamean <- apply(data,1,mean)
        datasd <- apply(data,1,sd)
        datarsd <- datasd/datamean
        if(type == 'e'){
                # generate the distribution from mean and rsd
                pdfmean <- ecdf(datamean)
                pdfrsd <- ecdf(datarsd)
                # simulate the sample mean and rsd from ecdf
                simmean <- as.numeric(quantile(pdfmean, runif(ncomp)))
                simrsd <- as.numeric(quantile(pdfrsd, runif(ncomp)))
        }else if(type == 'b'){
                # simulate the sample mean and rsd from bootstrap
                simmean <- sample(datamean,ncomp,replace = T)
                simrsd <- sample(datarsd,ncomp,replace = T)
        } else{
                stop("type should be 'e' or 'b'")
        }

        # simulate the data
        for (i in 1:ncomp) {
                samplei <- abs(stats::rnorm(ncol,mean = simmean[i], sd = simmean[i]*simrsd[i]))
                matrix[i,] <- samplei
        }

        # get the multi-peaks
        nmpeaks <- npeaks - ncomp
        indexm <- sample(1:ncomp,nmpeaks,replace = T)
        change <- exp(stats::rnorm(nmpeaks))
        mpeaks <- matrix[indexm,]*change

        # get the matrix
        matrix0 <- matrix <- rbind(matrix,mpeaks)

        # get the numbers of signal and batch peaks
        ncpeaks <- npeaks * ncpeaks
        nbpeaks <- npeaks * nbpeaks
        # simulation of condition
        index <- sample(1:npeaks, ncpeaks)
        matrixc <- matrix[index, ]
        changec <- NULL
        for (i in 1:ncond) {
                colindex <- condition == i
                change <- exp(stats::rnorm(ncpeaks))
                matrixc[, colindex] <- matrixc[, colindex] * change
                changec <- cbind(changec,change)
        }
        matrix[index, ] <- matrixc
        # simulation of batch
        indexb <- sample(1:npeaks, nbpeaks)
        matrixb <- matrix[indexb, ]
        matrixb0 <- matrix0[indexb,]
        changeb <- NULL
        for (i in 1:nbatch) {
                colindex <- batch == i
                change <- exp(stats::rnorm(nbpeaks))
                matrixb[, colindex] <-
                        matrixb[, colindex] * change
                matrixb0[, colindex] <-
                        matrixb0[, colindex] * change
                changeb <- cbind(changeb,change)
        }
        matrix[indexb, ] <- matrixb
        # add row names
        rownames(matrix) <- paste0('P', c(1:npeaks))
        return(list(
                data = matrix,
                conp = index,
                batchp = indexb,
                con = condition,
                batch = batch,
                cmatrix = matrixc,
                changec = changec,
                bmatrix = matrixb0,
                changeb = changeb,
                matrix = matrix0,
                compmatrix = matrix[1:ncomp,]
        ))
}
