#' Generate simulated count data with batch effects for npeaks
#'
#' @param npeaks Number of peaks to simulate
#' @param ncond Number of conditions to simulate
#' @param ncpeaks percentage of peaks influenced by conditions
#' @param nbatch Number of batches to simulate
#' @param nbpeaks percentage of peaks influenced by batchs
#' @param npercond Number of samples per condition to simulate
#' @param nperbatch Number of samples per batch to simulate
#' @param shape shape for Weibull distribution
#' @param scale scale for Weibull distribution
#' @param seed Random seed for reproducibility
#' @details the numbers of batch columns should be the same with the condition columns.
#' @return list with rtmz data matrix, row index of peaks influenced by conditions, row index of peaks influenced by batchs, column index of conditions, column of batchs, raw condition matrix, raw batch matrix
#' @export
#' @examples
#' sim <- mzrtsim()
mzrtsim <- function(npeaks = 1000,
                    ncond = 2,
                    ncpeaks = 0.05,
                    nbatch = 3,
                    nbpeaks = 0.1,
                    npercond = 10,
                    nperbatch = c(8,5,7),
                    shape = 8,
                    scale = 12,
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
        matrix <- matrix(0, nrow = npeaks, ncol = ncol)
        colnames(matrix) <- bc
        # generate the base peaks
        for (i in 1:ncol) {
                sampleio <- exp(stats::rweibull(npeaks, shape = shape, scale = scale))
                samplei <- sampleio[order(sampleio)]
                matrix[, i] <- samplei
        }
        # reorder the peaks
        matrix0 <- matrix <- matrix[sample(nrow(matrix)), ]
        # get the numbers of signal and batch peaks
        ncpeaks <- npeaks * ncpeaks
        nbpeaks <- npeaks * nbpeaks
        # simulation of condition
        index <- sample(1:npeaks, ncpeaks)
        matrixc <- matrix[index, ]
        for (i in 1:ncond) {
                colindex <- condition == i
                change <- exp(stats::rnorm(ncpeaks))
                matrixc[, colindex] <- matrixc[, colindex] * change
        }
        matrix[index, ] <- matrixc
        # simulation of batch
        indexb <- sample(1:npeaks, nbpeaks)
        matrixb <- matrix[indexb, ]
        matrixb0 <- matrix0[indexb,]
        for (i in 1:nbatch) {
                colindex <- batch == i
                change <- exp(stats::rnorm(nbpeaks))
                matrixb[, colindex] <-
                                matrixb[, colindex] * change
                matrixb0[, colindex] <-
                                matrixb0[, colindex] * change
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
                bmatrix = matrixb0,
                matrix = matrix0
        ))
}
