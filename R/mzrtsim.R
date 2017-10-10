#' Generate simulated count data with batch effects for npeaks
#'
#' @param npeaks Number of peaks to simulate
#' @param ncond Number of conditions to simulate
#' @param ncpeaks percentage of peaks influenced by conditions
#' @param upc up regulation proportion in peaks influenced by conditions
#' @param nbatch Number of batches to simulate
#' @param nbpeaks percentage of peaks influenced by batchs
#' @param upb up regulation proportion in peaks influenced by batchs
#' @param npercond Number of samples per condition per batch to simulate
#' @param shape shape for Weibull distribution
#' @param scale scale for Weibull distribution
#' @param seed Random seed for reproducibility
#' @return list with rtmz data matrix, row index of peaks influenced by conditions, row index of peaks influenced by batchs, column index of conditions, column of batchs, raw condition matrix, raw batch matrix
#' @export
#' @examples
#' sim <- mzrtsim()
mzrtsim <- function(npeaks = 1000,
                    ncond = 2,
                    ncpeaks = 0.05,
                    upc = 0.5,
                    nbatch = 3,
                    nbpeaks = 0.1,
                    upb = 0.3,
                    npercond = 10,
                    shape = 8,
                    scale = 12,
                    seed = 42) {
        set.seed(seed)
        ncol <- nbatch * ncond * npercond
        batch <- rep(1:nbatch, each = ncond * npercond)
        condition <- rep(rep(1:ncond, each = npercond), nbatch)
        bc <- paste0('C', condition, 'B', batch)
        matrix <- matrix(0, nrow = npeaks, ncol = ncol)
        colnames(matrix) <- bc

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
                if (i != 1) {
                        colindex <- condition == i
                        flips <-
                                stats::rbinom(ncpeaks,
                                       size = 1,
                                       prob = upc)
                        change <- stats::runif(ncpeaks, 1, 10)
                        flipschange <-
                                ifelse(flips == 0, change, 1 / change)
                        matrixc[, colindex] <-
                                matrixc[, colindex] * flipschange
                }
        }
        matrix[index, ] <- matrixc
        # simulation of batch
        indexb <- sample(1:npeaks, nbpeaks)
        matrixb <- matrix[indexb, ]
        matrixb0 <- matrix0[indexb,]
        for (i in 1:nbatch) {
                if (i != 1) {
                        colindex <- batch == i
                        flips <- stats::rbinom(nbpeaks,
                                        size = 1,
                                        prob = upb)
                        change <- stats::runif(nbpeaks, 1, 10)
                        flipschange <-
                                ifelse(flips == 0, change, 1 / change)
                        matrixb[, colindex] <-
                                matrixb[, colindex] * flipschange
                        matrixb0[, colindex] <-
                                matrixb0[, colindex] * flipschange
                }
        }
        matrix[indexb, ] <- matrixb
        # add row names
        rownames(matrix) <- paste0('P', c(1:npeaks))
        return(list(
                matrix = matrix,
                conp = index,
                batchp = indexb,
                con = condition,
                batch = batch,
                cmatrix = matrixc,
                bmatrix = matrixb0
        ))
}
