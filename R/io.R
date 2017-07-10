#' Write MSP files for NIST search
#' @param mz a intensity vector, who name is the mass in m/z
#' @param outfilename the name of the MSP file, default is 'unknown'
#' @return none a MSP file will be created at the subfolder working dictionary with name 'MSP'
#' @examples
#' \dontrun{
#' mz <- c(10000,20000,10000,30000,5000)
#' names(mz) <- c(101,143,189,221,234)
#' writeMSP(mz,'test')
#' }
#' @export
writeMSP <- function(mz, outfilename = "unknown") {
    mz <- paste(names(mz), round(mz))
    dir.create("MSP")
    zz <- file(file.path("MSP", paste(outfilename,
        ".msp", sep = "")), "w")
    nPeaks <- length(mz)
    cat("Name: unknown", paste("Num Peaks: ", nPeaks),
        file = zz, sep = "\n")
    while (length(mz) >= 5) {
        cat(paste(mz[1:5]), "", file = zz, sep = "; ")
        cat(paste("\n"), file = zz)
        mz <- mz[6:length(mz)]
    }
    if (!is.na(mz[1])) {
        cat(paste(mz), "", file = zz, sep = "; ")
        cat(paste("\n"), file = zz)
    }
    close(zz)
    print(paste("A data file", outfilename, ".MSP has been generated in the folder:",
        "MSP", cat("\n")))
}
