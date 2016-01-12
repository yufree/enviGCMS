#' import data and return the annotated matrix
#' @param data file type which xcmsRaw could handle
#' @return matrix with the row as increasing m/z in second and column as increasing scantime
getmd <- function(data){
        data <- xcmsRaw(data)
        z1 <- data@env$profile
        zf <- as.factor(round(data@scantime))
        df <- aggregate(t(z1), list(zf), sum)[-1]
        rownames(df) <- unique(round(data@scantime))
        colnames(df) <- seq(data@mzrange[1],data@mzrange[2])
        return(t(as.matrix(df)))
}

#' Combine two or more dataset
#'
#' @param data1 matrix lower mass range
#' @param data2 matrix higher mass range
#' @param ... matrix even higher mass range
#' @return matrix with the row as scantime in second and column as m/z
#'
conbinemd <- function(data1,data2,...){
        if(missing(...)){
                z1 <- getmd(data1)
                z2 <- getmd(data2)
                ind <- intersect(colnames(z1),colnames(z2))
                z <- rbind(,z1[as.character(ind)],z2[,as.character(ind)])
                rownames(z) <- c(seq(min(colnames(z1)),max(colnames(z2))),seq(min(colnames(z2)),max(colnames(z3))))
                colnames(z) <- ind
                return(z)
        }
        else{
                conbinedmd(data1, conbinedmd(data2, ...))
        }
}


#' Subset the data mass spectrum of certain retention time and plot them
#' @param data imported data matrix of GC-MS
#' @param rt vector range of the retention time
#' @param ms vector range of the m/z
#' @return data matrix
getsubmd <- function(data,rt,ms){
        mzindexstart <- as.numeric(head(rownames(data),1))
        rtindexstart <- as.numeric(head(colnames(data),1))
        rts <- rt*60-rtindexstart
        rt1 <- min(rts)
        rt2 <- max(rts)
        mzs <- ms-mzindexstart+1
        mz1 <- min(mzs)
        mz2 <- max(mzs)
        data <- data[mz1:mz2,]
        data <- t(data)[rt1:rt2,]
        return(t(data))
}

#' Write MSP files for NIST search
#' @param mz a intensity vector, who name is the mass in m/z
#' @param outfilename the name of the MSP file, default is 'unknown'
#' @return none a MSP file will be created at the subfolder working dictionary with name 'MSP'
writeMSP<-function(mz, outfilename="unknown"){
        mz <- paste(names(mz),round(mz))
        dir.create('MSP')
        zz <- file(file.path('MSP',paste(outfilename,".msp",sep="")), "w")
        nPeaks <- length(mz)
        cat("Name: unknown", paste("Num Peaks: ",nPeaks),  file = zz, sep = "\n")
        while (length(mz) >=5 ){
                cat(paste(mz[1:5]),"", file = zz, sep="; ")
                cat(paste("\n"), file = zz)
                mz<-mz[6:length(mz)]
        }
        if(!is.na(mz[1])){
                cat(paste(mz),"", file = zz, sep="; ")
                cat(paste("\n"), file = zz)
        }
        close(zz)
        print(paste("A data file",outfilename,".MSP has been generated in the folder:", 'MSP', cat("\n")))
}
