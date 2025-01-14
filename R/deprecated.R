#' Get the peak list with blank samples' peaks removed
#' @param xset the xcmsset object with blank and certain group samples' data
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 1000
#' @return diff report
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getbgremove(xset)
#' }
getbgremove <-
        function(xset,
                 method = "medret",
                 intensity = "into",
                 file = NULL,
                 rsdcf = 30,
                 inscf = 1000) {
                .Deprecated("getdoe")
                message("This function has been deprecated and you could use getdoe to remove background.")
        }

#' Get the report for biological replicates.
#' @param xset the xcmsset object which for all of your technique replicates for bio replicated sample in single group
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 0
#' @return dataframe with mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data
getbiotechrep <-
        function(xset,
                 method = "medret",
                 intensity = "into",
                 file = NULL,
                 rsdcf = 30,
                 inscf = 1000) {
                .Deprecated("getdoe")
                message("This function has been deprecated and you could use getdoe to process data.")
        }

#' Get the report for samples with biological and technique replicates in different groups
#' @param xset the xcmsset object all of samples with technique replicates
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for the peaklist to MetaboAnalyst
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 1000
#' @return dataframe with mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data in different groups if file are defaults NULL.
getgrouprep <-
        function(xset,
                 file = NULL,
                 method = "medret",
                 intensity = "into",
                 rsdcf = 30,
                 inscf = 1000) {
                .Deprecated("getdoe")
                message("This function has been deprecated and you could use getdoe to process data.")
        }

#' output the similarity of two dataset
#' @param xset1 the first dataset
#' @param xset2 the second dateset
#' @return similarity on retention time and rsd %
getsim <- function(xset1, xset2) {
        .Deprecated()
        message("This function has been deprecated.")
}

#' Get the report for technique replicates.
#' @param xset the xcmsset object which for all of your technique replicates for one sample
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for further annotation, default NULL
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 1000
#' @return dataframe with mean, standard deviation and RSD for those technique replicates combined with raw data
gettechrep <-
        function(xset,
                 method = "medret",
                 intensity = "into",
                 file = NULL,
                 rsdcf = 30,
                 inscf = 1000) {
                .Deprecated("getdoe")
                message("This function has been deprecated and you could use getdoe to process data.")
        }

#' Get the time series or two factor DoE report for samples with biological and technique replicates in different groups
#' @param xset the xcmsset object all of samples with technique replicates in time series or two factor DoE
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param file file name for the peaklist to MetaboAnalyst
#' @param rsdcf rsd cutoff for peaks, default 30
#' @param inscf intensity cutoff for peaks, default 1000
#' @return dataframe with time series or two factor DoE mean, standard deviation and RSD for those technique replicates & biological replicates combined with raw data in different groups if file are defaults NULL.
gettimegrouprep <-
        function(xset,
                 file = NULL,
                 method = "medret",
                 intensity = "into",
                 rsdcf = 30,
                 inscf = 1000) {
                .Deprecated("getdoe")
                message("This function has been deprecated and you could use getdoe to process data.")
        }

#' Plot the influences of DoE and Batch effects on each peaks
#' @param df data output from `svacor` function
#' @param dfsv data output from `svaplot` function for corrected data
#' @param dfanova data output from `svaplot` function for raw data
#' @return influences plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' dfsv <- svaplot(xset3)
#' dfanova <- svaplot(xset3, pqvalues = "anova")
#' svabatch(df,dfsv,dfanova)
#' }
#' @seealso \code{\link{svacor}}, \code{\link{svaplot}}, \code{\link{svapca}}
svabatch <- function(df, dfsv, dfanova) {
        .Deprecated()
        message(
                "This function has been deprecated and you could use mzrtsim package for batch effect correction."
        )
}

#' Surrogate variable analysis(SVA) to correct the unknown batch effects
#' @param xset xcmsset object
#' @param lv group information
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @details this is used for reviesed version of SVA to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, signal part, random errors part, batch part, p-values, q-values, mass, rt, Posterior Probabilities of Surrogate variables and Posterior Probabilities of Mod. If no surrogate variable found, corresponding part would miss.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' }
#' @seealso \code{\link{svapca}}, \code{\link{svaplot}}, \code{\link{svabatch}}
svacor <- function(xset,
                   lv = NULL,
                   method = "medret",
                   intensity = "into") {
        .Deprecated()
        message(
                "This function has been deprecated and you could use mzrtsim package for batch effect correction."
        )
}

#' Filter the data with p value and q value
#' @param list results from svacor function
#' @param pqvalues method for ANOVA or SVA
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @return data, corrected data, mz and retention for filerted data
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' svadata(df)
#' }
svadata <- function(list,
                    pqvalues = "sv",
                    pt = 0.05,
                    qt = 0.05) {
        .Deprecated()
        message(
                "This function has been deprecated and you could use mzrtsim package for batch effect correction."
        )
}

#' Principal component analysis(PCA) for SVA corrected data and raw data
#' @param list results from svacor function
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param lv group information
#' @return plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' svapca(df)
#' }
#' @seealso \code{\link{svacor}}, \code{\link{svaplot}}, \code{\link{svabatch}}
svapca <- function(list,
                   center = TRUE,
                   scale = TRUE,
                   lv = NULL) {
        .Deprecated()
        message(
                "This function has been deprecated and you could use mzrtsim package for batch effect correction."
        )
}

#' Filter the data with p value and q value and show them
#' @param list results from svacor function
#' @param pqvalues method for ANOVA or SVA
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @param lv group information
#' @param index index for selected peaks
#' @return heatmap for the data
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' df <- svacor(xset3)
#' svaplot(df)
#' }
#' @seealso \code{\link{svacor}}, \code{\link{svapca}}, \code{\link{svabatch}}
svaplot <- function(list,
                    pqvalues = "sv",
                    pt = 0.05,
                    qt = 0.05,
                    lv = NULL,
                    index = NULL) {
        .Deprecated()
        message(
                "This function has been deprecated and you could use mzrtsim package for batch effect correction."
        )
}

#' Get the corrected data after SVA for metabolanalyst
#' @param xset xcmsset object
#' @param lv group information
#' @return csv files for both raw and corrected data for metaboanalyst if SVA could be applied
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' xset <- xcmsSet(cdffiles)
#' xset <- group(xset)
#' xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
#' xset2 <- group(xset2, bw = 10)
#' xset3 <- fillPeaks(xset2)
#' svaupload(xset3)
#' }
svaupload <- function(xset, lv = NULL) {
        .Deprecated()
        message(
                "This function has been deprecated and you could use mzrtsim package for batch effect correction."
        )
}
#' Get the csv files from xcmsset/XCMSnExp/list object
#' @param xset the xcmsset/XCMSnExp/list object which you want to submitted to Metaboanalyst
#' @param method parameter for groupval function
#' @param value parameter for groupval function
#' @param name file name
#' @param type m means  Metaboanalyst, a means xMSannotator, o means full information csv
#' @param mzdigit m/z digits of row names of data frame
#' @param rtdigit retention time digits of row names of data frame
#' @return dataframe with data needed for Metaboanalyst/xMSannotator/pmd if your want to perform local analysis.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getupload(xset)
#' }
#' @seealso \code{\link{getdata}}, \code{\link{getmzrt}}
getupload <-
        function(xset,
                 method = "medret",
                 value = "into",
                 name = "Peaklist",
                 type = 'm',
                 mzdigit = 4,
                 rtdigit = 1) {
                .Deprecated()
                message("This function has been deprecated and you could use getmzrt to get csv file.")
        }
#' Get the csv files to be submitted to Metaboanalyst
#' @param xset a XCMSnExp object with processed data which you want to submitted to Metaboanalyst
#' @param value value for `xcms::featureValues`
#' @param name file name
#' @return dataframe with data needed for Metaboanalyst if your want to perform local analysis.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata2(cdfpath)
#' getupload2(xset)
#' }
#' @seealso \code{\link{getdata2}},\code{\link{getupload}}, \code{\link{getmzrt2}}

getupload2 <- function(xset, value = "into", name = "Peaklist") {
        .Deprecated()
        message("This function has been deprecated and you could use getupload to get csv file.")
}

#' Get the csv files to be submitted to Metaboanalyst
#' @param list list with data as peaks list, mz, rt and group information
#' @param name file name
#' @return dataframe with data needed for Metaboanalyst if your want to perform local analysis.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata2(cdfpath,
#' ppp = xcms::MatchedFilterParam(),
#' rtp = xcms::ObiwarpParam(),
#' gpp = xcms::PeakDensityParam())
#' xset <- enviGCMS::getmzrt2(xset)
#' getupload3(xset)
#' }
#' @seealso \code{\link{getmzrt}}, \code{\link{getmzrt2}}

getupload3 <- function(list, name = "Peaklist") {
        .Deprecated()
        message("This function has been deprecated and you could use getupload to get csv file.")
}
#' Get the mzrt profile and group information for batch correction and plot as a list for xcms 3 object
#' @param xset a XCMSnExp object with processed data
#' @param name file name for csv file, default NULL
#' @return list with rtmz profile and group information
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata2(cdfpath,
#' ppp = xcms::MatchedFilterParam(),
#' rtp = xcms::ObiwarpParam(),
#' gpp = xcms::PeakDensityParam())
#' getmzrt2(xset)
#' }
#' @seealso \code{\link{getdata2}},\code{\link{getupload2}}, \code{\link{getmzrt}}, \code{\link{getdoe}},\code{\link{getmzrtcsv}}

getmzrt2 <- function(xset, name = NULL) {
        .Deprecated()
        message(
                "This function has been deprecated and you could use getmzrt to get list object or csv file."
        )
}
#' plot EIC and boxplot for all peaks and return diffreport
#' @param xset xcmsset object
#' @param name filebase of the sub dir
#' @param test 't' means two-sample welch t-test, 't.equalvar' means two-sample welch t-test with equal variance, 'wilcoxon' means rank sum wilcoxon test, 'f' means F-test, 'pairt' means paired t test, 'blockf' means Two-way analysis of variance, default 't'
#' @param nonpara 'y' means using nonparametric ranked data, 'n' means original data
#' @param ... other parameters for `diffreport`
#' @return diffreport and pdf figure for EIC and boxplot
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' plote(xset)
#' }

plote <- function(xset,
                  name = "test",
                  test = "t",
                  nonpara = "n",
                  ...) {
        .Deprecated()
        message(
                "This function will be deprecated and you could use getmzrt to get related object to plot EIC."
        )
}

#' Get the features from t test, with p value, q value, rsd and power restriction
#' @param list list with data as peaks list, mz, rt and group information (two groups)
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param imputation parameters for `getimputation` function method
#' @return dataframe with peaks fit the setting above

getfeaturest <- function(list,
                         power = 0.8,
                         pt = 0.05,
                         qt = 0.05,
                         n = 3,
                         imputation = "l") {
        .Deprecated()
        message(
                "This function has been deprecated and you could use getpower to calculate post-hoc power for each peaks."
        )
}

#' Get the features from anova, with p value, q value, rsd and power restriction
#' @param list list with data as peaks list, mz, rt and group information (more than two groups)
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param ng group numbers
#' @param inscf Log intensity cutoff for peaks across samples. If any peaks show a intensity higher than the cutoff in any samples, this peaks would not be filtered. default 5
#' @param rsdcf the rsd cutoff of all peaks in all group
#' @param imputation parameters for `getimputation` function method
#' @param index the index of peaks considered, default NULL
#' @return dataframe with peaks fit the setting above
getfeaturesanova <- function(list,
                             power = 0.8,
                             pt = 0.05,
                             qt = 0.05,
                             n = 3,
                             ng = 3,
                             rsdcf = 100,
                             inscf = 5,
                             imputation = "l",
                             index = NULL) {
        .Deprecated()
        message(
                "This function has been deprecated and you could use getpower to calculate post-hoc power for each peaks."
        )
}

#' Get xcmsset object in one step with optimized methods.
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the reference
#' @param minfrac minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group, default 0.67
#' @param ... arguments for xcmsSet function
#' @details the parameters are extracted from the papers. If you use name other than the name above, you will use the default setting of XCMS. Also I suggest IPO packages or apLCMS packages to get reasonable data for your own instrumental. If you want to summit the results to a paper, remember to include those parameters.
#' @return a xcmsset object for that path or selected samples
#' @references Patti, G. J.; Tautenhahn, R.; Siuzdak, G. Nat. Protocols 2012, 7 (3), 508–516.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' }
#' @seealso \code{\link{getdata2}}, \code{\link{getmzrt}}
#' @export
getdata <-
        function(path,
                 index = FALSE,
                 BPPARAM = BiocParallel::SnowParam(),
                 pmethod = "hplcorbitrap",
                 minfrac = 0.67,
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
                if (index) {
                        cdffiles <- cdffiles[index]
                }
                .Deprecated()
                message(
                        "This function has been deprecated and you might check documents of xcms to process your data."
                )
        }
#' Get XCMSnExp object in one step from structured folder path for xcms 3.
#' @param path the path to your data
#' @param index the index of the files
#' @param snames sample names. By default the file name without extension is used
#' @param sclass sample classes.
#' @param phenoData data.frame or NAnnotatedDataFrame defining the sample names and classes and other sample related properties. If not provided, the argument sclass or the subdirectories in which the samples are stored will be used to specify sample grouping.
#' @param BPPARAM used for BiocParallel package
#' @param mode 'inMemory' or 'onDisk' see `?MSnbase::readMSData` for details, default 'onDisk'
#' @param ppp parameters for peaks picking, e.g. xcms::CentWaveParam()
#' @param rtp parameters for retention time correction, e.g. xcms::ObiwarpParam()
#' @param gpp parameters for peaks grouping, e.g. xcms::PeakDensityParam()
#' @param fpp parameters for peaks filling, e.g. xcms::FillChromPeaksParam(), PeakGroupsParam()
#' @details This is a wrap function for metabolomics data process for xcms 3.
#' @return a XCMSnExp object with processed data
#' @seealso \code{\link{getdata}},\code{\link{getmzrt}}
#' @export
getdata2 <- function(path,
                     index = FALSE,
                     snames = NULL,
                     sclass = NULL,
                     phenoData = NULL,
                     BPPARAM = BiocParallel::SnowParam(),
                     mode = "onDisk",
                     ppp,
                     rtp,
                     gpp,
                     fpp) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )

}
#' Get the mzrt profile and group information for batch correction and plot as a list directly from path with default setting
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the references
#' @param minfrac minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group, default 0.67
#' @param ... arguments for xcmsSet function
#' @return list with rtmz profile and group infomation
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' }
#' @seealso \code{\link{getdata}},\code{\link{getupload}}, \code{\link{getmzrt}}, \code{\link{getdoe}}
getmr <-
        function(path,
                 index = FALSE,
                 BPPARAM = BiocParallel::SnowParam(),
                 pmethod = "hplcorbitrap",
                 minfrac = 0.67,
                 ...) {
                .Deprecated()
                message(
                        "This function has been deprecated and you might check documents of xcms to process your data."
                )
        }

#' get the data of QC compound for a group of data
#' @param path data path for your QC samples
#' @param mzrange mass of the QC compound
#' @param rtrange retention time of the QC compound
#' @param index index of the files contained QC compounds, default is all of the compounds
#' @return number vector, each number indicate the peak area of that mass and retention time range
getQCraw <- function(path, mzrange, rtrange, index = NULL) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}
#' Import data and return the annotated matrix for GC/LC-MS by m/z range and retention time
#' @param data file type which xcmsRaw could handle
#' @param mzstep the m/z step for generating matrix data from raw mass spectral data
#' @param mzrange vector range of the m/z, default all
#' @param rtrange vector range of the retention time, default all
#' @return matrix with the row as increasing m/z second and column as increasing scantime
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' matrix <- getmd(cdffiles[1])
#' }
getmd <- function(data,
                  mzstep = 0.1,
                  mzrange = FALSE,
                  rtrange = FALSE) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}
#' Get the differences of two GC/LC-MS data
#'
#' @param data1 data file path of first data
#' @param data2 data file path of second data
#' @param mzstep the m/z step for generating matrix data from raw mass spectral data
#' @param rtstep the alignment accuracy of retention time, e.g. 0.01 means the retention times of combined data should be the same at the accuracy 0.01s. Higher rtstep would return less scans for combined data
#' @return list four matrix with the row as scantime in second and column as m/z, the first matrix refer to data 1, the second matrix refer to data 2, the third matrix refer to data1 - data2 while the fourth refer to data2 - data1, minus values are imputed by 0
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' matrix <- submd(cdffiles[1],cdffiles[7])
#' }
#' @export
submd <- function(data1,
                  data2,
                  mzstep = 0.1,
                  rtstep = 0.01) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )

}
#' Plot mass spectrum of certain retention time and return mass spectrum vector (MSP file) for NIST search
#' @param data imported data matrix of GC-MS
#' @param rt vector range of the retention time
#' @param ms vector range of the m/z
#' @param msp logical, return MSP files or not, default False
#' @return plot, vector and MSP files for NIST search
#' @examples
#' \dontrun{
#' plotrtms(matrix,rt = c(500,1000),ms = c(300,500))
#' }
#' @export
plotrtms <- function(data, rt, ms, msp = FALSE) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}
#' Plot EIC of certain m/z and return dataframe for integration
#' @param data imported data matrix of GC-MS
#' @param ms m/z to be extracted
#' @param rt vector range of the retention time
#' @param n logical smooth or not
#' @return dataframe with  with the first column RT and second column intensity of the SIM ions.
#' @examples
#' \dontrun{
#' matrix <- getmd(rawdata)
#' plotmsrt(matrix,rt = c(500,1000),ms = 300)
#' }
#' @export
plotmsrt <- function(data, ms, rt, n = FALSE) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}

#' Combine two data with similar retention time while different mass range
#'
#' @param data1 data file path of lower mass range
#' @param data2 data file path of higher mass range
#' @param mzstep the m/z step for generating matrix data from raw mass spectral data
#' @param rtstep the alignment accuracy of retention time, e.g. 0.01 means the retention times of combined data should be the same at the accuracy 0.01s. Higher rtstep would return less scans for combined data
#' @return matrix with the row as scantime in second and column as m/z
#' @examples
#' \dontrun{
#' # mz100_200 and mz201_300 were the path to the raw data
#' matrix <- getmd(mz100_200,mz201_300)
#' }
#' @export
cbmd <- function(data1,
                 data2,
                 mzstep = 0.1,
                 rtstep = 0.01) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}
#' Get the MIR and related information from the files
#' @param file data file, CDF or other format supportted by xcmsRaw
#' @param mz1 the lowest mass
#' @param mz2 the highest mass
#' @return Molecular isotope ratio
#' @examples
#' \dontrun{
#' mr <- batch(data,mz1 = 79, mz2 = 81)
#' }
#' @export
batch <- function(file, mz1, mz2) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}

#' Get the MIR from the file
#' @param file data file, CDF or other format supportted by xcmsRaw
#' @param mz1 the lowest mass
#' @param mz2 the highest mass
#' @param rt a rough RT range contained only one peak to get the area
#' @param brt a rough RT range contained only one peak and enough noises to get the area
#' @return arearatio
#' @examples
#' \dontrun{
#' arearatio <- qbatch(datafile)
#' }
#' @export
qbatch <- function(file,
                   mz1,
                   mz2,
                   rt = c(8.65, 8.74),
                   brt = c(8.74, 8.85)) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}

#' Convert an XCMSnExp object to an mzrt S3 object.
#'
#' @noRd
.XCMSnExp2mzrt <-
        function(XCMSnExp,
                 method = "medret",
                 value = "into",
                 mzdigit = 4,
                 rtdigit = 1)
        {
                .Deprecated()
                message(
                        "This function has been deprecated and you might check documents of xcms to process your data."
                )
        }
#' Convert an xcmsSet object to an mzrt S3 object.
#'
#' @noRd
.xcmsSet2mzrt <-
        function(xcmsSet,
                 method = "medret",
                 value = "into",
                 mzdigit = 4,
                 rtdigit = 1)
        {
                .Deprecated()
                message(
                        "This function has been deprecated and you might check documents of xcms to process your data."
                )

        }

#' Get the mzrt profile and group information as a mzrt list and/or save them as csv or rds for further analysis.
#' @param xset xcmsSet/XCMSnExp objects
#' @param name file name for csv and/or eic file, default NULL
#' @param mzdigit m/z digits of row names of data frame, default 4
#' @param rtdigit retention time digits of row names of data frame, default 1
#' @param method parameter for groupval or featureDefinitions function, default medret
#' @param value parameter for groupval or featureDefinitions function, default into
#' @param eic logical, save xcmsSet and xcmsEIC objects for further investigation with the same name of files, you will need raw files in the same directory as defined in xcmsSet to extract the EIC based on the binned data. You could use `plot` to plot EIC for specific peaks. For example, `plot(xcmsEIC,xcmsSet,groupidx = 'M123.4567T278.9')` could show the EIC for certain peaks with m/z 206 and retention time 2789. default F
#' @param type csv format for further analysis, m means  Metaboanalyst, a means xMSannotator, p means Mummichog(NA values are imputed by `getimputation`, and F test is used here to generate stats and p value), o means full information csv (for `pmd` package), default o. mapo could output all those format files.
#' @return mzrt object, a list with mzrt profile and group information
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getmzrt(xset, name = 'demo', type = 'mapo')
#' }
#' @seealso \code{\link{getdata}},\code{\link{getdata2}}, \code{\link{getdoe}}, \code{\link{getcsv}}, \code{\link{getfilter}}
#' @references
#' Smith, C.A., Want, E.J., O’Maille, G., Abagyan, R., Siuzdak, G., 2006. XCMS: Processing Mass Spectrometry Data for Metabolite Profiling Using Nonlinear Peak Alignment, Matching, and Identification. Anal. Chem. 78, 779–787.
#' @export

getmzrt <-
        function(xset,
                 name = NULL,
                 mzdigit = 4,
                 rtdigit = 1,
                 method = "medret",
                 value = "into",
                 eic = FALSE,
                 type = 'o') {
                .Deprecated()
                message(
                        "This function has been deprecated and you might check documents of xcms to process your data."
                )
        }
#' Perform MS/MS dot product annotation for mgf file
#' @param file mgf file generated from MS/MS data
#' @param db database could be list object from `getMSP`
#' @param ppm mass accuracy, default 10
#' @param prems precursor mass range, default 1.1 to include M+H or M-H
#' @param binstep bin step for consin similarity
#' @param consinc consin similarity cutoff for annotation. Default 0.6.
#' @return list with MSMS annotation results
#' @export
dotpanno <- function(file,
                     db = NULL,
                     ppm = 10,
                     prems = 1.1,
                     binstep = 1,
                     consinc = 0.6) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}

#' Perform MS/MS X rank annotation for mgf file
#' @param file mgf file generated from MS/MS data
#' @param db database could be list object from `getms2pmd`
#' @param ppm mass accuracy, default 10
#' @param prems precursor mass range, default 1.1 to include M+H or M-H
#' @param intc intensity cutoff for peaks. Default 0.1
#' @param quantile X rank quantiles cutoff for annotation. Default 0.75.
#' @return list with MSMS annotation results
#' @export
xrankanno <- function(file,
                      db = NULL,
                      ppm = 10,
                      prems = 1.1,
                      intc = 0.1,
                      quantile = 0.75) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}

#' Show MS/MS pmd annotation result
#' @param anno list from MSMS anno function
#' @param ... other parameter for plot function
#' @return NULL
#' @export
plotanno <- function(anno, ...) {
        .Deprecated()
        message(
                "This function has been deprecated and you might check documents of xcms to process your data."
        )
}
