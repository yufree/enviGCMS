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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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

#' Plot the influnces of DoE and Batch effects on each peaks
#' @param df data output from `svacor` function
#' @param dfsv data output from `svaplot` function for corrected data
#' @param dfanova data output from `svaplot` function for raw data
#' @return influnces plot
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
#' @export
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
#' @export
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
#' @return data, corrected data, mz and retention for fileted data
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
#' @export
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
#' @export
svapca <- function(list,
                   center = T,
                   scale = T,
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
#' @export
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
#' @return csv files for both raw and corrected data for metabolanalyst if SVA could be applied
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
#' @export
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
#' @param type m means  Metaboanalyst, a means xMSannotator, o means full infomation csv
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
#' @export
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
#' @export
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
#' @export
getupload3 <- function(list, name = "Peaklist") {
        .Deprecated()
        message("This function has been deprecated and you could use getupload to get csv file.")
}
#' Get the mzrt profile and group information for batch correction and plot as a list for xcms 3 object
#' @param xset a XCMSnExp object with processed data
#' @param name file name for csv file, default NULL
#' @return list with rtmz profile and group infomation
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
#' @export
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
#' @export
plote <- function(xset,
                  name = "test",
                  test = "t",
                  nonpara = "n",
                  ...) {
        .Deprecated()
        message(
                "This function will be deprecated and you could use getmzrt to get related object to plot EIC."
        )
        gt <- xcms::groups(xset)
        a <-
                xcms::diffreport(
                        xset,
                        filebase = name,
                        eicmax = nrow(gt),
                        nonpara = nonpara,
                        ...
                )
        return(a)
}

#' Get the features from t test, with p value, q value, rsd and power restriction
#' @param list list with data as peaks list, mz, rt and group information (two groups)
#' @param power defined power
#' @param pt p value threshold
#' @param qt q value threshold, BH adjust
#' @param n sample numbers in one group
#' @param imputation parameters for `getimputation` function method
#' @return dataframe with peaks fit the setting above
#' @export
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
#' @export
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
