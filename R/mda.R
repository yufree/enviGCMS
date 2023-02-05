#' Demo data for TBBPA metabolism in Pumpkin
#' @docType data
#' @usage data(TBBPA)
#' @format A list object with data, mass to charge ratio, retention time and group information. Three pumpkin seeding root samples' peaks list is extracted by xcms online.
#' @references Hou, X., Yu, M., Liu, A., Wang, X., Li, Y., Liu, J., Schnoor, J.L., Jiang, G., 2019. Glycosylation of Tetrabromobisphenol A in Pumpkin. Environ. Sci. Technol. https://doi.org/10.1021/acs.est.9b02122
"TBBPA"

#' Get mass defect with certain scaled factor
#' @param mass vector of mass
#' @param sf scaled factors
#' @return dataframe with mass, scaled mass and scaled mass defect
#' @examples
#' mass <- c(100.1022,245.2122,267.3144,400.1222,707.2294)
#' sf <- 0.9988
#' mf <- getmassdefect(mass,sf)
#' @seealso \code{\link{plotkms}}
#' @export

getmassdefect <- function(mass, sf) {
        sm <- mass * sf
        sd <- round(sm) - sm
        df <- as.data.frame(cbind(mass, sm, sd))
        graphics::plot(df$sd ~ df$sm, xlab = "m/z", ylab = "scaled MD")
        return(df)
}
#' plot the kendrick mass defect diagram
#' @param data vector with the name m/z
#' @param cutoff remove the low intensity
#' @return NULL
#' @seealso \code{\link{getmassdefect}}
#' @examples
#' \dontrun{
#' mz <- c(10000,5000,20000,100,40000)
#' names(mz) <- c(100.1022,245.2122,267.3144,400.1222,707.2294)
#' plotkms(mz)
#' }
#' @export
plotkms <- function(data, cutoff = 1000) {
        data <- data[data > cutoff]
        mz <- as.numeric(names(data))
        km <- mz * 14 / 14.01565
        kmd <- round(km) - km
        graphics::smoothScatter(kmd ~ round(km), xlab = "Kendrick nominal mass",
                                ylab = "Kendrick mass defect")
}
#' Get the exact mass of the isotopologues from a chemical formula or reaction's isotope patterns with the highest abundances
#' @param data a chemical formula or reaction e.g. 'Cl-H', 'C2H4'
#' @return numerical vector
#' @examples
#' getmass('CH2')
#' @export
getmass <- function(data) {
        if (grepl('-', data)) {
                name <- unlist(strsplit(data, '-'))
                if (name[1] == '') {
                        table2 <- Rdisop::getMolecule(name[2])$isotopes[[1]]
                        iso2 <-
                                table2[1, which.max(table2[2, ])]
                        iso <- 0 - iso2
                } else{
                        table1 <- Rdisop::getMolecule(name[1])$isotopes[[1]]
                        iso1 <-
                                table1[1, which.max(table1[2, ])]
                        table2 <-
                                Rdisop::getMolecule(name[2])$isotopes[[1]]
                        iso2 <-
                                table2[1, which.max(table2[2, ])]
                        iso <-
                                iso1 - iso2
                }

        } else{
                table <- Rdisop::getMolecule(data)$isotopes[[1]]
                iso <- table[1, which.max(table[2, ])]
        }
        return(iso)
}
#' Get the Relative Mass Defect
#' @param mz numeric vector for exact mass
#' @return Relative Mass Defect
#' @examples
#' getrmd(getmass('C2H4'))
#' @export
getrmd <- function(mz) {
        rmd <- round((round(mz) - mz) / mz * 10 ^ 6)
        return(rmd)
}
#' Get the raw Mass Defect
#' @param mz numeric vector for exact mass
#' @return raw Mass Defect
#' @examples
#' getmdr(getmass('C2H4'))
#' @export
getmdr <- function(mz) {
        md <- round((round(mz) - mz) * 10 ^ 3)
        return(md)
}
#' Get the high order unit based Mass Defect
#' @param mz numeric vector for exact mass
#' @param cus chemical formula or reaction
#' @param method you could use `round`, `floor` or `ceiling`
#' @return high order Mass Defect with details
#' @examples
#' getmdh(getmass('C2H4'))
#' @export
getmdh <- function(mz,
                   cus = c('CH2,H2'),
                   method = 'round') {
        getorder <- function(input) {
                if (grepl(',', input)) {
                        name <- unlist(strsplit(input, ','))
                } else{
                        name <- input
                }
                return(name)
        }
        temp <- getorder(cus)
        cus <- NULL
        for (i in seq_along(temp)) {
                cus <- c(cus, getmass(temp[i]))
        }
        if (length(cus) == 2) {
                omd <- mz * round(cus[1]) / cus[1]
                sumd <- cus[2] * round(cus[1]) / cus[1]

                if (method == 'round') {
                        MD1 <-
                                round(round(omd) - omd,
                                      digits = 6)
                        md2 <- round(round(sumd) - sumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        MD2 <-
                                round(round(smd) - smd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD2)

                } else if (method == 'floor') {
                        MD1 <-
                                round(floor(omd) - omd,
                                      digits = 6)
                        md2 <- round(floor(sumd) - sumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        MD2 <-
                                round(floor(smd) - smd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD2)

                } else if (method == 'ceiling') {
                        MD1 <-
                                round(ceiling(signif(omd)) - omd,
                                      digits = 6)
                        md2 <- round(ceiling(signif(sumd)) - sumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        MD2 <-
                                round(ceiling(signif(smd)) - smd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD1)
                }
        } else if (length(cus) == 3) {
                omd <- mz * round(cus[1]) / cus[1]
                sumd <- cus[2] * round(cus[1]) / cus[1]
                tumd <- cus[3] * round(cus[1]) / cus[1]

                if (method == 'round') {
                        MD1 <-
                                round(round(omd) - omd,
                                      digits = 6)
                        md2 <- round(round(sumd) - sumd,
                                     digits = 6)
                        md3 <- round(round(tumd) - tumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        tsmd <- md3 / md2
                        MD2 <-
                                round(round(smd) - smd,
                                      digits = 6)
                        md3 <- round(round(tsmd) - tsmd,
                                     digits = 6)
                        tmd <- MD2 / md3
                        MD3 <-
                                round(round(tmd) - tmd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD2, MD3)
                } else if (method == 'floor') {
                        MD1 <-
                                round(floor(omd) - omd,
                                      digits = 6)
                        md2 <- round(floor(sumd) - sumd,
                                     digits = 6)
                        md3 <- round(floor(tumd) - tumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        tsmd <- md3 / md2
                        MD2 <-
                                round(floor(smd) - smd,
                                      digits = 6)
                        md3 <- round(floor(tsmd) - tsmd,
                                     digits = 6)
                        tmd <- MD2 / md3
                        MD3 <-
                                round(floor(tmd) - tmd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD2, MD3)
                } else{
                        MD1 <-
                                round(ceiling(omd) - omd,
                                      digits = 6)
                        md2 <- round(ceiling(sumd) - sumd,
                                     digits = 6)
                        md3 <- round(ceiling(tumd) - tumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        tsmd <- md3 / md2
                        MD2 <-
                                round(ceiling(smd) - smd,
                                      digits = 6)
                        md3 <- round(ceiling(tsmd) - tsmd,
                                     digits = 6)
                        tmd <- MD2 / md3
                        MD3 <-
                                round(ceiling(tmd) - tmd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD2, MD3)
                }

        } else if (length(cus) > 3) {
                message("Sorry, only the first three unit would be used.")
                omd <- mz * round(cus[1]) / cus[1]
                sumd <- cus[2] * round(cus[1]) / cus[1]
                tumd <- cus[3] * round(cus[1]) / cus[1]

                if (method == 'round') {
                        MD1 <-
                                round(round(omd) - omd,
                                      digits = 6)
                        md2 <- round(round(sumd) - sumd,
                                     digits = 6)
                        md3 <- round(round(tumd) - tumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        tsmd <- md3 / md2
                        MD2 <-
                                round(round(smd) - smd,
                                      digits = 6)
                        md3 <- round(round(tsmd) - tsmd,
                                     digits = 6)
                        tmd <- MD2 / md3
                        MD3 <-
                                round(round(tmd) - tmd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD2, MD3)
                } else if (method == 'floor') {
                        MD1 <-
                                round(floor(omd) - omd,
                                      digits = 6)
                        md2 <- round(floor(sumd) - sumd,
                                     digits = 6)
                        md3 <- round(floor(tumd) - tumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        tsmd <- md3 / md2
                        MD2 <-
                                round(floor(smd) - smd,
                                      digits = 6)
                        md3 <- round(floor(tsmd) - tsmd,
                                     digits = 6)
                        tmd <- MD2 / md3
                        MD1_3 <-
                                round(floor(tmd) - tmd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD2, MD3)
                } else{
                        MD1 <-
                                round(ceiling(omd) - omd,
                                      digits = 6)
                        md2 <- round(ceiling(sumd) - sumd,
                                     digits = 6)
                        md3 <- round(ceiling(tumd) - tumd,
                                     digits = 6)
                        smd <-  MD1 / md2
                        tsmd <- md3 / md2
                        MD2 <-
                                round(ceiling(smd) - smd,
                                      digits = 6)
                        md3 <- round(ceiling(tsmd) - tsmd,
                                     digits = 6)
                        tmd <- MD2 / md3
                        MD3 <-
                                round(ceiling(tmd) - tmd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1, MD2, MD3)
                }
        } else{
                if (method == 'round') {
                        omd <- mz * round(cus) / cus
                        MD1 <-
                                round(round(omd) - omd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1)
                } else if (method == 'floor') {
                        omd <- mz * floor(cus) / cus
                        MD1 <-
                                round(floor(omd) - omd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1)
                } else{
                        omd <- mz * ceiling(cus) / cus
                        MD1 <-
                                round(ceiling(omd) - omd,
                                      digits = 6)
                        re <- cbind.data.frame(mz, MD1)
                }
        }
        return(re)
}

#' Screen organohalogen compounds by retention time, mass defect analysis and isotope relationship modified by literature report. Also support compounds with [M] and [M+2] ratio cutoff.
#' @param list list with data as peaks list, mz, rt and group information, retention time should be in seconds
#' @param sf scale factor, default 78/77.91051(Br)
#' @param step mass defect step, default 0.001
#' @param stepsd1 mass defect uncertainty for lower mass, default 0.003
#' @param stepsd2 mass defect uncertainty for higher mass, default 0.005
#' @param mzc threshold of lower mass and higher mass, default 700
#' @param cutoffint the cutoff of intensity, default 1000
#' @param cutoffr the cutoff of [M] and [M+2] ratio, default 0.4
#' @param clustercf the cutoff of cluster analysis to separate two different ions groups for retention time, default 10
#' @return list with filtered organohalogen compounds
#' @references Identification of Novel Brominated Compounds in Flame Retarded Plastics Containing TBBPA by Combining Isotope Pattern and Mass Defect Cluster Analysis Ana Ballesteros-Gómez, Joaquín Ballesteros, Xavier Ortiz, Willem Jonker, Rick Helmus, Karl J. Jobst, John R. Parsons, and Eric J. Reiner Environmental Science & Technology 2017 51 (3), 1518-1526 DOI: 10.1021/acs.est.6b03294
#' @export
findohc <-
        function(list,
                 sf = 78 / 77.91051,
                 step = 0.001,
                 stepsd1 = 0.003,
                 stepsd2 = 0.005,
                 mzc = 700,
                 cutoffint = 1000,
                 cutoffr = 0.4,
                 clustercf = 10) {
                mz <- list$mz
                ins <- apply(list$data, 1, mean, na.rm = TRUE)
                rt <- list$rt
                mzr <- round(mz)
                sm <- mz * sf
                sd <- ceiling(sm) - sm
                smsd <- ifelse(mz <= mzc, stepsd1, stepsd2)
                smstep <- seq(0, 1, step)
                data <-
                        cbind.data.frame(
                                mz = mz,
                                mzr = mzr,
                                sm = sm,
                                sd = sd,
                                ins = ins,
                                rt = rt
                        )

                result <- NULL
                for (i in seq_along(smstep)) {
                        mini <- smstep[i] - smsd
                        maxi <- smstep[i] + smsd
                        index <- sd < maxi & sd > mini

                        li <- data[index & ins > cutoffint, ]
                        mzt <- mzr[index & ins > cutoffint]
                        rtt <- rt[index & ins > cutoffint]
                        #dist(mzt) <-
                        if (length(mzt) >= 2) {
                                #c <- stats::cutree(stats::hclust(stats::dist(mzt)), h = clustercf)
                                t <-
                                        stats::cutree(stats::hclust(stats::dist(rtt)), h = clustercf)
                                # u <- paste0(c, t)
                                # cn <- length(unique(u))
                                # lit <- cbind.data.frame(li, u, i)
                                #
                                # for (j in 1:cn) {
                                lit <- cbind.data.frame(li, t, i)
                                for (j in seq_along(unique(t))) {
                                        li2 <- lit[lit[, 7] == j, ]
                                        mzt2 <-
                                                lit$mzr[lit[, 7] == j]
                                        if (length(mzt2) >= 2) {
                                                if (length(unique(li2$ins)) > 1) {
                                                        ratio <- max(li2$ins[li2$ins != max(li2$ins)]) / max(li2$ins)
                                                        diff <-
                                                                abs(li2$mzr[round(li2$ins) == round(max(li2$ins[li2$ins != max(li2$ins)]))] - li2$mzr[which.max(li2$ins)])
                                                } else{
                                                        ratio <- 1
                                                        diff <-
                                                                abs(li2$mzr[1] - li2$mzr[2])
                                                }

                                                if (ratio > cutoffr &
                                                    round(diff) == 2) {
                                                        li2 <- cbind.data.frame(li2, ratio)
                                                        result <-
                                                                as.data.frame(rbind(result, li2))
                                                }
                                        }
                                }
                        }
                }
                list$ohc <- result[!duplicated(result$mz),]
                return(list)
        }

#' Screen metabolites by Mass Defect
#' @param list list with data as peaks list, mz, rt and group information, retention time should be in seconds
#' @param mass mass to charge ratio of specific compounds
#' @param mdr mass defect range, default 50mDa
#' @return list with filtered metabolites mass to charge index of certain compound
#' @export
findmet <-
        function(list, mass, mdr = 50) {
                rmd <- (round(mass) - mass) * 1000
                rmdall <- (round(list$mz) - list$mz) * 1000

                if (length(mass) > 1) {
                        metindex <- NULL
                        for (i in seq_along(mass)) {
                                metindexi <- rmdall > rmd[i] - mdr & rmdall < rmd[i] + mdr
                                metindex <-
                                        cbind(metindex, metindexi)
                        }
                        colnames(metindex) <- mass
                        list$metindex <- metindex
                        return(list)
                } else{
                        list$metindex <- rmdall > rmd - mdr & rmdall < rmd + mdr
                        return(list)
                }

        }

#' Find lipid class of metabolites base on referenced Kendrick mass defect
#' @param list list with data as peaks list, mz, rt and group information, retention time should be in seconds
#' @param mode 'pos' for positive mode, 'neg' for negative mode and 'none' for neutral mass, only support [M+H] and [M-H] for each mode
#' @return list list with dataframe with the lipid referenced Kendrick mass defect(RKMD) and logical for class
#' @references Method for the Identification of Lipid Classes Based on Referenced Kendrick Mass Analysis. Lerno LA, German JB, Lebrilla CB. Anal Chem. 2010 May 15;82(10):4236–45.
#' @examples
#' data(list)
#' RKMD <- findlipid(list)
#' @export
findlipid <-
        function(list, mode = 'pos') {
                if (mode == 'pos' | mode == 'neg') {
                        km <-
                                (list$mz * 14 / 14.01565 - floor(list$mz * 14 / 14.01565)) / 0.0134
                } else{
                        adduct <- ifelse(mode == 'pos', 1.008, -1.008)
                        km <-
                                (((list$mz + adduct) * 14 / 14.01565) - floor((list$mz + adduct) * 14 / 14.01565)) /
                                0.0134
                }
                TAG <- 0.8355
                DAG <- 0.8719
                MAG <- 0.9082
                PC <- 0.7430
                PE <- 0.7430
                PS <- 0.6771
                PI <- 0.6141
                PtdG <- 0.6963
                PtdH <- 0.7422
                class <- c(TAG, DAG, MAG, PC, PE, PS, PI, PtdG, PtdH) /
                        0.0134
                t <- as.data.frame(outer(km, class, `-`))
                s <- t - floor(t)
                c <-
                        apply(s, 2, function(x)
                                round(x %% 1, digits = 2) == 0)
                colnames(t) <-
                        c(
                                'TAG_RKMD',
                                'DAG_RKMD',
                                'MAG_RKMD',
                                'PC_RKMD',
                                'PE_RKMD',
                                'PS_RKMD',
                                'PI_RKMD',
                                'PtdG_RKMD',
                                'PtdH_RKMD'
                        )
                colnames(c) <-
                        c('TAG',
                          'DAG',
                          'MAG',
                          'PC',
                          'PE',
                          'PS',
                          'PI',
                          'PtdG',
                          'PtdH')
                t$mz <- list$mz
                list$RKMD <- cbind.data.frame(t, c)
                return(list)
        }

#' Find PFCs based on mass defect analysis
#' @param list list with data as peaks list, mz, rt and group information, retention time should be in seconds
#' @return list list with potential PFCs compounds index
#' @references Liu, Y.; D’Agostino, L. A.; Qu, G.; Jiang, G.; Martin, J. W. High-Resolution Mass Spectrometry (HRMS) Methods for Nontarget Discovery and Characterization of Poly- and per-Fluoroalkyl Substances (PFASs) in Environmental and Human Samples. TrAC Trends in Analytical Chemistry 2019, 121, 115420.
#' @examples
#' data(list)
#' pfc <- findpfc(list)
#' @export
findpfc <- function(list) {
        md <- list$mz - round(list$mz)
        # mass defect -0.1 to 0.15
        list$pfc <- md > -0.1 & md < 0.15
        return(list)
}
