#' Get the exact mass of the isotopologues from a chemical formula or reaction's isotope patterns with the highest abundances
#' @param data a chemical formula or reaction e.g. 'Cl-H', 'C2H4'
#' @return numerical vector
#' @examples
#' \dontrun{
#' getmass('CH2')
#' }
#' @export
getmass <- function(data) {
        if (grepl('-', data)) {
                name <- unlist(strsplit(data, '-'))
                iso1 <-
                        rcdk::get.isotopes.pattern(rcdk::get.formula(name[1]))
                iso2 <-
                        rcdk::get.isotopes.pattern(rcdk::get.formula(name[2]))
                cus <-
                        as.numeric(iso1[max(iso1[, 2]), 1]) - as.numeric(iso2[max(iso2[, 2]), 1])
        } else{
                iso <- rcdk::get.isotopes.pattern(rcdk::get.formula(data))
                cus <-
                        as.numeric(iso[max(iso[, 2]), 1])
        }
        return(cus)
}
#' Get the Relative Mass Defect
#' @param mz numeric vector for exact mass
#' @return Relative Mass Defect
#' @examples
#' \dontrun{
#' getrmd(getmass('C2H4'))
#' }
#' @export
getrmd <- function(mz){
        rmd <- round((round(mz) - mz) / mz * 10 ^ 6)
        return(rmd)
}
#' Get the raw Mass Defect
#' @param mz numeric vector for exact mass
#' @return raw Mass Defect
#' @examples
#' \dontrun{
#' getmdr(getmass('C2H4'))
#' }
#' @export
getmdr <- function(mz){
        md <- round((round(mz) - mz) * 10 ^ 3)
        return(md)
}
#' Get the high order unit based Mass Defect
#' @param mz numeric vector for exact mass
#' @param cus chemical formula or reaction
#' @param method you could use `round`, `floor` or `ceiling`
#' @return high order Mass Defect with details
#' @examples
#' \dontrun{
#' getmdh(getmass('C2H4'))
#' }
#' @export
getmdh <- function(mz,cus = c('CH2,H2'), method = 'round'){
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
        for (i in 1:length(temp)) {
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
                re <- cbind.data.frame(mz,MD1,MD2)

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
                re <- cbind.data.frame(mz,MD1,MD2)

        } else if (method == 'ceiling'){
                MD1 <-
                        round(ceiling(signif(omd)) - omd,
                              digits = 6)
                md2 <- round(ceiling(signif(sumd))- sumd,
                               digits = 6)
                smd <-  MD1 / md2
                MD2 <-
                        round(ceiling(signif(smd)) - smd,
                              digits = 6)
                re <- cbind.data.frame(mz,MD1,MD1)
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
                        re <- cbind.data.frame(mz,MD1,MD2,MD3)
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
                        re <- cbind.data.frame(mz,MD1,MD2,MD3)
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
                        re <- cbind.data.frame(mz,MD1,MD2,MD3)
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
                        re <- cbind.data.frame(mz,MD1,MD2,MD3)
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
                        re <- cbind.data.frame(mz,MD1,MD2,MD3)
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
                        re <- cbind.data.frame(mz,MD1,MD2,MD3)
                }
        } else{

                if (method == 'round') {
                        omd <- mz * round(cus) / cus
                        MD1 <-
                                round(round(omd) - omd,
                                      digits = 6)
                        re <- cbind.data.frame(mz,MD1)
                } else if (method == 'floor') {
                        omd <- mz * floor(cus) / cus
                        MD1 <-
                                round(floor(omd) - omd,
                                      digits = 6)
                        re <- cbind.data.frame(mz,MD1)
                } else{
                        omd <- mz * ceiling(cus) / cus
                        MD1 <-
                                round(ceiling(omd) - omd,
                                      digits = 6)
                        re <- cbind.data.frame(mz,MD1)
                }
        }
        return(re)
}
