#' @export
runGloablStd <- function() {
    file <- system.file("shinyapp", "GlobalStd", "GlobalStd.Rmd", 
        package = "enviGCMS")
    if (file == "") {
        stop("Could not find directory. Try re-installing `enviGCMS`.", 
            call. = FALSE)
    }
    rmarkdown::run(file)
}
