#' Shiny application for Short-Chain Chlorinated Paraffins analysis
#' @export
runsccp <- function() {
        file <- system.file("shinyapps","sccp",
                            package = "enviGCMS")
        if (file == "") {
                stop("Could not find directory. Try re-installing `enviGCMS`.",
                     call. = FALSE)
        }
        shiny::runApp(file)
}
