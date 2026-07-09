#' Plot extracted ion chromatograms from raw mzML/mzXML data
#'
#' Reads mzML/mzXML raw data and plots extracted ion chromatograms (EIC) for a
#' list of target m/z values at MS1 or MS2 level. Uses \pkg{RaMS} for
#' lightweight data access and \pkg{plotly} for interactive visualization.
#'
#' @param filepath character. Path to an mzML or mzXML file.
#' @param featlist data.frame with columns: \code{name} (feature label),
#'   \code{mz} (target m/z), \code{ms_level} (either \code{"ms1"} or \code{"ms2"}).
#' @param diff numeric. Half-width of the m/z extraction window in Da. Default 0.005.
#'
#' @return A \pkg{plotly} object with overlaid extracted ion chromatograms.
#'
#' @examples
#' \dontrun{
#' featlist <- data.frame(
#'   name = c("target1", "target2"),
#'   mz = c(300.1234, 350.5678),
#'   ms_level = c("ms1", "ms2")
#' )
#' plotEIC("sample.mzML", featlist)
#' }
#' @export
plotEIC <- function(filepath, featlist, diff = 0.005) {
        if (!requireNamespace("RaMS", quietly = TRUE)) {
                stop("Package 'RaMS' is required for plotEIC. ",
                     "Install with: install.packages('RaMS')",
                     call. = FALSE)
        }
        if (!requireNamespace("plotly", quietly = TRUE)) {
                stop("Package 'plotly' is required for plotEIC.", call. = FALSE)
        }

        # Determine which MS levels are needed
        grab_what <- unique(ifelse(featlist$ms_level == "ms1", "MS1", "MS2"))
        msdata <- RaMS::grabMSdata(filepath, grab_what = grab_what, verbosity = 0)

        p <- plotly::plot_ly()

        for (i in seq_len(nrow(featlist))) {
                target_mz <- featlist$mz[i]
                mz_lo <- target_mz - diff
                mz_hi <- target_mz + diff

                if (featlist$ms_level[i] == "ms1" && !is.null(msdata$MS1)) {
                        dt <- msdata$MS1
                        sub <- dt[dt$mz >= mz_lo & dt$mz <= mz_hi, ]
                        # Sum intensities per RT scan within the mz window
                        eic <- stats::aggregate(int ~ rt, data = sub, FUN = sum)
                } else if (featlist$ms_level[i] == "ms2" && !is.null(msdata$MS2)) {
                        dt <- msdata$MS2
                        sub <- dt[dt$fragmz >= mz_lo & dt$fragmz <= mz_hi, ]
                        eic <- stats::aggregate(int ~ rt, data = sub, FUN = sum)
                } else {
                        next
                }

                if (nrow(eic) > 0) {
                        eic <- eic[order(eic$rt), ]
                        p <- plotly::add_lines(
                                p,
                                x = eic$rt,
                                y = eic$int,
                                name = paste0(round(target_mz, 4), "_",
                                              featlist$ms_level[i])
                        )
                }
        }
        p <- plotly::config(p, showTips = FALSE)
        return(p)
}


#' Extract top MS1 ions from MS2 EIC interactively
#'
#' A Shiny app that displays extracted ion chromatograms from MS2 data.
#' Click on a peak to select a retention time window, then extract the most
#' intense MS1 ions from that window and overlay their chromatograms.
#'
#' @param filepath character. Path to an mzML or mzXML file.
#' @param featlist data.frame with columns: \code{name}, \code{mz}, \code{ms_level}.
#' @param numTopIons integer. Number of most intense MS1 ions to extract. Default 10.
#' @param diff numeric. Half-width of m/z extraction window (Da). Default 0.01.
#' @param rtWindow numeric. Half-width of RT window (seconds) around clicked peak. Default 0.3.
#'
#' @return Opens an interactive Shiny app in the browser.
#'
#' @examples
#' \dontrun{
#' frags <- data.frame(
#'   name = c("Br79", "Br81"),
#'   mz = c(78.9183, 80.9163),
#'   ms_level = c("ms2", "ms2")
#' )
#' plotTopMS1Peaks("sample.mzML", frags, numTopIons = 3)
#' }
#' @export
plotTopMS1Peaks <- function(filepath, featlist, numTopIons = 10,
                            diff = 0.01, rtWindow = 0.3) {
        .check_eic_deps()
        .run_top_peaks_gadget(filepath, featlist, numTopIons, diff,
                              rtWindow, extract_level = 1)
}


#' Extract top MS2 ions from MS1 EIC interactively
#'
#' A Shiny app that displays extracted ion chromatograms from MS1 data.
#' Click on a peak to select a retention time window, then extract the most
#' intense MS2 ions from that window and overlay their chromatograms.
#'
#' @inheritParams plotTopMS1Peaks
#'
#' @return Opens an interactive Shiny app in the browser.
#'
#' @examples
#' \dontrun{
#' targets <- data.frame(
#'   name = c("precursor1"),
#'   mz = c(500.1234),
#'   ms_level = c("ms1")
#' )
#' plotTopMS2Peaks("sample.mzML", targets, numTopIons = 5)
#' }
#' @export
plotTopMS2Peaks <- function(filepath, featlist, numTopIons = 10,
                            diff = 0.01, rtWindow = 0.3) {
        .check_eic_deps()
        .run_top_peaks_gadget(filepath, featlist, numTopIons, diff,
                              rtWindow, extract_level = 2)
}


#' Check dependencies for EIC interactive functions
#' @noRd
.check_eic_deps <- function() {
        if (!requireNamespace("RaMS", quietly = TRUE)) {
                stop("Package 'RaMS' is required. ",
                     "Install with: install.packages('RaMS')",
                     call. = FALSE)
        }
        if (!requireNamespace("plotly", quietly = TRUE)) {
                stop("Package 'plotly' is required.", call. = FALSE)
        }
        if (!requireNamespace("shiny", quietly = TRUE)) {
                stop("Package 'shiny' is required.", call. = FALSE)
        }
}


#' Build EIC traces from RaMS data for a feature list
#' @noRd
.build_eic_plotly <- function(msdata, featlist, diff) {
        p <- plotly::plot_ly()
        for (i in seq_len(nrow(featlist))) {
                target_mz <- featlist$mz[i]
                mz_lo <- target_mz - diff
                mz_hi <- target_mz + diff

                if (featlist$ms_level[i] == "ms1" && !is.null(msdata$MS1)) {
                        dt <- msdata$MS1
                        sub <- dt[dt$mz >= mz_lo & dt$mz <= mz_hi, ]
                        eic <- stats::aggregate(int ~ rt, data = sub, FUN = sum)
                } else if (featlist$ms_level[i] == "ms2" && !is.null(msdata$MS2)) {
                        dt <- msdata$MS2
                        sub <- dt[dt$fragmz >= mz_lo & dt$fragmz <= mz_hi, ]
                        eic <- stats::aggregate(int ~ rt, data = sub, FUN = sum)
                } else {
                        next
                }

                if (nrow(eic) > 0) {
                        eic <- eic[order(eic$rt), ]
                        p <- plotly::add_lines(
                                p,
                                x = eic$rt,
                                y = eic$int,
                                name = paste0(round(target_mz, 4), "_",
                                              featlist$ms_level[i])
                        )
                }
        }
        p <- plotly::config(p, showTips = FALSE)
        plotly::layout(p, showlegend = TRUE)
}


#' Internal: Run the interactive peak extraction Shiny app
#'
#' Shared implementation for plotTopMS1Peaks and plotTopMS2Peaks.
#'
#' @param filepath,featlist,numTopIons,diff,rtWindow See public functions.
#' @param extract_level integer. 1 to extract MS1 from MS2 clicks,
#'   2 to extract MS2 from MS1 clicks.
#' @noRd
.run_top_peaks_gadget <- function(filepath, featlist, numTopIons, diff,
                                  rtWindow, extract_level) {

        extract_label <- paste0("MS", extract_level)

        ui <- shiny::fluidPage(
                shiny::titlePanel(
                        paste0("Click on peaks to select RT, then click Get ",
                               extract_label)
                ),
                shiny::sidebarLayout(
                        shiny::sidebarPanel(
                                width = 2,
                                shiny::textOutput("rtselect"),
                                shiny::br(),
                                shiny::textOutput("rtWindow"),
                                shiny::br(),
                                shiny::actionButton("sync", "Sync ranges"),
                                shiny::br(), shiny::br(),
                                shiny::actionButton("getMS",
                                                    paste0("Get ", extract_label)),
                                shiny::br(), shiny::br(),
                                shiny::actionButton("done", "Done")
                        ),
                        shiny::mainPanel(
                                width = 10,
                                plotly::plotlyOutput("plot1", height = "400px"),
                                plotly::plotlyOutput("plot2", height = "400px")
                        )
                )
        )

        server <- function(input, output, session) {

                # Read all MS data once at startup
                msdata <- RaMS::grabMSdata(filepath,
                                           grab_what = c("MS1", "MS2"),
                                           verbosity = 0)

                rtr <- NULL

                output$plot1 <- plotly::renderPlotly({
                        .build_eic_plotly(msdata, featlist, diff)
                })

                shiny::observeEvent(plotly::event_data("plotly_click"), {
                        click_rt <- as.data.frame(
                                plotly::event_data("plotly_click")
                        )[[3]]
                        rtr <<- c(round(click_rt - rtWindow, 1),
                                  round(click_rt + rtWindow, 1))
                        output$rtselect <- shiny::renderText(
                                paste0("RT:", rtr[1], "-", rtr[2])
                        )
                })

                shiny::observeEvent(input$getMS, {
                        if (is.null(rtr)) return()

                        # Get the target MS level data within the RT window
                        if (extract_level == 1 && !is.null(msdata$MS1)) {
                                dt <- msdata$MS1
                                dt_rt <- dt[dt$rt >= rtr[1] & dt$rt <= rtr[2], ]
                                top_ions <- dt_rt[order(dt_rt$int,
                                                        decreasing = TRUE), ]
                                top_ions <- utils::head(
                                        top_ions[!duplicated(round(top_ions$mz, 4)), ],
                                        numTopIons
                                )
                                ms_tag <- "ms1"
                        } else if (extract_level == 2 && !is.null(msdata$MS2)) {
                                dt <- msdata$MS2
                                dt_rt <- dt[dt$rt >= rtr[1] & dt$rt <= rtr[2], ]
                                top_ions <- dt_rt[order(dt_rt$int,
                                                        decreasing = TRUE), ]
                                top_ions <- utils::head(
                                        top_ions[!duplicated(
                                                round(top_ions$fragmz, 4)), ],
                                        numTopIons
                                )
                                ms_tag <- "ms2"
                        } else {
                                return()
                        }

                        if (nrow(top_ions) == 0) return()

                        mz_col <- if (extract_level == 1) "mz" else "fragmz"
                        featlist2 <- data.frame(
                                name = paste0("S", seq_len(nrow(top_ions))),
                                mz = top_ions[[mz_col]],
                                ms_level = ms_tag,
                                stringsAsFactors = FALSE
                        )
                        featlist2 <- rbind(featlist2, featlist)

                        output$plot2 <- plotly::renderPlotly({
                                .build_eic_plotly(msdata, featlist2, diff)
                        })
                })

                shiny::observeEvent(plotly::event_data("plotly_relayout"), {
                        rtsel <- plotly::event_data("plotly_relayout")
                        output$rtWindow <- shiny::renderText(
                                paste0("RT range: ", rtsel[[1]], " - ",
                                       rtsel[[2]])
                        )
                })

                shiny::observeEvent(input$sync, {
                        rtsel <- plotly::event_data("plotly_relayout")
                        output$plot1 <- plotly::renderPlotly({
                                p <- .build_eic_plotly(msdata, featlist, diff)
                                plotly::layout(p, xaxis = list(
                                        range = c(rtsel[[1]], rtsel[[2]])
                                ))
                        })
                })

                shiny::observeEvent(input$done, {
                        shiny::stopApp()
                })
        }

        app <- shiny::shinyApp(ui, server)
        shiny::runApp(app, launch.browser = TRUE)
}
