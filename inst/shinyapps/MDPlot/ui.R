library(shiny)
library(dplyr)
library(shinythemes)
library(DT)
library(plotly)
library(crosstalk)
library(shinyjs)
library(rcdk)

# UI function -------------------------------------------------------------

shinyUI(navbarPage(
        useShinyjs(),
        # Include shinyjs
        theme = shinytheme('spacelab'),
        tabPanel(
                "MDPlotR: interactive mass defect plots",
                fluidPage(sidebarLayout(
                        sidebarPanel(
                                fileInput(
                                        'file1',
                                        'Choose CSV File',
                                        accept = c('text/csv',
                                                   'text/comma-separated-values,text/plain',
                                                   '.csv')
                                ),
                                fluidRow(column(
                                        6,
                                        textInput("cus1", "MD formula 1", value = 'CH2,H2')
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = "mdr1",
                                                label = "Rounding 1",
                                                choices = c("round", "floor", "ceiling"),
                                                selected = "round"
                                        )
                                )),

                                fluidRow(column(
                                        6,
                                        textInput("cus2", "MD formula 2", value = 'Cl-H')
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = "mdr2",
                                                label = "Rounding 2",
                                                choices = c("round", "floor", "ceiling"),
                                                selected = "round"
                                        )
                                )),
                                actionButton('go', 'Plot', width = '100%'),
                                br(),
                                radioButtons(
                                        inputId = "single",
                                        label = "Single or Double plots",
                                        choices = c("Single", "Double"),
                                        selected = "Single",
                                        inline = TRUE
                                ),
                                uiOutput("plotctr"),
                                uiOutput("plotctr2"),
                                uiOutput("slide1"),
                                uiOutput("slide2"),
                                uiOutput("slide3"),
                                checkboxInput('ins', 'Show intensity as size', F),
                                width = 3
                        ),
                        mainPanel(
                                uiOutput("plot"),
                                DTOutput("x1"),
                                fluidRow(column(
                                        3, downloadButton("x3", "Download Filtered Data")
                                )),
                                tags$br(),

                                # Using Shinyjs to open websites
                                fluidRow(
                                        h4("Links to web tools for compound search"),
                                        column(
                                                3,
                                                align = "left",
                                                actionButton("open_1", "Chemistry Dashboard",
                                                             onclick =
                                                                     "window.open('https://comptox.epa.gov/dashboard/dsstoxdb/advanced_search')")
                                        ),
                                        column(
                                                3,
                                                align = "left",
                                                actionButton("open_2", "ChemSpider",
                                                             onclick =
                                                                     "window.open('http://www.chemspider.com/FullSearch.aspx')")
                                        ),
                                        column(
                                                2,
                                                align = "left",
                                                actionButton("open_3", "EnviPat",
                                                             onclick =
                                                                     "window.open('http://www.envipat.eawag.ch')")
                                        ),
                                        column(
                                                2,
                                                align = "left",
                                                actionButton("open_4", "EnviHomolog",
                                                             onclick =
                                                                     "window.open('http://www.envihomolog.eawag.ch')")
                                        ),
                                        column(
                                                2,
                                                align = "left",
                                                actionButton("open_5", "Norman Massbank",
                                                             onclick =
                                                                     "window.open('https://massbank.eu/MassBank/QuickSearch.html')")
                                        )


                                )
                        )
                ))
        ),
        tabPanel(
                "Instructions",
                h4("Data uploading"),
                p(
                        "Uploaded csv files should contain three columns with name 'mz','rt', 'intensity' and contain mass to charge(m/z), retention time and intensity data."
                ),
                p(
                        "After you uploaded the csv data, input your mass defect base in the input box(es) and click plot to show the MD plots. When you make changes on the left panel, you need to click plot to update the plot. However, you could explore interactively on the plot and table."
                ),
                h4('Equation'),
                p("Mass defect = measured mass - round(measured mass) "),
                p(
                        "Relative Mass defect = (measured mass - round(measured mass))/measured mass * 10^6 "
                ),
                h5('Unit based first order mass defect'),
                p(
                        "first-order mass = measured mass * round(first-order unit exact mass)/first order unit exact mass"
                ),
                br(),
                p(
                        "first order mass defect = first-order mass - round/floor/ceiling(first-order mass) "
                ),
                br(),
                p(
                        "second-order mass = first-order mass defect(unit 1)/first-order mass defect(unit 2)"
                ),
                br(),
                p(
                        "second-order mass defect = second-order mass - round/floor/ceiling(second-order mass)"
                ),
                br(),
                p(
                        "third-order mass = second-order mass defect(unit 1, unit 2)/second-order mass defect(unit 2,unit3) "
                ),
                br(),
                p(
                        "third-order mass defect = third-order mass - round/floor/ceiling(third-order mass)"
                ),

                h4('Chemical formula'),
                p("'CH2' means unit"),
                br(),
                p(
                        "'Br-H' means add Br atom and remove H atom, use minus sign to seperate them. This app only support two different unit."
                ),
                br(),
                p(
                        "'CH2,H2' means 'CH2' is the first-order mass defect unit and 'H2' is the second-order mass defect, use comma to seperate them. This app only support at most three-order mass defect."
                ),
                br(),
                p(
                        "'CH2,Br-H' means 'CH2' is the first-order mass defect unit and 'Br-H' is the second-order mass defect, use comma to seperate them."
                )
        )
)
)
