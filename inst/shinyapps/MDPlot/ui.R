library(shiny)
library(dplyr)
library(shinythemes)
library(DT)
library(plotly)
library(crosstalk)
library(Rdisop)

shinyUI(navbarPage(
        "MDPlotR: interactive mass defect plots",
        theme = shinytheme('spacelab'),
        tabPanel("Get Started",
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
                                        12,
                                        textInput("cus1", "MD formula 1", value = 'CH2,O')
                                )

                                ),

                                fluidRow(column(
                                        12,
                                        textInput("cus2", "MD formula 2", value = 'Cl-H')
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
                                checkboxInput('ins', 'Show intensity as size', F),
                                checkboxInput("show_leg", "Show plot legends", T),
                                uiOutput("plotctr"),
                                uiOutput("plotctr2"),
                                uiOutput("slide1"),
                                uiOutput("slide2"),
                                uiOutput("slide3"),
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
                                                shiny::a(h4("Chemistry Dashboard", class = "btn btn-default action-button" ,
                                                            style = "fontweight:600"), target = "_blank",
                                                         href = "https://comptox.epa.gov/dashboard/dsstoxdb/advanced_search")
                                        ),
                                        column(
                                                3,
                                                align = "left",
                                                shiny::a(h4("ChemSpider", class = "btn btn-default action-button" ,
                                                            style = "fontweight:600"), target = "_blank",
                                                         href = "http://www.chemspider.com/FullSearch.aspx")),
                                        column(
                                                        2,
                                                        align = "left",
                                                        shiny::a(h4("EnviPat", class = "btn btn-default action-button" ,
                                                                    style = "fontweight:600"), target = "_blank",
                                                                 href = "http://www.envipat.eawag.ch")),
                                        column(
                                                2,
                                                align = "left",
                                                shiny::a(h4("EnviHomolog", class = "btn btn-default action-button" ,
                                                            style = "fontweight:600"), target = "_blank",
                                                         href = "http://www.envihomolog.eawag.ch")),
                                        column(
                                                2,
                                                align = "left",
                                                shiny::a(h4("Norman Massbank", class = "btn btn-default action-button" ,
                                                            style = "fontweight:600"), target = "_blank",
                                                         href = "https://massbank.eu/MassBank/QuickSearch.html"))
                                        )


                                )
                        )
                )),
        tabPanel(
                "Instructions",
                sidebarLayout(
                        sidebarPanel(h3("Table of content"),
                                h4("File input"),
                                h4("Equation"),
                                     width = 3),
                        mainPanel(
                withMathJax(includeMarkdown("instructions.md"))
                        )
                )
        )
)
)
