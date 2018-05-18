library(shiny)
library(dplyr)
library(shinythemes)
library(DT)
library(plotly)
library(crosstalk)
library(shinyjs)
library(rcdk)

shinyServer(function(input, output, session) {
        MD_data <- reactive({
                #  require that the input is available
                req(input$file1)
                df <- read.csv(input$file1$datapath)
                df$RMD <- round((df$mz - round(df$mz)) / df$mz * 10 ^ 6)
                df$MD <- round((df$mz - round(df$mz)) * 10 ^ 3)
                # high order mass defect computation

                getorder <- function(input) {
                        if (grepl(',', input)) {
                                name <- unlist(strsplit(input, ','))
                        } else{
                                name <- input
                        }
                        return(name)
                }

                getnum <- function(data) {
                        if (grepl('-', data)) {
                                name <- unlist(strsplit(data, '-'))
                                iso1 <-
                                        rcdk::get.isotopes.pattern(rcdk::get.formula(name[1]))
                                iso2 <-
                                        rcdk::get.isotopes.pattern(rcdk::get.formula(name[1]))
                                cus <-
                                        as.numeric(iso1[max(iso1[, 2]), 1]) - as.numeric(iso2[max(iso2[, 2]), 1])
                        } else{
                                iso <- rcdk::get.isotopes.pattern(rcdk::get.formula(data))
                                cus <-
                                        as.numeric(iso[max(iso[, 2]), 1])
                        }
                        return(cus)
                }

                temp <- getorder(input$cus1)
                cus <- NULL
                for (i in 1:length(temp)) {
                        cus <- c(cus, getnum(temp[i]))
                }

                if (length(cus) == 2) {
                        omd <- df$mz * round(cus[1]) / cus[1]
                        sumd <- cus[2] * round(cus[1]) / cus[1]
                        if (input$mdr1 == 'round') {
                                df$MD1_1 <-
                                        round(omd - round(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - round(sumd),
                                               digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1_1 <-
                                        round(omd - floor(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - floor(sumd),
                                               digits = 6)
                        } else{
                                df$MD1_1 <-
                                        round(omd - ceiling(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - ceiling(sumd),
                                               digits = 6)
                        }
                        smd <-  df$MD1_1 / md1_2
                        if (input$mdr1 == 'round') {
                                df$MD1_2 <-
                                        round(smd - round(smd),
                                              digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1_2 <-
                                        round(smd - floor(smd),
                                              digits = 6)
                        } else{
                                df$MD1_2 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                        }
                } else if (length(cus) == 3) {
                        omd <- df$mz * round(cus[1]) / cus[1]
                        sumd <- cus[2] * round(cus[1]) / cus[1]
                        tumd <- cus[3] * round(cus[1]) / cus[1]

                        if (input$mdr1 == 'round') {
                                df$MD1_1 <-
                                        round(omd - round(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - round(sumd),
                                               digits = 6)
                                md1_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1_1 <-
                                        round(omd - floor(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - floor(sumd),
                                               digits = 6)
                                md1_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        } else{
                                df$MD1_1 <-
                                        round(omd - ceiling(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - ceiling(sumd),
                                               digits = 6)
                                md1_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        }
                        smd <-  df$MD1_1 / md1_2
                        tsmd <- md1_3 / md1_2

                        if (input$mdr1 == 'round') {
                                df$MD1_2 <-
                                        round(smd - round(smd),
                                              digits = 6)
                                md1_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1_2 <-
                                        round(smd - floor(smd),
                                              digits = 6)
                                md1_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        } else{
                                df$MD1_2 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                                md1_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        }
                        tmd <- df$MD1_2 / md1_3
                        if (input$mdr1 == 'round') {
                                df$MD1_3 <-
                                        round(tmd - round(tmd),
                                              digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1_3 <-
                                        round(tmd - floor(tmd),
                                              digits = 6)
                        } else{
                                df$MD1_3 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                        }
                } else if (length(cus) > 3) {
                        message("Sorry, only the first three unit would be used.")
                        omd <- df$mz * round(cus[1]) / cus[1]
                        sumd <- cus[2] * round(cus[1]) / cus[1]
                        tumd <- cus[3] * round(cus[1]) / cus[1]

                        if (input$mdr1 == 'round') {
                                df$MD1_1 <-
                                        round(omd - round(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - round(sumd),
                                               digits = 6)
                                md1_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1_1 <-
                                        round(omd - floor(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - floor(sumd),
                                               digits = 6)
                                md1_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        } else{
                                df$MD1_1 <-
                                        round(omd - ceiling(omd),
                                              digits = 6)
                                md1_2 <- round(sumd - ceiling(sumd),
                                               digits = 6)
                                md1_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        }
                        smd <-  df$MD1_1 / md1_2
                        tsmd <- md1_3 / md1_2

                        if (input$mdr1 == 'round') {
                                df$MD1_2 <-
                                        round(smd - round(smd),
                                              digits = 6)
                                md1_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1_2 <-
                                        round(smd - floor(smd),
                                              digits = 6)
                                md1_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        } else{
                                df$MD1_2 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                                md1_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        }
                        tmd <- df$MD1_2 / md1_3
                        if (input$mdr1 == 'round') {
                                df$MD1_3 <-
                                        round(tmd - round(tmd),
                                              digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1_3 <-
                                        round(tmd - floor(tmd),
                                              digits = 6)
                        } else{
                                df$MD1_3 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                        }
                } else{
                        omd <- df$mz * round(cus) / cus
                        if (input$mdr1 == 'round') {
                                df$MD1 <-
                                        round(omd - round(omd),
                                              digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD1 <-
                                        round(omd - floor(omd),
                                              digits = 6)
                        } else{
                                df$MD1 <-
                                        round(omd - ceiling(omd),
                                              digits = 6)
                        }
                }

                temp <- getorder(input$cus2)
                cus <- NULL
                for (i in 1:length(temp)) {
                        cus <- c(cus, getnum(temp[i]))
                }
                if (length(cus) == 2) {
                        omd <- df$mz * round(cus[1]) / cus[1]
                        sumd <- cus[2] * round(cus[1]) / cus[1]
                        if (input$mdr2 == 'round') {
                                df$MD2_1 <-
                                        round(omd - round(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - round(sumd),
                                               digits = 6)
                        } else if (input$mdr2 == 'floor') {
                                df$MD2_1 <-
                                        round(omd - floor(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - floor(sumd),
                                               digits = 6)
                        } else{
                                df$MD2_1 <-
                                        round(omd - ceiling(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - ceiling(sumd),
                                               digits = 6)
                        }
                        smd <-  df$MD2_1 / md2_2
                        if (input$mdr2 == 'round') {
                                df$MD2_2 <-
                                        round(smd - round(smd),
                                              digits = 6)
                        } else if (input$mdr2 == 'floor') {
                                df$MD2_2 <-
                                        round(smd - floor(smd),
                                              digits = 6)
                        } else{
                                df$MD2_2 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                        }
                } else if (length(cus) == 3) {
                        omd <- df$mz * round(cus[1]) / cus[1]
                        sumd <- cus[2] * round(cus[1]) / cus[1]
                        tumd <- cus[3] * round(cus[1]) / cus[1]

                        if (input$mdr2 == 'round') {
                                df$MD2_1 <-
                                        round(omd - round(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - round(sumd),
                                               digits = 6)
                                md2_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        } else if (input$mdr2 == 'floor') {
                                df$MD2_1 <-
                                        round(omd - floor(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - floor(sumd),
                                               digits = 6)
                                md2_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        } else{
                                df$MD2_1 <-
                                        round(omd - ceiling(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - ceiling(sumd),
                                               digits = 6)
                                md2_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        }
                        smd <-  df$MD2_1 / md2_2
                        tsmd <- md2_3 / md2_2

                        if (input$mdr2 == 'round') {
                                df$MD2_2 <-
                                        round(smd - round(smd),
                                              digits = 6)
                                md2_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        } else if (input$mdr2 == 'floor') {
                                df$MD2_2 <-
                                        round(smd - floor(smd),
                                              digits = 6)
                                md2_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        } else{
                                df$MD2_2 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                                md2_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        }
                        tmd <- df$MD2_2 / md2_3
                        if (input$mdr2 == 'round') {
                                df$MD2_3 <-
                                        round(tmd - round(tmd),
                                              digits = 6)
                        } else if (input$mdr2 == 'floor') {
                                df$MD2_3 <-
                                        round(tmd - floor(tmd),
                                              digits = 6)
                        } else{
                                df$MD2_3 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                        }
                } else if (length(cus) > 3) {
                        message("Sorry, only the first three unit would be used.")
                        omd <- df$mz * round(cus[1]) / cus[1]
                        sumd <- cus[2] * round(cus[1]) / cus[1]
                        tumd <- cus[3] * round(cus[1]) / cus[1]

                        if (input$mdr2 == 'round') {
                                df$MD2_1 <-
                                        round(omd - round(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - round(sumd),
                                               digits = 6)
                                md2_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        } else if (input$mdr2 == 'floor') {
                                df$MD2_1 <-
                                        round(omd - floor(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - floor(sumd),
                                               digits = 6)
                                md2_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        } else{
                                df$MD2_1 <-
                                        round(omd - ceiling(omd),
                                              digits = 6)
                                md2_2 <- round(sumd - ceiling(sumd),
                                               digits = 6)
                                md2_3 <- round(tumd - round(tumd),
                                               digits = 6)
                        }
                        smd <-  df$MD2_1 / md2_2
                        tsmd <- md2_3 / md2_2

                        if (input$mdr2 == 'round') {
                                df$MD2_2 <-
                                        round(smd - round(smd),
                                              digits = 6)
                                md2_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD2_2 <-
                                        round(smd - floor(smd),
                                              digits = 6)
                                md2_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        } else{
                                df$MD2_2 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                                md2_3 <- round(tsmd - ceiling(tsmd),
                                               digits = 6)
                        }
                        tmd <- df$MD2_2 / md2_3
                        if (input$mdr2 == 'round') {
                                df$MD2_3 <-
                                        round(tmd - round(tmd),
                                              digits = 6)
                        } else if (input$mdr1 == 'floor') {
                                df$MD2_3 <-
                                        round(tmd - floor(tmd),
                                              digits = 6)
                        } else{
                                df$MD2_3 <-
                                        round(smd - ceiling(smd),
                                              digits = 6)
                        }
                } else{
                        omd <- df$mz * round(cus) / cus
                        if (input$mdr2 == 'round') {
                                df$MD2 <-
                                        round(omd - round(omd),
                                              digits = 6)
                        } else if (input$mdr2 == 'floor') {
                                df$MD2 <-
                                        round(omd - floor(omd),
                                              digits = 6)
                        } else{
                                df$MD2 <-
                                        round(omd - ceiling(omd),
                                              digits = 6)
                        }
                }
                return(df)
        })

        # Filtering the intensity, mz, and rt
        output$slide1 <- renderUI({
                minZ <- min(MD_data()$intensity)
                maxZ <- max(MD_data()$intensity)

                sliderInput(
                        "slide1",
                        "Intensity range filter",
                        min = minZ,
                        max = maxZ,
                        value = c(minZ, maxZ)
                )
        })
        output$slide2 <- renderUI({
                minZ <- min(MD_data()$mz)
                maxZ <- max(MD_data()$mz)

                sliderInput(
                        "slide2",
                        "mass to charge ratio range",
                        min = minZ,
                        max = maxZ,
                        value = c(minZ, maxZ)
                )
        })
        output$slide3 <- renderUI({
                minZ <- min(MD_data()$rt)
                maxZ <- max(MD_data()$rt)

                sliderInput(
                        "slide3",
                        "retention time range",
                        min = minZ,
                        max = maxZ,
                        value = c(minZ, maxZ)
                )
        })


        ## for plot control ##
        output$plot <- renderUI({
                if (input$single == "Single") {
                        plotlyOutput("DTPlot1")
                } else{
                        fluidRow(column(6, plotlyOutput("DTPlot1")),
                                 column(6, plotlyOutput("DTPlot2")))
                }
        })
        output$plotctr <- renderUI({
                if (input$single == "Single") {
                        fluidRow(
                                h4("Plot controls"),
                                tags$br(),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'xvar1',
                                                label = 'X variable for plot',
                                                choices = names(MD_data()),
                                                selected = names(MD_data())[1]
                                        )
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'yvar1',
                                                label = 'Y variable for plot',
                                                choices = names(MD_data()),
                                                selected = names(MD_data())[4]
                                        )
                                )
                        )

                } else{
                        fluidRow(
                                h4("Plot controls"),
                                tags$br(),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'xvar1',
                                                label = 'X variable for Plot 1',
                                                choices = names(MD_data()),
                                                selected = names(MD_data())[1]
                                        )
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'yvar1',
                                                label = 'Y variable for Plot 1',
                                                choices = names(MD_data()),
                                                selected = names(MD_data())[4]
                                        )
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'xvar2',
                                                label = 'X variable for Plot 2',
                                                choices = names(MD_data()),
                                                selected = names(MD_data())[1]
                                        )
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'yvar2',
                                                label = 'Y variable for Plot 2',
                                                choices = names(MD_data()),
                                                selected = names(MD_data())[4]
                                        )
                                )
                        )
                }
        })
        output$plotctr2 <- renderUI({
                if (input$single == "Single") {
                        fluidRow(
                                tags$br(),
                                textInput('x1', 'x axis label', input$xvar1),
                                textInput('y1', 'y axis label', input$yvar1)
                        )
                } else{
                        fluidRow(
                                tags$br(),
                                textInput(
                                        'x1',
                                        'x axis label for plot 1',
                                        input$xvar1
                                ),
                                textInput(
                                        'y1',
                                        'y axis label for plot 1',
                                        input$yvar1
                                ),
                                textInput(
                                        'x2',
                                        'x axis label for plot 2',
                                        input$xvar2
                                ),
                                textInput(
                                        'y2',
                                        'y axis label for plot 2',
                                        input$yvar2
                                )
                        )
                }

        })
        #### For MD Plot Panel ####

        #OE#
        observeEvent(input$go, {
                m <- MD_data()
                m <-
                        m[m$intensity >= input$slide1[1] &
                                  m$intensity <= input$slide1[2] &
                                  m$mz >= input$slide2[1] &
                                  m$mz <= input$slide2[2] &
                                  m$rt >= input$slide3[1] &
                                  m$rt <= input$slide3[2],]
                d <- SharedData$new(m)

                MDplot_y1 <-
                        m[, input$yvar1]

                MDplot_x1 <-
                        m[, input$xvar1]

                # Checkbox option for size of markers by intensity
                if (input$ins) {
                        intensity <- m$intensity
                } else{
                        intensity <- NULL
                }

                if (input$single == "Double") {
                        MDplot_x2 <-
                                m[, input$xvar2]


                        MDplot_y2 <-
                                m[, input$yvar2]

                }

                # highlight selected rows in the scatterplot
                output$DTPlot1 <- renderPlotly({
                        s <- input$x1_rows_selected
                        if (!length(s)) {
                                p <- d %>%
                                        plot_ly(
                                                x = MDplot_x1,
                                                y = MDplot_y1,
                                                type = "scatter",
                                                size = intensity,
                                                mode = "markers",
                                                marker = list(
                                                        line = list(
                                                                width = 1,
                                                                color = '#FFFFFF'
                                                        )
                                                ),
                                                color = I('black'),
                                                name = 'Unfiltered'
                                        ) %>%
                                        layout(
                                                legend = list(
                                                        orientation = "h",
                                                        xanchor = "center",
                                                        x = 0.5,
                                                        y = 100
                                                ),
                                                showlegend = T,
                                                xaxis = list(title = input$x1),
                                                yaxis = list(title = input$y1)
                                        ) %>%
                                        highlight(
                                                "plotly_selected",
                                                color = I('red'),
                                                selected = attrs_selected(name = 'Filtered')
                                        )
                        } else if (length(s)) {
                                pp <- m %>%
                                        plot_ly() %>%
                                        add_trace(
                                                x = MDplot_x1,
                                                y = MDplot_y1,
                                                type = "scatter",
                                                size = intensity,
                                                mode = "markers",
                                                marker = list(
                                                        line = list(
                                                                width = 1,
                                                                color = '#FFFFFF'
                                                        )
                                                ),
                                                color = I('black'),
                                                name = 'Unfiltered'
                                        ) %>%
                                        layout(
                                                legend = list(
                                                        orientation = "h",
                                                        xanchor = "center",
                                                        x = 0.5,
                                                        y = 100
                                                ),
                                                showlegend = T,
                                                xaxis = list(title = input$x1),
                                                yaxis = list(title = input$y1)
                                        )

                                # selected data
                                pp <-
                                        add_trace(
                                                pp,
                                                data = m[s, , drop = F],
                                                x = MDplot_x1[s],
                                                y = MDplot_y1[s],
                                                type = "scatter",
                                                size = intensity[s],
                                                mode = "markers",
                                                marker = list(
                                                        line = list(
                                                                width = 1,
                                                                color = '#FFFFFF'
                                                        )
                                                ),
                                                color = I('red'),
                                                name = 'Filtered'
                                        )
                        }

                })

                # Plot 2
                if (input$single == "Double") {
                        output$DTPlot2 <- renderPlotly({
                                t <- input$x1_rows_selected

                                if (!length(t)) {
                                        p <- d %>%
                                                plot_ly(
                                                        x = MDplot_x2,
                                                        y = MDplot_y2,
                                                        type = "scatter",
                                                        size = intensity,
                                                        mode = "markers",
                                                        marker = list(
                                                                line = list(
                                                                        width = 1,
                                                                        color = '#FFFFFF'
                                                                )
                                                        ),
                                                        color = I('black'),
                                                        name = 'Unfiltered'
                                                ) %>%
                                                layout(
                                                        legend = list(
                                                                orientation = "h",
                                                                xanchor = "center",
                                                                x = 0.5,
                                                                y = 100
                                                        ),
                                                        showlegend = T,
                                                        xaxis = list(title = input$x2),
                                                        yaxis = list(title = input$y2)
                                                ) %>%
                                                highlight(
                                                        "plotly_selected",
                                                        color = I('red'),
                                                        selected = attrs_selected(name = 'Filtered')
                                                )
                                } else if (length(t)) {
                                        pp <- m %>%
                                                plot_ly() %>%
                                                add_trace(
                                                        x = MDplot_x2,
                                                        y = MDplot_y2,
                                                        type = "scatter",
                                                        size = intensity,
                                                        mode = "markers",
                                                        marker = list(
                                                                line = list(
                                                                        width = 1,
                                                                        color = '#FFFFFF'
                                                                )
                                                        ),
                                                        color = I('black'),
                                                        name = 'Unfiltered'
                                                ) %>%
                                                layout(
                                                        legend = list(
                                                                orientation = "h",
                                                                xanchor = "center",
                                                                x = 0.5,
                                                                y = 100
                                                        ),
                                                        showlegend = T,
                                                        xaxis = list(title = input$x2),
                                                        yaxis = list(title = input$y2)
                                                )

                                        # selected data
                                        pp <-
                                                add_trace(
                                                        pp,
                                                        data = m[t, , drop = F],
                                                        x = MDplot_x2[t],
                                                        y = MDplot_y2[t],
                                                        type = "scatter",
                                                        size = intensity[t],
                                                        mode = "markers",
                                                        marker = list(
                                                                line = list(
                                                                        width = 1,
                                                                        color = '#FFFFFF'
                                                                )
                                                        ),
                                                        color = I('red'),
                                                        name = 'Filtered'
                                                )
                                }

                        })
                }
                # highlight selected rows in the table
                output$x1 <- renderDT({
                        T_out1 <- m[d$selection(), ]
                        dt <-
                                DT::datatable(
                                        m,
                                        editable = TRUE,
                                        rownames = FALSE,
                                        filter = "top"
                                )
                        if (NROW(T_out1) == 0) {
                                dt
                        } else {
                                T_out1
                        }
                })

                # download the filtered data
                output$x3 = downloadHandler(
                        'MDplot-filtered.csv',
                        content = function(file) {
                                s <- input$x1_rows_selected
                                if (length(s)) {
                                        write.csv(m[s, , drop = FALSE], file)
                                } else if (!length(s)) {
                                        write.csv(m[d$selection(), ], file)
                                }
                        }
                )
        })
}
)
