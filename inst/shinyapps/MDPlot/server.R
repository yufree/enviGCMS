library(shiny)
library(DT)
library(plotly)
library(crosstalk)
library(enviGCMS)

shinyServer(function(input, output, session) {
        MD_data <- reactive({
                #  require that the input is available
                req(input$file1)
                df <- read.csv(input$file1$datapath)
                df$RMD <- round((round(df$mz) - df$mz) / df$mz * 10 ^ 6)
                df$OMD <- round((round(df$mz) - df$mz) * 10 ^ 3)
                # high order mass defect computation

                mdh1 <- getmdh(df$mz,cus = input$cus1)
                mdh2 <- getmdh(df$mz,cus = input$cus2)
                # change colname
                name1 <- paste0(colnames(mdh1),'_p1')
                name2 <- paste0(colnames(mdh2),'_p2')
                mdh <- cbind(mdh1[,-1],mdh2[,-1])
                colnames(mdh) <- c(name1[-1],name2[-1])

                df <- cbind(df,mdh)

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
                                                choices = names(MD_data())
                                        )
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'yvar1',
                                                label = 'Y variable for plot',
                                                choices = names(MD_data())
                                        )
                                ),
                                column(
                                        12,
                                        selectInput(
                                                inputId = 'zvar1',
                                                label = 'Symbol variable for plot',
                                                choices = list(`NULL` = 'NA',`Variable` = names(MD_data()))
                                                ,
                                                selected = 'NULL'
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
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'zvar1',
                                                label = 'Symbol variable for plot',
                                                choices = list(`NULL` = 'NA',`Variable` = names(MD_data())),
                                                selected = 'NULL'
                                        )
                                ),
                                column(
                                        6,
                                        selectInput(
                                                inputId = 'zvar2',
                                                label = 'Symbol variable for plot 2',
                                                choices = list(`NULL` = 'NA',`Variable` = names(MD_data())),
                                                selected = 'NULL'
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

                if(input$zvar1 == 'NA'){
                        MDplot_z1 <- 1
                }else{
                        MDplot_z1 <- m[, input$zvar1]
                }

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

                        if(input$zvar2 == 'NA'){
                                MDplot_z2 <- 1
                        }else{
                                MDplot_z2 <- m[, input$zvar2]
                        }

                }

                # highlight selected rows in the scatterplot
                output$DTPlot1 <- renderPlotly({
                        s <- input$x1_rows_selected
                        if (!length(s)) {
                                p <- d %>%
                                        plot_ly(
                                                x = MDplot_x1,
                                                y = MDplot_y1,
                                                symbol = MDplot_z1,
                                                showlegend = input$show_leg,
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
                                                symbol = MDplot_z1,
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
                                                        symbol = MDplot_z2,
                                                        showlegend = input$show_leg,
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
                                                        symbol = MDplot_z2,
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
