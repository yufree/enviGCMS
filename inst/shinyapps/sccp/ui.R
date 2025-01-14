library(shiny)

shinyUI(navbarPage(
        "MSCCP",
        tabPanel("Usage",
                 fluidPage(
                         titlePanel(
                                 "Online shiny apps for quantitative analysis of Short Chain Chlorinated Paraffins (SCCPs)"
                         ),
                         h3('Step 0'),
                         h4('Prepare standards with known %Cl.'),
                         h3('Step 1'),
                         h4('Covert you GC-HRMS data into mzxml format by',a('MsConverter', href ='http://proteowizard.sourceforge.net/tools.shtml')),
                         h3('Step 2'),
                         h4('Upload standards data in "Standards analysis" tab, set the parameters for internal standards and retention time range, clike "Go" to process the data. Record the intercept and slope of the model.'),
                         h3('Step 3'),
                         h4('Upload Samples data in "Sample analysis" tab, set the parameters for intercept and slope from the standards, internal standards and retention time range, clike "Go" to process the data.'),
                         h3('Step 4'),
                         h4('The concertration and compositon of your samples would be shown.'),
                         h3('Tips'),
                         h4('"Ions for analysis" show the ions used for this analysis in case you use other HRMS to analysis SCCPs.')
                 )),
        tabPanel(
                "Standards analysis",
                fluidPage(
                        titlePanel("Analysis for SCCPs standards"),
                        sidebarPanel(
                                h4('Uploading Standards File'),
                                fileInput('Standards',
                                          label = 'mzxml files',
                                          multiple = T,
                                          accept = c('.mzxml')),
                        sliderInput(
                                'ISmz',
                                'm/z',
                                min = 200,
                                max = 600,
                                value = 323,
                                step = 0.01,
                                round = 3
                        ),
                        sliderInput(
                                'ISrt',
                                'Retention Time range for IS',
                                min = 0,
                                max = 6000,
                                value = c(600,1200),
                                step = 50,
                                round = 0
                        ),
                        sliderInput(
                                'SCCPrt',
                                'Retention Time range for SCCPs',
                                min = 0,
                                max = 6000,
                                value = c(600,1200),
                                step = 50,
                                round = 0
                        ),
                        sliderInput(
                                'ppm',
                                'ppm',
                                min = 1,
                                max = 100,
                                value = 5,
                                step = 1,
                                round = 0
                        ),
                        sliderInput(
                                'con',
                                'Concertration in ppm',
                                min = 10,
                                max = 10000,
                                value = 2000,
                                step = 10,
                                round = 0
                        ),
                        actionButton("go", "Go")
                        ),

                        mainPanel(
                               plotOutput("plotstd") ,
                               plotOutput("plotcomp"),

                               h4('Parameters for linear regression'),
                               tableOutput("reg2"),

                               h4("Parameters for log transformation on response factor"),
                               tableOutput("reg")

                        )
                )
                ),
        tabPanel(
                "Sample analysis",
                fluidPage(
                        titlePanel("Analysis for SCCPs sample"),
                        sidebarPanel(
                                h4('Uploading Sample Files'),
                                fileInput('Samples',
                                          label = 'mzxml files',
                                          multiple = T,
                                          accept = c('.mzxml')),
                                checkboxInput('log','Log Trans', value = T),
                                sliderInput(
                                        'inc',
                                        'Intercept',
                                        min = -50,
                                        max = 0,
                                        value = -26.36,
                                        step = 0.01,
                                        round = 2
                                ),
                                sliderInput(
                                        'slope',
                                        'Slope',
                                        min = 0,
                                        max = 50,
                                        value = 34.66,
                                        step = 0.01,
                                        round = 2
                                ),
                                sliderInput(
                                'ISmz',
                                'm/z',
                                min = 200,
                                max = 600,
                                value = 323,
                                step = 0.01,
                                round = 3
                        ),
                        sliderInput(
                                'ISrt',
                                'Retention Time range for IS',
                                min = 0,
                                max = 6000,
                                value = c(600,1200),
                                step = 50,
                                round = 0
                        ),
                        sliderInput(
                                'SCCPrt',
                                'Retention Time range for SCCPs',
                                min = 0,
                                max = 6000,
                                value = c(600,1200),
                                step = 50,
                                round = 0
                        ),
                        sliderInput(
                                'ppm',
                                'ppm',
                                min = 1,
                                max = 100,
                                value = 5,
                                step = 1,
                                round = 0
                        ),
                        actionButton("go2", "Go")),

                mainPanel(plotOutput("plotcomps"),
                          h4("SCCPs Concertrations in sample(s)(ppm):"),
                          textOutput("results"))

        )),
        tabPanel("Ions for analysis",
                 dataTableOutput("data")),
        tabPanel("References",
                 p("Bogdal, C., Alsberg, T., Diefenbacher, P.S., MacLeod, M., Berger, U., 2015. Fast Quantification of Chlorinated Paraffins in Environmental Samples by Direct Injection High-Resolution Mass Spectrometry with Pattern Deconvolution. Anal. Chem. 87, 2852–2860. doi:10.1021/ac504444d"),
                 p("Gao, W., Wu, J., Wang, Y., Jiang, G., 2016. Quantification of short- and medium-chain chlorinated paraffins in environmental samples by gas chromatography quadrupole time-of-flight mass spectrometry. Journal of Chromatography A 1452, 98–106. doi:10.1016/j.chroma.2016.04.081"),
                 p("ISO 12010:2012 - Water quality -- Determination of short-chain polychlorinated alkanes (SCCPs) in water -- Method using gas chromatography-mass spectrometry (GC-MS) and negative-ion chemical ionization (NCI) [WWW Document], n.d. URL https://www.iso.org/standard/51124.html (accessed 10.25.17)."),
                 p("Nilsson, M.-L., Bengtsson, S., Kylin, H., 2012. Identification and determination of chlorinated paraffins using multivariate evaluation of gas chromatographic data. Environmental Pollution 163, 142–148. doi:10.1016/j.envpol.2011.12.010"),
                 p("Rusina, T., Korytar, P., Boer, J. de, 2011. Comparison of quantification methods for the analysis of polychlorinated alkanes using electron capture negative ionization mass spectrometry. Intern J Environ Anal Chem, INT J ENVIRON AN CH, Intern. J. Environ. Anal. Chem. 91, 319–332. doi:10.1080/03067311003602583"),
                 p("Tomy, G.T., Stern, G.A., Muir, D.C.G., Fisk, A.T., Cymbalisty, C.D., Westmore, J.B., 1997. Quantifying C10−C13 Polychloroalkanes in Environmental Samples by High-Resolution Gas Chromatography/Electron Capture Negative Ion High-Resolution Mass Spectrometry. Anal. Chem. 69, 2762–2771. doi:10.1021/ac961244y"),
                 p("van Mourik, L.M., Leonards, P.E.G., Gaus, C., de Boer, J., 2015. Recent developments in capabilities for analysing chlorinated paraffins in environmental matrices: A review. Chemosphere 136, 259–272. doi:10.1016/j.chemosphere.2015.05.045"),
                 p("Yuan, B., Bogdal, C., Berger, U., MacLeod, M., Gebbink, W.A., Alsberg, T., de Wit, C.A., 2017. Quantifying Short-Chain Chlorinated Paraffin Congener Groups. Environ. Sci. Technol. 51, 10633–10641. doi:10.1021/acs.est.7b02269")
)
))
