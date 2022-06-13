navbarPage("Systems Biology Ireland!",
           tabPanel("Cancer Statistics",
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons( "Stats","Statistics",
                                     c("Europe Cancer Incidence" = "europe", "Ireland Case Incidence 2020" ="ire")
                        )
                      ),
                      mainPanel(
                        plotOutput("statsPlot")
                      )
                    )
                    
           ),
           tabPanel("Cancer Mutations",
                             sidebarLayout(
                               sidebarPanel(
                                 radioButtons( "cancer","Mutations",
                                               c("Healthy Cells" = "hc", "Healthy Cells Division" = "hdiv",
                                                 "Cancer Cells with Mutation" ="cancer", 
                                                 "Cell Division Cancer Cells" = "cdiv")
                                 )
                               ),
                               mainPanel(
                                 plotOutput("cellshow")
                               )
                             )
                             
                    ),
           tabPanel("Mutation Frequencies",
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons( "cancer_type","Mutations in Common Cancer Types",
                                      c("Bowel Cancer" = "brc", "Bowel Cancer" = "bc",
                                        "Lung Cancer" ="lc", 
                                        "Prostate Cancer" = "pc", "Melanoma (Type of Skin Cancer)" ="sc" )
                        )
                      ),
                      mainPanel(
                        plotOutput("freqshow")
                      )
                    )
                    ),
           tabPanel(" Networks",
                    fluidRow(
                      column(width = 4, wellPanel(
                        radioButtons( "nets","Network Structure",
                                      c("Network 1" = "net1", "Network 2" = "net2", "Breast Cancer Network" = "net3",
                                        "Bowel Cancer Network" = "net4",
                                        "Lung Network" = "net5", "Melanoma Network" = "net6",
                                        "Prostate Network" = "net7")
                                  )
                             )
                      ),
                      column(width = 4, wellPanel(
                        plotOutput("netshow", 
                                   # Equivalent to: click = clickOpts(id = "plot_click")
                                   click = "plot1_click",dblclick = "plot1_dblclick",)
                                   )
                        
                        ),
                      column(width = 4, wellPanel(plotOutput("celldec"))
                    )
                  )
           ),
           tabPanel("Tumor Growth Simulation",
                    sidebarLayout(
                      # Sidebar with a slider input
                      sidebarPanel(
                        sliderInput("obs",
                                    "Days:",
                                    min = 0,
                                    max = 60,
                                    value = 0.1),
                        # Sidebar with a slider input
                        sliderInput("d1",
                                    "Drug 1:",
                                    min = 0,
                                    max = 1,
                                    value = 0.01),
                        # Sidebar with a slider input
                        sliderInput("d2",
                                    "Drug 2:",
                                    min = 0,
                                    max = 1,
                                    value = 0.01),
                        sliderInput("d3",
                                    "Drug 3:",
                                    min = 0,
                                    max = 1,
                                    value = 0.01)
                      ),
                      mainPanel(
                        plotOutput("distPlot")
                      )
                     )
           )
 
)