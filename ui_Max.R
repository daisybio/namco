library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(plotly)
library(networkD3)

ui <- dashboardPage(
  dashboardHeader(title="Microbiome Explorer"),
  dashboardSidebar(
    sidebarMenu(id="sidebar",
      br(),
      fluidRow(
        column(12,align="center",actionButton("upload","Load new dataset"))
      ),
      dataTableOutput("datasets"),
      br(),br(),
      menuItem("Quality Control",tabName="QC",icon=icon("search")),
      menuItem("Network Analysis",tabName="Network",icon=icon("project-diagram"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "QC",
        h4("Quality Control Report"),
        fluidRow(  
          tabBox(id="qualityPlots",width=12,
            tabPanel("Sample Information",
              tags$hr(),
              p("Sample meta data"),
              fluidRow(
                column(1),
                column(7,dataTableOutput("metaTable"))
            )),
            tabPanel("Rarefaction Curves",
              tags$hr(),
              fluidRow(
                column(1),
                column(6,plotlyOutput("rarefacCurve")),
                column(1),
                column(2,sliderInput("rareToShow","Number of samples to display:",10,100,50),
                       br(),
                       sliderInput("rareToHighlight","Quantile to highlight:",0,100,2))
            )),
            tabPanel("Taxa Distribution",
              tags$hr(),
              fluidRow(
              column(1),
                column(6,plotlyOutput("taxaDistribution")),
                column(1),
                column(2,sliderInput("otherCutoff","Lower percentage to bin:",1,20,1))
            )),
            tabPanel("PCA", 
              tags$hr(),
              fluidRow(
                column(1),
                column(6,plotlyOutput("pcaPlot")),
                column(1),
                column(2,br(),
                  radioButtons("pcaMode", "PCA Mode:",c("2D","3D")))
                  #  br(),
                  #  selectInput("pcaGroups","Group by:",""))
                ),br(),br(),br(),
              fluidRow(
                column(1),
                column(6,plotlyOutput("loadingsPlot")),
                column(1),
                column(2,selectInput("pcaLoading","Plot loadings on PC",1:3))
            )),
            tabPanel("Alpha Diversity",
              tags$hr(),
              fluidRow(
                column(1),
                column(7,plotlyOutput("alphaPlot")),
                column(3,
                  selectInput("alphaMethod","Method:",c("Shannon Entropy","effective Shannon Entropy","Simpson Index","effective Simpson Index")),
                  selectInput("alphaGroup","Group by:",c("bla","bli"))
                )
            )),
            tabPanel("Beta Diversity",
              tags$hr(),
              fluidRow(
              column(1),
                column(7,
                  plotOutput("betaTree"),
                  plotOutput("betaMDS"),
                  plotOutput("betaNMDS")
                ),
                column(3,
                  br(),
                  selectInput("betaMethod","Method:",c("Bray-Curtis Dissimilarity","Generalized UniFrac Distance")),
                  selectInput("betaGroup","Group by:",choices=""))
            ))
          )
        )
      ),
      tabItem(tabName = "Network",
        h4("Network Analysis"),
        fluidRow(  
          tabBox(id="netWorkPlots",width=12,
            tabPanel("Basic Approach",
              p("Co-occurrences are counted if OTU is present in both samples"),
              tags$hr(),
              fluidRow(
                column(4,sliderInput("binCutoff","Cutoff for Binarization",0,10,1,step = 0.01)),
                column(4,plotlyOutput("cutoffHist")),
                column(2,radioButtons("useFC","Calculation of Counts:",c("log2(fold-change)","difference"))),
                column(2,selectInput("groupCol","Select Column from META-file containing groupings:",c()))),
              fluidRow(
                column(2,actionButton("startCalc","Start Count Calculation & Reload Network!"))
              ),
              p("Network Plot"),
              fluidRow(
                column(1),
                column(9,simpleNetworkOutput("basicNetwork")),
                column(2,sliderInput("networkCutoff","Number of interactions to show:",5,500,30,5))
              )),
            tabPanel("Advanced Approach",
              tags$hr(),
              fluidRow(
                column(12)
            ))
          )
        )
      )
    )
  )
)