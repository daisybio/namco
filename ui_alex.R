library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(plotly)
library(visNetwork)


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
                column(7,dataTableOutput("betaTable")),
                column(3,selectInput("betaGroup","Group by:",c("bla","bli")))
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
                #,column(6,plotlyOutput("countDistr"))
                ),
              p("Network Plot"),
              fluidRow(
                column(1),
                column(10,visNetworkOutput("basicNetwork"))
            )),
            tabPanel("Advanced Approach",
              p("themetagenomics approach: Taxa Features"),       
              tags$hr(),
              fluidRow(
                column(3,sliderInput("K","Pick Number of Topics:", 10, 0, 100, step=1)),
                column(3,sliderInput("sigma_prior","Pick sigma_prior value:", 0, 0, 1, step=0.01)),
                column(3,actionButton("themeta","Start themetagenomics Calculation!"))
              ),
              tags$hr()
              ,
              fixedRow(
                column(1,radioButtons('dim',label=strong('Dim'),
                                      choices=list('2D'='2d','3D'='3d'),
                                      selected='2d')),
                column(3,selectInput('dist',label=strong('Method'),
                                     choices=list('Bray Curtis'='bray','Jaccard'='jaccard','Euclidean'='euclidean',
                                                  'Hellinger'='hellinger','Chi Squared'='chi2','Jensen Shannon'='jsd',
                                                  't-SNE'='tsne'),
                                     selected='jsd')),
                column(1,style='padding: 25px 0px;',actionButton('reset','Reset')),
                column(2,numericInput('k_in',label=strong('Topic Number'),value=0,min=0,max=K,step=1)),
                column(3,sliderInput('lambda',label=strong('Lambda'),min=0,max=1,value=1,step=lambda_step)),
                column(2,selectInput('taxon',label=strong('Taxon'),
                                     choices=list('Phylum'='Phylum','Class'='Class','Order'='Order',
                                                  'Family'='Family','Genus'='Genus')))
              ),
              
              fixedRow(
                column(1,tags$div('Number of components to plot.',class='capt')),
                column(3,tags$div('Type of distance and method for ordination.',class='capt')),
                column(1,tags$div('Reset topic selection.',class='capt')),
                column(2,tags$div('Current selected topic.',class='capt')),
                column(3,tags$div(paste0('Relative weighting that influences taxa shown in barplot.',
                                         ' If equal to 1, p(taxa|topic)l if 0, p(taxa|topic)/p(taxa).'),class='capt')),
                column(2,tags$div('Taxonomic group to dictate bar plot shading',class='capt'))
              ),
              fixedRow(
                column(6,offset=0,height='600px',plotlyOutput('ord')),
                column(6,offset=0,height='600px',plotOutput('bar'))
              ),
              
              fixedRow(
                column(6,tags$div(paste0('Ordination of the samples over topics distribution theta, colored according to',
                                         ' the weights shown in the scatter plot above. The radius of a given point',
                                         ' represents the marginal topic frequency. The amount of variation explained is',
                                         ' annotated on each axis.'),class='below')),
                column(6,tags$div(paste0('Bar plot representing the taxa frequencies. When no topic is selected, the overall',
                                         ' taxa frequencies are shown, colored based on the selected taxonomy and ordered in',
                                         ' in terms of saliency. When a topic is chosen, the red bars show the margina taxa',
                                         ' frequency within the selected topic, ordered in terms of relevency, which in turn',
                                         ' can be reweighted by adjusting the lambda slider.'),class='below'))
              ),
              
              br(),
              
              fixedRow(
                column(1,''),
                column(10,htmlOutput('text2')),
                column(1,'')
              ),
              networkD3::forceNetworkOutput('corr')
            )
          )
        )
      )
    )
  )
)