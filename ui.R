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
      menuItem("Welcome!",tabName="welcome",icon=icon("door-open")),
      menuItem("Basic Analyses",tabName="basics",icon=icon("search")),
      menuItem("Network Analysis",tabName="Network",icon=icon("project-diagram"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName="welcome",
        fluidRow(
          column(1),
          column(8,htmlOutput("welcome"))
        ),
        tags$hr(),
        br(),
        fluidRow(
          column(1),
          column(10,htmlOutput("authors"))
        ),
        tags$hr(),
        fluidRow(
          column(1),
          column(10,htmlOutput("welcome_ref"))
      )),
      tabItem(tabName = "basics",
        h4("Basic Analysis"),
        fluidRow(  
          tabBox(id="basicPlots",width=12,
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
                column(2,br(),
                  sliderInput("rareToShow","Number of samples to display:",min=1,max=1,step=1,value=1),
                  br(),
                  sliderInput("rareToHighlight","Quantile of most undersampled samples to highlight:",0,100,2))
              ),
              br(),br(),
              fluidRow(
                column(2),
                column(4,verbatimTextOutput("undersampled")),
                column(2,checkboxInput("excludeSamples","exclude undersampled samples"))
            )),
            tabPanel("Taxa Distribution",
              tags$hr(),
              fluidRow(
                column(1),
                column(10,plotlyOutput("taxaDistribution",height="auto"))
              ),br(),br(),
              fluidRow(
                column(4),
                column(4,
                  selectInput("taxLevel","Taxonomic Level",choices=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")),
                  sliderInput("otherCutoff","Bin organims below an average abundance of:",0,20,0))
              )
            ),
            tabPanel("Data Structure", 
              tags$hr(),
              fluidRow(
                column(1),
                column(6,plotlyOutput("structurePlot")),
                column(1),
                column(2,br(),
                  selectInput("structureMethod","",c("PCA","t-SNE","UMAP")),
                  br(),br(),
                  radioButtons("structureDim","Dimensions:",c("2D","3D")),
                  br(),
                  selectInput("structureGroup","Group by:",""))
                ),br(),br(),br(),
              conditionalPanel("input.structureMethod == 'PCA'",
                fluidRow(
                  column(1),
                  column(6,plotlyOutput("loadingsPlot")),
                  column(1),
                  column(2,selectInput("pcaLoading","Plot loadings on PC",1:3))
                )
              )
            ),
            tabPanel("Alpha Diversity",
              tags$hr(),
              fluidRow(
                column(1),
                column(7,plotlyOutput("alphaPlot")),
                column(3,
                  selectInput("alphaMethod","Method:",c("Shannon Entropy","effective Shannon Entropy","Simpson Index","effective Simpson Index","Richness")),
                  selectInput("alphaGroup","Group by:","")
                )
              ),
              br(),br(),
              fluidRow(
                column(1),
                column(7,dataTableOutput("alphaTable"))
              )
            ),
            tabPanel("Beta Diversity",
              tags$hr(),
              fluidRow(
                column(1),
                column(7,
                  plotOutput("betaTree")
                ),
                column(3,
                  br(),
                  selectInput("betaMethod","Method:",choices=""),
                  #conditionalPanel("length() == T",
                  #  infoBox("Remark","Bla")
                  #),
                  selectInput("betaGroup","Group by:",choices="")
                )
              ),
              fluidRow(
                column(1),
                column(5,plotOutput("betaMDS")),
                column(5,plotOutput("betaNMDS"))
              ),
              br(),br(),br(),
              fluidRow(
                column(1),
                column(7,dataTableOutput("betaTable"))
              )
            )
          )
        )
      ),
      tabItem(tabName = "Network",
        h4("Network Analysis"),
        fluidRow(  
          tabBox(id="netWorkPlots",width=12,
            tabPanel("Co-occurrence of OTUs",
              p("Co-occurrences are counted if OTU is present in both samples"),
              tags$hr(),
              htmlOutput("cutoff_title"),
              fluidRow(
                column(6,
                       fixedRow(
                         column(3,numericInput("binCutoff","Cutoff for Binarization",min=0.01,max=10,value=1,step = 0.01)),
                         column(3,htmlOutput("log_cutoff"))
                       ),
                  plotlyOutput("cutoffHist")),
                column(4,offset=1,
                  plotlyOutput("boolHeat"))),
              fluidRow(column(6,htmlOutput("cutoff_text")),
                column(1),
                column(4,htmlOutput("heatmap_text"))),
              tags$hr(),
              fluidRow(
                column(4,
                  htmlOutput("basic_calc_title"),
                  radioButtons("useFC","Calculation of Counts:",c("log2(fold-change)","difference"))
                ),
                column(4,selectInput("groupCol","Select Column from META-file containing groupings:",choices = c("Please Upload OTU & META file first!"),selected = "Please Upload OTU & META file first!")),
                column(3,actionButton("startCalc","Start Count Calculation & Reload Network!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))),
              tags$hr(),
              fluidRow(
                column(1),
                column(8,simpleNetworkOutput("basicNetwork")),
                column(3,sliderInput("networkCutoff","Number of edges to show:",50,500,50,10))
              ),
              fixedRow(
                column(1),
                column(10,
                  htmlOutput("basic_additional"),
                  htmlOutput("basic_calc_additional")),
                column(1)
              )
            ),
            tabPanel("Clustering based on functional Topics",
              fluidRow(column(12,htmlOutput("advanced_text"))),
              tags$hr(),
              fluidRow(
                column(3,
                  sliderInput("K","Pick Number of Topics:", 1, 150, 30, step=1),
                  htmlOutput("topic_text")),
                column(3,
                  sliderInput("sigma_prior","Pick Scalar between 0 and 1:", 0, 1, 0, step=0.01),
                  htmlOutput("sigma_text")),
                column(3,
                  selectInput("formula", label = "Formula for covariates of interest found in metadata:",choices="Please provide OTU-table & metadata first!"),
                  selectInput("refs", label = "Binary covariate in formula, indicating the reference level:",choices = "Please provide OTU-table & metadata first!")),
                column(3,
                  actionButton("themeta","Start functional Clustering!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                  #,downloadButton("downloadGeneTable","Download Gene-Table")
                  )
              ),
              tags$hr(),
              fluidRow(
                column(3,h4("Input Variables:")),
                column(10,htmlOutput("input_variables")),
                column(1)
              ),
              tags$hr(),
              fixedRow(
                column(1),
                column(10,htmlOutput('text1')),
                column(1)
              ),
              br(),
              fixedRow(
                column(4,selectInput('choose', label='Covariate',
                  choices="Please start calculation above first!"),
              fixedRow(column(1),
                column(11,tags$div(paste0('Choosing a covariate determines which weight estimates will shown',
                  ' The order of the topics will be adjusted accordingly. By clicking',
                  ' an estimate, all figures below will rerender.'),class='side')))),
                column(10,plotlyOutput('est',height='200px'))
              ),
              br(),
              fixedRow(
                column(1,radioButtons('dim',label=strong('Dim'),choices=list('2D'='2d','3D'='3d'),selected='2d')),
                column(3,selectInput('dist',label=strong('Method'),choices=list('Bray Curtis'='bray',
                  'Jaccard'='jaccard','Euclidean'='euclidean','Hellinger'='hellinger','Chi Squared'='chi2',
                  'Jensen Shannon'='jsd','t-SNE'='tsne'),selected='jsd')),
                column(1,style='padding: 25px 0px;',actionButton('reset','Reset')),
                column(2,numericInput('k_in',label=strong('Topic Number'),value=0,min=0,max=100,step=1)),
                column(3,sliderInput('lambda',label=strong('Lambda'),min=0,max=1,value=1,step=0.01)),
                column(2,selectInput('taxon',label=strong('Taxon'),choices=list('Phylum'='Phylum',
                  'Class'='Class','Order'='Order','Family'='Family','Genus'='Genus')))
              ),
              fixedRow(
                column(1,tags$div('Number of components to plot.',class='capt')),
                column(3,tags$div('Type of distance and method for ordination.',class='capt')),
                column(1,tags$div('Reset topic selection.',class='capt')),
                column(2,tags$div('Current selected topic.',class='capt')),
                column(3,tags$div(paste0('Relative weighting of selected topic that influences taxa shown in barplot.',
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
                column(1),
                column(10,htmlOutput('text2')),
                column(1)
              ),
              forceNetworkOutput('corr'),
              br(),
              fixedRow(
                column(1,''),
                column(10,htmlOutput('text3')),
                column(1)
              )
            )
          )
        )
      )
    )
  )
)