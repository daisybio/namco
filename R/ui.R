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
      menuItem("Network Analysis",tabName="Network",icon=icon("project-diagram")),
      menuItem("Info & Settings",tabName = "info",icon=icon("cogs"))
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
                column(10,dataTableOutput("metaTable"))
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
                  sliderInput("rareToHighlight","Quantile of most undersampled samples to highlight:",0,100,0))
              ),
              br(),br(),
              fluidRow(
                column(2),
                column(4,verbatimTextOutput("undersampled")),
                column(2,switchInput("excludeSamples","exclude undersampled samples",value=F))
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
                  sliderInput("otherCutoff","Bin taxa below an average abundance of:",0,20,0))
            )),
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
            )),
            tabPanel("OTU Correlation", 
              tags$hr(),
              fluidRow(
                column(1),
                column(10,plotlyOutput("OTUcorrPlot"))
              ),
              br(),br(),
              fluidRow(
                column(4),
                column(4,sliderInput("corrCluster","Number of Clusters:",1,20,2))
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
                  selectInput("betaGroup","Group by:",choices="")
                )
              ),
              fluidRow(
                column(1),
                column(5,plotOutput("betaMDS")),
                column(5,plotOutput("betaNMDS"))
              )
            ),
            tabPanel("Phylogenetic Tree",
              tags$hr(),
              fluidRow(
                column(4,sliderInput("phylo_prune","Number of OTUs to display (max of 200 is advised):",1,2,1,1))
              ),
              fluidRow(
                column(2,selectInput("phylo_method","Visualization Method:",choices = c("sampledodge","treeonly","-"))),
                column(2,selectInput("phylo_color","Colors (scroll down for taxonomic classes):",choices = c(""))),
                column(2,selectInput("phylo_shape","Shapes:",choices = c(""))),
                column(2,selectInput("phylo_size","Size:",choices = c(""))),
                column(2,selectInput("phylo_label.tips","Label tips:",choices = c("-","taxa_names"))),
                column(1,checkboxInput("phylo_ladderize","Ladderize Phylogenetic tree:",F))
              ),
              tags$hr(),
              fluidRow(
                column(1),
                column(10,
                  plotOutput("phyloTree")
                ),
                column(1)
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
                  plotlyOutput("cutoffHist")
                ),
                column(4,offset=1,plotlyOutput("boolHeat"))
              ),
              fluidRow(column(6,htmlOutput("cutoff_text")),
                column(1),
                column(4,htmlOutput("heatmap_text"))
              ),
              tags$hr(),
              fluidRow(
                column(4,
                  htmlOutput("basic_calc_title"),
                  radioButtons("useFC","Calculation of Counts:",c("log2(fold-change+1)","difference"))
                ),
                column(4,selectInput("groupCol","Select Column from META-file containing groupings:",choices = c("Please Upload OTU & META file first!"),selected = "Please Upload OTU & META file first!")),
                column(3,actionButton("startCalc","Start Count Calculation & Reload Network!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
              ),
              tags$hr(),
              fluidRow(
                column(1),
                column(8,forceNetworkOutput("basicNetwork")),
                column(3,sliderInput("networkCutoff","Number of edges to show:",50,500,50,10),
                  selectInput("netLevel","Taxonomic Level:",choices=c("-","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
                )
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
                  sliderInput("K","Pick Number of Topics:", 1, 150, 10, step=1),
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
            ),
            tabPanel("SPIEC-EASI",
              h4("This is SPIEC-EASI"),
              tags$hr(),
              fixedRow(
                column(6, p("Meinshausen-Buhlmann's (MB) neighborhood selection:"),
                       fixedRow(
                         column(3,numericInput('se_mb_lambda',label='nlambda',value=20,min=0,max=1000,step=1)),
                         column(3,numericInput('se_mb_lambda.min.ratio',label='lambda.min.ratio',value=0.01,min=0,max=1,step=0.001)),
                         column(3,numericInput("se_mb_repnumber",label="Number of repeates",value = 50,min=1,max=1000,step=1))
                       )),
                column(6, p("inverse covariance selection (glasso):"),
                       fixedRow(
                         column(3,numericInput('se_glasso_lambda',label='nlambda',value=20,min=0,max=1000,step=1)),
                         column(3,numericInput('se_glasso_lambda.min.ratio',label='lambda.min.ratio',value=0.01,min=0,max=1,step=0.001)),
                         column(3,numericInput("se_glasso_repnumber",label="Number of repeates",value = 50,min=1,max=1000,step=1))
                       ))
                # ,
                # column(4, p("SparCC:"),
                #        fixedRow(
                #          column(3,numericInput('se_sparcc_iter',label='iterations outer loop',value=20,min=0,max=1000,step=1)),
                #          column(3,numericInput('se_sparcc_iter_inner',label='iterations inner loop',value=20,min=0,max=1000,step=1)),
                #          column(3,numericInput("se_sparcc_threshold",label="correlation matrix threshold¹",value=0.3,min=0,max=1,step=0.01))
                #        )),
              ),
              tags$hr(),
              fixedRow(
                column(6,actionButton("se_mb_start", "Start MB")),
                column(6,actionButton("se_glasso_start", "Start glasso"))
                #,column(4,actionButton("se_sparcc_start", "Start SparCC"))
              ),
              tags$hr(),
              fixedRow(
                column(6,plotOutput("spiec_easi_mb_network")),
                column(6,plotOutput("spiec_easi_glasso_network"))
                #,column(4,plotOutput("spiec_easi_sparcc_network"))
              ),
              tags$hr(),
              fixedRow(
                column(6, 
                       fluidRow(
                         column(6,selectInput("mb_select_taxa","select taxa class for mb",choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))),
                         column(1,actionButton("mb_reload_plot", "Reload Plot"))
                       )),
                column(6, 
                       fixedRow(
                         column(6,selectInput("glasso_select_taxa","select taxa class",choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))),
                         column(1,actionButton("glasso_reload_plot", "Reload Plot"))
                       ))
                # ,
                # column(4, 
                #        fixedRow(
                #          column(6,selectInput("sparcc_select_taxa","select taxa class",choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))),
                #          column(1,actionButton("sparcc_reload_plot", "Reload Plot"))
                #        ))
              ),
              tags$hr(),
              fixedRow(
                column(1),
                #column(10,htmlOutput("spiec_easi_additional")),
                column(1)
              )
            )
          )
        )
      ),
      tabItem(tabName="info",
          h4("Information & global settings"),
          fluidRow(  
            tabBox(id="info",width=12, 
               tabPanel("Data Input Format",
                        tags$hr(),
                        fixedRow(
                          column(1,''),
                          column(10,htmlOutput('info_inputdata')),
                          column(1)
                        )
                        ),
               tabPanel("Global Settings",
                        tags$hr(),
                            fixedRow(
                              column(4,numericInput("ncores","Number of cores used for calculations:",min=1,max=500,value = 1,step=1))
                            )
                        )
                   )
            )   
      )
    )
  )
)