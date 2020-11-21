library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(DT)
library(plotly)
library(networkD3)

source("texts.R")
ui <- dashboardPage(
  dashboardHeader(title="Microbiome Explorer"),
  dashboardSidebar(
    sidebarMenu(id="sidebar",
      br(),
      fluidRow(
        column(12,align="center",actionButton("upload","Load new dataset"))
      ),
      fluidRow(
        column(12,align="center",actionButton("upload_testdata","Load sample dataset"))
      ),
      dataTableOutput("datasets"),
      br(),br(),
      menuItem("Welcome!",tabName="welcome",icon=icon("door-open")),
      menuItem("Basic Analyses",tabName="basics",icon=icon("search")),
      menuItem("Advanced Analyses",tabName="advanced",icon=icon("search")),
      menuItem("Network Analysis",tabName="Network",icon=icon("project-diagram")),
      menuItem("Info & Settings",tabName = "info",icon=icon("cogs"))
    )
  ),
  dashboardBody(
    useShinyjs(),
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
                column(9,dataTableOutput("metaTable")),
                column(3,wellPanel(
                  h4("Filter options for meta-data"),
                  selectInput("filterColumns","Meta-Variables",choices = ''),
                  selectInput("filterColumnValues","Variable values",choices = ''),
                  hr(),
                  selectInput("filterTaxa","Taxa-Groups",choices = c("NONE","Kingdom","Phylum","Class","Order","Family","Genus","Species")),
                  selectInput("filterTaxaValues","Group values",choices = ''),
                  hr(),
                  actionButton("filterApply","Apply Filter",style="background-color:blue; color:white"),
                  actionButton("filterReset","Reset all Filters", style="background-color:green; color:white")
                ))
            )),
            tabPanel("Rarefaction Curves",
              tags$hr(),
              fluidRow(
                column(1),
                column(6,plotlyOutput("rarefacCurve")),
                column(1),
                column(4,br(),
                  sliderInput("rareToHighlight","Number of samples with steepest rarefaction slope to be highlighted:",min=0,value=1,max=100,step=1))
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
                column(10,wellPanel(
                  plotlyOutput("taxaDistribution",height="auto") 
                )),
                column(2,wellPanel(
                  selectInput("taxLevel","Taxonomic Level",choices=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")),
                  sliderInput("taxCutoff","Group OTUs below average abundance of: ",0,20,0)
                ))
              )
            ),
            tabPanel("Data Structure", 
              tags$hr(),
              fluidRow(
                column(8,wellPanel(
                  p("Dimensionality reduction methods"),
                  plotlyOutput("structurePlot")
                )),
                column(4,wellPanel(
                  selectInput("structureMethod","",c("PCA","UMAP","t-SNE")),
                  selectInput("structureGroup","Group by:",""),
                  radioButtons("structureDim","Dimensions:",c("2D","3D")),
                  div(id="structureCompChoosing",
                      selectInput("structureCompOne","Choose component 1 to look at:",1),
                      selectInput("structureCompTwo","Choose component 2 to look at:",1),
                      selectInput("structureCompThree","Choose component 3 to look at:",1))
                ))
              ),
              hr(),
              conditionalPanel("input.structureMethod == 'PCA'",
                fluidRow(
                  column(8, wellPanel(
                    p("Top and Bottom Loadings (show those OTUs which have the most positive (blue) or negative (red) influence on the chosen principal component)"),
                    plotlyOutput("loadingsPlot")
                  )),
                  column(4, wellPanel(
                    selectInput("pcaLoading","Plot loadings on PC",1)
                  ))
                )
            )),
            tabPanel("Confounding Analysis & Explained Variation",
                     tags$hr(),
                     h4("Confounding Analysis:"),
                     fluidRow(
                       column(1),
                       column(3,selectInput("confounding_var","Choose variable to test for confounding (variables need to have at least 2 levels; for datasets with only a single meta variable, confounding analysis can not be performed):",choices = "")),
                       column(3,disabled(actionButton("confounding_start","Start calculation..")))
                     ),
                     tags$hr(),
                     fluidRow(
                       column(1),
                       column(3,htmlOutput("confounding_var_text"))
                     ),
                     tags$hr(),
                     fluidRow(
                       column(1),
                       column(5,tableOutput("confounding_table")),
                       column(2,downloadButton("confounding_table_download","Download Table")),
                       column(1)
                     ),
                     tags$hr(),
                     h4("Explained Variation:"),
                     fluidRow(
                       column(1),
                       p("The bars represent rsquare (4 digit rounded value is written over bars) and are colored by pvalue. The rsquare value corresponds to the explained variation a variable has"),
                       column(7, plotOutput("explainedVariationBar",height = "700px"))
                     ),
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
              h4("Raw values for alpha diversity scores:"),
              fluidRow(
                column(1),
                column(7,tableOutput("alphaTable"),
                       downloadButton("alphaTableDownload","Download Table"))
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
              h3("Phylogenetic Tree of OTU taxa"),
              tags$hr(),
              fixedRow(
                column(6,wellPanel(
                  h4("Basic tree visualization options:"),
                  div(id="phylo_basic",
                      sliderInput("phylo_prune","Number of OTUs to display (pick the x OTUs with the highest cumulative abundance):",2,2,1,1))
                )),
                column(6,wellPanel(
                  h4("Advanced tree visualization options:"),
                  actionButton("phylo_toggle_advanced","Show/hide advanced options"),
                  hidden(div(id="phylo_advanced",
                      selectInput("phylo_method","Visualization Method:",choices = c("sampledodge","treeonly")),
                      selectInput("phylo_color","Group OTUs by meta samples using: colors (scroll down for taxonomic classes)",choices = c("")),
                      selectInput("phylo_shape","Group OTUs by meta samples using: shapes",choices = c("")),
                      selectInput("phylo_size","Group OTUs by meta samples using: size",choices = c("")),
                      selectInput("phylo_tiplabels","Label tips:",choices = c("-","taxa_names")),
                      sliderInput("phylo_margin","Plotmargin (defines right-handed padding; 0.2 adds 20% extra space):",0,1,0.2,0.01),
                      checkboxInput("phylo_ladderize","Ladderize Phylogenetic tree",F),
                      checkboxInput("phylo_radial","Display radial tree",F)))
                ))
              ),
              hr(),
              fixedRow(
                column(12, wellPanel(
                  div(id="phylo_tree",
                      plotOutput("phyloTree"),style="height:600px")
                ))
              )
            )
          )
        )
      ),
      tabItem(tabName = "advanced",
        h4("Advanced Analysis"),
        fluidRow(
          tabBox(id="advancedPlots",width=12,
             # tabPanel("OTU Correlation", 
             #          tags$hr(),
             #          fluidRow(
             #            column(1),
             #            column(10,plotlyOutput("OTUcorrPlot"))
             #          ),
             #          br(),br(),
             #          fluidRow(
             #            column(4),
             #            column(4,sliderInput("corrCluster","Number of Clusters:",1,20,2))
             #          )
             # ),
             tabPanel("Abundance Heatmaps",
                      tags$hr(),
                      p(heatmapText),
                      tags$hr(),
                      fluidRow(
                        column(10,wellPanel(
                          plotlyOutput("abundanceHeatmap")
                        )),
                        column(2,wellPanel(
                          selectInput("heatmapDistance","Choose distance method",choices = c("bray","gunifrac","wunifrac","unifrac","dpcoa","jsd")),
                          selectInput("heatmapOrdination","Choose Orientation Method (Ordination)",choices = c("NMDS","MDS/PCoA","DPCoA","DCA","CCA","RDA")),
                          selectInput("heatmapSample","Choose labeling of X-axis",choices="")
                        ))
                      ),
                      fluidRow(
                        column(1),
                        column(10,
                               br(),
                               htmlOutput("heatmapOrdinationText"))
                      )
             ),
             tabPanel("Random Forests",
                tags$hr(),
                fixedRow(
                  column(6, wellPanel( 
                        h2("Options for building the model:"),
                        selectInput("forest_variable","Choose variable of meta file, for which a prediction model will be built:",choices = ""),
                        plotOutput("forest_continuous_density",height = "200px"),
                        sliderInput("forest_continuous_slider","If a numeric/continuous variable was chosen, select cutoff value to transform variable into 2 distinct groups:",0,1,0,.01),
                        selectInput("forest_type","Select mode of model calculation",choices=c("random forest"),selected = "randomForest"),
                        selectInput("forest_features","Select meta-features to build model",choices = "",multiple = T),
                        checkboxInput("forest_otu","Use OTU relative abundances to predict model",T)),
                        wellPanel(p("Model parameters:"),
                                  verbatimTextOutput("forest_model_parameters"))
                  ),
                  column(6, wellPanel(
                        actionButton("forest_toggle_advanced","Show/hide advanced options"),
                        tags$hr(),
                        hidden(
                          div(id="forest_advanced",
                              div(id="general_advanced",
                                  checkboxInput("forest_clr","Perform centered log ratio transformation on OTU abundace data",F),
                                  sliderInput("forest_partition","Choose ratio of dataset which will be used for training; the rest will be used for testing the model",0,1,.75,.01),
                                  selectInput("forest_exclude","Exclude OTUs from model calculation:",choices = "",selected = NULL,multiple = T),
                                  p("Resampling options:"),
                                  selectInput("forest_resampling_method","The resampling method",choices = c("boot","cv","LOOCV","LGOCV","repeatedcv"),selected = "repeatedcv"),
                                  numericInput("forest_cv_fold","Number of folds in K-fold cross validation/Number of resampling iterations for bootstrapping and leave-group-out cross-validation",min=1,max=100,step=1,value=10),
                                  #numericInput("forest_cv_repeats","Number of repeats (Only applied to repeatedcv)",min=1,max=100,value=3,step=1),
                                  numericInput("forest_seed","Set random seed:",min=-Inf,max=Inf,value=2020,step=.000001)
                              ),
                              tags$hr(),
                              div(id="ranger_advanced",
                                  numericInput("forest_ntrees","Number of decision trees to grow",value = 500,min=1,max=10000,step=1),
                                  textInput("forest_mtry","Number of variables to possibly split at in each node (multiple entries possible, seperate by comma)","1,2,3"),
                                  selectInput("forest_splitrule","Splitting rule",choices=c("gini","extratrees","hellinger"),selected = "gini",multiple = T),
                                  textInput("forest_min_node_size","Minimal node size (multiple entries possible, seperate by comma)","1,2,3"),
                                  selectInput("forest_importance","Variable importance mode",choices=c("impurity","impurity_corrected","permutation"),selected = "impurity"),
                              ),
                              #TODO: remove placeholders!! 
                              hidden(
                                div(id="gbm_advanced",
                                    textInput("gbm_ntrees","Number of decision trees to grow (multiple entries possible)",placeholder = "Enter mutiple numbers seperated by comma.."),
                                    textInput("gbm_interaction_depth","The maximum depth of variable interactions. (multiple entries possible)",placeholder = "Enter mutiple numbers seperated by comma.."),
                                    textInput("gbm_shrinkage","The shrinkage parameter applied to each tree in the expansion (multiple entries possible)",placeholder = "Enter mutiple numbers seperated by comma.."),
                                    textInput("gbm_n_minobsinoode","Integer specifying the minimum number of observations in the trees terminal nodes (multiple entries possible)",placeholder = "Enter mutiple numbers seperated by comma.."),
                                    
                                ) 
                              )
                          )
                        ),
                        actionButton("forest_start","Start model calculation!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        checkboxInput("forest_default","Use default parameters: (toggle advanced options for more flexibility)",T)
                        
                  ))
                ),
                hr(),
                fixedRow(
                  column(5),
                  column(2,h1("Results")),
                  column(5)
                ),
                fixedRow(
                  column(4, wellPanel(
                    p("Confusion Matrix for testing-dataset"),
                    plotOutput("forest_con_matrix"),
                    p("Confusion Matrix for full dataset"),
                    plotOutput("forest_con_matrix_full")
                  )),
                  column(4, wellPanel(
                    p("ROC-Plot: TP-rate vs. FP-rate including AUC for model"),
                    plotOutput("forest_roc"),
                    plotOutput("forest_roc_cv"),
                    p("The receiver operating characteristic (ROC) can show you how good the model can distuingish between sample-groups. A perfect ROC-Curve would go from (0,0) to (0,1) to (1,1). This means the model has a perfect measure of seperability. A ROC-Curve that goes diagonally from (0,0) to (1,1) tells you, that the model makes only random predictions."),
                    p("The AUC (area under the curve) is a good measure to compare multiple ROC curves and therefore models. Here a AUC of 1 tells you, that you have a perfect model, AUC of 0.5 is again only random.")
                  )),
                  column(4,wellPanel(
                    p("Show the top x most important features for building the model"),
                    sliderInput("top_x_features","Pick x",min=1,max=100,value=20,step=1),
                    plotOutput("forest_top_features")
                  ))
                ),
                hr(),
                fixedRow(
                  column(4),
                  column(4,h1("Apply the model")),
                  column(4)
                ),
                fixedRow(
                  column(5, wellPanel(
                    p("Use model classifier for new sample(s):"),
                    p("Upload file has to be table with same columns, used in model! Columns used: "), verbatimTextOutput("forest_model_variables"),
                    fileInput("forest_upload_file","Choose file for new sample",multiple = F,placeholder = "Choose file for new sample"),
                    actionButton("forest_upload","Upload new Sample",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                    tableOutput("forest_prediction")
                  ))
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
                    column(3,numericInput("binCutoff","Cutoff for Binarization (OTUs with a smaller value are considered as not present)",min=0.01,max=10,value=1,step = 0.01)),
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
                column(1),
                column(5,htmlOutput("basic_calc_title"),
                       wellPanel(radioButtons("useFC","Calculation of Counts:",c("log2(fold-change)","difference"))),
                       actionButton("startCalc","Start Count Calculation & Reload Network!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
                column(5,wellPanel(selectInput("groupCol","Select which sample group is to be compared:",choices = c("Please Upload OTU & META file first!"),selected = "Please Upload OTU & META file first!"),
                                   selectInput("groupVar1","Select variable of group to compare with",choices = c("Please Upload OTU & META file first!")),
                                   selectInput("groupVar2","Select variable of group to compare against (choose *all* to compaire against all other variables in group)", choices = c("Please Upload OTU & META file first!"))))
              ),
              tags$hr(),
              fluidRow(
                column(1),
                column(8,forceNetworkOutput("basicNetwork")),
                column(3,sliderInput("networkCutoff","Number of edges to show:",1,5000,100,10),
                  selectInput("netLevel","Taxonomic Level:",choices=c("-","Kingdom","Phylum","Class","Order","Family","Genus","Species")),
                  plotOutput("nodeDegree"),
                  sliderInput("nodeDegreeBins","Choose number of bins for plot:",10,200,50,1)
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
            tabPanel("Topic Modeling",
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
                  actionButton("themeta","Visualize topics!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
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
                         column(3,numericInput('se_mb_lambda',label='The number of regularization/thresholding parameters (=lambda)',value=10,min=0,max=100,step=1)),
                         column(3,numericInput('se_mb_lambda.min.ratio',label='Smallest value for lambda',value=0.1,min=0,max=1,step=0.001)),
                         column(3,numericInput("se_mb_repnumber",label="Number of subsamples for StARS",value = 20,min=1,max=1000,step=1))
                       )),
                column(6, p("inverse covariance selection (glasso):"),
                       fixedRow(
                         column(3,numericInput('se_glasso_lambda',label='The number of regularization/thresholding parameters (=lambda)',value=10,min=0,max=100,step=1)),
                         column(3,numericInput('se_glasso_lambda.min.ratio',label='Smallest value for lambda',value=0.1,min=0,max=1,step=0.001)),
                         column(3,numericInput("se_glasso_repnumber",label="Number of subsamples for StARS",value = 20,min=1,max=1000,step=1))
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
                column(6,
                       fluidRow(
                         column(2, actionButton("se_mb_start", "Start MB")),
                         column(2, radioButtons("se_mb_interactive",NULL,choices = c("fixed network","interactive network"))),
                         column(3,selectInput("mb_select_taxa","select taxa class",choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
                       )),
                column(6,
                       fluidRow(
                         column(2, actionButton("se_glasso_start", "Start glasso")),
                         column(2, radioButtons("se_glasso_interactive",NULL,choices = c("fixed network","interactive network"))),
                         column(3,selectInput("glasso_select_taxa","select taxa class",choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
                       ))
                # column(4, 
                #        fixedRow(
                #          column(6,selectInput("sparcc_select_taxa","select taxa class",choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))),
                #          column(1,actionButton("sparcc_reload_plot", "Reload Plot"))
                #        ))
                #,column(4,actionButton("se_sparcc_start", "Start SparCC"))
              ),
              tags$hr(),
              fixedRow(
                column(6,
                       conditionalPanel(
                         condition = "input.se_mb_interactive == 'interactive network'",
                         p("use your mouse to zoom in and out and move the network around"),
                         forceNetworkOutput("spiec_easi_mb_network_interactive")
                       ),
                       conditionalPanel(
                         condition = "input.se_mb_interactive == 'fixed network'",
                         plotOutput("spiec_easi_mb_network")
                       ),
                ),
                column(6,
                       conditionalPanel(
                         condition = "input.se_glasso_interactive == 'interactive network'",
                         p("use your mouse to zoom in and out and move the network around"),
                         forceNetworkOutput("spiec_easi_glasso_network_interactive")
                       ),
                       conditionalPanel(
                         condition = "input.se_glasso_interactive == 'fixed network'",
                         plotOutput("spiec_easi_glasso_network")
                       ),
                )
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
               tabPanel("Testdata",
                        tags$hr(),
                        fixedRow(
                          column(1,''),
                          column(10,htmlOutput('info_testdata')),
                          column(1)
                        )
               )
            )
         )   
      )
    )
  )
)