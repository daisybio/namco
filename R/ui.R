namco_packages <- c("DT","networkD3", "shiny", "shinyjs", "waiter", "plotly",
                    "fontawesome", "shinyWidgets","shinydashboard", "shinydashboardPlus")

suppressMessages(lapply(namco_packages, require, character.only=T, quietly=T, warn.conflicts=F))

source("texts.R")
ui <- dashboardPage(
  skin="blue",
  dashboardHeader(title="Microbiome Explorer", titleWidth = 300),
  dashboardSidebar(
    sidebarMenu(id="sidebar",
      br(),
      h2("DATA UPLOAD", style="text-align:center; font-weight:1000"),
      fluidRow(
        column(12,align="center",actionButton("upload_otu","Upload OTU/ASV table", icon = icon("table"), style="color:#3c8dbc"))
      ),
      fluidRow(
        column(12,align="center",actionButton("upload_fastq","Upload raw fastq files", icon = icon("dna"), style="color:#3c8dbc"))
      ),
      fluidRow(
        column(12,align="center",actionButton("upload_testdata","Load sample dataset", icon = icon("database"), style="color:#3c8dbc"))
      ),
      hr(), br(),
      tags$style(HTML("thead {
                    color: #3c8dbc;
                    }")),
      dataTableOutput("datasets"),
      hr(),br(),
      #h3("Analysis tabs:",style="text-align:center; font-weight:500"),
      menuItem("Welcome!",tabName="welcome",icon=icon("door-open")),
      menuItem("Data Overview & Filtering",tabName="overview",icon=icon("filter")),
      menuItemOutput("fastq_overview"),
      menuItem("Basic Analysis",tabName="basics",icon=icon("search")),
      menuItem("Advanced Analysis",tabName="advanced",icon=icon("search")),
      menuItem("Network Analysis",tabName="Network",icon=icon("project-diagram")),
      menuItem("Info & Settings",tabName = "info",icon=icon("info-circle")),
      hr(),br(),
      h4("Save or restore session:",style="text-align:center; font-weight:500"),
      fluidRow(
        column(6, align="center", downloadBttn("saveSession","Save Session",size = "xs",style="float")),
        column(6, align="center", actionBttn("loadSession","Restore Session",size = "xs",style="float"))
      )
    ),
    width = 300
  ),
  dashboardBody(
    use_waiter(),
    waiter_show_on_load(html = tagList(spin_rotating_plane(),"Loading necessary packages for NAMCO ...")),
    useShinyjs(),
    tabItems(
      tabItem(tabName="welcome",
        fluidRow(
          column(2,htmlOutput("startHere")),
          column(1),
          column(6,htmlOutput("welcome"))
        ),
        fluidRow(
          column(3),
          column(9, htmlOutput("contactText"))
        ),
        br(),
        fluidRow(
          column(3),
          column(9,box(title="Authors",htmlOutput("authors"), solidHeader = T, status="primary", collapsible = T, collapsed = T))
        ),
        fluidRow(
          column(3),
          column(9,box(title="References:",htmlOutput("welcome_ref"), solidHeader = T, status="primary", collapsible = T, collapsed = T))
      )),
      tabItem(tabName="overview",
              h4("Data Overview & Filtering"),
              fluidRow(
                valueBoxOutput("otus_box1"),
                valueBoxOutput("samples_box1"),
                valueBoxOutput("conditions_box1")
              ),
              tabBox(id="filters",width=12,
                       tabPanel("Filter Samples",
                                hr(),
                                p("Explore the meta-file you uploaded. Use the filtering options to use only specific groups of samples for your analysis."),
                                fluidRow(
                                  column(10,wellPanel(
                                    dataTableOutput("metaTable")
                                  )),
                                  column(2,wellPanel(
                                    h4("Filter options for samples"),
                                    selectInput("filterColumns","Meta-Variables",choices = ''),
                                    selectInput("filterColumnValues","Variable values",choices = ''),
                                    selectInput("filterSample", "Pick specific samples by name", choices = "", multiple = T),
                                    hr(),
                                    actionButton("filterApplySamples","Apply Filter",style="background-color:blue; color:white; display:inline-block"),
                                    actionButton("filterResetA","Restore original dataset", style="background-color:green; color:white")
                                  ))
                                )
                      ),
                      tabPanel("Filter Taxa",
                              hr(),
                              p("Explore the taxonomies in the samples. Use the filtering options to use only specific taxonomic groups for your analysis."),
                              fixedRow(
                                column(10, wellPanel(
                                  plotlyOutput("taxaDistribution",height="auto")
                                )),
                                column(2, 
                                  switchInput("taxaAbundanceType","Show relative or absolute abundance",onLabel = "relative",offLabel = "absolute",value = T,size = "mini"),
                                  wellPanel(
                                  h4("Filter options for taxa"),
                                  selectInput("filterTaxa","Taxa-Groups",choices = c("Kingdom","Phylum","Class","Order","Family","Genus")),
                                  selectInput("filterTaxaValues","Group values",choices = ''),
                                  hr(),
                                  actionButton("filterApplyTaxa","Apply Filter",style="background-color:blue; color:white"),
                                  actionButton("filterResetB","Restore original dataset", style="background-color:green; color:white")
                                ))
                              )
                      )
              )
      ),
      tabItem(tabName = "fastq_overview",
              h4("fastq Overview"),
              fluidRow(
                valueBoxOutput("otus_box2"),
                valueBoxOutput("samples_box2"),
                valueBoxOutput("conditions_box2")
              ),
              fluidRow(
                tabBox(id="fastq_dada2", width=12,
                       tabPanel("Quality and Filtering",
                         h3("Analysis of sequence quality for provided fastq files after filtering"),
                         htmlOutput("dada2SourceText"),
                         
                         fluidRow(column(12,
                                  selectizeInput("fastq_file_select_filtered", label="Select fastq-pair:", multiple = F, choices=c()), 
                                  fluidRow(
                                    column(6, wellPanel(
                                      h4("foreward"),
                                      div("",plotOutput("fastq_file_quality_fw_filtered"))
                                    )),
                                    column(6, wellPanel(
                                      h4("reverse"),
                                      div("",plotOutput("fastq_file_quality_rv_filtered"))
                                    ))
                                  )
                                  )
                         ),
                         fluidRow(column(10, htmlOutput("fastqQualityText"))),
                         hr(),
                         fluidRow(
                           column(12,h3("Number of reads after each step in the DADA2 pipeline"),
                                  fluidRow(
                                    column(2),
                                    column(8, wellPanel(
                                      plotlyOutput("fastq_pipeline_readloss")
                                    ))
                                  ))
                         )
                       ),
                       tabPanel("Downloads", 
                         h3("Download options"),
                         hr(),
                         fluidRow(
                           column(12, 
                             h3("Download the generated ASV-tables:"), wellPanel(
                             fixedRow(column(4, downloadBttn("download_asv_norm","Download normalized ASV table", style="float", size="sm")),
                                      column(4, downloadBttn("download_asv_raw", "Download unnormalized ASV table", style="float", size="sm")))),
                             h3("Download the ASV sequences:"),wellPanel(
                               fixedRow(column(6, downloadBttn("download_asv_fastq","Download fasta file of ASV sequences",style="float", size="sm")))),
                             h3("Download the taxonomic classification:"),wellPanel(
                               fixedRow(column(6, downloadBttn("download_taxonomy","Download taxonomic classification of ASVs",style="float", size="sm")))),
                             h3("Download a phyloseq R-object:"), wellPanel(
                               fixedRow(column(6, downloadBttn("download_phyloseq","Download phyloseq object", style="float", size="sm"))))
                           
                         )
                       ))
                ) 
              )
      ),
      tabItem(tabName = "basics",
        h4("Basic Analysis"),
        fluidRow(  
          tabBox(id="basicPlots",width=12,
                 
             tabPanel("Alpha Diversity",
                      h3("Analyse the diversity of species inside of samples"),
                      tags$hr(),
                      htmlOutput("alphaDivText"),
                      tags$hr(),
                      fluidRow(
                        column(8,wellPanel(
                          plotlyOutput("alphaPlot") 
                        )),
                        column(3,wellPanel(
                          selectInput("alphaMethod","Method:",c("Shannon_Entropy","effective_Shannon_Entropy","Simpson_Index","effective_Simpson_Index","Richness")),
                          selectInput("alphaGroup","Group by:",c("-"))
                        ))
                      ),
                      br(),br(),
                      
                      fluidRow(column(12,
                                      h4("Raw values for alpha diversity scores, including download:"),
                                      downloadButton("alphaTableDownload","Download Table"))
                      ),
                      fluidRow(
                        column(10,wellPanel(
                          tableOutput("alphaTable"))
                        )
                      ),
                      tags$hr(),
                      fluidRow(
                        column(1),
                        column(10, htmlOutput("alphaDivFormulas"))
                      )
             ),
             
             tabPanel("Beta Diversity",
                      h3("Analyse the diversity of species between samples"),
                      tags$hr(),
                      htmlOutput("betaDivText"),
                      hr(),
                      fluidRow(
                        column(1),
                        column(7,wellPanel(plotOutput("betaTree", width = "100%"))),
                        column(3,wellPanel(
                          selectInput("betaMethod","Method:",choices=""),
                          selectInput("betaGroup","Group by:",choices=""),
                          switchInput("betaShowLabels","Show label of samples",F)
                        ))
                      ),
                      fluidRow(
                        column(1),
                        column(10, wellPanel(
                          fluidRow(                
                            column(6,plotOutput("betaMDS")),
                            column(6,plotOutput("betaNMDS"))
                          )))
                      )
             ),
             
             tabPanel("Sample comparisons", 
                      h3("Compare samples using dimensionality reduction methods"),
                      tags$hr(),
                      htmlOutput("dimReductionInfoText"),
                      tags$hr(),
                      fluidRow(
                        column(9, wellPanel(
                          p("Dimensionality reduction methods"),
                          plotlyOutput("structurePlot", height = "500px")
                        )),
                        box(
                          width=3,
                          title="Options",
                          solidHeader = T, status = "primary",
                          selectInput("structureMethod","",c("PCA","UMAP","t-SNE")),
                          selectInput("structureGroup","Group by:",""),
                          radioButtons("structureDim","Dimensions:",c("2D","3D")),
                          div(id="structureCompChoosing",
                              selectInput("structureCompOne","Choose component 1 to look at:",1),
                              selectInput("structureCompTwo","Choose component 2 to look at:",1),
                              selectInput("structureCompThree","Choose component 3 to look at:",1))
                        )
                      ),
                      hr(),
                      conditionalPanel("input.structureMethod == 'PCA'",
                                       fluidRow(
                                         column(9, wellPanel(
                                           p("Top and Bottom Loadings (show those OTUs which have the most positive (blue) or negative (red) influence on the chosen principal component)"),
                                           plotlyOutput("loadingsPlot")
                                         )),
                                         box(
                                           width=3,
                                           title="Display loadings of one PC",
                                           solidHeader = T, status = "primary",
                                           selectInput("pcaLoading","Plot loadings on PC",1),
                                           sliderInput("pcaLoadingNumber", "Number of taxa, for which loadings are displayed", 0,1,1,1)
                                         )
                                       ),
                                       fluidRow(
                                         column(9, wellPanel(
                                           p("Scree-Plot: Shows the Fraction of explained Variation for each PC; can help to identify highly variant Principal Components."),
                                           plotOutput("screePlot")
                                         )),
                                         box(
                                           width = 3,
                                           title ="Number of PCs to display in the Screeplot",
                                           solidHeader = T, status="primary",
                                           sliderInput("screePCshow","Number of PCs:", 1,10,1,10)
                                         )
                                       )
                      )),
                 
            tabPanel("Rarefaction Curves",
              h3("Analysis of species richness with rarefaction curves"),
              tags$hr(),
              htmlOutput("rarefactionInfoText"),
              tags$hr(),
              fluidRow(
                column(8,wellPanel(plotlyOutput("rarefacCurve",width = "100%"), waiter_hide_on_render(id="rarefacCurve"))),
                column(4,wellPanel(
                  sliderInput("rareToHighlight","Number of samples with steepest rarefaction slope to be highlighted:",min=0,value=1,max=100,step=1),
                  br(),
                  verbatimTextOutput("undersampled"),
                  switchInput("excludeSamples","exclude undersampled samples",value=F)
                )),
                valueBoxOutput("samples_box3")
              )
            ),
            
            tabPanel("Confounding Analysis & Explained Variation",
                     h3("Analyse confounding factors"),
                     tags$hr(),
                     htmlOutput("confoundingInfoText"),
                     tags$hr(),
                     fluidRow(
                       column(1),
                       column(3,selectInput("confounding_var","Choose variable to test for confounding (only variables with at least 2 levels are displayed):",choices = "")),
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
            
            tabPanel("Associations",
              h3("Explore different measures of association between sample groups"),
              hr(),
              htmlOutput("associationsSourceText"),
              hr(),
              htmlOutput("associationsText"),
              hr(),
              fluidRow(
                column(9, plotOutput("associationsPlot", width="100%")),
                box(
                  width=3,
                  title="Options",
                  solidHeader = T, status = "primary",
                  sliderInput("associations_alpha","Significance level",0.00001,1,0.05,0.001),
                  selectInput("associations_label", "Select meta-label, for which associations are tested", c("")),
                  selectInput("associations_case", "Select, which value is considered case (will be compared against all other values in label)", c("")),
                  sliderInput("assiciation_show_numer", "How many significant features do you want to display?",1,100,25,1),
                  selectInput("associations_sort","Select how to sort the displayed features", choices=c("p-value","fold-change","prevalence shift")),
                  selectInput("associations_panels","Which additional values do you want to display?", choices=c("fold-change","AU-ROC","prevalence"), multiple = T),
                  hr(),
                  actionBttn("associations_start", "Generate Plot...", icon=icon("play"), style = "pill", color="primary", block=T, size="md")
                )
              )
              
            ),
            
            tabPanel("Phylogenetic Tree",
              h3("Phylogenetic Tree of OTU taxa", fontawesome::fa("tree", fill="red", height="1.5em")),
              tags$hr(),
              fixedRow(
                column(6,wellPanel(
                  h4("Basic tree visualization options:"),
                  div(id="phylo_basic",
                      sliderInput("phylo_prune","Number of OTUs to display (pick the x OTUs with the highest cumulative abundance):",0,1,1,1),
                      selectInput("phylo_tiplabels","Label tips (remove OTU labels by choosing \'-\'):",choices = c("taxa_names", "-")),
                      selectInput("phylo_method","Visualization Method (\'sampledodge\': display samples, in which an OTU is present as circles; or \'treeonly\'):",choices = c("sampledodge","treeonly")),
                      selectInput("phylo_color","Group OTUs by meta samples using: colors (open advanced options to add more than one grouping)",choices = c(""))
                  )
                )),
                column(6,wellPanel(
                  h4("Advanced tree visualization options:"),
                  actionButton("phylo_toggle_advanced","Show/hide advanced options"),
                  hidden(div(id="phylo_advanced",
                      p("Additional options to insert groupings into the tree (size is also able to display the \'abundance\' values of an OTU in each sample):"),
                      selectInput("phylo_shape","Group OTUs by meta samples using: shapes",choices = c("")),
                      selectInput("phylo_size","Group OTUs by meta samples using: size",choices = c("")),
                      hr(),
                      checkboxInput("phylo_ladderize","Ladderize Phylogenetic tree",F),
                      checkboxInput("phylo_radial","Display radial tree",F)))
                ))
              ),
              hr(),
              fixedRow(
                column(12, wellPanel(
                  div(id="phylo_tree",
                      plotOutput("phyloTree"),style="height:800px")
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
             tabPanel("Abundance Heatmaps",
                      h3("Generate ecologically-organized heatmaps"),
                      hr(),
                      htmlOutput("heatmapText"),
                      hr(),
                      fluidRow(
                        column(10, wellPanel(
                          plotOutput("abundanceHeatmap")
                        )),
                        column(2,
                          box(title="Options",
                              selectInput("heatmapDistance","Choose distance method",choices = c("bray","gunifrac","wunifrac","unifrac","jsd")),
                              selectInput("heatmapOrdination","Choose Orientation Method (Ordination)",choices = c("NMDS","MDS/PCoA","DPCoA","DCA","CCA","RDA")),
                              selectInput("heatmapSample","Choose labeling of X-axis",choices=""),
                              solidHeader = T, status = "primary", width=12)
                        )
                      ),
                      fluidRow(
                        column(10,
                               br(),
                               htmlOutput("heatmapOrdinationText"),
                               br(),
                               htmlOutput("heatmapSourceText"),
                               htmlOutput("neatmapSourceText"),
                        )
                      )
             ),
             tabPanel("Random Forests",
                h3("Build a random forest machine learning model"),
                tags$hr(),
                fixedRow(
                  column(6, wellPanel( 
                        h2("Options for building the model:"),
                        selectInput("forest_variable","Choose variable of meta file, for which a prediction model will be built:",choices = ""),
                        selectInput("forest_covariable","If chosen variable has more than 2 groups, choose level in variable which will be compared against the other levels:",choices=""),
                        plotOutput("forest_sample_preview",height = "200px"),
                        hidden(div(id="forest_continuous_options",
                            radioButtons("forest_continuous_radio","If a numeric/continuous variable was chosen, select cutoff value to transform variable into 2 distinct groups:",choices = c("Mean","Median","Custom (Use slider below)"),inline = T),
                            sliderInput("forest_continuous_slider","Custom split",0,1,0,.01),
                            plotOutput("forest_continuous_preview",height = "200px")
                        )),
                        hr(),
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
                                  selectizeInput("forest_exclude","Exclude OTUs from model calculation:",choices = "",selected = NULL,multiple = T),
                                  p("Resampling options:"),
                                  selectInput("forest_resampling_method","The resampling method",choices = c("boot","cv","LOOCV","LGOCV","repeatedcv"),selected = "repeatedcv"),
                                  numericInput("forest_cv_fold","Number of folds in K-fold cross validation/Number of resampling iterations for bootstrapping and leave-group-out cross-validation",min=1,max=100,step=1,value=10),
                                  numericInput("forest_cv_repeats","Number of repeats (Only applied to repeatedcv)",min=1,max=100,value=3,step=1)
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
                        actionBttn("forest_start","Start model calculation!", icon=icon("play"), style = "pill", color="primary", block=T, size="md"),
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
                    #plotOutput("forest_roc_cv"),
                    p("The receiver operating characteristic (ROC) can show you how good the model can distuingish between sample-groups. A perfect ROC-Curve would go from (0,0) to (0,1) to (1,1). This means the model has a perfect measure of seperability. A ROC-Curve that goes diagonally from (0,0) to (1,1) tells you, that the model makes only random predictions."),
                    p("The AUC (area under the curve) is a good measure to compare multiple ROC curves and therefore models. Here a AUC of 1 tells you, that you have a perfect model, AUC of 0.5 is again only random.")
                  )),
                  column(4,wellPanel(
                    p("Show the top x most important features for building the model. You can change how many features to display by moving the slider."),
                    sliderInput("top_x_features","Pick x",min=1,max=100,value=20,step=1),
                    plotOutput("forest_top_features")
                  )),
                  downloadButton("forest_save_model","Save model object as RDS file")
                )
             ),
             tabPanel("Functional prediction",
                      h3("Functional prediction using Picrust2"),
                      tags$hr(),
                      htmlOutput("picrust2Text"),
                      tags$br(),
                      htmlOutput("picrust2SourceText"),
                      tags$hr(),
                      fluidRow(
                        column(6, wellPanel(
                          h3("Parameters & Input Files"),
                          fileInput("picrustFastaFile","Upload .fasta file with sequences for your OTUs/ASVs:", accept = c()),
                          selectInput("picrust_test_condition", "Choose condition for which to test differential abundance", choices=c()),
                          numericInput("picrust_mc_samples","Choose number of MC iterations",min=4,max=1000, value=200, step=4),
                          p("A higher number of MC iterations will increase precision of estimating the sampling error but also increase runtime. For datasets with few samples a higher value can be chosen, with more samples a lower one should be used."),
                          hr(),
                          checkboxInput("picrust_copy_number_normalization","Normalize OTU abdunances by copy-number", value=T),
                          p("Next to the functional assignment of OTUs, Picrust2 also infers the copy numbers of each 16s-rRNA gene per OTU; you have the option to normalize your abundance values with the copy-numbers by selecting this checkbox.")
                        )),
                        column(6,
                               actionBttn("picrust2Start", "Start picrust2 & differential analysis!", icon=icon("play"), style = "pill", color="primary", block=T, size="lg"),
                               wellPanel(
                                 p("Download zip-archive with raw picrust2 results:"),
                                 h4("Please be aware:"),
                                 p("This will create a zip archive of all output files, so it might take a few seconds until the download window appears!"),
                                 p("This download window will not appear if you use a restored dataset!"),
                                 hidden(div(id="download_picrust_div",
                                            downloadButton("download_picrust_raw", "Download picrust2 results as zip archive:")
                               )
                               ))
                        )
                      ),
                      hr(),
                      h3("Differential functional analysis"),
                      htmlOutput("aldexSourceText"),
                      hr(),
                      fluidRow(column(12,
                        tabBox(
                          title="Relationships between effect size & p-value",
                          id="picrust_tabBox", width=12,
                          tabPanel("EC",
                                   fluidRow(
                                     column(6, div("",plotOutput("picrust_ec_effect_plot"))),
                                     column(6, div("",plotOutput("picrust_ec_vulcano_plot")))
                                   ),
                                   fluidRow(
                                     column(8,
                                       box(title="Names of significant functions",
                                           verbatimTextOutput("picrust_ec_effect_signif"),
                                           status = "info", width = 12, collapsible = T, collapsed = T)
                                     ),
                                     column(4, valueBoxOutput("picrust_ec_effect_signif_value"))
                                   ),
                                   h3("Details about significant functions:"),
                                   hr(),
                                   fluidRow(
                                     column(10, plotOutput("picrust_ec_signif_plot")),
                                     column(2, numericInput("picrust_ec_signif_plot_show", "Maximum number of displayed ECs", 20, min=1, max=100,step=1))
                                   ),
                                   p("Here the functions with BH adjusted P-value above the significance threshold are displayed; the boxplot shows the different abundance distributions of a function colored by each sample group. Also the BH-adjusted P-value and effect size is displayed as a barplot.")
                          ),
                          tabPanel("KO",
                                   fluidRow(
                                     column(6, div("",plotOutput("picrust_ko_effect_plot"))),
                                     column(6, div("",plotOutput("picrust_ko_vulcano_plot")))
                                   ),
                                   fluidRow(
                                     column(8,
                                            box(title="Names of significant functions",
                                                verbatimTextOutput("picrust_ko_effect_signif"),
                                                status = "info", width = 12, collapsible = T, collapsed = T)
                                     ),
                                     column(4, valueBoxOutput("picrust_ko_effect_signif_value"))
                                   ),
                                   h3("Details about significant functions:"),
                                   fluidRow(
                                     column(10, plotOutput("picrust_ko_signif_plot")),
                                     column(2, numericInput("picrust_ko_signif_plot_show", "Maximum number of displayed KOs", 20, min=1, max=100,step=1))
                                   ),
                                   p("Here the functions with BH adjusted P-value above the significance threshold are displayed; the boxplot shows the different abundance distributions of a function colored by each sample group. Also the BH-adjusted P-value and effect size is displayed as a barplot.")
                          ),
                          tabPanel("PW",
                                   fluidRow(
                                     column(6, div("",plotOutput("picrust_pw_effect_plot"))),
                                     column(6, div("",plotOutput("picrust_pw_vulcano_plot")))
                                   ),
                                   fluidRow(
                                     column(8,
                                            box(title="Names of significant functions",
                                                verbatimTextOutput("picrust_pw_effect_signif"),
                                                status = "info", width = 12, collapsible = T, collapsed = T)
                                     ),
                                     column(4, valueBoxOutput("picrust_pw_effect_signif_value"))
                                   ),
                                   h3("Details about significant functions:"),
                                   fluidRow(
                                     column(10, plotOutput("picrust_pw_signif_plot")),
                                     column(2, numericInput("picrust_pw_signif_plot_show", "Set max. number of displayed PWs", 20, min=1, max=100,step=1))
                                   ),
                                   p("Here the functions with BH adjusted P-value above the significance threshold are displayed; the boxplot shows the different abundance distributions of a function colored by each sample group. Also the BH-adjusted P-value and effect size is displayed as a barplot.")
                          ),
                          tabPanel("Information & Options",
                                   fluidRow(
                                     column(6, wellPanel(
                                       h3("Options for Visualization"),
                                       numericInput("picrust_signif_lvl","Change significance level (P-value)",min=0.01,max=1,value=0.05,step=0.01),
                                       p("Here you can set the significance cutoff for the BH adjusted P-value; functions with a p-value below it are considered significant."),
                                       numericInput("picrust_signif_lvl_effect","Change significance level (effect size)",min=-10,max=10,value=1,step=0.01),
                                       p("Here you can set the significance cutoff for the effect size; functions with a effect size greater than it it are considered significant. "),
                                       sliderInput("picrust_maxoverlaps", "Change number of overlaps for point labels",min=1, max=500, value=50,step=1),
                                       checkboxInput("picrust_signif_label","Label significant functions in scatterplots", value = F),
                                       p("If too many points in close proximity are considered significant, change the number of overlaps, to display more labels.")
                                     )),
                                     column(6, wellPanel(
                                       h3("Information"),
                                       htmlOutput("picrust_pval_info_text")
                                     ))
                                   )
                          ),
                          tabPanel("Downloads",
                                   fluidRow(
                                     column(5, 
                                       wellPanel(downloadBttn("picrust_download_ec","Download differential analysis of EC numbers (EC)",style="float", size="sm")),
                                       wellPanel(downloadBttn("picrust_download_ko","Download differential analysis of KEGG ortholog groups (KO)",style="float", size="sm")),
                                       wellPanel(downloadBttn("picrust_download_pw","Download differential analysis of pathways (PW)",style="float", size="sm"))
                                     )
                                   )
                          )
                        ),
                               
                      ))
                      )
          )
        )
        ),
      tabItem(tabName = "Network",
        h4("Network Analysis"),
        fluidRow(  
          tabBox(id="netWorkPlots",width=12,
            tabPanel("Co-occurrence of OTUs",
              h3("Co-occurrences network generation"),
              htmlOutput("basic_info"),
              tags$hr(),
              htmlOutput("cutoff_title"),
              fluidRow(
                column(1),
                column(5, numericInputIcon("binCutoff","Cutoff for Binarization (OTUs with a smaller value are considered as not present)",min=0.0000001,max=100,value=1,step = 0.01, icon=icon("cut"), width="400px")),
                column(6, box(title="Information",
                             htmlOutput("basic_additional"),
                             solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T))
              ),
              fixedRow(
                column(6,                
                       box(title="Cutoff-Barplot",
                           plotlyOutput("cutoffHist"),
                           htmlOutput("cutoff_text"),
                           solidHeader = T, status="primary", collapsible = T, collapsed = F, width=12)),
                column(6,
                       box(title="Cutoff-Heatmap",
                           plotlyOutput("boolHeat"),
                           htmlOutput("heatmap_text"),
                           solidHeader = T, status="primary", collapsible = T, collapsed = F, width=12))

              ),
              tags$hr(),
              htmlOutput("basic_calc_title"),
              fluidRow(
                column(6,
                       wellPanel(
                         radioButtons("useFC","Calculation of Counts:",c("log2(fold-change)","difference")),
                         selectInput("groupCol","Select which sample group is to be compared (minimum of 2 levels in group!):",choices = c("Please Upload OTU & META file first!"),selected = "Please Upload OTU & META file first!"),
                         selectInput("groupVar1","Select variable of group to compare with (reference group)",choices = c("Please Upload OTU & META file first!")),
                         selectInput("groupVar2","Select variable of group to compare against (choose *all* to compaire against all other variables in group)", choices = c("Please Upload OTU & META file first!"))
                       ),
                actionBttn("startCalc"," Start count-calculation & (re-)load network!",icon = icon("play"), style = "pill", color="primary", block=T, size="lg")
                ),
                column(6,
                       box(title="Information",
                           htmlOutput("basic_calc_additional"),
                           solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T))
              ),
              tags$hr(),
              htmlOutput("basic_network_title"),
              fluidRow(
                column(9,forceNetworkOutput("basicNetwork")),
                column(3,
                       box(title="Options",
                           sliderInput("networkCutoff","Number of edges to show (edges are sorted by most extreme values, positive and negative):",1,5000,100,10),
                           selectInput("netLevel","Color Taxonomic Level:",choices=c("-","Kingdom","Phylum","Class","Order","Family","Genus","Species")),
                           solidHeader = T, status = "primary", width=12
                       ),
                       box(title = "chosen parameters",
                           htmlOutput("chosen_network_params"),
                           solidHeader = F, status="info", width=12, collapsible = T, collapsed = T),
                       p("Green edges: OTU-pair more often occuring in selected reference group."),
                       p("Red edges: OTU-pair more often occuring in other sample group, which is compared against.")
                )
              )
            ),
            tabPanel("Topic Modeling",
              h3("Explore thematic structure using themetagenomics"),
              htmlOutput("advanced_text"),
              tags$hr(),
              fluidRow(
                column(3,
                  sliderInput("K","Pick Number of Topics:", 1, 150, 10, step=1),
                  htmlOutput("topic_text")),
                #column(3,
                #  sliderInput("sigma_prior","Pick Scalar between 0 and 1:", 0, 1, 0, step=0.01),
                #  htmlOutput("sigma_text")),
                column(3,
                  selectInput("formula", label = "Formula for covariates of interest found in metadata:",choices="Please provide OTU-table & metadata first!"),
                  selectInput("refs", label = "Binary covariate in formula, indicating the reference level:",choices = "Please provide OTU-table & metadata first!")),
                column(3,
                  #actionButton("themeta","Visualize topics!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                  actionBttn("themeta","Visualize topics!", icon=icon("play"), style="pill", size="md",color="primary")
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
              #  column(4,selectInput('choose', label='Covariate',
              #    choices="Please start calculation above first!"),
              #fixedRow(column(1),
              #  column(11,tags$div(paste0('Choosing a covariate determines which weight estimates will shown',
              #    ' The order of the topics will be adjusted accordingly. By clicking',
              #    ' an estimate, all figures below will rerender.'),class='side')))
              #),
                column(10,plotlyOutput('est',height='200px')),
                p("Topics colored red, have a strong association with the chosen reference level; 
                  the blue topics on the other hand are associated with the other level within the chosen covariate.
                  (Example: Chosen covariate is Gender and reference level is Female -> Female will be colored red, Male is blue)")
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
              forceNetworkOutput('corr')
            ),
            
            tabPanel("SPIEC-EASI",
              h3("SPIEC-EASI microbial ecological networks"),
              htmlOutput("spiecEasiInfoText"),
              br(),
              htmlOutput("spiecEasiSource"),
              br(),
              fluidRow(
                column(8, box(title="Parameter-information",
                              htmlOutput("spiecEasiParamsText"),
                              solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T))
              ),
              hr(),
              fluidRow(
                column(6, wellPanel(
                    h4("Meinshausen-Buhlmann's (MB) neighborhood selection:"),
                    fluidRow(
                      column(4,numericInput('se_mb_lambda',label='number of lambdas',value=10,min=0,max=100,step=1)),
                      column(4,numericInput('se_mb_lambda.min.ratio',label='Smallest value for lambda',value=0.1,min=0,max=1,step=0.001)),
                      column(4,numericInput("se_mb_repnumber",label="Number of subsamples for StARS",value = 20,min=1,max=1000,step=1))
                    ),
                    hr(),
                    fluidRow(
                      column(4, actionBttn("se_mb_start","Start MB", icon=icon("play"), style="pill", size="md",color="primary")),
                      column(4, radioButtons("se_mb_interactive",NULL,choices = c("fixed network","interactive network"))),
                      column(4, selectInput("mb_select_taxa","Color by taxonomic level",choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
                    ),
                    hr(),
                    fluidRow(
                      column(12, 
                            conditionalPanel(
                              condition = "input.se_mb_interactive == 'interactive network'",
                              p("use your mouse to zoom in and out and move the network around"),
                              forceNetworkOutput("spiec_easi_mb_network_interactive")
                            ),
                            conditionalPanel(
                              condition = "input.se_mb_interactive == 'fixed network'",
                              plotOutput("spiec_easi_mb_network")
                            ))
                    )
                  )),
                column(6,wellPanel(
                    h4("inverse covariance selection (glasso):"),
                    fluidRow(
                      column(4,numericInput('se_glasso_lambda',label='number of lambdas',value=10,min=0,max=100,step=1)),
                      column(4,numericInput('se_glasso_lambda.min.ratio',label='Smallest value for lambda',value=0.1,min=0,max=1,step=0.01)),
                      column(4,numericInput("se_glasso_repnumber",label="Number of subsamples for StARS",value = 20,min=1,max=1000,step=1))
                    ),
                    hr(),
                    fluidRow(
                      column(4, actionBttn("se_glasso_start","Start glasso", icon=icon("play"), style="pill", size="md",color="primary")),
                      column(4, radioButtons("se_glasso_interactive",NULL,choices = c("fixed network","interactive network"))),
                      column(4, selectInput("glasso_select_taxa","Color by taxonomic level",choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
                    ),
                    hr(),
                    fluidRow(
                      column(12,
                             conditionalPanel(
                               condition = "input.se_glasso_interactive == 'interactive network'",
                               p("use your mouse to zoom in and out and move the network around"),
                               forceNetworkOutput("spiec_easi_glasso_network_interactive")
                             ),
                             conditionalPanel(
                               condition = "input.se_glasso_interactive == 'fixed network'",
                               plotOutput("spiec_easi_glasso_network")
                             ))
                    )
                  ))
              ),
              tags$hr(),
              fixedRow(
                column(1),
                #column(10,htmlOutput("spiec_easi_additional")),
                
                column(1)
              )
            ),
            
            tabPanel("Taxonomic Rank Networks",
              h3("Explore network structures between the discovered taxonomic ranks"),
              htmlOutput("taxNetworkInfoText"),
              br(),
              htmlOutput("diffNetworkSourceCopy"),
              br(),
              fluidRow(
                column(8, box(title="Parameter-information",
                              htmlOutput("taxNetworkParamsText"),
                              solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T))
              ),
              hr(),
              fluidRow(column(12,
                wellPanel(
                  fluidRow(
                    column(3, 
                           selectInput("taxNetworkRank","Select taxonomic rank", choices=c("Kingdom","Phylum","Class","Order","Family","Genus")),
                           selectInput("taxNetworkMeasure","Choose the measure used for calculation of network", choices=c("spring","pearson","spearman","spieceasi","bicor","sparcc","euclidian","bray","jsd")),
                    ),
                    column(3,
                           selectInput("taxNetworkClustMethod","Choose method how to detect clusters in network",choices=c("cluster_fast_greedy", "hierarchical")),
                           selectInput("taxNetworkNormMethod","Choose normalization method (in order to make counts of different samples comparable)",choices=c("none","mclr","clr","rarefy","TSS")),
                           selectInput("taxNetworkzeroMethod", "Choose method how to replace zeros in data", choices = c("none","add pseudocount of 1 to data"="pseudo","multireplicative replacement"="multRepl"))
                    ),
                    box(width=4,
                        title="Additional Parameters",
                        solidHeader = T, status = "info", collapsed = T, collapsible = T,
                        hidden(div(id="taxNetworkAdditionalParamsSPRING.EASIDiv",
                                   numericInput("taxNetworkNlambda", "Number of lambdas", 10,1,100,1),
                                   numericInput("taxNetworkRepNum", "Number of subsamples for StARS", 20,1,100,1),
                                   numericInput("taxNetworkLambdaRatio", "Smallest value for lambda", 0.1,0,1,0.01)
                        )),
                        hidden(div(id="taxNetworkAdditionalParamsSPARCCdiv",
                                   numericInput("taxNetworkIter", "Number of iterations in outer loop", 20,1,100,1),
                                   numericInput("taxNetworkInnerIter", "Number of iterations in inner loop", 10,1,100,1),
                                   numericInput("taxNetworkTh", "Threshold for correlations", 0.1,0,1,0.01)
                        ))
                    ),
                    column(2,
                           actionBttn("taxNetworkCalculate","Start Calculation",style="pill", size="lg",color="primary")       
                    )
                  )
                ))),
              hr(),
              fluidRow(
                column(9, div(id="tax_network",
                              plotOutput("taxNetwork"), style="height:2000px")),
                box(width=3,
                    title="Display options",
                    solidHeader = T, status = "primary",
                    selectInput("taxNetworkLayout","Layout",choices = c("spring", "circle", "Fruchterman-Reingold"="layout_with_fr")),
                    hr(),
                    h4("Node options:"),
                    #selectInput("taxNetworkNodeColor","Choose how to color nodes", choices=c("by detected clusters"="cluster", "Kingdom", "Phylum")),
                    selectInput("taxNetworkNodeFilterMethod", "Choose method how to filter out nodes (keep top x nodes with ...)", choices=c("none","highestConnect","highestDegree", "highestBetween", "highestClose", "highestEigen")),
                    numericInput("taxNetworkNodeFilterValue", "Choose x for the node filtering method", value = 100,min = 1,max = 10000,step = 1),
                    selectInput("taxNetworkRmSingles","How to handle unconnected nodes (all: remove all; inboth: only if unconnected in both networks, none: remove no unconnected nodes)", choices=c("inboth","none","all")),
                    selectInput("taxNetworkNodeSize", "Choose value which indicates the size of nodes", choices = c("fix","degree","betweenness","closeness","eigenvector","counts","normCounts","clr","mclr","rarefy","TSS")),
                    hr(),
                    h4("Edge options:"),
                    selectInput("taxNetworkEdgeFilterMethod","Choose method how to filter out edges (threshold: keep edges with weight of at least x; highestWeight: keep first x edges with highest weight)", choices = c("none","threshold","highestWeight")),
                    numericInput("taxNetworkEdgeFilterValue","Choose x for edge filtering method",value=300,min=1,max=5000,step=1)
                )
              ),
              hr(),
              h3("Details about network difference"),
              #TODO: print output of net_compare
              
            ),
            
            tabPanel("Differential Networks",
              h3("Explore network structures in different sample groups"),
              htmlOutput("diffNetworkInfoText"),
              br(),
              htmlOutput("diffNetworkSource"),
              br(),
              fluidRow(
                column(8, box(title="Parameter-information",
                              htmlOutput("diffNetworkParamsText"),
                              solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T))
              ),
              hr(),
              fluidRow(column(12,
                wellPanel(
                  fluidRow(
                    column(3, 
                           selectInput("diffNetworkSplitVariable","Choose sample group you want to compare (only groups with 2 levels are shown)",choices=c()),
                           selectInput("diffNetworkMeasure","Choose the measure used for calculation of network", choices=c("spring","pearson","spearman","spieceasi","bicor","sparcc","euclidian","bray","jsd")),
                           #selectInput("diffNetworkSparsMethod","Choose method used for sparsification (how to select subset of edges that are connected in network)",choices=c("none","t-test",""))
                    ),
                    column(3, 
                           selectInput("diffNetworkClustMethod","Choose method how to detect clusters in network",choices=c("cluster_fast_greedy", "hierarchical")),
                           selectInput("diffNetworkNormMethod","Choose normalization method (in order to make counts of different samples comparable)",choices=c("none","mclr","clr","rarefy","TSS")),
                           selectInput("diffNetworkzeroMethod", "Choose method how to replace zeros in data", choices = c("none","add pseudocount of 1 to data"="pseudo","multireplicative replacement"="multRepl"))
                           #numericInput("diffNetworkSparsMethodParams","A Students t-test is used to select a subset of edges which are connected; choose significance level here",value = 0.05,min = 0.001,max=1,step = 0.001)
                    ),
                    box(width=4,
                      title="Additional Parameters",
                      solidHeader = T, status = "info", collapsed = T, collapsible = T,
                      hidden(div(id="diffNetworkAdditionalParamsSPRING.EASIDiv",
                                 numericInput("diffNetworkNlambda", "Number of lambdas", 10,1,100,1),
                                 numericInput("diffNetworkRepNum", "Number of subsamples for StARS", 20,1,100,1),
                                 numericInput("diffNetworkLambdaRatio", "Smallest value for lambda", 0.1,0,1,0.01)
                      )),
                      hidden(div(id="diffNetworkAdditionalParamsSPARCCdiv",
                                 numericInput("diffNetworkIter", "Number of iterations in outer loop", 20,1,100,1),
                                 numericInput("diffNetworkInnerIter", "Number of iterations in inner loop", 10,1,100,1),
                                 numericInput("diffNetworkTh", "Threshold for correlations", 0.1,0,1,0.01)
                      ))
                    ),
                    column(2,
                           actionBttn("diffNetworkCalculate","Start Calculation",style="pill", size="lg",color="primary")       
                    )
                  )
                ))
              ),
              hr(),
              fluidRow(
                column(9, div(id="diff_network",
                              plotOutput("diffNetwork"), style="height:2000px")),
                box(width=3,
                    title="Display options",
                    solidHeader = T, status = "primary",
                    selectInput("diffNetworkLayout","Layout",choices = c("spring", "circle",  "Fruchterman-Reingold"="layout_with_fr")),
                    hr(),
                    h4("Node options:"),
                    selectInput("diffNetworkNodeFilterMethod", "Choose method how to filter out nodes (keep top x nodes with ...)", choices=c("none","highestConnect","highestDegree", "highestBetween", "highestClose", "highestEigen")),
                    numericInput("diffNetworkNodeFilterValue", "Choose x for the node filtering method", value = 100,min = 1,max = 10000,step = 1),
                    selectInput("diffNetworkRmSingles","How to handle unconnected nodes (all: remove all; inboth: only if unconnected in both networks, none: remove no unconnected nodes)", choices=c("inboth","none","all")),
                    selectInput("diffNetworkNodeSize", "Choose value which indicates the size of nodes", choices = c("fix","degree","betweenness","closeness","eigenvector","counts","normCounts","clr","mclr","rarefy","TSS")),
                    hr(),
                    h4("Edge options:"),
                    selectInput("diffNetworkEdgeFilterMethod","Choose method how to filter out edges (threshold: keep edges with weight of at least x; highestWeight: keep first x edges with highest weight)", choices = c("none","threshold","highestWeight")),
                    numericInput("diffNetworkEdgeFilterValue","Choose x for edge filtering method",value=300,min=1,max=5000,step=1)
                )
              ),
              hr(),
              h3("Details about network difference"),
              #TODO: print output of net_compare
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