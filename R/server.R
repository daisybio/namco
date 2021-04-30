namco_packages <- c("ade4", "data.table", "cluster", "DT", "fpc", "GUniFrac",
                    "heatmaply", "networkD3", "klaR", "phangorn", "plotly",
                    "RColorBrewer", "reshape2", "Rtsne", "shiny", "textshape",
                    "tidyr", "umap", "themetagenomics", "igraph",
                    "Matrix", "phyloseq", "NbClust", "caret", "ranger", "gbm",
                    "shinyjs", "MLeval", "Rcpp", "MLmetrics", "mdine", "biomformat",
                    "waiter", "dada2", "Biostrings", "fontawesome", "shinyWidgets",
                    "shinydashboard", "shinydashboardPlus", "proxy", "parallel",
                    "DECIPHER", "SpiecEasi", "ALDEx2","ggrepel", "SIAMCAT",
                    "NetCoMi","metagMisc")
#renv::snapshot(packages= namco_packages)

#TODO:
#: phylo-tree not working after dada2.. [check]
#: https://github.com/rstudio/DT/issues/758
#: normalize OTU-table by copy numbers of 16S gene per OTU -> picrust2 output (marker_predicted_and_nsti.tsv.gz) [check]
#: SIAMCAT [check]
#: https://github.com/stefpeschel/NetCoMi
#: save/load namco sesseion [check]
#: allow different file-encodings [check]

suppressMessages(lapply(namco_packages, require, character.only=T, quietly=T, warn.conflicts=F))
overlay_color="rgb(51, 62, 72, .5)"
tree_logo <- fa("tree", fill="red", height="1.5em")  #indication logo where phylo-tree is needed

server <- function(input,output,session){
  waiter_hide()
  options(shiny.maxRequestSize=10000*1024^2,stringsAsFactors=F)  # upload up to 10GB of data
  source("algorithms.R")
  source("utils.R")
  source("texts.R")
  source("file_handlings.R")
  message(log_startText)
  
  vals = reactiveValues(datasets=list(),undersampled=c()) # reactiveValues is a container for variables that might change during runtime and that influence one or more outputs, e.g. the currently selected dataset
  currentSet = NULL # a pointer to the currently selected dataset
  ncores = 4  # number of cores used where it is possible to use multiple
  seed = 123 # Global variable to use as seed
  session$onSessionEnded(stopApp) #automatically stop app, if browser window is closed
  sample_column = "SampleID"    # the column with the sample IDs will be renamed to this 
  if(!interactive()){sink(stderr(), type="output")} #this makes it so that print statements and other stdOut is saved in log file

  
  #####################################
  #    save & restore session         #
  #####################################
  # Download current session
  output$saveSession <- downloadHandler(
    filename <- function(){
      paste("namco_session.RData")
    },
    
    content = function(file) {
      message(paste0(Sys.time(), " - Saving session ..."))
      session_lst <- vals$datasets[[currentSet()]]
      save(session_lst, file = file)
    }
  )
  
  # upload RData object of namco session
  uploadSessionModal <- function(failed=F,error_message=NULL) {
    modalDialog(
      title = "Restore previous namco session",
      h4("Select session-file:"),
      fluidRow(
        column(6,wellPanel(fileInput("sessionFile","Select file"), style="background:#3c8dbc"))
      ),
      if(failed) {
        div(tags$b(error_message,style="color:red;"))
      },
      footer = tagList(
        modalButton("Cancel", icon = icon("times-circle")),
        actionButton("upload_session_ok","OK",style="background-color:blue; color:white")
      ),
      easyClose = T, fade = T, size = "l"
    )
  }
  
  observeEvent(input$loadSession, {
    showModal(uploadSessionModal())
  })
  
  observeEvent(input$upload_session_ok, {
    message(paste0(Sys.time(), " - Restoring previous session ..."))
    tryCatch({
      load(input$sessionFile$datapath)
      session_name <- session_lst[["session_name"]]
      if(session_name%in%names(vals$datasets)){stop(duplicateSessionNameError,call. = F)}
      vals$datasets[[session_name]] <- session_lst
      vals$datasets[[session_name]]$is_restored <- T
      updateTabItems(session,"sidebar")
      removeModal()
    },error=function(e){
      print(e)
      showModal(uploadSessionModal(failed=T,error_message = e))
    })
  })
  
  #####################################
  #    menu items & info-box          #
  #####################################
  
  output$fastq_overview <- renderMenu({
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$is_fastq){
        menuItem("fastq Overview", tabName = "fastq_overview", icon=icon("dna"))
      }
    }
  })
  
  output$samples_box1 <- renderValueBox({
    samples <- 0
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        samples <- dim(vals$datasets[[currentSet()]]$metaData)[1]
      }else if(vals$datasets[[currentSet()]]$is_fastq){
        samples <- length(vals$datasets[[currentSet()]]$generated_files[["sample_names"]])
      }
    }
    valueBox(samples, "Samples",icon = icon("list"), color="blue")
  })
  
  output$conditions_box1 <- renderValueBox({
    groups <- 0
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        groups <- dim(vals$datasets[[currentSet()]]$metaData)[2]
      }
    }
    valueBox(groups, "Groups", icon = icon("list"), color="purple")
  })
  
  output$otus_box1 <- renderValueBox({
    groups <- 0
    if(!is.null(currentSet())){
      otus <- ntaxa(vals$datasets[[currentSet()]]$phylo)
    }
    valueBox(otus, "ASVs/OTUs", icon = icon("list"), color="orange")
  })
  
  output$samples_box2 <- renderValueBox({
    samples <- 0
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$is_fastq){
        samples <- length(vals$datasets[[currentSet()]]$generated_files[["sample_names"]])
      }
    }
    valueBox(samples, "Samples", icon = icon("list"), color="blue")
  })
  
  output$conditions_box2 <- renderValueBox({
    groups <- 0
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        groups <- dim(vals$datasets[[currentSet()]]$metaData)[2]
      }
    }
    valueBox(groups, "Groups", icon = icon("list"), color="purple")
  })
  
  output$otus_box2 <- renderValueBox({
    groups <- 0
    if(!is.null(currentSet())){
      otus <- ntaxa(vals$datasets[[currentSet()]]$phylo)
    }
    valueBox(otus, "ASVs", icon = icon("list"), color="orange")
  })
  
  output$samples_box3 <- renderValueBox({
    samples <- 0
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        samples <- dim(vals$datasets[[currentSet()]]$metaData)[1]
      }else if(vals$datasets[[currentSet()]]$is_fastq){
        samples <- length(vals$datasets[[currentSet()]]$generated_files[["sample_names"]])
      }
    }
    valueBox(subtitle="Samples", value=samples,icon = icon("list"), color="blue", width=6)
  })
  
  
  #####################################
  #    observers & datasets           #
  #####################################
  
    #observer to show only tabs with / without meta file
  observe({
    if(!is.null(currentSet())){
      if(!vals$datasets[[currentSet()]]$has_meta){
        hideTab(inputId = "filters", target = "Filter Samples")
        hideTab(inputId = "basicPlots", target = "Confounding Analysis & Explained Variation")
        hideTab(inputId = "basicPlots", target = "Associations")
        hideTab(inputId = "basicPlots", target = "Beta Diversity")
        hideTab(inputId = "advancedPlots", target = "Abundance Heatmaps")
        hideTab(inputId = "advancedPlots", target = "Random Forests")
        hideTab(inputId = "advancedPlots", target = "Functional prediction")
        hideTab(inputId = "netWorkPlots", target = "Co-occurrence of OTUs")
        hideTab(inputId = "netWorkPlots", target = "Topic Modeling")
        hideTab(inputId = "netWorkPlots", target = "SPIEC-EASI")
        hideTab(inputId = "netWorkPlots", target = "Co-occurrence of OTUs")
      }else{
        showTab(inputId = "filters", target = "Filter Samples")
        showTab(inputId = "basicPlots", target = "Confounding Analysis & Explained Variation")
        showTab(inputId = "basicPlots", target = "Associations")
        showTab(inputId = "basicPlots", target = "Beta Diversity")
        showTab(inputId = "advancedPlots", target = "Abundance Heatmaps")
        showTab(inputId = "advancedPlots", target = "Random Forests")
        showTab(inputId = "advancedPlots", target = "Functional prediction")
        showTab(inputId = "netWorkPlots", target = "Co-occurrence of OTUs")
        showTab(inputId = "netWorkPlots", target = "Topic Modeling")
        showTab(inputId = "netWorkPlots", target = "SPIEC-EASI")
        showTab(inputId = "netWorkPlots", target = "Co-occurrence of OTUs")
      }
    }
  })
  
  #observer for fastq-related stuff
  observe({
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$is_fastq){
        updateSelectInput(session, "fastq_file_select_raw", choices = vals$datasets[[currentSet()]]$generated_files$sample_names)
        updateSelectInput(session, "fastq_file_select_filtered", choices = vals$datasets[[currentSet()]]$generated_files$sample_names)
      }
    }
    #if(!is.null(fastqSampleNamesReact())){
    #  updateSelectInput()
    #}
  })
  
  # filter variables
  observe({
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        phylo <- vals$datasets[[currentSet()]]$phylo
        meta <- data.frame(sample_data(phylo))
        taxonomy <- data.frame(tax_table(phylo))
        
        #filter variables
        filterColumnValues <- unique(meta[[input$filterColumns]])
        updateSelectInput(session,"filterColumnValues",choices=filterColumnValues)
        
        filterTaxaValues <- unique(taxonomy[[input$filterTaxa]])
        updateSelectInput(session,"filterTaxaValues",choices=filterTaxaValues) 
      }
    }
  })
  
  #observer for inputs depending on choosing a meta-group first
  observe({
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        phylo <- vals$datasets[[currentSet()]]$phylo
        meta <- data.frame(sample_data(phylo))
        
        #display sample names which can be filtered
        if(input$filterColumns == "NONE"){
          samples_left <- meta[[sample_column]]
        }else if(input$filterColumns != "" && input$filterColumnValues != ""){
          samples_left <- meta[meta[eval(input$filterColumns)] == input$filterColumnValues,][[sample_column]]
        }else{
          samples_left <- NULL
        }
        updateSelectInput(session,"filterSample",choices=samples_left)
        
        #associations
        caseVariables <- unique(meta[[input$associations_label]])
        updateSelectInput(session, "associations_case", choices=c(caseVariables))
        
        #basic network variables
        groupVariables <- unique(meta[[input$groupCol]])
        updateSelectInput(session,"groupVar1",choices = c(groupVariables))
      }
    }
  })
  
  # update input selections
  observe({
    if(!is.null(currentSet())){ 
      phylo <- vals$datasets[[currentSet()]]$phylo
      otu <- otu_table(phylo)
      taxonomy <- data.frame(tax_table(phylo))
      
      updateSelectInput(session,"structureCompOne",choices=(1:nsamples(phylo)))
      updateSelectInput(session,"structureCompTwo",choices=(2:nsamples(phylo)))
      updateSelectInput(session,"structureCompThree",choices=(3:nsamples(phylo)))
      updateSelectInput(session,"pcaLoading",choices=(1:nsamples(phylo)))
      updateSliderInput(session, "pcaLoadingNumber", min=2, max=ntaxa(phylo)-1, value=ifelse(ntaxa(phylo)<10, ntaxa(phylo), 10), step=2)
      if(!vals$datasets[[currentSet()]]$is_fastq){updateSelectInput(session,"filterTaxa", choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))}
      
      if(vals$datasets[[currentSet()]]$has_meta){
        #get tables from phyloseq object
        meta <- data.frame(sample_data(phylo))
        if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
        phylo <- vals$datasets[[currentSet()]]$phylo
        
        updateSliderInput(session,"rareToShow",min=1,max=ncol(otu),value=min(50,ncol(otu)))
        updateSliderInput(session,"rareToHighlight",min=1,max=ncol(otu),value=round(ncol(otu)/10))
        updateSliderInput(session,"top_x_features",min=1,max=nrow(otu))
        if(ncol(meta)>2){enable("confounding_start")}
        if(is.null(vals$datasets[[currentSet()]]$tree)){disable("confounding_start")}
        
        #update silder for binarization cutoff dynamically based on normalized dataset
        min_value <- min(otu)
        max_value <- round(max(otu)/16)
        updateNumericInput(session,"binCutoff",min=min_value,max=max_value)
        
        updateNumericInput(session,"k_in",value=0,min=0,max=vals$datasets[[currentSet()]]$vis_out$K,step=1)
        
        ########updates based on meta info########
        covariates <- vals$datasets[[currentSet()]]$vis_out$covariates
        updateSelectInput(session,"choose",choices = covariates)
        
        #pick all column names, except the SampleID
        group_columns <- setdiff(colnames(meta),sample_column)
        updateSelectInput(session,"alphaGroup",choices = c("-",group_columns))
        updateSelectInput(session,"betaGroup",choices = group_columns)
        updateSelectInput(session,"structureGroup",choices = group_columns)
        updateSelectInput(session,"formula",choices = group_columns)
        updateSelectInput(session,"taxSample",choices=c("NULL",group_columns))
        updateSliderInput(session, "screePCshow", min=1, max=nsamples(phylo), step=1, value=20)
        
        #pick all categorical variables in meta dataframe (except SampleID)
        categorical_vars <- colnames(meta[,unlist(lapply(meta,is.character))])
        categorical_vars <- setdiff(categorical_vars,sample_column)
        updateSelectInput(session,"forest_variable",choices = group_columns)
        updateSelectInput(session,"heatmapSample",choices = c("SampleID",group_columns))
        updateSelectInput(session, "associations_label", choices=c(categorical_vars))
        
        if(is.null(access(phylo,"phy_tree"))) betaChoices="Bray-Curtis Dissimilarity" else betaChoices=c("Bray-Curtis Dissimilarity","Generalized UniFrac Distance", "Unweighted UniFrac Distance", "Weighted UniFrac Distance", "Variance adjusted weighted UniFrac Distance")
        updateSelectInput(session,"betaMethod",choices=betaChoices)
        
        updateSelectInput(session,"phylo_color",choices= c("-","abundance",group_columns,"Kingdom","Phylum","Class","Order","Family","Genus","Species"))
        updateSelectInput(session,"phylo_shape",choices = c("-",group_columns,"Kingdom","Phylum","Class","Order","Family","Genus","Species"))
        updateSelectInput(session,"phylo_size",choices = c("-","abundance",group_columns))
        updateSliderInput(session,"phylo_prune",min=2,max=ntaxa(phylo),value=round(ntaxa(phylo)/2),step=1)
        
        updateSelectInput(session,"filterColumns",choices = c("NONE",group_columns)) 
      }
    }
  })
  
  #observer for legit variables for confounding analysis & basic network
  # -> only categorical vars with more than 1 level
  observe({
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        #factorize meta data
        meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
        categorical_vars <- colnames(meta[,unlist(lapply(meta,is.character))])
        categorical_vars <- setdiff(categorical_vars,sample_column)
        meta <- meta[,c(categorical_vars)]
        meta[] <- lapply(meta,factor)
        #pick all variables which have 2 or more factors for possible variables for confounding!
        tmp <- names(sapply(meta,nlevels)[sapply(meta,nlevels)>1])
        tmp2 <- names(sapply(meta,nlevels)[sapply(meta,nlevels)==2])
        group_columns_no_single <- setdiff(tmp,sample_column)
        group_columns_only_two <- setdiff(tmp2,sample_column)
        updateSelectInput(session,"confounding_var",choices=group_columns_no_single) 
        updateSelectInput(session,"groupCol",choices = group_columns_no_single)
        updateSelectInput(session,"picrust_test_condition", choices = group_columns_only_two)
      }
    }
  })
  
  #this part needs to be in its own "observe" block
  #-> updates ref choice in section "functional topics"
  #-> also update var2 for basic network
  observe({
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        ref_choices <- unique(sample_data(vals$datasets[[currentSet()]]$phylo)[[input$formula]])
        updateSelectInput(session,"refs",choices=ref_choices)
        
        #do not display var chosen for var1 in var2 selection
        groupVariables <- unique(sample_data(vals$datasets[[currentSet()]]$phylo)[[input$groupCol]])
        remainingGroupVariables <- setdiff(groupVariables,input$groupVar1)
        updateSelectInput(session,"groupVar2",choices = c("all",remainingGroupVariables)) 
      }
    }
  })
  
  # check for update if undersampled columns are to be removed (rarefaction curves)
  observeEvent(input$excludeSamples, {
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_meta){
        if(!is.null(vals$undersampled) && input$excludeSamples == T){
          # remove undersampled columns from data
          sampledOTUData <- vals$datasets[[currentSet()]]$normalizedData[,!(colnames(vals$datasets[[currentSet()]]$normalizedData)%in%vals$undersampled)]
          sampledMetaData <- vals$datasets[[currentSet()]]$metaData[!(rownames(vals$datasets[[currentSet()]]$metaData)%in%vals$undersampled),]
          
          #save old (oversampled) data, in case the switch is turned OFF again
          vals$datasets[[currentSet()]]$old.normalizedData <- vals$datasets[[currentSet()]]$normalizedData
          vals$datasets[[currentSet()]]$old.metaData <- vals$datasets[[currentSet()]]$metaData
          old.phylo <- vals$datasets[[currentSet()]]$phylo
          vals$datasets[[currentSet()]]$old.phylo <- old.phylo
          
          #replace old (oversampled) data with new data
          vals$datasets[[currentSet()]]$normalizedData <- sampledOTUData
          vals$datasets[[currentSet()]]$metaData <- sampledMetaData
          
          #build new phyloseq object (with old tree & tax)
          py.otu <- otu_table(sampledOTUData,T)
          py.meta <- sample_data(sampledMetaData)
          #old.tree <- phy_tree(old.phylo)
          if(!is.null(access(old.phylo,"phy_tree"))) old.tree <- phy_tree(old.phylo) else old.tree <- NULL
          old.taxa <- tax_table(old.phylo)
          vals$datasets[[currentSet()]]$phylo <- merge_phyloseq(py.otu, py.meta, old.tree, old.taxa)
          
          #set global Set variable to TRUE, indicating, that undersampled data is already removed
          vals$datasets[[currentSet()]]$undersampled_removed <- T
          
        }else if (input$excludeSamples == F && vals$datasets[[currentSet()]]$undersampled_removed == T){
          #case: undersampled data was removed but shall be used again (switch turned OFF)
          #use old (oversampled) data again, which was saved 
          vals$datasets[[currentSet()]]$normalizedData <- vals$datasets[[currentSet()]]$old.normalizedData
          vals$datasets[[currentSet()]]$metaData <- vals$datasets[[currentSet()]]$old.metaData
          vals$datasets[[currentSet()]]$phylo <- vals$datasets[[currentSet()]]$old.phylo
          vals$datasets[[currentSet()]]$undersampled_removed = F
        } 
      }
    }
  })

  # choose current dataset; return NULL if no set is yet uploaded
  currentSet <- eventReactive(input$datasets_rows_selected, {
    if(length(vals$datasets) == 0){
      return (NULL)
    }
    return(input$datasets_rows_selected)
  })
  
  # update datatable holding currently loaded datasets
  output$datasets <- renderDataTable({
    if(length(vals$datasets)!=0){
      datatable(data.frame(Datasets=names(vals$datasets)),rownames=F,options=list(pageLength=10,dom="t"),selection=list(mode="single",selected=length(vals$datasets))) %>%
        formatStyle("Datasets",color="white",backgroundColor="#222D33")
    } else{
      datatable(data.frame(Datasets=""),rownames=F,options=list(pageLength=10,dom="t"),selection=list(mode="single")) %>%
        formatStyle("Datasets",color="white",backgroundColor="#222D33")
    }
  })
  #####################################
  #    otu upload                     #
  #####################################
  source(file.path("server", "upload_otu_server.R"),  local = TRUE)$value
  #####################################
  #    fastq upload                   #
  #####################################
  source(file.path("server", "upload_fastq_server.R"),  local = TRUE)$value
  #####################################
  #    sample upload                  #
  #####################################
  source(file.path("server", "upload_sample_server.R"),  local = TRUE)$value
  #####################################
  #    data filtering                 #
  #####################################
  source(file.path("server", "filtering_server.R"),  local = TRUE)$value
  #####################################
  #    Basic Analysis                 #
  #####################################
  source(file.path("server", "basic_server.R"),  local = TRUE)$value
  #####################################
  #    Advanced analysis              #
  #####################################
  source(file.path("server", "advanced_server.R"),  local = TRUE)$value
  #####################################
  #    Network analysis               #
  #####################################
  source(file.path("server", "network_server.R"),  local = TRUE)$value
  #####################################
  #     DADA2                         #
  #####################################
  source(file.path("server", "dada2_server.R"),  local = TRUE)$value
  #####################################
  #    Text fields                    #
  #####################################
  source(file.path("server", "texts_server.R"),  local = TRUE)$value
}