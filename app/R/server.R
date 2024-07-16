namco_packages <- c(
  "ade4", "data.table", "cluster", "DT", "fpc", "GUniFrac",
  "heatmaply", "networkD3", "klaR", "phangorn", "plotly",
  "RColorBrewer", "reshape2", "shiny", "textshape",
  "tidyr", "themetagenomics", "igraph", "grid", "dplyr",
  "Matrix", "phyloseq", "NbClust", "caret", "ranger", "gbm",
  "shinyjs", "MLeval", "Rcpp", "MLmetrics", "biomformat",
  "waiter", "dada2", "Biostrings", "fontawesome", "shinyWidgets",
  "shinydashboard", "shinydashboardPlus", "proxy", "parallel",
  "DECIPHER", "SpiecEasi", "ALDEx2", "ggrepel", "SIAMCAT", "gridExtra",
  "genefilter", "fastqcr", "NetCoMi", "metagMisc", "ggnewscale", "ggtree",
  "scales", "ggpubr", "ggsci", "Hmisc", "corrplot", "factoextra",
  "vegan", "decontam", "renv", "shinyBS", "R.utils", "MOFA2",
  "BiocVersion", "biomehorizon", "mixOmics", "ComplexHeatmap"
)

# renv::snapshot(packages= namco_packages, lockfile="app/renv.lock")

check_pkgs <- suppressMessages(lapply(namco_packages, require, character.only = T, quietly = T, warn.conflicts = F))
pkg_status <- ifelse(all(unlist(check_pkgs)), "Dependencies loaded successfully", "Not all dependencies could be loaded")
message(Sys.time(), " - ", pkg_status)
overlay_color <- "rgb(51, 62, 72, .5)"
tree_logo <- fa("tree", fill = "red") # indication logo where phylo-tree is needed
namco_version <- 'v1.1'
captured_messages <- character()  # collected console logs and prints

# Overload cat() to capture messages
cat <- function(...) {
  message_text <- paste(..., sep = " ")
  captured_messages <<- c(captured_messages, message_text)
  base::cat(...)  # Call the original cat function
}

# Overload print() to capture messages
print <- function(...) {
  message_text <- paste(..., sep = " ")
  captured_messages <<- c(captured_messages, message_text)
  base::print(...)  # Call the original print function
}

# Overload message() to capture messages
message <- function(...) {
  message_text <- paste(..., sep = " ")
  captured_messages <<- c(captured_messages, message_text)
  base::message(...)  # Call the original message function
}

server <- function(input, output, session) {
  waiter_hide()
  options(shiny.maxRequestSize = 10000 * 1024^2, stringsAsFactors = F) # upload up to 10GB of data
  source("algorithms.R")
  source("utils.R")
  source("texts.R")
  source("fastq_utils.R")
  message(log_startText)
  
  vals <- reactiveValues(datasets = list(), undersampled = c()) # reactiveValues is a container for variables that might change during runtime and that influence one or more outputs, e.g. the currently selected dataset
  currentSet <- NULL # a pointer to the currently selected dataset
  ncores <- 4 # number of cores used where it is possible to use multiple
  seed <- 123 # Global variable to use as seed
  sessionID <- paste(sample(1:9, 14, replace = T), collapse = "")
  session$onSessionEnded(stopApp) # automatically stop app, if browser window is closed
  sample_column <- "SampleID" # the column with the sample IDs will be renamed to this
  # session ID
  message(paste0("#############", " NAMCO_ID: ", sessionID, " #############"))
  output$sessionIdDiv <- renderText({
    paste0("Session ID: ", sessionID)
  })
  if (!interactive()) {
    sink(stdout(), type = "output")
  } # this makes it so that print statements and other stdOut are saved in log file
  
  
  # choose current dataset; return NULL if no set is yet uploaded
  currentSet <- eventReactive(input$datasets_rows_selected, {
    if (length(vals$datasets) == 0) {
      return(NULL)
    }
    return(input$datasets_rows_selected)
  })
  
  # display logs from all messages
  observeEvent(input$info, {
    output$consoleLogs <- renderText({
      paste(captured_messages, collapse = "\n")
    })
  })
  
  debugging <- F
  if(debugging){
    fastqc.path <- "/usr/bin/fastqc"
    namco_conda_env <- '/usr/local/bin/anaconda3/condabin/conda run -n namco_env'
    lotus2 <- paste0(namco_conda_env, ' lotus2')
    # do not build tree for dada2
  }else{
    fastqc.path <- "/opt/FastQC/fastqc"
    namco_conda_env <- '/opt/miniconda3/bin/conda run -n namco_env'
    lotus2 <- paste0(namco_conda_env, ' lotus2')
  }
  
  #####################################
  #    save & restore session         #
  #####################################
  # Download current session
  output$saveSession <- downloadHandler(
    filename = function() {
      paste("namco_session.RData")
    },
    content = function(file) {
      message(paste0(Sys.time(), " - Saving session ..."))
      session_lst <- vals$datasets[[currentSet()]]
      save(session_lst, file = file)
    }
  )
  
  # upload RData object of namco session
  uploadSessionModal <- function(failed = F, error_message = NULL) {
    modalDialog(
      title = "Restore previous namco session",
      h4("Select session-file:"),
      fluidRow(
        column(6, wellPanel(fileInput("sessionFile", "Select file"), style = "background:#3c8dbc"))
      ),
      if (failed) {
        div(tags$b(error_message, style = "color:red;"))
      },
      footer = tagList(
        modalButton("Cancel", icon = icon("times-circle")),
        actionButton("upload_session_ok", "OK", style = "background-color:blue; color:white")
      ),
      easyClose = T, fade = T, size = "l"
    )
  }
  
  observeEvent(input$loadSession, {
    showModal(uploadSessionModal())
  })
  
  observeEvent(input$upload_session_ok, {
    message(paste0(Sys.time(), " - Restoring previous session ..."))
    info_text <- NULL
    tryCatch(
      {
        base::load(input$sessionFile$datapath)
        session_name <- session_lst[["session_name"]]
        if (session_name %in% names(vals$datasets)) {
          stop(duplicateSessionNameError, call. = F)
        }
        vals$datasets[[session_name]] <- session_lst
        
        #warning if old version is uploaded
        if(vals$datasets[[session_name]]$namco_version != namco_version){
          info_text <- 'The uploaded session file is from an older version of namco. Please be aware that some features might have changed in the meantime.'
        }
        
        # a newly uploaded session is always "not" filtered
        vals$datasets[[session_name]]$filtered <- F
        
        # replace NA and NaN entries in OTU table
        phylo <- vals$datasets[[session_name]]$phylo
        otu_check <- data.frame(otu_table(phylo), check.names = F)
        otu_checked <- replace_NA_OTU(otu_check)
        otu_table(phylo) <- otu_table(otu_checked, taxa_are_rows = T)
        vals$datasets[[session_name]]$phylo <- phylo
        
        # calculate alpha diversity if not present
        if(is.null(vals$datasets[[session_name]]$alpha_diversity)){
          if(vals$datasets[[session_name]]$has_meta){
            alphaTabFull <- createAlphaTab(data.frame(phylo@otu_table, check.names=F), data.frame(phylo@sam_data, check.names = F))
          }else{
            alphaTabFull <- createAlphaTab(data.frame(phylo@otu_table, check.names=F))
          }
          vals$datasets[[session_name]]$alpha_diversity <- alphaTabFull
        }
        
        # store phylo.raw --> will not be subject to filtering
        phylo.raw <- merge_phyloseq(phylo@sam_data, phylo@tax_table, phylo@phy_tree, otu_table(vals$datasets[[session_name]]$rawData, T))
        vals$datasets[[session_name]]$phylo.raw <- phylo.raw
        
        if(vals$datasets[[session_name]]$has_picrust && is.null(vals$datasets[[session_name]]$picrust_results_list)){
          info_text<-"You are uploading an old namco_session file. Some features (picrust analysis) have been reworked in the meantime and you will have to re-calculate these results."
        }
        
        vals$datasets[[session_name]]$is_restored <- T
        updateTabItems(session,"sidebar", selected = "overview")
        removeModal()
        if(!is.null(info_text)){
          showModal(infoModal(info_text))
        }
        # check multi-omics flag for older versions
        if(!"has_omics"%in%names(vals$datasets[[session_name]])) {
          vals$datasets[[session_name]]$has_omics <- F
        }
      },
      error = function(e) {
        print(e$message)
        showModal(uploadSessionModal(failed = T, error_message = e$message))
      }
    )
  })
  
  #####################################
  #    menu items & info-box          #
  #####################################
  
  output$fastq_overview <- renderMenu({
    if (!is.null(currentSet())) {
      if (vals$datasets[[currentSet()]]$is_fastq) {
        menuItem("fastq Overview", tabName = "fastq_tab", icon = icon("dna"))
      }
    }
  })
  
  output$filtering_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Data Overview & Filtering", tabName = "overview", icon = icon("filter"))
    }
  })
  
  output$basic_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Basic Analysis", tabName = "basics", icon = icon("search"))
    }
  })
  
  output$differential_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Differential Analysis", tabName = "differential", icon = icon("object-group"))
    }
  })
  
  output$functional_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Functional Analysis", tabName = "functional", icon = icon("wrench"))
    }
  })
  
  output$phylo_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Phylogenetic Analysis", tabName = "phylogenetic", icon = icon("tree"))
    }
  })
  
  output$network_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Network Analysis", tabName = "network", icon = icon("project-diagram"))
    }
  })
  
  output$ml_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Machine Learning", tabName = "machineLearning", icon = icon("brain"))
    }
  })
  
  output$confounding_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Confounding Analysis", tabName = "confounding", icon = icon("bolt"))
    }
  })
  
  output$multiomics_menu <- renderMenu({
    if(!is.null(currentSet())){
      menuItem("Multi-omics Analysis", tabName = "multiomics", icon = icon("microscope"))
    }
  })
  
  output$samples_box1 <- renderValueBox({
    samples <- 0
    if (!is.null(currentSet())) {
      samples <- length(sample_names(vals$datasets[[currentSet()]]$phylo))
    }
    valueBox(samples, "Samples", icon = icon("list"), color = "blue")
  })
  
  output$conditions_box1 <- renderValueBox({
    groups <- 0
    if (!is.null(currentSet())) {
      if (vals$datasets[[currentSet()]]$has_meta) {
        groups <- dim(vals$datasets[[currentSet()]]$metaData)[2]
      }
    }
    valueBox(groups, "Groups", icon = icon("list"), color = "purple")
  })
  
  output$otus_box1 <- renderValueBox({
    groups <- 0
    if (!is.null(currentSet())) {
      otus <- ntaxa(vals$datasets[[currentSet()]]$phylo)
    }
    valueBox(otus, "ASVs/OTUs", icon = icon("list"), color = "orange")
  })
  
  output$samples_box2 <- renderValueBox({
    samples <- 0
    if (!is.null(currentSet())) {
      samples <- length(sample_names(vals$datasets[[currentSet()]]$phylo))
    }
    valueBox(samples, "Samples", icon = icon("list"), color = "blue")
  })
  
  output$conditions_box2 <- renderValueBox({
    groups <- 0
    if (!is.null(currentSet())) {
      if (vals$datasets[[currentSet()]]$has_meta) {
        groups <- dim(vals$datasets[[currentSet()]]$metaData)[2]
      }
    }
    valueBox(groups, "Groups", icon = icon("list"), color = "purple")
  })
  
  output$otus_box2 <- renderValueBox({
    groups <- 0
    if (!is.null(currentSet())) {
      otus <- ntaxa(vals$datasets[[currentSet()]]$phylo)
    }
    valueBox(otus, "ASVs", icon = icon("list"), color = "orange")
  })
  
  output$samples_box3 <- renderValueBox({
    samples <- 0
    if (!is.null(currentSet())) {
      samples <- length(sample_names(vals$datasets[[currentSet()]]$phylo))
    }
    valueBox(subtitle = "Samples", value = samples, icon = icon("list"), color = "blue", width = 6)
  })
  
  output$statTestSignifCount <- renderValueBox({
    count <- 0
    if(!is.null(currentSet())){
      if(!is.null(statTestReactive())){
        count <- length(statTestReactive())
      }
    }
    valueBox(subtitle = "Significant", value=count, icon=icon("star"), color="green")
  })
  
  
  ############################################
  #    observers & datasets & normalization  #
  ############################################
  
  # observer for normalization
  observeEvent(input$normalizationApply, {
    if (!is.null(currentSet())) {
      normMethod <- which(input$normalizationSelect == c("no Normalization", "by minimum Sampling Depth", "by Rarefaction", "centered log-ratio", "Total Sum Normalization (normalize to 10,000 reads)", "Spike-in Normalization")) - 1
      normalized_dat <- normalizeOTUTable(vals$datasets[[currentSet()]]$phylo, normMethod)
      vals$datasets[[currentSet()]]$normalizedData <- normalized_dat$norm_tab
      otu_table(vals$datasets[[currentSet()]]$phylo) <- otu_table(normalized_dat$norm_tab, T)
      vals$datasets[[currentSet()]]$normMethod <- normMethod
      vals$datasets[[currentSet()]]$normFunction <- normalized_dat$f
    }
  })
  
  # normalization info
  output$normalizationDropdown <- renderMenu({
    if(!is.null(currentSet())){
      methods <- c("no Normalization", "by minimum Sampling Depth", "by Rarefaction", "centered log-ratio", "Total Sum Normalization (normalize to 10,000 reads)", "Spike-in Normalization", "Copy Number (Picrust2)")
      normMethod <-  vals$datasets[[currentSet()]]$normMethod
      normMethodText <- methods[normMethod+1]
      
      f <- vals$datasets[[currentSet()]]$normFunction
      if (!is.null(f) && startsWith(f, "Error")) {
        dropdownMenu(type="notification", badgeStatus = "warning",
                     notificationItem(
                       text=f,
                       icon("exclamation-triangle"),
                       status="warning"
                     ))
      } else if (normMethod==0) {
        dropdownMenu(type="notification", badgeStatus = "warning",
                     notificationItem(
                       text="Currently no normalization applied!",
                       icon("exclamation-triangle"),
                       status="warning"
                     ))
      } else {
        dropdownMenu(type="notification", badgeStatus = "success",
                     notificationItem(
                       text=tags$div("Currently used normalization: ",
                                     tags$br(),
                                     normMethodText,
                                     style = "display: inline-block; vertical-align: middle;"),
                       icon("info"),
                       status="success"
                     ))
      }
    }
  })
  
  
  # observer to show only tabs with / without meta file
  observe({
    if (!is.null(currentSet())) {
      if (!vals$datasets[[currentSet()]]$has_meta) {
        hideTab(inputId = "filters", target = "Data Overview")
        hideTab(inputId = "basicPlots", target = "Beta Diversity")
        hideTab(inputId = "differentialPlots", target = "Associations")
        hideTab(inputId = "differentialPlots", target = "Correlations")
        hideTab(inputId = "differentialPlots", target = "Topic Modeling")
        hideTab(inputId = "differentialPlots", target = "Time-series analysis")
        hideTab(inputId = "differentialPlots", target = "Differential statistical Analysis")
        hideTab(inputId = "netWorkPlots", target = "Co-occurrence of OTUs/ASVs")
        hideTab(inputId = "netWorkPlots", target = "Network inference")
        hideTab(inputId = "netWorkPlots", target = "Taxonomic Rank Networks")
        hideTab(inputId = "netWorkPlots", target = "Differential Networks")
        hideTab(inputId = "machineLearning", target ="Random Forests")
        hideTab(inputId = "confoundingPlots", target = "Confounding Analysis & Explained Variation")
      } else {
        showTab(inputId = "filters", target = "Data Overview")
        showTab(inputId = "basicPlots", target = "Beta Diversity")
        showTab(inputId = "differentialPlots", target = "Associations")
        showTab(inputId = "differentialPlots", target = "Correlations")
        showTab(inputId = "differentialPlots", target = "Topic Modeling")
        showTab(inputId = "differentialPlots", target = "Time-series analysis")
        showTab(inputId = "differentialPlots", target = "Differential statistical Analysis")
        showTab(inputId = "netWorkPlots", target = "Co-occurrence of OTUs/ASVs")
        showTab(inputId = "netWorkPlots", target = "Network inference")
        showTab(inputId = "netWorkPlots", target = "Taxonomic Rank Networks")
        showTab(inputId = "netWorkPlots", target = "Differential Networks")
        showTab(inputId = "machineLearning", target ="Random Forests")
        showTab(inputId = "confoundingPlots", target = "Confounding Analysis & Explained Variation")
      }
    }
  })
  
  # observer for MSD 
  observe({
    if(!is.null(input$msdFile)){
      shinyjs::show("msdLinkSelect")
      msd_file <- read.table(input$msdFile$datapath, sep="\t", header=T)
      updateSelectInput(session, "msdLinkSelect", choices=msd_file$Link)
    }else{
      shinyjs::hide("msdLinkSelect")
    }
  })
  
  # observer for fastq-related stuff
  observe({
    if(input$clustering_lotus == 'dada2'){
      shinyjs::show('dada2_lotus2_warning', anim = T)
    }else{
      shinyjs::hide('dada2_lotus2_warning', anim = T)
    }
    if (!is.null(currentSet())) {
      if (vals$datasets[[currentSet()]]$is_fastq) {
        updateSelectInput(session, "fastq_file_select_pre", choices = vals$datasets[[currentSet()]]$generated_files$file_df$sample_names)
        updateSelectInput(session, "fastq_file_select_post", choices = vals$datasets[[currentSet()]]$generated_files$file_df$sample_names)
      }
    }
  })
  
  # observer for picrust2
  observe({
    if(!is.null(currentSet())){
      if(!is.null(vals$datasets[[currentSet()]]$picrust_analysis_list)){
        results <- vals$datasets[[currentSet()]]$picrust_analysis_list
        ec <- results$test_EC[order(results$test_EC$pval1),]
        ec_picks <- rownames(ec)
        names(ec_picks) <- paste0(rownames(ec)," (",round(ec$pval1, 5),")")
        updatePickerInput(session,"picrust_ec_select", choices=ec_picks, options = list(`liveSearch` = T))
        
        ko <- results$test_KO[order(results$test_KO$pval1),]
        ko_picks <- rownames(ko)
        names(ko_picks) <- paste0(rownames(ko)," (",round(ko$pval1, 5),")")
        updatePickerInput(session,"picrust_ko_select", choices=ko_picks, options = list(`liveSearch` = T))
        
        pw <- results$test_PW[order(results$test_PW$pval1),]
        pw_picks <- rownames(pw)
        names(pw_picks) <- paste0(rownames(pw)," (",round(pw$pval1, 5),")")
        updatePickerInput(session,"picrust_pw_select", choices=pw_picks, options = list(`liveSearch` = T))
      }
    }
  }, priority = 3)
  
  # filter variables
  observe({
    if (!is.null(currentSet())) {
      phylo <- vals$datasets[[currentSet()]]$phylo
      taxonomy <- data.frame(tax_table(phylo))
      
      filterTaxaValues <- unique(taxonomy[[input$filterTaxa]])
      updatePickerInput(session, "filterTaxaValues", choices = filterTaxaValues)
      
      if (vals$datasets[[currentSet()]]$has_meta) {
        meta <- data.frame(sample_data(phylo))
        
        shinyjs::show("taxBinningDiv")
        
        # filter variables
        filterColumnValues <- unique(meta[[input$filterColumns]])
        updateSelectInput(session, "filterColumnValues", choices = filterColumnValues)
        
      } else {
        shinyjs::hide("taxBinningDiv")
      }
    }
  }, priority = 3)
  
  # observer for inputs depending on choosing a meta-group first
  observe({
    if (!is.null(currentSet())) {
      phylo <- vals$datasets[[currentSet()]]$phylo
      methods <- c("no Normalization", "by minimum Sampling Depth", "by Rarefaction", "centered log-ratio", "Total Sum Normalization (normalize to 10,000 reads)")
      
      if (vals$datasets[[currentSet()]]$has_meta) {
        meta <- data.frame(sample_data(phylo), check.names = F)
        
        # display sample names which can be filtered
        if (input$filterColumns == "NONE") {
          samples_left <- sample_names(phylo)
        } else if (input$filterColumns != "" && !is.null(input$filterColumnValues)) {
          samples_left <- meta[which(unlist(meta[eval(input$filterColumns)]) %in% input$filterColumnValues),][[sample_column]]
        } else {
          samples_left <- NULL
        }
        updatePickerInput(session, "filterSample", choices = samples_left)
        
        # tax binning
        groupVariables <- unique(meta[[input$taxBinningYLabel]])
        updateSelectizeInput(session, 'taxBinningYOrder', choices = c(groupVariables), server=T,selected = sort(groupVariables)[1:min(100, length(groupVariables))])
        
        # associations
        caseVariables <- unique(meta[[input$associations_label]])
        updateSelectizeInput(session, "associations_case", choices = c(caseVariables), server=T)
        
        #decontamination
        groupVariables <- unique(meta[[input$controlSamplesColumn]])
        updateSelectizeInput(session, 'controlSamplesName', choices = c(groupVariables), server=T)
        
        # basic network variables
        groupVariables <- unique(meta[[input$groupCol]])
        updateSelectizeInput(session, "groupVar1", choices = c(groupVariables), server=T)
        
        # beta diversity
        groupVariables <- unique(meta[[input$betaGroup]])
        updateSelectizeInput(session, "betaLevel", choices = c("All", groupVariables), server=T)
        
        
        # alpha diversity -> select pairs for wilcoxon test
        if (input$alphaGroup != "-") {
          lev <- levels(as.factor(meta[[input$alphaGroup]]))
          if(length(lev) > 1){
            pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
            pairs_list <- sapply(pairs, paste, collapse = " vs. ")
            updatePickerInput(session, "alphaPairs", choices = pairs_list)  
          }
        }
        
        # time series
        timePoints <- unique(meta[[input$timeSeriesGroup]])
        updateSelectizeInput(session, "timeSeriesTimePointOrder", choices = c(timePoints), selected = c(timePoints), server=T)
        highlight_values <- unique(meta[[input$timeSeriesBackground]])
        updateSelectizeInput(session, "timeSeriesSampleHighlight", choices=c("NONE",highlight_values), server=T)
        
        #picrust
        groupVariables <- unique(meta[[input$picrust_test_condition]])
        updateSelectizeInput(session, "picrust_test_covariate", choices = c(groupVariables), server=T)
        
        # add spike-in normalization
        methods <- c(methods, "Spike-in Normalization")
      }else{
        updatePickerInput(session, "filterSample", choices = sample_names(phylo))
      }
      if(input$normalizationSelect!="") {
        updateSelectInput(session, "normalizationSelect", selected = input$normalizationSelect, choices = methods)
      } else {
        updateSelectInput(session, "normalizationSelect", choices = methods) 
      }
    }
  }, priority = 3)
  
  # update input selections
  observe({
    if (!is.null(currentSet())) {
      phylo <- vals$datasets[[currentSet()]]$phylo
      otu <- otu_table(phylo)
      taxonomy <- data.frame(tax_table(phylo))
      
      updateSelectInput(session, "structureCompOne", choices = (1:nsamples(phylo)))
      updateSelectInput(session, "structureCompTwo", choices = (2:nsamples(phylo)))
      updateSelectInput(session, "structureCompThree", choices = (3:nsamples(phylo)))
      updateSelectInput(session, "pcaLoading", choices = (1:nsamples(phylo)))
      updateSelectizeInput(session, "horizonTaxaSelect", choices = c("", rownames(taxonomy)), server=T)
      updateSliderInput(session, "pcaLoadingNumber", min = 2, max = ntaxa(phylo) - 1, value = ifelse(ntaxa(phylo) < 10, ntaxa(phylo), 10), step = 2)
      if (!vals$datasets[[currentSet()]]$is_fastq) {
        updateSelectInput(session, "filterTaxa", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
        updateSelectInput(session, "taxBinningLevel", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
        updateSelectInput(session, "horizonTaxaLevel", choices = rev(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
        updateSelectInput(session, "mofa2_taxa_level", choices = rev(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
      }
      updateSelectInput(session, "taxBinningYLabel", choices = c('None', sample_column))
      if(input$taxBinningOrderManually){
        shinyjs::show('taxBinningYOrder')
      }else{
        shinyjs::hide('taxBinningYOrder')
      }
      
      if(!is.null(phylo@phy_tree)){
        updateSelectInput(session, "confounding_distance", choices=c("Unifrac","Bray-Curtis"))
        updateSelectInput(session, 'heatmapDistance', choices=c("bray", "gunifrac", "wunifrac", "unifrac", "jsd"))
      }else{
        updateSelectInput(session, 'heatmapDistance', choices=c("bray", "jsd"))
      }
      
      if (vals$datasets[[currentSet()]]$has_meta) {
        # get tables from phyloseq object
        meta <- data.frame(sample_data(phylo), check.names = F)
        if (!is.null(access(vals$datasets[[currentSet()]]$phylo, "phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
        phylo <- vals$datasets[[currentSet()]]$phylo
        
        updateSliderInput(session, "rareToShow", min = 1, max = ncol(otu), value = min(50, ncol(otu)))
        updateSliderInput(session, "rareToHighlight", min = 1, max = ncol(otu), value = round(ncol(otu) / 10))
        updateSliderInput(session, "top_x_features", min = 1, max = nrow(otu))
        if (ncol(meta) > 2) {
          enable("confounding_start")
        }
        if (is.null(vals$datasets[[currentSet()]]$tree)) {
          disable("confounding_start")
        }
        
        # update silder for binarization cutoff dynamically based on normalized dataset
        min_value <- min(otu)
        max_value <- round(max(otu) / 16)
        updateNumericInput(session, "binCutoff", min = min_value, max = max_value)
        updateNumericInput(session, "k_in", value = 0, min = 0, max = vals$datasets[[currentSet()]]$vis_out$K, step = 1)
        
        ######## updates based on meta info########
        covariates <- vals$datasets[[currentSet()]]$vis_out$covariates
        updateSelectInput(session, "choose", choices = covariates)
        updateSelectInput(session, "taxBinningYLabel", choices = c("--Combined--", colnames(meta)), selected = sample_column)
        updateSelectInput(session, "alphaGroup", choices = c("-",sample_column, colnames(meta)), selected = "-")
        updateSelectInput(session, "horizonSubject", choices = colnames(meta))
        updateSelectInput(session, "horizonSample", choices = colnames(meta))
        updateSelectInput(session, "horizonCollectionDate", choices = colnames(meta))
        updateSelectInput(session, "horizonSubjectSelection", choices = colnames(meta))
        # add meta options to metabolomics tab
        updateSelectInput(session, "mofa2_sample_label", choices=colnames(meta))
        updateSelectInput(session, "mofa2_condition_label", choices=colnames(meta))
        updateSelectInput(session, "mofa2_group_label", choices=colnames(meta))
        
        # pick all column names, except the SampleID
        group_columns <- setdiff(colnames(meta), sample_column)
        updateSelectInput(session, "structureGroup", choices = group_columns)
        updateSelectInput(session, "formula", choices = group_columns)
        updateSelectInput(session, "taxSample", choices = c("NULL", group_columns))
        updateSliderInput(session, "screePCshow", min = 1, max = nsamples(phylo), step = 1, value = 20)
        updateSelectInput(session, "timeSeriesGroup", choices=c(group_columns))
        updateSelectInput(session, "DNAconcentrationColumn", choices=c('NULL',group_columns))
        updateSelectInput(session, "controlSamplesColumn", choices=c('NULL',group_columns))
        
        # pick all categorical variables in meta dataframe (except SampleID) == variables with re-appearing values
        # also do not show columns which have the same value for each entry
        categorical_vars <- names(which(sapply(meta, function(x) {
          length(unique(x))>1 && length(levels(as.factor(na.omit(x)))) < length(na.omit(x))
        })))
        categorical_vars <- setdiff(categorical_vars, sample_column)
        updateSelectInput(session, "forest_variable", choices = group_columns)
        updateSelectInput(session, "heatmapSample", choices = c("SampleID", categorical_vars))
        updateSelectInput(session, "associations_label", choices = c(categorical_vars))
        updateSelectInput(session, "betaGroup", choices = categorical_vars)
        updateSelectInput(session, "betaGroup2", choices = c("None",categorical_vars))
        updateSelectInput(session, "taxBinningGroup", choices = c("None", categorical_vars))
        updateSelectInput(session, "timeSeriesBackground", choices=c(sample_column, categorical_vars))
        updateSelectInput(session, "timeSeriesMeanLine", choices = c("NONE", categorical_vars))
        updateSelectInput(session, "statTestGroup", choices=c(categorical_vars))
        updateSelectInput(session, "picrust_test_condition", choices = c(categorical_vars))
        
        # pick all numerical/continuous variables in dataframe
        numerical_vars <- colnames(meta %>% select_if(is.numeric))
        updatePickerInput(session, "corrSelectGroups", choices = numerical_vars)
        
        if (is.null(access(phylo, "phy_tree"))) betaChoices <- "Bray-Curtis Dissimilarity" else betaChoices <- c("Bray-Curtis Dissimilarity", "Generalized UniFrac Distance", "Unweighted UniFrac Distance", "Weighted UniFrac Distance", "Variance adjusted weighted UniFrac Distance")
        updateSelectInput(session, "betaMethod", choices = betaChoices)
        
        updateSliderInput(session, "phylo_prune", min = 2, max = ntaxa(phylo), value = round(ntaxa(phylo) * 0.45), step = 1)
        updateSelectInput(session, "phylo_group", choices = c("NONE", categorical_vars))
        
        updateSelectInput(session, "filterColumns", choices = c("NONE", group_columns))
      }
    }
  }, priority = 3)
  
  # observer for legit variables for confounding analysis & basic network
  # -> only categorical vars with more than 1 level
  observe({
    if (!is.null(currentSet())) {
      if (vals$datasets[[currentSet()]]$has_meta) {
        # factorize meta data
        meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo), check.names=F)
        categorical_vars <- names(which(sapply(meta, function(x) {
          length(unique(x))>1 && length(levels(as.factor(na.omit(x)))) < length(na.omit(x))
        })))
        categorical_vars <- setdiff(categorical_vars, sample_column)
        meta <- meta[,which(colnames(meta) %in% categorical_vars)]
        meta[] <- lapply(meta, factor)
        
        tmp <- names(sapply(meta, nlevels)[sapply(meta, nlevels) > 1]) # pick all variables which have 2 or more factors for possible variables for confounding!
        tmp2 <- names(sapply(meta, nlevels)[sapply(meta, nlevels) == 2]) # variables with exactly 2 levels
        group_columns_no_single <- setdiff(tmp, sample_column)
        group_columns_only_two <- setdiff(tmp2, sample_column)
        updateSelectInput(session, "confounding_var", choices = group_columns_no_single)
        updateSelectInput(session, "groupCol", choices = group_columns_no_single)
        updateSelectInput(session, "diffNetworkSplitVariable", choices = group_columns_only_two)
        
        # own calculation for confounding analysis
        meta <- vals$datasets[[currentSet()]]$phylo@sam_data
        meta[[sample_column]]<-NULL
        variables <- colnames(meta[,apply(meta, 2, function(x){return(length(unique(x))!=1)})])
        updatePickerInput(session, "confounding_select_tested_variable", choices=variables, selected = variables)
        
      }
    }
  }, priority = 3)
  
  observe({
    if(!is.null(decontamReactive())){
      df <- decontamReactive()$contamdf
      contams <- df[which(df$contaminant),]$feature
      if(input$DNAconcentrationColumn != 'NULL') {
        updateSelectInput(session, 'contamCandidatesSelect', choices = c(contams))
      }else{
        updateSelectInput(session, 'contamCandidatesSelect', choices = c())
      }
      updatePickerInput(session, 'contamCandidatesSelectRemove', choices = c(contams))
    }
  })
  
  # this part needs to be in its own "observe" block
  #-> updates ref choice in section "functional topics"
  #-> also update var2 for basic network
  observe({
    if (!is.null(currentSet())) {
      if (vals$datasets[[currentSet()]]$has_meta) {
        ref_choices <- unique(sample_data(vals$datasets[[currentSet()]]$phylo)[[input$formula]])
        updateSelectizeInput(session, "refs", choices = ref_choices, server = T)
        
        # do not display var chosen for var1 in var2 selection
        groupVariables <- unique(sample_data(vals$datasets[[currentSet()]]$phylo)[[input$groupCol]])
        remainingGroupVariables <- setdiff(groupVariables, input$groupVar1)
        updateSelectizeInput(session, "groupVar2", choices = c("all", remainingGroupVariables), server = T)
      }
    }
  }, priority = 3)
  
  # update datatable holding currently loaded datasets
  output$datasets <- renderDataTable({
    if (length(vals$datasets) != 0) {
      datatable(data.frame(Datasets = names(vals$datasets)), rownames = F, options = list(pageLength = 10, dom = "t"), selection = list(mode = "single", selected = length(vals$datasets))) %>%
        formatStyle("Datasets", color = "white", backgroundColor = "#222D33")
    } else {
      datatable(data.frame(Datasets = ""), rownames = F, options = list(pageLength = 10, dom = "t"), selection = list(mode = "single")) %>%
        formatStyle("Datasets", color = "white", backgroundColor = "#222D33")
    }
  }, priority = 4)
  
  dataset_proxy <- dataTableProxy("datasets")
  #####################################
  #    otu upload                     #
  #####################################
  source(file.path("server", "upload_otu_server.R"), local = TRUE)$value
  #####################################
  #    fastq upload (DADA2)           #
  #####################################
  source(file.path("server", "upload_dada2_server.R"), local = TRUE)$value
  #####################################
  #    fastq upload (LotuS2)          #
  #####################################
  source(file.path("server", "upload_lotus2_server.R"), local = TRUE)$value
  #####################################
  #    msd upload                     #
  #####################################
  source(file.path("server", "upload_msd_server.R"), local = TRUE)$value
  #####################################
  #    sample upload                  #
  #####################################
  source(file.path("server", "upload_sample_server.R"), local = TRUE)$value
  #####################################
  #    omics upload                   #
  #####################################
  source(file.path("server", "upload_multiomics_server.R"), local = TRUE)$value
  #####################################
  #    data filtering                 #
  #####################################
  source(file.path("server", "filtering_server.R"), local = TRUE)$value
  #####################################
  #    Basic Analysis                 #
  #####################################
  source(file.path("server", "basic_server.R"), local = TRUE)$value
  #####################################
  #    Differential Analysis          #
  #####################################
  source(file.path("server", "differential_server.R"), local = TRUE)$value
  #####################################
  #    Functional Analysis            #
  #####################################
  source(file.path("server", "functional_server.R"), local = TRUE)$value
  #####################################
  #    Phylogenetic Analysis          #
  #####################################
  source(file.path("server", "phylogenetic_server.R"), local = TRUE)$value
  #####################################
  #    Network analysis               #
  #####################################
  source(file.path("server", "network_server.R"), local = TRUE)$value
  #####################################
  #    Confounding Analysis           #
  #####################################
  source(file.path("server", "confounding_server.R"), local = TRUE)$value
  #####################################
  #    Multi-omics Analysis           #
  #####################################
  source(file.path("server", "multiomics_server.R"), local = TRUE)$value
  #####################################
  #     fastq related things          #
  #####################################
  source(file.path("server", "fastq_server.R"), local = TRUE)$value
  #####################################
  #    Text fields                    #
  #####################################
  source(file.path("server", "texts_server.R"), local = TRUE)$value
}

