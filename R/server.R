packages <- c("ade4", "cluster", "data.table", "DT", "fpc", "GUniFrac", "heatmaply", "networkD3",
              "klaR", "phangorn", "plotly", "RColorBrewer", "reshape2", "Rtsne", "shiny", "textshape",
              "tidyr", "umap", "themetagenomics", "SpiecEasi", "igraph", "Matrix", "phyloseq", "NbClust", 
              "caret", "ranger", "gbm", "shinyjs", "MLeval", "Rcpp", "MLmetrics", "mdine")
suppressMessages(lapply(packages, require, character.only=T, quietly=T, warn.conflicts=F))

server <- function(input,output,session){
  options(shiny.maxRequestSize=1000*1024^2,stringsAsFactors=F)  # upload up to 1GB of data
  source("algorithms.R")
  source("utils.R")
  source("texts.R")
  
  vals = reactiveValues(datasets=list(),undersampled=c()) # reactiveValues is a container for variables that might change during runtime and that influence one or more outputs, e.g. the currently selected dataset
  currentSet = NULL # a pointer to the currently selected dataset
  ncores = 4  # number of cores used where it is possible to use multiple
  seed = 123 # Global variable to use as seed
  session$onSessionEnded(stopApp) #automatically stop app, if browser window is closed
  
  #####################################
  #    file upload                    #
  #####################################
  
  
  # Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
  uploadModal <- function(failed=F,error_message=NULL) {
    modalDialog(
      h4("Please provide pregenerated input files. For detailed information how the files have to look, check out the Info & Settings tab on the left!"),
      hr(),
      fluidRow(
        column(6,fileInput("otuFile","Select OTU table")),
        column(6,fileInput("metaFile","Select Metadata File"))
      ),
      fluidRow(
        column(6,checkboxInput("taxInOTU","Click here if the taxonomic classification is stored in a seperate file:",F))
      ),
      fixedRow(
        column(6,fileInput("taxFile","Select Taxonomic classification file")),
        column(6,fileInput("treeFile","Select Phylogenetic Tree File (optional)",width="100%"))
      ),
      hr(),
      fluidRow(
        column(10,radioButtons("normMethod","Normalization Method",c("no Normalization","by Sampling Depth","by Rarefaction"),inline=T))
      ),
      br(),
      textInput("dataName","Enter a project name:",placeholder="New_Project",value="New_Project"),
      if(failed) {
        #div(tags$b("The file you specified could not be loaded. Please check the Info tab and to confirm your data is in the correct format!",style="color: red;"),
            tags$p(error_message,style="color:red;")
      },
      footer = tagList(
        modalButton("Cancel"),
        actionButton("upload_ok","OK",style="background-color:blue; color:white")
      )
    )
  }
  #observer for taxInOTU checkbox
  observeEvent(input$taxInOTU,{
    if(!input$taxInOTU){shinyjs::hide("taxFile")} else {shinyjs::show("taxFile")}
  })
  
  uploadTestdataModal <- function(failed=F, error_message=NULL){
    modalDialog(
      h4("Choose a sample-Dataset (for details of each set, look into Info & Settings tab!)"),
      fluidRow(column(1),
               column(10, selectInput("selectTestdata", shiny::HTML("<p><span style='color: green'>Select Sample-Dataset</span></p>"),choices = c("Mueller et al. (Mice samples)"))),
               column(1)
      ),
      fluidRow(
        column(10,radioButtons("normMethod","Normalization Method",c("no Normalization","by Sampling Depth","by Rarefaction"),inline=T))
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("upload_testdata_ok","OK",style="background-color:blue; color:white")
      )
    )
  }
  
  errorModal <- function(error_message=NULL){
    modalDialog(
      p(error_message,style="color:red;"),
      easyClose = T,
      modalButton("Cancel")
    )
  }
  
  # launch upload dialog
  observeEvent(input$upload, {
    showModal(uploadModal())
  })
  
  # try to load the dataset specified in the dialog window
  observeEvent(input$upload_ok, {
    tryCatch({
      #withProgress()
      if(input$dataName%in%names(vals$datasets)){stop(duplicateSessionNameError,call. = F)}
      if(is.null(input$otuFile) || is.null(input$metaFile)){stop(otuOrMetaMissingError,call. = F)}
      if(!file.exists(input$otuFile$datapath)){stop(otuFileNotFoundError,call. = F)}
      if(!file.exists(input$metaFile$datapath)){stop(metaFileNotFoundError,call. = F)}
      
      ##read OTU (and taxa) file ##
      
      #case: taxonomy in otu-file -> no taxa file provided
      if(is.null(input$taxFile)){
        otu <- read.csv(input$otuFile$datapath,header=T,sep="\t",row.names=1,check.names=F) #load otu table -> rows are OTU, columns are samples
        #stop if no taxonomy column present
        if (! "taxonomy" %in% colnames(otu)){stop(noTaxaInOtuError,call. = F)}
        taxonomy = generateTaxonomyTable(otu) # generate taxonomy table from TAX column
        otu = otu[!apply(is.na(otu)|otu=="",1,all),-ncol(otu)] # remove "empty" rows
        otus <- row.names(otu) #save OTU names
        otu <- sapply(otu,as.numeric) #make OTU table numeric
        rownames(otu) <- otus
      }
      #case: taxonomy in seperate file -> taxonomic classification file
      else{
        otu <- read.csv(input$otuFile$datapath,header=T,sep="\t",row.names=1,check.names=F) #load otu table -> rows are OTU, columns are samples
        otus <- row.names(otu) #save OTU names
        taxonomy <- read.csv(input$taxFile$datapath,header=T,sep="\t",row.names=1,check.names=F) #load taxa file
        #check for consistent OTU naming in OTU and taxa file:
        if(!identical(rownames(otu),rownames(taxonomy))){stop(otuNoMatchTaxaError,call. = F)}
        taxonomy[is.na(taxonomy)] <- "NA"
        otu <- sapply(otu,as.numeric) #make OTU table numeric
        rownames(otu) <- otus
      }

      ## read meta file ##
      meta <- read.csv(input$metaFile$datapath,header=T,sep="\t")
      rownames(meta)=meta[,1]
      if(!setequal(colnames(otu),meta$SampleID)){stop(unequalSamplesError,call. = F)}
      meta = meta[match(colnames(otu),meta$SampleID),]
      #set SampleID column to be character, not numeric (in case the sample names are only numbers)
      meta$SampleID <- as.character(meta$SampleID)
      
      ## read phylo-tree file ##
      if(is.null(input$treeFile)) tree = NULL else {
        if(!file.exists(input$treeFile$datapath)){stop(treeFileNotFoundError,call. = F)}
        tree = read.tree(input$treeFile$datapath)
      }
      
      normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
      normalized_dat = normalizeOTUTable(otu, normMethod)
      #tax_binning = taxBinning(normalized_dat[[2]],taxonomy)
      #create phyloseq object from data (OTU, meta, taxonomic, tree)
      py.otu <- otu_table(normalized_dat$norm_tab,T)
      py.tax <- tax_table(as.matrix(taxonomy))
      py.meta <- sample_data(meta)
      
      #cannot build phyloseq object with NULL as tree input; have to check both cases:
      if (!is.null(tree)) phylo <- merge_phyloseq(py.otu,py.tax,py.meta, tree) else phylo <- merge_phyloseq(py.otu,py.tax,py.meta)
      
      #pre-build unifrac distance matrix
      if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalized_dat$norm_tab,meta,tree) else unifrac_dist <- NULL
      
      vals$datasets[[input$dataName]] <- list(rawData=otu,metaData=meta,taxonomy=taxonomy,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab,tree=tree,phylo=phylo,unifrac_dist=unifrac_dist,undersampled_removed=F, filtered=F, normMethod = normMethod)
      updateTabItems(session,"sidebar")
      removeModal()
      
    },error=function(e){
      print(e)
      showModal(uploadModal(failed=T,error_message = e))
    })

  })
  
  ##### sample datasets: ####
  
  #load specified testdataset
  observeEvent(input$upload_testdata, {
    showModal(uploadTestdataModal())
  })
  
  #load different testdata sets
  observeEvent(input$upload_testdata_ok, {
    if(input$selectTestdata == "Mueller et al. (Mice samples)"){
      dat <- read.csv("testdata/OTU_table.tab",header=T,sep="\t",row.names=1,check.names = F) # load otu table -> samples are columns
      taxonomy = generateTaxonomyTable(dat) # generate taxonomy table from TAX column
      dat = dat[!apply(is.na(dat)|dat=="",1,all),-ncol(dat)] # remove "empty" rows
      meta = read.csv("testdata/metafile.tab",header=T,sep="\t")
      rownames(meta) = meta[,1]
      meta = meta[match(colnames(dat),meta$SampleID),]
      tree = read.tree("testdata/tree.tre") # load phylogenetic tree
      
      normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
      normalized_dat = normalizeOTUTable(dat,normMethod)
      #tax_binning = taxBinning(normalized_dat[[2]],taxonomy)
      
      #create phyloseq object from data (OTU, meta, taxonomic, tree)
      py.otu <- otu_table(normalized_dat$norm_tab,T)
      py.tax <- tax_table(as.matrix(taxonomy))
      py.meta <- sample_data(meta)
      phylo <- merge_phyloseq(py.otu,py.tax,py.meta,tree)
      
      #pre-build unifrac distance matrix
      if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalized_dat$norm_tab,meta,tree) else unifrac_dist <- NULL
      
      #the final dataset
      dataset<- list(rawData=dat,metaData=meta,taxonomy=taxonomy,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab,tree=tree,phylo=phylo,unifrac_dist=unifrac_dist,undersampled_removed=F,filtered=F, normMethod = normMethod)
      
      vals$datasets[["Mueller et al."]] <- dataset
      updateTabItems(session,"sidebar")
      removeModal()
    }
    if(input$selectTestdata == "Global Patterns (environmental samples)"){
      gp <- readRDS("testdata/GlobalPatterns")
      dat <- as.data.frame(otu_table(gp))
      taxonomy <- as.data.frame(tax_table(gp))
      taxonomy <- addMissingTaxa(taxonomy)
      taxonomy[is.na(taxonomy)] <- "NA"
      meta <- as.data.frame(sample_data(gp))
      colnames(meta)[1]<-"SampleID"
      tree <- phy_tree(gp)
      
      normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
      normalized_dat = normalizeOTUTable(dat,normMethod)
      #tax_binning = taxBinning(normalized_dat[[2]],taxonomy)
      
      phylo_gp <- merge_phyloseq(otu_table(normalized_dat$norm_tab,T),tax_table(as.matrix(taxonomy)),sample_data(meta),tree)
      unifrac_dist <- readRDS("testdata/GlobalPatternsUnifrac")
      
      #the final dataset 
      dataset <- list(rawData=dat,metaData=meta,taxonomy=taxonomy,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab,tree=tree,phylo=phylo_gp,unifrac_dist=unifrac_dist,undersampled_removed=F,filtered=F, normMethod = normMethod)

      vals$datasets[["GlobalPatterns"]] <- dataset
      updateTabItems(session,"sidebar")
      removeModal()
    }
    if(input$selectTestdata == "Enterotype (facial samples)"){
      ent <- readRDS("testdata/enterotype")
      dat <- as.data.frame(otu_table(ent))
      rows <- rownames(dat)
      dat <- data.frame(lapply(dat, function(x){return(as.integer(x*10000))}))
      rownames(dat) <- rows
      taxonomy <- data.frame(tax_table(ent))
      taxonomy <- addMissingTaxa(taxonomy)
      taxonomy[is.na(taxonomy)] <- "NA"
      meta <- as.data.frame(sample_data(ent))
      names(meta)[names(meta) == "SampleID"] <- "SampleIDwNAs"
      names(meta)[names(meta) == "Sample_ID"] <- "SampleID"
      meta <- meta[,c(2,4,1,3,5,6,7,8,9)]
      
      tree <- NULL
      
      normalized_dat = normalizeOTUTable(dat,which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1)
      
      phylo_ent <- ent
      unifrac_dist <- NULL
      
      #the final datatset
      dataset <- list(rawData=dat,metaData=meta,taxonomy=taxonomy,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab,tree=tree,phylo=phylo_ent,unifrac_dist=unifrac_dist,undersampled_removed=F,filtered=F)
      
      vals$datasets[["Enterotype"]] <- dataset
      updateTabItems(session,"sidebar")
      removeModal()
    }
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
  
  # choose current dataset; return NULL if no set is yet uploaded
  currentSet <- eventReactive(input$datasets_rows_selected, {
    if(length(vals$datasets) == 0){
      return (NULL)
    }
    return(input$datasets_rows_selected)
  })
  
  ###obeservers####
  
  #observer for inputs depending on choosing a meta-group first
  observe({
    if(!is.null(currentSet())){
      phylo <- vals$datasets[[currentSet()]]$phylo
      meta <- data.frame(sample_data(phylo))
      taxonomy <- data.frame(tax_table(phylo))
      
      #filter variables
      filterColumnValues <- unique(meta[[input$filterColumns]])
      updateSelectInput(session,"filterColumnValues",choices=filterColumnValues)
      
      filterTaxaValues <- unique(taxonomy[[input$filterTaxa]])
      updateSelectInput(session,"filterTaxaValues",choices=filterTaxaValues)
      
      if(input$filterColumns == "NONE"){
        samples_left <- meta$SampleID
      }else if(input$filterColumns != "" && input$filterColumnValues != ""){
        samples_left <- meta[eval(input$filterColumns) == input$filterColumnValues,]$SampleID
      }else{
        samples_left <- NULL
      }
      updateSelectInput(session,"filterSample",choices=samples_left)
      
      #basic network variables
      groupVariables <- unique(meta[[input$groupCol]])
      updateSelectInput(session,"groupVar1",choices = groupVariables) 
    }
  })
  
  
  # update input selections
  observe({
    if(!is.null(currentSet())){  
      #get tables from phyloseq object
      phylo <- vals$datasets[[currentSet()]]$phylo
      otu <- otu_table(phylo)
      meta <- data.frame(sample_data(phylo))
      taxonomy <- data.frame(tax_table(phylo))
      #if(!is.null(phy_tree(vals$datasets[[currentSet()]]$phylo))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      phylo <- vals$datasets[[currentSet()]]$phylo
      
      updateSliderInput(session,"rareToShow",min=1,max=ncol(otu),value=min(50,ncol(otu)))
      updateSliderInput(session,"rareToHighlight",min=1,max=ncol(otu),value=round(ncol(otu)/10))
      updateSliderInput(session,"top_x_features",min=1,max=nrow(otu))
      if(ncol(meta)>2){enable("confounding_start")}
      
      updateSelectInput(session,"structureCompOne",choices=(1:nrow(meta)))
      updateSelectInput(session,"structureCompTwo",choices=(2:nrow(meta)))
      updateSelectInput(session,"structureCompThree",choices=(3:nrow(meta)))
      updateSelectInput(session,"pcaLoading",choices=(1:nrow(meta)))
      
      #update silder for binarization cutoff dynamically based on normalized dataset
      min_value <- min(otu)
      max_value <- round(max(otu)/16)
      updateNumericInput(session,"binCutoff",min=min_value,max=max_value)
      
      updateNumericInput(session,"k_in",value=0,min=0,max=vals$datasets[[currentSet()]]$vis_out$K,step=1)
      
      ########updates based on meta info########
      covariates <- vals$datasets[[currentSet()]]$vis_out$covariates
      updateSelectInput(session,"choose",choices = covariates)
      
      #pick all column names, except the SampleID
      group_columns <- setdiff(colnames(meta),"SampleID")
      updateSelectInput(session,"alphaGroup",choices = c("-",group_columns))
      updateSelectInput(session,"betaGroup",choices = group_columns)
      updateSelectInput(session,"structureGroup",choices = group_columns)
      updateSelectInput(session,"groupCol",choices = group_columns)
      updateSelectInput(session,"formula",choices = group_columns)
      updateSelectInput(session,"taxSample",choices=c("NULL",group_columns))
      
      #pick all categorical variables in meta dataframe (except SampleID)
      categorical_vars <- colnames(meta[,unlist(lapply(meta,is.character))])
      categorical_vars <- setdiff(categorical_vars,"SampleID")
      updateSelectInput(session,"forest_variable",choices = group_columns)
      updateSelectInput(session,"heatmapSample",choices = c("NULL",group_columns))
      
      if(is.null(access(phylo,"phy_tree"))) betaChoices="Bray-Curtis Dissimilarity" else betaChoices=c("Bray-Curtis Dissimilarity","Generalized UniFrac Distance")
      updateSelectInput(session,"betaMethod",choices=betaChoices)
      
      updateSelectInput(session,"phylo_color",choices= c("-","abundance",group_columns,"Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      updateSelectInput(session,"phylo_shape",choices = c("-",group_columns,"Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      updateSelectInput(session,"phylo_size",choices = c("-","abundance",group_columns))
      updateSliderInput(session,"phylo_prune",min=2,max=ntaxa(phylo),value=50,step=1)
      
      updateSelectInput(session,"filterColumns",choices = c("NONE",group_columns))
    }
  })
  
  #observer for legit variables for confounding analysis
  observe({
    if(!is.null(currentSet())){
      #factorize meta data
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      meta[] <- lapply(meta,factor)
      #pick all variables which have 2 or more factors for possible variables for confounding!
      tmp <- names(sapply(meta,nlevels)[sapply(meta,nlevels)>1])
      group_columns_no_single <- setdiff(tmp,"SampleID")
      updateSelectInput(session,"confounding_var",choices=group_columns_no_single)
    }
  })
  
  #this part needs to be in its own "observe" block
  #-> updates ref choice in section "functional topics"
  #-> also update var2 for basic network
  observe({
    if(!is.null(currentSet())){
      ref_choices <- unique(sample_data(vals$datasets[[currentSet()]]$phylo)[[input$formula]])
      updateSelectInput(session,"refs",choices=ref_choices)
      
      #do not display var chosen for var1 in var2 selection
      groupVariables <- unique(sample_data(vals$datasets[[currentSet()]]$phylo)[[input$groupCol]])
      remainingGroupVariables <- setdiff(groupVariables,input$groupVar1)
      updateSelectInput(session,"groupVar2",choices = c("all",remainingGroupVariables))
    }
  })
  
  # check for update if undersampled columns are to be removed (rarefation curves)
  observeEvent(input$excludeSamples, {
    if(!is.null(currentSet())){
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
  })
  
  
  ##### data filtering #####
  
  observeEvent(input$filterApplySamples, {
    if(!is.null(currentSet())){
      #save "old" dataset to reset filters later; only if there are no sample-filters applied to the current set
      if(!vals$datasets[[currentSet()]]$filtered){
        vals$datasets[[currentSet()]]$old.dataset <- vals$datasets[[currentSet()]]
      }

      if(input$filterColumns != "NONE"){
        #convert to datatable for filtering to work
        meta <- data.table(vals$datasets[[currentSet()]]$metaData, keep.rownames = F)
        
        #subset metatable by input 
        meta <- meta[get(input$filterColumns) == input$filterColumnValues,]
        if(!is.null(input$filterSample)){
          meta <- meta[meta$SampleID %in% input$filterSample,] 
        }
        
        #replace metaData
        vals$datasets[[currentSet()]]$metaData <- data.frame(meta, row.names = meta$SampleID)
        
        #remember that this dataset is filtered now
        vals$datasets[[currentSet()]]$filtered = T
        
        #build new dataset with filtered meta 
        filtered_samples <- as.vector(meta$SampleID)
        #adapt otu-tables to only have samples, which were not removed by filter
        vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[,filtered_samples]
        vals$datasets[[currentSet()]]$normalizedData <- vals$datasets[[currentSet()]]$normalizedData[,filtered_samples]
        vals$datasets[[currentSet()]]$relativeData <- vals$datasets[[currentSet()]]$relativeData[,filtered_samples]

        #build new phyloseq-object
        py.otu <- otu_table(vals$datasets[[currentSet()]]$normalizedData,T)
        py.tax <- tax_table(as.matrix(vals$datasets[[currentSet()]]$taxonomy))
        py.meta <- sample_data(data.frame(meta, row.names = meta$SampleID)) # with new meta df
        sample_names(py.meta) <- filtered_samples
        tree <- vals$datasets[[currentSet()]]$tree
        
        #cannot build phyloseq object with NULL as tree input; have to check both cases:
        if (!is.null(tree)) phylo <- merge_phyloseq(py.otu,py.tax,py.meta, tree) else phylo <- merge_phyloseq(py.otu,py.tax,py.meta)
        vals$datasets[[currentSet()]]$phylo <- phylo
        
        #re-calculate unifrac distance
        #pick correct subset of unifrac distance matrix, containing only the new filtered samples
        if(!is.null(tree)) unifrac_dist <- as.dist(as.matrix(vals$datasets[[currentSet()]]$unifrac_dist)[filtered_samples,filtered_samples]) else unifrac_dist <- NULL
        vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist
      }

    }
  })
  
  observeEvent(input$filterApplyTaxa,{
    if(!is.null(currentSet())){
      #save "old" dataset to reset filters later; only if there are no taxa-filters applied to the current set
      if(!vals$datasets[[currentSet()]]$filtered){
        vals$datasets[[currentSet()]]$old.dataset <- vals$datasets[[currentSet()]]
      }
      
      taxonomy <- data.table(vals$datasets[[currentSet()]]$taxonomy,keep.rownames = T)
      
      #subset taxonomy by input
      taxonomy <- taxonomy[get(input$filterTaxa) == input$filterTaxaValues,]
      #fix rownames and replace taxonomy
      taxonomy <- data.frame(taxonomy)
      rownames(taxonomy)<-taxonomy$rn
      remainingOTUs <- taxonomy$rn
      taxonomy$rn <- NULL
      vals$datasets[[currentSet()]]$taxonomy <- taxonomy
      vals$datasets[[currentSet()]]$filtered = T
      
      #adapt otu-tables to only have OTUs, which were not removed by filter
      vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[remainingOTUs,]

      #recalculate the relative abundances and normalize again 
      normalizedData <- normalizeOTUTable(vals$datasets[[currentSet()]]$rawData, vals$datasets[[currentSet()]]$normMethod)
      vals$datasets[[currentSet()]]$normalizedData <- normalizedData$norm_tab
      vals$datasets[[currentSet()]]$relativeData <- normalizedData$rel_tab
      
      #build new phyloseq-object
      py.otu <- otu_table(vals$datasets[[currentSet()]]$normalizedData,T)
      py.tax <- tax_table(as.matrix(taxonomy)) # with new taxonomy df
      py.meta <- sample_data(data.frame(vals$datasets[[currentSet()]]$metaData))
      tree <- vals$datasets[[currentSet()]]$tree
      
      #cannot build phyloseq object with NULL as tree input; have to check both cases:
      if (!is.null(tree)) phylo <- merge_phyloseq(py.otu,py.tax,py.meta, tree) else phylo <- merge_phyloseq(py.otu,py.tax,py.meta)
      vals$datasets[[currentSet()]]$phylo <- phylo
      
      #recalculate unifrac distance in this case
      if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalizedData$norm_tab, data.frame(vals$datasets[[currentSet()]]$metaData), tree) else unifrac_dist <- NULL
      vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist 
      
    }
  })
  
  observeEvent(input$filterResetA, {
    if(!is.null(currentSet())){
      #check if there filters applied to dataset
      if(vals$datasets[[currentSet()]]$filtered){
        print("reseting dataset")
        restored_dataset <- vals$datasets[[currentSet()]]$old.dataset
        vals$datasets[[currentSet()]] <- restored_dataset
        #remove old dataset from restored dataset
        vals$datasets[[currentSet()]]$old.dataset <- NULL
        #dataset is not filtered anymore
        vals$datasets[[currentSet()]]$filtered <- F
      }
    }
  })
  
  observeEvent(input$filterResetB, {
    if(!is.null(currentSet())){
      #check if there filters applied to dataset
      if(vals$datasets[[currentSet()]]$filtered){
        print("reseting dataset")
        restored_dataset <- vals$datasets[[currentSet()]]$old.dataset
        vals$datasets[[currentSet()]] <- restored_dataset
        #remove old dataset from restored dataset
        vals$datasets[[currentSet()]]$old.dataset <- NULL
        #dataset is not filtered anymore
        vals$datasets[[currentSet()]]$filtered <- F
      }
    }
  })
  
  #####################################
  #    Basic Analysis                 #
  #####################################
  # update targets table of the currently loaded dataset
  output$metaTable <- renderDataTable({
    if(!is.null(currentSet())){
      #datatable(sample_data(vals$datasets[[currentSet()]]$phylo),filter='top',options=list(searching=T,pageLength=20,dom="Blfrtip",scrollX=T),editable=F,rownames=F)
      meta_dt<-datatable(sample_data(vals$datasets[[currentSet()]]$phylo),filter='top',selection=list(mode='multiple'),options=list(pageLength=20,scrollX=T))
      meta_dt
    } 
    else datatable(data.frame(),options=list(dom="t"))
  },server=T)

  # Plot rarefaction curves
  output$rarefacCurve <- renderPlotly({
    if(!is.null(currentSet())){
      #needs integer values to work
      tab = as.matrix(vals$datasets[[currentSet()]]$rawData)
      class(tab)<-"integer"

      rarefactionCurve = lapply(1:ncol(tab),function(i){
        n = seq(1,colSums(tab)[i],by=20)
        if(n[length(n)]!=colSums(tab)[i]) n=c(n,colSums(tab)[i])
        drop(rarefy(t(tab[,i]),n))
      })
      
      #slope = apply(tab,2,function(x) rareslope(x,sum(x)-100))
      slopesDF = calcSlopes(rarefactionCurve,tab)
      slopes <- as.numeric(slopesDF[,2])
      #samples_order <- order(as.numeric(slopesDF[,2]), decreasing = TRUE)
      ordered_samples <- slopesDF[order(as.numeric(slopesDF[,2]),decreasing = T)]
      vals$undersampled = tail(ordered_samples,input$rareToHighlight)
      
      first = order(slopes,decreasing=T)[1]
      p <- plot_ly(x=attr(rarefactionCurve[[first]],"Subsample"),y=rarefactionCurve[[first]],text=paste0(colnames(tab)[first],"; slope: ",round(1e5*slopes[first],3),"e-5"),hoverinfo="text",color="high",type="scatter",mode="lines",colors=c("red","black"))
      for(i in order(slopes,decreasing=T)[2:ncol(tab)]){
        highslope = ifelse(i < input$rareToHighlight,2,1)
        p <- p %>% add_trace(x=attr(rarefactionCurve[[i]],"Subsample"),y=rarefactionCurve[[i]],text=paste0(colnames(tab)[i],"; slope: ",round(1e5*slopes[i],3),"e-5"),hoverinfo="text",color=c("low","high")[highslope],showlegend=F)
      }
      p %>% layout(title="Rarefaction Curves",xaxis=list(title="Number of Reads"),yaxis=list(title="Number of Species"))
      p
    }
  })
  
  # show undersampled samples
  output$undersampled <- renderText({
    paste0("The following samples might be undersampled:\n",paste0(vals$undersampled,collapse=", "))
  })
  
  #reactive for taxonomic binning
  taxBinningReact <- reactive({
    if(!is.null(currentSet())){
      phylo <- vals$datasets[[currentSet()]]$phylo
      rel_dat <- vals$datasets[[currentSet()]]$relativeData
      #tax_binning <- taxBinning(rel_dat,taxonomy)
      #create phyloseq-object with relative abundance data
      rel_phylo <- merge_phyloseq(otu_table(rel_dat,T),tax_table(phylo),sample_data(phylo))
      tax_binning <- taxBinningNew(if(input$taxaAbundanceType)rel_phylo else phylo)
      tax_binning
    }
  })
  
  # plot distribution of taxa
  output$taxaDistribution <- renderPlotly({
    if(!is.null(currentSet())){
      tab = taxBinningReact()[[which(c("Kingdom","Phylum","Class","Order","Family","Genus","Species")==input$filterTaxa)]]
      
      #taxa = ifelse(rowSums(tab)/ncol(tab)<(input$taxCutoff),"Other",rownames(tab))
      taxa = rownames(tab)
      if(any(taxa=="Other")){
        other <- tab[which(taxa=="Other"),]
        other <- colSums(other)
        
        tab <- tab[-which(taxa=="Other"),]
        tab <- rbind(tab,Other=other)
      }
      
      tab = melt(tab)
      
      plot_ly(tab,name=~Var1,x=~value,y=~Var2,type="bar",orientation="h") %>%
        layout(xaxis=list(title=ifelse(input$taxaAbundanceType,"Cumulative Relative Abundance (%)", "Absolute Abundance")),yaxis=list(title="Samples"),
               barmode="stack",showlegend=T) #,legend=list(orientation="h")
    } else plotly_empty()
    
  })
  
  
  #reactive for PCR, UMAP, tSNE
  structureReact <- reactive({
    if(!is.null(currentSet())){
      mat <- otu_table(vals$datasets[[currentSet()]]$phylo)
      mat_t <- t(mat)
      
      samples = colnames(mat)
      taxa = colnames(mat_t)
      
      #PCR-calculation:
      pca = prcomp(mat,center=T,scale=T)
      pca_t = prcomp(mat_t,center=T,scale=T)
      out_pca = data.frame(pca$rotation,txt=samples)
      percentage = signif(pca$sdev^2/sum(pca$sdev^2)*100,2)
      
      loadings = data.frame(Taxa=taxa,loading=pca_t$rotation[,as.numeric(input$pcaLoading)])
      loadings = loadings[order(loadings$loading,decreasing=T),][c(1:10,(nrow(loadings)-9):nrow(loadings)),]
      loadings$Taxa = factor(loadings$Taxa,levels=loadings$Taxa)
      
      #UMAP:
      UMAP = umap(t(mat),n_components=3)
      out_umap = data.frame(UMAP$layout,txt=samples)
      
      #tSNE:
      tsne = Rtsne(t(mat),dim=3,perplexity=min(floor((length(samples)-1)/3),30))
      out_tsne = data.frame(tsne$Y,txt=samples)
      
      l <- list(percentage=percentage, loadings=loadings, out_pca = out_pca, out_umap=out_umap, out_tsne=out_tsne)
      return(l)
    }
  })
  
  #observer if PCA is chosen:
  observeEvent(input$structureMethod,{
    if(input$structureMethod=="PCA"){
      shinyjs::show("structureCompChoosing",anim = T)
    }else{
      shinyjs::hide("structureCompChoosing",anim = T)
    }
  })

  # visualize data structure as PCA, t-SNE or UMAP plots
  output$structurePlot <- renderPlotly({
    if(!is.null(structureReact())){
      mode = input$structureDim
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      
      if(input$structureMethod=="PCA"){
        out = structureReact()$out_pca
        percentage = structureReact()$percentage

        xlab = paste0("PC1 (",percentage[as.numeric(input$structureCompOne)]," % variance)")
        ylab = paste0("PC2 (",percentage[as.numeric(input$structureCompTwo)]," % variance)")
        zlab = paste0("PC3 (",percentage[as.numeric(input$structureCompThree)]," % variance)")
      }
      else if(input$structureMethod=="t-SNE"){
        out = structureReact()$out_tsne
      }
      else if(input$structureMethod=="UMAP"){
        out = structureReact()$out_umap
      }
      
      xlab = input$structureCompOne
      ylab = input$structureCompTwo
      zlab = input$structureCompThree
      colnames(out) = c(paste0("Dim",1:(ncol(out)-1)), "txt")
      
      if(mode=="2D"){
        plot_ly(out,x=as.formula(paste0("~Dim",input$structureCompOne)),y=as.formula(paste0("~Dim",input$structureCompTwo)),color=meta[[input$structureGroup]],text=~txt,hoverinfo='text',type='scatter',mode="markers") %>%
          layout(xaxis=list(title=xlab),yaxis=list(title=ylab))
      } else{
        plot_ly(out,x=as.formula(paste0("~Dim",input$structureCompOne)),y=as.formula(paste0("~Dim",input$structureCompTwo)),z=as.formula(paste0("~Dim",input$structureCompThree)),color=meta[[input$structureGroup]],text=~txt,hoverinfo='text',type='scatter3d',mode="markers") %>%
          layout(scene=list(xaxis=list(title=xlab),yaxis=list(title=ylab),zaxis=list(title=zlab)))
      }
    } else{
      plotly_empty()
    }
  })
  
  # plot PCA loadings
  output$loadingsPlot <- renderPlotly({
    if(!is.null(structureReact())){
      loadings = structureReact()$loadings
      
      plot_ly(loadings,x=~Taxa,y=~loading,text=~Taxa,hoverinfo='text',type='bar',color=I(ifelse(loadings$loading>0,"blue","red"))) %>%
        layout(title="Loadings",xaxis=list(title="",zeroline=F,showline=F,showticklabels=F,showgrid=F),yaxis=list(title=paste0("loadings on PC",input$pcaLoading)),showlegend=F) %>% hide_colorbar()
    } else{
      plotly_empty()
    }
  })
  
  # plot correlations between OTUs
  output$OTUcorrPlot <- renderPlotly({
    if(!is.null(currentSet())){
      mat <- vals$datasets[[currentSet()]]$rawData
      
      corCluster = corclust(t(mat),method="ward.D2")
      cor_mat = cor(t(mat))[corCluster$cluster.numerics$order,corCluster$cluster.numerics$order]
      
      clusters = cutree(corCluster$cluster.numerics,k=input$corrCluster)[corCluster$cluster.numerics$order]
      cuts = c(0,cumsum(table(clusters)[unique(clusters)]))-0.5
      shapes=list()
      for(i in 1:input$corrCluster){
        shapes[[i]] = list(type="rect",
                           line=list(color="red",width=3),
                           x0=cuts[i],x1=cuts[i+1],
                           y0=cuts[i],y1=cuts[i+1])
      }
      
      plot_ly(x=colnames(cor_mat),y=rownames(cor_mat),z=cor_mat,type="heatmap") %>%
        layout(shapes=shapes) %>%
        colorbar(title="Pearson Correlation\n")
    } else{
      plotly_empty()
    }
  })
  
  #reactive alpha-diversity table; stores all measures for alpha-div of current set
  alphaReact <- reactive({
    if(!is.null(currentSet())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      
      alphaTab = data.frame(colnames(otu))
      for(i in c("Shannon Entropy","effective Shannon Entropy","Simpson Index","effective Simpson Index","Richness")){
        alphaTab = cbind(alphaTab,round(alphaDiv(otu,i),2))
      }
      colnames(alphaTab) = c("SampleID","Shannon Entropy","effective Shannon Entropy","Simpson Index","effective Simpson Index","Richness")
      metaColumn <- as.factor(meta[[input$alphaGroup]])
      alphaTab[,input$alphaMethod]
    }else{
      NULL
    }
  })
  
  #reactive table of explained variation; for each meta-variable calculate p-val and rsquare
  explVarReact <- reactive({
    if(!is.null(alphaReact())){
      OTUs <- data.frame(t(otu_table(vals$datasets[[currentSet()]]$phylo)))  #transposed otu-table (--> rows are samples, OTUs are columns)
      meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
      
      #alpha<-alphaReact()
      #alpha$SampleID <- NULL
      #variables <- cbind.as.data.frame(meta[rownames(OTUs),],alpha[rownames(OTUs),])
      variables <- data.frame(meta[rownames(OTUs),])
      variables$SampleID <- NULL

      
      plist <- vector()
      rlist <- vector()
      namelist <- vector()
      #iterate over all columns
      for (i in 1:dim(variables)[2]) {
        if(length(unique(variables[,i])) > 1){
          variables_nc <- completeFun(variables,i)
          #calculate distance matrix between OTUs (bray-curtis)
          BC <- vegdist(OTUs[which(row.names(OTUs) %in% row.names(variables_nc)),], method="bray")
          output <- adonis(BC ~ variables_nc[,i])
          pvalue <- output$aov.tab[1,"Pr(>F)"]
          rsquare <- output$aov.tab[1,"R2"]
          names <- names(variables_nc)[i]
          
          plist <- append(plist,pvalue)
          rlist <- append(rlist,rsquare)
          namelist <- append(namelist,names)
        }
      }
      
      df <- data.frame(Variable = namelist, pvalue = plist, rsquare = rlist)
      df
      
      
    }else{
      NULL
    }
  })
  
  # plot alpha diversity
  output$alphaPlot <- renderPlotly({
    if(!is.null(alphaReact())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      
      alphaTab = alphaReact()
      
      if(input$alphaGroup=="-") plot_ly(y=alphaTab,type='violin',box=list(visible=T),meanline=list(visible=T),x0=input$alphaMethod) %>% layout(yaxis=list(title="alpha Diversity",zeroline=F))
      else plot_ly(x=meta[[input$alphaGroup]],y=alphaTab,color=meta[[input$alphaGroup]],type='violin',box=list(visible=T),meanline=list(visible=T),x0=input$alphaMethod) %>% layout(yaxis=list(title="alpha Diversity",zeroline=F))
    }
    else plotly_empty()
  })
  
  # show alpha diversity table
  output$alphaTable <- renderTable({
    if(!is.null(alphaReact())){
      alphaReact()
    }
  })
  
  output$alphaTableDownload <- downloadHandler(
    filename=function(){paste("alpha_diversity.csv")},
    content = function(file){
      if(!is.null(alphaReact())){
        write.csv(alphaReact(),file,row.names = F)
      }
    }
  )
  
  output$explainedVariation <- renderTable({
    if(!is.null(explVarReact())){
      explVarReact()
    }
  })

  output$explainedVariationBar <- renderPlot({
    if(!is.null(explVarReact())){
      explVar <- explVarReact()
      explVar$Variable <- factor(explVar$Variable, levels = unique(explVar$Variable)[order(explVar$rsquare,decreasing = T)])
      ggplot(data=explVar,aes(x=Variable,y=rsquare))+
        geom_bar(stat = "identity",aes(fill=pvalue))+
        ggtitle("Explained Variation of meta variables")+
        theme(axis.text.x = element_text(angle = 90))+
        geom_text(aes(label=round(rsquare,digits = 4)),vjust=-1)
    }
  })
  
  #confounding analysis
  observeEvent(input$confounding_start,{
    if(!is.null(currentSet())){
      #only if pre-build distance matrix exists, this can be calcualted (depending on tree input)
      if (!is.null(vals$datasets[[currentSet()]]$unifrac_dist)){
        meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
        #remove first column --> SampleID column
        meta[,1]<-NULL
        
        #calulate confounding matrix
        withProgress(message="Calculating confounding factors...",value=0,{
          vals$datasets[[currentSet()]]$confounder_table <- calculateConfounderTable(var_to_test=input$confounding_var,
                                                                                     variables = meta,
                                                                                     distance=vals$datasets[[currentSet()]]$unifrac_dist,
                                                                                     seed=seed,
                                                                                     progress=T)
        })
      }
    }
  })
  
  output$confounding_table <- renderTable({
    if(!is.null(currentSet())){
      if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
        vals$datasets[[currentSet()]]$confounder_table$table
      }
    }
  })
  
  output$confounding_table_download <- downloadHandler(
    filename = function(){
      paste("confounding_factors.csv")
    },
    content = function(file){
      if(!is.null(currentSet())){
        if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
          write.csv(vals$datasets[[currentSet()]]$confounder_table$table,file,row.names = F)
        }
      }
    }
  )
  
  output$confounding_var_text <- renderUI({
    if(!is.null(currentSet())){
      if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
        HTML(paste0("<b> Chosen variable: </b> ",vals$datasets[[currentSet()]]$confounder_table$var))
      }
    }
  })
  
  #do calculation for beta diversity plots
  betaReactive <- reactive({
    if(!is.null(currentSet())){
      group <- input$betaGroup
      
      otu <- as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo))
      otu <- otu[,order(names(otu))]
      otu <- data.frame(t(otu))
      
      meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
      meta <- data.frame(meta[order(rownames(meta)),])
      meta_pos <- which(colnames(meta) == group)
      all_groups <- as.factor(meta[,meta_pos])
      
      if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      
      method <- ifelse(input$betaMethod=="Bray-Curtis Dissimilarity","brayCurtis","uniFrac")
      dist <- betaDiversity(otu=otu,meta=meta,tree=tree,method=method)
      
      all_fit <- hclust(dist,method="ward.D2")
      tree <- as.phylo(all_fit)
      
      adonis <- adonis(dist ~ all_groups)
      col = rainbow(length(levels(all_groups)))[all_groups]
      
      out <- list(dist=dist, col=col, all_groups=all_groups, adonis=adonis, tree=tree)
      return(out)
    }
  })
  
  # clustering tree of samples based on beta-diversity
  output$betaTree <- renderPlot({
    if(!is.null(betaReactive())){
      beta <- betaReactive()
      plot(beta$tree,type="phylogram",use.edge.length=T,tip.color=beta$col,label.offset=0.01)
      axisPhylo()
      tiplabels(pch=16,col=beta$col)
    }
  })
  
#  MDS plot based on beta-diversity
  output$betaMDS <- renderPlot({
    if(!is.null(betaReactive())){
      beta<-betaReactive()
      s.class(
        cmdscale(beta$dist,k=2),col=unique(beta$col),cpoint=2,fac=beta$all_groups,
        sub=paste("MDS plot of Microbial Profiles\n(p-value ",beta$adonis[[1]][6][[1]][1],")",sep="")
      )
    }
  })
  
  # NMDS plot based on beta-diversity
  output$betaNMDS <- renderPlot({
    if(!is.null(betaReactive())){
      beta<-betaReactive()
      meta_mds = metaMDS(beta$dist,k=2)
      s.class(
        meta_mds$points,col=unique(beta$col),cpoint=2,fac=beta$all_groups,
        sub=paste("MDS plot of Microbial Profiles\n(p-value ",beta$adonis[[1]][6][[1]][1],")",sep="")
      )
    }
  })
  
  #phylogenetic tree using phyloseq
  output$phyloTree <- renderPlot({
    if(!is.null(currentSet())){
      phylo <- vals$datasets[[currentSet()]]$phylo
      #need to look for tree object like this, did not find out better/cleaner way yet
      if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL

      #only can display a phylogenetic tree if tree object is present
      if(!is.null(tree)){
        #prune taxa to number the user sets with slider (changes size of tree)
        myTaxa = names(sort(taxa_sums(phylo), decreasing = TRUE)[1:input$phylo_prune])
        pruned_phylo <- prune_taxa(myTaxa, phylo)

        if(input$phylo_color == "-") phylo_color = NULL else phylo_color=input$phylo_color
        if(input$phylo_shape == "-") phylo_shape = NULL else phylo_shape=input$phylo_shape
        if(input$phylo_size == "-") phylo_size = NULL else phylo_size=input$phylo_size
        if(input$phylo_tiplabels == "-") phylo_label.tips = NULL else phylo_label.tips=input$phylo_tiplabels

        if(input$phylo_radial == T){
          plot_tree(pruned_phylo,method = input$phylo_method,color=phylo_color,shape = phylo_shape,size = phylo_size,label.tips = phylo_label.tips,ladderize = "left",plot.margin = input$phylo_margin)+coord_polar(theta = "y")
        }else{
          plot_tree(pruned_phylo,method = input$phylo_method,color=phylo_color,shape = phylo_shape,size = phylo_size,label.tips = phylo_label.tips,ladderize = input$phylo_ladderize,plot.margin = input$phylo_margin)
        }

      }
    }
  }, height=800)
  

  #javascript show/hide toggle for advanced options
  shinyjs::onclick("phylo_toggle_advanced",shinyjs::toggle(id="phylo_advanced",anim = T))
  
  #only use if otu sequenced together
  #only use if same ASVs were used (http://benjjneb.github.io/dada2/tutorial.html)
  
  #####################################
  #    Advanced analysis              #
  #####################################
  
  # plot heatmap of OTU abundances per sample
  output$abundanceHeatmap <- renderPlotly({
    if(!is.null(currentSet())){
      set.seed(seed)
      phylo <- vals$datasets[[currentSet()]]$phylo
      #save generalized unifrac distance as global variable to use it for heatmap
      gunifrac_heatmap <<- as.dist(vals$datasets[[currentSet()]]$unifrac_dist)
      hm_distance <- if(input$heatmapDistance == "gunifrac") "gunifrac_heatmap" else input$heatmapDistance
      if(input$heatmapSample != "NULL"){
        plot_heatmap(phylo,method = input$heatmapOrdination,distance = hm_distance, sample.label = input$heatmapSample)
      }else{
        plot_heatmap(phylo,method = input$heatmapOrdination,distance = hm_distance)
      }
    }
  })
  
  ####random forest models####
  
  #calculate confusion matrix using random forest aproach
  rForestDataReactive <- eventReactive(input$forest_start,{
    if(!is.null(currentSet())){
      withProgress(message="Predicting random forest model...",value=0,{
        set.seed(input$forest_seed)
        incProgress(1/4,message="preparing data...")
        meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
        otu_t <- as.data.frame(t(otu_table(vals$datasets[[currentSet()]]$phylo)))
        
        if(input$forest_clr){
          otu_t <- clr(otu_t)
        }
        
        combined_data <- buildForestDataset(meta, otu_t, input)
        class_labels <- as.factor(combined_data$variable)
        
        incProgress(1/4,message="partitioning dataset & resampling...")
        #partition dataset in training+testing; percentage can be set by user
        inTraining <- createDataPartition(combined_data$variable, p = input$forest_partition, list = FALSE)
        #create training & testing partitions
        training <- combined_data[inTraining,]
        testing <- combined_data[-inTraining,]
        #parameters for training resampling (cross-validation)
        fitControl <- trainControl(classProbs=T,
                                   savePredictions=T,
                                   summaryFunction=twoClassSummary,
                                   method=input$forest_resampling_method,
                                   number=input$forest_cv_fold,
                                   repeats=input$forest_cv_repeats)
        
        #train random forest with training partition
        incProgress(1/4,message="training model...")
        
        if(input$forest_type == "random forest"){
          #tuning parameters:
          tGrid <- expand.grid(
            .mtry=extract(input$forest_mtry),
            .splitrule=input$forest_splitrule,
            .min.node.size=extract(input$forest_min_node_size)
          )
          if(input$forest_default){
            #use default values
            #model<-train(x=training,y=class_labels[inTraining],method = "ranger",trControl=fitControl,metric="ROC")
            model<-train(variable~.,
                         data=training,
                         method = "ranger",
                         trControl=fitControl,
                         metric="ROC",
                         importance="impurity",
                         num.threads=ncores)
          }else{
            model <- train(variable~.,
                           data=training,
                           method="ranger",
                           tuneGrid = tGrid,
                           importance=input$forest_importance,
                           trControl=fitControl,
                           num.trees=input$forest_ntrees,
                           metric="ROC",
                           num.threads=ncores)
          }
          
        }else if (input$forest_type == "gradient boosted model"){
          if(input$forest_default) {
            tGrid<-expand.grid(n.trees=100,interaction.depth=1,shrinkage=.1,n.minobsinnode=10)
          }else{
            tGrid <- expand.grid(
              n.trees=extract(input$gbm_ntrees),
              interaction.depth=extract(input$gbm_interaction_depth),
              shrinkage=extract(input$gbm_shrinkage),
              n.minobsinnode = extract(input$gbm_n_minobsinoode)
            )
          }
          model <- train(variable~.,
                         data=training,
                         method="gbm",
                         tuneGrid = tGrid,
                         trControl=fitControl,
                         metric="ROC")
        }
        
        #test model with testing dataset
        incProgress(1/4,message="testing model on test dataset...")
        predictions_model <- predict(model, newdata=testing)
        predictions_model_full <- predict(model, newdata=combined_data)
        con_matrix<-confusionMatrix(data=predictions_model, reference= class_labels[-inTraining])
        con_matrix_full<-confusionMatrix(data=predictions_model_full, reference= class_labels)
        return(list(cmtrx=con_matrix,cmtrx_full=con_matrix_full,model=model,class_labels=class_labels))
      })
    }
  })
  

  #density plot for continous variables
  output$forest_continuous_density <- renderPlot({
    if(!is.null(currentSet())){
      meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
      if(is.numeric(meta[[input$forest_variable]])){
        g<-ggplot(data=meta,mapping=aes_string(x=input$forest_variable))+
          geom_density()+
          xlab(NULL)
        g
      }else{
        return (NULL)
      }
    }
  })
  
  #output plots using the reactive output of random forest calculation
  output$forest_con_matrix <- renderPlot({
    if(!is.null(rForestDataReactive())){
      draw_confusion_matrix(rForestDataReactive()$cmtrx)
    }
  })
  
  output$forest_con_matrix_full <- renderPlot({
    if(!is.null(rForestDataReactive())){
      draw_confusion_matrix(rForestDataReactive()$cmtrx_full)
    }
  })
  
  #ROC curve for model
  output$forest_roc <- renderPlot({
    if(!is.null(rForestDataReactive())){
      res<-evalm(rForestDataReactive()$model)
      res$roc
    }
  })
  


  output$forest_roc_cv <- renderPlot({
    if(!is.null(rForestDataReactive())){
      ldf <- lift_df(rForestDataReactive()$model,rForestDataReactive()$model$levels[1])
      ggplot(ldf) +
        geom_line(aes(1 - Sp, Sn, color = fold)) +
        scale_color_discrete(guide = guide_legend(title = "Fold"))+
        xlab("x")+ylab("y")
    }
  })
  
  output$forest_top_features <- renderPlot({
    if(!is.null(rForestDataReactive())){
      plot(varImp(rForestDataReactive()$model),top=input$top_x_features)
    }
  })

  #upload file to use model to predict variable for new sample(s)
  rForestPrediction <- eventReactive(input$forest_upload,{
    if(is.null(input$forest_upload_file)){showModal(errorModal(error_message = fileNotFoundError))}
    if(is.null(rForestDataReactive())){showModal(errorModal(error_message = "Please calculate the model first."))}
    else{
      new_sample <- read.csv(input$forest_upload_file$datapath,header=T,sep="\t",row.names=1,check.names = F)
      #transpose new_sample; same orientation as otu in model building 
      new_sample<-t(new_sample)
      if(is.null(new_sample)){showModal(errorModal(error_message = fileEmptyError))}
      else{
        model <- rForestDataReactive()$model
        #only use columns for prediction, which were used for model building
        model_predictors <- setdiff(colnames(model$trainingData),".outcome")
        if(!all(model_predictors %in% colnames(new_sample))){showModal(errorModal(error_message = inconsistentColumnsForest))}
        else{
          pred <- predict(model,newdata=new_sample)
          df <- data.frame(row.names = rownames(new_sample),prediction = pred)
        }
      }
    }
    return (df)
  }) 
  
  #output table for prediction of variables of new sample
  output$forest_prediction <- renderTable({
    if(!is.null(rForestPrediction())){
      rForestPrediction()
    }
  },rownames = T)
  
  #observer for random Forest menu items  
  observe({
    #switch between different versions of advanced options, depending on type of model
    if(input$forest_type == "random forest"){
      shinyjs::show("ranger_advanced",anim = T)
      shinyjs::hide("gbm_advanced")
    }else if (input$forest_type == "gradient boosted model"){
      shinyjs::hide("ranger_advanced")
      shinyjs::show("gbm_advanced",anim=T)
    }
    
    if(!is.null(currentSet())){
      meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
      group_columns <- setdiff(colnames(meta),"SampleID")
      #for continous_slider; update input based on varibale chosen in forest_variable; show/hide density plot
      if(is.numeric(meta[[input$forest_variable]])){
        updateSliderInput(session,"forest_continuous_slider",min=0,max=max(meta[[input$forest_variable]],na.rm = T),value = median(meta[[input$forest_variable]],na.rm = T))
        shinyjs::show("forest_continuous_density")
      }else{
        updateSliderInput(session,"forest_continuous_slider",min=0,max=1)
        shinyjs::hide("forest_continuous_density")
      }
      #for selecting OTUs which will not be used in rForest calculation
      otu_t<-as.data.frame(t(otu_table(vals$datasets[[currentSet()]]$phylo)))
      updateSelectInput(session, "forest_exclude",choices=colnames(otu_t))
      
      #for features to include in model calculation; remove feature from list, which will be predicted 
      features<-group_columns[group_columns != input$forest_variable]
      updateSelectInput(session,"forest_features",choices = features)
    }
  })
  
  #javascript show/hide toggle for advanced options
  shinyjs::onclick("forest_toggle_advanced",shinyjs::toggle(id="forest_advanced",anim = T))
  #show/hide exclude OTU-option if OTU abundances are to be used for model building
  shinyjs::onclick("forest_otu",shinyjs::toggle(id="forest_exclude",anim = T))
  
  #####################################
  #    Network analysis               #
  #####################################
  
  
  # show histogram of all OTU values -> user can pick cutoff for binarization here
  output$cutoffHist <- renderPlotly({
    if(!is.null(currentSet())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      dat <- log(as.data.frame(otu+1))
      
      plot_ly(x=unlist(dat),type="histogram") %>%
        layout(xaxis=list(title="log(OTU-values)"), yaxis = list(title="Frequency"),
               shapes=list(list(type="line",y0=0,y1=1,yref="paper",x0=log(input$binCutoff),
                                x1=log(input$binCutoff),line=list(color="black",width=2))))
    } else{
      plotly_empty()
    }
  })
  
  # cutoff heatmap
  output$boolHeat <- renderPlotly({
    if(!is.null(currentSet())){
      cutoff <- input$binCutoff
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      
      m <- as.matrix(otu)
      m <- apply(m,c(1,2),function(x){ifelse(x<cutoff,0,1)})
      
      melted_m<-melt(m)
      melted_m$Var1 <- as.factor(melted_m$Var1)
      melted_m$Var2 <- as.factor(melted_m$Var2)
      
      plot_ly(z=melted_m$value,x=melted_m$Var1,y=melted_m$Var2,type="heatmap",showscale=F)
      
    } else{
      plotly_empty()
    }
  })
  
  # check if button for new calculation of counts is clicked -> reload network with the new counts
  observeEvent(input$startCalc,{
    
    withProgress(message = 'Calculating Counts..', value = 0, {
      vals$datasets[[currentSet()]]$counts = generate_counts(OTU_table=data.frame(otu_table(vals$datasets[[currentSet()]]$phylo)),
                                                             meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo)),
                                                             group_column = input$groupCol,
                                                             cutoff = input$binCutoff,
                                                             fc = ifelse(input$useFC=="log2(fold-change)",T,F),
                                                             var1 = input$groupVar1,
                                                             var2 = input$groupVar2,
                                                             progress = T)
    })
  })
  
  #network reactive
  cooccurrenceReactive <- reactive({
    if(!is.null(currentSet())){
      if(!is.null(vals$datasets[[currentSet()]]$counts)){
        otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
        counts = vals$datasets[[currentSet()]]$counts
        tax = tax_table(vals$datasets[[currentSet()]]$phylo)
        
        #remove all rows of counts where value is 0 -> no edge existing
        counts<-counts[counts$value!=0,]
        
        #if cutoff value is higher than maximal number of count values -> display all
        if(input$networkCutoff >= length(counts$value)){networkCutoff <- length(counts$value)} else {networkCutoff <- input$networkCutoff}
        
        Links = counts[order(abs(counts$value),decreasing=T)[1:networkCutoff],]
        colnames(Links) = c("source","target","value")
        Links$source = as.character(Links$source); Links$target = as.character(Links$target); Links$valueToPlot = abs(Links$value)
        Links$valueToPlot = norm10(Links$valueToPlot)+0.1
        
        Nodes = data.frame(name=unique(c(Links$source,Links$target)),group="")
        Nodes$size = rowSums(otu[Nodes$name,])/1000
        
        if(input$netLevel!="-") Nodes$group = substring(tax[Nodes$name,input$netLevel],4)
        Nodes$group[Nodes$group==""] = "unknown"
        
        Links$source = match(Links$source,Nodes$name)-1
        Links$target = match(Links$target,Nodes$name)-1
        
        out<-list(Nodes=Nodes,Links=Links)
        out
      }
    }
  })
  
  output$nodeDegree <- renderPlot({
    if(!is.null(cooccurrenceReactive())){
      Links <- as.data.frame(cooccurrenceReactive()$Links) 
      dat<-c(Links$source,Links$target)
      ggplot(data=as.data.frame(dat),aes(x=dat))+
        geom_histogram(bins = input$nodeDegreeBins)+
        ggtitle("Node Degree Plot for current network \n (Degree = # of edges coming out of Node)")+
        xlab("Degree")+ylab("Density")
    }
  })
  
  
  # network plot
  output$basicNetwork <- renderForceNetwork({
    if(!is.null(cooccurrenceReactive())){
      Links <- as.data.frame(cooccurrenceReactive()$Links)  
      Nodes <- as.data.frame(cooccurrenceReactive()$Nodes)
      
      forceNetwork(Links,Nodes,Source="source",Target="target",Value="valueToPlot",NodeID="name",
                   Nodesize="size",Group="group",linkColour=c("red","green")[(Links$value>0)+1],zoom=T,legend=T,
                   bounded=T,fontSize=12,fontFamily='sans-serif',charge=-25,linkDistance=100,opacity = 1)
      
    }
  })
  
  
  #####################################
  #    themetagenomcis apporach       #
  #####################################
  
  #here all objects and values needed for the plots of themetagenomics are created and stored in vals$datasets[[currentSet()]]$vis_out
  #does not work with phyloseq-object (error if otu-table is of class phyloseq::otu_table)
  
  observeEvent(input$themeta,{
    withProgress(message='Calculating Topics..',value=0,{
      if(!is.null(currentSet())){
        #take otu table and meta file from user input
        otu <- data.frame(otu_table(vals$datasets[[currentSet()]]$phylo))
        #otu <- vals$datasets[[currentSet()]]$normalizedData
        meta <- vals$datasets[[currentSet()]]$metaData
        tax = vals$datasets[[currentSet()]]$taxonomy
        
        incProgress(1/7,message="preparing OTU data..")
        formula <- as.formula(paste0("~ ",input$formula))
        formula_char <- input$formula
        refs <- input$refs
        
        #create themetadata object
        obj <- prepare_data(otu_table = otu,
                            rows_are_taxa = T,
                            tax_table = tax,
                            metadata = meta,
                            formula=formula,
                            refs = refs,
                            cn_normalize = FALSE)
        
        incProgress(1/7,message = "finding topics..")
        K=input$K
        sigma_prior = input$sigma_prior
        #use themetadata object to find K topics
        topics_obj <- find_topics(themetadata_object=obj,
                                  K=K,
                                  sigma_prior = sigma_prior)
        
        incProgress(1/7,message="predicting topic functions..")
        #functions_obj <- predict.topics(topics_obj,reference_path = "themetagenomics/")
        
        incProgress(1/7,message = "generate gene.table")
        #topic_function_table <- gene.table(functions_obj)
        
        incProgress(1/7,message = "finding topic effects..")
        #measure relationship of covarite with samples over topics distribution from the STM
        topic_effects_obj <- est(topics_obj)
        
        
        #function_effects <- themetagenomics::est(functions_obj,topics_subset=3)
        
        class(topics_obj) <- "topics"
        incProgress(1/7,message = "preparing visualization..")
        vals$datasets[[currentSet()]]$vis_out <- prepare_vis(topic_effects_obj)
        vals$datasets[[currentSet()]]$vis_out$K <- K
        vals$datasets[[currentSet()]]$vis_out$sigma <- sigma_prior
        vals$datasets[[currentSet()]]$vis_out$formula <- formula_char
        vals$datasets[[currentSet()]]$vis_out$refs <- refs
        vals$datasets[[currentSet()]]$topic_effects <- topic_effects_obj
        #vals$datasets[[currentSet()]]$gene_table <- topic_function_table
        incProgress(1/7)
      }
    })
  })
  
  #make gene.table download-able
  # output$downloadGeneTable <- downloadHandler(
  #   filename = "gene_table.csv",
  #   content = function(file){
  #     write.csv(vals$datasets[[currentSet()]]$gene_table,file,row.names = T)
  #   }
  # )
  
  REL <- reactive({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      if (show_topic$k != 0){
        
        current_k <- paste0('T',show_topic$k)
        
        l <- input$lambda
        
        tinfo_k <- vis_out$tinfo[vis_out$tinfo$Category == current_k,]  #subset(tinfo,Category == current_k)
        rel_k <- l*tinfo_k$logprob + (1-l)*tinfo_k$loglift
        new_order <- tinfo_k[order(rel_k,decreasing=TRUE)[1:vis_out$taxa_bar_n],]
        new_order$Term <- as.character.factor(new_order$Term)
        new_order$Taxon <- vis_out$taxa[new_order$Term,input$taxon]
        new_order$Term <- factor(new_order$Term,levels=rev(new_order$Term),ordered=TRUE)
        
      }else{
        
        new_order <- vis_out$default
        new_order$Taxon <- vis_out$taxa[as.character.factor(new_order$Term),input$taxon]
        
      }
      
      new_order
    }
    
  })
  
  EST <- reactive({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    topic_effects <- vals$datasets[[currentSet()]]$topic_effects$topic_effects
    
    if(!is.null(vis_out) & input$choose != "Please start calculation above first!"){
      suppressWarnings({
        
        covariate <- input$choose
        
        est_mat <- topic_effects[[covariate]]$est
        
        df0 <- data.frame(topic=paste0('T',1:vis_out$K),
                          est=est_mat[,1],
                          lower=est_mat[,2],
                          upper=est_mat[,3],
                          sig=ifelse(1:vis_out$K %in% topic_effects[[covariate]]$sig,'1','0'))
        df0$sig <- factor(as.character(sign(df0$est) * as.numeric(as.character(df0$sig))),levels=c('0','1','-1'),ordered=TRUE)
        df <- df0[order(topic_effects[[covariate]][['rank']]),]
        df$topic <- factor(df$topic,levels=df$topic,ordered=TRUE)
        
        p_est <- ggplot(df,aes_(~topic,y=~est,ymin=~lower,ymax=~upper,color=~sig)) +
          geom_hline(yintercept=0,linetype=3)
        p_est <- p_est +
          geom_pointrange(size=2) +
          theme_minimal() +
          labs(x='',y='Estimate') +
          scale_color_manual(values=c('gray','indianred3','dodgerblue3'),drop=FALSE) +
          scale_fill_brewer(type='qual',palette=6,direction=-1) +
          theme(legend.position='none',
                axis.text.x=element_text(angle=-90,hjust=0,vjust=.5))
        
        list(p_est=p_est,k_levels=levels(df$topic),covariate=covariate,df0=df0)
        
      })
    }
  })
  
  output$est <- renderPlotly({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out)){
        suppressWarnings(ggplotly(EST()$p_est,source='est_hover',tooltip=c('topic','est','lower','upper'))) #Error in UseMethod: no applicable method for 'plotly_build' applied to an object of class "shiny.tag"
      }else{
        plotly_empty()
      }
    }

    
  })
  
  output$ord <- renderPlotly({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out) & !is.null(EST())){
        beta <- t(vis_out$beta)
        
        if (input$dist == 'hellinger'){
          
          d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'norm'),method='euclidean'),k=3,eig=TRUE)
          
        }else if (input$dist == 'chi2'){
          
          d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'chi.square'),method='euclidean'),k=3,eig=TRUE)
          
        }else if (input$dist == 'jsd'){
          
          d <- cmdscale(proxy::dist(beta,jsd),k=3,eig=TRUE)   #woher kommt jsd?
          
        } else if (input$dist == 'tsne'){
          p <- 30
          d <- try(Rtsne::Rtsne(beta,k=,theta=.5,perplexity=p),silent=TRUE)
          while(class(d) == 'try-error'){
            p <- p-1
            d <- try(Rtsne::Rtsne(beta,k=3,theta=.5,perplexity=p),silent=TRUE)
          }
          if (p < 30) cat(sprintf('Performed t-SNE with perplexity = %s.\n',p))
          d$points <- d$Y
          d$eig <- NULL
          
        }else{
          
          d <- cmdscale(vegan::vegdist(beta,method=input$dist),k=3,eig=TRUE)
          
        }
        
        eig <- d$eig[1:3]/sum(d$eig)
        colnames(d$points) <- c('Axis1','Axis2','Axis3')
        df <- data.frame(d$points,EST()$df0)
        df$marg <- vis_out$topic_marg
        
        df$colors <- vis_out$colors[as.character(df$sig)]
        
        if (input$dim == '2d'){
          
          p1 <- plot_ly(df,source='ord_click')
          p1 <- add_trace(p1,
                          x=~Axis1,y=~Axis2,size=~marg,
                          type='scatter',mode='markers',sizes=c(5,125),
                          color=I(df$colors),opacity=.5,
                          marker=list(symbol='circle',sizemode='diameter',line=list(width=3,color='#FFFFFF')),
                          text=~paste('<br>Topic:',topic),hoverinfo='text')
          p1 <- layout(p1,
                       showlegend=FALSE,
                       xaxis=list(title=sprintf('Axis 1 [%.02f%%]',eig[1]*100),
                                  showgrid=FALSE),
                       yaxis=list(title=sprintf('Axis 2 [%.02f%%]',eig[2]*100),
                                  showgrid=FALSE),
                       paper_bgcolor='rgb(243, 243, 243)',
                       plot_bgcolor='rgb(243, 243, 243)')
          p1 <- add_annotations(p1,x=df$Axis1,y=df$Axis2,text=df$topic,showarrow=FALSE,
                                font=list(size=10))
          
          h <- event_data('plotly_hover',source='est_hover')
          
          if (!is.null(h)){
            k <- EST()$k_levels[h[['x']]]
            
            if (length(k) > 0){
              df_update <- df[df$topic == k,]
              
              
              if (df_update$sig == '1') df_update$sig <- '2' else if(df_update$sig== '-1') df_update$sig <- '-2' else df_update$sig<- '00'
              df_update$colors <- vis_out$colors[df_update$sig]
              
              p1 <- add_markers(p1,
                                x=df_update$Axis1,y=df_update$Axis2,opacity=.8,color=I(df_update$color),
                                marker=list(size=150,symbol='circle',sizemode='diameter',line=list(width=3,color='#000000')))
            }
            
            p1
            
          }
          
        }
        
        if (input$dim == '3d'){
          
          p1 <- plot_ly(df,source='ord_click',
                        x=~Axis1,y=~Axis2,z=~Axis3,size=~marg,
                        type='scatter3d',mode='markers',sizes=c(5,125),
                        color=I(df$colors),opacity=.5,
                        marker=list(symbol='circle',sizemode='diameter'),
                        text=~paste('<br>Topic:',topic),hoverinfo='text')
          
          p1 <- layout(p1,
                       showlegend=FALSE,
                       scene=list(
                         xaxis=list(title=sprintf('Axis 1 [%.02f%%]',eig[1]*100),
                                    showgrid=FALSE),
                         yaxis=list(title=sprintf('Axis 2 [%.02f%%]',eig[2]*100),
                                    showgrid=FALSE),
                         zaxis=list(title=sprintf('Axis 3 [%.02f%%]',eig[3]*100),
                                    showgrid=FALSE)),
                       paper_bgcolor='rgb(243, 243, 243)',
                       plot_bgcolor='rgb(243, 243, 243)')
          
        }
        
        p1
      }
    }
  })
  
  show_topic <- reactiveValues(k=0)
  
  #maybe not working...
  # observeEvent(event_data('plotly_click',source='ord_click'),{
  #   
  #   s <- event_data('plotly_click',source='ord_click')
  #   
  #   if (is.null(s)){
  #     
  #     show_topic$k <- 0
  #     updateNumericInput(session,'k_in',value=0)
  #     
  #   }else{
  #     
  #     t_idx <- s$pointNumber + 1
  #     updateNumericInput(session,'k_in',value=t_idx)
  #     show_topic$k <- t_idx
  #     
  #   }
  #   
  # })
  
  observeEvent(input$k_in,{
    show_topic$k <- input$k_in
  })
  
  observeEvent(input$reset,{
    show_topic$k <- 0
    updateNumericInput(session,'k_in',value=0)
  })
  
  output$bar <- renderPlot({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out)){
        if (show_topic$k != 0){
          p_bar <- ggplot(data=REL()) +
            geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=.6) +
            geom_bar(aes_(~Term,~Freq),stat='identity',fill='darkred',color='white')
        } else{
          p_bar <- ggplot(data=REL()) +
            geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=1)
        }
        
        p_bar +
          coord_flip() +
          labs(x='',y='Frequency',fill='') +
          theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
                legend.position='bottom') +
          viridis::scale_fill_viridis(discrete=TRUE,drop=FALSE) +
          guides(fill=guide_legend(nrow=2))
      }
    }
  })
  
  output$corr <- renderForceNetwork({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      topic_effects <- vals$datasets[[currentSet()]]$topic_effects$topic_effects
      if(!is.null(vis_out)){
        suppressWarnings(effects_sig <- topic_effects[[EST()$covariate]][['sig']])  #Warning: Error in [[: attempt to select less than one element in get1index
        K <- nrow(vis_out$corr$posadj)
        
        suppressWarnings({suppressMessages({
          g <- igraph::graph.adjacency(vis_out$corr$posadj,mode='undirected',
                                       weighted=TRUE,diag=FALSE)
          
          wc <- igraph::cluster_walktrap(g)
          members <- igraph::membership(wc)
          
          g_d3 <- networkD3::igraph_to_networkD3(g,group=members)
          
          g_d3$links$edge_width <- 10*(.1+sapply(seq_len(nrow(g_d3$links)),function(r) vis_out$corr$poscor[g_d3$links$source[r]+1,g_d3$links$target[r]+1]))
          g_d3$nodes$color <- 25*ifelse(1:K %in% effects_sig,1,0)*sign(topic_effects[[EST()$covariate]]$est[,1])
          g_d3$nodes$node_size <- 10*(.5+norm10(c(0,abs(topic_effects[[EST()$covariate]]$est[,1])))[-1])
          g_d3$nodes$name <- paste0('T',g_d3$nodes$name)
          
          networkD3::forceNetwork(Links=g_d3$links,Nodes=g_d3$nodes,
                                  Source='source',Target='target',
                                  charge=-25,
                                  opacity=1,
                                  fontSize=12,
                                  zoom=TRUE,
                                  bounded=TRUE,
                                  NodeID='name',
                                  fontFamily='sans-serif',
                                  opacityNoHover=.7,
                                  Group='color',
                                  Value='edge_width',
                                  Nodesize='node_size',
                                  linkColour='#000000',
                                  linkWidth=networkD3::JS('function(d) {return d.value;}'),
                                  radiusCalculation=networkD3::JS('d.nodesize'),
                                  colourScale=networkD3::JS("color=d3.scaleLinear()\n.domain([-1,0,1])\n.range(['blue','gray','red']);"))
          
        })})
      }
    }
  })
  
  
  #####################################
  #    SPIEC-EASI                     #
  #####################################
  
  #WATCH OUT!! this tool needs the otu-table to be in format: OTU in column; sample in row
  #or simply the phyloseq class
  
  
  ## Meinshausen-Buhlmann's ##
  
  observeEvent(input$se_mb_start,{
    if(!is.null(currentSet())){
      withProgress(message = 'Calculating mb..', value = 0, {
        phylo <- vals$datasets[[currentSet()]]$phylo
        taxa <- tax_table(phylo)
        incProgress(1/2,message = "starting calculation..")
        se_mb <- spiec.easi(phylo, method = "mb", lambda.min.ratio = input$se_mb_lambda.min.ratio, nlambda = input$se_mb_lambda, pulsar.params = list(rep.num=input$se_mb_repnumber, ncores =ncores,seed=seed))
        incProgress(1/2,message = "building graph objects..")
        
        #pre-build graph object for phyloseq graph
        se_mb$ig <- adj2igraph(getRefit(se_mb), vertex.attr=list(name=taxa_names(phylo)))
        #pre-build grapg for interactive networkD3 graph
        nd3 <-igraph_to_networkD3(se_mb$ig, taxa)
        
        output$spiec_easi_mb_network <- renderPlot({
          plot_network(se_mb$ig,phylo,type = "taxa",color=as.character(input$mb_select_taxa))
        })
        output$spiec_easi_mb_network_interactive <- renderForceNetwork(forceNetwork(Links=nd3$links,Nodes=nd3$nodes,NodeID = "name",Group = as.character(input$mb_select_taxa),
                                                                                    zoom=T,legend=T,fontSize = 5,charge = -2,opacity = .9, height = 200,width = 100))
      })
    }
  })
  
  
  ## Glasso ## 
  
  observeEvent(input$se_glasso_start,{
    if(!is.null(currentSet())){
      withProgress(message = 'Calculating glasso..', value = 0, {
        py <- vals$datasets[[currentSet()]]$phylo
        taxa <- tax_table(py)
        incProgress(1/2,message = "starting calculation..")
        se_glasso <- spiec.easi(py, method = "glasso", lambda.min.ratio = input$glasso_mb_lambda.min.ratio, nlambda = input$glasso_mb_lambda, pulsar.params = list(rep.num=input$se_glasso_repnumber, ncores = ncores,seed=seed))
        incProgress(1/2,message = "building graph objects..")
        
        #pre-build graph object for phyloseq graph
        se_glasso$ig <- adj2igraph(getRefit(se_glasso), vertex.attr=list(name=taxa_names(py)))
        #pre-build grapg for interactive networkD3 graph
        nd3 <-igraph_to_networkD3(se_glasso$ig, taxa)
        
        output$spiec_easi_glasso_network <- renderPlot({
          plot_network(se_glasso$ig,py,type = "taxa",color=as.character(input$glasso_select_taxa))
        })
        output$spiec_easi_glasso_network_interactive <- renderForceNetwork(forceNetwork(Links=nd3$links,Nodes=nd3$nodes,NodeID = "name",Group = as.character(input$glasso_select_taxa),
                                                                                    zoom=T,legend=T,fontSize = 5,charge = -2,opacity = .9, height = 200,width = 100))
        
      })
    }
  })

  ## SparCC ##
  
  # observeEvent(input$se_sparcc_start,{
  #   if(!is.null(currentSet())){
  #     withProgress(message = 'Calculating SparCC..', value = 0, {
  #       otu <- vals$datasets[[currentSet()]]$normalizedDat
  #       otu.t <- t(otu)
  #       
  #       incProgress(1/2,message = "starting calculation..")
  #       se_sparcc <- sparcc(otu.t)
  #       
  #       #set threshold on matrix
  #       sparcc.graph <- abs(se_sparcc$Cor) >= input$se_sparcc_threshold
  #       diag(sparcc.graph) <- 0
  #       sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
  #       incProgress(1/2,message = "building graph..")
  #       vals$datasets[[currentSet()]]$se_sparcc$graph <- adj2igraph(sparcc.graph)
  #     })
  #   }
  # })
  # 
  # output$spiec_easi_sparcc_network <- renderPlot({
  #   if(!is.null(currentSet())){
  #     sparcc_ig <- vals$datasets[[currentSet()]]$se_sparcc$graph
  #     otu <- vals$datasets[[currentSet()]]$normalizedDat
  #     otu.t <- t(otu)
  #     if(!is.null(sparcc_ig) && !is.null(otu.t)){
  #       vsize  <- rowMeans(clr(otu.t, 1))+6
  #       am.coord <- layout.fruchterman.reingold(sparcc_ig)
  #       
  #       plot(sparcc_ig, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="SparCC")
  #     }
  #   }
  # }) 
  

  #####################################
  #    Text fields                    #
  #####################################
  output$welcome <- renderUI({
    welcomeText
  })
  
  output$authors <- renderUI({
    authorsText
  })
  
  output$welcome_ref <- renderUI({
    sourcesText
  })
  
  output$sourcePhyloseq <- renderUI({
    phyloseqSourceText
  })
  
  output$text1 <- renderUI({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out)){
        HTML(sprintf("Below are the results of a %s topic STM. The ordination of the topics over taxa distribution (left) and the frequencies of
                   the top %s taxa (in terms of saliency) across all topics. By selecting a topic, the relative
                   frequencies of the taxa within that topic are shown in red. The ordination figure can be shown in
                   either 2D or 3D and the ordination method can be adjusted. Lambda adjusts the relevance calculation.
                   Choosing the taxon adjusts the group coloring for the bar plot. Clicking Reset resets the topic selection.",
                     vis_out$K,vis_out$taxa_bar_n))
      }
    }
    
  })
  
  output$text2 <- renderUI({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out)){
        themetagenomicsText2
      }
    }
  })
  
  output$text3 <- renderUI({
    themetagenomicsText3
  })
  
  output$cutoff_title <- renderUI({
    coOcurrenceCutoffTitleText  })
  
  output$cutoff_text <- renderUI({
    coOcurrenceDistrText
  })
  
  output$heatmap_text <- renderUI({
    coOcurrenceHeatmapText
  })
  
  output$basic_additional <- renderUI({
    coOcurrenceCutoffText
  })
  
  output$basic_calc_title <- renderUI({
    coOcurrenceCountsTitleText
  })
  
  output$basic_calc_additional <- renderUI({
    coOcurrenceCountsText
  })
  
  output$input_variables <- renderUI({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out)){
        K <- vis_out$K
        sigma <- vis_out$sigma_prior
        formula<-vis_out$formula
        refs<-paste(vis_out$refs,collapse=", ")
        HTML(paste0("<b> Number of chosen topics (K): </b>",K,"<br>",
                    "<b> Value of sigma_prior: </b>",sigma,"<br>",
                    "<b> Group from META file: </b>",formula, "<br>",
                    "<b> Reference Level in this group: </b>",refs))
      }else{
        HTML(paste0("<b> Number of chosen topics (K): </b>","<br>",
                    "<b> Value of sigma_prior: </b>","<br>",
                    "<b> Group from META file: </b>","<br>",
                    "<b> Reference Level in this group: </b>"))
      }
    }
  })
  
  output$advanced_text <- renderUI({
    themetagenomicsTextTitle
  })
  
  output$topic_text <- renderUI({
    themetagenomicsTextTopics
  })
  
  output$sigma_text <- renderUI({
    themetagenomicsTextSigma
  })
  
  # output$spiec_easi_additional <- renderUI({
  #   HTML(paste0("<b>1:</b> absolute value of correlations below this threshold are considered zero by the inner SparCC loop."))
  # })
  
  output$info_inputdata <- renderUI({
    inputDataText
  })
  
  output$info_testdata <- renderUI({
    testdataText
  })
  
  output$forest_model_parameters <- renderPrint({
    if(!is.null(rForestDataReactive())){
      model<-rForestDataReactive()$model
      model$finalModel
    }
  })
  
  output$forest_model_variables <- renderPrint({
    if(!is.null(rForestDataReactive())){
      model<-rForestDataReactive()$model
      features<-setdiff(colnames(model$trainingData),".outcome")
      features
    }
  })
  
  output$heatmapOrdinationText <- renderUI({
    heatmapOrdinationText
  })
  
  output$heatmapSourceText <- renderUI({
    heatmapText2
  })
  
  output$spiecEasiSource <- renderUI({
    spiecEasiSourceText
  })
}