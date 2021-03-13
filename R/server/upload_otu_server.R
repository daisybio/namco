# Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
uploadOTUModal <- function(failed=F,error_message=NULL) {
  modalDialog(
    title = "UPLOAD OTU/ASV table!",
    HTML("<h5>[For detailed information on how the files have to look, check out the <b>Info & Settings</b> tab on the left!]</h5>"),
    hr(),
    h4("Files:"),
    fluidRow(
      column(6,wellPanel(fileInput("otuFile","Select OTU table"), style="background:#3c8dbc")),
      column(6,wellPanel(fileInput("metaFile","Select Metadata File"), style="background:#3c8dbc", 
                         textInput("metaSampleColumn", "Name of the sample-column:", value="SampleID")))
    ),
    fluidRow(
      column(6,wellPanel(checkboxInput("taxInOTU","Click here if the taxonomic classification is stored in a seperate file and not in the OTU-file:",F),
                         fileInput("taxFile","Select Taxonomic classification file"))),
      column(6,wellPanel(fileInput("treeFile","Select Phylogenetic Tree File (optional)",width="100%"), fontawesome::fa("tree", fill="red", height="1.5em")))
    ),
    hr(),
    h4("Additional parameters:"),
    fluidRow(
      column(10,wellPanel(radioGroupButtons("normMethod","Normalization Method",c("no Normalization","by Sampling Depth","by Rarefaction"), direction="horizontal")))
    ),
    br(),
    textInput("dataName","Enter a project name:",placeholder="New_Project",value="New_Project"),
    if(failed) {
      #div(tags$b("The file you specified could not be loaded. Please check the Info tab and to confirm your data is in the correct format!",style="color: red;"))
      div(tags$b(error_message,style="color:red;"))
    },
    footer = tagList(
      modalButton("Cancel", icon = icon("times-circle")),
      actionButton("upload_otu_ok","OK",style="background-color:blue; color:white")
    ),
    easyClose = T, fade = T, size = "l"
  )
}

# launch upload dialog
observeEvent(input$upload_otu, {
  showModal(uploadOTUModal())
})

#observer for taxInOTU checkbox
observeEvent(input$taxInOTU,{
  if(!input$taxInOTU){shinyjs::hide("taxFile")} else {shinyjs::show("taxFile")}
})

# try to load the dataset specified in the dialog window
observeEvent(input$upload_otu_ok, {
  
  cat("Starting OTU-table upload ...")
  
  tryCatch({

    if(input$dataName%in%names(vals$datasets)){stop(duplicateSessionNameError,call. = F)}
    if(is.null(input$otuFile) || is.null(input$metaFile)){stop(otuOrMetaMissingError,call. = F)}
    if(!file.exists(input$otuFile$datapath)){stop(otuFileNotFoundError,call. = F)}
    if(!file.exists(input$metaFile$datapath)){stop(metaFileNotFoundError,call. = F)}
    
    ##read OTU (and taxa) file ##
    
    #case: taxonomy in otu-file -> no taxa file provided
    if(is.null(input$taxFile)){
      otu <- read.csv(input$otuFile$datapath,header=T,sep="\t",row.names=1,check.names=F) #load otu table -> rows are OTU, columns are samples
      checked_taxa_column <- checkTaxonomyColumn(otu)
      if (checked_taxa_column[1] == F){
        wrong_rows = checked_taxa_column[3]
        err = checked_taxa_column[2]
        if (checked_taxa_column[2] == wrongTaxaColumnError){err = paste0(checked_taxa_column[2], wrong_rows)}
        stop(err, call. = F)
      }
      taxonomy = generateTaxonomyTable(otu) # generate taxonomy table from TAX column
      otu = otu[!apply(is.na(otu)|otu=="",1,all),-ncol(otu)] # remove "empty" rows
      otus <- row.names(otu) #save OTU names
      otu <- sapply(otu,as.numeric) #make OTU table numeric
      rownames(otu) <- otus
      cat(paste0(Sys.time()," - OTU-table loaded (with taxonomy column): ", dim(otu)))
      
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
      cat(paste0(Sys.time()," - OTU-table loaded (with taxonomy file): ", dim(otu)))
      
    }
    
    ## read meta file, replace sample-column with 'SampleID' and set it as first column in df ##
    meta <- read.csv(input$metaFile$datapath,header=T,sep="\t")
    if(!(input$metaSampleColumn %in% colnames(meta))){stop(didNotFindSampleColumnError, call. = F)}
    sample_column_idx <- which(colnames(meta)==input$metaSampleColumn)
    colnames(meta)[sample_column_idx] <- sample_column    # rename sample-column 
    if (sample_column_idx != 1) {meta <- meta[c(sample_column, setdiff(names(meta), sample_column))]}   # place sample-column at first position
    
    # remove columns with only NA values
    meta <- meta[, colSums(is.na(meta)) != nrow(meta)] 
    rownames(meta)=meta[[sample_column]]
    if(!setequal(colnames(otu),meta[[sample_column]])){stop(unequalSamplesError,call. = F)}
    meta = meta[match(colnames(otu),meta[[sample_column]]),]
    
    #set SampleID column to be character, not numeric (in case the sample names are only numbers)
    meta[[sample_column]] <- as.character(meta[[sample_column]])
    cat(paste0(Sys.time()," - Loaded meta file; colnames: ", colnames(meta)))
    
    
    ## read phylo-tree file ##
    if(is.null(input$treeFile)) tree = NULL else {
      if(!file.exists(input$treeFile$datapath)){stop(treeFileNotFoundError,call. = F)}
      tree = read.tree(input$treeFile$datapath)
      print(paste0(Sys.time()," - Loaded phylogenetic tree", colnames(meta)))
      
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
    
    cat(paste0(Sys.time()," - final phyloseq-object: "))
    cat(phylo)
    cat(paste0(Sys.time()," - Finished OTU-table data upload! "))
    
    vals$datasets[[input$dataName]] <- list(rawData=otu,
                                            metaData=meta,
                                            taxonomy=taxonomy,
                                            counts=NULL,
                                            normalizedData=normalized_dat$norm_tab,
                                            relativeData=normalized_dat$rel_tab,
                                            tree=tree,
                                            phylo=phylo,
                                            unifrac_dist=unifrac_dist,
                                            undersampled_removed=F,
                                            filtered=F, 
                                            normMethod = normMethod,
                                            is_fastq=F,
                                            has_meta=T)
    updateTabItems(session,"sidebar")
    removeModal()
    
  },error=function(e){
    print(e)
    showModal(uploadOTUModal(failed=T,error_message = e))
  })
  
})

errorModal <- function(error_message=NULL){
  modalDialog(
    p(error_message,style="color:red;"),
    easyClose = T,
    modalButton("Cancel")
  )
}