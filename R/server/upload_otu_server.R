# Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
uploadOTUModal <- function(failed=F,error_message=NULL) {
  modalDialog(
    h4("Please provide pregenerated input files. For detailed information on how the files have to look, check out the Info & Settings tab on the left!"),
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
      #div(tags$b("The file you specified could not be loaded. Please check the Info tab and to confirm your data is in the correct format!",style="color: red;"))
      div(tags$b(error_message,style="color:red;"))
    },
    footer = tagList(
      modalButton("Cancel"),
      actionButton("upload_otu_ok","OK",style="background-color:blue; color:white")
    )
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
    meta <- meta[, colSums(is.na(meta)) != nrow(meta)] # remove columns with only NA values
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