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
    #h4("Additional parameters:"),
    #fluidRow(
    #  column(10,wellPanel(
    #    numericInput("otu_abundance_cutoff", "ASVs with abundance over all samples below this value (in %) will be removed:", value=0.25, min=0, max=100, step=0.01),
    #  ))
    #),
    #br(),
    textInput("dataName","Enter a project name:",placeholder=paste0("Namco_project_",Sys.Date()),value=paste0("Namco_project_",Sys.Date())),
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

finishedOtuUploadModal <- modalDialog(
  title = "Success! Upload of your dataset is finished.",
  "Check out the \"Filter & Overview\" tab to get started or move on to the analysis tabs",
  easyClose = T, size="s"
)

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
  
  message("Starting OTU-table upload ...")
  waiter_show(html = tagList(spin_rotating_plane(),"Starting OTU-upload ..."),color=overlay_color)
  
  tryCatch({

    if(input$dataName%in%names(vals$datasets)){stop(duplicateSessionNameError,call. = F)}
    if(is.null(input$otuFile) || is.null(input$metaFile)){stop(otuOrMetaMissingError,call. = F)}
    if(!file.exists(input$otuFile$datapath)){stop(otuFileNotFoundError,call. = F)}
    if(!file.exists(input$metaFile$datapath)){stop(metaFileNotFoundError,call. = F)}
    
    ##read OTU (and taxa) file ##
    
    #case: taxonomy in otu-file -> no taxa file provided
    if(is.null(input$taxFile)){
      otu <- read_csv_custom(input$otuFile$datapath, "otu")
      if(is.null(otu)){stop(changeFileEncodingError, call. = F)}
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
      message(paste0(Sys.time()," - OTU-table loaded (with taxonomy column):"))
      message(paste0(dim(otu)[1], dim(otu)[2]))
    }
    #case: taxonomy in separate file -> taxonomic classification file
    else{
      otu <- read_csv_custom(input$otuFile$datapath, "otu")
      if(is.null(otu)){stop(changeFileEncodingError, call. = F)}
      otus <- row.names(otu) #save OTU names
      taxonomy <- read.csv(input$taxFile$datapath,header=T,sep="\t",row.names=1,check.names=F) #load taxa file
      if("taxonomy" %in% colnames(taxonomy)){taxonomy <- generateTaxonomyTable(taxonomy)} # if taxonomy needs to be split by ";"
      #check for consistent OTU naming in OTU and taxa file:
      if(!all(rownames(otu) %in% rownames(taxonomy))){stop(otuNoMatchTaxaError,call. = F)}
      taxonomy[is.na(taxonomy)] <- "NA"
      otu <- sapply(otu,as.numeric) #make OTU table numeric
      rownames(otu) <- otus
      message(paste0(Sys.time()," - OTU-table loaded (with taxonomy column): ", dim(otu)[1], " - ", dim(otu)[2]))
    }
    
    ## read meta file, replace sample-column with 'SampleID' and set it as first column in df ##
    #meta <- read.csv(input$metaFile$datapath,header=T,sep="\t", check.names=F)
    meta <- read_csv_custom(input$metaFile$datapath, file_type="meta")
    if(is.null(meta)){stop(changeFileEncodingError, call. = F)}
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
    message(paste0(Sys.time()," - Loaded meta file; colnames: "))
    message(paste(unlist(colnames(meta)), collapse = " "))
    
    
    ## read phylo-tree file ##
    if(is.null(input$treeFile)) tree = NULL else {
      if(!file.exists(input$treeFile$datapath)){stop(treeFileNotFoundError,call. = F)}
      tree = ape::read.tree(input$treeFile$datapath)
      message(paste0(Sys.time()," - Loaded phylogenetic tree"))
      
    }
    
    normalized_dat = list(norm_tab=otu, rel_tab = relAbundance(otu))

    #create phyloseq object from data (OTU, meta, taxonomic, tree)
    py.otu <- otu_table(normalized_dat$norm_tab,T)
    py.tax <- tax_table(as.matrix(taxonomy))
    py.meta <- sample_data(meta)
    
    #cannot build phyloseq object with NULL as tree input; have to check both cases:
    if (!is.null(tree)) phyloseq <- merge_phyloseq(py.otu,py.tax,py.meta, tree) else phyloseq <- merge_phyloseq(py.otu,py.tax,py.meta)
    
    #pre-build unifrac distance matrix
    if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalized_dat$norm_tab, tree) else unifrac_dist <- NULL
    
    message(paste0(Sys.time()," - final phyloseq-object: "))
    message(paste0("nTaxa: ", ntaxa(phyloseq)))
    message(paste0(Sys.time()," - Finished OTU-table data upload! "))
    
    vals$datasets[[input$dataName]] <- list(session_name=input$dataName,
                                            rawData=normalized_dat$norm_tab, # this is the raw data, since no normalization is applied during upload
                                            metaData=meta,
                                            taxonomy=taxonomy,
                                            counts=NULL,
                                            normalizedData=normalized_dat$norm_tab,
                                            relativeData=normalized_dat$rel_tab,
                                            tree=tree,
                                            phylo=phyloseq,
                                            unifrac_dist=unifrac_dist,
                                            undersampled_removed=F,
                                            filtered=F, 
                                            normMethod = 0,
                                            is_fastq=F,
                                            has_meta=T,
                                            has_picrust=F,
                                            is_sample_data=F,
                                            is_restored=F,
                                            has_rf=F,
                                            has_diff_nw=F,
                                            has_tax_nw=F)
    updateTabItems(session,"sidebar")
    removeModal()
    waiter_hide()
    showModal(finishedOtuUploadModal)
    
  },error=function(e){
    print(e)
    showModal(uploadOTUModal(failed=T,error_message = e))
  })
  
})

