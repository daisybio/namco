finishedOtuUploadModal <- function(message=NULL){
  
  if(is.null(message)){
    text <- HTML(paste0("Check out the \"Filter & Overview\" tab to get started or move on to the analysis tabs!"))
    s <- "s"
  }else{
    m <- paste(unlist(message), collapse=",<br>")
    text <- HTML(paste0("Check out the \"Filter & Overview\" tab to get started or move on to the analysis tabs<br>",
                        "Additional Info: The following samples were not present in both OTU table and meta-file and were removed: <br>",m))
    s <- "l"
  }
  modal<-modalDialog(
    title = "Success! Upload of your dataset is finished.",
    text,
    easyClose = T, size=s
  )
  showModal(modal)
}

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
      otu = otu[!apply(is.na(otu)|otu=="",1,all),-ncol(otu)] # remove "empty" rows & remove taxonomy column
      otus <- row.names(otu) #save OTU names
      otu <- sapply(otu,as.numeric) #make OTU table numeric
      rownames(otu) <- otus
      message(paste0(Sys.time()," - OTU-table loaded (with taxonomy column):"))
      message(paste0(dim(otu)[1],"-", dim(otu)[2]))
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
      taxonomy <- replaceNATaxonomy(taxonomy)
      otu <- as.data.frame(sapply(otu,as.numeric)) #make OTU table numeric
      rownames(otu) <- otus
      message(paste0(Sys.time()," - OTU-table loaded (with taxonomy column): ", dim(otu)[1], " - ", dim(otu)[2]))
    }
    
    ## read meta file, replace sample-column with 'SampleID' and set it as first column in df ##
    meta <- read_csv_custom(input$metaFile$datapath, file_type="meta")
    if(is.null(meta)){stop(changeFileEncodingError, call. = F)}
    if(!(input$metaSampleColumn %in% colnames(meta))){stop(didNotFindSampleColumnError, call. = F)}
    sample_column_idx <- which(colnames(meta)==input$metaSampleColumn)
    colnames(meta)[sample_column_idx] <- sample_column    # rename sample-column 
    if (sample_column_idx != 1) {meta <- meta[c(sample_column, setdiff(names(meta), sample_column))]} # place sample-column at first position
    
    #subset both meta & otu to only have same samples
    intersect_samples <- Reduce(intersect, list(v1=meta[[sample_column]], v2=colnames(otu)))
    missing_samples <- unique(c(meta[!meta[[sample_column]] %in% intersect_samples,][[sample_column]], setdiff(colnames(otu), intersect_samples)))
    meta <- meta[meta[[sample_column]] %in% intersect_samples,]
    otu <- otu[,c(intersect_samples)]
    colnames(meta) <- gsub("-","_",colnames(meta))
    
    # remove columns with only NA values
    meta <- meta[, colSums(is.na(meta)) != nrow(meta)] 
    rownames(meta)=meta[[sample_column]]
    
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
                                            has_tax_nw=F,
                                            has_comp_nw=F,
                                            filterHistory="")
    updateTabItems(session,"sidebar", selected = "overview")
    waiter_hide()
    finishedOtuUploadModal(missing_samples)
    
  },error=function(e){
    print(e)
    showModal(errorModal(error_message = e))
    waiter_hide()
  })
  
})

