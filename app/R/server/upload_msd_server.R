observeEvent(input$msdStart, {
  
  message("Starting MSD upload ...")
  waiter_show(html = tagList(spin_rotating_plane(),"Starting MSD upload..."),color=overlay_color)
  
  tryCatch({
    
    if(input$msdDataName%in%names(vals$datasets)){stop(duplicateSessionNameError,call. = F)}
    
    if(!is.null(input$msdLink)){
      link <- input$msdLink
    }else{
      link <- input$msdLinkSelect
    }
    
    download_file_name <- paste0(tempdir(),"/msd_data.zip")
    msd_data_dir <- paste0(tempdir(), "/msd_data_dir")
    # download zip file
    download.file(url = link, destfile = download_file_name, method="wget", extra = "--no-check-certificate")
    if(!file.exists(download_file_name)){stop("Could not download the MSD data. Please check that the link you are using is correct.",call. = F)}
    unzip(download_file_name, exdir=msd_data_dir)
    
    if(input$msdOTUType == "zOTUs"){
      otu <- read_csv_custom(paste0(msd_data_dir, "/ZOTUs-Table.tab"), "otu")
      sequences <- readDNAStringSet(paste0(msd_data_dir,"/ZOTUs-Seqs.fasta"))
      names(sequences) <- unlist(lapply(names(sequences), function(x){y<-strsplit(x, split=" "); return(unlist(y)[1])}))
      tree <- ape::read.tree(paste0(msd_data_dir,'/ZOTUs-Tree-nj.tre'))
    }else if(input$msdOTUType == "S-OTUs"){
      otu <- read_csv_custom(paste0(msd_data_dir, "/SOTUs-Table.tab"), "otu")
      sequences <- readDNAStringSet(paste0(msd_data_dir,"/SOTUs-Seqs.fasta"))
      names(sequences) <- unlist(lapply(names(sequences), function(x){y<-strsplit(x, split=" "); return(unlist(y)[1])}))
      tree <- ape::read.tree(paste0(msd_data_dir,'/SOTUs-Tree-nj.tre'))
    }
    
    # handle taxonomy
    taxonomy = generateTaxonomyTable(otu, 'Taxonomy') # generate taxonomy table from TAX column
    otu = otu[!apply(is.na(otu)|otu=="",1,all),-ncol(otu)] # remove "empty" rows & remove taxonomy column
    otus <- row.names(otu) #save OTU names
    otu <- sapply(otu,as.numeric) #make OTU table numeric
    rownames(otu) <- otus
    message(paste0(Sys.time()," - OTU-table loaded:"))
    message(paste0(dim(otu)[1],"-", dim(otu)[2]))
    
    # relative abundance
    normalized_dat <- list(norm_tab=otu, rel_tab = relAbundance(otu))
    
    # handle meta file
    meta <- read_csv_custom(paste0(msd_data_dir, "/mapping.tab"), file_type="meta")
    
    if('dataset_name_in_features_table' %in% colnames(meta)){
      colnames(meta)[which(colnames(meta)=="dataset_name_in_features_table")]<-"SampleID"
    }else if('##DatasetID' %in% colnames(meta)){
      colnames(meta)[which(colnames(meta)=="##DatasetID")]<-"SampleID"
    }else{
      stop('No matching sample ID columns detected in metadata. Please check your URL(s) or contact the Namco support.')
    }
    
    meta[["SampleID"]] <- gsub('-','_', meta[["SampleID"]])
    rownames(meta) <- meta[["SampleID"]]
    has_meta <- T
    
    # build phyloseq
    py.otu <- otu_table(normalized_dat$norm_tab,T)
    py.tax <- tax_table(as.matrix(taxonomy))
    py.meta <- sample_data(meta)
    py.refseq <- refseq(sequences)
    phyloseq <- merge_phyloseq(py.otu,py.tax,py.meta, tree, py.refseq)
    
    # rest
    waiter_update(html = tagList(spin_rotating_plane(),"Calculating UniFrac distances ..."))
    unifrac_dist <- buildGUniFracMatrix(normalized_dat$norm_tab, tree)
    alphaTabFull <- createAlphaTab(data.frame(phyloseq@otu_table, check.names=F), data.frame(phyloseq@sam_data, check.names = F))
    
    # finish up
    message(paste0(Sys.time()," - final phyloseq-object: "))
    message(paste0("nTaxa: ", ntaxa(phyloseq)))
    message(paste0(Sys.time()," - Finished OTU-table data upload! "))
    
    vals$datasets[[input$msdDataName]] <- list(session_name=input$msdDataName,
                                            rawData=normalized_dat$norm_tab, # this is the raw data, since no normalization is applied during upload
                                            metaData=meta,
                                            taxonomy=taxonomy,
                                            counts=NULL,
                                            normalizedData=normalized_dat$norm_tab,
                                            relativeData=normalized_dat$rel_tab,
                                            tree=tree,
                                            phylo=phyloseq,
                                            phylo.raw=phyloseq,             # this phylo object will not be normalized, only filtered
                                            unifrac_dist=unifrac_dist,
                                            alpha_diversity=alphaTabFull,
                                            undersampled_removed=F,
                                            filtered=F, 
                                            normMethod = 0,
                                            is_fastq=F,          
                                            has_meta=has_meta,
                                            has_picrust=F,
                                            is_sample_data=F,
                                            is_restored=F,
                                            has_rf=F,
                                            has_diff_nw=F,
                                            has_tax_nw=F,
                                            has_comp_nw=F,
                                            has_omics=F,
                                            filterHistory="",
                                            namco_version=namco_version)
    
    updateTabItems(session,"sidebar", selected = "overview")
    waiter_hide()
    finishedOtuUploadModal()
    
  }, error=function(e){
    showModal(errorModal(error_message=e$message))
    waiter_hide()
  })
})