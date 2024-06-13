#https://ofstack.com/Nginx/17072/nginx-upload-large-file-timeout-solution.html

observeEvent(input$upload_fastq_dada2, {
  
  message(paste0(Sys.time()," Starting DADA2 ... "))
  
  trim_primers <- ifelse(input$trim_primers_dada=="Yes", T, F)
  overlay_text <- "Starting DADA2 ..."
  is_paired <- input$fastqIsPaired
  if(is_paired) truncations <- c(input$truncFw, input$truncRv) else truncations <- c(input$truncFw)
  
  trim_primers <- c(input$trim_primers_fw_dada, input$trim_primers_rv_dada)
  
  waiter_show(html = tagList(spin_rotating_plane(),overlay_text),color=overlay_color)
  
  tryCatch({
    if(is.null(input$fastqFiles$datapath)){stop(noFileError, call. = F)}
    # store if fastQC has already been run
    fastqc_dir <- paste0(dirname(dirname(input$fastqFiles$datapath[1])),"/fastqc_out")
    if(dir.exists(fastqc_dir))fastqc_exists <- T else fastqc_exists <- F
  
    # load meta-file (..or not)
    m <- handleMetaFastqMode(meta_file = input$fastqMetaFile$datapath,
                             fastq_sample_column = input$metaSampleColumnFastq,
                             sample_column = sample_column)
    meta <- m$meta
    meta_file_path <- m$meta_file_path
    has_meta <- ifelse(is.null(input$fastqMetaFile$datapath), F, T)
    
    # check if fastQC folder exists
    fastqc_dir <- paste0(dirname(dirname(input$fastqFiles$datapath[1])),"/fastqc_out")
    if(dir.exists(fastqc_dir)) fastqc_exists <- T else fastqc_exists <- F
    
    # name of directory with fastq files
    waiter_update(html = tagList(spin_rotating_plane(),"Reading in files ..."))
    fastq_lst <- handle_fastqs(fastq_files = input$fastqFiles,
                               fastqc_exists = fastqc_exists,
                               sampleNameCutoff = input$sampleNameCutoff,
                               is_paired = is_paired)
    
    # collect results
    foreward_files <- fastq_lst$foreward_files
    reverse_files <- fastq_lst$reverse_files
    sample_names <- fastq_lst$sample_names
    dirname <- fastq_lst$dirname
   
    # check if all sample names are equal for mapping and fastq
    if (has_meta){if (!all(meta[[sample_column]] %in% sample_names)){stop(noEqualFastqPairsError, call.=F)}} 
    
    ##### starting DADA2 pipeline #####
    
    # trim&trunc files
    foreward_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
    if (is_paired) reverse_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
    names(foreward_files_filtered) <- sample_names
    if (is_paired) names(reverse_files_filtered) <- sample_names
    waiter_update(html = tagList(spin_rotating_plane(),"Filtering and trimming primers ..."))
    
    out_filter <- data.frame(filterAndTrim(fwd=foreward_files, 
                                           filt=foreward_files_filtered, 
                                           rev=if(is_paired) reverse_files else NULL, 
                                           filt.rev=if(is_paired) reverse_files_filtered else NULL, 
                                           truncLen=truncations, 
                                           trimLeft = trim_primers, 
                                           rm.phix=TRUE, 
                                           compress=TRUE, 
                                           multithread=TRUE, 
                                           maxEE = if(is_paired) c(2,2) else c(2)))
    
    files_filtered <- rownames(out_filter[out_filter$reads.out!=0,])        # get files(R1), which have more than 0 reads left after filtering
    samples_filtered <- sapply(strsplit(files_filtered, "_L001"), `[`, 1)       # get all samples, which have more than 0 reads left
    if(length(samples_filtered)==0){stop(noTaxaRemainingAfterFilterError, call.=F)}
    foreward_files_filtered <- foreward_files_filtered[samples_filtered]
    if (is_paired) reverse_files_filtered <- reverse_files_filtered[samples_filtered]
    message(paste0(Sys.time()," - Filtered fastqs: ", truncations,"; "))
    message(paste0(Sys.time(), " - Files with 0 reads after filtering: ", rownames(out_filter[out_filter$reads.out==0,])))
    
    # learn errors
    waiter_update(html = tagList(spin_rotating_plane(),"Learning Errors (foreward)..."))
    errF <- learnErrors(foreward_files_filtered, multithread=TRUE, nbases = 1e8, randomize = T)
    waiter_update(html = tagList(spin_rotating_plane(),"Learning Errors (reverse)..."))
    if (is_paired) errR <- learnErrors(reverse_files_filtered, multithread=TRUE, nbases = 1e8, randomize = T)
    message(paste0(Sys.time()," - Learned Errors. "))
    
    #dada2
    waiter_update(html = tagList(spin_rotating_plane(),"Sample inference ..."))
    if (is_paired) {
      dadaFs <- dada(foreward_files_filtered, err=errF, multithread=T)
      dadaRs <- dada(reverse_files_filtered, err=errR, multithread=T)
      dada_merged <- mergePairs(dadaFs, foreward_files_filtered, dadaRs, reverse_files_filtered)
    } else{
      dadaRs <- NULL
      dadaFs <- dada(foreward_files_filtered, err=errF, multithread=T)
      dada_merged <- dadaFs
    }

    message(paste0(Sys.time()," - Merged files. "))
    
    # create ASV table & removing chimeras
    waiter_update(html = tagList(spin_rotating_plane(),"Merging and removing chimeras ..."))
    seq_table <- makeSequenceTable(dada_merged)
    if(dim(seq_table)[2] == 0){stop(emptyASVTableError, call.=F)}
    seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=T)
    message(paste0(Sys.time()," - Created ASV table: ", dim(seq_table_nochim)[1], " - ", dim(seq_table_nochim)[2]))
    
    ##### done with DADA2 pipeline #####
    
    # calculate loss of reads during steps
    track <- calcReadLoss(out_filter, dadaFs, dadaRs, dada_merged, seq_table_nochim, sample_names, samples_filtered, is_paired)
    
    # assign taxonomy
    waiter_update(html = tagList(spin_rotating_plane(),"Assigning taxonomy ..."))
    taxa <- assignTaxonomy(seq_table_nochim, "data/taxonomy_annotation.fa.gz", multithread = T, tryRC = !is_paired)
    message(paste0(Sys.time()," - Assigned Taxonomy."))
    
    # build phylogenetic tree
    seqs <- getSequences(seq_table_nochim)
    if(input$buildPhyloTree=="Yes"){
      waiter_update(html = tagList(spin_rotating_plane(),"building phylogenetic tree ..."))
      tree <- buildPhyloTree(seqs, ncores)
    } else{tree<-NULL; unifrac_dist<-NULL}

    # combine results into phyloseq object
    waiter_update(html = tagList(spin_rotating_plane(),"Combining results & Normalizing ..."))
    cn_lst <- combineAndNormalize_dada2(seq_table = seq_table_nochim,
                                        taxonomy = taxa, 
                                        has_meta = has_meta, 
                                        meta = meta, 
                                        tree = tree, 
                                        sn = samples_filtered, 
                                        apply_filter = input$fastqApplyRelAbundanceFilter)
    
    if(!is.null(tree)){
      #pre-build unifrac distance matrix
      unifrac_dist <- buildGUniFracMatrix(otu_table(cn_lst$phylo), phy_tree(cn_lst$phylo))
    }else{unifrac_dist<-NULL}
    
    #pre-calculate alpha-diversity
    if(has_meta){
      alphaTabFull <- createAlphaTab(data.frame(cn_lst$phylo@otu_table, check.names=F), data.frame(cn_lst$phylo@sam_data, check.names = F))
    }else{
      alphaTabFull <- createAlphaTab(data.frame(cn_lst$phylo@otu_table, check.names=F))
    }
    
    # run FastQC for trimmed & filtered files
    if(!fastqc_exists){
      unlink(fastqc_dir, recursive = T)
      suppressMessages(fastqc(fq.dir = dirname(foreward_files_filtered)[1], qc.dir = fastqc_dir, threads = ncores, fastqc.path = fastqc.path))
      fastqc_fw <- list.files(fastqc_dir, pattern="F_filt_fastqc.zip", full.names = T)
      fastqc_rv <- list.files(fastqc_dir, pattern="R_filt_fastqc.zip", full.names = T)
    }else{
      fastqc_fw <- list.files(fastqc_dir, pattern="R1_001_fastqc.zip", full.names = T)
      fastqc_rv <- list.files(fastqc_dir, pattern="R2_001_fastqc.zip", full.names = T) 
    }
    
    # store all filepaths in one place
    if(!is_paired){
      reverse_files <- rep(NA, length(foreward_files))
      reverse_files_filtered <- rep(NA, length(foreward_files_filtered))
      fastqc_rv <- rep(NA, length(fastqc_fw))
    }
    raw_df <- data.frame(fw_files = foreward_files,
                         rv_files = reverse_files,
                         sample_names = sample_names)
    filtered_df <- data.frame(fw_files_filtered = foreward_files_filtered,
                              rv_files_filtered = reverse_files_filtered,
                              fastqc_fw = fastqc_fw,
                              fastqc_rv = fastqc_rv,
                              sample_names = samples_filtered)
    file_df <- merge(raw_df, filtered_df, all.x = T)
    log_files <- NULL
    
    
    message(paste0(Sys.time()," - Finished DADA2 upload!"))
    
    vals$datasets[[input$fastqDataName]] <- list(session_name=input$fastqDataName,
                                                 generated_files = list(file_df=file_df, log_files=log_files),
                                                 fastq_dir = dirname,
                                                 is_fastq = T,
                                                 track = track,
                                                 rawData=cn_lst$raw_asv,
                                                 metaData=cn_lst$meta,
                                                 taxonomy=cn_lst$taxonomy,
                                                 counts=NULL,
                                                 normalizedData=cn_lst$normalized_asv$norm_tab,
                                                 relativeData=cn_lst$normalized_asv$rel_tab,
                                                 tree=cn_lst$tree,
                                                 phylo=cn_lst$phylo,
                                                 phylo.raw=cn_lst$phylo,             # this phylo object will not be normalized, only filtered
                                                 unifrac_dist=unifrac_dist,
                                                 alpha_diversity=alphaTabFull,
                                                 undersampled_removed=F,
                                                 filtered=F,
                                                 normMethod = 4,                     # default to 10.000 reads
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
    
    updateBox('dada2_readloss_box', action = 'toggle')
    waiter_hide()
    showModal(finishedFastqUploadModal)
  },error=function(e){
    waiter_hide()
    message(e$message)
    showModal(errorModal(error_message = e$message))
  })
})


