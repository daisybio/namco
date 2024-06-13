
# run LotuS2 
observeEvent(input$upload_fastq_lotus2, {
  
  message(paste0(Sys.time()," Starting LotuS2 ... "))
  
  overlay_text <- "Starting LotuS2 ..."
  is_paired <- input$fastqIsPaired
  
  waiter_show(html = tagList(spin_rotating_plane(),overlay_text),color=overlay_color)
  
  tryCatch({
    # do not allow to change number of cores
    if(grepl('-t ', input$additional_params_lotus)){stop('You may not use the -t parameter. The number of cores are fixed by Namco.', call. = F)}
    
    if(is.null(input$fastqFiles$datapath)){stop(noFileError, call. = F)}
    
    # load meta-file (..or not)
    m <- handleMetaFastqMode(meta_file = input$fastqMetaFile$datapath,
                             fastq_sample_column = input$metaSampleColumnFastq,
                             sample_column = sample_column)
    meta <- m$meta
    meta_file_path <- m$meta_file_path
    has_meta <- ifelse(is.null(input$fastqMetaFile$datapath), F, T)
    
    # check if fastQC folder exists
    fastqc_dir <- paste0(dirname(dirname(input$fastqFiles$datapath[1])),"/fastqc_out")
    if(dir.exists(fastqc_dir) && length(list.files(fastqc_dir)) > 0) fastqc_exists <- T else fastqc_exists <- F
    
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
    fastq_dir <- fastq_lst$dirname
    
    # check if all sample names are equal for mapping and fastq
    if (has_meta){if (!all(meta[[sample_column]] %in% sample_names)){stop(noEqualFastqPairsError, call.=F)}} 
    
    ##### starting LotuS2 pipeline #####
    
    # create temporary mapping file for lotus to use demultiplexed files
    temp_meta <- tempfile(tmpdir = dirname(dirname(input$fastqFiles$datapath[1])))
    cmd <- paste0(lotus2, ' -create_map ', temp_meta, ' -i ', fastq_dir)
    system(cmd)
    
    # add path to usearch binary if it gets selected
    if(input$clustering_lotus == 'usearch'){
      if(file.exists('data/usearch')){
        usearch_bin <- normalizePath('data/usearch')
        cmd <- paste0(lotus2, ' -link_usearch ', usearch_bin)
        system(cmd)
      }else{
        stop(usearchBinaryNotFoundError, call. = F)
      }
    }
    
    # rename samples in this temp_meta file with correct sample_names
    temp_meta_df <- fread(temp_meta)
    temp_meta_df[['#SampleID']] <- sample_names
    # merge it with uploaded meta file (if available)
    if(has_meta){
      temp_meta_df <- merge(temp_meta_df, meta, by.x='#SampleID', by.y=sample_column)
    }
    write.table(temp_meta_df, temp_meta, quote = F, sep = '\t', row.names = F)
    
    # run main LotuS2 
    outdir <- paste0(dirname(temp_meta),'/lotus2_out')
    cmd <- paste0(lotus2, ' -m ',temp_meta ,' -i ', dirname(fastq_dir),' -o ', outdir, ' -t ', ncores, 
                  ' -forwardPrimer ', input$trim_primers_fw_lotus, 
                  ' -reversePrimer ', input$trim_primers_rv_lotus,
                  ' -CL ', input$clustering_lotus,
                  ' -tax4refDB data/taxonomy_annotation.fa.gz',
                  ' ', input$additional_params_lotus)
    waiter_update(html = tagList(spin_rotating_plane(),"LotuS2 pipeline running ..."))
    system(cmd)
    
    ##### reading results #####
    
    # load phyloseq biom file
    phylo <- import_biom(file.path(outdir, "OTU.biom"), 
                         treefilename = file.path(outdir, "OTUphylo.nwk"),
                         refseqfilename = file.path(outdir, "OTU.fna"))
    
    # get raw data (i.e. just take it from phylo object)
    cn_lst <- combineAndNormalize_lotus2(phylo = phylo,
                                         apply_filter = input$fastqApplyRelAbundanceFilter, 
                                         has_meta = has_meta,
                                         sample_names = sample_names, 
                                         sample_column = sample_column)
    
    # log files
    demulti.log <- paste0(outdir,'/LotuSLogS/demulti.log')
    run.log <- paste0(outdir,'/LotuSLogS/LotuS_run.log')
    
    ##### finalizing results #####
    
    # calculate Unifrac distance
    if(!is.null(cn_lst$tree)){unifrac_dist <- buildGUniFracMatrix(otu_table(cn_lst$phylo), phy_tree(cn_lst$phylo))}else{unifrac_dist <- NULL}
    
    #pre-calculate alpha-diversity
    if(has_meta){
      alphaTabFull <- createAlphaTab(data.frame(cn_lst$phylo@otu_table, check.names=F), data.frame(cn_lst$phylo@sam_data, check.names = F))
    }else{
      alphaTabFull <- createAlphaTab(data.frame(cn_lst$phylo@otu_table, check.names=F))
    }
    
    # run FastQC (place results next to meta file, output dir, fastq dir)
    if(!fastqc_exists){
      unlink(fastqc_dir, recursive = T)
      waiter_show(html = tagList(spin_rotating_plane(),"Generating FastQC plots ..."),color=overlay_color)
      suppressMessages(fastqc(fq.dir = fastq_dir, qc.dir = fastqc_dir, threads = ncores, fastqc.path = fastqc.path))
    }
    fastqc_fw <- list.files(fastqc_dir, pattern="R1_001_fastqc.zip", full.names = T)
    fastqc_rv <- list.files(fastqc_dir, pattern="R2_001_fastqc.zip", full.names = T) 
    
    # store all filepaths in one place
    if(!is_paired){
      reverse_files <- rep(NA, length(foreward_files))
      fastqc_rv <- rep(NA, length(fastqc_fw))
    }
    file_df <- data.frame(fw_files = foreward_files,
                         rv_files = reverse_files,
                         fastqc_fw = fastqc_fw,
                         fastqc_rv = fastqc_rv,
                         sample_names = sample_names)
    log_files <- list(demulti.log=demulti.log, run.log=run.log)
    
    message(paste0(Sys.time()," - Finished LotuS2 upload!"))
    
    vals$datasets[[input$fastqDataName]] <- list(session_name=input$fastqDataName,
                                                 generated_files = list(file_df=file_df, log_files=log_files),
                                                 fastq_dir = fastq_dir,
                                                 is_fastq = T,
                                                 track = NULL,
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
    
    waiter_hide()
    showModal(finishedFastqUploadModal)
    
  }, error = function(e){
    waiter_hide()
    message(e$message)
    showModal(errorModal(error_message = e$message))
  })
  
})