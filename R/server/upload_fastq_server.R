# Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
uploadFastqModal <- function(failed=F,error_message=NULL) {
  modalDialog(
    h4("Please provide the directory with demutiplexed fastq files and corresponding meta file. For detailed information on how the files have to look, check out the Info & Settings tab on the left!"),
    hr(),
    fluidRow(
      column(6,fileInput("fastqFiles","Select all fastq files", multiple = T, accept = ".fastq")),
      column(6,fileInput("metaFile","Select Metadata File"))
    ),
    fluidRow(
      column(12, wellPanel(
        
      )),
      column(6, radioButtons("normMethod","Normalization Method",c("no Normalization","by Sampling Depth","by Rarefaction"),inline=T)),
      column(6, radioButtons("rm_spikes", "Remove spikes", c("Yes","No"), inline=T))
    ),
    textInput("dataName","Enter a project name:",placeholder="New_Project",value="New_Project"),
    if(failed) {
      #div(tags$b("The file you specified could not be loaded. Please check the Info tab and to confirm your data is in the correct format!",style="color: red;"))
      div(tags$b(error_message,style="color:red;"))
    },
    footer = tagList(
      modalButton("Cancel"),
      actionButton("upload_fastq_ok","OK",style="background-color:blue; color:white")
    )
  )
}

# launch upload dialog
observeEvent(input$upload_fastq, {
  showModal(uploadFastqModal())
})

observeEvent(input$upload_fastq_ok, {
  
  rm_spikes <- ifelse(input$rm_spikes=="Yes", T, F)
  overlay_text <- ifelse(rm_spikes, "Starting DADA2 & spike removal ...", "Starting DADA2 ...")
  logfile <- tempfile(pattern="dada2_log",fileext = ".log")
  
  waiter_show(html = tagList(spin_rotating_plane(),overlay_text),color=overlay_color)
  Sys.sleep(1)
  
  tryCatch({
    # load meta-file
    meta_file <- read.csv(input$metaFile$datapath,header=T,sep="\t",check.names=F)
    meta <- meta_file[, colSums(is.na(meta_file)) != nrow(meta_file)] # remove columns with only NA values
    rownames(meta)=meta[["#SampleID"]]

    #files get "random" new filename in /tmp/ directory when uploaded in docker -> change filename to the upload-name
    dirname <- dirname(input$fastqFiles$datapath[1])  # this is the file-path if the fastq files
    file.rename(from=input$fastqFiles$datapath,to=paste0(dirname,"/",input$fastqFiles$name))
    #fastq_files <- list.files(dirname, pattern = c(".fastq", ".fastq.gz"))
    #cat("blub",  file = logfile)

    # remove spikes with python script 
    if(rm_spikes){
      waiter_update(html = tagList(spin_rotating_plane(),"Removing spikes ..."))
      rm_spikes_outdir <- paste0(dirname,"/rm_spikes_out")
      rm_spikes_outdir_woSpikes <- paste0(rm_spikes_outdir,"/omapping")
      rm_spikes_outdir_wSpikes <- paste0(rm_spikes_outdir,"/osamples")
      rm_spikes_spikes_file <- paste0(rm_spikes_outdir,"/ospikes.txt")
      rm_spikes_stats_file <- paste0(rm_spikes_outdir,"/ostats.txt")
      dir.create(rm_spikes_outdir)
      rm_spikes_command = paste0("python3 ../src/rm_spikes.py data/spikes.fasta ",
                                 input$metaFile$datapath,
                                 " ", dirname,
                                 " ", rm_spikes_outdir_woSpikes,
                                 " ", rm_spikes_outdir_wSpikes,
                                 " ", rm_spikes_spikes_file,
                                 " ", rm_spikes_stats_file,
                                 " ", ncores)
      print(rm_spikes_command)
      out <- system(rm_spikes_command, wait = T)
      #fastq_files <- list.files(rm_spikes_outdir_woSpikes, pattern=c(".fastq", ".fastq.gz"))
      dirname <- rm_spikes_outdir_woSpikes
      cat("removed spikes with command: ", rm_spikes_command, file = logfile)
    }
    
    # collect fw & rv files 
    foreward_files <- sort(list.files(dirname, pattern = "_R1_001.fastq", full.names = T))
    reverse_files <- sort(list.files(dirname, pattern = "_R2_001.fastq", full.names = T))
    sample_names <- sapply(strsplit(basename(foreward_files), "_"), `[`, 1) # sample name: everything until first "_"
    if (length(foreward_files) != length(reverse_files)){stop(noEqualFastqPairsError, call.=F)}
    if (!all(meta[["SampleID"]] == sample_names)){stop(noEqualFastqPairsError, call.=F)} # check if all sample names are equal for mapping and fastq
    
    ##### starting DADA2 pipeline #####
    # http://benjjneb.github.io/dada2/tutorial.html
    
    # trim&trunc files
    foreward_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
    reverse_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
    names(foreward_files_filtered) <- sample_names
    names(reverse_files_filtered) <- sample_names
    waiter_update(html = tagList(spin_rotating_plane(),"Filtering ..."))
    out <- filterAndTrim(foreward_files, foreward_files_filtered, reverse_files, reverse_files_filtered, truncLen=c(280,200),
                         rm.phix=TRUE, compress=TRUE, multithread=ncores) # maxEE filter not used
    cat("reads after filtering: ", out , file = logfile)
    # remove files, which have 0 reads after filtering
    
    
    #TODO: remove primers!
    #TODO: handle singletons!
    
    # learn errors
    waiter_update(html = tagList(spin_rotating_plane(),"Learning Errors ..."))
    errF <- learnErrors(foreward_files_filtered, multithread=ncores)
    cat("calculated fw errors ... \n ", file = logfile)
    errR <- learnErrors(reverse_files_filtered, multithread=ncores)
    cat("calculated rv errors ... \n ", file = logfile)
    
    #dada2 
    waiter_update(html = tagList(spin_rotating_plane(),"Sample inference ..."))
    dadaFs <- dada(foreward_files_filtered, err=errF, multithread=ncores)
    dadaRs <- dada(reverse_files_filtered, err=errR, multithread=ncores)
    dada_merged <- mergePairs(dadaFs, foreward_files_filtered, dadaRs, reverse_files_filtered)
    
    # create ASV table & removing chimeras
    waiter_update(html = tagList(spin_rotating_plane(),"Merging and removing chimeras ..."))
    seq_table <- makeSequenceTable(dada_merged)
    seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=ncores, verbose=TRUE)
    
    ##### done with DADA2 pipeline #####
    
    # calculate loss of reads during steps
    getN <- function(x) sum(getUniques(x))
    if (length(sample_names) > 1){
      track <- cbind(sample_names, out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(dada_merged, getN), rowSums(seq_table_nochim))
    }else{
      track <- cbind(sample_names, out, getN(dadaFs), getN(dadaRs), getN(dada_merged), rowSums(seq_table_nochim))
    }
    colnames(track) <- c("sample","input", "filtered", "denoisedF", "denoisedR", "merged", "non_chimera")
    
    # assign taxonomy
    waiter_update(html = tagList(spin_rotating_plane(),"Assigning taxonomy ..."))
    taxa <- assignTaxonomy(seq_table_nochim, "R/testdata/taxonomy_annotation.fa.gz", multithread = ncores)
    
    
    # build phylogenetic tree
    # https://f1000research.com/articles/5-1492/v1
    seqs <- getSequences(seq_table_nochim)
    
    # combine results into phyloseq object
    waiter_update(html = tagList(spin_rotating_plane(),"Combining results & Normalizing ..."))
    phylo <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows = F),
                      sample_data(meta),
                      tax_table(taxa))
    
    dna <- Biostrings::DNAStringSet(taxa_names(phylo))
    names(dna) <- taxa_names(dna)
    phylo <- merge_phyloseq(phylo, dna)
    taxa_names(phylo) <- paste0("ASV", seq(ntaxa(phylo)))
    
    # create final objects with "real" ASV names
    asv_table_final <- as.data.frame(otu_table(phylo))
    taxonomy_final <- as.data.frame(tax_table(phylo))
    meta_final <- as.data.frame(sample_data(phylo))
    
    # normalization
    normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
    normalized_asv <- normalizeOTUTable(otu, normMethod)
    
    # store all filepaths in one place
    file_df <- data.frame(fw_files = foreward_files, fw_files_filtered= foreward_files_filtered,
                          rv_files = reverse_files, rv_files_filtered = reverse_files_filtered,
                          sample_name = sample_names)
    
    vals$datasets[[input$dataName]] <- list(fastq_files = file_df,
                                            fastq_dir = dirname,
                                            is_fastq = T,
                                            rawData=asv_table_final,
                                            metaData=meta_final,
                                            taxonomy=taxonomy_final,
                                            counts=NULL,
                                            normalizedData=normalized_dat$norm_tab,
                                            relativeData=normalized_dat$rel_tab,
                                            tree=NULL,
                                            phylo=phylo,
                                            unifrac_dist=NULL,
                                            undersampled_removed=F,
                                            filtered=F, 
                                            normMethod = normMethod)
    
    updateTabItems(session,"sidebar")
    removeModal()
    
    waiter_hide()
  },error=function(e){
    waiter_hide()
    print(e)
    showModal(uploadFastqModal(failed=T,error_message = e))
  })
})