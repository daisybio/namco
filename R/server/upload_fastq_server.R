# Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
uploadFastqModal <- function(failed=F,error_message=NULL) {
  modalDialog(
    title = "UPLOAD fastq files:",
    HTML("<h5>[For detailed information on how the files have to look, check out the <b>Info & Settings</b> tab on the left!]</h5>"),
    hr(),
    h4("Files:"),
    fluidRow(
      column(6,wellPanel(fileInput("fastqFiles","Select all fastq files", multiple = T, accept = ".fastq"))),
      column(6,wellPanel(fileInput("metaFile","Select Metadata File"),
                         textInput("metaSampleColumn", "Name of the sample-column:", value="SampleID")))
    ),
    hr(),
    h4("Additional parameters:"),
    fluidRow(
      column(12, wellPanel(
        fluidRow(
          column(4, radioGroupButtons("rm_spikes", "Remove spikes", c("Yes","No"), direction="horizontal")),
          column(4, selectInput("trim_primers", "Trim Primers",choices = c("V3/V4")))
        ),
        div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("truncFw", "Truncation foreward:",value=280, min=1, max=500, step=1)),
        div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("truncRv", "Truncation reverse:",value=200, min=1, max=500, step=1)),
        radioGroupButtons("normMethod","Normalization Method",c("no Normalization","by Sampling Depth","by Rarefaction"), direction="horizontal")
      ))
    ),
    br(),
    textInput("dataName","Enter a project name:",placeholder="New_Project",value="New_Project"),
    if(failed) {
      #div(tags$b("The file you specified could not be loaded. Please check the Info tab and to confirm your data is in the correct format!",style="color: red;"))
      div(tags$b(error_message,style="color:red;"))
    },
    footer = tagList(
      modalButton("Cancel", icon = icon("times-circle")),
      actionButton("upload_fastq_ok","OK",style="background-color:blue; color:white")
    ),
    easyClose = T, fade = T, size = "l"
  )
}

# launch upload dialog
observeEvent(input$upload_fastq, {
  showModal(uploadFastqModal())
})

observeEvent(input$upload_fastq_ok, {
  
  print("Starting fastq data upload ...")
  
  rm_spikes <- ifelse(input$rm_spikes=="Yes", T, F)
  trim_primers <- ifelse(input$trim_primers=="Yes", T, F)
  overlay_text <- ifelse(rm_spikes, "Starting DADA2 & spike removal ...", "Starting DADA2 ...")
  trunc_fw <- as.numeric(input$truncFw)
  trunc_rv <- as.numeric(input$truncRv)
  
  waiter_show(html = tagList(spin_rotating_plane(),overlay_text),color=overlay_color)
  Sys.sleep(1)
  
  tryCatch({
    ## load meta-file and replace sample-column with 'SampleID' ##
    meta <- read.csv(input$metaFile$datapath,header=T,sep="\t",check.names=F)
    if(!(input$metaSampleColumn %in% colnames(meta))){stop(didNotFindSampleColumnError, call. = F)}
    sample_column_idx <- which(colnames(meta)==input$metaSampleColumn)
    colnames(meta)[sample_column_idx] <- sample_column           # rename sample-column 
    if (sample_column_idx != 1) {meta <- meta[c(sample_column, setdiff(names(meta), sample_column))]}   # place sample-column at first position
    
    meta <- meta[, colSums(is.na(meta)) != nrow(meta)] # remove columns with only NA values
    rownames(meta)=meta[[sample_column]]
    print(paste0(Sys.time()," - Loaded meta file; colnames: "))
    print(colnames(meta))

    #files get "random" new filename in /tmp/ directory when uploaded in docker -> change filename to the upload-name
    dirname <- dirname(input$fastqFiles$datapath[1])  # this is the file-path of the fastq files
    file.rename(from=input$fastqFiles$datapath,to=paste0(dirname,"/",input$fastqFiles$name))
    
    # trim primers of specified region(s)
    waiter_update(html = tagList(spin_rotating_plane(), "Trimming Primers ..."))
    if(input$trim_primers=="V3/V4"){
      print(paste0(Sys.time()," - Trimming V3/V4 Primers ... "))
      fw_primer <- "CCTACGGGNGGCWGCAG"
      rv_primer <- "GACTACHVGGGTATCTAATCC"
    }
    trimmed_dir <- paste0(dirname, "/trimmed")
    dir.create(trimmed_dir)
    trim_primers_command <- paste0("bash src/trim_primers.sh ", dirname, 
                                   " ", trimmed_dir,
                                   " ", fw_primer,
                                   " ", rv_primer)
    #out <- system(trim_primers_command, wait=T)
    #dirname <- trimmed_dir
    print(paste0(Sys.time()," - Trimmed Primers. "))
    

    # remove spikes with python script 
    if(rm_spikes){
      print("############ spike removal ############ \n")
      waiter_update(html = tagList(spin_rotating_plane(),"Removing spikes ..."))
      rm_spikes_outdir <- paste0(dirname,"/rm_spikes_out")
      rm_spikes_outdir_woSpikes <- paste0(rm_spikes_outdir,"/omapping")
      rm_spikes_outdir_wSpikes <- paste0(rm_spikes_outdir,"/osamples")
      rm_spikes_spikes_file <- paste0(rm_spikes_outdir,"/ospikes.txt")
      rm_spikes_stats_file <- paste0(rm_spikes_outdir,"/ostats.txt")
      dir.create(rm_spikes_outdir)
      rm_spikes_command <- paste0("python3 src/rm_spikes.py ../data/spikes.fasta ",
                                 input$metaFile$datapath,
                                 " ", dirname,
                                 " ", rm_spikes_outdir_woSpikes,
                                 " ", rm_spikes_outdir_wSpikes,
                                 " ", rm_spikes_spikes_file,
                                 " ", rm_spikes_stats_file,
                                 " ", ncores)
      out <- system(rm_spikes_command, wait = T)
      print(paste0(Sys.time()," - Removed spikes: ",rm_spikes_command))
      print ("############ spike removal finished ############ \n")
      dirname <- rm_spikes_outdir_woSpikes
    }
    
    # collect fw & rv files 
    foreward_files <- sort(list.files(dirname, pattern = "_R1_001.fastq", full.names = T))
    reverse_files <- sort(list.files(dirname, pattern = "_R2_001.fastq", full.names = T))
    sample_names <- sapply(strsplit(basename(foreward_files), "_"), `[`, 1) # sample name: everything until first "_"
    if (length(foreward_files) != length(reverse_files)){stop(noEqualFastqPairsError, call.=F)}
    if (!all(meta[[sample_column]] == sample_names)){stop(noEqualFastqPairsError, call.=F)} # check if all sample names are equal for mapping and fastq
    
    ##### starting DADA2 pipeline #####
    # http://benjjneb.github.io/dada2/tutorial.html
    
    # trim&trunc files
    foreward_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
    reverse_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
    names(foreward_files_filtered) <- sample_names
    names(reverse_files_filtered) <- sample_names
    waiter_update(html = tagList(spin_rotating_plane(),"Filtering ..."))
    out <- filterAndTrim(foreward_files, foreward_files_filtered, reverse_files, reverse_files_filtered, truncLen=c(trunc_fw,trunc_rv),
                         rm.phix=TRUE, compress=TRUE, multithread=ncores) # maxEE filter not used
    print(paste0(Sys.time()," - Filtered fastqs: "))
    print(c(trunc_fw, trunc_rv))
    # remove files, which have 0 reads after filtering

    #TODO: adapt meta-file header to have '#' at beginning!

    # learn errors
    waiter_update(html = tagList(spin_rotating_plane(),"Learning Errors ..."))
    errF <- learnErrors(foreward_files_filtered, multithread=ncores, nbases = 1e8, randomize = T)
    errR <- learnErrors(reverse_files_filtered, multithread=ncores, nbases = 1e8, randomize = T)
    print(paste0(Sys.time()," - Learned Errors. "))

    #dada2
    waiter_update(html = tagList(spin_rotating_plane(),"Sample inference ..."))
    dadaFs <- dada(foreward_files_filtered, err=errF, multithread=ncores)
    dadaRs <- dada(reverse_files_filtered, err=errR, multithread=ncores)
    dada_merged <- mergePairs(dadaFs, foreward_files_filtered, dadaRs, reverse_files_filtered)
    print(paste0(Sys.time()," - Merged files. "))

    # create ASV table & removing chimeras
    waiter_update(html = tagList(spin_rotating_plane(),"Merging and removing chimeras ..."))
    seq_table <- makeSequenceTable(dada_merged)
    seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=ncores)
    print(paste0(Sys.time()," - Created ASV table: "))
    print(dim(seq_table_nochim))

    ##### done with DADA2 pipeline #####

    # calculate loss of reads during steps
    getN <- function(x) sum(getUniques(x))
    if (length(sample_names) > 1){
      track <- cbind(sample_names, out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(dada_merged, getN), rowSums(seq_table_nochim))
    }else{
      track <- cbind(sample_names, out, getN(dadaFs), getN(dadaRs), getN(dada_merged), rowSums(seq_table_nochim))
    }
    rownames(track) <- sample_names
    colnames(track) <- c("sample","input_reads", "filtered", "denoisedF", "denoisedR", "merged", "non_chimera")

    # assign taxonomy
    waiter_update(html = tagList(spin_rotating_plane(),"Assigning taxonomy ..."))
    taxa <- assignTaxonomy(seq_table_nochim, "data/taxonomy_annotation.fa.gz", multithread = ncores)
    print(paste0(Sys.time()," - Assigned Taxonomy. "))

    # build phylogenetic tree
    # https://f1000research.com/articles/5-1492/v1
    seqs <- getSequences(seq_table_nochim)

    # combine results into phyloseq object
    waiter_update(html = tagList(spin_rotating_plane(),"Combining results & Normalizing ..."))
    phylo_unnormalized <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows = F),
                                   sample_data(meta),
                                   tax_table(as.matrix(taxa)))

    dna <- Biostrings::DNAStringSet(taxa_names(phylo_unnormalized))
    names(dna) <- taxa_names(phylo_unnormalized)
    phylo_unnormalized <- merge_phyloseq(phylo_unnormalized, dna)
    taxa_names(phylo_unnormalized) <- paste0("ASV", seq(ntaxa(phylo_unnormalized)))

    # create final objects with "real" ASV names
    asv_table_final <- t(as.data.frame(otu_table(phylo_unnormalized, F)))
    taxonomy_final <- as.data.frame(tax_table(phylo_unnormalized))
    meta_final <- as.data.frame(sample_data(phylo_unnormalized))
    refseq_final <- refseq(phylo_unnormalized)
    
    # create fasta file with sequences of ASVs
    asv_fasta <- paste0(dirname, "/asv_sequences.fasta")
    Biostrings::writeXStringSet(refseq_final, asv_fasta)

    # normalization
    normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
    normalized_asv <- normalizeOTUTable(asv_table_final, normMethod)

    # store all filepaths in one place
    file_df <- data.frame(fw_files = foreward_files, fw_files_filtered= foreward_files_filtered,
                          rv_files = reverse_files, rv_files_filtered = reverse_files_filtered,
                          sample_names = sample_names, asv_fasta = asv_fasta)

    phylo <- phyloseq(otu_table(normalized_asv$norm_tab, taxa_are_rows = T),
                      sample_data(meta_final),
                      tax_table(as.matrix(taxonomy_final)),
                      refseq_final)
    
    print(paste0(Sys.time()," - final phyloseq-object: \n"))
    print(phylo)
    print(paste0(Sys.time()," - Finished fastq data upload! \n"))

    vals$datasets[[input$dataName]] <- list(generated_files = file_df,
                                            fastq_dir = dirname,
                                            is_fastq = T,
                                            track = track,
                                            rawData=asv_table_final,
                                            metaData=meta_final,
                                            taxonomy=taxonomy_final,
                                            counts=NULL,
                                            normalizedData=normalized_asv$norm_tab,
                                            relativeData=normalized_asv$rel_tab,
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