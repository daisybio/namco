# Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
uploadFastqModal <- function(failed=F,error_message=NULL) {
  modalDialog(
    title = "UPLOAD fastq files:",
    HTML("<h5>[For detailed information on how the files have to look, check out the <b>Info & Settings</b> tab on the left!]</h5>"),
    hr(),
    div(id="test-div", p(123)),
    h4("Files:"),
    fluidRow(
      column(6,wellPanel(fileInput("fastqFiles","Select all fastq files", multiple = T, accept = c(".fastq", ".fastq.gz")), style="background:#3c8dbc")),
      column(6,wellPanel(fileInput("metaFile","Select Metadata File [optional]"),
                         textInput("metaSampleColumn", "Name of the sample-column:", value="SampleID")))
    ),
    hr(),
    h4("Additional parameters:"),
    fluidRow(
      column(12, wellPanel(
        fluidRow(
          column(4, radioGroupButtons("rm_spikes", "Remove spikes", c("Yes","No"), direction="horizontal")),
          column(4, selectInput("trim_primers", "Trim Primers",choices = c("V3/V4", "NONE")))
        ),
        div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("truncFw", "Truncation foreward:",value=280, min=1, max=500, step=1)),
        div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("truncRv", "Truncation reverse:",value=200, min=1, max=500, step=1)),
        numericInput("abundance_cutoff", "ASVs with abundance over all samples below this value (in %) will be removed:", value=0.25, min=0, max=100, step=0.01),
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

finishedFastqUploadModal <- modalDialog(
  title = "Success! Upload of your dataset is finished.",
  "Check out the fastq-overview tab on the left for your results and downloads.",
  easyClose = T, size="s"
)

# launch upload dialog
observeEvent(input$upload_fastq, {
  showModal(uploadFastqModal())
})

observeEvent(input$upload_fastq_ok, {
  
  message("Starting fastq data upload ... ")
  
  rm_spikes <- ifelse(input$rm_spikes=="Yes", T, F)
  trim_primers <- ifelse(input$trim_primers=="Yes", T, F)
  overlay_text <- ifelse(rm_spikes, "Starting DADA2 & spike removal ...", "Starting DADA2 ...")
  trunc_fw <- as.numeric(input$truncFw)
  trunc_rv <- as.numeric(input$truncRv)
  
  if(input$trim_primers == "V3/V4"){
    trim_primers <- c(17,21)   # "CCTACGGGNGGCWGCAG" & "GACTACHVGGGTATCTAATCC"
  }else if(input$trim_primers == "NONE"){
    trim_primers <- c(0,0)
  }
  
  waiter_show(html = tagList(spin_rotating_plane(),overlay_text),color=overlay_color)
  
  tryCatch({
    if(is.null(input$fastqFiles$datapath)){stop(noFileError, call. = F)}
    ## load meta-file (..or not)##
    m <- handleMetaFastqMode(input$metaFile$datapath, input$metaSampleColumn, rm_spikes)
    meta <- m$meta
    meta_file_path <- m$meta_file_path
    has_meta <- ifelse(is.null(meta), F, T)

    #files get "random" new filename in /tmp/ directory when uploaded in docker -> change filename to the upload-name
    dirname <- dirname(input$fastqFiles$datapath[1])  # this is the file-path of the fastq files
    file.rename(from=input$fastqFiles$datapath,to=paste0(dirname,"/",input$fastqFiles$name))
    
    # remove spikes with python script 
    if(rm_spikes){
      if(!has_meta){stop(rmSpikesNoMetaError, call. = F)}
      waiter_update(html = tagList(spin_rotating_plane(),"Removing spikes ..."))
      dirname <- removeSpikes(dirname, meta_file_path, ncores)
    }
    
    # collect fw & rv files 
    foreward_files <- sort(list.files(dirname, pattern = "_R1_001.fastq", full.names = T))
    reverse_files <- sort(list.files(dirname, pattern = "_R2_001.fastq", full.names = T))
    sample_names <- sapply(strsplit(basename(foreward_files), "_"), `[`, 1) # sample name: everything until first "_"
    if (length(foreward_files) != length(reverse_files)){stop(noEqualFastqPairsError, call.=F)}
    if (has_meta){if (!all(meta[[sample_column]] == sample_names)){stop(noEqualFastqPairsError, call.=F)}} # check if all sample names are equal for mapping and fastq
    
    ##### starting DADA2 pipeline #####
    # http://benjjneb.github.io/dada2/tutorial.html
    # trim&trunc files
    foreward_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
    reverse_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
    names(foreward_files_filtered) <- sample_names
    names(reverse_files_filtered) <- sample_names
    waiter_update(html = tagList(spin_rotating_plane(),"Filtering and trimming primers ..."))
    out <- filterAndTrim(foreward_files, 
                         foreward_files_filtered, 
                         reverse_files, 
                         reverse_files_filtered, 
                         truncLen=c(trunc_fw,trunc_rv), trimLeft = trim_primers, rm.phix=TRUE, compress=TRUE, multithread=ncores, maxEE = c(2,2)) 
    message(paste0(Sys.time()," - Filtered fastqs: ", trunc_fw, " - ", trunc_rv))

    # learn errors
    waiter_update(html = tagList(spin_rotating_plane(),"Learning Errors (foreward)..."))
    errF <- learnErrors(foreward_files_filtered, multithread=ncores, nbases = 1e8, randomize = T)
    waiter_update(html = tagList(spin_rotating_plane(),"Learning Errors (reverse)..."))
    errR <- learnErrors(reverse_files_filtered, multithread=ncores, nbases = 1e8, randomize = T)
    message(paste0(Sys.time()," - Learned Errors. "))

    #dada2
    waiter_update(html = tagList(spin_rotating_plane(),"Sample inference ..."))
    dadaFs <- dada(foreward_files_filtered, err=errF, multithread=ncores)
    dadaRs <- dada(reverse_files_filtered, err=errR, multithread=ncores)
    dada_merged <- mergePairs(dadaFs, foreward_files_filtered, dadaRs, reverse_files_filtered)
    message(paste0(Sys.time()," - Merged files. "))

    # create ASV table & removing chimeras
    waiter_update(html = tagList(spin_rotating_plane(),"Merging and removing chimeras ..."))
    seq_table <- makeSequenceTable(dada_merged)
    seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=ncores)
    message(paste0(Sys.time()," - Created ASV table: ", dim(seq_table_nochim)[1], " - ", dim(seq_table_nochim)[2]))

    ##### done with DADA2 pipeline #####
    
    # calculate loss of reads during steps
    getN <- function(x) sum(getUniques(x))
    if (length(sample_names) > 1){
      track <- cbind(sample_names, out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(dada_merged, getN), rowSums(seq_table_nochim))
    }else{
      track <- cbind(sample_names, out, getN(dadaFs), getN(dadaRs), getN(dada_merged), rowSums(seq_table_nochim))
    }
    rownames(track) <- sample_names
    colnames(track) <- c("sample","input_reads", "filtered & trimmed", "denoisedFW", "denoisedRV", "merged", "non_chimera")

    # assign taxonomy
    waiter_update(html = tagList(spin_rotating_plane(),"Assigning taxonomy ..."))
    taxa <- assignTaxonomy(seq_table_nochim, "data/taxonomy_annotation.fa.gz", multithread = ncores)
    message(paste0(Sys.time()," - Assigned Taxonomy. "))

    # TODO: build phylogenetic tree
    # https://f1000research.com/articles/5-1492/v1
    seqs <- getSequences(seq_table_nochim)

    # combine results into phyloseq object
    waiter_update(html = tagList(spin_rotating_plane(),"Combining results & Normalizing ..."))
    normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
    cn_lst <- combineAndNormalize(seq_table_nochim, taxa, has_meta, meta, sample_names, input$abundance_cutoff, normMethod)
    
    # store all filepaths in one place
    file_df <- data.frame(fw_files = foreward_files, fw_files_filtered= foreward_files_filtered,
                          rv_files = reverse_files, rv_files_filtered = reverse_files_filtered,
                          sample_names = sample_names)

    message(paste0(Sys.time()," - Finished fastq data upload!"))

    vals$datasets[[input$dataName]] <- list(generated_files = file_df,
                                            fastq_dir = dirname,
                                            is_fastq = T,
                                            track = track,
                                            rawData=cn_lst$raw_asv,
                                            metaData=cn_lst$meta,
                                            taxonomy=cn_lst$taxonomy,
                                            counts=NULL,
                                            normalizedData=cn_lst$normalized_asv$norm_tab,
                                            relativeData=cn_lst$normalized_asv$rel_tab,
                                            tree=NULL,
                                            phylo=cn_lst$phylo,
                                            unifrac_dist=NULL,
                                            undersampled_removed=F,
                                            filtered=F,
                                            normMethod = normMethod,
                                            has_meta=has_meta)
    
    updateTabItems(session,"sidebar")
    removeModal()
    
    showModal(finishedFastqUploadModal)
    
    waiter_hide()
  },error=function(e){
    waiter_hide()
    print(e)
    showModal(uploadFastqModal(failed=T,error_message = e))
  })
})