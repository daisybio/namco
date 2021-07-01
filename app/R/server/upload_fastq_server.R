#https://ofstack.com/Nginx/17072/nginx-upload-large-file-timeout-solution.html

# Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
uploadFastqModal <- function(failed=F,error_message=NULL) {
  modalDialog(
    title = "UPLOAD fastq files:",
    HTML("<h5>[For detailed information on how the files have to look, check out the <b>Info & Settings</b> tab on the left!]</h5>"),
    hr(),
    h4("Files:"),
    fluidRow(
      column(6,wellPanel(fileInput("fastqFiles","Select fastq-files or compressed folder", multiple = T, accept = c(".fastq", ".fastq.gz", ".tar", ".tar.gz", ".zip")), style="background:#3c8dbc")),
      column(6,wellPanel(fileInput("fastqMetaFile","Select Metadata File [optional]"),
                         textInput("metaSampleColumn", "Name of the sample-column:", value="SampleID")))
    ),
    hr(),
    fluidRow(
      column(1),
      column(4, actionButton("loadFastqc","Generate read quality profiles")),
      column(7, p("It is highly advised to first check the sequencing quality of your reads in order to set the filterin parameters below correctly. Please hit this button to generate quality plots for each file."))
    ),
    hidden(div(id="readQualityRaw",
        hr(),
        fluidRow(
          column(6, selectInput("qualityUploadSelectSample","Select Sample", choices=c("Waiting to finish file upload...")))
        ),
        fluidRow(
          tabBox(
            title="Sequencing quality of uploaded fastq-files",
            id="qualityUploadTabBox", width=12,
            tabPanel("Foreward",
                     fluidRow(column(12, plotOutput("fastq_file_quality_fw_raw")))
            ),
            tabPanel("Reverse",
                     fluidRow(column(12, plotOutput("fastq_file_quality_rv_raw")))
            )
          )
        )
    )),
    #fluidRow(column(10, htmlOutput("fastqQualityTextCopy"))),
    #hr(),
    h4("Additional parameters:"),
    fluidRow(
      column(12, wellPanel(
        fluidRow(
          column(4, radioGroupButtons("rm_spikes", "Remove spikes [needs meta file!]", c("Yes","No"), direction="horizontal", selected = "No")),
          column(4, selectInput("trim_primers", "Trim Primers",choices = c("V3/V4", "NONE")))
        ),
        div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("truncFw", "Truncation foreward:",value=280, min=1, max=500, step=1)),
        div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("truncRv", "Truncation reverse:",value=200, min=1, max=500, step=1)),
        #div(style="display: inline-block;vertical-align:top; width: 150px;",p("These two cutoff values are displayed as a vertical red line in the plots above.")),
        #numericInput("fastq_abundance_cutoff", "ASVs with abundance over all samples below this value (in %) will be removed:", value=0.25, min=0, max=100, step=0.01),
        radioGroupButtons("buildPhyloTree", "build phylogenetic tree [will increase runtime!]", c("Yes", "No"), direction = "horizontal", selected = "No")
      ))
    ),
    hr(),
    fluidRow(
      column(10, box(
               title="Parameter information",
               htmlOutput("dada2_filter_info"),
               solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
      ))
    ),
    br(),
    textInput("dataName","Enter a project name:",placeholder=paste0("Namco_project_",Sys.Date()),value=paste0("Namco_project_",Sys.Date())),
    if(failed) {
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

#TODO
# if cancel is pressed, reset path of meta file

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
    ## load meta-file (..or not)
    m <- handleMetaFastqMode(input$fastqMetaFile$datapath, input$metaSampleColumn, rm_spikes)
    meta <- m$meta
    meta_file_path <- m$meta_file_path
    has_meta <- ifelse(is.null(input$fastqMetaFile$datapath), F, T)
    
    #files get "random" new filename in /tmp/ directory when uploaded in docker -> change filename to the upload-name
    dirname <- dirname(input$fastqFiles$datapath[1])  # this is the file-path of the fastq files
    file.rename(from=input$fastqFiles$datapath,to=paste0(dirname,"/",input$fastqFiles$name))
    
    #check file-type: if compressed file or multiple fastq-files
    waiter_update(html = tagList(spin_rotating_plane(),"Reading in files ..."))
    outcome_decompress <- decompress(dirname)
    if(outcome_decompress == 1){stop(errorDuringDecompression, call. =F)}
    
    # remove spikes with python script 
    if(rm_spikes){
      if(!has_meta){stop(rmSpikesNoMetaError, call. = F)}
      waiter_update(html = tagList(spin_rotating_plane(),"Removing spikes ..."))
      dirname <- removeSpikes(dirname, meta_file_path, ncores)
    }
    
    # collect fw & rv files 
    foreward_files <- sort(list.files(dirname, pattern = "_R1_001.fastq", full.names = T))
    reverse_files <- sort(list.files(dirname, pattern = "_R2_001.fastq", full.names = T))
    sample_names <- sapply(strsplit(basename(foreward_files), "_L001"), `[`, 1) # sample name: everything until first "_L001"
    if (length(foreward_files) != length(reverse_files)){stop(noEqualFastqPairsError, call.=F)}
    if (has_meta){if (!all(meta[[sample_column]] == sample_names)){stop(noEqualFastqPairsError, call.=F)}} # check if all sample names are equal for mapping and fastq
    
    ##### starting DADA2 pipeline #####
    
    # trim&trunc files
    foreward_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
    reverse_files_filtered <- file.path(dirname, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
    names(foreward_files_filtered) <- sample_names
    names(reverse_files_filtered) <- sample_names
    waiter_update(html = tagList(spin_rotating_plane(),"Filtering and trimming primers ..."))
    out_filter <- data.frame(filterAndTrim(foreward_files, 
                                           foreward_files_filtered, 
                                           reverse_files, 
                                           reverse_files_filtered, 
                                           truncLen=c(trunc_fw,trunc_rv), 
                                           trimLeft = trim_primers, 
                                           rm.phix=TRUE, 
                                           compress=TRUE, 
                                           multithread=TRUE, 
                                           maxEE = c(2,2)))
    files_filtered <- rownames(out_filter[out_filter$reads.out!=0,])        # get files(R1), which have more than 0 reads left after filtering
    samples_filtered <- sapply(strsplit(files_filtered, "_"), `[`, 1)       # get all samples, which have more than 0 reads left
    if(length(samples_filtered)==0){stop(noTaxaRemainingAfterFilterError, call.=F)}
    foreward_files_filtered <- foreward_files_filtered[samples_filtered]
    reverse_files_filtered <- reverse_files_filtered[samples_filtered]
    message(paste0(Sys.time()," - Filtered fastqs: ", trunc_fw, " - ", trunc_rv))
    message(paste0(Sys.time(), " - Files with 0 reads after filtering: ", rownames(out_filter[out_filter$reads.out==0,])))
    
    # learn errors
    waiter_update(html = tagList(spin_rotating_plane(),"Learning Errors (foreward)..."))
    errF <- learnErrors(foreward_files_filtered, multithread=TRUE, nbases = 1e8, randomize = T)
    waiter_update(html = tagList(spin_rotating_plane(),"Learning Errors (reverse)..."))
    errR <- learnErrors(reverse_files_filtered, multithread=TRUE, nbases = 1e8, randomize = T)
    message(paste0(Sys.time()," - Learned Errors. "))
    
    #dada2
    waiter_update(html = tagList(spin_rotating_plane(),"Sample inference ..."))
    dadaFs <- dada(foreward_files_filtered, err=errF, multithread=T)
    dadaRs <- dada(reverse_files_filtered, err=errR, multithread=T)
    dada_merged <- mergePairs(dadaFs, foreward_files_filtered, dadaRs, reverse_files_filtered)
    message(paste0(Sys.time()," - Merged files. "))
    
    # create ASV table & removing chimeras
    waiter_update(html = tagList(spin_rotating_plane(),"Merging and removing chimeras ..."))
    seq_table <- makeSequenceTable(dada_merged)
    seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=T)
    message(paste0(Sys.time()," - Created ASV table: ", dim(seq_table_nochim)[1], " - ", dim(seq_table_nochim)[2]))
    
    ##### done with DADA2 pipeline #####
    
    # calculate loss of reads during steps
    track <- calcReadLoss(out_filter, dadaFs, dadaRs, dada_merged, seq_table_nochim, sample_names, samples_filtered)
    
    # assign taxonomy
    waiter_update(html = tagList(spin_rotating_plane(),"Assigning taxonomy ..."))
    taxa <- assignTaxonomy(seq_table_nochim, "data/taxonomy_annotation.fa.gz", multithread = T)
    message(paste0(Sys.time()," - Assigned Taxonomy. "))
    
    # build phylogenetic tree
    seqs <- getSequences(seq_table_nochim)
    if(input$buildPhyloTree=="Yes"){
      waiter_update(html = tagList(spin_rotating_plane(),"building phylogenetic tree ..."))
      tree<-buildPhyloTree(seqs, ncores)
    } else{tree<-NULL; unifrac_dist<-NULL}

    # combine results into phyloseq object
    waiter_update(html = tagList(spin_rotating_plane(),"Combining results & Normalizing ..."))
    cn_lst <- combineAndNormalize(seq_table_nochim, taxa, has_meta, meta, tree, samples_filtered, 0)
    
    if(!is.null(tree)){
      #pre-build unifrac distance matrix
      unifrac_dist <- buildGUniFracMatrix(otu_table(cn_lst$phylo), phy_tree(cn_lst$phylo))
    }else{unifrac_dist<-NULL}
    
    # run FastQC for trimmed & filtered files
    fastqc_dir <- paste0(dirname,"/fastqc_out")
    unlink(fastqc_dir)
    suppressMessages(fastqc(fq.dir = dirname(foreward_files_filtered)[1], qc.dir = fastqc_dir, threads = ncores, fastqc.path = "/opt/FastQC/fastqc"))
    fastqc_fw <- list.files(fastqc_dir, pattern="F_filt_fastqc.zip", full.names = T)
    fastqc_rv <- list.files(fastqc_dir, pattern="R_filt_fastqc.zip", full.names = T)
    
    # store all filepaths in one place
    raw_df <- data.frame(fw_files = foreward_files,
                         rv_files = reverse_files,
                         sample_names = sample_names)
    filtered_df <- data.frame(fw_files_filtered = foreward_files_filtered,
                              rv_files_filtered = reverse_files_filtered,
                              fastqc_fw = fastqc_fw,
                              fastqc_rv = fastqc_rv,
                              sample_names = samples_filtered)
    file_df <- merge(raw_df, filtered_df, all.x = T)
    
    message(paste0(Sys.time()," - Finished fastq data upload!"))
    
    vals$datasets[[input$dataName]] <- list(session_name=input$dataName,
                                            generated_files = file_df,
                                            fastq_dir = dirname,
                                            is_fastq = T,
                                            track = track,
                                            rawData=cn_lst$raw_asv,
                                            metaData=cn_lst$meta,
                                            taxonomy=cn_lst$taxonomy,
                                            counts=NULL,
                                            normalizedData=cn_lst$normalized_asv$norm_tab,
                                            relativeData=cn_lst$normalized_asv$rel_tab,
                                            tree=tree,
                                            phylo=cn_lst$phylo,
                                            unifrac_dist=unifrac_dist,
                                            undersampled_removed=F,
                                            filtered=F,
                                            normMethod = 0,
                                            has_meta=has_meta,
                                            has_picrust=F,
                                            is_sample_data=F,
                                            is_restored=F,
                                            has_rf=F,
                                            has_diff_nw=F,
                                            has_tax_nw=F,
                                            has_comp_nw=F,
                                            filterHistory="")
    
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

observeEvent(input$loadFastqc,{
  message("Generating FastQC files ...")
  waiter_show(html = tagList(spin_rotating_plane(),"Generating FastQC plots ..."),color=overlay_color)
  
  tryCatch({
    #files get "random" new filename in /tmp/ directory when uploaded in docker -> change filename to the upload-name
    dirname <- dirname(input$fastqFiles$datapath[1])  # this is the file-path of the fastq files
    file.rename(from=input$fastqFiles$datapath,to=paste0(dirname,"/",input$fastqFiles$name))
    
    #check file-type: if compressed file or multiple fastq-files
    outcome_decompress <- decompress(dirname)
    if(outcome_decompress == 1){stop(errorDuringDecompression, call. =F)}
    
    # collect fw & rv files (this is only to check for corrent fastq-pairs)
    foreward_files <- sort(list.files(dirname, pattern = "_R1_001.fastq", full.names = T))
    reverse_files <- sort(list.files(dirname, pattern = "_R2_001.fastq", full.names = T))
    if (length(foreward_files) != length(reverse_files)){stop(noEqualFastqPairsError, call.=F)}
    
    # create new folder for fastqc results 
    fastqc_dir <- paste0(dirname,"/fastqc_out")
    unlink(fastqc_dir)
    suppressWarnings(fastqc(fq.dir = dirname,qc.dir = fastqc_dir, threads = ncores, fastqc.path = "/opt/FastQC/fastqc"))
    shinyjs::show("readQualityRaw", anim = T)
    
    waiter_hide()
    
  }, error=function(e){
    waiter_hide()
    print(e$message)
    showModal(uploadFastqModal(failed=T,error_message = e$message))
  })

})
