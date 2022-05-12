finishedFastqUploadModal <- modalDialog(
  title = p("Success! Upload of your dataset is finished.", style="color:green; font-size:40px"),
  "Check out the fastq-overview tab on the left for your results and downloads. We advise you to download a namco_session file using the 'Save session' button in the side-menu; this allows you to save the generated ASV table and work with it at a later point again.",
  easyClose = T, size="l"
)

observeEvent(input$loadFastqc,{
  message("Generating FastQC files ...")
  is_paired <- input$fastqIsPaired
  waiter_show(html = tagList(spin_rotating_plane(),"Generating FastQC plots ..."),color=overlay_color)
  
  tryCatch({
    #files get "random" new filename in /tmp/ directory when uploaded in docker -> change filename to the upload-name
    dirname <- dirname(input$fastqFiles$datapath[1])  # this is the file-path of the fastq files
    file.rename(from=input$fastqFiles$datapath,to=paste0(dirname,"/",input$fastqFiles$name))
    
    #check file-type: if compressed file or multiple fastq-files
    outcome_decompress <- decompress(dirname, is_paired)
    if(outcome_decompress == 1){stop(errorDuringDecompression, call. =F)}
    
    # collect fw & rv files (this is only to check for current fastq-pairs)
    foreward_files <- sort(list.files(dirname, pattern = "_R1_001.fastq", full.names = T))
    reverse_files <- sort(list.files(dirname, pattern = "_R2_001.fastq", full.names = T))
    sample_names <- sapply(strsplit(basename(foreward_files), input$sampleNameCutoff), `[`, 1) # sample name: everything until first "_L001" (default)
    if (is_paired && (length(foreward_files) != length(reverse_files))){stop(noEqualFastqPairsError, call.=F)}
    
    # create new folder for fastqc results 
    fastqc_dir <- paste0(dirname(dirname(input$fastqFiles$datapath[1])),"/fastqc_out")
    unlink(fastqc_dir, recursive = T)
    suppressMessages(fastqc(fq.dir = dirname, qc.dir = fastqc_dir, threads = ncores, fastqc.path = fastqc.path))
    
    updateSelectInput(session, "qualityUploadSelectSample", choices = sample_names)
    updateBox('fastqcBox', action = 'toggle')
    
    waiter_hide()
    
  }, error=function(e){
    waiter_hide()
    print(e$message)
    showModal(errorModal(error_message = e$message))
  })
  
})

# this reactive looks for the fastQC output files after a pipeline has run
fastqSampleNamesReact <- reactive({
  if(!is.null(input$fastqFiles$datapath)){
    
    fastq_dir <- dirname(input$fastqFiles$datapath[1])
    fastqc_dir <- paste0(dirname(dirname(input$fastqFiles$datapath[1])),"/fastqc_out")
    
    #### get sample names for selectInput: qualityUploadSelectSample ####
    
    # only need to rename files if they have not been already (by fastqc button for example)
    if(!any(file.exists(paste0(fastq_dir,"/",input$fastqFiles$name)))){
      file.rename(from=input$fastqFiles$datapath,to=paste0(fastq_dir,"/",input$fastqFiles$name)) 
    }
    
    # forward files are enough
    forward_files <- sort(list.files(fastq_dir, pattern = "_R1_001.fastq", full.names = T))
    
    #get sample names
    sample_names <- sapply(strsplit(basename(forward_files), input$sampleNameCutoff), `[`, 1)
    if(length(sample_names)==0){
      showModal(errorModal('Could not extract sample names.'))
      return (NULL)
    }

    #### gather fastQC output files ####
    
    # check if fastQC result dir exists and the files in there are in line with the fastq files
    if(dir.exists(fastqc_dir)){
      fastqc_fw <- list.files(fastqc_dir, pattern="_R1_001_fastqc.zip", full.names = T)
      fastqc_rv <- list.files(fastqc_dir, pattern="_R2_001_fastqc.zip", full.names = T)
      
      if(length(fastqc_fw) == 0 | length(fastqc_rv) == 0){
        showModal(errorModal('No quality control files were found. Please re-run.'))
        return(NULL)
      }
      fastqc_df <- data.frame(sample_names=sample_names, fastqc_fw=fastqc_fw, fastqc_rv=fastqc_rv)
    }
    
    return(list(fastqc_df = fastqc_df))
    
  }else{
    NULL
  }
})

output$fastq_file_quality_fw_pre <- renderPlot({
  if(!is.null(fastqSampleNamesReact())){
    sample_df <- fastqSampleNamesReact()$fastqc_df
    if(dim(sample_df)[2]==3){
      selected_sample <- input$qualityUploadSelectSample
      if(selected_sample != "Press orange button first to generate quality profiles ..."){
        fw_file <- sample_df[["fastqc_fw"]][sample_df[["sample_names"]]==selected_sample]
        p<-qc_plot(fw_file,modules = "Per base sequence quality")
        p
      }
    }
  }
})

output$fastq_file_quality_rv_pre <- renderPlot({
  if(!is.null(fastqSampleNamesReact())){
    sample_df <- fastqSampleNamesReact()$fastqc_df
    if(dim(sample_df)[2]==3){
      selected_sample <- input$qualityUploadSelectSample
      if(selected_sample != "Press orange button first to generate quality profiles ..."){
        rv_file <- sample_df[["fastqc_rv"]][sample_df[["sample_names"]]==selected_sample]
        p<-qc_plot(rv_file,modules = "Per base sequence quality")
        p
      }
    }
  }
})


output$fastq_file_quality_fw_post <- renderPlot({
  if(!is.null(currentSet()) ){
    if(vals$datasets[[currentSet()]]$is_fastq && !vals$datasets[[currentSet()]]$is_restored){
      files <- vals$datasets[[currentSet()]]$generated_files$file_df
      fastq_pair = input$fastq_file_select_post
      fw_file <- files[["fastqc_fw"]][files[["sample_names"]]==fastq_pair]
      if(is.na(fw_file)){return(NULL)}
      
      p<-qc_plot(fw_file,modules = "Per base sequence quality")
      p
    }
  }
})

output$fastq_file_quality_rv_post <- renderPlot({
  if(!is.null(currentSet()) ){
    if(vals$datasets[[currentSet()]]$is_fastq && !vals$datasets[[currentSet()]]$is_restored){
      files <- vals$datasets[[currentSet()]]$generated_files$file_df
      fastq_pair = input$fastq_file_select_post
      rv_file <- files[["fastqc_rv"]][files[["sample_names"]]==fastq_pair]
      if(is.na(rv_file)){return(NULL)}
      
      p<-qc_plot(rv_file,modules = "Per base sequence quality")
      p
    }
  }
})

output$fastq_pipeline_readloss <- renderPlotly({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$is_fastq && !is.null(vals$datasets[[currentSet()]]$track)){
      track <- as.data.table(vals$datasets[[currentSet()]]$track)
      track_long <- gather(track, pipeline_step, n_reads, input_reads:non_chimera, factor_key = T)
      track_long[["n_reads"]] <- as.numeric(track_long[["n_reads"]])
      order_points <- c("sample","input_reads", "filtered", "denoisedF", "denoisedR", "merged", "non_chimera")
      plot_ly(track_long,
              type="scatter",
              x=~pipeline_step,
              y=~n_reads,
              color=~sample,
              mode="lines+markers")
    }
  }else{
    plotly_empty()
  }
})

#### show log file contents ####

# https://stackoverflow.com/questions/45395241/display-content-of-file-in-dashboard-shiny
output$demulti_log_output <- renderUI({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$is_fastq){
      if(!is.null(vals$datasets[[currentSet()]]$generated_files$log_files)){
        f <- vals$datasets[[currentSet()]]$generated_files$log_files$demulti.log
        rawText <- readLines(f)
        splitText <- stringi::stri_split(str = rawText, regex = '\\n')
        replacedText <- lapply(splitText, p)
        return(replacedText)
      }
    }
  }
})

output$run_log_output <- renderUI({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$is_fastq){
      if(!is.null(vals$datasets[[currentSet()]]$generated_files$log_files)){
        f <- vals$datasets[[currentSet()]]$generated_files$log_files$run.log
        rawText <- readLines(f)
        splitText <- stringi::stri_split(str = rawText, regex = '\\n')
        replacedText <- lapply(splitText, p)
        return(replacedText)
      }
    }
  }
})

#### download results ####

output$download_asv_norm <- downloadHandler(
  filename=function(){paste("asv_table_norm.tab")},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$is_fastq){
        write.table(vals$datasets[[currentSet()]]$normalizedData, file, row.names = T, sep="\t", quote=F)
      }
    }
  }
)

output$download_asv_raw <- downloadHandler(
  filename=function(){paste("asv_table_raw.tab")},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$is_fastq){
        write.table(vals$datasets[[currentSet()]]$rawData, file, row.names = T, sep="\t", quote=F)
      }
    }
  }
)

output$download_taxonomy <- downloadHandler(
  filename=function(){paste("taxonomy.tab")},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$is_fastq){
        write.table(as.data.frame(vals$datasets[[currentSet()]]$taxonomy), file, row.names = T, sep="\t", quote=F)
      }
    }
  }
)

output$download_phyloseq <- downloadHandler(
  filename=function(){paste("phyloseq.RDS")},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$is_fastq){
        saveRDS(vals$datasets[[currentSet()]]$phylo, file)
      }
    }
  }
)

output$download_asv_fastq <- downloadHandler(
  filename=function(){paste("asv_sequences.fasta")},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$is_fastq){
        dna <- refseq(vals$datasets[[currentSet()]]$phylo)
        writeXStringSet(dna, file)
      }
    }
  }
)
