fastqSampleNamesReact <- reactive({
  if(!is.null(input$fastqFiles$datapath)){
    dirname <- dirname(input$fastqFiles$datapath[1])
    file.rename(from=input$fastqFiles$datapath,to=paste0(dirname,"/",input$fastqFiles$name))
    foreward_files <- sort(list.files(dirname, pattern = "_R1_001.fastq", full.names = T))
    reverse_files <- sort(list.files(dirname, pattern = "_R2_001.fastq", full.names = T))
    sample_names <- sapply(strsplit(basename(foreward_files), "_"), `[`, 1) # sample name: everything until first "_"
    if(length(sample_names)==0){
      updateSelectInput(session, "qualityUploadSelectSample", choices = c("Could not extract sample names"))
      return (NULL)
    }
    updateSelectInput(session, "qualityUploadSelectSample", choices = sample_names)
    sample_df <- data.frame(sample_names=sample_names, fw_files=foreward_files, rv_files=reverse_files)
    return(list(sample_df=sample_df))
  }else{
    NULL
  }
})

output$fastq_file_quality_fw_raw <- renderPlot({
  if(!is.null(fastqSampleNamesReact())){
    sample_df <- fastqSampleNamesReact()$sample_df
    if(dim(sample_df)[2]>2){
      selected_sample <- input$qualityUploadSelectSample
      fw_file <- sample_df[["fw_files"]][sample_df[["sample_names"]]==selected_sample]
      p<-plotQualityProfile(fw_file)
      p + geom_vline(xintercept = as.numeric(input$truncFw), color = "red")
    }
  }
})

output$fastq_file_quality_rv_raw <- renderPlot({
  if(!is.null(fastqSampleNamesReact())){
    sample_df <- fastqSampleNamesReact()$sample_df
    if(dim(sample_df)[2]>2){
      selected_sample <- input$qualityUploadSelectSample
      rv_file <- sample_df[["rv_files"]][sample_df[["sample_names"]]==selected_sample]
      p<-plotQualityProfile(rv_file)
      p + geom_vline(xintercept = as.numeric(input$truncRv), color = "red")
    }
  }
})


output$fastq_file_quality_fw_filtered <- renderPlot({
  if(!is.null(currentSet()) ){
    if(vals$datasets[[currentSet()]]$is_fastq && !vals$dataset[[currentSet()]]$is_restored){
      files <- vals$datasets[[currentSet()]]$generated_files
      fastq_pair = input$fastq_file_select_filtered
      fw_file <- files[["fw_files_filtered"]][files[["sample_names"]]==fastq_pair]
      if(is.na(fw_file)){return(NULL)}
      
      p<-plotQualityProfile(fw_file)
      p + geom_vline(xintercept = as.numeric(input$truncFw), color = "red")
      p
    }
  }
})

output$fastq_file_quality_rv_filtered <- renderPlot({
  if(!is.null(currentSet()) ){
    if(vals$datasets[[currentSet()]]$is_fastq && !vals$dataset[[currentSet()]]$is_restored){
      files <- vals$datasets[[currentSet()]]$generated_files
      fastq_pair = input$fastq_file_select_filtered
      rv_file <- files[["rv_files_filtered"]][files[["sample_names"]]==fastq_pair]
      if(is.na(rv_file)){return(NULL)}
      
      p<-plotQualityProfile(rv_file)
      p + geom_vline(xintercept = as.numeric(input$truncRv), color = "red")
      p
    }
  }
})

output$fastq_pipeline_readloss <- renderPlotly({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$is_fastq){
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
