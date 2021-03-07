output$fastq_file_quality_fw <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$is_fastq){
      files <- vals$datasets[[currentSet()]]$fastq_files
      fastq_pair = input$fastq_file_select
      fw_file <- files[["fw_files"]][files[["sample_name"]]==fastq_pair]
      
      plotQualityProfile(fw_file)
    }
  }
})

output$fastq_file_quality_rv <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$is_fastq){
      files <- vals$datasets[[currentSet()]]$fastq_files
      fastq_pair = input$fastq_file_select
      rv_file <- files[["rv_files"]][files[["sample_name"]]==fastq_pair]
      
      plotQualityProfile(rv_file)
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