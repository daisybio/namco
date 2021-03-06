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
      track_long <- gather(track, pipeline_step, n_reads, input:non_chimera, factor_key = T)
      track_long[["n_reads"]] <- as.numeric(track_long[["n_reads"]])
      order_points <- c("sample","input", "filtered", "denoisedF", "denoisedR", "merged", "non_chimera")
      p<-ggplot(track_long, aes(x=pipeline_step, y = n_reads, group=sample, color=sample))+
        geom_line()+
        geom_point()+
        ggtitle("Number of reads after each processing steps in the DADA2 pipeline")
      ggplotly(p)
    }
  }
})