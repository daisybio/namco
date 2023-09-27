output$exampleInfoText <- renderUI({
  exampleInfoText
})

output$welcome <- renderUI({
  welcomeText
})

output$logo <- renderUI({
  logoText
})

output$workflow <- renderUI({
  workflowImage
})

output$biomedLogo <- renderUI({
  biomedLogo
})

output$contactText <- renderUI({
  contactText
})

output$newsText <- renderUI({
  newsText
})

output$documentation <- renderUI({
  documentationText
})

output$startHere <- renderUI({
  startHereText
})

output$authors <- renderUI({
  authorsText
})

output$welcome_ref <- renderUI({
  sourcesText
})

output$filterHistoryText <-renderUI({
  HTML(vals$datasets[[currentSet()]]$filterHistory)
})

output$taxBinningText <- renderUI({
  taxBinningText
})

output$rarefactionInfoText <- renderUI({
  rarefactionInfoText
})

output$dimReductionInfoText <- renderUI({
  dimReductionInfoText
})

output$confoundingInfoText <- renderUI({
  confoundingInfoText
})

output$randomForestText <- renderUI ({
  randomForestText
})

output$glmText <- renderUI ({
  glmTest
})

output$alphaDivText <- renderUI({
  alphaDivText
})

alphaParameterInfo <- renderUI({
  alphaParameterInfo
})

output$alphaDivFormulas <- renderUI({
  alphaDivFormulas
})

output$betaDivText <- renderUI({
  betaDivText
})

output$associationsText <- renderUI({
  associationsText
})

output$associationsSourceText <- renderUI({
  associationsSourceText
})

output$corrText <- renderUI({
  corrText
})

output$topicText <- renderUI({
  topicText
})

output$phyloTreeText <- renderUI({
  phyloTreeText
})

output$sourcePhyloseq <- renderUI({
  phyloseqSourceText
})

output$text1 <- renderUI({
  if(!is.null(currentSet())){
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      HTML(sprintf("Below are the results of a %s topic STM. The ordination of the topics over taxa distribution (left) and the frequencies of
                   the top %s taxa (in terms of saliency) across all topics. By selecting a topic, the relative
                   frequencies of the taxa within that topic are shown in red. The ordination figure can be shown in
                   either 2D or 3D and the ordination method can be adjusted. Lambda adjusts the relevance calculation.
                   Choosing the taxon adjusts the group coloring for the bar plot. Clicking Reset resets the topic selection.",
                   vis_out$K,vis_out$taxa_bar_n))
    }
  }
  
})

output$text2 <- renderUI({
  if(!is.null(currentSet())){
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      themetagenomicsText2
    }
  }
})

output$text3 <- renderUI({
  themetagenomicsText3
})

output$cutoff_title <- renderUI({
  coOcurrenceCutoffTitleText  })

output$cutoff_text <- renderUI({
  coOcurrenceDistrText
})

output$heatmap_text <- renderUI({
  coOcurrenceHeatmapText
})

output$neatmapSourceText <- renderUI({
  neatmapSourceText
})

output$basic_info <- renderUI({
  coOccurrenceInfoText
})

output$basic_additional <- renderUI({
  coOcurrenceCutoffText
})

output$basic_calc_title <- renderUI({
  coOcurrenceCountsTitleText
})

output$basic_calc_additional <- renderUI({
  coOcurrenceCountsText
})

output$basic_network_title <- renderUI({
  coOcurrenceNetworkTitleText
})

output$chosen_network_params <- renderUI({
  if(!is.null(currentSet())){
    if(!is.null(vals$datasets[[currentSet()]]$network_params)){
      params = vals$datasets[[currentSet()]]$network_params
      calc <- ifelse(params$fc,"fold-change","difference")
      HTML(paste0("sample-group: <b>", params$group_column, "</b><br>",
                  "comparing variable: <b>", params$var1, "</b> with <b>", params$var2, "</b><br>",
                  "chosen cutoff: <b>", params$cutoff, "</b><br>",
                  "counts calculated using: <b>", calc, "</b>"))
    }
  }
})

output$input_variables <- renderUI({
  if(!is.null(currentSet())){
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      K <- vis_out$K
      sigma <- vis_out$sigma
      formula<-vis_out$formula
      refs<-paste(vis_out$refs,collapse=", ")
      HTML(paste0("Number of chosen topics (K): <b>",K,"</b><br>",
                  "Value of sigma_prior: <b>",sigma,"</b><br>",
                  "Group from META file: <b>",formula, "</b><br>",
                  "Reference Level in this group: <b>",refs,"</b>"))
    }else{
      HTML(paste0("Number of chosen topics (K):","<br>",
                  "Value of sigma_prior:","<br>",
                  "Group from META file:","<br>",
                  "Reference Level in this group: "))
    }
  }
})

output$advanced_text <- renderUI({
  themetagenomicsTextTitle
})

output$topic_text <- renderUI({
  themetagenomicsTextTopics
})

output$sigma_text <- renderUI({
  themetagenomicsTextSigma
})

# output$spiec_easi_additional <- renderUI({
#   HTML(paste0("<b>1:</b> absolute value of correlations below this threshold are considered zero by the inner SparCC loop."))
# })

output$info_inputdata <- renderUI({
  inputDataText
})

output$info_testdata <- renderUI({
  testdataText
})

output$forest_model_parameters <- renderPrint({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_rf){
      model<-vals$datasets[[currentSet()]]$rf_lst$model
      model$finalModel
    }
  }
})

output$forest_model_variables <- renderPrint({
  if(!is.null(rForestDataReactive())){
    model<-rForestDataReactive()$model
    features<-setdiff(colnames(model$trainingData),".outcome")
    features
  }
})

output$heatmapOrdinationText <- renderUI({
  heatmapOrdinationText
})

output$heatmapText <- renderUI({
  heatmapText
})

output$heatmapSourceText <- renderUI({
  heatmapText2
})

output$timeSeriesText <- renderUI({
  timeSeriesText
})

output$biomehorizonText <- renderUI({
  biomehorizonText
})

output$statTestText <- renderUI({
  statTestText
})

output$networkNodeText <- renderUI({
  networkNodeText
})

output$networkNodeTextCopy <- renderUI({
  networkNodeText
})

output$networkNodeTextCopyCopy <- renderUI({
  networkNodeText
})

output$diffNetworkInfoText <- renderUI({
  diffNetworkInfoText
})

output$diffNetworkSource <- renderUI({
  diffNetworkSource
})

output$diffNetworkSourceCopy <- renderUI({
  diffNetworkSource
})

output$diffNetworkSourceCopyCopy <- renderUI({
  diffNetworkSource
})

output$taxNetworkInfoText <- renderUI({
  taxNetworkInfoText
})

output$diffNetworkParameterText <- renderUI({
  diffTaxNetworkParameterText
})

output$taxNetworkParamsText <- renderUI({
  diffTaxNetworkParameterText
})

output$compNetworkInfoText <- renderUI({
  compNetworkInfoText
})

output$compNetworkParamsText <- renderUI({
  diffTaxNetworkParameterText
})

output$picrust2Text <- renderUI({
  picrust2Text
})

output$picrust2SourceText <- renderUI({
  picrust2SourceText
})

output$picrust_pval_info_text <- renderUI({
  picrust_pval_info_text
})

output$aldexSourceText <- renderUI({
  aldexSourceText
})

output$dada2_filter_info <- renderUI({
  dada2_filter_info
})

output$dada2SourceText <- renderUI({
  dada2SourceText
})

output$fastqQualityText <- renderUI({
  fastqQualityText
})

output$fastqQualityTextCopy <- renderUI({
  fastqQualityText
})

output$decontamText <- renderUI({
  decontamText
})

output$dada2_lotus2_warning <- renderUI({
  HTML('<font color="orange">Warning: Due to differences in pre- and post-processing, the standalone DADA2 (left) pipeline and DADA2 integrated in LotuS2 will likely result in different results.<br></font>')
})