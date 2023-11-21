output$below_abundance_filter_features1 <- renderValueBox({
  otu_below <- 0
  cutoff <- 0.25
  color<-"red"
  if (!is.null(currentSet())) {
    phylo <- vals$datasets[[currentSet()]]$phylo
    rel_otu <- relAbundance(as.data.frame(otu_table(phylo)))
    min <- apply(rel_otu, 2, function(x) ifelse(x>cutoff,1,0))
    keep_taxa = names(which(rowSums(min)>0))
    otu_below <- ntaxa(phylo) - length(keep_taxa)
  }
  if(otu_below == 0){color="green"}
  valueBox(otu_below, "ASVs/OTUs below 0.25% rel. abundance", icon = icon("attention"), color = color)
})

output$below_prevalence_filter_features1 <- renderValueBox({
  otu_below <- 0
  cutoff <- 0.1
  color<-"red"
  if (!is.null(currentSet())) {
    phylo <- vals$datasets[[currentSet()]]$phylo
    rel_otu <- relAbundance(as.data.frame(otu_table(phylo)))
    min <- apply(rel_otu, 2, function(x) ifelse(x>cutoff,1,0))
    keep_taxa <- names(which(rowSums(min)/dim(rel_otu)[2] > cutoff))
    otu_below <- ntaxa(phylo) - length(keep_taxa)
  }
  if(otu_below == 0){color="green"}
  valueBox(otu_below, "ASVs/OTUs below 10% prevalence", icon = icon("attention"), color = color)
})

## rarefaction filtering ##
# check for update if undersampled columns are to be removed (rarefaction curves)
observeEvent(input$excludeSamples, {
  if (!is.null(currentSet())) {
    if (vals$datasets[[currentSet()]]$has_meta) {
      if (!is.null(vals$undersampled) && input$excludeSamples == T) {
        
        vals$datasets[[currentSet()]]$old.dataset.rarefaction <- vals$datasets[[currentSet()]]
        
        maintained_samples <- setdiff(sample_names(vals$datasets[[currentSet()]]$phylo), vals$undersampled)
        phylo.new <- prune_samples(maintained_samples, vals$datasets[[currentSet()]]$phylo) # keep all samples except undersampled ones
        phylo.raw.new <- prune_samples(maintained_samples, vals$datasets[[currentSet()]]$phylo.raw)
        filterMessage <- paste0("<br>",Sys.time()," - removed (undersampled) samples: ", paste(unlist(vals$undersampled), collapse = "; "))
        # add applied filter to history
        vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,filterMessage)
        
        # these are the new samples
        new_samples <- as.vector(sample_names(phylo.new))
        
        # set global variable to TRUE, indicating, that undersampled data is already removed
        vals$datasets[[currentSet()]]$undersampled_removed <- T
        
        if(vals$datasets[[currentSet()]]$has_meta){
          #replace metaData
          vals$datasets[[currentSet()]]$metaData <- data.frame(phylo.new@sam_data, check.names = F)
        }
        
        #adapt otu-tables to only have samples, which were not removed by filter
        vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[,new_samples,drop=F]
        vals$datasets[[currentSet()]]$normalizedData <- vals$datasets[[currentSet()]]$normalizedData[,new_samples,drop=F]
        vals$datasets[[currentSet()]]$relativeData <- vals$datasets[[currentSet()]]$relativeData[,new_samples,drop=F]
        
        # update phyloseq object
        vals$datasets[[currentSet()]]$phylo <- phylo.new
        vals$datasets[[currentSet()]]$phylo.raw <- phylo.raw.new
        message(paste0(Sys.time()," - filtered dataset: "))
        message(nsamples(phylo.new))
        
        #re-calculate unifrac distance
        #pick correct subset of unifrac distance matrix, containing only the new filtered samples
        tree <- phylo.new@phy_tree
        if(!is.null(tree)) unifrac_dist <- as.dist(as.matrix(vals$datasets[[currentSet()]]$unifrac_dist)[new_samples,new_samples]) else unifrac_dist <- NULL
        vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist 
        showModal(infoModal(paste0("Filtering successful. ", nsamples(phylo.new)," samples are remaining.")))
        
        # re-calculate alpha-diversity
        if(vals$datasets[[currentSet()]]$has_meta){
          alphaTabFull <- createAlphaTab(data.frame(phylo.new@otu_table, check.names=F), data.frame(phylo.new@sam_data, check.names = F))
        }else{
          alphaTabFull <- createAlphaTab(data.frame(phylo.new@otu_table, check.names=F))
        }
        vals$datasets[[currentSet()]]$alpha_diversity <- alphaTabFull
        
        # check if picrust results are stored; need to remove samples there as well
        if(!is.null(vals$datasets[[currentSet()]]$picrust_results_list)){
          picrust_results <-  vals$datasets[[currentSet()]]$picrust_results_list
          p2EC <- picrust_results$p2EC[,new_samples]
          p2KO <- picrust_results$p2KO[,new_samples]
          p2PW <- picrust_results$p2PW[,new_samples]
          vals$datasets[[currentSet()]]$picrust_results_list <- list(p2EC=p2EC,
                                                                     p2KO=p2KO,
                                                                     p2PW=p2PW,
                                                                     p2EC_descr=vals$datasets[[currentSet()]]$picrust_results_list$p2EC_descr,
                                                                     p2KO_descr=vals$datasets[[currentSet()]]$picrust_results_list$p2KO_descr,
                                                                     p2PW_descr=vals$datasets[[currentSet()]]$picrust_results_list$p2PW_descr)
          vals$datasets[[currentSet()]]$picrust_analysis_list <- NULL
        }
        
      } else if (input$excludeSamples == F && vals$datasets[[currentSet()]]$undersampled_removed == T) {
        # case: undersampled data was removed but shall be used again (switch turned OFF)
        # use old (oversampled) data again, which was saved
        message("reseting dataset ...")
        vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,"<br>",Sys.time()," - Restoring undersampled samples.")
        filterHistory <- vals$datasets[[currentSet()]]$filterHistory
        restored_dataset <- vals$datasets[[currentSet()]]$old.dataset.rarefaction
        vals$datasets[[currentSet()]] <- restored_dataset
        #remove old dataset from restored dataset
        vals$datasets[[currentSet()]]$old.dataset.rarefaction <- NULL
        vals$datasets[[currentSet()]]$undersampled_removed <- F
        # keep filter history
        vals$datasets[[currentSet()]]$filterHistory <- filterHistory
      }
    }
  }
})


observeEvent(input$filterApplySamples, {
  if(!is.null(currentSet())){
    message("Filtering Samples ...")
    #save "old" dataset to reset filters later; only if there are no sample-filters applied to the current set
    if(!vals$datasets[[currentSet()]]$filtered){
      vals$datasets[[currentSet()]]$old.dataset <- vals$datasets[[currentSet()]]
    }
    
    tryCatch({
      meta_changed = F 
      updateSwitchInput(session, "excludeSamples", value = F)
      
      # filter by specific sample names
      # this works without meta file
      if(!is.null(input$filterSample)){
        maintained_samples <- setdiff(sample_names(vals$datasets[[currentSet()]]$phylo), input$filterSample)
        phylo.new <- prune_samples(maintained_samples, vals$datasets[[currentSet()]]$phylo) # keep all samples except selected ones
        phylo.raw.new <- prune_samples(maintained_samples, vals$datasets[[currentSet()]]$phylo.raw)
        meta_changed = T
        filterMessage <- paste0("<br>",Sys.time()," - removed samples: ", paste(unlist(input$filterSample), collapse = "; "))
        message(filterMessage)
        # add applied filter to history
        vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,filterMessage)
      }
      
      # filter by samples groups
      # this does not work without meta file
      if(input$filterColumns != "NONE" && vals$datasets[[currentSet()]]$has_meta){
        #convert to datatable for filtering to work
        meta <- data.table(data.frame(vals$datasets[[currentSet()]]$phylo@sam_data, keep.rownames = F), check.names = F)
        
        #subset metatable by input 
        meta <- meta[!get(input$filterColumns) %in% input$filterColumnValues,]  # subset meta to have all samples except the ones in selected group
        keep_samples <- meta[[sample_column]]
        phylo.new <- prune_samples(keep_samples, vals$datasets[[currentSet()]]$phylo)
        phylo.raw.new <- prune_samples(keep_samples, vals$datasets[[currentSet()]]$phylo.raw)
        meta_changed = T
        filterMessage <- paste0("<br>",Sys.time()," - removed sample-group: ",input$filterColumns, "==",input$filterColumnValues, collapse = "; ")
        message(filterMessage)
        # add applied filter to history
        vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,filterMessage)
      }
      
      # these are the new samples
      new_samples <- as.vector(sample_names(phylo.new))
      
      if(meta_changed){
        #remember that this dataset is filtered now
        vals$datasets[[currentSet()]]$filtered = T
        
        if(vals$datasets[[currentSet()]]$has_meta){
          #replace metaData
          vals$datasets[[currentSet()]]$metaData <- data.frame(phylo.new@sam_data, check.names = F)
        }
        
        #adapt otu-tables to only have samples, which were not removed by filter
        vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[,new_samples,drop=F]
        vals$datasets[[currentSet()]]$normalizedData <- vals$datasets[[currentSet()]]$normalizedData[,new_samples,drop=F]
        vals$datasets[[currentSet()]]$relativeData <- vals$datasets[[currentSet()]]$relativeData[,new_samples,drop=F]
        
        # update phyloseq object
        vals$datasets[[currentSet()]]$phylo <- phylo.new
        vals$datasets[[currentSet()]]$phylo.raw <- phylo.raw.new
        message(paste0(Sys.time()," - filtered dataset: "))
        message(nsamples(phylo.new))
        
        #re-calculate unifrac distance
        #pick correct subset of unifrac distance matrix, containing only the new filtered samples
        tree <- phylo.new@phy_tree
        if(!is.null(tree)) unifrac_dist <- as.dist(as.matrix(vals$datasets[[currentSet()]]$unifrac_dist)[new_samples,new_samples]) else unifrac_dist <- NULL
        vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist 
        showModal(infoModal(paste0("Filtering successful. ", nsamples(phylo.new)," samples are remaining.")))
        
        # re-calculate alpha-diversity
        if(vals$datasets[[currentSet()]]$has_meta){
          alphaTabFull <- createAlphaTab(data.frame(phylo.new@otu_table, check.names=F), data.frame(phylo.new@sam_data, check.names = F))
        }else{
          alphaTabFull <- createAlphaTab(data.frame(phylo.new@otu_table, check.names=F))
        }
        vals$datasets[[currentSet()]]$alpha_diversity <- alphaTabFull
        
        # check if picrust results are stored; need to remove samples there as well
        if(!is.null(vals$datasets[[currentSet()]]$picrust_results_list)){
          picrust_results <-  vals$datasets[[currentSet()]]$picrust_results_list
          p2EC <- picrust_results$p2EC[,new_samples]
          p2KO <- picrust_results$p2KO[,new_samples]
          p2PW <- picrust_results$p2PW[,new_samples]
          vals$datasets[[currentSet()]]$picrust_results_list <- list(p2EC=p2EC,
                                                                     p2KO=p2KO,
                                                                     p2PW=p2PW,
                                                                     p2EC_descr=vals$datasets[[currentSet()]]$picrust_results_list$p2EC_descr,
                                                                     p2KO_descr=vals$datasets[[currentSet()]]$picrust_results_list$p2KO_descr,
                                                                     p2PW_descr=vals$datasets[[currentSet()]]$picrust_results_list$p2PW_descr)
          vals$datasets[[currentSet()]]$picrust_analysis_list <- NULL
        }
        
      }  
    }, error = function(e){
      vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,"<br>",Sys.time()," - error during sample filtering:", e$message)
      print(e$message)
      showModal(errorModal(e$message))
      return(NULL)
    })
  }
})

observeEvent(input$filterApplyTaxa,{
  if(!is.null(currentSet())){
    
    updateSwitchInput(session, "excludeSamples", value = F) #Why did I do this?
    
    taxonomy <- data.frame(vals$datasets[[currentSet()]]$phylo@tax_table)
    taxonomy <- taxonomy[taxonomy[[input$filterTaxa]] %in% input$filterTaxaValues,]
    remainingOTUs <- rownames(taxonomy)
    filterMessage <- paste0("<br>",Sys.time()," - keeping taxa: ", paste(unlist(input$filterTaxaValues), collapse = "; "))
    
    
    vals$datasets[[currentSet()]] <- filter_taxa_custom(otus_to_keep = remainingOTUs,
                                                        message = filterMessage,
                                                        current_dataset = vals$datasets[[currentSet()]])
  }
})

observeEvent(input$filterResetA, {
  if(!is.null(currentSet())){
    #check if there filters applied to dataset
    if(vals$datasets[[currentSet()]]$filtered){
      message("reseting dataset ...")
      vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,"<br>",Sys.time()," - Original dataset restored.")
      filterHistory <- vals$datasets[[currentSet()]]$filterHistory
      restored_dataset <- vals$datasets[[currentSet()]]$old.dataset
      vals$datasets[[currentSet()]] <- restored_dataset
      #remove old dataset from restored dataset
      vals$datasets[[currentSet()]]$old.dataset <- NULL
      #dataset is not filtered anymore
      vals$datasets[[currentSet()]]$filtered <- F
      # keep filter history
      vals$datasets[[currentSet()]]$filterHistory <- filterHistory
    }
  }
})

observeEvent(input$filterResetB, {
  if(!is.null(currentSet())){
    #check if there filters applied to dataset
    if(vals$datasets[[currentSet()]]$filtered){
      message("reseting dataset ...")
      vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,"<br>",Sys.time()," - Original dataset restored.")
      filterHistory <- vals$datasets[[currentSet()]]$filterHistory
      restored_dataset <- vals$datasets[[currentSet()]]$old.dataset
      vals$datasets[[currentSet()]] <- restored_dataset
      #remove old dataset from restored dataset
      vals$datasets[[currentSet()]]$old.dataset <- NULL
      #dataset is not filtered anymore
      vals$datasets[[currentSet()]]$filtered <- F
      # keep filter history
      vals$datasets[[currentSet()]]$filterHistory <- filterHistory
    }
  }
})

# observer to show tax binning by sample
observe({
  if(!is.null(currentSet())){
    if(!vals$datasets[[currentSet()]]$has_meta){
      shinyjs::show("taxBinningDiv",anim = T)
    }else{
      shinyjs::hide("taxBinningDiv",anim = T)
    }
  }
})

observe({
  if(input$advFilterMinAbundance){enable("advFilterMinAbundanceValue")}else{disable("advFilterMinAbundanceValue")}
  if(input$advFilterRelAbundance){enable("advFilterRelAbundanceValue")}else{disable("advFilterRelAbundanceValue")}
  if(input$advFilterNumSamples){enable("advFilterNumSamplesValue")}else{disable("advFilterNumSamplesValue")}
  if(input$advFilterMaxVariance){enable("advFilterMaxVarianceValue")}else{disable("advFilterMaxVarianceValue")}
  if(input$advFilterPrevalence){enable("advFilterPrevalenceValue")}else{disable("advFilterPrevalenceValue")}
})

observeEvent(input$filterApplyAdv, {
  if(!is.null(currentSet())){
    message("Filtering taxa (advanced) ...")
    #save "old" dataset to reset filters later; only if there are no taxa-filters applied to the current set
    if(!vals$datasets[[currentSet()]]$filtered){
      vals$datasets[[currentSet()]]$old.dataset <- vals$datasets[[currentSet()]]
    }
    
    # store phylo and taxa sums of that phylo
    f_list <- list(phylo = vals$datasets[[currentSet()]]$phylo, 
                   x = taxa_sums(vals$datasets[[currentSet()]]$phylo),
                   otu = as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo)),
                   rel_otu = relAbundance(as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo))))
    keep_taxa = NULL
    filterMessage = ""
    # apply different filtering functions
    if(input$advFilterMinAbundance){
      keep_taxa = names(which(f_list$x > input$advFilterMinAbundanceValue))
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa, f_list)
      filterMessage <- paste(filterMessage,Sys.time()," - Filtered by minimum Abundance:",input$advFilterMinAbundanceValue,"<br>")
    }
    if(input$advFilterRelAbundance){
      cutoff <- input$advFilterRelAbundanceValue
      min <- apply(f_list$rel_otu, 2, function(x) ifelse(x>cutoff,1,0))
      keep_taxa = names(which(rowSums(min)>0))
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa, f_list)
      filterMessage <- paste(filterMessage,Sys.time()," - Filtered by relative abundance:",input$advFilterRelAbundanceValue,"<br>")
    }
    if(input$advFilterNumSamples){
      keep_taxa = rownames(f_list$otu[rowSums(f_list$otu == 0) < input$advFilterNumSamplesValue, ])
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa, f_list)
      filterMessage <- paste(filterMessage,Sys.time()," - Filtered by occurrence in samples:",input$advFilterNumSamplesValue,"<br>")
    }
    if(input$advFilterMaxVariance){
      keep_taxa = names(sort(genefilter::rowVars(f_list$otu), decreasing = T)[1:input$advFilterMaxVarianceValue])
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa, f_list)
      filterMessage <- paste(filterMessage,Sys.time()," - Filtered by maximum variance:",input$advFilterMaxVarianceValue,"<br>")
    }
    if(input$advFilterPrevalence){
      cutoff <- input$advFilterPrevalenceValue/100
      min <- apply(f_list$rel_otu, 2, function(x) ifelse(x>cutoff,1,0))
      keep_taxa <- names(which(rowSums(min)/dim(f_list$rel_otu)[2] > cutoff))
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa, f_list)
      filterMessage <- paste(filterMessage,Sys.time()," - Filtered by prevalence:",input$advFilterPrevalenceValue,"<br>")
    }
    if(is.null(keep_taxa)){
      tmp <- applyFilterFunc(f_list$phylo, keep_taxa)
      return()
    }
    if(is.null(f_list)){
      return()
    }
    if(!is.null(f_list$message)){
      showModal(errorModal(f_list$message))
      return()
    }
    #adapt otu-tables to only have OTUs, which were not removed by filter
    vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[keep_taxa,]
    
    #recalculate the relative abundances and normalize again 
    if(vals$datasets[[currentSet()]]$has_picrust){
      normMethod<-0
    }else{
      normMethod <- vals$datasets[[currentSet()]]$normMethod
    }
    normalizedData <- normalizeOTUTable(f_list$phylo, normMethod)
    vals$datasets[[currentSet()]]$normalizedData <- normalizedData$norm_tab
    vals$datasets[[currentSet()]]$relativeData <- normalizedData$rel_tab
    
    vals$datasets[[currentSet()]]$phylo <- f_list$phylo
    tree <- vals$datasets[[currentSet()]]$tree 
    if(!is.null(tree)){vals$datasets[[currentSet()]]$tree <- phy_tree(f_list$phylo)}
  
    #recalculate unifrac distance in this case
    if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalizedData$norm_tab, vals$datasets[[currentSet()]]$tree) else unifrac_dist <- NULL
    vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist 
    
    vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,"<br>",filterMessage)
    message(paste0(Sys.time()," - filtered dataset: "))
    message(ntaxa(f_list$phylo))
    vals$datasets[[currentSet()]]$filtered = T
  }
})

observeEvent(input$filterResetC, {
  if(!is.null(currentSet())){
    #check if there filters applied to dataset
    if(vals$datasets[[currentSet()]]$filtered){
      message("reseting dataset ...")
      vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,Sys.time()," - Original dataset restored.<br>")
      filterHistory <- vals$datasets[[currentSet()]]$filterHistory
      restored_dataset <- vals$datasets[[currentSet()]]$old.dataset
      vals$datasets[[currentSet()]] <- restored_dataset
      #remove old dataset from restored dataset
      vals$datasets[[currentSet()]]$old.dataset <- NULL
      #dataset is not filtered anymore
      vals$datasets[[currentSet()]]$filtered <- F
      # keep filter history
      vals$datasets[[currentSet()]]$filterHistory <- filterHistory
    }
  }
})

output$advFilterMinAbundancePlot <- renderPlotly({
  if(!is.null(currentSet())){
    abundances <- data.frame(abundance=unlist(taxa_sums(vals$datasets[[currentSet()]]$phylo)))
    p <- plot_ly(x=abundances$abundance, type="histogram", nbinsx=ntaxa(vals$datasets[[currentSet()]]$phylo)/10)
    p %>% layout(shapes=list(vline(input$advFilterMinAbundanceValue)), xaxis=list(title="abundance (summed up over all samples)"),
                 yaxis=list(type="log"))
  }
})

output$advFilterRelAbundancePlot <- renderPlotly({
  if(!is.null(currentSet())){
    rel_abundances <- melt(relAbundance(data.frame(otu_table(vals$datasets[[currentSet()]]$phylo))))
    p <- plot_ly(x=rel_abundances$value, type="histogram", nbinsx=ntaxa(vals$datasets[[currentSet()]]$phylo)/10)
    p %>% layout(shapes=list(vline(input$advFilterRelAbundanceValue)), xaxis=list(title="relative abundances"),
                 yaxis=list(type="log"))
  }
})

output$advFilterNumSamplesPlot <- renderPlotly({
  if(!is.null(currentSet())){
    zero_counts <- data.frame(counts=unlist(rowSums(as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo)) == 0)))
    p <- plot_ly(x=zero_counts$counts, type="histogram")
    p %>% layout(shapes=list(vline(input$advFilterNumSamplesValue)), xaxis=list(title="number of OTUs with abundance of 0"))
  }
})

output$advFilterMaxVariancePlot <- renderPlotly({
  if(!is.null(currentSet())){
    vars <- sort(unlist(genefilter::rowVars(as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo)))), decreasing = T)
    cutoff_var <- vars[input$advFilterMaxVarianceValue]
    variance <- data.frame(var = vars)
    p <- plot_ly(x=variance$var, type="histogram", nbinsx=ntaxa(vals$datasets[[currentSet()]]$phylo)/10)
    p %>% layout(shapes=list(vline(cutoff_var)), xaxis=list(title="variance values of OTUs"),
                 yaxis=list(type="log"))
  }
})


### Decontamination ###

# plot to inspect library size of control and true samples

decontamPlotReactive <- reactive({
  if(!is.null(currentSet())){
    if(input$controlSamplesColumn != 'NULL'){
      phylo <- vals$datasets[[currentSet()]]$phylo
      meta <- as.data.frame(sample_data(phylo))
      meta$LibrarySize <- sample_sums(phylo)
      meta <- meta[order(meta$LibrarySize),]
      meta$Index <- seq(nrow(meta))
      meta$Sample_or_Control <- meta[[input$controlSamplesColumn]]
      
      p <- ggplot(data=meta, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + 
        geom_point()+
        theme_bw()+
        scale_color_manual(name = input$controlSamplesColumn, values = colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(meta$Sample_or_Control))))+
        theme(legend.position = 'top')
      
      return(list(p=p))
    }
  }
})

output$librarySizePlot <- renderPlot({
  if(!is.null(decontamPlotReactive)){
    decontamPlotReactive()$p
  }
})

output$librarySizePlotPDF <- downloadHandler(
  filename = function(){"decontamination_librarySize.pdf"},
  content = function(file){
    if(!is.null(decontamPlotReactive())){
      p <- decontamPlotReactive()$p
      ggsave(filename = file, plot=p, device="pdf", width = 10, height = 10)
    }
  }
)

decontamReactive <- eventReactive(input$startDecontam,{
  if(!is.null(currentSet())){
    
    waiter_show(html = tagList(spin_rotating_plane(),"Decontamination running ... "),color=overlay_color)
    
    tryCatch({
      if(input$DNAconcentrationColumn == 'NULL' && input$controlSamplesColumn == 'NULL'){
        stop(missingParametersDecontamError, call.=F)
      }
      phylo <- vals$datasets[[currentSet()]]$phylo
      if(any(is.na(sample_sums(phylo)))){
        zero_samples <- unlist(names(which(is.na(sample_sums(phylo)) | sample_sums(phylo) == 0)))
        zero_samples <- stringr::str_c(zero_samples, collapse = ', ')
        stop(paste0(sampleWithZeroReadsError, zero_samples), call.=F)
      }
      
      if(input$controlSamplesColumn != 'NULL'){
        sample_data(phylo)$is.neg <- sample_data(phylo)[[input$controlSamplesColumn]] == input$controlSamplesName
      }
      
      if(input$DNAconcentrationColumn != 'NULL'){
        if(any(is.na(sample_data(phylo)[[input$DNAconcentrationColumn]]))){stop(NAinDNAconcentrationError, call.=F)}
        if(any(sample_data(phylo)[[input$DNAconcentrationColumn]] < 0)){stop(negativesinDNAconcentrationError, call.=F)}
      }

      if(input$controlSamplesColumn != 'NULL' && input$DNAconcentrationColumn != 'NULL'){
        contamdf <- decontam::isContaminant(seqtab = phylo,
                                            conc = input$DNAconcentrationColumn,
                                            neg = 'is.neg',
                                            method = 'combined',
                                            batch = NULL, #TODO
                                            threshold = input$decontamThreshold, 
                                            normalize = FALSE)
        used_prev <- TRUE
        used_conc <- TRUE
      }
      else if(input$controlSamplesColumn != 'NULL'){
        contamdf <- decontam::isContaminant(seqtab = phylo,
                                            neg = 'is.neg',
                                            method = 'prevalence',
                                            batch = NULL, #TODO
                                            threshold = input$decontamThreshold, #TODO
                                            normalize = FALSE)
        used_prev <- TRUE
        used_conc <- FALSE
      }else if(input$DNAconcentrationColumn != 'NULL'){
        contamdf <- decontam::isContaminant(seqtab = phylo,
                                            conc = input$DNAconcentrationColumn,
                                            method = 'frequency',
                                            batch = NULL, #TODO
                                            threshold = input$decontamThreshold, #TODO
                                            normalize = FALSE)
        used_prev <- FALSE
        used_conc <- TRUE
      }
      
      contamdf$feature <-rownames(contamdf)
      
      n_contam <- length(which(contamdf$contaminant))
      waiter_hide()
      showModal(infoModal(info_message = paste0('Finished contaminants identification. A total of ', n_contam,' candidates were found. Scroll down to inspect and remove them.')))
      
      return(list(contamdf=contamdf, phylo=phylo, used_prev = used_prev, used_conc = used_conc))
      
    }, error=function(e){
      showModal(errorModal(e$message))
      return(NULL)
    })
  }
})

output$contamTableDownload <- downloadHandler(
  filename=function(){paste("decontamination_results.csv")},
  content = function(file){
    if(!is.null(decontamReactive())){
      write.csv(decontamReactive()$contamdf,file,row.names = F)
    }
  }
)

contamDiagnosticDNAReactive <- reactive({
  if(!is.null(decontamReactive())){
    
    if(!decontamReactive()$used_conc){
      return(NULL)
    }
    phylo <- decontamReactive()$phylo
    contamdf <- decontamReactive()$contamdf
    
    p <- plot_frequency(seqtab = phylo,
                        taxa = taxa_names(phylo)[which(taxa_names(phylo) == input$contamCandidatesSelect)],
                        conc = input$DNAconcentrationColumn)
    p <- p + 
      theme_bw()+
      ggtitle(input$contamCandidatesSelect)
    
    return(list(p=p))
    
  }
})

output$contamDiagnosticDNA <- renderPlot({
  if(!is.null(contamDiagnosticDNAReactive())){
    contamDiagnosticDNAReactive()$p
  }
})

output$contamDiagnosticDNA_PDF <- downloadHandler(
  filename = function(){"decontamination_diagnostics_DNAconcentration.pdf"},
  content = function(file){
    if(!is.null(contamDiagnosticDNAReactive())){
      p <- contamDiagnosticDNAReactive()$p
      ggsave(filename = file, plot=p, device="pdf", width = 10, height = 10)
    }
  }
)

contamDiagnosticPrevReactive <- reactive({
  if(!is.null(decontamReactive())){
    
    if(!decontamReactive()$used_prev){
      return(NULL)
    }
    
    phylo <- decontamReactive()$phylo
    contamdf <- decontamReactive()$contamdf
    
    # https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
    ps.pa <- transform_sample_counts(phylo, function(abund) 1*(abund>0))
    ps.pa.neg <- prune_samples(sample_data(ps.pa)$is.neg == TRUE, ps.pa)
    ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.neg == FALSE, ps.pa)
    
    df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), 
                        pa.neg=taxa_sums(ps.pa.neg),
                        contaminant=contamdf$contaminant)
    
    df.pa$feature <- rownames(df.pa)
    
    p <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant, group=feature)) + 
      geom_point(size = 2) +
      xlab("Prevalence (Control Samples)") + 
      ylab("Prevalence (True Samples)")+
      theme_bw()+
      scale_color_manual(values = colorRampPalette(brewer.pal(9, input$namco_pallete))(2))
    
    return(list(p=p))
  }
})

output$contamDiagnosticPrev <- renderPlotly({
  if(!is.null(contamDiagnosticPrevReactive())){
    ggplotly(contamDiagnosticPrevReactive()$p, tooltip=c('feature','contaminant')) 
  }
})

output$contamDiagnosticPrev_PDF <- downloadHandler(
  filename = function(){"decontamination_diagnostics_prevalence.pdf"},
  content = function(file){
    if(!is.null(contamDiagnosticPrevReactive())){
      p <- contamDiagnosticPrevReactive()$p
      ggsave(filename = file, plot=p, device="pdf", width = 10, height = 10)
    }
  }
)

observeEvent(input$removeContam, {
  if(!is.null(decontamReactive())){
    phylo <- vals$datasets[[currentSet()]]$phylo
    otus_to_keep <- taxa_names(phylo)[!taxa_names(phylo) %in% input$contamCandidatesSelectRemove]
    filterMessage <- paste0("<br>",Sys.time()," - Removed contaminants: ", length(input$contamCandidatesSelectRemove))

    vals$datasets[[currentSet()]] <- filter_taxa_custom(otus_to_keep = otus_to_keep,
                                                        message = filterMessage,
                                                        current_dataset = vals$datasets[[currentSet()]])
  }
})
