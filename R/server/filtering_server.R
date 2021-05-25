
observeEvent(input$filterApplySamples, {
  if(!is.null(currentSet())){
    message("Filtering Samples ...")
    #save "old" dataset to reset filters later; only if there are no sample-filters applied to the current set
    if(!vals$datasets[[currentSet()]]$filtered){
      vals$datasets[[currentSet()]]$old.dataset <- vals$datasets[[currentSet()]]
    }
    
    #convert to datatable for filtering to work
    meta <- data.table(vals$datasets[[currentSet()]]$metaData, keep.rownames = F)
    meta_changed = F 
    
    #filter by specific sample names
    if(!is.null(input$filterSample)){
      meta <- meta[meta[[sample_column]] %in% input$filterSample,]
      meta_changed = T
      message(paste0(Sys.time()," - filtered samples: ", input$filterSample))
    }
    
    
    #filter by samples groups
    if(input$filterColumns != "NONE" ){
      #subset metatable by input 
      meta <- meta[get(input$filterColumns) == input$filterColumnValues,]
      meta_changed = T
      message(paste0(Sys.time()," - filtered sample-groups: ",input$filterColumns, "==",input$filterColumnValues))
    }
    
    if(meta_changed){
      #replace metaData
      vals$datasets[[currentSet()]]$metaData <- data.frame(meta, row.names = meta[[sample_column]])
      
      #remember that this dataset is filtered now
      vals$datasets[[currentSet()]]$filtered = T
      
      #build new dataset with filtered meta 
      filtered_samples <- as.vector(meta[[sample_column]])
      #adapt otu-tables to only have samples, which were not removed by filter
      vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[,filtered_samples,drop=F]
      vals$datasets[[currentSet()]]$normalizedData <- vals$datasets[[currentSet()]]$normalizedData[,filtered_samples,drop=F]
      vals$datasets[[currentSet()]]$relativeData <- vals$datasets[[currentSet()]]$relativeData[,filtered_samples,drop=F]
      
      #build new phyloseq-object
      py.otu <- otu_table(vals$datasets[[currentSet()]]$normalizedData,T)
      py.tax <- tax_table(as.matrix(vals$datasets[[currentSet()]]$taxonomy))
      py.meta <- sample_data(data.frame(meta, row.names = meta[[sample_column]])) # with new meta df
      sample_names(py.meta) <- filtered_samples
      tree <- vals$datasets[[currentSet()]]$tree
      
      #cannot build phyloseq object with NULL as tree input; have to check both cases:
      if (!is.null(tree)) phylo <- merge_phyloseq(py.otu,py.tax,py.meta, tree) else phylo <- merge_phyloseq(py.otu,py.tax,py.meta)
      vals$datasets[[currentSet()]]$phylo <- phylo
      message(paste0(Sys.time()," - filtered dataset: "))
      message(nsamples(phylo))
            
      #re-calculate unifrac distance
      #pick correct subset of unifrac distance matrix, containing only the new filtered samples
      if(!is.null(tree)) unifrac_dist <- as.dist(as.matrix(vals$datasets[[currentSet()]]$unifrac_dist)[filtered_samples,filtered_samples]) else unifrac_dist <- NULL
      vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist 
    }
  }
})

observeEvent(input$filterApplyTaxa,{
  if(!is.null(currentSet())){
    message("Filtering taxa ...")
    #save "old" dataset to reset filters later; only if there are no taxa-filters applied to the current set
    if(!vals$datasets[[currentSet()]]$filtered){
      vals$datasets[[currentSet()]]$old.dataset <- vals$datasets[[currentSet()]]
    }
    
    taxonomy <- data.table(vals$datasets[[currentSet()]]$taxonomy,keep.rownames = T)
    
    #subset taxonomy by input
    taxonomy <- taxonomy[get(input$filterTaxa) == input$filterTaxaValues,]
    message(paste0(Sys.time()," - keeping taxa: ", input$filterTaxaValues))
    #fix rownames and replace taxonomy
    taxonomy <- data.frame(taxonomy)
    rownames(taxonomy)<-taxonomy$rn
    remainingOTUs <- taxonomy$rn
    taxonomy$rn <- NULL
    vals$datasets[[currentSet()]]$taxonomy <- taxonomy
    vals$datasets[[currentSet()]]$filtered = T
    
    #adapt otu-tables to only have OTUs, which were not removed by filter
    vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[remainingOTUs,]
    
    #recalculate the relative abundances and normalize again 
    normalizedData <- normalizeOTUTable(vals$datasets[[currentSet()]]$rawData, vals$datasets[[currentSet()]]$normMethod)
    vals$datasets[[currentSet()]]$normalizedData <- normalizedData$norm_tab
    vals$datasets[[currentSet()]]$relativeData <- normalizedData$rel_tab
    
    #build new phyloseq-object
    py.otu <- otu_table(vals$datasets[[currentSet()]]$normalizedData,T)
    py.tax <- tax_table(as.matrix(taxonomy)) # with new taxonomy df
    py.meta <- sample_data(data.frame(vals$datasets[[currentSet()]]$metaData))
    tree <- vals$datasets[[currentSet()]]$tree
    
    #cannot build phyloseq object with NULL as tree input; have to check both cases:
    if (!is.null(tree)) phylo <- merge_phyloseq(py.otu,py.tax,py.meta, tree) else phylo <- merge_phyloseq(py.otu,py.tax,py.meta)
    vals$datasets[[currentSet()]]$phylo <- phylo
    message(paste0(Sys.time()," - filtered dataset: "))
    message(ntaxa(phylo))
    
    #recalculate unifrac distance in this case
    if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalizedData$norm_tab, tree) else unifrac_dist <- NULL
    vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist 
    
  }
})

observeEvent(input$filterResetA, {
  if(!is.null(currentSet())){
    #check if there filters applied to dataset
    if(vals$datasets[[currentSet()]]$filtered){
      message("reseting dataset ...")
      restored_dataset <- vals$datasets[[currentSet()]]$old.dataset
      vals$datasets[[currentSet()]] <- restored_dataset
      #remove old dataset from restored dataset
      vals$datasets[[currentSet()]]$old.dataset <- NULL
      #dataset is not filtered anymore
      vals$datasets[[currentSet()]]$filtered <- F
    }
  }
})

observeEvent(input$filterResetB, {
  if(!is.null(currentSet())){
    #check if there filters applied to dataset
    if(vals$datasets[[currentSet()]]$filtered){
      message("reseting dataset ...")
      restored_dataset <- vals$datasets[[currentSet()]]$old.dataset
      vals$datasets[[currentSet()]] <- restored_dataset
      #remove old dataset from restored dataset
      vals$datasets[[currentSet()]]$old.dataset <- NULL
      #dataset is not filtered anymore
      vals$datasets[[currentSet()]]$filtered <- F
    }
  }
})

observe({
  if(input$advFilterMinAbundance){enable("advFilterMinAbundanceValue")}else{disable("advFilterMinAbundanceValue")}
  if(input$advFilterMaxAbundance){enable("advFilterMaxAbundanceValue")}else{disable("advFilterMaxAbundanceValue")}
  if(input$advFilterRelAbundance){enable("advFilterRelAbundanceValue")}else{disable("advFilterRelAbundanceValue")}
  if(input$advFilterNumSamples){enable("advFilterNumSamplesValue")}else{disable("advFilterNumSamplesValue")}
  if(input$advFilterMaxVariance){enable("advFilterMaxVarianceValue")}else{disable("advFilterMaxVarianceValue")}
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
                   otu = as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo)))
    keep_taxa = NULL
    # apply different filtering functions
    if(input$advFilterMinAbundance){
      keep_taxa = names(which(f_list$x > input$advFilterMinAbundanceValue))
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa)
    }
    if(input$advFilterMaxAbundance){
      keep_taxa = names(which(f_list$x < input$advFilterMaxAbundanceValue))
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa)
    }
    if(input$advFilterRelAbundance){
      keep_taxa = names(which((f_list$x / sum(f_list$x)) > input$advFilterRelAbundanceValue))
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa)
    }
    if(input$advFilterNumSamples){
      keep_taxa = rownames(f_list$otu[rowSums(f_list$otu == 0) < input$advFilterNumSamplesValue, ])
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa)
    }
    if(input$advFilterMaxVariance){
      keep_taxa = names(sort(genefilter::rowVars(f_list$otu), decreasing = T)[1:input$advFilterMaxVarianceValue])
      f_list <- applyFilterFunc(f_list$phylo, keep_taxa)
    }
    if(is.null(keep_taxa)){
      tmp <- applyFilterFunc(f_list$phylo, keep_taxa)
      return()
    }
    #adapt otu-tables to only have OTUs, which were not removed by filter
    vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[keep_taxa,]
    
    #recalculate the relative abundances and normalize again 
    normalizedData <- normalizeOTUTable(vals$datasets[[currentSet()]]$rawData, vals$datasets[[currentSet()]]$normMethod)
    vals$datasets[[currentSet()]]$normalizedData <- normalizedData$norm_tab
    vals$datasets[[currentSet()]]$relativeData <- normalizedData$rel_tab
    
    vals$datasets[[currentSet()]]$phylo <- f_list$phylo
    tree <- vals$datasets[[currentSet()]]$tree 
    if(!is.null(tree)){vals$datasets[[currentSet()]]$tree <- phy_tree(f_list$phylo)}
  
    #recalculate unifrac distance in this case
    if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalizedData$norm_tab, vals$datasets[[currentSet()]]$tree) else unifrac_dist <- NULL
    vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist 
    
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
      restored_dataset <- vals$datasets[[currentSet()]]$old.dataset
      vals$datasets[[currentSet()]] <- restored_dataset
      #remove old dataset from restored dataset
      vals$datasets[[currentSet()]]$old.dataset <- NULL
      #dataset is not filtered anymore
      vals$datasets[[currentSet()]]$filtered <- F
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

output$advFilterMaxAbundancePlot <- renderPlotly({
  if(!is.null(currentSet())){
    abundances <- data.frame(abundance=unlist(taxa_sums(vals$datasets[[currentSet()]]$phylo)))
    p <- plot_ly(x=abundances$abundance, type="histogram", nbinsx=ntaxa(vals$datasets[[currentSet()]]$phylo)/10)
    p %>% layout(shapes=list(vline(input$advFilterMaxAbundanceValue)), xaxis=list(title="abundance (summed up over all samples)"),
                 yaxis=list(type="log"))
  }
})

output$advFilterRelAbundancePlot <- renderPlotly({
  if(!is.null(currentSet())){
    x <- taxa_sums(vals$datasets[[currentSet()]]$phylo) 
    abundances <- data.frame(relative_abundance=unlist(x/sum(x))) # get summed up rel abundance value of taxa over all samples
    p <- plot_ly(x=abundances$relative_abundance, type="histogram", nbinsx=ntaxa(vals$datasets[[currentSet()]]$phylo)/10)
    p %>% layout(shapes=list(vline(input$advFilterRelAbundanceValue)), xaxis=list(title="relative abundance (summed up over all samples)"),
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
