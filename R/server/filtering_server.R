
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
    if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalizedData$norm_tab, data.frame(vals$datasets[[currentSet()]]$metaData), tree) else unifrac_dist <- NULL
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

