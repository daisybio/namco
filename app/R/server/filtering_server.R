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



observeEvent(input$filterApplySamples, {
  if(!is.null(currentSet())){
    message("Filtering Samples ...")
    #save "old" dataset to reset filters later; only if there are no sample-filters applied to the current set
    if(!vals$datasets[[currentSet()]]$filtered){
      vals$datasets[[currentSet()]]$old.dataset <- vals$datasets[[currentSet()]]
    }
    
    tryCatch({
      meta_changed = F 
      
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
    message("Filtering taxa ...")
    #save "old" dataset to reset filters later; only if there are no taxa-filters applied to the current set
    if(!vals$datasets[[currentSet()]]$filtered){
      vals$datasets[[currentSet()]]$old.dataset <- vals$datasets[[currentSet()]]
    }
    
    tryCatch({
      taxonomy <- data.frame(vals$datasets[[currentSet()]]$phylo@tax_table)
      
      #subset taxonomy by input
      taxonomy <- taxonomy[taxonomy[[input$filterTaxa]] %in% input$filterTaxaValues,]
      remainingOTUs <- rownames(taxonomy)
      filterMessage <- paste0("<br>",Sys.time()," - keeping taxa: ", paste(unlist(input$filterTaxaValues), collapse = "; "))
      message(filterMessage)
      vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,filterMessage)
      vals$datasets[[currentSet()]]$filtered = T
      
      #adapt otu-tables to only have OTUs, which were not removed by filter
      vals$datasets[[currentSet()]]$rawData <- vals$datasets[[currentSet()]]$rawData[remainingOTUs,]
      
      #recalculate the relative abundances and normalize again
      if(vals$datasets[[currentSet()]]$has_picrust){
        normMethod<-0
      }else{
        normMethod <- vals$datasets[[currentSet()]]$normMethod
      }
      normalizedData <- normalizeOTUTable(vals$datasets[[currentSet()]]$rawData, normMethod)
      vals$datasets[[currentSet()]]$normalizedData <- normalizedData$norm_tab
      vals$datasets[[currentSet()]]$relativeData <- normalizedData$rel_tab
      
      #adapt phyloseq-object
      vals$datasets[[currentSet()]]$phylo <- prune_taxa(remainingOTUs, vals$datasets[[currentSet()]]$phylo)
      vals$datasets[[currentSet()]]$phylo.raw <- prune_taxa(remainingOTUs, vals$datasets[[currentSet()]]$phylo.raw)
      message(paste0(Sys.time()," - filtered dataset: "))
      message(length(remainingOTUs))
      phylo_tree <- vals$datasets[[currentSet()]]$phylo@phy_tree
      
      #recalculate unifrac distance in this case
      if(!is.null(phylo_tree)) unifrac_dist <- buildGUniFracMatrix(normalizedData$norm_tab, phylo_tree) else unifrac_dist <- NULL
      vals$datasets[[currentSet()]]$unifrac_dist <- unifrac_dist  
      
      # re-calculate alpha-diversity
      phylo <- vals$datasets[[currentSet()]]$phylo
      if(vals$datasets[[currentSet()]]$has_meta){
        alphaTabFull <- createAlphaTab(data.frame(phylo@otu_table, check.names=F), data.frame(phylo@sam_data, check.names = F))
      }else{
        alphaTabFull <- createAlphaTab(data.frame(phylo@otu_table, check.names=F))
      }
      vals$datasets[[currentSet()]]$alpha_diversity <- alphaTabFull
      
      showModal(infoModal(paste0("Filtering successful. ", length(remainingOTUs)," OTUs are remaining.")))
    }, error=function(e){
      vals$datasets[[currentSet()]]$filterHistory <- paste(vals$datasets[[currentSet()]]$filterHistory,Sys.time()," - error during taxa filtering:", e$message)
      print(e$message)
      showModal(errorModal(e$message))
      return(NULL)
    })
    
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
    normalizedData <- normalizeOTUTable(vals$datasets[[currentSet()]]$rawData, normMethod)
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
