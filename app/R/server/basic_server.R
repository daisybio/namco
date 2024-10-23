# update targets table of the currently loaded dataset
output$metaTable <- renderDataTable({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_meta){
      meta_dt<-datatable(sample_data(vals$datasets[[currentSet()]]$phylo),filter='top',selection=list(mode='multiple'),options=list(pageLength=20,scrollX=T))
      meta_dt
    }
  } 
  else datatable(data.frame(),options=list(dom="t"))
},server=T)

# download button for combined download
output$downloadMetaOTU <- downloadHandler(
  filename = function(){"combined_data.tsv"},
  content = function(file){
    if(!is.null(currentSet())){
      if(input$downloadMetaOTUrelativeAbundance == 'relative'){
        phylo <- phyloseq::transform_sample_counts(vals$datasets[[currentSet()]]$phylo, function(x) x/sum(x))
      }else{
        phylo <- vals$datasets[[currentSet()]]$phylo
      }
      if(input$downloadMetaOTUTaxLevel != "OTU/ASV"){
        phylo <- glom_taxa_custom(phylo, input$downloadMetaOTUTaxLevel)$phylo_rank
      }

      d <- data.frame(cbind(phylo@sam_data, t(phylo@otu_table)))
      write.table(d, file, quote = F, sep="\t", row.names = F)
    }
  }
)

####rarefaction curves####
rarefactionReactive <- reactive({
  if(!is.null(currentSet())){
    waiter_show(html = tagList(spin_rotating_plane(),"Preparing plot ... "),color=overlay_color)
    #needs integer values to work
    tab = as.matrix(vals$datasets[[currentSet()]]$rawData)
    class(tab)<-"integer"
    
    rarefactionCurve <- lapply(1:ncol(tab),function(i){
      if(colSums(tab)[i]==0){
        c(0)
      }else{
        n = seq(1,colSums(tab)[i],by=30)
        if(n[length(n)]!=colSums(tab)[i]){
          n=c(n,colSums(tab)[i])
        }
        drop(rarefy(t(tab[,i]),n))  
      }
    })
    return(list(tab=tab, rarefactionCurve=rarefactionCurve))
  }
})

output$rarefacCurve <- renderPlotly({
  if(!is.null(rarefactionReactive())){
    #slope = apply(tab,2,function(x) rareslope(x,sum(x)-100))
    rarefactionCurve <- rarefactionReactive()$rarefactionCurve
    tab <- rarefactionReactive()$tab
    slopesDF = calcSlopes(rarefactionCurve,tab)
    slopes <- as.numeric(slopesDF[,2])
    #samples_order <- order(as.numeric(slopesDF[,2]), decreasing = TRUE)
    ordered_samples <- slopesDF[order(as.numeric(slopesDF[,2]),decreasing = F)] #order samples by slope; low slope == potentially undersampled
    vals$undersampled = tail(ordered_samples,input$rareToHighlight)
    
    first = order(slopes,decreasing=F)[1]
    p <- plot_ly(x=attr(rarefactionCurve[[first]],"Subsample"),y=rarefactionCurve[[first]],text=paste0(colnames(tab)[first],"; slope: ",round(1e5*slopes[first],3),"e-5"),hoverinfo="text",color="high",type="scatter",mode="lines",colors=c("black","red"))
    
    for(i in order(slopes,decreasing=F)[2:(ncol(tab)-input$rareToHighlight)]){
      highslope = ifelse(i == input$rareToHighlight,2,1) #highlight last x slopes, since they have lowest value
      p <- p %>% add_trace(x=attr(rarefactionCurve[[i]],"Subsample"),y=rarefactionCurve[[i]],text=paste0(colnames(tab)[i],"; slope: ",round(1e5*slopes[i],3),"e-5"),hoverinfo="text",color=c("high"),showlegend=F)
    }
    for(i in order(slopes,decreasing = F)[ncol(tab)-input$rareToHighlight+1:ncol(tab)]){
      p <- p %>% add_trace(x=attr(rarefactionCurve[[i]],"Subsample"),y=rarefactionCurve[[i]],text=paste0(colnames(tab)[i],"; slope: ",round(1e5*slopes[i],3),"e-5"),hoverinfo="text",color=c("low"),showlegend=F)
    }
    p %>% layout(title="Rarefaction Curves",xaxis=list(title="Number of Reads"),yaxis=list(title="Number of Species"))
    waiter_hide()
    p
    
  }
})

# show undersampled samples
output$undersampled <- renderText({
  paste0("The following samples might be undersampled:\n",paste0(vals$undersampled,collapse=", "))
})

####taxonomic binning####

taxBinningReact <- reactive({
  if(!is.null(currentSet())){
    phylo <- vals$datasets[[currentSet()]]$phylo
    rel_dat <- vals$datasets[[currentSet()]]$relativeData
    
    if(input$taxBinningGroup == input$taxBinningYLabel && vals$datasets[[currentSet()]]$has_meta){
      #stop("Cannot select same group variable for to split taxonomic binning and for the y-label.")
      return(NULL)
    }
    
    waiter_show(html = tagList(spin_rotating_plane(),"Generating data ... "),color=overlay_color)
    
    #create phyloseq-object with relative abundance data
    #otu_table(phylo) <- otu_table(rel_dat,T)
    rel_phylo <- merge_phyloseq(otu_table(rel_dat,T),tax_table(phylo))
    #TODO: move tax binning calculation during data upload, then no re-calculation has to be done each time..
    tax_binning <- taxBinningNew(if(input$taxaAbundanceType)rel_phylo else phylo, vals$datasets[[currentSet()]]$is_fastq, input$taxBinningLevel)
    
    if(vals$datasets[[currentSet()]]$is_fastq){
      binning <- tax_binning[[which(c("Kingdom","Phylum","Class","Order","Family","Genus")==input$taxBinningLevel)]]
    }else{
      binning <- tax_binning[[which(c("Kingdom","Phylum","Class","Order","Family","Genus","Species")==input$taxBinningLevel)]] 
    }
    
    if(input$taxBinningTop < nrow(binning) && !is.na(input$taxBinningTop)){
      top_taxa <- names(sort(rowSums(binning[-1]), decreasing = T)[1:input$taxBinningTop])
      taxa <- ifelse(rownames(binning) %in% top_taxa,rownames(binning),"Other")
      other <- data.frame(binning[which(taxa=="Other"),])
      other <- colSums(other[-1])
      other <- c('Other',other)
      
      binning <- binning[-which(taxa=="Other"),]
      binning <- rbind(binning, Other=other)
    }else{
      top_taxa <- rownames(binning)
    }
    
    if(vals$datasets[[currentSet()]]$has_meta){
      meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo), check.names = F)
      tab <- merge(melt(binning, id.vars = input$taxBinningLevel), meta, by.x = "variable", by.y=sample_column, all.x=T)
    }else{
      tab <- melt(binning, id.vars = input$taxBinningLevel)
    }
    colnames(tab)[which(colnames(tab)=="variable")] <- sample_column
    colnames(tab)[which(colnames(tab)==input$taxBinningLevel)] <- "custom_taxonomy_column"
    colnames(tab)[which(colnames(tab)==input$taxBinningYLabel)] <- "y_split"
    colnames(tab)[which(colnames(tab)==input$taxBinningGroup)] <- "facet_split"
    tab$value <- as.numeric(tab$value)
    
    if('y_split' %in% colnames(tab)) tab$y_split <- as.character(tab$y_split)
    if('facet_split' %in% colnames(tab)) tab$facet_split <- as.character(tab$facet_split)
    
    waiter_hide()
    return(tab)
  }
})

taxBinningPlotReact <- reactive({
  if(!is.null(taxBinningReact())){
    
    waiter_show(html = tagList(spin_rotating_plane(),"Generating plot ... "),color=overlay_color)
    
    tab <- taxBinningReact()
    
    if(!vals$datasets[[currentSet()]]$has_meta){
      #colnames(tab)[which(colnames(tab)==sample_column)] <- "SampleID"
      #case0: no meta file
      p <- ggplot(tab, aes(x=value, y=y_split, fill=custom_taxonomy_column))+
        geom_col()+
        xlab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance"))+
        ylab(sample_column)+
        scale_fill_manual(name=input$taxBinningLevel,values = colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(tab$custom_taxonomy_column))))+
        ggtitle(paste0("Taxonomic Binning of samples by ",sample_column))
    }
    
    else if(input$taxBinningYLabel != "--Combined--" && input$taxBinningGroup != "None"){
      
      # case1: split by variable and group y axis by other variable
      # subcase: use relative abundance -> calculate mean
      if(input$taxaAbundanceType){
        tab <- as.data.table(tab)
        tab <- tab[, mean(value), by=c("y_split","facet_split","custom_taxonomy_column")]
        colnames(tab)[which(colnames(tab)=="V1")] <- "value"
      }
      # order bars by selected reference and place reference at first position of stacked bar
      if(input$taxBinningOrderReference != "None"){
        tmp <- as.data.table(tab[which(tab$custom_taxonomy_column %in% c(input$taxBinningOrderReference)),])
        tmp[,value_group := sum(value), by=y_split] # sum of abundance by group
        tmp <- tmp[with(tmp, order(value_group)),]
        ordered_y <- unique(tmp$y_split)
        tab$y_split <- factor(tab$y_split, levels=ordered_y)
        reordered_taxa <- c(setdiff(unique(tab$custom_taxonomy_column), input$taxBinningOrderReference), input$taxBinningOrderReference)
        tab$custom_taxonomy_column <- factor(tab$custom_taxonomy_column, levels=reordered_taxa)
      }
      p <- ggplot(tab, aes(x=value, y=y_split, fill=custom_taxonomy_column))+
        geom_col()+
        facet_grid(~facet_split, scales='free_y')+
        xlab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance"))+
        ylab(input$taxBinningYLabel)+
        scale_fill_manual(name=input$taxBinningLevel,values = colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(tab$custom_taxonomy_column))))+
        ggtitle(paste0("Taxonomic Binning of samples by ",input$taxBinningGroup," and ", input$taxBinningYLabel))
      
    }else if(input$taxBinningYLabel == "--Combined--" && input$taxBinningGroup != "None"){
      
      # case2: combined y axis, split only by facet (one bar in each facet)
      # subcase: use relative abundance -> divide by number of groups in y axis to get to a total of 100
      if(input$taxaAbundanceType){
        tab <- as.data.table(tab)
        tab <- tab[,value/.N, by=c("facet_split","custom_taxonomy_column")]
        colnames(tab)[which(colnames(tab)=="V1")] <- "value"
      }
      
      p <- ggplot(tab, aes(x=facet_split, y=value, fill=custom_taxonomy_column))+
        geom_col()+
        facet_grid(~facet_split, scales="free")+
        ylab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance"))+
        xlab(input$taxBinningYLabel)+
        scale_fill_manual(name=input$taxBinningLevel,values = colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(tab$custom_taxonomy_column))))+
        ggtitle(paste0("Taxonomic Binning of samples by ",input$taxBinningGroup))
      
    }else if(input$taxBinningYLabel != "--Combined--" && input$taxBinningGroup == "None"){
      
      # case3: group y axis by variable, no facets
      # subcase: use relative abundance -> calculate mean
      if(input$taxaAbundanceType){
        tab <- as.data.table(tab)
        tab <- tab[,value/.N, by=c("y_split","custom_taxonomy_column")]
        colnames(tab)[which(colnames(tab)=="V1")] <- "value"
      }
      # order bars by selected reference and place reference at first position of stacked bar
      if(input$taxBinningOrderReference != "None"){
        tmp <- as.data.table(tab[which(tab$custom_taxonomy_column %in% c(input$taxBinningOrderReference)),])
        tmp[,value_group := sum(value), by=y_split] # sum of abundance by group
        tmp <- tmp[with(tmp, order(value_group)),]
        ordered_y <- unique(tmp$y_split)
        tab$y_split <- factor(tab$y_split, levels=ordered_y)
        reordered_taxa <- c(setdiff(unique(tab$custom_taxonomy_column), input$taxBinningOrderReference), input$taxBinningOrderReference)
        tab$custom_taxonomy_column <- factor(tab$custom_taxonomy_column, levels=reordered_taxa)
      }
      p <- ggplot(tab, aes(x=value, y=y_split, fill=custom_taxonomy_column))+
        geom_col()+
        xlab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance"))+
        ylab(input$taxBinningYLabel)+
        scale_fill_manual(name=input$taxBinningLevel,values = colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(tab$custom_taxonomy_column))))+
        ggtitle(paste0("Taxonomic Binning of samples by ",input$taxBinningYLabel))
      
    }else{
      #case4: no plot
      waiter_hide()
      return(NULL)
    }
    
    # manually change order of y axis labels
    if(input$taxBinningYLabel != '--Combined--' && input$taxBinningOrderManually){
      p <- p + scale_y_discrete(limits=input$taxBinningYOrder)
    }
    if(input$taxBinningRotate){
      p <- p + coord_flip() + theme(axis.text.x = element_text(angle = 90))
    }
    
    p <- p + theme_bw()
    
    # hide y labels
    if(!input$taxBinningShowY){
      if(input$taxBinningRotate){
        p <- p + theme(axis.text.x = element_blank())
      }else{
        p <- p + theme(axis.text.y = element_blank())
      }
    }
    
    waiter_hide()
    list(py=ggplotly(p, height = 800),gg=p)
  }
})

observe({
  if(!is.null(taxBinningReact())){
    tab <- taxBinningReact()
    taxa <- as.character(unique(tab$custom_taxonomy_column))
    updateSelectInput(session, "taxBinningOrderReference", choices = c("None", taxa))
  }
})

# plot distribution of taxa
output$taxaDistribution <- renderPlotly({
  if(!is.null(taxBinningPlotReact())){
    taxBinningPlotReact()$py
  } else plotly_empty()
  
})

#download as pdf
output$taxaPDF <- downloadHandler(
  filename = function(){"taxonomic_binning.pdf"},
  content = function(file){
    if(!is.null(taxBinningPlotReact())){
      ggsave(file, taxBinningPlotReact()$gg, device="pdf", width = 16, height = 10)
    }
  }
)


####alpha diversity####
#reactive alpha-diversity table; stores all measures for alpha-div of current set
alphaReact <- reactive({
  if(!is.null(currentSet())){
    otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
    
    alphaTabFull <- vals$datasets[[currentSet()]]$alpha_diversity

    alphaTab <- gather(alphaTabFull, measure, value, Shannon_Entropy:Richness)
    alphaTab <- alphaTab[alphaTab$measure %in% c(input$alphaMethod),]
    
    if(input$alphaGroup=="-") {
      p <- ggboxplot(alphaTab, x="measure", y="value", fill="#0072B2", facet.by = "measure",
                     scales=as.character(input$alphaScalesFree))+
        rremove("x.text")+
        ggtitle(paste("Alpha Diversity of all samples"))
    }else if(input$alphaGroup==sample_column){
      p <- ggdotplot(alphaTab, x=sample_column, y='value', size=.75, fill='black', ggtheme=theme_bw(), facet.by = 'measure', binwidth = 1)
      p <- p + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    }else{
      pairs <- sapply(input$alphaPairs, strsplit, split=" vs. ")
      p <- suppressWarnings(ggboxplot(alphaTab, x=input$alphaGroup, y="value", fill=input$alphaGroup, facet.by = "measure",
                                     palette = colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(alphaTab[[input$alphaGroup]]))), 
                                     scales = input$alphaScalesFree)+
                            rremove("x.text")+
                            ggtitle(paste("Alpha Diversity colored by",input$alphaGroup)))
      if(!is.null(pairs)){
        p <- p + stat_compare_means(comparisons = pairs,
                                    label=ifelse(input$alphaSignifView=="p-value","p.format","p.signif"))
      }
    } 
    
    list(gg=p, alphaTabFull=alphaTabFull)
  }else{
    NULL
  }
})

observe({
  if(input$alphaGroup == "-"){
    shinyjs::hide("alphaShowSamples", anim=T)
  }else{
    shinyjs::show("alphaShowSamples", anim=T)
  }
})

# plot alpha diversity
output$alphaPlot <- renderPlot({
  if(!is.null(alphaReact())){
    alphaReact()$gg
  }
})

# show alpha diversity table
output$alphaTable <- renderTable({
  if(!is.null(alphaReact())){
    alphaReact()$alphaTabFull[1:6]
  }
})

output$alphaTableDownload <- downloadHandler(
  filename=function(){paste("alpha_diversity.csv")},
  content = function(file){
    if(!is.null(alphaReact())){
      write.csv(alphaReact()$alphaTabFull,file,row.names = F)
    }
  }
)

#download as pdf
output$alphaPDF <- downloadHandler(
  filename = function(){"alpha_diversity.pdf"},
  content = function(file){
    if(!is.null(alphaReact())){
      ggsave(file, alphaReact()$gg, device="pdf", width = 15, height = 12)
    }
  }
)

####beta diversity####
#do calculation for beta diversity plots
betaReactive <- reactive({
  if(!is.null(currentSet())){
    waiter_show(html = tagList(spin_rotating_plane(),"Calculating distance ... "),color=overlay_color)
    
    group <- input$betaGroup
    if(input$betaGroup2 != "None") group2 <- input$betaGroup2 else group2 <- NULL
    phylo <- vals$datasets[[currentSet()]]$phylo
    
    meta <- data.frame(sample_data(phylo), check.names = F)
    meta <- data.frame(meta[order(rownames(meta)),], check.names = F)
    # remove all samples with NA in column of chosen group
    meta <- meta[!is.na(meta[[group]]),]
    if(nrow(meta) == 0){
      waiter_hide()
      return(NULL)
    }
    samples <- meta[[sample_column]]
    phylo <- prune_samples(samples, phylo)
    
    if(input$betaLevel != "All"){
      keep_samples <- meta[meta[[group]]==input$betaLevel,][["SampleID"]]
      phylo <- prune_samples(keep_samples, phylo)
      meta <- as.data.frame(sample_data(phylo))
    }
    
    group_vector <- as.factor(meta[[group]])
    group_vector2 <- if(!is.null(group2)) as.factor(meta[[group2]]) else NA
    
    if(!is.null(access(phylo,"phy_tree"))) tree <- phy_tree(phylo) else tree <- NULL
    
    method <- match(input$betaMethod,c("Bray-Curtis Dissimilarity","Generalized UniFrac Distance", "Unweighted UniFrac Distance", "Weighted UniFrac Distance", "Variance adjusted weighted UniFrac Distance"))
    my_dist <- betaDiversity(phylo, method=method)
    
    all_fit <- hclust(my_dist, method = 'ward.D2')
    tree <- as.phylo(all_fit)

    pval <- 0
    if(input$betaLevel == "All"){
      colnames(meta)[which(colnames(meta) == group)] <- "condition"
      tryCatch({
        adonis <- adonis2(my_dist ~ condition, data=meta, parallel = ncores)
        pval <- adonis[["Pr(>F)"]][1]
      }, error=function(e){
        print(e$message)
        showModal(errorModal(paste0("Error with adonis2 at beta-diversity: ", group)))
        pval <- 0
      })
    }

    pcoa <- data.frame(cmdscale(my_dist,k=2))
    pcoa$group <- group_vector
    pcoa$group2 <- group_vector2
    pcoa$method <- "Principal coordinates analysis (PCoA)"
    pcoa$sample <- meta[[sample_column]]
    
    
    nmds <- metaMDS(my_dist,k=2)
    nmds_stress <- round(nmds$stress, digits = 4)
    nmds_obj <- nmds
    nmds <- data.frame(nmds_obj$points)
    colnames(nmds) <- c("X1","X2")
    nmds$group <- group_vector
    nmds$group2 <- group_vector2
    nmds$method <- "non-metric multidimensional scaling (NMDS)"
    nmds$sample <- meta[[sample_column]]
    
    plot_df <- rbind(pcoa, nmds)
    show_ellipse <- (length(unique(plot_df[["group"]]))>1 & input$betaShowEllipse)
    centroids <- aggregate(cbind(X1,X2)~group+method, data=plot_df, mean)
    plot_df <- merge(plot_df, centroids, by=c("group", "method"), suffixes=c("",".centroid"))
    colors <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(plot_df$group)))
    
    beta_pcoa <- ggplot(plot_df[plot_df$method == 'Principal coordinates analysis (PCoA)',])+
      geom_point(data=centroids[centroids$method == 'Principal coordinates analysis (PCoA)',], aes(x=X1, y=X2, color=group), size=5)+
      geom_segment(aes(x=X1.centroid, y=X2.centroid, xend=X1, yend=X2, color=group))+
      theme_bw()+
      xlab("")+ylab("")+
      ggtitle(paste0("Principal Coordinates Analysis using ", input$betaMethod))+
      annotate("text",
               x=min(plot_df[plot_df$method == 'Principal coordinates analysis (PCoA)',]$X1)+input$betaDivTextSliderX, 
               y=min(plot_df[plot_df$method == 'Principal coordinates analysis (PCoA)',]$X2)+input$betaDivTextSliderY,
               label="")+
      labs(caption = paste0("p-value (PERMANOVA): ", pval))+
      theme(
        plot.caption = element_text(size = 12, hjust = 0)
      )

    beta_nmds <- ggplot(plot_df[plot_df$method == 'non-metric multidimensional scaling (NMDS)',])+
      geom_point(data=centroids[centroids$method == 'non-metric multidimensional scaling (NMDS)',], aes(x=X1, y=X2, color=group), size=5)+
      geom_segment(aes(x=X1.centroid, y=X2.centroid, xend=X1, yend=X2, color=group))+
      theme_bw()+
      xlab("")+ylab("")+
      ggtitle(paste0("Non-metric Multidimensional Scaling using ", input$betaMethod))+
      annotate("text",
               x=min(plot_df[plot_df$method == 'non-metric multidimensional scaling (NMDS)',]$X1)+input$betaDivTextSliderX, 
               y=min(plot_df[plot_df$method == 'non-metric multidimensional scaling (NMDS)',]$X2)+input$betaDivTextSliderY,
               label="") +
      labs(caption = paste0("p-value (PERMANOVA): ", pval, "; stress (NMDS): ", nmds_stress)) +
      theme(
        plot.caption = element_text(size = 12, hjust = 0)
      )
      
    if(input$betaShowLabels){
      beta_pcoa <- beta_pcoa + geom_label_repel(aes(x=X1, y=X2, color=group, label=sample))
      beta_nmds <- beta_nmds + geom_label_repel(aes(x=X1, y=X2, color=group, label=sample))
    }
    
    if(!is.null(group2)){
      beta_pcoa <- beta_pcoa + 
        geom_point(aes(x=X1, y=X2, color=group, shape=group2, group=sample), size=3)+
        labs(fill=input$betaGroup, color=input$betaGroup, shape=input$betaGroup2)
      
      beta_nmds <- beta_nmds + 
        geom_point(aes(x=X1, y=X2, color=group, shape=group2, group=sample), size=3)+
        labs(fill=input$betaGroup, color=input$betaGroup, shape=input$betaGroup2)
    }else{
      beta_pcoa <- beta_pcoa + 
        geom_point(aes(x=X1, y=X2, color=group, group=sample), size=3)+
        labs(fill=input$betaGroup, color=input$betaGroup)
      
      beta_nmds <- beta_nmds + 
        geom_point(aes(x=X1, y=X2, color=group, group=sample), size=3)+
        labs(fill=input$betaGroup, color=input$betaGroup)
    }
    if(show_ellipse){
      beta_pcoa <- beta_pcoa + stat_ellipse(aes(x=X1, y=X2, fill=group), geom = "polygon", alpha=.2)
      beta_nmds <- beta_nmds + stat_ellipse(aes(x=X1, y=X2, fill=group), geom = "polygon", alpha=.2)
    }
    
    beta_pcoa <- beta_pcoa +      
      scale_fill_manual(values = colors)+
      scale_color_manual(values = colors)
    
    beta_nmds <- beta_nmds +      
      scale_fill_manual(values = colors)+
      scale_color_manual(values = colors)
    
    
    out <- list(dist=my_dist, all_groups=group_vector, tree=tree, pval=pval, beta_pcoa=beta_pcoa, beta_nmds=beta_nmds, colors=colors, nmds_obj=nmds_obj)
    waiter_hide()
    return(out)
  }
})

# clustering tree of samples based on beta-diversity
output$betaTree <- renderPlot({
  if(!is.null(betaReactive())){
    beta <- betaReactive()
    plot(beta$tree,type="phylogram",use.edge.length=T,tip.color=betaReactive()$colors[beta$all_groups],label.offset=0.01)
    axisPhylo()
    tiplabels(pch=16,col=betaReactive()$colors[beta$all_groups])
    
    title(main = paste("Hierarchical Clustering using", input$betaMethod, "Distance"), font.main = 1)
  }
})

#download as pdf
output$betaTreePDF <- downloadHandler(
  filename = function(){"beta_diversity_clustering.pdf"},
  content = function(file){
    if(!is.null(betaReactive())){
      beta <- betaReactive()
      pdf(file, width=14, height=10)
      plot(beta$tree,type="phylogram",use.edge.length=T,tip.color=betaReactive()$colors[beta$all_groups],label.offset=0.01)
      axisPhylo()
      tiplabels(pch=16,col=betaReactive()$colors[beta$all_groups])
      dev.off()
    }
  }
)

output$betaDivStress <- renderPlot({
  if(!is.null(betaReactive())){
    nmds <- betaReactive()$nmds_obj
    vegan::stressplot(nmds, main='Shepard diagram (stressplot) of NMDS calculation')
  }
})

output$betaDivStressPDF <- downloadHandler(
  filename = function(){"beta_diversity_nmds_stressplot.pdf"},
  content = function(file){
    if(!is.null(betaReactive())){
      nmds <- betaReactive()$nmds_obj
      pdf(file, width=14, height=10)
      vegan::stressplot(nmds, main='Shepard diagram (stressplot) of NMDS calculation')
      dev.off()
    }
  }
)

output$betaDivPcoa <- renderPlot({
  if(!is.null(betaReactive())){
    betaReactive()$beta_pcoa
  }
})

output$betaDivNMDS <- renderPlot({
  if(!is.null(betaReactive())){
    betaReactive()$beta_nmds
  }
})

output$betaDivPocaPDF <- downloadHandler(
  filename = function(){"beta_diversity_pcoa.pdf"},
  content = function(file){
    if(!is.null(betaReactive())){
      p <- betaReactive()$beta_pcoa
      ggsave(filename = file, plot=p, device="pdf", width = 15, height = 15)
    }
  }
)

output$betaDivNMDSPDF <- downloadHandler(
  filename = function(){"beta_diversity_nmds.pdf"},
  content = function(file){
    if(!is.null(betaReactive())){
      p <- betaReactive()$beta_nmds
      ggsave(filename = file, plot=p, device="pdf", width = 15, height = 15)
    }
  }
)

output$betaDownloadDistance <- downloadHandler(
  filename = function(){"beta_diversity_distance.tab"},
  content = function(file){
    if(!is.null(betaReactive())){
      write.table(x=as.matrix(betaReactive()$dist),
                  file = file, quote = F, sep = '\t')
    }
  }
)


####heatmap#### 

# plot heatmap of OTU abundances per sample
abundanceHeatmapReact <- reactive({
  if(!is.null(currentSet())){
    set.seed(seed)
    phylo <- vals$datasets[[currentSet()]]$phylo
    phylo <- transform_sample_counts(phylo, function(x) x+1)  # pseudocount to not get -Inf values
    
    if(nsamples(phylo) <= 2){
      showModal(errorModal('The heatmap cannot be calculated with 2 or less samples.'))
      return(NULL)
    }

    #check for unifrac distance --> (needs phylo tree file):
    l <- list()
    if(!is.null(vals$datasets[[currentSet()]]$unifrac_dist)){
      #save generalized unifrac distance as global variable to use it for heatmap
      gunifrac_heatmap <<- as.dist(vals$datasets[[currentSet()]]$unifrac_dist)
    }

    hm_distance <- if(input$heatmapDistance == "gunifrac") "gunifrac_heatmap" else input$heatmapDistance
    # plot_heatmap only accepts meta data inputs which come from check.names=T dataframes:
    meta_T <- data.frame(sample_data(phylo))
    meta_F <- data.frame(sample_data(phylo), check.names = F)
    sample.label <- colnames(meta_T)[which(colnames(meta_F) == input$heatmapSample)] 
    
    # add taxa info of specified level to taxa name -> will show up in heatmap hover info
    if(input$heatmapOverlayTaxa != 'OTU/ASV'){
      taxa_names_level <- data.frame(phylo@tax_table[, input$heatmapOverlayTaxa])
      taxa_names_combined <- paste0(rownames(taxa_names_level),' (',taxa_names_level[,1],')') 
      taxa_names(phylo) <- taxa_names_combined  
    }
    
    if(input$heatmapOrderSamples){
      sample.order <- colnames(meta_T)[which(colnames(meta_F) == input$heatmapSample)] 
      p<-plot_heatmap(phylo,method = input$heatmapOrdination, distance = hm_distance, sample.label = sample.label, sample.order = sample.order, trans = 'log10')
    }else{
      p<-plot_heatmap(phylo,method = input$heatmapOrdination, distance = hm_distance, sample.label = sample.label, trans = 'log10')
    }
    l <- list(gg=p)
  }
})

output$abundanceHeatmap <- renderPlotly({
  if(!is.null(abundanceHeatmapReact())){
    ggplotly(abundanceHeatmapReact()$gg)
  }
})

#download as pdf
output$abundanceHeatmapPDF <- downloadHandler(
  filename = function(){"abundance_heatmap.pdf"},
  content = function(file){
    if(!is.null(abundanceHeatmapReact())){
      ggsave(file, abundanceHeatmapReact()$gg, device="pdf", width = 10, height = 7)
    }
  }
)

