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
      phylo <- vals$datasets[[currentSet()]]$phylo
      d <- data.frame(cbind(phylo@sam_data, t(phylo@otu_table)))
      write.table(d, file, quote = F, sep="\t")
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
      stop("Cannot select same group variable for to split taxonomic binning and for the y-label.")
      return(NULL)
    }
    
    waiter_show(html = tagList(spin_rotating_plane(),"Preparing plot ... "),color=overlay_color)
    
    #create phyloseq-object with relative abundance data
    #otu_table(phylo) <- otu_table(rel_dat,T)
    rel_phylo <- merge_phyloseq(otu_table(rel_dat,T),tax_table(phylo))
    tax_binning <- taxBinningNew(if(input$taxaAbundanceType)rel_phylo else phylo, vals$datasets[[currentSet()]]$is_fastq)

    if(vals$datasets[[currentSet()]]$is_fastq){
      binning = tax_binning[[which(c("Kingdom","Phylum","Class","Order","Family","Genus")==input$taxBinningLevel)]]
    }else{
      binning = tax_binning[[which(c("Kingdom","Phylum","Class","Order","Family","Genus","Species")==input$taxBinningLevel)]] 
    }
    if(input$taxBinningTop < nrow(binning)){
      top_taxa <- names(sort(rowSums(binning), decreasing = T)[1:input$taxBinningTop])
      taxa <- ifelse(rownames(binning) %in% top_taxa,rownames(binning),"Other")
      other <- binning[which(taxa=="Other"),]
      other <- colSums(other)
      
      binning <- binning[-which(taxa=="Other"),]
      binning <- rbind(binning, Other=other)
    }else{
      top_taxa <- rownames(binning)
    }
    
    if(vals$datasets[[currentSet()]]$has_meta){
      meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo), check.names = F)
      tab <- merge(melt(binning), meta, by.x = "Var2", by.y=sample_column, all.x=T)
    }else{
      tab <- melt(binning)
    }
    colnames(tab)[which(colnames(tab)=="Var2")] <- sample_column

    if(input$taxBinningYLabel != "--Combined--"){
      colnames(tab)[which(colnames(tab)==input$taxBinningYLabel)] <- "reference"
      if(input$taxBinningGroup == "None"){
        p <- ggplot(tab, aes(x=value, y=as.character(reference), fill=Var1))+
          geom_bar(stat="identity")+
          xlab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance"))+
          ylab(input$taxBinningYLabel)+
          scale_fill_discrete(name = input$taxBinningLevel)+
          ggtitle(paste0("Taxonomic Binning of samples"))
      }else{
        p <- ggplot(tab, aes(x=value, y=as.character(reference), fill=Var1))+
          geom_bar(stat="identity")+
          facet_wrap(as.formula(paste0("~",input$taxBinningGroup)), scales = "free")+
          xlab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance"))+
          ylab(input$taxBinningYLabel)+
          scale_fill_discrete(name = input$taxBinningLevel)+
          ggtitle(paste0("Taxonomic Binning, grouped by ", input$taxBinningGroup))
      }
    }else if(input$taxBinningGroup != "None"){
      if(vals$datasets[[currentSet()]]$has_meta && input$taxaAbundanceType){
        # if relative abundance is used + total binning per group
        colnames(tab)[which(colnames(tab)==input$taxBinningGroup)] <- "facet_split"
        tab <- as.data.table(tab)
        tab <- tab[, value/.N, by=c("facet_split", "Var1")]
        colnames(tab) <- c(as.character(input$taxBinningGroup), "Var1", "value")
        
        #tab <- merge(melt(binning), meta, by.x = "variable", by.y=sample_column, all.x=T)
      }
      p <- ggplot(tab, aes(x=input$taxBinningGroup, y=value, fill=Var1))+
        geom_bar(stat="identity")+
        facet_wrap(as.formula(paste0("~",input$taxBinningGroup)), scales = "free")+
        ylab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance"))+
        xlab("Group")+
        scale_fill_discrete(name = input$taxBinningLevel)+
        ggtitle(paste0("Taxonomic Binning, grouped by ", input$taxBinningGroup))+
        theme(axis.text.x = element_blank())
    }else{
      waiter_hide()
      return(NULL)
    }

    waiter_hide()
    list(py=ggplotly(p, height = 800),gg=p)
  }
})

# plot distribution of taxa
output$taxaDistribution <- renderPlotly({
  if(!is.null(taxBinningReact())){
    taxBinningReact()$py
  } else plotly_empty()
  
})

#download as pdf
output$taxaPDF <- downloadHandler(
  filename = function(){"taxonomic_binning.pdf"},
  content = function(file){
    if(!is.null(taxBinningReact())){
      ggsave(file, taxBinningReact()$gg, device="pdf", width = 16, height = 7)
    }
  }
)


####dimensionality reduction (PCA, UMAP, tSNE)####
structureReact <- reactive({
  if(!is.null(currentSet())){
    waiter_show(html = tagList(spin_rotating_plane(),"Preparing plots ... "),color=overlay_color)
    
    # need to remove OTUs and samples with variance of 0 --> PCA cannot rescale them
    mat <- as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo))
    mat <- data.frame(mat[,which(apply(mat, 2, var) != 0)])
    mat_t <- t(mat)
    mat_t <- data.frame(mat_t[,which(apply(mat_t, 2, var) != 0)])
    samples = colnames(mat)
    otus = colnames(mat_t)
    
    #PCR-calculation:
    pca = prcomp(mat,center=T,scale=T)
    pca_t = prcomp(mat_t,center=T,scale=T)
    out_pca = data.frame(pca$rotation,txt=samples)
    percentage = signif(pca$sdev^2/sum(pca$sdev^2)*100,2)
    
    loadings = data.frame(Taxa=otus,loading=pca_t$rotation[,as.numeric(input$pcaLoading)])
    show_loadings_pos = ceiling(input$pcaLoadingNumber/2)
    show_loadings_neg = ceiling(input$pcaLoadingNumber/2)
    loadings = loadings[order(loadings$loading,decreasing=T),][c(1:show_loadings_pos,(nrow(loadings)-show_loadings_pos+1):nrow(loadings)),]
    loadings$Taxa = factor(loadings$Taxa,levels=loadings$Taxa)
    
    #UMAP:
    if(length(samples) < 15){
      UMAP <- umap(mat_t, n_components=3, n_neighbors=length(samples))
    }else{
      UMAP = umap(mat_t,n_components=3)    #use default with 15 nearest neighbors
    }
    out_umap = data.frame(UMAP$layout,txt=samples)
    
    #tSNE:
    tsne = Rtsne(mat_t,dim=3,perplexity=min((length(samples)-1)/3,30))
    out_tsne = data.frame(tsne$Y,txt=samples)
    
    l <- list(percentage=percentage, loadings=loadings, out_pca = out_pca, raw_pca = pca, out_umap=out_umap, out_tsne=out_tsne)
    waiter_hide()
    return(l)
  }
})

#observer if PCA is chosen:
observeEvent(input$structureMethod,{
  if(input$structureMethod=="PCA"){
    shinyjs::show("structureCompChoosing",anim = T)
  }else{
    shinyjs::hide("structureCompChoosing",anim = T)
  }
})

# visualize data structure as PCA, t-SNE or UMAP plots
output$structurePlot <- renderPlotly({
  if(!is.null(structureReact())){
    mode = input$structureDim
    if(vals$datasets[[currentSet()]]$has_meta){meta<-sample_data(vals$datasets[[currentSet()]]$phylo)}else{meta<-NULL}
    
    if(input$structureMethod=="PCA"){
      out = structureReact()$out_pca
      percentage = structureReact()$percentage
      
      xlab = paste0("PC1 (",percentage[as.numeric(input$structureCompOne)]," % variance)")
      ylab = paste0("PC2 (",percentage[as.numeric(input$structureCompTwo)]," % variance)")
      zlab = paste0("PC3 (",percentage[as.numeric(input$structureCompThree)]," % variance)")
    }
    else if(input$structureMethod=="t-SNE"){
      out = structureReact()$out_tsne
    }
    else if(input$structureMethod=="UMAP"){
      out = structureReact()$out_umap
    }
    
    xlab = input$structureCompOne
    ylab = input$structureCompTwo
    zlab = input$structureCompThree
    colnames(out) = c(paste0("Dim",1:(ncol(out)-1)), "txt")
    if (is.null(meta)){color <- "Samples"}else{color<-meta[[input$structureGroup]]}
    
    if(mode=="2D"){
      plot_ly(out,x=as.formula(paste0("~Dim",input$structureCompOne)),y=as.formula(paste0("~Dim",input$structureCompTwo)),color=color,colors="Set1",text=~txt,hoverinfo='text',type='scatter',mode="markers", size=1) %>%
        layout(xaxis=list(title=xlab),yaxis=list(title=ylab))
    } else{
      plot_ly(out,x=as.formula(paste0("~Dim",input$structureCompOne)),y=as.formula(paste0("~Dim",input$structureCompTwo)),z=as.formula(paste0("~Dim",input$structureCompThree)),color=color,colors="Set1",text=~txt,hoverinfo='text',type='scatter3d',mode="markers") %>%
        layout(scene=list(xaxis=list(title=xlab),yaxis=list(title=ylab),zaxis=list(title=zlab)))
    }
  } else{
    plotly_empty()
  }
})

# plot PCA loadings
output$loadingsPlot <- renderPlotly({
  if(!is.null(structureReact())){
    loadings = structureReact()$loadings
    
    plot_ly(loadings,x=~Taxa,y=~loading,text=~Taxa,hoverinfo='text',type='bar',color=I(ifelse(loadings$loading>0,"blue","red"))) %>%
      layout(title="Loadings",xaxis=list(title="",zeroline=F,showline=F,showticklabels=F,showgrid=F),yaxis=list(title=paste0("loadings on PC",input$pcaLoading)),showlegend=F) %>% hide_colorbar()
  } else{
    plotly_empty()
  }
})

#screeplot for PCA
output$screePlot <- renderPlot({
  if(!is.null(structureReact())){
    pca <- structureReact()$raw_pca
    var_explained <- pca$sdev**2 / sum(pca$sdev**2)
    var_explained <- var_explained[1:input$screePCshow]
    p<-qplot(c(1:length(var_explained)), var_explained)+
      geom_line()+
      xlab("Principal Component")+
      ylab("Fraction of Variance explained")+
      ggtitle("Scree-plot")+
      ylim(0,1)
    p
  }
})

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
    }else{
      pairs <- sapply(input$alphaPairs, strsplit, split=" vs. ")
      p <- suppressWarnings(ggboxplot(alphaTab, x=input$alphaGroup, y="value", fill=input$alphaGroup, facet.by = "measure",
                                     palette = input$alphaPalette, 
                                     scales=input$alphaScalesFree)+
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
      ggsave(file, alphaReact()$gg, device="pdf", width = 10, height = 7)
    }
  }
)

####beta diversity####
#do calculation for beta diversity plots
betaReactive <- reactive({
  if(!is.null(currentSet())){
    waiter_show(html = tagList(spin_rotating_plane(),"Doing calculation ... "),color=overlay_color)
    
    group <- input$betaGroup
    phylo <- vals$datasets[[currentSet()]]$phylo
    
    meta <- as.data.frame(sample_data(phylo))
    meta <- data.frame(meta[order(rownames(meta)),])
    # remove all samples with NA in column of chosen group
    meta <- meta[!is.na(meta[[group]]),]
    samples <- meta[[sample_column]]
    phylo <- prune_samples(samples, phylo)
    
    if(input$betaLevel != "All"){
      keep_samples <- meta[meta[[group]]==input$betaLevel,][["SampleID"]]
      phylo <- prune_samples(keep_samples, phylo)
      meta <- as.data.frame(sample_data(phylo))
    }
    
    group_vector <- as.factor(meta[[group]])
    
    if(!is.null(access(phylo,"phy_tree"))) tree <- phy_tree(phylo) else tree <- NULL
    
    method <- match(input$betaMethod,c("Bray-Curtis Dissimilarity","Generalized UniFrac Distance", "Unweighted UniFrac Distance", "Weighted UniFrac Distance", "Variance adjusted weighted UniFrac Distance"))
    my_dist <- betaDiversity(phylo, method=method)
    
    all_fit <- hclust(my_dist,method="ward.D2")
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

    mds <- data.frame(cmdscale(my_dist,k=2))
    mds$group <- group_vector
    mds$method <- "multidimensional scaling (MDS)"
    mds$sample <- meta[[sample_column]]
    
    meta_mds <- data.frame(metaMDS(my_dist,k=2)[["points"]])
    colnames(meta_mds) <- c("X1","X2")
    meta_mds$group <- group_vector
    meta_mds$method <- "non-metric multidimensional scaling (NMDS)"
    meta_mds$sample <- meta[[sample_column]]
    
    plot_df <- rbind(mds, meta_mds)
    show_ellipse <- length(unique(plot_df[["group"]]))>1
    pval_label <- grobTree(textGrob(paste("p-value: ",pval), x=0.05,  y=0.05, hjust=0, gp=gpar(col="black", fontsize=10)))
    
    if(input$betaShowLabels){
      beta_scatter <- ggscatter(plot_df, x="X1",y="X2", color="group", show.legend.text = F, 
                                ellipse = show_ellipse, ellipse.type = "confidence",mean.point = T, palette = input$betaPalette,
                                star.plot = T,facet.by = "method", ggtheme = theme_bw(), xlab="",ylab="",
                                label = "sample") 
      beta_scatter <- beta_scatter+annotation_custom(pval_label)
    }else{
      beta_scatter <- ggscatter(plot_df, x="X1",y="X2", color="group", show.legend.text = F, 
                                ellipse = show_ellipse, ellipse.type = "confidence",mean.point = T, palette = input$betaPalette,
                                star.plot = T,facet.by = "method", ggtheme = theme_bw(), xlab="",ylab="")
      beta_scatter <- beta_scatter+annotation_custom(pval_label)
    }
    
    out <- list(dist=my_dist, all_groups=group_vector, tree=tree, pval=pval, beta_scatter=beta_scatter)
    waiter_hide()
    return(out)
  }
})

betaColReactive <- reactive({
  if(!is.null(betaReactive())){
    group_vector <- betaReactive()$all_groups
    col = switch (input$betaPalette,
                  "jco" = pal_jco("default")(length(levels(group_vector)))[group_vector],
                  "rainbow" = rainbow(length(levels(group_vector)))[group_vector],
                  "npg" = pal_npg("nrc")(length(levels(group_vector)))[group_vector],
                  "aaas" = pal_aaas("default")(length(levels(group_vector)))[group_vector],
                  "nejm" = pal_nejm("default")(length(levels(group_vector)))[group_vector],
                  "lancet" = pal_lancet("lanonc")(length(levels(group_vector)))[group_vector],
                  "jama" = pal_jama("default")(length(levels(group_vector)))[group_vector],
                  "ucscgb" = pal_ucscgb("default")(length(levels(group_vector)))[group_vector]
    )
    col
  }
})

# clustering tree of samples based on beta-diversity
output$betaTree <- renderPlot({
  if(!is.null(betaReactive())){
    beta <- betaReactive()
    plot(beta$tree,type="phylogram",use.edge.length=T,tip.color=betaColReactive(),label.offset=0.01)
    axisPhylo()
    tiplabels(pch=16,col=beta$col)
  }
})

#download as pdf
output$betaTreePDF <- downloadHandler(
  filename = function(){"beta_diversity_clustering.pdf"},
  content = function(file){
    if(!is.null(betaReactive())){
      pdf(file, width=8, height=10)
      plot(betaReactive()$tree,type="phylogram",use.edge.length=T,tip.color=betaColReactive(),label.offset=0.01, cex=0.7)
      axisPhylo()
      tiplabels(pch=16,col=betaReactive()$col)
      dev.off()
    }
  }
)

output$betaDivScatterInteractive <- renderPlotly({
  if(!is.null(betaReactive())){
    ggplotly(betaReactive()$beta_scatter,tooltip = c(),height = 600) 
  }
})

output$betaDivScatter <- renderPlot({
  if(!is.null(betaReactive())){
    betaReactive()$beta_scatter
  }
})

output$betaDivScatterPDF <- downloadHandler(
  filename = function(){"beta_diversity_scatter.pdf"},
  content=function(file){
    if(!is.null(betaReactive())){
      ggsave(file, betaReactive()$beta_scatter, device="pdf", width = 10, height = 10)
      #ggsave(file, taxBinningReact()$gg, device="pdf", width = 10, height = 7)
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
    #phylo <- transform_sample_counts(phylo, function(x) x/sum(x)*100) #relative abunance
    #check for unifrac distance --> (needs phylo tree file):
    l <- list()
    if(!is.null(vals$datasets[[currentSet()]]$unifrac_dist)){
      #save generalized unifrac distance as global variable to use it for heatmap
      gunifrac_heatmap <<- as.dist(vals$datasets[[currentSet()]]$unifrac_dist)
      hm_distance <- if(input$heatmapDistance == "gunifrac") "gunifrac_heatmap" else input$heatmapDistance
      if(input$heatmapOrderSamples){
        p<-plot_heatmap(phylo,method = input$heatmapOrdination, distance = hm_distance, sample.label = input$heatmapSample, sample.order = input$heatmapSample)
      }else{
        p<-plot_heatmap(phylo,method = input$heatmapOrdination, distance = hm_distance, sample.label = input$heatmapSample)
      }
      l <- list(gg=p)
    }
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

