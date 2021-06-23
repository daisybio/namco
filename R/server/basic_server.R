# update targets table of the currently loaded dataset
output$metaTable <- renderDataTable({
  if(!is.null(currentSet())){
    #datatable(sample_data(vals$datasets[[currentSet()]]$phylo),filter='top',options=list(searching=T,pageLength=20,dom="Blfrtip",scrollX=T),editable=F,rownames=F)
    meta_dt<-datatable(sample_data(vals$datasets[[currentSet()]]$phylo),filter='top',selection=list(mode='multiple'),options=list(pageLength=20,scrollX=T))
    meta_dt
  } 
  else datatable(data.frame(),options=list(dom="t"))
},server=T)

####rarefaction curves####
rarefactionReactive <- reactive({l
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
    
    waiter_show(html = tagList(spin_rotating_plane(),"Preparing plot ... "),color=overlay_color)
    #create phyloseq-object with relative abundance data
    #otu_table(phylo) <- otu_table(rel_dat,T)
    rel_phylo <- merge_phyloseq(otu_table(rel_dat,T),tax_table(phylo))
    tax_binning <- taxBinningNew(if(input$taxaAbundanceType)rel_phylo else phylo, vals$datasets[[currentSet()]]$is_fastq)

    if(vals$datasets[[currentSet()]]$is_fastq){
      tab= tax_binning[[which(c("Kingdom","Phylum","Class","Order","Family","Genus")==input$filterTaxa)]]
    }else{
      tab = tax_binning[[which(c("Kingdom","Phylum","Class","Order","Family","Genus","Species")==input$filterTaxa)]] 
    }
    
    if(vals$datasets[[currentSet()]]$has_meta){
      meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo), check.names = F)
      tab <- merge(melt(tab), meta, by.x = "Var2", by.y=sample_column, all.x=T)
    }else{
      tab <- melt(tab)
    }
    
    if(input$taxBinningGroup == "None"){
      p <- ggplot(tab, aes(x=value, y=Var2, fill=Var1))+
        geom_bar(stat="identity")+
        xlab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance / counts"))+
        ylab("Sample")+
        scale_fill_discrete(name = input$filterTaxa)+
        ggtitle(paste0("Taxonomic Binning of samples"))
      
    }else{
      p <- ggplot(tab, aes(x=value, y=Var2, fill=Var1))+
        geom_bar(stat="identity")+
        facet_wrap(as.formula(paste0("~",input$taxBinningGroup)), scales = "free")+
        xlab(ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance / counts"))+
        ylab("Sample")+
        scale_fill_discrete(name = input$filterTaxa)+
        ggtitle(paste0("Taxonomic Binning, grouped by ", input$taxBinningGroup))
    }
    if(input$taxBinningShowNames == "No"){
      p <- p + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
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
      ggsave(file, taxBinningReact()$gg, device="pdf", width = 10, height = 7)
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
    
    alphaTabFull = data.frame(colnames(otu))
    for(i in c("Shannon_Entropy","effective_Shannon_Entropy","Simpson_Index","effective_Simpson_Index","Richness")){
      alphaTabFull = cbind(alphaTabFull,round(alphaDiv(otu,i),2))
    }
    colnames(alphaTabFull) = c("SampleID","Shannon_Entropy","effective_Shannon_Entropy","Simpson_Index","effective_Simpson_Index","Richness")
    #metaColumn <- as.factor(meta[[input$alphaGroup]])
    
    if(vals$datasets[[currentSet()]]$has_meta){
      meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo), check.names = F)
      alphaTabFull <- merge(alphaTabFull, meta, by.x="SampleID", by.y="SampleID")
    } 
    
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
      write.csv(alphaReact(),file,row.names = F)
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

####explained variation####
#reactive table of explained variation; for each meta-variable calculate p-val and rsquare
explVarReact <- reactive({
  if(!is.null(currentSet())){
    
    waiter_show(html = tagList(spin_rotating_plane(),"Doing calculation ... "),color=overlay_color)
    
    OTUs <- data.frame(t(otu_table(vals$datasets[[currentSet()]]$phylo)))  #transposed otu-table (--> rows are samples, OTUs are columns)
    meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
    
    #alpha<-alphaReact()
    #alpha[[sample_column]] <- NULL
    #variables <- cbind.as.data.frame(meta[rownames(OTUs),],alpha[rownames(OTUs),])
    variables <- data.frame(meta[rownames(OTUs),])
    variables[[sample_column]] <- NULL
    
    
    plist <- vector()
    rlist <- vector()
    namelist <- vector()
    #iterate over all columns
    for (i in 1:dim(variables)[2]) {
      if(length(unique(variables[,i])) > 1){
        variables_nc <- completeFun(variables,i)
        colnames(variables_nc) <- colnames(variables)
        var_name <- colnames(variables)[i]
        #calculate distance matrix between OTUs (bray-curtis)
        BC <- vegdist(OTUs[which(row.names(OTUs) %in% row.names(variables_nc)),], method="bray")
        output <- adonis2(as.formula(paste0("BC ~ ",var_name)), data = variables_nc)
        pvalue <- output[["Pr(>F)"]][1]
        rsquare <- output[["R2"]][1]
        names <- names(variables_nc)[i]
        
        plist <- append(plist,pvalue)
        rlist <- append(rlist,rsquare)
        namelist <- append(namelist,names)
      }
    }
    
    df <- data.frame(Variable = namelist, pvalue = plist, rsquare = rlist)
    waiter_hide()
    df
    
    
  }else{
    NULL
  }
})


output$explainedVariation <- renderTable({
  if(!is.null(explVarReact())){
    explVarReact()
  }
})

output$explainedVariationBar <- renderPlot({
  if(!is.null(explVarReact())){
    explVar <- explVarReact()
    explVar$Variable <- factor(explVar$Variable, levels = unique(explVar$Variable)[order(explVar$rsquare,decreasing = T)])
    ggplot(data=explVar,aes(x=Variable,y=rsquare))+
      geom_bar(stat = "identity",aes(fill=pvalue))+
      ggtitle("Explained Variation of meta variables")+
      theme(axis.text.x = element_text(angle = 90))+
      geom_text(aes(label=round(rsquare,digits = 4)),vjust=-1)
  }
})

####confounding factors####
observeEvent(input$confounding_start,{
  if(!is.null(currentSet())){
    #only if pre-build distance matrix exists, this can be calcualted (depending on tree input)
    if (!is.null(vals$datasets[[currentSet()]]$unifrac_dist)){
      meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
      #remove first column --> SampleID column
      meta[[sample_column]]<-NULL
      
      #calulate confounding matrix
      waiter_show(html = tagList(spin_rotating_plane(),"Doing calculation ... "),color=overlay_color)
      vals$datasets[[currentSet()]]$confounder_table <- calculateConfounderTable(var_to_test=input$confounding_var,
                                                                                 variables = meta,
                                                                                 distance=vals$datasets[[currentSet()]]$unifrac_dist,
                                                                                 seed=seed,
                                                                                 progress=F)
      waiter_hide()
    }
  }
})

output$confounding_table <- renderTable({
  if(!is.null(currentSet())){
    if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
      vals$datasets[[currentSet()]]$confounder_table$table
    }
  }
})

output$confounding_table_download <- downloadHandler(
  filename = function(){
    paste("confounding_factors.csv")
  },
  content = function(file){
    if(!is.null(currentSet())){
      if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
        write.csv(vals$datasets[[currentSet()]]$confounder_table$table,file,row.names = F)
      }
    }
  }
)

output$confounding_var_text <- renderUI({
  if(!is.null(currentSet())){
    if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
      HTML(paste0("<b> Chosen variable: </b> ",vals$datasets[[currentSet()]]$confounder_table$var))
    }
  }
})

####beta diversity####
#do calculation for beta diversity plots
betaReactive <- reactive({
  if(!is.null(currentSet())){
    waiter_show(html = tagList(spin_rotating_plane(),"Doing calculation ... "),color=overlay_color)
    
    group <- input$betaGroup
    phylo <- vals$datasets[[currentSet()]]$phylo
    
    meta <- as.data.frame(sample_data(phylo))
    meta <- data.frame(meta[order(rownames(meta)),])
    
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

    if(input$betaLevel != "All"){
      pval <- 0
    }else{
      colnames(meta)[which(colnames(meta) == group)] <- "condition"
      tryCatch({
        adonis <- adonis2(my_dist ~ condition, data=meta, parallel = ncores)
        pval <- adonis[["Pr(>F)"]][1]
      }, error=function(e){
        print(e$message)
        message(paste0("Error with adonis2 at beta-diversity: ", group))
        pval <- NULL
      })
    }

    mds <- cmdscale(my_dist,k=2)
    meta_mds <- metaMDS(my_dist,k=2)
    
    out <- list(dist=my_dist, all_groups=group_vector, tree=tree, pval=pval, mds=mds, meta_mds=meta_mds)
    waiter_hide()
    return(out)
  }
})

betaColReactive <- reactive({
  if(!is.null(betaReactive())){
    group_vector <- betaReactive()$all_groups
    col = switch (input$betaPalette,
                  "JCO" = pal_jco("default")(length(levels(group_vector)))[group_vector],
                  "Rainbow" = rainbow(length(levels(group_vector)))[group_vector],
                  "NPG" = pal_npg("nrc")(length(levels(group_vector)))[group_vector],
                  "AAAS" = pal_aaas("default")(length(levels(group_vector)))[group_vector],
                  "NEJM" = pal_nejm("default")(length(levels(group_vector)))[group_vector],
                  "Lancet" = pal_lancet("lanonc")(length(levels(group_vector)))[group_vector],
                  "JAMA" = pal_jama("default")(length(levels(group_vector)))[group_vector],
                  "UCSCGB" = pal_ucscgb("default")(length(levels(group_vector)))[group_vector]
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
      pdf(file, width=8, height=6)
      plot(betaReactive()$tree,type="phylogram",use.edge.length=T,tip.color=betaColReactive(),label.offset=0.01)
      axisPhylo()
      tiplabels(pch=16,col=betaReactive()$col)
      dev.off()
    }
  }
)

# MDS plot based on beta-diversity
output$betaMDS <- renderPlot({
  if(!is.null(betaReactive())){
    beta <- betaReactive()
    mds <- beta$mds
    samples<-row.names(mds)
    s.class(
      mds,col=unique(betaColReactive()),cpoint=2,fac=beta$all_groups,
      sub=paste("MDS plot of Microbial Profiles; pvalue:", beta$pval)
    )
    if(input$betaShowLabels){
      text(mds,labels=samples,cex=0.7,adj = c(-.1,-.8))
    }
  }
})

#download as pdf
output$betaMDSPDF <- downloadHandler(
  filename = function(){"beta_diversity_MDS.pdf"},
  content = function(file){
    if(!is.null(betaReactive())){
      pdf(file, width=8, height=6)
      samples<-row.names(betaReactive()$mds)
      s.class(
        mds,col=unique(betaColReactive()),cpoint=2,fac=betaReactive()$all_groups,
        sub=paste("MDS plot of Microbial Profiles; pvalue:", betaReactive()$pval)
      )
      if(input$betaShowLabels){
        text(mds,labels=samples,cex=0.7,adj = c(-.1,-.8))
      }
      dev.off()
    }
  }
)

# NMDS plot based on beta-diversity
output$betaNMDS <- renderPlot({
  if(!is.null(betaReactive())){
    beta<-betaReactive()
    meta_mds = beta$meta_mds
    samples = row.names(meta_mds$points)
    s.class(
      meta_mds$points,col=unique(betaColReactive()),cpoint=2,fac=beta$all_groups,
      sub=paste("metaNMDS plot of Microbial Profiles; pvalue:", beta$pval)
    )
    if(input$betaShowLabels){
      text(meta_mds$points,labels=samples,cex=0.7,adj = c(-.1,-.8),offset = .1)
    }
  }
})

#download as pdf
output$betaNMDSPDF <- downloadHandler(
  filename = function(){"beta_diversity_metaNMDS.pdf"},
  content = function(file){
    if(!is.null(betaReactive())){
      pdf(file, width=8, height=6)
      meta_mds = betaReactive()$meta_mds
      samples = row.names(meta_mds$points)
      s.class(
        meta_mds$points,col=unique(betaColReactive()),cpoint=2,fac=betaReactive()$all_groups,
        sub=paste("metaNMDS plot of Microbial Profiles; pvalue:", betaReactive()$pval)
      )
      if(input$betaShowLabels){
        text(meta_mds$points,labels=samples,cex=0.7,adj = c(-.1,-.8),offset = .1)
      }
      dev.off()
    }
  }
)

####phylogenetic tree####

treeReactive <- reactive({
  if(!is.null(currentSet())){
    # prune taxa
    myTaxa = names(sort(taxa_sums(vals$datasets[[currentSet()]]$phylo), decreasing = TRUE)[1:input$phylo_prune])
    phy <- prune_taxa(myTaxa, vals$datasets[[currentSet()]]$phylo)
    
    #phy <- vals$datasets[[currentSet()]]$phylo
    meta <- as.data.frame(sample_data(phy))
    otu <- as.data.frame(otu_table(phy))
    taxonomy <- as.data.frame(tax_table(phy))
    if(!is.null(access(phy,"phy_tree"))) tree <- phy_tree(phy) else tree <- NULL
    if(!is.null(tree)){
      if(input$phylo_group != "NONE"){
        group <- input$phylo_group
        # count number of occurrences of the OTUs in each sample group
        l<-lapply(unique(meta[[group]]), function(x){
          samples_in_group <- meta[["SampleID"]][as.character(meta[[group]])==as.character(x)]
          d<-data.frame(otu[,samples_in_group])
          d<-data.frame(rowSums(apply(d,2,function(x) ifelse(x>0,1,0))))
          colnames(d) <- c(as.character(x))
          return(d)
        })
        info <- merge(data.frame(l), taxonomy,by.x=0, by.y=0)
        rownames(info) <- info$Row.names
        info$Row.names <- NULL
        group_cols <- suppressWarnings(which(colnames(info)==unique(meta[[group]])))
      }else{
        info <- taxonomy
        group_cols <- c()
      }
      if(input$phylo_taxonomy != "NONE"){
        # collect, which columns contain the info for the heatmap
        taxa_cols <- suppressWarnings(which(colnames(info)==input$phylo_taxonomy))
      }else{
        taxa_cols <- c()
      }
      
      tree_plot <- suppressWarnings(suppressMessages(ggtree::ggtree(tree, layout = input$phylo_method, branch.length = input$phylo_draw_clado) %<+% info +
                                                      geom_tiplab(size=input$phylo_size_tree,
                                                                  align = ifelse(input$phylo_edge_length=="Yes",T,F))))
      
      return(list(tree_plot = tree_plot,
                  info = info,
                  group_cols = group_cols,
                  taxa_cols = taxa_cols))
    }
  }
})


output$phyloTree <- renderPlot({
  if(!is.null(treeReactive())){
    h <- treeReactive()$tree_plot
    info <- treeReactive()$info
    if(!is.null(treeReactive()$group_cols)){
      h<-suppressWarnings(suppressMessages(ggtree::gheatmap(h, info[treeReactive()$group_cols], 
                                                            offset=input$phylo_offset,
                                                            color=NULL, 
                                                            width=input$phylo_width_meta,
                                                            colnames_position="top", 
                                                            colnames_angle=90, colnames_offset_y = 5, 
                                                            hjust=1, font.size=3,low="white")))
      h <- h + new_scale_fill()   
    }
    if(!is.null(treeReactive()$taxa_cols)){
      h<-suppressWarnings(suppressMessages(ggtree::gheatmap(h, info[treeReactive()$taxa_cols],
                                                            width=input$phylo_width_taxonomy,
                                                            offset=input$phylo_offset+5,
                                                            color="black",
                                                            colnames=T,
                                                            colnames_position="top",
                                                            colnames_angle=90, colnames_offset_y = 5,
                                                            hjust=1, font.size = 3)))  
    }
    h
  }
}, height = 1000)


#javascript show/hide toggle for advanced options
shinyjs::onclick("phylo_toggle_advanced",shinyjs::toggle(id="phylo_advanced",anim = T))

####associations####

observeEvent(input$associations_start,{
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_meta){
      message(paste0(Sys.time(), " - building SIAMCAT object ..."))
      waiter_show(html = tagList(spin_rotating_plane(),"Calculating differential associations ..."),color=overlay_color)
      
      phylo <- vals$datasets[[currentSet()]]$phylo
      
      if(input$associations_level != "OTU"){
        phylo_glom <- glom_taxa_custom(phylo, input$associations_level) #merge OTUs with same taxonomic level
        phylo <- phylo_glom$phylo_rank
        taxa_names(phylo) <- phylo_glom$taxtab[[input$associations_level]]
        rel_otu <- relAbundance(data.frame(otu_table(phylo), check.names=F))
      }else{
        rel_otu <- vals$datasets[[currentSet()]]$relativeData
      }
      
      meta <- data.frame(sample_data(phylo))
      meta <- data.frame(t(na.omit(t(meta))))
      
      tryCatch({
        # siamcat works only with 5 or more samples (https://git.embl.de/grp-zeller/SIAMCAT/-/blob/a1c662f343e99dabad4de024d4c993deba91bb0c/R/validate_data.r#L82)
        if(sum(meta[[input$associations_label]]==input$associations_case) <= 5){stop(siamcatNotEnoughSamplesError, call.=F)}
        
        s.obj <- siamcat(feat=rel_otu, meta=meta, label= input$associations_label, case= input$associations_case)
        s.obj <- filter.features(s.obj)
        vals$datasets[[currentSet()]]$siamcat <- s.obj
      },error = function(e){
        print(e$message)
        showModal(errorModal(e$message))
      })

      waiter_hide()
    }
  }
})

output$associationsPlot <- renderPlot({
  if(!is.null(currentSet())){
    if(!is.null(vals$datasets[[currentSet()]]$siamcat)){
      s.obj <- vals$datasets[[currentSet()]]$siamcat 
      sort.by <- c("p.val","fc","pr.shift")[which(input$associations_sort==c("p-value","fold-change","prevalence shift"))]
      panels <- c("fc","auroc","prevalence")[which(input$associations_panels==c("fold-change","AU-ROC","prevalence"))]
      suppressMessages(check.associations(s.obj, fn.plot = NULL, prompt=F, verbose=0,
                         alpha = input$associations_alpha, 
                         max.show = input$assiciation_show_numer, 
                         sort.by = sort.by,
                         panels = panels))
    }
  }
}, height=800)
