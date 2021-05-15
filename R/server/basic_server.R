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
output$rarefacCurve <- renderPlotly({
  if(!is.null(currentSet())){
    
    #needs integer values to work
    tab = as.matrix(vals$datasets[[currentSet()]]$rawData)
    class(tab)<-"integer"
    
    rarefactionCurve = lapply(1:ncol(tab),function(i){
      n = seq(1,colSums(tab)[i],by=20)
      if(n[length(n)]!=colSums(tab)[i]) n=c(n,colSums(tab)[i])
      drop(rarefy(t(tab[,i]),n))
    })
    
    #slope = apply(tab,2,function(x) rareslope(x,sum(x)-100))
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
    
    #create phyloseq-object with relative abundance data
    #otu_table(phylo) <- otu_table(rel_dat,T)
    rel_phylo <- merge_phyloseq(otu_table(rel_dat,T),tax_table(phylo))
    tax_binning <- taxBinningNew(if(input$taxaAbundanceType)rel_phylo else phylo, vals$datasets[[currentSet()]]$is_fastq)
    tax_binning
  }
})

# plot distribution of taxa
output$taxaDistribution <- renderPlotly({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$is_fastq){
      tab= taxBinningReact()[[which(c("Kingdom","Phylum","Class","Order","Family","Genus")==input$filterTaxa)]]
    }else{
      tab = taxBinningReact()[[which(c("Kingdom","Phylum","Class","Order","Family","Genus","Species")==input$filterTaxa)]] 
    }
    
    #taxa = ifelse(rowSums(tab)/ncol(tab)<(input$taxCutoff),"Other",rownames(tab))
    taxa = rownames(tab)
    if(any(taxa=="Other")){
      other <- tab[which(taxa=="Other"),]
      other <- colSums(other)
      
      tab <- tab[-which(taxa=="Other"),]
      tab <- rbind(tab,Other=other)
    }
    
    tab = melt(tab)
    tab$Var2 <- as.character(tab$Var2)
    
    plot_ly(tab,name=~Var1,x=~value,y=~Var2,type="bar",orientation="h") %>%
      layout(xaxis=list(title=ifelse(input$taxaAbundanceType,"Relative Abundance", "Absolute Abundance / counts")),yaxis=list(title="Samples"),
             barmode="stack",showlegend=T) #,legend=list(orientation="h")
  } else plotly_empty()
  
})

####dimensionality reduction (PCA, UMAP, tSNE)####
structureReact <- reactive({
  if(!is.null(currentSet())){
    mat <- as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo))
    mat_t <- t(mat)
    
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
    tsne = Rtsne(t(mat),dim=3,perplexity=min((length(samples)-1)/3,30))
    out_tsne = data.frame(tsne$Y,txt=samples)
    
    l <- list(percentage=percentage, loadings=loadings, out_pca = out_pca, raw_pca = pca, out_umap=out_umap, out_tsne=out_tsne)
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
    #meta <- ifelse(vals$datasets[[currentSet()]]$has_meta, sample_data(vals$datasets[[currentSet()]]$phylo), NULL) 
    
    alphaTab = data.frame(colnames(otu))
    for(i in c("Shannon_Entropy","effective_Shannon_Entropy","Simpson_Index","effective_Simpson_Index","Richness")){
      alphaTab = cbind(alphaTab,round(alphaDiv(otu,i),2))
    }
    colnames(alphaTab) = c("SampleID","Shannon_Entropy","effective_Shannon_Entropy","Simpson_Index","effective_Simpson_Index","Richness")
    #metaColumn <- as.factor(meta[[input$alphaGroup]])
    alphaTab
  }else{
    NULL
  }
})

# plot alpha diversity
output$alphaPlot <- renderPlotly({
  if(!is.null(alphaReact())){
    otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
    if(vals$datasets[[currentSet()]]$has_meta){meta<-sample_data(vals$datasets[[currentSet()]]$phylo)} 
    
    alphaTab = data.frame(alphaReact()[,c(input$alphaMethod, sample_column)])
    colnames(alphaTab) <- c(input$alphaMethod, sample_column)
    
    if(input$alphaGroup=="-") plot_ly(data=alphaTab, 
                                      y=as.formula(paste0("~ ",input$alphaMethod)),
                                      type='violin',
                                      box=list(visible=T),
                                      meanline=list(visible=T),
                                      x0=input$alphaMethod, 
                                      hoverinfo="text",
                                      text=paste0(as.character(input$alphaMethod),": ",alphaTab[[input$alphaMethod]]," Sample: ", alphaTab[[sample_column]]),
                                      points="all") %>% layout(yaxis=list(title="alpha Diversity",zeroline=F))
    else plot_ly(data=alphaTab, 
                 y=as.formula(paste0("~",input$alphaMethod)),
                 x=meta[[input$alphaGroup]],
                 color=meta[[input$alphaGroup]],
                 type='violin',
                 box=list(visible=T),
                 meanline=list(visible=T),
                 x0=input$alphaMethod, 
                 hoverinfo="text",
                 text=paste0(as.character(input$alphaMethod),": ",alphaTab[[input$alphaMethod]]," Sample: ", alphaTab[[sample_column]], " Group: ", meta[[input$alphaGroup]]),
                 points="all") %>% layout(yaxis=list(title="alpha Diversity",zeroline=F))
  }
  else plotly_empty()
})

# show alpha diversity table
output$alphaTable <- renderTable({
  if(!is.null(alphaReact())){
    alphaReact()
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

####explained variation####
#reactive table of explained variation; for each meta-variable calculate p-val and rsquare
explVarReact <- reactive({
  if(!is.null(currentSet())){
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
      withProgress(message="Calculating confounding factors...",value=0,{
        vals$datasets[[currentSet()]]$confounder_table <- calculateConfounderTable(var_to_test=input$confounding_var,
                                                                                   variables = meta,
                                                                                   distance=vals$datasets[[currentSet()]]$unifrac_dist,
                                                                                   seed=seed,
                                                                                   progress=T)
      })
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
    group <- input$betaGroup
    phylo <- vals$datasets[[currentSet()]]$phylo
    
    meta <- as.data.frame(sample_data(phylo))
    meta <- data.frame(meta[order(rownames(meta)),])
    meta_pos <- which(colnames(meta) == group)
    all_groups <- as.factor(meta[,meta_pos])
    
    if(!is.null(access(phylo,"phy_tree"))) tree <- phy_tree(phylo) else tree <- NULL
    
    method <- match(input$betaMethod,c("Bray-Curtis Dissimilarity","Generalized UniFrac Distance", "Unweighted UniFrac Distance", "Weighted UniFrac Distance", "Variance adjusted weighted UniFrac Distance"))
    my_dist <- betaDiversity(phylo, method=method)
    
    all_fit <- hclust(my_dist,method="ward.D2")
    tree <- as.phylo(all_fit)
    
    #adonis <- adonis2(as.formula(paste0("my_dist ~ ",eval("all_groups")), env = environment()))
    #pval <- adonis[["Pr(>F)"]][1]
    
    col = rainbow(length(levels(all_groups)))[all_groups]
    
    out <- list(dist=my_dist, col=col, all_groups=all_groups, tree=tree)
    return(out)
  }
})

# clustering tree of samples based on beta-diversity
output$betaTree <- renderPlot({
  if(!is.null(betaReactive())){
    beta <- betaReactive()
    plot(beta$tree,type="phylogram",use.edge.length=T,tip.color=beta$col,label.offset=0.01)
    axisPhylo()
    tiplabels(pch=16,col=beta$col)
  }
})

# MDS plot based on beta-diversity
output$betaMDS <- renderPlot({
  if(!is.null(betaReactive())){
    beta <- betaReactive()
    mds <- cmdscale(beta$dist,k=2)
    samples<-row.names(mds)
    s.class(
      mds,col=unique(beta$col),cpoint=2,fac=beta$all_groups,
      sub=paste("MDS plot of Microbial Profiles")
    )
    if(input$betaShowLabels){
      text(mds,labels=samples,cex=0.7,adj = c(-.1,-.8))
    }
  }
})

# NMDS plot based on beta-diversity
output$betaNMDS <- renderPlot({
  if(!is.null(betaReactive())){
    beta<-betaReactive()
    meta_mds = metaMDS(beta$dist,k=2)
    samples = row.names(meta_mds$points)
    s.class(
      meta_mds$points,col=unique(beta$col),cpoint=2,fac=beta$all_groups,
      sub=paste("metaNMDS plot of Microbial Profiles")
    )
    if(input$betaShowLabels){
      text(meta_mds$points,labels=samples,cex=0.7,adj = c(-.1,-.8),offset = .1)
    }
  }
})


####phylogenetic tree####
output$phyloTree <- renderPlot({
  if(!is.null(currentSet())){
    phylo <- vals$datasets[[currentSet()]]$phylo
    #need to look for tree object like this, did not find out better/cleaner way yet
    if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
    
    #only can display a phylogenetic tree if tree object is present
    if(!is.null(tree)){
      #prune taxa to number the user sets with slider (changes size of tree)
      myTaxa = names(sort(taxa_sums(phylo), decreasing = TRUE)[1:input$phylo_prune])
      pruned_phylo <- prune_taxa(myTaxa, phylo)
      
      if(input$phylo_color == "-") phylo_color = NULL else phylo_color=input$phylo_color
      if(input$phylo_shape == "-") phylo_shape = NULL else phylo_shape=input$phylo_shape
      if(input$phylo_size == "-") phylo_size = NULL else phylo_size=input$phylo_size
      if(input$phylo_tiplabels == "-") phylo_label.tips = NULL else phylo_label.tips=input$phylo_tiplabels
      
      if(input$phylo_radial == T){
        plot_tree(pruned_phylo,method = input$phylo_method,color=phylo_color,shape = phylo_shape,size = phylo_size,label.tips = phylo_label.tips,ladderize = "left", plot.margin = 0.1)+coord_polar(theta = "y")
      }else{
        plot_tree(pruned_phylo,method = input$phylo_method,color=phylo_color,shape = phylo_shape,size = phylo_size,label.tips = phylo_label.tips,ladderize = input$phylo_ladderize, plot.margin = 0.1)
      }
      
    }
  }
}, height=800)


#javascript show/hide toggle for advanced options
shinyjs::onclick("phylo_toggle_advanced",shinyjs::toggle(id="phylo_advanced",anim = T))

####associations####

observeEvent(input$associations_start,{
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_meta){
      message(paste0(Sys.time(), " - building SIAMCAT object ..."))
      meta <- vals$datasets[[currentSet()]]$metaData
      rel_otu <- vals$datasets[[currentSet()]]$relativeData
      meta <- data.frame(t(na.omit(t(meta))))
      
      s.obj <- siamcat(feat=rel_otu, meta=meta, label= input$associations_label, case= input$associations_case)
      s.obj <- filter.features(s.obj)
      vals$datasets[[currentSet()]]$siamcat <- s.obj
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
