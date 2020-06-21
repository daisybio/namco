library(ade4)
library(cluster)
library(data.table)
library(DT)
library(fpc)
library(GUniFrac)
library(heatmaply)
library(networkD3)
library(klaR)
library(phangorn)
library(plotly)
library(RColorBrewer)
library(reshape2)
library(Rtsne)
library(shiny)
library(textshape)
library(tidyr)
library(umap)
library(themetagenomics)
library(SpiecEasi)
library(igraph)
library(Matrix)
library(phyloseq)
library(NbClust)

server <- function(input,output,session){
  options(shiny.maxRequestSize=1000*1024^2,stringsAsFactors=F)  # upload up to 1GB of data
  source("algorithms.R")
  source("utils.R")
  
  vals = reactiveValues(datasets=list(),undersampled=c()) # reactiveValues is a container for variables that might change during runtime and that influence one or more outputs, e.g. the currently selected dataset
  currentSet = NULL # a pointer to the currently selected dataset
  ncores = 4  # number of cores used where it is possible to use multiple
  
  ################################################################################################################################################
  #
  # file upload
  #
  ################################################################################################################################################
  
  
  # Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
  uploadModal <- function(failed=F,error_message=NULL) {
    modalDialog(
      p("Please provide pregenerated input files. For detailed information how the files have to look, check out the Info & Settings tab on the left!"),
      fluidRow(
        column(6,fileInput("otuFile","Select OTU table")),
        column(6,fileInput("metaFile","Select Metadata File"))
      ),
      fluidRow(
        column(6,fileInput("treeFile","Select Phylogenetic Tree File (optional)",width="100%")),
        column(3,actionButton("testdata","Load Testdata",style="background-color:red"))
      ),
      br(),
      fluidRow(
        column(1),
        column(3,checkboxInput("inputNormalized","Input is already normalized")),
        column(6,radioButtons("normMethod","Normalization Method",c("by Sampling Depth","by Rarefaction"),inline=T))
      ),
      br(),
      textInput("dataName","Enter a project name:",placeholder="New_Project",value="New_Project"),
      
      if(failed) {
        div(tags$b("The file you specified could not be loaded. Please check the Info tab and to irm your data is in the correct format!",style="color: red;"),
            tags$p(error_message))
      },
      footer = tagList(
        modalButton("Cancel"),
        actionButton("upload_ok","OK",style="background-color:blue")
      )
    )
  }
  
  # launch upload dialog
  observeEvent(input$upload, {
    showModal(uploadModal())
  })
  
  # try to load the dataset specified in the dialog window
  observeEvent(input$upload_ok, {
    # Check that data object exists and is data frame.
    if(!is.null(input$otuFile) & !input$dataName%in%names(vals$datasets)){
      tryCatch({
        dat <- read.csv(input$otuFile$datapath,header=T,sep="\t",row.names=1,check.names=F) # load data table
        taxonomy = generateTaxonomyTable(dat) # generate taxonomy table from TAX column
        dat = dat[!apply(is.na(dat)|dat=="",1,all),-ncol(dat)] # remove "empty" rows
        #read meta file; set rownames to names of first column; 
        if(is.null(input$metaFile)) meta = NULL else {meta = read.csv(input$metaFile$datapath,header=T,sep="\t"); rownames(meta)=meta[,1]; meta = meta[match(colnames(dat),meta$SampleID),]}
        if(is.null(input$treeFile)) tree = NULL else tree = read.tree(input$treeFile$datapath)
        normalized_dat = normalizeOTUTable(dat,which(input$normMethod==c("by Sampling Depth","by Rarefaction"))-1,input$inputNormalized)
        tax_binning = taxBinning(normalized_dat[[2]],taxonomy)
        
        #create phyloseq object from data (OTU, meta, taxonomic, tree)
        py.otu <- otu_table(normalized_dat$norm_tab,T)
        py.tax <- tax_table(as.matrix(taxonomy))
        py.meta <- sample_data(meta)
        #cannot build phyloseq object with NULL as tree input; have to check both cases:
        if (!is.null(tree)) phylo <- merge_phyloseq(py.otu,py.tax,py.meta, tree) else phylo <- merge_phyloseq(py.otu,py.tax,py.meta)
        
        #pre-build unifrac distance matrix
        if(!is.null(tree)) unifrac_dist <- buildDistanceMatrix(normalized_dat$norm_tab,meta,tree) else unifrac_dist <- NULL
        
        vals$datasets[[input$dataName]] <- list(rawData=dat,metaData=meta,taxonomy=taxonomy,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab,tree=tree,tax_binning=tax_binning,phylo=phylo,unifrac_dist=unifrac_dist,undersampled_removed=F)
        updateTabItems(session,"sidebar",selected="basics")
        removeModal()
      },
      error = function(e){
        #print(e)
        showModal(uploadModal(failed=T,e))
      })
    } else{
      e="OTU-file empty or name for project already used in current session. Please change name and check your files."
      showModal(uploadModal(failed=T))
    }
  })
  
  
  # upload test data
  observeEvent(input$testdata, {
    dat <- read.csv("testdata/OTU_table.tab",header=T,sep="\t",row.names=1) # load data table
    taxonomy = generateTaxonomyTable(dat) # generate taxonomy table from TAX column
    dat = dat[!apply(is.na(dat)|dat=="",1,all),-ncol(dat)] # remove "empty" rows
    meta = read.csv("testdata/metafile.tab",header=T,sep="\t")
    rownames(meta) = meta[,1]
    meta = meta[match(colnames(dat),meta$SampleID),]
    tree = read.tree("testdata/tree.tre") # load phylogenetic tree
    
    normalized_dat = normalizeOTUTable(dat,which(input$normMethod==c("by Sampling Depth","by Rarefaction"))-1,F)
    tax_binning = taxBinning(normalized_dat[[2]],taxonomy)
    
    #create phyloseq object from data (OTU, meta, taxonomic, tree)
    py.otu <- otu_table(normalized_dat$norm_tab,T)
    py.tax <- tax_table(as.matrix(taxonomy))
    py.meta <- sample_data(meta)
    phylo <- merge_phyloseq(py.otu,py.tax,py.meta,tree)
    
    #pre-build unifrac distance matrix
    if(!is.null(tree)) unifrac_dist <- buildDistanceMatrix(normalized_dat$norm_tab,meta,tree) else unifrac_dist <- NULL
    
    vals$datasets[["Mueller et al."]] <- list(rawData=dat,metaData=meta,taxonomy=taxonomy,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab,tree=tree,tax_binning=tax_binning,phylo=phylo,unifrac_dist=unifrac_dist,undersampled_removed=F)
    updateTabItems(session,"sidebar",selected="basics")
    removeModal()
    
  })
  
  # update datatable holding currently loaded datasets
  output$datasets <- renderDataTable({
    if(length(vals$datasets)!=0){
      datatable(data.frame(Datasets=names(vals$datasets)),rownames=F,options=list(pageLength=10,dom="t"),selection=list(mode="single",selected=length(vals$datasets))) %>%
        formatStyle("Datasets",color="white",backgroundColor="#222D33")
    } else{
      datatable(data.frame(Datasets=""),rownames=F,options=list(pageLength=10,dom="t"),selection=list(mode="single")) %>%
        formatStyle("Datasets",color="white",backgroundColor="#222D33")
    }
  })
  
  # choose current dataset; return NULL if no set is yet uploaded
  currentSet <- eventReactive(input$datasets_rows_selected, {
    if(length(vals$datasets) == 0){
      return (NULL)
    }
    return(input$datasets_rows_selected)
  })
  
  # update input selections
  observe({
    if(!is.null(currentSet())){  
      #get tables from phyloseq object
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      #if(!is.null(phy_tree(vals$datasets[[currentSet()]]$phylo))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      phylo <- vals$datasets[[currentSet()]]$phylo
      
      updateSliderInput(session,"rareToShow",min=1,max=ncol(otu),value=min(50,ncol(otu)))
      
      #update silder for binarization cutoff dynamically based on normalized dataset
      min_value <- min(otu)
      max_value <- round(max(otu)/16)
      updateNumericInput(session,"binCutoff",min=min_value,max=max_value)
      updateNumericInput(session,"k_in",value=0,min=0,max=vals$datasets[[currentSet()]]$vis_out$K,step=1)
      
      ########updates based on meta info########
      covariates <- vals$datasets[[currentSet()]]$vis_out$covariates
      updateSelectInput(session,"choose",choices = covariates)
      
      #pick all column names, except the SampleID
      group_columns <- setdiff(colnames(meta),"SampleID")
      updateSelectInput(session,"alphaGroup",choices = c("-",group_columns))
      updateSelectInput(session,"betaGroup",choices = group_columns)
      updateSelectInput(session,"structureGroup",choices = group_columns)
      updateSelectInput(session,"groupCol",choices = group_columns)
      updateSelectInput(session,"formula",choices = group_columns)
      
      
      if(is.null(access(phylo,"phy_tree"))) betaChoices="Bray-Curtis Dissimilarity" else betaChoices=c("Bray-Curtis Dissimilarity","Generalized UniFrac Distance")
      updateSelectInput(session,"betaMethod",choices=betaChoices)
      
      
      updateSelectInput(session,"phylo_color",choices= c("-",group_columns,"Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      updateSelectInput(session,"phylo_shape",choices = c("-",group_columns))
      updateSelectInput(session,"phylo_size",choices = c("-","abundance",group_columns))
      updateSliderInput(session,"phylo_prune",min=1,max=ntaxa(phylo),value=50,step=1)
    }
    
  })
  
  observe({
    if(!is.null(currentSet())){
      #factorize meta data
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      meta[] <- lapply(meta,factor)
      #pick all variables which have 2 or more factors for possible variables for confounding!
      tmp <- names(sapply(meta,nlevels)[sapply(meta,nlevels)>1])
      group_columns_no_single <- setdiff(tmp,"SampleID")
      updateSelectInput(session,"confounding_var",choices=group_columns_no_single)
    }
  })
  
  #this part needs to be in its own "observe" block
  #-> updates ref choice in section "functional topics"
  observe({
    if(!is.null(currentSet())){
      ref_choices <- unique(sample_data(vals$datasets[[currentSet()]]$phylo)[[input$formula]])
      updateSelectInput(session,"refs",choices=ref_choices)
    }
  })
  
  # check for update if undersampled columns are to be removed (rarefation curves)
  observeEvent(input$excludeSamples, {
    if(!is.null(currentSet())){
      if(!is.null(vals$undersampled) && input$excludeSamples == T){
        # remove undersampled columns from data
        sampledOTUData <- vals$datasets[[currentSet()]]$normalizedData[,!(colnames(vals$datasets[[currentSet()]]$normalizedData)%in%vals$undersampled)]
        sampledMetaData <- vals$datasets[[currentSet()]]$metaData[!(rownames(vals$datasets[[currentSet()]]$metaData)%in%vals$undersampled),]
        
        #save old (oversampled) data, in case the switch is turned OFF again
        vals$datasets[[currentSet()]]$old.normalizedData <- vals$datasets[[currentSet()]]$normalizedData
        vals$datasets[[currentSet()]]$old.metaData <- vals$datasets[[currentSet()]]$metaData
        old.phylo <- vals$datasets[[currentSet()]]$phylo
        vals$datasets[[currentSet()]]$old.phylo <- old.phylo
        
        #replace old (oversampled) data with new data
        vals$datasets[[currentSet()]]$normalizedData <- sampledOTUData
        vals$datasets[[currentSet()]]$metaData <- sampledMetaData
        
        #build new phyloseq object (with old tree & tax)
        py.otu <- otu_table(sampledOTUData,T)
        py.meta <- sample_data(sampledMetaData)
        #old.tree <- phy_tree(old.phylo)
        if(!is.null(access(old.phylo,"phy_tree"))) old.tree <- phy_tree(old.phylo) else old.tree <- NULL
        old.taxa <- tax_table(old.phylo)
        vals$datasets[[currentSet()]]$phylo <- merge_phyloseq(py.otu, py.meta, old.tree, old.taxa)
        
        #set global Set variable to TRUE, indicating, that undersampled data is already removed
        vals$datasets[[currentSet()]]$undersampled_removed <- T

      }else if (input$excludeSamples == F && vals$datasets[[currentSet()]]$undersampled_removed == T){
        #case: undersampled data was removed but shall be used again (switch turned OFF)
        #use old (oversampled) data again, which was saved 
        vals$datasets[[currentSet()]]$normalizedData <- vals$datasets[[currentSet()]]$old.normalizedData
        vals$datasets[[currentSet()]]$metaData <- vals$datasets[[currentSet()]]$old.metaData
        vals$datasets[[currentSet()]]$phylo <- vals$datasets[[currentSet()]]$old.phylo
        vals$datasets[[currentSet()]]$undersampled_removed = F
      }
    }
  })
  
  ################################################################################################################################################
  #
  # QC plots
  #
  ################################################################################################################################################
  # update targets table of the currently loaded dataset
  output$metaTable <- renderDataTable({
    if(!is.null(currentSet())){
      datatable(sample_data(vals$datasets[[currentSet()]]$phylo),filter='bottom',options=list(searching=T,pageLength=20,dom="Blfrtip",scrollX=T),editable=T,rownames=F)
    } 
    else datatable(data.frame(),options=list(dom="t"))
  })
  
  # Plot rarefaction curves
  output$rarefacCurve <- renderPlotly({
    if(!is.null(currentSet())){
      #needs integer values to work
      tab = as.matrix(vals$datasets[[currentSet()]]$rawData)
      class(tab)<-"integer"
      
      # determine data points for rarefaction curve
      rarefactionCurve = lapply(1:ncol(tab),function(i){
        n = seq(1,colSums(tab)[i],by=5000)
        if(n[length(n)]!=colSums(tab)[i]) n=c(n,colSums(tab)[i])
        drop(rarefy(t(tab[,i]),n))
      })
      slope = apply(tab,2,function(x) rareslope(x,sum(x)-100))
      vals$undersampled = colnames(tab)[slope>=quantile(slope,1-input$rareToHighlight/100)]
      
      first = order(slope,decreasing=T)[1]
      p <- plot_ly(x=attr(rarefactionCurve[[first]],"Subsample"),y=rarefactionCurve[[first]],text=paste0(colnames(tab)[first],"; slope: ",round(1e5*slope[first],3),"e-5"),hoverinfo="text",color="high",type="scatter",mode="lines",colors=c("red","black"))
      for(i in order(slope,decreasing=T)[2:input$rareToShow]){
        highslope = as.numeric(slope[i]>=quantile(slope,1-input$rareToHighlight/100))+1
        p <- p %>% add_trace(x=attr(rarefactionCurve[[i]],"Subsample"),y=rarefactionCurve[[i]],text=paste0(colnames(tab)[i],"; slope: ",round(1e5*slope[i],3),"e-5"),hoverinfo="text",color=c("low","high")[highslope],showlegend=F)
      }
      p %>% layout(title="Rarefaction Curves",xaxis=list(title="Number of Reads"),yaxis=list(title="Number of Species"))
    }
  })
  
  # show undersampled samples
  output$undersampled <- renderText({
    paste0("The following samples might be undersampled:\n",paste0(vals$undersampled,collapse=", "))
  })
  
  # plot distribution of taxa
  output$taxaDistribution <- renderPlotly({
    if(!is.null(currentSet())){
      tab = vals$datasets[[currentSet()]]$tax_binning[[which(c("Kingdom","Phylum","Class","Order","Family","Genus","Species")==input$taxLevel)]]

      
      taxa = ifelse(rowSums(tab)/ncol(tab)<(input$otherCutoff),"Other",rownames(tab))
      if(any(taxa=="Other")){
        other <- tab[which(taxa=="Other"),]
        other <- colSums(other)
        
        tab <- tab[-which(taxa=="Other"),]
        tab <- rbind(tab,Other=other)
      }
      
      tab = melt(tab)
      
      plot_ly(tab,name=~Var1,x=~value,y=~Var2,type="bar",orientation="h") %>%
        layout(xaxis=list(title="Cumulative Relative Abundance (%)"),yaxis=list(title="Samples"),
               barmode="stack",showlegend=T) #,legend=list(orientation="h")
    } else plotly_empty()
  })
  
  # visualize data structure as PCA, t-SNE or UMAP plots
  output$structurePlot <- renderPlotly({
    if(!is.null(currentSet())){
      mat <- otu_table(vals$datasets[[currentSet()]]$phylo)
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)

      samples = colnames(mat)
      mode = input$structureDim
      
      if(input$structureMethod=="PCA"){
        pca = prcomp(mat,center=T,scale=T)
        out = data.frame(pca$rotation,txt=samples)
        percentage = signif(pca$sdev^2/sum(pca$sdev^2)*100,2)
        
        xlab = paste0("PC1 (",percentage[1]," % variance)")
        ylab = paste0("PC2 (",percentage[2]," % variance)")
        zlab = paste0("PC3 (",percentage[3]," % variance)")
      }
      else if(input$structureMethod=="t-SNE"){
        tsne = Rtsne(t(mat),dim=3,perplexity=min(floor((length(samples)-1)/3),30))
        out = data.frame(tsne$Y,txt=samples)
        
        xlab = "t-SNE Dimension 1"
        ylab = "t-SNE Dimension 2"
        zlab = "t-SNE Dimension 3"
      }
      else if(input$structureMethod=="UMAP"){
        UMAP = umap(t(mat),n_components=3)
        out = data.frame(UMAP$layout,txt=samples)
        
        xlab = "UMAP Dimension 1"
        ylab = "UMAP Dimension 2"
        zlab = "UMAP Dimension 3"
      }
      colnames(out)[1:3] = c("Dim1","Dim2","Dim3")
      
      if(mode=="2D"){
        plot_ly(out,x=~Dim1,y=~Dim2,color=meta[[input$structureGroup]],text=~txt,hoverinfo='text',type='scatter',mode="markers") %>%
          layout(xaxis=list(title=xlab),yaxis=list(title=ylab))
      } else{
        plot_ly(out,x=~Dim1,y=~Dim2,z=~Dim3,color=meta[[input$structureGroup]],text=~txt,hoverinfo='text',type='scatter3d',mode="markers") %>%
          layout(scene=list(xaxis=list(title=xlab),yaxis=list(title=ylab),zaxis=list(title=zlab)))
      }
    } else{
      plotly_empty()
    }
  })
  
  # plot PCA loadings
  output$loadingsPlot <- renderPlotly({
    if(!is.null(currentSet())){
      mat <- t(otu_table(vals$datasets[[currentSet()]]$phylo))

      taxa = colnames(mat)
      
      pca = prcomp(mat,center=T,scale=T)
      loadings = data.frame(Taxa=taxa,loading=pca$rotation[,as.numeric(input$pcaLoading)])
      loadings = loadings[order(loadings$loading,decreasing=T),][c(1:10,(nrow(loadings)-9):nrow(loadings)),]
      loadings$Taxa = factor(loadings$Taxa,levels=loadings$Taxa)
      
      plot_ly(loadings,x=~Taxa,y=~loading,text=~Taxa,hoverinfo='text',type='bar',color=I(ifelse(loadings$loading>0,"blue","red"))) %>%
        layout(title="Top and Bottom Loadings",xaxis=list(title="",zeroline=F,showline=F,showticklabels=F,showgrid=F),yaxis=list(title=paste0("loadings on PC",input$pcaLoading)),showlegend=F) %>% hide_colorbar()
    } else{
      plotly_empty()
    }
  })
  
  # plot correlations between OTUs
  output$OTUcorrPlot <- renderPlotly({
    if(!is.null(currentSet())){
      mat <- vals$datasets[[currentSet()]]$rawData
      
      corCluster = corclust(t(mat),method="ward.D2")
      cor_mat = cor(t(mat))[corCluster$cluster.numerics$order,corCluster$cluster.numerics$order]
      
      clusters = cutree(corCluster$cluster.numerics,k=input$corrCluster)[corCluster$cluster.numerics$order]
      cuts = c(0,cumsum(table(clusters)[unique(clusters)]))-0.5
      shapes=list()
      for(i in 1:input$corrCluster){
        shapes[[i]] = list(type="rect",
                           line=list(color="red",width=3),
                           x0=cuts[i],x1=cuts[i+1],
                           y0=cuts[i],y1=cuts[i+1])
      }
      
      plot_ly(x=colnames(cor_mat),y=rownames(cor_mat),z=cor_mat,type="heatmap") %>%
        layout(shapes=shapes) %>%
        colorbar(title="Pearson Correlation\n")
    } else{
      plotly_empty()
    }
  })
  
  
  # plot alpha diversity
  output$alphaPlot <- renderPlotly({
    if(!is.null(currentSet())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      
      alpha = alphaDiv(otu,input$alphaMethod)
      
      if(input$alphaGroup=="-") plot_ly(y=alpha,type='violin',box=list(visible=T),meanline=list(visible=T),x0=input$alphaMethod) %>% layout(yaxis=list(title="alpha Diversity",zeroline=F))
      else plot_ly(x=meta[[input$alphaGroup]],y=alpha,color=meta[[input$alphaGroup]],type='violin',box=list(visible=T),meanline=list(visible=T),x0=input$alphaMethod) %>% layout(yaxis=list(title="alpha Diversity",zeroline=F))
    }
    else plotly_empty()
  })
  
  # show alpha diversity table
  output$alphaTable <- renderDataTable({
    if(!is.null(currentSet())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)

      alphaTab = data.frame(colnames(otu))
      for(i in c("Shannon Entropy","effective Shannon Entropy","Simpson Index","effective Simpson Index","Richness")){
        alphaTab = cbind(alphaTab,round(alphaDiv(otu,i),2))
      }
      colnames(alphaTab) = c("SampleID","Shannon Entropy","effective Shannon Entropy","Simpson Index","effective Simpson Index","Richness")
      datatable(alphaTab,filter='bottom',options=list(searching=T,pageLength=20,dom="Blfrtip"),editable=T,rownames=F)
    }
    else datatable(data.frame(),options=list(dom="t"))
  })
  
  
  
  # clustering tree of samples based on beta-diversity
  output$betaTree <- renderPlot({
    if(!is.null(currentSet())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
      #need to flook for tree object like this, did not find out better/cleaner way yet
      if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      group = input$betaGroup
      
      method = ifelse(input$betaMethod=="Bray-Curtis Dissimilarity","brayCurtis","uniFrac")
      distMat = betaDiversity(otu=otu,meta=meta,tree=tree,group=group,method=method)
      
      all_fit = hclust(distMat,method="ward")
      tree = as.phylo(all_fit)
      all_groups = as.factor(meta[[group]])
      col = rainbow(length(levels(all_groups)))[all_groups]
      
      plot(tree,type="phylogram",use.edge.length=T,tip.color=col,label.offset=0.01)
      axisPhylo()
      tiplabels(pch=16,col=col)
    }
  })
  
#  MDS plot based on beta-diversity
  output$betaMDS <- renderPlot({
    if(!is.null(currentSet())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      #need to flook for tree object like this, did not find out better/cleaner way yet
      if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      group = input$betaGroup

      method = ifelse(input$betaMethod=="Bray-Curtis Dissimilarity","brayCurtis","uniFrac")
      dist = betaDiversity(otu=otu,meta=meta,tree=tree,group=group,method=method)

      #all_groups = as.factor(meta[[group]])
      all_groups = as.factor(meta[,group][[1]])
      adonis = adonis(dist ~ all_groups)
      all_groups = factor(all_groups,levels(all_groups)[unique(all_groups)])

      # Calculate and display the MDS plot (Multidimensional Scaling plot)
      col = rainbow(length(levels(all_groups)))
      s.class(
        cmdscale(dist,k=2),col=col,cpoint=2,fac=all_groups,
        sub=paste("MDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
      )
    }

  })
  

  
  # NMDS plot based on beta-diversity
  output$betaNMDS <- renderPlot({
    if(!is.null(currentSet())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      meta <- sample_data(vals$datasets[[currentSet()]]$phylo)
      #need to flook for tree object like this, did not find out better/cleaner way yet
      if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      group = input$betaGroup

      method = ifelse(input$betaMethod=="Bray-Curtis Dissimilarity","brayCurtis","uniFrac")
      dist = betaDiversity(otu=otu,meta=meta,tree=tree,group=group,method=method)

      all_groups = as.factor(meta[[group]])
      adonis = adonis(dist ~ all_groups)
      all_groups = factor(all_groups,levels(all_groups)[unique(all_groups)])

      # Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
      meta_mds = metaMDS(dist,k=2)
      col = rainbow(length(levels(all_groups)))
      s.class(
        meta_mds$points,col=col,cpoint=2,fac=all_groups,
        sub=paste("metaNMDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
      )
    }

  })
  
  #phylogenetic tree using phyloseq
  output$phyloTree <- renderPlot({
    if(!is.null(currentSet())){
      phylo <- vals$datasets[[currentSet()]]$phylo
      #need to flook for tree object like this, did not find out better/cleaner way yet
      if(!is.null(access(vals$datasets[[currentSet()]]$phylo,"phy_tree"))) tree <- phy_tree(vals$datasets[[currentSet()]]$phylo) else tree <- NULL
      
      #only can display a phylogenetic tree if tree object is present
      if(!is.null(tree)){
        #prune taxa to number the user sets with slider (changes size of tree)
        pruned_phylo <- prune_taxa(taxa_names(phylo)[1:input$phylo_prune], phylo)
        
        if(input$phylo_color == "-") phylo_color = NULL else phylo_color=input$phylo_color
        if(input$phylo_shape == "-") phylo_shape = NULL else phylo_shape=input$phylo_shape
        if(input$phylo_size == "-") phylo_size = NULL else phylo_size=input$phylo_size
        if(input$phylo_label.tips == "-") phylo_label.tips = NULL else phylo_label.tips=input$phylo_label.tips
        
        if(input$phylo_radial == T){
          plot_tree(pruned_phylo,method = input$phylo_method,color=phylo_color,shape = phylo_shape,size = phylo_size,label.tips = phylo_label.tips,ladderize = "left",plot.margin = input$phylo_margin)+coord_polar(theta = "y")
        }else{
          plot_tree(pruned_phylo,method = input$phylo_method,color=phylo_color,shape = phylo_shape,size = phylo_size,label.tips = phylo_label.tips,ladderize = input$phylo_ladderize,plot.margin = input$phylo_margin)
        }
      }
    }
  })
  
  #confounding analysis
  observeEvent(input$confounding_start,{
    if(!is.null(currentSet())){
      #only if pre-build distance matrix exists, this can be calcualted (depending on tree input)
      if (!is.null(vals$datasets[[currentSet()]]$unifrac_dist)){
        meta <- vals$datasets[[currentSet()]]$metaData
        #remove first column --> SampleID column
        meta[,1]<-NULL
        
        #calulate confounding matrix
        withProgress(message="Calculating confounding factors...",value=0,{
          vals$datasets[[currentSet()]]$confounder_table <- calculateConfounderTable(var_to_test=input$confounding_var,
                                                                                     variables = meta,
                                                                                     distance=vals$datasets[[currentSet()]]$unifrac_dist,
                                                                                     permutations = input$confounding_perm,
                                                                                     useSeed=input$confounding_seed,
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
  
  output$confounding_var_text <- renderUI({
    if(!is.null(currentSet())){
      if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
        HTML(paste0("<b> Chosen variable: </b> ",vals$datasets[[currentSet()]]$confounder_table$var))
      }
    }
  })
  
  
  #####################################
  #    Network analysis               #
  #####################################
  
  
  # show histogram of all OTU values -> user can pick cutoff for binarization here
  output$cutoffHist <- renderPlotly({
    if(!is.null(currentSet())){
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      dat <- log(as.data.frame(otu))
      
      plot_ly(x=unlist(dat),type="histogram") %>%
        layout(xaxis=list(title="log(normalized OTU-values)"), yaxis = list(title="Frequency"),
               shapes=list(list(type="line",y0=0,y1=1,yref="paper",x0=log(input$binCutoff),
                                x1=log(input$binCutoff),line=list(color="black",width=2))))
    } else{
      plotly_empty()
    }
  })
  
  # cutoff heatmap
  output$boolHeat <- renderPlotly({
    if(!is.null(currentSet())){
      cutoff <- input$binCutoff
      otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
      
      m <- as.matrix(otu)
      m <- apply(m,c(1,2),function(x){ifelse(x<cutoff,0,1)})
      
      melted_m<-melt(m)
      melted_m$Var1 <- as.factor(melted_m$Var1)
      melted_m$Var2 <- as.factor(melted_m$Var2)
      
      plot_ly(z=melted_m$value,x=melted_m$Var1,y=melted_m$Var2,type="heatmap",showscale=F)
      
    } else{
      plotly_empty()
    }
  })
  
  # check if button for new calculation of counts is clicked -> reload network with the new counts
  observeEvent(input$startCalc,{
    
    withProgress(message = 'Calculating Counts..', value = 0, {
      vals$datasets[[currentSet()]]$counts = generate_counts(OTU_table=vals$datasets[[currentSet()]]$normalizedData,
                                                             meta = vals$datasets[[currentSet()]]$metaData,
                                                             group_column = input$groupCol,
                                                             cutoff = input$binCutoff,
                                                             fc = ifelse(input$useFC=="log2(fold-change)",T,F),
                                                             progress = T)
    })
  })
  
  # network plot
  output$basicNetwork <- renderForceNetwork({
    if(!is.null(currentSet())){
      if(!is.null(vals$datasets[[currentSet()]]$counts)){
        
        otu <- otu_table(vals$datasets[[currentSet()]]$phylo)
        counts = vals$datasets[[currentSet()]]$counts
        tax = tax_table(vals$datasets[[currentSet()]]$phylo)
        
        
        Links = counts[order(abs(counts$value),decreasing=T)[1:input$networkCutoff],]
        colnames(Links) = c("source","target","value")
        Links$source = as.character(Links$source); Links$target = as.character(Links$target); Links$valueToPlot = abs(Links$value)
        Links$valueToPlot = (Links$valueToPlot-min(Links$valueToPlot))/(max(Links$valueToPlot)-min(Links$valueToPlot))*2
        
        Nodes = data.frame(name=unique(c(Links$source,Links$target)),group="")
        Nodes$size = rowSums(otu[Nodes$name,])/1000

        if(input$netLevel!="-") Nodes$group = substring(tax[Nodes$name,input$netLevel],4)
        Nodes$group[Nodes$group==""] = "unknown"
        
        Links$source = match(Links$source,Nodes$name)-1
        Links$target = match(Links$target,Nodes$name)-1
        
        forceNetwork(Links,Nodes,Source="source",Target="target",Value="valueToPlot",NodeID="name",
                     Nodesize="size",Group="group",linkColour=c("red","green")[(Links$value>0)+1],zoom=T,legend=T,
                     bounded=T,fontSize=12,fontFamily='sans-serif',charge=-25,linkDistance=100,opacity = 1)
        
      }
    }
  })
  
  
  #####################################
  #    themetagenomcis apporach       #
  #####################################
  
  #here all objects and values needed for the plots of themetagenomics are created and stored in vals$datasets[[currentSet()]]$vis_out
  #does not work with phyloseq-object (error if otu-table is of class phyloseq::otu_table)
  
  observeEvent(input$themeta,{
    withProgress(message='Calculating Topics..',value=0,{
      if(!is.null(currentSet())){
        #take otu table and meta file from user input
        otu <- vals$datasets[[currentSet()]]$normalizedData
        meta <- vals$datasets[[currentSet()]]$metaData
        tax = vals$datasets[[currentSet()]]$taxonomy
        
        incProgress(1/7,message="preparing OTU data..")
        formula <- as.formula(paste0("~ ",input$formula))
        formula_char <- input$formula
        refs <- input$refs
        
        #create themetadata object
        obj <- prepare_data(otu_table = otu,
                            rows_are_taxa = T,
                            tax_table = tax,
                            metadata = meta,
                            formula=formula,
                            refs = refs,
                            cn_normalize = FALSE)
        
        incProgress(1/7,message = "finding topics..")
        K=input$K
        sigma_prior = input$sigma_prior
        #use themetadata object to find K topics
        topics_obj <- find_topics(themetadata_object=obj,
                                  K=K,
                                  sigma_prior = sigma_prior)
        
        incProgress(1/7,message="predicting topic functions..")
        #functions_obj <- predict.topics(topics_obj,reference_path = "themetagenomics/")
        
        incProgress(1/7,message = "generate gene.table")
        #topic_function_table <- gene.table(functions_obj)
        
        incProgress(1/7,message = "finding topic effects..")
        #measure relationship of covarite with samples over topics distribution from the STM
        topic_effects_obj <- est(topics_obj)
        
        
        #function_effects <- themetagenomics::est(functions_obj,topics_subset=3)
        
        
        class(topics_obj) <- "topics"
        incProgress(1/7,message = "preparing visualization..")
        vals$datasets[[currentSet()]]$vis_out <- prepare_vis(topic_effects_obj)
        vals$datasets[[currentSet()]]$vis_out$K <- K
        vals$datasets[[currentSet()]]$vis_out$sigma <- sigma_prior
        vals$datasets[[currentSet()]]$vis_out$formula <- formula_char
        vals$datasets[[currentSet()]]$vis_out$refs <- refs
        vals$datasets[[currentSet()]]$topic_effects <- topic_effects_obj
        #vals$datasets[[currentSet()]]$gene_table <- topic_function_table
        incProgress(1/7)
      }
    })
  })
  
  #make gene.table download-able
  # output$downloadGeneTable <- downloadHandler(
  #   filename = "gene_table.csv",
  #   content = function(file){
  #     write.csv(vals$datasets[[currentSet()]]$gene_table,file,row.names = T)
  #   }
  # )
  
  REL <- reactive({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      if (show_topic$k != 0){
        
        current_k <- paste0('T',show_topic$k)
        
        l <- input$lambda
        
        tinfo_k <- vis_out$tinfo[vis_out$tinfo$Category == current_k,]  #subset(tinfo,Category == current_k)
        rel_k <- l*tinfo_k$logprob + (1-l)*tinfo_k$loglift
        new_order <- tinfo_k[order(rel_k,decreasing=TRUE)[1:vis_out$taxa_bar_n],]
        new_order$Term <- as.character.factor(new_order$Term)
        new_order$Taxon <- vis_out$taxa[new_order$Term,input$taxon]
        new_order$Term <- factor(new_order$Term,levels=rev(new_order$Term),ordered=TRUE)
        
      }else{
        
        new_order <- vis_out$default
        new_order$Taxon <- vis_out$taxa[as.character.factor(new_order$Term),input$taxon]
        
      }
      
      new_order
    }
    
  })
  
  EST <- reactive({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    topic_effects <- vals$datasets[[currentSet()]]$topic_effects$topic_effects
    
    if(!is.null(vis_out) & input$choose != "Please start calculation above first!"){
      suppressWarnings({
        
        covariate <- input$choose
        
        est_mat <- topic_effects[[covariate]]$est
        
        df0 <- data.frame(topic=paste0('T',1:vis_out$K),
                          est=est_mat[,1],
                          lower=est_mat[,2],
                          upper=est_mat[,3],
                          sig=ifelse(1:vis_out$K %in% topic_effects[[covariate]]$sig,'1','0'))
        df0$sig <- factor(as.character(sign(df0$est) * as.numeric(as.character(df0$sig))),levels=c('0','1','-1'),ordered=TRUE)
        df <- df0[order(topic_effects[[covariate]][['rank']]),]
        df$topic <- factor(df$topic,levels=df$topic,ordered=TRUE)
        
        p_est <- ggplot(df,aes_(~topic,y=~est,ymin=~lower,ymax=~upper,color=~sig)) +
          geom_hline(yintercept=0,linetype=3)
        p_est <- p_est +
          geom_pointrange(size=2) +
          theme_minimal() +
          labs(x='',y='Estimate') +
          scale_color_manual(values=c('gray','indianred3','dodgerblue3'),drop=FALSE) +
          scale_fill_brewer(type='qual',palette=6,direction=-1) +
          theme(legend.position='none',
                axis.text.x=element_text(angle=-90,hjust=0,vjust=.5))
        
        list(p_est=p_est,k_levels=levels(df$topic),covariate=covariate,df0=df0)
        
      })
    }
  })
  
  output$est <- renderPlotly({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out)){
        suppressWarnings(ggplotly(EST()$p_est,source='est_hover',tooltip=c('topic','est','lower','upper'))) #Error in UseMethod: no applicable method for 'plotly_build' applied to an object of class "shiny.tag"
      }else{
        plotly_empty()
      }
    }

    
  })
  
  output$ord <- renderPlotly({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out) & !is.null(EST())){
        beta <- t(vis_out$beta)
        
        if (input$dist == 'hellinger'){
          
          d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'norm'),method='euclidean'),3,eig=TRUE)
          
        }else if (input$dist == 'chi2'){
          
          d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'chi.square'),method='euclidean'),3,eig=TRUE)
          
        }else if (input$dist == 'jsd'){
          
          d <- cmdscale(proxy::dist(beta,jsd),3,eig=TRUE)   #woher kommt jsd?
          
        } else if (input$dist == 'tsne'){
          p <- 30
          d <- try(Rtsne::Rtsne(beta,3,theta=.5,perplexity=p),silent=TRUE)
          while(class(d) == 'try-error'){
            p <- p-1
            d <- try(Rtsne::Rtsne(beta,3,theta=.5,perplexity=p),silent=TRUE)
          }
          if (p < 30) cat(sprintf('Performed t-SNE with perplexity = %s.\n',p))
          d$points <- d$Y
          d$eig <- NULL
          
        }else{
          
          d <- cmdscale(vegan::vegdist(beta,method=input$dist),3,eig=TRUE)
          
        }
        
        eig <- d$eig[1:3]/sum(d$eig)
        colnames(d$points) <- c('Axis1','Axis2','Axis3')
        df <- data.frame(d$points,EST()$df0)
        df$marg <- vis_out$topic_marg
        
        df$colors <- vis_out$colors[as.character(df$sig)]
        
        if (input$dim == '2d'){
          
          p1 <- plot_ly(df,source='ord_click')
          p1 <- add_trace(p1,
                          x=~Axis1,y=~Axis2,size=~marg,
                          type='scatter',mode='markers',sizes=c(5,125),
                          color=I(df$colors),opacity=.5,
                          marker=list(symbol='circle',sizemode='diameter',line=list(width=3,color='#FFFFFF')),
                          text=~paste('<br>Topic:',topic),hoverinfo='text')
          p1 <- layout(p1,
                       showlegend=FALSE,
                       xaxis=list(title=sprintf('Axis 1 [%.02f%%]',eig[1]*100),
                                  showgrid=FALSE),
                       yaxis=list(title=sprintf('Axis 2 [%.02f%%]',eig[2]*100),
                                  showgrid=FALSE),
                       paper_bgcolor='rgb(243, 243, 243)',
                       plot_bgcolor='rgb(243, 243, 243)')
          p1 <- add_annotations(p1,x=df$Axis1,y=df$Axis2,text=df$topic,showarrow=FALSE,
                                font=list(size=10))
          
          h <- event_data('plotly_hover',source='est_hover')
          
          if (!is.null(h)){
            k <- EST()$k_levels[h[['x']]]
            
            if (length(k) > 0){
              df_update <- df[df$topic == k,]
              
              
              if (df_update$sig == '1') df_update$sig <- '2' else if(df_update$sig== '-1') df_update$sig <- '-2' else df_update$sig<- '00'
              df_update$colors <- vis_out$colors[df_update$sig]
              
              p1 <- add_markers(p1,
                                x=df_update$Axis1,y=df_update$Axis2,opacity=.8,color=I(df_update$color),
                                marker=list(size=150,symbol='circle',sizemode='diameter',line=list(width=3,color='#000000')))
            }
            
            p1
            
          }
          
        }
        
        if (input$dim == '3d'){
          
          p1 <- plot_ly(df,source='ord_click',
                        x=~Axis1,y=~Axis2,z=~Axis3,size=~marg,
                        type='scatter3d',mode='markers',sizes=c(5,125),
                        color=I(df$colors),opacity=.5,
                        marker=list(symbol='circle',sizemode='diameter'),
                        text=~paste('<br>Topic:',topic),hoverinfo='text')
          
          p1 <- layout(p1,
                       showlegend=FALSE,
                       scene=list(
                         xaxis=list(title=sprintf('Axis 1 [%.02f%%]',eig[1]*100),
                                    showgrid=FALSE),
                         yaxis=list(title=sprintf('Axis 2 [%.02f%%]',eig[2]*100),
                                    showgrid=FALSE),
                         zaxis=list(title=sprintf('Axis 3 [%.02f%%]',eig[3]*100),
                                    showgrid=FALSE)),
                       paper_bgcolor='rgb(243, 243, 243)',
                       plot_bgcolor='rgb(243, 243, 243)')
          
        }
        
        p1
      }
    }
  })
  
  show_topic <- reactiveValues(k=0)
  
  #maybe not working...
  observeEvent(event_data('plotly_click',source='ord_click'),{
    
    s <- event_data('plotly_click',source='ord_click')
    
    if (is.null(s)){
      
      show_topic$k <- 0
      updateNumericInput(session,'k_in',value=0)
      
    }else{
      
      t_idx <- s$pointNumber + 1
      updateNumericInput(session,'k_in',value=t_idx)
      show_topic$k <- t_idx
      
    }
    
  })
  
  observeEvent(input$k_in,{
    show_topic$k <- input$k_in
  })
  
  observeEvent(input$reset,{
    show_topic$k <- 0
    updateNumericInput(session,'k_in',value=0)
  })
  
  output$bar <- renderPlot({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out)){
        if (show_topic$k != 0){
          p_bar <- ggplot(data=REL()) +
            geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=.6) +
            geom_bar(aes_(~Term,~Freq),stat='identity',fill='darkred',color='white')
        } else{
          p_bar <- ggplot(data=REL()) +
            geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=1)
        }
        
        p_bar +
          coord_flip() +
          labs(x='',y='Frequency',fill='') +
          theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
                legend.position='bottom') +
          viridis::scale_fill_viridis(discrete=TRUE,drop=FALSE) +
          guides(fill=guide_legend(nrow=2))
      }
    }
  })
  
  output$corr <- renderForceNetwork({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      topic_effects <- vals$datasets[[currentSet()]]$topic_effects$topic_effects
      if(!is.null(vis_out)){
        suppressWarnings(effects_sig <- topic_effects[[EST()$covariate]][['sig']])  #Warning: Error in [[: attempt to select less than one element in get1index
        K <- nrow(vis_out$corr$posadj)
        
        suppressWarnings({suppressMessages({
          g <- igraph::graph.adjacency(vis_out$corr$posadj,mode='undirected',
                                       weighted=TRUE,diag=FALSE)
          
          wc <- igraph::cluster_walktrap(g)
          members <- igraph::membership(wc)
          
          g_d3 <- networkD3::igraph_to_networkD3(g,group=members)
          
          g_d3$links$edge_width <- 10*(.1+sapply(seq_len(nrow(g_d3$links)),function(r) vis_out$corr$poscor[g_d3$links$source[r]+1,g_d3$links$target[r]+1]))
          g_d3$nodes$color <- 25*ifelse(1:K %in% effects_sig,1,0)*sign(topic_effects[[EST()$covariate]]$est[,1])
          g_d3$nodes$node_size <- 10*(.5+norm10(c(0,abs(topic_effects[[EST()$covariate]]$est[,1])))[-1])
          g_d3$nodes$name <- paste0('T',g_d3$nodes$name)
          
          networkD3::forceNetwork(Links=g_d3$links,Nodes=g_d3$nodes,
                                  Source='source',Target='target',
                                  charge=-25,
                                  opacity=1,
                                  fontSize=12,
                                  zoom=TRUE,
                                  bounded=TRUE,
                                  NodeID='name',
                                  fontFamily='sans-serif',
                                  opacityNoHover=.7,
                                  Group='color',
                                  Value='edge_width',
                                  Nodesize='node_size',
                                  linkColour='#000000',
                                  linkWidth=networkD3::JS('function(d) {return d.value;}'),
                                  radiusCalculation=networkD3::JS('d.nodesize'),
                                  colourScale=networkD3::JS("color=d3.scaleLinear()\n.domain([-1,0,1])\n.range(['blue','gray','red']);"))
          
        })})
      }
    }
  })
  
  
  #####################################
  #    SPIEC-EASI                     #
  #####################################
  
  #WATCH OUT!! this tool needs the otu-table to be in format: OTU in column; sample in row
  #or simply the phyloseq class
  
  
  ## Meinshausen-Buhlmann's ##
  
  observeEvent(input$se_mb_start,{
    if(!is.null(currentSet())){
      withProgress(message = 'Calculating mb..', value = 0, {
        py <- vals$datasets[[currentSet()]]$phylo
        taxa <- tax_table(py)
        incProgress(1/2,message = "starting calculation..")
        se_mb <- spiec.easi(py, method = "mb", lambda.min.ratio = input$se_mb_lambda.min.ratio, nlambda = input$se_mb_lambda, pulsar.params = list(rep.num=input$se_mb_repnumber, ncores =ncores))
        incProgress(1/2,message = "building graph objects..")
        
        #pre-build graph object for phyloseq graph
        se_mb$ig <- adj2igraph(getRefit(se_mb), vertex.attr=list(name=taxa_names(py)))
        #pre-build grapg for interactive networkD3 graph
        nd3 <-igraph_to_networkD3(se_mb$ig, taxa)
        
        output$spiec_easi_mb_network <- renderPlot({
          plot_network(se_mb$ig,py,type = "taxa",color=as.character(input$mb_select_taxa))
        })
        output$spiec_easi_mb_network_interactive <- renderForceNetwork(forceNetwork(Links=nd3$links,Nodes=nd3$nodes,NodeID = "name",Group = as.character(input$mb_select_taxa),
                                                                                    zoom=T,legend=T,fontSize = 5,charge = -2,opacity = .9, height = 200,width = 100))
      })
    }
  })
  
  
  ## Glasso ## 
  
  observeEvent(input$se_glasso_start,{
    if(!is.null(currentSet())){
      withProgress(message = 'Calculating glasso..', value = 0, {
        py <- vals$datasets[[currentSet()]]$phylo
        taxa <- tax_table(py)
        incProgress(1/2,message = "starting calculation..")
        se_glasso <- spiec.easi(py, method = "glasso", lambda.min.ratio = input$glasso_mb_lambda.min.ratio, nlambda = input$glasso_mb_lambda, pulsar.params = list(rep.num=input$se_glasso_repnumber, ncores = ncores))
        incProgress(1/2,message = "building graph objects..")
        
        #pre-build graph object for phyloseq graph
        se_glasso$ig <- adj2igraph(getRefit(se_glasso), vertex.attr=list(name=taxa_names(py)))
        #pre-build grapg for interactive networkD3 graph
        nd3 <-igraph_to_networkD3(se_glasso$ig, taxa)
        
        output$spiec_easi_glasso_network <- renderPlot({
          plot_network(se_glasso$ig,py,type = "taxa",color=as.character(input$mb_select_taxa))
        })
        output$spiec_easi_glasso_network_interactive <- renderForceNetwork(forceNetwork(Links=nd3$links,Nodes=nd3$nodes,NodeID = "name",Group = as.character(input$glasso_select_taxa),
                                                                                    zoom=T,legend=T,fontSize = 5,charge = -2,opacity = .9, height = 200,width = 100))
        
      })
    }
  })

  ## SparCC ##
  
  # observeEvent(input$se_sparcc_start,{
  #   if(!is.null(currentSet())){
  #     withProgress(message = 'Calculating SparCC..', value = 0, {
  #       otu <- vals$datasets[[currentSet()]]$normalizedDat
  #       otu.t <- t(otu)
  #       
  #       incProgress(1/2,message = "starting calculation..")
  #       se_sparcc <- sparcc(otu.t)
  #       
  #       #set threshold on matrix
  #       sparcc.graph <- abs(se_sparcc$Cor) >= input$se_sparcc_threshold
  #       diag(sparcc.graph) <- 0
  #       sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
  #       incProgress(1/2,message = "building graph..")
  #       vals$datasets[[currentSet()]]$se_sparcc$graph <- adj2igraph(sparcc.graph)
  #     })
  #   }
  # })
  # 
  # output$spiec_easi_sparcc_network <- renderPlot({
  #   if(!is.null(currentSet())){
  #     sparcc_ig <- vals$datasets[[currentSet()]]$se_sparcc$graph
  #     otu <- vals$datasets[[currentSet()]]$normalizedDat
  #     otu.t <- t(otu)
  #     if(!is.null(sparcc_ig) && !is.null(otu.t)){
  #       vsize  <- rowMeans(clr(otu.t, 1))+6
  #       am.coord <- layout.fruchterman.reingold(sparcc_ig)
  #       
  #       plot(sparcc_ig, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="SparCC")
  #     }
  #   }
  # }) 
  

  #####################################
  #    Text fields                    #
  #####################################
  output$welcome <- renderUI({
    HTML(paste0("<h3> Welcome to <i>namco</i>, the free Microbiome Explorer!</h3>",
                "<img src=\"Logo.png\" alt=\"Logo\" width=400 height=400>"))
  })
  
  output$authors <- renderUI({
    HTML(paste0("<b>Authors of this tool:</b>",
                "Alexander Dietrich (Contact for Issues: ga89noj@mytum.de)",
                "Benjamin Ölke, ",
                "Maximilian Zwiebel <br> ",
                "Supervisor: Michael Lauber, Dr. Markus List, Prof. Dr. Jan Baumbach <br>",
                "Experimental Chair of Bioinformatics, TU München <br>"))
  })
  
  output$welcome_ref <- renderUI({
    HTML(paste0("This tool was built using source-code from <br> ",
                "<b>RHEA</b>: Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836 <br>",
                "<b>themetagenomics</b>: Stephen Woloszynek, Joshua Chang Mell, Gideon Simpson, and Gail Rosen. Exploring thematic structure in 16S rRNA marker gene surveys. 2017. bioRxiv 146126; doi: https://doi.org/10.1101/146126"))
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
        HTML(paste0('Below shows topic-to-topic correlations from the samples over topics distribution. The edges represent positive',
                    ' correlation between two topics, with the size of the edge reflecting to the magnitude of the correlation.',
                    ' The size of the nodes are consistent with the ordination figure, reflecting the marginal topic frequencies.'))
      }
    }
  })
  
  output$text3 <- renderUI({
    HTML(paste0("Integrates the samples over topics p(s|k) and topics over taxa p(k|t) distributions from the STM, the topic correlations from the p(s|k) component, the covariate effects from the p(s|k) component, and their relationship with the raw taxonomic abundances. The covariate effects for each topic are shown as a scatterplot of posterior weights with error bars corresponding the global approximation of uncertainty. If the covariate chosen is binary, then the posterior regression weights with uncertainty intervals are shown. This is analogous to the mean difference between factor levels in the posterior predictive distribution. For continuous covariates, the points again represent the mean regression weights (i.e., the posterior slope estimate of the covariate). If, however, a spline or polynomial expansion was used, then the figure shows the posterior estimates of the standard deviation of the predicted topic probabilities from the posterior predictive distribution. Colors indicate whether a given point was positive (red) or negative (blue) and did not enclose 0 at a user defined uncertainty interval.",
                "The ordination figure maintains the color coding just decribed. The ordination is performed on p(k|t) via either PCoA (using either Jensen-Shannon, Euclidean, Hellinger, Bray-Curtis, Jaccard, or Chi-squared distance) or t-SNE. The latter iterates through decreasing perplexity values (starting at 30) until the algorithm succeeds. The top 2 or 3 axes can be shown. The radius of the topic points corresponds to the topic frequencies marginalized over taxa.",
                "The bar plot behaves in accordance with LDAvis. When no topics are chosen, the overall taxa frequencies are shown. These frequencies do not equal the abundances found in the initial abundance table. Instead, they show p(k|t) multiplied by the marginal topic distribution (in counts). To determine the initial order in which taxa are shown, these two distributions are compared via Kullback-Liebler divergence and then weighted by the overall taxa frequency. The coloration of the bars indiciates the taxonomic group the inidividual taxa belong to. The groups shown are determined based on the abundance of that group in the raw abundance table. When a topic is selected, the relative frequency of a given taxa in that topic is shown in red.",
                "λ controls relevance of taxa within a topic, which in turn is used to adjust the order in which the taxa are shown when a topic is selected. Relevence is essentially a weighted sum between the probability of taxa in a given topic and the probability of taxa in a given topic relative to the overall frequency of that taxa. Adjusting λ influences the relative weighting such that",
                "r = λ x log p(t|k) + λ x log p(t|k)/p(x)",
                "The correlation graph shows the topic correlations from p(s|k) ~ MVN(mu,sigma). Again, the coloration described above is conserved. The size of the nodes reflects the magnitude of the covariate posterior regression weight, whereas the width of the edges represents the value of the positive correlation between the connected nodes. By default, the graph estimates are determined using the the huge package, which first performs a nonparanormal transformation of p(s|k), followed by a Meinhuasen and Buhlman procedure. Alternatively, by choosing the simple method, the correlations are simply a thresholded MAP estimate of p(s|k)."))
  })
  
  output$cutoff_title <- renderUI({
    HTML(paste0("<h4><b>Cutoff Assignment:<sup>1</sup></b></h4>"))
  })
  
  output$cutoff_text <- renderUI({
    HTML(paste0("This shows the logarithmic distribution in the normalized OTU table (black line is currently selected cutoff)"))
  })
  
  output$heatmap_text <- renderUI({
    HTML(paste0("Heatmap of cutoff-effect: dark fields are being set to 0 in co-occurrence calculation. Y are samples, X are OTUs"))
  })
  
  output$basic_additional <- renderUI({
    HTML(paste0("<sup>1</sup>: All values in the normalized OTU-table, which fall below the chosen cutoff are being set to 0. Those OTUs are then counted as <i> not present </i>. "))
  })
  
  output$basic_calc_title <- renderUI({
    HTML(paste0("<h4><b> igure Count Calculation:<sup>2</sup></b></h4>"))
  })
  
  output$basic_calc_additional <- renderUI({
    HTML(paste0("<sup>2</sup>: Two ways of calculating the counts are possible: <br>",
                "Both methods start by counting the co-occurrences of OTUs in all possible sample pairings. Those co-occurrences are then added up to generate a count table for each OTU-pair.",
                "By chosing a group from the meta file, this process is executed seperatly for all samples in the group corresponing to a unique covariate. The two tables are then compared: <br>",
                "<b> difference</b>: For each OTU pair x and y, calculate: counts(x) - counts(y), where x is the first occuring covariate. <br>",
                "<b> log <sub>2</sub> fold-change</b>: For each OTU pair x and y calculate: log<sub>2</sub>(counts(x)+0.001 / counts(y)+0.001), where x is the first occuring covariate."))
  })
  
  output$input_variables <- renderUI({
    if(!is.null(currentSet())){
      vis_out <- vals$datasets[[currentSet()]]$vis_out
      if(!is.null(vis_out)){
        K <- vis_out$K
        sigma <- vis_out$sigma_prior
        formula<-vis_out$formula
        refs<-paste(vis_out$refs,collapse=", ")
        HTML(paste0("<b> Number of chosen topics (K): </b>",K,"<br>",
                    "<b> Value of sigma_prior: </b>",sigma,"<br>",
                    "<b> Group from META file: </b>",formula, "<br>",
                    "<b> Reference Level in this group: </b>",refs))
      }else{
        HTML(paste0("<b> Number of chosen topics (K): </b>","<br>",
                    "<b> Value of sigma_prior: </b>","<br>",
                    "<b> Group from META file: </b>","<br>",
                    "<b> Reference Level in this group: </b>"))
      }
    }
  })
  
  output$advanced_text <- renderUI({
    HTML(paste0("<h5>Explore clustering by functional topics in your dataset! Choose covariate of interest to measure its relationship with the samples over topics distribution from the STM. </h5> <i>(for detailed explanation of tool scroll to bottom of page)</i>"))
  })
  
  
  output$topic_text <- renderUI({
    HTML(paste0("Pick the number of functional clusters you want to split your OTUs into"))
  })
  
  output$sigma_text <- renderUI({
    HTML(paste0("This sets the strength of regularization towards a diagonalized covariance matrix. Setting the value above 0 can be useful if topics are becoming too highly correlated. Default is 0"))
  })
  
  output$spiec_easi_additional <- renderUI({
    HTML(paste0("<b>1:</b> absolute value of correlations below this threshold are considered zero by the inner SparCC loop."))
  })
  
  output$info_inputdata <- renderUI({
    HTML(paste0("<p>Namco has 2 obligatory input files &amp; and one optional file:</p><p><span style='text-decoration: underline;'>1) OTU-Table:</span> tab-seperated table, where rows represent OTUs amd columns represent samples. Additionally the last column in the file needs to include the <strong>taxonomic information</strong> for the corresponding OTU of that row. This information must be in order: <em>Kingdom, Phylum, Class, Order, Family, Genus and Species</em>, seperated by semicolon. If information for any level is missing the space has to be kept empty, but still, the semicolon has to be present. For an OTU with only taxonomic information of the kingdom the entry would look like this: <code>Bacteria;;;;;;</code></p><p>Namco expects un-normalized input data and allows for normalization during file upload; this can be switched off in the upload window if the user has already normalized data.</p><p><span style='text-decoration: underline;'>2) Meta-file:</span> tab-seperated table with meta information for each sample, mapping at least one categorical variable to those. The first column has to be <b>identical</b> with the column names of the OTU file and has to named <b>SampleID</b></p><p><span style='text-decoration: underline;'>3) Phylogenetic tree:</span> To access the full functionalities provided by namco in addition to the OTU table and themapping file, we require a (rooted) phylogenetic tree of representative sequences from each OTU <strong>in Newick format</strong>. Providing an OTU tree is optional, however required to calculate certain measures of beta diversity for instance.</p>"))
  })
  
  output$info_testdata <- renderUI({
    HTML(paste0("<p>The testdata was taken from the following experiment: https://onlinelibrary.wiley.com/doi/abs/10.1002/mnfr.201500775. It investigates intestinal barrier integrity in diet induced obese
mice. </p>"))
  })
}