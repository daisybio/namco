
#####################################
#    co-occurrence                  #
#####################################
# show histogram of all OTU values -> user can pick cutoff for binarization here
output$cutoffHist <- renderPlotly({
  if(!is.null(currentSet())){
    otu <- vals$datasets[[currentSet()]]$rawData
    dat <- log(as.data.frame(otu))
    
    plot_ly(x=unlist(dat),type="histogram") %>%
      layout(xaxis=list(title="log(OTU abundance)"), yaxis = list(title="Frequency"),
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
    otu <- vals$datasets[[currentSet()]]$rawData
    
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
  
  waiter_show(html = tagList(spin_rotating_plane(),"Calculating Counts.."),color=overlay_color)
  
  counts = generate_counts(OTU_table = vals$datasets[[currentSet()]]$rawData,
                           meta = data.frame(sample_data(vals$datasets[[currentSet()]]$phylo)),
                           group_column = input$groupCol,
                           cutoff = input$binCutoff,
                           fc = ifelse(input$useFC=="log2(fold-change)",T,F),
                           var1 = input$groupVar1,
                           var2 = input$groupVar2,
                           progress = F)
  if (is.null(counts)){showModal(modalDialog(
    title="Warning!",
    "The column you chose for comparison only has one unique value! Please pick one with at least 2!",
    easyClose = T
  ))}
  vals$datasets[[currentSet()]]$counts <- counts
  vals$datasets[[currentSet()]]$network_params <- list(group_column = input$groupCol,
                                                       cutoff = input$binCutoff,
                                                       fc = ifelse(input$useFC=="log2(fold-change)",T,F),
                                                       var1 = input$groupVar1,
                                                       var2 = input$groupVar2)
  waiter_hide()
  
})

#network reactive
cooccurrenceReactive <- reactive({
  if(!is.null(currentSet())){
    if(!is.null(vals$datasets[[currentSet()]]$counts)){
      otu <- vals$datasets[[currentSet()]]$normalizedData
      counts = vals$datasets[[currentSet()]]$counts
      tax = tax_table(vals$datasets[[currentSet()]]$phylo)
      
      #remove all rows of counts where value is 0 -> no edge existing
      counts<-counts[counts$value!=0,]
      
      #if cutoff value is higher than maximal number of count values -> display all
      if(input$networkCutoff >= length(counts$value)){networkCutoff <- length(counts$value)} else {networkCutoff <- input$networkCutoff}
      
      Links = counts[order(abs(counts$value),decreasing=T)[1:networkCutoff],]
      colnames(Links) = c("source","target","value")
      Links$source = as.character(Links$source); Links$target = as.character(Links$target); Links$valueToPlot = abs(Links$value)
      Links$valueToPlot = norm10(Links$valueToPlot)+0.1
      
      Nodes = data.frame(name=unique(c(Links$source,Links$target)),group="")
      Nodes$size = rowSums(otu[Nodes$name,])/1000
      
      if(input$netLevel!="-") Nodes$group = substring(tax[Nodes$name,input$netLevel],4)
      Nodes$group[Nodes$group==""] = "unknown"
      
      Links$source = match(Links$source,Nodes$name)-1
      Links$target = match(Links$target,Nodes$name)-1
      
      network <- forceNetwork(Links,Nodes,Source="source",Target="target",Value="valueToPlot",NodeID="name",
                              Nodesize="size",Group="group",linkColour=c("red","green")[(Links$value>0)+1],zoom=T,legend=T,
                              bounded=T,fontSize=12,fontFamily='sans-serif',charge=-25,linkDistance=100,opacity = 1)
      out<-list(Nodes=Nodes,Links=Links,network=network)
      out
    }
  }
})

output$nodeDegree <- renderPlot({
  if(!is.null(cooccurrenceReactive())){
    Links <- as.data.frame(cooccurrenceReactive()$Links) 
    x<-table(c(Links$source,Links$target))
    names(x)<-NULL
    dat <- data.frame(degree=x)
    ggplot(data=dat,aes(x=log(degree.Freq)))+
      geom_histogram(bins = input$networkCutoff/25, stat="count")+
      ggtitle("Node Degree Plot for current network \n (Degree = # of edges coming out of a Node)")+
      xlab("log(Degree)")+ylab("Count")
  }
})


# network plot
output$basicNetwork <- renderForceNetwork({
  if(!is.null(cooccurrenceReactive())){
    cooccurrenceReactive()$network
  }
})

#download as pdf
output$basicNetworkPDF <- downloadHandler(
  filename = function(){"basic_network.html"},
  content = function(file){
    if(!is.null(cooccurrenceReactive())){
      saveNetwork(cooccurrenceReactive()$network, file, selfcontained = T)
    }
  }
)


#####################################
#    NetCoMi                        #
#####################################

##### observers for additional parameters ####
observe({
  if(input$diffNetworkMeasure == "spring"){
    shinyjs::show("diffNetworkAdditionalParamsSPRING.EASIDiv",anim = T)
    shinyjs::hide("diffNetworkAdditionalParamsSPARCCdiv")
  }else if(input$diffNetworkMeasure == "sparcc"){
    shinyjs::hide("diffNetworkAdditionalParamsSPRING.EASIDiv")
    shinyjs::show("diffNetworkAdditionalParamsSPARCCdiv",anim = T)
  }else{
    shinyjs::hide("diffNetworkAdditionalParamsSPRING.EASIDiv")
    shinyjs::hide("diffNetworkAdditionalParamsSPARCCdiv")
  }
})

observe({
  if(input$taxNetworkMeasure == "spring"){
    shinyjs::show("taxNetworkAdditionalParamsSPRING.EASIDiv",anim = T)
    shinyjs::hide("taxNetworkAdditionalParamsSPARCCdiv")
  }else if(input$taxNetworkMeasure == "sparcc"){
    shinyjs::hide("taxNetworkAdditionalParamsSPRING.EASIDiv")
    shinyjs::show("taxNetworkAdditionalParamsSPARCCdiv",anim = T)
  }else{
    shinyjs::hide("taxNetworkAdditionalParamsSPRING.EASIDiv")
    shinyjs::hide("taxNetworkAdditionalParamsSPARCCdiv")
  }
})

observe({
  if(input$compNetworkMeasure == "spring"){
    shinyjs::show("compNetworkAdditionalParamsSPRING.EASIDiv",anim = T)
    shinyjs::hide("compNetworkAdditionalParamsSPARCCdiv")
  }else if(input$compNetworkMeasure == "sparcc"){
    shinyjs::hide("compNetworkAdditionalParamsSPRING.EASIDiv")
    shinyjs::show("compNetworkAdditionalParamsSPARCCdiv",anim = T)
  }else{
    shinyjs::hide("compNetworkAdditionalParamsSPRING.EASIDiv")
    shinyjs::hide("compNetworkAdditionalParamsSPARCCdiv")
  }
})
##### single network #####

observeEvent(input$compNetworkCalculate, {
  if(!is.null(currentSet())){
    message(paste0(Sys.time()," - Starting single NetCoMi run ..."))
    waiter_show(html = tagList(spin_rotating_plane(),"Constructing network ..." ),color=overlay_color)
    phylo <- vals$datasets[[currentSet()]]$phylo

    measureParList = list()
    if(input$compNetworkMeasure == "sparcc"){
      measureParList<-append(measureParList, c(iter=input$compNetworkIter, inner_iter=compNetworkInnerIter, th=input$compNetworkTh))
    }
    if(input$compNetworkMeasure == "spring"){
      measureParList<-append(measureParList, c(nlambda=input$compNetworkNlambda, rep.num=input$compNetworkRepNum, lambda.min.ratio=input$compNetworkLambdaRatio))
    }
    
    tryCatch({
      net_con <- netConstruct(phylo,   
                              measure = input$compNetworkMeasure,
                              measurePar = measureParList,
                              normMethod = input$compNetworkNormMethod, 
                              zeroMethod = input$compNetworkzeroMethod,
                              sparsMethod = "none",
                              verbose = 0,
                              seed = seed,
                              cores=parallel::detectCores()) # use all of available cores
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      return()
    })
    
    waiter_update(html = tagList(spin_rotating_plane(),"Analyzing network ..."))
    tryCatch({
      net_ana <- netAnalyze(net_con, clustMethod = input$compNetworkClustMethod, weightDeg = FALSE, normDeg = FALSE,centrLCC = TRUE)
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      return(NULL)
    })
    
    if(input$compNetworkColor != "cluster"){
      tax_colors <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(tax_table(phylo)[,input$compNetworkColor])))
      featVecCol <- as.factor(tax_table(phylo)[,input$compNetworkColor])
      names(featVecCol) <- taxa_names(phylo)
      use_cluster_colors <- F
    }else{
      tax_colors <- NULL
      featVecCol <- NULL
      use_cluster_colors <- T
    }
    compNetworkList <- list(net_con=net_con, net_ana=net_ana, tax_colors=tax_colors, featVecCol=featVecCol, use_cluster_colors=use_cluster_colors)
    vals$datasets[[currentSet()]]$compNetworkList <- compNetworkList
    vals$datasets[[currentSet()]]$has_comp_nw <- T
    
    waiter_hide()
    message(paste0(Sys.time()," - Finished single NetCoMi run"))
  }
})


compNetworkPlotReactive <- reactive({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_comp_nw){
      tryCatch({
        p <- plot(vals$datasets[[currentSet()]]$compNetworkList$net_ana, 
                  sameLayout = T,
                  sameClustCol = T,
                  layout=input$compNetworkLayout,
                  layoutGroup = "union",
                  rmSingles = input$compNetworkRmSingles,
                  featVecCol = vals$datasets[[currentSet()]]$compNetworkList$featVecCol,
                  colorVec = vals$datasets[[currentSet()]]$compNetworkList$tax_colors,
                  nodeColor = ifelse(input$compNetworkColor == "cluster", "cluster", "feature"),
                  nodeTransp = 30,
                  nodeSize = input$compNetworkNodeSize,
                  nodeFilter = input$compNetworkNodeFilterMethod,
                  nodeFilterPar = input$compNetworkNodeFilterValue,
                  edgeInvisFilter = input$compNetworkEdgeFilterMethod,
                  edgeInvisPar =input$compNetworkEdgeFilterValue,
                  labelScale = T,
                  showTitle=T,
                  cexTitle=input$compNetworkTitleSize,
                  hubBorderCol  = "gray40",
                  title1=paste("Network on OTU level, edges calculated with",input$compNetworkMeasure))  
      }, error=function(e){
        p<-NULL
        showModal(errorModal(e$message))
      })
      return(p)
    }
  }
})


output$compNetwork <- renderPlot({
  if(!is.null(compNetworkPlotReactive())){
    compNetworkPlotReactive()
    if(input$compNetworkColor != "cluster" && !vals$datasets[[currentSet()]]$compNetworkList$use_cluster_colors){
      legend(x=-1.1,y=1.1, legend = levels(vals$datasets[[currentSet()]]$compNetworkList$featVecCol), 
             fill=NetCoMi:::colToTransp(vals$datasets[[currentSet()]]$compNetworkList$tax_colors, 30), 
             cex=1, bty="n",pt.cex=2.5, title=input$compNetworkColor, y.intersp = .7)  
    }
  }
}, height = 800)

output$compNetworkInteractive <- renderForceNetwork({
  if(!is.null(compNetworkPlotReactive())){
    p <- compNetworkPlotReactive()
    links <- data.frame(p$q1$Edgelist)
    nodes <- data.frame(node=rep(1:length(p$q1$Arguments$labels)), names = p$q1$Arguments$labels, size=p$q1$Arguments$vsize, color=p$q1$Arguments$color)
    # all nodes/edges indices need to start from 0
    links$from <- links$from-1
    links$to <- links$to-1
    nodes$node <- nodes$node-1
    
    forceNetwork(Links=links, Nodes=nodes, Source="from", Target="to", Value="weight", NodeID = "names", Nodesize = "size",
                 Group="color",zoom=T, bounded = T, opacity = 0.85, fontSize = 12, charge = -15)
  }
})

output$compNetworkSummary <- renderPrint(
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_comp_nw){
      summary(vals$datasets[[currentSet()]]$compNetworkList$net_ana)
    }
  }
)

#save plot as pdf
output$comp_networkPDF <- downloadHandler(
  filename = function(){"network.pdf"},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_comp_nw){
        pdf(file, width=9, height=7)
        plot(vals$datasets[[currentSet()]]$compNetworkList$net_ana, 
             sameLayout = T,
             sameClustCol = T,
             layout=input$compNetworkLayout,
             layoutGroup = "union",
             rmSingles = input$compNetworkRmSingles,
             featVecCol = vals$datasets[[currentSet()]]$compNetworkList$featVecCol,
             colorVec = vals$datasets[[currentSet()]]$compNetworkList$tax_colors,
             nodeColor = ifelse(input$compNetworkColor == "cluster", "cluster", "feature"),
             nodeTransp = 30,
             nodeSize = input$compNetworkNodeSize,
             nodeFilter = input$compNetworkNodeFilterMethod,
             nodeFilterPar = input$compNetworkNodeFilterValue,
             edgeInvisFilter = input$compNetworkEdgeFilterMethod,
             edgeInvisPar =input$compNetworkEdgeFilterValue,
             labelScale = T,
             showTitle=T,
             cexTitle=input$compNetworkTitleSize,
             hubBorderCol  = "gray40",
             title1=paste("Network on OTU level, edges calculated with",input$compNetworkMeasure))
        if(input$compNetworkColor != "cluster" && !vals$datasets[[currentSet()]]$compNetworkList$use_cluster_colors){
          legend(x=-1.1,y=1.1, legend = levels(vals$datasets[[currentSet()]]$compNetworkList$featVecCol), 
                 fill=NetCoMi:::colToTransp(vals$datasets[[currentSet()]]$compNetworkList$tax_colors, 30), 
                 cex=1, bty="n",pt.cex=2.5, title=input$compNetworkColor, y.intersp = .7)  
        }
        dev.off()
      }
    }
  }
)

##### differential network #####
observe({
  if(input$diffNetworkTaxaLevel=="OTU/ASV"){
    shinyjs::show("diffNetworkColor", anim=T)
  }else{
    shinyjs::hide("diffNetworkColor")
  }
})
observeEvent(input$diffNetworkCalculate, {
  if(!is.null(currentSet())){
    message(paste0(Sys.time()," - Starting differential NetCoMi run ..."))
    waiter_show(html = tagList(spin_rotating_plane(),"Constructing differential networks ..." ),color=overlay_color)
    phylo <- vals$datasets[[currentSet()]]$phylo
    if(input$diffNetworkTaxaLevel != "OTU/ASV"){
      tax_lst <- glom_taxa_custom(phylo, as.character(input$diffNetworkTaxaLevel))
      taxtable <- tax_lst$taxtab
      phylo <- tax_lst$phylo
      rownames(phylo@otu_table@.Data) <- taxtable[, as.character(input$diffNetworkTaxaLevel)]
    }
    phylo_split <- metagMisc::phyloseq_sep_variable(phylo, as.character(input$diffNetworkSplitVariable))
    
    measureParList = list()
    if(input$diffNetworkMeasure == "sparcc"){
      measureParList<-append(measureParList, c(iter=input$diffNetworkIter, inner_iter=diffNetworkInnerIter, th=input$diffNetworkTh))
    }
    if(input$diffNetworkMeasure == "spring"){
      measureParList<-append(measureParList, c(nlambda=input$diffNetworkNlambda, rep.num=input$diffNetworkRepNum, lambda.min.ratio=input$diffNetworkLambdaRatio))
    }
    phylo1 <- phylo_split[[1]]
    phylo2 <- phylo_split[[2]]
    
    tryCatch({
      net_con <- netConstruct(data = phylo1, 
                              data2 = phylo2,  
                              measure = input$diffNetworkMeasure,
                              measurePar = measureParList,
                              normMethod = input$diffNetworkNormMethod, 
                              zeroMethod = input$diffNetworkzeroMethod,
                              sparsMethod = "none",
                              verbose = 0,
                              seed = seed,
                              cores=round(parallel::detectCores())*0.5)
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      return()
    })

    
    waiter_update(html = tagList(spin_rotating_plane(),"Analyzing both networks ..."))
    tryCatch({
      net_ana <- netAnalyze(net_con, clustMethod = input$diffNetworkClustMethod, weightDeg = FALSE, normDeg = FALSE,centrLCC = TRUE)
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      net_ana <- NULL
    })
    
    waiter_update(html = tagList(spin_rotating_plane(),"Calculating differential networks ..."))
    tryCatch({
      diff_net <- diffnet(net_con, diffMethod = input$diffNetworkDiffMethod, cores = ncores)
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      diff_net <- NULL
    })
    
    if(input$diffNetworkColor != "cluster"){
      tax_colors1 <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(tax_table(phylo1)[,input$diffNetworkColor])))
      tax_colors2 <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(tax_table(phylo2)[,input$diffNetworkColor])))
      tax_colors <- list(tax_colors1, tax_colors2)
      featVecCol <- as.factor(tax_table(phylo1)[,input$diffNetworkColor])
      names(featVecCol) <- taxa_names(phylo1)
      use_cluster_colors <- F
    }else{
      tax_colors <- NULL
      featVecCol <- NULL
      use_cluster_colors <- T
    }
    
    waiter_update(html = tagList(spin_rotating_plane(),"Comparing differential networks ..."))
    net_comp <- netCompare(net_ana, permTest = FALSE, verbose = FALSE)
    
    diffNetworkList <- list(net_con=net_con, net_ana=net_ana, net_comp=net_comp, diff_net=diff_net, 
                            groups = names(phylo_split), tax_colors=tax_colors, featVecCol=featVecCol, 
                            use_cluster_colors=use_cluster_colors)
    vals$datasets[[currentSet()]]$diffNetworkList <- diffNetworkList
    vals$datasets[[currentSet()]]$has_diff_nw <- T
    
    waiter_hide()
    message(paste0(Sys.time()," - Finished differential NetCoMi run"))
  }
})

diffNetworkPlotReactive <- reactive({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_diff_nw){
      tryCatch({
        p<-plot(vals$datasets[[currentSet()]]$diffNetworkList$diff_net, 
                layout=ifelse(input$diffNetworkLayout=="Fruchterman-Reingold","spring",input$diffNetworkLayout),
                nodeTransp = 60,
                edgeInvisFilter = input$diffNetworkEdgeFilterMethod,
                edgeInvisPar = input$diffNetworkEdgeFilterValue,
                labelScale = T,
                cexLabels = input$diffNetworkLabelSize,
                hubBorderCol  = "gray40",
                legendGroupnames = c(vals$datasets[[currentSet()]]$diffNetworkList$groups[[1]],
                                     vals$datasets[[currentSet()]]$diffNetworkList$groups[[2]]),
                legendPos ="topright",
                cexTitle=input$diffNetworkTitleSize)  
      }, error=function(e){
        p<-NULL
        showModal(errorModal(e$message))
      })
      return(p)
    }
  }
})

output$diffNetwork <- renderPlot({
  if(!is.null(diffNetworkPlotReactive())){
    diffNetworkPlotReactive()
  }
}, height = 800)


#save plot as pdf
output$diff_networkPDF <- downloadHandler(
  filename = function(){"differential_network.pdf"},
  content = function(file){
    if(!is.null(diffNetworkPlotReactive())){
      pdf(file, width=9, height=7)
      plot(vals$datasets[[currentSet()]]$diffNetworkList$diff_net, 
              layout=ifelse(input$diffNetworkLayout=="Fruchterman-Reingold","spring",input$diffNetworkLayout),
              nodeTransp = 60,
              edgeInvisFilter = input$diffNetworkEdgeFilterMethod,
              edgeInvisPar = input$diffNetworkEdgeFilterValue,
              labelScale = T,
              cexLabels = input$diffNetworkLabelSize,
              hubBorderCol  = "gray40",
              legendGroupnames = c(vals$datasets[[currentSet()]]$diffNetworkList$groups[[1]],
                                   vals$datasets[[currentSet()]]$diffNetworkList$groups[[2]]),
              legendPos ="topright",
              cexTitle=input$diffNetworkTitleSize)
      dev.off()
    }
  }
)

##### 2 group network #####

observe({
  if(input$diffnet_tabs == "2 group network"){
    shinyjs::show("diffNetworkInteractiveSwitch")
    shinyjs::show("diffNetworkNodeFilterMethod")
    shinyjs::show("diffNetworkNodeFilterValue")
    shinyjs::show("diffNetworkRmSingles")
    shinyjs::show("diffNetworkNodeSize")
    updateSelectInput(session, label = "Choose method how to filter out edges (threshold: keep edges with weight of at least x; highestWeight: keep first x edges with highest weight)", "diffNetworkEdgeFilterMethod", choices = c("none", "threshold", "highestWeight"), selected = "highestWeight")
  }else{
    shinyjs::hide("diffNetworkInteractiveSwitch")
    shinyjs::hide("diffNetworkNodeFilterMethod")
    shinyjs::hide("diffNetworkNodeFilterValue")
    shinyjs::hide("diffNetworkRmSingles")
    shinyjs::hide("diffNetworkNodeSize")
    updateSelectInput(session, label ="Choose method how to filter out edges (highestDiff: the first x edges with highest absolute difference are shown)","diffNetworkEdgeFilterMethod", choices = c("none", "highestDiff"), selected="highestDiff")
  }
})

groupNetworkPlotReactive <- reactive({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_diff_nw){
      tryCatch({
        p<-plot(vals$datasets[[currentSet()]]$diffNetworkList$net_ana, 
                sameLayout = T,
                sameClustCol = T,
                layout=input$diffNetworkLayout,
                layoutGroup = "union",
                rmSingles = input$diffNetworkRmSingles,
                featVecCol = vals$datasets[[currentSet()]]$diffNetworkList$featVecCol,
                colorVec = vals$datasets[[currentSet()]]$diffNetworkList$tax_colors,
                nodeColor = ifelse(input$diffNetworkColor == "cluster", "cluster", "feature"),
                nodeTransp = 30,
                nodeSize = input$diffNetworkNodeSize,
                nodeFilter = input$diffNetworkNodeFilterMethod,
                nodeFilterPar = input$diffNetworkNodeFilterValue,
                edgeInvisFilter = input$diffNetworkEdgeFilterMethod,
                edgeInvisPar =input$diffNetworkEdgeFilterValue,
                labelScale = T,
                cexTitle=input$diffNetworkTitleSize,
                cexLabels = input$diffNetworkLabelSize,
                shortenLabels="none",
                groupNames = vals$datasets[[currentSet()]]$diffNetworkList$groups,
                showTitle=T,
                hubBorderCol  = "gray40")
      },error=function(e){
        p<-NULL
        showModal(errorModal(e$message))
      })

      return(p)
    }
  }
})

output$groupNetwork <- renderPlot({
  if(!is.null(groupNetworkPlotReactive())){
    groupNetworkPlotReactive()
    if(input$diffNetworkColor != "cluster" && !vals$datasets[[currentSet()]]$diffNetworkList$use_cluster_colors){
      legend(x=-1.1,y=1.1, legend = levels(vals$datasets[[currentSet()]]$diffNetworkList$featVecCol), 
             fill=NetCoMi:::colToTransp(vals$datasets[[currentSet()]]$diffNetworkList$tax_colors[[1]], 30), 
             cex=1, bty="n",pt.cex=2.5, title=input$diffNetworkColor, y.intersp = .7)  
    }
  }
}, height = 800)

output$groupNetworkInteractive1 <- renderForceNetwork({
  if(!is.null(groupNetworkPlotReactive())){
    p <- groupNetworkPlotReactive()
    links <- data.frame(p$q1$Edgelist)
    nodes <- data.frame(node=rep(1:length(p$q1$Arguments$labels)), names = p$q1$Arguments$labels, size=p$q1$Arguments$vsize, color=p$q1$Arguments$color)
    # all nodes/edges indices need to start from 0
    links$from <- links$from-1
    links$to <- links$to-1
    nodes$node <- nodes$node-1
    
    forceNetwork(Links=links, Nodes=nodes, Source="from", Target="to", Value="weight", NodeID = "names", Nodesize = "size",
                 Group="color",zoom=T, bounded = T, opacity = 0.85, fontSize = 12, charge = -5)
  }
})

output$groupNetworkInteractive2 <- renderForceNetwork({
  if(!is.null(groupNetworkPlotReactive())){
    p <- groupNetworkPlotReactive()
    links <- data.frame(p$q2$Edgelist)
    nodes <- data.frame(node=rep(1:length(p$q2$Arguments$labels)), names = p$q2$Arguments$labels, size=p$q2$Arguments$vsize, color=p$q2$Arguments$color)
    # all nodes/edges indices need to start from 0
    links$from <- links$from-1
    links$to <- links$to-1
    nodes$node <- nodes$node-1
    
    forceNetwork(Links=links, Nodes=nodes, Source="from", Target="to", Value="weight", NodeID = "names", Nodesize = "size",
                 Group="color",zoom=T, bounded = T, opacity = 0.85, fontSize = 12, charge = -5)
  }
})

output$groupNetworkSummary <- renderPrint(
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_diff_nw){
      summary(vals$datasets[[currentSet()]]$diffNetworkList$net_comp,
              groupNames=vals$datasets[[currentSet()]]$diffNetworkList$groups)
    }
  }
)

#save plot as pdf
output$group_networkPDF <- downloadHandler(
  filename = function(){"group_network.pdf"},
  content = function(file){
    if(!is.null(groupNetworkPlotReactive())){
      pdf(file, width=9, height=7)
      plot(vals$datasets[[currentSet()]]$diffNetworkList$net_ana, 
           sameLayout = T,
           sameClustCol = T,
           layout=input$diffNetworkLayout,
           layoutGroup = "union",
           rmSingles = input$diffNetworkRmSingles,
           featVecCol = vals$datasets[[currentSet()]]$diffNetworkList$featVecCol,
           colorVec = vals$datasets[[currentSet()]]$diffNetworkList$tax_colors,
           nodeColor = ifelse(input$diffNetworkColor == "cluster", "cluster", "feature"),
           nodeTransp = 30,
           nodeSize = input$diffNetworkNodeSize,
           nodeFilter = input$diffNetworkNodeFilterMethod,
           nodeFilterPar = input$diffNetworkNodeFilterValue,
           edgeInvisFilter = input$diffNetworkEdgeFilterMethod,
           edgeInvisPar = input$diffNetworkEdgeFilterValue,
           labelScale = T,
           shortenLabels="none",
           cexTitle=input$diffNetworkTitleSize,
           cexLabels = input$diffNetworkLabelSize,
           groupNames = vals$datasets[[currentSet()]]$diffNetworkList$groups,
           showTitle=T,
           hubBorderCol  = "gray40")
      if(input$diffNetworkColor != "cluster" && !vals$datasets[[currentSet()]]$diffNetworkList$use_cluster_colors){
        legend(x=-1.1,y=1.1, legend = levels(vals$datasets[[currentSet()]]$diffNetworkList$featVecCol), 
               fill=NetCoMi:::colToTransp(vals$datasets[[currentSet()]]$diffNetworkList$tax_colors[[1]], 30), 
               cex=1, bty="n",pt.cex=2.5, title=input$diffNetworkColor, y.intersp = .7)  
      }
      dev.off()
    }
  }
)


##### taxonomic network #####

observeEvent(input$taxNetworkCalculate, {
  if(!is.null(currentSet())){
    message(paste0(Sys.time()," - Starting taxonomic NetCoMi run ..."))
    waiter_show(html = tagList(spin_rotating_plane(),"Constructing taxonomic networks ..." ),color=overlay_color)
    
    phylo <- vals$datasets[[currentSet()]]$phylo
    tax_lst <- glom_taxa_custom(phylo, as.character(input$taxNetworkRank))
    taxtable <- tax_lst$taxtab
    phylo_rank <- tax_lst$phylo_rank
    rownames(phylo_rank@otu_table@.Data) <- taxtable[, as.character(input$taxNetworkRank)]
    
    measureParList = list()
    if(input$taxNetworkMeasure == "sparcc"){
      measureParList<-append(measureParList, c(iter=input$taxNetworkIter, inner_iter=taxNetworkInnerIter, th=input$taxNetworkTh))
    }
    if(input$taxNetworkMeasure == "spring"){
      measureParList<-append(measureParList, c(nlambda=input$taxNetworkNlambda, rep.num=input$taxNetworkRepNum, lambda.min.ratio=input$taxNetworkLambdaRatio))
    }
    
    tryCatch({
      net_con <- netConstruct(phylo_rank,  
                              measure = input$taxNetworkMeasure,
                              measurePar = measureParList,
                              normMethod = input$taxNetworkNormMethod, 
                              zeroMethod = input$taxNetworkzeroMethod,
                              sparsMethod = "none",
                              verbose = 0,
                              seed = seed,
                              cores=round(parallel::detectCores()*0.5))
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      return()
    })

    waiter_update(html = tagList(spin_rotating_plane(),"Analyzing taxonomic networks ..."))
    tryCatch({
      net_ana <- netAnalyze(net_con, clustMethod = input$taxNetworkClustMethod, weightDeg = FALSE, normDeg = FALSE,centrLCC = TRUE)
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      return(NULL)
    })
    
    taxNetworkList <- list(net_con=net_con, net_ana=net_ana, taxtable = taxtable, rank=input$taxNetworkRank, method=input$taxNetworkMeasure)
    vals$datasets[[currentSet()]]$taxNetworkList <- taxNetworkList
    vals$datasets[[currentSet()]]$has_tax_nw <- T
    
    waiter_hide()
    message(paste0(Sys.time()," - Finished taxonomic NetCoMi run"))
  }
})


taxNetworkPlotReactive <- reactive({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_tax_nw){
      tryCatch({
        rank <- as.character(vals$datasets[[currentSet()]]$taxNetworkList$rank)
        method <- as.character(vals$datasets[[currentSet()]]$taxNetworkList$method)
        p <-plot(vals$datasets[[currentSet()]]$taxNetworkList$net_ana, 
                 sameLayout = T, 
                 layout=input$taxNetworkLayout,
                 rmSingles = input$taxNetworkRmSingles,
                 nodeColor = "cluster",
                 nodeTransp = 60,
                 nodeSize = input$taxNetworkNodeSize,
                 nodeFilter = input$taxNetworkNodeFilterMethod,
                 nodeFilterPar = input$taxNetworkNodeFilterValue,
                 edgeInvisFilter = input$taxNetworkEdgeFilterMethod,
                 edgeInvisPar =input$taxNetworkEdgeFilterValue,
                 labelScale = T,
                 hubBorderCol  = "gray40",
                 shortenLabels = "none",
                 labelLength = 10,
                 cexNodes = 1.2,
                 cexLabels = 3.5,
                 cexHubLabels = 4,
                 cexTitle=input$taxNetworkTitleSize,
                 showTitle=T,
                 title1 = paste("Network on",rank," level, edges calculated with ", method))
      }, error=function(e){
        p<-NULL
        showModal(errorModal(e$message))
      })
      return(p)     
    }
  }
})

output$taxNetwork <- renderPlot({
  if(!is.null(taxNetworkPlotReactive())){
    taxNetworkPlotReactive()
  }
}, height = 800)

output$taxNetworkInteractive <- renderForceNetwork({
  p <- taxNetworkPlotReactive()
  links <- data.frame(p$q1$Edgelist)
  nodes <- data.frame(node=rep(1:length(p$q1$Arguments$labels)), names = p$q1$Arguments$labels, size=p$q1$Arguments$vsize, color=p$q1$Arguments$color)
  # all nodes/edges indices need to start from 0
  links$from <- links$from-1
  links$to <- links$to-1
  nodes$node <- nodes$node-1
  
  forceNetwork(Links=links, Nodes=nodes, Source="from", Target="to", Value="weight", NodeID = "names", Nodesize = "size",
               Group="color",zoom=T, bounded = T, opacity = 0.85, fontSize = 12, charge = -10)
})

output$taxNetworkSummary <- renderPrint(
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_tax_nw){
      summary(vals$datasets[[currentSet()]]$taxNetworkList$net_ana)
    }
  }
)

#save plot as pdf
output$tax_networkPDF <- downloadHandler(
  filename = function(){"taxonomic_network.pdf"},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_tax_nw){
        rank <- as.character(vals$datasets[[currentSet()]]$taxNetworkList$rank)
        method <- as.character(vals$datasets[[currentSet()]]$taxNetworkList$method)
        pdf(file, width = 9, height = 7)
        plot(vals$datasets[[currentSet()]]$taxNetworkList$net_ana, 
             sameLayout = T, 
             layout=input$taxNetworkLayout,
             rmSingles = input$taxNetworkRmSingles,
             nodeColor = "cluster",
             nodeTransp = 60,
             nodeSize = input$taxNetworkNodeSize,
             nodeFilter = input$taxNetworkNodeFilterMethod,
             nodeFilterPar = input$taxNetworkNodeFilterValue,
             edgeInvisFilter = input$taxNetworkEdgeFilterMethod,
             edgeInvisPar =input$taxNetworkEdgeFilterValue,
             labelScale = T,
             hubBorderCol  = "gray40",
             shortenLabels = "none",
             labelLength = 10,
             cexNodes = 1.2,
             cexLabels = 3.5,
             cexHubLabels = 4,
             cexTitle=input$taxNetworkTitleSize,
             showTitle=T,
             title1 = paste("Network on",rank," level, edges calculated with ", method)
        )
        dev.off()
      }
    }
  }
)
