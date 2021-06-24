
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
  
  withProgress(message = 'Calculating Counts..', value = 0, {
    counts = generate_counts(OTU_table = vals$datasets[[currentSet()]]$rawData,
                             meta = data.frame(sample_data(vals$datasets[[currentSet()]]$phylo)),
                             group_column = input$groupCol,
                             cutoff = input$binCutoff,
                             fc = ifelse(input$useFC=="log2(fold-change)",T,F),
                             var1 = input$groupVar1,
                             var2 = input$groupVar2,
                             progress = T)
  })
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

#download as html (does not work :( )
output$basicNetworkPDF <- downloadHandler(
  filename = function(){"basic_network.html"},
  content = function(file){
    if(!is.null(cooccurrenceReactive())){
      saveNetwork(cooccurrenceReactive()$network, file = file, selfcontained = T)
    }
  }
)

#####################################
#    themetagenomcis                #
#####################################

#here all objects and values needed for the plots of themetagenomics are created and stored in vals$datasets[[currentSet()]]$vis_out
observeEvent(input$themeta,{
  withProgress(message='Calculating Topics..',value=0,{
    if(!is.null(currentSet())){
      #take otu table and meta file from user input
      #otu <- as.data.frame(otu_table(vals$datasets[[currentSet()]]$phylo))
      otu <- round(vals$datasets[[currentSet()]]$normalizedData*100)
      meta <- vals$datasets[[currentSet()]]$metaData
      tax <- vals$datasets[[currentSet()]]$taxonomy
      
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
      sigma_prior = 0
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
      vals$datasets[[currentSet()]]$vis_out <- prepare_vis2(topic_effects_obj)
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
  
  if(!is.null(vis_out)){
    suppressWarnings({
      
      #covariate <- input$choose
      covariate<-vals$datasets[[currentSet()]]$vis_out$covariates[[1]]
      
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
        
        d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'norm'),method='euclidean'),k=3,eig=TRUE)
        
      }else if (input$dist == 'chi2'){
        
        d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'chi.square'),method='euclidean'),k=3,eig=TRUE)
        
      }else if (input$dist == 'jsd'){
        
        d <- cmdscale(proxy::dist(beta,jsd),k=3,eig=TRUE)   #woher kommt jsd?
        
      } else if (input$dist == 'tsne'){
        p <- 30
        d <- try(Rtsne::Rtsne(beta,k=,theta=.5,perplexity=p),silent=TRUE)
        while(class(d) == 'try-error'){
          p <- p-1
          d <- try(Rtsne::Rtsne(beta,k=3,theta=.5,perplexity=p),silent=TRUE)
        }
        if (p < 30) cat(sprintf('Performed t-SNE with perplexity = %s.\n',p))
        d$points <- d$Y
        d$eig <- NULL
        
      }else{
        
        d <- cmdscale(vegan::vegdist(beta,method=input$dist),k=3,eig=TRUE)
        
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
# observeEvent(event_data('plotly_click',source='ord_click'),{
#   
#   s <- event_data('plotly_click',source='ord_click')
#   
#   if (is.null(s)){
#     
#     show_topic$k <- 0
#     updateNumericInput(session,'k_in',value=0)
#     
#   }else{
#     
#     t_idx <- s$pointNumber + 1
#     updateNumericInput(session,'k_in',value=t_idx)
#     show_topic$k <- t_idx
#     
#   }
#   
# })

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
        
        # edge width (if edges are present)
        if(nrow(g_d3$links) == 0){
          return(NULL)
        }else{
          g_d3$links$edge_width <- 10*(.1+sapply(seq_len(nrow(g_d3$links)),function(r) vis_out$corr$poscor[g_d3$links$source[r]+1,g_d3$links$target[r]+1]))
        }
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
#    NetCoMi                        #
#####################################

##### observers for additional parameters ####
observe({
  if(input$diffNetworkMeasure == "spieceasi" || input$diffNetworkMeasure == "spring"){
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
  if(input$taxNetworkMeasure == "spieceasi" || input$taxNetworkMeasure == "spring"){
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
  if(input$compNetworkMeasure == "spieceasi" || input$compNetworkMeasure == "spring"){
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
    if(input$compNetworkMeasure == "spieceasi" || input$compNetworkMeasure == "spring"){
      measureParList<-append(measureParList, c(nlambda=input$compNetworkNlambda, rep.num=input$compNetworkRepNum, lambda.min.ratio=input$compNetworkLambdaRatio, seed=seed, ncores=ncores))
    }
    
    net_con <- netConstruct(phylo,   
                            measure = input$compNetworkMeasure,
                            measurePar = measureParList,
                            normMethod = input$compNetworkNormMethod, 
                            zeroMethod = input$compNetworkzeroMethod,
                            sparsMethod = "none",
                            verbose = 0,
                            seed = seed,
                            cores=round(parallel::detectCores())) # use 3/4 of available cores
    
    waiter_update(html = tagList(spin_rotating_plane(),"Analyzing network ..."))
    tryCatch({
      net_ana <- netAnalyze(net_con, clustMethod = input$compNetworkClustMethod, weightDeg = FALSE, normDeg = FALSE,centrLCC = TRUE)
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      return(NULL)
    })
    
    compNetworkList <- list(net_con=net_con, net_ana=net_ana)
    vals$datasets[[currentSet()]]$compNetworkList <- compNetworkList
    vals$datasets[[currentSet()]]$has_comp_nw <- T
    
    waiter_hide()
    message(paste0(Sys.time()," - Finished single NetCoMi run"))
  }
})


output$compNetwork <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_comp_nw){
      
      plot(vals$datasets[[currentSet()]]$compNetworkList$net_ana, 
           sameLayout = T,
           sameClustCol = T,
           layout=input$compNetworkLayout,
           layoutGroup = "union",
           rmSingles = input$compNetworkRmSingles,
           nodeColor = "cluster",
           nodeTransp = 60,
           nodeSize = input$compNetworkNodeSize,
           nodeFilter = input$compNetworkNodeFilterMethod,
           nodeFilterPar = input$compNetworkNodeFilterValue,
           edgeFilter = input$compNetworkEdgeFilterMethod,
           edgeFilterPar =input$compNetworkEdgeFilterValue,
           labelScale = T,
           showTitle=T,
           hubBorderCol  = "gray40",
           title1=paste("Network on OTU level, edges calculated with",input$compNetworkMeasure))
    }
  }
}, height = 800)

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
             nodeColor = "cluster",
             nodeTransp = 60,
             nodeSize = input$compNetworkNodeSize,
             nodeFilter = input$compNetworkNodeFilterMethod,
             nodeFilterPar = input$compNetworkNodeFilterValue,
             edgeFilter = input$compNetworkEdgeFilterMethod,
             edgeFilterPar =input$compNetworkEdgeFilterValue,
             labelScale = T,
             showTitle=T,
             hubBorderCol  = "gray40",
             title1=paste("Network on OTU level, edges calculated with",input$compNetworkMeasure))
        dev.off()
      }
    }
  }
)

##### differential network #####

observeEvent(input$diffNetworkCalculate, {
  if(!is.null(currentSet())){
    message(paste0(Sys.time()," - Starting differential NetCoMi run ..."))
    waiter_show(html = tagList(spin_rotating_plane(),"Constructing differential networks ..." ),color=overlay_color)
    phylo <- vals$datasets[[currentSet()]]$phylo
    phylo_split <- metagMisc::phyloseq_sep_variable(phylo, as.character(input$diffNetworkSplitVariable))
    
    measureParList = list()
    if(input$diffNetworkMeasure == "sparcc"){
      measureParList<-append(measureParList, c(iter=input$diffNetworkIter, inner_iter=diffNetworkInnerIter, th=input$diffNetworkTh))
    }
    if(input$diffNetworkMeasure == "spieceasi" || input$diffNetworkMeasure == "spring"){
      measureParList<-append(measureParList, c(nlambda=input$diffNetworkNlambda, rep.num=input$diffNetworkRepNum, lambda.min.ratio=input$diffNetworkLambdaRatio, seed=seed, ncores=ncores))
    }
    
    net_con <- netConstruct(data = phylo_split[[1]], 
                             data2 = phylo_split[[2]],  
                             measure = input$diffNetworkMeasure,
                             measurePar = measureParList,
                             normMethod = input$diffNetworkNormMethod, 
                             zeroMethod = input$diffNetworkzeroMethod,
                             sparsMethod = "none",
                             verbose = 0,
                             seed = seed,
                             cores=round(parallel::detectCores()*0.5))
    
    waiter_update(html = tagList(spin_rotating_plane(),"Analyzing differential networks ..."))
    tryCatch({
      net_ana <- netAnalyze(net_con, clustMethod = input$diffNetworkClustMethod, weightDeg = FALSE, normDeg = FALSE,centrLCC = TRUE)
    }, error=function(e){
      print(e$message)
      showModal(errorModal(e$message))
      return(NULL)
    })
    
    waiter_update(html = tagList(spin_rotating_plane(),"Comparing differential networks ..."))
    net_comp <- netCompare(net_ana, permTest = FALSE, verbose = FALSE)
    
    diffNetworkList <- list(net_con=net_con, net_ana=net_ana, net_comp=net_comp, groups = names(phylo_split))
    vals$datasets[[currentSet()]]$diffNetworkList <- diffNetworkList
    vals$datasets[[currentSet()]]$has_diff_nw <- T
    
    waiter_hide()
    message(paste0(Sys.time()," - Finished differential NetCoMi run"))
  }
})


output$diffNetwork <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_diff_nw){
      
      plot(vals$datasets[[currentSet()]]$diffNetworkList$net_ana, 
           sameLayout = T,
           sameClustCol = T,
           layout=input$diffNetworkLayout,
           layoutGroup = "union",
           rmSingles = input$diffNetworkRmSingles,
           nodeColor = "cluster",
           nodeTransp = 60,
           nodeSize = input$diffNetworkNodeSize,
           nodeFilter = input$diffNetworkNodeFilterMethod,
           nodeFilterPar = input$diffNetworkNodeFilterValue,
           edgeFilter = input$diffNetworkEdgeFilterMethod,
           edgeFilterPar =input$diffNetworkEdgeFilterValue,
           labelScale = T,
           groupNames = vals$datasets[[currentSet()]]$diffNetworkList$groups,
           showTitle=T,
           hubBorderCol  = "gray40")
           #title1=paste("Differential network on OTU level, edges calculated with",input$diffNetworkMeasure))
    }
  }
}, height = 800)

output$diffNetworkSummary <- renderPrint(
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_diff_nw){
      summary(vals$datasets[[currentSet()]]$diffNetworkList$net_comp,
              groupNames=vals$datasets[[currentSet()]]$diffNetworkList$groups)
    }
  }
)

#save plot as pdf
output$diff_networkPDF <- downloadHandler(
  filename = function(){"differential_network.pdf"},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_diff_nw){
        pdf(file, width=9, height=7)
        plot(vals$datasets[[currentSet()]]$diffNetworkList$net_ana, 
             sameLayout = T,
             sameClustCol = T,
             layout=input$diffNetworkLayout,
             layoutGroup = "union",
             rmSingles = input$diffNetworkRmSingles,
             nodeColor = "cluster",
             nodeTransp = 60,
             nodeSize = input$diffNetworkNodeSize,
             nodeFilter = input$diffNetworkNodeFilterMethod,
             nodeFilterPar = input$diffNetworkNodeFilterValue,
             edgeFilter = input$diffNetworkEdgeFilterMethod,
             edgeFilterPar =input$diffNetworkEdgeFilterValue,
             labelScale = T,
             groupNames = vals$datasets[[currentSet()]]$diffNetworkList$groups,
             showTitle=T,
             hubBorderCol  = "gray40")
             #paste("Differential network on OTU level, edges calculated with",input$diffNetworkMeasure))
        dev.off()
      }
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
    if(input$taxNetworkMeasure == "spieceasi" || input$diffNetworkMeasure == "spring"){
      measureParList<-append(measureParList, c(nlambda=input$taxNetworkNlambda, rep.num=input$taxNetworkRepNum, lambda.min.ratio=input$taxNetworkLambdaRatio, seed=seed, ncores=ncores))
    }
    
    net_con <- netConstruct(phylo_rank,  
                            measure = input$taxNetworkMeasure,
                            measurePar = measureParList,
                            normMethod = input$taxNetworkNormMethod, 
                            zeroMethod = input$taxNetworkzeroMethod,
                            sparsMethod = "none",
                            verbose = 0,
                            seed = seed,
                            cores=round(parallel::detectCores()*0.5))
    
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


output$taxNetwork <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_tax_nw){
      
      rank <- as.character(vals$datasets[[currentSet()]]$taxNetworkList$rank)
      method <- as.character(vals$datasets[[currentSet()]]$taxNetworkList$method)
      plot(vals$datasets[[currentSet()]]$taxNetworkList$net_ana, 
           sameLayout = T, 
           layout=input$taxNetworkLayout,
           rmSingles = input$taxNetworkRmSingles,
           nodeColor = "cluster",
           nodeTransp = 60,
           nodeSize = input$taxNetworkNodeSize,
           nodeFilter = input$taxNetworkNodeFilterMethod,
           nodeFilterPar = input$taxNetworkNodeFilterValue,
           edgeFilter = input$taxNetworkEdgeFilterMethod,
           edgeFilterPar =input$taxNetworkEdgeFilterValue,
           labelScale = T,
           hubBorderCol  = "gray40",
           shortenLabels = "none",
           labelLength = 10,
           cexNodes = 1.2,
           cexLabels = 3.5,
           cexHubLabels = 4,
           showTitle=T,
           title1 = paste("Network on",rank," level, edges calculated with ", method)
           )
    }
  }
}, height = 800)

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
             edgeFilter = input$taxNetworkEdgeFilterMethod,
             edgeFilterPar =input$taxNetworkEdgeFilterValue,
             labelScale = T,
             hubBorderCol  = "gray40",
             shortenLabels = "none",
             labelLength = 10,
             cexNodes = 1.2,
             cexLabels = 3.5,
             cexHubLabels = 4,
             showTitle=T,
             title1 = paste("Network on",rank," level, edges calculated with ", method)
        )
        dev.off()
      }
    }
  }
)
