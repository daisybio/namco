
#####################################
#    co-occurrence                  #
#####################################
# show histogram of all OTU values -> user can pick cutoff for binarization here
output$cutoffHist <- renderPlotly({
  if(!is.null(currentSet())){
    otu <- vals$datasets[[currentSet()]]$normalizedData
    dat <- log(as.data.frame(otu+1))
    
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
    otu <- vals$datasets[[currentSet()]]$normalizedData
    
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
  
  isolate(withProgress(message = 'Calculating Counts..', value = 0, {
    counts = generate_counts(OTU_table=otu <- vals$datasets[[currentSet()]]$normalizedData,
                             meta = data.frame(sample_data(vals$datasets[[currentSet()]]$phylo)),
                             group_column = input$groupCol,
                             cutoff = input$binCutoff,
                             fc = ifelse(input$useFC=="log2(fold-change)",T,F),
                             var1 = input$groupVar1,
                             var2 = input$groupVar2,
                             progress = T)
  }))
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
      
      out<-list(Nodes=Nodes,Links=Links)
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
    Links <- as.data.frame(cooccurrenceReactive()$Links)  
    Nodes <- as.data.frame(cooccurrenceReactive()$Nodes)
    
    forceNetwork(Links,Nodes,Source="source",Target="target",Value="valueToPlot",NodeID="name",
                 Nodesize="size",Group="group",linkColour=c("red","green")[(Links$value>0)+1],zoom=T,legend=T,
                 bounded=T,fontSize=12,fontFamily='sans-serif',charge=-25,linkDistance=100,opacity = 1)
    
  }
})

#####################################
#    themetagenomcis                #
#####################################

#here all objects and values needed for the plots of themetagenomics are created and stored in vals$datasets[[currentSet()]]$vis_out
#does not work with phyloseq-object (error if otu-table is of class phyloseq::otu_table)

observeEvent(input$themeta,{
  withProgress(message='Calculating Topics..',value=0,{
    if(!is.null(currentSet())){
      #take otu table and meta file from user input
      #otu <- data.frame(otu_table(vals$datasets[[currentSet()]]$phylo))
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
  
  if(!is.null(vis_out) & input$choose != "Please start calculation above first!"){
    suppressWarnings({
      
      #covariate <- input$choose
      covariate<-vals$datasets[[currentSet()]]$vis_out$covariates[1]
      
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


#Meinshausen-Buhlmann's #

observeEvent(input$se_mb_start,{
  
  if(!is.null(currentSet())){
    
    waiter_show(html = tagList(spin_rotating_plane(),"Running SPIEC-EASI (mb) ..." ),color=overlay_color)
    
    phylo <- vals$datasets[[currentSet()]]$phylo
    taxa <- tax_table(phylo)
    
    se_mb <- spiec.easi(phylo, method = "mb", lambda.min.ratio = input$se_mb_lambda.min.ratio, nlambda = input$se_mb_lambda, pulsar.params = list(rep.num=input$se_mb_repnumber, ncores =ncores,seed=seed))
    #pre-build graph object for phyloseq graph
    se_mb$ig <- adj2igraph(getRefit(se_mb), vertex.attr=list(name=taxa_names(phylo)))
    #pre-build grapg for interactive networkD3 graph
    nd3 <-igraph_to_networkD3(se_mb$ig, taxa)
    
    output$spiec_easi_mb_network              <- renderPlot({plot_network(se_mb$ig,phylo,type = "taxa",color=as.character(input$mb_select_taxa))})
    output$spiec_easi_mb_network_interactive  <- renderForceNetwork(forceNetwork(Links=nd3$links,Nodes=nd3$nodes,NodeID = "name",Group = as.character(input$mb_select_taxa),
                                                                                 zoom=T,legend=T,fontSize = 5,charge = -2,opacity = .9, height = 200,width = 100))
    
    waiter_hide()
  }
})

## Glasso ## 

observeEvent(input$se_glasso_start,{
  
  if(!is.null(currentSet())){
    
    waiter_show(html = tagList(spin_rotating_plane(),"Running SPIEC-EASI (glasso) ..." ),color=overlay_color)
    
    py <- vals$datasets[[currentSet()]]$phylo
    taxa <- tax_table(py)
    se_glasso <- isolate(spiec.easi(py, method = "glasso", lambda.min.ratio = input$glasso_mb_lambda.min.ratio, nlambda = input$glasso_mb_lambda, pulsar.params = list(rep.num=input$se_glasso_repnumber, ncores = ncores,seed=seed)))
    
    #pre-build graph object for phyloseq graph
    se_glasso$ig <- adj2igraph(getRefit(se_glasso), vertex.attr=list(name=taxa_names(py)))
    #pre-build grapg for interactive networkD3 graph
    nd3 <-igraph_to_networkD3(se_glasso$ig, taxa)
    
    output$spiec_easi_glasso_network              <- renderPlot({plot_network(se_glasso$ig,py,type = "taxa",color=as.character(input$glasso_select_taxa))})
    output$spiec_easi_glasso_network_interactive  <- renderForceNetwork(forceNetwork(Links=nd3$links,Nodes=nd3$nodes,NodeID = "name",Group = as.character(input$glasso_select_taxa),
                                                                                     zoom=T,legend=T,fontSize = 5,charge = -2,opacity = .9, height = 200,width = 100))
    
    
    waiter_hide()
    
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



