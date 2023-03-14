####associations####

observeEvent(input$associations_start,{
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_meta){
      message(paste0(Sys.time(), " - building SIAMCAT object ..."))
      waiter_show(html = tagList(spin_rotating_plane(),"Calculating differential associations ..."),color=overlay_color)
      
      phylo <- vals$datasets[[currentSet()]]$phylo
      
      if(input$associations_level != "OTU/ASV"){
        phylo_glom <- glom_taxa_custom(phylo, input$associations_level) #merge OTUs with same taxonomic level
        phylo <- phylo_glom$phylo_rank
        taxa_names(phylo) <- phylo_glom$taxtab[[input$associations_level]]
        rel_otu <- relAbundanceTo1(data.frame(otu_table(phylo), check.names=F))
      }else{
        rel_otu <- relAbundanceTo1(data.frame(otu_table(phylo), check.names=F))
      }
      
      # remove rows with NA in selected group
      meta <- data.frame(sample_data(phylo), check.names=F)
      meta <- meta[!is.na(meta[[input$associations_label]]),]
      # get check.name version of group name
      label <- make.names(colnames(data.frame(meta)))[which(colnames(meta) == input$associations_label)]
      
      tryCatch({
        # siamcat works only with 5 or more samples (https://git.embl.de/grp-zeller/SIAMCAT/-/blob/a1c662f343e99dabad4de024d4c993deba91bb0c/R/validate_data.r#L82)
        if(sum(meta[[input$associations_label]]==input$associations_case) <= 5){stop(siamcatNotEnoughSamplesError, call.=F)}
        
        s.obj <- siamcat(feat=rel_otu, meta=data.frame(meta), label= label, case= input$associations_case)
        s.obj <- filter.features(s.obj)
        s.obj <- check.associations(s.obj, alpha = input$associations_alpha)
        vals$datasets[[currentSet()]]$siamcat <- s.obj
        if(!any(associations(s.obj)[['p.adj']] < input$associations_alpha)){
          showModal(errorModal('Did not find any features with adjusted p-value below selected significance level. You might consider increasing it.'))
          waiter_hide()
        }
      },error = function(e){
        waiter_hide()
        print(e$message)
        showModal(errorModal(e$message))
      })
      
      waiter_hide()
    }
  }
})



output$associationsPlot <- renderPlot({
  if(!is.null(currentSet())){
    shiny::validate(
      need(!is.null(vals$datasets[[currentSet()]]$siamcat), "Please press 'Generate plot...' first to display the associations!")
    )
    sort.by <- c("p.val","fc","pr.shift")[which(input$associations_sort==c("p-value","fold-change","prevalence shift"))]
    panels <- c("fc","auroc","prevalence")[which(input$associations_panels==c("fold-change","AU-ROC","prevalence"))]
    shiny::validate(
      need(any(associations(vals$datasets[[currentSet()]]$siamcat)[["p.val"]]<input$associations_alpha), "Did not find any signficantly different features; try to adapt your significance level or change taxonomic level.")
    )
    # update stored siamcat object to extract significant features
    suppressMessages(SIAMCAT::association.plot(vals$datasets[[currentSet()]]$siamcat, fn.plot = NULL, prompt=F, verbose=0,
                                               max.show = input$assiciation_show_numer, 
                                               sort.by = sort.by,
                                               panels = panels,
                                               color.scheme = input$namco_pallete))
  }
}, height=800)

output$associationsPDF <- downloadHandler(
  filename=function(){"associations.pdf"},
  content = function(file){
    if(!is.null(vals$datasets[[currentSet()]]$siamcat)){
      s.obj <- vals$datasets[[currentSet()]]$siamcat
      sort.by <- c("p.val","fc","pr.shift")[which(input$associations_sort==c("p-value","fold-change","prevalence shift"))]
      panels <- c("fc","auroc","prevalence")[which(input$associations_panels==c("fold-change","AU-ROC","prevalence"))]
      
      suppressMessages(SIAMCAT::association.plot(vals$datasets[[currentSet()]]$siamcat, fn.plot = file, prompt=F, verbose=0,
                                                 max.show = input$assiciation_show_numer, 
                                                 sort.by = sort.by,
                                                 panels = panels,
                                                 color.scheme = input$namco_pallete))
      }
  }
)

output$associationsTable <- downloadHandler(
  filename = function(){'significant_features_associations.tab'},
  content = function(file){
    if(!is.null(vals$datasets[[currentSet()]]$siamcat)){
      df <- associations(vals$datasets[[currentSet()]]$siamcat)[,c('fc','p.val','auc','pr.shift','p.adj')]
      colnames(df) <- c('fold change','p-value','AUC','prevalence shift','FDR adjusted p-value')
      write.table(df, file=file, sep = '\t',row.names = T, quote = F)
    }
  }
)

####correlation####
corrReactive <- reactive({
  if(!is.null(currentSet())){
    
    waiter_show(html = tagList(spin_rotating_plane(),"Calculating pairwise correlations ..."),color=overlay_color)
    phylo <- vals$datasets[[currentSet()]]$phylo
    
    if(input$corrTaxLevel != "OTU/ASV"){
      phylo <- glom_taxa_custom(phylo, input$corrTaxLevel)$phylo_rank 
    }
    
    otu <- data.frame(phylo@otu_table, check.names = F)
    meta <- data.frame(phylo@sam_data, check.names = F)
    meta_numeric_all <- meta %>% select_if(is.numeric)
    meta_numeric <- as.data.frame(meta_numeric_all[,input$corrSelectGroups])
    colnames(meta_numeric) <- input$corrSelectGroups
    
    complete_data <- NULL
    tryCatch({
      if(input$corrIncludeTax){
        complete_data <- cbind(t(otu))
        meta_names <- NULL
      }
      if((length(meta_numeric) > 0 && input$corrIncludeTax) | length(meta_numeric) > 1){
        # fill NA entries 
        meta_fixed = as.data.frame(apply(meta_numeric, 2, fill_NA_INF.mean))
        
        # remove columns with only a single value
        meta_fixed[names(which(apply(meta_fixed,2,function(x){length(unique(x))})==1))] <- NULL
        if(is.null(complete_data)){
          complete_data <- cbind(meta_fixed)
        }else{
          complete_data <- cbind(complete_data, meta_fixed)
        }
        
        meta_names <- colnames(meta_fixed)
      }
      if(length(meta_numeric) <= 1 && !input$corrIncludeTax){
        waiter_hide()
        return(NULL)
      }
      
      # calculate correlation
      corr_df <- rcorr(as.matrix(complete_data, type = "pearson"))
      
      var_names <- row.names(corr_df$r)
      otu_names <- colnames(t(otu))
      
      waiter_update(html = tagList(spin_rotating_plane(),"Preparing plots ..."))
      
      corr_subset <- subsetCorrelation(input$corrIncludeTax, 
                                       ifelse(length(input$corrSelectGroups)>0,T,F),
                                       var_names, 
                                       otu_names, 
                                       meta_names, 
                                       corr_df,
                                       input$corrSignifCutoff,
                                       input$corrCorrelationCutoff)

      waiter_hide()
      if(is.null(dim(corr_subset$my_cor_matrix)) || dim(corr_subset$my_cor_matrix)[1] == 0){
        stop("No pairs found with the used correlation & p-value cutoffs!")
        return(NULL)
      }
      return(list(my_cor_matrix=corr_subset$my_cor_matrix,
                  my_pvl_matrix=corr_subset$my_pvl_matrix,
                  my_pairs=corr_subset$my_pairs))
    }, error = function(e){
      waiter_hide()
      print(e$message)
      showModal(errorModal(e$message))
    })
  }
})

output$corrPlot <- renderPlot({
  shiny::validate(
    need(!is.null(corrReactive()), "You need so select at least one meta-variable or include OTUs/ASVs in the options menu in order to display a plot!")
  )
  plot_correlation_custom(corrReactive()$my_cor_matrix, corrReactive()$my_pvl_matrix, input)
}, height=800)

#download as pdf
output$corrPlotPDF <- downloadHandler(
  filename = function(){"correlations.pdf"},
  content = function(file){
    if(!is.null(corrReactive())){
      pdf(file, width=12, height=12)
      plot_correlation_custom(corrReactive()$my_cor_matrix, corrReactive()$my_pvl_matrix, input)
      dev.off()
    }
  }
)


####topic modeling####

#here all objects and values needed for the plots of themetagenomics are created and stored in vals$datasets[[currentSet()]]$vis_out
observeEvent(input$themeta,{
  if(!is.null(currentSet())){
    #take otu table and meta file from user input
    phylo <- vals$datasets[[currentSet()]]$phylo
    otu <- round(data.frame(phylo@otu_table@.Data, check.names = F)*100)
    meta <- data.frame(phylo@sam_data, check.names = F)
    tax <- data.frame(phylo@tax_table, check.names = F)
    
    # check for "species" level column; this package does not work without it..
    if(dim(tax_table(phylo))[2]==6){
      message("No species level information in this dataset; adding a column with 'unkown' entries.")
      tax[["Species"]] <- "unkown"
      phylo <- merge_phyloseq(phylo, tax_table(as.matrix(tax)))
    }
    
    waiter_show(html = tagList(spin_rotating_plane(),"Finding Topics ..."),color=overlay_color)
    
    tryCatch({
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
      
      K=input$K
      sigma_prior = 0
      #use themetadata object to find K topics
      topics_obj <- find_topics(themetadata_object=obj,
                                K=K,
                                sigma_prior = sigma_prior)
      
      #measure relationship of covarite with samples over topics distribution from the STM
      topic_effects_obj <- est(topics_obj)
      
      #function_effects <- themetagenomics::est(functions_obj,topics_subset=3)
      
      class(topics_obj) <- "topics"
      vals$datasets[[currentSet()]]$vis_out <- prepare_vis2(topic_effects_obj)
      vals$datasets[[currentSet()]]$vis_out$K <- K
      vals$datasets[[currentSet()]]$vis_out$sigma <- sigma_prior
      vals$datasets[[currentSet()]]$vis_out$formula <- formula_char
      vals$datasets[[currentSet()]]$vis_out$refs <- refs
      vals$datasets[[currentSet()]]$topic_effects <- topic_effects_obj
      vals$datasets[[currentSet()]]$topics_table <- vals$datasets[[currentSet()]]$vis_out$tinfo
      waiter_hide()
    }, error =function(e){
      waiter_hide()
      print(e$message)
      showModal(errorModal(e$message))
    })
  }
})

output$downloadThemeta <- downloadHandler(
  filename = function(){"topics_assignment.csv"},
  content = function(file){
    if(!is.null(vals$datasets[[currentSet()]]$topics_table)){
      write.csv(vals$datasets[[currentSet()]]$topics_table, file = file, quote = F, sep = "\t")
    }
  } 
)


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
    shiny::validate(
      need(!is.null(vis_out), "Please press the 'Visualize topics!' button first to generate results!")
    )
    suppressWarnings(ggplotly(EST()$p_est,source='est_hover',tooltip=c('topic','est','lower','upper'))) #Error in UseMethod: no applicable method for 'plotly_build' applied to an object of class "shiny.tag"
  }
  
  
})

output$ord <- renderPlotly({
  if(!is.null(currentSet())){
    shiny::validate(
      need(!is.null(vals$datasets[[currentSet()]]$vis_out), "Please press the 'Visualize topics!' button first to generate results!"),
      need(!is.null(EST()), "")
    )
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    beta <- t(vis_out$beta)
    
    if (input$dist == 'hellinger'){
      d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'norm'),method='euclidean'),k=3,eig=TRUE)
    }else if (input$dist == 'chi2'){
      d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'chi.square'),method='euclidean'),k=3,eig=TRUE)
    }else if (input$dist == 'jsd'){
      d <- cmdscale(proxy::dist(beta,jsd),k=3,eig=TRUE)
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
})

show_topic <- reactiveValues(k=0)

observeEvent(input$k_in,{
  show_topic$k <- input$k_in
})

observeEvent(input$reset,{
  show_topic$k <- 0
  updateNumericInput(session,'k_in',value=0)
})

output$bar <- renderPlot({
  if(!is.null(currentSet())){
    shiny::validate(
      need(!is.null(vals$datasets[[currentSet()]]$vis_out), "Please press the 'Visualize topics!' button first to generate results!"),
      need(!is.null(REL()), "")
    )
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
})

output$themetaBarPDF <- downloadHandler(
  filename = function(){"themetagenomics_barplot.pdf"},
  content = function(file){
    if(!is.null(REL())){
      if (show_topic$k != 0){
        p_bar <- ggplot(data=REL()) +
          geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=.6) +
          geom_bar(aes_(~Term,~Freq),stat='identity',fill='darkred',color='white')
      } else{
        p_bar <- ggplot(data=REL()) +
          geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=1)
      }
      
      p_bar <- p_bar +
        coord_flip() +
        labs(x='',y='Frequency',fill='') +
        theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
              legend.position='bottom') +
        viridis::scale_fill_viridis(discrete=TRUE,drop=FALSE) +
        guides(fill=guide_legend(nrow=2))
      
      ggsave(file, p_bar, device="pdf", width = 10, height = 7)
    }
  }
)

output$corr <- renderForceNetwork({
  if(!is.null(currentSet())){
    shiny::validate(
      need(!is.null(vals$datasets[[currentSet()]]$vis_out), "Please press the 'Visualize topics!' button first to generate results!"),
      need(!is.null(vals$datasets[[currentSet()]]$topic_effects$topic_effects), "")
    )
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    topic_effects <- vals$datasets[[currentSet()]]$topic_effects$topic_effects
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
})


####time-series analysis####

timeSeriesReactive <- eventReactive(input$timeSeriesStart,{
  if(!is.null(currentSet())){
    
    phylo <- vals$datasets[[currentSet()]]$phylo
    waiter_show(html = tagList(spin_rotating_plane(),"Clustering and preparing plot ..."),color=overlay_color)
    
    if(vals$datasets[[currentSet()]]$has_meta){
      tryCatch({
        if(input$timeSeriesGroup != "" && input$timeSeriesGroup == input$timeSeriesBackground && input$timeSeriesClusterK == 0){stop(timeAndSampleGroupEqualError, call. = F)}
        if(input$timeSeriesGroup != "" && (input$timeSeriesGroup == input$timeSeriesMeanLine || input$timeSeriesBackground == input$timeSeriesMeanLine)){stop(timeSeriesEqualVariablesError, call.=F)}
        # apply k-means clustering to cluster OTUs
        # it makes no sense to cluster by diversity measures, since they are calculated per sample, not per OTU
        if(input$timeSeriesClusterK > 0 && input$timeSeriesMeasure %in% c("Abundance", "relative Abundance")){
          # use all OTU abundance values to cluster samples
          cluster_variables <- data.frame(t(phylo@otu_table@.Data), check.names = F)
          cluster_meta <- data.frame(phylo@sam_data, check.names = F)
          
          rownames(cluster_variables) <- sample_names(phylo)
          
          # run k-means
          km <- kmeans(cluster_variables, centers=input$timeSeriesClusterK, nstart=25)
          
          # add clusters as new variable to meta table
          clusters <- km$cluster
          cluster_meta[["cluster_group"]] <- as.character(clusters)
          
          # build new phyloseq object with new meta
          phylo <- merge_phyloseq(phylo, sample_data(cluster_meta))
          is_clustered <- T
        }else{
          cluster_variables <- NULL
          cluster_meta <- NULL
          is_clustered <- F
        }
        
        if(input$timeSeriesMeasure %in% c("Abundance", "relative Abundance")){
          # get sum of abundance per taxa level
          if(input$timeSeriesTaxa != "OTU/ASV"){
            phylo <- suppressMessages(glom_taxa_custom(phylo, input$timeSeriesTaxa)$phylo_rank) 
          }
          if(input$timeSeriesMeasure == "relative Abundance"){
            phylo <- transform_sample_counts(phylo, function(x) x/sum(x))
          }
          plot_df <- psmelt(phylo)
          # need to set the original column names to result of psmelt, as it adds . for whitespaces by default
          colnames(plot_df)[4:(3+ncol(sample_data(phylo)))] <- colnames(sample_data(phylo))
        }else{
          alphaTab <- vals$datasets[[currentSet()]]$alpha_diversity
          plot_df <- alphaTab
        }
        
        # calculate significant features 
        
        features_df <- over_time_serial_comparison(phylo, input$timeSeriesGroup, ifelse(input$timeSeriesClusterK == 0, input$timeSeriesBackground, "cluster_group"))
        waiter_hide()
        showModal(infoModal("Finished time-series analysis. Select one or more taxa to display the plot!"))
        return(list(plot_df=plot_df,
                    cluster_variables=cluster_variables,
                    cluster_meta=cluster_meta,
                    is_clustered=is_clustered,
                    features_df=features_df))  
      }, error=function(e){
        waiter_hide()
        print(e$message)
        showModal(errorModal(e$message))
      })
    }
  }
})

observe({
  if(!is.null(timeSeriesReactive())){
    possible_taxa <- sort(unique(timeSeriesReactive()$plot_df[["OTU"]]))
    updatePickerInput(session, "timeSeriesTaxaSelect", choices=possible_taxa)
  }
})

output$timeSeriesSignifFeatures <- renderPlot({
  req(timeSeriesReactive())
  df <- timeSeriesReactive()$features_df
  if(input$timeSeriesAdjPval == "default"){
    colnames(df)[which(colnames(df)=="pvalue_default")] <- "pvalue"
  }else{
    colnames(df)[which(colnames(df)=="pvalue_corrected")] <- "pvalue"
  }
  df$name <- factor(df$name, levels=df$name[order(df$pvalue)])
  ggplot(df, aes(x=-log10(pvalue), y=name))+
    geom_col(width = .75,na.rm = T)+
    ggtitle("Results of Friedman test for each taxonomic level & numeric meta-variable over the selected time-points and blocks.")+
    ylab("Name of taxa or meta-variable")+
    xlab("- log10(p-value)")
}, height=800)

timeSeriesPlotReactive <- reactive({
  req(timeSeriesReactive())
  plot_df <- as.data.table(timeSeriesReactive()$plot_df)
  colnames(plot_df)[which(colnames(plot_df)==input$timeSeriesGroup)] <- "reference" 
  colnames(plot_df)[which(colnames(plot_df)==input$timeSeriesBackground)] <- "sample_group"
  if(input$timeSeriesClusterK != 0 && !timeSeriesReactive()$is_clustered) {
    return(list(plot=NULL))
  }
  if(input$timeSeriesMeanLine == input$timeSeriesGroup){
    showModal(errorModal(error_message = "You cannot select the same group for the mean-line and as time-points!"))
    return(list(plot=NULL))
  }
  if(!input$timeSeriesMeanLine %in% c("NONE","")) colnames(plot_df)[which(colnames(plot_df)==input$timeSeriesMeanLine)] <- "time_series_mean"
  if(input$timeSeriesMeanLine == input$timeSeriesBackground) plot_df$time_series_mean <- plot_df[['sample_group']]
  p <- NULL
  title_text <- NULL
  vals$datasets[[currentSet()]]$has_ts_plot <- T
  # no need to select a taxa if diversity measure is chosen instead of abundance
  # --> no facet_wrap
  if (input$timeSeriesMeasure %in% c("Abundance", "relative Abundance")) {
    colnames(plot_df)[which(colnames(plot_df)=="Abundance")] <- "measure" 
    plot_df <- plot_df[plot_df[["OTU"]] %in% input$timeSeriesTaxaSelect,]
    shiny::validate(
      need(!is.null(input$timeSeriesTaxaSelect), "You need to select one or more taxa.")
    )

    p<-ggplot(plot_df, aes(x=as.character(reference), y=measure))+
      geom_line(aes(group=sample_group),alpha=0.35, color="black")+
      facet_wrap(~OTU, scales="free")+
      xlab(input$timeSeriesGroup)+
      ylab(input$timeSeriesMeasure)+
      labs(color=ifelse(input$timeSeriesClusterK > 0,"Cluster ID",input$timeSeriesMeanLine))+
      theme_bw()+
      ggtitle(paste0("Time-series analysis at ",input$timeSeriesTaxa," level; \n", 
                     input$timeSeriesBackground, " is displayed as small black lines in the back; \n",
                     "For ",input$timeSeriesMeanLine, " the mean ", input$timeSeriesMeasure, " over the time-points is displayed."))
    
  }else{
    colnames(plot_df)[which(colnames(plot_df)==input$timeSeriesMeasure)] <- "measure"
    p<-ggplot(plot_df, aes(x=as.character(reference), y=measure))+
      geom_line(aes(group=sample_group),alpha=0.35, color="black")+
      xlab(input$timeSeriesGroup)+
      ylab(input$timeSeriesMeasure)+
      labs(color=ifelse(input$timeSeriesClusterK > 0,"Cluster ID",input$timeSeriesMeanLine))+
      theme_bw()+
      ggtitle(paste0("Time-series analysis; \n", 
                     input$timeSeriesBackground, " is displayed as small black lines in the back; \n",
                     "For ",input$timeSeriesMeanLine, " the mean ", input$timeSeriesMeasure, " over the time-points is displayed."))
  }
  if(!is.null(p)){
    if(input$timeSeriesClusterK > 0){
      p <- p + 
        stat_summary(fun=mean, geom="line", size=input$timeSeriesLineSize, aes(group=cluster_group, color=cluster_group))+
        scale_color_manual(values=colorRampPalette(brewer.pal(9, input$namco_pallete))(input$timeSeriesClusterK))
    }
    if(input$timeSeriesMeanLine != "NONE"){
      p <- p + 
        stat_summary(fun=mean, geom="line", size=input$timeSeriesLineSize, aes(group=time_series_mean, color=as.character(time_series_mean)))+
        scale_color_manual(values=colorRampPalette(brewer.pal(9, input$namco_pallete))(length(unique(plot_df$time_series_mean))))
    }
    if(input$timeSeriesSampleHighlight != "NONE"){
      plot_df_filtered <- plot_df[plot_df$sample_group == input$timeSeriesSampleHighlight,]
      p <- p + 
        geom_line(aes(group=sample_group), data=plot_df_filtered, color=input$timeSeriesHighlightColor)
    }
    p <- p + scale_x_discrete(limits=input$timeSeriesTimePointOrder)
    vals$datasets[[currentSet()]]$has_ts_plot <- T
    return(list(plot=p))  
  }else{
    vals$datasets[[currentSet()]]$has_ts_plot <- F
    return(list(plot=NULL))
  }
  
})

output$timeSeriesPlot <- renderPlot({
  if(!is.null(timeSeriesPlotReactive()$plot)){
    timeSeriesPlotReactive()$plot  
  }
}, height = 800)


#download as pdf
output$timeSeriesPlotPDF <- downloadHandler(
  filename = function(){"time-series.pdf"},
  content = function(file){
    if(!is.null(timeSeriesPlotReactive())){
      ggsave(file, timeSeriesPlotReactive()$plot, device="pdf", width = 10, height = 7)
    }
  }
)

output$timeSeriesSignifTable <- downloadHandler(
  filename = function(){"significant_features.tsv"},
  content = function(file){
    if(!is.null(timeSeriesReactive())){
      write.table(timeSeriesReactive()$features_df, file = file, quote = F,sep = "\t", row.names = F)
    }
  }
)

output$timeSeriesOptimalClustersPlot <- renderPlot(
  if(!is.null(timeSeriesReactive())){
    if(input$timeSeriesClusterK > 0){
      df <- timeSeriesReactive()$cluster_variables
      p<-fviz_nbclust(df, kmeans, method="wss", k.max=20)
      p+geom_vline(xintercept=input$timeSeriesClusterK, color="red")
    }
  }
)

output$timeSeriesClusterSizePlot <- renderPlot(
  if(!is.null(timeSeriesReactive())){
    if(input$timeSeriesClusterK > 0 && timeSeriesReactive()$is_clustered){
      cluster_meta <- timeSeriesReactive()$cluster_meta
      colnames(cluster_meta)[which(colnames(cluster_meta)==input$timeSeriesBackground)] <- "time_points" 
      ggplot(cluster_meta, aes(y=cluster_group))+
        geom_bar(aes(fill=as.character(time_points)))+
        ggtitle("Composition and Size of individual clusters")+
        ylab("Cluster ID")+
        xlab("Number of samples in cluster")+
        scale_fill_discrete(name=input$timeSeriesBackground)
    }
  }
)

output$timeSeriesClusterContent <- renderDataTable({
  if(!is.null(timeSeriesReactive())){
    if(input$timeSeriesClusterK > 0 && timeSeriesReactive()$is_clustered){
      cluster_meta <- timeSeriesReactive()$cluster_meta
      dt <- cluster_meta[,c("SampleID","cluster_group")]
      colnames(dt) <- c("Sample ID","Cluster ID")
      dt
    }
  }
})

observe({
  x = input$horizonSubject
  if(input$horizonSubject!=""){
    updateSelectizeInput(session, 'horizonSubjectSelection', 
                         choices = c("", phylo@sam_data[,input$horizonSubject][[1]]), 
                         server=T)
  }
})

# look for biomehorizon data input
horizonData <- eventReactive(input$horizonStart, {
  waiter_show(html = tagList(spin_rotating_plane(), "Preparing horizon ..."),color=overlay_color)
  phylo <- vals$datasets[[currentSet()]]$phylo
  otu <- data.frame(OTU_ID=rownames(phylo@otu_table), phylo@otu_table, check.names = F)
  taxa <- data.frame(taxon_id=rownames(phylo@tax_table), phylo@tax_table, check.names = F)
  # prepare meta data
  meta_raw <- phylo@sam_data
  meta <- data.frame(meta_raw[,c(input$horizonSubject, input$horizonSample, input$horizonCollectionDate)], check.names = F)
  colnames(meta) <- c("subject", "sample", "time_point")
  # extract numbers from given time field
  meta$collection_date <- as.double(gsub("([0-9]+).*$", "\\1", meta$time_point))
  # remove NAs
  meta <- meta[complete.cases(meta),]
  waiter_hide()
  return(list(OTU=otu, META=meta, TAXA=taxa))
})

# biomehorizon plot
output$horizonPlot <- renderPlot({
  if(!is.null(horizonData())){
    hd <- horizonData()
    meta <- hd$META
    tps <- unique(meta[,c("collection_date", "time_point")])
    tps <- tps[order(tps$collection_date),]
    tax_data <- hd$TAXA[,input$horizonTaxaLevel]
    tax_data <- sapply(strsplit(tax_data, "__"), "[", 2)
    if(input$horizonTaxaSelect!=""){
      otulist <- input$horizonTaxaSelect
    } else {
      otulist <- NA
    }
    # select single subject
    if(input$horizonSubjectSelection!=""){
      subject <- input$horizonSubjectSelection
      subject_label <- paste0("Samples from Subject: ", subject)
      # samples per subject match number of time points selected
      if(max(table(meta$subject)) <= nrow(tps)) {
        data <- hd$OTU
      # more samples per subject then expected
      } else {
        # aggregate by input$horizonSubjectSelection
        meta$marker <- paste(meta$subject, meta$time_point, sep = "_")
        groups <- unique(meta$marker)
        data <- lapply(groups, function(g) rowSums(hd$OTU[,meta[meta$marker==g,"sample"]]))
        names(data) <- groups
        data <- data.frame(data, check.names = F)
        data <- data.frame(taxon_id=rownames(data), data, check.names = F)
        # collapse meta file to new samples
        meta$sample <- meta$marker
        meta$marker <- NULL
        meta <- unique(meta)
      }
    # use complete group if no specific patient is given
    } else {
      meta$subject <- "ALL"
      subject <- "ALL"
      subject_label <- "Samples from all subjects"
      # aggregate samples for each time point
      time_points <- unique(meta$time_point)
      data <- lapply(time_points, function(tp) rowSums(hd$OTU[,meta[meta$time_point==tp,"sample"]]))
      names(data) <- time_points
      data <- data.frame(data, check.names = F)
      data <- data.frame(taxon_id=rownames(data), data, check.names = F)
      meta <- unique(meta[,c("collection_date", "time_point", "subject")])
      meta$sample <- meta$time_point
    }

    paramList <- prepanel(otudata = data, metadata = meta, taxonomydata = tax_data,
                          subj = subject, facetLabelsByTaxonomy = input$horizonShowTaxa,
                          thresh_prevalence = as.numeric(input$horizonPrevalence), 
                          thresh_abundance = as.numeric(input$horizonAbundance), 
                          otulist = otulist, nbands = input$horizonNbands)
    horizonplot(paramList, aesthetics = horizonaes(title = "Microbiome Horizon Plot", 
                                                   xlabel = subject_label, 
                                                   ylabel = paste0("Taxa found in >", input$horizonPrevalence,"% of samples"), 
                                                   legendTitle = "Quartiles Relative to Taxon Median", 
                                                   legendPosition	= "bottom")) +
      scale_x_continuous(breaks = 1:length(unique(meta$collection_date)), labels = tps$time_point)
  }
})

observe({
  if(!is.null(currentSet())){
    if(input$timeSeriesClusterK > 0){
      updateSelectInput(session, "timeSeriesMeanLine", selected="NONE")
      shinyjs::hide("timeSeriesMeanLine")
    }else{
      shinyjs::show("timeSeriesMeanLine")
    }
    if(input$timeSeriesMeasure %in% c("Abundance", "relative Abundance")){
      shinyjs::show("timeSeriesTaxa")
      shinyjs::show("timeSeriesTaxaSelect")
      shinyjs::show("timeSeriesClusterK")
    }else{
      updateNumericInput(session, "timeSeriesClusterK", value = 0)
      shinyjs::hide("timeSeriesTaxa")
      shinyjs::hide("timeSeriesTaxaSelect")
      shinyjs::hide("timeSeriesClusterK")
    }
  }
})

####statistical tests#####
statTestReactive <- eventReactive(input$statTestStart, {
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_meta){
      
      if(input$statTestMethod == "Wilcoxon test"){vals$datasets[[currentSet()]]$has_wx_test<-F}
      if(input$statTestMethod == "Kruskal-Wallis test"){vals$datasets[[currentSet()]]$has_kw_test<-F}
      
      waiter_show(html = tagList(spin_rotating_plane(),"Finding significant taxa ..." ),color=overlay_color)
      
      phylo <- vals$datasets[[currentSet()]]$phylo
      phylo.rel <- transform_sample_counts(phylo, function(x) 100*x/sum(x))
      meta <- as.data.frame(phylo@sam_data, check.names=F)
      
      if(input$statTestcompLevel != "OTU/ASV"){
        phylo.rel <- glom_taxa_custom(phylo.rel, input$statTestcompLevel)$phylo_rank
      }
      all_taxa <- taxa_names(phylo.rel)
      # find all significant taxa for selected group
      tryCatch({
        signif <- lapply(all_taxa, function(i){
          group_vector <- meta[[input$statTestGroup]]
          abundance <- as.vector(otu_table(phylo.rel)[i,])
          df <- data.table(relative_abundance = abundance, group = group_vector, feature_name = i)
          # don't consider samples with NA as value for test group
          df <- df[complete.cases(df)]
          
          # perform wilcoxon test for all pairs of sample-groups; return if any pair is significantly different
          if(input$statTestMethod == "Wilcoxon test"){
            fit <- compare_means(relative_abundance~group, data=df, p.adjust.method = input$statTestPAdjust)  
            if(any(fit$p.adj < input$statTestCutoff)){
              # save pairs for which test was performed
              fit$pair <- paste0(fit$group1," vs. ", fit$group2)
              fit$pair_display <- paste0(fit$group1," vs. ", fit$group2, " (pval:", fit$p.adj,")")
              fit$tax_name <- i
              
              median_df <- df[, list(median=median(relative_abundance)), by=group]
              median_abundance <- median_df$median
              names(median_abundance) <- median_df$group
              
              mean_df <- df[, list(mean=mean(relative_abundance)), by=group]
              mean_abundance <- mean_df$mean
              names(mean_abundance) <- mean_df$group

              return(list(data=df,
                          fit_table=fit,
                          tax_name = i,
                          pval = fit$p.adj,
                          median_abundance = median_abundance,
                          mean_abundance = mean_abundance))
            }  
          }
          # perform KW test between all groups of selected meta variable
          else if(input$statTestMethod == "Kruskal-Wallis test"){
            fit <- kruskal.test(relative_abundance ~ group, data=df)
            fit$p.value <- p.adjust(fit$p.value, method=input$statTestPAdjust) 
            if(!is.na(fit$p.value)){
              if(fit$p.value < input$statTestCutoff){
                median_df <- df[, list(median=median(relative_abundance)), by=group]
                median_abundance <- median_df$median
                names(median_abundance) <- median_df$group
                
                mean_df <- df[, list(mean=mean(relative_abundance)), by=group]
                mean_abundance <- mean_df$mean
                names(mean_abundance) <- mean_df$group
                
                fit$tax_name <- i
                return(list(data=df,
                            fit=fit,
                            tax_name = i,
                            pval = fit$p.value,
                            median_abundance = median_abundance,
                            mean_abundance = mean_abundance,
                            test_method = input$statTestMethod))
              }  
            }
          }
        })
        names(signif) <- all_taxa
        signif <- Filter(Negate(is.null), signif)
        waiter_hide()
        if(length(signif) == 0){
          return(NULL)
        }
        if(input$statTestMethod == "Wilcoxon test"){vals$datasets[[currentSet()]]$has_wx_test<-T}
        if(input$statTestMethod == "Kruskal-Wallis test"){vals$datasets[[currentSet()]]$has_kw_test<-T}
        return(signif)  
      }, error=function(e){
        waiter_hide()
        print(e$message)
        showModal(errorModal(e$message))
        return(NULL)
      })
    }
  }
})

output$statTestDownloadTable <- downloadHandler(
  filename = function(){"significant_taxa_signifiance.tab"},
  content = function(file){
    if(!is.null(statTestReactive())){
      
      if(input$statTestMethod == "Wilcoxon test"){
        df <- rbindlist(lapply(statTestReactive(), function(x){return(x[["fit_table"]])}))
        df[,c('.y.','p.format','pair','pair_display')] <- NULL
      }
      if(input$statTestMethod == "Kruskal-Wallis test"){
        df <- rbindlist(lapply(statTestReactive(), function(x){return(x[["fit"]])}))
        df[,c('data.name')] <- NULL
        colnames(df) <- c('chi-squared','degrees.of.freedom','p.value','method','tax_name')
      }
      write.table(x = df,file = file, quote=F, sep="\t", row.names = F)
    }
  }
)

output$statTestAbundanceDownloadTable <- downloadHandler(
  filename = function(){"significant_taxa_abundance.tab"},
  content = function(file){
    if(!is.null(statTestReactive())){
      means <- data.frame(t(data.frame(lapply(statTestReactive(), function(x){return(x[["mean_abundance"]])}))))
      means$taxa <- rownames(means)
      means <- means %>% gather(group, mean_abundance, -taxa)
      
      medians <- data.frame(t(data.frame(lapply(statTestReactive(), function(x){return(x[["median_abundance"]])}))))
      medians$taxa <- rownames(medians)
      medians <- medians %>% gather(group, median_abundance, -taxa)
      
      df <- merge(means, medians, by=c('taxa', 'group'))
      write.table(x = df,file = file, quote=F, sep="\t", row.names = F)
    }
  }
)

# display significant taxa in select input 
observe({
  if(!is.null(statTestReactive())){
    signif_taxa <- names(statTestReactive())
    updateSelectInput(session, "statTestSignifPicker", choices = signif_taxa)
    if(input$statTestMethod == "Kruskal-Wallis test"){
      shinyjs::hide("statTestPairPicker")
    }else{
      shinyjs::show("statTestPairPicker")
    }
  }
})

# show pairs of comparison
observe({
  if(!is.null(statTestReactive())){
    pairs <- unlist(statTestReactive()[[input$statTestSignifPicker]]$fit_table$pair)
    if(!is.null(pairs)){
      names(pairs) = unlist(statTestReactive()[[input$statTestSignifPicker]]$fit_table$pair_display)
      order <- order(statTestReactive()[[input$statTestSignifPicker]]$fit_table$p)
      pairs <- pairs[order(order(order))]  # <-- maybe pick better variable name next time.. :p
      updatePickerInput(session, "statTestPairPicker", choices = pairs) 
    }
  }
})

statTestPlotReactive <- reactive({
  if(!is.null(statTestReactive())){
    selectedData <- statTestReactive()[[input$statTestSignifPicker]] # this is the taxa/OTU which was selected to plot
    if(!is.null(selectedData)){
      raw_data <- selectedData$data
      if(input$statTestMethod == "Wilcoxon test" && !is.null(vals$datasets[[currentSet()]]$has_wx_test)){
        test_data <- selectedData$fit_table
        selectedGroupsTable <- test_data[which(test_data$pair %in% input$statTestPairPicker),]
        selectedGroups <- unique(c(selectedGroupsTable$group1, selectedGroupsTable$group2)) # which groups will be displayed
        shiny::validate(
          need(length(selectedGroups)>0, "Please select at least one sup-group pair in the Options menu to display.")
        )
        if(!is.null(selectedGroups)){
          plot_data <- raw_data[raw_data$group %in% selectedGroups,]
          pairs <- sapply(selectedGroupsTable$pair, strsplit, split=" vs. ")
          p<-ggboxplot(plot_data, x="group", y="relative_abundance", 
                       title = paste0("Differential abundance for ",input$statTestSignifPicker))+
            stat_compare_means(comparisons = pairs) 
        }
      }else if(input$statTestMethod == "Kruskal-Wallis test" && !is.null(vals$datasets[[currentSet()]]$has_kw_test)){
        plot_data <- raw_data
        pval <- round(selectedData$fit$p.value,digits = 4)
        p<-ggboxplot(plot_data, x="group", y="relative_abundance",
                     title=paste0("Differential abundance for " ,input$statTestSignifPicker,"; p-value: ",pval))
      }
      p <- p+theme(axis.title=element_text(size=input$statTestLabelSize, face="bold"),
        axis.text=element_text(size=input$statTestLabelSize))

      return(list(plot=p))
    }
  }
})

statTestPlotReactiveAll <- reactive({
  if(!is.null(statTestReactive())){
    data <- statTestReactive()
    if(length(data) > 0){
      plot_data <- rbindlist(lapply(data, function(i){i[["data"]]}))
      p <- ggplot(plot_data, aes(x=group,y=relative_abundance))+
        geom_boxplot()+
        facet_wrap(~feature_name, scales="free")+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
        stat_compare_means()
      return(list(plot=p))
    }
  }
})

output$statTestPlot <- renderPlot({
  if(!is.null(statTestPlotReactive())){
    statTestPlotReactive()$plot
  }
})

#download as pdf
output$statTestPDF <- downloadHandler(
  filename = function(){"statistical_test.pdf"},
  content = function(file){
    if(!is.null(statTestPlotReactive())){
      ggsave(file, statTestPlotReactive()$plot, device="pdf", width = 10, height = 7)
    }
  }
)

output$statTestPDFall <- downloadHandler(
  filename = function(){"statistical_test_all.pdf"},
  content = function(file){
    if(!is.null(statTestPlotReactiveAll())){
      ggsave(file, statTestPlotReactiveAll()$plot, device = "pdf", width=20, height=20)
    }
  }
)
