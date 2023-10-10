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
    explVar$log_pval <- -log10(explVar$pvalue)
    ggplot(data=explVar,aes(x=log_pval,y=rsquare,color=Variable))+
      geom_point(size = input$variation_point_size)+
      geom_text_repel(aes(label=paste0(Variable, "\n",
                                       "pvalue: ", round(pvalue,digits = 4), "\n",
                                       "rsquare: ", round(rsquare,digits = 4)
      )), size = input$variation_label_size) +
      xlab("-log10(pvalue)") +
      scale_color_manual(values=colorRampPalette(brewer.pal(9, input$namco_pallete))(length(explVar$Variable)))
  }
})

####confounding factors####

observeEvent(input$confounding_start,{
  if(!is.null(currentSet())){
    
    phylo <- vals$datasets[[currentSet()]]$phylo
    if(input$confounding_distance == "Unifrac"){
      distance <- as.matrix(vals$datasets[[currentSet()]]$unifrac_dist)
    }else{
      distance <- as.matrix(betaDiversity(phylo, method=1))
    } 

    meta <- as.data.frame(sample_data(phylo))
    #remove first column --> SampleID column
    meta[[sample_column]]<-NULL
    
    #calulate confounding matrix
    waiter_show(html = tagList(spin_rotating_plane(),"Looking for confounding factors ... "),color=overlay_color)
    tryCatch({
      vals$datasets[[currentSet()]]$confounder_table <- calculateConfounderTableNew(variables = meta,
                                                                                    distance = distance,
                                                                                    seed = seed,
                                                                                    ncores = ncores)
      waiter_hide()
    }, error=function(e){
      waiter_hide()
      print(e$message)
      showModal(errorModal(e$message))
    })
    
  }
})

output$confounding_heatmap <- renderPlot({
  if(!is.null(currentSet())){
    if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
      data<-vals$datasets[[currentSet()]]$confounder_table$table
      colnames(data)[which(colnames(data)==input$confounding_heatmap_type)] <- "plot_variable"
      data <- data[which(tested_variable %in% input$confounding_select_tested_variable),]
      ggplot(data, aes(x=tested_variable, y=possible_confounder, fill=plot_variable))+
        geom_tile(color="black")+
        xlab("tested variable")+
        ylab("possible confounder")+
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              axis.title = element_text(size=15),
              axis.text = element_text(size=input$confounding_label_size))+
        ggtitle("Heatmap of confounding factors")+
        scale_fill_viridis(discrete=input$confounding_heatmap_type!="pvalue")
    }
  }
}, height = 800)

output$confounding_PDF_download <- downloadHandler(
  filename=function(){
    paste("confounding_heatmap.pdf")
  },
  content=function(file){
    if(!is.null(currentSet())){
      if(!is.null(vals$datasets[[currentSet()]]$confounder_table)){
        data<-vals$datasets[[currentSet()]]$confounder_table$table
        colnames(data)[which(colnames(data)==input$confounding_heatmap_type)] <- "plot_variable"
        data <- data[which(tested_variable %in% input$confounding_select_tested_variable),]
        p<-ggplot(data, aes(x=tested_variable, y=possible_confounder, fill=plot_variable))+
          geom_tile(color="black")+
          xlab("tested variable")+
          ylab("possible confounder")+
          theme(axis.text.x = element_text(angle = 45, hjust=1),
                axis.title = element_text(size=15),
                axis.text = element_text(size=input$confounding_label_size))+
          ggtitle("Heatmap of confounding factors")+
          scale_fill_viridis(discrete=input$confounding_heatmap_type!="pvalue")
        
        ggsave(file, plot=p, device="pdf", width=20, height=12)
      }
    }
  }
)

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
