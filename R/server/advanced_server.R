####heatmap#### 

# plot heatmap of OTU abundances per sample
output$abundanceHeatmap <- renderPlotly({
  if(!is.null(currentSet())){
    set.seed(seed)
    phylo <- vals$datasets[[currentSet()]]$phylo
    #check for unifrac distance --> (needs phylo tree file):
    if(!is.null(vals$datasets[[currentSet()]]$unifrac_dist)){
      #save generalized unifrac distance as global variable to use it for heatmap
      gunifrac_heatmap <<- as.dist(vals$datasets[[currentSet()]]$unifrac_dist)
      hm_distance <- if(input$heatmapDistance == "gunifrac") "gunifrac_heatmap" else input$heatmapDistance
      if(input$heatmapSample != "NULL"){
        plot_heatmap(phylo,method = input$heatmapOrdination,distance = hm_distance, sample.label = input$heatmapSample)
      }else{
        plot_heatmap(phylo,method = input$heatmapOrdination,distance = hm_distance)
      }
    }else{
      plotly_empty()
    }
    
  }
})

####random forest models####

#calculate confusion matrix using random forest aproach
rForestDataReactive <- eventReactive(input$forest_start,{
  if(!is.null(currentSet())){
    withProgress(message="Predicting random forest model...",value=0,{
      set.seed(input$forest_seed)
      incProgress(1/4,message="preparing data...")
      meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
      otu_t <- as.data.frame(t(otu_table(vals$datasets[[currentSet()]]$phylo)))
      
      if(input$forest_clr){
        otu_t <- clr(otu_t)
      }
      
      combined_data <- buildForestDataset(meta, otu_t, input)
      class_labels <- as.factor(combined_data$variable)
      
      incProgress(1/4,message="partitioning dataset & resampling...")
      #partition dataset in training+testing; percentage can be set by user
      inTraining <- createDataPartition(combined_data$variable, p = input$forest_partition, list = FALSE)
      #create training & testing partitions
      training <- combined_data[inTraining,]
      testing <- combined_data[-inTraining,]
      #parameters for training resampling (cross-validation)
      fitControl <- trainControl(classProbs=T,
                                 savePredictions=T,
                                 summaryFunction=twoClassSummary,
                                 method=input$forest_resampling_method,
                                 number=input$forest_cv_fold,
                                 repeats=input$forest_cv_repeats)
      
      #train random forest with training partition
      incProgress(1/4,message="training model...")
      
      if(input$forest_type == "random forest"){
        #tuning parameters:
        tGrid <- expand.grid(
          .mtry=extract(input$forest_mtry),
          .splitrule=input$forest_splitrule,
          .min.node.size=extract(input$forest_min_node_size)
        )
        if(input$forest_default){
          #use default values
          #model<-train(x=training,y=class_labels[inTraining],method = "ranger",trControl=fitControl,metric="ROC")
          model<-train(variable~.,
                       data=training,
                       method = "ranger",
                       trControl=fitControl,
                       metric="ROC",
                       importance="impurity",
                       num.threads=ncores)
        }else{
          model <- train(variable~.,
                         data=training,
                         method="ranger",
                         tuneGrid = tGrid,
                         importance=input$forest_importance,
                         trControl=fitControl,
                         num.trees=input$forest_ntrees,
                         metric="ROC",
                         num.threads=ncores)
        }
        
      }else if (input$forest_type == "gradient boosted model"){
        if(input$forest_default) {
          tGrid<-expand.grid(n.trees=100,interaction.depth=1,shrinkage=.1,n.minobsinnode=10)
        }else{
          tGrid <- expand.grid(
            n.trees=extract(input$gbm_ntrees),
            interaction.depth=extract(input$gbm_interaction_depth),
            shrinkage=extract(input$gbm_shrinkage),
            n.minobsinnode = extract(input$gbm_n_minobsinoode)
          )
        }
        model <- train(variable~.,
                       data=training,
                       method="gbm",
                       tuneGrid = tGrid,
                       trControl=fitControl,
                       metric="ROC")
      }
      
      #test model with testing dataset
      incProgress(1/4,message="testing model on test dataset...")
      predictions_model <- predict(model, newdata=testing)
      predictions_model_full <- predict(model, newdata=combined_data)
      con_matrix<-confusionMatrix(data=predictions_model, reference= class_labels[-inTraining])
      con_matrix_full<-confusionMatrix(data=predictions_model_full, reference= class_labels)
      return(list(cmtrx=con_matrix,cmtrx_full=con_matrix_full,model=model,class_labels=class_labels))
    })
  }
})


#density plot for continous variables
output$forest_sample_preview <- renderPlot({
  if(!is.null(currentSet())){
    meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
    if(is.numeric(meta[[input$forest_variable]])){
      g<-ggplot(data=meta,mapping=aes_string(x=input$forest_variable))+
        geom_density()+
        xlab(NULL)
      g
    }else{
      g<-ggplot(data=meta,mapping=aes_string(x=input$forest_variable))+
        geom_bar()+
        scale_fill_manual(values = )
      g
    }
  }
})

output$forest_continuous_preview<-renderPlot({
  if(!is.null(currentSet())){
    meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
    v = data.frame(cut(meta[[input$forest_variable]],breaks = c(-Inf,input$forest_continuous_slider,Inf),labels = c("low","high")))
    colnames(v)<-c("variable")
    g<-ggplot(data=v,aes(x=variable))+
      geom_bar()
    print(g)
  }
})

#output plots using the reactive output of random forest calculation
output$forest_con_matrix <- renderPlot({
  if(!is.null(rForestDataReactive())){
    draw_confusion_matrix(rForestDataReactive()$cmtrx)
  }
})

output$forest_con_matrix_full <- renderPlot({
  if(!is.null(rForestDataReactive())){
    draw_confusion_matrix(rForestDataReactive()$cmtrx_full)
  }
})

#ROC curve for model
output$forest_roc <- renderPlot({
  if(!is.null(rForestDataReactive())){
    res<-evalm(rForestDataReactive()$model)
    res$roc
  }
})


output$forest_roc_cv <- renderPlot({
  if(!is.null(rForestDataReactive())){
    ldf <- lift_df(rForestDataReactive()$model,rForestDataReactive()$model$levels[1])
    ggplot(ldf) +
      geom_line(aes(1 - Sp, Sn, color = fold)) +
      scale_color_discrete(guide = guide_legend(title = "Fold"))+
      xlab("x")+ylab("y")
  }
})

output$forest_top_features <- renderPlot({
  if(!is.null(rForestDataReactive())){
    plot(varImp(rForestDataReactive()$model),top=input$top_x_features)
  }
})

output$forest_save_model <- downloadHandler(
  filename=function(){paste("random_forest_model.rds")},
  content = function(file){
    if(!is.null(rForestDataReactive())){
      saveRDS(rForestDataReactive()$model,file)
    }
  }
)

#upload file to use model to predict variable for new sample(s)
rForestPrediction <- eventReactive(input$forest_upload,{
  if(is.null(input$forest_upload_file)){showModal(errorModal(error_message = fileNotFoundError))}
  if(is.null(rForestDataReactive())){showModal(errorModal(error_message = "Please calculate the model first."))}
  else{
    new_sample <- read.csv(input$forest_upload_file$datapath,header=T,sep="\t",row.names=1,check.names = F)
    #transpose new_sample; same orientation as otu in model building 
    new_sample<-t(new_sample)
    if(is.null(new_sample)){showModal(errorModal(error_message = fileEmptyError))}
    else{
      model <- rForestDataReactive()$model
      #only use columns for prediction, which were used for model building
      model_predictors <- setdiff(colnames(model$trainingData),".outcome")
      if(!all(model_predictors %in% colnames(new_sample))){showModal(errorModal(error_message = inconsistentColumnsForest))}
      else{
        pred <- predict(model,newdata=new_sample)
        df <- data.frame(row.names = rownames(new_sample),prediction = pred)
      }
    }
  }
  return (df)
}) 

#output table for prediction of variables of new sample
output$forest_prediction <- renderTable({
  if(!is.null(rForestPrediction())){
    rForestPrediction()
  }
},rownames = T)

#observer for random Forest menu items  
observe({
  #switch between different versions of advanced options, depending on type of model
  if(input$forest_type == "random forest"){
    shinyjs::show("ranger_advanced",anim = T)
    shinyjs::hide("gbm_advanced")
  }else if (input$forest_type == "gradient boosted model"){
    shinyjs::hide("ranger_advanced")
    shinyjs::show("gbm_advanced",anim=T)
  }
  
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_meta){
      meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
      group_columns <- setdiff(colnames(meta),"SampleID")
      #for continous_slider; update input based on varibale chosen in forest_variable
      if(is.numeric(meta[[input$forest_variable]])){
        cut<-switch (input$forest_continuous_radio,
                     Median = median(meta[[input$forest_variable]],na.rm=T),
                     Mean = mean(meta[[input$forest_variable]],na.rm=T)
        )
        updateSliderInput(session,"forest_continuous_slider",min=0,max=max(meta[[input$forest_variable]],na.rm = T),value = cut)
        shinyjs::show("forest_continuous_options",anim=T)
      }else{
        updateSliderInput(session,"forest_continuous_slider",min=0,max=1)
        shinyjs::hide("forest_continuous_options")
        #test for more than 2 groups in variable
        if(length(unique(meta[[input$forest_variable]])) > 2){
          shinyjs::show("forest_covariable",anim = T)
        }else{
          shinyjs::hide("forest_covariable")
        }
      }
      #for selecting OTUs which will not be used in rForest calculation
      otu_t<-as.data.frame(t(otu_table(vals$datasets[[currentSet()]]$phylo)))
      updateSelectizeInput(session, "forest_exclude",choices=sort(colnames(otu_t)), server = T)
      updateSelectInput(session, "forest_covariable",choices=unique(meta[[input$forest_variable]]))
      
      #for features to include in model calculation; remove feature from list, which will be predicted 
      features<-group_columns[group_columns != input$forest_variable]
      updateSelectInput(session,"forest_features",choices = features) 
    }
  }
})

#javascript show/hide toggle for advanced options
shinyjs::onclick("forest_toggle_advanced",shinyjs::toggle(id="forest_advanced",anim = T))
#show/hide exclude OTU-option if OTU abundances are to be used for model building
shinyjs::onclick("forest_otu",shinyjs::toggle(id="forest_exclude",anim = T))

####picrust2 ####

observe({
  if (!is.null(currentSet())){
    if (vals$datasets[[currentSet()]]$is_fastq){
      shinyjs::hide("fastaFile")
    }else{
      shinyjs::show("fastaFile")
    }
  }
})

observeEvent(input$picrust2Start,{
  if(!is.null(currentSet())){
    print(paste0(Sys.time(), " - Starting picrust2 analysis ..."))
    
    waiter_show(html = tagList(spin_rotating_plane(),"Running Picrust2 ..." ),color=overlay_color)
    
    phylo <- vals$datasets[[currentSet()]]$phylo
    vals$datasets[[currentSet()]]$picrust_output <- NULL    # reset old picrust-output variable
    shinyjs::hide("download_picrust_div")
    
    foldername <- sprintf("/picrust2_%s", digest::digest(phylo))  # unique folder name for this output
    outdir <- paste0(tempdir(),foldername)
    if (dir.exists(outdir)){unlink(outdir, recursive = T)}    # remove output directory if it exists
    dir.create(outdir)
    
    biom_file = paste0(outdir,"/biom_picrust.biom")
    biom <- make_biom(data=otu_table(phylo))
    write_biom(biom, biom_file)
    print(paste0(Sys.time(), " - Wrote biom-file: ", biom_file))
    
    # use fasta file of dada2 pipeline if available
    if (vals$datasets[[currentSet()]]$is_fastq){
      fasta_file <- paste0(outdir,"/seqs.fasta")
      writeXStringSet(refseq(phylo), fasta_file)
      print(paste0(Sys.time(), " - Using dada2-generated fasta file:", fasta_file))
    }else{
      fasta_file <- input$fastaFile$datapath
      print(paste0(Sys.time(), " - Using user-uploaded fasta file: ", fasta_file))
    }
    
    picrust_outdir <- paste0(outdir,"/picrust2out")       # this is the name of the final output directory of this picrust run
    command = paste0("/opt/anaconda3/bin/conda run -n picrust2 picrust2_pipeline.py -s ",fasta_file," -i ",biom_file, " -o ", picrust_outdir, " -p", ncores)
    #command = paste0("/home/alex/anaconda3/bin/conda run -n picrust2 picrust2_pipeline.py -s ",fasta_file," -i ",biom_file, " -o ", picrust_outdir, " -p", ncores)
    out <- system(command, wait = TRUE)
    shinyjs::show("download_picrust_div", anim = T)
    vals$datasets[[currentSet()]]$picrust_output <- picrust_outdir 
    
    print(paste0(Sys.time(), " - Finished picrust2 run; output in: ", picrust_outdir))
    waiter_hide()
  }
})

output$download_picrust <- downloadHandler(
  filename = function(){
    paste("picrust2_output.zip")
  },
  content = function(file){
    if(!is.null(currentSet())){
      if (!is.null(vals$datasets[[currentSet()]]$picrust_output)){
        withProgress(message("Creating zip archive of picrust result-files...", value=0), {
          zip(zipfile = file, files=list.files(vals$datasets[[currentSet()]]$picrust_output, full.names = T))
        })
      }
    }
  }, 
  contentType = "application/zip"
)

