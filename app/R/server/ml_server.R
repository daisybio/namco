####random forest models####

#calculate confusion matrix using random forest aproach
observeEvent(input$forest_start,{
  if(!is.null(currentSet())){
    
    set.seed(seed)
    message(paste0(Sys.time(), " - Starting random forest ..."))
    waiter_show(html = tagList(spin_rotating_plane(),"Preparing Data ..." ),color=overlay_color)
    meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
    otu_t <- as.data.frame(t(otu_table(vals$datasets[[currentSet()]]$phylo)))
    
    combined_data <- buildForestDataset(meta, otu_t, input)
    class_labels <- as.factor(combined_data$variable)
    
    waiter_update(html = tagList(spin_rotating_plane(),"Partitioning Dataset & Resampling ..." ))
    
    # create initial split
    split <- rsample::initial_split(combined_data,
                                    prop = input$forest_partition,
                                    strata = variable)
    
    # define training and testint
    train_data <- rsample::training(split)
    test_data <- rsample::testing(split)
    
    # predict variable with all available features
    train_recipe <- recipe(variable ~ ., data = train_data)
    
    # setup RF model
    tuning_specs <- parsnip::rand_forest(
      mtry = tune(), 
      trees = tune(), 
      min_n = tune()
    ) |>
      set_mode('classification') |>
      set_engine('ranger', importance = input$forest_importance)
    
    # complete workflow
    tuning_workflow <- workflow() |>
      add_recipe(train_recipe) |>
      add_model(tuning_specs)
    
    # setup CV
    training_folds <- rsample::vfold_cv(data = train_data, v = input$forest_cv_fold, repeats = input$forest_cv_repeats, strata = variable)
    
    # tuning grid
    tuning_grid <- NULL
    if(input$forest_tuning_type == 'random'){
      tuning_grid <- dials::grid_random(
        trees(range = input$forest_ntrees), 
        min_n(range = input$forest_min_node_size), 
        mtry(range = input$forest_mtry), 
        size = input$forest_random_parameters
      )
    }else if(input$forest_tuning_type == 'regular'){
      tuning_grid <- dials::grid_regular(
        trees(range = input$forest_ntrees), 
        min_n(range = input$forest_min_node_size), 
        mtry(range = input$forest_mtry), 
        levels = input$forest_regular_parameters
      )
    }
    
    library(future)
    future::plan(future::multisession(workers = ncores))
    
    waiter_update(html = tagList(spin_rotating_plane(),"Trying out different hyperparameter using cross-validation ..." ))
    # hyperparameter tuning using CV
    tuning_results <- tune_grid(
      tuning_workflow, 
      resamples = training_folds,
      grid=tuning_grid,
      metrics = metric_set(f_meas, accuracy),
      control = control_grid(save_pred = TRUE, verbose = TRUE)
    )
    
    # collect best hyperparameters
    tuning_metrics <- tuning_results |>
      collect_metrics()
    

    
    
    
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
    waiter_update(html = tagList(spin_rotating_plane(),"Training Model ..." ))
    
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
        model<-caret::train(variable~.,
                            data=training,
                            method = "ranger",
                            trControl=fitControl,
                            metric="ROC",
                            importance="impurity",
                            num.threads=ncores)
      }else{
        model <- caret::train(variable~.,
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
      model <- caret::train(variable~.,
                            data=training,
                            method="gbm",
                            tuneGrid = tGrid,
                            trControl=fitControl,
                            metric="ROC")
    }
    
    #test model with testing dataset
    waiter_update(html = tagList(spin_rotating_plane(),"Testing Model ..." ))
    predictions_model <- caret::predict.train(model, newdata=testing)
    predictions_model_full <- caret::predict.train(model, newdata=combined_data)
    con_matrix<-confusionMatrix(data=predictions_model, reference= class_labels[-inTraining])
    con_matrix_full<-confusionMatrix(data=predictions_model_full, reference= class_labels)
    
    message(paste0(Sys.time(), " - Finished random forest"))
    waiter_hide()
    showModal(forestReadyModal)
    rf_lst <- list(cmtrx=con_matrix,cmtrx_full=con_matrix_full,model=model,class_labels=class_labels)
    vals$datasets[[currentSet()]]$rf_lst <- rf_lst
    vals$datasets[[currentSet()]]$has_rf <- T
    
  }
})

forestReadyModal <- modalDialog(
  title = "Finished Random Forest!",
  "Scroll down to check the performance of your model!",
  easyClose = T, size="s"
)

#density plot for continuous variables
output$forest_sample_preview <- renderPlot({
  if(!is.null(currentSet())){
    meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
    if(is.numeric(meta[[input$forest_variable]])){
      g<-ggplot(data=meta,mapping=aes_string(x=input$forest_variable))+
        geom_density()+
        xlab(NULL)+
        theme_bw()
      g
    }else{
      g<-ggplot(data=meta,mapping=aes_string(x=input$forest_variable, fill=input$forest_variable))+
        geom_bar()+
        theme_bw()+
        scale_fill_brewer(palette = input$namco_pallete)
      g
    }
  }
})

output$forest_continuous_preview<-renderPlot({
  if(!is.null(currentSet())){
    meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
    v = data.frame(cut(meta[[input$forest_variable]],breaks = c(-Inf,input$forest_continuous_slider,Inf),labels = c("low","high")))
    colnames(v)<-c("variable")
    g<-ggplot(data=v,aes(x=variable, fill=variable))+
      geom_bar()+
      theme_bw()+
      scale_fill_brewer(palette = input$namco_pallete)
    print(g)
  }
})

#output plots using the reactive output of random forest calculation
output$forest_con_matrix <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_rf){
      draw_confusion_matrix(vals$datasets[[currentSet()]]$rf_lst$cmtrx) 
    }
  }
})

#download as pdf
output$forest_con_matrixPDF <- downloadHandler(
  filename = function(){"randomForest_confusion_matrix.pdf"},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_rf){
        pdf(file, width=8, height=6)
        draw_confusion_matrix(vals$datasets[[currentSet()]]$rf_lst$cmtrx) 
        dev.off()
      }
    }
  }
)

output$forest_con_matrix_full <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_rf){
      draw_confusion_matrix(vals$datasets[[currentSet()]]$rf_lst$cmtrx_full) 
    }
  }
})

#download as pdf
output$forest_con_matrix_fullPDF <- downloadHandler(
  filename = function(){"randomForest_confusion_matrix_full.pdf"},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_rf){
        pdf(file, width=8, height=6)
        draw_confusion_matrix(vals$datasets[[currentSet()]]$rf_lst$cmtrx_full) 
        dev.off()
      }
    }
  }
)

#ROC curve for model
output$forest_roc <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_rf){
      res<-evalm(vals$datasets[[currentSet()]]$rf_lst$model)
      res$roc 
    }
  }
})

#download as pdf
output$forest_rocPDF <- downloadHandler(
  filename = function(){"randomForest_roc.pdf"},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_rf){
        pdf(file, width=8, height=6)
        res<-evalm(vals$datasets[[currentSet()]]$rf_lst$model)
        res$roc 
        dev.off()
      }
    }
  }
)

output$forest_roc_cv <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_rf){
      ldf <- lift_df(vals$datasets[[currentSet()]]$rf_lst$model,vals$datasets[[currentSet()]]$rf_lst$model$levels[1])
      ggplot(ldf) +
        geom_line(aes(1 - Sp, Sn, color = fold)) +
        scale_color_discrete(guide = guide_legend(title = "Fold"))+
        xlab("x")+ylab("y") 
    }
  }
})

output$forest_top_features <- renderPlot({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_rf){
      plot(varImp(vals$datasets[[currentSet()]]$rf_lst$model),top=input$top_x_features) 
    }
  }
})

#download as pdf
output$forest_top_featuresPDF <- downloadHandler(
  filename = function(){"randomForest_features.pdf"},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_rf){
        pdf(file, width=8, height=6)
        plot(varImp(vals$datasets[[currentSet()]]$rf_lst$model)) 
        dev.off()
      }
    }
  }
)


output$forest_save_model <- downloadHandler(
  filename=function(){paste("random_forest_model.rds")},
  content = function(file){
    if(!is.null(currentSet())){
      if(vals$datasets[[currentSet()]]$has_rf){
        saveRDS(vals$datasets[[currentSet()]]$rf_lst$model,file) 
      }
    }
  }
)


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
  
  if(input$forest_tuning_type == "random"){
    shinyjs::show("forest_tuning_random",anim = F)
    shinyjs::hide("forest_tuning_regular")
  }else if (input$forest_tuning_type == "regular"){
    shinyjs::hide("forest_tuning_random")
    shinyjs::show("forest_tuning_regular",anim=F)
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
      updateSelectizeInput(session, "forest_covariable",choices=unique(meta[[input$forest_variable]]), server=T)
      
      #for features to include in model calculation; remove feature from list, which will be predicted 
      features<-group_columns[group_columns != input$forest_variable]
      updateSelectInput(session,"forest_features",choices = features) 
      
      # for metric to be chosen for best hyperparameters
      updateSelectInput(session, 'forest_select_best_metric', choices=input$forest_metrics)
    }
  }
})

#javascript show/hide toggle for advanced options
shinyjs::onclick("forest_toggle_advanced",shinyjs::toggle(id="forest_advanced",anim = T))
#show/hide exclude OTU-option if OTU abundances are to be used for model building
shinyjs::onclick("forest_otu",shinyjs::toggle(id="forest_exclude",anim = T))

