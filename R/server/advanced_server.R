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

    set.seed(input$forest_seed)
    message(paste0(Sys.time(), " - Starting random forest ..."))
    waiter_show(html = tagList(spin_rotating_plane(),"Preparing Data ..." ),color=overlay_color)
    meta <- as.data.frame(sample_data(vals$datasets[[currentSet()]]$phylo))
    otu_t <- as.data.frame(t(otu_table(vals$datasets[[currentSet()]]$phylo)))
    
    if(input$forest_clr){
      otu_t <- clr(otu_t)
    }
    
    combined_data <- buildForestDataset(meta, otu_t, input)
    class_labels <- as.factor(combined_data$variable)
    
    waiter_update(html = tagList(spin_rotating_plane(),"Partitioning Dataset & Resampling ..." ))
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
    waiter_update(html = tagList(spin_rotating_plane(),"Testing Model ..." ))
    predictions_model <- predict(model, newdata=testing)
    predictions_model_full <- predict(model, newdata=combined_data)
    con_matrix<-confusionMatrix(data=predictions_model, reference= class_labels[-inTraining])
    con_matrix_full<-confusionMatrix(data=predictions_model_full, reference= class_labels)
    
    message(paste0(Sys.time(), " - Starting random forest ..."))
    waiter_hide()
    showModal(forestReadyModal)
    
    return(list(cmtrx=con_matrix,cmtrx_full=con_matrix_full,model=model,class_labels=class_labels))
  }
})

forestReadyModal <- modalDialog(
  title = "Finished Random Forest!",
  "Scroll down to check the performance of your model!",
  easyClose = T, size="s"
)

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
    if (vals$datasets[[currentSet()]]$is_fastq || vals$datasets[[currentSet()]]$is_sample_data){
      shinyjs::hide("picrustFastaFile")
    }else{
      shinyjs::show("picrustFastaFile")
    }
  }
})

observeEvent(input$picrust2Start,{
  if(!is.null(currentSet())){
    message(paste0(Sys.time(), " - Starting picrust2 analysis ..."))
    
    waiter_show(html = tagList(spin_rotating_plane(),"Running Picrust2 ..." ),color=overlay_color)
    
    phylo <- vals$datasets[[currentSet()]]$phylo
    vals$datasets[[currentSet()]]$picrust_output <- NULL    # reset old picrust-output variable
    vals$datasets[[currentSet()]]$aldex_list <- NULL
    shinyjs::hide("download_picrust_div")
    
    foldername <- sprintf("/picrust2_%s", digest::digest(phylo))  # unique folder name for this output
    outdir <- paste0(tempdir(),foldername)
    if (dir.exists(outdir)){unlink(outdir, recursive = T)}    # remove output directory if it exists
    dir.create(outdir)
    
    biom_file = paste0(outdir,"/biom_picrust.biom")
    biom <- make_biom(data=otu_table(phylo))
    write_biom(biom, biom_file)
    message(paste0(Sys.time(), " - Wrote biom-file: ", biom_file))
    
    # use fasta file of dada2 pipeline if available
    if (vals$datasets[[currentSet()]]$is_fastq){
      fasta_file <- paste0(outdir,"/seqs.fasta")
      writeXStringSet(refseq(phylo), fasta_file)
      message(paste0(Sys.time(), " - Using dada2-generated fasta file:", fasta_file))
    }else{
      fasta_file <- input$picrustFastaFile$datapath
      message(paste0(Sys.time(), " - Using user-uploaded fasta file: ", fasta_file))
    }
    
    if(!vals$datasets[[currentSet()]]$is_sample_data){
      picrust_outdir <- paste0(outdir,"/picrust2out")       # this is the name of the final output directory of this picrust run
      command = paste0("/opt/anaconda3/bin/conda run -n picrust2 picrust2_pipeline.py --remove_intermediate -s ",fasta_file," -i ",biom_file, " -o ", picrust_outdir, " -p", ncores)
      message(paste0(Sys.time(), " - picrust2-command:"))
      message(command)
      #command = paste0("/home/alex/anaconda3/bin/conda run -n picrust2 picrust2_pipeline.py -s ",fasta_file," -i ",biom_file, " -o ", picrust_outdir, " -p", ncores)
      out <- system(command, wait = TRUE)
      shinyjs::show("download_picrust_div", anim = T)
      vals$datasets[[currentSet()]]$picrust_output <- picrust_outdir 
      message(paste0(Sys.time(), " - Finished picrust2 run; output in: ", picrust_outdir))
      outfiles <- list.files(picrust_outdir)
      message(outfiles)
      
      # generate tables
      p2_EC <- paste0(picrust_outdir, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
      p2_KO <- paste0(picrust_outdir, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
      p2_PW <- paste0(picrust_outdir, "/pathways_out/path_abun_unstrat.tsv.gz") 
    }else{
      p2_EC <- "testdata/picrust2/ec_pred_metagenome_unstrat.tsv.gz"
      p2_KO <- "testdata/picrust2/ko_pred_metagenome_unstrat.tsv.gz"
      p2_PW <- "testdata/picrust2/path_abun_unstrat.tsv.gz"
    }
    
    if(all(file.exists(c(p2_EC, p2_KO, p2_PW)))){
      vals$datasets[[currentSet()]]$has_picrust <- T
      
      message(paste0(Sys.time(), " - Starting differential analysis with ALDEx2 ... "))
      waiter_update(html = tagList(spin_rotating_plane(),"Differential analysis ..."))
      p2EC = as.data.frame(fread(p2_EC))
      rownames(p2EC) = p2EC$"function"
      p2EC = as.matrix(p2EC[,-1])
      p2EC = round(p2EC)
      
      p2KO = as.data.frame(fread(p2_KO))
      rownames(p2KO) = p2KO$"function"
      p2KO = as.matrix(p2KO[,-1])
      p2KO = round(p2KO)
      
      p2PW = as.data.frame(fread(p2_PW))
      rownames(p2PW) = p2PW$"pathway"
      p2PW = as.matrix(p2PW[,-1])
      p2PW = round(p2PW)
      
      # run ALDEx2 to perform differential abundance testing between 2(!) conditions 
      test <- "t"
      waiter_update(html = tagList(spin_rotating_plane(),"Differential analysis (EC)..."))
      aldex2_EC <- aldex(p2EC, sample_data(phylo)[[input$picrust_test_condition]], mc.samples = input$picrust_mc_samples, test=test, effect=T,denom="iqlr")
      message(paste0(Sys.time(), " - finished EC ... "))
      waiter_update(html = tagList(spin_rotating_plane(),"Differential analysis (KO)..."))
      aldex2_KO <- aldex(p2KO, sample_data(phylo)[[input$picrust_test_condition]], mc.samples = input$picrust_mc_samples, test=test, effect=T,denom="iqlr")
      message(paste0(Sys.time(), " - finished KO ... "))
      waiter_update(html = tagList(spin_rotating_plane(),"Differential analysis (PW)..."))
      aldex2_PW <- aldex(p2PW, sample_data(phylo)[[input$picrust_test_condition]], mc.samples = input$picrust_mc_samples, test=test, effect=T,denom="iqlr")
      message(paste0(Sys.time(), " - finished PW ... "))
      vals$datasets[[currentSet()]]$aldex_list <- list(aldex2_EC=aldex2_EC,
                                                       aldex2_KO=aldex2_KO,
                                                       aldex2_PW=aldex2_PW)
      message(paste0(Sys.time(), " - Finished differential analysis with ALDEx2 "))
      
    }else{
      message(paste0(Sys.time(), " - Did not find all files for differential analysis with ALDEx2; stopping ... "))
      message(paste0(c(p2_EC, p2_KO, p2_PW)))
      message(paste0(file.exists(c(p2_EC, p2_KO, p2_PW))))
      vals$datasets[[currentSet()]]$has_picrust <- F
      showModal(errorPicrustModal)
    }
    
    waiter_hide()
  }
})

errorPicrustModal <- modalDialog(
  title = "Error with picrust2",
  "Something went wrong with picrust2, not all files were created. Did your fastq-files have the correct OTU/ASV-names? 
  If you feel like you did nothing wrong, please contanct the author of namco.",
  easyClose = T, size="s"
)

#download results
output$download_picrust_raw <- downloadHandler(
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

output$picrust_download_ec <- downloadHandler(
  filename = function(){
    paste("picrust2_analysis_ec.tab")
  },
  content = function(file){
    if(!is.null(currentSet())){
      if (!is.null(vals$datasets[[currentSet()]]$aldex_list)){
        write.table(vals$datasets[[currentSet()]]$aldex_list$aldex2_EC, file=file, quote = F, sep="\t")
      }
    }
  }
)

output$picrust_download_ko <- downloadHandler(
  filename = function(){
    paste("picrust2_analysis_ko.tab")
  },
  content = function(file){
    if(!is.null(currentSet())){
      if (!is.null(vals$datasets[[currentSet()]]$aldex_list)){
        write.table(vals$datasets[[currentSet()]]$aldex_list$aldex2_KO, file=file, quote = F, sep="\t")
      }
    }
  }
)

output$picrust_download_pw <- downloadHandler(
  filename = function(){
    paste("picrust2_analysis_pw.tab")
  },
  content = function(file){
    if(!is.null(currentSet())){
      if (!is.null(vals$datasets[[currentSet()]]$aldex_list)){
        write.table(vals$datasets[[currentSet()]]$aldex_list$aldex2_PW, file=file, quote = F, sep="\t")
      }
    }
  }
)

# reactive data for long dataframes with sigificance column
aldex_reactive <- reactive({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_picrust){
    aldex2_EC <- vals$datasets[[currentSet()]]$aldex_list$aldex2_EC
    aldex2_KO <- vals$datasets[[currentSet()]]$aldex_list$aldex2_KO
    aldex2_PW <- vals$datasets[[currentSet()]]$aldex_list$aldex2_PW
    
    list<-lapply(list(aldex2_EC,aldex2_KO,aldex2_PW), function(x){
      x[["func"]] <- rownames(x)
      a.long <- gather(x[,c(12,4,6,8,9)], pvalue_type, pvalue, 4:5)
      colnames(a.long) <- c("func","difference","effect","pvalue_type","pvalue")
      a.long$significant <- ifelse(a.long$pvalue<input$picrust_signif_lvl,T,F) 
      return(a.long)
    })
    
    out <- list(EC_long=list[[1]], KO_long=list[[2]], PW_long=list[[3]])
    return(out)
    }
  }
})

output$picrust_ec_effect_signif <- renderPrint({
  if(!is.null(aldex_reactive())){
    aldex_reactive()$EC_long[aldex_reactive()$EC_long[["significant"]]==T,][["func"]]
  }
})
  

output$picrust_ko_effect_signif <- renderPrint({
  if(!is.null(aldex_reactive())){
    aldex_reactive()$KO_long[aldex_reactive()$KO_long[["significant"]]==T,][["func"]]
  }
})


output$picrust_pw_effect_signif <- renderPrint({
  if(!is.null(aldex_reactive())){
    aldex_reactive()$PW_long[aldex_reactive()$PW_long[["significant"]]==T,][["func"]]
  }
})

output$picrust_ec_effect_signif_value <- renderValueBox({
  if(!is.null(aldex_reactive())){
    val <- length(aldex_reactive()$EC_long[aldex_reactive()$EC_long[["significant"]]==T,][["func"]])
    valueBox(val, "Sinificant ECs",icon = icon("arrow-up"), color="olive")
  }
})

output$picrust_ko_effect_signif_value <- renderValueBox({
  if(!is.null(aldex_reactive())){
    val <- length(aldex_reactive()$KO_long[aldex_reactive()$KO_long[["significant"]]==T,][["func"]])
    valueBox(val, "Sinificant KOs",icon = icon("arrow-up"), color="olive")
  }
})

output$picrust_pw_effect_signif_value <- renderValueBox({
  if(!is.null(aldex_reactive())){
    val <- length(aldex_reactive()$PW_long[aldex_reactive()$PW_long[["significant"]]==T,][["func"]])
    valueBox(val, "Sinificant PWs",icon = icon("arrow-up"), color="olive")
  }
})

#### pvalue plots ####
output$picrust_ec_effect_plot <- renderPlot({
  if(!is.null(aldex_reactive())){
    a.long <- aldex_reactive()$EC_long
    p<-ggplot(data=a.long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Effect size vs P-value")+
      scale_color_manual(labels=c("P-value","BH-adjusted"), values=c("#0072B2", "#D55E00"))+
      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
      theme_minimal()
    
    if(input$picrust_signif_label){
      p<-p+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    p
  }
})

output$picrust_ec_vulcano_plot <- renderPlot({
  if(!is.null(aldex_reactive())){
    a.long <- aldex_reactive()$EC_long
    p<-ggplot(data=a.long,aes(x=difference,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Difference vs P-value")+
      scale_color_manual(labels=c("P-value","BH-adjusted"), values=c("#0072B2", "#D55E00"))+
      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
      theme_minimal()
    
    if(input$picrust_signif_label){
      p<-p+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    p
  }
})

output$picrust_ko_effect_plot <- renderPlot({
  if(!is.null(aldex_reactive())){
    a.long <- aldex_reactive()$KO_long
    p<-ggplot(data=a.long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Effect size vs P-value")+
      scale_color_manual(labels=c("P-value","BH-adjusted"), values=c("#0072B2", "#D55E00"))+
      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
      theme_minimal()
    
    if(input$picrust_signif_label){
      p<-p+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    p
  }
})

output$picrust_ko_vulcano_plot <- renderPlot({
  if(!is.null(aldex_reactive())){
    a.long <- aldex_reactive()$KO_long
    p<-ggplot(data=a.long,aes(x=difference,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Difference vs P-value")+
      scale_color_manual(labels=c("P-value","BH-adjusted"), values=c("#0072B2", "#D55E00"))+
      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
      theme_minimal()
    
    if(input$picrust_signif_label){
      p<-p+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    p
  }
})

output$picrust_pw_effect_plot <- renderPlot({
  if(!is.null(aldex_reactive())){
    a.long <- aldex_reactive()$PW_long
    p<-ggplot(data=a.long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Effect size vs P-value")+
      scale_color_manual(labels=c("P-value","BH-adjusted"), values=c("#0072B2", "#D55E00"))+
      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
      theme_minimal()
    
    if(input$picrust_signif_label){
      p<-p+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    p
  }
})

output$picrust_pw_vulcano_plot <- renderPlot({
  if(!is.null(aldex_reactive())){
    a.long <- aldex_reactive()$PW_long
    p<-ggplot(data=a.long,aes(x=difference,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Difference vs P-value")+
      scale_color_manual(labels=c("P-value","BH-adjusted"), values=c("#0072B2", "#D55E00"))+
      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
      theme_minimal()
    
    if(input$picrust_signif_label){
      p<-p+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    p
    
  }
})


