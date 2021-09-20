####heatmap#### 

# plot heatmap of OTU abundances per sample
abundanceHeatmapReact <- reactive({
  if(!is.null(currentSet())){
    set.seed(seed)
    phylo <- vals$datasets[[currentSet()]]$phylo
    phylo <- transform_sample_counts(phylo, function(x) x+1)  # pseudocount to not get -Inf values
    #phylo <- transform_sample_counts(phylo, function(x) x/sum(x)*100) #relative abunance
    #check for unifrac distance --> (needs phylo tree file):
    l <- list()
    if(!is.null(vals$datasets[[currentSet()]]$unifrac_dist)){
      #save generalized unifrac distance as global variable to use it for heatmap
      gunifrac_heatmap <<- as.dist(vals$datasets[[currentSet()]]$unifrac_dist)
      hm_distance <- if(input$heatmapDistance == "gunifrac") "gunifrac_heatmap" else input$heatmapDistance
      if(input$heatmapOrderSamples){
        p<-plot_heatmap(phylo,method = input$heatmapOrdination, distance = hm_distance, sample.label = input$heatmapSample, sample.order = input$heatmapSample)
      }else{
        p<-plot_heatmap(phylo,method = input$heatmapOrdination, distance = hm_distance, sample.label = input$heatmapSample)
      }
      l <- list(gg=p)
    }
  }
})

output$abundanceHeatmap <- renderPlotly({
  if(!is.null(abundanceHeatmapReact())){
    ggplotly(abundanceHeatmapReact()$gg)
  }
})

#download as pdf
output$abundanceHeatmapPDF <- downloadHandler(
  filename = function(){"abundance_heatmap.pdf"},
  content = function(file){
    if(!is.null(abundanceHeatmapReact())){
      ggsave(file, abundanceHeatmapReact()$gg, device="pdf", width = 10, height = 7)
    }
  }
)

####associations####

observeEvent(input$associations_start,{
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_meta){
      message(paste0(Sys.time(), " - building SIAMCAT object ..."))
      waiter_show(html = tagList(spin_rotating_plane(),"Calculating differential associations ..."),color=overlay_color)
      
      phylo <- vals$datasets[[currentSet()]]$phylo
      
      if(input$associations_level != "OTU"){
        phylo_glom <- glom_taxa_custom(phylo, input$associations_level) #merge OTUs with same taxonomic level
        phylo <- phylo_glom$phylo_rank
        taxa_names(phylo) <- phylo_glom$taxtab[[input$associations_level]]
        rel_otu <- relAbundanceTo1(data.frame(otu_table(phylo), check.names=F))
      }else{
        rel_otu <- relAbundanceTo1(data.frame(otu_table(phylo), check.names=F))
      }
      
      meta <- data.frame(sample_data(phylo))
      meta <- data.frame(t(na.omit(t(meta))))
      
      tryCatch({
        # siamcat works only with 5 or more samples (https://git.embl.de/grp-zeller/SIAMCAT/-/blob/a1c662f343e99dabad4de024d4c993deba91bb0c/R/validate_data.r#L82)
        if(sum(meta[[input$associations_label]]==input$associations_case) <= 5){stop(siamcatNotEnoughSamplesError, call.=F)}
        
        s.obj <- siamcat(feat=rel_otu, meta=meta, label= input$associations_label, case= input$associations_case)
        s.obj <- filter.features(s.obj)
        vals$datasets[[currentSet()]]$siamcat <- s.obj
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

output$associationsPDF <- downloadHandler(
  filename=function(){"associations.pdf"},
  content = function(file){
    if(!is.null(vals$datasets[[currentSet()]]$siamcat_plot)){
      s.obj <- vals$datasets[[currentSet()]]$siamcat
      sort.by <- c("p.val","fc","pr.shift")[which(input$associations_sort==c("p-value","fold-change","prevalence shift"))]
      panels <- c("fc","auroc","prevalence")[which(input$associations_panels==c("fold-change","AU-ROC","prevalence"))]
      print("hello")
      suppressMessages(check.associations(s.obj, fn.plot = file, prompt=F, verbose=0,
                                          alpha = input$associations_alpha, 
                                          max.show = input$assiciation_show_numer, 
                                          sort.by = sort.by,
                                          panels = panels))
    }
  }
)

####correlation####
corrReactive <- reactive({
  if(!is.null(currentSet())){
    
    waiter_show(html = tagList(spin_rotating_plane(),"Calculating pairwise correlations ..."),color=overlay_color)
    
    phylo <- vals$datasets[[currentSet()]]$phylo 
    otu <- data.frame(phylo@otu_table, check.names = F)
    meta <- data.frame(phylo@sam_data, check.names = F)
    meta_numeric_all <- meta[,unlist(lapply(meta, is.numeric))]
    meta_numeric <- as.data.frame(meta_numeric_all[,input$corrSelectGroups])
    colnames(meta_numeric) <- input$corrSelectGroups
    
    tryCatch({
      if(length(meta_numeric) > 0){
        # fill NA entries 
        meta_fixed = as.data.frame(apply(meta_numeric, 2, fill_NA_INF.mean))
        
        # remove columns with only a single value
        meta_fixed[names(which(apply(meta_fixed,2,function(x){length(unique(x))})==1))] <- NULL
        
        complete_data <- cbind(meta_fixed, t(otu))
        meta_names <- colnames(meta_fixed)
      }else if(input$corrIncludeTax){
        complete_data <- cbind(t(otu))
        meta_names <- NULL
      }else{
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
                                       input$corrSignifCutoff)
      waiter_hide()
      return(list(my_cor_matrix=corr_subset$my_cor_matrix,
                  my_pvl_matrix=corr_subset$my_pvl_matrix,
                  my_pairs=corr_subset$my_pairs,
                  diagonale=corr_subset$diagonale))
    }, error = function(e){
      waiter_hide()
      print(e$message)
      showModal(errorModal(e$message))
    })
  }
})

output$corrPlot <- renderPlot({
  if(!is.null(corrReactive())){
    plot_correlation_custom(corrReactive()$my_cor_matrix, corrReactive()$my_pvl_matrix, input)
  }
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


####random forest models####

#calculate confusion matrix using random forest aproach
observeEvent(input$forest_start,{
  if(!is.null(currentSet())){

    set.seed(seed)
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
    predictions_model <- predict(model, newdata=testing)
    predictions_model_full <- predict(model, newdata=combined_data)
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
    
    tryCatch({
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
        #TODO: check if all OTUs have a sequence in fasta file
        fasta_file <- input$picrustFastaFile$datapath
        if(vals$datasets[[currentSet()]]$is_sample_data) fasta_file <- "testdata/seqs.fasta"
        fasta <- readDNAStringSet(fasta_file, format="fasta")
        all_taxa <- taxa_names(phylo)
        fasta <- fasta[fasta@ranges@NAMES %in% all_taxa]
        writeXStringSet(fasta, fasta_file)
        message(paste0(Sys.time(), " - Using user-uploaded fasta file: ", fasta_file))
      }
      
      if(!vals$datasets[[currentSet()]]$is_sample_data){
        picrust_outdir <- paste0(outdir,"/picrust2out")       # this is the name of the final output directory of this picrust run
        command = paste0("/opt/anaconda3/bin/conda run -n picrust2 picrust2_pipeline.py --remove_intermediate -s ",fasta_file," -i ",biom_file, " -o ", picrust_outdir, " -p", ncores*2)
        #command = paste0("/home/alex/anaconda3/bin/conda run -n picrust2 picrust2_pipeline.py --remove_intermediate -s ",fasta_file," -i ",biom_file, " -o ", picrust_outdir, " -p", ncores)
        message(paste0(Sys.time(), " - picrust2-command:"))
        message(command)
        # here picrust2 is started:
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
        marker_nsti <- paste0(picrust_outdir, "/marker_predicted_and_nsti.tsv.gz")
      }else{
        p2_EC <- "testdata/picrust2/ec_pred_metagenome_unstrat.tsv.gz"
        p2_KO <- "testdata/picrust2/ko_pred_metagenome_unstrat.tsv.gz"
        p2_PW <- "testdata/picrust2/path_abun_unstrat.tsv.gz"
        marker_nsti <- "testdata/picrust2/marker_predicted_and_nsti.tsv.gz"
      }
      
      if(all(file.exists(c(p2_EC, p2_KO, p2_PW, marker_nsti)))){
        vals$datasets[[currentSet()]]$has_picrust <- T
        
        # normalize OTU table by copy-number
        if(input$picrust_copy_number_normalization){
          message(paste0(Sys.time(), " - Normalizing OTU-table by copy-numbers ..."))
          otu <- vals$datasets[[currentSet()]]$rawData
          copy_number_df <- as.data.frame(fread(marker_nsti))
          copy_numbers <- copy_number_df[order(as.character(copy_number_df$sequence)),][["16S_rRNA_Count"]]
          normalized_dat <- otu[order(as.character(rownames(otu))),] / copy_numbers
          vals$datasets[[currentSet()]]$normalized_dat$norm_tab <- normalized_dat
          #update phylo-object
          otu_table(vals$datasets[[currentSet()]]$phylo) <- otu_table(normalized_dat,T)
        }
        
        # do not run differential analysis if selected group does not exist
        if(is.null(sample_data(phylo)[[input$picrust_test_condition]])){
          message(paste0(Sys.time(), " - Skipping differential analysis; sample group not found ... "))
          stop(picrustDifferentialGroupNotFoundError, call. = F)
        }
        
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
        label <- sample_data(vals$datasets[[currentSet()]]$phylo)[[input$picrust_test_condition]]
        names(label) <- sample_names(vals$datasets[[currentSet()]]$phylo) 
        vals$datasets[[currentSet()]]$aldex_list <- list(aldex2_EC=aldex2_EC,
                                                         aldex2_KO=aldex2_KO,
                                                         aldex2_PW=aldex2_PW,
                                                         p2EC = p2EC,
                                                         p2KO = p2KO,
                                                         p2PW = p2PW,
                                                         label = label)
        message(paste0(Sys.time(), " - Finished differential analysis with ALDEx2 "))
        waiter_hide()
        
      }else{
        message(paste0(Sys.time(), " - Did not find all files for differential analysis with ALDEx2; stopping ... "))
        message(paste0(c(p2_EC, p2_KO, p2_PW, marker_nsti)))
        message(paste0(file.exists(c(p2_EC, p2_KO, p2_PW, marker_nsti))))
        vals$datasets[[currentSet()]]$has_picrust <- F
        stop(picrustFilesMissingError, call. = F)
      }
      
    }, error=function(e){
      waiter_hide()
      print(e$message)
      showModal(errorModal(e$message))
    })
  }
})

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

# reactive data for long dataframes with significance column
aldex_reactive <- reactive({
  if(!is.null(currentSet())){
    if(!is.null(vals$datasets[[currentSet()]]$aldex_list)){
      message(Sys.time(), " - generating picrust analysis tables ...")
      abundances <- list(vals$datasets[[currentSet()]]$aldex_list$p2EC,
                         vals$datasets[[currentSet()]]$aldex_list$p2KO,
                         vals$datasets[[currentSet()]]$aldex_list$p2PW)
      aldex_results <- list(vals$datasets[[currentSet()]]$aldex_list$aldex2_EC,
                            vals$datasets[[currentSet()]]$aldex_list$aldex2_KO,
                            vals$datasets[[currentSet()]]$aldex_list$aldex2_PW)
 
      nsamples <- nsamples(vals$datasets[[currentSet()]]$phylo)
      label <- vals$datasets[[currentSet()]]$aldex_list$label
      if(input$picrustTest=="Welch's t-test"){
        pval_test <- "we.eBH"
        aldex_columns <- c("significant","func","diff.btw","effect","we.ep","we.eBH")
      }else{
        pval_test <- "wi.eBH"
        aldex_columns <- c("significant","func","diff.btw","effect","wi.ep","wi.eBH")
      }
      out_list<-lapply(c(1,2,3), function(x){
        aldex.tab <- aldex_results[[x]]
        abundances.tab <- abundances[[x]]
        aldex.tab[["func"]] <- rownames(aldex.tab)
        aldex.tab$significant <- ifelse((aldex.tab[[pval_test]]<input$picrust_signif_lvl & aldex.tab[["effect"]]>input$picrust_signif_lvl_effect),T,F)
        a.long <- gather(aldex.tab[,aldex_columns], pvalue_type, pvalue, 5:6)
        colnames(a.long) <- c("significant","func","difference","effect","pvalue_type","pvalue")
        
        x.merged <- merge(abundances.tab, aldex.tab, by=0)
        rownames(x.merged)<-x.merged$Row.names
        x.merged$func <- rownames(x.merged)
        x.merged[["Row.names"]] <- NULL
        x.merged.long <- gather(x.merged,SampleID,abundance,1:nsamples)
        x.merged.long$label<-unlist(lapply(x.merged.long$SampleID, function(y){return(label[[y]])}))
        colnames(x.merged.long)[colnames(x.merged.long) == pval_test] <- "p_value"
        signif_df <- data.frame(x.merged.long[x.merged.long$significant==T,])
        signif_df[["p_value"]] <- -log10(signif_df[["p_value"]])
        signif_df <- signif_df[order(signif_df[["p_value"]], decreasing = F),]
        
        return(list(a.long, signif_df))
      })
      
      out <- list(EC_long=out_list[[1]][[1]], KO_long=out_list[[2]][[1]], PW_long=out_list[[3]][[1]],
                  EC_signif=out_list[[1]][[2]], KO_signif=out_list[[2]][[2]], PW_signif=out_list[[3]][[2]])
      return(out)
    }
  }
})


####picrust2 plots ####

picrust_plots_reactive <- reactive({
  if(!is.null(aldex_reactive())){
    #### pvalue plots ####
    EC_long <- aldex_reactive()$EC_long
    ec_effect_plot<-ggplot(data=EC_long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
                    geom_point(alpha=0.8)+
                    ggtitle("Effect size vs P-value")+
                    scale_color_manual(labels=c("BH-adjusted","P-value"), values=c("#0072B2", "#E18E4C"))+
                    geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
                    geom_vline(xintercept = input$picrust_signif_lvl_effect, color="black", linetype="dashed", alpha=0.8)+
                    theme_minimal()
    if(input$picrust_signif_label){
      ec_effect_plot<-ec_effect_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    
    ec_vulcano_plot<-ggplot(data=EC_long,aes(x=difference,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
                      geom_point(alpha=0.8)+
                      ggtitle("Difference vs P-value")+
                      scale_color_manual(labels=c("BH-adjusted","P-value"), values=c("#0072B2", "#E18E4C"))+
                      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
                      theme_minimal()
    if(input$picrust_signif_label){
      ec_vulcano_plot<-ec_vulcano_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    
    
    KO_long <- aldex_reactive()$KO_long
    ko_effect_plot<-ggplot(data=KO_long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
                    geom_point(alpha=0.8)+
                    ggtitle("Effect size vs P-value")+
                    scale_color_manual(labels=c("BH-adjusted","P-value"), values=c("#0072B2", "#E18E4C"))+
                    geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
                    geom_vline(xintercept = input$picrust_signif_lvl_effect, color="black", linetype="dashed", alpha=0.8)+
                    theme_minimal()
    if(input$picrust_signif_label){
      ko_effect_plot<-ko_effect_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    
    ko_vulcano_plot<-ggplot(data=KO_long,aes(x=difference,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
                            geom_point(alpha=0.8)+
                            ggtitle("Difference vs P-value")+
                            scale_color_manual(labels=c("BH-adjusted","P-value"), values=c("#0072B2", "#E18E4C"))+
                            geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
                            theme_minimal()
    if(input$picrust_signif_label){
      ko_vulcano_plot<-ko_vulcano_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    
    
    PW_long <- aldex_reactive()$PW_long
    pw_effect_plot<-ggplot(data=PW_long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
                    geom_point(alpha=0.8)+
                    ggtitle("Effect size vs P-value")+
                    scale_color_manual(labels=c("BH-adjusted","P-value"), values=c("#0072B2", "#E18E4C"))+
                    geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
                    geom_vline(xintercept = input$picrust_signif_lvl_effect, color="black", linetype="dashed", alpha=0.8)+
                    theme_minimal()
    if(input$picrust_signif_label){
      pw_effect_plot<-pw_effect_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    
    pw_vulcano_plot <-ggplot(data=PW_long,aes(x=difference,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
                      geom_point(alpha=0.8)+
                      ggtitle("Difference vs P-value")+
                      scale_color_manual(labels=c("BH-adjusted","P-value"), values=c("#0072B2", "#E18E4C"))+
                      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
                      theme_minimal()
    if(input$picrust_signif_label){
      pw_vulcano_plot<-pw_vulcano_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    
    #### signficance boxplots ####
    EC_signif <- aldex_reactive()$EC_signif
    nsamples <- nsamples(vals$datasets[[currentSet()]]$phylo)
    if(dim(EC_signif)[1] > 0){
      if(dim(EC_signif)[1] > input$picrust_ec_signif_plot_show){
        rows_to_keep <- input$picrust_ec_signif_plot_show*nsamples
        EC_signif <- na.omit(EC_signif[1:rows_to_keep,])
      }
      
      p1 <- ggplot(data=EC_signif, aes(x=abundance,y=func,fill=as.factor(label)))+
        geom_boxplot()+
        theme_minimal()+
        theme(legend.position = c(0.9,0.9),
              legend.title = element_blank())+
        scale_fill_brewer(palette = "Set1")+
        ylab("Function")
      
      p2 <- ggplot(data=EC_signif,aes(x=p_value/nsamples,y=func))+
        geom_bar(stat="identity", width=0.35, aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none")+
        xlab("-log10(P-value)")+
        geom_vline(xintercept = -log10(input$picrust_signif_lvl), color="red", linetype="dashed")
      
      p3 <- ggplot(data=EC_signif,aes(x=effect/nsamples,y=func))+
        geom_bar(stat="identity", width=0.35,aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(), 
              legend.position = "none")+
        xlab("Effect size")+
        geom_vline(xintercept = input$picrust_signif_lvl_effect, color="red", linetype="dashed")
      
      ec_signif_plot_list <- list(p1=p1,p2=p2,p3=p3)#grid.arrange(p1,p2,p3,ncol=3,widths=c(3,1,1))
    }else{
      ec_signif_plot_list <- NULL
    }
    
    KO_signif <- aldex_reactive()$KO_signif
    if(dim(KO_signif)[1] > 0){
      if(dim(KO_signif)[1] > input$picrust_ko_signif_plot_show){
        rows_to_keep <- input$picrust_ko_signif_plot_show*nsamples
        KO_signif <- na.omit(KO_signif[1:rows_to_keep,])
      }    
      p1 <- ggplot(data=KO_signif, aes(x=abundance,y=func,fill=as.factor(label)))+
        geom_boxplot()+
        theme_minimal()+
        theme(legend.position = c(0.9,0.9),
              legend.title = element_blank())+
        scale_fill_brewer(palette = "Set1")+
        ylab("Function")
      
      p2 <- ggplot(data=KO_signif,aes(x=p_value/nsamples,y=func))+
        geom_bar(stat="identity", width=0.35, aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none")+
        xlab("-log10(P-value)")+
        geom_vline(xintercept = -log10(input$picrust_signif_lvl), color="red", linetype="dashed")
      
      p3 <- ggplot(data=KO_signif,aes(x=effect/nsamples,y=func))+
        geom_bar(stat="identity", width=0.35,aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(), 
              legend.position = "none")+
        xlab("Effect size")+
        geom_vline(xintercept = input$picrust_signif_lvl_effect, color="red", linetype="dashed")
      
      ko_signif_plot_list <- list(p1=p1,p2=p2,p3=p3)#grid.arrange(p1,p2,p3,ncol=3,widths=c(3,1,1))
    }else{
      ko_signif_plot_list <- NULL
    }
    
    PW_signif <- aldex_reactive()$PW_signif
    if(dim(PW_signif)[1] > 0){
      if(dim(PW_signif)[1] > input$picrust_pw_signif_plot_show){
        rows_to_keep <- input$picrust_pw_signif_plot_show*nsamples
        PW_signif <- na.omit(PW_signif[1:rows_to_keep,])
      }     
      p1 <- ggplot(data=PW_signif, aes(x=abundance,y=func,fill=as.factor(label)))+
        geom_boxplot()+
        theme_minimal()+
        theme(legend.position = c(0.9,0.9),
              legend.title = element_blank())+
        scale_fill_brewer(palette = "Set1")+
        ylab("Function")
      
      p2 <- ggplot(data=PW_signif,aes(x=p_value/nsamples,y=func))+
        geom_bar(stat="identity", width=0.35, aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none")+
        xlab("-log10(P-value)")+
        geom_vline(xintercept = -log10(input$picrust_signif_lvl), color="red", linetype="dashed")
      
      p3 <- ggplot(data=PW_signif,aes(x=effect/nsamples,y=func))+
        geom_bar(stat="identity", width=0.35,aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(), 
              legend.position = "none")+
        xlab("Effect size")+
        geom_vline(xintercept = input$picrust_signif_lvl_effect, color="red", linetype="dashed")
      
      pw_signif_plot_list <- list(p1=p1,p2=p2,p3=p3)#grid.arrange(p1,p2,p3,ncol=3,widths=c(3,1,1))
    }else{
      pw_signif_plot_list <- NULL
    }
    
    plots <- list(ec_effect_plot=ec_effect_plot,
                  ec_vulcano_plot=ec_vulcano_plot,
                  ko_effect_plot=ko_effect_plot,
                  ko_vulcano_plot=ko_vulcano_plot,
                  pw_effect_plot=pw_effect_plot,
                  pw_vulcano_plot=pw_vulcano_plot,
                  ec_signif_plot_list=ec_signif_plot_list,
                  ko_signif_plot_list=ko_signif_plot_list,
                  pw_signif_plot_list=pw_signif_plot_list)
  }
})


output$picrust_ec_effect_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    picrust_plots_reactive()$ec_effect_plot
  }
})

output$picrust_ec_effectPDF <- downloadHandler(
  filename = function(){"EC_effect.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      ggsave(file, picrust_plots_reactive()$ec_effect_plot, device="pdf", width = 10, height = 7)
    }
  }
)

output$picrust_ec_vulcano_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    picrust_plots_reactive()$ec_vulcano_plot
  }
})

output$picrust_ec_vulcanoPDF <- downloadHandler(
  filename = function(){"EC_vulcano.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      ggsave(file, picrust_plots_reactive()$ec_vulcano_plot, device="pdf", width = 10, height = 7)
    }
  }
)

output$picrust_ko_effect_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    picrust_plots_reactive()$ko_effect_plot
  }
})

output$picrust_ko_effectPDF <- downloadHandler(
  filename = function(){"KO_effect.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      ggsave(file, picrust_plots_reactive()$ko_effect_plot, device="pdf", width = 10, height = 7)
    }
  }
)

output$picrust_ko_vulcano_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    picrust_plots_reactive()$ko_vulcano_plot
  }
})

output$picrust_ko_vulcanoPDF <- downloadHandler(
  filename = function(){"KO_vulcano.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      ggsave(file, picrust_plots_reactive()$ko_vulcano_plot, device="pdf", width = 10, height = 7)
    }
  }
)

output$picrust_pw_effect_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    picrust_plots_reactive()$pw_effect_plot
  }
})

output$picrust_pw_effectPDF <- downloadHandler(
  filename = function(){"PW_effect.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      ggsave(file, picrust_plots_reactive()$pw_effect_plot, device="pdf", width = 10, height = 7)
    }
  }
)

output$picrust_pw_vulcano_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    picrust_plots_reactive()$pw_vulcano_plot
  }
})

output$picrust_pw_vulcanoPDF <- downloadHandler(
  filename = function(){"PW_vulcano.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      ggsave(file, picrust_plots_reactive()$pw_vulcano_plot, device="pdf", width = 10, height = 7)
    }
  }
)

output$picrust_ec_signif_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    if(!is.null(picrust_plots_reactive()$ec_signif_plot_list)){
      l <- picrust_plots_reactive()$ec_signif_plot_list
      grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1))
    }
  }
})

output$picrust_ec_signifPDF <- downloadHandler(
  filename = function(){"EC_significant_boxplots.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      l <- picrust_plots_reactive()$ec_signif_plot_list
      ggsave(file, grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1)), device="pdf", width = 12, height = 8)
    }
  }
)

output$picrust_ko_signif_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    if(!is.null(picrust_plots_reactive()$ko_signif_plot_list)){
      l <- picrust_plots_reactive()$ko_signif_plot_list
      grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1))
    }
  }
})

output$picrust_ko_signifPDF <- downloadHandler(
  filename = function(){"KO_significant_boxplots.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      l <- picrust_plots_reactive()$ko_signif_plot_list
      ggsave(file, grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1)), device="pdf", width = 12, height = 8)
    }
  }
)

output$picrust_pw_signif_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    if(!is.null(picrust_plots_reactive()$pw_signif_plot_list)){
      l <- picrust_plots_reactive()$pw_signif_plot_list
      grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1))
    }
  }
})

output$picrust_pw_signifPDF <- downloadHandler(
  filename = function(){"PW_significant_boxplots.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      l <- picrust_plots_reactive()$pw_signif_plot_list
      ggsave(file, grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1)), device="pdf", width = 12, height = 8)
    }
  }
)

output$picrust_ec_effect_signif <- renderPrint({
  if(!is.null(aldex_reactive())){
    unique(aldex_reactive()$EC_long[aldex_reactive()$EC_long[["significant"]]==T,][["func"]])
  }
})


output$picrust_ko_effect_signif <- renderPrint({
  if(!is.null(aldex_reactive())){
    unique(aldex_reactive()$KO_long[aldex_reactive()$KO_long[["significant"]]==T,][["func"]])
  }
})


output$picrust_pw_effect_signif <- renderPrint({
  if(!is.null(aldex_reactive())){
    unique(aldex_reactive()$PW_long[aldex_reactive()$PW_long[["significant"]]==T,][["func"]])
  }
})

output$picrust_ec_effect_signif_value <- renderValueBox({
  val <- 0
  if(!is.null(aldex_reactive())){
    val <- length(unique(aldex_reactive()$EC_long[aldex_reactive()$EC_long[["significant"]]==T,][["func"]]))
  }
  valueBox(val, "Sinificant ECs",icon = icon("arrow-up"), color="olive")
})

output$picrust_ko_effect_signif_value <- renderValueBox({
  val <- 0
  if(!is.null(aldex_reactive())){
    val <- length(unique(aldex_reactive()$KO_long[aldex_reactive()$KO_long[["significant"]]==T,][["func"]]))
  }
  valueBox(val, "Sinificant KOs",icon = icon("arrow-up"), color="olive")
})

output$picrust_pw_effect_signif_value <- renderValueBox({
  val <- 0
  if(!is.null(aldex_reactive())){
    val <- length(unique(aldex_reactive()$PW_long[aldex_reactive()$PW_long[["significant"]]==T,][["func"]]))
  }
  valueBox(val, "Sinificant PWs",icon = icon("arrow-up"), color="olive")
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
          cluster_meta[["sample_group"]] <- as.character(clusters)
          
          # build new phyloseq object with new meta
          phylo <- merge_phyloseq(phylo, sample_data(cluster_meta))
        }else{
          cluster_variables <- NULL
          cluster_meta <- NULL
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
        }else{
          alphaTab <- vals$datasets[[currentSet()]]$alpha_diversity
          plot_df <- alphaTab
        }
        
        # calculate significant features 
        
        features_df <- suppressWarnings(over_time_serial_comparison(phylo, input$timeSeriesGroup, ifelse(input$timeSeriesClusterK == 0, input$timeSeriesBackground, "sample_group")))
        waiter_hide()
        showModal(infoModal("Finished time-series analysis. Select one or more taxa to display the plot!"))
        return(list(plot_df=plot_df,
                    cluster_variables=cluster_variables,
                    cluster_meta=cluster_meta,
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
  if(!is.null(timeSeriesReactive())){
    df <- timeSeriesReactive()$features_df
    df$name <- factor(df$name, levels=df$name[order(df$pvalue)])
    ggplot(df, aes(x=pvalue, y=name))+
      geom_col(width = .75,na.rm = T)+
      ggtitle("Results of Friedman test for each taxonomic level & numeric meta-variable over the selected time-points and blocks.")+
      ylab("Name of taxa or meta-variable")+
      xlab("p-value")
  }
}, height=800)

timeSeriesPlotReactive <- reactive({
  if(!is.null(timeSeriesReactive())){
    plot_df <- as.data.table(timeSeriesReactive()$plot_df)
    colnames(plot_df)[which(colnames(plot_df)==input$timeSeriesGroup)] <- "reference" 
    if(input$timeSeriesClusterK == 0) colnames(plot_df)[which(colnames(plot_df)==input$timeSeriesBackground)] <- "sample_group" 
    if(!input$timeSeriesMeanLine %in% c("NONE","")) colnames(plot_df)[which(colnames(plot_df)==input$timeSeriesMeanLine)] <- "time_series_mean"
    p <- NULL
    title_text <- NULL
    
    # no need to select a taxa if diversity measure is chosen instead of abundance
    # --> no facet_wrap
    if (input$timeSeriesMeasure %in% c("Abundance", "relative Abundance")) {
      colnames(plot_df)[which(colnames(plot_df)=="Abundance")] <- "measure" 
      plot_df <- plot_df[plot_df[["OTU"]] %in% input$timeSeriesTaxaSelect,]
      
      if(!is.null(input$timeSeriesTaxaSelect)){
        p<-ggplot(plot_df, aes(x=reference, y=measure))+
          geom_line(aes(group=sample_group),alpha=0.5, color="grey")+
          facet_wrap(~OTU, scales="free")+
          xlab(input$timeSeriesGroup)+
          ylab(input$timeSeriesMeasure)+
          labs(color=ifelse(input$timeSeriesClusterK > 0,"Cluster ID",input$timeSeriesMeanLine))+
          theme_bw()+
          ggtitle(paste0("Time-series analysis at ",input$timeSeriesTaxa," level. \n", 
                         input$timeSeriesBackground, " is displayed as small grey lines in the back; \n",
                         "For ",input$timeSeriesMeanLine, " the mean ", input$timeSeriesMeasure, " over the time-points is displayed."))
      }
    }else{
      colnames(plot_df)[which(colnames(plot_df)==input$timeSeriesMeasure)] <- "measure"
      p<-ggplot(plot_df, aes(x=reference, y=measure))+
        geom_line(aes(group=sample_group),alpha=0.5, color="grey")+
        xlab(input$timeSeriesGroup)+
        ylab(input$timeSeriesMeasure)+
        labs(color=ifelse(input$timeSeriesClusterK > 0,"Cluster ID",input$timeSeriesMeanLine))+
        theme_bw()+
        ggtitle(paste0("Time-series analysis \n", 
                       input$timeSeriesBackground, " is displayed as small grey lines in the back; \n",
                       "For ",input$timeSeriesMeanLine, " the mean ", input$timeSeriesMeasure, " over the time-points is displayed."))
    }
    if(!is.null(p)){
      if(input$timeSeriesClusterK > 0){
        p <- p + stat_summary(fun=mean, geom="line", size=input$timeSeriesLineSize, aes(group=sample_group, color=sample_group))
      }
      if(input$timeSeriesMeanLine != "NONE"){
        p <- p + stat_summary(fun=mean, geom="line", size=input$timeSeriesLineSize, aes(group=time_series_mean, color=as.character(time_series_mean)))
      }
      return(list(plot=p))  
    }
  }
})

output$timeSeriesPlot <- renderPlot({
  if(!is.null(timeSeriesPlotReactive())){
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
  filename = function(){"sigificant_featires.tsv"},
  content = function(file){
    if(!is.null(timeSeriesReactive())){
      write.table(timeSeriesReactive()$features_df, file = file, quote = F,sep = "\t")
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
    if(input$timeSeriesClusterK > 0){
      cluster_meta <- timeSeriesReactive()$cluster_meta
      colnames(cluster_meta)[which(colnames(cluster_meta)==input$timeSeriesBackground)] <- "time_points" 
      ggplot(cluster_meta, aes(y=sample_group))+
        geom_bar(aes(fill=as.character(time_points)))+
        ggtitle("Composition and Size of individual clusters")+
        ylab("Cluster ID")+
        xlab("Number of samples in cluster")
    }
  }
)

output$timeSeriesClusterContent <- renderDataTable({
  if(!is.null(timeSeriesReactive())){
    if(input$timeSeriesClusterK > 0){
      cluster_meta <- timeSeriesReactive()$cluster_meta
      dt <- cluster_meta[,c("SampleID","sample_group")]
      colnames(dt) <- c("Sample ID","Cluster ID")
      dt
    }
  }
})

observe({
  if(!is.null(currentSet())){
    if(input$timeSeriesClusterK > 0){
      updateSelectInput(session, "timeSeriesMeanLine", selected="NULL")
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
          df <- data.frame(relative_abundance = abundance, group = group_vector)
  
          # perform wilcoxon test for all pairs of sample-groups; return if any pair is significantly different
          if(input$statTestMethod == "Wilcoxon test"){
            fit <- compare_means(relative_abundance~group, data=df, p.adjust.method = input$statTestPAdjust)  
            if(any(fit$p.adj < input$statTestCutoff)){
              # save pairs for which test was performed
              fit$pair <- paste0(fit$group1," vs. ", fit$group2)
              fit$pair_display <- paste0(fit$group1," vs. ", fit$group2, " (pval:", fit$p.adj,")")
              return(list(data=df,
                          fit_table=fit,
                          tax_name = i))
            }  
          }
          # perform KW test between all groups of selected meta variable
          else if(input$statTestMethod == "Kruskal-Wallis test"){
            fit <- kruskal.test(relative_abundance ~ group, data=df)
            fit$p.value <- p.adjust(fit$p.value, method=input$statTestPAdjust)  # this is technically useless, since correcting a single p-value makes no sense..
            if(fit$p.value < input$statTestCutoff){
              return(list(data=df,
                          fit=fit,
                          tax_name = i))
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
