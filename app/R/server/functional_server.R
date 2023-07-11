####picrust2 ####

observe({
  if(input$picrustTestNormalization == "clr"){
    shinyjs::show("aldex2Additional")
  }else{
    shinyjs::hide("aldex2Additional")
  }
  if (!is.null(currentSet())){
    phylo <- vals$datasets[[currentSet()]]$phylo 
    if (vals$datasets[[currentSet()]]$is_fastq || vals$datasets[[currentSet()]]$is_sample_data || !is.null(phylo@refseq)){
      shinyjs::hide("picrustFastaFile")
    }else{
      shinyjs::show("picrustFastaFile")
    }
    if(vals$datasets[[currentSet()]]$has_picrust){
      shinyjs::enable("picrustDiffStart")
    }else{
      shinyjs::disable("picrustDiffStart")
    }
  }
})

output$hasPicrustInfoBox <- renderInfoBox({
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_picrust){
      infoBox("Picrust2 ready", icon = icon("thumbs-up", lib = "glyphicon"), color="green", width=10)
    }else{
      infoBox("Picrust2 not ready", icon = icon("thumbs-down", lib = "glyphicon"), color="red", width=10)
    }
  }
})


#### picrust2 run ####

observeEvent(input$picrust2Start,{
  if(!is.null(currentSet())){
    message(paste0(Sys.time(), " - Starting picrust2 analysis ..."))
    waiter_show(html = tagList(spin_rotating_plane(),"Running Picrust2 ..." ),color=overlay_color)
    
    tryCatch({
      phylo <- vals$datasets[[currentSet()]]$phylo
      vals$datasets[[currentSet()]]$picrust_output <- NULL    # reset old picrust-output variable
      shinyjs::hide("download_picrust_div")
      
      foldername <- sprintf("/picrust2_%s", digest::digest(phylo))  # unique folder name for this output
      outdir <- paste0(tempdir(),foldername)
      if (dir.exists(outdir)){unlink(outdir, recursive = T)}    # remove output directory if it exists
      dir.create(outdir)
      
      biom_file = paste0(outdir,"/biom_picrust.biom")
      biom <- make_biom(data=otu_table(vals$datasets[[currentSet()]]$phylo.raw)) # Picrust2 requires absolute abundances
      write_biom(biom, biom_file)
      message(paste0(Sys.time(), " - Wrote biom-file: ", biom_file))
      
      # use fasta file of dada2 pipeline if available or if sequence data is in phylo object
      if (vals$datasets[[currentSet()]]$is_fastq || !is.null(phylo@refseq)){
        fasta_file <- paste0(outdir,"/seqs.fasta")
        writeXStringSet(refseq(phylo), fasta_file)
        message(paste0(Sys.time(), " - Using dada2-generated fasta file:", fasta_file))
      }else{
        fasta_file <- input$picrustFastaFile$datapath
        if(vals$datasets[[currentSet()]]$is_sample_data) fasta_file <- "testdata/seqs.fasta"
        fasta <- readDNAStringSet(fasta_file, format="fasta")
        all_taxa <- taxa_names(phylo)
        missing_otus <- setdiff(all_taxa, names(fasta))
        # check if all OTUs have a sequence in fasta file
        if(length(missing_otus)>0) {
          warn <- paste0("Some OTU IDs do not have a record in the provided fasta file: ",
                         paste(missing_otus, collapse = ", "), "\n",
                         "These OTUs will be excluded from the picrust2 analyses.")
          modal<-modalDialog(
            title = p("Warning:", style="color:orange; font-size:40px"),
            warn,
            easyClose = T, size=s
          )
          showModal(modal)
        }
        fasta <- fasta[fasta@ranges@NAMES %in% all_taxa]
        writeXStringSet(fasta, fasta_file)
        message(paste0(Sys.time(), " - Using user-uploaded fasta file: ", fasta_file))
      }
      
      picrust_outdir <- paste0(outdir,"/picrust2out")       # this is the name of the final output directory of this picrust run
      
      if(!vals$datasets[[currentSet()]]$is_sample_data){
        command = paste0(namco_conda_env, " picrust2_pipeline.py --remove_intermediate -s ",fasta_file," -i ",biom_file, " -o ", picrust_outdir, " -p", ncores)
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
        p2_EC_tmp <- paste0(picrust_outdir, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
        p2_KO_tmp <- paste0(picrust_outdir, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
        p2_PW_tmp <- paste0(picrust_outdir, "/pathways_out/path_abun_unstrat.tsv.gz") 
        marker_nsti <- paste0(picrust_outdir, "/marker_predicted_and_nsti.tsv.gz")
        
        #generate descriptions
        p2_EC <- paste0(picrust_outdir,"/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz")
        command_EC <- paste0(namco_conda_env, " add_descriptions.py -i ",p2_EC_tmp, " -m EC -o ",p2_EC)
        out_EC <- system(command_EC, wait=T) 
        p2_KO <- paste0(picrust_outdir,"/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz")
        command_KO <- paste0(namco_conda_env, " add_descriptions.py -i ",p2_KO_tmp, " -m KO -o ",p2_KO)
        out_KO <- system(command_KO, wait=T) 
        p2_PW<- paste0(picrust_outdir,"/pathways_out/path_abun_unstrat_descrip.tsv.gz")
        command_PW <- paste0(namco_conda_env, " add_descriptions.py -i ",p2_PW_tmp, " -m METACYC -o ",p2_PW)
        out_PW <- system(command_PW, wait=T) 
        
      }else{
        p2_EC <- "testdata/picrust2/ec_pred_metagenome_unstrat_descrip.tsv.gz"
        p2_KO <- "testdata/picrust2/ko_pred_metagenome_unstrat_descrip.tsv.gz"
        p2_PW <- "testdata/picrust2/path_abun_unstrat_descrip.tsv.gz"
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
          vals$datasets[[currentSet()]]$normMethod <- 5
        }
        
        p2EC = as.data.frame(fread(p2_EC))
        rownames(p2EC) = p2EC[["function"]]
        p2EC_descr <- p2EC[,c("function","description")]
        p2EC[["description"]]<-NULL
        p2EC = as.matrix(p2EC[,-1])
        p2EC = round(p2EC)
        
        p2KO = as.data.frame(fread(p2_KO))
        rownames(p2KO) = p2KO$"function"
        p2KO_descr <- p2KO[,c("function","description")]
        p2KO[["description"]]<-NULL
        p2KO = as.matrix(p2KO[,-1])
        p2KO = round(p2KO)
        
        p2PW = as.data.frame(fread(p2_PW))
        rownames(p2PW) = p2PW$"pathway"
        p2PW_descr <- p2PW[,c("pathway","description")]
        colnames(p2PW_descr)<-c("function","description")
        p2PW[["description"]]<-NULL
        p2PW = as.matrix(p2PW[,-1])
        p2PW = round(p2PW)
        
        vals$datasets[[currentSet()]]$picrust_results_list <- list(p2EC=p2EC,
                                                                   p2KO=p2KO,
                                                                   p2PW=p2PW, 
                                                                   p2EC_rel=relAbundance(p2EC),
                                                                   p2KO_rel=relAbundance(p2KO),
                                                                   p2PW_rel=relAbundance(p2PW),
                                                                   p2EC_descr=p2EC_descr,
                                                                   p2KO_descr=p2KO_descr,
                                                                   p2PW_descr=p2PW_descr)

        message(paste0(Sys.time(), " - Finished picrust2"))
        waiter_hide()
        
      }else{
        message(paste0(Sys.time(), " - Did not find all files for differential analysis; stopping ... "))
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


#### picrust2 differential analysis ####

observeEvent(input$picrustDiffStart, {
  if(!is.null(currentSet())){
    if(vals$datasets[[currentSet()]]$has_picrust){
      message(paste0(Sys.time(), " - Starting picrust2 differential analysis ..."))
      waiter_show(html = tagList(spin_rotating_plane(),"Running Picrust2 differential analysis ..." ),color=overlay_color)
      
      tryCatch({
        phylo <- vals$datasets[[currentSet()]]$phylo
        picrust_results <- vals$datasets[[currentSet()]]$picrust_results_list
        p2EC <- picrust_results$p2EC
        p2PW <- picrust_results$p2PW
        p2KO <- picrust_results$p2KO
        
        # do not run differential analysis if selected group does not exist
        if(is.null(sample_data(phylo)[[input$picrust_test_condition]])){
          message(paste0(Sys.time(), " - Skipping differential analysis; sample group not found ... "))
          stop(picrustDifferentialGroupNotFoundError, call. = F)
        }
        
        # run ALDEx2 to perform differential abundance testing between conditions 
        # for t & wilcox: need to select test condition and covariate, against which to compare all others
        
        meta <- as.data.frame(phylo@sam_data, check.names=F)
        sample_vector <- meta[[input$picrust_test_condition]]
        names(sample_vector) <- sample_names(phylo)
        if(input$picrustTest != "kw"){
          sample_vector[which(sample_vector != input$picrust_test_covariate)] <- "other"
        }
        
        waiter_update(html = tagList(spin_rotating_plane(),"Analyzing EC ..."))
        test_EC <- picrust2_statistical_analysis(abundances = p2EC, 
                                                 normalization =  input$picrustTestNormalization, 
                                                 test = input$picrustTest, 
                                                 sample_vector = sample_vector,
                                                 test_covariate = input$picrust_test_covariate, 
                                                 mc.samples = input$picrust_mc_samples)
        waiter_update(html = tagList(spin_rotating_plane(),"Analyzing PW ..."))
        test_PW <- picrust2_statistical_analysis(abundances = p2PW, 
                                                 normalization =  input$picrustTestNormalization, 
                                                 test = input$picrustTest, 
                                                 sample_vector = sample_vector,
                                                 test_covariate = input$picrust_test_covariate, 
                                                 mc.samples = input$picrust_mc_samples)
        waiter_update(html = tagList(spin_rotating_plane(),"Analyzing KO ..."))
        test_KO <- picrust2_statistical_analysis(abundances = p2KO, 
                                                 normalization =  input$picrustTestNormalization, 
                                                 test = input$picrustTest, 
                                                 sample_vector = sample_vector,
                                                 test_covariate = input$picrust_test_covariate, 
                                                 mc.samples = input$picrust_mc_samples)
        
        vals$datasets[[currentSet()]]$picrust_analysis_list <- list(test_EC=test_EC$test,
                                                                    test_KO=test_KO$test,
                                                                    test_PW=test_PW$test,
                                                                    abundances_EC=test_EC$abundances,
                                                                    abundances_KO=test_KO$abundances,
                                                                    abundances_PW=test_PW$abundances,
                                                                    label = sample_vector)
        waiter_hide()
        
      }, error=function(e){
        waiter_hide()
        print(e$message)
        showModal(errorModal(e$message))
      })
    }
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
    if(!is.null(vals$datasets[[currentSet()]]$picrust_analysis_list) && vals$datasets[[currentSet()]]$has_picrust){
      message(Sys.time(), " - generating picrust analysis tables ...")
      abundances <- list(vals$datasets[[currentSet()]]$picrust_results_list$p2EC_rel,
                         vals$datasets[[currentSet()]]$picrust_results_list$p2KO_rel,
                         vals$datasets[[currentSet()]]$picrust_results_list$p2PW_rel)
      descriptions <- list(vals$datasets[[currentSet()]]$picrust_results_list$p2EC_descr,
                           vals$datasets[[currentSet()]]$picrust_results_list$p2KO_descr,
                           vals$datasets[[currentSet()]]$picrust_results_list$p2PW_descr)
      test_results <- list(vals$datasets[[currentSet()]]$picrust_analysis_list$test_EC,
                            vals$datasets[[currentSet()]]$picrust_analysis_list$test_KO,
                            vals$datasets[[currentSet()]]$picrust_analysis_list$test_PW)
      
      nsamples <- nsamples(vals$datasets[[currentSet()]]$phylo)
      label <- vals$datasets[[currentSet()]]$picrust_analysis_list$label
      has_effect_measure <- F
      pval_var <- "pval1"  # can change to "pvalAdj1" to use adjusted p-values for significance test

      # order: 1=EC, 2=KO, 3=PW
      out_list<-lapply(c(1,2,3), function(x){
        test.tab <- test_results[[x]]
        abundances.tab <- abundances[[x]]
        description <- descriptions[[x]]
        
        # create df with all functions and two entries per func for the 2 pvalues
        if(input$picrustTestNormalization == "clr" && input$picrustTest %in% c("t","wilcox")){
          test.tab$significant <- ifelse((test.tab[[pval_var]]<input$picrust_signif_lvl & test.tab[["effect"]]>input$picrust_signif_lvl_effect),T,F)
          has_effect_measure <<- T
        }else{
          test.tab$significant <- ifelse((test.tab[[pval_var]]<input$picrust_signif_lvl),T,F)
        }
        test.tab <- merge(test.tab, description, by.x="func",by.y="function")
        tab.long <- gather(test.tab, pvalue_type, pvalue, 4:5)

        # create df with only significant functions
        x.merged <- merge(abundances.tab, test.tab, by.x=0, by.y="func")
        x.merged$func<-x.merged$Row.names
        x.merged[["Row.names"]] <- NULL
        x.merged.long <- gather(x.merged,SampleID,abundance,1:nsamples)
        x.merged.long$label<-unlist(lapply(x.merged.long$SampleID, function(y){return(label[[y]])}))
        
        # use only features which are selected, otherwise use all signif features
        if(x==1){
          feature_to_add <- input$picrust_ec_select
        }else if(x==2){
          feature_to_add <- input$picrust_ko_select
        }else if(x==3){
          feature_to_add <- input$picrust_pw_select
        }
        if(!is.null(feature_to_add)){
          signif_df <- x.merged.long[x.merged.long$func %in% feature_to_add,]
        }else{
          signif_df <- data.frame(x.merged.long[x.merged.long$significant==T,])
        }
        signif_df[["pval1"]] <- -log10(signif_df[["pval1"]])
        signif_df <- signif_df[order(signif_df[["pval1"]], decreasing = F),]
        
        return(list(tab.long, signif_df))
      })
      
      out <- list(EC_long=out_list[[1]][[1]], KO_long=out_list[[2]][[1]], PW_long=out_list[[3]][[1]],
                  EC_signif=out_list[[1]][[2]], KO_signif=out_list[[2]][[2]], PW_signif=out_list[[3]][[2]],
                  has_effect_measure = has_effect_measure)
      return(out)
    }
  }
})

output$picrustDiffDownloadEC <- downloadHandler(
  filename = function(){
    paste("picrust2_differential_analysis_EC.tab")
  },
  content = function(file){
    if(!is.null(aldex_reactive())){
      write.table(aldex_reactive()$EC_long, file = file, quote = F, sep = "\t")
    }
  }
)

output$picrustDiffDownloadKO <- downloadHandler(
  filename = function(){
    paste("picrust2_differential_analysis_KO.tab")
  },
  content = function(file){
    if(!is.null(aldex_reactive())){
      write.table(aldex_reactive()$KO_long, file = file, quote = F, sep = "\t")
    }
  }
)

output$picrustDiffDownloadPW <- downloadHandler(
  filename = function(){
    paste("picrust2_differential_analysis_PW.tab")
  },
  content = function(file){
    if(!is.null(aldex_reactive())){
      write.table(aldex_reactive()$PW_long, file = file, quote = F, sep = "\t")
    }
  }
)


####picrust2 plots ####

picrust_plots_reactive <- reactive({
  if(!is.null(aldex_reactive())){
    #### pvalue plots ####
    EC_long <- aldex_reactive()$EC_long
    KO_long <- aldex_reactive()$KO_long
    PW_long <- aldex_reactive()$PW_long
    
    if(aldex_reactive()$has_effect_measure){
      ec_effect_plot<-ggplot(data=EC_long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
        geom_point(alpha=0.8)+
        ggtitle("Effect size vs P-value")+
        scale_color_manual(labels=c("pvalAdj1"="BH-adjusted","pval1"="P-value"), values=c("#0072B2", "#E18E4C"))+
        geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
        geom_vline(xintercept = input$picrust_signif_lvl_effect, color="black", linetype="dashed", alpha=0.8)+
        theme_minimal()
      if(input$picrust_signif_label){
        ec_effect_plot<-ec_effect_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
      }
      ko_effect_plot<-ggplot(data=KO_long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
        geom_point(alpha=0.8)+
        ggtitle("Effect size vs P-value")+
        scale_color_manual(labels=c("pvalAdj1"="BH-adjusted","pval1"="P-value"), values=c("#0072B2", "#E18E4C"))+
        geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
        geom_vline(xintercept = input$picrust_signif_lvl_effect, color="black", linetype="dashed", alpha=0.8)+
        theme_minimal()
      if(input$picrust_signif_label){
        ko_effect_plot<-ko_effect_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
      }
      pw_effect_plot<-ggplot(data=PW_long,aes(x=effect,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
        geom_point(alpha=0.8)+
        ggtitle("Effect size vs P-value")+
        scale_color_manual(labels=c("pvalAdj1"="BH-adjusted","pval1"="P-value"), values=c("#0072B2", "#E18E4C"))+
        geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
        geom_vline(xintercept = input$picrust_signif_lvl_effect, color="black", linetype="dashed", alpha=0.8)+
        theme_minimal()
      if(input$picrust_signif_label){
        pw_effect_plot<-pw_effect_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
      }
    }else{
      ec_effect_plot <- NULL
      ko_effect_plot <- NULL
      pw_effect_plot <- NULL
    }
    
    ec_vulcano_plot<-ggplot(data=EC_long,aes(x=diff,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Difference vs P-value")+
      scale_color_manual(labels=c("pvalAdj1"="BH-adjusted","pval1"="P-value"), values=c("#0072B2", "#E18E4C"))+
      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
      theme_minimal()
    if(input$picrust_signif_label){
      ec_vulcano_plot<-ec_vulcano_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    ko_vulcano_plot<-ggplot(data=KO_long,aes(x=diff,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Difference vs P-value")+
      scale_color_manual(labels=c("pvalAdj1"="BH-adjusted","pval1"="P-value"), values=c("#0072B2", "#E18E4C"))+
      geom_hline(yintercept = -log10(input$picrust_signif_lvl), color="black", linetype="dashed",alpha=0.8)+
      theme_minimal()
    if(input$picrust_signif_label){
      ko_vulcano_plot<-ko_vulcano_plot+geom_label_repel(aes(label=ifelse(significant,as.character(func),"")), box.padding = 0.3,point.padding = 0.2,color="black",max.overlaps = input$picrust_maxoverlaps)
    }
    pw_vulcano_plot <-ggplot(data=PW_long,aes(x=diff,y=-log10(pvalue),color=pvalue_type, label=as.character(func)))+
      geom_point(alpha=0.8)+
      ggtitle("Difference vs P-value")+
      scale_color_manual(labels=c("pvalAdj1"="BH-adjusted","pval1"="P-value"), values=c("#0072B2", "#E18E4C"))+
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
      if(input$picrust_show_descripton_ec){
        colnames(EC_signif)[which(colnames(EC_signif)=="description")]<-"func_label"
      }else{
        colnames(EC_signif)[which(colnames(EC_signif)=="func")]<-"func_label"
      }
      EC_signif[["func_label"]]<-as.character(EC_signif[["func_label"]])

      p1 <- ggplot(data=EC_signif, aes(x=abundance,y=func_label,fill=as.factor(label)))+
        geom_boxplot()+
        theme_minimal()+
        theme(legend.position = c(0.9,0.9),
              legend.title = element_blank(),
              axis.text=element_text(size=input$picrust_ylab_size_ec))+
        scale_fill_brewer(palette = "Set1")+
        ylab("Function")+xlab('relative abundance')
      
      p2 <- ggplot(data=unique(EC_signif[,c("func_label","pval1", "label")]),aes(x=pval1,y=func_label))+
        geom_bar(stat="identity", width=0.35, aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none")+
        xlab("-log10(P-value)")+
        geom_vline(xintercept = -log10(input$picrust_signif_lvl), color="red", linetype="dashed")
      
      if(aldex_reactive()$has_effect_measure){
        p3 <- ggplot(data=unique(EC_signif[,c("func_label","effect", "label")]),aes(x=effect/nsamples,y=func_label))+
          geom_bar(stat="identity", width=0.35,aes(alpha=0.8))+
          theme_bw()+
          theme(axis.title.y=element_blank(),
                axis.ticks.y =element_blank(),
                axis.text.y = element_blank(), 
                legend.position = "none")+
          xlab("Effect size")+
          geom_vline(xintercept = input$picrust_signif_lvl_effect, color="red", linetype="dashed")
        ec_signif_plot_list <- list(p1=p1,p2=p2,p3=p3)
      }else{
        ec_signif_plot_list <- list(p1=p1,p2=p2) 
      }
    }else{
      ec_signif_plot_list <- NULL
    }
    
    KO_signif <- aldex_reactive()$KO_signif
    if(dim(KO_signif)[1] > 0){
      if(dim(KO_signif)[1] > input$picrust_ko_signif_plot_show){
        rows_to_keep <- input$picrust_ko_signif_plot_show*nsamples
        KO_signif <- na.omit(KO_signif[1:rows_to_keep,])
      }    
      if(input$picrust_show_descripton_ko){
        colnames(KO_signif)[which(colnames(KO_signif)=="description")]<-"func_label"
      }else{
        colnames(KO_signif)[which(colnames(KO_signif)=="func")]<-"func_label"
      }
      KO_signif[["func_label"]]<-as.character(KO_signif[["func_label"]])
      
      p1 <- ggplot(data=KO_signif, aes(x=abundance,y=func_label,fill=as.factor(label)))+
        geom_boxplot()+
        theme_minimal()+
        theme(legend.position = c(0.9,0.9),
              legend.title = element_blank(),
              axis.text=element_text(size=input$picrust_ylab_size_ko))+
        scale_fill_brewer(palette = "Set1")+
        ylab("Function")+xlab('relative abundance')
      
      p2 <- ggplot(data=unique(KO_signif[,c("func_label","pval1", "label")]),aes(x=pval1,y=func_label))+
        geom_bar(stat="identity", width=0.35, aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none")+
        xlab("-log10(P-value)")+
        geom_vline(xintercept = -log10(input$picrust_signif_lvl), color="red", linetype="dashed")
      
      if(aldex_reactive()$has_effect_measure){
        p3 <- ggplot(data=unique(KO_signif[,c("func_label","effect", "label")]),aes(x=effect,y=func_label))+
          geom_bar(stat="identity", width=0.35,aes(alpha=0.8))+
          theme_bw()+
          theme(axis.title.y=element_blank(),
                axis.ticks.y =element_blank(),
                axis.text.y = element_blank(), 
                legend.position = "none")+
          xlab("Effect size")+
          geom_vline(xintercept = input$picrust_signif_lvl_effect, color="red", linetype="dashed")
        ko_signif_plot_list <- list(p1=p1,p2=p2,p3=p3)
      }else{
        ko_signif_plot_list <- list(p1=p1,p2=p2) 
      }
    }else{
      ko_signif_plot_list <- NULL
    }
    
    PW_signif <- aldex_reactive()$PW_signif
    if(dim(PW_signif)[1] > 0){
      if(dim(PW_signif)[1] > input$picrust_pw_signif_plot_show){
        rows_to_keep <- input$picrust_pw_signif_plot_show*nsamples
        PW_signif <- na.omit(PW_signif[1:rows_to_keep,])
      }     
      if(input$picrust_show_descripton_pw){
        colnames(PW_signif)[which(colnames(PW_signif)=="description")]<-"func_label"
      }else{
        colnames(PW_signif)[which(colnames(PW_signif)=="func")]<-"func_label"
      }
      PW_signif[["func_label"]]<-as.character(PW_signif[["func_label"]])
      
      p1 <- ggplot(data=PW_signif, aes(x=abundance,y=func_label,fill=as.factor(label)))+
        geom_boxplot()+
        theme_minimal()+
        theme(legend.position = c(0.9,0.9),
              legend.title = element_blank(),
              axis.text=element_text(size=input$picrust_ylab_size_pw))+
        scale_fill_brewer(palette = "Set1")+
        ylab("Function")+xlab('relative abundance')
      
      p2 <- ggplot(data=unique(PW_signif[,c("func_label","pval1", "label")]),aes(x=pval1,y=func_label))+
        geom_bar(stat="identity", width=0.35, aes(alpha=0.8))+
        theme_bw()+
        theme(axis.title.y=element_blank(),
              axis.ticks.y =element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none")+
        xlab("-log10(P-value)")+
        geom_vline(xintercept = -log10(input$picrust_signif_lvl), color="red", linetype="dashed")
      
      if(aldex_reactive()$has_effect_measure){
        p3 <- ggplot(data=unique(PW_signif[,c("func_label","effect", "label")]),aes(x=effect,y=func_label))+
          geom_bar(stat="identity", width=0.35,aes(alpha=0.8))+
          theme_bw()+
          theme(axis.title.y=element_blank(),
                axis.ticks.y =element_blank(),
                axis.text.y = element_blank(), 
                legend.position = "none")+
          xlab("Effect size")+
          geom_vline(xintercept = input$picrust_signif_lvl_effect, color="red", linetype="dashed")
        pw_signif_plot_list <- list(p1=p1,p2=p2,p3=p3)
      }else{
        pw_signif_plot_list <- list(p1=p1,p2=p2) 
      }
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
      if(aldex_reactive()$has_effect_measure){
        grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1))
      }else{
        grid.arrange(l$p1,l$p2,ncol=2,widths=c(3,1))
      }
    }
  }
})

output$picrust_ec_signifPDF <- downloadHandler(
  filename = function(){"EC_significant_boxplots.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      l <- picrust_plots_reactive()$ec_signif_plot_list
      if(aldex_reactive()$has_effect_measure){
        ggsave(file, grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1)), device="pdf", width = 12, height = 8)
      }else{
        ggsave(file, grid.arrange(l$p1,l$p2,ncol=2,widths=c(3,1)), device="pdf", width = 12, height = 8)
      }
    }
  }
)

output$picrust_ko_signif_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    if(!is.null(picrust_plots_reactive()$ko_signif_plot_list)){
      l <- picrust_plots_reactive()$ko_signif_plot_list
      if(aldex_reactive()$has_effect_measure){
        grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1))
      }else{
        grid.arrange(l$p1,l$p2,ncol=2,widths=c(3,1))
      }
    }
  }
})

output$picrust_ko_signifPDF <- downloadHandler(
  filename = function(){"KO_significant_boxplots.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      l <- picrust_plots_reactive()$ko_signif_plot_list
      if(aldex_reactive()$has_effect_measure){
        ggsave(file, grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1)), device="pdf", width = 12, height = 8)
      }else{
        ggsave(file, grid.arrange(l$p1,l$p2,ncol=2,widths=c(3,1)), device="pdf", width = 12, height = 8)
      }
    }
  }
)

output$picrust_pw_signif_plot <- renderPlot({
  if(!is.null(picrust_plots_reactive())){
    if(!is.null(picrust_plots_reactive()$pw_signif_plot_list)){
      l <- picrust_plots_reactive()$pw_signif_plot_list
      if(aldex_reactive()$has_effect_measure){
        grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1))
      }else{
        grid.arrange(l$p1,l$p2,ncol=2,widths=c(3,1))
      }
    }
  }
})

output$picrust_pw_signifPDF <- downloadHandler(
  filename = function(){"PW_significant_boxplots.pdf"},
  content = function(file){
    if(!is.null(picrust_plots_reactive())){
      l <- picrust_plots_reactive()$pw_signif_plot_list
      if(aldex_reactive()$has_effect_measure){
        ggsave(file, grid.arrange(l$p1,l$p2,l$p3,ncol=3,widths=c(3,1,1)), device="pdf", width = 12, height = 8)
      }else{
        ggsave(file, grid.arrange(l$p1,l$p2,ncol=2,widths=c(3,1)), device="pdf", width = 12, height = 8)
      }
    }
  }
)


#### significant features ####

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
  valueBox(val, "Significant ECs",icon = icon("arrow-up"), color="olive")
})

output$picrust_ko_effect_signif_value <- renderValueBox({
  val <- 0
  if(!is.null(aldex_reactive())){
    val <- length(unique(aldex_reactive()$KO_long[aldex_reactive()$KO_long[["significant"]]==T,][["func"]]))
  }
  valueBox(val, "Significant KOs",icon = icon("arrow-up"), color="olive")
})

output$picrust_pw_effect_signif_value <- renderValueBox({
  val <- 0
  if(!is.null(aldex_reactive())){
    val <- length(unique(aldex_reactive()$PW_long[aldex_reactive()$PW_long[["significant"]]==T,][["func"]]))
  }
  valueBox(val, "Significant PWs",icon = icon("arrow-up"), color="olive")
})


