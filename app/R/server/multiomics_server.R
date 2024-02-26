### multi-omics methods ###

read_metabolomics_expr <- function(file){
  # read expression file, can be read as an OTU table
  expr <- read_csv_custom(file, file_type = "otu")
  # check for 15 or more samples
  if(ncol(expr)<15) stop("Metabolomics expression table has to contain a minimum of 15 samples.")
  # check common sample IDs
  otu_samples <- colnames(vals$datasets[[currentSet()]]$phylo@otu_table)
  expr_samples <- colnames(expr)
  if(length(intersect(otu_samples, expr_samples))<15) stop("Less than 15 shared samples between OTU table and metabolomics expression found.")
  return(expr)
}
# manage warning
observe({
  if(!is.null(currentSet())) {
    if(vals$datasets[[currentSet()]]$normMethod == 0) {
      show("norm_warning")
    } else {
      hide("norm_warning")
    }
  }
})

observeEvent(input$upload_metabolomics,{
  if(is.null(currentSet())) {
    showModal(errorModal(error_message = "Please first create a dataset by uploading data to add the metabolomics expression to."))
  }
  else if(!is.null(input$metabolomicsExpressionFile$datapath)){
    tryCatch({
      metabolomicsExpression <- read_metabolomics_expr(input$metabolomicsExpressionFile$datapath[1])
      # assign expression to slot in current dataset
      vals$datasets[[currentSet()]]$metabolomicsExpression <- metabolomicsExpression
      vals$datasets[[currentSet()]]$has_metabolomics <- T
      message(paste0(Sys.time(), " - Successfully loaded metabolomics expression data with dimensions: ", dim(metabolomicsExpression)))
    },
    error=function(e){
      message(e$message)
      showModal(errorModal(error_message = e$message))
    })
  } else {
    showModal(errorModal(error_message = "Please first upload a metabolomics expression file to continue."))
  }
})

#### MOFA2 ####

run_namco_mofa2 <- function(){
  waiter_show(html = tagList(spin_rotating_plane(), "Preparing MOFA2 data ..."),color=overlay_color)
  # align samples of otu table and metabolomics expression 
  otu <- as.matrix(as.data.frame(vals$datasets[[currentSet()]]$phylo@otu_table))
  expr <- as.matrix(vals$datasets[[currentSet()]]$metabolomicsExpression)
  common_samples <- sort(intersect(colnames(otu), colnames(expr)))
  # build matrix list
  matrix_list <- list(
    "microbiome"=otu[,common_samples],
    "metabolites"=expr[,common_samples]
  )
  MOFA_object <- NULL
  condition_labels <- NULL
  group_labels <- NULL
  # let user assign conditions and groups if that option is enabled
  if(input$mofa2_sample_label!=""){
    meta <- vals$datasets[[currentSet()]]$metaData
    rownames(meta) <- meta[[input$mofa2_sample_label]]
    
    # set conditions
    if(input$mofa2_condition_label!=""){
      condition_labels <- meta[common_samples,input$mofa2_condition_label]
    }
    # set groups
    if(input$mofa2_grouped_run && input$mofa2_group_label!=""){
      group_labels <- meta[common_samples,input$mofa2_group_label]
    }
  }
  
  # create initial MOFA object
  MOFA_object <- create_mofa(matrix_list, groups = group_labels)
  # save overview plot
  data_overview <- plot_data_overview(MOFA_object)
  #### save data options
  data_opts <- get_default_data_options(MOFA_object)
  # user input/ default values
  data_opts$scale_views <- input$mofa2_scale_views
  data_opts$scale_groups <- input$mofa2_scale_groups
  #### save model options
  model_opts <- get_default_model_options(MOFA_object)
  # user input/ default values
  model_opts$num_factors <- input$mofa2_num_factors
  model_opts$likelihoods <- setNames(rep(input$mofa2_likelihood, length(model_opts$likelihoods)),
                                     names(model_opts$likelihoods))
  model_opts$spikeslab_factors <- input$mofa2_spikeslab_factors
  model_opts$spikeslab_weights <- input$mofa2_spikeslab_weights
  model_opts$ard_factors <- input$mofa2_ard_factors
  model_opts$ard_weights <- input$mofa2_ard_weights
  #### save train options
  train_opts <- get_default_training_options(MOFA_object)
  train_opts$verbose <- T   # always use verbose mode for debugging
  train_opts$seed <- seed   # global seed
  # user input/ default values
  train_opts$maxiter <- input$mofa2_maxiter
  train_opts$convergence_mode <- input$mofa2_convergence_mode
  train_opts$stochastic <- input$mofa2_stochastic
  ### add options to MOFA object
  MOFA_object <- prepare_mofa(
    object = MOFA_object,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  # build session and run specific output dir
  mofa2_run_id <- paste(sample(1:9, 4, replace = T), collapse = "")
  mofa2_dir <- paste0("mofa2_", sessionID, "_", mofa2_run_id)
  mofa2_tmp_dir <- file.path(tempdir(), mofa2_dir)
  # build hdf5 output file
  mofa2_tmp_out_file <- file.path(mofa2_tmp_dir, "model.hdf5")
  
  waiter_update(html = tagList(spin_rotating_plane(), "Setting up environment ..."))
  # set conda env with mofapy2
  reticulate::use_condaenv(condaenv = "namco_env")
  ### run MOFA2
  waiter_update(html = tagList(spin_rotating_plane(), "Running MOFA2 ..."))
  MOFA_object_trained <- run_mofa(MOFA_object, mofa2_tmp_out_file, use_basilisk = F)
  
  if(is.null(condition_labels)) condition_labels <- "none"
  if(is.null(group_labels)) group_labels <- "group1"
  # add meta data
  mofa_meta <- data.frame(
    sample = common_samples,
    condition = condition_labels,
    group = group_labels
  )
  
  samples_metadata(MOFA_object_trained) <- mofa_meta
  
  # explained variation plots
  explained_variance_plot <- plot_variance_explained(MOFA_object_trained, x="group")
  
  waiter_hide()
  # build result object
  return(list(
    # plots
    "data_overview"=data_overview,
    "explained_variance_plot"=explained_variance_plot,
    # objects
    "mofa_model"=MOFA_object_trained
  ))
}

mofa2_data <- eventReactive(input$mofa2_start,{
  if(vals$datasets[[currentSet()]]$has_metabolomics) {
    tryCatch({
      data <- run_namco_mofa2()
      modal <- modalDialog(
        title = p("Success!", style="color:green; font-size:40px"),
        HTML(paste0("MOFA2 finished successfully")),
        easyClose = T, size = "l"
      )
      showModal(modal)
      return(data)
    },error=function(e){
      waiter_hide()
      print(e$message)
      showModal(errorModal(e$message))
    })
  }
})

# observe mofa2 output and update some options
observe({
  if(!is.null(mofa2_data())){
    shinyjs::show("mofa2_object_overview")
    shinyjs::show("mofa2_results_div")
    m <- mofa2_data()$mofa_model@dimensions$K
    updateSliderInput(session, "mofa2_selected_factors", min = 1, max = m, step = 1, value = c(1,3))
    updateSliderInput(session, "mofa2_selected_top_factors", min = 1, max = m, step = 1, value = c(1,3))
    updateSliderInput(session, "mofa2_scatter_selected_factors", min = 1, max = m, step = 1, value = c(1,3))
    updateSliderInput(session, "mofa2_selected_impact_factors", min = 1, max = m , step = 1, value = 1)
  }
})

observe({
  if(input$mofa2_grouped_run){
    enable("mofa2_group_label")
  } else {
    disable("mofa2_group_label")
  }
})

# data overview plot
output$mofa2_data_overview <- renderPlot({
  if(!is.null(mofa2_data()$data_overview)) {
    mofa2_data()$data_overview
  }
})

# explained variance plot
output$mofa2_explained_variance_group <- renderPlot({
  if(!is.null(mofa2_data()$explained_variance_plot)) {
    mofa2_data()$explained_variance_plot +
      theme(text = element_text(size=input$mofa2_variance_text_size))
  }
})

check_factors <- function(factors, model){
  n_factors <- model@dimensions$K
  factors <- factors[1]:factors[2]
  factors <- intersect(factors, 1:n_factors)
  if(length(factors)==0) factors <- 1:3
  return(factors)
}

# plot factors
output$mofa2_factor_plot <- renderPlot({
  if(!is.null(mofa2_data()$mofa_model)){
    # collect inputs
    model <- mofa2_data()$mofa_model
    factors <- check_factors(input$mofa2_selected_factors, model)
    # call plot function
    p <- plot_factor(model,
                     factors = factors,
                     add_violin = input$mofa2_factor_violin,
                     violin_alpha = input$mofa2_factor_violin_alpha/100,
                     legend = input$mofa2_factor_legend,
                     dodge = T,
                     color_by = "condition",
                     dot_size = input$mofa2_factor_dot_size)
    # colors
    conditions <- unique(model@samples_metadata$condition)
    condition_cols <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(conditions))
    names(condition_cols) <- conditions
    p <- p + scale_color_manual(values=condition_cols) +
      scale_fill_manual(values=condition_cols) +
      theme(text = element_text(size = input$mofa2_factor_text_size))
    p
  }
})

# factor scatter plot
output$mofa2_factor_scatters <- renderPlot({
  if(!is.null(mofa2_data()$mofa_model)){
    model <- mofa2_data()$mofa_model
    factors <- check_factors(input$mofa2_scatter_selected_factors, model)
    dot_size <- input$mofa2_factor_scatter_dot_size
    object <- model
    p <- MOFA2::plot_factors(object,
                             factors = factors,
                             color_by = "condition", 
                             dot_size = dot_size)
    # add correct colors
    conditions <- unique(model@samples_metadata$condition)
    condition_cols <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(conditions))
    names(condition_cols) <- conditions
    p <- p +
      scale_fill_manual(values=condition_cols) +
      scale_color_manual(values=condition_cols) +
      theme(text = element_text(size = input$mofa2_factor_scatter_text_size))
    p
  }
})

output$mofa2_top_weights <- renderPlot({
  if(!is.null(mofa2_data()$mofa_model)){
    waiter_show(html = tagList(spin_rotating_plane(), "Preparing result plots ..."),color=overlay_color)
    # collect model
    model <- mofa2_data()$mofa_model
    factors <- check_factors(input$mofa2_selected_top_factors, model)
    views <- c()
    if(input$mofa2_top_weights_show_microbiome) views <- c(views,"microbiome")
    if(input$mofa2_top_weights_show_metabolomics) views <- c(views,"metabolites")
    # extract weights from model
    W <- get_weights(model, factors = factors, views = views, 
                     as.data.frame = TRUE)
    
    # merge otus with the same taxa level
    if (input$mofa2_top_weights_taxa){
      waiter_update(html = tagList(spin_rotating_plane(), "Aggregating taxa with the same rank ..."))
      # convert OTUs to taxa
      phylo <- vals$datasets[[currentSet()]]$phylo
      tax_info <- as.data.frame(phylo@tax_table)
      taxa_level <- input$mofa2_taxa_level
      # mean weights
      Ws <- split(W, W$view)
      Ws$microbiome <- as.data.table(merge(Ws$microbiome, tax_info, by.x = "feature", by.y = 0))
      Ws$microbiome[,value:=mean(value,na.rm=T),by=c(taxa_level, "factor")]
      Ws$microbiome <- as.data.frame(Ws$microbiome)
      Ws$microbiome <- unique(Ws$microbiome[,c("factor", "view", taxa_level, "value")])
      colnames(Ws$microbiome)[3] <- "feature"
      W <- rbind(Ws$microbiome[,colnames(Ws$metabolites)], Ws$metabolites)
    }
    
    # scale weights
    W$value <- W$value/max(abs(W$value))
    # select top features
    W <- W[W$value != 0, ]
    W$sign <- ifelse(W$value > 0, "+", "-")
    max_w <- W[W$sign=="+",] %>% slice_max(order_by = value, n = input$mofa2_top_n_weights)
    min_w <- W[W$sign=="-",] %>% slice_min(order_by = value, n = input$mofa2_top_n_weights)
    W <- rbind(max_w, arrange(min_w, desc(value)))
    # order y bars
    W$feature <- factor(W$feature, levels = unique(W$feature))
    # get correct colors
    namco_colors <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(factors))
    names(namco_colors) <- paste0("Factor",factors)
    waiter_hide()
    # plot
    p <- ggplot(W, aes(x=value, y=feature, fill=factor)) +
      geom_col(position="identity") +
      xlab("scaled weights") + ylab("Feature") +
      geom_vline(xintercept = 0) +
      theme(text = element_text(size=input$mofa2_weight_text_size)) +
      scale_fill_manual(values=namco_colors) +
      facet_wrap(~view)
    if(input$mofa2_top_weights_fix_x_axis) p <- p + xlim(c(-1,1))
    p
  }
})


output$mofa2_microbiome_scatter <- renderPlot({
  if(!is.null(mofa2_data()$mofa_model)){
    model <- mofa2_data()$mofa_model
    factors <- input$mofa2_selected_impact_factors
    
    # get colors
    conditions <- unique(model@samples_metadata$condition)
    condition_cols <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(conditions))
    names(condition_cols) <- conditions
    # plot
    plot_data_scatter(model,
                      view = 1,
                      factor = factors,
                      features = input$mofa2_microbiome_scatter_features,
                      add_lm = T,
                      color_by = "condition"
    ) + theme(text = element_text(size=input$mofa2_microbiome_scatter_text_size)) +
      scale_color_manual(values=condition_cols)
  }
})

