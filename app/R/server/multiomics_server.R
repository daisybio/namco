### multi-omics methods ###

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

#### MOFA2 ####
tidy_expr <- function(mat) {
  m <- as.matrix(mat)
  return(m[,!is.na(colSums(m))])
}

run_namco_mofa2 <- function(){
  waiter_show(html = tagList(spin_rotating_plane(), "Preparing MOFA2 data ..."),color=overlay_color)
  # align samples of otu table
  otu <- as.matrix(as.data.frame(vals$datasets[[currentSet()]]$phylo@otu_table))
  omics <- lapply(vals$datasets[[currentSet()]]$omics, tidy_expr)
  omics_samples <- unique(unlist(lapply(omics, colnames)))
  
  common_samples <- sort(intersect(colnames(otu), omics_samples))
  
  # build matrix list
  matrix_list <- lapply(omics, function(expr) as.matrix(expr[,common_samples]))
  matrix_list[["Microbiome"]] <- otu[,common_samples]
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
  MOFA_object <- MOFA2::create_mofa(matrix_list, groups = group_labels)
  # save overview plot
  data_overview <- MOFA2::plot_data_overview(MOFA_object)
  #### save data options
  data_opts <- MOFA2::get_default_data_options(MOFA_object)
  # user input/ default values
  data_opts$scale_views <- input$mofa2_scale_views
  data_opts$scale_groups <- input$mofa2_scale_groups
  #### save model options
  model_opts <- MOFA2::get_default_model_options(MOFA_object)
  # user input/ default values
  model_opts$num_factors <- input$mofa2_num_factors
  model_opts$likelihoods <- setNames(rep(input$mofa2_likelihood, length(model_opts$likelihoods)),
                                     names(model_opts$likelihoods))
  model_opts$spikeslab_factors <- input$mofa2_spikeslab_factors
  model_opts$spikeslab_weights <- input$mofa2_spikeslab_weights
  model_opts$ard_factors <- input$mofa2_ard_factors
  model_opts$ard_weights <- input$mofa2_ard_weights
  #### save train options
  train_opts <- MOFA2::get_default_training_options(MOFA_object)
  train_opts$verbose <- T   # always use verbose mode for debugging
  train_opts$seed <- seed   # global seed
  # user input/ default values
  train_opts$maxiter <- input$mofa2_maxiter
  train_opts$convergence_mode <- input$mofa2_convergence_mode
  train_opts$stochastic <- input$mofa2_stochastic
  ### add options to MOFA object
  MOFA_object <- MOFA2::prepare_mofa(
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
  MOFA_object_trained <- MOFA2::run_mofa(MOFA_object, mofa2_tmp_out_file, use_basilisk = F)
  
  if(is.null(condition_labels)) condition_labels <- "none"
  if(is.null(group_labels)) group_labels <- "group1"
  # add meta data
  mofa_meta <- data.frame(
    sample = common_samples,
    condition = condition_labels,
    group = group_labels
  )
  
  MOFA2::samples_metadata(MOFA_object_trained) <- mofa_meta
  
  # explained variation plots
  explained_variance_plot <- MOFA2::plot_variance_explained(MOFA_object_trained, x="group")
  
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
  if(vals$datasets[[currentSet()]]$has_omics) {
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
    data <- mofa2_data()
    # extract model
    m <- data$mofa_model
    # build a heatmap for each data view
    conditions <- unique(m@samples_metadata$condition)
    condition_cols <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(conditions))
    names(condition_cols) <- conditions
    
    meta <- data.frame(m@samples_metadata %>% select(-group), row.names = 1)
    hma <- HeatmapAnnotation(df = meta, col = list(condition = condition_cols),
                             annotation_legend_param = list(title_gp = gpar(fontsize = 20),
                                                            labels_gp = gpar(18))
    )
    p <- NULL
    for (g in names(m@data)) {
      d <- m@data[[g]]
      mat <- as.matrix(unlist(d[[1]]))
      mat <- mat[,rownames(meta)]
      hm <- Heatmap(mat, name = g, na_col = "grey", top_annotation = hma,
                    cluster_rows = T, cluster_columns = T,
                    show_row_names = F, show_column_names = F,
                    show_row_dend = F, row_title = g,
                    row_title_gp = gpar(fontsize = 22), 
                    heatmap_legend_param = list(
                      title_gp = gpar(fontsize = 20),
                      labels_gp = gpar(18)
                    ))
      if(is.null(p)) {
        p <- hm
      } else {
        p <- p %v% hm
      }
    }
    draw(p, ht_gap = unit(2, "cm"), padding = unit(c(2, 0.25, 0.25, 0.25), "cm"),
         column_title = paste0("# shared samples: ", m@dimensions$N), 
         column_title_gp = gpar(fontsize = 22))
  }
})

# explained variance plot
output$mofa2_explained_variance_group <- renderPlot({
  if(!is.null(mofa2_data()$explained_variance_plot)) {
    mofa2_data()$explained_variance_plot +
      theme(text = element_text(size=input$mofa2_variance_text_size)) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
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
    p <- MOFA2::plot_factor(model,
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

find_indices_below_diagonal <- function(n) {
  positions <- vector()
  # Loop over each possible row and column index
  for (i in 1:n) {
    for (j in 1:n) {
      if (i > j) {  # Condition for below the diagonal
        # Calculate linear position
        k <- (i-1) * n + j
        positions <- c(positions, k)
      }
    }
  }
  return(positions)
}

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
    # remove half of plot due to redundancy
    kill <- find_indices_below_diagonal(max(factors))
    for (k in kill) {
      p$plots[[k]] <- ggplot() + theme_void()
    }
    
    p <- p +
      scale_fill_manual(values=condition_cols) +
      scale_color_manual(values=condition_cols) +
      theme(text = element_text(size = input$mofa2_factor_scatter_text_size)) +
      theme(legend.position = "none")
    suppressWarnings(p)
  }
})

# update select boxes
observe({
  if(!is.null(mofa2_data()$mofa_model)) {
    updateSelectInput(session, "mofa2_top_weights_omics_selection", choices = names(mofa2_data()$mofa_model@data),
                      selected = names(mofa2_data()$mofa_model@data))
    updateSelectInput(session, "mofa2_microbiome_scatter_view", choices = names(mofa2_data()$mofa_model@data),
                      selected = names(mofa2_data()$mofa_model@data)[1])
  }
})

# (top) factor weights
output$mofa2_top_weights <- renderPlot({
  if(!is.null(mofa2_data()$mofa_model)){
    waiter_show(html = tagList(spin_rotating_plane(), "Preparing result plots ..."),color=overlay_color)
    # collect model
    model <- mofa2_data()$mofa_model
    factors <- check_factors(input$mofa2_selected_top_factors, model)
    views <- input$mofa2_top_weights_omics_selection
    
    # extract weights from model
    W <- MOFA2::get_weights(model, factors = factors, views = views, 
                            as.data.frame = TRUE)
    
    # merge otus with the same taxa level
    if (input$mofa2_top_weights_taxa){
      waiter_update(html = tagList(spin_rotating_plane(), "Aggregating taxa with the same rank ..."))
      # convert OTUs to taxa
      phylo <- vals$datasets[[currentSet()]]$phylo
      tax_info <- as.data.frame(phylo@tax_table)
      taxa_level <- input$mofa2_taxa_level
      # mean weights
      Ws = W[W$view=="Microbiome",]
      Ws <- as.data.table(merge(Ws, tax_info, by.x = "feature", by.y = 0))
      Ws <- Ws[,value:=mean(value,na.rm=T),by=c(taxa_level, "factor")]
      Ws <- as.data.frame(Ws)
      Ws <- unique(Ws[,c("factor", "view", taxa_level, "value")])
      colnames(Ws)[3] <- "feature"
      W <- rbind(Ws[,colnames(W)], W[W$view!="Microbiome",])
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
    if(input$mofa2_top_weights_sign=="positive") W <- W[W$sign=="+",]
    if(input$mofa2_top_weights_sign=="negative") W <- W[W$sign=="-",]
    # plot
    p <- ggplot(W, aes(x=value, y=feature, fill=factor)) +
      geom_col(position="identity") +
      xlab("scaled weights") + ylab("Feature") +
      geom_vline(xintercept = 0) +
      theme(text = element_text(size=input$mofa2_weight_text_size)) +
      scale_fill_manual(values=namco_colors, name="Factor") +
      facet_wrap(~view, scales = "free_y")
    if(input$mofa2_top_weights_fix_x_axis) p <- p + xlim(c(-1,1))
    p
  }
})

# scatter of microbiome data
output$mofa2_microbiome_scatter <- renderPlot({
  if(!is.null(mofa2_data()$mofa_model)){
    model <- mofa2_data()$mofa_model
    factors <- input$mofa2_selected_impact_factors
    
    # get colors
    conditions <- unique(model@samples_metadata$condition)
    condition_cols <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(conditions))
    names(condition_cols) <- conditions
    # plot
    MOFA2::plot_data_scatter(model,
                             view = input$mofa2_microbiome_scatter_view,
                             factor = factors,
                             features = input$mofa2_microbiome_scatter_features,
                             add_lm = T,
                             color_by = "condition"
    ) + theme(text = element_text(size=input$mofa2_microbiome_scatter_text_size)) +
      ylab("Observations")
  }
})

#### DIABLO ####

run_diablo <- function() {
  waiter_show(html = tagList(spin_rotating_plane(), "Preparing DIABLO run ..."),color=overlay_color)
  # get OTU data
  otu <- as.matrix(as.data.frame(vals$datasets[[currentSet()]]$phylo@otu_table))
  expr <- as.matrix(vals$datasets[[currentSet()]]$metabolomicsExpression)
  tmp <- t(expr)
  tmp <- t(tmp[complete.cases(tmp),])
  common_samples <- sort(intersect(colnames(otu), colnames(expr)))
  # get omics data
  X <- list(
    "microbiome"=t(otu[,common_samples]),
    "metabolites"=t(expr[,common_samples])
  )
  # get meta data
  meta <- vals$datasets[[currentSet()]]$metaData
  rownames(meta) <- meta[[input$diablo_sample_label]]
  # define outcome variable, e.g. subtype etc.
  Y <- meta[common_samples,input$diablo_condition_label]
  # run either Projection to Latent Structure Discriminant Analysis (PLS-DA) or sparse PLS (sPLS-DA)
  run <- NULL
  if(input$diablo_method=="PLS"){
    run <- block.plsda(X, Y)
  } else if(input$diablo_method=="sPLS") {
    run <- block.splsda(X, Y)
  } else {
    stop("DIABLO method has to be one of PLS, sPLS")
  }
  waiter_hide()
  return(run)
}

diablo_data <- eventReactive(input$diablo_start,{
  if(vals$datasets[[currentSet()]]$has_metabolomics) {
    tryCatch({
      data <- run_diablo()
      modal <- modalDialog(
        title = p("Success!", style="color:green; font-size:40px"),
        HTML(paste0("DIABLO finished successfully")),
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

# plot diablo samples
output$diablo_samples_plot <- renderPlot({
  if(!is.null(diablo_data())){
    p <- plotIndiv(diablo_data())
    p <- p + theme(text = element_text(size=input$diablo_samples_plot_ts))
    p
  }
})

# plot diablo variables
output$diablo_var_plot <- renderPlot({
  if(!is.null(diablo_data())){
    p <- plotVar(diablo_data())
    p <- p + theme(text = element_text(size=input$diablo_samples_plot_ts))
    p
  }
})



