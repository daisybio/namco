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

output$multi_omics_overview <- renderUI({
  if(vals$datasets[[currentSet()]]$has_omics) {
    m <- length(vals$datasets[[currentSet()]]$omics)
    w <- round(12/m)
    boxes <- lapply(names(vals$datasets[[currentSet()]]$omics), function(n) {
      valueBox(
        subtitle = n,
        value = nrow(vals$datasets[[currentSet()]]$omics[[n]]),
        color = "green",
        width = w
      )
    })
    outUI <- do.call(tagList, boxes)
  } else {
    outUI <- valueBox(
      subtitle = "Upload to continue",
      value = 0,
      color = "green"
    )
  }
  return(outUI)
})

#### MOFA2 ####
tidy_expr <- function(mat) {
  m <- as.matrix(mat)
  return(m[,!is.na(colSums(m))])
}

# filter OTUs based on linear regression
regression_filter <- function(phylo, label_col){
  # prepare data
  x = t(as.data.frame(phylo@otu_table))
  y = as.matrix(data.frame(phylo@sam_data[,label_col]))
  # convert to factors
  y = apply(y, 2, function(j){
    if(class(j)=="character") {
      as.double(factor(j))
    } else {
      as.double(j)
    }
  })
  # outline numbers
  n_observations <- nrow(x)
  n_predictors <- ncol(x)
  real_p <- ncol(y)
  # fit optimal model
  alphas = seq(0, 1, 0.1)
  grid <- glmnetUtils::cva.glmnet(x, y,
                                  type.measure="mse",
                                  alpha=alphas,
                                  family="gaussian")
  # get best run
  best_alpha_idx = which.min(sapply(grid$modlist, "[[", "lambda.min"))
  best_alpha = alphas[[best_alpha_idx]]
  model = grid$modlist[[best_alpha_idx]]
  b_l = model$lambda.min
  b_m = model$glmnet.fit
  # extract coefficients
  coefficients <- coef(b_m, s = b_l)
  # Convert coefficients to a data frame
  coefficients_df <- as.data.frame(as.matrix(coefficients))
  coefficients_df$feature <- rownames(coefficients_df)
  # sort and exclude 0s
  coefficients_df <- coefficients_df[order(-abs(coefficients_df$s1)),]
  coefficients_df <- coefficients_df[coefficients_df$s1 != 0,]
  # extract the important features
  important_features <- coefficients_df[coefficients_df$feature!=("(Intercept)"),"feature"]
  return(list(coefs=important_features, data=coefficients_df, lambda = b_l, alpha = best_alpha))
}

run_namco_mofa2 <- function(){
  waiter_show(html = tagList(spin_rotating_plane(), "Preparing MOFA2 data ..."),color=overlay_color)
  # align samples of otu table
  otu <- as.matrix(as.data.frame(vals$datasets[[currentSet()]]$phylo@otu_table))
  omics <- lapply(vals$datasets[[currentSet()]]$omics, tidy_expr)
  meta <- vals$datasets[[currentSet()]]$metaData
  # normalize omics if argument is given
  if(input$mofa2_norm_omics){
    omics <- lapply(omics, normalizeExpression, meta, input$mofa2_condition_label)
  }
  # apply taxonomic regression filter if specified
  if(input$mofa2_regression_filter) {
    waiter_update(html = tagList(spin_rotating_plane(), "Applying regression filter ..."))
    reg_filter = regression_filter(phylo = vals$datasets[[currentSet()]]$phylo, label_col = input$mofa2_condition_label)
    otu = otu[reg_filter$coefs,]
  } else if(input$mofa2_taxa_reduce_level!="Do not reduce"){
    waiter_update(html = tagList(spin_rotating_plane(), paste0("Collapsing to ", input$mofa2_taxa_reduce_level, " level ...")))
    tax = glom_taxa_custom(vals$datasets[[currentSet()]]$phylo, rank = input$mofa2_taxa_reduce_level)
    otu = as.matrix(as.data.frame(tax$phylo_rank@otu_table))
  }
  waiter_update(html = tagList(spin_rotating_plane(), "Still preparing MOFA2 data ..."))
  # check sample consistency
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
      modal_text = "MOFA2 finished successfully"
      if(length(warnings())>0) {
        # modal_text = paste0(modal_text)
      }
      modal <- modalDialog(
        title = p("Success!", style="color:green; font-size:40px"),
        HTML(modal_text),
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
    updateSelectInput(session, "mofa2_selected_impact_factors", choices = 1:m, selected = 1)
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
    # specify color scheme
    conditions <- unique(m@samples_metadata$condition)
    condition_cols <- colorRampPalette(brewer.pal(9, input$namco_pallete))(length(conditions))
    names(condition_cols) <- conditions
    # build an overview for the omic dimensions
    overview <- melt(m@dimensions$D)
    overview$Omic <- rownames(overview)
    meta <- data.frame(m@samples_metadata %>% select(-group), row.names = 1)
    # plot as bar plot
    bars <- ggplot(overview, aes(x=Omic, y=value, fill=Omic)) +
      geom_bar(stat="identity") +
      ylab("Number of features") + xlab("") +
      theme(legend.position="none", text = element_text(size=16)) +
      scale_fill_brewer(palette = input$namco_pallete) +
      ggtitle("")
    # show shared samples and conditions in pie plots
    pieData <- melt(table(meta$condition))
    colnames(pieData)[1] <- "Condition"
    # pieData <- rbind(pieData, data.frame(Var1="samples not shared", value=2))
    pie <- ggplot(pieData, aes(x="", y=value, fill=Condition)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() +
      geom_text(aes(y = value/nrow(pieData) + c(0, cumsum(value)[-length(value)]), label = value), size = 8) +
      theme(text = element_text(size = 16)) +
      scale_fill_manual(values = condition_cols)
    
    # plot next to each other
    ggarrange(bars, pie, 
              labels = c("Overview of number of features for each dataset", "Shared samples grouped by condition"),
              nrow = 1, ncol = 2)
  }
})

# explained variance plot
output$mofa2_explained_variance_group <- renderPlot({
  if(!is.null(mofa2_data()$explained_variance_plot)) {
    e <- mofa2_data()$explained_variance_plot
    d <- e$data
    d$color <- ifelse(d$value > max(d$value)/2, "white", "black")
    fs <- sapply(strsplit(as.character(d$factor), "Factor"), paste, collapse = "Factor ")
    d$factor <- factor(fs, levels = unique(fs))
    # plot variance
    ggplot(d, aes(x = view, y = factor, fill = value)) +
      geom_tile() +
      geom_text(aes(label = paste0(round(value, 2), "%")), color = d$color) +
      scale_fill_gradient(low = "white", high = "blue", name = "Var. (%)") +
      theme_minimal() +
      labs(title = "", x = "", y = "") +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size = input$mofa2_variance_text_size)
      )
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
output$mofa2_factor_plot <- renderPlotly({
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
    ggplotly(p) %>%
      style(
        text = paste("Sample ID:", p$data$sample, "<br>Value:", round(p$data$value, 2)),
        hoverinfo = "text"
      )
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
    p
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
    factors <- paste0("Factor", input$mofa2_selected_impact_factors)
    
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

### normalizing functions ###
normalizeExpression <- function(data, meta, design_col, force=T, pc=1) {
  # check for raw counts
  if (round(data[1,1])!=data[1,1]) {
    warning("Data was not raw counts (not integers)")
    if (force) {
      warning("Forcing normalization on non-integer counts (rounding)")
      data <- round(data)
    }
  }
  # check for 0s, add pseudocount of pc
  if (any(data==0)) {
    data = data + pc
  }
  ## Create DESeq2Dataset object
  meta$design <- meta[,design_col]
  m <- meta[colnames(data),]
  dds <- DESeqDataSetFromMatrix(countData = data, colData = m, design = ~ design)
  # normalize
  dds <- DESeq(dds)
  # extract counts
  normalized_counts <- counts(dds, normalized=T)
  # Variance Stabilizing Transformation
  vst_data <- varianceStabilizingTransformation(dds, blind = FALSE)
  return(assay(vst_data))
}

output$mofa2_download_results <- downloadHandler(
  filename = function() {
    paste("mofa2_results", ".rds", sep = "")
  },
  content = function(file) {
    if(!is.null(mofa2_data())){
      tempFile <- tempfile(fileext = ".rds")
      saveRDS(diablo_data(), file = tempFile)
      # Copy the RDS file to the file path provided by Shiny
      file.copy(tempFile, file)
    }
  },
  contentType = "application/octet-stream"
)

#### DIABLO ####

run_diablo <- function(otu) {
  waiter_show(html = tagList(spin_rotating_plane(), "Preparing DIABLO run ..."),color=overlay_color)
  # get OTU data
  otu <- as.matrix(as.data.frame(otu))
  omics <- lapply(vals$datasets[[currentSet()]]$omics, tidy_expr)
  # get meta data
  meta <- vals$datasets[[currentSet()]]$metaData
  rownames(meta) <- meta[[input$diablo_sample_label]]
  
  if(input$diablo_norm_omics){
    waiter_update(html = tagList(spin_rotating_plane(), "normalizing omics ..."))
    omics <- lapply(omics, normalizeExpression, meta, input$diablo_condition_label)
    waiter_update(html = tagList(spin_rotating_plane(), "Preparing DIABLO run ..."))
  }
  # check sample consistency
  omics_samples <- unique(unlist(lapply(omics, colnames)))
  common_samples <- sort(intersect(colnames(otu), omics_samples))
  
  # build matrix list
  matrix_list <- lapply(omics, function(expr) as.matrix(expr[,common_samples]))
  matrix_list[["Microbiome"]] <- otu[,common_samples]
  # get omics data
  X <- lapply(matrix_list, t)
  # define outcome variable, e.g. subtype etc.
  Y <- as.factor(meta[common_samples,input$diablo_condition_label])
  if (input$diablo_det_opt_ncomp) {
    run <- block.plsda(X, Y, ncomp=min(sapply(X, dim)[2,]))
    ncomp <- sum(run$prop_expl_var$microbiome > .1)
  } else {
    ncomp <- 2
  }
  # run either Projection to Latent Structure Discriminant Analysis (PLS-DA) or sparse PLS (sPLS-DA)
  waiter_update(html = tagList(spin_rotating_plane(), "Running DIABLO ..."))
  run <- NULL
  if(input$diablo_method=="PLS"){
    run <- block.plsda(X, Y, ncomp = ncomp)
  } else if(input$diablo_method=="sPLS") {
    run <- block.splsda(X, Y, ncomp = ncomp)
  } else {
    stop("DIABLO method has to be one of PLS, sPLS")
  }
  waiter_hide()
  return(run)
}

diablo_data <- eventReactive(input$diablo_start,{
  if(vals$datasets[[currentSet()]]$has_omics) {
    tryCatch({
      data <- run_diablo(vals$datasets[[currentSet()]]$phylo@otu_table)
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

# plotDiablo() method from mixOmics
output$diablo_plot <- renderPlot({
  if(!is.null(diablo_data())){
    p <- mixOmics::plotDiablo(diablo_data())
    p <- p + theme(text = element_text(size=input$diablo_samples_plot_ts))
    p
  }
})

agg_func <- function(x, func) {
  if(func=="mean") return(mean(x))
  else if(func=="median") return(median(x))
  else if(func=="sum") return(sum(x))
}

# plot loadings of each block/omic
output$diablo_loadings <- renderPlot({
  if(!is.null(diablo_data())){
    run <- diablo_data()
    d <- run$loadings
    # collapse to taxa level if specified
    if(input$diablo_taxa_lvl!="OTU/ASV") {
      phylo <- vals$datasets[[currentSet()]]$phylo
      tax_info <- as.data.frame(phylo@tax_table)
      tax_lvl <- input$diablo_taxa_lvl
      col_tax <- merge(d$Microbiome, tax_info %>% select(all_of(tax_lvl)), by=0)
      col_tax <- col_tax %>% group_by_at(tax_lvl) %>%
        summarise(comp1=agg_func(comp1, input$diablo_taxa_agg),
                  comp2=agg_func(comp2, input$diablo_taxa_agg))
      d$Microbiome <- as.matrix(data.frame(col_tax, row.names = 1))
    }
    d <- melt(d)
    # remove predictor variable
    d <- d[d$L1!="Y",]
    colnames(d) <- c("Feature", "Component", "weight", "omic")
    # choose which Components to show
    d <- d[d$Component%in%paste0("comp",input$diablo_comp_select),]
    # sort by importance over all components
    imp <- d %>% group_by(Feature) %>% summarise(importance=sum(abs(weight)))
    d <- merge(d, imp)
    d %>% arrange(importance) %>%
      mutate(Feature=factor(Feature, levels=unique(Feature))) %>%
      ggplot(aes(x=weight, y=Feature, fill=Component)) +
      geom_bar(stat="identity", position = "dodge") +
      scale_fill_brewer(palette = input$namco_pallete) +
      facet_wrap(~omic, scales = "free") +
      theme(text = element_text(size=input$diablo_samples_plot_ts))
  }
})

# download diablo results
output$diablo_download_results <- downloadHandler(
  filename = function() {
    paste("diablo_results", ".rds", sep = "")
  },
  content = function(file) {
    if(!is.null(diablo_data())){
      tempFile <- tempfile(fileext = ".rds")
      saveRDS(diablo_data(), file = tempFile)
      # Copy the RDS file to the file path provided by Shiny
      file.copy(tempFile, file)
    }
  },
  contentType = "application/octet-stream"
)

# plot diablo samples
output$diablo_samples_arrow_plot <- renderPlot({
  if(!is.null(diablo_data())){
    if(length(input$diablo_arrow_comp)==2) {
      p <- plotArrow(diablo_data(), 
                     comp = as.numeric(input$diablo_arrow_comp),
                     arrow.size = input$diablo_arrow_size,
                     pch.size = input$diablo_arrow_pch, 
                     ind.names.size = input$diablo_arrow_label_size)
    } else {
      p <- ggplot() +
        theme_void() +
        annotate("text", 
                 x = 0.5, y = 0.5,
                 label = "Please provide two components to plot", 
                 size = 8,
                 color = "grey")
    }
    p
  }
})

agg_diablo_run <- reactiveValues(cache=list(run=NULL,lvl=NULL))

# plot diablo variables
output$diablo_circos <- renderPlot({
  if(!is.null(diablo_data())){
    run <- NULL
    if (input$diablo_cir_taxa_lvl!="OTU/ASV") {
      # if new run, cache it
      if (is.null(agg_diablo_run$cache$lvl) || input$diablo_cir_taxa_lvl!=agg_diablo_run$cache$lvl) {
        waiter_update(html = tagList(spin_rotating_plane(), "Aggregating taxa ..."))
        # Aggregate taxa and run DIABLO
        tax = glom_taxa_custom(vals$datasets[[currentSet()]]$phylo, rank = input$diablo_cir_taxa_lvl)
        otu = as.matrix(as.data.frame(tax$phylo_rank@otu_table))
        run = run_diablo(otu)
        # update aggregation run
        agg_diablo_run$cache <- list(run=run,lvl=input$diablo_cir_taxa_lvl)
      } else {
        # recover aggregation run if same
        run <- agg_diablo_run$cache$run
      }
    } else {
      run <- diablo_data()
    }
    # plot circos
    circosPlot(run, cutoff=input$diablo_cir_cutoff,
               group=run$Y,
               line=input$diablo_cir_line, 
               var.adj = input$diablo_cir_var_adj, 
               blocks.labels.adj = input$diablo_cir_lab_adj,
               size.labels = input$diablo_cir_label_ps, 
               size.variables = input$diablo_cir_var_ps)
  }
})

observe({
  if(input$diablo_taxa_lvl!="OTU/ASV"){
    enable("diablo_taxa_agg")
  }
})

observe({
  if(!is.null(diablo_data())){
    shinyjs::show("diablo_results_div")
    comps <- 1:max(diablo_data()$ncomp)
    updateSelectInput(session, "diablo_arrow_comp", choices = comps, selected = comps)
    updateSelectInput(session, "diablo_comp_select", choices = comps, selected = comps)
    updateSelectInput(session, "diablo_cir_comp_select", choices = comps, selected = comps)
  }
})


