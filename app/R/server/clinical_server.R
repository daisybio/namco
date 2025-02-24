output$selected_number <- renderText({
  paste("You selected:", input$number)
})

### Alpha Diversity ###

# Fetch JSON data from the API
fetch_json_data_alfa <- function() {
  library(httr)
  library(jsonlite)
  
  url <- "https://www.misigdb.org/api/difuture/div_metrics/"
  headers <- add_headers(
    Authorization = "Token ZNXrbp13Y7IgHPU9Qdsjmnu44D24sXpkqW6tb8KSDmcb0uX0OcwRNWke5uM2la3YH8rTHcW0z14LjVbZU2Vr1w",
    Accept = "application/json",
    `User-Agent` = "R-httr"
  )
  
  response <- GET(url, headers)
  content_type <- headers(response)[["content-type"]]
  
  if (status_code(response) != 200) {
    return(NULL)
  }
  
  if (!grepl("application/json", content_type)) {
    print("Received non-JSON response (API may be down).")
    return(NULL)
  }
  
  data <- rawToChar(response$content)
  json_data_alfa <- fromJSON(data)
  return(json_data_alfa)
}

# Reactive function to store JSON data
json_data_alfa <- reactiveVal(NULL)

# Load JSON data when the app starts
observe({
  json_data_alfa(fetch_json_data_alfa())
})

otu <- reactive({
  req(currentSet())
  phylo <- vals$datasets[[currentSet()]]$phylo
  as.matrix(phylo@otu_table)  # Convert OTU table to a matrix
})

# Observe changes and update the sample dropdown
observe({
  req(otu())  
  sample_names <- colnames(otu())[-ncol(otu())]  # Exclude taxonomy column
  updateSelectizeInput(session, "sampleSelector_Alpha", choices = sample_names, server = TRUE)
})

# Store the selected samples when the button is clicked
observeEvent(input$ClinicalReport_generateAlpha, {
  selected_samples <- input$sampleSelector_Alpha 
  req(selected_samples)  
  req(json_data_alfa())  # Ensure JSON data is available
  
  waiter_show(html = tagList(spin_rotating_plane(), "Generating Plots ... "), color = overlay_color)
  otu_table <- otu()  
  
  # Calculate richness values for selected samples (from OTU table)
  richness_values_selected <- sapply(selected_samples, function(sample) {
    sum(otu_table[, sample] > 0)  
  })
  
  # Extract all richness values from JSON data (for boxplot)
  #richness_values_json <- sapply(json_data_alfa(), function(sample) sample$richness)
  richness_values_json <- json_data_alfa()$richness
  
  # Remove NA values (to prevent errors)
  richness_values_json <- na.omit(richness_values_json)
  
  if (length(richness_values_json) == 0) {
    showNotification("No richness data available from JSON.", type = "warning")
    return()
  }
  
  # Render the boxplot when the button is clicked
  output$boxplot <- renderPlot({
    ylim <- range(c(richness_values_json, richness_values_selected), na.rm = TRUE) + c(-10, 10) 
    
    num_samples <- length(selected_samples)
    num_cols <- 7
    num_rows <- 4
    
    par(mfrow = c(num_rows, num_cols), mar = c(4, 4, 4, 1))
    
    # Create the boxplots for each selected sample
    for (i in seq_along(selected_samples)) {
      sample <- selected_samples[i]
      selected_richness <- richness_values_selected[i]
      
      # Create the boxplot using JSON richness values
      boxplot(richness_values_json,
              main = paste("Boxplot for", sample),
              ylab = "Richness",
              col = "lightblue",
              border = "darkblue",
              ylim = ylim)  # Ensure the y-axis includes the point
      
      # Add the selected sample's richness value as a red dot
      points(1, selected_richness, col = "red", pch = 19, cex = 1.5)
    }
  })
  waiter_hide()
})


### Stacked Bar Charts ###

# Observe event for sample selection
observe({
  req(otu(), ncol(otu()) > 1) 
  sample_names <- colnames(otu())[-ncol(otu())]  # Exclude last column (metadata)
  updateSelectizeInput(session, "sampleSelector_Krona", choices = sample_names, server = TRUE)
})

# Event to generate stacked bar data when the button is clicked
stacked_bar_data <- eventReactive(input$ClinicalReport_generateKrona, {
  req(currentSet())  # Ensure dataset is selected
  req(input$sampleSelector_Krona)  # Ensure sample selector is used
  req(input$ClinicalReport_phylolevel)  # Ensure phylogenetic level is selected
  waiter_show(html = tagList(spin_rotating_plane(), "Generating plot ... "), color = overlay_color)
  
  # Prune most abundant taxa
  myTaxa <- names(sort(taxa_sums(vals$datasets[[currentSet()]]$phylo), decreasing = TRUE)[1:input$phylo_prune])
  phy <- prune_taxa(myTaxa, vals$datasets[[currentSet()]]$phylo)
  
  # Extract OTU table and taxonomy data
  otu <- vals$datasets[[currentSet()]]$normalizedData[myTaxa, input$sampleSelector_Krona]  # Subset based on selected samples
  taxonomy <- as.data.frame(tax_table(phy))
  
  # Convert row names to a proper column
  taxonomy$OTU <- rownames(taxonomy)
  rownames(taxonomy) <- NULL
  
  # Check if the "Class" column exists and use it directly
  if (!"Class" %in% colnames(taxonomy)) {
    taxonomy$Class <- "Unknown"
  }
  
  # Merge OTU data with taxonomy
  otu_table <- as.data.frame(t(otu))  # Transpose OTU table so samples are rows
  otu_table$Sample <- rownames(otu_table)
  otu_table <- pivot_longer(otu_table, cols = -Sample, names_to = "OTU", values_to = "Abundance")
  
  # Merge taxonomy and OTU table by "OTU"
  otu_table <- left_join(otu_table, taxonomy, by = "OTU")
  
  # Get the selected taxonomic level (to group by it)
  tax_selected_level <- input$ClinicalReport_phylolevel  # Phylogenetic level selected by user
  
  # Summarize data by Sample and selected taxonomic level to calculate relative abundance
  otu_summarized <- otu_table %>%
    group_by(Sample, !!sym(tax_selected_level)) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(Relative_Abundance = Abundance / sum(Abundance))  # Calculate relative abundance
  waiter_hide()
  return(otu_summarized)
})

# Phylo Level
taxBinningReact_ClinicalReport <- reactive({
  req(currentSet())  # Ensure dataset is selected
  phylo <- vals$datasets[[currentSet()]]$phylo
  rel_dat <- vals$datasets[[currentSet()]]$relativeData
  
  # Get the selected taxonomic level from UI input
  tax_selected_level <- input$ClinicalReport_phylolevel  # Updated input name
  
  # Ensure the level is valid (non-empty and one of the taxonomic levels)
  if (is.null(tax_selected_level) || tax_selected_level == "") {
    return(NULL)  # Return NULL if no level is selected or it's an invalid level
  }
  
  waiter_show(html = tagList(spin_rotating_plane(), "Generating data ... "), color = overlay_color)
  
  # Create the phyloseq object with the relative abundance data
  rel_phylo <- merge_phyloseq(otu_table(rel_dat, TRUE), tax_table(phylo))
  
  # Perform the taxonomic binning based on the selected level
  tax_binning_data <- taxBinningNew(rel_phylo, vals$datasets[[currentSet()]]$is_fastq, tax_selected_level)
  
  # Adjust the taxonomic binning level (e.g., Kingdom, Phylum, etc.)
  if (vals$datasets[[currentSet()]]$is_fastq) {
    tax_binning_table <- tax_binning_data[[which(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus") == tax_selected_level)]]
  } else {
    tax_binning_table <- tax_binning_data[[which(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") == tax_selected_level)]]
  }
  
  # Processing: Top N taxa, Metadata merging
  if (input$taxBinningTop < nrow(tax_binning_table) && !is.na(input$taxBinningTop)) {
    top_taxa <- names(sort(rowSums(tax_binning_table[-1]), decreasing = TRUE)[1:input$taxBinningTop])
    taxa <- ifelse(rownames(tax_binning_table) %in% top_taxa, rownames(tax_binning_table), "Other")
    other <- data.frame(tax_binning_table[which(taxa == "Other"),])
    other <- colSums(other[-1])
    other <- c("Other", other)
    
    tax_binning_table <- tax_binning_table[-which(taxa == "Other"), ]
    tax_binning_table <- rbind(tax_binning_table, Other = other)
  } else {
    top_taxa <- rownames(tax_binning_table)
  }
  
  # Merge with metadata (if present)
  if (vals$datasets[[currentSet()]]$has_meta) {
    meta <- data.frame(sample_data(vals$datasets[[currentSet()]]$phylo), check.names = FALSE)
    final_data <- merge(melt(tax_binning_table, id.vars = tax_selected_level), meta, by.x = "variable", by.y = sample_column, all.x = TRUE)
  } else {
    final_data <- melt(tax_binning_table, id.vars = tax_selected_level)
  }
  
  # Final data preparation (renaming columns to avoid conflicts)
  colnames(final_data)[which(colnames(final_data) == "variable")] <- sample_column
  colnames(final_data)[which(colnames(final_data) == tax_selected_level)] <- "custom_taxonomy_column"
  colnames(final_data)[which(colnames(final_data) == input$taxBinningYLabel)] <- "y_split"
  colnames(final_data)[which(colnames(final_data) == input$taxBinningGroup)] <- "facet_split"
  final_data$value <- as.numeric(final_data$value)
  
  if ('y_split' %in% colnames(final_data)) final_data$y_split <- as.character(final_data$y_split)
  if ('facet_split' %in% colnames(final_data)) final_data$facet_split <- as.character(final_data$facet_split)
  
  waiter_hide()
  return(final_data)
})

# Render Stacked Bar Plot
output$stacked_barplot <- renderPlot({
  req(stacked_bar_data())  # Ensure the stacked bar data is available
  
  # Extract the final data
  plot_data <- stacked_bar_data()
  
  # Ensure the plot_data has the necessary columns
  if (nrow(plot_data) == 0) {
    return(NULL)  # Return nothing if no data
  }
  
  # Get the selected phylogenetic level from input
  tax_selected_level <- input$ClinicalReport_phylolevel
  
  # Dynamically set the fill column based on the selected taxonomic level
  fill_column <- tax_selected_level
  
  # Check if the required column exists
  if (!fill_column %in% colnames(plot_data)) {
    return(NULL)  # Return nothing if the required column is missing
  }
  
  # Plot the stacked bar plot (horizontal)
  ggplot(plot_data, aes(x = Relative_Abundance, y = Sample, fill = !!sym(fill_column))) +
    geom_bar(stat = "identity", position = "stack") +  # Stack bars
    scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(plot_data[[fill_column]])))) +  # Custom colors
    coord_flip() +  # Make the bars horizontal
    labs(y = "Samples", x = "Relative Abundance", fill = tax_selected_level) +  # Use tax_selected_level here
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed
})


### BETA DIVERSITY ###

library(httr)
library(jsonlite)
library(vegan)
library(shiny)
library(DT)
library(ggplot2)
library(gridExtra)

# --- Standardization Function ---

# Function to remove taxonomic prefixes (for uploaded data)
remove_prefixes <- function(taxonomy) {
  gsub("(k__|p__|c__|o__|f__|g__|s__)", "", taxonomy)
}

# --- Processing Functions ---

# Process the uploaded OTU table
process_uploaded_otu <- function(phylo) {
  otu_data <- as.data.frame(as.matrix(phylo@otu_table))
  taxonomy_matrix <- phylo@tax_table@.Data
  taxonomy <- apply(taxonomy_matrix, 1, function(x) {
    x[x == ""] <- ""
    tax <- paste(remove_prefixes(x), collapse = ";")
    return(tax)
  })
  otu_data$taxonomy <- taxonomy
  otu_data <- aggregate(. ~ taxonomy, data = otu_data, FUN = sum)
  rownames(otu_data) <- otu_data$taxonomy
  otu_data <- otu_data[, -1, drop = FALSE]
  otu_data[otu_data > 0] <- 1  # Convert counts to binary
  return(otu_data)
}

# Process JSON data to build a presence/absence table
process_json_pa <- function(json_data) {
  if (is.null(names(json_data)) || any(names(json_data) == "")) {
    names(json_data) <- paste0("JSON_Sample_", seq_along(json_data))
  }
  
  json_std <- lapply(json_data, function(sample) {
    sample[sapply(sample, function(x) grepl("^Bacteria;", x))]
  })
  
  otu_vec <- unlist(json_std)
  sample_ids <- rep(names(json_data), times = sapply(json_std, length))
  
  pa_table <- table(otu_vec, sample_ids)
  pa_table <- (pa_table > 0) * 1
  return(as.data.frame.matrix(pa_table))
}

# Fetch JSON Data Function
fetch_json_data_beta <- function() {
  url <- "https://www.misigdb.org/api/beta_diversity/taxas/"
  headers <- add_headers(
    Authorization = "Token ZNXrbp13Y7IgHPU9Qdsjmnu44D24sXpkqW6tb8KSDmcb0uX0OcwRNWke5uM2la3YH8rTHcW0z14LjVbZU2Vr1w",
    Accept = "application/json",
    `User-Agent` = "R-httr"
  )
  response <- GET(url, headers)
  if (status_code(response) == 200) {
    json_data <- fromJSON(rawToChar(response$content))
    names(json_data) <- paste0("JSON_Sample_", seq_along(json_data))
    return(json_data)
  } else {
    stop("Failed to fetch JSON data from API.")
  }
}

# Final Process: Build combined presence/absence table and compute Jaccard dissimilarity
final_beta_diversity <- reactive({
  req(currentSet())
  waiter_show(html = tagList(spin_rotating_plane(),"Calculating distance ... "),color=overlay_color)
  phylo <- vals$datasets[[currentSet()]]$phylo  
  uploaded_pa <- process_uploaded_otu(phylo)
  json_data <- fetch_json_data_beta()
  json_pa <- process_json_pa(json_data)
  
  # Select only the first 10% of JSON samples
  json_sample_count <- ncol(json_pa)
  selected_json_samples <- colnames(json_pa)[seq_len(ceiling(json_sample_count * 0.1))]
  json_pa <- json_pa[, selected_json_samples, drop = FALSE]
  
  all_otus <- union(rownames(uploaded_pa), rownames(json_pa))
  
  missing_uploaded <- setdiff(all_otus, rownames(uploaded_pa))
  if (length(missing_uploaded) > 0) {
    add_mat <- matrix(0, nrow = length(missing_uploaded), ncol = ncol(uploaded_pa),
                      dimnames = list(missing_uploaded, colnames(uploaded_pa)))
    uploaded_pa <- rbind(uploaded_pa, add_mat)
  }
  
  missing_json <- setdiff(all_otus, rownames(json_pa))
  if (length(missing_json) > 0) {
    add_mat <- matrix(0, nrow = length(missing_json), ncol = ncol(json_pa),
                      dimnames = list(missing_json, colnames(json_pa)))
    json_pa <- rbind(json_pa, add_mat)
  }
  
  uploaded_pa <- uploaded_pa[all_otus, , drop = FALSE]
  json_pa <- json_pa[all_otus, , drop = FALSE]
  
  combined_pa <- cbind(uploaded_pa, json_pa)
  jaccard_matrix <- as.matrix(vegdist(t(combined_pa), method = "jaccard"))
  waiter_hide()
  return(jaccard_matrix)
})

# NMDS Plot
output$nmds_plot <- renderPlot({
  req(final_beta_diversity())
  dist_matrix <- final_beta_diversity()
  
  nmds_result <- metaMDS(dist_matrix, trymax = 50)  
  sample_names <- rownames(dist_matrix)
  
  json_samples <- grep("^JSON_Sample_", sample_names, value = TRUE)
  uploaded_samples <- setdiff(sample_names, json_samples)
  
  plot_data <- data.frame(
    Dimension1 = nmds_result$points[, 1],
    Dimension2 = nmds_result$points[, 2],
    SampleType = ifelse(rownames(nmds_result$points) %in% json_samples, "JSON Data", "OTU Table"),
    Label = ifelse(rownames(nmds_result$points) %in% uploaded_samples, rownames(nmds_result$points), NA)
  )
  
  ggplot(plot_data, aes(x = Dimension1, y = Dimension2, color = SampleType)) +
    geom_point(size = 2) +
    geom_text(aes(label = Label), vjust = -1, hjust = 0.5, color = "red", size = 3, na.rm = TRUE) +
    scale_color_manual(values = c("OTU Table" = "red", "JSON Data" = "blue")) +
    labs(title = "NMDS Plot", x = "NMDS 1", y = "NMDS 2", color = "Sample Type") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.justification = "center",
          legend.box.margin = margin(0, 15, 0, 0))
})

# MDS Plot (Classical MDS)
output$mds_plot <- renderPlot({
  req(final_beta_diversity())
  dist_matrix <- final_beta_diversity()
  
  mds_result <- cmdscale(dist_matrix, k = 2)
  sample_names <- rownames(dist_matrix)
  
  json_samples <- grep("^JSON_Sample_", sample_names, value = TRUE)
  uploaded_samples <- setdiff(sample_names, json_samples)
  
  plot_data <- data.frame(
    Dimension1 = mds_result[, 1],
    Dimension2 = mds_result[, 2],
    SampleType = ifelse(rownames(mds_result) %in% json_samples, "JSON Data", "OTU Table"),
    Label = ifelse(rownames(mds_result) %in% uploaded_samples, rownames(mds_result), NA)
  )
  
  ggplot(plot_data, aes(x = Dimension1, y = Dimension2, color = SampleType)) +
    geom_point(size = 2) +
    geom_text(aes(label = Label), vjust = -1, hjust = 0.5, color = "red", size = 3, na.rm = TRUE) +
    scale_color_manual(values = c("OTU Table" = "red", "JSON Data" = "blue")) +
    labs(title = "Classical MDS Plot", x = "Dimension 1", y = "Dimension 2", color = "Sample Type") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.justification = "center",
          legend.box.margin = margin(0, 15, 0, 0))
})

## UI Output: Display the Jaccard Distance Matrix Table
#output$beta_diversity_table <- DT::renderDataTable({
#  req(final_beta_diversity())
#  round(as.data.frame(final_beta_diversity()), 2) 
#}, options = list(pageLength = 50, scrollX = TRUE))
