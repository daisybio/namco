### PARAMETERS
omics = c("Metabolomics", "Proteomics")

observe({
  updateSelectInput(session, "omicsSelection", choices = omics)
})

# read expression
read_omics_expr <- function(file){
  # read expression file, can be read as an OTU table
  expr <- read_csv_custom(file, file_type = "otu")
  # check for 15 or more samples
  if(ncol(expr)<15) stop("Multi-omics expression table has to contain a minimum of 15 samples.")
  # check common sample IDs
  otu_samples <- colnames(vals$datasets[[currentSet()]]$phylo@otu_table)
  expr_samples <- colnames(expr)
  cs <- intersect(otu_samples, expr_samples)
  if(length(cs)<15) stop("Less than 15 shared samples between OTU table and expression found.")
  return(expr[,cs])
}

### upload omics
observeEvent(input$upload_omics,{
  if(is.null(currentSet())) {
    showModal(errorModal(error_message = "Please first create a dataset by uploading data to add the metabolomics expression to."))
  }
  else if(!is.null(input$omicsExpressionFile$datapath)){
    tryCatch({
      omicsExpression <- read_omics_expr(input$omicsExpressionFile$datapath[1])
      # assign expression to slot in current dataset
      vals$datasets[[currentSet()]]$has_omics <- T
      if(!"omics"%in%names(vals$datasets[[currentSet()]])) vals$datasets[[currentSet()]]$omics = list()
      vals$datasets[[currentSet()]]$omics[[input$omicsSelection]] = omicsExpression
      message(paste0(Sys.time(), pasteß(" - Successfully loaded ", input$omicsSelection, " expression data with dimensions: ", dim(metabolomicsExpression))))
    },
    error=function(e){
      message(e$message)
      showModal(errorModal(error_message = e$message))
    })
  } else {
    showModal(errorModal(error_message = "Please first upload a metabolomics expression file to continue."))
  }
})

