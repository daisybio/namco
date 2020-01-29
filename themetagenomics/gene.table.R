gene.table <- function(functions_obj,level=2){
  fxn_table <- functions_obj$fxn_table
  fxn_meta <- functions_obj$fxn_meta
  
  fxn_table <- fxn_table[,colSums(fxn_table) > 0]
  fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])
  
  functions_obj$fxn_table <- fxn_table
  functions_obj$fxn_meta <- fxn_meta
  
  gene_table <- list(gene_table=format_gene_table(functions_obj,level=level))
  
  return(gene_table)
}