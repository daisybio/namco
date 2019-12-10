############
# store functions for different methods of calculating co-occurence
############

library(data.table)
library(tidyr)

basic_approach <- function(file){
  #TODO: test file to be a valid OTU table
  
  
  table <- fread(file = file,
                 sep = "\t",
                 dec = ".")
  
  #getting table into right format (row == sample, column == OTU)
    #as.data.table(t(as.matrix(read.table(file = file,
    #                      sep = "\t",
    #                      dec = ".",
    #                      header = T,
    #                      row.names = 1))))
  
  #save OTU names
  OTUs <- as.list(table[,1])
  table[,V1 := NULL]  
  
  
  cols <- colnames(table)
  sample_size <- nrow(table)
  
  #binarization & normalization
  table[,(cols) := lapply(.SD,function(x){ ifelse(x > 0, 1, 0)}),.SDcols = cols][
    ,(cols) := lapply(.SD,function(x) {return (x*100/sample_size)}),.SDcols = cols]
  
  
}


basic_approach_flipped <- function(file){
  table1 <- as.data.table(t(as.matrix(read.table(file = file,
                                      sep = "\t",
                                      dec = ".",
                                      header = T,
                                      row.names = 1))))
  
  OTUs <- colnames(table1)
  sample_size <- ncol(table1)
  
  table1[, (OTUs) := lapply(.SD,function(x){ ifelse(x>0,1,0)}),.SDcols =OTUs][
    ,(OTUs) := lapply(.SD,function(x){return(x/sample_size)}),.SDcols =OTUs]
  
  sums <- colSums(table1)
  
  m <- as.matrix(table1)
  result <- as.data.table(expand.grid(OTUs,OTUs))
  counts <- list() 
  
  j=1
  for(i in seq(1,nrow(result))){
    
    otu1 <- result[i,1]
    otu2 <- result[i,2]
    value <- sums[otu1]+sums[otu2]
    print(value)
    append(counts,value)
    
    if(j==130){
      j=1
    }else{
      j=j+1
    }
  }
}






