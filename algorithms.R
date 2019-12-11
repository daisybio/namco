############
# store functions for different methods of calculating co-occurence
############

#file = "Documents/Studium_Bioinformatik/Master/Semester1/SystemsBioMedicine/Project/namco/namco/OTUs_Table-norm.tab"

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
  
  #binarization & normalization
  table1[, (OTUs) := lapply(.SD,function(x){ ifelse(x>0,1,0)}),.SDcols =OTUs][
    ,(OTUs) := lapply(.SD,function(x){return(x/sample_size)}),.SDcols =OTUs]
  
  sums <- colSums(table1)
  
  m <- as.matrix(table1)
  result <- as.data.table(expand.grid(OTU1=OTUs,OTU2=OTUs))
  counts <- list() 
  
  counts <- list(lapply(seq(1,nrow(result)), function(x){
    otu1 <- result$OTU1[x]
    otu2 <- result$OTU2[x]
    value <- sums[otu1]+sums[otu2]
    return(value)
  }))
  
  
  # for(i in seq(1,nrow(result))){
  #   
  #   otu1 <- result$OTU1[i]
  #   otu2 <- result$OTU2[i]
  #   value <- sums[otu1]+sums[otu2]
  #   #print(value)
  #   counts<-append(counts,value)
  #   
  # }
  
  result[,counts:=counts]
}






