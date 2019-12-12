############
# store functions for different methods of calculating co-occurence
############

file = "testdata/OTUs_Table-norm.tab"
file2 = "testdata/gut_16s_abundance.txt"

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
  print(Sys.time())
  options(stringsAsFactors = F)
  #read in file (flip and change first column to column names -> this needs to be adressed differently for type of OTU-table)
  table1 <- as.data.table(t(as.matrix(read.table(file = file,
                                      sep = "\t",
                                      dec = ".",
                                      header = T,
                                      row.names = 1))))
  
  #save OTU names
  OTUs <- colnames(table1)
  sample_size <- nrow(table1)
  
  #cutoff: take 0.1 percentile of each sample and then mean over all samples to find cutoff; values below are considered to be 0
  cutoff <- mean(t(table1[, apply(table1,1,quantile,probs =.1,na.rm=TRUE)] ))
  
  #binarization & normalization
  #cutoff nicht bei 1, sondern percentile wise...
  table1[, (OTUs) := lapply(.SD,function(x){ ifelse(x>cutoff,1,0)}),.SDcols =OTUs][
    ,(OTUs) := lapply(.SD,function(x){return(x/sample_size)}),.SDcols =OTUs]
  
  result <- unique_combinations(OTUs,OTUs)
  
  counts <- unlist(lapply(seq(1,nrow(result)), function(x){
    #find 2 OTUs, for which value is calculated
    otu1 <- unlist(result[x,1])
    otu2 <- unlist(result[x,2])
    
    #skip identical OTUs
    # -> can't skip, otherwise there are not enough values in counts list in order to merge with result dt
    if(otu1==otu2){
      #result[OTU1 == otu1 & OTU2 == otu2] <- NULL
      return (0)
    }
    
    #add OTU values for 2 OTUs for each sample, where they both appear
    l<-unlist(sapply(seq(1,nrow(table1)),function(y){
      otu1_value <- table1[[otu1]][y]
      otu2_value <- table1[[otu2]][y]
      if(otu1_value > 0 & otu2_value > 0){
        s <- otu1_value + otu2_value
        return(s)
      }else{
        return (NULL)
      }
    }))
    
    total <- Reduce(`+`,l)
    #debugging stuff
    if(x%%1000==0){
      print(x)
    }
    #return count value
    return(total)

  }))

  #merge count list to data table; order is consistent
  result[,counts:=counts]
  print(Sys.time())
  return(result)
}

#create data table with all unique combinations, where (A,B) == (B,A)
unique_combinations <- function(l1, l2){
  t<-as.data.table(expand.grid(l1,l2))
  t$key<-apply(t, 1, function(x)paste(sort(x), collapse=''))
  t<-subset(t, !duplicated(t$key))
  t$key<-NULL
  names(t)<-c("OTU1","OTU2")
  return(t)
}


basic_approach_flipped(file2)


