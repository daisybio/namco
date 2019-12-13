############
# store functions for different methods of calculating co-occurence
############

file = "testdata/OTUs_Table-norm.tab"
file2 = "testdata/gut_16s_abundance.txt"

library(data.table)
library(tidyr)
library(reshape2)

basic_approach <- function(file){
  print(Sys.time())
  table1 <- as.data.table(t(as.matrix(read.table(file = file,
                                                 sep = "\t",
                                                 dec = ".",
                                                 header = T,
                                                 row.names = 1))))
  
  #save OTU names
  OTUs <- colnames(table1)
  n_otus <- length(OTUs)
  sample_size <- nrow(table1)
  
  #cutoff: take 0.1 percentile of each sample and then mean over all samples to find cutoff; values below are considered to be 0
  cutoff <- mean(t(table1[, apply(table1,1,quantile,probs =.1,na.rm=TRUE)] ))
  
  #binarization & normalization
  #cutoff nicht bei 1, sondern percentile wise...
  table1[, (OTUs) := lapply(.SD,function(x){ ifelse(x>cutoff,1,0)}),.SDcols =OTUs][
    ,(OTUs) := lapply(.SD,function(x){return(x/sample_size)}),.SDcols =OTUs]
  
  
  mat <- matrix(,nrow=n_otus,ncol=n_otus)
  table1<-as.data.frame(table1)
  
  for(i in 1:n_otus){
    rowsToCount <- which(table1[,i] > 0)
    for(j in 1:n_otus){
      #this skips upper triangle
      if(i<j){
        next()
      }else{
        mat[i,j]<-sum(table1[rowsToCount,j])
      }
    }
    #set diagonale to NA
    mat[i,i]<-NA
  }
  rownames(mat)<-OTUs
  colnames(mat)<-OTUs
  
  counts<-setNames(melt(mat,na.rm = T),c("OTU1","OTU2","value"))
  print(Sys.time())
  return(counts)
}

####deprecated#####
#DO NOT USE!!
#..only if you have a lot of time 
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
        #s <- otu1_value + otu2_value
        return(1)
      }else{
        return (0)
      }
    }))
    
    total <- Reduce(`+`,l)
    #debugging stuff
    if(x%%1000==0){
      print(x)
    }
    if(is.null(total)){
      print(paste("Error! ",x," ",otu1," ",otu2))
      return(NA)
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


out<-basic_approach(file2)


