############
# store functions for different methods of calculating co-occurence
############

file = "testdata/OTUs_Table-norm.tab"
file = "testdata/gut_16s_abundance.txt"

library(data.table)
library(tidyr)
library(reshape2)

table2 <- as.data.table(t(as.matrix(read.table(file = file,
                                               sep = "\t",
                                               dec = ".",
                                               header = T,
                                               row.names = 1))))


#wrapper function to calculate counts values from input files
generate_counts<-function(OTU_table, meta, group_column, cutoff, fc){
  meta<-as.data.table(meta)
  
  OTUs<-rownames(OTU_table)
  groups <- unique(meta[[group_column]])
  
  #pick OTU tables for each group
  otus_by_group<-lapply(groups, function(x){
    samples <- meta[meta[[group_column]] == x]$SampleID
    return(OTU_table[,samples])
  })
  
  counts_by_group <- lapply(otus_by_group, function(x){
    return(basic_approach(x,OTUs=OTUs,cutoff=cutoff))
  })
  
  counts <- compare_counts(counts_by_group[[1]],counts_by_group[[2]],fc)
  return(counts)
}

#rows are OTU; column is sample
basic_approach <- function(table, OTUs, cutoff){
  print(Sys.time())

  #save OTU names
  #OTUs are in rows!
  samples <- colnames(table)
  n_otus <- length(OTUs)
  sample_size <- ncol(table)
  
  #binarization & normalization
    
  if(!is.data.table(table)){
    table <- as.data.table(table)
  }
  
  #cutoff: take 0.1 percentile of each sample and then mean over all samples to find cutoff; values below are considered to be 0
  #cutoff <- mean(t(table[, apply(table,1,quantile,probs =.1,na.rm=TRUE)]))
  
  table[, (samples) := lapply(.SD,function(x){ ifelse(x>cutoff,1,0)}),.SDcols =samples][
    ,(samples) := lapply(.SD,function(x){return(x/1)}),.SDcols =samples]
  
  
  #TODO: normalizing after calculation of counts!
  
  mat <- matrix(,nrow=n_otus,ncol=n_otus)
  table<-as.data.frame(table)
  
  start_time <- Sys.time()
  for(i in 1:n_otus){
    colsToCount <- which(table[i,] > 0)
    for(j in 1:n_otus){
      #this skips upper triangle of matrix
      if(i<j){
        next()
      }else{
        mat[i,j]<-sum(table[colsToCount,j])
      }
    }
    #set diagonale to NA
    mat[i,i]<-NA
  }
  end_time <- Sys.time()
  print(paste("Count calculation took",(end_time-start_time),"seconds."))
  
  rownames(mat)<-OTUs
  colnames(mat)<-OTUs
  
  counts<-setNames(melt(mat,na.rm = T),c("OTU1","OTU2","value"))
  
  print(Sys.time())
  return(counts)
}

compare_counts <- function(tab1, tab2, fc){
  
  #test if both tables are equal in size and order
  if((nrow(tab1) != nrow(tab2)) | sum(tab1$OTU1 != tab2$OTU1) != 0 | sum(tab1$OTU2 != tab2$OTU2) != 0){
    print("This should not happen..")
    return()
  }
  
  #both tables have same amount of rows and same order if calculated with "basic_approach()"
  out_tab <- data.table(OTU1 = tab1$OTU1, OTU2 = tab1$OTU2)
  
  if(fc){
    out_tab$value <- (log2(tab1$value + 1e-20) - log2(tab2$value + 1e-20))
  }else{
    out_tab$value <- (tab1$value - tab2$value)
  }
  
  return(out_tab)
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


  



