############
# store functions for different methods of calculating co-occurence
############
library(data.table)
library(tidyr)
library(reshape2)

#file = "testdata/OTUs_Table-norm.tab"
#file = "testdata/gut_16s_abundance.txt"
#table2 <- as.data.table(t(as.matrix(read.table(file = file,
#                                               sep = "\t",
#                                               dec = ".",
#                                               header = T,
#                                               row.names = 1))))

# wrapper function to calculate counts values from input files
generate_counts <- function(OTU_table,meta,group_column,cutoff,fc,progress=F){
  meta <- as.data.table(meta)
  
  OTUs <- rownames(OTU_table)
  groups <- unique(meta[[group_column]])
  
  #pick OTU tables for each group
  otus_by_group <- lapply(groups, function(x){
    samples <- meta[meta[[group_column]] == x]$SampleID
    return(OTU_table[,samples])
  })
  if(progress){
    incProgress(1/4)
  }
  
  counts_by_group <- lapply(otus_by_group, function(x){
    return(basic_approach(x,OTUs=OTUs,cutoff=cutoff))
  })
  if(progress){
    incProgress(1/2)
  }
  
  counts <- compare_counts(counts_by_group[[1]],counts_by_group[[2]],fc)
  if(progress){
    incProgress(1/4)
  }
  return(counts)
}

# rows are OTU; column is sample
basic_approach <- function(table,OTUs,cutoff){
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
  
  table[, (samples) := lapply(.SD,function(x){ ifelse(x>cutoff,1,0)}),.SDcols=samples][
    ,(samples) := lapply(.SD,function(x){return(x/1)}),.SDcols=samples]
  
  
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
  return(counts)
}

# 
compare_counts <- function(tab1, tab2, fc){
  #test if both tables are equal in size and order
  if((nrow(tab1) != nrow(tab2)) | sum(tab1$OTU1 != tab2$OTU1) != 0 | sum(tab1$OTU2 != tab2$OTU2) != 0){
    print("This should not happen..")
    return()
  }
  
  #both tables have same amount of rows and same order if calculated with "basic_approach()"
  out_tab <- data.table(OTU1 = tab1$OTU1, OTU2 = tab1$OTU2)
  
  if(fc){
    out_tab$value <- log2(((tab2$value - tab1$value)/tab1$value)+1)
      #log2((tab1$value/tab2$value)+1e-20)#(log2(tab1$value + 1e-20) - log2(tab2$value + 1e-20))
  }else{
    out_tab$value <- (tab1$value - tab2$value)
  }
  
  return(out_tab)
}

# create data table with all unique combinations, where (A,B) == (B,A)
unique_combinations <- function(l1, l2){
  t<-as.data.table(expand.grid(l1,l2))
  t$key<-apply(t, 1, function(x)paste(sort(x), collapse=''))
  t<-subset(t, !duplicated(t$key))
  t$key<-NULL
  names(t)<-c("OTU1","OTU2")
  return(t)
}
