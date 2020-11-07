############
# store functions for different methods of calculating co-occurence
############
library(data.table)
library(tidyr)
library(reshape2)

# wrapper function to calculate counts values from input files
# parameters:
#@OTU_table: table with OTU abundances, rows are taxa, cols are samples
#@meta: meta data table
#@group_column: group of meta table for which to build network
#@cutoff: abundance values below cutoff will be set to 0
#var1: variable of group
#var2: variable of group
generate_counts <- function(OTU_table,meta,group_column,cutoff,fc,var1,var2,progress=F){
  
  OTUs <- rownames(OTU_table)
  #groups: all variables in the group_column
  groups <- na.exclude(unique(meta[[group_column]]))
  
  #pick OTU tables for each sample-group -> skip entries if group-label is NA there
  otus_by_group <- lapply(groups, function(x){
    samples <- na.exclude(meta[meta[[group_column]] == x,]$SampleID)
    return(OTU_table[,samples])
  })
  names(otus_by_group) <- groups
  if(progress){
    incProgress(1/4)
  }
  
  counts_by_group <- lapply(otus_by_group, function(x){
    return(basic_approach(x,OTUs=OTUs,cutoff=cutoff))
  })
  if(progress){
    incProgress(1/2)
  }
  
  #add option to compare one group variable against all other
  if (var2 == "all"){
    #index of var1 inside the counts_by_group list
    idx1 <- which(names(counts_by_group)==var1)
    #loop over remaining tables in list, do "compare_counts()" and combine into one data frame
    #--> this creates x edges between OTUs, with x== the number of remaining tables
    counts<-lapply(counts_by_group[c(-idx1)], function(x){
      counts_x <- compare_counts(counts_by_group[[var1]], x,fc)
      return(counts_x)
    })
    counts<-rbindlist(counts)
  }else{
    counts <- compare_counts(counts_by_group[[var1]],counts_by_group[[var2]],fc)
  }
  if(progress){
    incProgress(1/4)
  }
  #sort counts (edges) by absolute value, to display most "extreme" edges
  counts<-counts[order(abs(counts$value),decreasing = T),]
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
compare_counts <- function(tab1, tab2, type){
  #both tables have same amount of rows and same order if calculated with basic_approach()
  out_tab <- data.table(OTU1 = tab1$OTU1, OTU2 = tab1$OTU2)
  
  
  if(type==0){
    out_tab$value <- log2((tab2$value+0.001) / (tab1$value+0.001))
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


###############################
#     topological sorting     #
###############################


# runTopologicalSorting<-function(otu,cutoff=0.1){
#   sourceCpp('src/topological_sorting.cpp')
#   otu_names <- rownames(otu)
#   out<-lapply(otu, function(x){
#       calculate_topological_sorting(x,cutoff)
#   })
#   return(out)
# }
# 
# ggplot(data=out_long,aes(x=value,fill=key))+
#   geom_histogram(position = "dodge")


