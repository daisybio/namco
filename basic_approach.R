basic_approach <- function(tab){
  #save OTU names
  OTUs = colnames(tab)
  n_otus = length(OTUs)
  sample_size = nrow(tab)
  
  #cutoff: take 0.1 percentile of each sample and then mean over all samples to find cutoff; values below are considered to be 0
  cutoff = mean(apply(tab,1,quantile,probs =.1,na.rm=TRUE))
  
  #binarization & normalization
  #cutoff nicht bei 1, sondern percentile wise...
  tab = ifelse(tab>cutoff,1,0)/sample_size
  
  #TODO: normalizing after calculation of counts!
  
  mat <- matrix(NA,nrow=n_otus,ncol=n_otus)
  for(i in 1:(n_otus-1)){
    for(j in (i+1):n_otus){
      mat[i,j] <- sum(tab[which(tab[,i]>0),j])
      #TODO: test this!
      #mat[,j]<-colSums(tab[rowsToCount,])
    }
  }
  rownames(mat) = OTUs
  colnames(mat) = OTUs
  
  counts <- setNames(reshape2::melt(mat,na.rm=T),c("OTU1","OTU2","value"))
  return(counts)
}
