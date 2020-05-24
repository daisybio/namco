# Jensen Shannon distance
jsd <- function(p,q) {
  
  m <- .5*(p + q)
  div <- .5*sum(p*log(p/m)) + .5*sum(q*log(q/m))
  dis <- sqrt(div)
  
  return(dis)
}

# Normalize to 0-1 range
norm10 <- function(x){
  return((x-min(x))/(max(x)-min(x)))
} 


# normalize input data (Rhea)
normalizeOTUTable <- function(tab,method=0,normalized=F){
  min_sum <- min(colSums(tab))
  
  if(!normalized){
    if(method==0){
      # Divide each value by the sum of the sample and multiply by the minimal sample sum
      norm_tab <- t(min_sum*t(tab)/colSums(tab))
    } else{
      # Rarefy the OTU table to an equal sequencing depth
      norm_tab <- Rarefy(t(tab),depth=min_sum)
      norm_tab <- t(as.data.frame(norm_tab$otu.tab.rff))
    }
  } else norm_tab = tab
  
  # Calculate relative abundances for all OTUs over all samples
  # Divide each value by the sum of the sample and multiply by 100
  rel_tab <- t(100*t(tab)/colSums(tab))
  
  return(list(norm_tab=norm_tab,rel_tab=rel_tab))
} 

# generates a taxonomy table using the taxonomy string of an otu
generateTaxonomyTable <- function(otu){
  taxonomy = otu[,ncol(otu)]
  
  splitTax = strsplit(x = as.character(taxonomy),";")
  splitTax=lapply(splitTax,function(x) if(length(x)==6) x=c(x,"") else x=x)
  taxonomy_new = data.frame(matrix(unlist(splitTax),ncol=7,byrow=T))
  colnames(taxonomy_new) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  rownames(taxonomy_new) = rownames(otu)
  
  # Add level information to all taxonomies
  taxonomy_new[,1] <- gsub("^","k__",taxonomy_new[,1]) # For taxonomies related to kingdom level
  taxonomy_new[,2] <- sub("^","p__",taxonomy_new[,2]) # For taxonomies related to phylum level
  taxonomy_new[,3] <- sub("^","c__",taxonomy_new[,3]) # For taxonomies related to class level
  taxonomy_new[,4] <- sub("^","o__",taxonomy_new[,4]) # For taxonomies related to order level
  taxonomy_new[,5] <- sub("^","f__",taxonomy_new[,5]) # For taxonomies related to family level
  taxonomy_new[,6] <- sub("^","g__",taxonomy_new[,6]) # For taxonomies related to genus level
  taxonomy_new[,7] <- sub("^","s__",taxonomy_new[,7]) # For taxonomies related to species level
  
  return(taxonomy_new)
}

# return a list of tables with taxonomically binned samples
taxBinning <- function(otuFile,taxonomy){
  # Create empty list for further processing
  sample_list <- vector(mode="list",length=7)
  list_length <- NULL
  
  # Preallocate list with right dimensions
  for(i in 1:7){
    list_length[i] <- num_taxa <- length(unique(taxonomy[,i]))
    
    for(j in 1:num_taxa){
      # Initialize list with the value zero for all taxonomies
      sample_list[[i]][[j]] <- list(rep.int(0,ncol(otuFile)))
    }
  }
  
  # Save relative abundances of all samples for each taxonomy
  for(i in 1:nrow(otuFile)){
    for(j in 1:7){
      taxa_in_list <- taxonomy[i,j]
      position <- which(unique(taxonomy[,j]) == taxa_in_list) # Record the current position
      
      sub_sample_tax <- otuFile[taxonomy[,j] == taxa_in_list,] # Filter rows with given taxonomic identifier
      if(length(dim(sub_sample_tax))>1) temp = colSums(sub_sample_tax) else temp = sub_sample_tax # Calculate the summed up relative abundances for the particular taxonomic class for n-th sample
      
      # Replace values by new summed values
      sample_list[[j]][[position]] <- list(temp)
    }
  }
  
  # Generate tables for each taxonomic class
  out_list <- vector(mode="list",length=7)
  for(i in 1:7){
    mat = matrix(unlist(sample_list[[i]]),nrow=list_length[i],ncol=ncol(otuFile),byrow=T,dimnames=list(unique(taxonomy[,i]),colnames(otuFile)))
    rownames(mat) = sapply(rownames(mat),substring,4)
    rownames(mat)[rownames(mat)==""] = "unknown"
    if(nrow(mat)>1) mat = mat[order(rownames(mat)),]
    out_list[[i]] = mat
  }
  
  return(out_list)
}

# calculate various measures of beta diversity
betaDiversity <- function(otu,meta,tree,group,method){
  otu = t(otu[,order(colnames(otu))])
  meta = meta[order(rownames(meta)),]

  all_groups = as.factor(meta[[group]])
  
  switch(method,
         uniFrac = {
           # rooting the tree
           rooted_tree = midpoint(tree)
           unifracs = GUniFrac(otu,rooted_tree,alpha=c(0.0,0.5,1.0))$unifracs
           unifract_dist = unifracs[,,"d_0.5"]
           return(as.dist(unifract_dist))
         },
         brayCurtis = {
           return(vegdist(otu))
         }
  )
}

# calculate various measures of alpha diversity
alphaDiv <- function(tab,method){
  alpha = apply(tab,2,function(x) switch(method,
                                         "Shannon Entropy" = {-sum(x[x>0]/sum(x)*log(x[x>0]/sum(x)))},
                                         "effective Shannon Entropy" = {round(exp(-sum(x[x>0]/sum(x)*log(x[x>0]/sum(x)))),digits=2)},
                                         "Simpson Index" = {sum((x[x>0]/sum(x))^2)},
                                         "effective Simpson Index" = {round(1/sum((x[x>0]/sum(x))^2),digits =2)},
                                         "Richness" = {sum(x>0)}
  ))
  
  return(alpha)
}