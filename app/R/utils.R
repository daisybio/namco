calcReadLoss <- function(out, dadaFs, dadaRs, dada_merged, seq_table_nochim, samples, samples_filtered, is_paired){
  out_0_reads <- out[out$reads.out==0,]
  samples_0_reads <- samples[!samples %in% samples_filtered]
  out_more_reads <- out[out$reads.out!=0,]
  getN <- function(x) sum(getUniques(x))
  if (length(samples_filtered) > 1){
    if (is_paired) {
      track <- cbind(samples_filtered, out_more_reads, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(dada_merged, getN), rowSums(seq_table_nochim))
    }else{
      track <- cbind(samples_filtered, out_more_reads, sapply(dadaFs, getN), rowSums(seq_table_nochim))
    }
  }else{
    if (is_paired){
      track <- cbind(samples_filtered, out_more_reads, getN(dadaFs), getN(dadaRs), getN(dada_merged), rowSums(seq_table_nochim))
    } else{
      track <- cbind(samples_filtered, out_more_reads, getN(dadaFs), rowSums(seq_table_nochim))
    }
  }
  rownames(track) <- samples_filtered
  if(is_paired){
    colnames(track) <- c("sample","input_reads", "filtered & trimmed", "denoisedFW", "denoisedRV", "merged", "non_chimera")
  }else{
    colnames(track) <- c("sample","input_reads", "filtered & trimmed", "denoised", "non_chimera")
  }
  
  track2 <- cbind(samples_0_reads, out_0_reads, rep(0, dim(out_0_reads)[1]), rep(0, dim(out_0_reads)[1]), rep(0, dim(out_0_reads)[1]), rep(0, dim(out_0_reads)[1]))
  rownames(track2) <- samples_0_reads
  colnames(track2) <- c("sample","input_reads", "filtered & trimmed", "denoisedFW", "denoisedRV", "merged", "non_chimera")
  
  track <- rbind(track, track2)
  
  return(track)
}

# build phylogenetic tree
# https://f1000research.com/articles/5-1492/v2
buildPhyloTree <- function(seqs, ncores){
  names(seqs) <- seqs
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, processors = 30, verbose = T)
  message(paste0(Sys.time()," - phylogenetic tree: finished alignment step. "))
  
  # write alignment to file
  writeXStringSet(DNAStringSet(alignment), "/tmp/otus.fasta")
  
  # use FastTree to build tree
  command <- paste0("cd /opt && ./FastTreeMP -gtr -nt /tmp/otus.fasta > /tmp/tree.nwk")
  system(command)
  message(paste0(Sys.time()," - phylogenetic tree: finished tree building with FastTree. "))
  
  # load tree back into R and root it
  tree <- ape::read.tree("/tmp/tree.nwk")
  tree <- midpoint(tree)
  
  return(tree)
}


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

#calculate slope values for rarefaction curve
calcSlopes <- function(rarefactionDf,otu){
  slope = vector()
  SampleID = vector()
  otu<-t(otu)
  for (i in seq_along(rarefactionDf)) {
    # If the sequencing depth is greater than 100, the difference between the last and last-100 richness is calculated
    richness <- ifelse(length(rarefactionDf[[i]]) >= 100, rarefactionDf[[i]][length(rarefactionDf[[i]])] - rarefactionDf[[i]][length(rarefactionDf[[i]]) - 100], 1000)
    slope <- c(slope,richness)
    SampleID <- c(SampleID,as.character(rownames(otu)[i]))
  }
  return(cbind(SampleID,slope))
}

# normalize input data (Rhea)
normalizeOTUTable <- function(tab,method=0){
  min_sum <- min(colSums(tab))
  
  if(method==0){
    #no normalization
    norm_tab <- tab
  } else if(method==1){
    # Rarefy the OTU table to an equal sequencing depth
    py.tab <- otu_table(tab,T) #create phyloseq otu_table object
    py.norm_tab <- rarefy_even_depth(py.tab,min_sum,rngseed = 711) #use phyloseq rarefy function
    norm_tab <- as.data.frame(otu_table(py.norm_tab)) 
    
    #norm_tab <- Rarefy(t(tab),depth=min_sum)
    #norm_tab <- t(as.data.frame(norm_tab$otu.tab.rff))
  } else if (method ==2){
    # Divide each value by the sum of the sample and multiply by the minimal sample sum
    norm_tab <- t(min_sum * t(tab) / colSums(tab))
  } else if (method == 3){
    #use centered log-ratio normalization:
    #It is based on dividing each sample by the geometric mean of its values, and taking the logarithm; inlcuding 1 pseudocount to not get negative values
    norm_tab <- log1p(tab/colMeans(tab))
  } else if (method == 4){ 
    # normalize to 10.000 reads per sample
    norm_tab <- t(10000 * t(tab) / colSums(tab))
  }
  
  # Calculate relative abundances for all OTUs over all samples
  # Divide each value by the sum of the sample and multiply by 100
  rel_tab <- relAbundance(tab)
  
  message(paste0(Sys.time()," - Normalized OTU/ASV table with method: ", method))
  return(list(norm_tab=norm_tab,rel_tab=rel_tab))
} 

relAbundance <-function(otu){
  ret <- t(100*t(otu)/colSums(otu))
  ret[is.nan(ret)] <- 0    
  return (ret)
}

relAbundanceTo1 <- function(otu){
  ret <- t(t(otu)/colSums(otu))
  ret[is.nan(ret)] <- 0    
  return (ret)
}

checkTaxonomyColumn <- function(otu){
  
  #stop if no taxonomy column present
  if (! "taxonomy" %in% colnames(otu)){
    return (c(FALSE, noTaxaInOtuError, 0))
  }
  
  taxonomy_col = otu$taxonomy
  col_length = lapply(strsplit(x=as.character(taxonomy_col), ";"), length)
  if (length(unique(col_length)) != 1){
    wrong_rows = paste(unlist(which(col_length != 6)), collapse = ", ")
    return (c(FALSE, wrongTaxaColumnError, wrong_rows))
  }
  return (c(TRUE, 0, 0))
}

# generates a taxonomy table using the taxonomy string of an otu
generateTaxonomyTable <- function(otu){
  taxonomy = otu$taxonomy
  
  splitTax = strsplit(x = as.character(taxonomy),";")
  splitTax=lapply(splitTax,function(x) if(length(x)==6) x=c(x,"") else x=x)
  taxonomy_new = data.frame(matrix(unlist(splitTax),ncol=7,byrow=T))
  colnames(taxonomy_new) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  rownames(taxonomy_new) = rownames(otu)
  
  # Add level information to all taxonomies
  taxonomy_new[,1] <- gsub("^","k__",taxonomy_new[,1]) # For taxonomies related to kingdom level
  taxonomy_new[,2] <- sub("^","p__",taxonomy_new[,2])  # For taxonomies related to phylum level
  taxonomy_new[,3] <- sub("^","c__",taxonomy_new[,3])  # For taxonomies related to class level
  taxonomy_new[,4] <- sub("^","o__",taxonomy_new[,4])  # For taxonomies related to order level
  taxonomy_new[,5] <- sub("^","f__",taxonomy_new[,5])  # For taxonomies related to family level
  taxonomy_new[,6] <- sub("^","g__",taxonomy_new[,6])  # For taxonomies related to genus level
  taxonomy_new[,7] <- sub("^","s__",taxonomy_new[,7])  # For taxonomies related to species level
  
  taxonomy_new <- replaceNATaxonomy(taxonomy_new)
  
  return(taxonomy_new)
}

# replace possible NA entries with *__
replaceNATaxonomy <- function(tax_table){
  
  tax_table["Kingdom"][is.na(tax_table["Kingdom"])]<-"k__"
  tax_table["Phylum"][is.na(tax_table["Phylum"])]<-"p__"
  tax_table["Class"][is.na(tax_table["Class"])]<-"c__"
  tax_table["Order"][is.na(tax_table["Order"])]<-"o__"
  tax_table["Family"][is.na(tax_table["Family"])]<-"f__"
  tax_table["Genus"][is.na(tax_table["Genus"])]<-"g__"
  tax_table["Species"][is.na(tax_table["Species"])]<-"s__"
  
  return(tax_table)
  
}

#run this function to add missing columns to tax-table
addMissingTaxa <- function(taxonomy){
  if(!("Kingdom" %in% colnames(taxonomy))){
    taxonomy$Kingdom <- rep("k__", nrow(taxonomy))
  }
  if(!("Phylum" %in% colnames(taxonomy))){
    taxonomy$Phylum <- rep("p__", nrow(taxonomy))
  }
  if(!("Class" %in% colnames(taxonomy))){
    taxonomy$Class <- rep("c__", nrow(taxonomy))
  }
  if(!("Order" %in% colnames(taxonomy))){
    taxonomy$Order <- rep("o__", nrow(taxonomy))
  }
  if(!("Family" %in% colnames(taxonomy))){
    taxonomy$Family <- rep("f__", nrow(taxonomy))
  }
  if(!("Genus" %in% colnames(taxonomy))){
    taxonomy$Genus <- rep("g__", nrow(taxonomy))
  }
  if(!("Species" %in% colnames(taxonomy))){
    taxonomy$Species <- rep("s__", nrow(taxonomy))
  }
  return(taxonomy)
}

taxBinningNew <- function(phylo, is_fastq){
  mdf <- as.data.table(psmelt(phylo))
  if (is_fastq){
    taxas <- c("Kingdom","Phylum","Class","Order","Family","Genus")
  }else{
    taxas <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  }
  mdf[Kingdom == "k__", Kingdom:="unknown"]
  mdf[Phylum == "p__", Phylum:="unknown"]
  mdf[Class == "c__", Class:="unknown"]
  mdf[Order == "o__", Order:="unknown"]
  mdf[Family == "f__", Family:="unknown"]
  mdf[Genus == "g__", Genus:="unknown"]
  if(!is_fastq){mdf[Species == "s__", Species:="unknown"]}
  
  out_l<-lapply(taxas, function(x){
    return(acast(mdf[,sum(Abundance),by=c(x,"Sample")],as.formula(paste(x,"~ Sample")),value.var = "V1"))
  })
  return(out_l)
}

# pool taxa and fill empty labels with higher-rank labels
# top_k: output only the top k taxa in the chosen level (by total abundance)
glom_taxa_custom <- function(phylo, rank, top_k = NULL){
  phylo_rank <- phyloseq::tax_glom(phylo, taxrank = rank)
  
  if(!is.null(top_k)){
    top_taxa <- names(sort(rowSums(phylo_rank@otu_table@.Data), decreasing = T)[1:top_k])
    phylo_rank <- prune_taxa(top_taxa, phylo_rank)
  }
  
  taxtab <- phylo_rank@tax_table@.Data
  
  miss_k <- which(taxtab[, "Kingdom"] == "k__")
  miss_p <- which(taxtab[, "Phylum"] == "p__")
  miss_c <- which(taxtab[, "Class"] == "c__")
  miss_o <- which(taxtab[, "Order"] == "o__")
  miss_f <- which(taxtab[, "Family"] == "f__")
  miss_g <- which(taxtab[, "Genus"] == "g__")
  
  taxtab[miss_k, "Kingdom"] <- paste0("k__", 1:length(miss_k))
  taxtab[miss_p, "Phylum"] <- paste0("p__", 1:length(miss_p))
  taxtab[miss_c, "Class"] <- paste0("c__", 1:length(miss_c))
  taxtab[miss_o, "Order"] <- paste0("o__", 1:length(miss_o))
  taxtab[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
  taxtab[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))
  
  
  for(i in seq(taxtab)){
   # The next higher non-missing rank is assigned to unspecified genera
   if(i %in% miss_f && i %in% miss_g && i %in% miss_o && i %in% miss_c && i %in% miss_p && i %in% miss_k){
     taxtab[i, rank] <- paste0(taxtab[i, rank], "(no_rank)")
   } else if(i %in% miss_p){
     taxtab[i, rank] <- paste0(taxtab[i, rank], "(", taxtab[i, "Kingdom"], ")")
   } else if(i %in% miss_c){
     taxtab[i, rank] <- paste0(taxtab[i, rank], "(", taxtab[i, "Phylum"], ")")
   } else if( i %in% miss_o){
     taxtab[i, rank] <- paste0(taxtab[i, rank], "(", taxtab[i, "Class"], ")")
   } else if( i %in% miss_f){
     taxtab[i, rank] <- paste0(taxtab[i, rank], "(", taxtab[i, "Order"], ")")
   } else if(i %in% miss_g ){
     taxtab[i, rank] <- paste0(taxtab[i, rank], "(", taxtab[i, "Family"], ")")
   }
  }
  
  # it can happen that the same taxonomic level appears more than once in this phylo_rank table
  # https://github.com/joey711/phyloseq/issues/927
  # i will add a running index to those entries so that are unique for later on
  #taxtab_tmp <- data.frame(phylo_rank@tax_table)
  #colnames(taxtab) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  taxtab <- as.data.frame(taxtab)
  non_unique <- c(names(which(table(taxtab[[rank]]) > 1)))
  if(length(non_unique)==1){
    ids <- which(taxtab[rank] == non_unique)
    taxtab[rank][ids,] <- paste(non_unique,rep(1:length(ids)))

    taxtab_tmp <- tax_table(taxtab)
    taxa_names(taxtab_tmp) <- taxa_names(phylo_rank)
    tax_table(phylo_rank) <- taxtab_tmp
  }
  
  taxa_names(phylo_rank) <- taxtab[[rank]]
  rownames(taxtab) <- taxtab[[rank]]
  
  return(list(taxtab=data.frame(taxtab), phylo_rank=phylo_rank))
}


# calculate alpha diversity scores and merge with meta table if there is one
createAlphaTab <- function(otu, meta=NULL){
  
  alphaTab = data.frame(colnames(otu))
  for(i in c("Shannon_Entropy","effective_Shannon_Entropy","Simpson_Index","effective_Simpson_Index","Richness")){
    alphaTab = cbind(alphaTab,round(alphaDiv(otu,i),2))
  }
  colnames(alphaTab) = c("SampleID","Shannon_Entropy","effective_Shannon_Entropy","Simpson_Index","effective_Simpson_Index","Richness")
  
  if(!is.null(meta)){
    # need to check if there are already columns with these scores -> rename them
    columns <- which(colnames(meta) %in% c("Shannon_Entropy","effective_Shannon_Entropy","Simpson_Index","effective_Simpson_Index","Richness"))
    if(!is.null(columns)){
      new_names <- paste0("provided.", colnames(meta)[columns])
      colnames(meta)[columns] <- new_names
    }
    
    alphaTab <- merge(alphaTab, meta, by.x="SampleID", by.y="SampleID")
  }
  return(alphaTab)

}

# calculate various measures of beta diversity
betaDiversity <- function(phylo,method){
  
  otu <- as.data.frame(otu_table(phylo))
  otu <- t(otu[,order(names(otu))])
  
  meta <- as.data.frame(sample_data(phylo))
  meta <- data.frame(meta[order(rownames(meta)),])
  
  if(!is.null(access(phylo,"phy_tree"))) tree <- phy_tree(phylo) else tree <- NULL
  
  if(method == 1){ #Bray-Curtis Dissimilarity
    return(vegdist(otu, method="bray"))
  }else{
    rooted_tree = midpoint(tree)
    unifracs = GUniFrac(otu,rooted_tree,alpha=c(0.0,0.5,1.0))$unifracs
    if(method == 2){ #Generalized UniFrac Distance
      unifract_dist = unifracs[,,"d_0.5"]
      return(as.dist(unifract_dist))
    }
    else if(method == 3){ #Unweighted UniFrac Distance
      unifract_dist = unifracs[,,"d_UW"]
      return(as.dist(unifract_dist))
    }
    else if(method == 4){ #Weighted UniFrac Distance
      unifract_dist = unifracs[,,"d_1"]
      return(as.dist(unifract_dist))
    }
    else if(method == 5){ #Variance adjusted weighted UniFrac Distance
      unifract_dist = unifracs[,,"d_VAW"]
      return(as.dist(unifract_dist))
    }
  }
}

# remove rows with missing information from a subset of selected columns
completeFun <- function(data, desiredCols) {
  names<-rownames(data)
  completeVec <- complete.cases(data[, desiredCols])
  return(data.frame(data[completeVec, ],row.names = names[completeVec]))
}

#build unifrac distance matrix
buildGUniFracMatrix <- function(otu,tree){
  
  message("Calculating unifrac distance matrix...")
  # Order the OTU-table by sample names (ascending)
  otu <- otu[,order(colnames(otu))]
  # Transpose OTU-table and convert format to a data frame
  otu <- data.frame(t(otu), check.names = F)
  # Root the OTU tree at midpoint 
  if(!is.rooted(tree)){
    tree <- midpoint(tree)
  }
  # Calculate the UniFrac distance matrix for comparing microbial communities
  unifracs <- GUniFrac(otu, tree, alpha = c(0.0,0.5,1.0))$unifracs
  # Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
  unifract_dist <- unifracs[, , "d_0.5"]
  
  return (unifract_dist)
}

#calculate confounding factors given a single variable to test
calculateConfounderTable <- function(var_to_test,variables,distance,seed,progress=T){
  
  set.seed(seed)
  
  namelist <- vector()
  confounderlist <-vector()
  directionList <- vector()
  pvalList <- vector()
  #variables[is.na(variables)]<-"none"
  loops <- dim(variables)[2]
  for (i in 1:loops) {
    if (dim(unique(variables[, i]))[1] > 1) {
      variables_nc <- completeFun(variables, i)
      position <- which(row.names(distance) %in% row.names(variables_nc))
      dist <- distance[position, position]
      
      # Test outcome without variables
      without <- adonis2(as.formula(paste0("dist ~ ", var_to_test)), data = variables_nc)
      #Test outcome with variable
      with <- adonis2(as.formula(paste0("dist ~ ",var_to_test," + ", colnames(variables_nc)[i])), data = variables_nc)

      names <- names(variables_nc)[i]
      namelist <- append(namelist, names)
      pval_without <- without[["Pr(>F)"]][1]
      pval_with <- with[["Pr(>F)"]][1]
      if (pval_without <= 0.05) {
        # variable is significant
        if (pval_with <= 0.05) {
          # no confounder
          confounder <- "NO"
          direction <- "signficant"
        }
        else {
          # confounder
          confounder <- "YES"
          direction <- "not_signficant"
        }
      } else {
        # variable is not signficant
        if (pval_with <= 0.05) {
          #  confounder
          confounder <- "YES"
          direction <- "signficant"
        }
        else {
          # no confounder
          confounder <- "NO"
          direction <- "not_signficant"
        }
      }
      
      confounderlist <- append(confounderlist, confounder)
      directionList <- append(directionList, direction)
      pvalList <- append(pvalList, pval_without)
      incProgress(amount=1/loops)
    }
  }
  
  df <- data.frame(name = namelist, confounder = confounderlist, value = directionList)
  #remove var-to-test from output dataframe; makes no sense that variable is confounding factor for itself
  df <- df [!(df$name == var_to_test),]
  #writeOutputHeader(df,paste(directory,"confounder.tab",sep=""))
  return(list(table=df,var=var_to_test))
}


#for all possible variables to test (all columns in meta-data) calculate the confounding factors for each and build matrix 
#matrix-rows: variable which is compared to the other variables(columns) -> diagonal should always be 0/NO
calculateConfounderMatrix <- function(meta, distance_matrix,progress=T){
  #look at all vairables except first (SampleID)
  vars_to_test <- unique(colnames(meta))[-1]
  meta[is.na(meta)]<-0
  meta[,1]<-NULL
  
  #create matrix; row order corresponds to vars_to_test order
  #Var1 = variable that was tested; Var2 the other variables, to test against
  m<-do.call(rbind, lapply(vars_to_test,function(x){
    incProgress(1/length(vars_to_test))
    return(calculateConfounderTable(x,meta,distance_matrix)$confounder)
  }))
  rownames(m)<-colnames(m)<-vars_to_test
  return(m)
}


# calculate various measures of alpha diversity
alphaDiv <- function(tab,method){
  alpha = apply(tab,2,function(x) switch(method,
                                         "Shannon_Entropy" = {-sum(x[x>0]/sum(x)*log(x[x>0]/sum(x)))},
                                         "effective_Shannon_Entropy" = {round(exp(-sum(x[x>0]/sum(x)*log(x[x>0]/sum(x)))),digits=2)},
                                         "Simpson_Index" = {sum((x[x>0]/sum(x))^2)},
                                         "effective_Simpson_Index" = {round(1/sum((x[x>0]/sum(x))^2),digits =2)},
                                         "Richness" = {sum(x>0)}
  ))
  
  return(alpha)
}

## Plotting of confusion matrix
# Code found on StackOverFlow (https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-matrix-using-the-caret-package)
draw_confusion_matrix <- function(cmtrx) {
  
  total <- sum(cmtrx$table)
  res <- as.numeric(cmtrx$table)
  
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
    
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix
  classes = colnames(cmtrx$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cmtrx results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cmtrx$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cmtrx$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cmtrx$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cmtrx$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cmtrx$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cmtrx$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cmtrx$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cmtrx$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cmtrx$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cmtrx$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information
  text(30, 35, names(cmtrx$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cmtrx$overall[1]), 3), cex=1.4)
  text(70, 35, names(cmtrx$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cmtrx$overall[2]), 3), cex=1.4)
  
}


#extract all numbers as list from comma seperated text
extract <- function(text){
  text <- gsub(" ", "", text)
  split <- strsplit(text, ",", fixed = FALSE)[[1]]
  as.numeric(split)
}

#build dataframe for ML model
#parameters: meta table, otu table (rows=samples,col=otu),input of shiny ui
buildForestDataset <- function(meta, otu, input){
  
  #use cutoff slider for numeric/continuous variables to transform it into 2 categories
  if(is.numeric(meta[[input$forest_variable]])){
    variable = cut(meta[[input$forest_variable]],breaks = c(-Inf,input$forest_continuous_slider,Inf),labels = c("low","high"))
    #add variable vector to OTU dataframe
    tmp <- data.frame(variable=as.factor(variable))
  }else{
    #use covariable if more than 2 groups
    if(length(unique(meta[[input$forest_variable]])) > 2 && !is.null(input$forest_covariable)){
      c <- replace(meta[[input$forest_variable]], which(meta[[input$forest_variable]] != input$forest_covariable), "rest")
      tmp <- data.frame(variable = as.factor(c))
    }else{
      #add variable of interest to otu dataframe 
      tmp <- data.frame(variable=as.factor(meta[[input$forest_variable]]))
    }
  }
  
  features <- input$forest_features
  #add meta features to model
  if(!is.null(features)){
    tmp <- cbind.data.frame(tmp, meta[,features],stringsAsFactors=T)
  }
  #add otu abundances to model
  if(input$forest_otu){
    tmp <- cbind.data.frame(tmp,otu,stringsAsFactors=T)
  }
  
  #remove rows, where variable is NA
  combined_data <- tmp[complete.cases(tmp),]
  
  #remove OTUs which the user wants to exclude from model
  if(!is.null(input$forest_exclude)){
    combined_data <- combined_data[, -which(names(combined_data) %in% input$forest_exclude)] 
  }
  
  return(combined_data)
  
}


#build dataframe to plot roc for every CV fold
# from here: https://stackoverflow.com/questions/51904700/how-to-plot-roc-curves-for-every-cross-validations-using-caret
lift_df <- function(model,parameter){
  for_lift <- data.frame(Class = model$pred$obs, rf = model$pred[[eval(parameter)]], resample = model$pred$Resample)
  lift_df <-  data.frame()
  for (fold in unique(for_lift$resample)) {
    fold_df <- dplyr::filter(for_lift, resample == fold)
    lift_obj_data <- lift(Class ~ rf, data = fold_df, class = eval(parameter))$data
    lift_obj_data$fold = fold
    lift_df = rbind(lift_df, lift_obj_data)
  }
  #lift_obj <- lift(Class ~ rf, data = for_lift, class = parameter)
  return(lift_df)
}

errorModal <- function(error_message=NULL){
  modalDialog(
    p(error_message,style="color:red;"),
    easyClose = T,
    modalButton("Cancel")
  )
}

# filtering function
# message object in f_list tells you if filtering was effective; if it is NULL, everything went OK
applyFilterFunc <- function(phylo, keep_taxa, f_list_old){
  message <- NULL
  if(is.null(keep_taxa)){message <- noTaxaRemainingAfterFilterError}
  if(length(keep_taxa) == 0){message <- noTaxaRemainingAfterFilterError}
  if(length(keep_taxa)==ntaxa(phylo)){message <- filteringHadNoEffectError}
  if(!is.null(message)){
    f_list_old$message <- message
    return(f_list_old)
  }else{
    tryCatch({
      #filter is applied here 
      phylo_pruned <- prune_taxa(keep_taxa, phylo)
      x <- taxa_sums(phylo_pruned)
      return(list(phylo=phylo_pruned, 
                  x=x, 
                  otu=as.data.frame(otu_table(phylo_pruned)), 
                  rel_otu= relAbundance(as.data.frame(otu_table(phylo_pruned))),
                  message=message))
    },error=function(e){
      f_list_old$message <- e$message
      return(f_list_old)
    })
  }
}

# add vertical line to plotly
# https://stackoverflow.com/questions/34093169/horizontal-vertical-line-in-plotly/34097929#34097929
vline <- function(x = 0, color = "red") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color)
  )
}



#### functions from Rhea ####
fill_NA_INF.mean <- function(vec)
{
  # Calculate mean value of each column (excluding missing values & Inf values)
  m <- mean(vec[!is.infinite(vec)], na.rm = TRUE)
  # Replace missing values with mean
  vec[is.na(vec)] <- m
  # Replace Inf values with mean
  vec[is.infinite(vec)] <- m
  # Return the new input data frame
  return(vec)
}

# Depending on which parameters were set at the beginning, one query type is selected
# In each query type, three matrices are generated: p-value matrix, correlation matrix, support matrix
# All possibles pairs are saved in a vector (pairs)
subsetCorrelation <- function(includeTax, includeMeta, var_names, otu_names, meta_names, my_rcorr, signif_cutoff){
  if(includeTax & !includeMeta){
    
    # Correlation among OTUs and NO correlation among meta-variables
    row_names <- otu_names
    col_names <- var_names
    pairs <-expand.grid(row_names, col_names)
    my_cor_matrix <- my_rcorr$r[otu_names,]
    my_pvl_matrix <-my_rcorr$P[otu_names,]
    my_num_matrix <- my_rcorr$n[otu_names,]
    
    # Set variable for plotting
    diagonale=0
    
  } else if(includeTax & includeMeta){
    
    # Correlation among OTUs and correlation among meta-variables
    row_names <-var_names
    col_names <- var_names
    pairs <-expand.grid(row_names, col_names)
    my_cor_matrix <- my_rcorr$r
    my_pvl_matrix <-my_rcorr$P
    my_num_matrix <- my_rcorr$n
    # Set variable for plotting
    diagonale=0
    
  } else if (!includeTax & includeMeta) {
    # NO correlation among OTUs and correlation among meta-variables
    
    row_names <- meta_names
    col_names <- var_names
    pairs <-expand.grid(row_names, col_names)
    my_cor_matrix <- my_rcorr$r[meta_names,]
    my_pvl_matrix <-my_rcorr$P[meta_names,]
    my_num_matrix <- my_rcorr$n[meta_names,]
    # Set variable for plotting
    diagonale=1
    
  } else {
    # NO correlation among OTUs and NO correlation among meta-variables
    
    row_names <- meta_names
    col_names <- otu_names
    pairs <-expand.grid(row_names, col_names)
    my_cor_matrix <- my_rcorr$r[meta_names,otu_names]
    my_pvl_matrix <-my_rcorr$P[meta_names,otu_names]
    my_num_matrix <- my_rcorr$n[meta_names,otu_names]
    # Set variable for plotting
    diagonale=1
    
  }
  
  # Select the corresponding p-value for each pair 
  p_vector <- as.vector(my_pvl_matrix)
  
  # Select the corresponding correlation coefficient for each pair 
  c_vector <- as.vector(my_cor_matrix)
  
  # Select the corresponding number of observations for each pair 
  n_vector <- as.vector(my_num_matrix)
  
  # Generate matrix with the pairwise comparisons
  my_pairs <-
    matrix(ncol = 5,
           c(
             as.character(pairs[, 2]),
             as.character(pairs[, 1]),
             c_vector,
             p_vector,
             n_vector
           ))
  
  # Delete all pairs with insufficient number of pairs 
  #TODO: default is 4; can be added as parameter later..
  my_pairs <- subset(my_pairs, as.numeric(my_pairs[,5]) > 4)
  
  # Adjust p-value for multiple testing using the Benjamin-Hochberg method
  pVal_BH <- round(p.adjust(my_pairs[,4], method = "BH"), 4)
  
  # Add the corrected p-value in the table 
  my_pairs <- cbind(my_pairs,as.numeric(pVal_BH))
  
  # Remove similar pairs (values along the diagonal)
  my_pairs <- my_pairs[!as.character(my_pairs[, 1]) == as.character(my_pairs[, 2]),]
  
  # Remove duplicate pairs
  my_pairs <- my_pairs[!duplicated(my_pairs[, 3]), ]
  
  colnames(my_pairs) <- c(
    "variable1",
    "variable2",
    "correlation",
    "pValue",
    "support",
    "Corrected"
  )
  
  ## remove non-significant pairs from correlation & pvalue matrix
  # Create subset of pairs with significant p-values
  my_pairs_cutoff <- my_pairs[as.numeric(my_pairs[, "pValue"]) <= signif_cutoff, ]
  
  # Missing values in the correlation matrix are set to zero
  my_cor_matrix[is.na(my_cor_matrix)] <- 0
  
  # only consider those pairs for plotting, which are significant
  if(length(meta_names)==1 && !includeTax){
    my_cor_matrix <- as.matrix(my_cor_matrix[which(names(my_cor_matrix) %in% unique(my_pairs_cutoff[,1]))])
    colnames(my_cor_matrix)<-meta_names
    my_cor_matrix <- t(my_cor_matrix)
    my_pvl_matrix <- as.matrix(my_pvl_matrix[which(names(my_pvl_matrix) %in% unique(my_pairs_cutoff[,1]))])
    colnames(my_pvl_matrix)<-meta_names
    my_pvl_matrix<-t(my_pvl_matrix)
  }else{
    my_cor_matrix <- my_cor_matrix[,colnames(my_cor_matrix) %in% unique(my_pairs_cutoff[,1])]
    my_pvl_matrix <- my_pvl_matrix[,colnames(my_pvl_matrix) %in% unique(my_pairs_cutoff[,1])]
  }
  
  return(list(diagonale=diagonale, 
              my_pairs=my_pairs,
              my_cor_matrix=my_cor_matrix,
              my_pvl_matrix=my_pvl_matrix))
}

plot_correlation_custom <- function(my_cor_matrix, my_pvl_matrix, input){
  if(input$corrPval=="highlight"){
    corrplot(my_cor_matrix, tl.col="black", tl.srt = 65, tl.cex = 0.4, cl.cex = 0.4, p.mat=my_pvl_matrix,
             sig.level = input$corrSignifCutoff)
  }else if(input$corrPval=="blank"){
    corrplot(my_cor_matrix, tl.col="black", tl.srt = 65, tl.cex = 0.4, cl.cex = 0.4, p.mat=my_pvl_matrix,
             sig.level = input$corrSignifCutoff, insig = "blank")
  }else{
    corrplot(my_cor_matrix, tl.col="black", tl.srt = 65, tl.cex = 0.4, cl.cex = 0.4)
  }
}
