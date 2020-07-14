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
      norm_tab <- t(min_sum * t(tab) / colSums(tab))
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
  
  # Generate tables for each taxonomic class (relative values)
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

# remove rows with missing information from a subset of selected columns
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data.frame(data[completeVec, ]))
}

#build unifrac distance matrix
buildDistanceMatrix <- function(otu,meta,tree){
  
  print("Calculating unifrac distance matrix...")
  # Order the OTU-table by sample names (ascending)
  otu <- otu[,order(colnames(otu))]
  # Transpose OTU-table and convert format to a data frame
  otu <- data.frame(t(otu))
  # Root the OTU tree at midpoint 
  rooted_tree <- midpoint(tree)
  # Order the mapping file by sample names (ascending)
  meta <- data.frame(meta[order(row.names(meta)),])
  # Calculate the UniFrac distance matrix for comparing microbial communities
  unifracs <- GUniFrac(otu, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
  # Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
  unifract_dist <- unifracs[, , "d_0.5"]
  
  return (unifract_dist)
}

#calculate confounding factors given a single variable to test
calculateConfounderTable <- function(var_to_test,variables,distance,useSeed,progress=T){
  
  if(useSeed) set.seed(123)

  
  
  namelist <- vector()
  confounderlist <-vector()
  directionList <- vector()
  variables[is.na(variables)]<-"none"
  loops <- dim(variables)[2]
  for (i in 1:loops) {
    if (length(unique(variables[, i])) > 1) {
      variables_nc <- completeFun(variables, i)
      position <- which(row.names(distance) %in% row.names(variables_nc))
      
      # Test outcome with variables
      without <-
        adonis(distance[position, position] ~ variables_nc[, var_to_test])
      #Test outcome without variable
      with <-
        adonis(distance[position, position] ~ variables_nc[, var_to_test] + variables_nc[, i])
      
      names <- names(variables_nc)[i]
      namelist <- append(namelist, names)
      if (without$aov.tab[1, "Pr(>F)"] <= 0.05) {
        # variable is significant
        if (with$aov.tab[1, "Pr(>F)"] <= 0.05) {
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
        if (with$aov.tab[1, "Pr(>F)"] <= 0.05) {
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
      incProgress(amount=1/loops)
    }
  }
  
  df <- data.frame(name = namelist, confounder = confounderlist, value = directionList)
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
                                         "Shannon Entropy" = {-sum(x[x>0]/sum(x)*log(x[x>0]/sum(x)))},
                                         "effective Shannon Entropy" = {round(exp(-sum(x[x>0]/sum(x)*log(x[x>0]/sum(x)))),digits=2)},
                                         "Simpson Index" = {sum((x[x>0]/sum(x))^2)},
                                         "effective Simpson Index" = {round(1/sum((x[x>0]/sum(x))^2),digits =2)},
                                         "Richness" = {sum(x>0)}
  ))
  
  return(alpha)
}

## Plotting of confusion matrix
#Code found on StackOverFlow (https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-matrix-using-the-caret-package)
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