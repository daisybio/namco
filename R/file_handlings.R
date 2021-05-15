##### dada2 - related functions #####

combineAndNormalize <- function(seq_table, taxonomy, has_meta, meta, tree, sn, abundance_cutoff){
  # generate unnormalized object to get unified ASV names
  if(has_meta){
    phylo_unnormalized <- phyloseq(otu_table(seq_table, taxa_are_rows = F),
                                   sample_data(meta),
                                   tax_table(as.matrix(taxonomy)))
  }else{
    phylo_unnormalized <- phyloseq(otu_table(seq_table, taxa_are_rows = F),
                                   tax_table(as.matrix(taxonomy)))
    sample_names(phylo_unnormalized) <- sn
  } 
  
  # add tree to phyloseq object
  if(!is.null(tree)){
    phylo_unnormalized <- merge_phyloseq(phylo_unnormalized, tree)
  }
  
  # store sequences of ASV in object
  dna <- Biostrings::DNAStringSet(taxa_names(phylo_unnormalized))
  names(dna) <- taxa_names(phylo_unnormalized)
  phylo_unnormalized <- merge_phyloseq(phylo_unnormalized, dna)
  
  # give "normal" names to ASVs
  taxa_names(phylo_unnormalized) <- paste0("ASV", seq(ntaxa(phylo_unnormalized)))
  message(paste0("First phyloseq-object: ", ntaxa(phylo_unnormalized)))
  
  # abundance filtering
  message(paste0("Filtering out ASVs with total abundance below ", abundance_cutoff, "% abundance"))
  abundance_cutoff <- abundance_cutoff/100
  phylo_unnormalized <- removeLowAbundantOTUs(phylo_unnormalized, abundance_cutoff, "fastq")
  
  # create final objects with "real" ASV names
  asv_table_final <- t(as.data.frame(otu_table(phylo_unnormalized, F)))
  taxonomy_final <- as.data.frame(tax_table(phylo_unnormalized))
  if(has_meta){meta_final <- as.data.frame(sample_data(phylo_unnormalized))}else{meta_final <- NULL}
  refseq_final <- refseq(phylo_unnormalized)
  
  # normalization (none applied during first data upload)
  normalized_asv = list(norm_tab=asv_table_final, rel_tab = relAbundance(asv_table_final))
  
  # final object
  if(has_meta){
    phylo <- phyloseq(otu_table(normalized_asv$norm_tab, taxa_are_rows = T),
                      sample_data(meta_final),
                      tax_table(as.matrix(taxonomy_final)),
                      refseq_final)
  }else{
    phylo <- phyloseq(otu_table(normalized_asv$norm_tab, taxa_are_rows = T),
                      tax_table(as.matrix(taxonomy_final)),
                      refseq_final)
  }
  if(!is.null(tree)){phylo <- merge_phyloseq(phylo, phy_tree(phylo_unnormalized))}
  
  message(paste0(Sys.time()," - final phyloseq-object: ", ntaxa(phylo)))
  print(phylo)
  
  return(list(phylo=phylo,
              normalized_asv=normalized_asv,
              raw_asv=asv_table_final,
              meta=meta_final,
              taxonomy=taxonomy_final))
}


removeLowAbundantOTUs <- function(phy, cutoff, mode){
  if(mode=="fastq"){otu_tab <- t(as.data.frame(otu_table(phy)))}
  if(mode=="otu"){otu_tab <- as.data.frame(otu_table(phy))}
  
  keep_otus<-do.call(rbind, lapply((1:nrow(otu_tab)), function(x){
    row <- otu_tab[x,]
    asv <- rownames(otu_tab)[x]
    keep = F
    for(i in (1:ncol(otu_tab))){
      perc <- (row[i])/(colSums(otu_tab)[i])
      if (perc > 0.0025){
        keep = T
        break
      }
    }
    if (keep){
      return(asv)
    }
  }))
  
  filtered_phylo <- prune_taxa(keep_otus[,1], phy)
  # differently detailed outputs depending on mode
  if(mode=="fastq"){
    return(filtered_phylo)
  }else if(mode=="otu"){
    if(!is.null(access(filtered_phylo,"phy_tree"))) tree <- phy_tree(filtered_phylo) else tree <- NULL
    out_lst <- list(phylo=filtered_phylo, 
                    otu=as.data.frame(otu_table(filtered_phylo, T)), 
                    taxonomy=as.data.frame(tax_table(filtered_phylo)),
                    tree=tree)
    return(out_lst)
  }else{return(NULL)}
  
  
}

# load meta file and rename sample-column to 'SampleID' in df and to '#SampleID' in file
handleMetaFastqMode <- function(meta_file, fastq_sample_column, rm_spikes, sample_names){
  if (is.null(meta_file)){
    message("No meta-file uploaded. Using only sample names as meta-group")
    return (list(NULL, NULL))
  }
  meta <- read.csv(meta_file, header=T, sep="\t", check.names=F)
  
  # check for correct column names
  if (rm_spikes){
    if(!("total_weight_in_g" %in% colnames(meta))){stop(didNotFindWeightColumnError, call. = F)}
    if(!("amount_spike" %in% colnames(meta))){stop(didNotFindSpikeColumnError, call. = F)} 
  }
  if(!(fastq_sample_column %in% colnames(meta))){stop(didNotFindSampleColumnError, call. = F)}
  sample_column_idx <- which(colnames(meta)==fastq_sample_column)
  colnames(meta)[sample_column_idx] <- fastq_sample_column           # rename sample-column 
  if (sample_column_idx != 1) {meta <- meta[c(fastq_sample_column, setdiff(names(meta), fastq_sample_column))]}   # place sample-column at first position
  
  meta <- meta[, colSums(is.na(meta)) != nrow(meta)] # remove columns with only NA values
  rownames(meta)=meta[[fastq_sample_column]]
  
  # create second version of meta-file for rm_spikes.py -> needs first column to start with #
  if (rm_spikes){
    if(startsWith(colnames(meta)[1], "#")){
      meta_file_path = meta_file
    }else{
      colnames(meta)[1] <- "#SampleID"
      meta_file_path <- paste0(dirname(meta_file), "/meta_file_rm_spikes.tab")
      write.table(meta, meta_file_path, quote = F, sep="\t", row.names = F)
    }
  }else{
    meta_file_path = meta_file
  }
  
  message(paste0(Sys.time()," - Loaded meta file; colnames: "))
  message(paste(unlist(colnames(meta)), collapse = " "))
  return(list(meta=meta, meta_file_path=meta_file_path))
}

removeSpikes <- function(fastq_dir, meta_filepath, ncores){
  message("############ spike removal ############")
  
  #creating output dirs and files
  rm_spikes_outdir <- paste0(fastq_dir,"/rm_spikes_out")
  rm_spikes_outdir_woSpikes <- paste0(rm_spikes_outdir,"/omapping")
  rm_spikes_outdir_wSpikes <- paste0(rm_spikes_outdir,"/osamples")
  rm_spikes_spikes_file <- paste0(rm_spikes_outdir,"/ospikes.txt")
  rm_spikes_stats_file <- paste0(rm_spikes_outdir,"/ostats.txt")
  dir.create(rm_spikes_outdir)
  
  # running python script
  rm_spikes_command <- paste0("python3 src/rm_spikes.py data/spikes.fasta ",
                              meta_filepath,
                              " ", fastq_dir,
                              " ", rm_spikes_outdir_woSpikes,
                              " ", rm_spikes_outdir_wSpikes,
                              " ", rm_spikes_spikes_file,
                              " ", rm_spikes_stats_file,
                              " ", ncores)
  out<-system(rm_spikes_command, wait = T)
  message(paste0(Sys.time()," - ",out))
  
  message(paste0(Sys.time()," - Removed spikes: ",rm_spikes_command))
  message ("############ spike removal finished ############")
  return(rm_spikes_outdir_woSpikes)
}


##### write output-files/Rdata objects #####

writephyloseq<-function(phylo,path,fileprefix){
  otu<-as.data.frame(otu_table(phylo))
  meta<-as.data.frame(sample_data(phylo))
  taxa<-as.data.frame(tax_table(phylo))
  tree<-phy_tree(phylo)
  
  write.table(otu,paste0(path,fileprefix,"_otu.tsv"),quote = F,sep="\t")
  write.table(meta,paste0(path,fileprefix,"_meta.tsv"),quote = F,sep="\t")
  write.table(taxa,paste0(path,fileprefix,"_taxa.tsv"),quote = F,sep="\t")
  write.tree(tree,paste0(path,fileprefix,"_tree.tre"))
  
}

save_session <- function(dataset, name, filename){
  
}


##### other stuff #####

# handle differently encoded files
reading_makes_sense <- function(content_read) {
  out <- 
    (
      is.data.frame(content_read) &&
        nrow(content_read) > 0 &&
        ncol(content_read) > 0
    )
  
  return(out)
}

read_csv_custom <- function(file, file_type){
  try_encodings <- c("UTF-8","UTF-16LE")
  #testing the different encodings:
  out_lst <- lapply(try_encodings, function(x){
    message(paste0("Trying to read file with encoding: ", x))
    out <- NULL
    if(file_type=="meta"){out<-suppressWarnings(read.csv(file, header=TRUE, sep="\t", fileEncoding=x, check.names = F))}
    if(file_type=="otu"){out<-suppressWarnings(read.csv(file, header=TRUE, sep="\t", fileEncoding=x, check.names = F,row.names=1))}
    if(reading_makes_sense(out)){
      message(paste0(x,"-encoding resulted in useful output!"))
      return(out)
    }
  })
  out_tab <- NULL
  for (x in out_lst) {
    if(is.null(x)){next}else{out_tab<-x}
  }
  return(out_tab)
}