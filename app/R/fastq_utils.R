##### fastq related functions #####

combineAndNormalize_dada2 <- function(seq_table, taxonomy, has_meta, meta, tree, sn, abundance_cutoff, apply_filter){
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
  
  # 0.25% relative abundance filtering
  # create relative abundance table first to find taxa to keep
  if(apply_filter){
    rel_otu_tmp <- relAbundance(t(as.data.frame(otu_table(phylo_unnormalized, F))))
    min <- apply(rel_otu_tmp, 2, function(x) ifelse(x>0.25, 1, 0))
    keep_taxa = names(which(rowSums(min)>0))
    phylo_unnormalized <- prune_taxa(keep_taxa, phylo_unnormalized)    
  }
  
  # give "normal" names to ASVs
  taxa_names(phylo_unnormalized) <- paste0("ASV", seq(ntaxa(phylo_unnormalized)))
  message(paste0("First phyloseq-object: ", ntaxa(phylo_unnormalized)))

  # create final objects with "real" ASV names
  asv_table_final <- t(as.data.frame(otu_table(phylo_unnormalized, F)))
  taxonomy_final <- as.data.frame(tax_table(phylo_unnormalized))
  if(has_meta){meta_final <- as.data.frame(sample_data(phylo_unnormalized))}else{meta_final <- NULL}
  refseq_final <- refseq(phylo_unnormalized)
  
  # normalization: per default to 10.000 reads
  # default normalization to 10.000 reads
  normMethod <- 4
  normalized_asv = normalizeOTUTable(asv_table_final, normMethod)
  
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
  if(!is.null(tree)){
    phylo <- merge_phyloseq(phylo, phy_tree(phylo_unnormalized))
    tree <- phy_tree(phylo)
  }else{
    tree <- NULL
  }
  
  message(paste0(Sys.time()," - final phyloseq-object: ", ntaxa(phylo)))
  print(phylo)
  
  return(list(phylo=phylo,
              normalized_asv=normalized_asv,
              raw_asv=asv_table_final,
              meta=meta_final,
              taxonomy=taxonomy_final, 
              tree=tree))
}

combineAndNormalize_lotus2 <- function(phylo, apply_filter, has_meta, sample_names, sample_column){
  
  # correct sample names (they are in reversed order in phyloseq object)
  sample_data(phylo)[[sample_column]] <- rev(sample_names)
  
  # 0.25% relative abundance filtering
  # create relative abundance table first to find taxa to keep
  if(apply_filter){
    rel_otu_tmp <- relAbundance(as.data.frame(otu_table(phylo)))
    min <- apply(rel_otu_tmp, 2, function(x) ifelse(x>0.25, 1, 0))
    keep_taxa = names(which(rowSums(min)>0))
    phylo <- prune_taxa(keep_taxa, phylo)    
  }
  
  raw_otu <- as.data.frame(otu_table(phylo))
  raw_meta <- as.data.frame(sample_data(phylo))
  raw_taxonomy <- as.data.frame(tax_table(phylo))
  raw_tree <- phy_tree(phylo)
  raw_refseq <- refseq(phylo)
  
  colnames(raw_taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # replace '?' in taxonomy assignment
  raw_taxonomy[,1] <- gsub("\\?","k__",raw_taxonomy[,1])  # For taxonomies related to kingdom level
  raw_taxonomy[,2] <- gsub("\\?","p__",raw_taxonomy[,2])  # For taxonomies related to phylum level
  raw_taxonomy[,3] <- gsub("\\?","c__",raw_taxonomy[,3])  # For taxonomies related to class level
  raw_taxonomy[,4] <- gsub("\\?","o__",raw_taxonomy[,4])  # For taxonomies related to order level
  raw_taxonomy[,5] <- gsub("\\?","f__",raw_taxonomy[,5])  # For taxonomies related to family level
  raw_taxonomy[,6] <- gsub("\\?","g__",raw_taxonomy[,6])  # For taxonomies related to genus level
  raw_taxonomy[,7] <- gsub("\\?","s__",raw_taxonomy[,7])  # For taxonomies related to species level
  
  normMethod <- 4
  normalized_asv = normalizeOTUTable(raw_otu, normMethod)
  
  # final object
  if(has_meta){
    phylo <- phyloseq(otu_table(normalized_asv$norm_tab, taxa_are_rows = T),
                      sample_data(raw_meta),
                      tax_table(as.matrix(raw_taxonomy)),
                      phy_tree(raw_tree),
                      raw_refseq)
  }else{
    phylo <- phyloseq(otu_table(normalized_asv$norm_tab, taxa_are_rows = T),
                      tax_table(as.matrix(raw_taxonomy)),
                      phy_tree(raw_tree),
                      raw_refseq)
    raw_meta <- NULL
  }
  
  message(paste0(Sys.time()," - final phyloseq-object: ", ntaxa(phylo)))
  print(phylo)
  
  return(list(phylo=phylo,
              normalized_asv=normalized_asv,
              raw_asv=raw_otu,
              meta=raw_meta,
              taxonomy=raw_taxonomy, 
              tree=raw_tree))
  
}

# load meta file and rename sample-column to 'SampleID' in df and to '#SampleID' in file
handleMetaFastqMode <- function(meta_file, fastq_sample_column, sample_column){
  
  if (is.null(meta_file)){
    message("No meta-file uploaded. Using only sample names as meta-group")
    return (list(NULL, NULL))
  }

  meta <- read_csv_custom(meta_file, "meta")
  
  # check for correct column names
  # if (rm_spikes){
  #   if(!("total_weight_in_g" %in% colnames(meta))){stop(didNotFindWeightColumnError, call. = F)}
  #   if(!("amount_spike" %in% colnames(meta))){stop(didNotFindSpikeColumnError, call. = F)} 
  # }
  
  if(!(fastq_sample_column %in% colnames(meta))){stop(didNotFindSampleColumnError, call. = F)}
  sample_column_idx <- which(colnames(meta)==fastq_sample_column)
  colnames(meta)[sample_column_idx] <- sample_column           # rename sample-column 
  if (sample_column_idx != 1) {meta <- meta[c(sample_column, setdiff(names(meta), sample_column))]}   # place sample-column at first position
  
  meta <- meta[, colSums(is.na(meta)) != nrow(meta)] # remove columns with only NA values
  rownames(meta)=meta[[sample_column]]
  
  # create second version of meta-file for rm_spikes.py -> needs first column to start with #
  # if (rm_spikes){
  #   if(startsWith(colnames(meta)[1], "#")){
  #     meta_file_path = meta_file
  #   }else{
  #     colnames(meta)[1] <- "#SampleID"
  #     meta_file_path <- paste0(dirname(meta_file), "/meta_file_rm_spikes.tab")
  #     write.table(meta, meta_file_path, quote = F, sep="\t", row.names = F)
  #   }
  # }else{
  #   meta_file_path = meta_file
  # }
  meta_file_path = meta_file
  
  message(paste0(Sys.time()," - Loaded meta file; colnames: "))
  message(paste(unlist(colnames(meta)), collapse = " "))
  return(list(meta=meta, meta_file_path=meta_file_path))
}

# not in use 
# TODO 
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

#function to unzip/untar files 
decompress <- function(dirname, is_paired){
  # get all files in dir, but no folders! (this excludes the fastqc_out folder)
  files <- setdiff(list.files(dirname, full.names = T), list.dirs(dirname, recursive = F, full.names = T))
  # only one file present -> handle it as compressed file
  if(length(files) == 1){
    compressed_file <- files[1]
    file_exts <- strsplit(compressed_file, split="\\.")[[1]]
    if(c("tar") %in% file_exts){
      fastq_files <- untar(compressed_file, list=T)
      untar(compressed_file, exdir=dirname)
    }else if(c("zip") %in% file_exts){
      fastq_files <- unzip(compressed_file, list=T)[["Name"]]
      unzip(compressed_file, exdir=dirname)
    }else{
      return(1)# no valid file extension/compression
    }
    if(length(fastq_files) == 0){
      return(1)# no fastq-files or no even number of fastq files in compressed file -> error
    }
    if(is_paired && (length(fastq_files) %% 2 != 0)){
      return(1) # no even amount of files uploaded for paired end experiment
    }
    unlink(compressed_file)
    return(0)
  # more than one file present -> no need for decompression
  }else if(length(files) > 1){
    if(is_paired && (length(files) %% 2 != 0)){
      return(1) 
    }else{
      return(0)  
    }
  }else{
    return(1) 
  }
}


# fastq_files: value from dataUpload 
# check fastq files for correct file name
# change temporary names to correct names
# return directory path with fastq files, fw files, rv files, sample names
handle_fastqs <- function(fastq_files, fastqc_exists, sampleNameCutoff, is_paired){

  dirname <- dirname(fastq_files$datapath[1]) 
  # if no fastQC has been run -> change filenames 
  if(!any(file.exists(paste0(dirname,"/",fastq_files$name)))){
    #files get "random" new filename in /tmp/ directory when uploaded in shiny -> change filename to the upload-name
    file.rename(from=fastq_files$datapath,to=paste0(dirname,"/",fastq_files$name)) 
  }
  
  #check file-type: if compressed file or multiple fastq-files
  outcome_decompress <- decompress(dirname, is_paired)
  if(outcome_decompress == 1){stop(errorDuringDecompression, call. =F)}
  
  # collect fw & rv files 
  foreward_files <- sort(list.files(dirname, pattern = "_R1_001.fastq", full.names = T))
  reverse_files <- sort(list.files(dirname, pattern = "_R2_001.fastq", full.names = T))
  
  if(!all(grepl(sampleNameCutoff, basename(foreward_files)))){stop(sampleNameCutoffNotPresent, call. =F)}
  # get correct sample names
  sample_names <- sapply(strsplit(basename(foreward_files), sampleNameCutoff), `[`, 1) # sample name: everything until cutoff
  
  # checks of fastq files
  if((length(foreward_files) == 0) || length(reverse_files) == 0){stop(noFilesWithCorrectExtensionFoundError, call.=F)}
  if (is_paired && (length(foreward_files) != length(reverse_files))){stop(noEqualFastqPairsError, call.=F)}

  
  return(list(dirname=dirname, 
              foreward_files=foreward_files,
              reverse_files=reverse_files,
              sample_names=sample_names))
}
