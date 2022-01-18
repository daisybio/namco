

##### sample dataset: ####

#load different testdata sets
observeEvent(input$upload_testdata, {
  dat <- read.csv("testdata/OTU_table.tab",header=T,sep="\t",row.names=1,check.names = F) # load otu table -> samples are columns
  taxonomy = generateTaxonomyTable(dat) # generate taxonomy table from TAX column
  dat = dat[!apply(is.na(dat)|dat=="",1,all),-ncol(dat)] # remove "empty" rows
  meta = read.csv("testdata/metafile.tab",header=T,sep="\t")
  rownames(meta) = meta[,1]
  meta = meta[match(colnames(dat),meta[[sample_column]]),]
  tree = ape::read.tree("testdata/tree.tre") # load phylogenetic tree
  
  #normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
  #normalized_dat = normalizeOTUTable(dat,normMethod)
  #tax_binning = taxBinning(normalized_dat[[2]],taxonomy)
  
  normalized_dat = list(norm_tab=dat, rel_tab = relAbundance(dat))
  
  #create phyloseq object from data (OTU, meta, taxonomic, tree)
  py.otu <- otu_table(normalized_dat$norm_tab,T)
  py.tax <- tax_table(as.matrix(taxonomy))
  py.meta <- sample_data(meta)
  phyloseq <- merge_phyloseq(py.otu,py.tax,py.meta,tree)
  
  #pre-build unifrac distance matrix
  if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalized_dat$norm_tab, tree) else unifrac_dist <- NULL
  
  # pre-calculate alpha-diversity
  alphaTabFull <- createAlphaTab(data.frame(phyloseq@otu_table, check.names=F), data.frame(phyloseq@sam_data, check.names = F))
  
  message(paste0(Sys.time()," - using Mueller sampledata "))
  
  #the final dataset
  dataset<- list(session_name="Mueller et al.",
                 rawData=dat,
                 metaData=meta,
                 taxonomy=taxonomy,
                 counts=NULL,
                 normalizedData=normalized_dat$norm_tab,
                 relativeData=normalized_dat$rel_tab,
                 tree=tree,
                 phylo=phyloseq,
                 phylo.raw=phyloseq,
                 unifrac_dist=unifrac_dist,
                 alpha_diversity=alphaTabFull,
                 undersampled_removed=F,
                 filtered=F, 
                 normMethod = 0,
                 is_fastq=F,
                 has_meta=T,
                 has_picrust=F,
                 is_sample_data=T,
                 is_restored=F,
                 has_rf=F,
                 has_diff_nw=F,
                 has_tax_nw=F,
                 has_comp_nw=F,
                 filterHistory="")
  
  vals$datasets[["Mueller et al."]] <- dataset
  updateTabItems(session,"sidebar", selected = "overview")
})
