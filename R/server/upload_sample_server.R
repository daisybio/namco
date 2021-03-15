uploadTestdataModal <- function(failed=F, error_message=NULL){
  modalDialog(
    title="Choose a sample-Dataset",
    HTML("<h5>[For details of each set, look into the <b>Info & Settings</b> tab!]</h5>"),
    hr(),
    fluidRow(column(1),
             column(10, selectInput("selectTestdata", shiny::HTML("<p><span style='color: green'>Select Sample-Dataset</span></p>"),choices = c("Mueller et al. (Mice samples)"))),
             column(1)
    ),
    fluidRow(
      column(10,radioButtons("normMethod","Normalization Method",c("no Normalization","by Sampling Depth","by Rarefaction"),inline=T))
    ),
    footer = tagList(
      modalButton("Cancel", icon = icon("times-circle")),
      actionButton("upload_testdata_ok","OK",style="background-color:blue; color:white")
    ),
    easyClose = T, fade = T, size = "l",
  )
}




##### sample datasets: ####

#load specified testdataset
observeEvent(input$upload_testdata, {
  showModal(uploadTestdataModal())
})

#load different testdata sets
observeEvent(input$upload_testdata_ok, {
  if(input$selectTestdata == "Mueller et al. (Mice samples)"){
    dat <- read.csv("testdata/OTU_table.tab",header=T,sep="\t",row.names=1,check.names = F) # load otu table -> samples are columns
    taxonomy = generateTaxonomyTable(dat) # generate taxonomy table from TAX column
    dat = dat[!apply(is.na(dat)|dat=="",1,all),-ncol(dat)] # remove "empty" rows
    meta = read.csv("testdata/metafile.tab",header=T,sep="\t")
    rownames(meta) = meta[,1]
    meta = meta[match(colnames(dat),meta[[sample_column]]),]
    tree = read.tree("testdata/tree.tre") # load phylogenetic tree
    
    normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
    normalized_dat = normalizeOTUTable(dat,normMethod)
    #tax_binning = taxBinning(normalized_dat[[2]],taxonomy)
    
    #create phyloseq object from data (OTU, meta, taxonomic, tree)
    py.otu <- otu_table(normalized_dat$norm_tab,T)
    py.tax <- tax_table(as.matrix(taxonomy))
    py.meta <- sample_data(meta)
    phylo <- merge_phyloseq(py.otu,py.tax,py.meta,tree)
    
    #pre-build unifrac distance matrix
    if(!is.null(tree)) unifrac_dist <- buildGUniFracMatrix(normalized_dat$norm_tab,meta,tree) else unifrac_dist <- NULL
    
    message(paste0(Sys.time()," - using Mueller sampledata "))
    
    #the final dataset
    dataset<- list(rawData=dat,
                   metaData=meta,
                   taxonomy=taxonomy,
                   counts=NULL,
                   normalizedData=normalized_dat$norm_tab,
                   relativeData=normalized_dat$rel_tab,
                   tree=tree,
                   phylo=phylo,
                   unifrac_dist=unifrac_dist,
                   undersampled_removed=F,
                   filtered=F, 
                   normMethod = normMethod,
                   is_fastq=F,
                   has_meta=T)
    
    vals$datasets[["Mueller et al."]] <- dataset
    updateTabItems(session,"sidebar")
    removeModal()
    showModal(finishedFastqUploadModal)
  }
  if(input$selectTestdata == "Global Patterns (environmental samples)"){
    gp <- readRDS("testdata/GlobalPatterns")
    dat <- as.data.frame(otu_table(gp))
    taxonomy <- as.data.frame(tax_table(gp))
    taxonomy <- addMissingTaxa(taxonomy)
    taxonomy[is.na(taxonomy)] <- "NA"
    meta <- as.data.frame(sample_data(gp))
    colnames(meta)[1]<-"SampleID"
    tree <- phy_tree(gp)
    
    normMethod = which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1
    normalized_dat = normalizeOTUTable(dat,normMethod)
    #tax_binning = taxBinning(normalized_dat[[2]],taxonomy)
    
    phylo_gp <- merge_phyloseq(otu_table(normalized_dat$norm_tab,T),tax_table(as.matrix(taxonomy)),sample_data(meta),tree)
    unifrac_dist <- readRDS("testdata/GlobalPatternsUnifrac")
    
    #the final dataset 
    dataset <- list(rawData=dat,metaData=meta,taxonomy=taxonomy,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab,tree=tree,phylo=phylo_gp,unifrac_dist=unifrac_dist,undersampled_removed=F,filtered=F, normMethod = normMethod)
    
    vals$datasets[["GlobalPatterns"]] <- dataset
    updateTabItems(session,"sidebar")
    removeModal()
  }
  if(input$selectTestdata == "Enterotype (facial samples)"){
    ent <- readRDS("testdata/enterotype")
    dat <- as.data.frame(otu_table(ent))
    rows <- rownames(dat)
    dat <- data.frame(lapply(dat, function(x){return(as.integer(x*10000))}))
    rownames(dat) <- rows
    taxonomy <- data.frame(tax_table(ent))
    taxonomy <- addMissingTaxa(taxonomy)
    taxonomy[is.na(taxonomy)] <- "NA"
    meta <- as.data.frame(sample_data(ent))
    names(meta)[names(meta) == "SampleID"] <- "SampleIDwNAs"
    names(meta)[names(meta) == "Sample_ID"] <- "SampleID"
    meta <- meta[,c(2,4,1,3,5,6,7,8,9)]
    
    tree <- NULL
    
    normalized_dat = normalizeOTUTable(dat,which(input$normMethod==c("no Normalization","by Sampling Depth","by Rarefaction","centered log-ratio"))-1)
    
    phylo_ent <- ent
    unifrac_dist <- NULL
    
    #the final datatset
    dataset <- list(rawData=dat,metaData=meta,taxonomy=taxonomy,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab,tree=tree,phylo=phylo_ent,unifrac_dist=unifrac_dist,undersampled_removed=F,filtered=F)
    
    vals$datasets[["Enterotype"]] <- dataset
    updateTabItems(session,"sidebar")
    removeModal()
  }
})
