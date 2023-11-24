####phylogenetic tree####

treeReactive <- reactive({
  if(!is.null(currentSet())){
    # prune taxa
    myTaxa = names(sort(taxa_sums(vals$datasets[[currentSet()]]$phylo), decreasing = TRUE)[1:input$phylo_prune])
    phy <- prune_taxa(myTaxa, vals$datasets[[currentSet()]]$phylo)
    
    meta <- as.data.frame(phy@sam_data)
    otu <- as.data.frame(otu_table(phy))
    taxonomy <- as.data.frame(tax_table(phy))
    if(!is.null(access(phy,"phy_tree"))) tree <- phy_tree(phy) else tree <- NULL
    if(!is.null(tree)){
      if(input$phylo_group != "NONE"){
        group <- input$phylo_group
        # count number of occurrences of the OTUs in each sample group
        l<-lapply(na.omit(unique(meta[[group]])), function(x){
          samples_in_group <- na.omit(meta[["SampleID"]][as.character(meta[[group]])==as.character(x)])
          d<-data.frame(otu[,samples_in_group])
          d<-data.frame(rowSums(apply(d,2,function(y) ifelse(y>0,1,0))))
          colnames(d) <- c(as.character(x))
          return(d)
        })
        info <- merge(data.frame(l, check.names = F), taxonomy,by.x=0, by.y=0)
        rownames(info) <- info$Row.names
        info$Row.names <- NULL
        group_cols <- suppressWarnings(which(colnames(info)%in%unique(meta[[group]])))
      }else{
        info <- taxonomy
        group_cols <- c()
      }
      if(input$phylo_taxonomy != "NONE"){
        # collect, which columns contain the info for the heatmap
        taxa_cols <- suppressWarnings(which(colnames(info)==input$phylo_taxonomy))
      }else{
        taxa_cols <- c()
      }
      
      tree_plot <- suppressWarnings(suppressMessages(ggtree::ggtree(tree, layout = input$phylo_method, branch.length = input$phylo_draw_clado) %<+% info +
                                                       geom_tiplab(size=input$phylo_size_tree,
                                                                   align = ifelse(input$phylo_edge_length=="Yes",T,F))))
      
      return(list(tree_plot = tree_plot,
                  info = info,
                  group_cols = group_cols,
                  taxa_cols = taxa_cols))
    }
  }
})


output$phyloTree <- renderPlot({
  if(!is.null(treeReactive())){
    h <- treeReactive()$tree_plot
    info <- treeReactive()$info
    if(!is.null(treeReactive()$group_cols)){
      h<-suppressWarnings(suppressMessages(ggtree::gheatmap(h, info[treeReactive()$group_cols], 
                                                            offset=input$phylo_offset,
                                                            color=NULL, 
                                                            width=input$phylo_width_meta,
                                                            colnames_position="top", 
                                                            colnames_angle=90, colnames_offset_y = 5, 
                                                            hjust=1, font.size=3,low="white")))
      h <- h + new_scale_fill()   
    }
    if(!is.null(treeReactive()$taxa_cols)){
      h<-suppressWarnings(suppressMessages(ggtree::gheatmap(h, info[treeReactive()$taxa_cols],
                                                            width=input$phylo_width_taxonomy,
                                                            offset=input$phylo_offset+5,
                                                            color="black",
                                                            colnames=T,
                                                            colnames_position="top",
                                                            colnames_angle=90, colnames_offset_y = 5,
                                                            hjust=1, font.size = 3))) 
      h <- h + scale_fill_manual(values=colorRampPalette(brewer.pal(9, input$namco_pallete))(nrow(unique(info[treeReactive()$taxa_cols]))))
    }
    h
  }
}, height = 1000)


#javascript show/hide toggle for advanced options
shinyjs::onclick("phylo_toggle_advanced",shinyjs::toggle(id="phylo_advanced",anim = T))





