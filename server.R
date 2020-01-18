library(data.table)
library(DT)
library(GUniFrac)
library(plotly)
library(RColorBrewer)
library(reshape2)
library(shiny)
library(textshape)
library(tidyr)
library(visNetwork)
packages <-c("ade4","GUniFrac","phangorn","cluster","fpc") 

server <- function(input,output,session){
  options(shiny.maxRequestSize=1000*1024^2,stringsAsFactors=F)  # upload up to 1GB of data
  vals = reactiveValues(datasets=list()) # reactiveValues is a container for variables that might change during runtime and that influence one or more outputs, e.g. the currently selected dataset
  currentSet = NULL # a pointer to the currently selected dataset
  minLibSize = 1000 # a lower boundary for the expected library size; if in a bait less than this number of proteins are quantified, a warning message will be shown
  lastplot = NULL
  source("algorithms.R")
  
  ################################################################################################################################################
  #
  # file upload
  #
  ################################################################################################################################################
  normalizeOTUTable <- function(tab,method=0){
    min_sum <- min(colSums(tab))
    
    if(method==0){
      # Divide each value by the sum of the sample and multiply by the minimal sample sum
      norm_tab <- t(min_sum*t(tab)/colSums(tab))
    } else{
      # Rarefy the OTU table to an equal sequencing depth
      norm_tab <- Rarefy(t(tab),depth=min_sum)
      norm_tab <- t(as.data.frame(norm_tab$otu.tab.rff))
    }
    
    # Calculate relative abundances for all OTUs over all samples
    # Divide each value by the sum of the sample and multiply by 100
    rel_tab <- t(100*t(tab)/colSums(tab))
    
    return(list(norm_tab=norm_tab,rel_tab=rel_tab))
  } 
  
  # Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
  uploadModal <- function(failed=F) {
    modalDialog(
      p("Please provide an OTU table."),
      fileInput("otuFile","Select OTU table",width="100%"),
      fileInput("metaFile","Select Metadata File",width="100%"),
      textInput("dataName","Enter a project name:",placeholder="New_Project",value="New_Project"),
      
      if(failed) div(tags$b("The file you specified could not be loaded.",style="color: red;")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("upload_ok","OK")
      )
    )
  }
  
  # launch upload dialog
  observeEvent(input$upload, {
    showModal(uploadModal())
  })
  
  # try to load the dataset specified in the dialog window
  observeEvent(input$upload_ok, {
    # Check that data object exists and is data frame.
    if(!is.null(input$otuFile)&!input$dataName%in%names(vals$datasets)){
      tryCatch({
        dat <- read.csv(input$otuFile$datapath,header=T,sep="\t",row.names=1)
        dat = dat[!apply(is.na(dat)|dat=="",1,all),]
        if(is.null(input$metaFile)) meta = NULL else meta = read.csv(input$metaFile$datapath,header=T,sep="\t",row.names=1)
        normalized_dat = normalizeOTUTable(dat,0)
        
        vals$datasets[[input$dataName]] <- list(rawData=dat,metaData=meta,counts=NULL,normalizedData=normalized_dat$norm_tab,relativeData=normalized_dat$rel_tab)
        removeModal()
      },
      error = function(e){
        showModal(uploadModal(failed=T))
      })
    } else{
      showModal(uploadModal(failed=T))
    }
  })
  
  # update datatable holding currently loaded datasets
  output$datasets <- renderDataTable({
    if(length(vals$datasets)!=0){
      datatable(data.frame(Datasets=names(vals$datasets)),rownames=F,options=list(pageLength=10,dom="t"),selection=list(mode="single",selected=length(vals$datasets))) %>%
        formatStyle("Datasets",color="white",backgroundColor="#222D33")
    } else{
      datatable(data.frame(Datasets=""),rownames=F,options=list(pageLength=10,dom="t"),selection=list(mode="single")) %>%
        formatStyle("Datasets",color="white",backgroundColor="#222D33")
    }
  })
  currentSet <- eventReactive(input$datasets_rows_selected, {return(input$datasets_rows_selected)})
  
  # update input selections
  observe({
    
  })
  
  
  ################################################################################################################################################
  #
  # QC plots
  #
  ################################################################################################################################################
  # update targets table of the currently loaded dataset
  output$metaTable <- renderDataTable({
    if(!is.null(currentSet())) datatable(vals$datasets[[currentSet()]]$metaData,filter='bottom',options=list(searching=T,pageLength=20,dom="Blfrtip"),editable=T)
    else datatable(data.frame(),options=list(dom="t"))
  })
  
  # Plot rarefaction curves
  output$rarefacCurve <- renderPlotly({
    tab = as.matrix(vals$datasets[[currentSet()]]$rawData)
    
    # determine data points for rarefaction curve
    rarefactionCurve = lapply(1:ncol(tab),function(i){
      n = seq(1,colSums(tab)[i],by=5000)
      if(n[length(n)]!=colSums(tab)[i]) n=c(n,colSums(tab)[i])
      drop(rarefy(t(tab[i,]),n))
    })
    slope = apply(tab,2,function(x) rareslope(x,sum(x)-100))
    
    first = order(slope,decreasing=T)[1]
    p <- plot_ly(x=attr(rarefactionCurve[[first]],"Subsample"),y=rarefactionCurve[[first]],text=paste0(rownames(tab)[first],"; slope: ",round(1e5*slope[first],3),"e-5"),hoverinfo="text",color="high",type="scatter",mode="lines",colors=c("red","black"))
    for(i in order(slope,decreasing=T)[2:input$rareToShow]){
      highslope = as.numeric(slope[i]>quantile(slope,1-input$rareToHighlight/100))+1
      p <- p %>% add_trace(x=attr(rarefactionCurve[[i]],"Subsample"),y=rarefactionCurve[[i]],text=paste0(rownames(tab)[i],"; slope: ",round(1e5*slope[i],3),"e-5"),hoverinfo="text",color=c("low","high")[highslope],showlegend=F)
    }
    p %>% layout(title="Rarefaction Curves",xaxis=list(title="Number of Reads"),yaxis=list(title="Number of Species"))
  })
  
  # calculate various measures of alpha diversity
  alphaDiv <- function(x,method){
    switch(method,
      richness = {
        # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
        count = sum(x>0)
        return(count)},
      shannon = {
        se = -sum(x[x>0]/sum(x)*log(x[x>0]/sum(x)))
        return(se)
      },
      shannon_eff = {
        se = round(exp(-sum(x[x>0]/sum(x)*log(x[x>0]/sum(x)))),digits =2)
        return(se)
      },
      simpson = {
        si = sum((x[x>0]/sum(x))^2)
        return(si)
      },
      simpson_eff = {
        si = round(1/sum((x[x>0]/sum(x))^2),digits =2)
        return(si)
      }
    )
  }
  
  # plot distribution of taxa
  output$taxaDistribution <- renderPlotly({
    if(!is.null(currentSet())){
      data = vals$datasets[[currentSet()]]$rawData
      data = apply(data,2,function(x) x/sum(x))
      
      taxa <- ifelse(rowSums(data)/ncol(data)<(input$otherCutoff/100),"Other",rownames(data))
      other <- data[which(taxa=="Other"),]
      other <- colSums(other)
      
      data <- data[-which(taxa=="Other"),]
      data <- rbind(data,Other=other)
      data = 100*data
      
      data = data[,1:10]
      
      p <- plot_ly(x=colnames(data),y=data[1,],name=rownames(data)[1],type="bar") #%>%
      for(i in 2:nrow(data)) p <- p %>% add_trace(y=data[i,],name=rownames(data)[i])
      p %>% layout(yaxis=list(title='Count [%]'),barmode='stack')
    } else plotly_empty()
  })
  
  # make PCA plots
  output$pcaPlot <- renderPlotly({
    if(!is.null(currentSet())){
      samples = rownames(vals$datasets[[currentSet()]]$rawData)
      mat = t(vals$datasets[[currentSet()]]$rawData)
      mode = input$pcaMode
      
      pca = prcomp(mat,center=T,scale=T)
      out = data.frame(pca$rotation,txt=samples)
      percentage = signif(pca$sdev^2/sum(pca$sdev^2)*100,2)
      
      if(mode=="2D"){
        plot_ly(out,x=~PC1,y=~PC2,text=~txt,hoverinfo='text',type='scatter',mode="markers") %>%
          layout(title="PCA (2D)",xaxis=list(title=paste0("PC1 (",percentage[1]," % variance)")),yaxis=list(title=paste0("PC2 (",percentage[2]," % variance)")))
      } else{
        plot_ly(out,x=~PC1,y=~PC2,z=~PC3,text=~txt,hoverinfo='text',type='scatter3d',mode="markers") %>%
          layout(title="PCA (3D)",scene=list(xaxis=list(title=paste0("PC1 (",percentage[1]," % variance)")),yaxis=list(title=paste0("PC2 (",percentage[2]," % variance)")),zaxis=list(title=paste0("PC3 (",percentage[3]," % variance)"))))
      }
    } else{
      plotly_empty()
    }
  })
  
  # plot PCA loadings
  output$loadingsPlot <- renderPlotly({
    if(!is.null(currentSet())){
      phyla = colnames(vals$datasets[[currentSet()]]$rawData)
      mat = vals$datasets[[currentSet()]]$rawData
      
      pca = prcomp(mat,center=T,scale=T)
      loadings = data.frame(Phyla=phyla,loading=pca$rotation[,as.numeric(input$pcaLoading)])
      loadings = loadings[order(loadings$loading,decreasing=T),][c(1:10,(nrow(loadings)-9):nrow(loadings)),]
      loadings$Phyla = factor(loadings$Phyla,levels =loadings$Phyla)
      
      plot_ly(loadings,x=~Phyla,y=~loading,text=~Phyla,hoverinfo='text',type='bar',color=I(ifelse(loadings$loading>0,"blue","red"))) %>%
        layout(title="Top and Bottom Loadings",xaxis=list(title="",zeroline=F,showline=F,showticklabels=F,showgrid=F),yaxis=list(title=paste0("loadings on PC",input$pcaLoading)),showlegend=F) %>% hide_colorbar()
    } else{
      plotly_empty()
    }
  })
  
  output$alphaPlot <- renderPlotly({
    if(!is.null(currentSet())){
      tab = vals$datasets[[currentSet()]]$rawData
      out = data.frame(Sample = rownames(tab))
      for(i in c("richness","shannon","shannon_eff","simpson","simpson_eff")){
        div = apply(tab,1,function(x) alphaDiv(x,i))
        out = cbind(out,div)
      }
      colnames(out) = c("Sample","richness","shannon","shannon_eff","simpson","simpson_eff")
      
      
    }
    else plotly_empty()
  })
  
  
  output$betaTable <- renderDataTable({
      
  })
  
  
  ################################################################################################################################################
  #
  # Network analysis
  #
  ################################################################################################################################################
  
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
  
  output$cutoffHist <- renderPlotly({
    
  })
  
  observe({
    if(is.null(vals$datasets[[currentSet()]]$counts)) vals$datasets[[currentSet()]]$counts = basic_approach(vals$datasets[[currentSet()]]$rawData)
  })
  
  # visualize network for basic approach
  output$basicNetwork <- renderVisNetwork({
    if(!is.null(currentSet())&!is.null(vals$datasets[[currentSet()]]$counts)){
      data = vals$datasets[[currentSet()]]$rawData
      counts = vals$datasets[[currentSet()]]$counts
      
      nodes <- data.frame(id=colnames(data),label=paste0("OTU_",1:ncol(data)))[1:10,]
      edges <- setNames(counts,c("from","to","label"))
      edges$label = round(edges$label,2)
      edges$color = ifelse(edges$label>0,"green","red")
      edges$width = round((edges$label-min(edges$label))/(max(edges$label)-min(edges$label))*9)
      edges = edges[edges$from%in%nodes$id|edges$to%in%nodes$id,]
      
      v = visNetwork(nodes,edges)
    } else v = visNetwork(data.frame(id="",label=""))
    v %>% visEdges(color=list(color="lightblue"))
  })
}