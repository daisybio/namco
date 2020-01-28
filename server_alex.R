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

#TODO: check for all the libraries -> themetagenomics!

server <- function(input,output,session){
  options(shiny.maxRequestSize=1000*1024^2,stringsAsFactors=F)  # upload up to 1GB of data
  vals = reactiveValues(datasets=list()) # reactiveValues is a container for variables that might change during runtime and that influence one or more outputs, e.g. the currently selected dataset
  currentSet = NULL # a pointer to the currently selected dataset
  minLibSize = 1000 # a lower boundary for the expected library size; if in a bait less than this number of proteins are quantified, a warning message will be shown
  lastplot = NULL
  source("algorithms.R")
  source("themetagenomics/themetagenomics.R")
  source("themetagenomics/formatting.R")
  source("themetagenomics/utils.R")
  source("themetagenomics/predict.R")
  source("themetagenomics/picrust.R")
  source("themetagenomics/RcppExports.R")
  source("gene.table.R")
  
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
  
  # update input selections based on dataset
  observe({
    
    #update silder for binarization cutoff dynamically based on normalized dataset
    min_value <- min(vals$datasets[[currentSet()]]$normalizedData)
    max_value <- round(max(vals$datasets[[currentSet()]]$normalizedData)/16)
    updateNumericInput(session,"binCutoff",value = 1,min=min_value,max=max_value)
    updateNumericInput(session,"k_in",value=0,min=0,max=vals$datasets[[currentSet()]]$vis_out$K,step=1)
    
    ########updates based on meta info########
    
    covariates <- vals$datasets[[currentSet()]]$vis_out$covariates
    updateSelectInput(session, "choose",choices = covariates)
    
    #pick all column names, exepct the SampleID
    group_columns <- colnames(subset(vals$datasets[[currentSet()]]$metaData, select= -c(SampleID)))
    updateSelectInput(session,"groupCol",choices = group_columns)
    
    updateSelectInput(session,"formula",choices = group_columns)
    
  })
  observe({
    ref_choices <- unique(vals$datasets[[currentSet()]]$metaData[[input$formula]])
    updateSelectInput(session, "refs",choices=ref_choices)
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
      samples =colnames(vals$datasets[[currentSet()]]$rawData)
      mat = vals$datasets[[currentSet()]]$rawData
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
      phyla = rownames(vals$datasets[[currentSet()]]$rawData)
      mat = t(vals$datasets[[currentSet()]]$rawData)
      
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
  
  
  
  
  
  
  
  ################################################################################################################################################
  #
  # Network analysis
  #
  ################################################################################################################################################
  
  #show histogram of all OTU values -> user can pick cutoff for binarization here
  output$cutoffHist <- renderPlotly({
    if(!is.null(currentSet())){
      dat <- log(as.data.frame(vals$datasets[[currentSet()]]$normalizedData))
      
      plot_ly(x=unlist(dat),type="histogram")%>%
        layout(xaxis = list(title="log(normalized OTU-values)"), yaxis = list(title="Frequency"),
               shapes=list(list(type="line",y0=0,y1=1,yref="paper",x0=log(input$binCutoff),x1=log(input$binCutoff),line=list(color="black",width=2))))
    }else{
      plotly_empty()
    }
  })
  
  output$boolHeat <- renderPlotly({
    if(!is.null(currentSet())){
      
      cutoff <- input$binCutoff
      df <- vals$datasets[[currentSet()]]$normalizedData
      
      m <- as.matrix(df)
      m<-apply(m,c(1,2),function(x){ifelse(x<cutoff,0,1)})
      
      melted_m<-melt(m)
      melted_m$Var1 <- as.factor(melted_m$Var1)
      melted_m$Var2 <- as.factor(melted_m$Var2)
      
      plotly::plot_ly(z=melted_m$value,x=melted_m$Var1,y=melted_m$Var2,type="heatmap")
      
    }else{
      plotly_empty()
    }
  })
  
  #check if button for new calculation of counts is clicked -> reload network with the new counts
  observeEvent(input$startCalc,{
    withProgress(message = 'Calculating Counts..', value = 0, {
      vals$datasets[[currentSet()]]$counts = generate_counts(OTU_table=vals$datasets[[currentSet()]]$normalizedData,
                                                             meta = vals$datasets[[currentSet()]]$metaData,
                                                             group_column = input$groupCol,
                                                             cutoff = input$binCutoff,
                                                             fc = ifelse(input$useFC=="log2(fold-change)",TRUE,FALSE),
                                                             progress = T)
      View(vals$datasets[[currentSet()]]$counts)
    })
    
    # visualize network for basic approach
    output$basicNetwork <- renderVisNetwork({
      #if(!is.null(currentSet())&!is.null(vals$datasets[[currentSet()]]$counts)){
      
      data = vals$datasets[[currentSet()]]$rawData
      counts = vals$datasets[[currentSet()]]$counts
      
      nodes <- data.frame(id=colnames(data),label=paste0("OTU_",1:ncol(data)))[1:50,]
      edges <- setNames(counts,c("from","to","label"))
      edges$label = round(edges$label,2)
      edges$color = ifelse(edges$label>0,"green","red")
      edges$width = round((edges$label-min(edges$label))/(max(edges$label)-min(edges$label))*9)
      edges = edges[edges$from%in%nodes$id|edges$to%in%nodes$id,]
      
      v = visNetwork(nodes,edges)
      
      #else v = visNetwork(data.frame(id="",label=""))
      #v %>% visEdges(color=list(color="lightblue"))
    })
  })
  
  output$countDistr <- renderPlotly({
    if(!is.null(vals$datasets[[currentSet()]]$counts)){
      plot_ly(x=vals$datasets[[currentSet()]]$counts$value,type="histogram")%>%
        layout(xaxis = list(title="log(count-values)"), yaxis = list(title="Frequency"))
    }else{
      plotly_empty()
    }
  })
  
  ##############################
  #  themetagenomcis apporach  #
  ##############################
  
  observeEvent(input$themeta,{
    
    withProgress(message = 'Calculating Topics..', value = 0,{
      if(!is.null(currentSet())){
        otu <- vals$datasets[[currentSet()]]$rawData
        meta <- vals$datasets[[currentSet()]]$metaData
        #tax <- vals$datasets[[currentSet()]]$taxonomy
        tax <- GEVERS$TAX
        
        incProgress(1/7,message = "preparing OTU data..")
        formula <- as.formula(paste0("~ ",input$formula))
        formula_char <- input$formula
        refs <- input$refs
        
        obj <- themetagenomics::prepare_data(otu_table = otu,
                                             rows_are_taxa = T,
                                             tax_table = tax,
                                             metadata = meta,
                                             formula=formula,
                                             refs = refs)
        
        incProgress(1/7,message = "finding topics..")
        K=input$K
        sigma_prior = input$sigma_prior
        topics_obj <- themetagenomics::find_topics(themetadata_object=obj,
                                                   K=K,
                                                   sigma_prior = sigma_prior)
        
        incProgress(1/7,message="predicting topic functions..")
        functions_obj <- predict.topics(topics_obj,reference_path = "themetagenomics/")
        
        incProgress(1/7,message = "generate gene.table")
        topic_function_table <- gene.table(functions_obj)
        
        incProgress(1/7,message = "finding topic effects..")
        #measure relationship of covarite with samples over topics distribution from the STM
        topic_effects_obj <- themetagenomics::est(topics_obj)
        
        
        #function_effects <- themetagenomics::est(functions_obj,topics_subset=3)
        
        class(topics_obj) <- "topics"
        incProgress(1/7,message = "preparing visualization..")
        vals$datasets[[currentSet()]]$vis_out <- prepare_vis(topic_effects_obj)
        vals$datasets[[currentSet()]]$vis_out$K <- K
        vals$datasets[[currentSet()]]$vis_out$sigma <- sigma_prior
        vals$datasets[[currentSet()]]$vis_out$formula <- formula_char
        vals$datasets[[currentSet()]]$vis_out$refs <- refs
        vals$datasets[[currentSet()]]$topic_effects <- topic_effects_obj
        vals$datasets[[currentSet()]]$gene_table <- topic_function_table
        incProgress(1/7)
      }
    })
  })
  
  #make gene.table download-able
  output$downloadGeneTable <- downloadHandler(
    filename = "gene_table.csv",
    content = function(file){
      write.csv(vals$datasets[[currentSet()]]$gene_table,file,row.names = T)
    }
  )
  
  
  REL <- reactive({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      if (show_topic$k != 0){
        
        current_k <- paste0('T',show_topic$k)
        
        l <- input$lambda
        
        tinfo_k <- vis_out$tinfo[vis_out$tinfo$Category == current_k,]  #subset(tinfo,Category == current_k)
        rel_k <- l*tinfo_k$logprob + (1-l)*tinfo_k$loglift
        new_order <- tinfo_k[order(rel_k,decreasing=TRUE)[1:vis_out$taxa_bar_n],]
        new_order$Term <- as.character.factor(new_order$Term)
        new_order$Taxon <- vis_out$taxa[new_order$Term,input$taxon]
        new_order$Term <- factor(new_order$Term,levels=rev(new_order$Term),ordered=TRUE)
        
      }else{
        
        new_order <- vis_out$default
        new_order$Taxon <- vis_out$taxa[as.character.factor(new_order$Term),input$taxon]
        
      }
      
      new_order
    }
    
  })
  
  EST <- reactive({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    topic_effects <- vals$datasets[[currentSet()]]$topic_effects$topic_effects
    
    if(!is.null(vis_out) & input$choose != "Please start calculation above first!"){
      suppressWarnings({
        
        covariate <- input$choose
        
        est_mat <- topic_effects[[covariate]]$est
        
        df0 <- data.frame(topic=paste0('T',1:vis_out$K),
                          est=est_mat[,1],
                          lower=est_mat[,2],
                          upper=est_mat[,3],
                          sig=ifelse(1:vis_out$K %in% topic_effects[[covariate]]$sig,'1','0'))
        df0$sig <- factor(as.character(sign(df0$est) * as.numeric(as.character(df0$sig))),levels=c('0','1','-1'),ordered=TRUE)
        df <- df0[order(topic_effects[[covariate]][['rank']]),]
        df$topic <- factor(df$topic,levels=df$topic,ordered=TRUE)
        
        p_est <- ggplot(df,aes_(~topic,y=~est,ymin=~lower,ymax=~upper,color=~sig)) +
          geom_hline(yintercept=0,linetype=3)
        p_est <- p_est +
          geom_pointrange(size=2) +
          theme_minimal() +
          labs(x='',y='Estimate') +
          scale_color_manual(values=c('gray','indianred3','dodgerblue3'),drop=FALSE) +
          scale_fill_brewer(type='qual',palette=6,direction=-1) +
          theme(legend.position='none',
                axis.text.x=element_text(angle=-90,hjust=0,vjust=.5))
        
        list(p_est=p_est,k_levels=levels(df$topic),covariate=covariate,df0=df0)
        
      })
    }
  })
  
  output$est <- renderPlotly({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      ggplotly(EST()$p_est,source='est_hover',tooltip=c('topic','est','lower','upper'))
    }else{
      plotly_empty()
    }
    
    
  })
  
  output$text1 <- renderUI({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      HTML(sprintf("Below are the results of a %s topic STM. The ordination of the topics over taxa distribution (left) and the frequencies of
                    the top %s taxa (in terms of saliency) across all topics. By selecting a topic, the relative
                    frequencies of the taxa within that topic are shown in red. The ordination figure can be shown in
                    either 2D or 3D and the ordination method can be adjusted. Lambda adjusts the relevance calculation.
                    Choosing the taxon adjusts the group coloring for the bar plot. Clicking Reset resets the topic selection.",
                   vis_out$K,vis_out$taxa_bar_n))
    }
  })
  
  output$text2 <- renderUI({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      HTML(paste0('Below shows topic-to-topic correlations from the samples over topics distribution. The edges represent positive',
                  ' correlation between two topics, with the size of the edge reflecting to the magnitude of the correlation.',
                  ' The size of the nodes are consistent with the ordination figure, reflecting the marginal topic frequencies.'))
    }
  })
  
  output$text3 <- renderUI({
    HTML(paste0("Integrates the samples over topics p(s|k) and topics over taxa p(k|t) distributions from the STM, the topic correlations from the p(s|k) component, the covariate effects from the p(s|k) component, and their relationship with the raw taxonomic abundances. The covariate effects for each topic are shown as a scatterplot of posterior weights with error bars corresponding the global approximation of uncertainty. If the covariate chosen is binary, then the posterior regression weights with uncertainty intervals are shown. This is analogous to the mean difference between factor levels in the posterior predictive distribution. For continuous covariates, the points again represent the mean regression weights (i.e., the posterior slope estimate of the covariate). If, however, a spline or polynomial expansion was used, then the figure shows the posterior estimates of the standard deviation of the predicted topic probabilities from the posterior predictive distribution. Colors indicate whether a given point was positive (red) or negative (blue) and did not enclose 0 at a user defined uncertainty interval.",
                "The ordination figure maintains the color coding just decribed. The ordination is performed on p(k|t) via either PCoA (using either Jensen-Shannon, Euclidean, Hellinger, Bray-Curtis, Jaccard, or Chi-squared distance) or t-SNE. The latter iterates through decreasing perplexity values (starting at 30) until the algorithm succeeds. The top 2 or 3 axes can be shown. The radius of the topic points corresponds to the topic frequencies marginalized over taxa.",
                "The bar plot behaves in accordance with LDAvis. When no topics are chosen, the overall taxa frequencies are shown. These frequencies do not equal the abundances found in the initial abundance table. Instead, they show p(k|t) multiplied by the marginal topic distribution (in counts). To determine the initial order in which taxa are shown, these two distributions are compared via Kullback-Liebler divergence and then weighted by the overall taxa frequency. The coloration of the bars indiciates the taxonomic group the inidividual taxa belong to. The groups shown are determined based on the abundance of that group in the raw abundance table. When a topic is selected, the relative frequency of a given taxa in that topic is shown in red.",
                "λ controls relevance of taxa within a topic, which in turn is used to adjust the order in which the taxa are shown when a topic is selected. Relevence is essentially a weighted sum between the probability of taxa in a given topic and the probability of taxa in a given topic relative to the overall frequency of that taxa. Adjusting λ influences the relative weighting such that",
                "r = λ x log p(t|k) + λ x log p(t|k)/p(x)",
                "The correlation graph shows the topic correlations from p(s|k) ~ MVN(mu,sigma). Again, the coloration described above is conserved. The size of the nodes reflects the magnitude of the covariate posterior regression weight, whereas the width of the edges represents the value of the positive correlation between the connected nodes. By default, the graph estimates are determined using the the huge package, which first performs a nonparanormal transformation of p(s|k), followed by a Meinhuasen and Buhlman procedure. Alternatively, by choosing the simple method, the correlations are simply a thresholded MAP estimate of p(s|k)."))
  })
  
  output$input_variables <- renderUI({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      K <- vis_out$K
      sigma <- vis_out$sigma_prior
      formula<-vis_out$formula
      refs<-paste(vis_out$refs,collapse=", ")
      HTML(paste0("Number of chosen topics (K): ",K,"  Value of sigma_prior: ",sigma,"  Formula: ",formula,"  Reference Level: ",refs))
    }
    
  })
  
  output$ord <- renderPlotly({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out) & !is.null(EST())){
      beta <- t(vis_out$beta)
      
      if (input$dist == 'hellinger'){
        
        d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'norm'),method='euclidean'),3,eig=TRUE)
        
      }else if (input$dist == 'chi2'){
        
        d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'chi.square'),method='euclidean'),3,eig=TRUE)
        
      }else if (input$dist == 'jsd'){
        
        d <- cmdscale(proxy::dist(beta,jsd),3,eig=TRUE)   #woher kommt jsd?
        
      }else if (input$dist == 'tsne'){
        
        p <- 30
        d <- try(Rtsne::Rtsne(beta,3,theta=.5,perplexity=p),silent=TRUE)
        while(class(d) == 'try-error'){
          p <- p-1
          d <- try(Rtsne::Rtsne(beta,3,theta=.5,perplexity=p),silent=TRUE)
        }
        if (p < 30) cat(sprintf('Performed t-SNE with perplexity = %s.\n',p))
        d$points <- d$Y
        d$eig <- NULL
        
      }else{
        
        d <- cmdscale(vegan::vegdist(beta,method=input$dist),3,eig=TRUE)
        
      }
      
      eig <- d$eig[1:3]/sum(d$eig)
      colnames(d$points) <- c('Axis1','Axis2','Axis3')
      df <- data.frame(d$points,EST()$df0)
      df$marg <- vis_out$topic_marg
      
      df$colors <- vis_out$colors[as.character(df$sig)]
      
      if (input$dim == '2d'){
        
        p1 <- plot_ly(df,source='ord_click')
        p1 <- add_trace(p1,
                        x=~Axis1,y=~Axis2,size=~marg,
                        type='scatter',mode='markers',sizes=c(5,125),
                        color=I(df$colors),opacity=.5,
                        marker=list(symbol='circle',sizemode='diameter',line=list(width=3,color='#FFFFFF')),
                        text=~paste('<br>Topic:',topic),hoverinfo='text')
        p1 <- layout(p1,
                     showlegend=FALSE,
                     xaxis=list(title=sprintf('Axis 1 [%.02f%%]',eig[1]*100),
                                showgrid=FALSE),
                     yaxis=list(title=sprintf('Axis 2 [%.02f%%]',eig[2]*100),
                                showgrid=FALSE),
                     paper_bgcolor='rgb(243, 243, 243)',
                     plot_bgcolor='rgb(243, 243, 243)')
        p1 <- add_annotations(p1,x=df$Axis1,y=df$Axis2,text=df$topic,showarrow=FALSE,
                              font=list(size=10))
        
        h <- event_data('plotly_hover',source='est_hover')
        
        if (!is.null(h)){
          k <- EST()$k_levels[h[['x']]]
          
          if (length(k) > 0){
            df_update <- df[df$topic == k,]
            
            
            if (df_update$sig == '1') df_update$sig <- '2' else if(df_update$sig== '-1') df_update$sig <- '-2' else df_update$sig<- '00'
            df_update$colors <- vis_out$colors[df_update$sig]
            
            p1 <- add_markers(p1,
                              x=df_update$Axis1,y=df_update$Axis2,opacity=.8,color=I(df_update$color),
                              marker=list(size=150,symbol='circle',sizemode='diameter',line=list(width=3,color='#000000')))
          }
          
          p1
          
        }
        
      }
      
      if (input$dim == '3d'){
        
        p1 <- plot_ly(df,source='ord_click',
                      x=~Axis1,y=~Axis2,z=~Axis3,size=~marg,
                      type='scatter3d',mode='markers',sizes=c(5,125),
                      color=I(df$colors),opacity=.5,
                      marker=list(symbol='circle',sizemode='diameter'),
                      text=~paste('<br>Topic:',topic),hoverinfo='text')
        
        p1 <- layout(p1,
                     showlegend=FALSE,
                     scene=list(
                       xaxis=list(title=sprintf('Axis 1 [%.02f%%]',eig[1]*100),
                                  showgrid=FALSE),
                       yaxis=list(title=sprintf('Axis 2 [%.02f%%]',eig[2]*100),
                                  showgrid=FALSE),
                       zaxis=list(title=sprintf('Axis 3 [%.02f%%]',eig[3]*100),
                                  showgrid=FALSE)),
                     paper_bgcolor='rgb(243, 243, 243)',
                     plot_bgcolor='rgb(243, 243, 243)')
        
      }
      
      p1
    }
  })
  
  show_topic <- reactiveValues(k=0)
  
  #maybe not working...
  observeEvent(event_data('plotly_click',source='ord_click'),{
    
    s <- event_data('plotly_click',source='ord_click')
    
    if (is.null(s)){
      
      show_topic$k <- 0
      updateNumericInput(session,'k_in',value=0)
      
    }else{
      
      t_idx <- s$pointNumber + 1
      updateNumericInput(session,'k_in',value=t_idx)
      show_topic$k <- t_idx
      
    }
    
  })
  
  observeEvent(input$k_in,{
    
    show_topic$k <- input$k_in
    
  })
  
  observeEvent(input$reset,{
    
    show_topic$k <- 0
    updateNumericInput(session,'k_in',value=0)
    
  })
  
  output$bar <- renderPlot({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      if (show_topic$k != 0){
        
        p_bar <- ggplot(data=REL()) +
          geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=.6) +
          geom_bar(aes_(~Term,~Freq),stat='identity',fill='darkred',color='white')
        
      }else{
        
        p_bar <- ggplot(data=REL()) +
          geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=1)
        
      }
      
      p_bar +
        coord_flip() +
        labs(x='',y='Frequency',fill='') +
        theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
              legend.position='bottom') +
        viridis::scale_fill_viridis(discrete=TRUE,drop=FALSE) +
        guides(fill=guide_legend(nrow=2))
    }
  })
  
  output$corr <- networkD3::renderForceNetwork({
    vis_out <- vals$datasets[[currentSet()]]$vis_out
    if(!is.null(vis_out)){
      K <- nrow(vis_out$corr$posadj)
      
      suppressWarnings({suppressMessages({
        g <- igraph::graph.adjacency(vis_out$corr$posadj,mode='undirected',
                                     weighted=TRUE,diag=FALSE)
        
        wc <- igraph::cluster_walktrap(g)
        members <- igraph::membership(wc)
        
        g_d3 <- networkD3::igraph_to_networkD3(g,group=members)
        
        g_d3$links$edge_width <- 10*(.1+sapply(seq_len(nrow(g_d3$links)),function(r) vis_out$corr$poscor[g_d3$links$source[r]+1,g_d3$links$target[r]+1]))
        g_d3$nodes$color <- 1
        g_d3$nodes$node_size <- 25*norm10(c(0,vis_out$topic_marg))[-1]
        g_d3$nodes$name <- paste0('T',g_d3$nodes$name)
        
        networkD3::forceNetwork(Links=g_d3$links,Nodes=g_d3$nodes,
                                Source='source',Target='target',
                                charge=-25,
                                opacity=.7,
                                fontSize=12,
                                zoom=TRUE,
                                bounded=TRUE,
                                NodeID='name',
                                fontFamily='sans-serif',
                                opacityNoHover=.7,
                                Group='color',
                                Value='edge_width',
                                Nodesize='node_size',
                                linkColour='#000000',
                                linkWidth=networkD3::JS('function(d) {return d.value;}'),
                                radiusCalculation=networkD3::JS('d.nodesize'),
                                colourScale=networkD3::JS("color=d3.scaleLinear()\n.domain([1,1])\n.range(['red','red']);"))
        
      })})
    }
  })
  
}
