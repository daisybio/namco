library(biomaRt)
#library(data.table)
#library(dplyr)
library(DT)
#library(factoextra)
#library(future)
#library(ggplot2)
#library(ggpubr)
#library(ggrepel)
library(plotly)
#library(plyr)
#library(pryr)
library(RColorBrewer)
library(shiny)
#library(splitstackshape)
#library(tidyr)
library(textshape)
library(visNetwork)


server <- function(input,output,session){
  options(shiny.maxRequestSize=1000*1024^2,stringsAsFactors=F)  # upload up to 1GB of data
  vals = reactiveValues(datasets=list()) # reactiveValues is a container for variables that might change during runtime and that influence one or more outputs, e.g. the currently selected dataset
  currentSet = NULL # a pointer to the currently selected dataset
  minLibSize = 1000 # a lower boundary for the expected library size; if in a bait less than this number of proteins are quantified, a warning message will be shown
  lastplot = NULL
  
  ################################################################################################################################################
  #
  # file upload
  #
  ################################################################################################################################################
  # Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
  uploadModal <- function(failed=F) {
    modalDialog(
      p("Please provide an OTU table."),
      fileInput("analysisFile","Select OTU table",width="100%"),
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
    if(!is.null(input$analysisFile)&!input$dataName%in%names(vals$datasets)){
      tryCatch({
        dat <- read.csv(input$analysisFile$datapath,header=T,sep="\t")
        colnames(dat)[1] = "taxa"
        vals$datasets[[input$dataName]] <- dat
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
  output$otuTable <- renderDataTable({
    if(!is.null(currentSet())) datatable(vals$datasets[[currentSet()]],filter='bottom',options=list(searching=T,pageLength=20,dom="Blfrtip"),editable=T)
    else datatable(data.frame(),options=list(dom="t"))
  })
  
  output$distribution <- renderPlotly({
    if(!is.null(currentSet())){
      data = vals$datasets[[currentSet()]]
      #data=dat
      rownames(data)=data$taxa
      data = data[,-1]
      
      data = apply(data,2,function(x) x/sum(x))
      
      taxa <- ifelse(rowSums(data)/ncol(data)<0.05,"Other",rownames(data))
      other <- data[which(taxa=="Other"),]
      other <- colSums(other)
      
      data <- data[-which(taxa=="Other"),]
      data <- rbind(data,Other=as.data.frame(t(other)))
      data = 100*data
      
      data = data[,1:10]
      
      p <- plot_ly(x=colnames(data),y=t(data[1,]),name=rownames(data)[1],type="bar") #%>%
      for(i in 2:nrow(data)) p <- p %>% add_trace(y=t(data[i,]),name=rownames(data)[i])
      p %>% layout(yaxis=list(title='Count'),barmode='stack')
    } else plotly_empty()
  })
  
  
  ################################################################################################################################################
  #
  # Network analysis
  #
  ################################################################################################################################################
}

