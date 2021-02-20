# Return a dialog window for dataset selection and upload. If 'failed' is TRUE, then display a message that the previous value was invalid.
uploadFastqModal <- function(failed=F,error_message=NULL) {
  modalDialog(
    h4("Please provide the directory with fastq files and corresponding meta file. For detailed information on how the files have to look, check out the Info & Settings tab on the left!"),
    hr(),
    fluidRow(
      column(6,fileInput("fastqFiles","Select all fastq files", multiple = T, accept = ".fastq")),
      column(6,fileInput("metaFile","Select Metadata File"))
    ),
    
    textInput("dataName","Enter a project name:",placeholder="New_Project",value="New_Project"),
    if(failed) {
      #div(tags$b("The file you specified could not be loaded. Please check the Info tab and to confirm your data is in the correct format!",style="color: red;"))
      div(tags$b(error_message,style="color:red;"))
    },
    footer = tagList(
      modalButton("Cancel"),
      actionButton("upload_fastq_ok","OK",style="background-color:blue; color:white")
    )
  )
}

# launch upload dialog
observeEvent(input$upload_fastq, {
  showModal(uploadFastqModal())
})

observeEvent(input$upload_fastq_ok, {
  tryCatch({
    rm_spikes_command = "python3.8 ../src/rm_spikes.py -h"
    out <- system(rm_spikes_command, wait = T)
  },error=function(e){
    print(e)
    showModal(uploadFastqModal(failed=T,error_message = e))
  })
})