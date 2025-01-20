output$selected_number <- renderText({
  paste("You selected:", input$number)
})