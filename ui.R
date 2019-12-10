library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(plotly)
library(visNetwork)


ui <- dashboardPage(
  dashboardHeader(title="Microbiome Explorer"),
  dashboardSidebar(
    sidebarMenu(id="sidebar",
                br(),
                fluidRow(
                  column(12,align="center",
                         actionButton("upload","Load new dataset")
                  )
                ),
                dataTableOutput("datasets"),
                br(),br(),
                menuItem("Quality Control",tabName="QC",icon=icon("search")),
                menuItem("Network Analysis",tabName="Network",icon=icon("project-diagram"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "QC",
        h4("Quality Control Report"),
        fluidRow(  
          tabBox(id="outputPlots",width=12,
            tabPanel("Sample Information",
              tags$hr(),
              p("OTU table"),
              fluidRow(
                column(1),
                column(7,dataTableOutput("otuTable"))
              )),
              tabPanel("Taxa Distribution",
                tags$hr(),
                fluidRow(
                  column(1),
                  column(7,plotlyOutput("distribution")),
                  column(1),
                  column(2)
              ))
          )
        )
      ),
      tabItem(tabName = "Network",
              h4("Network Analysis")
      )
    )
  )
)