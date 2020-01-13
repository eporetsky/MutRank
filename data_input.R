
dataInputUI <- function(id) {
  ns <- NS(id)
  useShinyjs()
  tabPanel("Data Input",
    tags$head(tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
    sidebarPanel(width = 12,
      selectInput(ns("expression"), "Choose an expression dataset:", 
                  choices = c("mini.csv", list.files("data/", pattern='*.csv'))
      ),
      fileInput(ns("submit_expression"), "Choose CSV File", multiple = FALSE,
                accept = c("text/csv","text/comma-separated-values,text/plain",".csv")
      ), hr(),
      actionButton(ns('reset_submit'), 'Reset'),
      tags$hr()
    ),
    mainPanel(
      tableOutput({ns("df_input")})
    )
  )
}

dataInput <- function(input, output, session) {
  ns <- session$ns
  # Different input functions, should try to de-clatter it in the near future              ###
  
  default <- reactiveVal()
  observeEvent(input$expression,{default(input$expression)})
  submitted <- reactiveVal()
  observeEvent(input$submit_expression,{submitted(input$submit_expression)})
  
  observeEvent(input$reset_submit, {
    submitted(NULL)
    default(NULL)
  })
  
  expression <- reactive({
    if(is.null(submitted())){
      data_loader(read.csv(paste("data/", default(), sep=""), header=T))}
    else{
      data_loader(read.csv(submitted()$datapath, header = T))
    }
  })
  
  metabolites <- reactive({
    NA
    #  data_loader(input$metabolites)
  })

  output$df_input <- renderTable({
    expression()[1:5,1:5]
  })
  
  return(list(expression,metabolites))
}

data_loader <- function(temp_df){
  # User has not uploaded a file yet
  duplicated <- duplicated(temp_df[,1])
  deduped <- temp_df[!duplicated,-1]
  row.names(deduped) <- temp_df[!duplicated,1]
  new_data <- as.matrix(deduped)
  return(new_data)
}