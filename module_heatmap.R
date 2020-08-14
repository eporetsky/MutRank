heatmapUI <- function(id) {
  ns <- NS(id)
  tabPanel("Heatmap",
    sidebarPanel(
      width = 4,
      numericInput(ns("n_rows"), "Number of rows to include:",25),
      numericInput(ns("text_threshold"), "Text MR threshold:",10),
      numericInput(ns("text_size"), "Text font size:",10,step=1),
      checkboxInput(ns("convert_symbols"), "Convert gene names to symbols", value = T, width = NULL),
      p(),
      downloadButton(ns('heatmap_download'),"Download Heatmap"),
     ),
    mainPanel(
      plotOutput(ns("heatmap_plot"),width="1000px", height="600px")
    )
  )
}

heatmap <- function(input, output, session, 
                    coexpression, symbols) {
  ns <- session$ns
  
  # Convert the gene names in the Mutual Rank table into symbols
  row_names <- reactive({ifelse(input$convert_symbols, 
                                return(symbol_converter(symbols(),rownames(final_cormat()))), 
                                return(rownames(final_cormat())))
  })
  
  # Shiny reactive function to reduce the amount of rows in the heatmap based on user specification
  final_cormat <- reactive({
    coexpression <- coexpression()
    if(input$n_rows<length(rownames(coexpression))){
      coexpression <- coexpression[1:input$n_rows,1:input$n_rows]}
    return(coexpression)
  })
  
  # Main reactive function to convert the Mutual Rank data frame into an adjescency function and melt for ggplot
  melted_cormat <- reactive({
    adj_matrix <- as.matrix(final_cormat())              # Convert the Mutual Rank data.frame into a matrix
    rownames(adj_matrix)<-row_names()                    # Use the row_names() reactive function to set matrix rownames
    colnames(adj_matrix)<-row_names()                    # Use the row_names() reactive function to set matrix colnames
    diag(adj_matrix) <- rep(-1,length(diag(adj_matrix))) # Change the diagonal of the matrix to -1, helpful later
    get_tri <- get_upper_tri(adj_matrix)                 # Convert the lower triangle into NAs
    rotate <- function(x) t(apply(x, 2, rev))            # The function that rotates the matrix
    get_tri <- rotate(rotate(get_tri))                   # Rotate the matrix twice to make first gene top-left
    # Melt the matrix into Var1, Var2 and value table to be used in ggplot
    return(reshape2::melt(as.matrix(get_tri), id=c("rows", "cols"), na.rm = TRUE))
  })
  
  # Call the ggplot function that creates that heatmap
  output$heatmap_plot <- renderPlot({
    draw_heatmap(melted_cormat(), input$text_threshold,input$text_size)
  })
  
  # Shiny handle for downloading the heatmap as a png
  output$heatmap_download <- downloadHandler(
    filename ="coexpression_heatmap.png",
    content = function(file) {ggsave(file, plot = draw_heatmap(melted_cormat()), device = "png")}
  )
}

# Get upper triangle of the correlation matrix (convert the lower triangle to NAs)
# Based on https://stackoverflow.com/questions/59917970/how-do-i-create-a-heat-map-in-r
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
