heatmapUI <- function(id) {
  ns <- NS(id)
  ?numericInput
  tabPanel("Heatmap",
    sidebarPanel(
      width = 4,
      numericInput(ns("n_rows"), "Number of rows to include:",25),
      checkboxInput(ns("order_coexpression"), "Order by MR?", value = T, width = NULL),
      numericInput(ns("text_threshold"), "Text MR threshold:",10),
      numericInput(ns("text_size"), "Text font size:",10,step=1),
      p(),
      downloadButton(ns('heatmap_download'),"Download Heatmap"),
     ),
    mainPanel(
      plotOutput(ns("heatmap_plot"),width="1000px", height="600px")
    )
  )
}



heatmap <- function(input, output, session, coexpression, symbols) {
  ns <- session$ns
  
  row_names <- reactive({symbol_converter(symbols(),rownames(final_cormat()))})
  
  final_cormat <- reactive({
    coexpression <- coexpression()
    if(input$order_coexpression){coexpression <- order_coexpression_table(coexpression)}
    if(input$n_rows<length(rownames(coexpression))){coexpression <- coexpression[1:input$n_rows,1:input$n_rows]}
    return(coexpression)
  })
  
  melted_cormat <- reactive({
    adj_matrix <- as.matrix(final_cormat())
    rownames(adj_matrix)<-row_names()
    colnames(adj_matrix)<-row_names()
    adj_matrix <- round(adj_matrix)
    diag(adj_matrix) <- rep(-1,length(diag(adj_matrix)))
    get_tri <- get_upper_tri(adj_matrix)
    rotate <- function(x) t(apply(x, 2, rev)) # rotates matrix
    get_tri <- rotate(rotate(get_tri)) # rotate twice to make first gene top-left
    return(reshape2::melt(as.matrix(get_tri), id=c("rows", "cols"), na.rm = TRUE))
  })
  
  output$heatmap_plot <- renderPlot({
    draw_heatmap(melted_cormat(), input$text_threshold,input$text_size)
  })
  
  output$heatmap_download <- downloadHandler(
    filename ="coexpression_heatmap.png",
    content = function(file) {ggsave(file, plot = draw_heatmap(melted_cormat()), device = "png")}
  )
}

draw_heatmap <- function(melted_cormat, text_threshold, text_size){
  num <- length(rownames(melted_cormat))
  # Change all the large values to 100 so they appear white instead of grey
  if(sum(melted_cormat$value>100)>0){melted_cormat[melted_cormat$value>100,]$value=100}
  # We can start making the heatmap using ggplot
  ggplot(data = melted_cormat, aes(Var1,Var2,fill=value))+
    geom_tile(color = "black",size=1,aes(fill = value))+
    geom_text(data=subset(melted_cormat,(melted_cormat$value>0 & melted_cormat$value<=text_threshold)),
              aes(label=value),size=log(text_size*10000/num))+
    scale_fill_gradient2(low = "red", high = "white",# mid = "white", 
                         limit = c(0,100), space = "Lab", midpoint = 95, 
                         breaks=c(0,20,40,60,80,100),
                         labels=c("0","20","40","60","80",">100"),
                         name="Mutual rank scores", guide="colourbar"
                         ) +
    theme(# Hide panel borders and remove grid lines
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Remove axis titles
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      # Change axis line
      axis.line = element_line(colour = "white"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, hjust = 0, colour="black",face="bold"),
      axis.text.y = element_text(angle = 0, vjust = 0.5, size = 12, hjust = 0, colour="black",face="bold"),
      axis.ticks = element_blank(),
      # legend stuff
      legend.title = element_text(size=13,face="bold"),
      legend.text = element_text(size=13,face="bold"),
      legend.justification = c(1, 0),
      legend.position = c(1, 0.7),
      legend.direction = "horizontal",
      ) +
    
    guides(fill=guide_colorbar(barwidth=15,barheight=2,title.position="top",title.hjust=0.5,
                               frame.colour = "black", frame.linewidth=2, ticks = F))+
    coord_fixed() +
    scale_x_discrete(limits = rev(levels(melted_cormat$Var1))) # reverses the x-axis
}

# Get upper triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Not used, but to get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}