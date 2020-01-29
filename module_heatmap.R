heatmapUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Heat Map",
    sidebarPanel(
      width = 2,
      actionButton(ns("update_heatmap"), "Update Heatmap"),
      p(),
      downloadButton(ns('heatmap_download'),"Download Heatmap"),
     ),
    mainPanel(
      plotOutput(ns("heatmap_plot"),width="1000px", height="600px")
    )
  )
}

heatmap <- function(input, output, session, coexpression) {
  ns <- session$ns
  
  melted_cormat <- reactive({
    input$update_heatmap
    adj_matrix <- as.matrix(coexpression())
    adj_matrix <- round(adj_matrix)
    diag(adj_matrix) <- rep(-1,length(diag(adj_matrix)))
    get_tri <- get_lower_tri(adj_matrix)
    # Melt the correlation matrix
    #print(melt(as.matrix(get_tri), id=c("rows", "cols"), na.rm = TRUE))
    melted_cormat <- melt(as.matrix(get_tri), id=c("rows", "cols"), na.rm = TRUE)
    melt(as.matrix(get_tri), id=c("rows", "cols"), na.rm = TRUE)
  })

  output$heatmap_plot <- renderPlot({
    draw_heatmap(melted_cormat())
  })
  output$heatmap_download <- downloadHandler(
    #filename = function() {paste("coexpression_heatmap.", input$image_type, sep = "")},
    filename ="coexpression_heatmap.png",
    content = function(file) {ggsave(file, plot = draw_heatmap(melted_cormat()), device = "png")}
  )
}


# https://github.com/wilkox/ggfittext - cool, maybe future useful
#library(ggrepel) # don't use, but useful if you don't want many labels to overlap

draw_heatmap <- function(melted_cormat){
  num <- length(rownames(melted_cormat))
  # Change all the large values to 100 so they appear white instead of grey
  if(sum(melted_cormat$value>100)>0){melted_cormat[melted_cormat$value>100,]$value=100}
  ggplot(data = melted_cormat, aes(Var1,Var2,fill=value))+
    # create the heatmap
    geom_tile(color = "black",size=1,aes(fill = value))+
    geom_text(data=subset(melted_cormat,(melted_cormat$value>0 & melted_cormat$value<=10)),aes(label=value),size=log(3000/num))+
    
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
  
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}