citation("ggplot2")

######################################################
## MutRank Version 1.0
## Written by Elly Poretsky
## Alisa Huffaker Lab
## University of California, San Diego
##
## mutRank is Shiny web application for Mutual Rank calculations
## More information at: https://github.com/eporetsky/MutRank
######################################################

# Install packages that are not yet installed and loads them
# Based on - https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
package.check <- lapply(
  c("shiny","shinythemes", "igraph","reshape2","visNetwork", "ggplot2",
    "data.table", "RColorBrewer", "ontologyIndex","hypergea"),
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

source("global.R")             # Loads default_files.csv to specifit which files to load when MutRank starts
source("module_datainput.R")   # Shiny module that handles loading user specified data for the analysis
source("module_mutualrank.R")  # Shiny module that handles the caluclation and printing of the Mutual Rank table 
source("module_heatmap.R")     # Shiny module that handles the plotting part of the Mutual Rank table as a heatmap
source("module_network.R")     # Shiny module that handles plotting the Mutual Rank table as a network using JS
source("module_enrichment.R")  # Shiny module that handles calculation of the GO term enrichment based on the Mutual Rank table
source("mutrank_functions.R")  # R script file that contains the functions used by the different modules

# Shiny UI function
ui <- fluidPage(
  navbarPage("MutRank v1.0",
    theme = shinytheme("flatly"),
    dataInputUI("dataInputNS"),
    mutualRankUI("mutualRankNS"),
    heatmapUI("heatmapNS"),
    networkUI("networkNS"),
    enrichmentUI("enrichmentNS")
  ),
  tags$head(tags$style('body {color:black; font-size: 175%;}')),
  tags$head(tags$style(HTML('.shiny-bound-input {color:black; font-size: 100%;}')))
)

# Shiny server function
server <- function(input, output, session){
  # Stops the running Shiny app when the app window is closed
  session$onSessionEnded(stopApp)
  
  # Calls module_datainput.R and returns a list containing all the loaded data
  data <- callModule(dataInput,"dataInputNS")
  
  # Calls module_mutualrank.R and returns the Mutual Rank table in  mutuak_rank_table
  mutuak_rank_table <- callModule(mutualRank,"mutualRankNS", 
                                  data$expression,
                                  data$annotations,
                                  data$symbols,
                                  data$categories,
                                  data$foldchange,
                                  data$GO_genes,
                                  data$domains)
  
  # Calls module_heatmap.R
  callModule(heatmap,"heatmapNS",
             mutuak_rank_table, 
             data$symbols)
  
  # Calls module_network.R
  callModule(network,"networkNS", 
             mutuak_rank_table, 
             data$annotations, 
             data$symbols, 
             data$foldchange, 
             data$categories,
             data$GO_genes,
             data$domains)
  
  # Calls module_enrichment.R
  callModule(enrichment,"enrichmentNS", 
             mutuak_rank_table, 
             data$GO_db, 
             data$GO_genes)
}

# Create and runs the Shiny app
shinyApp(ui, server)