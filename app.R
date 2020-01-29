######################################################
## mutRank Version 0.8
## Written by Elly Poretsky
## Alisa Huffaker Lab
## University of California, San Diego
##
## mutRank is Shiny web application
## Other dependencies are cited when possible,
## or mentioned below
######################################################

## Load necessary packages
if (!require("phytools")) {
  install.packages("shiny", dependencies = TRUE)
  library(shiny)
}

library(shinyjs)
library(shiny)
library(shinycssloaders)
library(shinythemes)
library(igraph)
library(visNetwork)
library(ggplot2)
library(RColorBrewer) # for network node colors
library(reshape2) # used for the heatmap
library(ontologyIndex)
library(hypergea)
library(data.table)

source("data_input.R")
source("module_mutualrank.R")
source("module_network.R")
source("module_heatmap.R")
source("module_enrichment.R")

ui <- fluidPage(
  navbarPage("mutRank v0.8",
    theme = shinytheme("flatly"),
    dataInputUI("dataInputNS"),
    mutualRankUI("mutualRankNS"),
    heatmapUI("heatmapNS"),
    MRnetworkUI("MRnetworkNS"),
    enrichmentUI("enrichmentNS")
    
  )
)

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  data <- callModule(dataInput,"dataInputNS")
  coexpression <- callModule(mutualRank,"mutualRankNS", data$expression,data$annotations)
  callModule(heatmap,"heatmapNS",coexpression)
  callModule(MRnetwork,"MRnetworkNS", coexpression, data$annotations, data$foldchange, data$association)
  callModule(enrichment,"enrichmentNS", coexpression, data$GO_db, data$GO_genes)
}

shinyApp(ui, server)