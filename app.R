######################################################
## mutRank Version 0.7
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
library(shinythemes)
library(igraph)
library(visNetwork)
library(ggplot2)
library(RColorBrewer) # for network node colors
library(reshape2) # used for the heatmap

source("data_input.R")
source("module_mutualrank.R")
source("module_network.R")
source("module_heatmap.R")

ui <- fluidPage(
  navbarPage("mutRank 0.7",
  theme = shinytheme("flatly"),
  dataInputUI("dataInputNS"),
  mutualRankUI("mutualRankNS"),
  heatmapUI("heatmapNS"),
  MRnetworkUI("MRnetworkNS")
  )
)

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  df_loaded <- callModule(dataInput,"dataInputNS")
  df_coexpression <- callModule(mutualRank,"mutualRankNS", df_loaded[[1]])
  callModule(heatmap,"heatmapNS",df_coexpression)
  callModule(MRnetwork,"MRnetworkNS", df_coexpression)
}

shinyApp(ui, server)