######################################################
## MutRank Version 0.9
## Written by Elly Poretsky
## Alisa Huffaker Lab
## University of California, San Diego
##
## mutRank is Shiny web application
## Other dependencies are cited when possible,
## or mentioned below
######################################################

#install.packages("rmarkdown")
#install.packages("pdflatex")
#install.packages("rmarkdown")
#rmarkdown::render("README.md", "pdf_document")

# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
# Now load or install&load all

package.check <- lapply(
  c("shiny","shinythemes", "igraph","visNetwork", "ggplot2",
    "data.table", "RColorBrewer","reshape2", "ontologyIndex","hypergea"),
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }})


source("global.R")
source("data_input.R")
source("module_mutualrank.R")
source("module_network.R")
source("module_heatmap.R")
source("module_enrichment.R")
source("mutrank_functions.R")

ui <- fluidPage(
  
  navbarPage("MutRank v0.10",
    theme = shinytheme("flatly"),
    dataInputUI("dataInputNS"),
    mutualRankUI("mutualRankNS"),
    heatmapUI("heatmapNS"),
    MRnetworkUI("MRnetworkNS"),
    enrichmentUI("enrichmentNS")
  ),
  tags$head(tags$style('body {color:black; font-size: 175%;}')),
  tags$head(tags$style(HTML('.shiny-bound-input {color:black; font-size: 100%;}')))
)

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  #reactive({print(session$clientData$url_pathname())})
  data <- callModule(dataInput,"dataInputNS")
  coexpression <- callModule(mutualRank,"mutualRankNS", 
                             data$expression,
                             data$annotations,
                             data$symbols,
                             data$categories,
                             data$foldchange,
                             data$GO_genes,
                             data$domains)
  callModule(heatmap,"heatmapNS",
             coexpression, 
             data$symbols)
  callModule(MRnetwork,"MRnetworkNS", 
             coexpression, 
             data$annotations, 
             data$symbols, 
             data$foldchange, 
             data$categories,
             data$GO_genes,
             data$domains)
  callModule(enrichment,"enrichmentNS", 
             coexpression, 
             data$GO_db, 
             data$GO_genes)
}

shinyApp(ui, server)