enrichmentUI <- function(id) {
  ns <- NS(id)
  useShinyjs()
  tabPanel("Enrichment",
    tags$head(tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
    sidebarPanel(width = 4,
      selectInput(ns("threshold"), "MR Threshold:", selected=12, choices = c(5,10,12,24,50,58,100,116,231,Inf)),
      actionButton(ns("update_enrichment"), "Update Network"),
      textOutput({ns("included_elements")}),
    ),
    mainPanel(
      tableOutput(ns("enrichment_df"))
    )
  )
}

enrichment <- function(input, output, session, coexpression, GO_db, GO_genes) {
  ns <- session$ns
  
  # Returns a vector of gene names that are below the MR threshold
  gene_list <- reactive({
    input$update_enrichment
    coexpression <- coexpression()
    coexpression <- coexpression[coexpression[,1]<=as.numeric(input$threshold)]
    return(rownames(coexpression))
  })
  
  # Returns how many genes are included in the enrichment calculations
  output$included_elements <- renderText(paste("Number of elements included: ", toString(length(gene_list()))))
  
  # Returns the enrichment table
  output$enrichment_df <- renderTable({
    GO_enrichment(gene_list(),GO_db(),GO_genes())},
    class = 'cell-border stripe nowrap compact dt-center',
    bordered = TRUE,
    rownames = TRUE
  )
}

GO_enrichment <- function(gene_set, GO_db, GO_genes){
  ### http://geneontology.org/page/download-ontology
  # Prepare all the data for a hypergeomteric test
  # num_genes - total number of genes
  # num_genes_go - total number of genes with GO term
  # num_genes_go_set - number of genes in the set with GO term
  # num_set - number of genes in the set
  terms <- row.names(GO_db)
  filtered <- GO_genes[GO_genes[,2] %in% terms,] # filters GOs when Slim is selected
  filtered[,2] <- factor(filtered[,2]) # removes factors
  count_go <- as.data.frame(table(filtered[,2])) # returns frequency table for each GO term
  num_genes <- length(unique(GO_genes[,1]))
  num_set <- length(gene_set)
  filtered_set <- GO_genes[GO_genes[,1] %in% gene_set,]
  filtered_set_slim <- filtered_set[filtered_set[,2] %in% terms,]
  filtered_set_slim[,2] <- factor(filtered_set_slim[,2]) # removes factors
  count_set <- as.data.frame(table(filtered_set_slim[,2])) # returns frequency table
  go_list <- unique(count_set$Var1)
  
  p.list <- c()  # list of p.values
  t.list <- c()  # list of GO term description
  fc.list <- c() # list of fold-enrichment
  for(GO in go_list){
    num_genes_go_set <- count_set[count_set$Var1==GO,]$Freq
    num_genes_go <- count_go[count_go$Var1==GO,]$Freq
    dmn <- list(set=c('y', 'n'), total=c('y', 'n'))
    CT <- matrix(c(num_genes_go_set,num_set,num_genes_go,num_genes), nrow=2, dimnames=dmn)
    p.val <- hypergeom.test(CT)$p.value
    p.list <- append(p.list,p.val)
    t.list <- append(t.list, as.vector(GO_db[GO,]))
    fc_enrichment <- (num_genes_go_set/num_set)/(num_genes_go/num_genes)
    fc.list <- append(fc.list,fc_enrichment)
  }
  # Adjust p. value using Benjamini & Yekutieli method
  p.list.adj <- p.adjust(p.list,method="BY")
  p.df <- data.frame(p.val=p.list.adj, fc=fc.list, GO=go_list,description=t.list)
  return(p.df[order(p.df$p.val),])
}