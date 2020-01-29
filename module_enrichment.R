#install.packages("ontologyIndex")
#install.packages("ontologySimilarity")
##############

enrichmentUI <- function(id) {
  ns <- NS(id)
  useShinyjs()
  tabPanel("Enrichment",
    tags$head(tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
    sidebarPanel(width = 4,
      selectInput(ns("threshold"), "MR Threshold:", selected=12, choices = c(5,10,12,24,50,58,100,116,231)),
      actionButton(ns("update_enrichment"), "Update Network"),
      textOutput({ns("included_elements")}),
    ),
    mainPanel(
      tags$head(tags$style(HTML("#enrichmentNS-enrichment_df table tr td {white-space: nowrap;}"))),
      tableOutput(ns("enrichment_df"))
    )
  )
}
#mr_exponential_decay(231)#100
#mr_exponential_decay(116)#50
#mr_exponential_decay(58) #25
#mr_exponential_decay(24) #10
#mr_exponential_decay(12) #5


enrichment <- function(input, output, session, coexpression, GO_db, GO_genes) {
  ns <- session$ns
  # Different input functions, should try to de-clatter it in the near future ###
  
  gene_list <- reactive({
    input$update_enrichment
    coexpression <- coexpression()
    coexpression <- coexpression[,1][coexpression[,1]<=as.numeric(input$threshold)]
    return(as.data.frame(coexpression))
  })
  
  #output$included_elements <- renderText(paste("Number of elements included: ", toString(length(gene_list))))
  output$included_elements <- renderText(paste("Number of elements included: ", toString(length(rownames(gene_list())))))
  
  output$enrichment_df <- renderTable({
    GO_enrichment(gene_list(),GO_db,GO_genes)},
    class = 'cell-border stripe nowrap compact dt-center',
    bordered = TRUE,
    rownames = TRUE
  )
}

GO_enrichment <- function(coexpression, GO_db, GO_genes){
  GO_genes <- GO_genes()
  terms <- row.names(GO_db())
  filtered <- GO_genes[GO_genes$term_accession %in% terms,] # filters GOs not in slim list
  filtered$term_accession <- factor(filtered$term_accession) # removes factors
  count_go <- as.data.frame(table(filtered$term_accession)) # returns frequency table for each GO term
  num_genes <- length(unique(GO_genes$db_object_id))
  gene_set <- row.names(coexpression)
  num_set <- length(gene_set)
  filtered_set <- GO_genes[GO_genes$db_object_id %in% gene_set,] # filters GOs not in slim list
  filtered_set_slim <- filtered_set[filtered_set$term_accession %in% terms,] # filters GOs not in slim list
  filtered_set_slim$term_accession <- factor(filtered_set_slim$term_accession) # removes factors
  count_set <- as.data.frame(table(filtered_set_slim$term_accession)) # returns frequency table for each GO term
  go_list <- unique(count_set$Var1)
  
  p.list <- c()
  t.list <- c()
  for(GO in go_list){
    num_genes_go_set <- count_set[count_set$Var1==GO,]$Freq
    num_genes_go <- count_go[count_go$Var1==GO,]$Freq
    dmn <- list(set=c('y', 'n'), total=c('y', 'n'))
    CT <- matrix(c(num_genes_go_set,num_set,num_genes_go,num_genes), nrow=2, dimnames=dmn)
    p.val <- hypergeom.test(CT)$p.value
    p.list <- append(p.list,p.val)
    t.list <- append(t.list, as.vector(GO_db()[GO,]))
  }
  p.list.adj <- p.adjust(p.list,method="BY")
  p.df <- data.frame(p.val=p.list.adj, description=t.list)
  #row.names(p.df) <- t.list
  #p.df <- p.df[p.df$p.val<0.05,]
  return(p.df[order(p.df$p.val),])
}

### http://geneontology.org/page/download-ontology
# num_genes - total number of genes
# num_genes_go - total number of genes with GO term
# num_genes_go_set - number of genes in the set with GO term
# num_set - number of genes in the set