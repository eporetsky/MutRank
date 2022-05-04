enrichmentUI <- function(id) {
  ns <- NS(id)
  tabPanel("Enrichment",
    tags$head(tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
    sidebarPanel(width = 4,
      numericInput(ns("mr_cutoff"), "Include genes with MR below:", 100),
      textOutput({ns("included_elements")}),
      actionButton(ns("update_enrichment"), "Update results"),
      p(), hr(), p(),
      checkboxInput(ns("include_pvals"), "Include non-adjusted p values?", value = T, width = NULL),
      checkboxInput(ns("include_values"), "Include values used for enrichment calculations?", value = F, width = NULL),
      checkboxInput(ns("include_genes"), "Include genes in each GO term?", value = F, width = NULL),
      selectInput(ns("fdr_method"), "Choose FDR method:", choices =c("holm","hochberg","hommel","bonferroni","BH","BY"), selected="BH"),
      p(), hr(), p(),
      downloadButton(ns('download_enrichment'),"Download table")
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
    mr_cutoff <- isolate(input$mr_cutoff)
    coexpression <- coexpression[(coexpression[,1]<=mr_cutoff),]
    return(rownames(coexpression))
  })
  
  # Returns how many genes are included in the enrichment calculations
  output$included_elements <- renderText(paste("Number of genes included: ", toString(length(gene_list()))))
  
  GO_enrichment_table <- reactive({
    GO_enrichment(gene_list(),GO_db(),GO_genes(),input$fdr_method,
                  input$include_pvals, input$include_values, input$include_genes)
  })
  
  # Returns the enrichment table
  output$enrichment_df <- renderTable({
    GO_enrichment_table()},
    class = 'cell-border stripe nowrap compact dt-center',
    bordered = TRUE,
    rownames = TRUE
  )
  
  output$download_enrichment <- downloadHandler(
    filename = function(){"mutRank_GO_enrichment.tsv"}, 
    content = function(fname){write.table(GO_enrichment_table(), fname, sep="\t", quote=F, col.names=NA)}
  )
}

GO_enrichment <- function(gene_set, GO_db, GO_genes, fdr_method,
                          include_pvals,include_values,include_genes){
  ## Gene Ontology OBO files downloaded from GO Consortium
  ## http://geneontology.org/page/download-ontology
  
  # Prepare all the data for a hypergeomteric test:
  # 1. num_genes -        Total number of genes in the GO database
  # 2. num_genes_go -     Total number of genes with GO term
  # 3. num_set -          Number of genes in the selected set
  # 4. num_genes_go_set - Number of genes in the selected set annotated with GO term
  
  terms <- row.names(GO_db)                                       # The list of GO terms in the OBO file
  filtered <- GO_genes[GO_genes[,2] %in% terms,]                  # GO Slim has less terms so I filter based on it
  filtered[,2] <- factor(filtered[,2])                            # Removes R factors
  count_go <- as.data.frame(table(filtered[,2]))                  # Returns frequency table for each GO term
  num_genes <- length(unique(GO_genes[,1]))                       ### Returns 1. Total number of genes in the GO database
  num_set <- length(gene_set)                                     ### Returns 3. Number of genes in the selected set
  filtered_set <- GO_genes[GO_genes[,1] %in% gene_set,]           # Filter the GO_genes included in genes set 
  filtered_set_slim <- filtered_set[filtered_set[,2] %in% terms,] # GO Slim has less terms so I filter based on it
  filtered_set_slim[,2] <- factor(filtered_set_slim[,2])          # Removes factor from filtered_set_slim
  count_set <- as.data.frame(table(filtered_set_slim[,2]))        # Create a frequency table
  go_list <- unique(count_set$Var1)                               # Get a list of unique GO terms for the analysis
  
  # These values are not used for calculations, but to be shown in the final Shiny table
  N.list <-  c() # 1. Total number of genes in the GO database
  M.list <-  c() # 2. Total number of genes with GO term
  n.list <-  c() # 3. Number of genes in the selected set
  m.list <-  c() # 4. Number of genes in the selected set annotated with GO term
  p.list <-  c() # List of p.values
  t.list <-  c() # List of GO term description
  fc.list <- c() # List of fold-enrichment
  
  for(GO in go_list){
    num_genes_go <- count_go[count_go$Var1==GO,]$Freq       ### 3. Number of genes in the selected set
    num_genes_go_set <- count_set[count_set$Var1==GO,]$Freq ### 4. Number of genes in the selected set annotated with GO term 
    N.list <- append(N.list,num_genes)                      # Add N to the list for each GO term
    M.list <- append(M.list,num_set)                        # Add M to the list for each GO term
    n.list <- append(n.list,num_genes_go)                   # Add n to the list for each GO term
    m.list <- append(m.list,num_genes_go_set)               # Add m to the list for each GO term
    
    #dmn <- list(set=c('y', 'n'), total=c('y', 'n'))         # A list to specifiy the matrix dimnames needed for the hypergeom.test
    #CT <- matrix(c(num_genes_go_set,num_set,num_genes_go,   # The matrix as specified for the hypergeom.test
    #               num_genes), nrow=2, dimnames=dmn)
    
    
    p.val <- dhyper(num_genes_go_set,num_set,
                    num_genes,num_genes_go)                 # Get the p.value after performing the hypergeom.test
    p.list <- append(p.list,p.val)                          # Add the p.value to the p.list
    t.list <- append(t.list, as.vector(GO_db[GO,]))         # Add the GO term description to the list 
    fc_enrichment <- (num_genes_go_set/num_set) /           # Calculate the GO term fold-change for the gene-set
                     (num_genes_go/num_genes)               
    fc.list <- append(fc.list,fc_enrichment)                # Add the GO term fold-change enrichment to the list
  }
  
  fdr_name <- paste("p.adj",fdr_method,sep="_")             # The name of the p.value adjustment method
  p.list.adj <- p.adjust(p.list,method=fdr_method)          # Adjust the p.value using the specified method
  p.df <- data.frame(GO=go_list,                            # Create the data frame that will be returned to the Shiny app
                     N=N.list,
                     M=M.list,
                     n=n.list,
                     m=m.list,
                     p.val=p.list,
                     p.adj=p.list.adj, 
                     fc=fc.list, 
                     description=t.list)
  colnames(p.df)[7] <- fdr_name                             # Change the column name to speficiy which FDR was used
  
  # Delete the p.value column if user specified not to keep it
  if(!include_pvals){p.df["p.val"]<-NULL}
  # Delete these columns if user specified not to keep them
  if(!include_values){
    p.df["N"]<-NULL
    p.df["M"]<-NULL
    p.df["n"]<-NULL
    p.df["m"]<-NULL}
  # Add a column contianing the set of genes included in each GO term
  if(include_genes){
    coln <- colnames(filtered_set_slim)
    genes=aggregate(filtered_set_slim,
                    by=list(filtered_set_slim[,coln[2]]),
                    FUN="paste",collapse=",")[c("Group.1",coln[1])]
    p.df <- merge(p.df,genes,by.x="GO",by.y="Group.1")
  }
  return(p.df[order(p.df[fdr_name]),])
}
