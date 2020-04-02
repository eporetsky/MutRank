enrichmentUI <- function(id) {
  ns <- NS(id)
  useShinyjs()
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
  
  N.list <- c()
  M.list <- c()
  n.list <- c()
  m.list <- c()
  
  p.list <- c()  # list of p.values
  t.list <- c()  # list of GO term description
  fc.list <- c() # list of fold-enrichment
  for(GO in go_list){
    num_genes_go_set <- count_set[count_set$Var1==GO,]$Freq
    num_genes_go <- count_go[count_go$Var1==GO,]$Freq
    dmn <- list(set=c('y', 'n'), total=c('y', 'n'))
    
    N.list <- append(N.list,num_genes)
    M.list <- append(M.list,num_set)
    n.list <- append(n.list,num_genes_go)
    m.list <- append(m.list,num_genes_go_set)
    
    CT <- matrix(c(num_genes_go_set,num_set,num_genes_go,num_genes), nrow=2, dimnames=dmn)
    p.val <- hypergeom.test(CT)$p.value
    p.list <- append(p.list,p.val)
    t.list <- append(t.list, as.vector(GO_db[GO,]))
    fc_enrichment <- (num_genes_go_set/num_set)/(num_genes_go/num_genes)
    fc.list <- append(fc.list,fc_enrichment)
  }
  # Adjust p. value using Benjamini & Yekutieli method
  fdr_name <- paste("p.adj",fdr_method,sep="_")
  p.list.adj <- p.adjust(p.list,method=fdr_method)
  p.df <- data.frame(
                     GO=go_list,
                     N=N.list,
                     M=M.list,
                     n=n.list,
                     m=m.list,
                     p.val=p.list,
                     p.adj=p.list.adj, 
                     fc=fc.list, 
                     description=t.list
  )
  colnames(p.df)[7]<-fdr_name
  
  if(!include_pvals){p.df["p.val"]<-NULL}
  if(!include_values){
    p.df["N"]<-NULL
    p.df["M"]<-NULL
    p.df["n"]<-NULL
    p.df["m"]<-NULL}
  if(include_genes){
    coln <- colnames(filtered_set_slim)
    genes=aggregate(filtered_set_slim,by=list(filtered_set_slim[,coln[2]]),FUN="paste",collapse=",")[c("Group.1",coln[1])]
    p.df <- merge(p.df,genes,by.x="GO",by.y="Group.1")
  }
  return(p.df[order(p.df[fdr_name]),])
}

