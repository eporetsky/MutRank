mutualRankUI <- function(id){
  ns <- NS(id)
  
  tabPanel("Mutual Rank",
     sidebarPanel(
       width = 4,
       textInput(ns('gene'), "Gene Name","GRMZM2G075459"),
       numericInput(ns("num"), "Number of Genes to correlate:", 5),
       actionButton(ns("update_button"), "Update"),
       selectInput(ns("mrVScor"), "MR or Correlation:", choices = c("MR", "PCC")),
       checkboxInput(ns("firstColumn"), "Only show first column?", value = T, width = NULL),
       checkboxInput(ns("order"), "Order table by MR?", value = T, width = NULL),
       checkboxInput(ns("round"), "Round to nearest integer?", value = T, width = NULL),
       checkboxInput(ns("annotate"), "Annotate Genes", value = T, width = NULL),
       hr()
     ),
     mainPanel(
      tags$head(tags$style(HTML("#mutualRankNS-coexpression_table table tr td {white-space: nowrap;}"))),
      tableOutput(ns("coexpression_table"))
      )
    )
}

mutualRank <- function(input, output, session, data, annotations) {
  coexpression <- reactive({
    input$update_button
    gene_list <- isolate(input$gene)
    num <- isolate(input$num)
    mrVScor <- isolate(input$mrVScor)
    coexpression_table(data(), gene_list, num, mrVScor)
  })
  
  coexpression_df_prefs <- reactive({
    df_output_editor(coexpression(),input$firstColumn,input$mrVScor,input$order,input$round)
  })
  
  output$coexpression_table <- renderTable({
    if(input$annotate){output_df <- df_annotator(coexpression_df_prefs(),annotations())} 
    else{output_df <- coexpression_df_prefs()}},
    class = 'cell-border stripe nowrap compact dt-center',
    bordered = TRUE,
    rownames = TRUE
  )
  return(coexpression)
}

coexpression_table <- function(datm, genes, n, mrVScor){
  #### Working Stand-Alone Version For Calculating MR of Top Correlating Genes ####
  #### Written by Elly Poretsky on 12.27.19, for a more efficient mutRank      ####
  # datm <- read.table("data/NAM_qTeller_ALL_4.csv", header=T,  row.names=1, sep=",")
  # genes <- "GRMZM2G090245 GRMZM2G048073 GRMZM2G169261"
  # Remove all Zero-Sum rows, they interfere with cor, but mostly as a warning
  datm <- as.data.frame(datm)
  datm <- datm[rowSums(datm[, -1])>0, ]
  gene_list <- unlist(strsplit(genes, "[[:space:]]+"))
  # Find the index row of the gene, use 'tolower' incase user misspelled
  # gene_index <- match(tolower(gene_list), tolower(row.names(datm)))
  # 
  if(length(gene_list)==1){
    # Calculate the correlation of the gene with all other genes
    gene_cors <- cor(t(datm[gene_list,]), t(datm))^2
    # Reorder the correlation values so highest one is first
    gene_cors_ordered <- order(gene_cors, decreasing=T)
    # Select the top n correlating genes to calculate MR on
    genes_for_mr <- datm[gene_cors_ordered[1:n],]
  } else {
    genes_for_mr <- datm[gene_list,]
  }
  #gene_index <- match(tolower(gene_list), tolower(row.names(datm)))
  # Calculate the correlation between selected genes and complete data
  cor_for_mr <- cor(t(genes_for_mr),t(datm))
  # Fast rank all the columns that contain the selected genes
  rank_for_mr <- apply(cor_for_mr, 1, frankv, order=-1) # maybe replace frankv from final version?
  # Add row names since it is lost in the 'cor' function
  row.names(rank_for_mr) <- row.names(datm)
  # Remove all the rows except for the selected genes
  rank_for_mr <- rank_for_mr[row.names(genes_for_mr),]
  # Calculate the MR values between selected genes
  mr <- sqrt(rank_for_mr*t(rank_for_mr))
  ifelse(mrVScor=="Correlation", return(cor_for_mr[,row.names(genes_for_mr)]), return(mr))
}

df_annotator <- function(df, annotations){
  annotations <- as.vector(annotations[, "annotation"][match(tolower(rownames(df)), tolower(rownames(annotations)))])
  df <- cbind(df,annotations)
  return(df)
}

df_output_editor <- function(df,firstColumn,mrVScor,order,round){
  if(order & mrVScor=="MR"){df<-df[order(df[,1]),]}
  if(firstColumn){df<-df[,1, drop=FALSE]} 
  if(round){df<-round(df, 0)}
  #colnames(df) <- sprintf('<div style="transform:rotate(90deg);margin-top:10px;">%s</div>', colnames(df))
  return(df)
}