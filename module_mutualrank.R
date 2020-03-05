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
       checkboxInput(ns("symbols"), "Add gene symbols", value = T, width = NULL),
       hr()
     ),
     mainPanel(
      tags$head(tags$style(HTML("#mutualRankNS-coexpression_table table tr td {white-space: nowrap;}"))),
      tableOutput(ns("coexpression_table"))
      )
    )
}

mutualRank <- function(input, output, session, data, annotations, symbols, categories, go_mapping, domain_mapping) {
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
  
  coexpression_df_prefs_annotations <- reactive({
    output_df <- coexpression_df_prefs()
    if(input$symbols){output_df <- df_add_symbols(output_df,symbols())}
    if(TRUE){output_df <- df_add_categories(output_df, categories(),go_mapping(),domain_mapping())}
    if(input$annotate){output_df <- df_annotator(output_df,annotations())}
    return(output_df)
  })
  
  output$coexpression_table <- renderTable(
    output_df <- coexpression_df_prefs_annotations(),
    class = 'cell-border stripe nowrap compact dt-center',
    bordered = TRUE,
    rownames = TRUE
  )
  
  coexpression_return <- reactive({
    df <- as.data.frame(coexpression()) # column reordering doesn't work on matrix
    if(input$order){
      df <- df[order(df[,1]),] # order rows from lowest MR values
      df <- df[rownames(df)]   # reorder column with rowname order
    }
    return(df)
  })
  
  return(coexpression_return)
}

coexpression_table <- function(datm, genes, n, mrVScor){
  #### Working Stand-Alone Version For Calculating MR of Top Correlating Genes ####
  #### Written by Elly Poretsky on 12.27.19, for a more efficient mutRank      ####
  # Remove all Zero-Sum rows, they interfere with cor, but mostly as a warning
  datm <- as.data.frame(datm)
  datm <- datm[rowSums(datm[, -1])>0, ]
  gene_list <- unlist(strsplit(genes, "[[:space:]]+"))
  # Find the index row of the gene, use 'tolower' incase user misspelled
  # mutRank has an option to calculate MR between a list of genes as input
  # If gene_list==1 it first finds top coexpressed genes based on PCC and
  # uses them as the list of genes to calculate MR for
  if(length(gene_list)==1){
    # Calculate the correlation of the gene with all other genes
    gene_cors <- cor(t(datm[gene_list,]), t(datm))^2
    # Reorder the correlation values so highest one is first
    gene_cors_ordered <- order(gene_cors, decreasing=T)
    # Select the top n correlating genes to calculate MR on
    genes_for_mr <- datm[gene_cors_ordered[1:n],]
  } else {
    genes_for_mr <- datm[gene_list,]
  } # Used to automatically support tolower, worth returning: match(tolower(gene_list), tolower(row.names(datm)))
  
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

# Add an annotation column to the MR table
df_annotator <- function(df, annotations){
  annotations <- as.vector(annotations[,1][match(tolower(rownames(df)), tolower(rownames(annotations)))])
  return(cbind(df,annotations=annotations))
}

# Add a symbol coumn to the MR table
df_add_symbols <- function(df, symbols){
  symbols <- as.vector(symbols[,1][match(tolower(rownames(df)), tolower(rownames(symbols)))])
  return(cbind(df,symbols=symbols))
}

# Add a column for each category in the MR table
df_add_categories <- function(df, categories, go_mapping, domain_mapping){
  genes <- tolower(rownames(df))
  for(category in colnames(categories)){
    matched_gos <- match(genes, tolower(unique(as.vector(go_mapping[go_mapping[,2] %in% categories[,category],][,1]))),nomatch=0)
    matched_domains <-  match(genes, tolower(unique(as.vector(domain_mapping[domain_mapping[,2] %in% categories[,category],][,1]))),nomatch=0)
    matched_categories <- match(genes, tolower(categories[,category]),nomatch=0)
    all_matched <- matched_gos + matched_domains + matched_categories
    all_matched[all_matched==0] <- NA
    all_matched[!is.na(all_matched)] <- "Y"
    df <- cbind(df,all_matched)
    colnames(df)[length(colnames(df))]<-category
  }
  return(df)
}

# Edit the MR table output based on Shiny inputs
df_output_editor <- function(df,firstColumn,mrVScor,order,round){
  if(order & mrVScor=="MR"){df<-df[order(df[,1]),]}
  if(firstColumn){df<-df[,1, drop=FALSE]} 
  if(round){df<-round(df, 0)}
  return(df)
}
