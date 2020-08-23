?selectInput
mutualRankUI <- function(id){
  ns <- NS(id)
  
  tabPanel("Mutual Rank",
     sidebarPanel(
       width = 4,
       selectInput(ns("reference_gene_method"), "Select reference gene method:", selectize=F, width="100%",
                   choices = c("Single reference gene", "Compound reference gene", "Reference gene list"),
                   selected="Single reference gene"),
       p(),
       uiOutput(ns("reference_method")),p(),
       actionButton(ns("update_button"), "Calculate MR Values"),p(),
       textOutput(ns("missing_genes")),
       p(), hr(),
       checkboxInput(ns("firstColumn"), "Only show first column?", value = T, width = NULL),
       checkboxInput(ns("round"), "Round to nearest integer", value = T, width = NULL),
       checkboxInput(ns("annotate"), "Add gene annotations", value = T, width = NULL),
       checkboxInput(ns("symbols"), "Add gene symbols", value = T, width = NULL),
       checkboxInput(ns("categories"), "Add custom categories", value = T, width= NULL),
       checkboxInput(ns("foldchange"), "Add foldchange values", value = T, width= NULL),
       p(), hr(), p(),
       downloadButton(ns('download_table'),"Download table")
     ),
     mainPanel(
      tags$head(tags$style(HTML("#mutualRankNS-coexpression_table table tr td {white-space: nowrap;}"))),
      # Suppress error message
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      tableOutput(ns("coexpression_table"))
      )
    )
}

mutualRank <- function(input, output, session, 
                       expression, annotations, symbols, categories, foldchange, go_mapping, domain_mapping){
  
  ns <- session$ns
  
  
  ref_genes <- reactive({
    input$update_button
    # Tries to convert the input gene(s) to a vector of gene_list and remove duplicates
    gene_list <- isolate(unique(unlist(strsplit(input$reference_gene, "[[:space:],]+"))))
    # If the selected method is a single reference gene it returns the first gene even if given a gene list
    ifelse(input$reference_gene_method=="Single reference gene",return(gene_list[1]),return(gene_list))
  })
  reference_genes <- reactive({ref_genes()[ref_genes() %in% row.names(expression())]}) # Returns the genes present in the expression data
  missing_genes <- reactive({ref_genes()[!ref_genes() %in% row.names(expression())]})  # Returns the genes not present in the expression data
    
  # The main reactive function to calculate to Mutual Rank table
  coexpression <- reactive({
    input$update_button
    # Set the maximum number of of top PCC coexpressed genes to 500
    num_top_pcc <- input$num_top_pcc
    if(num_top_pcc>500){num_top_pcc<-500}
    # Calculatue the MR coexpression table using the input parameters
    coexpression_table(expression(), 
                       # Values below are isolated so it is only reactive to input$update_button 
                       isolate(reference_genes()),
                       isolate(input$reference_gene_method),
                       isolate(num_top_pcc),
                       isolate(input$compound_method))
  })
  
  output$reference_method <- renderUI({
    if(input$reference_gene_method=="Single reference gene"){
      list(textInput(ns('reference_gene'), "Reference gene ID:","GRMZM2G085381"),# Bx1 for reference
      numericInput(ns("num_top_pcc"), "Number of genes for coexpression (Max 500):", 200))
    } else{
    if(input$reference_gene_method=="Compound reference gene"){
      list(selectInput(ns("compound_method"), "Choose compounding method:", choices=c("Sum","Average","Max","Min"), selected="Sum"),
      textInput(ns('reference_gene'), "Reference gene IDs to compound:","GRMZM2G085381\tGRMZM2G085054"), # Bx1 and Bx8 for reference
      numericInput(ns("num_top_pcc"), "Number of genes for coexpression (Max 500):", 200))
    } else{
    if(input$reference_gene_method=="Reference gene list"){
      textInput(ns('reference_gene'), "List of reference genes:","GRMZM2G085381\tGRMZM2G085054") # Bx1 and Bx8 for reference
    }}}
  })

  coexpression_df_prefs <- reactive({
    # Calls the Mutual Rank coexpression reative table to filter all columns except first and/or round MR values
    df_output_editor(coexpression(),input$firstColumn,input$round)
  })
  
  coexpression_df_prefs_annotations <- reactive({
    output_df <- coexpression_df_prefs()
    if(input$symbols)   {output_df <- df_add_symbols(output_df,symbols())}
    if(input$categories){output_df <- df_add_categories(output_df, categories(),go_mapping(),domain_mapping())}
    if(input$foldchange){output_df <- df_add_foldchange(output_df,foldchange())}
    if(input$annotate)  {output_df <- df_annotator(output_df,annotations())}
    return(output_df)
  })
  
  output$coexpression_table <- renderTable(
    output_df <- coexpression_df_prefs_annotations(),
    class = 'cell-border stripe nowrap compact dt-center',
    bordered = TRUE,
    rownames = TRUE
  )
  
  # Prints the list of genes not found in the expression table
  output$missing_genes <- renderText(paste("Genes not found: ", toString(missing_genes()))) 

  output$download_table <- downloadHandler(
    filename = function(){"mutRank_table.tsv"}, 
    content = function(fname){write.table(coexpression_df_prefs_annotations(), fname, sep="\t", quote=F, col.names=NA)}
  )
  
  return(coexpression)
}
