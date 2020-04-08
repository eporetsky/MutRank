?selectInput
mutualRankUI <- function(id){
  ns <- NS(id)
  
  tabPanel("Mutual Rank",
     sidebarPanel(
       width = 4,
       selectInput(ns("reference_gene_method"), "Select reference gene method:", selectize=F, width="100%",
                   choices = c("Single reference gene", "Compound reference gene", "Reference gene list"),
                   selected="Single reference gene"),
       p(), hr(), p(),
       uiOutput(ns("reference_method")),
       p(),
       actionButton(ns("update_button"), "Calculate MR Values"),
       p(), hr(),
       checkboxInput(ns("firstColumn"), "Only show first column?", value = T, width = NULL),
       checkboxInput(ns("order_coexpression"), "Order by MR?", value = T, width = NULL),
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

mutualRank <- function(input, output, session, expression, annotations, symbols, 
                       categories, foldchange, go_mapping, domain_mapping){
  
  ns <- session$ns
  
  coexpression <- reactive({
    input$update_button
    coexpression_table(expression(), 
                       isolate(input$reference_gene),
                       isolate(input$reference_gene_method),
                       isolate(input$num_top_pcc),
                       isolate(input$compound_method))
  })
  
  output$reference_method <- renderUI({
    if(input$reference_gene_method=="Single reference gene"){
      list(textInput(ns('reference_gene'), "Reference gene ID:","GRMZM2G085381"),# Bx1 for reference
      numericInput(ns("num_top_pcc"), "Number of genes for coexpression:", 200))
    } else{
    if(input$reference_gene_method=="Compound reference gene"){
      list(selectInput(ns("compound_method"), "Choose compounding method:", choices=c("Sum","Average","Max","Min"), selected="Sum"),
      textInput(ns('reference_gene'), "Reference gene IDs to compound:","GRMZM2G085381\tGRMZM2G085054"),
      numericInput(ns("num_top_pcc"), "Number of genes for coexpression:", 200)) # Bx1 for reference
    } else{
    if(input$reference_gene_method=="Reference gene list"){
      textInput(ns('reference_gene'), "List of reference genes:","GRMZM2G085381\tGRMZM2G085054") # Bx1 and Bx8 for reference
    }}}
  })

  coexpression_df_prefs <- reactive({
    coexpression <- coexpression()
    if(input$order_coexpression){coexpression <- order_coexpression_table(coexpression)}
    df_output_editor(coexpression,input$firstColumn,input$round)
  })
  
  coexpression_df_prefs_annotations <- reactive({
    output_df <- coexpression_df_prefs()
    if(input$symbols){output_df <- df_add_symbols(output_df,symbols())}
    if(input$categories){output_df <- df_add_categories(output_df, categories(),go_mapping(),domain_mapping())}
    if(input$foldchange){output_df <- df_add_foldchange(output_df,foldchange())}
    if(input$annotate){output_df <- df_annotator(output_df,annotations())}
    
    
    #output_df <- output_df[, !apply(is.na(output_df), 2, all),drop=FALSE] # remove empty columns
    return(output_df)
  })
  
  output$coexpression_table <- renderTable(
    output_df <- coexpression_df_prefs_annotations(),
    class = 'cell-border stripe nowrap compact dt-center',
    bordered = TRUE,
    rownames = TRUE
  )
  
  output$download_table <- downloadHandler(
    filename = function(){"mutRank_table.tsv"}, 
    content = function(fname){write.table(coexpression_df_prefs_annotations(), fname, sep="\t", quote=F, col.names=NA)}
  )
  
  return(coexpression)
}
