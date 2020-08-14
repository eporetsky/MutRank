dataInputUI <- function(id) {
  ns <- NS(id)
  tabPanel("Data Input",
    sidebarLayout(position="left", fluid=F,
      sidebarPanel(width = 4,
        p("Use the dropdown menu to select the expression data and supporting data located in the home 
          directory where the MutRank app is located. The selected files will be loaded from the appropriate subfolders."),
        selectInput(ns("expression_list"), "Load expression data (.csv):", selectize=F, width="100%",
          choices = c("",list.files("data/", pattern='*.csv')),selected=defaults["data",]),
        textOutput({ns("expression_name")}),
      hr(),
        selectInput(ns("annotation_list"), "Load gene annotations (.tsv):", selectize=F, width="100%",
                    choices = c("",list.files("annotations/", pattern='*.tsv')),selected=defaults["annotations",]),
        selectInput(ns("symbol_list"), "Load gene symbols (.tsv):", selectize=F, width="100%",
                    choices = c("",list.files("symbols/", pattern='*.tsv')), selected=defaults["symbols",]),
      p(),hr(),p(),
        actionButton(ns("save_default"), "Remember Selected Files")
      ),
       
      mainPanel(width=8,
       fluidRow(
         column(6, style = "",
          selectInput(ns("foldchange"), "Load DEG File (.csv):", selectize=F, width="100%",
            choices = c("",list.files("foldchange/", pattern='*.csv')), selected=defaults["foldchange",]),
          p(),hr(),p(),
          selectInput(ns("category_list"), "Load Custom Categories (.tsv):", selectize=F, width="100%",
            choices = c("",list.files("categories/", pattern='*.tsv')), selected=defaults["categories",]),
          p(),hr(),p(),
          selectInput(ns("domain_list"), "Load Protein Domains (.tsv):", selectize=F, width="100%",
                      choices = c("",list.files("domains/", pattern='*.tsv')),selected=defaults["domains",]),
          p(),hr(),p(),
          selectInput(ns("GO_db"), "Load GO Database (.obo):", selectize=F, width="100%",
                      choices = c("",list.files("GO/", pattern='*.obo')), selected=defaults["GO_db",]),
          selectInput(ns("GO_genes"), "Load GO Annotations (.tsv):", selectize=F, width="100%",
                      choices = c("",list.files("GO/", pattern='*.tsv')),selected=defaults["GO_genes",]),
         )
       )
     )
))}

dataInput <- function(input, output, session) {
  ns <- session$ns
  
  # Saves the selected files in all input fields as the default
  observeEvent(input$save_default,{
    defaults["data",]<-input$expression_list
    defaults["annotations",]<-input$annotation_list
    defaults["symbols",]<-input$symbol_list
    defaults["foldchange",]<-input$foldchange
    defaults["categories",]<-input$category_list
    defaults["GO_db",]<-input$GO_db
    defaults["GO_genes",]<-input$GO_genes
    defaults["domains",]<-input$domain_list
    write.table(defaults,"default_files.csv",quote=F,sep=",", col.names=NA)
  })
  
  # Below are the functions to load the different tables into memory
  # The tables are loaded into a Shiny reactiveVal and changes in the selected file are monitored by observe()
  # If no file is specified the reactiveVal becomes an empty data.frame 
  ### Support for proper file format testing is not implemented - Follow example file format ###
  expression <- reactiveVal()
  observe({if(input$expression_list==""){return(data.frame())}
    data <- data_loader(read.csv(paste("data/", input$expression_list, sep=""), header=T))    # Load the expression data from CSV file
    output$expression_name <- renderText(paste("Selected table size: ", toString(dim(data)))) # Print the expression data dimensions
    expression(data[rowSums(data[, -1])>0, ])})                                               # Removes row with 0 sum and save into reactiveVal
  
  annotations <- reactiveVal()
  observe({if(input$annotation_list==""){return(annotations(data.frame(annotation=c(NA))))}
    annotations(read.table(paste("annotations/", input$annotation_list, sep=""), sep="\t", header=T, row.names=1, quote=""))})
  
  symbols <- reactiveVal()
  observe({if(input$symbol_list==""){return(symbols(data.frame(symbol=c(NA))))}
    symbols(read.table(paste("symbols/", input$symbol_list, sep=""), sep="\t", header=T, row.names=1, quote=""))})
  
  foldchange <- reactiveVal()
  observe({if(input$foldchange==""){foldchange(data.frame(FC=c(1),row.names=c("a")))} else{
    foldchange(read.table(paste("foldchange/", input$foldchange, sep=""), row.names=1, header=T,sep=","))}})
  
  GO_db <- reactiveVal()
  observe({if(input$GO_db==""){return(GO_db(data.frame(onthology=c(NA))))}
    onthology <- get_ontology(paste("GO/", input$GO_db,sep=""))$name
    GO_db(as.data.frame(onthology))})
  
  GO_genes <- reactiveVal()
  observe({if(input$GO_genes==""){return(GO_genes(data.frame(genes=c(NA),term_accession=c(NA))))}
    GO_genes(wide2long(read.table(paste("GO/",input$GO_genes,sep=""),header=T, sep="\t")))})
  
  domains <- reactiveVal()
  observe({if(input$domain_list==""){return(domains(data.frame(genes=c(NA),domains=c(NA))))}
    domains(wide2long(read.table(paste("domains/",input$domain_list,sep=""),header=T, sep="\t")))})
  
  categories <- reactiveVal()
  observe({if(input$category_list==""){return(categories(data.frame(symbol=c(NA))))}
    categories(read.table(paste("categories/", input$category_list, sep=""), sep="\t", row.names=NULL, header=T, quote=""))})
  
  return(list(expression=expression,
              annotations=annotations,
              symbols=symbols,
              GO_db=GO_db,
              GO_genes=GO_genes,
              domains=domains,
              foldchange=foldchange, 
              categories=categories))
}