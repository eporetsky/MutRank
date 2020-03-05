dataInputUI <- function(id) {
  ns <- NS(id)
  useShinyjs()
  tabPanel("Data Input",
    tags$head(tags$style(HTML("hr {border-top: 1px solid #95a5a6; margin-bottom:0px; margin-top:0px}
                               .form-group {margin-bottom: 2px; margin-top: 2px;}"))),
    sidebarLayout(position="left", fluid=F,
      sidebarPanel(width = 4,
        p('To begin the coexpression analysis load a gene expression data set using one of the options below.
          Data can be loaded from a list of ".tsv" files in the /data folder or directly uploaded from your computer.
          A text output showing the number of rows and columns will updated when expression data is loaded.'),
        selectInput(ns("expression_list"), "Load expression data:", selectize=F, width="100%",
          choices = c("",list.files("data/", pattern='*.csv'))),
        fileInput(ns("expression_submit"), NULL, multiple = FALSE, width="100%",
          accept = c("text/csv","text/comma-separated-values,text/plain",".csv")), 
        textOutput({ns("expression_name")}),
      hr(),
        selectInput(ns("annotation_list"), "Load gene annotations:", selectize=F, width="100%",
                    choices = c("",list.files("annotations/", pattern='*.tsv'))),
        fileInput(ns("annotation_submit"), NULL, multiple = FALSE, width="100%",
                  accept = c("text/csv","text/comma-separated-values,text/plain",".tsv")), 
        hr(),
        selectInput(ns("symbol_list"), "Load gene symbols:", selectize=F, width="100%",
                    choices = c("",list.files("symbols/", pattern='*.tsv'))),
        fileInput(ns("symbol_submit"), NULL, multiple = FALSE, width="100%",
                  accept = c("text/tsv","text/tab-separated-values,text/plain",".tsv")), 
        ),
     
    mainPanel(width=8,
     fluidRow(
       column(6, style = "",
        selectInput(ns("foldchange"), "Load fold-change data file:", selectize=F, width="100%",
          choices = c("",list.files("foldchange/", pattern='*.csv'))),
        fileInput(ns("foldchange_submit"), NULL, multiple = FALSE, width="100%",
          c("text/csv","text/comma-separated-values,text/plain",".csv"))
       ),
       column(6, style = "",
        selectInput(ns("association"), paste("Load association data file:"), selectize=F, width="100%",
          choices = c("",list.files("association/", pattern='*.csv'))),
        fileInput(ns("association_submit"), NULL, multiple = FALSE, width="100%",
          c("text/csv","text/comma-separated-values,text/plain",".csv"))
       )
     ),
     hr(),
     fluidRow(
       column(6, style = "",
          selectInput(ns("GO_db"), "Load GO database file:", selectize=F, width="100%",
            choices = c("",list.files("GO/", pattern='*.obo'))),
          fileInput(ns("GO_db_submit"), NULL, multiple = FALSE, width="100%",
            c("text/csv","text/comma-separated-values,text/plain",".csv"))
       ),
       column(6, style = "",
        selectInput(ns("GO_genes"), "Load GO for genes:", selectize=F, width="100%",
          choices = c("",list.files("GO/", pattern='*.tsv'))),
        fileInput(ns("GO_genes_submit"), NULL, multiple = FALSE, width="100%",
          c("text/csv","text/comma-separated-values,text/plain",".tsv"))
       )
     ),
     hr(),
     fluidRow(
       column(6, style = "",
          selectInput(ns("category_list"), "Annotate using custom categories:", selectize=F, width="100%",
            choices = c("",list.files("categories/", pattern='*.tsv'))),
          fileInput(ns("category_submit"), NULL, multiple = FALSE, width="100%",
            c("text/csv","text/comma-separated-values,text/plain",".tsv")),
          #p("", style = "font-weight: bold;"),
          p("The presence of every coexpressed genes, and their respective GO terms and Pfam domains, will be examined in each category.
            "),
          p("Genes with the current 1. Diamond, 2. Star, 3. Triangle, 4. Down Triangle, 5. Square. Column 6 and above won't 
            change the network vertices shapes but will appear in the table.")
          
       ),

       column(6, style = "",
          selectInput(ns("domain_list"), "Load a gene-specific domain file:", selectize=F, width="100%",
              choices = c("",list.files("domains/", pattern='*.tsv'))),
          fileInput(ns("domain_submit"), NULL, multiple = FALSE, width="100%",
              c("text/csv","text/comma-separated-values,text/plain",".csv")),
          fileInput(ns("submit_protein_fasta"), "Submit full protein FASTA file (.fa)", multiple = FALSE, width="100%",
              c("text/csv","text/comma-separated-values,text/plain",".csv")),
          actionButton(ns("download_gene_domains"), "Download predicted domain file")
              
       ),
     )
     
    )
  )
)}

dataInput <- function(input, output, session) {
  ns <- session$ns
  
  expression <- reactiveVal(data.frame())
  observe({if(input$expression_list==""){return(data.frame())}
    data <- data_loader(read.csv(paste("data/", input$expression_list, sep=""), header=T))
    output$expression_name <- renderText(paste("Selected table size: ", toString(dim(data))))
    expression(data)})
  observe({if(is.null(input$expression_submit)){return(data.frame())}
    data <- data_loader(read.csv(input$expression_submit$datapath, header = T))
    output$expression_name <- renderText(paste("Uploaded table size: ", toString(dim(data))))
    expression(data)})
  
  annotations <- reactiveVal()
  observe({if(input$annotation_list==""){return(annotations(data.frame(annotation=c(NA))))}
    annotations(read.table(paste("annotations/", input$annotation_list, sep=""), sep="\t", header=T, row.names=1, quote=""))})
  observe({if(is.null(input$annotation_submit)){return(data.frame())}
    annotations(read.table(paste(input$annotation_submit$datapath, sep=""), sep="\t", header=T, row.names=1, quote=""))})
  
  symbols <- reactiveVal()
  observe({if(input$symbol_list==""){return(symbols(data.frame(symbol=c(NA))))}
    symbols(read.table(paste("symbols/", input$symbol_list, sep=""), sep="\t", header=T, row.names=1, quote=""))})
  observe({if(is.null(input$symbol_submit)){return(data.frame())}
    symbols(read.table(paste(input$symbol_submit$datapath, sep=""), sep="\t", header=T, row.names=1, quote=""))})
  
  foldchange <- reactiveVal()
  observe({if(input$foldchange==""){return(data.frame(GeneID=c(1),ZmPep3=c(1)))}
    foldchange(read.table(paste("foldchange/", input$foldchange, sep=""), header=T,sep=","))})
  observe({if(is.null(input$foldchange_submit)){return(data.frame(GeneID=c(2),ZmPep3=c(2)))}
    foldchange(read.table(input$foldchange_submit$datapath, header = T,sep=",",row.names=1))})
  
  association <- reactiveVal(data.frame())
  observe({if(input$association==""){return(data.frame())}
    association(read.table(paste("association/", input$association, sep=""), header=T,sep=",",row.names=1))})
  observe({if(is.null(input$association_submit)){return(data.frame())}
    association(read.table(input$association_submit$datapath, header = T,sep=",",row.names=1))})
  
  GO_db <- reactiveVal()
  observe({if(input$GO_db==""){return(GO_db(data.frame(onthology=c(NA))))}
    onthology <- get_ontology(paste("GO/", input$GO_db,sep=""))$name
    GO_db(as.data.frame(onthology))})
  observe({if(is.null(input$GO_db_submit)){return(data.frame())}
    onthology <- get_ontology(input$GO_db_submit$datapath)$name
    GO_db(as.data.frame(onthology))})
  
  GO_genes <- reactiveVal()
  observe({if(input$GO_genes==""){return(GO_genes(data.frame(genes=c(NA),term_accession=c(NA))))}
    GO_genes(read.table(paste("GO/",input$GO_genes,sep=""),header=T))})
  observe({if(is.null(input$GO_genes_submit)){return(data.frame())}
    GO_genes(read.table(input$GO_genes_submit$datapath, header=T))})
  
  domains <- reactiveVal()
  observe({if(input$domain_list==""){return(domains(data.frame(genes=c(NA),domains=c(NA))))}
    domains(read.table(paste("domains/",input$domain_list,sep=""),header=T, sep="\t"))})
  observe({if(is.null(input$domain_submit)){return(data.frame(domains=c(NA)))}
    domains(read.table(input$domain_submit$datapath, header=T))})
  
  categories <- reactiveVal()
  observe({if(input$category_list==""){return(categories(data.frame(symbol=c(NA))))}
    categories(read.table(paste("categories/", input$category_list, sep=""), sep="\t", row.names=NULL, header=T, quote=""))})
  observe({if(is.null(input$category_submit)){return(data.frame())}
    categories(read.table(paste(input$category_submit$datapath, sep=""), sep="\t", row.names=NULL, header=T, quote=""))})
  
  return(list(expression=expression,
              annotations=annotations,
              symbols=symbols,
              GO_db=GO_db,
              GO_genes=GO_genes,
              domains=domains,
              foldchange=foldchange, 
              association=association,
              categories=categories))
}

data_loader <- function(temp_df){
  # Quick function to load expression data
  duplicated <- duplicated(temp_df[,1])
  deduped <- temp_df[!duplicated,-1]
  row.names(deduped) <- temp_df[!duplicated,1]
  new_data <- as.matrix(deduped)
  return(new_data)
}