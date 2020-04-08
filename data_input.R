dataInputUI <- function(id) {
  ns <- NS(id)
  tabPanel("Data Input",
    tags$head(tags$style(HTML("hr {border-top: 1px solid #95a5a6; margin-bottom:0px; margin-top:0px}
                               .form-group {margin-bottom: 2px; margin-top: 2px;}"))),
    sidebarLayout(position="left", fluid=F,
      sidebarPanel(width = 4,
        p('Load the expression data and support data to start using mutRank. You can select the files 
          located in the app folder from the dropdown menu.'),
        
        selectInput(ns("expression_list"), "Load expression data:", selectize=F, width="100%",
          choices = c("",list.files("data/", pattern='*.csv')),selected=defaults["data",]),
        fileInput(ns("expression_submit"), NULL, multiple = FALSE, width="100%",
          accept = c("text/csv","text/comma-separated-values,text/plain",".csv")), 
        textOutput({ns("expression_name")}),
      hr(),
        selectInput(ns("annotation_list"), "Load gene annotations:", selectize=F, width="100%",
                    choices = c("",list.files("annotations/", pattern='*.tsv')),selected=defaults["annotations",]),
        fileInput(ns("annotation_submit"), NULL, multiple = FALSE, width="100%",
                  accept = c("text/csv","text/comma-separated-values,text/plain",".tsv")), 
      
        selectInput(ns("symbol_list"), "Load gene symbols:", selectize=F, width="100%",
                    choices = c("",list.files("symbols/", pattern='*.tsv')), selected=defaults["symbols",]),
        fileInput(ns("symbol_submit"), NULL, multiple = FALSE, width="100%",
                  accept = c("text/tsv","text/tab-separated-values,text/plain",".tsv")), 
        ),
     
    mainPanel(width=8,
     fluidRow(
       column(6, style = "",
        selectInput(ns("foldchange"), "Load fold-change data file:", selectize=F, width="100%",
          choices = c("",list.files("foldchange/", pattern='*.csv')), selected=defaults["foldchange",]),
        fileInput(ns("foldchange_submit"), NULL, multiple = FALSE, width="100%",
          c("text/csv","text/comma-separated-values,text/plain",".csv"))
       ),
       column(6, style = "",
              selectInput(ns("category_list"), "Annotate using custom categories:", selectize=F, width="100%",
                          choices = c("",list.files("categories/", pattern='*.tsv')), selected=defaults["categories",]),
              fileInput(ns("category_submit"), NULL, multiple = FALSE, width="100%",
                        c("text/csv","text/comma-separated-values,text/plain",".tsv"))
       )
     ),
     hr(), p(),
     fluidRow(
       column(6, style = "",
          selectInput(ns("GO_db"), "Load GO database file:", selectize=F, width="100%",
            choices = c("",list.files("GO/", pattern='*.obo')), selected=defaults["GO_db",]),
          fileInput(ns("GO_db_submit"), NULL, multiple = FALSE, width="100%",
            c("text/csv","text/comma-separated-values,text/plain",".obo"))
       ),
       column(6, style = "",
        selectInput(ns("GO_genes"), "Load GO for genes:", selectize=F, width="100%",
          choices = c("",list.files("GO/", pattern='*.tsv')),selected=defaults["GO_genes",]),
        fileInput(ns("GO_genes_submit"), NULL, multiple = FALSE, width="100%",
          c("text/csv","text/comma-separated-values,text/plain",".tsv"))
       )
     ),
     hr(), p(),
     fluidRow(
       column(6, style = "",
        selectInput(ns("domain_list"), "Load a gene-specific domain file:", selectize=F, width="100%",
          choices = c("",list.files("domains/", pattern='*.tsv')),selected=defaults["domains",]),
        fileInput(ns("domain_submit"), NULL, multiple = FALSE, width="100%",
          c("text/csv","text/comma-separated-values,text/plain",".csv"))
       ),
       column(6, style = "",
              p('Press the button below to save the selected files for the next time you run mutRank', style = "font-weight: bold;"),
              actionButton(ns("save_default"), "Save Default")
       )
     )
    )
  )
)}

dataInput <- function(input, output, session) {
  ns <- session$ns
  
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
  observe({if(input$foldchange==""){foldchange(data.frame(FC=c(1),row.names=c("a")))} else{
    foldchange(read.table(paste("foldchange/", input$foldchange, sep=""), row.names=1, header=T,sep=","))}})
  observe({if(is.null(input$foldchange_submit)){return(data.frame(FC=c(1),row.names=c("a")))}
    foldchange(read.table(input$foldchange_submit$datapath, row.names=1, header = T,sep=",",row.names=1))})
  
  GO_db <- reactiveVal()
  observe({if(input$GO_db==""){return(GO_db(data.frame(onthology=c(NA))))}
    onthology <- get_ontology(paste("GO/", input$GO_db,sep=""))$name
    GO_db(as.data.frame(onthology))})
  observe({if(is.null(input$GO_db_submit)){return(data.frame())}
    onthology <- get_ontology(input$GO_db_submit$datapath)$name
    GO_db(as.data.frame(onthology))})
  
  GO_genes <- reactiveVal()
  observe({if(input$GO_genes==""){return(GO_genes(data.frame(genes=c(NA),term_accession=c(NA))))}
    GO_genes(wide2long(read.table(paste("GO/",input$GO_genes,sep=""),header=T, sep="\t")))})
  observe({if(is.null(input$GO_genes_submit)){return(data.frame())}
    GO_genes(read.table(input$GO_genes_submit$datapath, header=T))})
  
  domains <- reactiveVal()
  observe({if(input$domain_list==""){return(domains(data.frame(genes=c(NA),domains=c(NA))))}
    domains(wide2long(read.table(paste("domains/",input$domain_list,sep=""),header=T, sep="\t")))})
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