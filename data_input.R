?column
?tabPanel
?fixedPanel
?sidebarPanel
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
                    choices = c("",list.files("annotations/", pattern='*.csv'))),
        fileInput(ns("annotation_submit"), NULL, multiple = FALSE, width="100%",
                  accept = c("text/csv","text/comma-separated-values,text/plain",".csv")), 
        hr(),
        selectInput(ns("symbol_list"), "Load gene symbols:", selectize=F, width="100%",
                    choices = c("",list.files("symbols/", pattern='*.csv'))),
        fileInput(ns("symbol_submit"), NULL, multiple = FALSE, width="100%",
                  accept = c("text/csv","text/comma-separated-values,text/plain",".csv")), 
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
        selectInput(ns("gene_domains"), "Load a gene-specific domain file:", selectize=F, width="100%",
          choices = c("",list.files("domains/", pattern='*.csv'))),
        fileInput(ns("submit_quantitative"), NULL, multiple = FALSE, width="100%",
          c("text/csv","text/comma-separated-values,text/plain",".csv")),
        fileInput(ns("submit_protein_fasta"), "Submit full protein FASTA file (.fa)", multiple = FALSE, width="100%",
          c("text/csv","text/comma-separated-values,text/plain",".csv")),
        actionButton(ns("download_gene_domains"), "Download predicted domain file")
        
       ),
       column(6, style = "",
        p("Don't have a gene-specific domain file?", style = "font-weight: bold;"),
        p("You can upload a protein fasta file and mutRank will generate a gene-specific domain file 
          that you can download for future use. We use hmmscan with default settings on the full Pfam
          database of hmm profiles to predict domains in your fasta file. The longest protein sequence
          for each gene used for domain predictions if the protein names in the fasta file contain _P001.")
       )
     )
    )
  )
)}

dataInput <- function(input, output, session) {
  ns <- session$ns
  # Different input functions, should try to de-clatter it in the near future              ###
  
  expression <- reactiveVal(data.frame())
  observe({if(input$expression_list==""){return(data.frame())}
    data <- data_loader(read.csv(paste("data/", input$expression_list, sep=""), header=T))
    output$expression_name <- renderText(paste("Selected table size: ", toString(dim(data))))
    expression(data)})
  observe({if(is.null(input$expression_submit)){return(data.frame())}
    data <- data_loader(read.csv(input$expression_submit$datapath, header = T))
    output$expression_name <- renderText(paste("Uploaded table size: ", toString(dim(data))))
    expression(data)})
  
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
  
  annotations <- reactiveVal()
  observe({if(input$annotation_list==""){return(annotations(data.frame(annotation=c(NA))))}
    annotations(read.table(paste("annotations/", input$annotation_list, sep=""), sep="\t", header=T, row.names=1, quote=""))})
  observe({if(is.null(input$annotation_submit)){return(data.frame())}
    annotations(read.table(paste(input$annotation_submit$datapath, sep=""), sep="\t", header=T, row.names=1, quote=""))})
  
  GO_db <- reactiveVal()
  observe({if(input$GO_db==""){return(GO_db(data.frame(onthology=c(NA))))}
    onthology <- get_ontology(paste("GO/", input$GO_db,sep=""))$name
    GO_db(as.data.frame(onthology))})
  observe({if(is.null(input$GO_db_submit)){return(data.frame())}
    onthology <- get_ontology(input$GO_db_submit$datapath)$name
    GO_db(as.data.frame(onthology))})
  
  GO_genes <- reactiveVal()
  observe({if(input$GO_genes==""){return(GO_genes(data.frame(term_accession=c(NA))))}
    GO_genes(read.table(paste("GO/",input$GO_genes,sep=""),header=T))})
  observe({if(is.null(input$GO_genes_submit)){return(data.frame())}
    GO_genes(read.table(input$GO_genes_submit$datapath, header=T))})
  
  return(list(expression=expression,
              annotations=annotations,
              GO_db=GO_db,
              GO_genes=GO_genes,
              foldchange=foldchange, 
              association=association))
  
  
  #return(list(expression=expression,annotations=annotations,symbols=symbols,))
}

data_loader <- function(temp_df){
  # User has not uploaded a file yet
  duplicated <- duplicated(temp_df[,1])
  deduped <- temp_df[!duplicated,-1]
  row.names(deduped) <- temp_df[!duplicated,1]
  new_data <- as.matrix(deduped)
  return(new_data)
}


#p('Gene expression data is the only necessary component for the coexpression analysis. 
#          Mutual-rank (MR) is used as a measure of coexpression instead of Pearson correlation coeeficient (PCC).
#          After the initial MR analysis, additional types of information can be overlayed for more informative results.'),
