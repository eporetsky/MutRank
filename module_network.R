MRnetworkUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Network",
    sidebarPanel(
      width = 4,
      numericInput(ns("threshold"), "Select MR Threshold:",10),
      checkboxInput(ns("star_bait"), "Change bait gene to star shape?", value = T, width = NULL),
      checkboxInput(ns("quant_vertices"), "Color vertices with foldchange data?", value=T, width=NULL),
      uiOutput(ns("foldchange_columns")),
      hr(),
      uiOutput(ns("category_diamond")),
      uiOutput(ns("category_star")),
      uiOutput(ns("category_triangle")),
      uiOutput(ns("category_downtriangle")),
      uiOutput(ns("category_square"))
     ),
    mainPanel(
      visNetworkOutput(ns("network_plot"),width="1000px", height="600px")
    )
  )
}

MRnetwork <- function(input, output, session, coexpression, annotations, 
                      symbols, foldchange, categories, go_mapping, domain_mapping) {
  ns <- session$ns

  # This is the gene that will be colored red
  reference_gene <- reactive({row.names(coexpression())[1]})
  # To make a network a 
  
  raw_network <- reactive({
    # Use input coexpression table as an adjacency matrix
    adj_matrix <- as.matrix(coexpression())
    adj_matrix[adj_matrix > as.integer(input$threshold)] <- 0
    # Create the igraph network from the adjecency matrix
    igraph_network <- graph_from_adjacency_matrix(adj_matrix, mode="undirected", weighted = T,
                                                  diag = F, add.colnames = NA, add.rownames = NULL)
    igraph_network <- set_edge_attr(igraph_network,"width",value=1)          # keep uniform
    igraph_network <- set_edge_attr(igraph_network,"color", value="gray")   # change edge color based on cluster
    igraph_network <- set_vertex_attr(igraph_network,"color",value="gray")  # change vertex color based on fc
  })
  
  # Select a shape for any of the Category columns
  category_diamond <- reactive({append("None",colnames(categories()))})
  output$category_diamond <- renderUI({selectInput(ns("category_diamond"), "Category for diamond shape:", selected="None", category_diamond()) })
  category_star <- reactive({append("None",colnames(categories()))})
  output$category_star <- renderUI({selectInput(ns("category_star"), "Category for star shape:", selected="None", category_star()) })
  category_triangle <- reactive({append("None",colnames(categories()))})
  output$category_triangle <- renderUI({selectInput(ns("category_triangle"), "Category for triangle shape:", selected="None", category_triangle()) })
  category_downtriangle <- reactive({append("None",colnames(categories()))})
  output$category_downtriangle <- renderUI({selectInput(ns("category_downtriangle"), "Category for down-triangle shape:", selected="None", category_downtriangle()) })
  category_square <- reactive({append("None",colnames(categories()))})
  output$category_square <- renderUI({selectInput(ns("category_square"), "Category for square shape:", selected="None", category_square()) })
  # Change the shape of the vertices based on selected Category-shape
  reactive_category_shape <- reactive({
    vertices_names <- get.data.frame(raw_network(), what= c("vertices"))[,1]
    shapes <- rep("dot", length(vertices_names))
    if(input$category_diamond!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_diamond, "diamond", go_mapping(),domain_mapping())}
    if(input$category_star!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_star, "star", go_mapping(),domain_mapping())}
    if(input$category_triangle!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_triangle, "triangle", go_mapping(),domain_mapping())}
    if(input$category_downtriangle!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_downtriangle, "triangleDown", go_mapping(),domain_mapping())}
    if(input$category_square!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_square, "square",go_mapping(),domain_mapping())}
    return(shapes)
  })
  
  # Change the colors of vertices based on log2 fold-change data. Select column to use
  foldchange_columns <- reactive({colnames(foldchange())[1:length(colnames(foldchange()))]})
  output$foldchange_columns <- renderUI({selectInput(ns("foldchange_column"), "Choose Foldchange Column:", 
                                                     selected=foldchange_columns()[1], foldchange_columns()) })
  reactive_foldchange <- reactive({
    if(input$quant_vertices){return(foldchange_vertices(raw_network(), foldchange(), input$foldchange_column))}  
    else{return("gray")}
  })
  
  # Take all the reactive igraph parameters and apply them to the network attributes
  final_network <- reactive({
    igraph_network <- raw_network()
    # Get the name of all the vertices in the network
    vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
    # Set name of vertices to black
    annot <- as.vector(annotations()[,1][match(tolower(vertices_names), tolower(rownames(annotations())))])
    igraph_network <- set_vertex_attr(igraph_network,"annotation",value=(annot))
    if(length(E(igraph_network))>0){igraph_network <- remove.edge.attribute(igraph_network, "weight")}
    igraph_network <- set_vertex_attr(igraph_network,"color",value=reactive_foldchange())           # change vertex color based on fc
    igraph_network <- set_vertex_attr(igraph_network,"shape",value=reactive_category_shape())     # change shape for pfam/GO
    if(input$star_bait){
      vertex_color <- list(background="red",border="red")
      igraph_network <- set_vertex_attr(igraph_network,"shape",index=V(igraph_network)[1], value="star")#list(background="gray",border="red"))
    }
    new_symbols <- symbol_converter(symbols(), vertex_attr(igraph_network,'name'))
    igraph_network <- set_vertex_attr(igraph_network,"name",value=new_symbols)           # change vertex color based on fc
    return(igraph_network)
  })
  
  # Convert the igraph network instance into a visNetwork
  network_plot <- reactive({
    data <- toVisNetworkData(final_network())
    visNetwork(nodes = data$nodes, edges = data$edges, randomSeed=1) %>% 
      visIgraphLayout() %>%
      #visEdges(smooth = FALSE) %>% 
      visEvents(selectNode = "function(properties) {alert(this.body.data.nodes._data[properties.nodes[0]].annotation); }")
  })
  output$network_plot <- renderVisNetwork({network_plot()})
}

check_for_pfams <- function(all_pfams, pfam_list, gene_list){
  shapes <- rep(NA,length(gene_list))
  for(ix in 1:length(gene_list)){
    pfams_of_gene <- all_pfams[which(all_pfams$target_name == gene_list[ix]),]$query_accession 
    if(TRUE %in% (pfams_of_gene %in% pfam_list$query_accession)){shapes[ix]="diamond"}
  }
  return(shapes)
}

foldchange_vertices <- function(igraph_network, foldchange, column){
  foldchange[,"NULL"]<-NA # work-around because it returns "incorrect number of dimensions" if single column
  # Get the name of all the vertices in the network
  vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
  # Set default fold-change values to 1 on all vertices
  fc <- rep(0, length(vertices_names))
  foldchange_genes <- rownames(foldchange)
  for(name in vertices_names){
    if(name %in% foldchange_genes){fc[match(name,vertices_names)] <- foldchange[foldchange_genes==name,][,column]}}
  br <- rev(brewer.pal(n = 7, name = 'RdBu')) # brewer colors, reversed
  br[4] <- "gray"
  bn <- .bincode(fc,breaks=c(-Inf,-3,-2,-1,1,2,3,Inf)) # bins of values
  attribute <- br[bn] # assign color to FC values
  return(attribute)
}

category_shape <- function(vertices_names, shapes, categories, column, shape, go_mapping, domain_mapping){
  vertices_names <- tolower(vertices_names)
  matched_gos <- match(vertices_names, tolower(unique(as.vector(go_mapping[go_mapping[,2] %in% categories[,column],][,1]))),nomatch=0)
  matched_domains <-  match(vertices_names, tolower(unique(as.vector(domain_mapping[domain_mapping[,2] %in% categories[,column],][,1]))),nomatch=0)
  matched_categories <- match(vertices_names, tolower(categories[,column]),nomatch=0)
  all_matched <- matched_gos + matched_domains + matched_categories
  all_matched[all_matched==0] <- NA
  all_matched[!is.na(all_matched)] <- shape
  shapes[!is.na(all_matched)] <- all_matched[!is.na(all_matched)]
  return(shapes)
}