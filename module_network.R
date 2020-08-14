networkUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Network",
    sidebarPanel(
      width = 4,
      numericInput(ns("n_rows"), "Number of rows to include:",25),
      numericInput(ns("threshold"), "Select MR Threshold:",10),
      numericInput(ns("font_size"), "Select label size:",15),
      checkboxInput(ns("star_bait"), "Change bait gene to star shape?", value = T, width = NULL),
      checkboxInput(ns("convert_symbols"), "Convert gene names to symbols", value = T, width = NULL),
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

network <- function(input, output, session, coexpression, annotations, 
                      symbols, foldchange, categories, go_mapping, domain_mapping) {
  ns <- session$ns

  # The first gene in the Mutual Rank table will be the reference gene
  reference_gene <- reactive({row.names(coexpression())[1]})
  
  # Shiny reactive function to reduce the amount of rows in the heatmap based on user specification
  nrow_coexpression <- reactive({
    coexpression <- coexpression()
    if(input$n_rows<length(rownames(coexpression))){
      coexpression <- coexpression[1:input$n_rows,1:input$n_rows]}
    return(coexpression)
  })
  
  raw_network <- reactive({
    # Convert the Mutual Rank coexpression table to an adjacency matrix
    adj_matrix <- as.matrix(nrow_coexpression())
    # Set all MR values larger than threshold to 0
    adj_matrix[adj_matrix > as.integer(input$threshold)] <- 0
    # Create the igraph network from the adjecency matrix
    igraph_network <- graph_from_adjacency_matrix(adj_matrix, mode="undirected", weighted = T,
                                                  diag = F, add.colnames = NA, add.rownames = NULL)
    igraph_network <- set_edge_attr(igraph_network,"width",value=1)         # Set the edge width to a constant 1
    igraph_network <- set_edge_attr(igraph_network,"color", value="gray")   # Set edge color attribute to gray
    igraph_network <- set_vertex_attr(igraph_network,"color",value="gray")  # Set vertex color attribute to gray
    return(igraph_network)
  })
  
  # Makes a reactive list of categories and renders as a select input in the UI side panel
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
    # Get a list of the names of the vertices in the network
    vertices_names <- get.data.frame(raw_network(), what= c("vertices"))[,1]
    # Creates a list of shapes starting with the default "dot" shape for each vertex
    shapes <- rep("dot", length(vertices_names))
    # Call the category_shape function to edit the list of shapes, keep dots as default
    if(input$category_diamond!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_diamond, "diamond", go_mapping(),domain_mapping())}
    if(input$category_star!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_star, "star", go_mapping(),domain_mapping())}
    if(input$category_triangle!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_triangle, "triangle", go_mapping(),domain_mapping())}
    if(input$category_downtriangle!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_downtriangle, "triangleDown", go_mapping(),domain_mapping())}
    if(input$category_square!="None"){shapes<-category_shape(vertices_names, shapes, categories(), input$category_square, "square",go_mapping(),domain_mapping())}
    return(shapes)
  })
  
  # Gets the list of columns in the foldchange data.frame
  foldchange_columns <- reactive({colnames(foldchange())[1:length(colnames(foldchange()))]})
  # Creates a UI selectInput object to select one of the columns from the foldchange data
  output$foldchange_columns <- renderUI({selectInput(ns("foldchange_column"), "Choose Foldchange Column:", 
                                                     selected=foldchange_columns()[1], foldchange_columns()) })
  # Change the colors of vertices based on fold-change data in the selected column
  reactive_foldchange <- reactive({
    if(input$quant_vertices){
      return(foldchange_vertices(raw_network(), foldchange(), input$foldchange_column))}  
      else{return("gray")}
  })
  
  # Take all the reactive igraph parameters and apply them to the network attributes
  final_network <- reactive({
    igraph_network <- raw_network()
    # Get a list of the names of the vertices in the network
    vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
    # Get the annotations for all of the vertices and save as vector of the same order
    annot <- as.vector(annotations()[,1][match(tolower(vertices_names), tolower(rownames(annotations())))])
    # Set the vertex annocation attribute in the igraph_network object
    igraph_network <- set_vertex_attr(igraph_network,"annotation",value=(annot))
    # If there are edges in the network remove the weight attribute
    if(length(E(igraph_network))>0){igraph_network <- remove.edge.attribute(igraph_network, "weight")}
    # Change vertex attributes using reactive specified functions
    igraph_network <- set_vertex_attr(igraph_network,"color",value=reactive_foldchange())           
    igraph_network <- set_vertex_attr(igraph_network,"shape",value=reactive_category_shape())
    igraph_network <- set_vertex_attr(igraph_network,"font",value=paste(input$font_size,"px arial black",sep=""))
    # By default set the reference gene as a star
    if(input$star_bait){igraph_network <- set_vertex_attr(igraph_network,"shape",index=V(igraph_network)[1], value="star")}
    # Convert gene names to symbols if specified
    if(input$convert_symbols){
      new_symbols <- symbol_converter(symbols(), vertex_attr(igraph_network,'name'))
      igraph_network <- set_vertex_attr(igraph_network,"name",value=new_symbols)
    }
    return(igraph_network)
  })
  
  # Convert the igraph network instance into a visNetwork
  network_plot <- reactive({
    data <- toVisNetworkData(final_network())
    visNetwork(nodes = data$nodes, edges = data$edges, randomSeed=1) %>% 
      visIgraphLayout() %>%
      # visEvents that calls a pop-up alert window containing the gene annotation when a vertex is pressed
      visEvents(selectNode = "function(properties) {alert(this.body.data.nodes._data[properties.nodes[0]].annotation); }")
  })
  output$network_plot <- renderVisNetwork({network_plot()})
}