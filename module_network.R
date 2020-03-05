MRnetworkUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Network",
    sidebarPanel(
      width = 4,
      actionButton(ns("update_network"), "Update Network"),
      checkboxInput(ns("mr_or_ed"), "Exponential decay instead of MR threshold:", value = F, width = NULL),
      uiOutput(ns("threshold")),
      checkboxInput(ns("star_bait"), "Change bait gene to star shape?", value = T, width = NULL),
      checkboxInput(ns("draw_clusters"), "Find Clusters?", value = F, width = NULL),
      checkboxInput(ns("quant_vertices"), "Set size of vertices?", value=T, width=NULL),
      #checkboxInput(ns("size_or_color"), "FC by size or color?", value=T, width=NULL),
      numericInput(ns("min_overlap"), "Minimum overlapping genes in cluster:", 1),
      uiOutput(ns("foldchange_columns")),
      uiOutput(ns("select_clusters")),
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
                      symbols, foldchange, association, categories, go_mapping, domain_mapping) {
  ns <- session$ns

  # This is the gene that will be colored red
  reference_gene <- reactive({row.names(coexpression())[1]})
  # To make a network a 
  output$threshold <- renderUI({
    if(input$mr_or_ed){selectInput(ns("threshold"), "Exponential decay threshold:", selected="all", choices = c(1,2,3,4,5))}
    else{numericInput(ns("threshold"), "MR Threshold:",10)}
  })
  
  raw_network <- reactive({
    # Use input coexpression table as an adjacency matrix
    adj_matrix <- as.matrix(coexpression())
    if(input$mr_or_ed){
      # transform data using Wisecavers exponential deccay parameters
      adj_matrix <- apply(adj_matrix, 2, mr_exponential_decay)
      # Filter out vertices with higher MR values than specified
      adj_matrix[adj_matrix < 0.01] <- 0
    }
    else{
      adj_matrix[adj_matrix > as.integer(input$threshold)] <- Inf
      adj_matrix <- apply(adj_matrix, 2, mr_exponential_decay)
    }
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
  
  # Apply ClusterONE to the network and change edge colors based on detected clusters
  reactive_clusterONEs <- reactive({clusterONEs(raw_network())})
  select_clusters <- reactive({append("all",1:length(reactive_clusterONEs()))})
  output$select_clusters <- renderUI({selectInput(ns("select_clusters"), "Choose Cluster:", selected="all", select_clusters()) })
  reactive_cluster_final <- reactive({
    if(input$draw_clusters & !is.na(reactive_clusterONEs())){
      final_clusters(raw_network(),reactive_clusterONEs(),input$select_clusters,input$min_overlap)
    } else{return(list(c(list(E(raw_network())),"gray")))}
  })
  
  # Change the colors of vertices based on log2 fold-change data. Select column to use
  foldchange_columns <- reactive({colnames(foldchange())[2:length(colnames(foldchange()))]})
  output$foldchange_columns <- renderUI({selectInput(ns("foldchange_column"), "Choose Foldchange Column:", 
                                                     selected=foldchange_columns()[1], foldchange_columns()) })
  reactive_foldchange <- reactive({
    if(input$quant_vertices){return(foldchange_vertices(raw_network(), foldchange, input$foldchange_column))}  
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
    for(i in 1:(length(reactive_cluster_final()))){
      igraph_network <- set_edge_attr(igraph_network,"color",reactive_cluster_final()[[i]][[1]],reactive_cluster_final()[[i]][[2]])
      igraph_network <- set_edge_attr(igraph_network,"width",reactive_cluster_final()[[i]][[1]], value=5)
    }
    new_symbols <- symbol_converter(symbols(), vertex_attr(igraph_network,'name'))
    igraph_network <- set_vertex_attr(igraph_network,"name",value=new_symbols)           # change vertex color based on fc
    return(igraph_network)
  })
  
  # Convert the igraph network instance into a visNetwork
  network_plot <- reactive({
    input$update_network
    data <- toVisNetworkData(final_network())
    visNetwork(nodes = data$nodes, edges = data$edges, randomSeed=1) %>% 
      visIgraphLayout() %>%
      #visEdges(smooth = FALSE) %>% 
      visEvents(selectNode = "function(properties) {alert(this.body.data.nodes._data[properties.nodes[0]].annotation); }")
  })
  
  output$network_plot <- renderVisNetwork({network_plot()})

  
}

mr_exponential_decay <- function(mr){exp(-(mr-1)/5)}

union_overlapping_clusters <- function(clusters,min_overlap){
  changes_were_made = TRUE
  while(changes_were_made){
    changes_were_made = FALSE
    if(length(clusters)==1){break}
    for(cmb in combn(1:length(clusters),2,simplify=F)){
      if(length(intersect(clusters[[cmb[1]]], clusters[[cmb[2]]]))>=min_overlap){
        union_cluster <- list(union(clusters[[cmb[1]]], clusters[[cmb[2]]]))
        clusters <- clusters[-cmb]
        clusters <- append(clusters, union_cluster)
        changes_were_made = TRUE
        break # Break cause list order just changed
      }
    }
  }
  return(clusters)
}

color_clusters <- function(clusters, vertices_names){
    vertices_colors <- rep("", length(vertices_names))
    # color pellete and color selection assuming non-overlapping clusters
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    for(i in 1:length(clusters)){
      for(c in 1:length(vertices_names)){
        if(vertices_names[c] %in% clusters[[i]]){
          vertices_colors[c] <- col_vector[i]
        }
      }
    }
  return(vertices_colors)
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
  foldchange <- foldchange()
  # Get the name of all the vertices in the network
  vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
  # Set default fold-change values to 1 on all vertices
  fc <- rep(0, length(vertices_names))
  gene_col <- colnames(foldchange[1])
  for(name in vertices_names){
    if(name %in% foldchange[,gene_col]){fc[match(name,vertices_names)] <- foldchange[foldchange[,gene_col]==name,][,column]}}
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

clusterONEs <- function(igraph_network){
  # Get the name of all the vertices in the network
  vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
  # Set name of vertices to black, ended up not using now for the time being
  #vertices_colors <- rep("black", length(vertices_names))
  if(length(E(igraph_network))==0){return(NA)}
  # Instead of writing the actual igraph network table to file, capture the output to variable
  clusterONE_IO_edges <- capture.output(write.table(get.data.frame(igraph_network, what= c("edges")), row.names=FALSE,col.names=F,quote=F,sep="\t"))
  # Run the ClusterONE algorithm on the edges tsv file that's captured as a variable and save to variable with intern=TRUE
  clusterONE_IO_clusters <- system("java -jar cluster_one-1.0.jar /dev/stdin -F 'csv'", intern=TRUE, input=clusterONE_IO_edges)
  # Read the ClusterONE output from a temporary file 
  clusterONEs <- read.table(text=clusterONE_IO_clusters, header=T, sep=",")
  # Filter the predicted clusters based on p-value<0.1, as used in Wisecave paper
  clusters <- strsplit(as.vector(clusterONEs[clusterONEs[,"P.value"]<0.1,]$Members), " ")
  # If significant clusters were detected, color them
  return(clusters)}


final_clusters <- function(igraph_network,clusters,selected_cluster,min_overlap){
  vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
  if(selected_cluster!="all"){
    clusters <- clusters[as.integer(selected_cluster)]
  }
  if(length(clusters)>1){
    # If more than one cluster was detected, check for overlapping nodes and combine them
    clusters <- union_overlapping_clusters(clusters,min_overlap)
    vertices_colors <- color_clusters(clusters, vertices_names)
  } else{vertices_colors <- color_clusters(clusters, vertices_names)}
  return_list <- list()
  for(cluster in clusters){
    # get all pair-wise combinations of vertices, unlist to fit next function
    vertex_pairs <- unlist(combn(cluster,2,simplify=F))
    # get all existing edges between pairs of vertices (equal zero=no edge)
    existing_edges <- get.edge.ids(igraph_network, vertex_pairs)
    # get the actual edges from the graph based on existing edge pairs
    return_list <- append(return_list, list(c(
      list(E(igraph_network, P=vertex_pairs[rep(existing_edges!=0,each=2)])),
      list(vertices_colors[match(cluster[[1]], vertices_names)]))))
  }
  return(return_list)
}