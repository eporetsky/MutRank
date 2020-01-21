MRnetworkUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Network",
    sidebarPanel(
      width = 2,
      actionButton(ns("update_network"), "Update Network"),
      checkboxInput(ns("red_bait"), "Color Bait Gene?", value = F, width = NULL),
      checkboxInput(ns("find_clusterONEs"), "Find Clusters?", value = F, width = NULL),
      checkboxInput(ns("quant_vertices"), "Set size of vertices?", value=T, width=NULL),
      #checkboxInput(ns("size_or_color"), "FC by size or color?", value=T, width=NULL),
      numericInput(ns("min_overlap"), "Minimum overlapping genes in cluster:", 1)
      
      #sliderInput(ns("range"),"Range:",min = 0, max = 100,value = c(2,10))
     ),
    mainPanel(
      #tableOutput({ns("df_mr")})
      visNetworkOutput(ns("network_plot"),width="1000px", height="600px")
      #
    )
  )
}

MRnetwork <- function(input, output, session, data) {
  ns <- session$ns
  
  output$df_mr <- renderTable({
    data[1:5,1:5]
  })
  
  slider_values <- reactive({
    input$range
  })
  
  igraph_network <- reactive({
    # For the purpose of testing, I'm using a small premade network, remove later
    # Don't forget to add () to all data() variables, because it's a reactive objecgt
    # igraph needs to work with as.matrix, fyi
    #data <- read.table("mr_adjacency_matrix.csv", header=T, sep="\t")
    #adj_matrix <- as.matrix(data[,-1])
    #rownames(adj_matrix) <- data[,1]
    # Use input coexpression table as an adjacency matrix
    adj_matrix <- as.matrix(data())
    # This is the gene that will be colored red
    selected_gene <- rownames(adj_matrix)[1]
    # transform data using Wisecavers exponential deccay parameters
    adj_matrix <- apply(adj_matrix, 2, mr_exponential_decay )
    # Filter out vertices with higher MR values than specified
    adj_matrix[adj_matrix <= 0.01] <- 0
    # Create the igraph network from the adjecency matrix
    igraph_network <- graph_from_adjacency_matrix(adj_matrix, mode="undirected", weighted = T,
                                                  diag = F, add.colnames = NA, add.rownames = NULL)
    
    # Get the name of all the vertices in the network
    vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
    # Set name of vertices to black
    vertices_colors <- rep("black", length(vertices_names))
    annot <- read.table("annotation/ZmV3.csv", sep="\t", header=T, row.names=1, quote="")
    annot <- as.vector(annot[, "annotation"][match(tolower(rownames(adj_matrix)), tolower(rownames(annot)))])
    igraph_network <- set_vertex_attr(igraph_network,"annotation",value=(annot))
    igraph_network <- remove.edge.attribute(igraph_network, "weight")
    igraph_network <- set_vertex_attr(igraph_network,"color",value=vertices_colors)
    #igraph_network<-clusterONEs(igraph_network,input$min_overlap)
    if(input$find_clusterONEs){igraph_network<-clusterONEs(igraph_network,input$min_overlap)}
    quant_data <- read.table("zm3_fc.csv",sep=",", header=T)
    if(input$quant_vertices){igraph_network<-quantVertices(igraph_network, quant_data)}
    if(input$red_bait){
      vertices_colors <- get.vertex.attribute(igraph_network,"color")
      vertices_colors[match(selected_gene, vertices_names)] <- "red"
      igraph_network <- set_vertex_attr(igraph_network,"color",value=vertices_colors)
      }
    return(igraph_network)
    })
  
  output$network_plot <- renderVisNetwork({network_plot()})
  
  network_plot <- reactive({
    input$update_network
    data <- toVisNetworkData(igraph_network())
    visNetwork(nodes = data$nodes, edges = data$edges) %>% 
       visIgraphLayout() %>%
       #visEdges(smooth = FALSE) %>% 
       visEvents(selectNode = "function(properties) {alert(this.body.data.nodes._data[properties.nodes[0]].annotation); }")
  })
}

# exponential decay function
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


color_clusters <- function(clusters, vertices_names, vertices_colors){
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
    if(TRUE %in% (pfams_of_gene %in% pfam_list$query_accession)){shapes[ix]="star"}
  }
  return(shapes)
}

quantVertices <- function(igraph_network, data){
  # Get the name of all the vertices in the network
  vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
  # Set default fold-change values to 1 on all vertices
  fc <- rep(0, length(vertices_names))
  for(name in vertices_names){
    if(name %in% data[,"GeneID"]){fc[match(name,vertices_names)] <- data[data[,"GeneID"]==name,][,"ZmPep3"]}}
  #if(size){attribute <- fc**3+10
  #  igraph_network <- set_vertex_attr(igraph_network,"value",value=attribute)
  #} else{#display.brewer.pal(n = 7, name = 'RdBu')
  br <- rev(brewer.pal(n = 7, name = 'RdBu')) # brewer colors, reversed
  bn <- .bincode(fc,breaks=c(-Inf,-3,-2,-1,1,2,3,Inf)) # bins of values
  attribute <- br[bn] # assign color to FC values
  igraph_network <- set_vertex_attr(igraph_network,"color",value=attribute)
  return(igraph_network)
}

clusterONEs <- function(igraph_network,min_overlap){
  # Get the name of all the vertices in the network
  vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
  # Set name of vertices to black
  vertices_colors <- rep("black", length(vertices_names))
  # Instead of writing the actual igraph network table to file, capture the output to variable
  clusterONE_IO_edges <- capture.output(write.table(get.data.frame(igraph_network, what= c("edges")), row.names=FALSE,col.names=F,quote=F,sep="\t"))
  # Run the ClusterONE algorithm on the edges tsv file that's captured as a variable and save to variable with intern=TRUE
  clusterONE_IO_clusters <- system("java -jar cluster_one-1.0.jar /dev/stdin -F 'csv'", intern=TRUE, input=clusterONE_IO_edges)
  # Read the ClusterONE output from a temporary file 
  clusterONEs <- read.table(text=clusterONE_IO_clusters, header=T, sep=",")
  # Filter the predicted clusters based on p-value<0.1, as used in Wisecave paper
  clusters <- strsplit(as.vector(clusterONEs[clusterONEs[,"P.value"]<0.1,]$Members), " ")
  # If significant clusters were detected, color them
  if(length(clusters)>1){
    # If more than one cluster was detected, check for overlapping nodes and combine them
    clusters <- union_overlapping_clusters(clusters,min_overlap)
    vertices_colors <- color_clusters(clusters, vertices_names, vertices_colors)
  } else{if(length(clusters)==1){
    # If only one cluster was detected, just color it
    vertices_colors <- color_clusters(clusters, vertices_names, vertices_colors)
  }}
  if(TRUE & length(clusters)>0){
    for(cluster in clusters){
      # get all pair-wise combinations of vertices, unlist to fit next function
      vertex_pairs <- unlist(combn(cluster,2,simplify=F))
      # get all existing edges between pairs of vertices (equal zero=no edge)
      existing_edges <- get.edge.ids(igraph_network, vertex_pairs)
      # get the actual edges from the graph based on existing edge pairs
      modify_edges <- E(igraph_network, P=vertex_pairs[rep(existing_edges!=0,each=2)])
      print(modify_edges)
      modify_color <- vertices_colors[match(cluster[[1]], vertices_names)]
      igraph_network <- set_edge_attr(igraph_network,"color",modify_edges,modify_color)
      igraph_network <- set_edge_attr(igraph_network,"width",modify_edges,5)
    }
  } else{igraph_network <- set_vertex_attr(igraph_network,"color",value=vertices_colors)}
  
  proteins_pfams <- read.table("domains/Zea_mays.AGPv3.21.pep.longest.parsed.split.tsv", header=T, sep="\t")
  selected_pfams <- read.table("domains/wisecaver_pfams_sm_split.tsv", header=T, sep="\t")
  shapes <- check_for_pfams(proteins_pfams,selected_pfams,V(igraph_network)$name)
  igraph_network <- set_vertex_attr(igraph_network,"shape",value=shapes)
  return(igraph_network)
}

# it's annotated so gotta remove the last column
#adj_matrix <- adj_matrix[,1:(dim(adj_matrix)[2]-1)]
# I think this is an old implementation for filtering self-vertices. Not necessary but nice
#for(col in names(adj_matrix)) set(adj_matrix, i=which(adj_matrix[[col]]<=1), j=col, value=NA)
# I used this incase the matrix is not symmetrical, but no longer relevant
#adj_matrix <- adj_matrix[order(row.names(adj_matrix)),]
#adj_matrix <- adj_matrix[,order(colnames(adj_matrix))]
# I don't think I need this anymore. Only relevant if network was edited in previous steps
#igraph_network <- graph_from_adjacency_matrix(adj_matrix_original, mode="undirected", weighted = T,
#                                              diag = F, add.colnames = NA, add.rownames = NULL)
#igraph_network <- delete.edges(igraph_network, which(E(igraph_network)$weight>10))
#igraph_network <- delete.edges(igraph_network, which(E(igraph_network)$weight<2))
#igraph_network <- delete.edges(igraph_network, which(E(igraph_network)$weight<slider_values()[1]))
#igraph_network <- delete.edges(igraph_network, which(E(igraph_network)$weight>slider_values()[2]))

