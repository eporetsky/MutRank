if("test"=="te"){print("Test")}
# Returns a Mutual Rank table based on user specified parameters
coexpression_table <- function(datm, gene_list, reference_method, num_top_pcc, compound_method=NA, order_ref_list){
  datm <- as.data.frame(datm)
  
  # keep order_ref_list TRUE unless user is using a list of reference genes
  if(is.null(order_ref_list)){order_ref_list<-T}
  if(reference_method!="Reference gene list"){order_ref_list<-T}
  
  # Handles the Compound reference gene option according to selected method
  if(reference_method=="Compound reference gene"){
    if(compound_method=="Sum"){compound_gene <- colSums(datm[gene_list,])}
    if(compound_method=="Average"){compound_gene <- colMeans(datm[gene_list,])}
    if(compound_method=="Min"){compound_gene <- apply(datm[gene_list,],2,min)}
    if(compound_method=="Max"){compound_gene <- apply(datm[gene_list,],2,max)}
    datm <- rbind(datm,"Compound_Gene"=compound_gene) # Adds the compound gene as a new row
    gene_list <- "Compound_Gene"  # Replaces the reference gene to Compound_Gene
  }
  
  # If gene_list==1 it first finds top coexpressed genes based on PCC values
  if(length(gene_list)==1){
    gene_cors <- cor(t(datm[gene_list,]), t(datm))          # Calculate PCC between the reference gene with all other genes
    gene_cors_ordered <- order(gene_cors, decreasing=T)     # Reorder the correlation values so highest one is first
    genes_for_mr <- datm[gene_cors_ordered[1:num_top_pcc],] # Select the top n correlating genes to calculate MR on
  } else{
    genes_for_mr <- datm[gene_list,]
  }
  
  cor_for_mr <- cor(t(genes_for_mr),t(datm))            # Calculate PCC between genes_for_mr vector and complete data
  rank_for_mr <- apply(cor_for_mr, 1, frankv, order=-1) # Fast rank all the columns that contain the selected genes
  row.names(rank_for_mr) <- row.names(datm)             # Add back row names since it is lost in the 'cor' function
  rank_for_mr <- rank_for_mr[row.names(genes_for_mr),]  # Remove all rows except for the genes_for_mr vector
  mr <- sqrt(rank_for_mr*t(rank_for_mr))                # Calculate all Mutual Rank values between selected genes
  mr <- as.data.frame(mr)                               # Column reordering doesn't work on matrix so convert to data.frame
  if(order_ref_list){mr<-mr[order(mr[,1]),]}            # Order rows from lowest MR values to highest based on first column
  mr <- mr[, row.names(mr)]                             # Reorder columns using the order of row names for symmetrical table
  return(mr)
}

# Returns a list of gene names with some converted to symbols when possible
symbol_converter <- function(symbols,gene_names){
  symbols <- as.vector(symbols[,1][match(tolower(gene_names), tolower(rownames(symbols)))])
  no_symbols <- attr(na.omit(symbols),"na.action")
  symbols[no_symbols] <- gene_names[no_symbols]
  return(symbols)
}

# Add an annotation column to the Mutual Rank table
df_annotator <- function(df, annotations){
  annotations <- as.vector(annotations[,1][match(tolower(rownames(df)), tolower(rownames(annotations)))])
  return(cbind(df,annotations=annotations))
}
# Add a symbol coumn to the Mutural Rank table
df_add_symbols <- function(df, symbols){
  symbols <- as.vector(symbols[,1][match(tolower(rownames(df)), tolower(rownames(symbols)))])
  return(cbind(df,symbols=symbols))
}

# Add a fold-change column to the Mutual Rank table
df_add_foldchange <- function(df, foldchange){
  FC <- foldchange[match(tolower(rownames(df)), tolower(rownames(foldchange))),]
  return(cbind(df,FC))
}

# Add a column for each category in the MR table
df_add_categories <- function(df, categories, go_mapping, domain_mapping){
  genes <- tolower(rownames(df))
  for(category in colnames(categories)){
    matched_gos <- match(genes, tolower(unique(as.vector(go_mapping[go_mapping[,2] %in% categories[,category],][,1]))),nomatch=0)
    matched_domains <-  match(genes, tolower(unique(as.vector(domain_mapping[domain_mapping[,2] %in% categories[,category],][,1]))),nomatch=0)
    matched_categories <- match(genes, tolower(categories[,category]),nomatch=0)
    all_matched <- matched_gos + matched_domains + matched_categories
    all_matched[all_matched==0] <- NA
    all_matched[!is.na(all_matched)] <- "Y"
    df <- cbind(df,all_matched)
    colnames(df)[length(colnames(df))]<-category
  }
  return(df)
}

# Edit the MR table output based on Shiny inputs
df_output_editor <- function(df,firstColumn,round){
  if(firstColumn){df<-df[,1, drop=FALSE]}  # Filter all columns except the first one
  if(round){df<-round(df, digits=0)}                 # Round to lower integer using floor
  return(df)
}

# https://stackoverflow.com/questions/13773770/split-comma-separated-strings-in-a-column-into-separate-rows
wide2long <- function(df){
  dt <- setDT(df)[, lapply(.SD, function(x) unlist(tstrsplit(x, ",", fixed=TRUE))), 
                  by = genes][!is.na(names(df)[2])]
  return(setDF(dt)) # convert data.table back to data.frame
}

# The main function for drawing the Mutual Rank heatmap using ggplot
draw_heatmap <- function(melted_cormat, text_threshold, text_size){
  num <- length(rownames(melted_cormat))
  # Check if there is atleast one value>100, otherwise function would return an error
  # Change all values higher than 100 to 101 so they appear white instead of grey
  if(sum(melted_cormat$value>100)>0){melted_cormat[melted_cormat$value>100,]$value=100}
  # If text_threshold is larger than 100 change it to 100 so larger values don't show
  if(text_threshold>100){text_threshold<-100}
  # Create the heatmap using ggplot
  ggplot(data = melted_cormat, aes(Var1,Var2,fill=value))+
    geom_tile(color = "black",size=1,aes(fill = value))+
    geom_text(data=subset(melted_cormat,(melted_cormat$value>0 & melted_cormat$value<text_threshold)),
              aes(label=floor(value)), size=log(text_size**2), fontface="bold", colour="black")+
    scale_fill_gradient2(low = "red", high = "white",# mid = "white", 
                         limit=c(0,100), space = "Lab", midpoint = 95, 
                         breaks=c(0,20,40,60,80,100),
                         labels=c("0","20","40","60","80",">100"),
                         name="Mutual rank scores", guide="colourbar"
    ) +
    
    theme(# Hide panel borders and remove grid lines
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Remove axis titles
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      # Change axis line
      axis.line = element_line(colour = "white"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = text_size, colour="black",hjust=1,face="bold"),
      axis.text.y = element_text(angle = 0, vjust = 0.5, size = text_size, colour="black",face="bold"),
      axis.ticks = element_blank(),
      # Edit the legend
      legend.title = element_text(size=13,face="bold"),
      legend.text = element_text(size=13,face="bold"),
      legend.justification = c(1, 0),
      legend.position = c(1, 0.7),
      legend.direction = "horizontal",
    ) +
    guides(fill=guide_colorbar(barwidth=15,barheight=2,title.position="top",title.hjust=0.5,
                               frame.colour = "black", frame.linewidth=2, ticks = F))+
    coord_fixed() +
    # Reverse the x-axis
    scale_x_discrete(limits = rev(levels(melted_cormat$Var1))) 
}

# Get upper triangle of the correlation matrix (convert the lower triangle to NAs)
# Based on https://stackoverflow.com/questions/59917970/how-do-i-create-a-heat-map-in-r
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Function to remove duplicated rows in the expression data
data_loader <- function(temp_df){
  duplicated <- duplicated(temp_df[,1])
  deduped <- temp_df[!duplicated,-1]
  row.names(deduped) <- temp_df[!duplicated,1]
  return(deduped)
}

# Returns a list of colors for each vertex in the network based on fold-change data
foldchange_vertices <- function(igraph_network, foldchange, column){
  # Work-around because it returns "incorrect number of dimensions" if single column
  foldchange[,"NULL"]<-NA 
  # Get the name of all the vertices in the network
  vertices_names <- get.data.frame(igraph_network, what= c("vertices"))[,1]
  # Set default log fold-change values to 0 on all vertices
  fc <- rep(0, length(vertices_names))
  foldchange_genes <- rownames(foldchange)
  # For every gene in the network, check if it's in the FC data and add it to fc base on index
  for(name in vertices_names){
    if(name %in% foldchange_genes){
      fc[match(name,vertices_names)] <- foldchange[foldchange_genes==name,][,column]
    }
  }
  # Used reversed brewer colors to color the vertices based on fold change
  br <- rev(brewer.pal(n = 7, name = 'RdBu')) 
  # Set the middle color, between -1 and 1, to gray
  br[4] <- "gray"
  # Set FC value bins manually
  bn <- .bincode(fc,breaks=c(-Inf,-3,-2,-1,1,2,3,Inf)) 
  # Assign color to FC values based on appropriate bins
  attribute <- br[bn] 
  return(attribute)
}

# Returns a list of shapes for each vertex in the network based on given categories
category_shape <- function(vertices_names, shapes, categories, column, shape, go_mapping, domain_mapping){
  ### Note: This function will overwrite previously set "shapes" and it depends on the order this function is called ###
  vertices_names <- tolower(vertices_names)
  # To find if a vertex is associated with a given category use match. If there is a match it returns a number, else returns 0
  matched_gos <- match(vertices_names, tolower(unique(as.vector(go_mapping[go_mapping[,2] %in% categories[,column],][,1]))),nomatch=0)
  matched_domains <-  match(vertices_names, tolower(unique(as.vector(domain_mapping[domain_mapping[,2] %in% categories[,column],][,1]))),nomatch=0)
  matched_categories <- match(vertices_names, tolower(categories[,column]),nomatch=0)
  # Add the three vectors checking the specified category with all 3 categories
  all_matched <- matched_gos + matched_domains + matched_categories
  # If summed value is 0 it means that the gene vertex didn't match any of the categories
  all_matched[all_matched==0] <- NA
  # Set all the non-NA values to the specified cells
  all_matched[!is.na(all_matched)] <- shape
  # Change the the shapes that agreed with the currently specified category
  shapes[!is.na(all_matched)] <- all_matched[!is.na(all_matched)]
  return(shapes)
}

