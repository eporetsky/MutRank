##############
############## Functions used in module_mutualrank ##############
coexpression_table <- function(datm, reference_gene, reference_method, num_top_pcc, compound_method=NA){
  #### Working Stand-Alone Version For Calculating MR of Top Correlating Genes ####
  #### Written by Elly Poretsky on 12.27.19, for a more efficient mutRank      ####
  if(is.null(reference_gene)){return(data.frame())}
  datm <- as.data.frame(datm)
  datm <- datm[rowSums(datm[, -1])>0, ] # removes zero rows to avoid warning 
  gene_list <- unlist(strsplit(reference_gene, "[[:space:],]+")) # Space regex characters: tab, newline, vertical tab, form feed, carriage return, space
  
  if(reference_method=="Compound reference gene"){
    if(compound_method=="Sum"){compound_gene <- colSums(datm[gene_list,])}
    if(compound_method=="Average"){compound_gene <- colMeans(datm[gene_list,])}
    if(compound_method=="Min"){compound_gene <- apply(datm[gene_list,],2,min)}
    if(compound_method=="Max"){compound_gene <- apply(datm[gene_list,],2,max)}
    datm <- rbind(datm,"Compound_Gene"=compound_gene)
    gene_list <- "Compound_Gene"
  }
  
  # If gene_list==1 it first finds top coexpressed genes based on PCC and
  if(length(gene_list)==1){
    gene_cors <- cor(t(datm[gene_list,]), t(datm))^2 # Calculate the correlation of the gene with all other genes
    gene_cors_ordered <- order(gene_cors, decreasing=T) # Reorder the correlation values so highest one is first
    genes_for_mr <- datm[gene_cors_ordered[1:num_top_pcc],] # Select the top n correlating genes to calculate MR on
  } else{
    genes_for_mr <- datm[gene_list,]
  }
  
  # Used to automatically support tolower, worth returning: match(tolower(gene_list), tolower(row.names(datm)))
  cor_for_mr <- cor(t(genes_for_mr),t(datm)) # Calculate the correlation between selected genes and complete data
  rank_for_mr <- apply(cor_for_mr, 1, frankv, order=-1) # Fast rank all the columns that contain the selected genes
  row.names(rank_for_mr) <- row.names(datm) # Add row names since it is lost in the 'cor' function
  rank_for_mr <- rank_for_mr[row.names(genes_for_mr),] # Remove all the rows except for the selected genes
  mr <- sqrt(rank_for_mr*t(rank_for_mr)) # Calculate the MR values between selected genes
  mr <- as.data.frame(mr)  # column reordering doesn't work on matrix
  return(mr)
}

order_coexpression_table <- function(mr){
  mr <- mr[order(mr[,1]),] # order rows from lowest MR values
  mr <- mr[rownames(mr)]   # reorder column with rowname orde
  return(mr)
}

symbol_converter <- function(symbols,gene_names){
  symbols <- as.vector(symbols[,1][match(tolower(gene_names), tolower(rownames(symbols)))])
  no_symbols <- attr(na.omit(symbols),"na.action")
  symbols[no_symbols] <- gene_names[no_symbols]
  return(symbols)
}

# Add an annotation column to the MR table
df_annotator <- function(df, annotations){
  annotations <- as.vector(annotations[,1][match(tolower(rownames(df)), tolower(rownames(annotations)))])
  return(cbind(df,annotations=annotations))
}
# Add a symbol coumn to the MR table
df_add_symbols <- function(df, symbols){
  symbols <- as.vector(symbols[,1][match(tolower(rownames(df)), tolower(rownames(symbols)))])
  return(cbind(df,symbols=symbols))
}

df_add_foldchange <- function(df, foldchange){
  # if foldchange contains a single column the cbind changes its name to the variable name...
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
  if(firstColumn){df<-df[,1, drop=FALSE]} 
  if(round){df<-round(df, 0)}
  return(df)
}

# https://stackoverflow.com/questions/13773770/split-comma-separated-strings-in-a-column-into-separate-rows
# if 
wide2long <- function(df){
  dt <- setDT(df)[, lapply(.SD, function(x) unlist(tstrsplit(x, ",", fixed=TRUE))), 
                         by = genes][!is.na(names(df)[2])]
  return(setDF(dt)) # convert data.table back to data.frame
}



