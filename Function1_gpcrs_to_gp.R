library(igraph) # network analysis library
library(tidyverse)
library(dplyr)
library(data.table)
## Create a function that takes the name of the GPCRs and return the small G protein 
# Define the funtion
find_the_g_protein <- function(
    gpcrs_gene,
    matrix,
    G_protein,
    threshold= 0.05
) {

# Convert the gpcrs_gene into a vector so it can allow single and multiple selection of GPCRs
  if (!is.vector(gpcrs_gene)) {
    gpcrs_gene <- as.vector(gpcrs_gene)
  }
  
# Check if the genes are in the matrix
  if (!all(gpcrs_gene %in% matrix$Target)) {
    stop("GPCRs gene not found")
  }
  
  gpcrs_gene_row <- matrix %>%
    filter(Target %in% gpcrs_gene,
           Correlation %in% G_protein,
           !is.na(r),
           r < threshold)
  
  #[matrix$Target == gpcrs_gene, ] # [ ] this is for indexing - write before the 
  #comma if you want to select rows - and after the comma if you want to select columns

  # Find associated small G proteins
  associated_gprotein <- gpcrs_gene_row$Correlation[!is.na(gpcrs_gene_row$r) & gpcrs_gene_row$r < threshold & gpcrs_gene_row$Correlation %in% G_protein]
  
  # Create a data frame edges
  edges <- gpcrs_gene_row %>%    # if only single parameter is passed -> data.frame(from = gpcrs_gene, to = associated_gprotein )
    select(from = Target, to = Correlation)
   
  
  # Create the network graph
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Assign colour
  V(g)$color <- ifelse(V(g)$name %in% gpcrs_gene, "red", "lightblue")
  # Plot the graph
  plot(g, vertex.label.color = "black", vertex.label.cex = 0.8, vertex.size = 30, edge.arrow.size = 5.0)
  
  return(g)
}

#Example

gpcrs_gene <- c("OR56B4", "OR10AG1")



associated_gprotein <- find_the_g_protein(
  gpcrs_gene,
  cor_results_gpcrs_gp_tfs,
  G_protein_with_significant_impact
)

print(associated_gprotein)
