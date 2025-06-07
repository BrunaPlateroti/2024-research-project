library(igraph) # network analysis library
library(dplyr)

### Create network for visualisation
# Filter the data for r < 0.05
cor_results_htr_gp_tf_05 <- cor_results_htr_gp_tf[r<0.05]
# Filter further the data for q < 0.05
cor_results_htr_gp_tf_q05 <- cor_results_htr_gp_tf[q<0.05]
# Filter further the data for q < 0.05
cor_results_htr_gp_tf_q01 <- cor_results_htr_gp_tf[q<0.01]

## You can also create a subgroup of the data to display only specific target with its correlation
filtered_results_htr_gp_tf <- cor_results_htr_gp_tf_q05 %>%
  filter(Target == "HTR4" | Correlation == "HTR4" )

# Create a data frame edges to include all the data
edges <- filtered_results_htr_gp_tf %>%    
  select(from = Target, to = Correlation)

# Create the graph
g <- graph_from_data_frame(edges, directed = FALSE)

# Find all vertices directly connected to HTR4
#connected_vertices <- ego(g, order = 1, nodes = "HTR4", mode = "all")[[1]]

## If you want to include only edges with higher degree
#vertex_degrees <- degree(g) # Calculate degree of each vertex

#degree_threshold <- 6 # Adjust as needed

# Filter based on the degree threshold
#vertices_to_keep <- names(vertex_degrees[vertex_degrees > degree_threshold]) 

# Create a subgraph with only the verteces that meet the criteria
#subgraph <- induced_subgraph(g, vids = connected_vertices)

# Use this one for subgraph
#V(subgraph)$color <- ifelse(V(subgraph)$name == "HTR4", "red", "lightblue")

# Assign colour
#V(g)$color <- ifelse(degree(g) > 2, "red", "lightblue") # choose the vertices degree
# depending on the number of edges connected to it and how many you want to highlight.

## The colour assignment must be different if specific target are selected
V(g)$color <- ifelse(V(g)$name == "HTR4", "red", "lightblue")
layout <- layout_with_fr(g, niter = 500)

# Plot the graph
plot(g,
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.size = 30, 
     edge.arrow.size = 5.0)



# Isolate metadata from the DepMap_RNAi original data
CCLE_meta <- DepMap_RNAi$CCLE_Name

# Identify tumour lineage
tumour_lineages <- DepMap_RNAi %>%
  mutate(tumour_lineages = sub("^[^_]*_", "", CCLE_Name)) %>%
  # Extract after the first underscore
  distinct(tumour_lineages) %>%
  # Extract as a vector
  pull(tumour_lineages)
  
  
  




