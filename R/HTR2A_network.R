## You can also create a subgroup of the data to display only specific target with its correlation
filtered_results_htr_gp_tf_HTR2A <- cor_results_htr_gp_tf_q05 %>%
  filter(Target == "HTR2A" | Correlation == "HTR2A" )

# Create a data frame edges to include all the data
edges <- filtered_results_htr_gp_tf_HTR2A %>%    
  select(from = Target, to = Correlation)

# Create the graph
g <- graph_from_data_frame(edges, directed = FALSE)

# Find all vertices directly connected to HTR4
#connected_vertices <- ego(g, order = 1, nodes = "HTR2A", mode = "all")[[1]]

## If you want to include only edges with higher degree
#vertex_degrees <- degree(g) # Calculate degree of each vertex

#degree_threshold <- 6 # Adjust as needed

# Filter based on the degree threshold
#vertices_to_keep <- names(vertex_degrees[vertex_degrees > degree_threshold]) 

# Create a subgraph with only the verteces that meet the criteria
#subgraph <- induced_subgraph(g, vids = connected_vertices)

# Use this one for subgraph
#V(subgraph)$color <- ifelse(V(subgraph)$name == "HTR2A", "red", "lightblue")

# Assign colour
#V(g)$color <- ifelse(degree(g) > 2, "red", "lightblue") # choose the vertices degree
# depending on the number of edges connected to it and how many you want to highlight.

## The colour assignment must be different if specific target are selected
V(g)$color <- ifelse(V(g)$name == "HTR2A", "red", "lightblue")
layout <- layout_with_fr(g, niter = 500)

# Plot the graph
plot(g,
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.size = 30, 
     edge.arrow.size = 5.0)
