# Load libraries
library(igraph) # network analysis library
library(tidyverse)
library(dplyr)
library(knitr)


## FUNCTION_2 -> Take the name of the GPCRs and return the TFs associated with it
# Define the function

find_the_tf <- function(
    gpcrs_gene,
    matrix,
    tf_gene,
    threshold = 0.05
) {
  
  # Convert the gpcrs_gene into a vector so it can allow single and multiple selection of GPCRs
  if(!is.vector(gpcrs_gene)) {
    gpcrs_gene <- as.vector(gpcrs_gene)
  }  
  
  # Check if the genes are in the matrix
  if(!all(gpcrs_gene %in% matrix$Target)) {
    stop("GPCRs gene not found ")
  }
  
  gpcrs_gene_row <- matrix %>%
    filter(Target %in% gpcrs_gene,
           Correlation %in% tf_gene,
           !is.na(r),
           r < threshold) %>%
    select(GPCRs = Target, TFs = Correlation)
  
  return(gpcrs_gene_row)
  
}


## Example
# Select single or multiple GPCRs gene of interest
gpcrs_gene <- c("ADGRG4", "GRM3")

# Call the function
associated_tf <- find_the_tf(
  gpcrs_gene,
  cor_results_gpcrs_gp_tfs,
  tfs_gene_with_significant_impact
)

# Visualize the summary table 
summary_tf_table <- find_the_tf(
  gpcrs_gene,
  cor_results_gpcrs_gp_tfs,
  tfs_gene_with_significant_impact
)

# Optional - to investigate                     
library(ggplot2)
ggplot(summary_table, aes(x = GPCR, y = Count, fill = GPCR)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Associated Small G Proteins",
       x = "GPCR Gene", y = "Count") +
  theme_minimal()
