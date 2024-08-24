# Load the libraries
library(igraph) # network analysis library
library(tidyverse)
library(dplyr)
library(knitr)


## FUNCTION_3 -> Take the name of the TFs and return the GPCRs associated with it
# Define the function
find_the_gpcrs <- function(
    tf_gene,
    matrix,
    gpcrs_gene,
    threshold = 0.05
) {
  
 # Convert the tf_gene into a vector so it can allow single and multiple selection of TFs
  if(!is.vector(tf_gene)) {
    tf_gene <- as.vector(tf_gene)
  }
  # Check if the genes are in the matrix
  if(!all(tf_gene %in% matrix$Target)) {
    stop("TFs gene not found")
  }
  
  tf_gene_row <- matrix %>%
    filter(Target %in% tf_gene,
           Correlation %in% gpcrs_gene,
           !is.na(r),
           r < threshold) %>%
    select(TFs = Target, GPCRs = Correlation)
  
  return(tf_gene_row)
  
}
  
## Example
tf_gene <- "HDX"

associated_gpcrs <- find_the_gpcrs(
  tf_gene,
  cor_results_gpcrs_gp_tfs,
  gpcrs_gene_with_significant_impact
)

summary_gpcrs_table <- find_the_gpcrs(
  tf_gene,
  cor_results_gpcrs_gp_tfs,
  gpcrs_gene_with_significant_impact
)






