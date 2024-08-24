### Analysis of HTR receptor genes in details
## From the serotonin receptor family the following genes are choose, HTR2A, HTR2B,
## HTR1A, HTR1B, HTR1D, HTR4

# Load packages
require(magrittr)
require(tidyr)
require(purrr)
require(collapse)
require(data.table)
library(cordial)

# Choose from the HTR receptor the key genes to analyse in depth
subset_htr <- c("HTR2A", "HTR2B", "HTR1A", "HTR1B", "HTR1D", "HTR4")

# Select genes that contain data from the DepMap_RNAi table
htr_with_data <- intersect(colnames(DepMap_RNAi), subset_htr)

# Create a table with genes from the htr_with_data
data_rnai_htr<- DepMap_RNAi[,..htr_with_data]

#HTR genes important for the survival of some cells
data_rnai_htr[is.na(data_rnai_htr)] <- 0

# cells that are kill/grow slow when the gene is deleted
data_rnai_impact_htr <- data.frame(below.minus.0.5= colSums(data_rnai_htr<(-cutoff)),
                                     average.score=apply(data_rnai_htr,2,mean))
# colSums(data_rnai_htr<(-cutoff) - calculates the count of elements in each columns 
# that meet the conditions
data_rnai_impact_htr[order(-data_rnai_impact_htr$below.minus.0.5),]

#add a column with all the values above the cutoff - Cells that can grow when a gene is deleted
data_rnai_impact_htr$above.plus.0.5 <- colSums(data_rnai_htr>cutoff)

data_rnai_impact_htr <- data_rnai_impact_htr[order(-data_rnai_impact_htr$below.minus.0.5),]

# genes with significant impact
htr_gene_with_significant_impact <- 
  rownames(data_rnai_impact_htr[data_rnai_impact_htr$below.minus.0.5>10 | 
                                    data_rnai_impact_htr$above.plus.0.5>10,])

# Select only htr, gp and tf with significant impact 
subset_htr_gp_tf <- c(htr_gene_with_significant_impact, G_protein_with_significant_impact, tfs_gene_with_significant_impact)

# Subsetting
df_filtered_data_for_htr <- DepMap_RNAi[,..subset_htr_gp_tf]

htr_gp_tf_df_DT <- data.table::as.data.table(df_filtered_data_for_htr)  # GPCRs vs Gp
# Perform cor_target_map function from the cordial package to perform correlation only on selected targets in parallel for GPCRs vs Gp.
cor_results_htr_gp_tf <- cor_target_map(
  target = htr_gene_with_significant_impact,
  dataset = htr_gp_tf_df_DT,
  select_cols = colnames(htr_gp_tf_df_DT),
  self = "no") # to remove self-correlation

write.csv(cor_results_htr_gp_tf, file = "cor_results_htr_gp_tf.csv")
