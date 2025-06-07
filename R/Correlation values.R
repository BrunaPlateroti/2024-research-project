# Load packages
require(magrittr)
require(tidyr)
require(purrr)
require(future)
require(furrr)
require(collapse)
require(data.table)
library(cordial)

### TASK 2 - Perform correlation analysis on gene dataset
## Pearson correlation analysis for an entire dataset using cordial package function cor_map()

#Selected genes with data in dep map are combined - all data combined
combined_gpcrs_gp_tfs <- c(gpcrs_gene_with_significant_impact, G_protein_with_significant_impact, tfs_gene_with_significant_impact)

# Select only GPCRs and G protein
combined_gpcrs_gp <- c(gpcrs_gene_with_significant_impact, G_protein_with_significant_impact)

# Select only GPCRs and TFs
combined_gpcrs_tf <- c(gpcrs_gene_with_significant_impact, tfs_gene_with_significant_impact)

## Subsetting
# Subset the DepMap_RNAi data table to include only the columns specified in the combined_gpcrs_gp_tfs vector.
# The resulting data table, df_gpcrs_gp_tfs, will contain all rows but only the selected columns taken from combined_gpcrs_gp_tfs
df_gpcrs_gp_tfs <- DepMap_RNAi[,..combined_gpcrs_gp_tfs]

df_gpcrs_gp <- DepMap_RNAi[,..combined_gpcrs_gp] # GPCRs vs Gp

df_gpcrs_tfs <- DepMap_RNAi[,..combined_gpcrs_tf] # GPCRs vs TFs

# Convert the df_gpcrs_gp_tfs object into a data table to leverage the efficient 
# data manipulation capabilities of the data.table package.

gpcrs_gp_tfs_df_DT <- data.table::as.data.table(df_gpcrs_gp_tfs)
cor_results_gpcrs_gp_tfs <- cor_map(
  gpcrs_gp_tfs_df_DT,
  select_cols = colnames(gpcrs_gp_tfs_df_DT),
  self = "no")# to remove self-correlation

gpcrs_gp_df_DT <- data.table::as.data.table(df_gpcrs_gp)  # GPCRs vs Gp
# Perform cor_target_map function from the cordial package to perform correlation only on selected targets in parallel for GPCRs vs Gp.
cor_results_gpcrs_gp <- cor_target_map(
  target = gpcrs_gene_with_significant_impact,
  dataset = gpcrs_gp_df_DT,
  select_cols = colnames(gpcrs_gp_df_DT),
  self = "no") # to remove self-correlation

gpcrs_tf_df_DT <- data.table::as.data.table(df_gpcrs_tfs) # GPCRs vs TFs
# Perform cor_target_map function from the cordial package to perform correlation only on selected targets in parallel for GPCRs vs TFs.
cor_results_gpcrs_tfs <- cor_target_map(
  target = gpcrs_gene_with_significant_impact,
  dataset = gpcrs_tf_df_DT,
  select_cols = colnames(gpcrs_tf_df_DT),
  self = "no") # to remove self-correlation

#cor_results_gpcrs_gp_tfs <- cor_results_gpcrs_gp_tfs[Target != "Correlation"] # to remove self correlation (to adjust) - it doesn't work
write.csv(cor_results_gpcrs_gp_tfs, file = "results_gpcrs_gp_tfs.csv") # correlation of the combined data

write.csv(cor_results_gpcrs_gp, file = "results_gpcrs_gp.csv")

write.csv(cor_results_gpcrs_tfs, file = "results_gpcrs_tfs.csv")

## To subset further the heatmap 
# Subset only key genes for TFs, GPCRs and Gp 
# Check before if those are part of the DepMap data

subset_tfs <- c("PML", "RUNX1","TP53", "ATF1", "ADCYAP1R1") #(choose 5 key genes for each subset)
#subset_tfs <- subset_tfs %in% combined_gpcrs_gp_tfs

subset_gps <- c("HRAS", "NRAS", "KRAS", "MRAS", "CDC42") 
subset_gpcrs <- c("LPAR1", "LPAR2", "LPAR5", "LPAR6", "S1PR1")

# to check if gene are part of the data
#"AML"%in%cor_results_gpcrs_gp_tfs[,Target] 

# combined the subsets
subset_combined <- c(subset_tfs, subset_gps, subset_gpcrs)

df_gpcrs_gp_tfs <- DepMap_RNAi[,..combined_gpcrs_gp_tfs]
df_subset_combined <- DepMap_RNAi[,..subset_combined] # get the genes from the original data


# Convert df_subset_combined into data.table
subset_combined_DT <- data.table::as.data.table(df_subset_combined)
# Perform correlation on the subset
cor_results_subset_combined <- cor_map(
  subset_combined_DT,
  select_cols = colnames(subset_combined_DT),
  self = "no") # to remove self-correlation

#cor_results_subset_combined <- cor_results_subset_combined[Target != "Correlation"]#
write.csv(cor_results_subset_combined, file = "results_subset_combined.csv")




 


