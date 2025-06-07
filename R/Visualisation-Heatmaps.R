## LOAD PACKAGES 
library(data.table)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(grDevices)
library(ggplot2)
library(tidytext)
library(viridis)
library(htmltools)
        
## Create heatmap input - subset the data to have better visualization
# get significant results - choose q<0.05 and q<0.01 if you want to subset further
# create subset for cor_results_gpcrs_gp_tfs
cor_combined_q05 <- cor_results_gpcrs_gp_tfs[q<0.05] # select only q values less then 
cor_combined_q01 <- cor_results_gpcrs_gp_tfs[q<0.01] # to subset further select q less then

#cor_combined_top <- cor_combined_q05[, .SD[1:10], by = Target] # use for bar plot - this code has not been run
# Select only the r value < 0.05 for the cor_results_gpcrs_gp
cor_gpcrs_gp_q05 <- cor_results_gpcrs_gp[q<0.05]

# Filter for p values - significant values
cor_gpcrs_gp_p05 <- cor_results_gpcrs_gp[p<0.05]
# Order the significant values in descending order
cor_gpcrs_gp_p05 <- cor_gpcrs_gp_p05[order(cor_gpcrs_gp_p05$Correlation), ]
# Choose top 20 significant
top_20_gpcrs_gp <- head(cor_gpcrs_gp_p05, 100)

# Select only the r value < 0.05 for the cor_results_gpcrs_tf
cor_gpcrs_tfs_q05 <- cor_results_gpcrs_tfs[q<0.05]

### Prepare the heatmap
# convert wide format to long format - this data is used cor_combined_q01 because the data is very big
dt_wide <- dcast(cor_combined_q01, Correlation ~ Target, value.var = "r")
# For the combined subset data
dt_wide_subset_combined <- dcast(cor_results_subset_combined, Correlation ~ Target, value.var = "r")

# Convert wide format to long format - GPCRs vs Gp
dt_gpcrs_gp <- dcast(cor_gpcrs_gp_q05, Correlation ~ Target, value.var = "r")

# Convert wide format to long format - GPCRs vs TFs
dt_gpcrs_tfs <- dcast(cor_gpcrs_tfs_q05, Correlation ~ Target, value.var = "r")

## Convert data.table into data frame to allow rownames 
# All data
dt_wide <- as.data.frame(dt_wide)
rownames(dt_wide) <- dt_wide$Correlation
# Combined subset data
dt_wide_subset_combined <- as.data.frame(dt_wide_subset_combined) 
rownames(dt_wide_subset_combined) <- dt_wide_subset_combined$Correlation 
# GPCRs vs Gp
dt_gpcrs_gp <- as.data.frame(dt_gpcrs_gp)
rownames(dt_gpcrs_gp) <- dt_gpcrs_gp$Correlation
# GPCRs vs TFs
dt_gpcrs_tfs <- as.data.frame(dt_gpcrs_tfs)
rownames(dt_gpcrs_tfs) <- dt_gpcrs_tfs$Correlation

### Convert data frame to matrix 
# All data
heatmap_matrix <- as.matrix(dt_wide[,-1]) # remove the first column 
heatmap_matrix[is.na(heatmap_matrix)] <- 0 # replace N/a with zero 
# For combined data
heatmap_matrix_subset_combined <- as.matrix(dt_wide_subset_combined[,-1]) # remove the first column 
heatmap_matrix_subset_combined[is.na(heatmap_matrix_subset_combined)] <- 0 # replace N/a with zero 
# GPCRs vs Gp
heatmap_matrix_gpcrs_gp <- as.matrix(dt_gpcrs_gp[,-1])
heatmap_matrix_gpcrs_gp[is.na(heatmap_matrix_gpcrs_gp)] <- 0 #Clustered heatmap doesn't work with N/a values
# GPCRs vs TFs
heatmap_matrix_gpcrs_tfs <- as.matrix(dt_gpcrs_tfs[,-1])
heatmap_matrix_gpcrs_tfs[is.na(heatmap_matrix_gpcrs_tfs)] <- 0

### First step
# Calculate dynamic figure height according to number of rows
# This is optional, so can comment out
fig_height_unit <- nrow(heatmap_matrix_gpcrs_tfs) * grid::unit(0.1, "inches")
fig_width_unit <- ncol(heatmap_matrix_gpcrs_tfs) * grid::unit(1, "cm")
fig_height <- grid::convertY(
  fig_height_unit,
  unitTo = "inch",
  valueOnly = TRUE
) + 2 # margin
fig_width <- grid::convertX(
  fig_width_unit,
  unitTo = "inch",
  valueOnly = TRUE
) + 3 # margin

# Create colour palette
color_palette <- viridisLite::viridis(n = 3, alpha = 1, begin = 0.1, end = 0.9, option = "turbo")

# Create colour mapping function
color_mapping <- circlize::colorRamp2(
  breaks = c(min(heatmap_matrix_gpcrs_tfs), 0, max(heatmap_matrix_gpcrs_tfs)),
  colors = color_palette
)

### Second STEP
svglite::svglite(  # can change device to pdf() instead of svglite::svglite()
  "heatmap_matrix_gpcrs_tfs.svg",
  width = fig_width,
  height = fig_height,
  # paper = "special"  # for use in pdf()
)

### Third STEP
# Heatmap
# Use the same code to run the all data
ComplexHeatmap::Heatmap(
  matrix = as.matrix(heatmap_matrix_gpcrs_tfs),
  col = color_mapping,
  name = "r values", # legend_title,
  # Rows #
  row_title = "", # row_title,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 9, fontface = "plain"),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 4),
  show_row_dend = TRUE,
  row_dend_side = "right",
  row_dend_width = grid::unit(1.5, "cm"),
  cluster_rows = TRUE, # for dendogram on rows
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  # Columns #
  column_title = "", # column_title,
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 9, fontface = "plain"),
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 7),
  show_column_dend = FALSE,
  # column_dend_side = "top",
  # column_dend_height = unit(2, "cm"),
  cluster_columns = TRUE, # for dendogram on columns
  # clustering_distance_columns = "euclidean",
  # clustering_method_columns = "ward.D2",
  # Body #
  width = fig_width_unit, # ncol(hm.zcr)*unit(1, "cm"), # unit(fig_width - 1, "inches"),  # Heatmap body
  height = fig_height_unit, # nrow(hm.zcr)*unit(0.1, "inches"), # unit(fig_height - 0.5, "inches"),  # Heatmap body
  # heatmap_width = unit(fig_width, "inches"),  # Whole figure
  # heatmap_height = unit(fig_height, "inches"),  # Whole figure
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    legend_direction = "vertical", # "horizontal",
    legend_width = grid::unit(2, "cm"),
    title_position = "topcenter",  # "leftcenter", # "leftcenter-rot",
    # title = "Clustered heatmap",
    title_gp = gpar(fontsize = 5, fontface = "plain"),
    labels_gp = gpar(fontsize = 5)
  ),
  border = TRUE,
  border_gp = gpar(col = "black", alpha = 0.2),
  # Cell annotation #
  # cell_fun = function(j, i, x, y, width, height, fill) {
  #   grid.text(
  #     df.ss[i, j], x, y,
  #     gp = gpar(fontsize = 4)
  #   )
  # }
)
### Close the PDF file
# To view plot instead of saving, comment this out
dev.off()

### you can use this one to view the map using an external browser
browseURL("heatmap_matrix_gpcrs_tfs.svg")

### Read the SVG content
# you can use htmltools to render it and open it on Rstudio
svg_content <- paste(readLines("heatmap_matrix_gpcrs_tfs.svg"), collapse = "\n")
# Display in RStudio viewer
html_print(HTML(svg_content))









