library(data.table)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(grDevices)
library(ggplot2)
library(tidytext)
library(viridis)
library(htmltools)
library(igraph) # network analysis library

### Prepare the heatmap
# convert wide format to long format
dt_htr_gp_tf <- dcast(cor_results_htr_gp_tf, Correlation ~ Target, value.var = "r")

## Convert data.table into data frame to allow rownames 
dt_htr_gp_tf <- as.data.frame(dt_htr_gp_tf)
rownames(dt_htr_gp_tf) <- dt_htr_gp_tf$Correlation

### Convert data frame to matrix 
# All data
heatmap_matrix_htr_gp_tf <- as.matrix(dt_htr_gp_tf[,-1]) # remove the first column 
heatmap_matrix_htr_gp_tf[is.na(heatmap_matrix_htr_gp_tf)] <- 0 # replace N/a with zero 

### First step
# Calculate dynamic figure height according to number of rows
# This is optional, so can comment out
fig_height_unit <- nrow(heatmap_matrix_htr_gp_tf) * grid::unit(0.1, "inches")
fig_width_unit <- ncol(heatmap_matrix_htr_gp_tf) * grid::unit(1, "cm")
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
  breaks = c(min(heatmap_matrix_htr_gp_tf), 0, max(heatmap_matrix_htr_gp_tf)),
  colors = color_palette
)

### Second STEP
svglite::svglite(  # can change device to pdf() instead of svglite::svglite()
  "heatmap_matrix_htr_gp_tf.svg",
  width = fig_width,
  height = fig_height,
  # paper = "special"  # for use in pdf()
)

### Third STEP
# Heatmap
# Use the same code to run the all data
ComplexHeatmap::Heatmap(
  matrix = as.matrix(heatmap_matrix_htr_gp_tf),
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
browseURL("heatmap_matrix_htr_gp_tf.svg")

### Read the SVG content
# you can use htmltools to render it and open it on Rstudio
svg_content <- paste(readLines("heatmap_matrix_htr_gp_tf.svg"), collapse = "\n")
# Display in RStudio viewer
html_print(HTML(svg_content))


