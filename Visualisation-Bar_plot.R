## LOAD PACKAGES
library(ggplot2)
library(tidytext)
library(viridis)

### Bar plot to illustrate example of Gprotein such as NRAS taking key genes of the NRAS pathway. 
# Create a subset of cor_results_subset_combined data table containing rows where the 'Target' column is either "HRAS", "NRAS", "KRAS", "MRAS"
#key_combined_cor_q1 <- cor_combined_q01[Target %in% c("HRAS", "NRAS", "KRAS", "MRAS")]
key_combined_cor <- cor_results_subset_combined[q<0.05] # if you want narrow down further

# Select the top 10 rows within each group defined by the 'Target' column,
# based on the values in the 'r' column in descending order.
top10_each_group <- cor_results_subset_combined[, .SD[order(-r)[1:10]], by = Target]
# Select the bottom 10 rows within each group defined by the 'Target' column,
# based on the values in the 'r' column in ascending order.
bottom10_each_group <- cor_results_subset_combined[, .SD[order(r)[1:10]], by = Target]

# Combine the top 10 and bottom 10 rows from each group (defined by 'Target') into a single data table.
combined_each_group <- rbindlist(list(top10_each_group, bottom10_each_group))
# Add a new column 'Direction' to the combined_each_group data table.
# The 'Direction' column indicates whether the value in the 'r' column is positive or negative.
results <- combined_each_group[, `:=`( Direction = data.table::fifelse(r > 0, "Positive", "Negative") )]

# Lets plot the all data without taking the top and bottom 10
results_all <- cor_results_subset_combined[, `:=`( Direction = data.table::fifelse(r > 0, "Positive", "Negative") )]

top_bar_plot <- function(
    x, 
    save_fig, 
    save_name, 
    facet_fig = "wrap", 
    width = 10, 
    height = 3
) {
  
  
## Bar plots
  combine_all_col_top <- ggplot(
    x,
    aes(
      x = r,
      #y = reorder_within(
      # Correlation, by = Index, within = Target
      #),
      y = Correlation,
      fill = Direction
    )
  ) +
    geom_col() +
    geom_text(
      aes(
        label = format(q, digits=3, nsmall=2, scientific=TRUE),
        x = fifelse(Direction == "Positive", -0.4, 0.4) #if_else
      ), 
      colour = "DimGrey",
      size = rel(2)
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(
      yintercept = seq(0, length(unique(x$Correlation))) + 0.5, 
      alpha = 0.1
    ) +
    scale_x_continuous(limits = c(-0.6, 0.6)) +
    scale_y_reordered(limits = rev) +
    scale_fill_manual(
      values = c(
        viridis_pal(begin = 0.1, end = 0.9, option = "D")(3)[3], 
        viridis_pal(begin = 0.1, end = 0.9, option = "D")(3)[1]
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(size = rel(0.5)),
      axis.title = element_text(size = rel(0.7)),
      strip.text = element_text(size = rel(0.5))
    ) +
    labs(
      x = "r-value",
      y = "Correlation"
    )
  
  if(facet_fig == "wrap") {
    # Facet according to Gene
    combine_all_col_top_facet <- combine_all_col_top + 
      facet_wrap(
        . ~ Target, scales = "free_y", ncol = 5
      )
    # Save
    if(isTRUE(save_fig)) {
      ggsave(
        file = save_name, 
        plot = combine_all_col_top_facet,
        width = width, height = height, units = "in", dpi = "retina"
      )
    }
    # View plot
    combine_all_col_top_facet
    
  } else if(facet_fig == "grid") {
    combine_all_col_top_facet <- combine_all_col_top + 
      facet_grid(
        Direction ~ Target, scales = "free_y", space = "free_y"
      )
    # Save
    if(isTRUE(save_fig)) {
      ggsave(
        file = save_name, 
        plot = combine_all_col_top_facet,
        width = width, height = height, units = "in", dpi = "retina"
      )
    }
    # View plot
    print(combine_all_col_top_facet)
    
  } else if(facet_fig == "none") {
    # Save
    if(isTRUE(save_fig)) {
      ggsave(
        file = save_name, 
        plot = combine_all_col_top,
        width = width, height = height, units = "in", dpi = "retina"
      )
    }
    # View plot
    print(combine_all_col_top)
  }
}

bar_plot_res_all <- top_bar_plot(
  results_all,
  save_fig = TRUE, 
  save_name = "All_bar_plot.pdf", 
  facet_fig = "wrap",
  width = 9, 
  height = 7.8
)

print(bar_plot_res)
