library(ggplot2)


gpcrs_gp_scaled <- scale(heatmap_matrix_gpcrs_gp)

# Perform PCA
pca_result_gpcrs_gp <- prcomp(gpcrs_gp_scaled, center = TRUE, scale. = TRUE)

# View the proportion of variance explained by each principal component
summary(pca_result_gpcrs_gp)

# Plot a Scree Plot to visualize how much variance each principal component explains
plot(pca_result_gpcrs_gp, type = "l", main = "Scree Plot")# tell how many components you should capture

# PCA Plot of the first two principal components
pca_result_gpcrs_gp_df <- as.data.frame(pca_result_gpcrs_gp$x)  # Convert PCA results into a data frame

# Access proportion of variance directly from the PCA result
# Proportion of variance explained is stored in pca_result$sdev (standard deviations of PCs)

# Calculate variance explained by each principal component
pca_variance_gpcrs_gp <- pca_result_gpcrs_gp$sdev^2

# Calculate proportion of variance
pca_variance_proportion_gpcrs_gp <- pca_variance_gpcrs_gp / sum(pca_variance_gpcrs_gp)

# Now you can use these proportions for the axis labels in your ggplot
pc1_variance_gpcrs_gp <- round(pca_variance_proportion_gpcrs_gp[1] * 100, 1)
pc2_variance_gpcrs_gp <- round(pca_variance_proportion_gpcrs_gp[2] * 100, 1)

# Proceed with the ggplot for PCA
ggplot(pca_result_gpcrs_gp_df, aes(x = PC1, y = PC2)) +
  geom_point() +
  ggtitle("PCA Plot: PC1 vs PC2") +
  xlab(paste0("PC1 (", pc1_variance_gpcrs_gp, "% Variance)")) +
  ylab(paste0("PC2 (", pc2_variance_gpcrs_gp, "% Variance)"))




## For Biplot - optional
biplot(pca_result_gpcrs_gp, scale = 0)
# Loadings (PC loadings for variables)
loadings <- as.data.frame(pca_result_gpcrs_gp$rotation)  

# Plot PCA points (observations, e.g., GPCRs or genes)
ggplot(pca_result_gpcrs_gp_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = rownames(pca_result_gpcrs_gp_df)), size = 2) +  # Plot points and color by observation names
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +  # Plot loadings as arrows
  geom_text(data = loadings, aes(x = PC1 * 5, y = PC2 * 5, label = rownames(loadings)), vjust = 1.5) +  # Loadings labels
  ggtitle("PCA Biplot: PC1 vs PC2") +
  xlab(paste0("PC1 (", pc1_variance_gpcrs_gp, "% Variance)")) +  # Adjusted x-axis label with variance explained
  ylab(paste0("PC2 (", pc2_variance_gpcrs_gp, "% Variance)")) +  # Adjusted y-axis label with variance explained
  theme_minimal()  # A cleaner theme
