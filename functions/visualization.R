# Create a standardized theme
standard_theme <- function() {
  theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      plot.title = element_text(size = 16)
    )
}



pca_plots <- function(explained_variance, nA, nB, pca_res, meta_subset, hc, 
                      k, color_var) {
 
  # Create a data frame with PC1/PC2 and your sample metadata
  df_pca <- data.frame(
    PCa = pca_res$x[, nA],
    PCb = pca_res$x[, nB],
    meta_subset)
  
  # Cluster samples (columns) using hierarchical clustering
  clusters <- cutree(hc, k)  # set k as desired (here: 3)
  
  # Add cluster info to PCA data frame
  df_pca$Cluster <- factor(clusters)
  # Plot PCa vs PCb, coloring points by cluster, label axes with
  pcA <- paste0("PC", nA)
  pcB <- paste0("PC", nB)
  # explained variance
  p <- ggplot(df_pca, aes_string(x = "PCa", y = "PCb", color = color_var)) +
    geom_point(size = 3, alpha = 0.5) +
    standard_theme() +
    theme_classic(base_size = 16) +
    labs(
      x = paste0(pcA, ": ", 
                 round(explained_variance[nA], 1), "% variance"),
      y = paste0(pcB, ": ", 
                 round(explained_variance[nB], 1), "% variance")
    )

  return(p)
}

elbow_plot <- function(explained_variance) {
  # First 25 principal components and their corresponding explained 
  # variance percentages
  df_variance <- data.frame(
    PC = seq_len(25),                 # Principal Component numbers 1 to 25
    # Variance explained by each PC (from earlier PCA)
    Variance = explained_variance[1:25]  
  )
  
  # Plot an "elbow plot" of explained variance for the top 
  # 25 principal components
  p <- ggplot(df_variance, aes(x = PC, y = Variance)) +
    geom_point(size = 3) +                # Add points for each pc
    standard_theme() +
    theme_minimal(base_size = 16) +                     # Use a minimal theme for a clean look
    labs(
      title = "Elbow Plot (Top 25 PCs)",  # Title for the plot
      x = "Principal Component",           # X-axis label
      y = "Explained Variance (%)"         # Y-axis label
    )

  return(p)
}

plot_hierarchical <- function(vsd, hc, k, sample_annotation_df, tcga_code) {
  
  
  # Select top 500 most variable genes
  high_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
  high_var_mat <- assay(vsd)[high_var_genes, ]
  
  # Perform hierarchical clustering of samples (columns)
  hc <- hclust(dist(t(high_var_mat)), method = "complete")
  clusters <- cutree(hc, k)
  
  # Generate colors automatically
  sample_annotation_colors <- list(
    Pathology = c("primary" = "black", "metastasis" = "grey",
                  "benign_neoplasia" = "red" ),
    Gender = c("male" = "#009E73", "female" = "#F0E442", 
               "missing_info" = "#D55E00"),
    tcga_code
  )
  
  column_annotation <- HeatmapAnnotation(
    df = sample_annotation_df,
    col = sample_annotation_colors,
    annotation_label = c("Pathology", "Gender", "tcga_code")  
  )
  
  # Draw the final heatmap
  p <- Heatmap(
    high_var_mat,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    show_row_names = TRUE,
    show_column_names = FALSE,
    name = "VST_Counts",
    column_split = clusters,
    top_annotation = column_annotation
  ) 
  
  return(p)
}

plot_genes_pca <- function(explained_variance, pca_res, vsd, meta_subset, pc) {
  # plot the 6 genes with larger 
  pc_loadings <- pca_res$rotation[, pc]
  pc_sorted_loadings <- sort(pc_loadings, decreasing = TRUE)
  pc_top_genes <- names(head(pc_sorted_loadings, n = 6))
  pc_top_genes_vsd <- data.frame(t(assay(vsd)[pc_top_genes, ]))
  
  df_pca <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    meta_subset)
  
  df_genes_long <- data.frame(
    pc_top_genes_vsd,  
    df_pca[, c("PC1", "PC2", "Pathology")]
  )
  
  df_long <- pivot_longer(
    df_genes_long,
    cols = all_of(pc_top_genes),
    names_to = "Gene",
    values_to = "Expression"
  )
  
  ggplot(df_long, aes(x = PC1, y = PC2, color = Expression)) +
    geom_point(size = 3, alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red") +
    standard_theme() +
    theme_classic(base_size = 16) +
    labs(title = paste0("Genes from PC", pc), 
         x = paste0("PC1: ", round(explained_variance[1], 1), "% variance"), 
         y = paste0("PC2: ", round(explained_variance[2], 1), "% variance"),
         color = "Expression (VST)") +
    facet_wrap( ~ Gene)
}

# Plot histogram of log2FoldChange
plot_2fold <- function(res) {
  p <- hist(res$log2FoldChange,
       breaks = 50,
       xlab = "log2 Fold Change",
       main = "Distribution of log2 Fold Changes",
       col = "lightblue")
  # Add a vertical line at zero
  abline(v = 0, col = "red", lty = 2)
  return(p)
}

# Plot 2fold genes
plot_2fold_genes <-function(explained_variance, pca_res, res, vsd, top_genes, 
                            meta_subset) {
  
  pc2_top_genes_vsd <- data.frame(t(assay(vsd)[top_genes, ]))
  
  df_pca <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    meta_subset)
  
  df_genes_long <- data.frame(
    pc2_top_genes_vsd,  
    df_pca[, c("PC1", "PC2", "Pathology")]
  )
  
  df_long <- pivot_longer(
    df_genes_long,
    cols = all_of(top_genes),
    names_to = "Gene",
    values_to = "Expression"
  )
  
  # enforce facet order as in top_6
  df_long$Gene <- factor(df_long$Gene, levels = top_genes)
  
  p <- ggplot(df_long, aes(x = PC1, y = PC2, color = Expression)) +
    geom_point(size = 3, alpha = 0.5) +
    scale_color_gradient(low = "blue", high = "red") +
    standard_theme() +
    theme_classic(base_size = 16) +
    labs(title = "Genes with larger differential expression",
         x = paste0("PC1: ", round(explained_variance[1], 1), "% variance"), 
         y = paste0("PC2: ", round(explained_variance[2], 1), "% variance"),
         color = "Expression (VST)") +
    facet_wrap( ~ Gene) 
    
  return(p)
}

plot_differential <-function(top_genes, res_reordered) {
  
  # Check dimensions
  print("inside plot_differential")
  print(top_genes)
  print(nrow(top_genes))        
  print(nrow(res_reordered))
  
  p <- EnhancedVolcano(
    res_reordered, # Input data frame with differential expression results
    x = "log2FoldChange",        # Column for log2 fold changes
    y = "padj",                  # Column for adjusted p-values
    #lab = res_reordered$gene,    # Label the genes based on their gene names
    lab = rownames(res_reordered), # Use rownames instead
    selectLab = intersect(rownames(res_reordered), top_genes),
    pCutoff = 0.00001,           # Significance threshold for p-value
    FCcutoff = 2,                # Threshold for log2 fold-change
    labSize = 5,                 # Size of the gene labels
    pointSize = 3,     # smaller points (default is ~2)
    colAlpha = 0.3,
    drawConnectors = TRUE,
    widthConnectors = 1,
    max.overlaps = Inf,
    title = "Cluster 1 vs Cluster 2 DEGs", # Title of the plot
    subtitle = "",               # No subtitle
    caption = "",                 # No caption
  ) +
    standard_theme() +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") 
  
  return(p)
}