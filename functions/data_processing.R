# Function to query Perplexity API
query_perplexity <- function(top_genes, selected_organ, api_key) {
  
  # Build query string (respecting 80-character limit)
  top_genes_string <- paste(top_genes, collapse = ",")
  query <- paste0(
    "Explain what could be the role of the genes ", 
    top_genes_string, 
    " in the ", 
    selected_organ, 
    " cancer. These genes were the ones with the largest differential ",
    "expression from the bulk RNA-seq tcga database from the ", 
    selected_organ, 
    " primary site. Max 200 words."
  )
  
  # API endpoint
  url <- "https://api.perplexity.ai/chat/completions"
  
  # Request headers
  headers <- add_headers(
    "Authorization" = paste("Bearer", api_key),
    "Content-Type" = "application/json"
  )
  
  # Request body
  body <- list(
    model = "sonar",
    messages = list(
      list(role = "user", content = query)
    ),
    search_mode = "academic",
    max_tokens = 500,
    temperature = 0.2
  )
  
  # Make API request
  response <- POST(
    url = url,
    body = body,
    encode = "json",
    headers
  )
  
  print(response)
  
  # Parse and return response
  if (status_code(response) == 200) {
    
    # Get response content as text first
    response_text <- content(response, "text", encoding = "UTF-8")
    
  } else {
    stop("API request failed with status: ", status_code(response))
  }
  
  return(response_text)
}

# Safe function to extract content
extract_content <- function(response_string) {
  
  # Find "content": " and get everything after it
  start_marker <- '"content": "'
  start_pos <- regexpr(start_marker, response_string, fixed = TRUE)
  
  if (start_pos[1] == -1) return("Content not found")
  
  # Extract from after the marker
  content_start <- start_pos[1] + nchar(start_marker)
  substring_after <- substr(response_string, content_start, nchar(response_string))
  
  # Find the closing pattern and extract content
  end_pos <- regexpr('"}', substring_after, fixed = TRUE)
  
  if (end_pos[1] == -1) return("End not found")
  
  # Get the content
  content <- substr(substring_after, 1, end_pos[1] - 1)
  
  # Clean up escaped characters
  content <- gsub('\\\\n', '\n', content)
  content <- gsub('\\\\"', '"', content)
  
  return(content)


}

extract_citations <- function(response) {
  # Find the start pattern "citations": [
  start_pattern <- '"citations": ['
  start_pos <- regexpr(start_pattern, response, fixed = TRUE)[1]
  
  if (start_pos == -1) {
    return("Citations not found")
  }
  
  # Move to after the start pattern
  content_start <- start_pos + nchar(start_pattern)
  remaining <- substr(response, content_start, nchar(response))
  
  # Find the closing bracket ]
  end_pos <- regexpr(']', remaining, fixed = TRUE)[1]
  
  if (end_pos == -1) {
    return("Citations end not found")
  }
  
  # Extract everything between [ and ]
  citations_text <- substr(remaining, 1, end_pos - 1)
  
  return(citations_text)
}


clean_citations <- function(citations_text) {
  # Split by comma to get individual URLs
  urls <- unlist(strsplit(citations_text, ', ', fixed = TRUE))
  
  # Remove all quotes and escape characters
  urls <- gsub('[\\"\\\\]', '', urls)
  
  # Remove any leading/trailing whitespace
  urls <- trimws(urls)
  
  return(urls)
}

create_sample_annotation_colors <- function(sample_annotation_df) {
  # Get unique TCGA codes from the dataframe
  unique_tcga_codes <- unique(sample_annotation_df$tcga_code)
  unique_tcga_codes <- unique_tcga_codes[!is.na(unique_tcga_codes)]  # Remove NAs
  
  # Number of unique codes
  n_codes <- length(unique_tcga_codes)
  
  # Generate visually distinct colors
  if (n_codes <= 10) {
    # Use Set3 palette for up to 10 categories (good contrast and visibility)
    colors <- RColorBrewer::brewer.pal(max(3, n_codes), "Set3")[1:n_codes]
  } else if (n_codes <= 20) {
    # Use a combination of Set3 and Paired palettes for more categories
    colors1 <- RColorBrewer::brewer.pal(12, "Set3")
    colors2 <- RColorBrewer::brewer.pal(8, "Set2")
    colors <- c(colors1, colors2)[1:n_codes]
  } else {
    # For many categories, use rainbow with good separation
    colors <- rainbow(n_codes, s = 0.8, v = 0.8)
  }
  
  # Create named vector
  names(colors) <- unique_tcga_codes
  
  return(colors)
}


filter_ccle_data_by_organ <- function(selected_organ = NULL, k = 2,
                                      verbose = FALSE) {
  
  print(selected_organ)
  # Filter the ccle_meta dataframe to keep only rows where CCLE_ID with "CENTRAL"
  #ccle_meta_subset <- ccle_meta[grepl(selected_organ, ccle_meta$Site_Primary), ]
  
  # Further filter the subset to remove rows where the "Pathology" column is
  # NA or an empty string
  #ccle_meta_subset <- ccle_meta_subset[!(is.na(ccle_meta_subset$Pathology) |
  #                                         ccle_meta_subset$Pathology == ""), ]
  
  # Subset the ccle_counts matrix to keep only columns that match 
  # the remaining CCLE_IDs
  #ccle_counts_subset <- ccle_counts[, colnames(ccle_counts) %in% 
  #                                    ccle_meta_subset$CCLE_ID]
  
  # CRITICAL: Ensure exact matching
  #common_samples <- intersect(ccle_meta_subset$CCLE_ID, 
  #                            colnames(ccle_counts_subset))
  
  #ccle_meta_subset <- ccle_meta_subset[ccle_meta_subset$CCLE_ID %in% 
  #                                       common_samples, ]
  #ccle_counts_subset <- ccle_counts_subset[, common_samples]
  
  # ccle_meta_subset <- "none"
  
  index_a <- ccle_meta$Site_Primary %in% c(selected_organ)
  index_b = ccle_meta$Pathology != ''
  ccle_counts_subset <- ccle_counts[,index_a & index_b]
  ccle_meta_subset <- ccle_meta[index_a & index_b,]
  
  print("creating the subsets")
  
  # Check number of unique values in Pathology column
  unique_pathology_count <- length(unique(ccle_meta_subset$Pathology))
  
  if (unique_pathology_count == 1) {
    # All samples have the same pathology value - use intercept-only design
    dds <- DESeqDataSetFromMatrix(countData = data.frame(ccle_counts_subset),
                                  colData = ccle_meta_subset,
                                  design = ~ 1)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = data.frame(ccle_counts_subset),
                                  colData = ccle_meta_subset,
                                  design = ~ Pathology)
  }
  
  # filter only those genes with total more than 100 counts
  dds_100 <- rowSums(counts(dds)) >= 100
  
  # keep only those genes with more than 100 counts
  dds <- dds[dds_100,]
  
  vsd <- vst(dds, blind = FALSE)
  
  print("variance analysis")
  
  
  # Compute variance for each gene (row) across all samples
  variances <- rowVars(assay(vsd))
  
  # Get indices of top 500 most variable genes (highest variance first)
  top_var_genes <- order(variances, decreasing = TRUE)[1:500]
  
  # Extract the expression matrix for these genes only
  high_var_mat <- assay(vsd)[top_var_genes, ]
  
  # check that the matrix is full numeric
  high_var_mat[is.na(high_var_mat)] <- 0
  

  df <- as.data.frame(high_var_mat)
  # Drop non-numeric columns
  df_numeric <- df[, sapply(df, is.numeric), drop=FALSE]
  df_numeric <- df_numeric[, unlist(lapply(df_numeric, function(x) is.numeric(x) & !is.factor(x))), drop=FALSE]
  
  high_var_mat <- as.matrix(df_numeric)
  storage.mode(high_var_mat) <- "numeric"
  
  print("high_var_mat information")
  str(high_var_mat)
  summary(high_var_mat)

  # Perform PCA on transposed matrix (samples as rows, genes as columns)
  pca_res <- prcomp(t(high_var_mat), center = TRUE, scale. = FALSE)
  
  # Calculate the percent variance explained by each principal component
  explained_variance <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
  
  # do hierarchical clustering
  hc <- hclust(dist(t(high_var_mat)), method = "complete")
  
  print("pca by expression level")
  
  
  # Check dimensions
  print(nrow(pca_res$x))        
  print(nrow(ccle_meta_subset))
  
  ##########################################@
  # for the differential expression
  #row variance from the variance-stabilized data (vsd)
  high_var_genes <- head(order(rowVars(assay(vsd)), 
                               decreasing = TRUE), 50)
  
  # Extract the expression matrix of the high variance genes
  high_var_mat <- assay(vsd)[high_var_genes,]
  
  
  # Cut the dendrogram into two clusters
  clusters <- cutree(hc, k)
  
  # Add the cluster assignment to the metadata of the CCLE dataset
  ccle_meta_subset$Cluster <- factor(clusters)
  
  # Create a DESeqDataSet object from the CCLE count data and
  # metadata with a design based on cluster assignment
  dds <- DESeqDataSetFromMatrix(countData = data.frame(ccle_counts_subset),
                                colData = ccle_meta_subset,
                                design = ~ Cluster)
  
  # Relevel the 'Cluster' factor so that the reference level is cluster "2"
  dds$Cluster <- relevel(dds$Cluster, ref = "2")
  
  # Filter out genes with low counts (those with fewer than 
  # 100 reads across all samples)
  dds_low_counts <- dds[rowSums(counts(dds)) >= 100,]
  
  # Perform differential expression analysis using DESeq
  dds_deseq <- DESeq(dds_low_counts)
  
  res <- results(dds_deseq)
  
  sample_annotation_df <- ccle_meta_subset
  sample_annotation_df <- sample_annotation_df %>% 
    dplyr::select(Pathology, Gender, tcga_code)
  
  tcga_code <- create_sample_annotation_colors(sample_annotation_df)
  
  print("differential expression analysis")
  
  
  ########## get top genes
  
  # Remove rows with NA in padj column
  res <- res[!is.na(res$padj), ]
  
  res_filter <- res[res$padj < 0.001, , drop = FALSE]
  
  print("res_filter")
  
  # filter by the genes with p-value adjusted less than value
  res_filter <- res[res$padj < 0.001, , drop = FALSE]
  
  print("res_reordered")
  
  res_reordered <- res_filter[order(-abs(res_filter$log2FoldChange)), , 
                              drop = FALSE]
  print("res_reordered")
  
  
  # lowest differential expression gen
  idx_low <- order(res_filter$log2FoldChange, na.last = NA)[1:3]
  names_low <- rownames(res_filter)[idx_low]
  
  # Highest differential expression gen
  idx_high <- order(-res_filter$log2FoldChange, na.last = NA)[1:3]
  names_high <- rownames(res_filter)[idx_high]
  
  top_genes <- c(names_low, names_high)
  
  print("top genes by differential expression analysis")
  
  # Get interpretation by GenAI
  # Example usage
  api_key <- Sys.getenv("PERPLEXITY_API_KEY")

  # Execute query
  query_response <- query_perplexity(top_genes, selected_organ, api_key)
  
  #print(query_response)
  
  content_only <- extract_content(query_response)
  
  #print(content_only)
  
  citations_dirty <- extract_citations(query_response)
  
  citations <- clean_citations(citations_dirty)
  
  #content_only <- "Interpretation done by genAI"
  
  #citations <- c("https://pmc.ncbi.nlm.nih.gov/articles/PMC19327/",   
                 #"https://pmc.ncbi.nlm.nih.gov/articles/PMC5505673/")
  
  print("All the analysis done without any problem")
  
  return(list(ccle_meta_subset, ccle_counts_subset, dds, vsd, 
              pca_res, explained_variance, hc, 
              dds_low_counts, res, res_reordered, top_genes, k,
              content_only, citations, sample_annotation_df, tcga_code))
}


dds_subset_data <- function(meta_subset, counts_subset, 
                            verbose = False) {
  
  dds <- DESeqDataSetFromMatrix(countData = data.frame(counts_subset),
                                colData = meta_subset,
                                design = ~ Pathology)
  
  # filter only those genes with total more than 100 counts
  dds_100 <- rowSums(counts(dds)) >= 100
  
  # keep only those genes with more than 100 counts
  dds <- dds[dds_100,]
  
  vsd <- vst(dds)
  
  #print(head(data.frame(counts(dds))[,1:5], 10))
  
  return(list(dds, vsd))
}

# do PCA from the data subset
pca_subset <-function(vsd, verbose = False) {
  # Compute variance for each gene (row) across all samples
  variances <- rowVars(assay(vsd))
  
  # Get indices of top 500 most variable genes (highest variance first)
  top_var_genes <- order(variances, decreasing = TRUE)[1:500]
  
  # Extract the expression matrix for these genes only
  high_var_mat <- assay(vsd)[top_var_genes, ]
  
  
  
  
  
  
  high_var_mat <- as.matrix(df_numeric)
  storage.mode(high_var_mat) <- "numeric"
  
  
  
  
  # Perform PCA on transposed matrix (samples as rows, genes as columns)
  pca_res <- prcomp(t(high_var_mat), center = TRUE, scale. = FALSE)
  
  # Calculate the percent variance explained by each principal component
  explained_variance <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
  
  # do hierarchical clustering
  hc <- hclust(dist(t(high_var_mat)), method = "complete")
  
  return(list(pca_res, explained_variance, hc))
}

dea <- function(vsd, hc, meta_subset, counts_subset) {
  #row variance from the variance-stabilized data (vsd)
  high_var_genes <- head(order(rowVars(assay(vsd)), 
                               decreasing = TRUE), 50)
  
  # Extract the expression matrix of the high variance genes
  high_var_mat <- assay(vsd)[high_var_genes,]
  
  
  # Cut the dendrogram into two clusters
  clusters <- cutree(hc, k = 2)
  
  # Add the cluster assignment to the metadata of the CCLE dataset
  meta_subset$Cluster <- factor(clusters)
  
  # Create a DESeqDataSet object from the CCLE count data and
  # metadata with a design based on cluster assignment
  dds <- DESeqDataSetFromMatrix(countData = data.frame(counts_subset),
                                colData = meta_subset,
                                design = ~ Cluster)
  
  # Relevel the 'Cluster' factor so that the reference level is cluster "2"
  dds$Cluster <- relevel(dds$Cluster, ref = "2")
  
  # Filter out genes with low counts (those with fewer than 
  # 100 reads across all samples)
  dds <- dds[rowSums(counts(dds)) >= 100,]
  
  # Perform differential expression analysis using DESeq
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  return(list(dds, res))
}

# find top genes
top_genes_2fold <- function(res) {
  
  # Remove rows with NA in padj column
  res_clean <- res[!is.na(res$padj), ]
  
  res_filter <- res_clean[res_clean$padj < 0.001, , drop = FALSE]
  
  # filter by the genes with p-value adjusted less than value
  res_filter <- res[res$padj < 0.001, , drop = FALSE]
  
  res_reordered <- res_filter[order(-abs(res_filter$log2FoldChange)), , 
                              drop = FALSE]
  
  # lowest differential expression gen
  idx_low <- order(res_filter$log2FoldChange, na.last = NA)[1:3]
  names_low <- rownames(res_filter)[idx_low]
  
  # Highest differential expression gen
  idx_high <- order(-res_filter$log2FoldChange, na.last = NA)[1:3]
  names_high <- rownames(res_filter)[idx_high]
  
  top_genes <- c(names_low, names_high)
  
  print("top_genes_2fold")
  print(top_genes)
  print(nrow(top_genes))        
  print(nrow(res_reordered))
  
  return(list(top_genes, res_reordered))
}
