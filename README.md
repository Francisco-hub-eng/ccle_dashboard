# CCLE Dashboard: RNA-seq Cancer Genome Atlas Analysis

An interactive R Shiny application for analyzing RNA-seq data from the Cancer Cell Line Encyclopedia (CCLE) using advanced bioinformatics techniques and AI-powered biological interpretation.

## Project Overview

This project provides a comprehensive analysis platform for cancer genomics data, specifically focusing on differential gene expression analysis across different cancer types and primary sites. The application combines traditional bioinformatics approaches with modern AI capabilities to provide biological insights into cancer-related gene expression patterns.

## Key Features

- **Interactive Cancer Type Selection**: Analyze RNA-seq data from multiple primary sites available in the CCLE dataset
- **Dynamic Clustering Analysis**: Perform hierarchical clustering with user-defined number of clusters (2-10)
- **Comprehensive Visualization Suite**: 
  - Principal Component Analysis (PCA) plots with explained variance
  - Hierarchical clustering heatmaps
  - Volcano plots for differential expression
  - Gene expression overlays on PCA space
- **Statistical Analysis**: Differential expression analysis using DESeq2
- **AI-Powered Interpretation**: Integration with Perplexity AI to provide biological context for top differentially expressed genes

## Scientific Methodology

### Data Processing Pipeline

1. **Data Filtering**: 
   - Subset CCLE data by selected primary site
   - Remove samples with missing pathology information
   - Filter genes with low expression (< 100 total counts)

2. **Normalization**: 
   - Variance Stabilizing Transformation (VST) using DESeq2
   - Selection of top 500 most variable genes for dimensionality reduction

3. **Dimensionality Reduction**:
   - Principal Component Analysis (PCA) on variance-stabilized data
   - Hierarchical clustering using Euclidean distance and complete linkage

4. **Differential Expression Analysis**:
   - DESeq2-based analysis comparing expression between hierarchical clusters
   - Multiple testing correction using Benjamini-Hochberg method
   - Identification of significantly differentially expressed genes (padj < 0.001)

### Statistical Techniques Employed

- **Variance Stabilization**: VST transformation to account for mean-variance relationship in RNA-seq data
- **Hierarchical Clustering**: Complete linkage clustering for sample grouping
- **Principal Component Analysis**: Linear dimensionality reduction for data visualization and exploration
- **Differential Expression**: Negative binomial generalized linear models via DESeq2
- **Multiple Testing Correction**: Benjamini-Hochberg procedure for FDR control

## AI Integration

The application leverages **Perplexity AI** to provide biological interpretation of the most significantly differentially expressed genes. This integration:

- Queries the AI with top differentially expressed genes and the specific cancer type
- Provides contextual biological information about gene functions in cancer
- Includes academic citations for further reading
- Delivers interpretations in a user-friendly format within the dashboard

The AI component enhances the analytical pipeline by bridging the gap between computational results and biological understanding, making the findings more accessible to researchers with varying levels of bioinformatics expertise.

## Technical Implementation

### Architecture
- **Frontend**: R Shiny with responsive UI components
- **Backend**: R with Bioconductor packages for genomics analysis
- **Visualization**: ggplot2, ComplexHeatmap, and EnhancedVolcano
- **Statistical Computing**: DESeq2, biomaRt integration
- **AI Integration**: RESTful API calls to Perplexity AI service

### Key Dependencies
```r
library(shiny)
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(biomaRt)
library(httr)
library(jsonlite)
```

### File Structure
```
├── app.R                 # Main Shiny application
├── functions/
│   ├── data_processing.R # Data filtering, DESeq2 analysis, AI queries
│   └── visualization.R   # Plotting functions and themes
└── data/
    ├── 1_ccle_counts.csv # RNA-seq count matrix
    └── 1_ccle_meta.csv   # Sample metadata
```

## Usage

1. Select a cancer primary site from the dropdown menu
2. Choose the number of clusters for hierarchical clustering
3. Explore the various visualization panels:
   - Sample composition (pathology, gender, TCGA codes)
   - PCA plots with explained variance
   - Gene expression patterns contributing to principal components
   - Hierarchical clustering heatmap
   - Differential expression results with volcano plots
4. Review AI-generated biological interpretations of top differentially expressed genes

## Data Source

This analysis utilizes data from the **Cancer Cell Line Encyclopedia (CCLE)**, a comprehensive resource containing genomic and transcriptomic profiles of human cancer cell lines across diverse tissue types and molecular subtypes.

## Future Directions

- Integration of additional omics data types (mutations, copy number variations)
- Extended AI capabilities for pathway analysis
- Batch effect correction for multi-study integration
- Enhanced visualization with interactive plots

## Contributing

This project demonstrates the integration of traditional bioinformatics workflows with modern AI capabilities, providing a template for similar genomics analysis applications.

---

*This project showcases advanced bioinformatics analysis techniques combined with AI-powered interpretation, representing a modern approach to cancer genomics research.*