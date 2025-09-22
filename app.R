#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(DT)
library(tidyr)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(DESeq2)
library(httr)
library(jsonlite)
library(rmarkdown)

# source the function files
source("functions/data_processing.R")
source("functions/visualization.R")

# load files
ccle_counts <- read.csv("data/1_ccle_counts.csv", row.names = 1)
ccle_meta <- read.csv("data/1_ccle_meta.csv")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    div(align = "center",
      titlePanel("RNA-seq Cancer Genome Atlas Analysis"),
    ),
    
    # Sidebar with a slider input for number of bins 
    div(align = "center",
      selectInput("organ",
                    "Select primary site:",
                    choices = sort(unique(ccle_meta$Site_Primary)))
        ),
    
    div(align = "center",
      sliderInput("num_clusters", 
                  label = "Number of Clusters:",
                  min = 2, 
                  max = 10, 
                  value = 2,
                  step = 1),
      ),
    
    div(align = "center",
      tabPanel("Data Sampple",
               tableOutput("countsSample")),
    ),
    
    fluidRow(
      column(4, plotOutput("pie_1", height = "300px")),
      column(4, plotOutput("pie_2", height = "300px")),
      column(4, plotOutput("pie_3", height = "300px")),
    ),
    
    
    titlePanel("Principal Component Analysis"),
    
    fluidRow(
      column(6, plotOutput("pca_plot1", height = "300px")),
      column(6, plotOutput("pca_plot2", height = "300px")),
    ),
    
    div(align = "center",
      plotOutput("elbow_plot", height = "300px", width = "600px"),
    ),
    
    fluidRow(
      column(6, plotOutput("pca_genes_plot1", height = "300px")),
      column(6, plotOutput("pca_genes_plot2", height = "300px")),
    ),
    
    titlePanel("Hierarchical Clustering"),
    
    div(align = "center",
      plotOutput("hierarchical_plot", height = "900px"),
      plotOutput("pca_genes_plot_cluster", height = "300px",
                 width = "600px")
    ),
    
    titlePanel("Differential Expression Analysis"),
    
    div(align = "center",
        plotOutput("plot_2fold", height = "300px", width = "600px"),
        plotOutput("plot_2fold_genes", height = "300px", width = "600px"),
        plotOutput("plot_volcano", height = "300px", width = "600px")
        
    ),
    
    titlePanel("Top genes by differential expression"),
    
    div(align = "center",
        tableOutput("top_genes_table")
    ),
    
    titlePanel("Interpretation (done by GenAI)"),
    
    htmlOutput("query_response"),

    htmlOutput("links_output")
    
)

# Define server logic required to do a cluster analysis
server <- function(input, output) {
  
  # create the subsets
  subsets <- eventReactive(list(input$organ, input$num_clusters), {
    result <- filter_ccle_data_by_organ(ccle_counts, ccle_meta, 
                                        input$organ, input$num_clusters)
    return(result)
  })

  meta_subset <- reactive({subsets()[[1]]})
  counts_subset <- reactive({subsets()[[2]]})
  dds <- reactive({subsets()[[3]]})
  vsd <- reactive({subsets()[[4]]})
  pca_res <- reactive({subsets()[[5]]})
  explained_variance <- reactive({subsets()[[6]]})
  hc <- reactive({subsets()[[7]]})
  dds_low_counts<- reactive({subsets()[[8]]}) 
  res <- reactive({subsets()[[9]]}) 
  res_reordered <- reactive({subsets()[[10]]}) 
  top_genes <- reactive({subsets()[[11]]})
  k <- reactive({subsets()[[12]]})
  content_query <- reactive({subsets()[[13]]})
  query_citations <- reactive({subsets()[[14]]})
  sample_annotation_df <- reactive({subsets()[[15]]})
  tcga_code <- reactive({subsets()[[16]]})
  
  output$countsSample <- renderTable({
    req(subsets())
    head(meta_subset(),4)[1:3]
  })
  
  output$pie_1 <- renderPlot({
    req(subsets())
    counts <- table(meta_subset()$Pathology)
    pie(counts, labels = names(counts), main = paste0("Pathology in "
    , input$organ, "\nprimary site"))
  })
  
  output$pie_2 <- renderPlot({
    req(subsets())
    counts <- table(meta_subset()$Gender)
    pie(counts, labels = names(counts), main = paste0("Gender in "
                  , input$organ, "\nprimary site"))
  })
  
  output$pie_3 <- renderPlot({
    req(subsets())
    counts <- table(meta_subset()$tcga_code)
    pie(counts, labels = names(counts), main = paste0("TCGA code in "
                    , input$organ, "\nprimary site"))
  })
  
  #results_dds <- reactive({
  #  req(subsets())
  #  dds_subset_data(ccle_meta_subset(), 
  #                  ccle_counts_subset())
  #})
  
  #dds <- reactive({results_dds()[[1]]})
  #vsd <- reactive({results_dds()[[2]]})
  
  #observe({
  #print(head(data.frame(counts(dds()))[,1:5], 10))
  #})
  
  #observe({
  #  print(head(ccle_meta_subset(),4)[1:3])
  #})
  
  # PCA values
  #results_pca <- reactive({
  #  req(subsets())
  #  pca_subset(vsd())
  #})
  
  #pca_res <- reactive({results_pca()[[1]]})
  #explained_variance <- reactive({results_pca()[[2]]})
  #hc <- reactive({results_pca()[[3]]})
  
  #observe({
  #  print(explained_variance())
  #})
  
  
  pca_results <- reactive({
    req(subsets())
    pca_plots(explained_variance(), 1, 2, pca_res(), meta_subset(), hc(), 
              k(), "Pathology")
  })
  
  output$pca_plot1 <- renderPlot({
    req(subsets())
    pca_results()
  })
  
  pca_results2 <- reactive({
    req(subsets())
    pca_plots(explained_variance(), 3, 4, pca_res(), meta_subset(), hc(), 
              k(), "Pathology")
  })
  
  output$pca_plot2 <- renderPlot({
    req(pca_results2())
    pca_results2()
  })
  
  elbow_results <- reactive({
    req(subsets())
    elbow_plot(explained_variance())
  })
  
  output$elbow_plot <- renderPlot({
    req(elbow_results())
    elbow_results()
  })
  
  pca_genes_pc1_results <- reactive ({
    plot_genes_pca(explained_variance(), pca_res(), vsd(), meta_subset(), 1)
  })
  
  output$pca_genes_plot1 <- renderPlot({
    pca_genes_pc1_results()
  })
  
  pca_genes_pc2_results <- reactive ({
    plot_genes_pca(explained_variance(), pca_res(), vsd(), meta_subset(), 2)
  })
  
  output$pca_genes_plot2 <- renderPlot({
    pca_genes_pc2_results()
  })
  
  
  
  hiearchical_results <- reactive ({
    plot_hierarchical(vsd(), hc(), k(), sample_annotation_df(), tcga_code())
  })
  
  output$hierarchical_plot <- renderPlot({
    hiearchical_results()
  })
  
  # clusters
  pca_results_cluster <- reactive({
    req(subsets())
    pca_plots(explained_variance(), 1, 2, pca_res(), meta_subset(), hc(), 
              k(), "Cluster")
  })
  
  output$pca_genes_plot_cluster <- renderPlot({
    pca_results_cluster()
  })
  
  pca_genes_results <- reactive ({
    plot_genes_pca(explained_variance(), pca_res(), vsd(), meta_subset(), 1)
  })
  
  output$pca_genes_plot <- renderPlot({
    pca_genes_results()
  })
  
  
  #dds_cluster_result <- reactive({
  #  dea(vsd(), hc(), meta_subset(), 
        #counts_subset())
  #})
  
  # res <- reactive({dds_cluster_result()[[2]]})
  
  # plot 2fold differential expression histogram
  plot_2fold_results <- reactive ({
    plot_2fold(res())
  })
  
  output$plot_2fold <- renderPlot({
    plot_2fold_results()
  })
  
  # find the top genes with more differential expression
  # and order res by differential expression
  top_genes_2fold_result <- reactive({
    top_genes_2fold(res())
  })
  
  #top_genes <- reactive({top_genes_2fold_result()[[1]]})
  #res_reordered <- reactive({top_genes_2fold_result()[[2]]})
  
  
  # plot genes with more differential expression
  plot_2fold_genes_result <- reactive({
    req(top_genes_2fold_result())
    plot_2fold_genes(explained_variance(), pca_res(), res(), vsd(), 
                     top_genes(), meta_subset())
  })
  
  output$plot_2fold_genes <- renderPlot({
    req(plot_2fold_genes_result())
    plot_2fold_genes_result()
  })
  
  
  # plot volcano differential expression
  plot_differential_result <- reactive({
    plot_differential(top_genes(), res_reordered())
  })
  
  output$plot_volcano <- renderPlot({
    plot_differential_result()
  })
  
  output$query_response <- renderUI({
    
    # Clean brackets first
    content_cleaned <- gsub("\\[", "&#91;", content_query())
    content_cleaned <- gsub("\\]", "&#93;", content_cleaned)
    
    # Use HTML with no fragment limits
    HTML(markdown::markdownToHTML(text = content_cleaned, 
                                  fragment.only = TRUE,
                                  options = c("use_xhtml", 
                                              "smartypants", "base64_images")))
  })
  
  
  output$links_output <- renderUI({
    numbered_links <- paste0("[", 1:length(query_citations()), "] ", 
                             '<a href="', query_citations(), 
                             '" target="_blank">', query_citations(), '</a>')
    HTML(paste(numbered_links, collapse = '<br/>'))
  })
  
  # Define columns you want to show
  columns_to_show <- c("log2FoldChange", "pvalue", "padj")
  
  output$top_genes_table <-renderTable({
    subset_data <- res_reordered()[top_genes(), columns_to_show, 
                                 drop = FALSE]
    
    # Add gene names as a column
    subset_data$Gene <- rownames(subset_data)
    
    new_order <- c("Gene", "log2FoldChange", "pvalue", "padj")
    subset_data <- subset_data[,new_order]
    
    # Return the final data frame
    subset_data
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
