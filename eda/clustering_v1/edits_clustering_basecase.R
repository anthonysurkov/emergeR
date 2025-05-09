library(tidyverse)
library(here)
library(umap)
library(stats)
library(ggseqlogo)
library(dbscan)


setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("modeling/linear_pairwise_glms.R")

infile_model <- "R255X_E488QD_modeling.csv"
infile_clean <- "R255X_E488QD_clean.csv"

lower_edit_threshold <- 0.01
upper_edit_threshold <- 0.10
step_size            <- 0.01
reads_threshold      <- 10


# ============================================================================ #


main <- function() {
  # Load data
  model_data <- read_csv(here("data", infile_model)) %>%
    filter(n >= reads_cutoff) %>%
    filter(edit >= edits_cutoff)
  
  # Generate UMAP in edit>0.01 UMAP space
  umap_data <- generate_umap_representation(model_data=model_data)
  
  # Filter editing according to current iteration
  umap_data %>% 
  
  # Cluster with dbscan
  message("Clustering started...")
  db <- dbscan(as.matrix(umap_data %>% select(X1, X2)), eps=0.5, minPts=5)
  
  axes <- calculate_plot_axes(umap_data)
  
  plot <- plot_umap_clusters(umap_data, clusters=db$cluster, 
    x_limits=axes$x_limits, y_limits=axes$y_limits
  )
  print(plot)
  
  umap_data <- umap_data %>% mutate(cluster = db$cluster)
  
  valid_clusters <- sort(unique(umap_data$cluster))
  valid_clusters <- valid_clusters[valid_clusters != 0]
  
  for (k in valid_clusters) {
    s <- umap_data %>%
      filter(cluster == k) %>%
      pull(N10) %>%
      as.character()
    
    if (length(s) > 0 && all(nchar(s) == nchar(s[1]))) {
      logo <- plot_seq_logo(s, k)
      print(logo)
    } else {
      message("Skipping cluster ", k, ": empty or invalid sequences.")
    }
  }
}


generate_umap_representation <- function(model_data) {
  feature_data <- model_data %>% select(-n_other, -N10, -GGA, -GAA, -n, -edit) 
  
  # PCA to reduce feature dimensions
  pca_result <- prcomp(feature_data, center=TRUE, scale.=TRUE, rank.=10)
  pca_matrix <- pca_result$x
  
  # Generate UMAP reduction from PCA results
  message("UMAP reduction started...")
  feature_umap <- umap(pca_matrix)
  edit_data    <- tibble(
    N10  = model_data$N10,
    X1   = feature_umap$layout[,1],
    X2   = feature_umap$layout[,2],
    edit = model_data$edit
  )
  
  return(edit_data)
}

calculate_plot_axes <- function(umap_data) {
  x_limits <- range(umap_data$X1)
  y_limits <- range(umap_data$X2)
  x_pad    <- 0.05 * diff(x_limits)
  y_pad    <- 0.05 * diff(y_limits)
  x_limits <- c(x_limits[1] - x_pad, x_limits[2] + x_pad)
  y_limits <- c(y_limits[1] - y_pad, y_limits[2] + y_pad)
  
  return(list(x_limits=x_limits, y_limits=y_limits))
}

plot_seq_logo <- function(s, k) {
  logo <- ggseqlogo(s, method="bits") +
    ggtitle(paste("Sequence Logo - Cluster", k, "for editing >", target_edit)) +
    theme_minimal() +
    labs(
      x="Position",
      y="Information (bits)"
    )
  
  return(logo)
}

plot_umap_clusters <- function(umap_data, clusters, x_limits, y_limits) {
  umap_data <- umap_data %>% mutate(cluster = clusters)
  
  edit_plot <- ggplot(umap_data, aes(x=X1, y=X2, color=as.factor(cluster))) +
    geom_point(size=1, alpha=0.5) +
    scale_color_viridis_d() +
    xlim(x_limits) +
    ylim(y_limits) +
    labs(
      title="Local Sequence Landscape by Editing Percentage",
      subtitle=paste0("edit > ", target_edit),
      x="UMAP Dimension 1",
      y="UMAP Dimension 2"
    )
  
  return(edit_plot)
}


main()
