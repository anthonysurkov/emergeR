library(tidyverse)
library(here)
library(stats)
library(umap)
library(dbscan)
library(ggseqlogo)

setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("eda/edits_clustering.R")


# ============================================================================ #
# Parameters

infile_model <- "R255X_E488QD_modeling.csv"
infile_clean <- "R255X_E488QD_clean.csv"
infile_dir   <- "data"
outfile_dir  <- "data/motifs"

# Controls UMAP space. E.g. a value of 0.01 means all subsequent clusters will
# exist in the UMAP coordinate space defined by the subset of sequences with
# 1% editing or higher.
lower_edit_threshold <- 0.01

# Highest editing threshold the program will explore.
upper_edit_threshold <- 0.30

# Step size of the program. E.g. a value of 0.01 means the program will explore
# all values (lower_edit_threshold, lower_edit_threshold + 0.01, + 0.02, ...)
step_size <- 0.01

# Cutoff for reads. E.g. a value of 10 means the program will only evaluate
# sequences with 10 or more reads associated.
reads_threshold <- 10


# ============================================================================ #
# Main pipeline logic

main <- function() {
  model_data     <- load_model_data()
  umap_data_full <- generate_umap_representation(model_data)
  
  for (threshold in seq(lower_edit_threshold,
                        upper_edit_threshold,
                        by = step_size)) {
    
    # Subset UMAP data by editing cutoff; 
    umap_data <- filter_umap_data(umap_data_full, threshold)
    
    # Create directory;
    edit_dir  <- file.path(outfile_dir,
                           paste0("edit_", sprintf("%.2f", threshold)))
    dir.create(edit_dir, recursive=TRUE, showWarnings=FALSE)
  
    # Cluster UMAP;
    umap_clustered <- cluster_umap_data(umap_data)
    
    # Plot and save UMAP;
    plot_axes <- calculate_plot_axes(umap_clustered) 
    umap_plot <- plot_umap_clustered(
      umap_data = umap_clustered,
      x_limits  = plot_axes$x_limits,
      y_limits  = plot_axes$y_limits
    )
    plot_outfile  = paste0(
      "edit_", 
      sprintf("%.2f", threshold), 
      "_umap_clusters.png"
    )
    ggsave(
      filename = file.path(edit_dir, plot_outfile),
      plot     = umap_plot,
      width    = 6,
      height   = 4,
      dpi      = 300
    )
  
    # Summarize clusters and save data;
    cluster_info <- summarize_clusters(umap_clustered, edit_dir)
    
    # Build and save sequence logos;
    seq_logos <- plot_seq_logos(cluster_info, threshold)
    save_seq_logos(seq_logos, cluster_info$cluster_ids, edit_dir)
  }
}


# ============================================================================ #
# Helper functions

load_model_data <- function() {
  model_data <- read_csv(here(infile_dir, infile_model)) %>%
    filter(n >= reads_threshold) %>%
    filter(edit >= lower_edit_threshold)
  return(model_data)
}

generate_umap_representation <- function(model_data) {
  feature_data <- model_data %>% select(-N10, -GGA, -GAA, -n, -n_other, -edit)
  
  pca_result <- prcomp(feature_data, center=TRUE, scale.=TRUE, rank.=10)
  pca_matrix <- pca_result$x
  
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

filter_umap_data <- function(umap_data, threshold) {
   umap_data <- umap_data %>% filter(edit >= threshold)
   return(umap_data)
}

cluster_umap_data <- function(umap_data) {
   db <- dbscan(as.matrix(umap_data %>% select(X1, X2)), eps=0.5, minPts=5)
   umap_data <- umap_data %>% mutate(cluster = db$cluster)
   return(umap_data)
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

plot_umap_clustered <- function(umap_data, x_limits, y_limits) {
  umap_plot <- ggplot(umap_data, aes(x=X1, y=X2, color=as.factor(cluster))) +
    geom_point(size=1, alpha=0.5) +
    scale_color_viridis_d(name="Cluster") +
    xlim(x_limits) +
    ylim(y_limits) +
    labs(
      title="Local Sequence Landscape by Editing Percentage",
      subtitle=paste0("For editing fraction greater than or equal to ", 
                      target_edit),
      x="UMAP Dimension 1",
      y="UMAP Dimension 2"
    )
  
    return(umap_plot) 
}

summarize_clusters <- function(umap_data, outdir) {
  umap_split <- umap_data %>% filter(cluster > 0) %>% group_split(cluster)
  cluster_ids <- umap_data %>% filter(cluster > 0) %>% 
    distinct(cluster) %>% pull(cluster)
  
  lapply(seq_along(cluster_ids), function(i) {
    cid   <- cluster_ids[i]
    cdata <- umap_split[[i]]
    cdir  <- file.path(outdir, paste0("cluster_", cid))
    
    dir.create(cdir, showWarnings=FALSE, recursive=TRUE)
    write_csv(cdata, file.path(cdir, paste0("cluster_", cid, "_data.csv")))
  })
  
  return(list(
    cluster_ids    = cluster_ids, 
    clustered_data = umap_split 
  ))
}

plot_seq_logo <- function(seqs, k, threshold) {
  logo <- ggseqlogo(seqs, method="bits") +
    ggtitle(paste("Sequence Logo - Cluster", k, "(edit >=", threshold, ")")) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold"),
      panel.background = element_rect(fill="white", colour=NA),
      plot.background  = element_rect(fill="white", colour=NA)
    ) +
    labs(
      x="Position",
      y="Information (bits)"
    )
  
  return(logo)
}

plot_seq_logos <- function(cluster_info, threshold) {
  cluster_ids    <- cluster_info$cluster_ids
  clustered_data <- cluster_info$clustered_data
  
  seq_logos <- lapply(seq_along(cluster_ids), function(i) {
    seqs <- clustered_data[[i]] %>%
      pull(N10) %>%
      as.character()
    
    if (length(seqs) > 0 && all(nchar(seqs) == nchar(seqs[1]))) {
      plot_seq_logo(seqs, cluster_ids[i], threshold)
    } else{
      NULL
    }
  })
  
  return(seq_logos)
}

save_seq_logos <- function(seq_logos, cluster_ids, outdir) {
  for (i in seq_along(seq_logos)) {
    logo <- seq_logos[[i]]
    if (is.null(logo)) next
    
    this_outfile <- here(
      paste0(outdir, "/cluster_", cluster_ids[i]), 
      paste0("seqlogo_cluster_", cluster_ids[i], ".png"))
    ggsave(
      filename = this_outfile, 
      plot     = logo, 
      width    = 6,
      height   = 4, 
      dpi      = 300
    )
  }
}

main()
