library(tidyverse)
library(here)
library(dbscan)
library(pracma)

setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("eda/hdbscan_diagnostics.R")


# ============================================================================ #
# Parameters

infile <- "R255X_E488QD_umap.csv"
indir  <- "data"
outdir <- "data/hdbscan"

minPts_min     <- 2
minPts_max     <- 100
edit_threshold <- 0.05


# ============================================================================ #

main <- function() {
  set.seed(0)
  
  umap_data     <- read_csv(here(indir, infile))
  umap_features <- umap_data %>% 
    filter(edit > edit_threshold) %>% 
    select(X1, X2)
  
  message("Now building hdbscan sweep plot")
  sweep_results <- hdbscan_sweep_plot(umap_features)
  
  x <- sweep_results$minPts
  noise <- sweep_results$noise_prop
  prob  <- sweep_results$mean_membership_prob
  
  maxima_idx <- findpeaks(prob, sortstr=TRUE)[,2]
  minima_idx <- findpeaks(-noise, sortstr=TRUE)[,2]
  
  peak_pts <- x[maxima_idx]
  dip_pts  <- x[minima_idx]
  
  print(sweep_results[maxima_idx[1:10], ])
  print(sweep_results[minima_idx[1:10], ])
}


hdbscan_sweep_plot <- function(umap_data, minPts_range = minPts_min:minPts_max){
  umap_matrix <- as.matrix(umap_data)
  
  sweep_results <- map_dfr(minPts_range, function(minPts) {
    message(paste0("Now starting minPts=", minPts))
    
    hdb <- hdbscan(umap_matrix, minPts = minPts)
    
    cluster_labels <- hdb$cluster
    membership_prob <- hdb$membership_prob
    
    n_clusters <- length(unique(cluster_labels[cluster_labels != 0]))
    noise_prop <- mean(cluster_labels == 0)
    mean_prob <- mean(membership_prob[cluster_labels != 0])
    
    tibble(
      minPts = minPts,
      n_clusters = n_clusters,
      noise_prop = noise_prop,
      mean_membership_prob = mean_prob
    )
  })
  
  # Plot number of clusters
  p1 <- ggplot(sweep_results, aes(x = minPts, y = n_clusters)) +
    geom_line(color = "steelblue", size = 1) +
    labs(
      title = "Number of Clusters vs minPts", 
      y = "# Clusters", 
      x = "minPts"
    )
  
  # Plot noise proportion
  p2 <- ggplot(sweep_results, aes(x = minPts, y = noise_prop)) +
    geom_line(color = "darkred", size = 1) +
    labs(
      title = "Noise Proportion vs minPts", 
      y = "Proportion Noise", 
      x = "minPts"
    )
  
  # Plot mean membership probability
  p3 <- ggplot(sweep_results, aes(x = minPts, y = mean_membership_prob)) +
    geom_line(color = "darkgreen", size = 1) +
    labs(
      title = "Mean Membership Probability vs minPts", 
      y = "Mean Prob (non-noise)", 
      x = "minPts"
    )
  
  print(p1)
  print(p2)
  print(p3)
  
  return(sweep_results)
}


main()
