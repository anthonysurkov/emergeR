library(tidyverse)
library(dbscan)
library(stringdist)
library(here)
library(ggprism)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR")
here::i_am("eda/edits/bayes/MAP_dbscan.R")

infile  <- "R255X_E488QD_onehot_dummy_demeaned.csv"
outfile <- "MAP_clustered.csv"
indir   <- "data"
outdir  <- "eda/edits/bayes/data"

# hyperparams
pcs_to_use <- 2

# ============================================================================ #
# Clusters top sequences (by MAP) with dbscan for a conservative approach to
# high-editing motifs. Visualizes clusters in principal component space and as
# sequence logos as qualitative check.

pca_vis <- function(
    model_data, 
    minPts        = 20, 
    map_threshold = 0.05, 
    pcs           = c("PC1", "PC2"),
    center        = FALSE, 
    scale         = FALSE, 
    save_csv      = NULL
) {
  pca_mat <- model_data %>%
    select(
      -N10, -n, -k,
      -map, -upper_cred, -lower_cred,
      -mle, -upper_ci, -lower_ci,
      -alpha, -beta
    ) %>%
    as.matrix()
  
  pca_fit <- prcomp(pca_mat, center = center, scale = scale)
  scores <- as_tibble(pca_fit$x) %>%
    mutate(N10 = model_data$N10) %>%
    left_join(model_data, by="N10") %>%
    filter(map > map_threshold)
  
  # scree plot
  var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  cum_explained <- cumsum(var_explained)
  elbow_df      <- tibble(
    PC = seq_along(var_explained),
    Variance = var_explained,
    Cumulative = cum_explained
  )
  elbow_plot    <- ggplot(elbow_df, aes(x=PC, y=Variance)) +
    geom_point(size = 2) +
    geom_line() +
    scale_x_continuous(breaks = seq_along(var_explained)) +
    labs(
      title = "PCA Variance Explained per Principal Component",
      x = "Principal Component",
      y = "Proportion of Variance Explained"
    ) +
    theme_minimal()
  print(elbow_plot)
  
  pc_coords <- scores[, pcs]
  hd        <- hdbscan(pc_coords, minPts=minPts)
  
  scores <- scores %>%
    mutate(
      hd_cluster_raw = hd$cluster,
      hd_cluster     = factor(
        ifelse(hd$cluster == 0, "Noise", as.character(hd$cluster)),
               levels = c("Noise", sort(unique(hd$cluster[hd$cluster != 0])))
      )
    )
  
  p_clusters <- ggplot(
      scores, 
      aes_string(x=pcs[1], y=pcs[2], colour="hd_cluster")
    ) +
    geom_point(alpha=0.7, size=1.6) +
    coord_fixed() +
    theme_prism(
      palette="winter_bright",
      base_size=14,
      base_family="sans",
      base_fontface="bold",
      base_line_size=1,
      base_rect_size=1,
      axis_text_angle=0,
      border=FALSE
    ) +
    scale_color_manual(
      values = c(
        "Noise" = "#D4D4D4",  # special color for noise
        setNames(
          viridis::viridis(length(unique(scores$hd_cluster)) - 1),
          setdiff(unique(scores$hd_cluster), "Noise"))
      ),
      aesthetics = "colour",
    ) + 
    theme(
      plot.subtitle = element_text(size=12, margin=margin(t=-10)),
      legend.title  = element_text(size=12),
      plot.title    = element_text(size = 16),
      legend.margin = margin(2,2,2,2),
      legend.box.margin = margin(0,0,0,0)
    ) +
    labs(
      x = pcs[1],
      y = pcs[2],
      title = "PCA of Editing-Capable Guides",
      subtitle = paste0("(MAP editing > ", signif(map_threshold, 2), ")")
    ) +
    scale_y_continuous(limits = c(-0.77, 1.2), guide = "prism_offset_minor") +
    scale_x_continuous(limits = c(-0.77, 1.2), guide = "prism_offset_minor")
  
  print(p_clusters)
}

main <- function() {
  data <- read_csv(here(indir, infile))
  print(head(data))
  
  top_data <- data %>%
    filter(map > quantile(map, 0.999, na.rm=TRUE))
  top_seqs <- top_data %>%
    pull(N10)
  top_edit <- top_data %>%
    pull(map)
  
  pcs = paste0("PC", 1:pcs_to_use)
  pca_vis(data, pcs=pcs, center=TRUE)
}

main()
