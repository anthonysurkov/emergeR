library(tidyverse)
library(ggprism)
library(scales)
library(here)
source("pipeline/5_mca/mca_vis.R")

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("pipeline/7_hdbscan/tsne_clusters_vis.R")

indir  <- "data"
outdir <- "data"
infile <- "R255X_E488QD_7_tSNEclustered.csv"
outpdf <- "R255X_E488QD_7_report.pdf"


# ============================================================================ #
# Generates a report visualizing a clustered t-SNE data set.

# Visualizes editing in MCA space
make_editing_vis <- function(X_emb) {
  p <- ggplot(X_emb, aes(tsne1, tsne2, colour = map)) +
    geom_point(size = 1.0, alpha = 0.40) +
    scale_colour_gradientn(
      colours = c("#b0f3ff", "#003e4c", "#ff0000"),
      name    = "Editing Rate",
      limits  = c(0, 0.4),
      oob     = scales::squish
    ) +
    guides(
      colour = guide_legend(position="right")
    ) +
    theme_prism(
      palette = "winter_bright", 
      base_size = 14,
      base_family = "sans",
      base_fontface = "bold",
      base_line_size = 1,
      base_rect_size = 1,
      axis_text_angle = 0,
      border = FALSE
    ) +
    theme(
      plot.subtitle = element_text(size=12, margin=margin(t=-10)),
      legend.title = element_text(size=12),
      plot.title = element_text(size=16),
      legend.position = c(0.05, 0.95),
      legend.margin = margin(2, 2, 2, 2),
      legend.box.margin = margin(0, 0, 0, 0)
    ) +
    labs(
      x = "MC1",
      y = "MC2",
      title = "MCA of Editing-Capable Guides",
      subtitle = "99.9th percentile Maximum A Posteriori Editors"
    ) +
    scale_y_continuous(
      limits = c(-21, 21),
      guide = "prism_offset_minor"
    ) +
    scale_x_continuous(
      limits = c(-21, 21),
      guide = "prism_offset_minor"
    ) +
    coord_fixed()
}

# Visualizes MCA clustering via colored MC1/2 plot
make_cluster_vis <- function(X_emb) {
  p <- ggplot(X_emb, aes(tsne1, tsne2, colour = cluster)) +
    geom_point(size = 1.0, alpha = 0.40) +
    scale_colour_discrete(
      name = "Cluster ID",
    ) +
    guides(
      colour = guide_legend(position="right")
    ) +
    theme_prism(
      palette = "winter_bright", 
      base_size = 14,
      base_family = "sans",
      base_fontface = "bold",
      base_line_size = 1,
      base_rect_size = 1,
      axis_text_angle = 0,
      border = FALSE
    ) +
    theme(
      plot.subtitle = element_text(size=12, margin=margin(t=-10)),
      legend.title = element_text(size=12),
      plot.title = element_text(size=16),
      legend.position = c(0.05, 0.95),
      legend.margin = margin(2, 2, 2, 2),
      legend.box.margin = margin(0, 0, 0, 0)
    ) +
    labs(
      x = "MC1",
      y = "MC2",
      title = "MCA of Editing-Capable Guides",
      subtitle = "99.9th percentile Maximum A Posteriori Editors"
    ) +
    scale_y_continuous(
      limits = c(-21, 21),
      guide = "prism_offset_minor"
    ) +
    scale_x_continuous(
      limits = c(-21, 21),
      guide = "prism_offset_minor"
    ) +
    coord_fixed()
}

# Prints pdf of clustered MC1/2 plot along with sequence logo of all clusters.
tsne_tri_vis <- function(X, outpdf) {
  seq_logos   <- make_seq_logos(X) # from mca_vis.R
  
  cluster_vis <- make_cluster_vis(X)
  editing_vis <- make_editing_vis(X)
  combined    <- c(
    list(cluster_vis, editing_vis), seq_logos
  )
  
  cairo_pdf(here(outdir, outpdf), width = 8, height = 10, onefile = TRUE)
  walk(combined, print)
  dev.off()
}

main <- function() {
  X <- read_csv(here(indir, infile))
  print(head(X))
  tsne_tri_vis(X, outpdf = outpdf)
}

main()
