library(tidyverse)
library(ggseqlogo)
library(ggprism)
library(scales)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/mca/mca_vis.R")

indir    <- "data"
outdir   <- "data"
infile_a <- "R255X_E488QD_3a_topeditors.csv"
outpdf   <- "R255X_E488QD_report.pdf"


# ============================================================================ #
# Generates a report visualizing a clustered MCA data set.

# Visualizes MCA clusters via a list of sequence logos
make_seq_logos <- function(X_emb) {
  seq_grouped <- X_emb %>%
    group_by(cluster) %>%
    group_split(.keep = FALSE) %>%
    { .[ order(map_int(., nrow), decreasing = TRUE) ]}
  
  logo_list <- imap(seq_grouped, ~ {
    ggseqlogo(
      pull(.x, N10),
      method = "bits",
      seq_type = "dna",
      col_scheme = "nucleotide"
    ) +
      ggtitle(
        paste0("99.9th Percentile Editors (MAP), Cluster ", .y, 
               " (n = ", nrow(.x), ")"
        )
      ) +
      guides(colour = "none") +
      theme_prism(
        palette         = "winter_bright",
        base_size       = 14,
        base_family     = "sans",
        base_fontface   = "bold",
        base_line_size  = 1,
        base_rect_size  = 1,
        axis_text_angle = 0,
        border          = FALSE
      ) +
      theme(
        plot.title   = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_text(face = "bold"),
        axis.text.y  = element_text(face = "bold")
      )
    
  })
}

# Visualizes editing in MCA space
make_MCA_editing_plot <- function(X_emb) {
  p <- ggplot(X_emb, aes(Axis1, Axis2, colour = map)) +
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
      limits = c(-0.87, 0.87),
      guide = "prism_offset_minor"
    ) +
    scale_x_continuous(
      limits = c(-0.87, 0.87),
      guide = "prism_offset_minor"
    ) +
    coord_fixed()
}

# Visualizes MCA clustering via colored MC1/2 plot
make_MCA_clustered_plot <- function(X_emb) {
  p <- ggplot(X_emb, aes(Axis1, Axis2, colour = cluster)) +
    geom_point(size = 1.0, alpha = 0.40) +
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
      limits = c(-0.87, 0.87),
      guide = "prism_offset_minor"
    ) +
    scale_x_continuous(
      limits = c(-0.87, 0.87),
      guide = "prism_offset_minor"
    ) +
    coord_fixed()
}

# Prints pdf of clustered MC1/2 plot along with sequence logo of all clusters.
tri_vis <- function(X) {
  editing_plot   <- make_MCA_editing_plot(X)
  clustered_plot <- make_MCA_clustered_plot(X) 
  seq_logos      <- make_seq_logos(X)
  combined       <- c(
    list(editing_plot, clustered_plot), 
    seq_logos
  ) 
  
  cairo_pdf(here(outdir, outpdf), width = 8, height = 10, onefile = TRUE)
  walk(combined, print)
  dev.off()
}

main <- function() {
  X <- read_csv(here(indir, infile))
  tri_vis(X)
}

main()
