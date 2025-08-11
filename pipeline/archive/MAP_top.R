library(tidyverse)
library(cowplot)
library(stringdist)
library(ggseqlogo)
library(ggprism)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR")
here::i_am("eda/edits/bayes/archive/MAP_top.R")

infile  <- "R255X_E488QD_CV.csv"
outcsv  <- "MAP_top_editors.csv"
outpdf  <- "MAP_motifs.pdf"
indir   <- "data"
outdir  <- "eda/edits/bayes"
outfile <- "top_editors_by_cluster.csv"

# Motif hyperparameters
max_diff  <- 6
motif_min <- 10


# ============================================================================ #
# Clustering top sequences (by MAP) for a conservative look at high-editing 
# motifs.

main <- function() {
  data <- read_csv(here(indir, infile))
  
  # Find 99.9th percentile by MAP.
  top_data <- data %>%
    filter(map > quantile(map, 0.999, na.rm=TRUE))
  print(top_data)
  write_csv(top_data, here(outdir, outcsv))
  message("Number of sequences matching 99.9th percentile: ", nrow(top_data))
  
  # Print 99.9th percentile.
  edit_long <- top_data %>%
    pivot_longer(cols=c(map, mle), names_to="estimator", values_to="value")
  both_plot <- ggplot(edit_long, aes(x=value, fill=estimator)) +
    geom_histogram(position="identity", alpha=0.5, bins=30) +
    labs(
      title="Distribution of MAP and MLE over the top 0.1% of MAP editors", 
      x="Estimate", 
      y="Count") +
    theme_minimal()
  print(both_plot)
  
  # Cluster 99.9th percentile.
  top_seqs <- top_data %>% pull(N10)
  top_edit <- top_data %>% pull(map)
  dmat     <- stringdistmatrix(top_seqs, method="hamming")
  hc       <- hclust(dmat, method="average")
  groups   <- cutree(hc, h=max_diff)
  
  seq_df   <- tibble(seq=top_seqs, edit=top_edit, grp=groups) %>%
    group_by(grp) %>%
    filter(n() >= motif_min)
  message("Number of motifs: ", length(unique(seq_df$grp)))
  

  
  seq_grouped <- group_split(seq_df, .keep=FALSE) %>%
    { .[ order(map_int(., nrow), decreasing = TRUE) ] }
  avg_map <- map_dbl(seq_grouped, ~ mean(.x$edit))
  names(avg_map) <- paste0("MAP AS", seq_along(avg_map))
  print(avg_map)
  
  # Fancy sequence logos.
  logo_list <-imap(seq_grouped, ~ {
    
    ggseqlogo(
      pull(.x, seq), 
      method="bits",
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
      axis.text.x  = element_blank(),
      axis.text.y  = element_text(face = "bold")
    )
    
  })
  
  n_logos  <- length(logo_list)
  halfway  <- ceiling(n_logos / 2)
  logo_pg1 <- logo_list[1:halfway]
  logo_pg2 <- logo_list[(halfway+1):n_logos]
  
  cairo_pdf(here(outdir, outpdf), width=8, height=10, onefile=TRUE)
  walk(logo_list, print)
  dev.off()
  
  seq_df_for_xlsx <- seq_df %>%
    arrange(desc(grp))
  write_csv(seq_df_for_xlsx, here(outdir, outfile))
}

main()
