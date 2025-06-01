library(tidyverse)
library(stringr)
library(purrr)
library(here)
library(ggseqlogo)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/edits/tcc_motifs/tcc_6to8/tcc_6to8_slidingkmer_logos.R")

infile <- "R255X_E488QD_clean.csv"
indir  <- "data"
outdir <- "D:/ngs storage/Natalie/R255X E488QD/emergeR/eda/edits/tcc_motifs/tcc_6to8"

# Animation hyperparameters
window_size <- 50
step_size   <- 5
img_width   <- 8
img_height  <- 4
img_dpi     <- 100


# ============================================================================ #
# Intuition-building for important submotifs given TCC @ positions 6~8

main <- function() {
  data <- read_csv(here(indir, infile))
  
  tcc_6to8 <- data %>%
    mutate(
      tcc_6to8 = (str_sub(N10, 6, 8) == "TCC")
    ) %>%
    filter(tcc_6to8 == TRUE, map > 0.05)
  
  total <- nrow(tcc_6to8)
  start_indices <- seq(1, total - window_size + 1, by = step_size)
  
  for (i in seq_along(start_indices)) {
    start <- start_indices[i]
    end   <- start + window_size - 1
    window_df <- tcc_6to8[start:end, ]
    logo_seqs <- window_df$N10
  
    p <- ggseqlogo(logo_seqs, method="bits") +
      ggtitle(paste0("Sequences ", start, "-", end, " (ranked by MAP)")) +
      theme_classic()
    
    outfile <- file.path(outdir, sprintf("logo_%03d_%03d.png", start, end))
    ggsave(outfile, plot=p, width=img_width, height=img_height, dpi=img_dpi)
  }  
  
}

main()
