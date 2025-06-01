library(tidyverse)
library(purrr)
library(here)
library(cowplot)
library(stringdist)
library(ggseqlogo)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/edits/frequentist/top_pwms/mle_top_pwms.R")

infile_clean <- "R255X_E488QD_clean.csv"
indir        <- "data"

outfile      <- "seed_logos.pdf"
outdir       <- "eda/seed_pwms"

# Define editing threshold above which sequences are considered for use as seeds
edit_min <- 0.30

# Define reads threshold above which sequences are considered for use as seeds
read_min <- 10

# Define minimum number of sequences present in each motif
motif_min <- 10


# ============================================================================ #
# Clustering top sequences (by MLE) to find well-editing motifs

main <- function() {
  clean_data <- read_csv(here(indir, infile_clean))
  clean_seqs <- clean_data %>%
    filter(edit >= edit_min, n >= read_min) %>%
    pull(N10)
  message("Number of sequences considered: ", length(clean_seqs))
  
  dmat   <- stringdistmatrix(clean_seqs, method="hamming")
  hc     <- hclust(dmat, method="average")
  groups <- cutree(hc, h=5)
  
  seq_df <- tibble(seq=clean_seqs, grp=groups) %>%
    group_by(grp) %>% 
    filter(n() >= motif_min)
  message("Number of motifs: ", length(unique(seq_df$grp)))
  
  seq_grouped <- group_split(seq_df, .keep=FALSE) %>%
    { .[ order(map_int(., nrow), decreasing = TRUE) ] }
  
  logo_list <-imap(seq_grouped, ~ {
    ggseqlogo(pull(.x, seq), method="bits") +
      ggtitle(paste0("R255X E488QD Motif AS", .y, "\n( >30% editing, n = ", 
                    nrow(.x), " )")) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size=10))
  })
  
  n_logos  <- length(logo_list)
  halfway  <- ceiling(n_logos / 2)
  logo_pg1 <- logo_list[1:halfway]
  logo_pg2 <- logo_list[(halfway+1):n_logos]
  
  cairo_pdf(here(outdir, outfile), width=8, height=10, onefile=TRUE)
  print(plot_grid(plotlist=logo_pg1, ncol=2))
  print(plot_grid(plotlist=logo_pg2, ncol=2))
  dev.off()
}


main()
