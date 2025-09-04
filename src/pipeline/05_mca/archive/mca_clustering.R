library(tidyverse)
library(dbscan)
library(purrr)
library(ggseqlogo)
library(ggprism)
library(scales)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/mca/mca_clustering.R")

indir    <- "data"
infile_a <- "R255X_E488QD_3a_topeditors.csv"
outfile  <- "R255X_E488QD_4_topclustered.csv"
outpdf   <- "R255X_E488QD_4_clustervis.pdf"


# ============================================================================ #
# Clusters editors according to multiple correlate components.

main <- function() {
  X_embed <- read_csv(here(indir, infile_a))
   
  n         <- nrow(X_embed)
  m         <- 10
  cols      <- paste0("Axis", 1:m)
  minPts    <- 15
  
  X_hdbscan <- as.matrix(X_embed[, cols, drop = FALSE])
  hd        <- hdbscan(X_hdbscan, minPts = minPts)
  nn        <- hd$cluster > 0
  
  X <- X_embed %>%
    mutate(
      cluster         = hd$cluster,
      membership_prob = hd$membership_prob,
      outlier_score   = hd$outlier_scores 
    )
  print(head(X))

  write_csv(X, here(indir, outfile))
}

main()
