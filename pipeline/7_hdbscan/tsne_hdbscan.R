library(tidyverse)
library(dbscan)
library(here)
source("pipeline/5_mca/mca_vis.R")

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("pipeline/7_hdbscan/tsne_hdbscan.R")

indir   <- "data"
infile  <- "R255X_E488QD_6_tSNE.csv"
outfile <- "R255X_E488QD_7_tSNEclustered.csv"


# ============================================================================ #
# Cluster t-SNE embedded MCA data.

main <- function() {
  X      <- read_csv(here(indir, infile))
  X_tsne <- X %>%
    select(tsne1, tsne2)
  
  minPts <- 10
  hd     <- hdbscan(X_tsne, minPts = minPts)
  
  X_clustered <- X %>%
    mutate(
      cluster_raw = hd$cluster,
      cluster = factor(
        ifelse(hd$cluster == 0, "Noise", as.character(hd$cluster)),
        levels = c("Noise", sort(unique(hd$cluster[hd$cluster != 0])))
      ) 
    )
  
  write_csv(X_clustered, here(indir, outfile))
}

main()
