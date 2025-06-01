library(tidyverse)
library(here)
library(stringdist)
library(igraph)
library(Matrix)
library(pheatmap)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("integration/hamming_distances.R")

indir <- "data"
infile <- "R255X_E488QD_clean.csv"
set.seed(0)


# ============================================================================ #
# Takes a random sample of good editors, computes pairwise distances, and plots
# as a heat map

seqs <- read_csv(here(indir, infile)) %>% 
  filter(edit > 0.10) %>%
  slice_sample(n = 2000) %>%
  pull(N10)

message("Now starting distance matrix calculations...")
dmat <- stringdistmatrix(seqs, seqs, method="hamming", useBytes=TRUE)
smat <- 10 - dmat

thr <- which(smat >= 4 & upper.tri(smat), arr.ind=TRUE)
sparse_sim <- as.matrix(sparseMatrix(
  i = thr[, 1],
  j = thr[, 2],
  x = smat[thr],
  dims = c(length(seqs), length(seqs)),
  symmetric = TRUE
))

g <- graph_from_adjacency_matrix(sparse_sim, mode="undirected", weighted=TRUE)
comm <- cluster_louvain(g)

order <- order(membership(comm))
ord_mat <- as.matrix(sparse_sim)[order, order]

message("Processing heatmap...")
heat <- pheatmap(
  ord_mat,
  cluster_cols=TRUE,
  na_col="white",
  color=colorRampPalette(c("navy","white","firebrick3"))(50),
  show_rownames=FALSE,
  show_colnames=FALSE,
  main = "Sequence similarity (similar=red, dissimilar=blue)"
)
print(heat)
