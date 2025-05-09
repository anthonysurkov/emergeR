library(tidyverse)
library(here)
library(dbscan)
library(ggseqlogo)
library(patchwork)

setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("eda/hdbscan_main.R")


# ============================================================================ #
# Parameters

infile <- "R255X_E488QD_umap.csv"
indir  <- "data"
outdir <- "data"

edit_threshold <- 0.05
read_threshold <- 10
my_minPts      <- 8

# ============================================================================ #

main <- function() {
  set.seed(0)
  
  umap_data <- read_csv(here(indir, infile))
  umap_filtered <- umap_data %>% filter(edit > edit_threshold)
  cat("Umap filtered nrow: ", nrow(umap_filtered), '\n')
  
  umap_features <- umap_filtered %>% select(X1, X2)
  
  hdb <- hdbscan(as.matrix(umap_features), minPts=my_minPts)
  #print(hdb)
  #plot(
  #  hdb$hc, labels=FALSE, 
  #  main=paste0("HDBSCAN condensed tree, minPts=", my_minPts)
  #)
  
  clusters_df <- tibble(
    cluster = hdb$cluster,
    prob    = hdb$membership_prob
  )
  
  # Identify 4 earliest/root/main families of sequences
  tree4 <- cutree(hdb$hc, k=4)
  families4 <- umap_filtered %>%
    mutate(family = factor(tree4)) %>%
    group_by(family) %>%
    summarise(seqs = list(N10), .groups="drop")
  families4 %>% mutate(size=map_int(seqs, length))
  logos4 <- families4 %>%
    pmap(function(family, seqs) {
      ggseqlogo(seqs, method="bits") +
        ggtitle(paste("Main family:", family))
    })
  #print(wrap_plots(logos4, ncol=2))
  
  # Identify 4 most stable families of sequences
  top_clusters <- clusters_df %>%
    filter(cluster != 0) %>%
    count(cluster, sort=TRUE) %>%
    slice_head(n=4) %>%
    pull(cluster)
  
  umap_clustered <- umap_filtered %>% bind_cols(clusters_df) %>%
    group_by(cluster) %>%
    mutate(cluster_median = median(edit, na.rm=TRUE)) %>%
    ungroup() %>%
    arrange(desc(edit))
  write_csv(umap_clustered, here(outdir, "R255X_E488QD_clustered.csv"))
  
  seqs_by_cluster <- umap_clustered %>%
    filter(cluster %in% top_clusters) %>%
    group_by(cluster) %>%
    summarise(seqs=list(N10), .groups="drop")
  
  logos <- seqs_by_cluster %>%
    pmap(function(cluster, seqs) {
      ggseqlogo(seqs, method="bits") +
        ggtitle(paste("Family", cluster))
    })
    
  #print(wrap_plots(logos, ncol=2))
  
  # Identify top100 families 
  top100 <- read_csv(here("nd/data", "R255X_E488QD_top_100.csv"))
  top100_clustered <- umap_clustered %>% 
    filter(N10 %in% top100$N10) %>%
    arrange(desc(edit))
  print(top100_clustered, n=Inf)
  
  print(umap_clustered %>% 
          distinct(cluster,cluster_median) %>% 
          arrange(desc(cluster_median)), n=100)
  
  old_top_10 = c("TTTGTTCCAC", "AGAGTTCCGC", "TTCTTTCGAC",
                 "TTTTCCTATT", "AGAGTTCCAA", "ATTTCCGATT",
                 "TGAATTCCGC", "AGGGTTCCAT", "TCATTCCGGT",
                 "CGAATTCCGC")
  
  print(umap_clustered %>% filter(N10 %in% old_top_10))
  
  print(
    umap_clustered %>% filter(cluster==68)
  )
  print(
    umap_clustered %>% filter(cluster==86)
  )
}




main()