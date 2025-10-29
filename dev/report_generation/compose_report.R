library(emergeR)
library(ggprism)
library(tidyverse)

compose_report <- function(X, name, dataset) {
  tsne_edits <- emergeR::plot_tsne_edits(X = X)
  tsne_clust <- emergeR::plot_tsne_clusters(X = X)
  seq_logos  <- emergeR::plot_motif_logos(X = X)
  
  rmarkdown::render(
    input       = "report_template.Rmd",
    output_file = paste0("report_", name, ".html"),
    params  = list(
      title        = paste("Analysis:", name),
      date         = Sys.Date(),
      dataset_name = dataset,
      tsne_edits   = tsne_edits,
      tsne_clust   = tsne_clust,
      seq_logos    = seq_logos
    ),
    envir = new.env(parent = globalenv())
  )
}

run_pipeline <- function(X, template, flanks) {
  
  message("Cleaning...")
  X_clean <- emergeR::clean(
    X,
    template = template,
    flanks = flanks
  )
  write_csv(X_clean, "R270X_MOD_clean.csv")
  
  message("Computing statistics...")
  X_stats <- emergeR::append_stats(X_clean)
  write_csv("R270X_MOD_stats.csv")
  
  #message("Computing MCA...")
  #X_mca   <- emergeR::append_mca(X_stats, top_quantile=0.999)
  
  #message("Computing t-SNE...")
  #X_tsne  <- emergeR::append_tsne(X_mca, verbose=TRUE)
  
  #message("Computing clusters...")
  #X_clust <- emergeR::append_clusters(X_tsne, minPts=10)
  
  #message("Composing report...")
  #compose_report(
  #  X_clust,
  #  name = "R255X",
  #  dataset = X
  #)
  
  message("Done!")
}

X_in <- "R270X_mod.fastq"
template = paste0(
  "GGGAATGGCGCTTCCCGGCAGGAAGTGAAAAGCTGAGGCCGACTTTCTTTCTTTCGTCGGCCTCA",
  "NNNNNNNNNNTCCTGCCGGGTTGTAGGGTC"
)
flanks = c("GGAAGT", "AAGCTG", "GCCTCA", "TCCTGC")
run_pipeline(X_in, template, flanks)