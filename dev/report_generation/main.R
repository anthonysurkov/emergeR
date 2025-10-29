library(targets)
library(tarchetypes)

setwd("D:/ngs storage/Natalie/emergeR/")


# ============================================================================ #

# Take input arguments from exec.sh
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < EXPECTED_ARGS) {
  stop("Usage: Rscript main.R <working_directory> <path_to_file>")
}
fastq_path <- args[1]
if (!file.exists(fastq_path)) {
  stop(paste("File does not exist:", fastq_path))
}
outdir <- args[2]
if (!dir.exists(outdir)) {
  stop(paste("Directory does not exist:", outdir))
}

# Activate pipeline scripts
lapply(list.files("pipeline", full.name=TRUE, source))
set.seed(1)

# Set up dependencies
tar_option_set(packages = c(
  "tidyverse",
  "ShortRead",
  "binom",
  "ade4",
  "dbscan",
  "Rtsne",
  "ggseqlogo",
))

# Execute pipeline
list(
  # 01: Read ht_stream preprocessed data, clean into tidy format
  tar_target(data_clean, emerge_clean(fastq_path)),
  # 02: Append statistics, calculate EB editing rate with k-fold CV, save data
  tar_target(data02_processed, emerge_stats(data01_clean)),
  tar_target(file02_processed, {
    outfile <- "mid/processed.csv"
    dir.create(dirname(outfile), showWarnings=FALSE, recursive=TRUE)
    readr::write_csv(data02_processed, outfile)
    outfile
  }, format = "file"),
  # 03: Transform data with MCA, save data
  tar_target(data03_mca, emerge_mca(data02_processed)),
  tar_target(file03_mca, {
    outfile <- "mid/mca_data.csv"
    readr::write_csv(data03_mca, outfile)
    outfile
  }, format = "file"),
  # 04: Reduce dimensionality with t-SNE, save data and 2D plot colored by
  #     editing
  tar_target(data04_tsne, emerge_tsne(data03_mca)),
  tar_target(plot04_tsne, emerge_tsne_plot(data02_processed, data04_tsne)),
  tar_target(file04_tsne, {
    outfile <- "img/tsne_unclustered.png"
    dir.create(dirname(outfile), showWarnings=FALSE, recursive=TRUE)
    ggplot2::ggsave(outfile, plot04_tsne, width=6, height=4)
    outfile
  }, format = "file"),
  # 05: Cluster reduced data with HDBSCAN, save data and 2D plot colored by
  #     cluster
  tar_target(data05_hdbscan, emerge_hdbscan(data04_tsne)),
  tar_target(plot05_hdbscan, emerge_hdbscan_plot(data04_tsne, data05_hdbscan)),
  tar_target(file05_hdbscan, {
    outfile <- "img/tsne_clustered.png"
    ggplot2::ggsave(outfile, plot05_hdbscan, width=6, height=4)
    outfile
  }, format = "file"),
  # 06: Identify statistically significant motifs, generate sequence logos
  tar_target(data06_motifs, emerge_motif_stats(data05_hdbscan)),
  tar_target(plots06_motifs, emerge_motif_logos(data05_hdbscan)),
  tar_target(files06_motifs, {
    outdir <- "img/"
    files <- purrr::imap(plots06_motifs, function(plot, i) {
      outfile <- file.path(outdir, paste0("clusterlogo_", i, ".png"))
      ggplot2::ggsave(outfile, plot, width=6, height=4)
      outfile
    })
    unlist(files)
  }, format = "file"),
  # 07: Unify outputs into a report
  tar_target(
    img_all,
    c(file04_tsne, file05_hdbscan, files06_motifs),
    format = "file"
  ),
  tar_render(
    report_html,
    "report.Rmd",
    params = list(
      processed_csv = file02_processed,
      mca_csv       = file03_mca,
      images        = img_all
    ),
    output_file = file.path(outdir, "report.html")
  )
)