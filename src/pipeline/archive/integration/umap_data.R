library(tidyverse)
library(here)
library(stats)
library(umap)

setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("integration/umap_data.R")

infile  <- "R255X_E488QD_modeling.csv"
indir   <- "data"
outfile <- "R255X_E488QD_umap.csv"

edit_threshold <- 0.01


# ============================================================================ #
# Process data into a UMAP representation

main <- function() {
  model_data <- load_model_data()
  umap_data  <- generate_umap_representation(model_data)
  
  write_csv(umap_data, here(indir, outfile))
}


load_model_data <- function() {
  model_data <- read_csv(here(indir, infile)) %>%
    filter(n >= read_threshold) %>%
    filter(edit >= edit_threshold)
  return(model_data)
}

generate_umap_representation <- function(model_data) {
  feature_data <- model_data %>% select(-N10, -GGA, -GAA, -n, -n_other, -edit)
  
  pca_result <- prcomp(feature_data, center=TRUE, scale.=TRUE, rank.=10)
  pca_matrix <- pca_result$x
  
  message("UMAP reduction started...")
  feature_umap <- umap(pca_matrix)
  edit_data    <- tibble(
    N10  = model_data$N10,
    X1   = feature_umap$layout[,1],
    X2   = feature_umap$layout[,2],
    edit = model_data$edit
  )
  
  return(edit_data)
}


main()
