library(tidyverse)
library(Rtsne)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("pipeline/6_tsne/tsne_embedding.R")

indir   <- "data"
infile  <- "R255X_E488QD_5a_MCAtopeditors.csv"
outfile <- "R255X_E488QD_6_tSNE.csv"


# ============================================================================ #
# Dimensional reduction of MCA data via t-SNE, in order to preserve local
# neighborhood structure for clustering algorithm.

main <- function() {
  k     <- 30
  
  X     <- read_csv(here(indir, infile))
  X_mca <- as.matrix(X[ , paste0("Axis", 1:k) ])
  
  n     <- nrow(X_mca)
  alpha <- n / 10      # exaggeration factor from Linderman & Steinerberger
  h     <- 1           # learning rate from Linderman & Steinerberger
  
  set.seed(42)
  tsne_fit <- Rtsne(
    as.matrix(X_mca),
    exaggeration_factor = alpha,
    eta                 = h,
    dims                = 2,
    perplexity          = 30,
    verbose             = TRUE,
    max_iter            = 2000,
    pca                 = FALSE
  )
  
  X_tsne <- X %>%
    mutate(
      tsne1 = tsne_fit$Y[, 1],
      tsne2 = tsne_fit$Y[, 2]
    )
  
  # update this image later for acs
  p <- ggplot(X_tsne, aes(x = tsne1, y = tsne2, colour = map)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_colour_gradient(low = "skyblue", high = "red") +
    theme_minimal()
  print(p)
  
  write_csv(X_tsne, here(indir, outfile))
}

main()
