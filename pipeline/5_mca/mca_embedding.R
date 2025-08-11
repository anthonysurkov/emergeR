library(tidyverse)
library(ade4)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/mca/mca_embedding.R")

indir     <- "data"
infile    <- "R255X_E488QD_2_clean.csv"
outfile_a <- "R255X_E488QD_3a_topeditors.csv"
outfile_b <- "R255X_E488QD_3b_alleditors.csv"

# Globals
L <- 10 # sequence length

# ============================================================================ #
# Multiple correspondence analysis (MCA) of top-editing sequences.

# Factors sequences into columns of characters.
factor_sequences <- function(X) {
  Xfac <- X$N10 %>%
    strsplit("") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(Xfac) <- paste0("p", 1:L)
  Xfac[] <- lapply(Xfac, function(x) factor(x, levels = c("A", "G", "C", "T")))
}

# Organizes MCA data attributed to all editors above "edit" cutoff.
get_MCs <- function(X, mca, edit) {
  scores <- as_tibble(mca$li) %>%
    mutate(N10 = X$N10) %>%
    left_join(X, by = "N10") %>%
    filter(map >= edit)
}

main <- function() {
  X    <- read_csv(here(indir, infile))                # read data
  Xfac <- factor_sequences(X = X)                      # factor into characters
  mca  <- dudi.acm(df = Xfac, nf = 50, scannf = FALSE) # conduct MCA
  
  top_cutoff <- quantile(X$map, 0.999)                       # 99.9th %tile
  edit_0.02  <- get_MCs(X = X, mca = mca, edit = 0.02)       # get >2% editors
  edit_top   <- get_MCs(X = X, mca = mca, edit = top_cutoff) # get top editors
  
  print(mca$eig)        # observe eigenvalues
  print(head(edit_top)) # sanity check
  
  write_csv(edit_top,  here(indir, outfile_a)) # save 99.9th %tile editors
  write_csv(edit_0.02, here(indir, outfile_b)) # save >2% editors
}

main()
