library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("adsurkov/ND1to14_MAPs.R")

indir  <- "data"
infile <- "R255X_E488QD_clean.csv"


# ============================================================================ #
# Quick-n-dirty requisitions of a few ND propositions

main <- function() {
  data <- read_csv(here(indir, infile))
  
  seqs <- c(
    ND1  = "TTTGTTCCAC",
    ND2  = "AGAGTTCCGC",
    ND3  = "TTCTTTCGAC",
    ND4  = "TTTTCCTATT",
    ND5  = "AGAGTTCCAA",
    ND6  = "ATTTCCGATT",
    ND7  = "TGAATTCCGC",
    ND8  = "AGGGTTCCAT",
    ND9  = "TCATTCCGGT",
    ND10 = "AGATTTTGAC",
    ND11 = "CGAATTCCGC",
    ND12 = "CCTGTTGATT",
    ND13 = "AGTTTTTGAC",
    ND14 = "CGCTCTCGAC"
  )
  seq_lookup <- tibble(
    id = names(seqs),
    N10 = unname(seqs)
  )
  
  seq_data <- data %>%
    filter(N10 %in% seqs) %>%
    left_join(seq_lookup, by="N10") %>%
    arrange(desc(id))
  
  print(seq_data, n=Inf)
  
}

main()
