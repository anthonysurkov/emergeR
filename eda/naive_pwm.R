library(tidyverse)
library(here)
library(ggseqlogo)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/naive_pwm.R")

infile <- "R255X_E488QD_clean.csv"
indir  <- "data"


# ============================================================================ #

main <- function() {
  data <- read_csv(here(indir, infile))
  
  all_seqs   <- data %>% pull(N10)
  edit_seqs  <- data %>% filter(edit != 0) %>% pull(N10)
  other_seqs <- data %>% filter(edit == 0) %>% pull(N10)
  
  all_logo <- ggseqlogo(all_seqs, method="bits") +
    ggtitle("R255X E488Qd: All Sequences")
  edit_logo <- ggseqlogo(edit_seqs, method="bits") +
    ggtitle("R255X E488Qd: All Editing-Capable Sequences")
  other_logo <- ggseqlogo(other_seqs, method="bits") +
    ggtitle("R255X E488Qd: All Editing-Incapable Sequences")
  
  all_logo   <- all_logo + ylim(0, 2)
  edit_logo  <- edit_logo + ylim(0, 2)
  other_logo <- other_logo + ylim(0, 2)
  
  print(all_logo)
  print(edit_logo)
  print(other_logo)
}

main()
