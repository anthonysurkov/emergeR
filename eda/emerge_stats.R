library(tidyverse)
library(here)


setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/emerge_stats.R")

clean_data <- read_csv(here("data", "R255X_E488QD_clean.csv"))


# ============================================================================ #


main <- function() {
  print_basic_stats(clean_data)
  print_read_depths(clean_data)
  
  #highlight_nas(clean_data) 
  #graph_editing_frequencies(clean_data)
  #print_old_top_10(clean_data)
}


print_basic_stats <- function(my_data) {
  seqs <- my_data %>% pull(N10)
  good_seqs <- my_data %>% filter(n >= 10) 
  editing_seqs <- my_data %>% filter(edit > 0) %>% filter(n >= 10)
  
  edit1_seqs  <- editing_seqs %>% filter(edit >= 0.01) %>% pull(N10)
  edit5_seqs  <- editing_seqs %>% filter(edit >= 0.05) %>% pull(N10)
  edit10_seqs <- editing_seqs %>% filter(edit >= 0.10) %>% pull(N10)
  edit20_seqs <- editing_seqs %>% filter(edit >= 0.20) %>% pull(N10)
  edit30_seqs <- editing_seqs %>% filter(edit >= 0.30) %>% pull(N10)
  edit40_seqs <- editing_seqs %>% filter(edit >= 0.40) %>% pull(N10)
  edit50_seqs <- editing_seqs %>% filter(edit >= 0.50) %>% pull(N10)
  edit60_seqs <- editing_seqs %>% filter(edit >= 0.60) %>% pull(N10)
  edit70_seqs <- editing_seqs %>% filter(edit >= 0.70) %>% pull(N10)
  
  good_seqs    <- good_seqs %>% pull(N10)
  editing_seqs <- editing_seqs %>% pull(N10)
  
  cat("\nR255X E488Qd: Basic Statistics\n\n")
  cat("Total sequences screened:", length(seqs), '\n\n')
  cat("Num. sequences with 10 or more reads:", length(good_seqs), '\n')
  cat("Total editors (n>=10):", length(editing_seqs), '\n')
  cat("Total editors (n>=10) by edit% cutoffs::\n")
  cat(">= 1%:",  length(edit1_seqs),  '\n')
  cat(">= 5%:",  length(edit5_seqs),  '\n')
  cat(">= 10%:", length(edit10_seqs), '\n')
  cat(">= 20%:", length(edit20_seqs), '\n')
  cat(">= 30%:", length(edit30_seqs), '\n')
  cat(">= 40%:", length(edit40_seqs), '\n')
  cat(">= 50%:", length(edit50_seqs), '\n')
  cat(">= 60%:", length(edit60_seqs), '\n')
  cat(">= 70%:", length(edit70_seqs), '\n')
}

# Find NAs
highlight_nas <- function(my_data) {
  nas <- my_data %>% filter(is.na(edit))
  print(nas, n=Inf)
}


# Graph editing frequencies
graph_editing_frequencies <- function(my_data) {
  all_edits <- my_data %>%
    select(edit)
  nz_edits <- all_edits %>%
    filter(edit > 0)
  high_edits <- nz_edits %>%
    filter(edit > 0.05)
  print("High edit count:")
  print(nrow(high_edits))
  
  all_hist <- ggplot(data=all_edits, aes(x=edit)) +
    geom_histogram(binwidth=0.01, fill="#3B4252") +
    labs(
      title="Histogram of Edit Counts (binwidth=0.01)",
      x="Fraction of Reads Edited",
      y="Frequency"
    ) +
    theme_minimal()
  
  nz_hist <- ggplot(data=nz_edits, aes(x=edit)) +
    geom_histogram(binwidth=0.01, fill="#3B4252") +
    labs(
      title="Histogram of Nonzero Edit Counts (binwidth=0.01)",
      x="Fraction of Reads Edited",
      y="Frequency"
    ) +
    theme_minimal()
  
  high_hist <- ggplot(data=high_edits, aes(x=edit)) +
    geom_histogram(binwidth=0.01, fill="#3B4252") +
    labs(
      title="Histogram of Editing > 0.05 Frequencies (binwidth=0.01)",
      x="Fraction of Reads Edited",
      y="Frequency"
    ) +
    theme_minimal()
  
  print(all_hist)
  print(nz_hist)
  print(high_hist)
}


# Read depths
print_read_depths <- function(my_data) {
  my_data_nonzero <- my_data %>% filter(edit > 0)
  
  n    <- sum(my_data$n)
  n_nz <- sum(my_data_nonzero$n)
  
  read_depth    <- n / nrow(my_data)
  read_depth_nz <- n_nz / nrow(my_data_nonzero)
  
  cat("\nRead depth:", read_depth, '\n')
  cat("Editor read depth:", read_depth_nz, '\n')
}


# Find old top-10 sequences in data
print_old_top_10 <- function(my_data) {
  nd1  <- my_data %>% filter(N10 == "TTTGTTCCAC")
  nd2  <- my_data %>% filter(N10 == "AGAGTTCCGC")
  nd3  <- my_data %>% filter(N10 == "TTCTTTCGAC")
  nd4  <- my_data %>% filter(N10 == "TTTTCCTATT")
  nd5  <- my_data %>% filter(N10 == "AGAGTTCCAA")
  nd6  <- my_data %>% filter(N10 == "ATTTCCGATT")
  nd7  <- my_data %>% filter(N10 == "TGAATTCCGC")
  nd8  <- my_data %>% filter(N10 == "AGGGTTCCAT")
  nd9  <- my_data %>% filter(N10 == "TCATTCCGGT")
  nd10 <- my_data %>% filter(N10 == "CGAATTCCGC")

  cat("ND1  ", nd1$N10,  sprintf(" Edit: %.4f", nd1$edit),  
      " Reads: ", nd1$n, " Other reads: ", nd1$n_other, '\n')
  cat("ND2  ", nd2$N10,  sprintf(" Edit: %.4f", nd2$edit),
      " Reads: ", nd2$n, " Other reads: ", nd2$n_other, '\n')
  cat("ND3  ", nd3$N10,  sprintf(" Edit: %.4f", nd3$edit),
      " Reads: ", nd3$n, " Other reads: ", nd3$n_other, '\n')
  cat("ND4  ", nd4$N10,  sprintf(" Edit: %.4f", nd4$edit),
      " Reads: ", nd4$n, " Other reads: ", nd4$n_other, '\n')
  cat("ND5  ", nd5$N10,  sprintf(" Edit: %.4f", nd5$edit),
      " Reads: ", nd5$n, " Other reads: ", nd5$n_other, '\n')
  cat("ND6  ", nd6$N10,  sprintf(" Edit: %.4f", nd6$edit),
      " Reads: ", nd6$n, " Other reads: ", nd6$n_other, '\n')
  cat("ND7  ", nd7$N10,  sprintf(" Edit: %.4f", nd7$edit),
      " Reads: ", nd7$n, " Other reads: ", nd7$n_other, '\n')
  cat("ND8  ", nd8$N10,  sprintf(" Edit: %.4f", nd8$edit),
      " Reads: ", nd8$n, " Other reads: ", nd8$n_other, '\n')
  cat("ND9  ", nd9$N10,  sprintf(" Edit: %.4f", nd9$edit),
      " Reads: ", nd9$n, " Other reads: ", nd9$n_other, '\n')
  cat("ND10 ", nd10$N10, sprintf(" Edit: %.4f", nd10$edit),
      " Reads: ", nd10$n," Other reads: ", nd10$n_other,'\n')
}


main()
