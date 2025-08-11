library(tidyverse)
library(here)
library(stringr)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/edits/tcc_motifs/tcc_6to8/tcc_6to8_model.R")

infile  <- "tcc_6to8.csv"
outfile <- "tcc_6to8_model.csv"
indir   <- "eda/edits/tcc_motifs/tcc_6to8"

nucleotides <- c("A", "C", "T", "G")
L <- 10 # length


# ============================================================================ #
# Modeling of submotifs in a TCC @ nt 6~8 context

main <- function() {
  clean_data <- read_csv(here(indir, infile))
  
  long_data <- clean_data %>%
    select(N10) %>%
    mutate(chars = strsplit(N10, "")) %>%
    unnest(chars) %>%
    group_by(N10) %>%
    mutate(position = row_number(), value = 1L) %>%
    ungroup() %>%
    complete(N10, position = seq_len(L), chars = nucleotides,
             fill = list(value = 0L))

  wide_data <- long_data %>%
    pivot_wider(
      id_cols     = N10,
      names_from  = c(position, chars),
      names_sep   = "",
      values_from = value,
      values_fill = list(value = 0L)
    )

  X_main <- as.matrix(wide_data %>% select(-N10))
  feature_names <- colnames(X_main)

  interaction_names <- unique(unlist(
    apply(X_main, 1, function(row) {
      active <- feature_names[row == 1]
      unlist(lapply(2:L, function(k) combn(active, k,
        FUN = function(x) paste(x, collapse = ":"), simplify = TRUE)))
  })
))

  X_int <- matrix(0L, nrow = nrow(X_main), ncol = length(interaction_names),
                  dimnames = list(NULL, interaction_names))

  for (i in seq_len(nrow(X_main))) {
    message(i, " of ", nrow(X_main))
    active <- feature_names[X_main[i, ] == 1]
    combos_i <- unlist(lapply(2:L, function(k) combn(active, k,
      FUN = function(x) paste(x, collapse = ":"), simplify = TRUE)))
    X_int[i, combos_i] <- 1L
  }

  features_pruned <- cbind(
    wide_data,
    as.data.frame(X_int, check.names = FALSE)
  )
  
  message("cooked")

  response_cols <- setdiff(names(clean_data), "N10")
  model_data <- clean_data %>%
    distinct(N10, across(all_of(response_cols))) %>%
    right_join(features_pruned, by = "N10") %>%
    arrange(N10)
  
  message("cooketh")

  write.csv(model_data, out_file, row.names = FALSE)
  message("Saved encoded model data to: ", out_file)
}

main()
