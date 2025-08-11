library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR")
here::i_am("modeling/structure/mfe_analysis.R")

indir   <- "data"
outdir  <- "modeling/structure"
infile  <- "R255X_E488QD_mfe_paired.csv"
outfile <- "R255X_E488QD_base_pair_freqs.csv"


# ============================================================================ #
# Analyzes patterns in minimum free energy (MFE) predicted secondary structures
# of best-editing guides.

main <- function() {
  mfe <- read_csv(here(indir, infile))
  
  # n1  = n66
  # n2  = n67
  # n3  = n68
  # n4  = n69
  # n5  = n70
  # n6  = n71
  # n7  = n72
  # n8  = n73
  # n9  = n74
  # n10 = n75
  # at  = target A
  
  grps <- mfe %>%
    group_by(group) %>%
    group_split()
  
  grp_counts <- map(
    grps, ~ .x %>%
    select(group, seq, edit, mfe, at, n66:n75) %>%
    pivot_longer(
      cols = c(at, n66:n75),
      names_to = "NX",
      values_to = "value"
    ) %>%
    count(NX, value, name="count") %>%
    group_by(NX) %>%
    mutate(prop = count / sum(count)) %>%
    arrange(NX) %>%
    ungroup()
  )
  
  for (i in seq_along(grp_counts)) {
    cat("=== Group", i, "===\n")
    print(grp_counts[[i]], n=Inf)
  }
  
  grp_counts <- bind_rows(grp_counts, .id="group")
  write_csv(grp_counts, here(outdir, outfile))
}

main()
