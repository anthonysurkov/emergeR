library(tidyverse)
library(stringr)
library(purrr)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/edits/tcc_motifs/seqs_by_TCC.R")

infile <- "R255X_E488QD_clean.csv"
outdir <- "eda/edits/tcc_motifs"
indir  <- "data"


# ============================================================================ #
# An exploration of sequences sharing position-variant TCC/TTG motifs

convert_to_pyrpur <- function(s) {
  s <- str_split(s, "", simplify=TRUE)
  s[s == "A" | s == "G"] <- "R" # purines
  s[s == "C" | s == "T"] <- "Y" # pyrimidines
  apply(s, 1, paste0, collapse="")
}

main <- function() {
  data <- read_csv(here(indir, infile))
  
  tcc_grps <- data %>%
    filter(map > quantile(map, 0.999, na.rm=TRUE)) %>%
    mutate(
      tcc_4to6 = str_sub(N10, 4, 6) == "TCC",
      tcc_5to7 = str_sub(N10, 5, 7) == "TCC",
      tcc_6to8 = str_sub(N10, 6, 8) == "TCC",
      tcc_7to9 = str_sub(N10, 7, 9) == "TCC",
      tcc_group = case_when(
        tcc_4to6 ~ "tcc_4to6",
        tcc_5to7 ~ "tcc_5to7",
        tcc_6to8 ~ "tcc_6to8",
        tcc_7to9 ~ "tcc_7to9",
        TRUE     ~ "None"
      )
    ) %>%
    group_by(tcc_group)
  prop_none <- sum(tcc_grps$tcc_group == "None") / nrow(tcc_grps)
  print(prop_none)
  
  tcc_list  <- group_split(tcc_grps, .keep=TRUE)
  tcc_names <- group_keys(tcc_grps)$tcc_group
  names(tcc_list) <- tcc_names
  tcc_4to6  <- tcc_list[["tcc_4to6"]] %>% select(N10, map)
  tcc_5to7  <- tcc_list[["tcc_5to7"]] %>% select(N10, map)
  tcc_6to8  <- tcc_list[["tcc_6to8"]] %>% select(N10, map)
  tcc_7to9  <- tcc_list[["tcc_7to9"]] %>% select(N10, map)
  
  print(tcc_4to6, n=Inf, width=Inf)
  print(tcc_5to7, n=Inf, width=Inf)
  print(tcc_6to8, n=Inf, width=Inf)
  print(tcc_7to9, n=Inf, width=Inf)
  write_csv(tcc_4to6, here(outdir, "tcc_4to6.csv"))
  write_csv(tcc_5to7, here(outdir, "tcc_5to7.csv"))
  write_csv(tcc_6to8, here(outdir, "tcc_6to8.csv"))
  write_csv(tcc_7to9, here(outdir, "tcc_7to9.csv"))
  
  groupwise_summary <- tcc_grps %>% summarise(avg_map = mean(map), count=n())
  print(groupwise_summary)
  
  # PYR/PUR EXPERIMENT
  lapply(tcc_list, function(tcc) {
    seqs   <- tcc %>% pull(N10)
    pyrpur <- map_chr(seqs, convert_to_pyrpur)
    print(ggseqlogo(pyrpur, method="bits"))
    print(ggseqlogo(seqs,   method="bits"))
  })

}

main()
