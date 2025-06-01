library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("integration/one_hot_encoding.R")

indir  <- "data"
infile <- "R255X_E488QD_clean.csv"


# ============================================================================ #
# One-hot encode sequence data to prepare it for a GLM
# (Deprecated; kept alive as a reference for future TODO)

clean_data <- read_csv(here(indir, infile))

nucleotides <- c("A", "T", "C", "G")
sequences   <- clean_data$N10
L           <- nchar(sequences[1])

long_data <- clean_data %>%
  select(N10) %>%
  mutate(chars = str_split(N10, "")) %>%
  unnest(chars) %>%
  group_by(N10) %>%
  mutate(
    position = row_number(),
    values   = 1L,
  ) %>%
  ungroup() %>%
  complete(
    N10,
    position = seq_len(L),
    chars    = nucleotides,
    fill     = list(values=0L)
  )

wide_data <- long_data %>%
  pivot_wider(
    id_cols     = N10,
    names_from  = c(position, chars),
    names_sep   = "",
    values_from = values,
    values_fill = list(values = 0L)
  )

model_data <- clean_data %>%
  distinct(N10, GGA, GAA, n, n_other, edit) %>%
  right_join(wide_data, by="N10") %>%
  arrange(N10)


print(model_data)
write_csv(model_data, here("data", "R255X_E488QD_modeling.csv"))
