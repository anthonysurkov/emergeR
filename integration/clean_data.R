library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("integration/clean_data.R")


# Cleans Casey pipeline output into CSV with...
# N10 sequence, GAA count, GGA count, n=GAA+GGA, editing %, n_other=anything not
# GAA or GGA;
# ...sorted by decreasing editing %.

input_file <- "R255X_E488QD_casey_output.csv"
output_file <- "R255X_E488QD_clean.csv"


# ============================================================================ #


full_data <- read_csv(here("data", input_file))
colnames(full_data)[1] <- "N10"

#R255X_E488QD <- full_data %>%
#  mutate(n = as.numeric(GAA + GGA)) %>%
#  mutate(n_other = rowSums(
#    select(., where(is.numeric), -GAA, -GGA, -n),
#    na.rm=TRUE
#  )) %>%
#  mutate(edit = if_else((GAA + GGA) == 0, 0, as.numeric(GGA / (GAA + GGA)))) %>%
#  select(N10, GAA, GGA, n, edit, n_other) %>%
#  arrange(desc(edit))

R255X_E488QD <- full_data %>%
  mutate(n = as.numeric(rowSums(
    select(., where(is.numeric)),
    na.rm=TRUE
  ))) %>%
  mutate(edit = if_else(n == 0, 0, as.numeric(GGA / n))) %>%
  select(N10, GAA, GGA, n, edit) %>%
  arrange(desc(edit))

write_csv(R255X_E488QD, here("data", output_file))
