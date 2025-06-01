library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("nd/nd_validation.R")

# ============================================================================ #

old_ac <- read_csv(here("nd/data", "old_AC_data.csv"))
old_ac <- old_ac %>%
  mutate(n = Reads) %>%
  select(N10, GAA, GGA, edit, n)

clean_data <- read_csv(here("data", "R255X_E488QD_clean.csv"))
clean_data <- clean_data #%>% select(-n_other)


joined <- inner_join(clean_data, old_ac, by="N10", suffix=c("_1", "_2"))

comparison <- joined %>%
  select(N10, edit_1, edit_2, n_1, n_2) %>%
  mutate(
    edit_diff = edit_1 - edit_2,
    n_diff = n_1 - n_2
  )

comparison_edits <- comparison %>%
  filter(abs(edit_diff) > 1e-6)

comparison_n <- comparison %>%
  filter(n_diff != 0)

print(comparison_edits, n=100)
print(comparison_n)
