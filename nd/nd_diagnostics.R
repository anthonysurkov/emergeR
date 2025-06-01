library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("nd/nd_diagnostics.R")

input_file <- "R255X_E488QD_clean.csv"
output_file <- "R255X_E488QD_nd.csv"

top_35_file <- "old_top_35_editors.csv"

# ============================================================================ #
# Attempt to diagnose differences between 2022 AC csv and 2025 full csv

new_data <- read_csv(here("data", input_file))
top_35 <- read_csv(here("nd/data", top_35_file), col_names=FALSE)
colnames(top_35)[1] <- "N10"

filtered <- new_data %>%
  filter(edit > 0) %>%
  filter(n >= 10)
#  %>% select(-n_other)


# VALIDATION

# See where ~5 top-35 sequences are lost:
top_35_in_filtered <- top_35 %>%
  filter(N10 %in% filtered$N10)
top_35_in_new <- top_35 %>%
  filter(N10 %in% new_data$N10)
print(nrow(top_35))
print(nrow(top_35_in_filtered))
print(nrow(top_35_in_new))
# found that n + n_other >= 10 but n < 10 for 5 top-35s in Excel
# write_csv(top_35_in_filtered, here("data/old_nd", "top_35_in_filtered.csv"))
# write_csv(top_35_in_new, here("data/old_nd", "top_35_in_new.csv"))

# See if %edits diverge
new_top35_merged <- top_35 %>%
  inner_join(new_data, by="N10")
print(new_top35_merged, n=Inf)
# write_csv(new_top35_merged, here("data/old_nd", "new_top35_merged.csv"))
# They do, along with the read counts :(
# Must be HTStream
