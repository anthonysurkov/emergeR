library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("nd/nd_data_gen.R")

indir  <- "data"
infile <- "R255X_E488QD_clean.csv"


# ============================================================================ #
# Generate top-100 MLE sequences for ND's R255X E488QD dataset

R255X_clean <- read_csv(here(indir, infile))
R255X_clean <- R255X_clean %>% filter(edit < 1) %>% filter(n >= 10)
R255X_top100 <- R255X_clean[1:100, ]

#print(R255X_top100)
write_csv(R255X_clean, here("nd/data", "R255X_E488QD_fulldata.csv"))
write_csv(R255X_top100, here("nd/data","R255X_E488QD_top_100.csv"))
