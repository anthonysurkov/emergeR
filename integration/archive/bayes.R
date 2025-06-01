library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("integration/bayes.R")

infile <- "R255X_E488QD_clean.csv"
indir  <- "data"


# ============================================================================ #
# Introducing Bayesian estimator

main <- function() {
  clean_data <- read_csv(here(indir, infile))
  
  alpha <- 1.37949
  beta  <- 89.8415
  
  data <- read_csv(here(indir, infile))
  data <- data %>%
    mutate(
      bayes = (alpha + GGA) / (alpha + beta + n),
      lower_credible = qbeta(0.025, alpha + GGA, beta + n - GGA),
      upper_credible = qbeta(0.975, alpha + GGA, beta + n - GGA),
      mle = edit,
      k = GGA
    ) %>%
    select(-GAA, -GGA)
  
  data_bayesian <- data %>% 
    arrange(desc(bayesian_estimator))
  print(data_bayesian, n=100)
  
  print(colnames(data_bayesian))
  
}

main()
