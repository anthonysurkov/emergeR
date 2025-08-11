library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/seed_pwms/seed_pwms.R")

infile <- "R255X_E488QD_clean.csv"
indir  <- "data"


# ============================================================================ #
# Is the MLE or MAP preferable for our data? (Deprecated - faulty logic)

main <- function() {
  data <- read_csv(here(indir, infile)) 
  
  alpha <- 1.37949
  beta  <- 89.84158 
  A     <- alpha + beta
  
  data <- data %>%
    mutate(D = (alpha - A * edit)^2 - (2 * A * edit * (1-edit))) %>%
    mutate(N_not = (A^2 * edit * (1-edit) / D)) %>%
    mutate(EB_flag = as.integer(n < N_not))
  
  EB_prop <- data %>% summarise(prop = mean(EB_flag)) %>% pull(prop)
  EB_edit <- data %>% filter(EB_flag == 1) %>% pull(edit)
  EB_n    <- data %>% filter(EB_flag == 1) %>% pull(n)
  
  EB_edit_mean <- mean(EB_edit)
  EB_edit_var  <- var(EB_edit)
  EB_n_mean    <- mean(EB_n)
  EB_n_var     <- var(EB_n)
  
  cat("Proportion of data benefitting from EB estimator:", EB_prop, "\n")
  cat("Average editor benefitting from EB estimator:", EB_edit_mean, "\n")
  cat("Variance of said editors:", EB_edit_var, "\n")
  cat("Average overall read count:", mean(data$n), "\n")
  cat("Average read count of editors benefitting from EB estimator:")
  cat(EB_n_mean, "\n")
  cat("Stdev of said read counts:", sqrt(EB_n_var), "\n")
}

main()
