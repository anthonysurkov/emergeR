library(tidyverse)
library(here)
library(binom)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("integration/clean_data.R")

infile  <- "R255X_E488QD_casey_output.csv"
indir   <- "data"
outfile <- "R255X_E488QD_clean.csv"
outdir  <- "data"


# ============================================================================ #
# Processes Casey pipeline output into clean tibbles with statistics attached.
# All Bayesian statistics are not cross-validated. Please see k_fold_eb.R.
# Treat them as placeholders for now.

main <- function() {
  alpha <- 1.37949
  beta  <- 89.8415
  
  preprocessed <- read_csv(here(indir, infile)) 
  colnames(preprocessed)[1] <- "N10"
  
  clean_data <- preprocessed %>%
    mutate(
      n = as.numeric(rowSums(
        select(., where(is.numeric)),
        na.rm=TRUE)
      ),
      k = GGA
    ) %>%
    mutate(
      mle = if_else(n == 0, 0, as.numeric(k / n)),
      map = (alpha + k) / (alpha + beta + n),
      cis = binom.confint(k, n, method="wilson"),
      lower_ci = cis$lower,
      upper_ci = cis$upper,
      lower_cred = qbeta(0.025, alpha + k, beta + n - k),
      upper_cred = qbeta(0.975, alpha + k, beta + n - k)
    ) %>%
    select(
      N10, n, k, 
      map, lower_cred, upper_cred, 
      mle, lower_ci, upper_ci
    ) %>%
    arrange(desc(map))
  
  write_csv(clean_data, here(outdir, outfile))
}

main()
