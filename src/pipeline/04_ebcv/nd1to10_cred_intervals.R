library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/edits/bayes/credible_intervals.R")

infile <- "R255X_E488QD_clean.csv"
indir  <- "data"


# ============================================================================ #
# What are the credible intervals for ND1~10?

main <- function() {
  alpha <- 1.37949
  beta  <- 89.8415
  
  data <- read_csv(here(indir, infile))
  data <- data %>%
    mutate(
      bayesian_estimator = (alpha + GGA) / (alpha + beta + n),
      lower_credible = qbeta(0.025, alpha + GGA, beta + n - GGA),
      upper_credible = qbeta(0.975, alpha + GGA, beta + n - GGA)
    )
  
  top_n10s <- c(
    "TTTGTTCCAC", "AGAGTTCCGC", "TTCTTTCGAC",
    "TTTTCCTATT", "AGAGTTCCAA", "ATTTCCGATT",
    "TGAATTCCGC", "AGGGTTCCAT", "TCATTCCGGT",
    "CGAATTCCGC")
  top_n10_names <- c(
    "ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "ND7", "ND8", "ND9", "ND10"
  )
  
  data_bayesian_top <- data %>% 
    arrange(desc(bayesian_estimator)) %>%
    dplyr::select(N10, bayesian_estimator)
  print(data_bayesian_top, n=100)
  
  data_top <- data %>%
    filter(N10 %in% top_n10s) %>%
    mutate(N10 = factor(N10, levels=top_n10s)) %>%
    dplyr::select(N10, bayesian_estimator, lower_credible, upper_credible)
  
  data_frequentist_top <- data %>%
    filter(edit != 1) %>%
    arrange(desc(edit)) %>%
    dplyr::select(N10, edit)
  print(data_frequentist_top, n=10)
  
  data_top_plot <- ggplot(data_top, aes(x=N10, y=bayesian_estimator)) +
    geom_boxplot(outlier.shape=NA, fill="gray90") +
    geom_errorbar(
      aes(ymin=lower_credible, ymax=upper_credible),
      width=0.2,
      color="steelblue",
      size=0.8
    ) +
    geom_point(
      aes(y=bayesian_estimator),
      shape=21, fill="white",
      color="steelblue",
      size=2
    ) +
    theme_minimal() +
    scale_x_discrete(labels = top_n10_names) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    labs(
      x = "Sequence",
      y = "Editing rate",
      title = "Editing Rate of Top Guides With 95% Credible Whiskers"
    )
  print(data_top_plot)
}


main()
