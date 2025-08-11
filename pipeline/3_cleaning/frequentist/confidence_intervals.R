library(tidyverse)
library(here)
library(binom)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/edits/frequentist/confidence_intervals.R")

infile <- "R255X_E488QD_clean.csv"
indir  <- "data"


# ============================================================================ # 
# Frequentist confidence intervals for ND1~10 guides

main <- function() {
  Z <- 1.960 # for two-sided 95% CI
  
  data <- read_csv(here(indir, infile))
  
  top_n10s <- c(
    "TTTGTTCCAC", "AGAGTTCCGC", "TTCTTTCGAC",
    "TTTTCCTATT", "AGAGTTCCAA", "ATTTCCGATT",
    "TGAATTCCGC", "AGGGTTCCAT", "TCATTCCGGT",
    "CGAATTCCGC")
  top_n10_names <- c(
    "ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "ND7", "ND8", "ND9", "ND10"
  )
  
  data_top <- data %>%
    filter(N10 %in% top_n10s) %>%
    mutate(N10 = factor(N10, levels=top_n10s))
  
  print(data_top)
    
  data_top_plot <- ggplot(data_top, aes(x=N10, y=mle)) +
    geom_boxplot(outlier.shape=NA, fill="gray90") +
    geom_errorbar(
      aes(ymin=lower_CI, ymax=upper_CI),
      width=0.2,
      color="steelblue",
      size=0.8
    ) +
    geom_point(
      aes(y=edit),
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
      title = "Editing Rate of Top Guides With 95% Wilson CIs"
    )
  print(data_top_plot)
}

main()
