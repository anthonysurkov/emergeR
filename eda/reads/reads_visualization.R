library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/reads/reads_visualization.R")

infile <- "R255X_E488QD_clean.csv"
indir  <- "data"


# ============================================================================ #
# How many reads do we get? What does their distribution look like?

R255X_E488QD <- read_csv(here("data", input_file))
R255X_E488QD_nz <- R255X_E488QD %>% # NZ = nonzero
  filter(edit > 0)

main <- function() {
  print(n_hist)
  print(n_hist_nz)
  print(n_hist_zoomed)
}


# Visualize distributions of reads
# Full data
n_hist <- ggplot(data=R255X_E488QD, aes(x=n)) +
  geom_histogram(binwidth=1, fill="#3B4252") +
  coord_cartesian(xlim = c(-1, 500), ylim= c(0, 60000)) +
  labs(
    title="Histogram of Read Counts (binwidth=1)",
    x="Read Count (n)",
    y="Frequency"
  ) +
  theme_minimal()

# Full data \ any zero-editing sequences
n_hist_nz <- ggplot(data=R255X_E488QD_nz, aes(x=n)) +
  geom_histogram(binwidth=10, fill="#3B4252") +
  coord_cartesian(xlim = c(-1, 1500)) +
  labs(
    title="Distribution of Nonzero Read Counts",
    x="Read Count (n)",
    y="Frequency"
  ) +
  theme_minimal()

# Zoom in on first 20 entries of full data
n_hist_zoomed <- ggplot(data=R255X_E488QD, aes(x=n)) +
  geom_histogram(
    binwidth=1, fill="#3B4252", color="white", boundary=0, closed="left"
  ) +
  coord_cartesian(xlim = c(0, 20)) +
  labs(
    title="Distribution of Read Counts (x=0,1,...,10)",
    x="Read Count (n)",
    y="Frequency"
  ) +
  theme_minimal()


main()
