library(tidyverse)
library(here)
library(MASS)
library(transport)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/reads/reads_diagnostics.R")

infile <- "R255X_E488QD_clean.csv"
indir  <- "data"


# ============================================================================ #
# Does a negative binomial describe the distribution of reads well?

main <- function(print_plots=FALSE) {
  clean_data <- read_csv(here(indir, infile))
  
  N <- clean_data %>% pull(n)
  
  nb_fit <- fitdistr(N, densfun="Negative Binomial")
  print(summary(nb_fit))
  
  mu    <- nb_fit$estimate["mu"]
  size  <- nb_fit$estimate["size"]
  N_sim <- rnbinom(length(N), mu=mu, size=size)
  
  obs_cdf <- ecdf(N)
  syn_cdf <- ecdf(N_sim)
  xs <- seq(0, max(N), length.out=1000)
  dif <- obs_cdf(xs) - syn_cdf(xs)
  rmse <- sqrt(mean(dif^2))
  mae <- mean(abs(dif))
  wass <- wasserstein1d(N, N_sim) / mean(N)
  
  cat("Params:\n")
  cat("Mu:", mu, "\n")
  cat("Size:", size, "\n")
  cat("Diagnostics:\n")
  cat("RMSE:", rmse, "\n")
  cat("MAE:", mae, "\n")
  cat("Wasserstein distance:", wass, "\n") 
  
  df_plot <- data.frame(
    Reads = c(N, N_sim),
    Source = rep(c("Observed", "Simulated"), each = length(N))
  )
  dens_plot <- ggplot(df_plot, aes(x = Reads, fill = Source)) +
    geom_density(alpha = 0.4) +
    ggtitle("Negative Binomial Fit to Read Counts") +
    theme_minimal()
  print(dens_plot)
}


main()
