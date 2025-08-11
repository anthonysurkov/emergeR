library(tidyverse)
library(here)
library(bbmle)
library(transport)
library(scales)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/edits/edits_distribution.R")

infile <- "R255X_E488QD_3_cleaned.csv"
indir  <- "data"


# ============================================================================ #
# What is the distribution of our edits as described by the MLE estimator?
# Can we successfully fit a beta-binomial distribution to the data?

main <- function() {
  # Data
  data_all <- read_csv(here(indir, infile))
  
  k <- data_all %>% pull(k)
  N <- data_all %>% pull(n)
  epsilon  <- 0.01
  
  # Fitting beta-bino
  nll_bb <- function(log_alpha, log_beta) {
    alpha <- exp(log_alpha)
    beta  <- exp(log_beta)
    -sum(lchoose(N, k) + lbeta(k + alpha, N - k + beta) - lbeta(alpha, beta))
  }
  
  mle_fit <- mle2(nll_bb,
                  start = list(log_alpha = 0, log_beta = 0),
                  method = "Nelder-Mead")

  print(summary(mle_fit))
  params    <- coef(mle_fit)
  alpha_hat <- exp(params["log_alpha"])
  beta_hat  <- exp(params["log_beta"])
  
  simulate_beta_binomial <- function(n_trials, alpha, beta) {
    p <- rbeta(length(n_trials), alpha, beta)
    rbinom(length(n_trials), n_trials, p)
  }
  
  k_sim <- simulate_beta_binomial(N, alpha_hat, beta_hat)
  
  # Diagnostics
  edit_rate_obs <- k / N
  edit_rate_sim <- k_sim / N
  obs_cdf <- ecdf(edit_rate_obs)
  syn_cdf <- ecdf(edit_rate_sim)
  xs <- seq(0, 1, length.out=1000)
  dif <- obs_cdf(xs) - syn_cdf(xs)
  
  rmse <- sqrt(mean(dif^2))
  mae <- mean(abs(dif))
  wass <- wasserstein1d(edit_rate_obs, edit_rate_sim)
  
  cat("Diagnostics:\n")
  cat("Alpha_hat:", alpha_hat, "\n")
  cat("Beta_hat:", beta_hat, "\n")
  cat("RMSE:", rmse, "\n")
  cat("MAE:", mae, "\n")
  cat("Wasserstein distance:", wass, "\n") 
  
  df_plot <- data.frame(
    EditingRate = c(edit_rate_obs, edit_rate_sim),
    Source = rep(c("Observed", "Simulated"), each = length(N))
  )
  dens_plot <- ggplot(df_plot, aes(x = EditingRate, fill = Source)) +
    geom_density(alpha = 0.2) +
    ggtitle("Beta-Binomial Fit to Editing Rates") +
    theme_minimal()
  print(dens_plot)
}


main()
