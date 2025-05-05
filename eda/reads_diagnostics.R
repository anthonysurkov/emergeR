library(tidyverse)
library(here)
library(MASS)
library(dgof) # for a discrete-aware KS test


# Diagnostics on number of reads for the R255X_E488QD dataset.
# Temporary: uncomment necessary print statements.
input_file <- "R255X_E488QD_clean.csv"


# ============================================================================ #


R255X_E488QD <- read_csv(here("data", input_file))
R255X_E488QD_nz <- R255X_E488QD %>% # NZ = nonzero
  filter(edit > 0)

# Fit NB parameters to data --> generate real NB with same parameters -->
# evaluate goodness-of-fit.
# obs=observed, syn=synthetic


diagnostics_reads <- function(clean_data, print_plots=TRUE, plot_name="") {
  dataset <- deparse(substitute(clean_data))
  cat("\nDataset name: ", dataset, "\n")
  
  obs_n     <- clean_data$n
  obs_model <- fitdistr(obs_n, densfun = "Negative Binomial")
  syn_n     <- rnbinom(
    obs_n, size=obs_model$estimate["size"], mu=obs_model$estimate["mu"]
  )
  cat("Parameter estimates:\n")
  print(obs_model$estimate)

  # Use Kolmogorov-Smirnov and CDF RMSE tests to evaluate goodness-of-fit.
  # (p value is not very meaningful at this dataset size)
  ks <- dgof::ks.test(
    obs_n, 
    "pnbinom", 
    size=obs_model$estimate["size"], 
    mu=obs_model$estimate["mu"],
    simulate.p.value=TRUE,
    B=5000
  )
  cat("KS: ", ks$statistic, '\n')

  obs_ecdf <- ecdf(obs_n)
  syn_ecdf <- ecdf(syn_n)
  xs <- seq(0, max(obs_n), length.out=1000)
  dif <- obs_ecdf(xs) - syn_ecdf(xs)
  rmse <- sqrt(mean(dif^2))
  mae <- mean(abs(dif))
  cat("RMSE: ", rmse, '\n')
  cat("MAE: ", mae, '\n')


  # eCDF plots to visually compare theoretical vs empirical distributions.
  df_cdf <- tibble(
    value = c(obs_n, syn_n),
    type = rep(c("Observed reads", "NB fitted reads"), each = length(obs_n))
  )
  plot_title <- paste0(
    "(", plot_name, 
    ") Empirical vs Theoretical Negative Binomial Distribution (CDF)"
  )
  ecdf_plot <- ggplot(df_cdf, aes(x=value, color=type)) +
    stat_ecdf(geom="step") +
    coord_cartesian(xlim=c(0,1000)) +
    labs(
      title=str_wrap(plot_title, width=60),
      x="Reads (n)",
      y="Cumulative Probability",
      color="Read Type"
    ) +
    theme_minimal()
  
  if (print_plots) print(ecdf_plot)

  
  # Find where KS difference is coming from. Find max and plot graph of CDF gap.
  i_max <- which.max(abs(dif))
  x_at_max <- xs[i_max]
  gap_at_max <- dif[i_max]

  cat("X at max KS difference: ", x_at_max, '\n')
  cat("Gap at max difference: ", gap_at_max, '\n')

  cat("Data quantiles:\n")
  print(quantile(obs_n, probs=c(.5, .9, .95, .99, 1)))

  df_gap <- data.frame(x=xs, dif=dif)
  gap_plot <- ggplot(df_gap, aes(x=x, y=dif)) +
    geom_line() +
    geom_vline(xintercept=x_at_max, linetype="dashed") +
    labs(
      title=paste0("(",plot_name,") CDF Gap (Empirical vs. Theoretical)"),
      x="Reads (n)",
      y="eCDF Difference"
    ) +
    theme_minimal()
  
  if (print_plots) print(gap_plot)
}


# ============================================================================ #


main <- function() {
  diagnostics_reads(
    R255X_E488QD, print_plots=TRUE, plot_name="Full dataset"
  ) 
  diagnostics_reads(
    R255X_E488QD_nz, print_plots=TRUE, plot_name="Zeroes removed"
  )
}

main()
