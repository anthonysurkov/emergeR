library(tidyverse)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("integration/k_fold_eb.R")

indir   <- "data"
outdir  <- "data"
infile  <- "R255X_E488QD_clean.csv"
outfile <- "R255X_E488QD_CV.csv"

num_folds <- 5


# ============================================================================ #
# K-fold empirical Bayesian procedure. Estimates editing rate based on a biased
# estimator derived from a fitted beta-binomial.

# Helper functions:
# Fit MLE parameters for beta-binomial/
fit_bb_mle <- function(k, n, start = c(log(1), log(1))) {
  bb_loglik <- function(par, k, n) {
    a <- exp(par[1]);
    b <- exp(par[2]);
    sum(
      lchoose(n, k) +
      lbeta(k + a, n - k + b) - 
      lbeta(a, b)
    )
  }
  opt <- optim(
    start, bb_loglik, k=k, n=n, 
    control=list(fnscale=-1), method="L-BFGS-B"
  )
  list(alpha = exp(opt$par[1]), beta = exp(opt$par[2]))
}

# Compute posterior mean.
post_mean <- function(k, n, alpha, beta) (k + alpha) / (n + alpha + beta)

# Compute log-pmf of beta-binomial.
dbb_logpmf <- function(k, n, alpha, beta) {
  lchoose(n, k) +
    lbeta(k + alpha, n - k + beta) -
    lbeta(alpha, beta)
}


main <- function() {
  clean_data <- read_csv(here(indir, infile))
  
  # Split data into 5 folds. Ensure proportion of editors is same per fold.
  model_data <- clean_data %>%
    select(N10, k, n, mle) %>%
    mutate(fold = sample(rep(1:5, length.out = n()))) %>%
    ungroup()
  
  # K-fold EB
  results <- model_data %>%
    group_nest(fold) %>%
    mutate(
      train = map2(data, fold, ~filter(model_data, fold != .y)),
      test  = data,
      fit   = map(train, ~fit_bb_mle(.x$k, .x$n)),
      
      stats = pmap(
        list(test=test, fit=fit),
        function(test, fit) {
          k <- test$k
          n <- test$n
          a <- fit$alpha
          b <- fit$beta
          N <- length(k)
          
          ll   <- sum(dbb_logpmf(k, n, a, b)) # test log-likelihood
          mnll <- -ll / N                     # mean neg-Log-lik
          aic  <- -2*ll + 2*2                 # AIC (p=2 params)
          
          # Pearson chi-square GOF
          mu_k    <- n*(a/(a+b))
          var_k   <- n*(a/(a+b))*(b/(a+b)) * ((a+b+n)/(a+b+1))
          chi2    <- sum((k - mu_k)^2 / var_k)
          
          # RMSE, MAE
          rmse_k  <- sqrt(mean((k - mu_k)^2))
          mae_k   <- mean(abs(k - mu_k))
          
          tibble(
            alpha       = a,
            beta        = b,
            test_n      = N,
            logLik      = ll,
            meanNegLL   = mnll,
            AIC         = aic,
            PearsonChi2 = chi2,
            RMSE_counts = rmse_k,
            MAE_counts  = mae_k
          )
        }),
      
        predictions = map2(test, fit, function(test, fit) {
          a <- fit$alpha
          b <- fit$beta
          map <- (test$k + a) / (test$n + a + b)
          lower_cred <- qbeta(0.025, test$k + a, test$n - test$k + b)
          upper_cred <- qbeta(0.975, test$k + a, test$n - test$k + b)
          
          test %>%
            select(N10) %>%
            mutate(
              map        = map,
              lower_cred = lower_cred,
              upper_cred = upper_cred,
              alpha      = a,
              beta       = b
            )
        })
    )
  
  posterior_updates <- bind_rows(results$predictions)
  clean_data <- clean_data %>%
    select(-map, -lower_cred, -upper_cred) %>%
    left_join(posterior_updates, by="N10") %>%
    arrange(
      N10, n, k, 
      map, upper_cred, lower_cred, 
      mle, upper_ci, lower_ci, 
      alpha, beta
    )
  
  diagnostics <- results %>%
    select(fold, fit, stats) %>%
    tidyr::unnest_wider(stats) %>%
    tidyr::unnest(fit)
  
  print(diagnostics)
  print(clean_data)
  
  write_csv(clean_data, here(outdir, outfile))
}

main()
