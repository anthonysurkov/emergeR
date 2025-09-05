#' Calculate statistics on a cleaned EMERGe dataset.
#' 
#' Calculates editing MLE, Wilson CIs (binomial), posterior mean under an
#' empirical Bayes Beta prior (with k-fold CV), and Bayesian credible intervals.
#' 
#' @param X_clean (tibble) Cleaned data from emergeR::clean.R.
#' @param num_folds (integer) Default 5, number of folds to use for k-fold CV.
#' @param seed (float) Default 1, controls RNG for reproducibility.
#' @param get_cv_stats (logical) Default FALSE, returns CV diagnostics if TRUE.
#' @return If get_cv_stats is TRUE, a list with two components:
#'   \describe{
#'     \item{X_stats}{A tibble containing EMERGe data: n10, mle, upper_ci,
#'     lower_ci, postmean, upper_cred, lower_cred, alpha, beta.}
#'     \item{cv_stats}{A tibble containing alpha MLE, beta MLE, number of
#'           samples in the held-in fold (test_n), log-likelihood of the
#'           test data (logLik), mean negative log-likelihood per observation
#'           (meanNegLL), Akaike Information Criterion for the fit (p = 2)
#'           (AIC), Pearson chi-square statistic comparing observed vs.
#'           expected counts (PearsonChi2), Root mean squared error between
#'           observed and expected counts (RMSE_counts), Mean absolute
#'           error between observed and expected counts (MAE_counts),
#'           Optimizer convergence code (0 indicates success) (convergence),
#'           and the Optimizer convergence message (message).}
#'   }
#'   Otherwise, only X_stats is returned.
#' @export
#' 
#' @examples
#' \dontrun{
#' emergeR::stats(X_clean = clean_data)
#' emergeR::stats(X_clean = clean_data, num_folds = 4)
#' emergeR::stats(X_clean = clean_data, num_folds = 4, seed = 42)
#' }
#'
#' @importFrom dplyr mutate if_else select left_join 
#' @importFrom dplyr arrange group_by summarise n
#' @importFrom purrr map map2 pmap
#' @importFrom tidyr group_nest
#' @importFrom tibble tibble
#' @importFrom binom binom.confint
#' @importFrom stats optim qbeta

append_stats <- function(
    X_clean, num_folds = 5, seed = 1, get_cv_stats = FALSE
) {
  stopifnot(all(c("n10", "n", "k") %in% names(X_clean)))
  if (!is.null(seed)) set.seed(seed)
  
  X_mle <- X_clean %>% mutate(mle = if_else(n == 0, 0, as.numeric(k / n)))
  cis <- binom.confint(x = X_mle$k, n = X_mle$n, method = "wilson")
  X_mle <- X_mle %>%
    mutate(
      lower_ci = cis$lower,
      upper_ci = cis$upper
    )
  
  cv_eb_out <- cv_eb(data = X_clean, num_folds = num_folds, seed = seed)
  X_bayes <- cv_eb_out$X_bayes
  
  X_out <- X_mle %>%
    left_join(X_bayes$eb_data, by = "n10") %>%
    arrange(
      n10, n, k,
      mle, lower_ci, upper_ci,
      postmean, upper_cred, lower_cred,
      alpha, beta
    )
  
  if (get_cv_stats == TRUE) {
    return(list(X_stats = X_out, cv_stats = cv_eb_out$cv_stats))
  }
  return(X_out)
}

#' Fit a beta-binomial to data via maximum likelihood optimization.
#' 
#' @param k (integer) number of edits for a gRNA.
#' @param n (integer) number of reads for a gRNA.
#' @param start (float, float) default c(0, 0), initial guesses for alpha, 
#'              beta parameters of beta-binomial.
#' @param method (string) default "L-BFGS-B", optimizer method used.
#' @return parameters for fitted beta-binomial
#' 
#' @keywords internal
#' @noRd
#'
#' @importFrom stats lchoose lbeta optim
fit_bb_mle <- function(k, n, start = c(0, 0), method = "L-BFGS-B") {
  bb_loglik <- function(par, k, n) {
    a <- exp(par[1])
    b <- exp(par[2])
    sum(
      lchoose(n, k) +
      lbeta(k + a, n - k + b) -
      lbeta(a, b)
    )
  }
  opt <- optim(
    start, bb_loglik, k = k, n = n,
    control = list(fnscale = -1), method = "L-BFGS-B"
  )
  list(
    alpha = exp(opt$par[1]), 
    beta = exp(opt$par[2]),
    convergence = opt$convergence,
    message = opt$message
  )
}

#' Compute posterior mean for data given k, n, and beta binomial prior
#' parameters.
#' 
#' @param k (integer) number of edits for a gRNA.
#' @param n (integer) number of reads for a gRNA.
#' @param alpha (float) alpha parameter for a beta-binomial prior.
#' @param beta (float) beta parameter for a beta-binomial prior.
#' @return posterior expectation for the beta-binomial given n and k.
#' 
#' @keywords internal
#' @noRd
post_mean <- function(k, n, alpha, beta) (k + alpha) / (n + alpha + beta)

#' Compute log-pmf of beta-binomial.
#'
#' @param k (integer) number of edits for a gRNA.
#' @param n (integer) number of reads for a gRNA.
#' @param alpha (float) alpha parameter for a beta-binomial prior.
#' @param beta (float) beta parameter for a beta-binomial prior.
#' @return log-pmf of a beta-binomial.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats lchoose lbeta
dbb_logpmf <- function(k, n, alpha, beta) {
  lchoose(n, k) +
    lbeta(k + alpha, n - k + beta) -
    lbeta(alpha, beta)
}

#' Computes empirical Bayesian statistics with cross-validation for cleaned
#' EMERGe data.
#' 
#' Fits a beta-binomial prior to EMERGe screening data using k-fold
#' cross-validation, and uses the fitted prior to compute posterior expectations
#' of editing rates with associated Bayesian credible intervals. For each fold,
#' the model is trained on the heldin data and evaluated on the held-out data
#' to assess fit quality.
#' 
#' The function returns both (i) per-guide posterior summaries (editing
#' posterior mean, lower and upper credible intervals, fitted prior parameters),
#' and (ii) cross-validation diagnostics (log-likelihood, mean negative
#' log-likelihood, AIC, Pearson chi-square, RMSE, MAE, and optimize convergence
#' status).
#' 
#' @param data (tibble) cleaned data with n10, n, and k columns.
#' @param num_folds (integer) default 5, number of folds to split data into
#'                  for cross-validation.
#' @param seed (float) default NULL.
#' @return A list with two components:
#'   \describe{
#'     \item{X_bayes}{A tibble containing n10, postmean, upper_cred, 
#'           lower_cred, alpha, beta}
#'     \item{cv_stats}{A tibble containing alpha MLE, beta MLE, number of
#'           samples in the held-in fold (test_n), log-likelihood of the
#'           test data (logLik), mean negative log-likelihood per observation
#'           (meanNegLL), Akaike Information Criterion for the fit (p = 2)
#'           (AIC), Pearson chi-square statistic comparing observed vs.
#'           expected counts (PearsonChi2), Root mean squared error between
#'           observed and expected counts (RMSE_counts), Mean absolute
#'           error between observed and expected counts (MAE_counts),
#'           Optimizer convergence code (0 indicates success) (convergence),
#'           and the Optimizer convergence message (message).}
#'   }
#' @return posterior expectation of editing, credible intervals about each
#'         expectation, and beta-binomial fit statistics.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr mutate select arrange filter group_nest left_join
#' @importFrom purrr map map2 pmap
#' @importFrom tibble tibble
#' @importFrom tidy group_nest
#' @importFrom stats qbeta
#' @importFrom dplyr any_of bind_rows
cv_eb <- function(X_clean, num_folds, seed = NULL) {
  stopifnot(all(c("n10", "k", "n") %in% names(X_clean)))
  if (!is.null(seed)) set.seed(seed)
  
  data <- data %>%
    mutate(fold = sample(rep(1:num_folds, length.out = nrow(.))))
  
  results <- X_clean %>%
    group_nest(fold, .key = "df") %>%
    mutate(
      train = map(fold, ~ filter(X_clean, fold != .x)),
      test  = map(fold, ~ filter(X_clean, fold == .x)),
      fit   = map(train, ~ fit_bb_mle(.x$k, .x$n)),
      
      stats = pmap(
        list(test=test, fit=fit),
        function(test, fit) {
          k <- test$k
          n <- test$n
          a <- fit$alpha
          b <- fit$beta
          N <- length(k)
          
          ll   <- sum(dbb_logpmf(k, n, a, b)) # test log-likelihood
          mnll <- -ll / max(1L, N)            # mean neg-Log-lik
          aic  <- -2*ll + 2*2                 # AIC (p=2 params)
          
          mu_k  <- n * (a/(a+b))
          var_k <- n * (a * b) / (a + b)^2 * ((a+b+n)/(a+b+1))
          var_k[var_k <= 0] <- 1e-12
          chi2  <- sum((k - mu_k)^2 / var_k)
          
          rmse  <- sqrt(mean((k - mu_k)^2))
          mae   <- mean(abs(k - mu_k))
          
          convergence <- fit$convergence
          message     <- fit$message
          
          tibble(
            alpha       = a,
            beta        = b,
            test_n      = N,
            logLik      = ll,
            meanNegLL   = mnll,
            AIC         = aic,
            PearsonChi2 = chi2,
            RMSE_counts = rmse,
            MAE_counts  = mae,
            convergence = convergence,
            message     = message
          )
        }),
      
      predictions = map2(test, fit, function(test, fit) {
        a <- fit$alpha
        b <- fit$beta
        tibble(
          n10 = test$n10,
          postmean = (test$k + a) / (test$n + a + b),
          lower_cred = qbeta(0.025, test$k + a, test$n - test$k + b),
          upper_cred = qbeta(0.975, test$k + a, test$n - test$k + b),
          alpha = a,
          beta = b
        )
      })
    )
  
  preds_all <- bind_rows(results$predictions, .id = "fold_id") %>%
    mutate(fold = as.integer(fold_id)) %>% select(-fold_id)
  
  stats_all <- bind_rows(results$stats, .id = "fold_id") %>%
    mutate(fold = as.integer(fold_id)) %>% select(-fold_id)
  
  X_bayes <- X_clean %>%
    select(-any_of(c(n, k, mle, lower_ci, upper_ci))) %>%
    left_join(preds_all, by="n10") %>%
    arrange(
      n10, postmean, upper_cred, lower_cred, alpha, beta
    )
  
  return(
    list(
      X_bayes = X_bayes,
      cv_stats = stats_all
    )
  )
}