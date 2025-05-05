library(tidyverse)
library(here)


clean_data <- read_csv(here("data","R255X_E488QD_clean.csv"))
p_guess    <- sum(clean_data$GGA) / sum(clean_data$n)

alpha <- 0.05
z     <- qnorm(1 - alpha/2)

epsilons <- c(0.10, 0.08, 0.06, 0.05)  

r_hat  <- 0.6699527
mu_hat <- 62.9877483

tradeoff <- tibble(epsilon = epsilons) %>%
  mutate(
    n0         = ceiling(z^2 * p_guess * (1-p_guess) / epsilon^2),
    keep_frac  = 1 - pnbinom(n0-1, size = r_hat, mu = mu_hat),
    keep_pct   = keep_frac * 100
  )

print(tradeoff)
