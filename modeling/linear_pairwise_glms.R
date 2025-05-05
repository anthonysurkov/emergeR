library(tidyverse)
library(here)

library(pscl)
library(pROC)
library(glmnet)


setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("modeling/linear_pairwise_glms.R")

model_data <- read_csv(here("data", "R255X_E488QD_modeling.csv"))


# ============================================================================ #
# Try a normal linear fit                                                      #

model_data <- model_data %>%
  mutate(
    success = GGA,
    failure = n - GGA
  ) %>%
  select(-N10, -GAA, -GGA, -n, -n_other, -edit)

model_data[is.na(model_data)] <- 0 # some n = 0 while success > 0

bino_fit <- glm(
  cbind(success, failure) ~ .,
  data   = model_data,
  family = binomial(link = "logit")
)
print(summary(bino_fit))

calib_df <- model_data %>%
  mutate(
    obs_prop = success / (success + failure),
    pred_prop = predict(bino_fit, type="response")
  )

calib_df <- calib_df %>% drop_na()

cat("correlation:", cor(calib_df$obs_prop, calib_df$pred_prop), '\n',
    "Brier MSE :", mean((calib_df$pred_prop - calib_df$obs_prop)^2), '\n') 


# ============================================================================ #
# Identify important linear and pairwise features.                             #

X_main <- model_data %>%
  select(-success, -failure)
X_int <- model.matrix(
  ~ (.)^2 - 1,
  data = X_main
)

y <- cbind(model_data$success, model_data$failure)

set.seed(0)
cvfit <- cv.glmnet(
  x = X_int,
  y = y,
  family = "binomial",
  alpha = 0.5,
  standardize = TRUE,
  nlambda = 20,
  nfolds = 5,
  parallel = TRUE
)

fit_lasso <- cvfit$glmnet.fit
print(fit_lasso)

coefs <- coef(cvfit, s="lambda.1se")
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "coefficient"
features <- coefs_df %>%
  filter(abs(coefficient) > 1e-6) %>%
  arrange(desc(abs(coefficient)))
print(features)
