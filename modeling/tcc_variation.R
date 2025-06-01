library(tidyverse)
library(here)
library(Matrix)
library(glmnet)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("modeling/tcc_variation.R")

indir  <- "data"
infile <- "R255X_E488QD_modeling.csv"


# ============================================================================ #
# TODO: Make this not terrible, I think I literally tried to guess motifs
# to see how much of the variation they explain. iirc something like 15% overall
# but that's weak

model_data <- read_csv(here(indir, infile))
model_data <- model_data %>% mutate(
  motif_AS1 = as.integer(substr(N10,5,10) == "TTGATT"),
  motif_AS2 = as.integer(substr(N10,6,10) == "TCCGC"),
  motif_AS3 = as.integer(substr(N10,6,10) == "TTGAC"),
  motif_AS4 = as.integer(substr(N10,6,9)  == "TCCA"),
  motif_AS5 = as.integer(substr(N10,5,7) == "TCC")
) %>% mutate(
  success = GGA,
  failure = n - GGA
) %>% 
  dplyr::select(-N10, -GAA, -GGA, -n, -n_other, -edit)

y_mat <- as.matrix(model_data %>% dplyr::select(success, failure))
w_vec <- with(model_data, success + failure)

X_df <- model_data %>% dplyr::select(-success, -failure)
X_mat <- model.matrix(~ ., data=X_df)[, -1]

X_sp <- as(Matrix(X_mat, sparse=TRUE), "dgCMatrix")

set.seed(1)
cvfit <- cv.glmnet(
  x = X_sp,
  y = y_mat,
  family = "binomial",
  weights = w_vec,
  alpha = 0.5,
  nfolds = 10,
  type.measure = "deviance"
)

plot(cvfit)
opt_lambda <- cvfit$lambda.min

nz_coef <- coef(cvfit, s="lambda.min")
nz_coef <- nz_coef[which(nz_coef != 0), , drop=FALSE]
print(nz_coef)

opt_idx   <- which(cvfit$lambda == cvfit$lambda.min)
pseudoR2  <- cvfit$glmnet.fit$dev.ratio[opt_idx]        # McFadden-style

cat("Pseudo-R² (glmnet dev.ratio):", round(pseudoR2, 4), "\n")
