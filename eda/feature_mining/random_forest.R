library(tidyverse)
library(here)
library(ranger)

set.seed(0)
setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("modeling/random_forest.R")

model_data <- read_csv(here("data", "R255X_E488QD_modeling.csv"))


# ============================================================================ #
# Random forest (rf) to identify important interactions/motifs

model_data <- model_data %>% drop_na()

X <- model_data %>%
  select(-N10, -GGA, -GAA, -n, -n_other, -edit)
y <- model_data$edit

rf_fit <- ranger(
  dependent.variable.name = "edit_rate",
  data = cbind(edit_rate = y, X),
  num.trees = 2000,
  mtry = floor(sqrt(ncol(X))),
  importance = "impurity",
  probability = FALSE,
)

importance_df <- data.frame(
  feature = names(rf_fit$variable.importance),
  importance = rf_fit$variable.importance
) %>%
  arrange(desc(importance))

print(importance_df)
