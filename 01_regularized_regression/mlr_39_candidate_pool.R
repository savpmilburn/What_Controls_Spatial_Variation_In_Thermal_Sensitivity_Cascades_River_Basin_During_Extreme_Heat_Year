# mlr_39_candidate_pool.R
# Ridge, LASSO, and elastic net regularized regression on original 39-covariate pool implemented via
# geographic 5-fold cross-validation.
# Tests whether continuous shrinkage (rather than explicit inclusion/exclusion) handles collinearity differently
# and helps validate MLR variable selection procedure.
# Not used for prediction: pooled R^2 is a diagnostic for whether unregularized OLS at p~n is unstable.

library(tidyverse)
library(glmnet)

set.seed(2021)  # cv.glmnet's internal CV fold assignment

data_2021 <- read.csv("data/clackamas_thermal_sensitivity_covariates_2021.csv", stringsAsFactors = FALSE)
nrow(data_2021)  # should be 72

mlr_candidate_covariates_2021 <- c(
  "channel_slope", "elevation_m", "base_flow_index", "summer_mean_max_air_temp_c",
  "summer_max_air_temp_c", "annual_precip_mm", "summer_precip_mm", "summer_mean_air_temp_c",
  "wet_season_precip_mm", "developed_upstream_pct", "lakes_upstream_pct",
  "agricultural_upstream_pct", "burned_upstream_pct", "road_density_upstream",
  "high_cascades_upstream_pct", "wetlands_upstream_pct", "veg_cover_upstream_pct",
  "veg_height_upstream_m", "forest_upstream_pct", "shrub_upstream_pct", "upstream_area_km2",
  "burned_reach_pct", "agricultural_reach_pct", "wetlands_reach_pct", "lakes_reach_pct",
  "high_cascades_reach_pct", "developed_reach_pct", "road_density_reach",
  "veg_cover_reach_pct", "veg_height_reach_m", "developed_buffer_pct",
  "agricultural_buffer_pct", "burned_buffer_pct", "wetlands_buffer_pct", "lakes_buffer_pct",
  "high_cascades_buffer_pct", "road_density_buffer", "veg_height_buffer_m", "veg_cover_buffer_pct"
)
length(mlr_candidate_covariates_2021)  # should be 39 (solar_exposure excluded due to instability)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Assign indices for sites ordered north to south and assign interleaved geographic folds
data_2021 <- data_2021 %>% arrange(index) %>% mutate(fold_number = paste0("Fold", ((index - 1) %% 5) + 1))
table(data_2021$fold_number)

drop_zero_variance <- function(train_data, covariates) {
  has_variance <- sapply(covariates, function(v) length(unique(na.omit(train_data[[v]]))) > 1)
  covariates[has_variance]
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Single parameterized fold-runner for ridge/LASSO/elastic net (alpha = 0/1/between).
run_glmnet_fold <- function(train_data, test_data, covariates, alpha_mix) {
  covariates_usable <- drop_zero_variance(train_data, covariates)
  x_train <- as.matrix(train_data[, covariates_usable]); y_train <- train_data$thermal_sensitivity
  x_test <- as.matrix(test_data[, covariates_usable]); y_test <- test_data$thermal_sensitivity

  cv_fit <- cv.glmnet(x_train, y_train, alpha = alpha_mix, standardize = TRUE)

  extract_coefs <- function(lambda_s) {
    coefs <- as.matrix(coef(cv_fit, s = lambda_s))
    setNames(as.numeric(coefs[-1, 1]), rownames(coefs)[-1])
  }
  coefs_min <- extract_coefs("lambda.min")
  coefs_1se <- extract_coefs("lambda.1se")

  # Standardized coefficients (train-fold SD, not full-data, to avoid leakage) 
  predictor_sds <- apply(x_train, 2, sd); response_sd <- sd(y_train)
  coefs_min_std <- coefs_min * predictor_sds / response_sd
  coefs_1se_std <- coefs_1se * predictor_sds / response_sd

  pred_min <- as.numeric(predict(cv_fit, newx = x_test, s = "lambda.min"))
  pred_1se <- as.numeric(predict(cv_fit, newx = x_test, s = "lambda.1se"))
  resid_min <- y_test - pred_min
  resid_1se <- y_test - pred_1se

  list(
    lambda_min = cv_fit$lambda.min, lambda_1se = cv_fit$lambda.1se,
    coefs_min = coefs_min, coefs_1se = coefs_1se,
    coefs_min_std = coefs_min_std, coefs_1se_std = coefs_1se_std,
    resid_min = resid_min, resid_1se = resid_1se,
    rmse_min = sqrt(mean(resid_min^2)), rmse_1se = sqrt(mean(resid_1se^2))
  )
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
run_ols_fold <- function(train_data, test_data, covariates) {
  covariates_usable <- drop_zero_variance(train_data, covariates)
  model <- lm(as.formula(paste("thermal_sensitivity ~", paste(covariates_usable, collapse = " + "))), data = train_data)
  residuals <- test_data$thermal_sensitivity - predict(model, newdata = test_data)
  list(residuals = residuals, rmse = sqrt(mean(residuals^2)))
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run ridge, LASSO, elastic net (alpha=0.5), and OLS diagnostic across all 5 folds.
fold_ids <- paste0("Fold", 1:5)
methods <- list(ridge = 0, lasso = 1, elastic_net = 0.5)
results_by_method <- list()

for (method_name in names(methods)) {
  cat(sprintf("=== %s ===\n", method_name))
  fold_results <- list()
  for (fold_id in fold_ids) {
    train_fold <- filter(data_2021, fold_number != fold_id)
    test_fold <- filter(data_2021, fold_number == fold_id)
    fold_results[[fold_id]] <- run_glmnet_fold(train_fold, test_fold, mlr_candidate_covariates_2021, methods[[method_name]])
  }
  results_by_method[[method_name]] <- fold_results
}

cat("=== OLS diagnostic (unregularized, all 39 covariates) ===\n")
ols_resid_all <- c()
for (fold_id in fold_ids) {
  train_fold <- filter(data_2021, fold_number != fold_id)
  test_fold <- filter(data_2021, fold_number == fold_id)
  ols_resid_all <- c(ols_resid_all, run_ols_fold(train_fold, test_fold, mlr_candidate_covariates_2021)$residuals)  
}
sst <- sum((data_2021$thermal_sensitivity - mean(data_2021$thermal_sensitivity))^2)
ols_pooled_r2 <- 1 - sum(ols_resid_all^2) / sst
cat(sprintf("OLS pooled R2: %.4f\n\n", ols_pooled_r2))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Output long-format master coefficient table + performance summary.
build_long_table <- function(fold_results, method_name) {
  do.call(rbind, lapply(names(fold_results), function(f) {
    r <- fold_results[[f]]
    rbind(
      data.frame(fold = f, method = method_name, lambda_type = "min", lambda_value = r$lambda_min,
                 covariate = names(r$coefs_min), coefficient = as.numeric(r$coefs_min),
                 coefficient_standardized = as.numeric(r$coefs_min_std),
                 is_nonzero = as.numeric(r$coefs_min) != 0, rmse_fold = r$rmse_min),
      data.frame(fold = f, method = method_name, lambda_type = "1se", lambda_value = r$lambda_1se,
                 covariate = names(r$coefs_1se), coefficient = as.numeric(r$coefs_1se),
                 coefficient_standardized = as.numeric(r$coefs_1se_std),
                 is_nonzero = as.numeric(r$coefs_1se) != 0, rmse_fold = r$rmse_1se)
    )
  }))
}

master_coefficients <- do.call(rbind, lapply(names(results_by_method), function(m) build_long_table(results_by_method[[m]], m)))
master_coefficients$pool <- "39_covariate"
write_csv(master_coefficients, "results_2021/regularized_regression/master_coefficients.csv")

performance_rows <- lapply(names(results_by_method), function(m) {
  r <- results_by_method[[m]]
  resid_min_all <- unlist(lapply(r, function(x) x$resid_min))
  resid_1se_all <- unlist(lapply(r, function(x) x$resid_1se))
  data.frame(
    method = m,
    lambda_type = c("min", "1se"),
    pooled_r2 = c(
      1 - sum(resid_min_all^2) / sst,
      1 - sum(resid_1se_all^2) / sst
    )
  )
})

performance_summary <- bind_rows(performance_rows) %>%
  bind_rows(data.frame(method = "ols_unregularized", lambda_type = NA, pooled_r2 = ols_pooled_r2))
performance_summary$pool <- "39_covariate"
print(performance_summary)
write_csv(performance_summary, "results_2021/regularized_regression/performance_summary.csv")

retention_table <- master_coefficients %>%
  filter(lambda_type == "min") %>%
  group_by(method, covariate) %>%
  summarise(pct_retained = 100 * sum(is_nonzero) / length(fold_ids), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = pct_retained, values_fill = 0) %>%
  arrange(desc(elastic_net), desc(lasso))

write_csv(retention_table, "results_2021/regularized_regression/retention_table_S1_2.csv")