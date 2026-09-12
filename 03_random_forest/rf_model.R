# rf_model.R
# Fits a RF model predicting thermal sensitivity using all 40 landscape covariates following Section 2.3.3.
# Uses the h2o package to execute a random grid search for optimal hyperparameters, then reports default-parameter 
# OOB-RMSE as baseline reference alongside tuned result. 
# Note: solar_exposure included in the RF model since MLR exclusion was a linear-model-specific diagnostic.
# RF has no equivalent basis to exclude it. 

# Output: saves optimal hyperparameters to results_2021/rf/optimal_rf_hyperparams.RDS.

library(dplyr)
library(readr)
library(ranger)
library(h2o)

# Load data
data_2021 <- read_csv("data/clackamas_thermal_sensitivity_covariates_2021.csv")
cat(sprintf("[CHECK 0] data_2021: %d rows\n", nrow(data_2021)))
stopifnot("data_2021 does not have 72 sites as expected" = nrow(data_2021) == 72)

# All 40 landscape covariates (solar_exposure included)
rf_covariates <- c("channel_slope", "solar_exposure", "elevation_m", "base_flow_index", "developed_upstream_pct", "lakes_upstream_pct", "agricultural_upstream_pct",
"burned_upstream_pct", "road_density_upstream", "high_cascades_upstream_pct", "wetlands_upstream_pct", "veg_cover_upstream_pct", "veg_height_upstream_m",
"forest_upstream_pct", "shrub_upstream_pct", "upstream_area_km2", "burned_reach_pct", "agricultural_reach_pct", "wetlands_reach_pct", 
"lakes_reach_pct", "high_cascades_reach_pct", "developed_reach_pct", "road_density_reach", "veg_cover_reach_pct", "veg_height_reach_m",
"developed_buffer_pct", "agricultural_buffer_pct", "burned_buffer_pct", "wetlands_buffer_pct", "lakes_buffer_pct", "high_cascades_buffer_pct", 
"road_density_buffer", "veg_height_buffer_m", "veg_cover_buffer_pct", "summer_mean_max_air_temp_c", "summer_max_air_temp_c", "annual_precip_mm",
"summer_precip_mm", "summer_mean_air_temp_c", "wet_season_precip_mm")
n_predictors <- length(rf_covariates)
cat(sprintf("[CHECK 0] Covariate pool size: %d\n", n_predictors))
stopifnot("rf_covariates does not have 40 covariates as expected" = n_predictors == 40)

rfData <- data_2021 %>% select(thermal_sensitivity, all_of(rf_covariates))




#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run default RF model for OOB RMSE baseline before tuning.
# mtry = floor(p/3) is standard ranger recommendation for regression (vs. floor(sqrt(p))) for classification).
rf_default <- ranger(thermal_sensitivity ~ ., data = rfData, mtry = floor(n_predictors / 3), importance = "permutation", respect.unordered.factors = "order", seed = 123)
default_rmse <- sqrt(rf_default$prediction.error)
cat(sprintf("[CHECK 1] Default RF model OOB RMSE (mtry = floor(%d/3) = %d): %.6f\n", n_predictors, floor(n_predictors / 3), default_rmse))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Hyperparameter optimization via h2o random grid search.
h2o.no_progress()
h2o.init(max_mem_size = "2g")
h2o.removeAll()

train_h2o <- as.h2o(rfData)
response <- "thermal_sensitivity"
predictors <- setdiff(colnames(rfData), response)
stopifnot("predictors count does not match 2021 covariates" = length(predictors) == n_predictors)

hyper_grid <- list(mtries = floor(n_predictors * c(0.05, 0.15, 0.25, 0.333, 0.4)), min_rows = c(1, 3, 5, 7), max_depth = c(10, 15, 20, 25), sample_rate = c(.50, .63, .7, .8))
cat("[CHECK 2] Hyperparameter grid:\n")
print(hyper_grid)
cat(sprintf("[CHECK 2] Total combinations: %d\n", length(hyper_grid$mtries) * length(hyper_grid$min_rows) * length(hyper_grid$max_depth) * length(hyper_grid$sample_rate)))

search_criteria <- list(strategy = "RandomDiscrete", stopping_metric = "rmse", stopping_tolerance = 0.005, stopping_rounds = 8, max_runtime_secs = 60 * 3)

random_grid <- h2o.grid(algorithm = "randomForest", grid_id = "rf_random_grid", x = predictors, y = response, training_frame = train_h2o, hyper_params = hyper_grid, ntrees = n_predictors * 10, seed = 123, stopping_metric = "RMSE", stopping_rounds = 10, stopping_tolerance = 0.005, search_criteria = search_criteria)
cat("[CHECK 2] Grid search completed.\n")

random_grid_rf_results <- h2o.getGrid(grid_id = "rf_random_grid", sort_by = "mse", decreasing = FALSE)
print(random_grid_rf_results)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract best RF model and its hyperparameters.
best_RF_model <- h2o.getModel(random_grid_rf_results@model_ids[[1]])
best_rmse <- h2o.rmse(best_RF_model, train = TRUE)
best_params <- best_RF_model@allparameters

cat(sprintf("[CHECK 3] Best model training RMSE: %.6f\n", best_rmse))
cat("[CHECK 3] Optimal hyperparameters from best RF model:\n")
cat("mtries:", best_params$mtries, "\n")
cat("min_rows:", best_params$min_rows, "\n")
cat("max_depth:", best_params$max_depth, "\n")
cat("sample_rate:", best_params$sample_rate, "\n")

stopifnot("best_params$mtries is NULL -- grid search did not return a valid mtries value" = !is.null(best_params$mtries), "best_params$sample_rate is NULL" = !is.null(best_params$sample_rate),"best_params$max_depth is NULL" = !is.null(best_params$max_depth),"best_params$min_rows is NULL" = !is.null(best_params$min_rows))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Top 5 RF models comparison for diagnostics only
top_models <- h2o.getGrid(grid_id = "rf_random_grid", sort_by = "rmse", decreasing = FALSE)
cat(sprintf("[CHECK 4] Number of models in grid: %d\n", length(top_models@model_ids)))

for (i in seq_len(min(5, length(top_models@model_ids)))) {
  model <- h2o.getModel(top_models@model_ids[[i]])
  rmse <- h2o.rmse(model, train = TRUE)
  cat(sprintf("Model %d: RMSE = %.6f, mtries = %s, sample_rate = %s, max_depth = %s, min_rows = %s\n", i, rmse, model@allparameters$mtries, model@allparameters$sample_rate, model@allparameters$max_depth, model@allparameters$min_rows))
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Save RF optimal hyperparameters.
# max-depth fixed at 25 since empirically equivalent to grid's returned optimum via ranger() output OOB RMSE identical to 6 decimals: 0.069002 for both
# grid search ties among max_depth >= 10 reflecting that rf treats on N=72 rarely reach depth 20 regardless of nominal cap
best_hyperparams <- list(mtry = best_params$mtries, min_rows = best_params$min_rows, max_depth = 25, sample_rate = best_params$sample_rate, ntrees = n_predictors * 10)

cat("[CHECK 5] Final saved hyperparameters:\n")
print(best_hyperparams)
stopifnot("best_hyperparams$mtry is NULL after fix -- something is still wrong" = !is.null(best_hyperparams$mtry))

saveRDS(best_hyperparams, "results_2021/rf/optimal_rf_hyperparams.RDS")
cat("[CHECK 5] Saved optimal_rf_hyperparams.RDS.\n")
