# rf_loocv.R
# Leave-one-out CV of the fixed, final RF model: fixed hyperparameters from
# optimal_rf_hyperparams.RDS on all 40 covariates becuase RF has no MLR-style collinearity
# exclusion basis. To enable direct comparison with SSN LOOCV.
# Output: results_2021/rf/rf_loocv_results.csv

library(dplyr)
library(readr)
library(ranger)

data_2021 <- read_csv("data/clackamas_thermal_sensitivity_covariates_2021.csv")
cat(sprintf("[CHECK 0] data_2021: %d rows\n", nrow(data_2021)))
stopifnot("data_2021 does not have 72 sites as expected" = nrow(data_2021) == 72)

optimal_hyperparams <- readRDS("results_2021/rf/optimal_rf_hyperparams.RDS")
cat("[CHECK 0] Loaded fixed (non-refit) hyperparameters for LOOCV:\n")
print(optimal_hyperparams)
stopifnot(
  "optimal_hyperparams$mtry is NULL" = !is.null(optimal_hyperparams$mtry),
  "optimal_hyperparams$max_depth is NULL" = !is.null(optimal_hyperparams$max_depth),
  "optimal_hyperparams$sample_rate is NULL" = !is.null(optimal_hyperparams$sample_rate),
  "optimal_hyperparams$min_rows is NULL" = !is.null(optimal_hyperparams$min_rows)
)

# Same 40-covariate pool as rf_model.R / optimal_rf_model.R (RF's full model: no MLR-style reduction applies)
covariates_40 <- c("channel_slope", "solar_exposure", "elevation_m", "base_flow_index",
             "developed_upstream_pct", "lakes_upstream_pct", "agricultural_upstream_pct",
             "burned_upstream_pct", "road_density_upstream", "high_cascades_upstream_pct",
             "wetlands_upstream_pct", "veg_cover_upstream_pct", "veg_height_upstream_m",
             "forest_upstream_pct", "shrub_upstream_pct", "upstream_area_km2",
             "burned_reach_pct", "agricultural_reach_pct", "wetlands_reach_pct",
             "lakes_reach_pct", "high_cascades_reach_pct", "developed_reach_pct",
             "road_density_reach", "veg_cover_reach_pct", "veg_height_reach_m",
             "developed_buffer_pct", "agricultural_buffer_pct", "burned_buffer_pct",
             "wetlands_buffer_pct", "lakes_buffer_pct", "high_cascades_buffer_pct",
             "road_density_buffer", "veg_height_buffer_m", "veg_cover_buffer_pct",
             "summer_mean_max_air_temp_c", "summer_max_air_temp_c", "annual_precip_mm",
             "summer_precip_mm", "summer_mean_air_temp_c", "wet_season_precip_mm")
n_predictors <- length(covariates_40)
stopifnot("covariates_40 does not have 40 covariates as expected" = n_predictors == 40)

rfData <- data_2021 %>% select(thermal_sensitivity, all_of(covariates_40))
n <- nrow(rfData)

# ----
# Full-sample fit kept only for an OOB RMSE/R2 reference row alongside the LOOCV numbers, no AIC computed
rf_full_model <- ranger(
  thermal_sensitivity ~ ., data = rfData,
  mtry = optimal_hyperparams$mtry,
  num.trees = optimal_hyperparams$ntrees,
  sample.fraction = optimal_hyperparams$sample_rate,
  min.node.size = optimal_hyperparams$min_rows,
  max.depth = optimal_hyperparams$max_depth,
  importance = "none",
  respect.unordered.factors = "order",
  seed = 123
)
cat(sprintf("[CHECK 1] RF OOB RMSE (full-sample fixed model): %.4f\n", sqrt(rf_full_model$prediction.error)))
cat(sprintf("[CHECK 1] RF OOB R-squared (full-sample fixed model): %.4f\n", rf_full_model$r.squared))

# LOOCV: 72 refits with fixed hyperparameters (no retuning per fold)
loo_residuals <- numeric(n)

cat("[CHECK 2] Running 72 leave-one-out RF refits (fixed hyperparameters)...\n")
for (i in seq_len(n)) {
  train_i <- rfData[-i, ]
  test_i <- rfData[i, , drop = FALSE]

  model_i <- ranger(
    thermal_sensitivity ~ ., data = train_i,
    mtry = optimal_hyperparams$mtry,
    num.trees = optimal_hyperparams$ntrees,
    sample.fraction = optimal_hyperparams$sample_rate,
    min.node.size = optimal_hyperparams$min_rows,
    max.depth = optimal_hyperparams$max_depth,
    importance = "none",
    respect.unordered.factors = "order",
    seed = 123
  )

  pred_i <- predict(model_i, data = test_i)$predictions
  loo_residuals[i] <- test_i$thermal_sensitivity - pred_i

  if (i %% 10 == 0) cat(sprintf("  ...%d/%d complete\n", i, n))
}

rf_loocv_rmse <- sqrt(mean(loo_residuals^2))
rf_loocv_ssr <- sum(loo_residuals^2)
rf_loocv_sst <- sum((rfData$thermal_sensitivity - mean(rfData$thermal_sensitivity))^2)
rf_loocv_r2 <- 1 - rf_loocv_ssr / rf_loocv_sst

cat(sprintf("[CHECK 3] RF LOOCV RMSE: %.4f\n", rf_loocv_rmse))
cat(sprintf("[CHECK 3] RF LOOCV R-squared: %.4f\n", rf_loocv_r2))

#-------
# Save results for Table 3 (no AIC column)
rf_loocv_results <- data.frame(
  model = "RF",
  n_predictors = n_predictors,
  training_rmse = round(sqrt(rf_full_model$prediction.error), 4),  # OOB, not in-sample -- see note below
  training_r2 = round(rf_full_model$r.squared, 4),                  # OOB R2, not in-sample -- see note below
  loocv_rmse = round(rf_loocv_rmse, 4),
  loocv_r2 = round(rf_loocv_r2, 4)
)
# Note on training_rmse/training_r2 for RF: these are OOB estimates from the full-sample
# fixed model, NOT in-sample fit -- in-sample RF fit is near-perfect by construction
# (min.node.size = 1) and would not be a meaningful "training" reference alongside MLR's
# true in-sample training RMSE/R2.

cat("[CHECK 4] RF LOOCV summary row:\n")
print(rf_loocv_results)

write_csv(rf_loocv_results, "results_2021/rf/loocv_results.csv")
cat("\nSaved rf/loocv_results.csv.\n")