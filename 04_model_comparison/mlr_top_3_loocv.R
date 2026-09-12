# mlr_top_3_loocv.R
# Leave-one-out CV for top-performing MLR model and 2 competing models from ArcGIS Exploratory Regression. 
#
# Output: results_2021/exploratory_mlr/mlr_top_3_loocv_results.csv

library(dplyr)
library(readr)

data_2021 <- read_csv("data/clackamas_thermal_sensitivity_covariates_2021.csv", show_col_types = FALSE)
cat(sprintf("[CHECK 0] data_2021: %d rows\n", nrow(data_2021)))
stopifnot("data_2021 does not have 72 sites as expected" = nrow(data_2021) == 72)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Three MLR model formulas for top-performing model and two competing models from ArcGIS Exploratory Regression
model_formulas <- list(
  Model1 = thermal_sensitivity ~ elevation_m + high_cascades_upstream_pct + shrub_upstream_pct + upstream_area_km2 + burned_buffer_pct,
  Model2 = thermal_sensitivity ~ elevation_m + base_flow_index + developed_upstream_pct + shrub_upstream_pct + upstream_area_km2,
  Model3 = thermal_sensitivity ~ channel_slope + base_flow_index + high_cascades_upstream_pct + shrub_upstream_pct + upstream_area_km2
)

results_list <- list()

for (model_name in names(model_formulas)) {
  formula_i <- model_formulas[[model_name]]
  covariates_i <- all.vars(formula_i)[-1]  # drop response variable

  cat(sprintf("\n[CHECK 1] Fitting %s: %s\n", model_name, deparse(formula_i)))

  modelData <- data_2021 %>% select(thermal_sensitivity, all_of(covariates_i))

  # Full-sample fit used for AIC and for training RMSE and R-squared references
  full_model <- lm(formula_i, data = modelData)
  model_aic <- AIC(full_model)
  train_rmse <- sqrt(mean(residuals(full_model)^2))
  train_r2 <- summary(full_model)$r.squared

  cat(sprintf("[CHECK 1] %s full-sample AIC: %.4f | training RMSE: %.4f | training R2: %.4f\n", model_name, model_aic, train_rmse, train_r2))

  # LOOCV loop 
  n <- nrow(modelData)
  loo_residuals <- numeric(n)

  for (i in seq_len(n)) {
    train_i <- modelData[-i, ]
    test_i <- modelData[i, , drop = FALSE]

    model_i <- lm(formula_i, data = train_i)
    pred_i <- predict(model_i, newdata = test_i)

    loo_residuals[i] <- test_i$thermal_sensitivity - pred_i
  }

  loocv_rmse <- sqrt(mean(loo_residuals^2))
  loocv_ssr <- sum(loo_residuals^2)
  loocv_sst <- sum((modelData$thermal_sensitivity - mean(modelData$thermal_sensitivity))^2)
  loocv_r2 <- 1 - loocv_ssr / loocv_sst

  cat(sprintf("[CHECK 2] %s LOOCV RMSE: %.4f | LOOCV R2: %.4f\n", model_name, loocv_rmse, loocv_r2))

  results_list[[model_name]] <- data.frame(
    model = paste0("MLR_", model_name),
    n_predictors = length(covariates_i),
    aic = round(model_aic, 4),
    aic_type = "Exact Gaussian AIC (stats::AIC on full-sample lm fit)",
    training_rmse = round(train_rmse, 4),
    training_r2 = round(train_r2, 4),
    loocv_rmse = round(loocv_rmse, 4),
    loocv_r2 = round(loocv_r2, 4)
  )
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Combined summary table 
mlr_three_models_summary <- bind_rows(results_list)
cat("\n[CHECK 3] Combined MLR LOOCV summary (Table 1 format, 3 competing models):\n")
print(mlr_three_models_summary)

write_csv(mlr_three_models_summary, "results_2021/exploratory_mlr/mlr_top_3_loocv_results.csv")
cat("\nSaved mlr_top_3_loocv_results.csv.\n")