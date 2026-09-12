# optimal_rf_model.R
# Fits the final RF model on the complete dataset using the optimal hyperparameters found
# in rf_model.R and computes permutation importance (Section 2.3.4's reported importance ranking, 

library(dplyr)
library(readr)
library(ranger)

data_2021 <- read_csv("data/clackamas_thermal_sensitivity_covariates_2021.csv")
cat(sprintf("[CHECK 0] data_2021: %d rows\n", nrow(data_2021)))
stopifnot("data_2021 does not have 72 sites as expected" = nrow(data_2021) == 72)

optimal_hyperparams <- readRDS("results_2021/rf/optimal_rf_hyperparams.RDS")
cat("[CHECK 0] Loaded optimal hyperparameters:\n")
print(optimal_hyperparams)
stopifnot(
  "optimal_hyperparams$mtry is NULL" = !is.null(optimal_hyperparams$mtry),
  "optimal_hyperparams$max_depth is NULL" = !is.null(optimal_hyperparams$max_depth),
  "optimal_hyperparams$sample_rate is NULL" = !is.null(optimal_hyperparams$sample_rate),
  "optimal_hyperparams$min_rows is NULL" = !is.null(optimal_hyperparams$min_rows)
)

# Same 40-covariate pool as rf_model.R (solar_exposure included -- no MLR-style exclusion
# basis applies to RF; see rf_model.R header note).
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
stopifnot("40 covariates does not have 40 covariates as expected" = length(covariates_40) == 40)

rfData <- data_2021 %>% select(thermal_sensitivity, all_of(covariates_40))

# ============================================================
# Fit the final model on the complete dataset (all 72 sites) using the optimal
# hyperparameters found via HPO in rf_model.R.
# ============================================================
rf_final_model <- ranger(
  thermal_sensitivity ~ ., data = rfData,
  mtry = optimal_hyperparams$mtry,
  num.trees = optimal_hyperparams$ntrees,
  sample.fraction = optimal_hyperparams$sample_rate,
  min.node.size = optimal_hyperparams$min_rows,
  max.depth = optimal_hyperparams$max_depth,
  importance = "permutation",
  respect.unordered.factors = "order",
  seed = 123
)

cat(sprintf("[CHECK 1] Final RF model OOB RMSE: %.6f\n", sqrt(rf_final_model$prediction.error)))
cat(sprintf("[CHECK 1] Final RF model OOB R-squared: %.4f\n", rf_final_model$r.squared))

# ============================================================
# Permutation importance -- full ranking (all 40 covariates), plus top 5.
# ============================================================
importance_scores <- ranger::importance(rf_final_model)
importance_df <- data.frame(covariate = names(importance_scores), importance = as.numeric(importance_scores)) %>%
  arrange(desc(importance))

cat("[CHECK 2] Full permutation importance ranking:\n")
print(importance_df)

top5_covariates <- importance_df$covariate[1:5]
cat(sprintf("\n[CHECK 2] Top 5 most important covariates: %s\n", paste(top5_covariates, collapse = ", ")))

dir.create("results_2021/rf", recursive = TRUE, showWarnings = FALSE)
write_csv(importance_df, "results_2021/rf/rf_final_model_importance.csv")
saveRDS(rf_final_model, "results_2021/rf/rf_final_model.RDS")
cat("\nSaved rf_final_model.RDS and rf_final_model_importance.csv.\n")