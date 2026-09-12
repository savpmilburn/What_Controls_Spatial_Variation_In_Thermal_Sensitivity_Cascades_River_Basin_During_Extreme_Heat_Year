# ssn_covariance_investigate.R
# Re-investigates the covariance structure selection in the SSN models separately for each
# of the 3 ArcGIS Exploratory Regression MLR models; previous covariance structure is justified. 
# Previously: selected linear tail-up, no tail-down, jbessel euclidean, nugget components after
# rigorous testing for originally submitted SSN model using 5 covariates from MLR model. 
# Output: results_2021/exploratory_ssn/ssn_covariance_model_comparison_Model1.csv
#         results_2021/exploratory_ssn/ssn_covariance_model_comparison_Model2.csv
#         results_2021/exploratory_ssn/ssn_covariance_model_comparison_Model3.csv
#         results_2021/exploratory_ssn/ssn_covariance_top_structure_by_model.csv (summary: best structure per model)

library(SSN2)
library(dplyr)
library(readr)

CRB_SSN <- readRDS("results_2021/exploratory_ssn/CRB_SSN_cleaned_top_3.RDS")
cat(sprintf("[CHECK 0] Loaded CRB_SSN_cleaned_top_3.RDS: %d sites\n", nrow(CRB_SSN$obs)))
stopifnot("CRB_SSN_cleaned_top_3.RDS does not have 72 sites as expected" = nrow(CRB_SSN$obs) == 72)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3 ArcGIS Exploratory Regression model formulas
model_formulas <- list(
  Model1 = thermal_sensitivity ~ elevation_m + high_cascades_upstream_pct + shrub_upstream_pct + upstream_area_km2 + burned_buffer_pct,
  Model2 = thermal_sensitivity ~ elevation_m + base_flow_index + developed_upstream_pct + shrub_upstream_pct + upstream_area_km2,
  Model3 = thermal_sensitivity ~ channel_slope + base_flow_index + high_cascades_upstream_pct + shrub_upstream_pct + upstream_area_km2
)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Using previously tested components restrictions in originally submitted manuscript"
tailup_types <- c("mariah", "linear", "epa")
taildown_types <- c("mariah", "exponential", "linear", "none")
euclid_types <- c("jbessel", "cubic", "gaussian", "none")

covariance_grid <- expand.grid(tailup = tailup_types, taildown = taildown_types, euclid = euclid_types, stringsAsFactors = FALSE)
cat(sprintf("[CHECK 1] Testing %d covariance combinations x %d formulas = up to %d total fits\n", nrow(covariance_grid), length(model_formulas), nrow(covariance_grid) * length(model_formulas)))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Fit one covariance combination for one formula.
# Returns NULL on convergence failure
fit_one_combination <- function(formula, ssn_object, tu, td, eu) {
  args <- list(formula = formula, ssn.object = ssn_object, tailup_type = tu, additive = "afvArea")
  if (td != "none") args$taildown_type <- td
  if (eu != "none") args$euclid_type <- eu
  tryCatch(do.call(ssn_lm, args), error = function(e) NULL)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run full grid search separately for each of the 3 Exploratory Regression MLR formulas
all_results <- list()
top_structure_by_model <- data.frame(model = character(), best_structure = character(), AIC = numeric(), tailup = character(), taildown = character(), euclid = character(), stringsAsFactors = FALSE)

for (model_name in names(model_formulas)) {
  cat(sprintf("\n[CHECK 2] Running covariance grid search for %s...\n", model_name))
  formula <- model_formulas[[model_name]]

  model_comparison <- data.frame(model = character(), tailup = character(), taildown = character(),
                                  euclid = character(), AIC = numeric(), AICc = numeric(),
                                  pseudo_r2 = numeric(), tailup_de = numeric(), tailup_range = numeric(),
                                  taildown_de = numeric(), taildown_range = numeric(), euclid_de = numeric(),
                                  euclid_range = numeric(), nugget = numeric(), tailup_prop = numeric(),
                                  taildown_prop = numeric(), euclid_prop = numeric(), nugget_prop = numeric(),
                                  stringsAsFactors = FALSE)

  for (i in seq_len(nrow(covariance_grid))) {
    tu <- covariance_grid$tailup[i]
    td <- covariance_grid$taildown[i]
    eu <- covariance_grid$euclid[i]

    fit <- fit_one_combination(formula, CRB_SSN, tu, td, eu)
    if (is.null(fit)) {
      cat(sprintf("  Skipped (fit failed): tailup=%s, taildown=%s, euclid=%s\n", tu, td, eu))
      next
    }

    vc <- varcomp(fit)$proportion
    g <- glance(fit)

    tu_coef <- coef(fit, type = "tailup")
    td_coef <- if (td != "none") coef(fit, type = "taildown") else c(de = NA, range = NA)
    eu_coef <- if (eu != "none") coef(fit, type = "euclid") else c(de = NA, range = NA)
    nug_coef <- coef(fit, type = "nugget")

    model_comparison <- rbind(model_comparison, data.frame(
      model = paste(tu, td, eu, sep = "_"), tailup = tu, taildown = td, euclid = eu,
      AIC = g$AIC, AICc = g$AICc, pseudo_r2 = g$pseudo.r.squared,
      tailup_de = tu_coef[["de"]], tailup_range = tu_coef[["range"]],
      taildown_de = td_coef[["de"]], taildown_range = td_coef[["range"]],
      euclid_de = eu_coef[["de"]], euclid_range = eu_coef[["range"]],
      nugget = nug_coef[["nugget"]],
      tailup_prop = vc[2], taildown_prop = vc[3], euclid_prop = vc[4], nugget_prop = vc[5]
    ))
  }

  # Ranking selection criteria metrics for SSN models with corroborating metrics for each model
  # rank_aic: 1 = lowest AIC (best)
  # rank_pseudo_r2: 1 = highest pseudo R2 (best "overall variance explained")
  # rank_total_spatial_de: 1 = highest sum of tailup_de + taildown_de + euclid_de (best: "amount of variability explained by each specific covariance component")
  model_comparison$total_spatial_de <- rowSums(
    model_comparison[, c("tailup_de", "taildown_de", "euclid_de")], na.rm = TRUE
  )

  model_comparison$rank_aic <- rank(model_comparison$AIC, ties.method = "min")
  model_comparison$rank_pseudo_r2 <- rank(-model_comparison$pseudo_r2, ties.method = "min")
  model_comparison$rank_total_spatial_de <- rank(-model_comparison$total_spatial_de, ties.method = "min")

  # Flag implausible ranges where ranges are much larger than any plausible in-basin distance suggests
  # Threshold set at 1,000,000 map units
  RANGE_IMPLAUSIBLE_THRESHOLD <- 1e6
  model_comparison$range_flag <- with(model_comparison,
    (!is.na(tailup_range) & tailup_range > RANGE_IMPLAUSIBLE_THRESHOLD) |
    (!is.na(taildown_range) & taildown_range > RANGE_IMPLAUSIBLE_THRESHOLD) |
    (!is.na(euclid_range) & euclid_range > RANGE_IMPLAUSIBLE_THRESHOLD)
  )

  model_comparison <- model_comparison[order(model_comparison$rank_aic), ]
  rownames(model_comparison) <- NULL

  cat(sprintf("[CHECK 3] Top 5 covariance structures for %s (sorted by AIC, with corroborating ranking for other metrics):\n", model_name))
  print(head(model_comparison[, c("model", "AIC", "pseudo_r2", "total_spatial_de", "rank_aic", "rank_pseudo_r2", "rank_total_spatial_de","range_flag")], 5))

  write_csv(model_comparison, sprintf("results_2021/exploratory_ssn/ssn_covariance_model_comparison_%s.csv", model_name))
  all_results[[model_name]] <- model_comparison

  best_row <- model_comparison[1, ]
  top_structure_by_model <- rbind(top_structure_by_model, data.frame(
    model = model_name,
    best_structure = best_row$model,
    AIC = best_row$AIC,
    tailup = best_row$tailup,
    taildown = best_row$taildown,
    euclid = best_row$euclid,
    rank_pseudo_r2 = best_row$rank_pseudo_r2,
    rank_total_spatial_de = best_row$rank_total_spatial_de,
    range_flag = best_row$range_flag
  ))
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Summary of best covariance structure for 3 formulas specified by Exploratory Regression MLR models
cat("\n[CHECK 4] Best covariance structure per model:\n")
print(top_structure_by_model)

previously_selected <- c(tailup = "linear", taildown = "none", euclid = "jbessel")
top_structure_by_model$matches_previous_selection <- with(top_structure_by_model,
  tailup == previously_selected["tailup"] &
  taildown == previously_selected["taildown"] &
  euclid == previously_selected["euclid"]
)

cat("\n[CHECK 5] Does each model's best structure match the previously-selected structure",
    "(linear tailup, no taildown, jbessel euclid)?\n")
print(top_structure_by_model %>% select(model, best_structure, matches_previous_selection))

write_csv(top_structure_by_model, "results_2021/exploratory_ssn/ssn_covariance_top_structure_by_model.csv")
cat("\nSaved per-model comparison CSVs and ssn_covariance_top_structure_by_model.csv.\n")