# ssn_fit_top_3.R
# Fits the 3 ArcGIS Exploratory Regression competing models (Delta AICc <= 2) as SSN models,
# holding the covariance structure fixed at the previously-selected qne justified specification:
# linear tail-up, no tail-down, jbessel euclid, nugget components: see ssn_covariance_investigate.R
# Only the fixed-effects formula differs across the 3 fits.

library(SSN2)
library(readr)
library(dplyr)

CRB_SSN <- readRDS("results_2021/exploratory_ssn/CRB_SSN_cleaned_top_3.RDS")
cat(sprintf("[CHECK 0] Loaded CRB_SSN_cleaned_top_3.RDS: %d sites\n", nrow(CRB_SSN$obs)))
stopifnot("CRB_SSN_cleaned_top_3.RDS does not have 72 sites as expected" = nrow(CRB_SSN$obs) == 72)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Three ArcGIS Exploratory Regression MLR model formulas
model_formulas <- list(
  Model1 = thermal_sensitivity ~ elevation_m + high_cascades_upstream_pct + shrub_upstream_pct + upstream_area_km2 + burned_buffer_pct,
  Model2 = thermal_sensitivity ~ elevation_m + base_flow_index + developed_upstream_pct + shrub_upstream_pct + upstream_area_km2,
  Model3 = thermal_sensitivity ~ channel_slope + base_flow_index + high_cascades_upstream_pct + shrub_upstream_pct + upstream_area_km2
)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Fixed covariance structure: identical across all 3 SSN fits matching the previously selected specification
fit_ssn_model <- function(formula, ssn_object) {
  ssn_lm(
    formula = formula,
    ssn.object = ssn_object,
    tailup_type = "linear",
    euclid_type = "jbessel",
    additive = "afvArea"
  )
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Fit all 3 MLR models to SSN, print convergence and covariance structure diagnostics for each
ssn_models <- list()

for (model_name in names(model_formulas)) {
  cat(sprintf("\n[CHECK 1] Fitting %s: %s\n", model_name, deparse(model_formulas[[model_name]])))

  fit <- tryCatch(
    fit_ssn_model(model_formulas[[model_name]], CRB_SSN),
    error = function(e) {
      cat(sprintf("  FAILED to converge: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) next
  ssn_models[[model_name]] <- fit

  g <- glance(fit)
  vc <- varcomp(fit)$proportion
  tu_coef <- coef(fit, type = "tailup")
  eu_coef <- coef(fit, type = "euclid")
  nug_coef <- coef(fit, type = "nugget")

  cat(sprintf("  Converged. AIC = %.4f | pseudo R2 = %.4f\n", g$AIC, g$pseudo.r.squared))
  cat(sprintf("  Tailup (linear): de = %.4f, range = %.3e\n", tu_coef[["de"]], tu_coef[["range"]]))
  cat(sprintf("  Euclid (jbessel): de = %.4f, range = %.3e\n", eu_coef[["de"]], eu_coef[["range"]]))
  cat(sprintf("  Nugget: %.4f\n", nug_coef[["nugget"]]))
  vc_named <- data.frame(
    component = c("covariates", "tailup", "taildown", "euclid", "nugget"),
    proportion = round(as.numeric(vc), 4)
  )
  cat("  Variance proportions:\n")
  print(vc_named)

  cat("  Fixed-effects coefficients:\n")
  print(summary(fit)$coefficients$fixed)

  saveRDS(fit, sprintf("results_2021/exploratory_ssn/ssn_%s.RDS", model_name))
}

cat(sprintf("\nDone. %d of %d models fit successfully and saved.\n", length(ssn_models), length(model_formulas)))