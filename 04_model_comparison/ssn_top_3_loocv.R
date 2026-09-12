# ssn_top_3_loocv.R
# Runs LOOCV on the 3 SSN model refits (see ssn_fit_top_3.R)

library(SSN2)
library(readr)
library(dplyr)

model_names <- c("Model1", "Model2", "Model3")
n_predictors_by_model <- c(Model1 = 5, Model2 = 5, Model3 = 5)

results_list <- list()

for (model_name in model_names) {
  cat(sprintf("\n[CHECK 0] Loading ssn_%s.RDS...\n", model_name))
  fit <- readRDS(sprintf("results_2021/exploratory_ssn/ssn_%s.RDS", model_name))
  cat(sprintf("[CHECK 0] Loaded (n = %d)\n", fit$n))

  cat(sprintf("[CHECK 1] Running LOOCV for %s...\n", model_name))
  loocv_results <- loocv(fit)
  print(loocv_results)

  ssn_aic <- AIC(fit)

  ssn_resid <- residuals(fit)
  ssn_train_rmse <- sqrt(mean(ssn_resid^2))
  ssn_fitted <- fitted(fit)
  ssn_observed <- ssn_fitted + ssn_resid
  ssn_train_r2 <- 1 - sum(ssn_resid^2) / sum((ssn_observed - mean(ssn_observed))^2)

  cat(sprintf("[CHECK 2] %s -- AIC: %.4f | training RMSE: %.4f | training R2: %.4f | LOOCV RMSE: %.4f | LOOCV R2: %.4f\n",
              model_name, ssn_aic, ssn_train_rmse, ssn_train_r2, loocv_results$RMSPE, loocv_results$cor2))

  results_list[[model_name]] <- data.frame(
    model = paste0("SSN_", model_name),
    n_predictors = n_predictors_by_model[[model_name]],
    aic = round(ssn_aic, 4),
    aic_type = "Exact Gaussian AIC (stats::AIC on full-sample ssn_lm ML/REML fit)",
    training_rmse = round(ssn_train_rmse, 4),
    training_r2 = round(ssn_train_r2, 4),
    loocv_rmse = round(loocv_results$RMSPE, 4),
    loocv_r2 = round(loocv_results$cor2, 4)
  )
}

ssn_three_models_summary <- bind_rows(results_list)
cat("\n[CHECK 3] Combined SSN LOOCV summary (Table 1 format, 3 competing models):\n")
print(ssn_three_models_summary)

write_csv(ssn_three_models_summary, "results_2021/exploratory_ssn/ssn_top_3_loocv_results.csv")
cat("\nSaved ssn_top_3_loocv_results.csv.\n")