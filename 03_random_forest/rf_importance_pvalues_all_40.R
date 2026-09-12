# rf_importance_pvalues_all_40.R
# Altmann et al. (2010) permutation-based significance test for RF permutation importance,
# tested across ALL 40 covariates with Benjamini-Hochberg FDR correction applied across all 40.
#
# DESIGN CHANGE FROM PRIOR VERSIONS: earlier drafts corrected significance only within a
# pre-selected subset (the 6 MLR-selected covariates, or briefly, RF's own top-5). Correcting
# across all 40 is the only defensible choice when reporting RF's own top-ranked
# covariates: those covariates are selected BECAUSE they ranked highest on the same importance
# metric being tested, so correcting only within that pre-selected family understates the
# number of comparisons effectively made when the top set was chosen from all 40.
#
# PERMUTATIONS increased from 500 to 1999. A permutation test's minimum resolvable p-value is
# 1/(B+1) (Phipson & Smyth, 2010); at B=500 this floor is ~0.002, which sits ABOVE the nominal p-value 
# several borderline covariates need to survive 40-way BH-FDR correction (~0.00125 under a conservative 
# Bonferroni-like bound). B=1999 lowers the floor to ~0.0005 and separates covariates previously tied at the B=500 floor.
#
# IMPORTANCE SCALING: sum-normalized importance added (each covariate's raw importance divided
# by the sum of all 40 raw importances) so values are interpretable as each covariate's share of
# total importance and sum to 1 across the full set. Negative raw importances (possible with
# permutation importance) are left as negative shares rather than floored at zero, so the
# column still sums to exactly 1.
#
# Output: results_2021/rf/significance_testing/rf_importance_pvalues_all_40.csv

library(dplyr)
library(readr)
library(ranger)

data_2021 <- read_csv("data/clackamas_thermal_sensitivity_covariates_2021.csv", show_col_types = FALSE)
cat(sprintf("[CHECK 0] data_2021: %d rows\n", nrow(data_2021)))
stopifnot("data_2021 does not have 72 sites as expected" = nrow(data_2021) == 72)

optimal_hyperparams <- readRDS("results_2021/rf/optimal_rf_hyperparams.RDS")
cat("[CHECK 0] Loaded fixed hyperparameters:\n")
print(optimal_hyperparams)

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
stopifnot("covariates_40 does not have 40 covariates as expected" = length(covariates_40) == 40)

rfData <- data_2021 %>% select(thermal_sensitivity, all_of(covariates_40))
rf_formula <- as.formula(paste("thermal_sensitivity ~", paste(covariates_40, collapse = " + ")))

# ============================================================
# Altmann permutation p-values across all 40 covariates. Refits the forest once per
# permutation: runtime scales with num.permutations (~4x the earlier 500-permutation run).
# ============================================================
n_permutations <- 1999
cat(sprintf("[CHECK 1] Running Altmann permutation importance test (%d permutations)...\n", n_permutations))
cat("[CHECK 1] This refits the RF model once per permutation, so expect substantially longer\n")
cat("[CHECK 1] runtime than the earlier 500-permutation version. Consider running as a\n")
cat("[CHECK 1] background job rather than waiting interactively.\n")

set.seed(123)
altmann_result <- importance_pvalues(
  x = ranger(
    rf_formula, data = rfData,
    mtry = optimal_hyperparams$mtry,
    num.trees = optimal_hyperparams$ntrees,
    sample.fraction = optimal_hyperparams$sample_rate,
    min.node.size = optimal_hyperparams$min_rows,
    max.depth = optimal_hyperparams$max_depth,
    importance = "permutation",
    respect.unordered.factors = "order",
    seed = 123
  ),
  method = "altmann",
  num.permutations = n_permutations,
  formula = rf_formula,
  data = rfData
)

altmann_df <- as.data.frame(altmann_result)
altmann_df$covariate <- rownames(altmann_df)

cat("[CHECK 2] Raw Altmann importance/p-value results (all 40 covariates):\n")
print(altmann_df)

# ============================================================
# Sum-normalized importance + BH-FDR correction, both across all 40 covariates.
# ============================================================
altmann_df <- altmann_df %>%
  mutate(
    importance_normalized = importance / sum(importance),
    pvalue_fdr_bh = p.adjust(pvalue, method = "BH")
  ) %>%
  arrange(desc(importance))

cat("\n[CHECK 3] Full results: raw importance, normalized importance, raw p-value, and\n")
cat("[CHECK 3] BH-FDR corrected p-value (corrected across all 40 covariates):\n")
print(altmann_df %>% select(covariate, importance, importance_normalized, pvalue, pvalue_fdr_bh))

n_significant <- sum(altmann_df$pvalue_fdr_bh < 0.05)
cat(sprintf("\n[CHECK 4] %d of 40 covariates significant after BH-FDR correction (alpha = 0.05):\n", n_significant))
print(altmann_df %>% filter(pvalue_fdr_bh < 0.05) %>%
        select(covariate, importance, importance_normalized, pvalue_fdr_bh))

dir.create("results_2021/rf/significance_testing", recursive = TRUE, showWarnings = FALSE)
write_csv(altmann_df %>% select(covariate, importance, importance_normalized, pvalue, pvalue_fdr_bh),
          "results_2021/rf/significance_testing/rf_importance_pvalues_all_40.csv")

cat("\nSaved results_2021/rf/significance_testing/rf_importance_pvalues_all_40.csv.\n")