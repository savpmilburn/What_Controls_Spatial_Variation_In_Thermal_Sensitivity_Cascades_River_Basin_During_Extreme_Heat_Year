# rf_ale_plots_rf_significant.R
# Accumulated Local Effects (ALE) plots for the 6 RF covariates significant after
# Benjamini-Hochberg FDR correction across all 40 candidates (Altmann et al., 2010
# permutation test at B=1999; see rf_importance_pvalues_all_40.R).
#
# After the all-40 significance test was run, agricultural_upstream_pct
# (pvalue_fdr_bh = 0.0067) was added here as the 6th significant covariate; the other 5
# (base_flow_index, high_cascades_upstream_pct, veg_height_upstream_m, elevation_m,
# developed_upstream_pct) are the other 5 significant covariates.
#
# ALE is used instead of PDPs because PDPs sweep a covariate across its full range while
# holding all other covariates at observed values, which can evaluate the model on covariate
# combinations with little to no real data support (confirmed with lakes_upstream_pct, where
# the tail beyond the 65th percentile was driven by only 2/72 sites). ALE effects are computed
# from differences in predictions within local, quantile-based neighborhoods of the observed
# data, rather than by substituting a swept value into every observation regardless of whether
# that combination occurs in the data (Apley & Zhu, 2020, JRSS-B).

library(dplyr)
library(readr)
library(ranger)
library(iml)
library(ggplot2)
library(patchwork)
library(extrafont)
loadfonts(quiet = TRUE)

data_2021 <- read_csv("data/clackamas_thermal_sensitivity_covariates_2021.csv", show_col_types = FALSE)
cat(sprintf("Check 0: data_2021: %d rows\n", nrow(data_2021)))
stopifnot("data_2021 does not have 72 sites as expected" = nrow(data_2021) == 72)

rf_final_model <- readRDS("results_2021/rf/rf_final_model.RDS")

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

# The 6 covariates significant after BH-FDR correction across all 40 (alpha = 0.05).
rf_significant_covariates <- c("base_flow_index", "high_cascades_upstream_pct", "veg_height_upstream_m", "elevation_m","developed_upstream_pct", "agricultural_upstream_pct")
stopifnot("rf_significant_covariates does not have 6 covariates as expected" = length(rf_significant_covariates) == 6)

# Full display names for axis/title labels (not raw variable names).
covariate_full_names <- c(
  base_flow_index            = "Base flow index (BFI)",
  high_cascades_upstream_pct = "High Cascades geology % (upstream)",
  veg_height_upstream_m      = "Vegetation height (upstream)",
  elevation_m                = "Elevation",
  developed_upstream_pct     = "Developed land area % (upstream)",
  agricultural_upstream_pct  = "Agricultural land area % (upstream)"
)

# Background colors for each covariate.
covariate_colors <- c(
  base_flow_index             = "#DCEEF2",
  high_cascades_upstream_pct  = "#F5E0E3",
  veg_height_upstream_m       = "#F5F0DC",
  elevation_m                 = "#E8E4F0",
  developed_upstream_pct      = "#EAEAEA",
  agricultural_upstream_pct   = "#F0EAD6"
)

# Scale percent-type covariates to percentage (0-100) for display.
pct_covariates <- c("high_cascades_upstream_pct", "developed_upstream_pct", "agricultural_upstream_pct")

GRID_SIZE <- 20
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# iml Predictor wrapper for ranger. ranger's predict() returns a list with $predictions, not a plain vector.
# iml needs a function returning a plain numeric vector/data.frame.
predict_function <- function(model, newdata) {
  predict(model, data = newdata)$predictions
}

predictor <- Predictor$new(
  model = rf_final_model,
  data = rfData[, covariates_40],
  y = rfData$thermal_sensitivity,
  predict.function = predict_function
)

# ============================================================
# Compute ALE for each of the 6 significant covariates.
# ============================================================
cat("[CHECK 1] Computing ALE effects for 6 significant covariates...\n")

ale_results_list <- list()

for (cov in rf_significant_covariates) {
  cat("  ...", cov, "\n")
  ale_effect <- FeatureEffect$new(predictor, feature = cov, method = "ale", grid.size = GRID_SIZE)

  df <- ale_effect$results
  df <- df %>% rename(grid_value_raw = all_of(cov), ale_value = .value) %>%
    select(grid_value_raw, ale_value)
  df$covariate <- cov

  df$grid_value_display <- if (cov %in% pct_covariates) df$grid_value_raw * 100 else df$grid_value_raw

  ale_results_list[[cov]] <- df
}

ale_results_significant <- bind_rows(ale_results_list)
cat("[CHECK 2] ALE results (6 significant covariates):\n")
print(ale_results_significant)

dir.create("results_2021/rf/ale", recursive = TRUE, showWarnings = FALSE)
write_csv(ale_results_significant, "results_2021/rf/ale/ale_results_significant.csv")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting: ALE curve plus a data-density rug (full, untrimmed observed values, scaled to match the display axis).
cat("[CHECK 3] Generating individual panels...\n")

panel_plots <- list()

for (cov in rf_significant_covariates) {
  df <- ale_results_list[[cov]]
  bg_color <- covariate_colors[[cov]]
  full_name <- covariate_full_names[[cov]]

  rug_values <- if (cov %in% pct_covariates) rfData[[cov]] * 100 else rfData[[cov]]
  rug_df <- data.frame(display_value = rug_values)

  tag_pos <- if (cov == "agricultural_upstream_pct") c(0.05, 0.94) else c(0.05, 0.97)
  
  p <- ggplot(df, aes(x = grid_value_display, y = ale_value)) +
    geom_line(color = "black", linewidth = 1) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dashed") +
    geom_rug(data = rug_df, aes(x = display_value), inherit.aes = FALSE, sides = "b", color = "grey50", alpha = 0.5) +
    labs(
      x = full_name,
      y = "ALE (centered effect)"
    ) +
    theme_minimal(base_family = "Verdana", base_size = 12) +
    theme(
      panel.background = element_rect(fill = bg_color, color = NA),
      plot.background = element_rect(fill = bg_color, color = "white", linewidth = 6),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 10, b = 8, l = 8),
      axis.title = element_text(face = "bold"),
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = tag_pos
    )

  panel_letter <- letters[which(rf_significant_covariates == cov)]

  p <- p + labs(tag = paste0("(", panel_letter, ")"))

  panel_plots[[cov]] <- p

  ggsave(paste0("results_2021/rf/ale/ale_", cov, ".png"),
         plot = p, width = 6, height = 4.5, dpi = 300, bg = bg_color)
}

cat(sprintf("Done. Saved %d individual panel PNGs (300 DPI, backup only).\n", length(rf_significant_covariates)))
cat("Saved ale_results_significant6.csv.\n")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Figure 3: ALE plots for 6 significant covariates (a-f), combined into a single figure (3 rows x 2 columns).
# (a) base_flow_index, (b) high_cascades_upstream_pct, (c) veg_height_upstream_m,
# (d) elevation_m, (e) developed_upstream_pct, (f) agricultural_upstream_pct.
cat("[CHECK 4] Assembling combined Figure 3 (3x2, panels a-f)...\n")

figure3_combined <- wrap_plots(panel_plots, ncol = 2, nrow = 3) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(
  "results_2021/rf/ale/Figure_3_ALE_combined.png",
  plot = figure3_combined,
  width = 12, height = 13.5, dpi = 600, bg = "white"
)
cat("Saved Figure_3_ALE_combined.png (600 DPI, 12 x 13.5 in).\n")