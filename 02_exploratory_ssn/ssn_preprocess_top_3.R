# ssn_preprocess_top_3.R
# Joins the union of covariates needed across the 3 ArcGIS Exploratory Regression competing models (Delta AICc <= 2). 

library(readr)
library(dplyr)
library(sf)
library(SSN2)

rm(list = ls())

CRB_SSN <- ssn_import("data/ClackHist.ssn", predpts = "preds", overwrite = TRUE)
data_2021 <- read_csv("data/clackamas_thermal_sensitivity_covariates_2021.csv")
n_raw_ssn <- nrow(CRB_SSN$obs)
n_data <- nrow(data_2021)
cat(sprintf("[CHECK 0] Raw SSN obs: %d rows | data_2021: %d rows\n", n_raw_ssn, n_data))
stopifnot("data_2021 does not have 72 sites as expected" = n_data == 72)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
gridcodes_missing_from_ssn <- setdiff(data_2021$gridcode, CRB_SSN$obs$GRIDCODE)
cat(sprintf("[CHECK 1] GRIDCODE(s) in data_2021 with no match in SSN object: %s\n", ifelse(length(gridcodes_missing_from_ssn) == 0, "none", paste(gridcodes_missing_from_ssn, collapse = ", "))))
excluded_sites <- data_2021 %>% filter(gridcode %in% gridcodes_missing_from_ssn) %>% select(index, site, stream_name, gridcode)
cat("Site(s) that will be EXCLUDED from the SSN model (no stream network match) -- empty table below means none:\n")
print(excluded_sites)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
ssn_pid_counts <- CRB_SSN$obs %>% st_drop_geometry() %>% count(GRIDCODE, name = "n_ssn_pids")
data_row_counts <- data_2021 %>% count(gridcode, name = "n_data_rows") %>% rename(GRIDCODE = gridcode)
gridcode_mismatch <- ssn_pid_counts %>% inner_join(data_row_counts, by = "GRIDCODE") %>% filter(n_ssn_pids > n_data_rows)
cat("[CHECK 2] GRIDCODE(s) with more SSN pids than data_2021 rows -- empty table below means none:\n")
print(gridcode_mismatch)

CRB_SSN$obs <- CRB_SSN$obs %>% filter(!(GRIDCODE == 479590 & pid == 503))
cat("[CHECK 2] Resolved: excluded pid 503 for GRIDCODE 479590 (395.2m from site 26's recorded location; matched pid 602 was 30.9m away).\n")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Join the union of covariates needed across all 3 ArcGIS Exploratory Regression models:
# 1) MLR Model 1: elevation_m, high_cascades_upstream_pct, shrub_upstream_pct, upstream_area_km2, burned_buffer_pct
# 2) MLR Model 2: elevation_m, base_flow_index, developed_upstream_pct, shrub_upstream_pct, upstream_area_km2
# 3) MLR Model 3: channel_slope, base_flow_index, high_cascades_upstream_pct, shrub_upstream_pct, upstream_area_km2
model_covariates_union <- c("channel_slope", "elevation_m", "high_cascades_upstream_pct", "shrub_upstream_pct", "upstream_area_km2", "burned_buffer_pct", "base_flow_index", "developed_upstream_pct")

ssn_covariates <- data_2021 %>% select(index, gridcode, thermal_sensitivity, all_of(model_covariates_union))
CRB_SSN$obs <- CRB_SSN$obs %>% left_join(ssn_covariates, by = join_by(GRIDCODE == gridcode), relationship = "many-to-many")
n_after_join <- nrow(CRB_SSN$obs)
cat(sprintf("[CHECK 3] Rows after join: %d\n", n_after_join))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# NA filter across the 8-covariate union to ensure stability across models
na_counts_by_covariate <- sapply(model_covariates_union, function(v) sum(is.na(CRB_SSN$obs[[v]])))
cat("[CHECK 4] NA count per covariate (union filter applied across all 8):\n")
print(na_counts_by_covariate)

complete_case_mask <- complete.cases(CRB_SSN$obs %>% st_drop_geometry() %>% select(thermal_sensitivity, all_of(model_covariates_union)))
CRB_SSN$obs <- CRB_SSN$obs %>% filter(complete_case_mask)
n_after_na_filter <- nrow(CRB_SSN$obs)
cat(sprintf("[CHECK 4] Rows after NA filter: %d (dropped %d)\n", n_after_na_filter, n_after_join - n_after_na_filter))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
n_before_dedup <- nrow(CRB_SSN$obs)
CRB_SSN$obs <- CRB_SSN$obs %>% filter(!(pid == 53 & index == 24)) %>% filter(!(pid == 56 & index == 25)) %>% filter(!(pid == 167 & index == 66)) %>% filter(!(pid == 383 & index == 65)) %>% filter(!(pid == 167 & index == 69)) %>% filter(!(pid == 383 & index == 68)) %>% filter(!(pid == 242 & index == 20)) %>% filter(!(pid == 332 & index == 15)) %>% filter(!(pid == 644 & index == 45)) %>% filter(!(pid == 294 & index == 44))
n_after_dedup <- nrow(CRB_SSN$obs)
cat(sprintf("[CHECK 5] Rows removed by duplicate filter: %d (expected 10)\n", n_before_dedup - n_after_dedup))
stopifnot("Duplicate filter removed an unexpected number of rows -- verify pid/index pairs still apply to this join" = (n_before_dedup - n_after_dedup) == 10)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
max_pid <- max(CRB_SSN$obs$pid, na.rm = TRUE)
max_locID <- max(CRB_SSN$obs$locID, na.rm = TRUE)
CRB_SSN$obs <- CRB_SSN$obs %>% mutate(pid = case_when(pid == 77 & index == 42 ~ max_pid + 1, pid == 167 & index == 68 ~ max_pid + 2, pid == 383 & index == 69 ~ max_pid + 3, pid == 400 & index == 50 ~ max_pid + 4, TRUE ~ pid), locID = case_when(locID == 158 & index == 42 ~ max_locID + 1, locID == 141 & index == 68 ~ max_locID + 2, locID == 74 & index == 69 ~ max_locID + 3, locID == 68 & index == 50 ~ max_locID + 4, TRUE ~ locID))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
CRB_SSN$obs <- CRB_SSN$obs %>% mutate(netgeom = case_when(pid == max_pid + 1 ~ gsub(" 77 ", paste0(" ", max_pid + 1, " "), netgeom), pid == max_pid + 2 ~ gsub(" 167 ", paste0(" ", max_pid + 2, " "), netgeom), pid == max_pid + 3 ~ gsub(" 383 ", paste0(" ", max_pid + 3, " "), netgeom), pid == max_pid + 4 ~ gsub(" 400 ", paste0(" ", max_pid + 4, " "), netgeom), TRUE ~ netgeom))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
final_pid_dupes <- CRB_SSN$obs %>% st_drop_geometry() %>% count(pid) %>% filter(n > 1)
final_locID_dupes <- CRB_SSN$obs %>% st_drop_geometry() %>% count(locID) %>% filter(n > 1)
final_index_dupes <- CRB_SSN$obs %>% st_drop_geometry() %>% count(index) %>% filter(n > 1)
cat(sprintf("[CHECK 8] Final pid duplicates: %d | locID duplicates: %d | index duplicates: %d\n", nrow(final_pid_dupes), nrow(final_locID_dupes), nrow(final_index_dupes)))
stopifnot("Duplicate pid remaining after cleanup" = nrow(final_pid_dupes) == 0, "Duplicate locID remaining after cleanup" = nrow(final_locID_dupes) == 0, "Duplicate index remaining after cleanup" = nrow(final_index_dupes) == 0)
cat("[CHECK 8] PASSED: object is valid.\n")
cat(sprintf("[CHECK 8] NOTE: final N = %d. Expected 72 (all data_2021 sites; site 38's GRIDCODE\n", nrow(CRB_SSN$obs)))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
row.names(CRB_SSN$obs) <- seq_len(nrow(CRB_SSN$obs))
ssn_create_distmat(ssn.object = CRB_SSN, overwrite = TRUE)
saveRDS(CRB_SSN, "results_2021/exploratory_ssn/CRB_SSN_cleaned_top_3.RDS")
cat(sprintf("\nSaved CRB_SSN_cleaned_top_3.RDS with %d validated sites.\n", nrow(CRB_SSN$obs)))