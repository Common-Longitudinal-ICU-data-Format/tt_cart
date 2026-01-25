# ==============================================================================
# 06_export.R
# Final assembly and export of poolable outputs for CAR-T project
# Outputs: summary tables for pooling (n's, sums, sum_sq by group)
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("save_output")) {
  stop("Run 00_setup.R through 05_easix.R first.", call. = FALSE)
}

message("Starting 06_export.R")


# save cohort data frame (intermediate, not for sharing) =======================

save_output(df, "cohort_analytic", shareable = FALSE)


# trajectory distribution ======================================================

traj_dist = 
  fgroup_by(df, trajectory) |>
  fsummarize(n = fnobs(patient_id)) |>
  fungroup()

save_output(traj_dist, "trajectory_distribution", shareable = TRUE)

# summary statistics for pooling ===============================================

## continuous variables: output n, sum, sum_sq for each group ------------------

continuous_vars = c(
  "age", "vw", "los_hospital_d", "time_to_death",
  "h_until_toci", "h_until_anakinra", "h_until_steroids",
  "apap_day_0", "apap_day_1", "apap_day_2",
  "easix", "log2_easix", "wbc_infusion", "lymph_abs_infusion",
  "cart_hour_precise"
)

continuous_summary = 
  lapply(continuous_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    
    x = df[[v]]
    
    fselect(df, trajectory) |>
      fmutate(
        x    = x,
        x_sq = x^2
      ) |>
      fgroup_by(trajectory) |>
      fsummarize(
        n_obs    = fnobs(x),
        sum_x    = fsum(x, na.rm = TRUE),
        sum_x_sq = fsum(x_sq, na.rm = TRUE)
      ) |>
      fungroup() |>
      fmutate(variable = v)
  }) |>
  rbindlist()

save_output(continuous_summary, "continuous_vars_by_trajectory", shareable = TRUE)

## binary variables: output n and sum (count of 1s) for each group -------------

binary_vars = c(
  "dc_dead_01", "dc_hospice_01", 
  "dead_30_01", "dead_90_01", "dead_180_01", "dead_365_01",
  "icu_30_01", "resp_30_01", "imv_30_01", "va_30_01",
  "bsi_30_01", "abx_30_01",
  "steroid_any_72h", "hct_01", "miss_vw"
)

binary_summary = 
  lapply(binary_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    
    x = df[[v]]
    
    fselect(df, trajectory) |>
      fmutate(x = x) |>
      fgroup_by(trajectory) |>
      fsummarize(
        n_obs    = fnobs(x),
        n_events = fsum(x, na.rm = TRUE)
      ) |>
      fungroup() |>
      fmutate(variable = v)
  }) |>
  rbindlist()

save_output(binary_summary, "binary_vars_by_trajectory", shareable = TRUE)

## categorical variables: output n per category per group ----------------------

categorical_vars = c(
  "sex_category", "race_category", "ethnicity_category",
  "cancer_cat", "cart_name", "cart_group",
  "easix_f_cat", "easix_fc_cat", "fit_flag"
)

categorical_summary = 
  lapply(categorical_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    
    tmp = tidytable(
      trajectory = df$trajectory,
      cat_val    = df[[v]]
    )
    
    fgroup_by(tmp, trajectory, cat_val) |>
      fsummarize(n = fnobs(trajectory)) |>
      fungroup() |>
      fmutate(variable = v) |>
      frename(category = cat_val)
  }) |>
  rbindlist()

save_output(categorical_summary, "categorical_vars_by_trajectory", shareable = TRUE)

# overall summary (not stratified) =============================================

## for comparison and total N verification

overall_continuous = 
  lapply(continuous_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    
    tidytable(
      variable = v,
      n_obs    = fnobs(df[[v]]),
      sum_x    = fsum(df[[v]], na.rm = TRUE),
      sum_x_sq = fsum(df[[v]]^2, na.rm = TRUE)
    )
  }) |>
  rbindlist()

overall_binary = 
  lapply(binary_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    
    tidytable(
      variable = v,
      n_obs    = fnobs(df[[v]]),
      n_events = fsum(df[[v]], na.rm = TRUE)
    )
  }) |>
  rbindlist()

save_output(overall_continuous, "continuous_vars_overall", shareable = TRUE)
save_output(overall_binary,     "binary_vars_overall",     shareable = TRUE)

# circadian variables (carried forward for future analysis) ====================

## summary by hour of infusion -------------------------------------------------

circadian_by_hour = 
  fgroup_by(df, cart_hour_floor) |>
  fsummarize(
    n               = fnobs(patient_id),
    dead_30_n       = fsum(dead_30_01, na.rm = TRUE),
    icu_30_n        = fsum(icu_30_01, na.rm = TRUE),
    toci_given_n    = fsum(toci_n_doses_30d > 0, na.rm = TRUE),
    steroid_72h_n   = fsum(steroid_any_72h, na.rm = TRUE)
  ) |>
  fungroup()

save_output(circadian_by_hour, "outcomes_by_infusion_hour", shareable = TRUE)


# summary of what was exported =================================================

message("\n=== Export Summary ===")
message("Shareable files saved to: ", file.path(project_location, "upload_to_box"))

export_files = list.files(
  file.path(project_location, "upload_to_box"), 
  pattern = paste0("_", site_lowercase, "_", today),
  full.names = FALSE
)

message("Files created:")
invisible(lapply(export_files, \(f) message("  - ", f)))

message("\nTotal cohort N: ", fnrow(df))
message("Trajectories assigned: ", fnrow(fsubset(df, !is.na(trajectory))))

message("\n06_export.R complete.")