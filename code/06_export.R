# ==============================================================================
# 06_export.R
# Final assembly and export of poolable outputs for CAR-T project
# All outputs saved to project_location/upload_to_box/{trajectory|circadian}/
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("save_output")) {
  stop("Run 00_setup.R through 05_abxs_cultures.R first.", call. = FALSE)
}

message("Starting 06_export.R")


# derived variables ============================================================

data_end_date = ymd("2024-12-31")

df = df |>
  fmutate(
    # everyone observed until death or 12/31/24
    followup_end = fifelse(!is.na(death_dttm), death_dttm, as.POSIXct(data_end_date)),
    followup_d   = as.numeric(difftime(followup_end, cart_t0, units = "days")),
    died         = as.integer(!is.na(time_to_death)),
    
    # calendar year
    infusion_year = year(cart_t0),
    
    # CAR-T target (CD19 vs BCMA)
    cart_target = case_when(
      cart_group %in% c(3L, 4L)         ~ "BCMA",
      cart_group %in% c(1L, 2L, 5L, 6L) ~ "CD19",
      TRUE                              ~ NA_character_
    ),
    
    # race/ethnicity collapsed (NHW vs other)
    race_eth = case_when(
      race_category == "white" & ethnicity_category == "non-hispanic" ~ "nhw",
      TRUE ~ "other"
    ),
    
    # treatment indicators
    toci_any_30d     = as.integer(toci_n_doses_30d > 0),
    anakinra_any_30d = as.integer(anakinra_n_doses_30d > 0)
  )


# save analytic cohort (intermediate) ==========================================

save_output(df, "cohort_analytic", shareable = FALSE)


# ==============================================================================
# TRAJECTORY PROJECT EXPORTS
# ==============================================================================

message("\n--- Trajectory Project Exports ---")

# table 1: patient characteristics =============================================

## continuous (including baseline labs) ----------------------------------------

t1_cont_vars = c(
  "age", "vw", "easix", "log2_easix", "wbc_infusion", "lymph_abs_infusion",
  "baseline_ldh", "baseline_crp", "baseline_ferritin", 
  "baseline_wbc", "baseline_lymphocytes_absolute"
)

t1_continuous = 
  lapply(t1_cont_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    x = df[[v]]
    fselect(df, trajectory) |>
      fmutate(x = x, x_sq = x^2) |>
      fgroup_by(trajectory) |>
      fsummarize(
        n      = fnobs(x), 
        n_miss = fsum(is.na(x)), 
        sum_x  = fsum(x, na.rm = TRUE), 
        sum_x_sq = fsum(x_sq, na.rm = TRUE)
      ) |>
      fungroup() |>
      fmutate(variable = v)
  }) |>
  rbindlist()

save_output(t1_continuous, "table1_continuous", shareable = TRUE, subfolder = "trajectory")

## categorical -----------------------------------------------------------------

t1_cat_vars = c(
  "sex_category", "race_eth", "cancer_cat", "cart_target", 
  "easix_fc_cat", "hct_01", "infusion_year"
)

t1_categorical = 
  lapply(t1_cat_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    tmp = tidytable(trajectory = df$trajectory, cat_val = as.character(df[[v]]))
    fgroup_by(tmp, trajectory, cat_val) |>
      fsummarize(n = fnobs(trajectory)) |>
      fungroup() |>
      fmutate(variable = v) |>
      frename(category = cat_val)
  }) |>
  rbindlist()

save_output(t1_categorical, "table1_categorical", shareable = TRUE, subfolder = "trajectory")


# table 2: outcomes ============================================================

## survival data at q30d intervals for 1 year ----------------------------------

time_points = seq(0, 365, by = 30)

surv_data = 
  lapply(time_points, \(t) {
    prev_t = ifelse(t == 0, 0, t - 30)
    
    fgroup_by(df, trajectory) |>
      fsummarize(
        n_at_risk  = fsum(followup_d >= prev_t, na.rm = TRUE),
        n_events   = fsum(died == 1 & time_to_death > prev_t & time_to_death <= t, na.rm = TRUE),
        n_censored = fsum(died == 0 & followup_d > prev_t & followup_d <= t, na.rm = TRUE)
      ) |>
      fungroup() |>
      fmutate(time_point = t)
  }) |>
  rbindlist()

save_output(surv_data, "survival_data_q30d", shareable = TRUE, subfolder = "trajectory")

## acute outcomes --------------------------------------------------------------

acute_vars = c(
  "icu_30_01", "imv_30_01", "va_30_01", "bsi_30_01", "abx_30_01",
  "toci_any_30d", "anakinra_any_30d", "steroid_any_72h",
  "severe_icaht_01", "hyperferritin_5k_01", "hyperferritin_10k_01"
)

acute_summary = 
  lapply(acute_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    x = df[[v]]
    fselect(df, trajectory) |>
      fmutate(x = x) |>
      fgroup_by(trajectory) |>
      fsummarize(n = fnobs(x), n_events = fsum(x, na.rm = TRUE)) |>
      fungroup() |>
      fmutate(variable = v)
  }) |>
  rbindlist()

save_output(acute_summary, "table2_acute", shareable = TRUE, subfolder = "trajectory")

## hematotoxicity detail -------------------------------------------------------

hematox_vars = c("anc_nadir", "plt_nadir", "hgb_nadir", "ferritin_max_30d")

hematox_summary = 
  lapply(hematox_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    x = df[[v]]
    fselect(df, trajectory) |>
      fmutate(x = x, x_sq = x^2) |>
      fgroup_by(trajectory) |>
      fsummarize(
        n        = fnobs(x), 
        sum_x    = fsum(x, na.rm = TRUE), 
        sum_x_sq = fsum(x_sq, na.rm = TRUE)
      ) |>
      fungroup() |>
      fmutate(variable = v)
  }) |>
  rbindlist()

save_output(hematox_summary, "table2_hematox", shareable = TRUE, subfolder = "trajectory")

## ICAHT grade distribution ----------------------------------------------------

icaht_dist = 
  tidytable(trajectory = df$trajectory, icaht_grade = as.character(df$icaht_grade)) |>
  fgroup_by(trajectory, icaht_grade) |>
  fsummarize(n = fnobs(trajectory)) |>
  fungroup()

save_output(icaht_dist, "table2_icaht_grades", shareable = TRUE, subfolder = "trajectory")

## time-to-event / LOS ---------------------------------------------------------

time_vars = c("los_hospital_d", "h_until_toci", "h_until_anakinra", "h_until_steroids")

time_summary = 
  lapply(time_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    x = df[[v]]
    fselect(df, trajectory) |>
      fmutate(x = x, x_sq = x^2) |>
      fgroup_by(trajectory) |>
      fsummarize(
        n        = fnobs(x), 
        sum_x    = fsum(x, na.rm = TRUE), 
        sum_x_sq = fsum(x_sq, na.rm = TRUE)
      ) |>
      fungroup() |>
      fmutate(variable = v)
  }) |>
  rbindlist()

save_output(time_summary, "table2_time", shareable = TRUE, subfolder = "trajectory")


# RMST pseudo-observations =====================================================

if (requireNamespace("pseudo", quietly = TRUE)) {
  
  message("Calculating RMST pseudo-observations at 365 days...")
  
  surv_time  = pmin(df$followup_d, 365, na.rm = TRUE)
  surv_event = fifelse(df$died == 1 & df$time_to_death <= 365, 1L, 0L)
  
  valid_idx = !is.na(surv_time) & !is.na(surv_event)
  pseudo_vals = rep(NA_real_, nrow(df))
  
  if (sum(valid_idx) > 10) {
    pseudo_vals[valid_idx] = pseudo::pseudomean(
      surv_time[valid_idx], 
      surv_event[valid_idx], 
      tmax = 365
    )
  }
  
  df$pseudo_rmst_365 = pseudo_vals
  
  # patient-level data for local analysis only
  pseudo_export = fselect(df, 
    patient_id, trajectory, fit_flag,
    age, sex_category, race_eth, cancer_cat, cart_target, vw, easix_fc_cat, infusion_year,
    cart_hour_precise,
    pseudo_rmst_365
  )
  
  save_output(pseudo_export, "rmst_pseudo_365", shareable = FALSE)
  
  message("Pseudo-observations saved to intermediate_output/ for use in 07/08 analysis scripts.")
  
} else {
  message("Package 'pseudo' not installed - skipping RMST.")
  message("Install with: install.packages('pseudo')")
}


# misclassification sensitivity ================================================

misclass = 
  fgroup_by(df, trajectory, second_best, fit_flag) |>
  fsummarize(
    n               = fnobs(patient_id),
    dead_90_n       = fsum(dead_90_01, na.rm = TRUE),
    icu_30_n        = fsum(icu_30_01, na.rm = TRUE),
    imv_30_n        = fsum(imv_30_01, na.rm = TRUE),
    toci_n          = fsum(toci_any_30d, na.rm = TRUE),
    anakinra_n      = fsum(anakinra_any_30d, na.rm = TRUE),
    hyperferr_10k_n = fsum(hyperferritin_10k_01, na.rm = TRUE)
  ) |>
  fungroup()

save_output(misclass, "sensitivity_misclassification", shareable = TRUE, subfolder = "trajectory")

margin_dist = 
  fgroup_by(df, trajectory) |>
  fsummarize(
    n             = fnobs(margin_rel),
    sum_margin    = fsum(margin_rel, na.rm = TRUE),
    sum_margin_sq = fsum(margin_rel^2, na.rm = TRUE),
    n_borderline  = fsum(fit_flag == "borderline", na.rm = TRUE)
  ) |>
  fungroup()

save_output(margin_dist, "sensitivity_margin", shareable = TRUE, subfolder = "trajectory")


# model data for local adjusted analysis =======================================

model_vars = c(
  "patient_id", "trajectory", "fit_flag",
  "dead_30_01", "dead_90_01", "dead_180_01", "dead_365_01",
  "icu_30_01", "imv_30_01", "va_30_01", "bsi_30_01",
  "toci_any_30d", "anakinra_any_30d",
  "severe_icaht_01", "hyperferritin_5k_01", "hyperferritin_10k_01",
  "followup_d", "died", "time_to_death",
  "age", "sex_category", "race_eth", "cart_target", "vw", "easix_fc_cat", "infusion_year",
  "cart_hour_precise", "cart_hour_floor"
)

model_vars_present = intersect(model_vars, names(df))
model_data = fselect(df, model_vars_present)
save_output(model_data, "model_data", shareable = FALSE)


# site summary =================================================================

site_summary = tidytable(
  site            = site_lowercase,
  total_n         = fnrow(df),
  n_with_traj     = fsum(!is.na(df$trajectory)),
  n_ht            = fsum(df$trajectory == "HT", na.rm = TRUE),
  n_nt            = fsum(df$trajectory == "NT", na.rm = TRUE),
  n_hfr           = fsum(df$trajectory == "HFR", na.rm = TRUE),
  n_hsr           = fsum(df$trajectory == "HSR", na.rm = TRUE),
  n_borderline    = fsum(df$fit_flag == "borderline", na.rm = TRUE),
  n_cd19          = fsum(df$cart_target == "CD19", na.rm = TRUE),
  n_bcma          = fsum(df$cart_target == "BCMA", na.rm = TRUE),
  n_severe_icaht  = fsum(df$severe_icaht_01, na.rm = TRUE),
  n_hyperferr_10k = fsum(df$hyperferritin_10k_01, na.rm = TRUE),
  date_min        = as.character(min(df$admin_dttm, na.rm = TRUE)),
  date_max        = as.character(max(df$admin_dttm, na.rm = TRUE))
)

save_output(site_summary, "site_summary", shareable = TRUE, subfolder = "trajectory")


# ==============================================================================
# CIRCADIAN PROJECT EXPORTS
# ==============================================================================

message("\n--- Circadian Project Exports ---")

## outcomes by hour ------------------------------------------------------------

circ_by_hour = 
  fgroup_by(df, cart_hour_floor) |>
  fsummarize(
    n               = fnobs(patient_id),
    dead_30_n       = fsum(dead_30_01, na.rm = TRUE),
    dead_90_n       = fsum(dead_90_01, na.rm = TRUE),
    icu_30_n        = fsum(icu_30_01, na.rm = TRUE),
    imv_30_n        = fsum(imv_30_01, na.rm = TRUE),
    va_30_n         = fsum(va_30_01, na.rm = TRUE),
    toci_n          = fsum(toci_any_30d, na.rm = TRUE),
    anakinra_n      = fsum(anakinra_any_30d, na.rm = TRUE),
    severe_icaht_n  = fsum(severe_icaht_01, na.rm = TRUE),
    hyperferr_10k_n = fsum(hyperferritin_10k_01, na.rm = TRUE)
  ) |>
  fungroup()

save_output(circ_by_hour, "outcomes_by_hour", shareable = TRUE, subfolder = "circadian")

## outcomes by hour and trajectory (for interaction) ---------------------------

circ_by_hour_traj = 
  fgroup_by(df, cart_hour_floor, trajectory) |>
  fsummarize(
    n         = fnobs(patient_id),
    dead_90_n = fsum(dead_90_01, na.rm = TRUE),
    icu_30_n  = fsum(icu_30_01, na.rm = TRUE),
    toci_n    = fsum(toci_any_30d, na.rm = TRUE)
  ) |>
  fungroup()

save_output(circ_by_hour_traj, "outcomes_by_hour_trajectory", shareable = TRUE, subfolder = "circadian")

## covariates by hour ----------------------------------------------------------

circ_covars = 
  fgroup_by(df, cart_hour_floor) |>
  fsummarize(
    n          = fnobs(patient_id),
    age_sum    = fsum(age, na.rm = TRUE),
    age_sum_sq = fsum(age^2, na.rm = TRUE),
    n_male     = fsum(sex_category == "male", na.rm = TRUE),
    n_nhw      = fsum(race_eth == "nhw", na.rm = TRUE),
    n_cd19     = fsum(cart_target == "CD19", na.rm = TRUE),
    vw_sum     = fsum(vw, na.rm = TRUE),
    vw_sum_sq  = fsum(vw^2, na.rm = TRUE)
  ) |>
  fungroup()

save_output(circ_covars, "covariates_by_hour", shareable = TRUE, subfolder = "circadian")


# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n=== Export Complete ===")
message("Trajectory outputs: ", file.path(project_location, "upload_to_box", "trajectory"))
message("Circadian outputs:  ", file.path(project_location, "upload_to_box", "circadian"))
message("\nTotal N: ", fnrow(df))
message("With trajectory: ", fsum(!is.na(df$trajectory)))

message("\n06_export.R complete.")
