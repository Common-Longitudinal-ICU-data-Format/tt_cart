# ==============================================================================
# 04_outcomes.R
# Outcomes and treatment variables for CAR-T project
# Requires: df, hid_jid_crosswalk, anchor_df from prior scripts
# Outputs: outcomes added to df, blood culture table by trajectory, lab trends
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("hid_jid_crosswalk") || !exists("anchor_df") || !exists("save_output")) {
  stop("Run 00_setup.R through 03_covariates.R first.", call. = FALSE)
}

message("Starting 04_outcomes.R")


# death outcomes ===============================================================

df$death_dttm = case_when(
  !is.na(df$death_dttm)                       ~ df$death_dttm, 
  tolower(df$discharge_category) == "expired" ~ df$discharge_dttm,
  TRUE                                        ~ as.POSIXct(NA)
)

df$discharge_category_clean = case_when(
  !is.na(df$death_dttm) & df$death_dttm <= df$discharge_dttm ~ "expired",
  TRUE                                                       ~ tolower(df$discharge_category)
)

df$dc_dead_01    = fifelse(df$discharge_category_clean == "expired", 1L, 0L)
df$dc_hospice_01 = fifelse(df$discharge_category_clean == "hospice", 1L, 0L)
df$time_to_death = as.numeric(difftime(df$death_dttm, df$cart_t0, units = "days"))
df$dead_30_01    = as.integer(!is.na(df$time_to_death) & df$time_to_death >= 0 & df$time_to_death <= 30)
df$dead_90_01    = as.integer(!is.na(df$time_to_death) & df$time_to_death >= 0 & df$time_to_death <= 90)
df$dead_180_01   = as.integer(!is.na(df$time_to_death) & df$time_to_death >= 0 & df$time_to_death <= 180)
df$dead_365_01   = as.integer(!is.na(df$time_to_death) & df$time_to_death >= 0 & df$time_to_death <= 365)


# index hospital length of stay ===============================================================

df$los_hospital_d = as.numeric(difftime(df$discharge_dttm, df$admission_dttm, units = "days"))


# create crosswalk for 30 day hospitalizations ================================================
all_hid_crosswalk = 
  data_list$hospitalization |>
  dplyr::filter(patient_id %in% cart_pids) |>
  dplyr::filter(admission_dttm >= start_date) |>
  dplyr::select(patient_id, hospitalization_id) |>
  dplyr::collect()

# ICU outcomes =================================================================

## ICU admission within 30 days ------------------------------------------------

icu_flag = 
  data_list$adt |>
  dplyr::filter(tolower(location_category) == "icu") |>
  dplyr::filter(in_dttm >= start_date) |>
  dplyr::distinct(hospitalization_id, in_dttm) |>
  dplyr::collect() 

icu_flag = 
  join(icu_flag, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, icu_dttm = in_dttm)

icu_flag = 
  join(icu_flag, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>  
  fsubset(icu_dttm >= cart_t0 & icu_dttm <= cart_t0 + ddays(30)) |>
  pull(patient_id) |>
  funique()

## vasopressors within 30 days -------------------------------------------------

vaso_flag = 
  data_list$medication_admin_continuous |>
  dplyr::filter(tolower(med_group) == "vasoactives") |> 
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::distinct(hospitalization_id, admin_dttm) |>
  dplyr::select(hospitalization_id, admin_dttm) |>
  dplyr::collect() 

vaso_flag = 
  join(vaso_flag, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, va_dttm = admin_dttm)

vaso_flag = 
  join(vaso_flag, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fsubset(va_dttm >= cart_t0 & va_dttm <= cart_t0 + ddays(30)) |>
  pull(patient_id) |>
  funique()

## respiratory support within 30 days ------------------------------------------

resp_flag = 
  data_list$respiratory_support |>
  dplyr::filter(tolower(device_category) %in% c("imv", "nippv", "high flow nc")) |>
  dplyr::filter(recorded_dttm >= start_date) |>
  dplyr::distinct(hospitalization_id, recorded_dttm, device_category) |>
  dplyr::select(hospitalization_id, recorded_dttm, device_category) |>
  dplyr::collect()

resp_flag = 
  join(resp_flag, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, recorded_dttm, device_category)

resp_flag = 
  join(resp_flag, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fsubset(recorded_dttm >= cart_t0 & recorded_dttm <= cart_t0 + ddays(30))

imv_flag =
  fsubset(resp_flag, tolower(device_category) == "imv") |>
  pull(patient_id) |>
  funique()

resp_flag =
  pull(resp_flag, patient_id) |>
  funique()

## add ICU flags to df ---------------------------------------------------------

df$resp_30_01 = as.integer(df$patient_id %in% resp_flag)
df$imv_30_01  = as.integer(df$patient_id %in% imv_flag)
df$va_30_01   = as.integer(df$patient_id %in% vaso_flag)

df$icu_30_01 = case_when(
  df$patient_id %in% icu_flag ~ 1L,
  df$resp_30_01 == 1L         ~ 1L,
  df$va_30_01   == 1L         ~ 1L,
  TRUE                        ~ 0L
)

rm(vaso_flag, resp_flag, imv_flag, icu_flag)
gc()


# lab outcomes =================================================================

## define labs of interest -----------------------------------------------------

lab_categories = c(
  "ldh", "crp", "ferritin", "wbc", "lymphocytes_absolute",
  "hemoglobin", "platelet_count", "neutrophils_absolute", "creatinine", 
  "alt", "ast", "bilirubin_total"
)

## load all relevant labs ------------------------------------------------------

labs_raw = 
  data_list$labs |>
  dplyr::filter(tolower(lab_category) %in% lab_categories) |>
  dplyr::filter(lab_result_dttm >= start_date) |>
  dplyr::filter(!is.na(lab_value_numeric)) |>
  dplyr::select(hospitalization_id, lab_result_dttm, lab_category, lab_value_numeric) |>
  dplyr::collect() |>
  distinct()

labs_raw = 
  join(labs_raw, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, lab_result_dttm, lab_category, lab_value_numeric) |>
  distinct()

labs_raw = 
  join(labs_raw, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fmutate(
    d_since_cart = as.integer(floor(as.numeric(difftime(lab_result_dttm, cart_t0, units = "days"))))
  ) |>
  frename(lab_dttm = lab_result_dttm, lab_value = lab_value_numeric)

##Baseline labs (day -1 to day 0, closest to infusion) ------------------------

baseline_labs = 
  fsubset(labs_raw, lab_dttm >= cart_t0 - dhours(24) & lab_dttm <= cart_t0) |>
  fmutate(hours_before_cart = as.numeric(difftime(cart_t0, lab_dttm, units = "hours"))) |>
  roworder(patient_id, lab_category, hours_before_cart) |>
  fgroup_by(patient_id, lab_category) |>
  ffirst() |>
  fselect(patient_id, lab_category, lab_value) |>
  pivot_wider(
    names_from   = lab_category,
    values_from  = lab_value,
    names_prefix = "baseline_"
  )

df = join(df, baseline_labs, how = "left", multiple = FALSE)

rm(baseline_labs)
gc()

##Nadir labs--------------------------------------------------------------------
# Filter to hematology labs for days 0-30
hematox_labs = 
  fsubset(labs_raw, d_since_cart >= 0 & d_since_cart <= 30) |>
  fsubset(lab_category %in% c("neutrophils_absolute", "hemoglobin", "lymphocytes_absolute", "platelet_count"))

# Calculate daily nadirs (minimum value per patient per day per lab)
daily_nadirs = 
  fgroup_by(hematox_labs, patient_id, lab_category, d_since_cart) |>
  fsummarize(nadir_value = fmin(lab_value, na.rm = TRUE)) |>
  fungroup()

# Also calculate overall nadirs (minimum across all 30 days)
overall_nadirs = 
  fgroup_by(hematox_labs, patient_id, lab_category) |>
  fsummarize(nadir_value = fmin(lab_value, na.rm = TRUE)) |>
  fungroup() |>
  pivot_wider(
    names_from   = lab_category,
    values_from  = nadir_value)

# Join overall nadirs to df
df = join(df, overall_nadirs, how = "left", multiple = FALSE)

# Rename for clarity
if ("neutrophils_absolute" %in% names(df)) {
  df = frename(df, anc_nadir = neutrophils_absolute)
}
if ("platelet_count" %in% names(df)) {
  df = frename(df, plt_nadir = platelet_count)
}
if ("hemoglobin" %in% names(df)) {
  df = frename(df, hgb_nadir = hemoglobin)
}
if ("lymphocytes_absolute" %in% names(df)) {
  df = frename(df, alc_nadir = lymphocytes_absolute)
}

rm(overall_nadirs, hematox_labs)
gc()

## hematotoxicity (ICAHT) - nadirs in day 0-30 ---------------------------------

# Get ANC daily nadirs and fill gaps (assume state persists between measurements)
anc_daily = 
  fsubset(daily_nadirs, lab_category == "neutrophils_absolute") |>
  fselect(patient_id, d_since_cart, anc = nadir_value) |>
  roworder(patient_id, d_since_cart)

# Fill gaps: for each patient, create complete day sequence 0-30 and fill forward
anc_complete = 
  anc_daily |>
  as.data.table() |>
  setorder(patient_id, d_since_cart)

# Create complete day grid for each patient
anc_complete = 
  anc_complete[, .(d_since_cart = 0:30), by = patient_id] |>
  merge(anc_daily, by = c("patient_id", "d_since_cart"), all.x = TRUE) |>
  setorder(patient_id, d_since_cart)

# Forward fill ANC values (LOCF - last observation carried forward)
anc_complete[, anc_filled := nafill(anc, type = "locf"), by = patient_id]

# Calculate consecutive days with ANC <500 and <100
# Function to find max consecutive days below threshold
find_max_consecutive = function(values, threshold) {
  below = values < threshold
  if (!any(below, na.rm = TRUE)) return(0L)
  
# Create run-length encoding
  rle_result = rle(below)
  max_consecutive = max(rle_result$lengths[rle_result$values], na.rm = TRUE)
  return(as.integer(max_consecutive))
}

# Calculate for each patient
icaht_calc = 
  anc_complete[!is.na(anc_filled), .(
    max_consec_lt500 = find_max_consecutive(anc_filled, 0.5),
    max_consec_lt100 = find_max_consecutive(anc_filled, 0.1),
    max_anc          = max(anc_filled, na.rm = TRUE),
    min_anc          = min(anc_filled, na.rm = TRUE),
    n_days_with_data = sum(!is.na(anc))
  ), by = patient_id]

# Apply ICAHT grading criteria (Early: 0-30 days)
# Grade 1: <500 for 1-6 consecutive days
# Grade 2: <500 for 7-13 consecutive days
# Grade 3: <500 for ≥14 consecutive days OR <100 for ≥7 consecutive days
# Grade 4: Never above 500 OR <100 for ≥14 consecutive days

icaht_calc = 
  icaht_calc |>
  as_tidytable() |>
  fmutate(
    icaht_grade = case_when(
      max_anc <= 0.5                                      ~ 4L,  # Never above 500
      max_consec_lt100 >= 14                              ~ 4L,  # <100 for ≥14 consecutive days
      max_consec_lt100 >= 7                               ~ 3L,  # <100 for ≥7 consecutive days
      max_consec_lt500 >= 14                              ~ 3L,  # <500 for ≥14 consecutive days
      max_consec_lt500 >= 7 & max_consec_lt500 < 14       ~ 2L,  # <500 for 7-13 consecutive days
      max_consec_lt500 >= 1 & max_consec_lt500 < 7        ~ 1L,  # <500 for 1-6 consecutive days
      TRUE                                                ~ 0L   # No significant neutropenia
    ),
    severe_icaht_01 = as.integer(icaht_grade >= 3)
  ) |>
  fselect(patient_id, icaht_grade, severe_icaht_01)

# Join to df
df = join(df, icaht_calc, how = "left", multiple = FALSE)

# Fill missing ICAHT grades with 0
df$icaht_grade       = replace_na(df$icaht_grade, 0L)
df$severe_icaht_01   = replace_na(df$severe_icaht_01, 0L)

## CTCAE grading for hematologic toxicities ------------------------------------
# Using overall nadirs, calculate worst CTCAE grades

df = df |>
  fmutate(
    neutropenia_grade = case_when(
      anc_nadir < 0.1  ~ 4L,
      anc_nadir < 0.5  ~ 3L,
      anc_nadir < 1.0 ~ 2L,
      anc_nadir < 1.5 ~ 1L,
      anc_nadir >= 1.5 ~ 0L,
      TRUE ~ NA_integer_),

    anemia_grade = case_when(
      hgb_nadir < 8.0  ~ 3L,
      hgb_nadir < 10.0 ~ 2L,
      hgb_nadir < 12.0 ~ 1L,
      hgb_nadir >= 12.0 ~ 0L,
      TRUE ~ NA_integer_),

    thrombocytopenia_grade = case_when(
      plt_nadir < 25  ~ 4L,
      plt_nadir < 50  ~ 3L,
      plt_nadir < 75  ~ 2L,
      plt_nadir < 150 ~ 1L,
      plt_nadir >= 150 ~ 0L,
      TRUE ~ NA_integer_)
  )

## IEC-HLH -------------------------------------------

# Step 1: Identify patients with ferritin >10,000
ferritin_hlh = 
  fsubset(labs_raw, lab_category == "ferritin" & d_since_cart >= 0 & d_since_cart <= 30) |>
  fgroup_by(patient_id) |>
  fsummarize(
    ferritin_max_30d = fmax(lab_value, na.rm = TRUE)) |>
  fungroup() |>
  fmutate(ferritin_10k_01 = as.integer(ferritin_max_30d >= 10000))

# Only evaluate HLH criteria for patients with ferritin >10,000
hlh_eligible = fsubset(ferritin_hlh, ferritin_10k_01 == 1)$patient_id

if (length(hlh_eligible) > 0) {
  
  # Get baseline values for comparison
  baseline_comparators = 
    df |>
    fselect(patient_id, 
            baseline_alt, baseline_ast, baseline_bilirubin_total, baseline_creatinine)
  
  # Prepare lab data for HLH criteria (days 0-30)
  hlh_labs = 
    fsubset(labs_raw, patient_id %in% hlh_eligible & d_since_cart >= 0 & d_since_cart <= 30) |>
    join(baseline_comparators, on = "patient_id", how = "left")
  
  # Define upper limits of normal (ULN)
  uln_alt = 40  
  uln_ast = 40  
  uln_bili = 1.2  
  uln_cr = 1.3   
  
  # Calculate criterion flags for each lab measurement
  hlh_labs = hlh_labs |>
    fmutate(
      # Liver dysfunction: >5x ULN or >5x baseline
      liver_alt = case_when(
        lab_category == "alt" & lab_value > 5 * uln_alt ~ 1L,
        lab_category == "alt" & !is.na(baseline_alt) & lab_value > 5 * baseline_alt ~ 1L,
        TRUE ~ 0L
      ),
      liver_ast = case_when(
        lab_category == "ast" & lab_value > 5 * uln_ast ~ 1L,
        lab_category == "ast" & !is.na(baseline_ast) & lab_value > 5 * baseline_ast ~ 1L,
        TRUE ~ 0L
      ),
      liver_bili = case_when(
        lab_category == "bilirubin_total" & lab_value > 5 * uln_bili ~ 1L,
        lab_category == "bilirubin_total" & !is.na(baseline_bilirubin_total) & 
          lab_value > 5 * baseline_bilirubin_total ~ 1L,
        TRUE ~ 0L
      ),
      liver_flag = pmax(liver_alt, liver_ast, liver_bili, na.rm = TRUE),
      
      # Renal dysfunction: >3x ULN or >3x baseline
      renal_flag = case_when(
        lab_category == "creatinine" & lab_value > 3 * uln_cr ~ 1L,
        lab_category == "creatinine" & !is.na(baseline_creatinine) & 
          lab_value > 3 * baseline_creatinine ~ 1L,
        TRUE ~ 0L
      ),
      
      # Cytopenia flags (will aggregate later)
      hgb_low   = fifelse(lab_category == "hemoglobin" & lab_value < 8, 1L, 0L),
      plt_low   = fifelse(lab_category == "platelet_count" & lab_value < 50, 1L, 0L),
      anc_low   = fifelse(lab_category == "neutrophils_absolute" & lab_value < 500, 1L, 0L)
    )
  
  # For each patient-day, determine if criteria are met
  hlh_daily = 
    hlh_labs |>
    fgroup_by(patient_id, d_since_cart) |>
    fsummarize(
      liver_crit = fmax(liver_flag, na.rm = TRUE),
      renal_crit = fmax(renal_flag, na.rm = TRUE),
      hgb_crit   = fmax(hgb_low, na.rm = TRUE),
      plt_crit   = fmax(plt_low, na.rm = TRUE),
      anc_crit   = fmax(anc_low, na.rm = TRUE)
    ) |>
    fungroup() |>
    fmutate(
      # Cytopenia in ≥2 lineages
      cytopenia_lineages = hgb_crit + plt_crit + anc_crit,
      cytopenia_crit = as.integer(cytopenia_lineages >= 2))
  
  # Look for concurrent criteria 
  # For each day, check if ≥2 criteria met within that day ± 3 days
  hlh_windows = 
    hlh_daily |>
    as.data.table() |>
    setorder(patient_id, d_since_cart)
  
  # For each patient, find if any 3-day window has ≥2 criteria
  hlh_result = hlh_windows[, {
    max_criteria_3d = 0L
    hlh_met = 0L
    
    if (.N > 0) {
      for (i in seq_len(.N)) {
        day_i = d_since_cart[i]
        window = d_since_cart >= (day_i - 1) & d_since_cart <= (day_i + 1)
        
        # Check if criteria are met in this window
        liver_in_window = any(liver_crit[window] == 1, na.rm = TRUE)
        renal_in_window = any(renal_crit[window] == 1, na.rm = TRUE)
        cyto_in_window  = any(cytopenia_crit[window] == 1, na.rm = TRUE)
        
        criteria_in_window = sum(liver_in_window, renal_in_window, cyto_in_window)
        max_criteria_3d = max(max_criteria_3d, criteria_in_window)
      }
      hlh_met = as.integer(max_criteria_3d >= 2)
    }
    
    .(hlh_criteria_met = hlh_met, max_concurrent_criteria = max_criteria_3d)
  }, by = patient_id]
  
  # Combine with ferritin data
  hlh_final = 
    join(ferritin_hlh, as_tidytable(hlh_result), on = "patient_id", how = "left") |>
    fmutate(
      hlh_criteria_met = replace_na(hlh_criteria_met, 0L),
      max_concurrent_criteria = replace_na(max_concurrent_criteria, 0L),
      hlh_01 = as.integer(ferritin_10k_01 == 1 & hlh_criteria_met == 1)
    )
  
  rm(hlh_labs, hlh_daily, hlh_windows, hlh_result, baseline_comparators)
  gc()
  
} else {
  # No patients with ferritin >10,000
  hlh_final = 
    ferritin_hlh |>
    fmutate(
      hlh_criteria_met = 0L,
      max_concurrent_criteria = 0L,
      hlh_01 = 0L
    )
}

# Join to main df
df = join(df, hlh_final |> fselect(patient_id, hlh_01), on = "patient_id", how = "left", multiple = FALSE)

# Fill missing values
df$hlh_01 = replace_na(df$hlh_01, 0L)


rm(ferritin_hlh, hlh_eligible, hlh_final, icaht_calc)
gc()



## lab trends by day (for export) ----------------------------------------------

labs_daily_mean = 
  fsubset(labs_raw, d_since_cart >= 0 & d_since_cart <= 30) |>
  fgroup_by(patient_id, lab_category, d_since_cart) |>
  fsummarize(lab_value = fmean(lab_value, na.rm = TRUE)) |>
  fungroup()

### add trajectory for grouping
labs_daily_mean = join(labs_daily_mean, fselect(df, patient_id, trajectory), how = "inner")

daily_nadirs = join(daily_nadirs, fselect(df, patient_id, trajectory), how = "inner")

### aggregate by trajectory and day
labs_trend_mean_traj = 
  fgroup_by(labs_daily_mean, trajectory, lab_category, d_since_cart) |>
  fsummarize(
    n       = fnobs(lab_value),
    sum_x   = fsum(lab_value, na.rm = TRUE),
    sum_x_sq = fsum(lab_value^2, na.rm = TRUE)
  ) |>
  fungroup() |>
  roworder(lab_category, trajectory, d_since_cart)

labs_trend_nadir_traj = 
  fgroup_by(daily_nadirs, trajectory, lab_category, d_since_cart) |>
  fsummarize(
    n       = fnobs(nadir_value),
    sum_x   = fsum(nadir_value, na.rm = TRUE),
    sum_x_sq = fsum(nadir_value^2, na.rm = TRUE)
  ) |>
  fungroup() |>
  roworder(lab_category, trajectory, d_since_cart)

### save for export (will be picked up in 06_export.R)
labs_trend_mean_export = labs_trend_mean_traj
labs_trend_nadir_export = labs_trend_nadir_traj

rm(labs_daily_mean, labs_trend_mean_traj, labs_raw, daily_nadirs, labs_trend_nadir_traj)
gc()


# medication outcomes ==========================================================

## tocilizumab -----------------------------------------------------------------

toci = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(tolower(med_category) == "tocilizumab") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm) |>
  dplyr::collect()

toci = 
  join(toci, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, toci_dttm = admin_dttm)

toci = 
  join(toci, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fsubset(toci_dttm >= cart_t0 & toci_dttm <= cart_t0 + ddays(30)) |>
  fmutate(d_since_cart = as.integer(floor(as.numeric(difftime(toci_dttm, cart_t0, units = "days"))))) |>
  funique()

# Daily dose counts
toci_daily = 
  fgroup_by(toci, patient_id, d_since_cart) |>
  fsummarize(n_doses = fnobs(toci_dttm)) |>
  fungroup()

toci_patient_daily = 
  toci_daily |>
  join(fselect(df, patient_id, trajectory), on = "patient_id", how = "inner") |>
  fselect(patient_id, trajectory, d_since_cart, toci_n_doses = n_doses)

meds_trend_toci = 
  fgroup_by(toci_patient_daily, trajectory, d_since_cart) |>
  fsummarize(
    n_patients  = fnobs(patient_id),
    sum_doses   = fsum(toci_n_doses, na.rm = TRUE),
    sum_doses_sq = fsum(toci_n_doses^2, na.rm = TRUE)
  ) |>
  fungroup() |>
  fmutate(medication = "tocilizumab") |>
  roworder(trajectory, d_since_cart)

# Overall summary
toci_summary = 
  roworder(toci, patient_id, toci_dttm) |>
  fgroup_by(patient_id) |>
  fsummarize(
    toci_first_dttm  = ffirst(toci_dttm),
    toci_n_doses_30d = fnobs(toci_dttm),
    cart_t0          = ffirst(cart_t0)
  ) |>
  fmutate(h_until_toci = as.numeric(difftime(toci_first_dttm, cart_t0, units = "hours"))) |>
  fselect(patient_id, toci_n_doses_30d, h_until_toci)

rm(toci)
gc()

## anakinra --------------------------------------------------------------------

anakinra = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(tolower(med_category) == "anakinra") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm) |>
  dplyr::collect()

anakinra = 
  join(anakinra, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, ana_dttm = admin_dttm)

anakinra = 
  join(anakinra, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fsubset(ana_dttm >= cart_t0 & ana_dttm <= cart_t0 + ddays(30)) |>
  fmutate(d_since_cart = as.integer(floor(as.numeric(difftime(ana_dttm, cart_t0, units = "days"))))) |>
  funique()

# Daily dose counts
anakinra_daily = 
  fgroup_by(anakinra, patient_id, d_since_cart) |>
  fsummarize(n_doses = fnobs(ana_dttm)) |>
  fungroup()

anakinra_patient_daily = 
  anakinra_daily |>
  join(fselect(df, patient_id, trajectory), on = "patient_id", how = "inner") |>
  fselect(patient_id, trajectory, d_since_cart, anakinra_n_doses = n_doses)

meds_trend_anakinra = 
  fgroup_by(anakinra_patient_daily, trajectory, d_since_cart) |>
  fsummarize(
    n_patients  = fnobs(patient_id),
    sum_doses   = fsum(anakinra_n_doses, na.rm = TRUE),
    sum_doses_sq = fsum(anakinra_n_doses^2, na.rm = TRUE)
  ) |>
  fungroup() |>
  fmutate(medication = "anakinra") |>
  roworder(trajectory, d_since_cart)

# Overall summary
anakinra_summary = 
  roworder(anakinra, patient_id, ana_dttm) |>
  fgroup_by(patient_id) |>
  fsummarize(
    anakinra_first_dttm  = ffirst(ana_dttm),
    anakinra_n_doses_30d = fnobs(ana_dttm),
    cart_t0              = ffirst(cart_t0)
  ) |>
  fmutate(h_until_anakinra = as.numeric(difftime(anakinra_first_dttm, cart_t0, units = "hours"))) |>
  fselect(patient_id, anakinra_n_doses_30d, h_until_anakinra)

rm(anakinra)

## acetaminophen ---------------------------------------------------------------

apap = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(tolower(med_category) == "acetaminophen") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm, med_dose) |>
  dplyr::collect()

apap = 
  join(apap, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, apap_dttm = admin_dttm, med_dose)

apap = 
  join(apap, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fsubset(apap_dttm >= cart_t0 & apap_dttm <= cart_t0 + ddays(3)) |>
  fmutate(d_since_cart = as.integer(floor(as.numeric(difftime(apap_dttm, cart_t0, units = "days"))))) |>
  funique()

# Daily cumulative doses (days 0-2)
apap_summary = 
  fgroup_by(apap, patient_id, d_since_cart) |>
  fsummarize(apap_cum = fsum(med_dose, na.rm = TRUE)) |>
  fungroup() |>
  pivot_wider(
    names_from   = d_since_cart,
    values_from  = apap_cum,
    names_prefix = "apap_day_"
  ) |>
  fmutate(
    apap_cum_0_2 = replace_na(apap_day_0, 0) +
                   replace_na(apap_day_1, 0) +
                   replace_na(apap_day_2, 0))

rm(apap)
gc()

## steroids --------------------------------------------------------------------

steroid_map = tidytable(
  med_category = c(
    "dexamethasone", "betamethasone", "methylprednisolone", 
    "prednisone", "prednisolone", "hydrocortisone"
  ),
  dex_factor = c(1.0, 1.25, 0.1875, 0.15, 0.15, 0.0375)
)

weight = 
  data_list$vitals |>
  dplyr::filter(vital_category == "weight_kg") |>
  dplyr::filter(recorded_dttm >= start_date) |>
  dplyr::select(hospitalization_id, recorded_dttm, weight_kg = vital_value) |>
  dplyr::collect() 

weight = 
  join(weight, hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
  fselect(patient_id, recorded_dttm, weight_kg)

weight = 
  join(weight, anchor_df, how = "inner", multiple = TRUE) |>
  fsubset(recorded_dttm >= cart_t0 - ddays(10)) |>
  fsubset(recorded_dttm <= cart_t0 + ddays(30)) |>
  roworder(recorded_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(weight_kg = ffirst(weight_kg))

weight = 
  join(weight, df, how = "right", multiple = FALSE) |>
  fmutate(weight_kg = fifelse(is.na(weight_kg), fmedian(weight_kg, na.rm = TRUE), weight_kg)) |>
  fselect(patient_id, weight_kg) 

steroid_admin = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(!is.na(med_dose) & med_group == "steroid") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm, med_category, med_dose, med_dose_unit) |>
  dplyr::collect() 

steroid_admin = 
  join(steroid_admin, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(-hospitalization_id) |>
  funique()

steroid_baseline = 
  join(steroid_admin, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fsubset(admin_dttm >= cart_t0 - dhours(24) & admin_dttm < cart_t0)

steroid_baseline =
  join(steroid_baseline, weight, on = "patient_id", how = "left", multiple = FALSE) |>
  fmutate(u = fifelse(med_dose_unit == "mg/kg/day", "mg/kg", med_dose_unit)) |>
  fmutate(dose_mg = fifelse(u == "mg", med_dose, med_dose*weight_kg)) |>
  fsubset(!is.na(dose_mg)) |>
  join(steroid_map, on = "med_category", how = "inner") |>
  fmutate(dex_mg_equivalent = dose_mg * dex_factor) |>
  fgroup_by(patient_id) |>
  fsummarize(baseline_steroid = fsum(dex_mg_equivalent, na.rm = TRUE)) |>
  fungroup()

steroid_post = 
  join(steroid_admin, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fsubset(admin_dttm >= cart_t0) |>
  fsubset(admin_dttm <= cart_t0 + ddays(30))

steroid_post = 
  join(steroid_post, weight, how = "left", multiple = FALSE) |>
  fmutate(u = fifelse(med_dose_unit == "mg/kg/day", "mg/kg", med_dose_unit)) |>
  fmutate(dose_mg = fifelse(u == "mg", med_dose, med_dose*weight_kg)) |>
  fsubset(!is.na(dose_mg)) |>
  ftransform(h_since = as.numeric(difftime(admin_dttm, cart_t0), "hours")) |>
  ftransform(d_since_cart = as.integer(floor(h_since/24))) |>
  fselect(patient_id, d_since_cart, cart_t0, admin_dttm, med_category, dose_mg) |>
  funique()

steroid_daily = 
  join(steroid_post, steroid_map, how = "inner") |>
  ftransform(dex_mg_equivalent = dose_mg*dex_factor) |>
  fgroup_by(patient_id, d_since) |>
  fsummarise(dex_mg_equivalent = fsum(dex_mg_equivalent, na.rm = TRUE)) |>
  fungroup() |>
  fsubset(dex_mg_equivalent >= 10)

steroid_patient_daily = 
  steroid_daily |>
  join(fselect(df, patient_id, trajectory), on = "patient_id", how = "inner") |>
  fselect(patient_id, trajectory, d_since_cart, steroid_dex_eq = dex_mg_equivalent)

meds_trend_steroid = 
  fgroup_by(steroid_patient_daily, trajectory, d_since_cart) |>
  fsummarize(
    n_patients  = fnobs(patient_id),
    sum_doses   = fsum(steroid_dex_eq, na.rm = TRUE),
    sum_doses_sq = fsum(steroid_dex_eq^2, na.rm = TRUE)
  ) |>
  fungroup() |>
  fmutate(medication = "steroid_dex_eq") |>
  roworder(trajectory, d_since_cart)

steroid_summary = 
  fgroup_by(steroid_daily, patient_id) |>
  fsummarize(
    steroid_any_72h  = as.integer(any(d_since_cart <= 2)),
    steroid_days_30d = fnobs(d_since_cart)
  )

##Medication trends day 0-30 by trajectory---------------------------------------
meds_trend_by_traj = 
  rbind(
    fselect(meds_trend_toci, trajectory, medication, d_since_cart, n_patients, sum_doses, sum_doses_sq),
    fselect(meds_trend_anakinra, trajectory, medication, d_since_cart, n_patients, sum_doses, sum_doses_sq),
    fselect(meds_trend_steroid, trajectory, medication, d_since_cart, n_patients, sum_doses, sum_doses_sq)
  )

## add meds to data frame ------------------------------------------------------

df = join(df, toci_summary,              how = "left", multiple = FALSE)
df = join(df, anakinra_summary,          how = "left", multiple = FALSE)
df = join(df, apap_summary,              how = "left", multiple = FALSE)
df = join(df, steroid_baseline,          how = "left", multiple = FALSE)
df = join(df, steroid_summary,                   how = "left", multiple = FALSE)

df$toci_n_doses_30d     = replace_na(df$toci_n_doses_30d,     0L)
df$anakinra_n_doses_30d = replace_na(df$anakinra_n_doses_30d, 0L)
df$steroid_days_30d     = replace_na(df$steroid_days_30d,     0L)
df$steroid_any_72h      = replace_na(df$steroid_any_72h,      0L)
df$baseline_steroid     = replace_na(df$baseline_steroid,     0L)
df$apap_day_0           = replace_na(df$apap_day_0,           0)
df$apap_day_1           = replace_na(df$apap_day_1,           0)
df$apap_day_2           = replace_na(df$apap_day_2,           0)
df$apap_cum_0_2         = replace_na(df$apap_cum_0_2,         0)

rm(weight, steroid_post, steroid_admin, steroid_map, toci_summary, anakinra_summary, apap_summary, steroid_baseline, steroid_summary, toci_daily, 
   toci_patient_daily, meds_trend_toci, anakinra_daily, anakinra_patient_daily, meds_trend_anakinra, steroid_daily, steroid_patient_daily, meds_trend_steroid)
gc()

# clean up =====================================================================

df = fselect(df, -discharge_category, -discharge_category_clean)

save_output(labs_trend_mean_export, "lab_trends_mean_by_trajectory", shareable = TRUE, subfolder = "trajectory")
save_output(labs_trend_nadir_export, "lab_trends_nadir_by_trajectory", shareable = TRUE, subfolder = "trajectory")
save_output(meds_trend_by_traj, "meds_daily_by_trajectory", shareable = TRUE, subfolder = "trajectory")

rm(labs_trend_mean_export,  labs_trend_nadir_export, meds_trend_by_traj)
gc()

message("04_outcomes.R complete.")
