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
  "hemoglobin", "platelet_count", "neutrophils_absolute"
)

## load all relevant labs ------------------------------------------------------

message("Loading lab data...")

labs_raw = 
  data_list$labs |>
  dplyr::filter(tolower(lab_category) %in% lab_categories) |>
  dplyr::filter(lab_result_dttm >= start_date) |>
  dplyr::filter(!is.na(lab_value_numeric)) |>
  dplyr::select(hospitalization_id, lab_result_dttm, lab_category, lab_value_numeric) |>
  dplyr::collect() |>
  distinct()

labs_raw = 
  join(labs_raw, hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
  fselect(-hospitalization_id, -joined_hosp_id) |>
  distinct()

labs_raw = 
  join(labs_raw, anchor_df, how = "inner", multiple = TRUE) |>
  fmutate(
    d_since_cart = as.integer(floor(as.numeric(difftime(lab_result_dttm, cart_t0, units = "days"))))
  ) |>
  frename(lab_dttm = lab_result_dttm, lab_value = lab_value_numeric)


## baseline labs (day -1 to day 0, closest to infusion) ------------------------

message("Calculating baseline labs...")

baseline_labs = 
  fsubset(labs_raw, d_since_cart >= -1 & d_since_cart <= 0) |>
  fmutate(time_to_cart = abs(as.numeric(difftime(lab_dttm, cart_t0, units = "hours")))) |>
  roworder(patient_id, lab_category, time_to_cart) |>
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


## hematotoxicity (ICAHT) - nadirs in day 0-30 ---------------------------------

message("Calculating hematotoxicity nadirs...")

hematox_labs = fsubset(labs_raw, d_since_cart >= 0 & d_since_cart <= 30)

### ANC nadir
anc_nadir = 
  fsubset(hematox_labs, lab_category == "neutrophils_absolute") |>
  fgroup_by(patient_id) |>
  fsummarize(anc_nadir = fmin(lab_value, na.rm = TRUE)) |>
  fungroup()

### platelet nadir
plt_nadir = 
  fsubset(hematox_labs, lab_category == "platelet_count") |>
  fgroup_by(patient_id) |>
  fsummarize(plt_nadir = fmin(lab_value, na.rm = TRUE)) |>
  fungroup()

### hemoglobin nadir
hgb_nadir = 
  fsubset(hematox_labs, lab_category == "hemoglobin") |>
  fgroup_by(patient_id) |>
  fsummarize(hgb_nadir = fmin(lab_value, na.rm = TRUE)) |>
  fungroup()

### combine nadirs
df = join(df, anc_nadir, how = "left", multiple = FALSE)
df = join(df, plt_nadir, how = "left", multiple = FALSE)
df = join(df, hgb_nadir, how = "left", multiple = FALSE)

### ICAHT severity (simplified - based on ANC and platelets)
### Grade 1: ANC <1.0 or Plt <75
### Grade 2: ANC <0.5 or Plt <50
### Grade 3: ANC <0.2 or Plt <25
### Grade 4: ANC <0.1 or Plt <10

df = df |>
  fmutate(
    # ANC in 10^9/L, platelets in 10^9/L (CLIF standard)
    icaht_grade = case_when(
      anc_nadir < 0.1 | plt_nadir < 10   ~ 4L,
      anc_nadir < 0.2 | plt_nadir < 25   ~ 3L,
      anc_nadir < 0.5 | plt_nadir < 50   ~ 2L,
      anc_nadir < 1.0 | plt_nadir < 75   ~ 1L,
      TRUE                               ~ 0L
    ),
    severe_icaht_01 = as.integer(icaht_grade >= 3)
  )

rm(anc_nadir, plt_nadir, hgb_nadir, hematox_labs)
gc()


## hyperferritinemia (HLH surrogate) -------------------------------------------

message("Calculating hyperferritinemia...")

ferritin_max = 
  fsubset(labs_raw, lab_category == "ferritin" & d_since_cart >= 0 & d_since_cart <= 30) |>
  fgroup_by(patient_id) |>
  fsummarize(ferritin_max_30d = fmax(lab_value, na.rm = TRUE)) |>
  fungroup()

df = join(df, ferritin_max, how = "left", multiple = FALSE)

df = df |>
  fmutate(
    hyperferritin_5k_01  = as.integer(!is.na(ferritin_max_30d) & ferritin_max_30d >= 5000),
    hyperferritin_10k_01 = as.integer(!is.na(ferritin_max_30d) & ferritin_max_30d >= 10000)
  )

rm(ferritin_max)
gc()


## lab trends by day (for export) ----------------------------------------------

message("Calculating daily lab trends...")

labs_daily = 
  fsubset(labs_raw, d_since_cart >= 0 & d_since_cart <= 30) |>
  fgroup_by(patient_id, lab_category, d_since_cart) |>
  fsummarize(lab_value = fmean(lab_value, na.rm = TRUE)) |>
  fungroup()

### add trajectory for grouping
labs_daily = join(labs_daily, fselect(df, patient_id, trajectory), how = "inner")

### aggregate by trajectory and day
labs_trend_by_traj = 
  fgroup_by(labs_daily, trajectory, lab_category, d_since_cart) |>
  fsummarize(
    n       = fnobs(lab_value),
    sum_x   = fsum(lab_value, na.rm = TRUE),
    sum_x_sq = fsum(lab_value^2, na.rm = TRUE)
  ) |>
  fungroup() |>
  roworder(lab_category, trajectory, d_since_cart)

### save for export (will be picked up in 06_export.R)
labs_trend_export = labs_trend_by_traj

rm(labs_daily, labs_trend_by_traj, labs_raw)
gc()


# medication outcomes ==========================================================

## tocilizumab -----------------------------------------------------------------

toci = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(tolower(med_category) == "tocilizumab") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm) |>
  dplyr::collect() |>
  join(hid_jid_crosswalk, how = "inner", multiple = FALSE) |>
  fselect(patient_id, toci_dttm = admin_dttm)

toci = 
  join(toci, anchor_df, how = "inner", multiple = TRUE) |>
  fsubset(toci_dttm >= cart_t0 & toci_dttm <= cart_t0 + ddays(30)) |>
  ftransform(one = 1L) |>
  funique() |>
  roworder(toci_dttm) |>
  fgroup_by(patient_id) |>
  fsummarize(
    tfirst           = ffirst(toci_dttm),
    toci_n_doses_30d = fsum(one),
    cart_t0          = ffirst(cart_t0)
  ) |>
  fmutate(h_until_toci = as.numeric(difftime(tfirst, cart_t0, units = "hours"))) |>
  fselect(patient_id, toci_n_doses_30d, h_until_toci)

## anakinra --------------------------------------------------------------------

anakinra = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(tolower(med_category) == "anakinra") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm) |>
  dplyr::collect() |>
  join(hid_jid_crosswalk, how = "inner", multiple = FALSE) |>
  fselect(patient_id, ana_dttm = admin_dttm)

anakinra = 
  join(anakinra, anchor_df, how = "inner", multiple = TRUE) |>
  fsubset(ana_dttm >= cart_t0 & ana_dttm <= cart_t0 + ddays(30)) |>
  ftransform(one = 1L) |>
  funique() |>
  roworder(ana_dttm) |>
  fgroup_by(patient_id) |>
  fsummarize(
    afirst               = ffirst(ana_dttm),
    anakinra_n_doses_30d = fsum(one),
    cart_t0              = ffirst(cart_t0)
  ) |>
  fmutate(h_until_anakinra = as.numeric(difftime(afirst, cart_t0, units = "hours"))) |>
  fselect(patient_id, anakinra_n_doses_30d, h_until_anakinra)

## acetaminophen ---------------------------------------------------------------

apap = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(tolower(med_category) == "acetaminophen") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm, med_dose) |>
  dplyr::collect() |>
  join(hid_jid_crosswalk, how = "inner", multiple = FALSE) |>
  fselect(patient_id, apap_dttm = admin_dttm, med_dose)

apap = 
  join(apap, anchor_df, how = "inner", multiple = TRUE) |>
  fsubset(apap_dttm >= cart_t0 & apap_dttm <= cart_t0 + ddays(3)) |>
  funique() |>
  roworder(apap_dttm) |>
  fmutate(h_since = as.numeric(difftime(apap_dttm, cart_t0), "hours")) |>
  fmutate(d_since = as.integer(floor(h_since/24))) |>
  fgroup_by(patient_id, d_since) |>
  fsummarize(apap_cum = fsum(med_dose, na.rm = TRUE)) |>
  pivot_wider(
    names_from   = d_since,
    values_from  = apap_cum,
    names_prefix = "apap_day_"
  ) 

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
  join(steroid_admin, hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
  fselect(-hospitalization_id, -joined_hosp_id) |>
  funique()

steroid_admin = 
  join(steroid_admin, anchor_df, how = "inner", multiple = TRUE) |>
  fsubset(admin_dttm >= cart_t0) |>
  fsubset(admin_dttm <= cart_t0 + ddays(30))

steroid_admin = 
  join(steroid_admin, weight, how = "left", multiple = FALSE) |>
  fmutate(u = fifelse(med_dose_unit == "mg/kg/day", "mg/kg", med_dose_unit)) |>
  fmutate(dose_mg = fifelse(u == "mg", med_dose, med_dose*weight_kg)) |>
  fsubset(!is.na(dose_mg)) |>
  ftransform(h_since = as.numeric(difftime(admin_dttm, cart_t0), "hours")) |>
  ftransform(d_since = as.integer(floor(h_since/24))) |>
  fselect(patient_id, d_since, cart_t0, admin_dttm, med_category, dose_mg) |>
  funique()

steroid_daily = 
  join(steroid_admin, steroid_map, how = "inner") |>
  ftransform(dex_mg_equivalent = dose_mg*dex_factor) |>
  fgroup_by(patient_id, d_since) |>
  fsummarise(dex_mg_equivalent = fsum(dex_mg_equivalent, na.rm = TRUE)) |>
  fungroup() |>
  fsubset(dex_mg_equivalent >= 10)

steroid_daily = 
  fgroup_by(steroid_daily, patient_id) |>
  fsummarize(
    steroid_any_72h  = as.integer(any(d_since <= 2)),
    steroid_days_30d = fnobs(d_since)  
  )

steroid = 
  roworder(steroid_admin, admin_dttm) |>
  fgroup_by(patient_id) |>
  fsummarize(steroid_first = ffirst(admin_dttm)) |>
  join(steroid_daily, how = "full", multiple = FALSE)

## add meds to data frame ------------------------------------------------------

df = join(df, toci,     how = "left", multiple = FALSE)
df = join(df, anakinra, how = "left", multiple = FALSE)
df = join(df, apap,     how = "left", multiple = FALSE)
df = join(df, steroid,  how = "left", multiple = FALSE)

df$toci_n_doses_30d     = replace_na(df$toci_n_doses_30d,     0L)
df$anakinra_n_doses_30d = replace_na(df$anakinra_n_doses_30d, 0L)
df$steroid_days_30d     = replace_na(df$steroid_days_30d,     0L)
df$apap_day_0           = replace_na(df$apap_day_0,           0)
df$apap_day_1           = replace_na(df$apap_day_1,           0)
df$apap_day_2           = replace_na(df$apap_day_2,           0)
df$steroid_any_72h      = replace_na(df$steroid_any_72h,      0L)
df$h_until_steroids     = as.numeric(difftime(df$steroid_first, df$admin_dttm), "hours")

df = fselect(df, -steroid_first)

rm(toci, anakinra, apap, weight, steroid, steroid_admin, steroid_daily, steroid_map)
gc()


# antibiotics and cultures =====================================================

abx_list = c(
  "amikacin", "aztreonam", "cefepime", "cefiderocol", "ceftaroline", "ceftazidime",
  "ceftazidime_avibactam", "ceftazidime_clavulanate", "ceftolozane_tazobactam",
  "colistin", "dalbavancin", "daptomycin", "doripenem", "ertapenem", "fidaxomicin",
  "gentamicin", "imipenem", "imipenem_relebactam", "linezolid", "meropenem",
  "meropenem_vaborbactam", "metronidazole", "piperacillin_tazobactam", "tedizolid",
  "telavancin", "tigecycline", "tobramycin", "vancomycin"
)

abx = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(tolower(med_category) %in% abx_list) |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, med_category, admin_dttm) |>
  dplyr::collect() 

abx = 
  join(abx, hid_jid_crosswalk, how = "inner") |>
  fselect(joined_hosp_id, admin_dttm, med_category) |>
  funique()

abx = 
  join(abx, anchor_df, how = "inner") |>
  fsubset(admin_dttm >= cart_t0 & admin_dttm <= cart_t0 + ddays(30)) |>
  roworder(admin_dttm) |>
  fgroup_by(patient_id, med_category) |>
  fsummarize(admin_dttm = ffirst(admin_dttm))

df$abx_30_01 = fifelse(df$patient_id %in% abx$patient_id, 1L, 0L)

rm(abx, abx_list)
gc()

## blood cultures --------------------------------------------------------------

blood_cx = 
  data_list$microbiology_culture |>
  dplyr::filter(tolower(fluid_category) == "blood_buffy") |>
  dplyr::filter(tolower(organism_category) != "no_growth") |>
  dplyr::select(hospitalization_id, result_dttm, organism_category) |>
  dplyr::collect() 

blood_cx = 
  join(blood_cx, hid_jid_crosswalk, how = "inner") |>
  fselect(joined_hosp_id, patient_id, result_dttm, organism_category) |>
  funique()

blood_cx = 
  join(blood_cx, anchor_df, how = "inner") |>
  fsubset(result_dttm >= cart_t0 & result_dttm <= cart_t0 + ddays(30)) |>
  roworder(result_dttm) |>
  fgroup_by(patient_id, organism_category) |>
  fsummarize(result_dttm = ffirst(result_dttm))

df$bsi_30_01 = fifelse(df$patient_id %in% blood_cx$patient_id, 1L, 0L)

blood_cx_traj = 
  join(blood_cx, fselect(df, patient_id, trajectory), how = "inner") |>
  fgroup_by(trajectory, organism_category) |>
  fsummarize(n = fnobs(patient_id)) |>
  fungroup() |>
  roworder(trajectory, -n)

save_output(blood_cx_traj, "blood_culture_bugs_by_trajectory", shareable = TRUE, subfolder = "trajectory")

rm(blood_cx, blood_cx_traj)
gc()


# clean up =====================================================================

df = fselect(df, -discharge_category, -discharge_category_clean)

save_output(labs_trend_export, "lab_trends_by_trajectory", shareable = TRUE, subfolder = "trajectory")

rm(labs_trend_export)
gc()

message("04_outcomes.R complete.")
