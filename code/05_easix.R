# ==============================================================================
# 05_easix.R
# EASIX and related severity scores for CAR-T project
# Requires: df, anchor_df, hid_jid_crosswalk from prior scripts
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("anchor_df") || !exists("save_output")) {
  stop("Run 00_setup.R through 04_outcomes.R first.", call. = FALSE)
}

message("Starting 05_easix.R")

# define thresholds ============================================================

thr_easix    = 4.6
thr_ferritin = 321L
thr_crp      = 21L

# load lab data ================================================================

easix_labs = c(
  "ldh", 
  "creatinine",
  "platelet_count",
  "crp", 
  "ferritin",
  "lymphocytes_absolute", 
  "wbc"
)

## load labs -------------------------------------------------------------------

lab_window = 
  data_list$labs |>
  dplyr::filter(tolower(lab_category) %in% easix_labs) |>
  dplyr::filter(lab_result_dttm >= start_date) |>
  dplyr::filter(tolower(lab_specimen_category) == "blood/plasma/serum") |>
  dplyr::filter(!is.na(lab_value_numeric)) |>
  dplyr::select(
    hospitalization_id, 
    lab_time = lab_result_dttm,
    lab_category,
    lab_value_numeric
  ) |>
  dplyr::collect() |>
  distinct()

lab_window = 
  join(lab_window, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fselect(-hospitalization_id, -joined_hosp_id) |>
  distinct()

lab_window = 
  join(lab_window, anchor_df, how = "inner", multiple = T) |>
  fsubset(lab_time >= cart_t0 - dhours(12)) |>
  fsubset(lab_time <= cart_t0 + dhours(12)) |>
  roworder(lab_time) |>
  fgroup_by(patient_id, lab_category) |>
  ffirst()

# calculate scores =============================================================

scores = 
  fselect(lab_window, patient_id, lab_category, lab_value_numeric) |>
  pivot_wider(names_from = lab_category, values_from = lab_value_numeric) |>
  fmutate(
    easix        = safe_div(creatinine*ldh, platelet_count),
    m_easix      = safe_div(crp*ldh, platelet_count),
    s_easix      = safe_div(ldh, platelet_count),
    log2_easix   = safe_log2(easix),
    log2_m_easix = safe_log2(m_easix),
    log2_s_easix = safe_log2(s_easix),
    easix_f_cat  = case_when(
      !is.na(easix) & easix > thr_easix                                                ~ "High",
      !is.na(easix) & easix <= thr_easix & !is.na(ferritin) & ferritin > thr_ferritin  ~ "Intermediate",
      !is.na(easix) & easix <= thr_easix & !is.na(ferritin) & ferritin <= thr_ferritin ~ "Low",
      !is.na(easix) & easix <= thr_easix & is.na(ferritin)                             ~ NA_character_,
      TRUE ~ NA_character_
    ),
    easix_fc_cat = case_when(
      !is.na(ferritin) & ferritin > thr_ferritin                                      ~ "High",
      !is.na(ferritin) & ferritin <= thr_ferritin & !is.na(easix) & easix > thr_easix ~ "Intermediate",
      !is.na(ferritin) & ferritin <= thr_ferritin & !is.na(crp) & crp > thr_crp       ~ "Intermediate",
      !is.na(ferritin) & ferritin <= thr_ferritin & !is.na(easix) & easix <= thr_easix & 
        !is.na(crp) & crp <= thr_crp                                                  ~ "Low",
      TRUE                                                                            ~ NA_character_
    ),
    easix_f_cat  = factor(easix_f_cat,  levels = c("Low", "Intermediate", "High"), ordered = TRUE),
    easix_fc_cat = factor(easix_fc_cat, levels = c("Low", "Intermediate", "High"), ordered = TRUE)
  ) |>
  fselect(
    patient_id,
    easix, 
    log2_easix,
    easix_f_cat, 
    easix_fc_cat,
    wbc_infusion       = wbc, 
    lymph_abs_infusion = lymphocytes_absolute
  ) |>
  distinct()

# add to cohort df =============================================================

df = join(df, scores, how = "left", multiple = FALSE)

rm(thr_easix, thr_ferritin, thr_crp, lab_window, scores, easix_labs); gc()

message("05_easix.R complete.")
