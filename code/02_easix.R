


# EASIX severity scores --------------------------------------------------------

safe_div  = function(num, den) ifelse(!is.na(num) & !is.na(den) & den > 0, num / den, NA_real_)
safe_log2 = function(x) ifelse(!is.na(x) & x > 0, log2(x), NA_real_)

thr_easix    = 4.6
thr_ferritin = 321
thr_crp      = 21

easix_labs = c(
  "ldh", 
  "creatinine",
  "platelet_count",
  "crp", 
  "ferritin",
  "lymphocytes_absolute", 
  "wbc"
)

lab_window = 
  data_list$labs |>
  dplyr::filter(lab_category %in% easix_labs) |>
  dplyr::filter(lab_specimen_category == "blood/plasma/serum") |>
  dplyr::filter(hospitalization_id %in% admissions$hospitalization_id) |>
  dplyr::select(
    hospitalization_id, 
    lab_time = lab_result_dttm,
    lab_category,
    lab_value_numeric
  ) |>
  dplyr::collect() |>
  distinct()

scores = 
  join(lab_window, admissions, how = "inner", multiple = T) |>
  fsubset(lab_time >= cart_t0 - lubridate::dhours(12)) |>
  fsubset(lab_time <= cart_t0 + lubridate::dhours(12)) |>
  roworder(lab_time) |>
  fgroup_by(hospitalization_id, lab_category) |>
  ffirst() |>
  fselect(patient_id, lab_category, lab_value_numeric) |>
  pivot_wider(names_from = lab_category, values_from = lab_value_numeric) |>
  fmutate(
    easix        = safe_div(creatinine*ldh, platelet_count),
    m_easix      = safe_div(crp*ldh, platelet_count),
    s_easix      = safe_div(ldh, platelet_count),
    log2_easix   = safe_log2(easix),
    log2_m_easix = safe_log2(m_easix),
    log2_s_easix = safe_log2(s_easix),
    easix_f_cat = case_when(
      !is.na(easix) & easix > thr_easix ~ "High",
      !is.na(easix) & easix <= thr_easix & !is.na(ferritin) & ferritin > thr_ferritin ~ "Intermediate",
      !is.na(easix) & easix <= thr_easix & !is.na(ferritin) & ferritin <= thr_ferritin ~ "Low",
      !is.na(easix) & easix <= thr_easix & is.na(ferritin) ~ NA_character_,  # explicit: can't classify without ferritin
      TRUE ~ NA_character_
    ),
    easix_fc_cat = case_when(
      !is.na(ferritin) & ferritin > thr_ferritin ~ "High",
      !is.na(ferritin) & ferritin <= thr_ferritin & !is.na(easix) & easix > thr_easix ~ "Intermediate",
      !is.na(ferritin) & ferritin <= thr_ferritin & !is.na(crp) & crp > thr_crp ~ "Intermediate",
      !is.na(ferritin) & ferritin <= thr_ferritin & 
        !is.na(easix) & easix <= thr_easix & 
        !is.na(crp) & crp <= thr_crp ~ "Low",
      TRUE ~ NA_character_
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
  funique()

cohort = left_join(cohort, scores, by = "patient_id")

rm(safe_div, safe_log2, thr_easix, thr_ferritin, thr_crp, lab_window, scores, easix_labs); gc()
