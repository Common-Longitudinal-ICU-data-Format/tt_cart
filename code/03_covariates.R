# ==============================================================================
# 03_covariates.R
# Patient demographics, diagnoses, and comorbidities for CAR-T project
# Requires: df, hid_jid_crosswalk, anchor_df from prior scripts
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("hid_jid_crosswalk") || !exists("save_output")) {
  stop("Run 00_setup.R through 02_traj.R first.", call. = FALSE)
}

message("Starting 03_covariates.R")

# patient demographics =========================================================

patient = 
  data_list$patient |>
  dplyr::filter(patient_id %in% df$patient_id) |>
  dplyr::select(
    patient_id, 
    race_category, 
    ethnicity_category, 
    sex_category, 
    death_dttm
  ) |>
  dplyr::collect()

df = 
  join(df, patient, how = "left", multiple = TRUE) |>
  ftransform(
    race_category = case_when(
      tolower(race_category) %in% c("white", "asian")       ~ tolower(race_category),
      tolower(race_category) == "black or african american" ~ "black or african american",
      TRUE                                                  ~ "other/unknown"
    ),
    ethnicity_category = case_when(
      tolower(ethnicity_category) %in% c("hispanic", "non-hispanic") ~ tolower(ethnicity_category),
      TRUE                                                           ~ "unknown"
    ),
    sex_category = case_when(
      tolower(sex_category) %in% c("male", "female") ~ tolower(sex_category),
      TRUE                                           ~ "unknown"
    )
  )

rm(patient); gc()

# age at admission =============================================================

enc = 
  data_list$hospitalization |>
  dplyr::filter(patient_id %in% df$patient_id) |>
  dplyr::select(hospitalization_id, age = age_at_admission) |>
  dplyr::collect()

enc = 
  join(enc, hid_jid_crosswalk, how = "inner") |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(age = fmin(age))

df = join(df, enc, how = "left", multiple = FALSE)

rm(enc); gc()

# diagnoses ====================================================================

hids_of_interest = 
  fselect(df, joined_hosp_id) |>
  join(hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
  pull(hospitalization_id)

dx = 
  data_list$hospital_diagnosis |>
  dplyr::filter(hospitalization_id %in% hids_of_interest) |>
  dplyr::select(hospitalization_id, diagnosis_code) |>
  dplyr::collect()

dx = 
  join(dx, hid_jid_crosswalk, how = "inner") |> 
  fselect(-hospitalization_id, -joined_hosp_id) |>
  funique()

# cancer diagnoses =============================================================

expand_icd = function(prefix) {
  paste0(prefix, 0:9)
}

expand_icd_range = function(base, subtypes) {
  unlist(lapply(paste0(base, ".", subtypes), expand_icd))
}

dx_groups <- list(
  follicular         = expand_icd_range("C82", c("0","1","2","3","4","5","6","8","9")),
  scbcl              = expand_icd("C83.0"),                
  mantle_cell        = expand_icd("C83.1"),
  dlbcl              = expand_icd("C83.3"),
  lymphoblastic      = expand_icd("C83.5"),
  burkitt            = expand_icd("C83.7"),
  non_follicular_nos = expand_icd(c("C83.8","C83.9")),
  nhl_nos            = c(expand_icd("C88.4"), expand_icd_range("C85", c("1","2","8","9"))),
  myeloma            = expand_icd_range("C90", c("0","1","2","3")),
  b_all              = c("C91.00","C91.01","C91.02"),
  cll                = expand_icd(c("C91.10","C91.11","C91.12")),
  hct                = c("Z94.81","Z94.84", paste0("T86.0", 0:3), "T86.09","T86.5")
)

dx_lookup = 
  lapply(names(dx_groups), \(g) {
    tidytable(
      cancer_cat     = g,
      diagnosis_code = dx_groups[[g]]
    )
  }) |>
  rbindlist()

dx_cancer = join(dx, dx_lookup, how = "inner") 

## HCT flag --------------------------------------------------------------------

hct_vec = 
  fsubset(dx_cancer, cancer_cat == "hct") |> 
  pull(patient_id)

## determine primary cancer category -------------------------------------------

lymphoma_priority = tidytable(
  cancer_cat = c(
    "dlbcl",
    "mantle_cell",
    "burkitt",
    "lymphoblastic",
    "follicular",
    "scbcl",
    "cll",
    "non_follicular_nos",
    "nhl_nos"
  ),
  lym_rank = 1L:9L
)

cart_indication = tidytable(
  cart_name = c(
    "axicabtagene_ciloleucel",
    "brexucabtagene_autoleucel",
    "ciltacabtagene_autoleucel",
    "idecabtagene_vicleucel",
    "lisocabtagene_maraleucel",
    "tisagenlecleucel"
  ),
  expected_dx = list(
    c("lymphoma"),
    c("lymphoma", "b_all"),
    c("myeloma"),
    c("myeloma"),
    c("lymphoma"),
    c("lymphoma", "b_all")
  )
)

dx_cancer = 
  fsubset(dx_cancer, cancer_cat != "hct") |>
  fselect(-diagnosis_code) |>
  funique() |>
  fmutate(
    dx_broad = case_when(
      cancer_cat %in% c("myeloma", "b_all") ~ cancer_cat,
      TRUE                                  ~ "lymphoma"
    )
  )

dx_with_cart = 
  join(dx_cancer, fselect(df, patient_id, cart_name), how = "inner") |>
  join(cart_indication, on = "cart_name", how = "left") |>
  rowwise() |>
  mutate(indication_match = dx_broad %in% unlist(expected_dx)) |>
  ungroup() |>
  join(lymphoma_priority, on = "cancer_cat", how = "left") |>
  fmutate(lym_rank = fifelse(is.na(lym_rank), 0L, lym_rank))

dx_primary = 
  roworder(dx_with_cart, patient_id, -indication_match, lym_rank) |>
  fgroup_by(patient_id) |>
  ffirst() |>
  fsubset(indication_match) |>
  fselect(patient_id, cancer_cat)

df = join(df, dx_primary, how = "left")


df$hct_01 = as.integer(df$patient_id %in% hct_vec)

## fill in missing cancer_cat based on cart product ----------------------------

df$cancer_cat = case_when(
  !is.na(df$cancer_cat)                          ~ df$cancer_cat,
  df$cart_name == "ciltacabtagene_autoleucel"    ~ "myeloma",
  df$cart_name == "idecabtagene_vicleucel"       ~ "myeloma",
  TRUE                                           ~ "unknown"
)

df$dx_broad_final = case_when(
  df$cancer_cat %in% c("myeloma", "b_all", "unknown") ~ df$cancer_cat,
  TRUE                                     ~ "lymphoma"
)
  
rm(dx_primary, dx_with_cart, cart_indication)
rm(hct_vec, dx_cancer, dx_lookup, dx_groups, lymphoma_priority)
rm(expand_icd, expand_icd_range, hids_of_interest)
gc()

# elixhauser comorbidities =====================================================

elix = 
  comorbidity::comorbidity(
    dx, 
    id      = "patient_id", 
    code    = "diagnosis_code", 
    map     = "elixhauser_icd10_quan", 
    assign0 = TRUE
  )

## calculate van walraven scores -----------------------------------------------

vw   = comorbidity::score(elix, weights = "vw", assign0 = TRUE)
elix = cbind(elix, vw = vw) |> fselect(patient_id, vw)
df   = join(df, elix, how = "left", multiple = FALSE)

df$miss_vw = as.integer(is.na(df$vw))
df$vw      = fifelse(is.na(df$vw), 0L, df$vw)

rm(elix, vw, dx); gc()

message("03_covariates.R complete.")
