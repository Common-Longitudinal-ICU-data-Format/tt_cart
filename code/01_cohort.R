# ==============================================================================
# 01_cohort.R
# Cohort identification and exclusions for CAR-T temperature trajectory project
# Outputs: exclusion flow diagram data, cohort df, temp df, crosswalks
# ==============================================================================

# check setup ran ==============================================================

if (!exists("data_list") || !exists("save_output")) {
  stop("Run 00_setup.R first.", call. = FALSE)
}

message("Starting 01_cohort.R")


# link contiguous hospitalizations =============================================

## encounters with age >= 18 and dates within study window ---------------------

hosp_blocks = 
  dplyr::filter(data_list$hospitalization, age_at_admission >= 18) |>
  dplyr::filter(admission_dttm >= start_date & admission_dttm <= end_date) |> 
  dplyr::filter(admission_dttm < discharge_dttm & !is.na(discharge_dttm)) |>
  dplyr::select(patient_id, hospitalization_id, admission_dttm, discharge_dttm) |>
  dplyr::collect() |>
  roworder(patient_id, admission_dttm)

## use data.table to find joined hospitalizations with <= 6h gaps --------------

link_hours = 6L
linked     = as.data.table(hosp_blocks)
setorder(linked, patient_id, admission_dttm)

linked[, next_admit := shift(admission_dttm, type = "lead"), by = patient_id]
linked[, next_gap   := as.numeric(difftime(next_admit, discharge_dttm, units = "hours"))]
linked[, prev_dc    := shift(discharge_dttm, type = "lag"),  by = patient_id]  
linked[, prev_gap   := as.numeric(difftime(admission_dttm, prev_dc, units = "hours"))]
linked[, link_flag  := (next_gap < link_hours | prev_gap < link_hours)]
linked[is.na(link_flag), link_flag := FALSE]

linked[, new_group     := is.na(prev_gap) | prev_gap >= link_hours]
linked[, joined_hosp_id := .GRP, by = .(patient_id, cumsum(new_group))]

## create crosswalk ------------------------------------------------------------

hid_jid_crosswalk = 
  fselect(linked, patient_id, hospitalization_id, joined_hosp_id) |> 
  as_tidytable()

rm(hosp_blocks, linked, link_hours)
gc()


# identify CAR-T infusions =====================================================

cart = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(med_group == "car_t") |>
  dplyr::filter(admin_dttm >= start_date & admin_dttm <= end_date) |>
  dplyr::mutate(
    cart_name = tolower(med_category),
    cart_group = dplyr::case_when(
      cart_name == "axicabtagene_ciloleucel"   ~ 1L,
      cart_name == "brexucabtagene_autoleucel" ~ 2L,
      cart_name == "ciltacabtagene_autoleucel" ~ 3L,
      cart_name == "idecabtagene_vicleucel"    ~ 4L,
      cart_name == "lisocabtagene_maraleucel"  ~ 5L,
      cart_name == "tisagenlecleucel"          ~ 6L,
      TRUE                                     ~ NA_integer_
    )
  ) |>
  dplyr::distinct(hospitalization_id, admin_dttm, cart_name, cart_group) |>
  dplyr::collect()

cart = join(hid_jid_crosswalk, cart, how = "inner") 

cart_pids = funique(cart$patient_id)

cart =
  fselect(cart, joined_hosp_id, admin_dttm, cart_name, cart_group) |>
  funique()


# build encounter-level data ===================================================

enc = 
  data_list$hospitalization |>
  dplyr::filter(patient_id %in% cart_pids & admission_dttm >= start_date) |>
  dplyr::filter(age_at_admission >= 18) |>
  dplyr::select(
    ends_with("_id"), 
    ends_with("dttm"), 
    age = age_at_admission, 
    discharge_category
  ) |>
  dplyr::collect() 

enc = 
  join(enc, hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
  roworder(admission_dttm, discharge_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    patient_id         = ffirst(patient_id),
    admission_dttm     = ffirst(admission_dttm),
    discharge_dttm     = flast(discharge_dttm),
    discharge_category = tolower(flast(discharge_category))
  )

## join to first CAR-T per patient ---------------------------------------------

df =
  join(enc, cart, how = "inner") |>
  roworder(patient_id, admin_dttm) |>
  fmutate(cart_t0 = admin_dttm) |>
  fgroup_by(patient_id) |>
  ffirst() |>
  fungroup() |>
  fselect(-admin_dttm)

rm(cart, enc)
gc()


# exclusion tracking ===========================================================

exclusions = tidytable(
  step        = character(),
  reason      = character(),
  n_excluded  = integer(),
  n_remaining = integer()
)

n_start = fnrow(df)

exclusions = rbind(
  tidytable(
    step        = "0",
    reason      = "adult car-t infusions in date range",
    n_excluded  = 0L,
    n_remaining = n_start
  ),
  exclusions
)


# exclusion 1: < 24h in-hospital post-cart =====================================

df$time_to_dc = as.numeric(difftime(df$discharge_dttm, df$cart_t0), "hours")

short_stay    = fsubset(df, time_to_dc < 24)
n_short_death = fnrow(fsubset(short_stay, discharge_category %in% c("expired", "hospice")))
n_short_other = fnrow(short_stay) - n_short_death

exclusions = rbind(
  exclusions,
  tidytable(
    step        = "1a",
    reason      = "< 24h post-infusion: death or hospice",
    n_excluded  = n_short_death,
    n_remaining = n_start - n_short_death
  ),
  tidytable(
    step        = "1b",
    reason      = "< 24h post-infusion: early discharge",
    n_excluded  = n_short_other,
    n_remaining = NA_integer_
  )
)

df = 
  fsubset(df, time_to_dc >= 24) |>
  fselect(-time_to_dc)

exclusions[step == "1b", n_remaining := fnrow(df)]

rm(short_stay, n_short_death, n_short_other)
gc()


# load temperature data (needed for exclusions 1c-d and 2) =====================

cart_pats  = funique(df$patient_id)
all_c_hids = fsubset(hid_jid_crosswalk, patient_id %in% cart_pats) |> pull(hospitalization_id)

temp = 
  data_list$vitals |>
  dplyr::filter(vital_category == "temp_c") |>
  dplyr::filter(hospitalization_id %in% all_c_hids) |>
  dplyr::filter(recorded_dttm >= start_date) |>
  dplyr::filter(vital_value > 30 & vital_value < 45) |>
  dplyr::select(hospitalization_id, recorded_dttm, temp_c = vital_value) |>
  dplyr::collect()

temp = 
  join(temp, hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
  fselect(-hospitalization_id) |>
  join(df, how = "inner") |>
  fsubset(recorded_dttm >= admission_dttm) |>
  fsubset(recorded_dttm <= discharge_dttm) |>
  fselect(joined_hosp_id, patient_id, cart_t0, recorded_dttm, temp_c) |>
  funique()

rm(all_c_hids, cart_pats)
gc()


# exclusion 1c-d: < 24h temperature data post-infusion =========================

has_temp_24h = 
  fsubset(temp, recorded_dttm >= cart_t0 + dhours(24)) |>
  fselect(joined_hosp_id) |>
  funique()

n_before_temp = fnrow(df)
no_temp_24h   = fsubset(df, !joined_hosp_id %in% has_temp_24h$joined_hosp_id)
no_temp_death = fnrow(fsubset(no_temp_24h, discharge_category %in% c("expired", "hospice")))
no_temp_other = fnrow(no_temp_24h) - no_temp_death

df = fsubset(df, joined_hosp_id %in% has_temp_24h$joined_hosp_id)

exclusions = rbind(
  exclusions,
  tidytable(
    step        = "1c",
    reason      = "< 24h temp data: death or hospice",
    n_excluded  = no_temp_death,
    n_remaining = n_before_temp - no_temp_death
  ),
  tidytable(
    step        = "1d",
    reason      = "< 24h temp data: other/unknown",
    n_excluded  = no_temp_other,
    n_remaining = fnrow(df)
  )
)

rm(no_temp_24h, no_temp_death, no_temp_other, n_before_temp, has_temp_24h)
gc()

# finalize exclusions and filter temp ==========================================
## keep only temp for included patients ----------------------------------------

temp = fsubset(temp, joined_hosp_id %in% df$joined_hosp_id)

## print and save exclusion summary --------------------------------------------

print(exclusions)

save_output(exclusions, "exclusion_flow", shareable = TRUE, subfolder = "trajectory")


# add time of day (local) ======================================================

df$admin_localtime   = with_tz(df$cart_t0, tz = tz_local)
df$cart_hour_precise = hour(df$admin_localtime) + minute(df$admin_localtime)/60
df$cart_hour_floor   = floor(df$cart_hour_precise)
df$cart_hour_round   = round(df$cart_hour_precise)

df = fselect(df, -admin_localtime)

# create anchor df for downstream joins ========================================

anchor_df = fselect(df, patient_id, joined_hosp_id, cart_t0)

# clean up =====================================================================

rm(n_start, exclusions); gc()

message("01_cohort.R complete. Objects: df, temp, hid_jid_crosswalk, anchor_df")
