
# Cohort script for CLIF project validating risk tools in oncology.
# Requires data_list to be loaded/validated from 00_*

# check/create output directories ---------------------------------------------

if (!exists("project_location")) {
  stop("project_location not found. Did you run 00_setup first?", call. = FALSE)
}

if (!dir.exists(paste0(project_location, "/upload_to_box"))) {
  dir.create(paste0(project_location, "/upload_to_box"), recursive = TRUE)
}

if (!dir.exists(paste0(project_location, "/proj_tables"))) {
  dir.create(paste0(project_location, "/proj_tables"), recursive = TRUE)
}

# resources for RAM heavy wrangling --------------------------------------------

## cores & RAM (reuse from 00 if available) ------------------------------------

os_type   = if (exists("os_type"))   os_type   else Sys.info()[["sysname"]]
all_cores = if (exists("all_cores")) all_cores else {
  x = parallel::detectCores(logical = TRUE); if (is.na(x)) 1L else as.integer(x)
}

if (!exists("avail_ram_gb") || !is.finite(avail_ram_gb)) {
  get_ram_gb = function() {
    tryCatch({
      if (Sys.info()[["sysname"]] == "Darwin") {
        bytes = suppressWarnings(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)))
        if (length(bytes) > 0 && !is.na(bytes)) bytes / 1024^3 else NA_real_
      } else if (file.exists("/proc/meminfo")) {
        kb = suppressWarnings(as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)))
        if (length(kb) > 0 && !is.na(kb)) kb / 1024^2 else NA_real_
      } else if (requireNamespace("ps", quietly = TRUE)) {
        ps::ps_system_memory()[["available"]] / 1024^3
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  avail_ram_gb = get_ram_gb()
}

## choose threads (2 GB per thread; leave 1 core free) -------------------------

reserve_cores  = 1L
gb_per_thread  = 2.0
max_by_cores   = max(1L, all_cores - reserve_cores)
max_by_memory  = if (is.finite(avail_ram_gb)) max(1L, floor(avail_ram_gb / gb_per_thread)) else max_by_cores
n_threads      = as.integer(max(1L, min(max_by_cores, max_by_memory)))
n_math_threads = as.integer(max(1L, min(n_threads, 8L)))

## apply thread settings -------------------------------------------------------

data.table::setDTthreads(threads = n_threads)
collapse::set_collapse(nthreads  = n_threads)
options(arrow.use_threads        = TRUE)
Sys.setenv(ARROW_NUM_THREADS     = n_threads)
options(mc.cores                 = n_threads)

message(
  sprintf("01 resources | OS=%s | Cores=%d | Threads=%d | MathThreads=%d | Avail RAM≈%s GB (rule: 2 GB/core)",
          os_type, all_cores, n_threads, n_math_threads,
          ifelse(is.finite(avail_ram_gb), round(avail_ram_gb, 1), "NA"))
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# cohort identification --------------------------------------------------------

## start by linking contiguous hospitalizations --------------------------------

### encounters with age >= 18 and dates within study window --------------------

hosp_blocks = 
  dplyr::filter(data_list$hospitalization, age_at_admission >= 18) |>
  dplyr::filter(admission_dttm >= start_date & admission_dttm <= end_date) |> 
  dplyr::filter(admission_dttm < discharge_dttm & !is.na(discharge_dttm)) |>
  dplyr::select(patient_id, hospitalization_id, admission_dttm, discharge_dttm) |>
  dplyr::collect() |>
  roworder(patient_id, admission_dttm)

### use data.table to find joined hospitalizations with <= 6h gaps -------------

link_hours = 6L
linked     = as.data.table(hosp_blocks)
setorder(linked, patient_id, admission_dttm)

#### calculate gaps between encounters and mark encounters that should be linked
linked[, next_admit := shift(admission_dttm, type = "lead"),           by    = patient_id]
linked[, next_gap   := as.numeric(difftime(next_admit, discharge_dttm, units = "hours"))]
linked[, prev_dc    := shift(discharge_dttm, type = "lag"),            by    = patient_id]  
linked[, prev_gap   := as.numeric(difftime(admission_dttm, prev_dc,    units = "hours"))]
linked[, link_flag  := (next_gap < link_hours | prev_gap < link_hours)]
linked[is.na(link_flag), link_flag := FALSE]

#### create unique joined hospitalization ID, new group whenever gap > link_hours from  previous discharge
linked[, new_group := is.na(prev_gap) | prev_gap >= link_hours]
linked[, joined_hosp_id := .GRP, by = .(patient_id, cumsum(new_group))]

#### create hid_jid_crosswalk --------------------------------------------------
hid_jid_crosswalk = select(linked, ends_with("id")) |> as_tidytable()

# inclusion --------------------------------------------------------------------

## all adult CAR-T infusions ---------------------------------------------------

cart = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(med_group == "car_t") |>
  dplyr::filter(admin_dttm >= start_date & admin_dttm <= end_date) |>
  dplyr::mutate(
    cart_name = tolower(med_category),
    cart_group = dplyr::case_when(
      cart_name == "axicabtagene_ciloleucel"      ~ 1L,
      cart_name == "brexucabtagene_autoleucel"    ~ 2L,
      cart_name == "ciltacabtagene_autoleucel"    ~ 3L,
      cart_name == "idecabtagene_vicleucel"       ~ 4L,
      cart_name == "lisocabtagene_maraleucel"     ~ 5L,
      cart_name == "tisagenlecleucel"             ~ 6L,
      TRUE                                        ~ NA_integer_
    )) |>
  dplyr::distinct(hospitalization_id, admin_dttm, cart_name, cart_group) |>
  dplyr::collect()

cart = join(hid_jid_crosswalk, cart, how = "inner") 

cart_pids = funique(cart$patient_id)
cart_hids = funique(cart$hospitalization_id)
cart_jids = funique(cart$joined_hosp_id)

cart =
  fselect(cart, joined_hosp_id, admin_dttm, cart_name, cart_group) |>
  funique()

## link to hospitalizations for these encounters -------------------------------

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
  join(enc, hid_jid_crosswalk, how = "inner", multiple = T) |>
  roworder(admission_dttm, discharge_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    patient_id         = ffirst(patient_id),
    admission_dttm     = ffirst(admission_dttm),
    discharge_dttm     = flast(discharge_dttm),
    discharge_category = tolower(flast(discharge_category))
  )

df = 
  join(enc, cart, how = "inner") |>
  roworder(admin_dttm) |>
  fgroup_by(patient_id) |>
  ffirst()

rm(cart, cart_vec, enc, hosp_blocks, linked); gc()

df |> janitor::get_dupes(joined_hosp_id)

# exclusion --------------------------------------------------------------------

## set up tracker --------------------------------------------------------------

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

## < 24h in-hospital post-cart (death/hospice vs early discharge) --------------

df$time_to_dc = as.numeric(difftime(df$discharge_dttm, df$admin_dttm), "hours")

short_stay    = fsubset(df, time_to_dc < 24)
n_short_death = fnrow(fsubset(short_stay, tolower(discharge_category) %in% c("expired", "hospice")))
n_short_other = fnrow(short_stay) - n_short_death

exclusions = rbind(
  exclusions,
  data.table(
    step        = "1a",
    reason      = "< 24h post-infusion: death or hospice",
    n_excluded  = n_short_death,
    n_remaining = n_start - n_short_death
  ),
  data.table(
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

rm(short_stay, n_short_death, n_short_other); gc()

cart_pats  = funique(df$patient_id)
all_c_hids = fsubset(hid_jid_crosswalk, patient_id %in% cart_pats) |> pull(hospitalization_id)
  
## at least 24h temperature data post-infusion ---------------------------------

### load temperature data ------------------------------------------------------

temp = 
  data_list$vitals |>
  dplyr::filter(vital_category == "temp_c") |>
  dplyr::filter(hospitalization_id %in% all_c_hids) |>
  dplyr::filter(recorded_dttm >= start_date) |>
  dplyr::filter(vital_value > 30 & vital_value < 45) |>
  dplyr::select(hospitalization_id, recorded_dttm, temp_c = vital_value) |>
  dplyr::collect()

temp = 
  join(temp, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fselect(-hospitalization_id) |>
  join(df, how = "inner") |>
  fsubset(recorded_dttm >= admission_dttm) |>
  fsubset(recorded_dttm <= discharge_dttm) |>
  fselect(joined_hosp_id, cart_t0 = admin_dttm, recorded_dttm, temp_c) |>
  distinct()

has_temp_24h = 
  fsubset(temp, recorded_dttm >= cart_t0 + lubridate::dhours(24)) |>
  fselect(joined_hosp_id) |>
  funique()

n_before_temp = fnrow(df)
no_temp_24h   = fsubset(df, !joined_hosp_id %in% has_temp_24h$joined_hosp_id)
no_temp_death = fnrow(fsubset(no_temp_24h, tolower(discharge_category) %in% c("expired", "hospice")))
no_temp_other = fnrow(no_temp_24h) - no_temp_death

df = fsubset(df, joined_hosp_id %in% has_temp_24h$joined_hosp_id)

### update tracker -------------------------------------------------------------

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

rm(no_temp_24h, no_temp_death, no_temp_other, n_before_temp, has_temp_24h); gc()

## not febrile in 12h prior to infusion ----------------------------------------

febrile_before = 
  fsubset(temp, temp_c >= 38) |>
  fsubset(recorded_dttm < cart_t0) |>
  fsubset(recorded_dttm >= cart_t0 - lubridate::dhours(12)) |>
  fselect(joined_hosp_id) |>
  funique()

n_before_febrile = fnrow(df)
n_febrile = sum(df$joined_hosp_id %in% febrile_before$joined_hosp_id)

df = fsubset(df, !joined_hosp_id %in% febrile_before$joined_hosp_id)

exclusions = rbind(
  exclusions,
  tidytable(
    step        = "2",
    reason      = "febrile within 12h pre-infusion",
    n_excluded  = n_febrile,
    n_remaining = fnrow(df)
  )
)

## finalize exclusion summary --------------------------------------------------

exclusions = rbind(
  data.table(
    step        = "0",
    reason      = "adult car-t infusions in date range",
    n_excluded  = 0L,
    n_remaining = n_start
  ),
  exclusions
)

exclusions[is.na(n_remaining), n_remaining := shift(n_remaining, type = "lead"), by = .I]

rm(n_start, febrile_before, n_febrile, n_before_febrile); gc()

print(exclusions)

## save exclusion summary for flow diagram -------------------------------------

excl_summary_path = file.path(
  project_location, "final_output",
  paste0("exclusions_", site_lowercase, "_", today, ".csv")
)

fwrite(exclusions, excl_summary_path)

message("Exclusion summary saved to: ", excl_summary_path)

# add patient details ----------------------------------------------------------

## patient table ---------------------------------------------------------------

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
  join(df, patient, how = "left", multiple = T) |>
  ftransform(
    race_category = case_when(
      tolower(race_category) %in% c("white", "asian")                ~ tolower(race_category),
      tolower(race_category) == "black or african american"          ~ tolower(race_category),
      TRUE                                                           ~ "other/unknown"
    ),
    ethnicity_category = case_when(
      tolower(ethnicity_category) %in% c("hispanic", "non-hispanic") ~ tolower(ethnicity_category),
      TRUE                                                           ~ "unknown"
    ),
    sex_category = case_when(
      tolower(sex_category) %in% c("male", "female")                 ~ tolower(sex_category),
      TRUE                                                           ~ "unknown"
    )
  )

## encounter table -------------------------------------------------------------

enc = 
  data_list$hospitalization |>
  dplyr::filter(patient_id %in% df$patient_id) |>
  dplyr::select(hospitalization_id, age = age_at_admission) |>
  dplyr::collect()

enc = 
  join(enc, hid_jid_crosswalk, how = "inner") |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(age = fmin(age))

df = join(df, enc, how = "left", multiple = T)

anchor_df = fselect(df, patient_id, joined_hosp_id, cart_t0 = admin_dttm)

rm(patient, exclusions, enc); gc()

# time of day, local -----------------------------------------------------------

df$admin_localtime   = with_tz(df$admin_dttm, tz = tz_local)
df$cart_hour_precise = hour(df$admin_localtime) + minute(df$admin_localtime)/60L
df$cart_hour_floor   = floor(df$cart_hour_precise)
df$cart_hour_round   = round(df$cart_hour_precise)

df = fselect(df, -admin_localtime)

# diagnoses --------------------------------------------------------------------

hids_of_interest = 
  fselect(df, joined_hosp_id) |>
  join(hid_jid_crosswalk, how = "inner", multiple = T) |>
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

## cancer diagnoses ------------------------------------------------------------

expand_icd = function(prefix) {
  paste0(prefix, 0:9)
}

expand_icd_range = function(base, subtypes) {
  unlist(lapply(paste0(base, ".", subtypes), expand_icd))
}

dx_groups = list(
  myeloma      = expand_icd_range("C90", c("0", "1", "2", "3")),
  dlbcl        = expand_icd("C83.3"),
  follicular   = expand_icd_range("C82", c("0", "1", "2", "3", "4", "5", "6", "8", "9")),
  mantle_cell  = expand_icd("C83.1"),
  marginal     = expand_icd("C83.0"),
  burkitt      = expand_icd("C83.7"),
  other_nhl    = expand_icd("C83.5"),
  nhl_nos      = expand_icd("C83.9"),
  lymphoma_nos = expand_icd_range("C85", c("1", "2", "8", "9")),
  cll_sll      = expand_icd("C91.1"),
  b_all        = c("C91.00", "C91.01", "C91.02"),
  hct          = c("Z94.81", "Z94.84", paste0("T86.0", 0:3), "T86.09", "T86.5")
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

hct_vec = 
  fsubset(dx_cancer, cancer_cat == "hct") |> 
  pull(patient_id)

lymphoma_priority = tidytable(
  cancer_cat = c("dlbcl", "mantle_cell", "follicular", "marginal", 
                 "burkitt", "cll_sll", "other_nhl", "nhl_nos", "lymphoma_nos"),
  lym_rank   = 1L:9L
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

df$cancer_cat = case_when(
  !is.na(df$cancer_cat)                        ~ df$cancer_cat,
  df$cancer_cat == "ciltacabtagene_autoleucel" ~ "myeloma",
  df$cancer_cat == "idecabtagene_vicleucel"    ~ "myelomaa",
  TRUE                                         ~ "unknown"
)

rm(dx_matched, dx_primary, dx_with_cart, cart_indication)
rm(hct_vec, dx_cancer, dx_lookup, dx_groups, lymphoma_priority); gc()

## elixhauser comorbidities ----------------------------------------------------

elix = 
  comorbidity::comorbidity(
    dx, 
    id      = "patient_id", 
    code    = "diagnosis_code", 
    map     = "elixhauser_icd10_quan", 
    assign0 = T
  )

### calculate vw scores --------------------------------------------------------

vw   = comorbidity::score(elix, weights = "vw", assign0 = T)
elix = cbind(elix, vw = vw) |> fselect(patient_id, vw)
df   = join(df, elix, how = "left", multiple = F)

df$miss_vw = as.integer(is.na(df$vw))
df$vw      = if_else(is.na(df$vw), 0L, df$vw)

rm(elix, vw, dx); gc()

# outcomes ---------------------------------------------------------------------

## death -----------------------------------------------------------------------

df$death_dttm = case_when(
  !is.na(df$death_dttm)                       ~ df$death_dttm, 
  tolower(df$discharge_category) == "expired" ~ df$discharge_dttm,
  TRUE                                        ~ df$death_dttm
)

df$discharge_category = case_when(
  df$death_dttm <= df$discharge_dttm          ~ "expired",
  TRUE                                        ~ tolower(df$discharge_dttm)
)

df$dc_dead_01    = if_else(df$discharge_category == "expired", 1L, 0L)
df$dc_hospice_01 = if_else(df$discharge_category == "hospice", 1L, 0L)
df$time_to_death = as.numeric(difftime(df$death_dttm, df$admin_dttm, units = "days"))
df$dead_30_01    = as.integer(!is.na(df$time_to_death) & df$time_to_death <= 30)
df$dead_90_01    = as.integer(!is.na(df$time_to_death) & df$time_to_death <= 90)
df$dead_180_01   = as.integer(!is.na(df$time_to_death) & df$time_to_death <= 180)
df$dead_365_01   = as.integer(!is.na(df$time_to_death) & df$time_to_death <= 365)

## length of stay --------------------------------------------------------------

df$los_hospital_d = as.numeric(difftime(df$discharge_dttm, df$admission_dttm, units = "days"))

## icu outcomes ----------------------------------------------------------------

### icu admission within 30 days -----------------------------------------------

icu_flag = 
  data_list$adt |>
  dplyr::filter(tolower(location_category) == "icu") |>
  dplyr::filter(in_dttm >= start_date) |>
  dplyr::distinct(hospitalization_id, in_dttm) |>
  dplyr::collect() 

icu_flag = 
  join(icu_flag, hid_jid_crosswalk, how = "inner") |>
  fselect(patient_id, icu_dttm = in_dttm)

icu_flag = 
  join(icu_flag, anchor_df, how = "inner", multiple = T) |>
  fsubset(icu_dttm >= cart_t0 & icu_dttm <= cart_t0 + ddays(30)) |>
  pull(patient_id) |>
  funique()

### vasopressors within 30 days ------------------------------------------------

vaso_flag = 
  data_list$medication_admin_continuous |>
  dplyr::filter(tolower(med_group) == "vasoactives") |> 
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm) |>
  dplyr::collect() 

vaso_flag = 
  join(vaso_flag, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fsubset(patient_id %in% df$patient_id) |>
  fselect(patient_id, va_dttm = admin_dttm)

vaso_flag = 
  join(vaso_flag, anchor_df, how = "inner", multiple = T) |>
  fsubset(va_dttm >= cart_t0 & va_dttm <= cart_t0 + ddays(30)) |>
  pull(patient_id) |>
  funique()

### respiratory support within 30 days -----------------------------------------

resp_flag = 
  data_list$respiratory_support |>
  dplyr::filter(tolower(device_category) %in% c("imv", "nippv", "high flow nc")) |>
  dplyr::filter(recorded_dttm >= start_date) |>
  dplyr::select(hospitalization_id, recorded_dttm, device_category) |>
  dplyr::collect()

resp_flag = 
  join(resp_flag, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fsubset(patient_id %in% df$patient_id) |>
  fselect(patient_id, recorded_dttm, device_category)

imv_flag = 
  join(resp_flag, anchor_df, how = "inner", multiple = T) |>
  fsubset(recorded_dttm >= cart_t0 & recorded_dttm <= cart_t0 + ddays(30)) |>
  fsubset(tolower(device_category) == "imv") |>
  pull(patient_id) |>
  funique()

resp_flag = 
  join(resp_flag, anchor_df, how = "inner", multiple = T) |>
  fsubset(recorded_dttm >= cart_t0 & recorded_dttm <= cart_t0 + ddays(30)) |>
  pull(patient_id) |>
  funique()

### add icu flags to cohort data frame -----------------------------------------

df$resp_30_01 = as.integer(df$patient_id %in% resp_flag)
df$imv_30_01  = as.integer(df$patient_id %in% imv_flag)
df$va_30_01   = as.integer(df$patient_id %in% vaso_flag)

df$icu_30_01 = case_when(
  df$patient_id %in% icu_flag ~ 1L,
  df$resp_30_01 == 1L         ~ 1L,
  df$va_30_01   == 1L         ~ 1L,
  TRUE                        ~ 0L
)

df = fselect(df, -discharge_category)

rm(vaso_flag, resp_flag, imv_flag, icu_flag); gc()

# medication-based outcomes ----------------------------------------------------

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
  join(toci, anchor_df, how = "inner", multiple = T) |>
  fsubset(toci_dttm >= cart_t0 & toci_dttm <= cart_t0 + lubridate::ddays(30)) |>
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
  dplyr::filter(med_category == "anakinra") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm) |>
  dplyr::collect() |>
  join(hid_jid_crosswalk, how = "inner", multiple = FALSE) |>
  fselect(patient_id, ana_dttm = admin_dttm)

anakinra = 
  join(anakinra, anchor_df, how = "inner", multiple = T) |>
  fsubset(ana_dttm >= cart_t0 & ana_dttm <= cart_t0 + lubridate::ddays(30)) |>
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

## apap ------------------------------------------------------------------------

apap = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(med_category == "acetaminophen") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm, med_dose) |>
  dplyr::collect() |>
  join(hid_jid_crosswalk, how = "inner", multiple = FALSE) |>
  fselect(patient_id, apap_dttm = admin_dttm, med_dose)

apap = 
  join(apap, anchor_df, how = "inner", multiple = T) |>
  fsubset(apap_dttm >= cart_t0 & apap_dttm <= cart_t0 + lubridate::ddays(3)) |>
  funique() |>
  roworder(apap_dttm) |>
  fmutate(h_since = as.numeric(difftime(apap_dttm, cart_t0), "hours")) |>
  fmutate(d_since = as.integer(floor(h_since/24))) |>
  fgroup_by(patient_id, d_since) |>
  fsummarize(apap_cum = fsum(med_dose)) |>
  pivot_wider(
    names_from   = d_since,
    values_from  = apap_cum,
    names_prefix = "apap_day_"
  ) 

## steroids --------------------------------------------------------------------

### prep conversion data frame -------------------------------------------------

steroid_map = tidytable(
  med_category = c(
    "dexamethasone", 
    "betamethasone", 
    "methylprednisolone", 
    "prednisone", 
    "prednisolone", 
    "hydrocortisone"
  ),
  dex_factor = c(1.0, 1.25, 0.1875, 0.15, 0.15, 0.0375)
)

### prep patient weights -------------------------------------------------------

weight = 
  data_list$vitals |>
  dplyr::filter(vital_category == "weight_kg") |>
  dplyr::filter(recorded_dttm >= start_date) |>
  dplyr::select(hospitalization_id, recorded_dttm, weight_kg = vital_value) |>
  dplyr::collect() 

weight = 
  join(weight, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fselect(patient_id, recorded_dttm, weight_kg)

weight = 
  join(weight, anchor_df, how = "inner", multiple = T) |>
  fsubset(recorded_dttm >= cart_t0 - ddays(10)) |>
  fsubset(recorded_dttm <= cart_t0 + ddays(30)) |>
  roworder(recorded_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(weight_kg = ffirst(weight_kg))

weight = 
  join(weight, df, how = "right", multiple = F) |>
  fmutate(weight_kg = if_else(is.na(weight_kg), fmedian(weight_kg), weight_kg)) |>
  fselect(patient_id, weight_kg) 

### steroid table --------------------------------------------------------------

steroid_admin = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(!is.na(med_dose) & med_group == "steroid") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(
    hospitalization_id, 
    admin_dttm, 
    med_category, 
    med_dose, 
    med_dose_unit
  ) |>
  dplyr::collect() 

steroid_admin = 
  join(steroid_admin, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fselect(-hospitalization_id, -joined_hosp_id) |>
  funique()

steroid_admin = 
  join(steroid_admin, anchor_df, how = "inner", multiple = T) |>
  fsubset(admin_dttm >= cart_t0) |>
  fsubset(admin_dttm <= cart_t0 + ddays(30))

### convert steroids to dexamethasone equivalents per day ----------------------

steroid_admin = 
  join(steroid_admin, weight, how = "left", multiple = F) |>
  fmutate(u = if_else(med_dose_unit == "mg/kg/day", "mg/kg", med_dose_unit)) |>
  fmutate(dose_mg = if_else(u == "mg", med_dose, med_dose*weight_kg)) |>
  fsubset(!is.na(dose_mg)) |>
  ftransform(h_since = as.numeric(difftime(admin_dttm, cart_t0), "hours")) |>
  ftransform(d_since = as.integer(floor(h_since/24))) |>
  fselect(patient_id, d_since, cart_t0, admin_dttm, med_category, dose_mg) |>
  funique()

steroid_daily = 
  join(steroid_admin, steroid_map, how = "inner") |>
  ftransform(dex_mg_equivalent = dose_mg*dex_factor) |>
  fgroup_by(patient_id, d_since) |>
  fsummarise(dex_mg_equivalent = fsum(dex_mg_equivalent)) |>
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
  join(steroid_daily, how = "full", multiple = F)

## add meds to data frame ------------------------------------------------------

df = join(df, toci,     how = "left", multiple = F)
df = join(df, anakinra, how = "left", multiple = F)
df = join(df, apap,     how = "left", multiple = F)
df = join(df, steroid,  how = "left", multiple = F)

df$toci_n_doses_30d     = replace_na(df$toci_n_doses_30d,     0L)
df$anakinra_n_doses_30d = replace_na(df$anakinra_n_doses_30d, 0L)
df$steroid_days_30d     = replace_na(df$steroid_days_30d,     0L)
df$apap_day_0           = replace_na(df$apap_day_0,           0L)
df$apap_day_1           = replace_na(df$apap_day_1,           0L)
df$apap_day_2           = replace_na(df$apap_day_2,           0L)
df$steroid_any_72h      = replace_na(df$steroid_any_72h,      0L)
df$steroid_days_30d     = replace_na(df$steroid_days_30d,     0L)
df$h_until_steroids     = as.numeric(difftime(df$steroid_first, df$admin_dttm), "hours")
  
df = fselect(df, -steroid_first)

rm(toci, anakinra, apap, weight, steroid, steroid_admin, steroid_daily, steroid_map); gc()

## antibiotics and cultures ----------------------------------------------------

### antibiotics ----------------------------------------------------------------

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
  fsubset(admin_dttm >= cart_t0 & admin_dttm <= cart_t0 + lubridate::ddays(30)) |>
  roworder(admin_dttm) |>
  fgroup_by(patient_id, med_category) |>
  fsummarize(admin_dttm = ffirst(admin_dttm))

df$abx_30_01 = if_else(df$patient_id %in% abx$patient_id, 1L, 0L)

### blood cultures -------------------------------------------------------------

blood_cx = 
  data_list$microbiology_culture |>
  dplyr::filter(fluid_category == "blood_buffy") |>
  dplyr::filter(organism_category != "no_growth") |>
  dplyr::select(hospitalization_id, result_dttm, organism_category) |>
  dplyr::collect() 

blood_cx = 
  join(blood_cx, hid_jid_crosswalk, how = "inner") |>
  fselect(joined_hosp_id, result_dttm, organism_category) |>
  funique()

blood_cx = 
  join(blood_cx, anchor_df, how = "inner") |>
  fsubset(result_dttm >= cart_t0 & result_dttm <= cart_t0 + lubridate::ddays(30)) |>
  roworder(result_dttm) |>
  fgroup_by(patient_id, organism_category) |>
  fsummarize(result_dttm = ffirst(result_dttm))
  
 df$bsi_30_01 = if_else(df$patient_id %in% blood_cx$patient_id, 1L, 0L)
   
 ### save blood culture table for traj link - need to know what bugs for which traj

rm(blood_cx, abx_list, abx); gc()

# save -------------------------------------------------------------------------


