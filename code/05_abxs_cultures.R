# ==============================================================================
# 05_abxs_cultures.R
# Infection outcome variables for CAR-T project
# Requires: df, all_hid_crosswalk, anchor_df from prior scripts
# Outputs: outcomes added to df, blood culture table by trajectory
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("all_hid_crosswalk") || !exists("anchor_df") || !exists("save_output")) {
  stop("Run 00_setup.R through 04_outcomes.R first.", call. = FALSE)
}
# Antibiotic exposure and infections ===========================================
## Organism classification -----------------------------------------------------

viral_organisms = c(
  "adenovirus", "cytomegalovirus", "enterovirus", "epstein_barr_virus", "hepatitis_a",
  "hepatitis_b", "hepatitis_c", "herpes_simplex", "herpes_zoster", "hhv_6", "hiv_htlv",
  "influenza", "measles", "mumps", "papovavirus", "parainfluenza", "polyomavirus",
  "respiratory_syncytial_virus", "rhinovirus", "rotavirus", "rubella", "sars_cov2",
  "viral_other")

fungal_organisms = c(
  "acremonium_sp", "alternaria_alternata", "alternaria_sp", "arthrographis_sp",
  "aspergillus_flavus", "aspergillus_fumigatus", "aspergillus_nidulans",
  "aspergillus_niger", "aspergillus_ochraceus", "aspergillus_sp",
  "aspergillus_sydowii", "aspergillus_terreus", "aspergillus_versicolor",
  "aureobasidium_sp", "blastomyces_dermatitidis", "blastomyces_sp",
  "candida_albicans", "candida_auris", "candida_colliculosa", "candida_dubliniensis",
  "candida_famata", "candida_fermentati", "candida_glabrata", "candida_guilliermondii",
  "candida_haemulonii", "candida_intermedia", "candida_kefyr", "candida_krusei",
  "candida_lambica", "candida_lipolytica", "candida_lusitaniae", "candida_nivariensis",
  "candida_orthopsilosis", "candida_parapsilosis", "candida_pararugosa", "candida_pelliculosa",
  "candida_pulcherrima", "candida_rugosa", "candida_sp", "candida_tropicalis",
  "candida_utilis", "chaetomium_sp", "cladosporium_sp", "coccidioides_immitis",
  "coccidioides_posadasii", "coccidioides_sp", "cokeromyces_recurvatus",
  "cryptococcus_albidus", "cryptococcus_neoformans", "cunninghamella_sp",
  "curvularia_sp", "dematiaceous_sp", "exophiala_dermatitidis", "exophiala_sp",
  "fonsecaea_sp", "fusarium_oxysporum", "fusarium_solani", "fusarium_sp",
  "geotrichum_sp", "histoplasma_capsulatum", "histoplasma_sp", "lichtheimia_sp",
  "lomentospora_prolificans", "malassezia_furfur", "malassezia_sp", "metarhizium_sp",
  "microsporum_canis", "microsporum_sp", "mucor_sp", "paecilomyces_sp",
  "penicillium_sp", "phialophora_sp", "piedraia_hortae", "pneumocystis_jirovecii",
  "purpureocillium_lilacinum", "ramichloridium_sp", "rhinocladiella_sp",
  "rhizomucor_sp", "rhizopus_arrhizus", "rhizopus_microsporus", "rhizopus_sp",
  "rhodotorula_mucilaginosa", "rhodotorula_sp", "saccharomyces_cerevisiae",
  "saccharomyces_sp", "sarocladium_sp", "scedosporium_apiospermum",
  "scedosporium_sp", "scopulariopsis_sp", "sporothrix_schenckii",  "sporothrix_sp",
  "syncephalastrum_sp", "talaromyces_marneffei", "talaromyces_sp", "trichophyton_mentagrophytes",
  "trichophyton_rubrum", "trichophyton_sp", "trichophyton_tonsurans", "trichosporon_asahii",
  "trichosporon_sp", "yeast")

exclude_results = c("not_detected", "indeterminate")

## Load culture data -----------------------------------------------------------

cultures = 
  data_list$microbiology_culture |>
  dplyr::filter(!is.na(collect_dttm)) |>
  dplyr::filter(!is.na(organism_group)) |>
  dplyr::filter(organism_group != "no_growth") |>
  dplyr::select(hospitalization_id, organism_id, collect_dttm, result_dttm, fluid_category, organism_category) |>
  dplyr::collect()

cultures = 
  join(cultures, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, collect_dttm, result_dttm, fluid_category, organism_category, organism_id) |>
  funique()

cultures =
  join(cultures, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fmutate(
    h_since_cart = as.numeric(difftime(collect_dttm, cart_t0, units = "hours")),
    d_since_cart = as.integer(floor(h_since_cart / 24)))

## Load PCR data ---------------------------------------------------------------

pcr = 
  data_list$microbiology_nonculture |>
  dplyr::filter(!is.na(collect_dttm)) |>
  dplyr::filter(!is.na(organism_group)) |>
  dplyr::filter(organism_group != "no_growth") |>
  dplyr::filter(tolower(result_category) %!in% exclude_results) |>
  dplyr::select(hospitalization_id, collect_dttm, result_dttm, fluid_category, organism_group, organism_category) |>
  dplyr::collect()

pcr = 
  join(pcr, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, collect_dttm, result_dttm, fluid_category, organism_category) |>
  funique()

pcr =
  join(pcr, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fmutate(
    h_since_cart = as.numeric(difftime(collect_dttm, cart_t0, units = "hours")),
    d_since_cart = as.integer(floor(h_since_cart / 24)))

infection <- dplyr::bind_rows(cultures, pcr)

rm(cultures, pcr)
gc()
    
##Classify organism type and infection location

infection <- infection |>
  dplyr::mutate(
    organism_type = dplyr::case_when(
      tolower(organism_category) %in% viral_organisms  ~ "viral",
      tolower(organism_category) %in% fungal_organisms ~ "fungal",
      TRUE                                             ~ "bacterial"),
    infection_source = dplyr::case_when(
      tolower(fluid_category) %in% c("blood", "blood_buffy", "catheter_tip")                         ~ "blood",
      tolower(fluid_category) %in% c("genito_urinary_tract", "kidneys_renal_pelvis_ureters_bladder") ~ "urine",
      tolower(fluid_category) %in% c("resp", "respiratory_tract", "respiratory_tract_lower",
                                     "nasopharynx_upperairway", "oropharynx_tongue_oralcavity")      ~ "respiratory",
      tolower(fluid_category) %in% c("feces_stool", "gastrointestinal_tract", "stomach",
                                     "small_intestine", "large_intestine", "woundsite", 
                                     "skin_unspecified", "skin_rash_pustules_abscesses", 
                                     "skin_disseminated_multiple_sites", "central_nervous_system", 
                                     "brain", "meninges_csf", "spinal_cord","other_unspecified")        ~ "other",
                                     TRUE                                                               ~ "other"))
                                                            
#Deduplicate infections
infection <- infection |> dplyr::distinct(patient_id, d_since_cart, fluid_category, organism_category, .keep_all = TRUE)

#Create binary flags and join back to df
infection_base_flags =
  infection |>
  fsubset(collect_dttm >= cart_t0 - ddays(1) & collect_dttm <= cart_t0) |>
  fmutate(
    bacterial_base_01 = as.integer(organism_type == "bacterial"),
    viral_base_01     = as.integer(organism_type == "viral"),
    fungal_base_01    = as.integer(organism_type ==  "fungal"),
    infection_base_01 = as.integer(!is.na(organism_type)),
    blood_base_01     = as.integer(infection_source == "blood"),
    urine_base_01     = as.integer(infection_source == "urine"),
    resp_base_01      = as.integer(infection_source == "respiratory"),
    other_inf_base_01 = as.integer(infection_source == "other")
  ) |>
  fselect(patient_id, bacterial_base_01, viral_base_01, fungal_base_01, 
         infection_base_01, blood_base_01, urine_base_01, resp_base_01, 
         other_inf_base_01)

infection_72h_flags =
  infection |>
  fsubset(collect_dttm >= cart_t0 & collect_dttm <= cart_t0 + ddays(2)) |>
  fmutate(
    bacterial_72h_01 = as.integer(organism_type    == "bacterial"),
    viral_72h_01     = as.integer(organism_type    == "viral"),
    fungal_72h_01    = as.integer(organism_type    == "fungal"),
    infection_72h_01 = as.integer(!is.na(organism_type)),
    blood_72h_01     = as.integer(infection_source == "blood"),
    urine_72h_01     = as.integer(infection_source == "urine"),
    resp_72h_01      = as.integer(infection_source == "respiratory"),
    other_inf_72h_01 = as.integer(infection_source == "other")
  ) |>
  fgroup_by(patient_id) |>
  fsummarize(
    bacterial_72h_01 = max(bacterial_72h_01, na.rm = TRUE),
    viral_72h_01     = max(viral_72h_01, na.rm = TRUE),
    fungal_72h_01    = max(fungal_72h_01, na.rm = TRUE),
    infection_72h_01 = max(infection_72h_01, na.rm = TRUE),
    blood_72h_01     = max(blood_72h_01, na.rm = TRUE),
    urine_72h_01     = max(urine_72h_01, na.rm = TRUE),
    resp_72h_01      = max(resp_72h_01, na.rm = TRUE),
    other_inf_72h_01 = max(other_inf_72h_01, na.rm = TRUE)
  ) |>
  fselect(patient_id, bacterial_72h_01, viral_72h_01, fungal_72h_01, 
          infection_72h_01, blood_72h_01, urine_72h_01, resp_72h_01, 
          other_inf_72h_01)

infection_3_30d_flags =
  infection |>
  fsubset(collect_dttm >= cart_t0 + ddays(3) & collect_dttm <= cart_t0 + ddays(30)) |>
  fmutate(
    bacterial_3_30d_01 = as.integer(organism_type    == "bacterial"),
    viral_3_30d_01     = as.integer(organism_type    == "viral"),
    fungal_3_30d_01    = as.integer(organism_type    == "fungal"),
    infection_3_30d_01 = as.integer(!is.na(organism_type)),
    blood_3_30d_01     = as.integer(infection_source == "blood"),
    urine_3_30d_01     = as.integer(infection_source == "urine"),
    resp_3_30d_01      = as.integer(infection_source == "respiratory"),
    other_inf_3_30d_01 = as.integer(infection_source == "other")
  ) |>
  fgroup_by(patient_id) |>
  fsummarize(
    bacterial_3_30d_01 = max(bacterial_3_30d_01, na.rm = TRUE),
    viral_3_30d_01     = max(viral_3_30d_01, na.rm = TRUE),
    fungal_3_30d_01    = max(fungal_3_30d_01, na.rm = TRUE),
    infection_3_30d_01 = max(infection_3_30d_01, na.rm = TRUE),
    blood_3_30d_01     = max(blood_3_30d_01, na.rm = TRUE),
    urine_3_30d_01     = max(urine_3_30d_01, na.rm = TRUE),
    resp_3_30d_01      = max(resp_3_30d_01, na.rm = TRUE),
    other_inf_3_30d_01 = max(other_inf_3_30d_01, na.rm = TRUE)
  ) |>
  fselect(patient_id, bacterial_3_30d_01, viral_3_30d_01, fungal_3_30d_01, 
          infection_3_30d_01, blood_3_30d_01, urine_3_30d_01, resp_3_30d_01, 
          other_inf_3_30d_01)


df =
  join(df, infection_base_flags,    on = "patient_id", how = "left", multiple = FALSE) |>
  join(infection_72h_flags,         on = "patient_id", how = "left", multiple = FALSE) |>
  join(infection_3_30d_flags,         on = "patient_id", how = "left", multiple = FALSE)

df$infection_base_01        = replace_na(df$infection_base_01, 0L)
df$bacterial_base_01        = replace_na(df$bacterial_base_01, 0L)
df$viral_base_01            = replace_na(df$viral_base_01, 0L)
df$fungal_base_01           = replace_na(df$fungal_base_01, 0L)
df$blood_base_01            = replace_na(df$blood_base_01, 0L)
df$urine_base_01            = replace_na(df$urine_base_01, 0L)
df$resp_base_01             = replace_na(df$resp_base_01, 0L)
df$other_inf_base_01        = replace_na(df$other_inf_base_01, 0L)

df$infection_72h_01         = replace_na(df$infection_72h_01, 0L)
df$bacterial_72h_01         = replace_na(df$bacterial_72h_01, 0L)
df$viral_72h_01             = replace_na(df$viral_72h_01, 0L)
df$fungal_72h_01            = replace_na(df$fungal_72h_01, 0L)
df$blood_72h_01             = replace_na(df$blood_72h_01, 0L)
df$urine_72h_01             = replace_na(df$urine_72h_01, 0L)
df$resp_72h_01              = replace_na(df$resp_72h_01, 0L)
df$other_inf_72h_01         = replace_na(df$other_inf_72h_01, 0L)

df$infection_3_30d_01         = replace_na(df$infection_3_30d_01, 0L)
df$bacterial_3_30d_01         = replace_na(df$bacterial_3_30d_01, 0L)
df$viral_3_30d_01             = replace_na(df$viral_3_30d_01, 0L)
df$fungal_3_30d_01            = replace_na(df$fungal_3_30d_01, 0L)
df$blood_3_30d_01             = replace_na(df$blood_3_30d_01, 0L)
df$urine_3_30d_01             = replace_na(df$urine_3_30d_01, 0L)
df$resp_3_30d_01              = replace_na(df$resp_3_30d_01, 0L)
df$other_inf_3_30d_01         = replace_na(df$other_inf_3_30d_01, 0L)


rm(infection_base_flags, infection_72h_flags, infection_3_30d_flags)
gc()

## Infection summaries by trajectory =====================================================================
#Overall prevalence of infection per trajectory by type and source

infections_base_prevalence =
  df |>
  fgroup_by(trajectory) |>
  fsummarize(
    n_patients    = fnobs(patient_id),
    n_infected    = fsum(infection_base_01, na.rm = TRUE),
    n_bacterial   = fsum(bacterial_base_01, na.rm = TRUE),
    n_viral       = fsum(viral_base_01, na.rm = TRUE),
    n_fungal      = fsum(fungal_base_01, na.rm = TRUE),
    n_blood       = fsum(blood_base_01, na.rm = TRUE),
    n_urine       = fsum(urine_base_01, na.rm = TRUE),
    n_resp        = fsum(resp_base_01, na.rm = TRUE),
    n_other_source   = fsum(other_inf_base_01, na.rm = TRUE)
  ) |>
  fmutate(
    pct_infected  = fifelse(n_infected > 0, round(100 * n_infected / n_patients, 1), NA_real_),
    pct_bacterial = fifelse(n_infected > 0, round(100 * n_bacterial / n_infected, 1), NA_real_),
    pct_viral     = fifelse(n_infected > 0, round(100 * n_viral / n_infected, 1), NA_real_),
    pct_fungal    = fifelse(n_infected > 0, round(100 * n_fungal / n_infected, 1), NA_real_),
    pct_blood     = fifelse(n_infected > 0, round(100 * n_blood / n_infected, 1), NA_real_),
    pct_urine     = fifelse(n_infected > 0, round(100 * n_urine / n_infected, 1), NA_real_),
    pct_resp     = fifelse(n_infected > 0, round(100 * n_resp / n_infected, 1), NA_real_),
    pct_other_source = fifelse(n_infected > 0, round(100 * n_other_source / n_infected, 1), NA_real_)
  ) |>
  fungroup() |>
  roworder(trajectory)

infections_72h_prevalence =
  df |>
  fgroup_by(trajectory) |>
  fsummarize(
    n_patients    = fnobs(patient_id),
    n_infected    = fsum(infection_72h_01, na.rm = TRUE),
    n_bacterial   = fsum(bacterial_72h_01, na.rm = TRUE),
    n_viral       = fsum(viral_72h_01, na.rm = TRUE),
    n_fungal      = fsum(fungal_72h_01, na.rm = TRUE),
    n_blood       = fsum(blood_72h_01, na.rm = TRUE),
    n_urine       = fsum(urine_72h_01, na.rm = TRUE),
    n_resp        = fsum(resp_72h_01, na.rm = TRUE),
    n_other_source   = fsum(other_inf_72h_01, na.rm = TRUE)
  ) |>
  fmutate(
    pct_infected  = fifelse(n_infected > 0, round(100 * n_infected / n_patients, 1), NA_real_),
    pct_bacterial = fifelse(n_infected > 0, round(100 * n_bacterial / n_infected, 1), NA_real_),
    pct_viral     = fifelse(n_infected > 0, round(100 * n_viral / n_infected, 1), NA_real_),
    pct_fungal    = fifelse(n_infected > 0, round(100 * n_fungal / n_infected, 1), NA_real_),
    pct_blood     = fifelse(n_infected > 0, round(100 * n_blood / n_infected, 1), NA_real_),
    pct_urine     = fifelse(n_infected > 0, round(100 * n_urine / n_infected, 1), NA_real_),
    pct_resp     = fifelse(n_infected > 0, round(100 * n_resp / n_infected, 1), NA_real_),
    pct_other_source = fifelse(n_infected > 0, round(100 * n_other_source / n_infected, 1), NA_real_)
  ) |>
  fungroup() |>
  roworder(trajectory)

infections_3_30d_prevalence =
  df |>
  fgroup_by(trajectory) |>
  fsummarize(
    n_patients    = fnobs(patient_id),
    n_infected    = fsum(infection_3_30d_01, na.rm = TRUE),
    n_bacterial   = fsum(bacterial_3_30d_01, na.rm = TRUE),
    n_viral       = fsum(viral_3_30d_01, na.rm = TRUE),
    n_fungal      = fsum(fungal_3_30d_01, na.rm = TRUE),
    n_blood       = fsum(blood_3_30d_01, na.rm = TRUE),
    n_urine       = fsum(urine_3_30d_01, na.rm = TRUE),
    n_resp        = fsum(resp_3_30d_01, na.rm = TRUE),
    n_other_source   = fsum(other_inf_3_30d_01, na.rm = TRUE)
  ) |>
  fmutate(
    pct_infected  = fifelse(n_infected > 0, round(100 * n_infected / n_patients, 1), NA_real_),
    pct_bacterial = fifelse(n_infected > 0, round(100 * n_bacterial / n_infected, 1), NA_real_),
    pct_viral     = fifelse(n_infected > 0, round(100 * n_viral / n_infected, 1), NA_real_),
    pct_fungal    = fifelse(n_infected > 0, round(100 * n_fungal / n_infected, 1), NA_real_),
    pct_blood     = fifelse(n_infected > 0, round(100 * n_blood / n_infected, 1), NA_real_),
    pct_urine     = fifelse(n_infected > 0, round(100 * n_urine / n_infected, 1), NA_real_),
    pct_resp     = fifelse(n_infected > 0, round(100 * n_resp / n_infected, 1), NA_real_),
    pct_other_source = fifelse(n_infected > 0, round(100 * n_other_source / n_infected, 1), NA_real_)
  ) |>
  fungroup() |>
  roworder(trajectory)


# Save summaries
save_output(infections_base_prevalence, "infections_72h_prevalence_by_trajectory", shareable = TRUE, subfolder = "trajectory")
save_output(infections_72h_prevalence, "infections_72h_time_by_trajectory", shareable = TRUE, subfolder = "trajectory")
save_output(infections_3_30d_prevalence, "infections_72h_time_by_trajectory", shareable = TRUE, subfolder = "trajectory")

rm(infections_base_prevalence, infections_72h_prevalence, infections_3_30d_prevalence)
gc()

## Daily infection summaries by trajectory -----------------------------------
#Overall prevalence of infection at 3-30 days by type and source
infection_bacterial_daily <-
  infection |>
  join(df |> fselect(patient_id, trajectory), on = "patient_id", how = "inner", multiple = FALSE) |>
  fsubset(collect_dttm >= cart_t0 & collect_dttm <= cart_t0 + ddays(30)) |>
  dplyr::filter(organism_type == "bacterial") |>
  dplyr::group_by(trajectory, d_since_cart) |>
  dplyr::summarise(
    n_patients_bacterial_infection = dplyr::n_distinct(patient_id),
    .groups = "drop"
  ) |>
  dplyr::arrange(trajectory, d_since_cart)


infection_viral_daily <-
  infection |>
  join(df |> fselect(patient_id, trajectory), on = "patient_id", how = "inner", multiple = FALSE) |>
  fsubset(collect_dttm >= cart_t0 & collect_dttm <= cart_t0 + ddays(30)) |>
  dplyr::filter(organism_type == "viral") |>
  dplyr::group_by(trajectory, d_since_cart) |>
  dplyr::summarise(
    n_patients_viral_infection = dplyr::n_distinct(patient_id),
    .groups = "drop"
  ) |>
  dplyr::arrange(trajectory, d_since_cart)

infection_fungal_daily <-
  infection |>
  join(df |> fselect(patient_id, trajectory), on = "patient_id", how = "inner", multiple = FALSE) |>
  fsubset(collect_dttm >= cart_t0 & collect_dttm <= cart_t0 + ddays(30)) |>
  dplyr::filter(organism_type == "fungal") |>
  dplyr::group_by(trajectory, d_since_cart) |>
  dplyr::summarise(
    n_patients_fungal_infection = dplyr::n_distinct(patient_id),
    .groups = "drop"
  ) |>
  dplyr::arrange(trajectory, d_since_cart)

infection_daily <-
  join(infection_bacterial_daily, infection_viral_daily, on = c("trajectory", "d_since_cart"),
       how = "left", multiple = FALSE) |>
  join(infection_fungal_daily, on = c("trajectory", "d_since_cart"), how = "left", multiple = FALSE)

save_output(infection_daily, "infection_daily_by_trajectory", shareable = TRUE, subfolder = "trajectory")


rm(infection, infection_bacterial_daily, infection_viral_daily, infection_fungal_daily, fungal_organisms, viral_organisms)
gc()

##Antimicrobial exposure
## broad-spectrum antibacterials -----------------------------------------------

broad_spectrum_abx = c(
  "amikacin", "aztreonam", "cefepime", "cefiderocol", "ceftaroline", 
  "ceftazidime", "ceftazidime_avibactam", "ceftazidime_clavulanate", 
  "ceftolozane_tazobactam", "colistin", "dalbavancin", "daptomycin", 
  "doripenem", "eravacycline", "ertapenem", "fidaxomicin", "gentamicin", 
  "imipenem", "imipenem_relebactam", "linezolid", "meropenem", 
  "meropenem_vaborbactam", "metronidazole", "piperacillin_tazobactam", 
  "tedizolid", "telavancin", "tigecycline", "tobramycin", "vancomycin"
)

##anaerobic antibacterials -----------------------------------------------------

anaerobic_abx <- c(
  "metronidazole", "piperacillin_tazobactam", "tigecycline", "doripenem",
  "eravacycline", "ertapenem", "imipenem", "imipenem_relebactam", "meropenem",
  "meropenem_vaborbactam"
)

## all antibacterials ----------------------------------------------------------

antibacterials = c(
  broad_spectrum_abx,
  "amoxicillin", "amoxicillin_clavulanate", "ampicillin", "ampicillin_sulbactam",
  "azithromycin", "cefaclor", "cefadroxil", "cefamandole", "cefazolin", "cefdinir",
  "cefixime", "cefmetazole", "cefonicid", "cefoperazone", "cefotaxime", "cefotetan",
  "cefoxitin", "cefpodoxime", "cefprozil", "ceftibuten", "ceftizoxime", "ceftriaxone",
  "cefuroxime", "cephalexin", "cephalothin", "cephapirin", "cephradine", "chloramphenicol",
  "cinoxacin", "ciprofloxacin", "clarithromycin", "clindamycin", "cloxacillin",
  "delafloxacin", "dicloxacillin", "doxycycline", "erythromycin",
  "fosfomycin", "gatifloxacin", "kanamycin", "lefamulin", "levofloxacin", "lincomycin",
  "methicillin", "mezlocillin", "minocycline", "moxifloxacin", "nafcillin",
  "nitrofurantoin", "norfloxacin", "ofloxacin", "omadacycline", "oritavancin",
  "oxacillin", "penicillin", "pivampicillin", "plazomicin", "polymyxin_b",
  "quinupristin_dalfopristin", "rifampin", "streptomycin", "sulbactam_durlobactam",
  "sulfadiazine", "sulfadiazine_trimethoprim", "sulfisoxazole", "telithromycin",
  "tetracycline", "ticarcillin", "ticarcillin_clavulanate", "trimethoprim",
  "trimethoprim_sulfamethoxazole"
)

## antifungals -----------------------------------------------------------------

antifungals = c(
  "amphotericin_b", "anidulafungin", "caspofungin", "fluconazole", "isavuconazole",
  "isavuconazonium", "itraconazole", "micafungin", "posaconazole", "rezafungin", "voriconazole"
)

## antivirals ------------------------------------------------------------------

antivirals = c(
  "acyclovir", "cidofovir", "foscarnet", "ganciclovir", "maribavir", "oseltamivir",
  "peramivir", "ribavirin", "valacyclovir", "valganciclovir"
)

## prophylactic antimicrobials (PCP and bacterial ppx) -------------------------

prophylactic_abx = c(
  "ciprofloxacin", "levofloxacin", "moxifloxacin", "cefpodoxime",                                                 
  "trimethoprim_sulfamethoxazole", "dapsone", "atovaquone"                                                   
)

# load antimicrobial administration data =======================================

antimicrobials_raw = 
  data_list$medication_admin_intermittent |>
  dplyr::filter(med_group == "CMS_sepsis_qualifying_antibiotics") |>
  dplyr::filter(admin_dttm >= start_date) |>
  dplyr::select(hospitalization_id, admin_dttm, med_category) |>
  dplyr::collect()

antimicrobials_raw = 
  join(antimicrobials_raw, all_hid_crosswalk, on = "hospitalization_id", how = "inner", multiple = TRUE) |>
  fselect(patient_id, admin_dttm, med_category) |>
  funique()

antimicrobials_raw = 
  join(antimicrobials_raw, anchor_df, on = "patient_id", how = "inner", multiple = TRUE) |>
  fmutate(
    med_class = case_when(
      tolower(med_category) %in% antibacterials ~ "antibacterial",
      tolower(med_category) %in% antifungals    ~ "antifungal",
      tolower(med_category) %in% antivirals     ~ "antiviral",
      TRUE                                      ~ "other"
    ),
    is_broad_spectrum = as.integer(tolower(med_category) %in% broad_spectrum_abx),
    is_prophylactic   = as.integer(tolower(med_category) %in% prophylactic_abx),
    is_broad_anaerobe = as.integer(tolower(med_category) %in% anaerobic_abx)
  ) |>
  fsubset(med_class != "other")


# baseline antimicrobials (24h before CAR-T) ===================================

baseline_antimicrobials = 
  fsubset(antimicrobials_raw, admin_dttm >= cart_t0 - dhours(24) & admin_dttm < cart_t0)

##identify patients on any antimicrobial class before CAR-T -----------------------------

baseline_by_class = 
  fgroup_by(baseline_antimicrobials, patient_id, med_class) |>
  fsummarize(
    n_doses = fnobs(admin_dttm),
    meds    = paste(unique(tolower(med_category)), collapse = "|")
  ) |>
  fungroup() |>
  pivot_wider(
    names_from   = med_class,
    values_from  = c(n_doses, meds),
    names_glue   = "baseline_{med_class}_{.value}"
  )

baseline_flags = 
  baseline_by_class |>
  fmutate(
    baseline_antibacterial_01 = as.integer(!is.na(baseline_antibacterial_meds)),
    baseline_antifungal_01    = as.integer(!is.na(baseline_antifungal_meds)),
    baseline_antiviral_01     = as.integer(!is.na(baseline_antiviral_meds))
  )

##flag if baseline antibacterial is broad-spectrum --------------------------------------

baseline_broad =
  baseline_antimicrobials |>
  fgroup_by(patient_id) |>
  fsummarize(baseline_broad_spectrum_01 = as.integer(any(med_class == "antibacterial" & is_broad_spectrum == 1L)),
             baseline_anaerobe_01       = as.integer(any(med_class == "antibacterial" & is_broad_anaerobe == 1L))
             ) |>
  fungroup()

baseline_flags = 
  join(baseline_flags, baseline_broad, on = "patient_id", how = "left", multiple = FALSE) |>
  fmutate(baseline_broad_spectrum_01 = replace_na(baseline_broad_spectrum_01, 0L),
          baseline_anaerobe_01       = replace_na(baseline_anaerobe_01, 0L))

##flag if baseline is likely prophylactic --------------------------------------------------

baseline_ppx =
  baseline_antimicrobials |>
  fgroup_by(patient_id) |>
  fsummarize(baseline_prophylactic_01 = as.integer(any(med_class == "antibacterial" & is_prophylactic == 1L))) |>
  fungroup()

baseline_flags = 
  join(baseline_flags, baseline_ppx, on = "patient_id", how = "left", multiple = FALSE) |>
  fmutate(baseline_prophylactic_01 = replace_na(baseline_prophylactic_01, 0L))

rm(baseline_antimicrobials, baseline_by_class, baseline_broad, baseline_ppx)
gc()

# post-CAR-T antimicrobials (0-30 days) ========================================

post_antibacterial =
  antimicrobials_raw |>
  dplyr::filter(admin_dttm >= cart_t0, admin_dttm <= cart_t0 + ddays(30)) |>
  dplyr::filter(med_class == "antibacterial") |>
  fmutate(
    med_lc       = tolower(med_category),
    h_since_cart = as.numeric(difftime(admin_dttm, cart_t0, units = "hours")),
    d_since_cart = as.integer(floor(h_since_cart / 24))
  )

post_antiviral =
  antimicrobials_raw |>
  dplyr::filter(admin_dttm >= cart_t0, admin_dttm <= cart_t0 + ddays(30)) |>
  dplyr::filter(med_class == "antiviral") |>
  fmutate(
    med_lc       = tolower(med_category),
    h_since_cart = as.numeric(difftime(admin_dttm, cart_t0, units = "hours")),
    d_since_cart = as.integer(floor(h_since_cart / 24))
  )

post_antifungal =
  antimicrobials_raw |>
  dplyr::filter(admin_dttm >= cart_t0, admin_dttm <= cart_t0 + ddays(30)) |>
  dplyr::filter(med_class == "antifungal") |>
  fmutate(
    med_lc       = tolower(med_category),
    h_since_cart = as.numeric(difftime(admin_dttm, cart_t0, units = "hours")),
    d_since_cart = as.integer(floor(h_since_cart / 24))
  )

## identify NEW antimicrobials (different from baseline) -----------------------

# create baseline medication lists per patient
baseline_meds_by_patient =
  baseline_flags |>
  fselect(patient_id, baseline_antibacterial_meds) |>
  fmutate(baseline_antibacterial_meds = tolower(replace_na(baseline_antibacterial_meds, "")))


# flag new medications
post_antibacterial_flag =
  post_antibacterial |>
  join(baseline_meds_by_patient, on = "patient_id", how = "left", multiple = FALSE) |>
  fmutate(
    baseline_antibacterial_meds = tolower(replace_na(baseline_antibacterial_meds, "")),
    baseline_abx_list = strsplit(baseline_antibacterial_meds, "\\|"),
    abx_in_baseline = mapply(function(m, lst) m %in% lst, med_lc, baseline_abx_list),
    new_abx = as.integer(baseline_antibacterial_meds == "" | !abx_in_baseline),
    abx_01 = as.integer(med_class == "antibacterial")
    )

##identify total days on antimicrobial----------------------------------

abx_patient_day =
  post_antibacterial_flag |>
  dplyr::group_by(patient_id, d_since_cart) |>
  dplyr::summarise(
    broad_spectrum_01 = as.integer(any(is_broad_spectrum == 1L, na.rm = TRUE)),
    prophylactic_01   = as.integer(any(is_prophylactic   == 1L, na.rm = TRUE)),
    broad_anaerobe_01 = as.integer(any(is_broad_anaerobe == 1L, na.rm = TRUE)),
    new_abx_01        = as.integer(any(new_abx           == 1L, na.rm = TRUE)),
    abx_01            = as.integer(any(abx_01            == 1L, na.rm = TRUE)),
    .groups = "drop"
  )

antiviral_patient_day =
  post_antiviral |>
  dplyr::group_by(patient_id, d_since_cart) |>
  dplyr::summarise(
    antiviral_01 = as.integer(any(med_class == "antiviral", na.rm = TRUE)),
    .groups = "drop"
  )

antifungal_patient_day =
  post_antifungal |>
  dplyr::group_by(patient_id, d_since_cart) |>
  dplyr::summarise(
    antifungal_01 = as.integer(any(med_class == "antifungal", na.rm = TRUE)),
    .groups = "drop"
  )

# Count days on broad spectrum abxs---------------------------------------------
broad_spectrum_days <-
  abx_patient_day |>
  fsubset(broad_spectrum_01 == 1) |>
  fgroup_by(patient_id) |>
  fsummarize(broad_spectrum_days = fnobs(d_since_cart)) |>
  fungroup()

broad_anaerobe_days <-
  abx_patient_day |>
  fsubset(broad_anaerobe_01 == 1) |>
  fgroup_by(patient_id) |>
  fsummarize(broad_anaerobe_days = fnobs(d_since_cart)) |>
  fungroup()

df = join(df, broad_spectrum_days, on = "patient_id", how = "left", multiple = FALSE)
df = join(df, broad_anaerobe_days, on = "patient_id", how = "left", multiple = FALSE)
df$broad_spectrum_days = replace_na(df$broad_spectrum_days, 0L)
df$broad_anaerobe_days = replace_na(df$broad_anaerobe_days, 0L)

rm(broad_spectrum_days, broad_anaerobe_days)
gc()

##72 hour antimicrobial exposure------------------------------------------------

abx_patient_day_72h =
  abx_patient_day |>
  dplyr::filter(d_since_cart >= 0, d_since_cart <= 2) |>
  join(baseline_flags, on = "patient_id", how = "left", multiple = FALSE)

abx_patient_day_72h =
  abx_patient_day_72h |>
  fmutate(
    abx_escalation_01 = case_when(
      new_abx_01 == 1 & baseline_broad_spectrum_01 == 0 & broad_spectrum_01 == 1 ~ 1L,
      TRUE ~ 0L
    ),
    abx_deescalation_01 = case_when(
      new_abx_01 == 1 & baseline_broad_spectrum_01 == 1 & broad_spectrum_01 == 0 ~ 1L,
      TRUE ~ 0L))

antiviral_patient_day_72h =
  antiviral_patient_day |>
  dplyr::filter(d_since_cart >= 0, d_since_cart <= 2) |>
  join(baseline_flags, on = "patient_id", how = "left", multiple = FALSE)

antifungal_patient_day_72h =
  antifungal_patient_day |>
  dplyr::filter(d_since_cart >= 0, d_since_cart <= 2) |>
  join(baseline_flags, on = "patient_id", how = "left", multiple = FALSE)

combined_antimicrobial_day_72h <- 
  abx_patient_day_72h |>
  dplyr::left_join(antiviral_patient_day_72h, by = c("patient_id", "d_since_cart"), relationship = "one-to-one") |>
  dplyr::left_join(antifungal_patient_day_72h, by = c("patient_id", "d_since_cart"), relationship = "one-to-one") |>
  dplyr::mutate(
    dplyr::across(
      c(antiviral_01, antifungal_01),
      ~ tidyr::replace_na(.x, 0L)))


##3-30 day antimicrobial exposure-----------------------------------------------

abx_patient_day_3_30d =
  abx_patient_day |>
  dplyr::filter(d_since_cart >= 3, d_since_cart <= 30) |>
  dplyr::select(-new_abx_01)

antiviral_patient_day_3_30d =
  antiviral_patient_day |>
  dplyr::filter(d_since_cart >= 3, d_since_cart <= 30)

antifungal_patient_day_3_30d =
  antifungal_patient_day |>
  dplyr::filter(d_since_cart >= 3, d_since_cart <= 30)

combined_antimicrobial_day <- 
  abx_patient_day_3_30d |>
  dplyr::left_join(antiviral_patient_day_3_30d, by = c("patient_id", "d_since_cart"), relationship = "one-to-one") |>
  dplyr::left_join(antifungal_patient_day_3_30d, by = c("patient_id", "d_since_cart"), relationship = "one-to-one") |>
  dplyr::mutate(
    dplyr::across(
      c(antiviral_01, antifungal_01),
      ~ tidyr::replace_na(.x, 0L)))

rm(abx_patient_day, abx_patient_day_3_30d, abx_patient_day_72h, antifungal_patient_day, antifungal_patient_day_3_30d, 
antifungal_patient_day_72h, antiviral_patient_day, antiviral_patient_day_3_30d, antiviral_patient_day_72h, post_antibacterial, 
post_antifungal, post_antiviral)
gc()

# antimicrobial trends by trajectory ===========================================
## daily antimicrobial flags by trajectory ---------------------------------

antimicrobial_day_traj_72h = 
  combined_antimicrobial_day_72h |>
  join(fselect(df, patient_id, trajectory), on = "patient_id", how = "inner")

antimicrobial_baseline_traj_summary <-
  antimicrobial_day_traj_72h |>
  fgroup_by(trajectory) |>
  fsummarize(
    n_patients_base_abx      = fsum(baseline_antibacterial_01,  na.rm = TRUE),
    n_patients_base_broad    = fsum(baseline_broad_spectrum_01, na.rm = TRUE),
    n_patients_base_anaerobe = fsum(baseline_anaerobe_01,       na.rm = TRUE),
    n_patients_base_ppx      = fsum(baseline_prophylactic_01,   na.rm = TRUE),
    n_patients_base_af       = fsum(baseline_antifungal_01,     na.rm = TRUE),
    n_patients_base_av       = fsum(baseline_antiviral_01,      na.rm = TRUE)
  ) |>
  fungroup() |>
  roworder(trajectory)

antimicrobial_day_traj_72h_summary <-
  antimicrobial_day_traj_72h |>
  fgroup_by(trajectory, d_since_cart) |>
  fsummarize(
    n_patients_abx         = fsum(abx_01, na.rm = TRUE),
    n_patients_broad       = fsum(broad_spectrum_01, na.rm = TRUE),
    n_patients_anaerobe    = fsum(broad_anaerobe_01, na.rm = TRUE),
    n_patients_ppx         = fsum(prophylactic_01, na.rm = TRUE),
    n_patients_new_abx     = fsum(new_abx_01, na.rm = TRUE),
  
    n_patients_af          = fsum(antifungal_01, na.rm = TRUE),
    n_patients_av          = fsum(antiviral_01, na.rm = TRUE),
    n_patients_abx_esc     = fsum(abx_escalation_01, na.rm = TRUE),
    n_patients_abx_de_esc  = fsum(abx_deescalation_01, na.rm = TRUE)) |>
  fungroup() |>
  roworder(trajectory, d_since_cart)

antimicrobial_day_traj_3_30d = 
  combined_antimicrobial_day |>
  join(fselect(df, patient_id, trajectory), on = "patient_id", how = "inner")

antimicrobial_day_traj_3_30d_summary <-
  antimicrobial_day_traj_3_30d |>
  fgroup_by(trajectory, d_since_cart) |>
  fsummarize(
    n_patients_abx = fsum(abx_01, na.rm = TRUE),
    n_patients_af         = fsum(antifungal_01, na.rm = TRUE),
    n_patients_av         = fsum(antiviral_01, na.rm = TRUE),
    n_patients_broad      = fsum(broad_spectrum_01, na.rm = TRUE),
    n_patients_anaerobe   = fsum(broad_anaerobe_01, na.rm = TRUE),
    n_patients_ppx        = fsum(prophylactic_01, na.rm = TRUE)
  ) |>
  fungroup() |>
  roworder(trajectory, d_since_cart)

save_output(antimicrobial_baseline_traj_summary,  "antimicrobials_baseline_by_traj", shareable = TRUE, subfolder = "trajectory")
save_output(antimicrobial_day_traj_72h_summary,   "antimicrobials_72h_by_traj",      shareable = TRUE, subfolder = "trajectory")
save_output(antimicrobial_day_traj_3_30d_summary, "antimicrobials_3_30d_by_traj",    shareable = TRUE, subfolder = "trajectory")

rm(antimicrobial_day_traj_3_30d, antimicrobial_day_traj_3_30d_summary, antimicrobial_day_traj_72h, antimicrobial_day_traj_72h_summary, 
post_antibacterial_flag, combined_antimicrobial_day, combined_antimicrobial_day_72h, antimicrobials_raw, baseline_flags, 
baseline_meds_by_patient, infection_fungal_daily, infection_bacterial_daily, infection_viral_daily)

rm(anaerobic_abx, antibacterials, antifungals, antivirals, broad_spectrum_abx, lab_categories, 
  prophylactic_abx, req_labs, req_vitals, required_tables, uln_alt, uln_ast, uln_bili, uln_cr, exclude_results, find_max_consecutive)
    
gc()


