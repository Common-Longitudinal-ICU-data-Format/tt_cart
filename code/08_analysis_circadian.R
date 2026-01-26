# ==============================================================================
# 08_analysis_circadian.R
# Circadian analysis: time of day → CAR-T outcomes
# Exports shareable coefficient tables for meta-analysis
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("save_output")) {
  stop("Run 00_setup.R through 06_export.R first.", call. = FALSE)
}

# ensure site_lowercase is set
if (!exists("site_lowercase") || is.null(site_lowercase) || site_lowercase == "") {
  site_lowercase = "site"
  message("Warning: site_lowercase not set, using 'site' as default")
}

message("Starting 08_analysis_circadian.R")

# packages =====================================================================

if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
library(broom)

# constants ====================================================================

EARLY_LATE_CUTOFF = 15  # 15:00 (3 PM) as in original paper

.allowed = list(
  table1      = c("continuous", "categorical"),
  table2      = c("binary"),
  survival    = c("earlylate"),
  models      = c("logistic_coefs", "rmst_coefs"),
  mediation   = c("logistic_coefs", "rmst_coefs"),
  subgroup    = c("logistic_coefs", "rmst_coefs"),
  sensitivity = c("logistic_coefs", "rmst_coefs"),
  diagnostics = c("overall", "by_hour")
)

# helper functions =============================================================

write_artifact = function(dt, analysis, artifact, site,
                          strata = NULL, variant = NULL,
                          root = file.path(project_location, "upload_to_box", "circadian")) {
  
  stopifnot(analysis %in% names(.allowed), artifact %in% .allowed[[analysis]])
  
  if (!dir.exists(root)) dir.create(root, recursive = TRUE, showWarnings = FALSE)
  
  parts = paste0(analysis, "-", artifact)
  if (!is.null(strata) && nzchar(strata)) parts = paste0(parts, "-", strata)
  if (!is.null(variant) && nzchar(variant)) parts = paste0(parts, "-", variant)
  fn = paste0(parts, "-", site, ".csv")
  
  path = file.path(root, fn)
  fwrite(dt, path)
  message("  Wrote: ", fn)
  invisible(path)
}

# prepare analysis data ========================================================

message("\n== Preparing analysis data ==")

## factors ---------------------------------------------------------------------

df$trajectory_f  = factor(df$trajectory, levels = c("NT", "HT", "HFR", "HSR"))
df$sex_f         = factor(df$sex_category)
df$race_eth_f    = factor(df$race_eth, levels = c("nhw", "other"))
df$cart_target_f = factor(df$cart_target, levels = c("CD19", "BCMA"))
df$easix_fc_f    = factor(df$easix_fc_cat)
df$year_f        = factor(df$infusion_year)

## early/late dichotomization --------------------------------------------------

df$early_late = fifelse(df$cart_hour_precise < EARLY_LATE_CUTOFF, "early", "late")
df$early_late_f = factor(df$early_late, levels = c("early", "late"))

## placeholder for sunrise-adjusted time ---------------------------------------
## TODO: Calculate cart_hour_sunrise from:
##   - latitude (from config$latitude)
##   - date of infusion (admin_dttm)
##   - Use suncalc::getSunlightTimes() or similar
## df$cart_hour_sunrise = ... 

message("  N = ", fnrow(df))
message("  Early (<", EARLY_LATE_CUTOFF, ":00): ", fsum(df$early_late == "early", na.rm = TRUE))
message("  Late (>=", EARLY_LATE_CUTOFF, ":00): ", fsum(df$early_late == "late", na.rm = TRUE))
message("  Mean infusion hour: ", round(fmean(df$cart_hour_precise, na.rm = TRUE), 2))

# ==============================================================================
# DIAGNOSTICS
# ==============================================================================

message("\n== Diagnostics ==")

diag_overall = tidytable(
  site           = site_lowercase,
  n_total        = fnrow(df),
  n_early        = fsum(df$early_late == "early", na.rm = TRUE),
  n_late         = fsum(df$early_late == "late", na.rm = TRUE),
  hour_mean      = fmean(df$cart_hour_precise, na.rm = TRUE),
  hour_sd        = fsd(df$cart_hour_precise, na.rm = TRUE),
  hour_min       = fmin(df$cart_hour_precise, na.rm = TRUE),
  hour_max       = fmax(df$cart_hour_precise, na.rm = TRUE),
  n_dead_90      = fsum(df$dead_90_01, na.rm = TRUE),
  rate_dead_90   = fmean(df$dead_90_01, na.rm = TRUE)
)

write_artifact(diag_overall, "diagnostics", "overall", site_lowercase)

diag_by_hour = 
  fgroup_by(df, cart_hour_floor) |>
  fsummarize(
    n          = fnobs(patient_id),
    n_dead_90  = fsum(dead_90_01, na.rm = TRUE),
    n_icu      = fsum(icu_30_01, na.rm = TRUE),
    n_toci     = fsum(toci_any_30d, na.rm = TRUE),
    n_anakinra = fsum(anakinra_any_30d, na.rm = TRUE)
  ) |>
  fungroup()

diag_by_hour$site = site_lowercase

write_artifact(diag_by_hour, "diagnostics", "by_hour", site_lowercase)

# ==============================================================================
# TABLE 1: PATIENT CHARACTERISTICS BY EARLY/LATE
# ==============================================================================

message("\n== Table 1 (early vs late) ==")

## continuous variables --------------------------------------------------------

t1_cont_vars = c("age", "vw", "easix", "log2_easix", "cart_hour_precise")

t1_continuous = 
  lapply(t1_cont_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    x = df[[v]]
    out = fselect(df, early_late) |>
      fmutate(x = x, x_sq = x^2) |>
      fgroup_by(early_late) |>
      fsummarize(
        n        = fnobs(x),
        n_miss   = fsum(is.na(x)),
        sum_x    = fsum(x, na.rm = TRUE),
        sum_x_sq = fsum(x_sq, na.rm = TRUE)
      ) |>
      fungroup()
    out$variable = v
    out$site = site_lowercase
    out
  }) |>
  rbindlist()

write_artifact(t1_continuous, "table1", "continuous", site_lowercase)

## categorical variables -------------------------------------------------------

t1_cat_vars = c("sex_category", "race_eth", "cart_target", "trajectory", "infusion_year")

t1_categorical = 
  lapply(t1_cat_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    tmp = tidytable(early_late = df$early_late, cat_val = as.character(df[[v]]))
    out = fgroup_by(tmp, early_late, cat_val) |>
      fsummarize(n = fnobs(early_late)) |>
      fungroup()
    out$variable = v
    out$site = site_lowercase
    setnames(out, "cat_val", "category")
    out
  }) |>
  rbindlist()

write_artifact(t1_categorical, "table1", "categorical", site_lowercase)

# ==============================================================================
# TABLE 2: OUTCOMES BY EARLY/LATE
# ==============================================================================

message("\n== Table 2 (early vs late) ==")

t2_binary_vars = c("dead_30_01", "dead_90_01", "dead_365_01",
                   "icu_30_01", "imv_30_01", "va_30_01", "bsi_30_01",
                   "toci_any_30d", "anakinra_any_30d")

t2_binary = 
  lapply(t2_binary_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    x = df[[v]]
    out = fselect(df, early_late) |>
      fmutate(x = x) |>
      fgroup_by(early_late) |>
      fsummarize(n = fnobs(x), n_events = fsum(x, na.rm = TRUE)) |>
      fungroup()
    out$variable = v
    out$site = site_lowercase
    out
  }) |>
  rbindlist()

write_artifact(t2_binary, "table2", "binary", site_lowercase)

# ==============================================================================
# SURVIVAL DATA (for KM curves by early/late)
# ==============================================================================

message("\n== Survival data (early vs late) ==")

time_points = seq(0, 365, by = 30)

surv_earlylate = 
  lapply(time_points, \(t) {
    prev_t = ifelse(t == 0, 0, t - 30)
    out = fgroup_by(df, early_late) |>
      fsummarize(
        n_at_risk  = fsum(followup_d >= prev_t, na.rm = TRUE),
        n_events   = fsum(died == 1 & time_to_death > prev_t & time_to_death <= t, na.rm = TRUE),
        n_censored = fsum(died == 0 & followup_d > prev_t & followup_d <= t, na.rm = TRUE)
      ) |>
      fungroup()
    out$time_point = t
    out$site = site_lowercase
    out
  }) |>
  rbindlist()

write_artifact(surv_earlylate, "survival", "earlylate", site_lowercase)

# ==============================================================================
# MAIN MODELS: TIME OF DAY → OUTCOMES (CONTINUOUS)
# ==============================================================================

message("\n== Main models (continuous time) ==")

logistic_outcomes = list(
  list(name = "dead_90",   var = "dead_90_01"),
  list(name = "dead_365",  var = "dead_365_01"),
  list(name = "icu_30",    var = "icu_30_01"),
  list(name = "imv_30",    var = "imv_30_01"),
  list(name = "bsi_30",    var = "bsi_30_01"),
  list(name = "toci",      var = "toci_any_30d"),
  list(name = "anakinra",  var = "anakinra_any_30d")
)

## logistic: time (continuous) → outcomes --------------------------------------

logistic_results = lapply(logistic_outcomes, \(o) {
  
  if (!o$var %in% names(df)) return(NULL)
  
  formula = as.formula(paste(
    o$var, "~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  model = tryCatch(glm(formula, data = df, family = binomial), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  
  tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    as.data.table() |>
    fmutate(outcome = o$name, site = site_lowercase, n = nobs(model), time_var = "clock")
})

logistic_coefs = rbindlist(logistic_results, fill = TRUE)
write_artifact(logistic_coefs, "models", "logistic_coefs", site_lowercase)

## RMST: time (continuous) → 365-day survival ----------------------------------

if ("pseudo_rmst_365" %in% names(df)) {
  
  rmst_model = tryCatch({
    lm(pseudo_rmst_365 ~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  if (!is.null(rmst_model)) {
    rmst_coefs = tidy(rmst_model, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(outcome = "rmst_365", site = site_lowercase, n = nobs(rmst_model), time_var = "clock")
    
    write_artifact(rmst_coefs, "models", "rmst_coefs", site_lowercase)
  }
} else {
  message("  pseudo_rmst_365 not found - skipping RMST")
}

# ==============================================================================
# MEDIATION: DOES TRAJECTORY MEDIATE TIME → OUTCOMES?
# ==============================================================================

message("\n== Mediation analysis ==")

## logistic: total effect vs direct effect -------------------------------------

mediation_outcomes = list(
  list(name = "dead_90",  var = "dead_90_01"),
  list(name = "dead_365", var = "dead_365_01"),
  list(name = "anakinra", var = "anakinra_any_30d")
)

mediation_logistic = lapply(mediation_outcomes, \(o) {
  
  if (!o$var %in% names(df)) return(NULL)
  
  # total effect (no trajectory)
  formula_total = as.formula(paste(
    o$var, "~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  # direct effect (with trajectory)
  formula_direct = as.formula(paste(
    o$var, "~ cart_hour_precise + trajectory_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  model_total = tryCatch(glm(formula_total, data = df, family = binomial), error = function(e) NULL)
  model_direct = tryCatch(glm(formula_direct, data = df, family = binomial), error = function(e) NULL)
  
  if (is.null(model_total) || is.null(model_direct)) return(NULL)
  
  total_coefs = tidy(model_total, exponentiate = TRUE, conf.int = TRUE) |>
    as.data.table() |>
    fmutate(outcome = o$name, site = site_lowercase, n = nobs(model_total), effect = "total")
  
  direct_coefs = tidy(model_direct, exponentiate = TRUE, conf.int = TRUE) |>
    as.data.table() |>
    fmutate(outcome = o$name, site = site_lowercase, n = nobs(model_direct), effect = "direct")
  
  rbind(total_coefs, direct_coefs)
})

mediation_logistic_coefs = rbindlist(mediation_logistic, fill = TRUE)
write_artifact(mediation_logistic_coefs, "mediation", "logistic_coefs", site_lowercase)

## RMST: total vs direct -------------------------------------------------------

if ("pseudo_rmst_365" %in% names(df)) {
  
  rmst_total = tryCatch({
    lm(pseudo_rmst_365 ~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  rmst_direct = tryCatch({
    lm(pseudo_rmst_365 ~ cart_hour_precise + trajectory_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  if (!is.null(rmst_total) && !is.null(rmst_direct)) {
    
    total_coefs = tidy(rmst_total, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(outcome = "rmst_365", site = site_lowercase, n = nobs(rmst_total), effect = "total")
    
    direct_coefs = tidy(rmst_direct, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(outcome = "rmst_365", site = site_lowercase, n = nobs(rmst_direct), effect = "direct")
    
    mediation_rmst_coefs = rbind(total_coefs, direct_coefs)
    write_artifact(mediation_rmst_coefs, "mediation", "rmst_coefs", site_lowercase)
    
    # print mediation summary
    total_est = total_coefs[term == "cart_hour_precise", estimate]
    direct_est = direct_coefs[term == "cart_hour_precise", estimate]
    pct_med = 100 * (total_est - direct_est) / total_est
    message("  RMST mediation by trajectory: ", round(pct_med, 1), "% of effect explained")
  }
}

# ==============================================================================
# SUBGROUP: CD19 vs BCMA
# ==============================================================================

message("\n== Subgroup: CD19 vs BCMA ==")

for (target in c("CD19", "BCMA")) {
  
  df_sub = fsubset(df, cart_target == target)
  n_sub = fnrow(df_sub)
  
  if (n_sub < 30) {
    message("  ", target, ": insufficient N (", n_sub, ")")
    next
  }
  
  message("  ", target, ": N = ", n_sub)
  
  # logistic
  sub_logistic = lapply(logistic_outcomes, \(o) {
    if (!o$var %in% names(df_sub)) return(NULL)
    
    formula = as.formula(paste(
      o$var, "~ cart_hour_precise + age + sex_f + race_eth_f + vw + easix_fc_f + year_f"
    ))
    
    model = tryCatch(glm(formula, data = df_sub, family = binomial), error = function(e) NULL)
    if (is.null(model)) return(NULL)
    
    tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(outcome = o$name, site = site_lowercase, n = nobs(model), subgroup = target)
  })
  
  sub_coefs = rbindlist(sub_logistic, fill = TRUE)
  if (nrow(sub_coefs) > 0) {
    write_artifact(sub_coefs, "subgroup", "logistic_coefs", site_lowercase, strata = tolower(target))
  }
  
  # RMST
  if ("pseudo_rmst_365" %in% names(df_sub)) {
    rmst_sub = tryCatch({
      lm(pseudo_rmst_365 ~ cart_hour_precise + age + sex_f + race_eth_f + vw + easix_fc_f + year_f,
         data = df_sub)
    }, error = function(e) NULL)
    
    if (!is.null(rmst_sub)) {
      rmst_sub_coefs = tidy(rmst_sub, conf.int = TRUE) |>
        as.data.table() |>
        fmutate(outcome = "rmst_365", site = site_lowercase, n = nobs(rmst_sub), subgroup = target)
      
      write_artifact(rmst_sub_coefs, "subgroup", "rmst_coefs", site_lowercase, strata = tolower(target))
    }
  }
}

# ==============================================================================
# SUBGROUP: SEX
# ==============================================================================

message("\n== Subgroup: Sex ==")

for (sex in c("male", "female")) {
  
  df_sub = fsubset(df, tolower(sex_category) == sex)
  n_sub = fnrow(df_sub)
  
  if (n_sub < 30) {
    message("  ", sex, ": insufficient N (", n_sub, ")")
    next
  }
  
  message("  ", sex, ": N = ", n_sub)
  
  # logistic for key outcomes
  sub_logistic = lapply(logistic_outcomes[1:3], \(o) {  # dead_90, dead_365, icu
    if (!o$var %in% names(df_sub)) return(NULL)
    
    formula = as.formula(paste(
      o$var, "~ cart_hour_precise + age + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
    ))
    
    model = tryCatch(glm(formula, data = df_sub, family = binomial), error = function(e) NULL)
    if (is.null(model)) return(NULL)
    
    tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(outcome = o$name, site = site_lowercase, n = nobs(model), subgroup = sex)
  })
  
  sub_coefs = rbindlist(sub_logistic, fill = TRUE)
  if (nrow(sub_coefs) > 0) {
    write_artifact(sub_coefs, "subgroup", "logistic_coefs", site_lowercase, strata = sex)
  }
}

# ==============================================================================
# SENSITIVITY: SUNRISE-ADJUSTED TIME (PLACEHOLDER)
# ==============================================================================

message("\n== Sensitivity: sunrise-adjusted time ==")

## TODO: Implement when latitude available in config
## 
## if (!is.null(config$latitude) && "cart_hour_sunrise" %in% names(df)) {
##   
##   # repeat main models with cart_hour_sunrise instead of cart_hour_precise
##   
##   logistic_sunrise = lapply(logistic_outcomes, \(o) { ... })
##   write_artifact(..., variant = "sunrise")
##   
##   rmst_sunrise = lm(pseudo_rmst_365 ~ cart_hour_sunrise + ..., data = df)
##   write_artifact(..., variant = "sunrise")
## }

message("  Sunrise adjustment not implemented - requires latitude in config")

# ==============================================================================
# SENSITIVITY: VARYING EARLY/LATE CUTPOINT
# ==============================================================================

message("\n== Sensitivity: varying cutpoint ==")

cutpoints = c(12, 13, 14, 15, 16, 17)

cutpoint_results = lapply(cutpoints, \(cp) {
  
  df$early_cp = fifelse(df$cart_hour_precise < cp, 1L, 0L)
  
  # 90-day mortality
  model = tryCatch({
    glm(dead_90_01 ~ early_cp + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
        data = df, family = binomial)
  }, error = function(e) NULL)
  
  if (is.null(model)) return(NULL)
  
  coef_row = tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    as.data.table() |>
    fsubset(term == "early_cp")
  
  coef_row$cutpoint = cp
  coef_row$n_early = fsum(df$early_cp == 1, na.rm = TRUE)
  coef_row$n_late = fsum(df$early_cp == 0, na.rm = TRUE)
  coef_row$site = site_lowercase
  coef_row
})

cutpoint_coefs = rbindlist(cutpoint_results, fill = TRUE)
if (nrow(cutpoint_coefs) > 0) {
  write_artifact(cutpoint_coefs, "sensitivity", "logistic_coefs", site_lowercase, variant = "cutpoint")
}

# ==============================================================================
# SENSITIVITY: HFR COLLAPSED INTO NT (FOR MEDIATION)
# ==============================================================================

message("\n== Sensitivity: HFR collapsed into NT (mediation) ==")

n_hfr = fsum(df$trajectory == "HFR", na.rm = TRUE)
message("  HFR patients: ", n_hfr)

# collapse HFR into NT
df$trajectory_collapsed = fifelse(df$trajectory == "HFR", "NT", df$trajectory)
df$trajectory_collapsed_f = factor(df$trajectory_collapsed, levels = c("NT", "HT", "HSR"))

## logistic mediation with collapsed trajectory --------------------------------

mediation_collapsed_logistic = lapply(mediation_outcomes, \(o) {
  
  if (!o$var %in% names(df)) return(NULL)
  
  # total effect (no trajectory) - same as before
  formula_total = as.formula(paste(
    o$var, "~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  # direct effect with COLLAPSED trajectory
  formula_direct = as.formula(paste(
    o$var, "~ cart_hour_precise + trajectory_collapsed_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  model_total = tryCatch(glm(formula_total, data = df, family = binomial), error = function(e) NULL)
  model_direct = tryCatch(glm(formula_direct, data = df, family = binomial), error = function(e) NULL)
  
  if (is.null(model_total) || is.null(model_direct)) return(NULL)
  
  total_coefs = tidy(model_total, exponentiate = TRUE, conf.int = TRUE) |>
    as.data.table() |>
    fmutate(outcome = o$name, site = site_lowercase, n = nobs(model_total), effect = "total", variant = "hfr_collapsed")
  
  direct_coefs = tidy(model_direct, exponentiate = TRUE, conf.int = TRUE) |>
    as.data.table() |>
    fmutate(
      term = gsub("trajectory_collapsed_f", "trajectory_f", term),
      outcome = o$name, 
      site = site_lowercase, 
      n = nobs(model_direct), 
      effect = "direct",
      variant = "hfr_collapsed"
    )
  
  rbind(total_coefs, direct_coefs)
})

mediation_collapsed_coefs = rbindlist(mediation_collapsed_logistic, fill = TRUE)
write_artifact(mediation_collapsed_coefs, "sensitivity", "logistic_coefs", site_lowercase, variant = "hfr_collapsed")

## RMST mediation with collapsed trajectory ------------------------------------

if ("pseudo_rmst_365" %in% names(df)) {
  
  rmst_total = tryCatch({
    lm(pseudo_rmst_365 ~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  rmst_direct_collapsed = tryCatch({
    lm(pseudo_rmst_365 ~ cart_hour_precise + trajectory_collapsed_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  if (!is.null(rmst_total) && !is.null(rmst_direct_collapsed)) {
    
    total_coefs = tidy(rmst_total, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(outcome = "rmst_365", site = site_lowercase, n = nobs(rmst_total), effect = "total", variant = "hfr_collapsed")
    
    direct_coefs = tidy(rmst_direct_collapsed, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(
        term = gsub("trajectory_collapsed_f", "trajectory_f", term),
        outcome = "rmst_365",
        site = site_lowercase,
        n = nobs(rmst_direct_collapsed),
        effect = "direct",
        variant = "hfr_collapsed"
      )
    
    rmst_collapsed_coefs = rbind(total_coefs, direct_coefs)
    write_artifact(rmst_collapsed_coefs, "sensitivity", "rmst_coefs", site_lowercase, variant = "hfr_collapsed")
    
    # print mediation summary
    total_est = total_coefs[term == "cart_hour_precise", estimate]
    direct_est = direct_coefs[term == "cart_hour_precise", estimate]
    pct_med = 100 * (total_est - direct_est) / total_est
    message("  RMST mediation (HFR collapsed): ", round(pct_med, 1), "% of effect explained")
  }
}

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

message("\n== Summary ==")

cat("\nInfusion time distribution:\n")
cat("  Mean:", round(fmean(df$cart_hour_precise, na.rm = TRUE), 2), "\n")
cat("  SD:", round(fsd(df$cart_hour_precise, na.rm = TRUE), 2), "\n")
cat("  Range:", round(fmin(df$cart_hour_precise, na.rm = TRUE), 2), "-",
    round(fmax(df$cart_hour_precise, na.rm = TRUE), 2), "\n")

cat("\nEarly vs Late (cutoff =", EARLY_LATE_CUTOFF, ":00):\n")
print(table(df$early_late))

if (exists("logistic_coefs") && nrow(logistic_coefs) > 0) {
  cat("\nAdjusted ORs per hour later (clock time):\n")
  time_ors = logistic_coefs[term == "cart_hour_precise", 
                            .(outcome, estimate, conf.low, conf.high, p.value)]
  print(time_ors)
}

message("\n08_analysis_circadian.R complete.")
message("Artifacts written to: ", file.path(project_location, "upload_to_box", "circadian"))
