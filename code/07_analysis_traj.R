# ==============================================================================
# 07_analysis_traj.R
# Trajectory analysis: temperature trajectory → CAR-T outcomes
# Exports shareable coefficient tables for meta-analysis
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("save_output")) {
  stop("Run 00_setup.R through 06_export.R first.", call. = FALSE)
}

message("Starting 07_analysis_traj.R")

# packages =====================================================================

if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
library(broom)

# constants ====================================================================

VARIANTS = c("main", "se_borderline_reassign")

.allowed = list(
  table1      = c("continuous", "categorical"),
  table2      = c("binary", "time"),
  survival    = c("q30d"),
  models      = c("logistic_coefs", "rmst_coefs"),
  subgroup    = c("logistic_coefs", "rmst_coefs"),
  sensitivity = c("logistic_coefs", "rmst_coefs"),
  diagnostics = c("overall", "by_trajectory")
)

# helper functions =============================================================

write_artifact = function(dt, analysis, artifact, site,
                          strata = NULL, variant = NULL,
                          root = file.path(project_location, "upload_to_box", "trajectory")) {
  
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

## sensitivity: reassign borderline to second_best -----------------------------

df$trajectory_reassign = fifelse(
  df$fit_flag == "borderline" & !is.na(df$second_best),
  df$second_best,
  df$trajectory
)
df$trajectory_reassign_f = factor(df$trajectory_reassign, levels = c("NT", "HT", "HFR", "HSR"))

message("  N = ", fnrow(df))
message("  Borderline fits: ", fsum(df$fit_flag == "borderline", na.rm = TRUE))

# ==============================================================================
# DIAGNOSTICS
# ==============================================================================

message("\n== Diagnostics ==")

diag_overall = tidytable(
  site        = site_lowercase,
  n_total     = fnrow(df),
  n_cd19      = fsum(df$cart_target == "CD19", na.rm = TRUE),
  n_bcma      = fsum(df$cart_target == "BCMA", na.rm = TRUE),
  n_dead_90   = fsum(df$dead_90_01, na.rm = TRUE),
  n_dead_365  = fsum(df$dead_365_01, na.rm = TRUE),
  n_icu       = fsum(df$icu_30_01, na.rm = TRUE),
  n_imv       = fsum(df$imv_30_01, na.rm = TRUE),
  n_toci      = fsum(df$toci_any_30d, na.rm = TRUE),
  n_anakinra  = fsum(df$anakinra_any_30d, na.rm = TRUE),
  rate_dead_90 = fmean(df$dead_90_01, na.rm = TRUE),
  n_borderline = fsum(df$fit_flag == "borderline", na.rm = TRUE)
)

write_artifact(diag_overall, "diagnostics", "overall", site_lowercase)

diag_by_traj = 
  fgroup_by(df, trajectory) |>
  fsummarize(
    n           = fnobs(patient_id),
    n_dead_90   = fsum(dead_90_01, na.rm = TRUE),
    n_dead_365  = fsum(dead_365_01, na.rm = TRUE),
    n_icu       = fsum(icu_30_01, na.rm = TRUE),
    n_toci      = fsum(toci_any_30d, na.rm = TRUE),
    n_anakinra  = fsum(anakinra_any_30d, na.rm = TRUE),
    rate_dead_90 = fmean(dead_90_01, na.rm = TRUE)
  ) |>
  fungroup() |>
  fmutate(site = site_lowercase)

write_artifact(diag_by_traj, "diagnostics", "by_trajectory", site_lowercase)

# ==============================================================================
# TABLE 1: PATIENT CHARACTERISTICS BY TRAJECTORY
# ==============================================================================

message("\n== Table 1 ==")

## continuous variables --------------------------------------------------------

t1_cont_vars = c("age", "vw", "easix", "log2_easix", "wbc_infusion", "lymph_abs_infusion",
                 "baseline_ldh", "baseline_crp", "baseline_ferritin")

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
      fmutate(variable = v, site = site_lowercase)
  }) |>
  rbindlist()

write_artifact(t1_continuous, "table1", "continuous", site_lowercase)

## categorical variables -------------------------------------------------------

t1_cat_vars = c("sex_category", "race_eth", "cancer_cat", "cart_target",
                "easix_fc_cat", "hct_01", "infusion_year")

t1_categorical = 
  lapply(t1_cat_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    tmp = tidytable(trajectory = df$trajectory, cat_val = as.character(df[[v]]))
    fgroup_by(tmp, trajectory, cat_val) |>
      fsummarize(n = fnobs(trajectory)) |>
      fungroup() |>
      fmutate(variable = v, site = site_lowercase) |>
      frename(category = cat_val)
  }) |>
  rbindlist()

write_artifact(t1_categorical, "table1", "categorical", site_lowercase)

# ==============================================================================
# TABLE 2: OUTCOMES BY TRAJECTORY
# ==============================================================================

message("\n== Table 2 ==")

## binary outcomes -------------------------------------------------------------

t2_binary_vars = c("dead_30_01", "dead_90_01", "dead_180_01", "dead_365_01",
                   "icu_30_01", "imv_30_01", "va_30_01", "bsi_30_01",
                   "toci_any_30d", "anakinra_any_30d", "steroid_any_72h",
                   "severe_icaht_01", "hyperferritin_5k_01", "hyperferritin_10k_01")

t2_binary = 
  lapply(t2_binary_vars, \(v) {
    if (!v %in% names(df)) return(NULL)
    x = df[[v]]
    fselect(df, trajectory) |>
      fmutate(x = x) |>
      fgroup_by(trajectory) |>
      fsummarize(n = fnobs(x), n_events = fsum(x, na.rm = TRUE)) |>
      fungroup() |>
      fmutate(variable = v, site = site_lowercase)
  }) |>
  rbindlist()

write_artifact(t2_binary, "table2", "binary", site_lowercase)

## time variables --------------------------------------------------------------

t2_time_vars = c("los_hospital_d", "h_until_toci", "h_until_anakinra", "h_until_steroids")

t2_time = 
  lapply(t2_time_vars, \(v) {
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
      fmutate(variable = v, site = site_lowercase)
  }) |>
  rbindlist()

write_artifact(t2_time, "table2", "time", site_lowercase)

# ==============================================================================
# SURVIVAL DATA (for KM curves)
# ==============================================================================

message("\n== Survival data ==")

time_points = seq(0, 365, by = 30)

surv_q30d = 
  lapply(time_points, \(t) {
    prev_t = ifelse(t == 0, 0, t - 30)
    fgroup_by(df, trajectory) |>
      fsummarize(
        n_at_risk  = fsum(followup_d >= prev_t, na.rm = TRUE),
        n_events   = fsum(died == 1 & time_to_death > prev_t & time_to_death <= t, na.rm = TRUE),
        n_censored = fsum(died == 0 & followup_d > prev_t & followup_d <= t, na.rm = TRUE)
      ) |>
      fungroup() |>
      fmutate(time_point = t, site = site_lowercase)
  }) |>
  rbindlist()

write_artifact(surv_q30d, "survival", "q30d", site_lowercase)

# ==============================================================================
# MAIN MODELS: TRAJECTORY → OUTCOMES
# ==============================================================================

message("\n== Main models ==")

## logistic regression for binary outcomes -------------------------------------

logistic_outcomes = list(
  list(name = "dead_90",   var = "dead_90_01"),
  list(name = "dead_365",  var = "dead_365_01"),
  list(name = "icu_30",    var = "icu_30_01"),
  list(name = "imv_30",    var = "imv_30_01"),
  list(name = "va_30",     var = "va_30_01"),
  list(name = "bsi_30",    var = "bsi_30_01"),
  list(name = "toci",      var = "toci_any_30d"),
  list(name = "anakinra",  var = "anakinra_any_30d")
)

logistic_results = lapply(logistic_outcomes, \(o) {
  
  if (!o$var %in% names(df)) return(NULL)
  
  formula = as.formula(paste(
    o$var, "~ trajectory_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  model = tryCatch(glm(formula, data = df, family = binomial), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  
  tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    as.data.table() |>
    fmutate(outcome = o$name, site = site_lowercase, n = nobs(model))
})

logistic_coefs = rbindlist(logistic_results, fill = TRUE)
write_artifact(logistic_coefs, "models", "logistic_coefs", site_lowercase)

## RMST pseudo-regression ------------------------------------------------------

if ("pseudo_rmst_365" %in% names(df)) {
  
  rmst_model = tryCatch({
    lm(pseudo_rmst_365 ~ trajectory_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  if (!is.null(rmst_model)) {
    rmst_coefs = tidy(rmst_model, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(outcome = "rmst_365", site = site_lowercase, n = nobs(rmst_model))
    
    write_artifact(rmst_coefs, "models", "rmst_coefs", site_lowercase)
  }
} else {
  message("  pseudo_rmst_365 not found - skipping RMST")
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
    
    # remove cart_target from model since we're stratified
    formula = as.formula(paste(
      o$var, "~ trajectory_f + age + sex_f + race_eth_f + vw + easix_fc_f + year_f"
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
      lm(pseudo_rmst_365 ~ trajectory_f + age + sex_f + race_eth_f + vw + easix_fc_f + year_f,
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
# SENSITIVITY: BORDERLINE REASSIGNMENT
# ==============================================================================

message("\n== Sensitivity: borderline reassignment ==")

n_reassigned = fsum(df$trajectory != df$trajectory_reassign, na.rm = TRUE)
message("  Reassigned: ", n_reassigned, " patients")

if (n_reassigned > 0) {
  
  # logistic with reassigned trajectory
  sens_logistic = lapply(logistic_outcomes, \(o) {
    if (!o$var %in% names(df)) return(NULL)
    
    formula = as.formula(paste(
      o$var, "~ trajectory_reassign_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
    ))
    
    model = tryCatch(glm(formula, data = df, family = binomial), error = function(e) NULL)
    if (is.null(model)) return(NULL)
    
    tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(
        term = gsub("trajectory_reassign_f", "trajectory_f", term),
        outcome = o$name, 
        site = site_lowercase, 
        n = nobs(model),
        variant = "borderline_reassign"
      )
  })
  
  sens_coefs = rbindlist(sens_logistic, fill = TRUE)
  write_artifact(sens_coefs, "sensitivity", "logistic_coefs", site_lowercase, variant = "borderline_reassign")
  
  # RMST with reassigned
  if ("pseudo_rmst_365" %in% names(df)) {
    rmst_sens = tryCatch({
      lm(pseudo_rmst_365 ~ trajectory_reassign_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
         data = df)
    }, error = function(e) NULL)
    
    if (!is.null(rmst_sens)) {
      rmst_sens_coefs = tidy(rmst_sens, conf.int = TRUE) |>
        as.data.table() |>
        fmutate(
          term = gsub("trajectory_reassign_f", "trajectory_f", term),
          outcome = "rmst_365",
          site = site_lowercase,
          n = nobs(rmst_sens),
          variant = "borderline_reassign"
        )
      
      write_artifact(rmst_sens_coefs, "sensitivity", "rmst_coefs", site_lowercase, variant = "borderline_reassign")
    }
  }
}

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

message("\n== Summary ==")

cat("\nTrajectory distribution:\n")
print(table(df$trajectory))

cat("\n90-day mortality by trajectory:\n")
mort_tbl = df |>
  fgroup_by(trajectory) |>
  fsummarize(n = fnobs(dead_90_01), deaths = fsum(dead_90_01, na.rm = TRUE)) |>
  fmutate(pct = round(100 * deaths / n, 1))
print(mort_tbl)

if (exists("logistic_coefs") && nrow(logistic_coefs) > 0) {
  cat("\nAdjusted ORs for 90-day mortality (trajectory effects):\n")
  traj_ors = logistic_coefs[outcome == "dead_90" & grepl("trajectory", term), 
                            .(term, estimate, conf.low, conf.high, p.value)]
  print(traj_ors)
}

message("\n07_analysis_traj.R complete.")
message("Artifacts written to: ", file.path(project_location, "upload_to_box", "trajectory"))
