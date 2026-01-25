# ==============================================================================
# 07_analysis_traj.R
# Analysis script for trajectory project
# Can run on single-site data OR pooled multi-site data
# ==============================================================================

# This script can be run in two modes:
# 1. SINGLE SITE: After running 00-06, using df in memory
# 2. POOLED: Loading CSVs from multiple sites

# Set mode here:
analysis_mode = "single"  # "single" or "pooled"
pooled_data_dir = "path/to/upload_to_box/trajectory"  # only used if mode = "pooled"


# ==============================================================================
# SETUP
# ==============================================================================

library(data.table)
library(tidytable)
library(collapse)

if (analysis_mode == "pooled") {
  
  # function to load and stack site files
  pool_files = function(pattern, data_dir) {
    files = list.files(data_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) return(NULL)
    
    lapply(files, \(f) {
      dt = fread(f)
      # extract site from filename (assumes: name_SITE_date.csv)
      parts = strsplit(basename(f), "_")[[1]]
      dt$site = parts[length(parts) - 1]
      dt
    }) |> rbindlist(fill = TRUE)
  }
  
  # load pooled data
  t1_cont     = pool_files("^table1_continuous", pooled_data_dir)
  t1_cat      = pool_files("^table1_categorical", pooled_data_dir)
  t2_mort     = pool_files("^table2_mortality", pooled_data_dir)
  t2_acute    = pool_files("^table2_acute", pooled_data_dir)
  t2_time     = pool_files("^table2_time", pooled_data_dir)
  surv_data   = pool_files("^survival_curve_data", pooled_data_dir)
  pseudo_data = pool_files("^rmst_pseudo", pooled_data_dir)
  misclass    = pool_files("^sensitivity_misclass", pooled_data_dir)
  site_summ   = pool_files("^site_summary", pooled_data_dir)
  
} else {
  # single site mode - check df exists
  if (!exists("df")) stop("Run scripts 00-06 first, or set mode to 'pooled'")
}


# ==============================================================================
# TABLE 1: PATIENT CHARACTERISTICS
# ==============================================================================

if (analysis_mode == "pooled") {
  
  ## pool continuous variables -------------------------------------------------
  
  table1_cont_pooled = t1_cont |>
    fgroup_by(trajectory, variable) |>
    fsummarize(
      n      = fsum(n),
      n_miss = fsum(n_miss),
      sum_x  = fsum(sum_x),
      sum_sq = fsum(sum_x_sq)
    ) |>
    fungroup() |>
    fmutate(
      mean = sum_x / n,
      var  = (sum_sq - (sum_x^2 / n)) / (n - 1),
      sd   = sqrt(var)
    ) |>
    fselect(trajectory, variable, n, n_miss, mean, sd)
  
  ## pool categorical variables ------------------------------------------------
  
  table1_cat_pooled = t1_cat |>
    fgroup_by(trajectory, variable, category) |>
    fsummarize(n = fsum(n)) |>
    fungroup() |>
    fgroup_by(trajectory, variable) |>
    fmutate(pct = 100 * n / fsum(n)) |>
    fungroup()
  
} else {
  
  ## single site - calculate directly ------------------------------------------
  
  t1_cont_vars = c("age", "vw", "easix", "log2_easix")
  
  table1_cont_pooled = 
    lapply(t1_cont_vars, \(v) {
      x = df[[v]]
      fselect(df, trajectory) |>
        fmutate(x = x) |>
        fgroup_by(trajectory) |>
        fsummarize(
          n    = fnobs(x),
          mean = fmean(x, na.rm = TRUE),
          sd   = fsd(x, na.rm = TRUE)
        ) |>
        fungroup() |>
        fmutate(variable = v)
    }) |>
    rbindlist()
  
  t1_cat_vars = c("sex_category", "race_eth", "cart_target", "infusion_year")
  
  table1_cat_pooled = 
    lapply(t1_cat_vars, \(v) {
      tmp = tidytable(trajectory = df$trajectory, category = as.character(df[[v]]))
      fgroup_by(tmp, trajectory, category) |>
        fsummarize(n = fnobs(trajectory)) |>
        fungroup() |>
        fgroup_by(trajectory) |>
        fmutate(pct = 100 * n / fsum(n)) |>
        fungroup() |>
        fmutate(variable = v)
    }) |>
    rbindlist()
}

print(table1_cont_pooled)
print(table1_cat_pooled)


# ==============================================================================
# TABLE 2: OUTCOMES
# ==============================================================================

if (analysis_mode == "pooled") {
  
  ## mortality -----------------------------------------------------------------
  
  mortality_pooled = t2_mort |>
    fgroup_by(trajectory) |>
    fsummarize(
      n_30  = fsum(n_30d_denom),  dead_30  = fsum(n_dead_30),
      n_90  = fsum(n_90d_denom),  dead_90  = fsum(n_dead_90),
      n_180 = fsum(n_180d_denom), dead_180 = fsum(n_dead_180),
      n_365 = fsum(n_365d_denom), dead_365 = fsum(n_dead_365)
    ) |>
    fungroup() |>
    fmutate(
      pct_30  = 100 * dead_30  / n_30,
      pct_90  = 100 * dead_90  / n_90,
      pct_180 = 100 * dead_180 / n_180,
      pct_365 = 100 * dead_365 / n_365
    )
  
  ## acute outcomes ------------------------------------------------------------
  
  acute_pooled = t2_acute |>
    fgroup_by(trajectory, variable) |>
    fsummarize(n = fsum(n), n_events = fsum(n_events)) |>
    fungroup() |>
    fmutate(pct = 100 * n_events / n)
  
} else {
  
  ## single site mortality -----------------------------------------------------
  
  mortality_pooled = 
    fgroup_by(df, trajectory) |>
    fsummarize(
      n_30  = fnobs(dead_30_01),
      dead_30 = fsum(dead_30_01, na.rm = TRUE),
      n_90  = fnobs(dead_90_01),
      dead_90 = fsum(dead_90_01, na.rm = TRUE),
      n_180 = fnobs(dead_180_01),
      dead_180 = fsum(dead_180_01, na.rm = TRUE),
      n_365 = fnobs(dead_365_01),
      dead_365 = fsum(dead_365_01, na.rm = TRUE)
    ) |>
    fungroup() |>
    fmutate(
      pct_30  = 100 * dead_30  / n_30,
      pct_90  = 100 * dead_90  / n_90,
      pct_180 = 100 * dead_180 / n_180,
      pct_365 = 100 * dead_365 / n_365
    )
  
  ## single site acute outcomes ------------------------------------------------
  
  acute_vars = c("icu_30_01", "imv_30_01", "va_30_01", "bsi_30_01", 
                 "toci_any_30d", "anakinra_any_30d")
  
  acute_pooled = 
    lapply(acute_vars, \(v) {
      if (!v %in% names(df)) return(NULL)
      x = df[[v]]
      fselect(df, trajectory) |>
        fmutate(x = x) |>
        fgroup_by(trajectory) |>
        fsummarize(n = fnobs(x), n_events = fsum(x, na.rm = TRUE)) |>
        fungroup() |>
        fmutate(variable = v, pct = 100 * n_events / n)
    }) |>
    rbindlist()
}

print(mortality_pooled)
print(acute_pooled)


# ==============================================================================
# SURVIVAL CURVES (POOLED KM RECONSTRUCTION)
# ==============================================================================

if (analysis_mode == "pooled" && !is.null(surv_data)) {
  
  # This is approximate - for proper pooled KM, use IPD meta-analysis
  # But this gives reasonable visualization
  
  surv_pooled = surv_data |>
    fgroup_by(trajectory, time_point) |>
    fsummarize(
      n_at_risk  = fsum(n_at_risk),
      n_events   = fsum(n_events),
      n_censored = fsum(n_censored)
    ) |>
    fungroup() |>
    roworder(trajectory, time_point)
  
  # calculate hazard and survival
  surv_pooled = surv_pooled |>
    fgroup_by(trajectory) |>
    fmutate(
      hazard   = n_events / n_at_risk,
      surv_prob = cumprod(1 - hazard)
    ) |>
    fungroup()
  
  print(surv_pooled)
  
  # example plot (requires ggplot2)
  # library(ggplot2)
  # ggplot(surv_pooled, aes(x = time_point, y = surv_prob, color = trajectory)) +
  #   geom_step() +
  #   labs(x = "Days since CAR-T", y = "Survival probability") +
  #   theme_minimal()
}


# ==============================================================================
# RMST PSEUDO-REGRESSION
# ==============================================================================

if (analysis_mode == "pooled" && !is.null(pseudo_data) && requireNamespace("geepack", quietly = TRUE)) {
  
  library(geepack)
  
  # prep data - trajectory as factor with NT as reference
  pseudo_data$trajectory = factor(pseudo_data$trajectory, levels = c("NT", "HT", "HFR", "HSR"))
  pseudo_data$easix_fc_cat = factor(pseudo_data$easix_fc_cat)
  
  # 90-day RMST model
  model_rmst90 = geeglm(
    pseudo_rmst_90 ~ trajectory + age + sex_category + race_eth + 
                     cart_target + vw + easix_fc_cat + infusion_year,
    data    = pseudo_data,
    id      = site,
    family  = gaussian,
    corstr  = "exchangeable"
  )
  
  cat("\n=== 90-day RMST Model ===\n")
  print(summary(model_rmst90))
  
  # 365-day RMST model (restricted to those with follow-up)
  pseudo_365 = fsubset(pseudo_data, !is.na(pseudo_rmst_365))
  
  if (fnrow(pseudo_365) > 50) {
    model_rmst365 = geeglm(
      pseudo_rmst_365 ~ trajectory + age + sex_category + race_eth + 
                        cart_target + vw + easix_fc_cat + infusion_year,
      data    = pseudo_365,
      id      = site,
      family  = gaussian,
      corstr  = "exchangeable"
    )
    
    cat("\n=== 365-day RMST Model ===\n")
    print(summary(model_rmst365))
  }
  
} else if (analysis_mode == "single") {
  
  # single site logistic regression
  if (requireNamespace("broom", quietly = TRUE)) {
    
    df$trajectory_f = factor(df$trajectory, levels = c("NT", "HT", "HFR", "HSR"))
    
    # 90-day mortality
    model_90d = glm(
      dead_90_01 ~ trajectory_f + age + sex_category + race_eth + 
                   cart_target + vw + easix_fc_cat + infusion_year,
      data   = fsubset(df, has_90d_fu == 1),
      family = binomial
    )
    
    cat("\n=== 90-day Mortality (Logistic) ===\n")
    print(broom::tidy(model_90d, exponentiate = TRUE, conf.int = TRUE))
    
    # ICU admission
    model_icu = glm(
      icu_30_01 ~ trajectory_f + age + sex_category + race_eth + 
                  cart_target + vw + easix_fc_cat + infusion_year,
      data   = df,
      family = binomial
    )
    
    cat("\n=== ICU Admission (Logistic) ===\n")
    print(broom::tidy(model_icu, exponentiate = TRUE, conf.int = TRUE))
  }
}


# ==============================================================================
# SENSITIVITY: MISCLASSIFICATION
# ==============================================================================

if (analysis_mode == "pooled" && !is.null(misclass)) {
  
  misclass_pooled = misclass |>
    fgroup_by(trajectory, second_best, fit_flag) |>
    fsummarize(
      n         = fsum(n),
      dead_90   = fsum(dead_90_n),
      icu_30    = fsum(icu_30_n),
      toci      = fsum(toci_n),
      anakinra  = fsum(anakinra_n)
    ) |>
    fungroup() |>
    fmutate(
      pct_dead_90 = 100 * dead_90 / n,
      pct_icu     = 100 * icu_30 / n,
      pct_toci    = 100 * toci / n
    )
  
  cat("\n=== Misclassification Sensitivity ===\n")
  print(misclass_pooled)
  
  # compare clear vs borderline within each trajectory
  sens_comparison = misclass_pooled |>
    fgroup_by(trajectory, fit_flag) |>
    fsummarize(
      n         = fsum(n),
      dead_90   = fsum(dead_90),
      pct_dead  = 100 * fsum(dead_90) / fsum(n)
    ) |>
    fungroup()
  
  print(sens_comparison)
}


# ==============================================================================
# STRATIFIED ANALYSIS: CD19 vs BCMA
# ==============================================================================

if (analysis_mode == "single") {
  
  cat("\n=== Stratified by CAR-T Target ===\n")
  
  # CD19 cohort
  df_cd19 = fsubset(df, cart_target == "CD19")
  cat("\nCD19 cohort N:", fnrow(df_cd19), "\n")
  
  traj_by_target = df |>
    fgroup_by(cart_target, trajectory) |>
    fsummarize(
      n        = fnobs(patient_id),
      dead_90  = fsum(dead_90_01, na.rm = TRUE),
      icu_30   = fsum(icu_30_01, na.rm = TRUE)
    ) |>
    fungroup() |>
    fmutate(
      pct_dead = 100 * dead_90 / n,
      pct_icu  = 100 * icu_30 / n
    )
  
  print(traj_by_target)
}


message("\n07_analysis_traj.R complete.")
