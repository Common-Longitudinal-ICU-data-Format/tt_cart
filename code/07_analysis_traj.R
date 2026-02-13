# ==============================================================================
# 07_analysis_traj.R
# Trajectory analysis: temperature trajectory → CAR-T outcomes
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

# safe tidy for glm - falls back to Wald CIs if profiling fails
safe_tidy_glm = function(model, exponentiate = TRUE) {
  tryCatch(
    tidy(model, exponentiate = exponentiate, conf.int = TRUE),
    error = function(e) {
      out = tidy(model, exponentiate = exponentiate, conf.int = FALSE)
      if (exponentiate) {
        out$conf.low = exp(log(out$estimate) - 1.96 * out$std.error)
        out$conf.high = exp(log(out$estimate) + 1.96 * out$std.error)
      } else {
        out$conf.low = out$estimate - 1.96 * out$std.error
        out$conf.high = out$estimate + 1.96 * out$std.error
      }
      out
    }
  )
}

# prepare analysis data ========================================================

message("\n== Preparing analysis data ==")

## factors ---------------------------------------------------------------------

df$trajectory_f  = factor(df$trajectory, levels = c("NT", "HT", "HFR", "HSR"))
df$sex_f         = factor(df$sex_category)
df$race_eth_f    = factor(df$race_eth, levels = c("nhw", "other"))
df$cart_target_f = factor(df$cart_target, levels = c("CD19", "BCMA"))
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
# MAIN MODELS: TRAJECTORY → OUTCOMES
# ==============================================================================

message("\n== Main models ==")

## logistic regression for binary outcomes -------------------------------------

logistic_outcomes = list(
  list(name = "dead_90",         var = "dead_90_01"),
  list(name = "dead_365",        var = "dead_365_01"),
  list(name = "icu_30",          var = "icu_30_01"),
  list(name = "resp_30",         var = "resp_30_01"),
  list(name = "va_30",           var = "va_30_01"),
  list(name = "infection_30",    var = "infection_overall_01"),
  list(name = "toci",            var = "toci_any_30d"),
  list(name = "anakinra",        var = "anakinra_any_30d")
)

logistic_results = lapply(logistic_outcomes, \(o) {
  
  if (!o$var %in% names(df)) return(NULL)
  
  formula = as.formula(paste(
    o$var, "~ trajectory_f + age + sex_f + race_eth_f + cart_target_f + vw + log2_m_easix + year_f"
  ))
  
  model = tryCatch(glm(formula, data = df, family = binomial), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  
  safe_tidy_glm(model) |>
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
    
    safe_tidy_glm(model) |>
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
    
    safe_tidy_glm(model) |>
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
# SENSITIVITY: HFR COLLAPSED INTO NT
# ==============================================================================

message("\n== Sensitivity: HFR collapsed into NT ==")

n_hfr = fsum(df$trajectory == "HFR", na.rm = TRUE)
message("  HFR patients: ", n_hfr)

# collapse HFR into NT
df$trajectory_collapsed = fifelse(df$trajectory == "HFR", "NT", df$trajectory)
df$trajectory_collapsed_f = factor(df$trajectory_collapsed, levels = c("NT", "HT", "HSR"))

message("  After collapse:")
message("    NT (includes HFR): ", fsum(df$trajectory_collapsed == "NT", na.rm = TRUE))
message("    HT: ", fsum(df$trajectory_collapsed == "HT", na.rm = TRUE))
message("    HSR: ", fsum(df$trajectory_collapsed == "HSR", na.rm = TRUE))

# logistic with collapsed trajectory
sens_collapsed_logistic = lapply(logistic_outcomes, \(o) {
  if (!o$var %in% names(df)) return(NULL)
  
  formula = as.formula(paste(
    o$var, "~ trajectory_collapsed_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  model = tryCatch(glm(formula, data = df, family = binomial), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  
  safe_tidy_glm(model) |>
    as.data.table() |>
    fmutate(
      term = gsub("trajectory_collapsed_f", "trajectory_f", term),
      outcome = o$name, 
      site = site_lowercase, 
      n = nobs(model),
      variant = "hfr_collapsed"
    )
})

sens_collapsed_coefs = rbindlist(sens_collapsed_logistic, fill = TRUE)
write_artifact(sens_collapsed_coefs, "sensitivity", "logistic_coefs", site_lowercase, variant = "hfr_collapsed")

# RMST with collapsed
if ("pseudo_rmst_365" %in% names(df)) {
  rmst_collapsed = tryCatch({
    lm(pseudo_rmst_365 ~ trajectory_collapsed_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  if (!is.null(rmst_collapsed)) {
    rmst_collapsed_coefs = tidy(rmst_collapsed, conf.int = TRUE) |>
      as.data.table() |>
      fmutate(
        term = gsub("trajectory_collapsed_f", "trajectory_f", term),
        outcome = "rmst_365",
        site = site_lowercase,
        n = nobs(rmst_collapsed),
        variant = "hfr_collapsed"
      )
    
    write_artifact(rmst_collapsed_coefs, "sensitivity", "rmst_coefs", site_lowercase, variant = "hfr_collapsed")
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
