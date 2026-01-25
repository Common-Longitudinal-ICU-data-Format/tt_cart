# ==============================================================================
# 08_analysis_circadian.R
# Circadian analysis: time of day → CAR-T outcomes
# Requires: df from 00-06 scripts
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df")) {
  stop("Run 00_setup.R through 06_export.R first.", call. = FALSE)
}

message("Starting 08_analysis_circadian.R")

# packages for analysis ========================================================

if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
library(broom)


# prepare analysis dataset =====================================================

## ensure factor levels --------------------------------------------------------

df$trajectory_f = factor(df$trajectory, levels = c("NT", "HT", "HFR", "HSR"))
df$sex_f        = factor(df$sex_category)
df$race_eth_f   = factor(df$race_eth, levels = c("nhw", "other"))
df$cart_target_f = factor(df$cart_target, levels = c("CD19", "BCMA"))
df$easix_fc_f   = factor(df$easix_fc_cat)
df$year_f       = factor(df$infusion_year)


# ==============================================================================
# DESCRIPTIVE: INFUSION TIME DISTRIBUTION
# ==============================================================================

cat("\n=== Infusion Time Distribution ===\n")
cat("Mean hour:", round(mean(df$cart_hour_precise, na.rm = TRUE), 2), "\n")
cat("SD:", round(sd(df$cart_hour_precise, na.rm = TRUE), 2), "\n")
cat("Range:", round(min(df$cart_hour_precise, na.rm = TRUE), 2), "-", 
    round(max(df$cart_hour_precise, na.rm = TRUE), 2), "\n")

## histogram by hour
hour_dist = table(df$cart_hour_floor)
cat("\nInfusions by hour:\n")
print(hour_dist)


# ==============================================================================
# PRIMARY ANALYSIS: TIME OF DAY → OUTCOMES (CONTINUOUS TIME)
# ==============================================================================

## define outcomes -------------------------------------------------------------

outcomes = list(
  list(name = "90-day mortality",  var = "dead_90_01"),
  list(name = "365-day mortality", var = "dead_365_01"),
  list(name = "ICU admission",     var = "icu_30_01"),
  list(name = "Mechanical vent",   var = "imv_30_01"),
  list(name = "BSI (infection)",   var = "bsi_30_01"),
  list(name = "Tocilizumab (CRS)", var = "toci_any_30d"),
  list(name = "Anakinra (ICANS)",  var = "anakinra_any_30d")
)

## run unadjusted models -------------------------------------------------------

cat("\n=== UNADJUSTED: Time of Day → Outcomes ===\n")
cat("(OR per 1-hour later infusion)\n\n")

unadj_results = lapply(outcomes, \(o) {
  
  if (!o$var %in% names(df)) return(NULL)
  
  formula = as.formula(paste(o$var, "~ cart_hour_precise"))
  
  model = tryCatch(
    glm(formula, data = df, family = binomial),
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NULL)
  
  tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "cart_hour_precise") |>
    mutate(outcome = o$name)
})

unadj_df = do.call(rbind, unadj_results)
print(unadj_df[, c("outcome", "estimate", "conf.low", "conf.high", "p.value")])


## run adjusted models ---------------------------------------------------------

cat("\n=== ADJUSTED: Time of Day → Outcomes ===\n")
cat("(OR per 1-hour later infusion, adjusted for age, sex, race, target, comorbidity, EASIX, year)\n\n")

adj_results = lapply(outcomes, \(o) {
  
  if (!o$var %in% names(df)) return(NULL)
  
  formula = as.formula(paste(
    o$var, 
    "~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  model = tryCatch(
    glm(formula, data = df, family = binomial),
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NULL)
  
  tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "cart_hour_precise") |>
    mutate(outcome = o$name)
})

adj_df = do.call(rbind, adj_results)
print(adj_df[, c("outcome", "estimate", "conf.low", "conf.high", "p.value")])


# ==============================================================================
# SECONDARY: DOES TIME OF DAY PREDICT TRAJECTORY?
# ==============================================================================

cat("\n=== Time of Day → Temperature Trajectory ===\n")

## descriptive: mean infusion hour by trajectory -------------------------------

traj_time = df |>
  fgroup_by(trajectory) |>
  fsummarize(
    n = fnobs(cart_hour_precise),
    mean_hour = fmean(cart_hour_precise, na.rm = TRUE),
    sd_hour = fsd(cart_hour_precise, na.rm = TRUE)
  ) |>
  fungroup()

cat("\nMean infusion hour by trajectory:\n")
print(traj_time)

## multinomial logistic: time → trajectory -------------------------------------

if (requireNamespace("nnet", quietly = TRUE)) {
  
  library(nnet)
  
  # unadjusted
  multi_unadj = multinom(trajectory_f ~ cart_hour_precise, data = df, trace = FALSE)
  
  cat("\nUnadjusted multinomial (ref = NT):\n")
  coefs = summary(multi_unadj)$coefficients
  ses   = summary(multi_unadj)$standard.errors
  z     = coefs / ses
  p     = 2 * (1 - pnorm(abs(z)))
  
  cat("OR (95% CI) per hour later:\n")
  for (traj in rownames(coefs)) {
    or = exp(coefs[traj, "cart_hour_precise"])
    se = ses[traj, "cart_hour_precise"]
    ci_low = exp(coefs[traj, "cart_hour_precise"] - 1.96 * se)
    ci_high = exp(coefs[traj, "cart_hour_precise"] + 1.96 * se)
    pval = p[traj, "cart_hour_precise"]
    cat(sprintf("  %s vs NT: OR=%.2f (%.2f-%.2f), p=%.3f\n", traj, or, ci_low, ci_high, pval))
  }
  
  # adjusted
  multi_adj = multinom(
    trajectory_f ~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + year_f, 
    data = df, 
    trace = FALSE
  )
  
  cat("\nAdjusted multinomial (ref = NT):\n")
  coefs = summary(multi_adj)$coefficients
  ses   = summary(multi_adj)$standard.errors
  z     = coefs / ses
  p     = 2 * (1 - pnorm(abs(z)))
  
  cat("OR (95% CI) per hour later:\n")
  for (traj in rownames(coefs)) {
    or = exp(coefs[traj, "cart_hour_precise"])
    se = ses[traj, "cart_hour_precise"]
    ci_low = exp(coefs[traj, "cart_hour_precise"] - 1.96 * se)
    ci_high = exp(coefs[traj, "cart_hour_precise"] + 1.96 * se)
    pval = p[traj, "cart_hour_precise"]
    cat(sprintf("  %s vs NT: OR=%.2f (%.2f-%.2f), p=%.3f\n", traj, or, ci_low, ci_high, pval))
  }
  
} else {
  message("Install 'nnet' package for multinomial regression: install.packages('nnet')")
}


# ==============================================================================
# MEDIATION: DOES TRAJECTORY MEDIATE TIME → OUTCOME?
# ==============================================================================

cat("\n=== Mediation Analysis: Does Trajectory Mediate Time → Outcomes? ===\n")
cat("Comparing models with and without trajectory adjustment\n\n")

mediation_outcomes = list(
  list(name = "90-day mortality", var = "dead_90_01"),
  list(name = "Anakinra (ICANS)", var = "anakinra_any_30d")
)

for (o in mediation_outcomes) {
  
  if (!o$var %in% names(df)) next
  
  cat("---", o$name, "---\n")
  
  # model without trajectory (total effect)
  formula_total = as.formula(paste(
    o$var, 
    "~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  model_total = tryCatch(
    glm(formula_total, data = df, family = binomial),
    error = function(e) NULL
  )
  
  # model with trajectory (direct effect)
  formula_direct = as.formula(paste(
    o$var, 
    "~ cart_hour_precise + trajectory_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f"
  ))
  
  model_direct = tryCatch(
    glm(formula_direct, data = df, family = binomial),
    error = function(e) NULL
  )
  
  if (!is.null(model_total) && !is.null(model_direct)) {
    
    total_effect = tidy(model_total, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == "cart_hour_precise")
    
    direct_effect = tidy(model_direct, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == "cart_hour_precise")
    
    cat(sprintf("  Total effect (no trajectory):  OR=%.2f (%.2f-%.2f), p=%.3f\n",
                total_effect$estimate, total_effect$conf.low, total_effect$conf.high, total_effect$p.value))
    cat(sprintf("  Direct effect (with trajectory): OR=%.2f (%.2f-%.2f), p=%.3f\n",
                direct_effect$estimate, direct_effect$conf.low, direct_effect$conf.high, direct_effect$p.value))
    
    # percent attenuation
    log_or_total = log(total_effect$estimate)
    log_or_direct = log(direct_effect$estimate)
    pct_mediated = 100 * (log_or_total - log_or_direct) / log_or_total
    cat(sprintf("  Percent mediated by trajectory: %.1f%%\n\n", pct_mediated))
  }
}


# ==============================================================================
# SUBGROUP: BY CAR-T TARGET (CD19 vs BCMA)
# ==============================================================================

cat("\n=== Subgroup Analysis: CD19 vs BCMA ===\n")

for (target in c("CD19", "BCMA")) {
  
  df_sub = fsubset(df, cart_target == target)
  n_sub = fnrow(df_sub)
  
  cat(sprintf("\n--- %s (n=%d) ---\n", target, n_sub))
  
  if (n_sub < 30) {
    cat("  Insufficient sample size\n")
    next
  }
  
  # 90-day mortality
  model = tryCatch(
    glm(dead_90_01 ~ cart_hour_precise + age + sex_f + race_eth_f + vw + easix_fc_f + year_f, 
        data = df_sub, family = binomial),
    error = function(e) NULL
  )
  
  if (!is.null(model)) {
    result = tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == "cart_hour_precise")
    cat(sprintf("  90-day mortality: OR=%.2f (%.2f-%.2f), p=%.3f\n",
                result$estimate, result$conf.low, result$conf.high, result$p.value))
  }
}


# ==============================================================================
# SUMMARY TABLE FOR EXPORT
# ==============================================================================

## combine adjusted results for export -----------------------------------------

if (exists("adj_df") && nrow(adj_df) > 0) {
  
  adj_df$site = site_lowercase
  
  save_output(adj_df, "circadian_logistic_coefs", shareable = TRUE, subfolder = "circadian")
  
  cat("\n=== Logistic regression coefficients saved ===\n")
}


# ==============================================================================
# RMST PSEUDO-REGRESSION (CIRCADIAN → SURVIVAL)
# ==============================================================================

cat("\n=== RMST Pseudo-Regression: Time of Day → 365-day Survival ===\n")

if ("pseudo_rmst_365" %in% names(df)) {
  
  ## Model 1: Total effect (time → RMST, no trajectory) ------------------------
  
  rmst_total = tryCatch({
    lm(pseudo_rmst_365 ~ cart_hour_precise + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  if (!is.null(rmst_total)) {
    
    cat("\nTotal effect (without trajectory adjustment):\n")
    total_results = broom::tidy(rmst_total, conf.int = TRUE)
    time_row = total_results[total_results$term == "cart_hour_precise", ]
    cat(sprintf("  Per hour later: %.1f days (95%% CI: %.1f to %.1f), p=%.3f\n",
                time_row$estimate, time_row$conf.low, time_row$conf.high, time_row$p.value))
    
    # export
    total_coefs = total_results |>
      as.data.table() |>
      fmutate(model = "rmst365_circadian_total", site = site_lowercase, n = nobs(rmst_total))
    
    save_output(total_coefs, "rmst_coefs_circadian_total", shareable = TRUE, subfolder = "circadian")
  }
  
  ## Model 2: Direct effect (time → RMST, with trajectory) ---------------------
  
  rmst_direct = tryCatch({
    lm(pseudo_rmst_365 ~ cart_hour_precise + trajectory_f + age + sex_f + race_eth_f + cart_target_f + vw + easix_fc_f + year_f,
       data = df)
  }, error = function(e) NULL)
  
  if (!is.null(rmst_direct)) {
    
    cat("\nDirect effect (with trajectory adjustment):\n")
    direct_results = broom::tidy(rmst_direct, conf.int = TRUE)
    time_row = direct_results[direct_results$term == "cart_hour_precise", ]
    cat(sprintf("  Per hour later: %.1f days (95%% CI: %.1f to %.1f), p=%.3f\n",
                time_row$estimate, time_row$conf.low, time_row$conf.high, time_row$p.value))
    
    # export
    direct_coefs = direct_results |>
      as.data.table() |>
      fmutate(model = "rmst365_circadian_direct", site = site_lowercase, n = nobs(rmst_direct))
    
    save_output(direct_coefs, "rmst_coefs_circadian_direct", shareable = TRUE, subfolder = "circadian")
  }
  
  ## Mediation summary ---------------------------------------------------------
  
  if (!is.null(rmst_total) && !is.null(rmst_direct)) {
    
    total_est  = total_results[total_results$term == "cart_hour_precise", "estimate"]
    direct_est = direct_results[direct_results$term == "cart_hour_precise", "estimate"]
    
    pct_mediated = 100 * (total_est - direct_est) / total_est
    
    cat(sprintf("\nMediation by trajectory: %.1f%% of effect explained\n", pct_mediated))
  }
  
} else {
  cat("pseudo_rmst_365 not found - run 06_export.R with 'pseudo' package installed.\n")
}


message("\n08_analysis_circadian.R complete.")
