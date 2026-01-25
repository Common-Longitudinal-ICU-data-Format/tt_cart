# ==============================================================================
# 02_traj.R
# Temperature trajectory assignment for CAR-T project
# Must run while temp df is still in memory from 01_cohort.R
# Outputs: trajectory assignments added to df, figure data for pooling
# ==============================================================================

# check dependencies ===========================================================

if (!exists("df") || !exists("temp") || !exists("save_output")) {
 stop("Run 00_setup.R and 01_cohort.R first.", call. = FALSE
)
}

message("Starting 02_traj.R")

# prepare hourly temperature data ==============================================

## filter to 0-72h post-infusion -----------------------------------------------

temp = fsubset(temp, recorded_dttm >= cart_t0)
temp = fsubset(temp, recorded_dttm <= cart_t0 + dhours(72))
temp = funique(temp)

## get discharge times for censoring -------------------------------------------

frame = fselect(df, joined_hosp_id, discharge_dttm)

## create hourly grid and fill -------------------------------------------------

temp_hourly = 
  fselect(temp, joined_hosp_id, patient_id, cart_t0, temp_time = recorded_dttm, temp_c) |>
  roworder(joined_hosp_id, temp_time) |>
  fgroup_by(joined_hosp_id) |>
  fmutate(t0 = ffirst(cart_t0)) |>
  fungroup()

temp_hourly = 
  group_by(temp_hourly, joined_hosp_id) |>
  complete(temp_time = seq.POSIXt(from = t0[1], to = t0[1] + dhours(72), by = "1 hour")) |>
  fill(everything(), .direction = "down") |>
  fill(everything(), .direction = "up") |>
  ungroup()

temp_hourly = 
  join(temp_hourly, frame, how = "inner") |>
  fsubset(temp_time <= discharge_dttm) |>
  fmutate(hour = as.integer(round(as.numeric(difftime(temp_time, cart_t0, units = "hours"))))) |>
  fmutate(temp_z = standardize(temp_c))

rm(frame); gc()

# fit to pre-specified trajectory curves =======================================

## trajectory polynomial coefficients (from prior study) -----------------------
## traj1 = HT  (hypothermic)
## traj2 = NT  (normothermic)
## traj3 = HFR (hyperthermic fast recovery)
## traj4 = HSR (hyperthermic slow recovery)

temp_model = 
  fselect(temp_hourly, joined_hosp_id, hour, temp_z) |>
  fmutate(
    traj_ht  = -0.89548 - 0.00298*hour + 0.00010*(hour^2),
    traj_nt  = -0.00667 + 0.00050*hour - 0.00001*(hour^2),
    traj_hfr =  1.35157 - 0.06946*hour + 0.00065*(hour^2),
    traj_hsr =  1.22203 - 0.00590*hour - 0.00007*(hour^2),
    error_ht  = (temp_z - traj_ht)^2,
    error_nt  = (temp_z - traj_nt)^2,
    error_hfr = (temp_z - traj_hfr)^2,
    error_hsr = (temp_z - traj_hsr)^2
  ) 

## calculate SSE for each trajectory -------------------------------------------

temp_model_fit = 
  fgroup_by(temp_model, joined_hosp_id) |>
  fsummarise(
    sse_ht  = fsum(error_ht,  na.rm = TRUE),
    sse_nt  = fsum(error_nt,  na.rm = TRUE),
    sse_hfr = fsum(error_hfr, na.rm = TRUE),
    sse_hsr = fsum(error_hsr, na.rm = TRUE)
  ) |>
  fmutate(
    sse_best = pmin(sse_ht, sse_nt, sse_hfr, sse_hsr),
    sse_second = pmap_dbl(
      list(sse_ht, sse_nt, sse_hfr, sse_hsr),
      ~ sort(c(...))[2]
    ),
    margin_abs = sse_second - sse_best,
    margin_rel = margin_abs / sse_best,
    trajectory = case_when(
      sse_best == sse_hsr ~ "HSR",
      sse_best == sse_hfr ~ "HFR",
      sse_best == sse_nt  ~ "NT",
      sse_best == sse_ht  ~ "HT",
      TRUE                ~ NA_character_
    ),
    second_best = case_when(
      sse_second == sse_hsr ~ "HSR",
      sse_second == sse_hfr ~ "HFR",
      sse_second == sse_nt  ~ "NT",
      sse_second == sse_ht  ~ "HT",
      TRUE                  ~ NA_character_
    ),
    fit_flag = fifelse(margin_rel < 0.10, "borderline", "clear")
  ) |>
  fselect(joined_hosp_id, trajectory, second_best, margin_rel, fit_flag)

## add trajectory to main df ---------------------------------------------------

df = join(df, temp_model_fit, how = "left")

rm(temp_model, temp_model_fit); gc()

# create figure data for pooling ===============================================

## add trajectory to temp_hourly for grouping ----------------------------------

temp_hourly = 
  join(temp_hourly, fselect(df, joined_hosp_id, trajectory, cart_group), how = "left")

## figure data by trajectory group ---------------------------------------------

fig_data_by_traj = 
  fgroup_by(temp_hourly, trajectory, hour) |>
  fsummarize(
    n_obs       = fnobs(temp_c),
    sum_temp    = fsum(temp_c, na.rm = TRUE),
    sum_temp_sq = fsum(temp_c^2, na.rm = TRUE)
  ) |>
  fungroup() |>
  roworder(trajectory, hour)

save_output(fig_data_by_traj, "fig_temp_by_trajectory", shareable = TRUE)

## figure data by CAR-T product ------------------------------------------------

fig_data_by_product = 
  fgroup_by(temp_hourly, cart_group, hour) |>
  fsummarize(
    n_obs       = fnobs(temp_c),
    sum_temp    = fsum(temp_c, na.rm = TRUE),
    sum_temp_sq = fsum(temp_c^2, na.rm = TRUE)
  ) |>
  fungroup() |>
  roworder(cart_group, hour)

save_output(fig_data_by_product, "fig_temp_by_product", shareable = TRUE)

# clean up large temp objects ==================================================

rm(temp, temp_hourly, fig_data_by_traj, fig_data_by_product); gc()

message("02_traj.R complete. Trajectory assigned to df, temp removed from memory.")
