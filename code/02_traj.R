
# assign temperature trajectory ------------------------------------------------

## hourly temperature data -----------------------------------------------------

standardize <- function(x) {(x - fmean(x, na.rm = TRUE)) / fsd(x, na.rm = TRUE)}

frame = fselect(df, joined_hosp_id, discharge_dttm)

temp = fsubset(temp, recorded_dttm >= cart_t0)
temp = fsubset(temp, recorded_dttm <= cart_t0 + dhours(72))
temp = distinct(temp)

temp_hourly = 
  rename(temp, temp_time = recorded_dttm) %>%
  roworder(joined_hosp_id, temp_time) |> 
  fgroup_by(joined_hosp_id) %>%
  fmutate(t0 = ffirst(cart_t0)) %>%
  fungroup() 

temp_hourly = 
  group_by(temp_hourly, joined_hosp_id) |>
  complete(temp_time = seq.POSIXt(from = t0[1], to = t0[1] + dhours(72), by = "1 hour")) %>%
  fill(everything(), .direction = "down") %>%
  fill(everything(), .direction = "up") %>%
  ungroup() 

temp_hourly = 
  join(temp_hourly, frame, how = "inner") |>
  fsubset(temp_time <= discharge_dttm) %>%
  fmutate(hour   = as.numeric(difftime(temp_time, cart_t0, units = "hours"))) |>
  fmutate(temp_z = standardize(temp_c))

## fit curves and calculate errors ---------------------------------------------

temp_model = 
  fselect(temp_hourly, joined_hosp_id, hour, temp_z) |>
  fmutate(
    traj1  = -0.89548 - 0.00298*hour + 0.00010*(hour^2),
    traj2  = -0.00667 + 0.00050*hour - 0.00001*(hour^2),
    traj3  =  1.35157 - 0.06946*hour + 0.00065*(hour^2),
    traj4  =  1.22203 - 0.00590*hour - 0.00007*(hour^2),
    error1 = (temp_z - traj1)^2,
    error2 = (temp_z - traj2)^2,
    error3 = (temp_z - traj3)^2,
    error4 = (temp_z - traj4)^2
  ) 

temp_model_fit = 
  fgroup_by(temp_model, joined_hosp_id) |>
  fsummarise(
    sse_ht  = fsum(error1, na.rm = TRUE),
    sse_nt  = fsum(error2, na.rm = TRUE),
    sse_hfr = fsum(error3, na.rm = TRUE),
    sse_hsr = fsum(error4, na.rm = TRUE)
  ) |>
  fmutate(
    sse_best   = pmin(sse_ht, sse_nt, sse_hfr, sse_hsr),
    sse_second = pmap_dbl(
      list(sse_ht, sse_nt, sse_hfr, sse_hsr),
      ~ sort(c(...))[2]
    ),
    margin_abs = sse_second-sse_best,
    margin_rel = margin_abs/sse_best,
    trajectory = case_when(
      sse_best   == sse_hsr ~ "HSR",
      sse_best   == sse_hfr ~ "HFR",
      sse_best   == sse_nt  ~ "NT",
      sse_best   == sse_ht  ~ "HT",
      TRUE                  ~ NA_character_
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






temp_plot = 
  fgroup_by(temp_hourly, hour) |>
  fsummarize(
    n       = fnobs(temp_c),
    mean    = fmean(temp_c),
    sd      = fsd(temp_c)
  ) |>
  fmutate(
    se      = sd / sqrt(n),
    ci_low  = mean - 1.96 * se,
    ci_high = mean + 1.96 * se
  )

temp_plot = 
  fgroup_by(temp_hourly, hour) |>
  fsummarize(
    n      = fnobs(temp_c),
    sum_x  = fsum(temp_c),
    sum_x2 = fsum(temp_c^2)
  )
