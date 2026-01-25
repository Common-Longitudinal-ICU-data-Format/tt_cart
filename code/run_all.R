# ==============================================================================
# run_all.R
# Master script to run CAR-T temperature trajectory project
# ==============================================================================

# clear environment
rm(list = ls())
gc()

# record start time
start_time = Sys.time()

message("\n", paste(rep("=", 60), collapse = ""))
message("CAR-T Temperature Trajectory Project")
message("Run started: ", start_time)
message(paste(rep("=", 60), collapse = ""), "\n")

# source scripts in order
source("00_setup.R")
message("\n", paste(rep("-", 40), collapse = ""), "\n")

source("01_cohort.R")
message("\n", paste(rep("-", 40), collapse = ""), "\n")

source("02_traj.R")
message("\n", paste(rep("-", 40), collapse = ""), "\n")

source("03_covariates.R")
message("\n", paste(rep("-", 40), collapse = ""), "\n")

source("04_outcomes.R")
message("\n", paste(rep("-", 40), collapse = ""), "\n")

source("05_easix.R")
message("\n", paste(rep("-", 40), collapse = ""), "\n")

source("06_export.R")

# report completion
end_time = Sys.time()
elapsed  = round(difftime(end_time, start_time, units = "mins"), 2)

message("\n", paste(rep("=", 60), collapse = ""))
message("Project complete!")
message("Elapsed time: ", elapsed, " minutes")
message("Final cohort N: ", fnrow(df))
message(paste(rep("=", 60), collapse = ""), "\n")
