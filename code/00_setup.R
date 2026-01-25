# ==============================================================================
# 00_setup.R
# Setup, validation, and helper functions for CAR-T temperature trajectory project
# ==============================================================================

# install and load packages ====================================================

packages_to_install =
  c(
    "data.table",   
    "tidyverse", 
    "lubridate",
    "rprojroot",
    "tidytable",    
    "collapse",     
    "janitor",      
    "readxl",       
    "arrow",        
    "readr",        
    "yaml",         
    "here",         
    "fst",          
    "ps",
    "comorbidity"
  )

packages_to_load = 
  c(
    "rprojroot",
    "data.table",
    "lubridate",
    "tidytable",    
    "collapse", 
    "stringr",
    "arrow",        
    "here",         
    "ps"
  )

fn_install_if_missing = function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {              
    try(install.packages(p, dependencies = TRUE), silent = TRUE)  
  }
}

fn_load_quiet = function(p) {
  suppressPackageStartupMessages(                     
    library(p, character.only = TRUE)             
  )
}

invisible(lapply(packages_to_install, fn_install_if_missing))
invisible(lapply(packages_to_load, fn_load_quiet))

options(dplyr.summarise.inform = FALSE) 

rm(packages_to_install, packages_to_load, fn_install_if_missing, fn_load_quiet)  
gc()                          


# threads and ram ==============================================================

os_type   = Sys.info()[["sysname"]]
all_cores = parallel::detectCores(logical = TRUE)
all_cores = if (is.na(all_cores)) 1L else as.integer(all_cores)

get_ram_gb = function() {
  tryCatch({
    if (os_type == "Darwin") {
      bytes = suppressWarnings(
        as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
      )
      if (length(bytes) > 0 && !is.na(bytes)) bytes / 1024^3 else NA_real_
    } else {
      mem_info = ps::ps_system_memory()
      key = intersect(c("available", "avail"), names(mem_info))
      val = if (length(key) > 0) mem_info[[key]] / 1024^3 else NA_real_
      if (is.finite(val)) val else {
        if (file.exists("/proc/meminfo")) {
          kb = suppressWarnings(
            as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
          )
          if (length(kb) > 0 && !is.na(kb)) kb / 1024^2 else NA_real_
        } else {
          NA_real_
        }
      }
    }
  }, error = function(e) NA_real_)
}

avail_ram_gb = get_ram_gb()

## choose threads conservatively -----------------------------------------------

reserve_cores   = 1L
gb_per_thread   = 0.50
max_by_cores    = max(1L, all_cores - reserve_cores)
max_by_memory   = if (is.finite(avail_ram_gb)) max(1L, floor(avail_ram_gb / gb_per_thread)) else max_by_cores
n_threads       = as.integer(max(1L, min(max_by_cores, max_by_memory, 8L)))
n_math_threads  = as.integer(max(1L, min(n_threads, 8L)))

## apply thread settings -------------------------------------------------------

data.table::setDTthreads(threads = n_threads)
collapse::set_collapse(nthreads  = n_threads)
options(arrow.use_threads        = TRUE)
Sys.setenv(ARROW_NUM_THREADS     = n_threads)
options(mc.cores                 = n_threads)

message(
  sprintf("Env OK | OS=%s | Cores=%d | Threads=%d | MathThreads=%d | Avail RAM≈%s GB",
          os_type, all_cores, n_threads, n_math_threads,
          ifelse(is.finite(avail_ram_gb), round(avail_ram_gb, 1), "NA"))
)


# site details =================================================================

project_root = find_root(rprojroot::has_dir("config"))
config       = yaml::read_yaml(file.path(project_root, "config", "config.yaml"))

allowed_files = c("parquet", "csv", "fst")

## user enters site details ----------------------------------------------------

site_lowercase   = tolower(config$site_lowercase)
tz_local         = config$time_zone
file_type        = tolower(config$file_type)  
tables_location  = config$tables_location
project_location = config$project_location

## file_type must be one of c(parquet, csv, fst) -------------------------------

if (!(file_type %in% allowed_files)) {
  stop(
    paste0(
      "Invalid '", file_type,
      "'. Expected one of: ", 
      paste(allowed_files, collapse = ", "), call. = FALSE
    )
  )
}

## create output directories ---------------------------------------------------

output_dirs = c(
  "intermediate_output",
  "final_output",
  "upload_to_box",
  "upload_to_box/trajectory",
  "upload_to_box/circadian",
  "proj_tables"
)

walk(output_dirs, ~ {
  p = file.path(project_location, .x)
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
})

rm(output_dirs)

## dates -----------------------------------------------------------------------

start_date = ymd_hms("2018-01-01 00:00:00", tz = tz_local) 
end_date   = ymd_hms("2024-12-31 23:59:59", tz = tz_local)
start_date = with_tz(start_date, tz = "UTC") 
end_date   = with_tz(end_date,   tz = "UTC")
today      = format(Sys.Date(), "%y%m%d")


# helper functions =============================================================

## standardized save function --------------------------------------------------

save_output = function(obj, name, shareable = FALSE, subfolder = NULL, format = "csv") {
  
  # base folder
  folder = if (shareable) "upload_to_box" else "intermediate_output"
  
  # add subfolder if specified (e.g., "trajectory" or "circadian")
  if (!is.null(subfolder) && shareable) {
    folder = file.path(folder, subfolder)
  }
  
  path = file.path(
    project_location, 
    folder, 
    paste0(name, "_", site_lowercase, "_", today, ".", format)
  )
  
  if (format == "csv") {
    fwrite(obj, path)
  } else if (format == "fst") {
    fst::write_fst(obj, path)
  } else if (format == "parquet") {
    arrow::write_parquet(obj, path)
  } else {
    stop("Unsupported format: ", format)
  }
  
  message("Saved: ", path)
  invisible(path)
}

## safe division and log -------------------------------------------------------

safe_div  = function(num, den) fifelse(!is.na(num) & !is.na(den) & den > 0, num / den, NA_real_)
safe_log2 = function(x) fifelse(!is.na(x) & x > 0, log2(x), NA_real_)

## standardize for trajectory fitting ------------------------------------------

standardize = function(x) {
  (x - fmean(x, na.rm = TRUE)) / fsd(x, na.rm = TRUE)
}


# clif tables needed for project ===============================================

required_tables = 
  c(
    "patient",  
    "hospitalization", 
    "adt", 
    "vitals", 
    "labs", 
    "hospital_diagnosis",
    "respiratory_support", 
    "microbiology_culture",
    "medication_admin_continuous", 
    "medication_admin_intermittent"
  )

## check required tables against available tables ------------------------------

clif_table_filenames = 
  list.files(
    path       = tables_location, 
    pattern    = paste0("^clif_.*\\.", file_type, "$"), 
    full.names = TRUE
  )

clif_table_basenames = 
  basename(clif_table_filenames) |>
  str_remove(paste0("\\.", file_type, "$")) |>
  str_remove("^clif_") |>
  str_remove("(_\\d{4}(_\\d{4})*)$")

table_file_map     = setNames(clif_table_filenames, clif_table_basenames)
missing_tables     = setdiff(required_tables, clif_table_basenames)
required_filenames = table_file_map[required_tables]

if (length(missing_tables) > 0) {
  stop(paste("Error: Missing required tables:", paste(missing_tables, collapse = ", ")))
} else {
  message("All required tables are present.")
}

rm(allowed_files, missing_tables)
gc()


# load tables ==================================================================

if (file_type == "parquet") {
  data_list = lapply(required_filenames, arrow::open_dataset)
} else if (file_type == "csv") {
  data_list = lapply(required_filenames, \(f) arrow::read_csv_arrow(f))
} else if (file_type == "fst") {
  data_list = lapply(required_filenames, \(f) {
    tmp = fst::read_fst(f, as.data.table = TRUE)
    arrow::arrow_table(tmp)
  })
} else {
  stop("Unsupported file format")
}

names(data_list) = names(required_filenames)


# validate table contents ======================================================

validate_table = function(tbl, table_name, req_vars = NULL, req_values = list()) {
  
  if (is.null(req_vars)) req_vars = character(0)
  
  problems     = character()
  missing_vars = setdiff(req_vars, names(tbl))
  
  if (length(missing_vars)) {
    problems = c(problems, sprintf("Missing required vars: %s", paste(missing_vars, collapse = ", ")))
  }
  
  for (var in names(req_values)) {
    if (!var %in% names(tbl)) {
      problems = c(problems, sprintf("Missing '%s' needed for value checks.", var))
      next
    }
    
    if (inherits(tbl, "FileSystemDataset")) {
      present_vals = character(0)
      tryCatch({
        value_counts =
          dplyr::select(tbl, !!rlang::sym(var)) |>
          dplyr::group_by(!!rlang::sym(var)) |>
          dplyr::summarize(n = dplyr::n()) |>
          dplyr::arrange(dplyr::desc(n)) |>
          dplyr::collect()
        if (nrow(value_counts) > 0) {
          present_vals = na.omit(as.character(value_counts[[var]]))
        }
      }, error = function(e) { })
      
      if (length(present_vals) == 0) {
        tryCatch({
          sample_data =
            dplyr::select(tbl, !!rlang::sym(var)) |>
            utils::head(1000) |>
            dplyr::collect()
          if (nrow(sample_data) > 0) {
            present_vals = na.omit(unique(as.character(sample_data[[var]])))
          }
        }, error = function(e) { })
      }
    } else if (!is.data.frame(tbl) && inherits(tbl, "Table")) {
      present_vals = character(0)
      tryCatch({
        distinct_vals = tbl |>
          dplyr::select(!!rlang::sym(var)) |>
          dplyr::distinct() |>
          dplyr::collect()
        if (nrow(distinct_vals) > 0) {
          present_vals = na.omit(as.character(distinct_vals[[var]]))
        }
      }, error = function(e) { })
    } else {
      present_vals = na.omit(unique(as.character(tbl[[var]])))
    }
    
    present_vals  = tolower(trimws(as.character(present_vals)))
    expected_vals = unique(tolower(trimws(req_values[[var]])))
    missing_vals  = setdiff(expected_vals, present_vals)
    
    if (length(missing_vals)) {
      problems = c(
        problems,
        sprintf("Variable '%s' is missing expected values: %s", var, paste(missing_vals, collapse = ", "))
      )
    }
  }
  
  if (length(problems)) {
    return(sprintf("Table '%s':\n- %s", table_name, paste(problems, collapse = "\n- ")))
  }
  
  invisible(NULL)
}

validate_all_tables = function(data_list, validation_specs) {
  all_problems = character()
  
  for (spec in validation_specs) {
    tbl_name = spec$table_name
    if (!tbl_name %in% names(data_list)) {
      all_problems = c(all_problems, sprintf("Table '%s' is missing entirely.", tbl_name))
      next
    }
    
    tbl = data_list[[tbl_name]]
    
    tbl_problems = 
      validate_table(
        tbl        = tbl,
        table_name = tbl_name,
        req_vars   = spec$req_vars,
        req_values = spec$req_values
      )
    
    if (!is.null(tbl_problems)) all_problems = c(all_problems, tbl_problems)
  }
  
  if (length(all_problems)) {
    stop("Validation errors found:\n", paste(all_problems, collapse = "\n\n"), call. = FALSE)
  }
  
  message("Validation passed: all required tables present with needed values.")
}

## validation specs ------------------------------------------------------------

req_vitals = c(
  "heart_rate", "respiratory_rate", "map", "spo2", "temp_c",
  "height_cm", "weight_kg"
)

req_labs = c(
  "crp", "wbc", "creatinine", "platelet_count",
  "lymphocytes_absolute", "ldh", "ferritin",
  "hemoglobin", "neutrophils_absolute"
)

validation_specs = list(
  list(
    table_name = "patient",
    req_vars   = c("patient_id", "race_category", "ethnicity_category", 
                   "sex_category", "death_dttm")
  ),
  list(
    table_name = "hospitalization",
    req_vars   = c("patient_id", "hospitalization_id", "age_at_admission", 
                   "admission_dttm", "discharge_dttm")
  ),
  list(
    table_name = "adt",
    req_vars   = c("hospitalization_id", "location_category", "in_dttm", "out_dttm"),
    req_values = list(location_category = c("icu", "ward"))
  ),
  list(
    table_name = "medication_admin_continuous",
    req_vars   = c("hospitalization_id", "med_group", "med_category"),
    req_values = list(
      med_group    = c("vasoactives"),
      med_category = c("norepinephrine", "vasopressin", "epinephrine", "phenylephrine")
    )
  ),
  list(
    table_name = "medication_admin_intermittent",
    req_vars   = c("hospitalization_id", "med_group"),
    req_values = list(med_group = c("steroid", "analgesia", "car_t"))
  ),
  list(
    table_name = "respiratory_support",
    req_vars   = c("hospitalization_id", "device_category", "recorded_dttm"),
    req_values = list(device_category = c("IMV", "NIPPV", "High Flow NC"))
  ),
  list(
    table_name = "vitals",
    req_vars   = c("hospitalization_id", "vital_category", "vital_value", "recorded_dttm"),
    req_values = list(vital_category = req_vitals)
  ),
  list(
    table_name = "labs",
    req_vars   = c("hospitalization_id", "lab_category", "lab_value", "lab_result_dttm"),
    req_values = list(lab_category = req_labs)
  ),
  list(
    table_name = "microbiology_culture",
    req_vars   = c("hospitalization_id", "result_dttm", "fluid_category", "organism_category")
  ),
  list(
    table_name = "hospital_diagnosis",
    req_vars   = "hospitalization_id"
  )
)

## run validation --------------------------------------------------------------

validate_all_tables(data_list, validation_specs)

rm(validation_specs, req_vitals, req_labs)
rm(clif_table_basenames, clif_table_filenames, required_filenames, required_tables, 
   table_file_map, validate_table, validate_all_tables)
gc()

message("00_setup.R complete.")
