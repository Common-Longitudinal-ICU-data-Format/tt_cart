# CLIF: CAR-T Temperature Trajectory & Circadian Analysis

**CLIF VERSION 2.1**

## Overview

This project provides an R-based pipeline for analyzing temperature trajectories and circadian infusion timing in patients undergoing chimeric antigen receptor T-cell (CAR-T) therapy. The pipeline supports two related analyses:

1. **Trajectory Analysis**: Testing whether post-infusion temperature trajectories predict CAR-T complications, including cytokine release syndrome, infections, and ICU care.

2. **Circadian Analysis**: Testing whether time of day of CAR-T infusion predicts outcomes, and whether temperature trajectory mediates this relationship.

Both analyses are designed for federated multi-site analysis. Each site runs the pipeline locally and exports aggregated summary statistics (no patient-level data) for pooled meta-analysis.

---

## Quick Start

### 1. Clone the Repository and Open the R Project

```bash
git clone [repository-url]
```

**Important**: Open the `.Rproj` file in RStudio before running any scripts. The pipeline uses `rprojroot` to locate the config file.

### 2. Configure `config/config.yaml`

```yaml
# Your institution name (used in output filenames)
institution: "your_site"

# Your local time zone
# Options: 'America/New_York', 'America/Chicago', 'America/Los_Angeles', etc.
time_zone: "America/Chicago"

# Path to your CLIF 2.1 tables
tables_location: "/path/to/your/clif/tables"

# File format of your CLIF tables
file_type: "parquet"  # Options: parquet, csv, fst

# Project output location (leave blank to use project root)
project_location: ""
```

### 3. Run the Pipeline

Run scripts in order from within the RStudio project. Scripts must run sequentially in a single R session (objects are passed in memory).

```r
source("code/00_setup.R")               # packages, validation, helpers
source("code/01_cohort.R")              # cohort identification + exclusions
source("code/02_traj.R")                # temperature trajectory assignment
source("code/03_covariates.R")          # demographics, diagnoses, comorbidities
source("code/04_outcomes.R")            # death, ICU, labs, ICAHT, IEC-HLH, EASIX, meds
source("code/05_abxs_cultures.R")       # infections + antimicrobial exposure
source("code/06_export.R")              # final assembly + Table 1 / Table 2 exports

# Analysis scripts (run after 00-06)
source("code/07_analysis_traj.R")      # Trajectory paper
source("code/08_analysis_circadian.R") # Circadian paper
```
---

## Required CLIF Tables and Fields

This project requires the following CLIF 2.1 tables. Ensure all listed fields are present.

### Core Tables

| Table | Required Fields |
|-------|-----------------|
| **patient** | `patient_id`, `race_category`, `ethnicity_category`, `sex_category`, `death_dttm` |
| **hospitalization** | `patient_id`, `hospitalization_id`, `age_at_admission`, `admission_dttm`, `discharge_dttm`, `discharge_category` |
| **adt** | `hospitalization_id`, `location_category`, `in_dttm`, `out_dttm` |
| **vitals** | `hospitalization_id`, `vital_category`, `vital_value`, `recorded_dttm` |
| **labs** | `hospitalization_id`, `lab_category`, `lab_value`, `lab_value_numeric`, `lab_result_dttm` |
| **medication_admin_continuous** | `hospitalization_id`, `med_group`, `med_category`, `admin_dttm` |
| **medication_admin_intermittent** | `hospitalization_id`, `med_group`, `med_category`, `admin_dttm`, `med_dose`, `med_dose_unit` |
| **respiratory_support** | `hospitalization_id`, `device_category`, `recorded_dttm` |
| **microbiology_culture** | `hospitalization_id`, `organism_id`, `collect_dttm`, `result_dttm`, `fluid_category`, `organism_category`, `organism_group` |
| **microbiology_nonculture** | `hospitalization_id`, `collect_dttm`, `result_dttm`, `method_category`, `result_category`, `organism_category`, `organism_group` |
| **hospital_diagnosis** | `hospitalization_id`, `diagnosis_code`|

### Required Values Within Tables

**vitals.vital_category**:
- `temp_c` (required for trajectory assignment)
- `heart_rate`, `respiratory_rate`, `map`, `spo2`, `weight_kg`, `height_cm`

**labs.lab_category**:
- `creatinine`, `platelet_count`, `lymphocytes_absolute`, `wbc` (required for EASIX scores)
- `ldh`, `crp`, `ferritin` (required for baseline labs and EASIX scores)
- `hemoglobin`, `neutrophils_absolute` (required for hematotoxicity/ICAHT)
- `ast`, `alt`, `bilirubin_total` (required for IEC-HLH organ-dysfunction criteria)

**medication_admin_continuous.med_group**:
- `vasoactives` (with med_category: norepinephrine, vasopressin, epinephrine, phenylephrine, etc.)

**medication_admin_intermittent.med_group**:
- `car_t` (required for cohort identification)
- `analgesia`, `steroid` (required for outcomes)

Note: `med_group` value used downstream but not strictly validated: `CMS_sepsis_qualifying_antibiotics` (used in 05_abxs_cultures.R)
      `med_category` values used downstream but not strictly validated: `tocilizumab`, `anakinra` (used for toxicity outcomes)  
      
**respiratory_support.device_category**:
- `IMV`, `NIPPV`, `High Flow NC`, etc.

**adt.location_category**:
- `icu`, `ward`

---

## Cohort Definition

### Inclusion Criteria
- Adult patients (≥18 years) admitted to the hospital between 2018–2025 (can include narrower range if certain years unavailable)
- Received CD19 or BCMA targeted CAR-T therapy infusion (identified via `medication_admin_intermittent` where `med_group = 'car_t'`), restricted to recognized CAR-T agents:
`axicabtagene ciloleucel`, `brexucabtagene autoleucel`, `ciltacabtagene autoleucel`, `idecabtagene vicleucel`, `lisocabtagene maraleucel`, `tisagenlecleucel`
- At least 24 hours of in-hospital follow-up and recorded temperatures after infusion

### Exclusion Criteria
- Patients discharged within 24 hours of CAR-T infusion
- Patients who died within 24 hours of CAR-T infusion
- Insufficient temperature data for trajectory assignment

### Hospitalization Linking
Contiguous hospitalization_id records for the same patient separated by ≤ 6 hours (discharge-to-readmission gap) are linked into a single joined_hosp_id for analysis. This handles transfers and brief inter-encounter gaps without losing post-infusion follow-up.

---
## Expected Outputs

After running the full pipeline, you will find outputs in:

### `upload_to_box/trajectory/` (Shareable - for Trajectory Paper)

| File Pattern | Description |
|--------------|-------------|
| `exclusion_flow-{site}.csv` | Cohort exclusion flow counts |
| `concordance_24h_72h_summary_{site}_{date}.csv` | 24h-vs-72h trajectory concordance summary |
| `concordance_24h_72h_crosstab_{site}_{date}.csv` | 24h × 72h trajectory cross-tabulation |
| `concordance_24h_72h_discordant_detail_{site}_{date}.csv` | Discordant-case detail (when n > 0) |
| `fig_temp_by_trajectory_{site}_{date}.csv` | Hourly temperature aggregates by trajectory |
| `fig_temp_by_product_{site}_{date}.csv` | Hourly temperature aggregates by CAR-T product |
| `lab_trends_mean_by_trajectory_{site}__{date}.csv`| Daily mean lab aggregates by trajectory (days 0–30) |
| `lab_trends_nadir_by_trajectory_{site}_{date}.csv` | Daily nadir lab aggregates by trajectory (days 0–30) |
| `meds_daily_by_trajectory_{site}_{date}.csv` | Daily tocilizumab / anakinra / steroid (dex-eq) by trajectory |
| `infection_daily_by_trajectory_{site}_{date}.csv` | Daily counts of bacterial / viral / fungal infections by trajectory |
| `bugs_by_trajectory_{site}_{date}.csv` | Distinct organisms cultured per trajectory (0–30 d) |
| `antimicrobials_baseline_by_traj_{site}_{date}.csv` | Baseline antimicrobial exposure by trajectory |
| `antimicrobials_72h_by_traj_{site}_{date}.csv` | Daily antimicrobial flags days 0–2 by trajectory |
| `antimicrobials_3_30d_by_traj_{site}_{date}.csv` | Daily antimicrobial flags days 3–30 by trajectory |
| `diagnostics-overall-{site}.csv` | Site-level summary statistics |
| `diagnostics-by_trajectory-{site}.csv` | N and outcome counts by trajectory |
| `site_summary_{site}_{date}.csv` | Site-level diagnostics and missingness|
| `table1-continuous-{site}.csv` | Continuous covariates by trajectory (n, sum, sum_sq) |
| `table1-categorical-{site}.csv` | Categorical covariates by trajectory (counts) |
| `table2_acute_{site}_{date}.csv` | Acute binary outcomes by trajectory |
| `table2_medication_outcomes_{site}_{date}.csv` | Medication dose / count outcomes by trajectory |
| `table2_nadir_{site}_{date}.csv` | ANC / Plt / Hgb / ALC nadir aggregates by trajectory |
| `table2_grade_distributions_{site}_{date}.csv` | ICAHT, neutropenia, anemia, thrombocytopenia grade distributions |
| `table2_time_{site}_{date}.csv` | LOS, hours-to-toci, hours-to-anakinra by trajectory |
| `table2_infection_outcomes_{site}_{date}.csv` | Infection outcomes by trajectory (baseline / 72 h / 3–30 d windows) |
| `survival-q30d-{site}.csv` | Survival data at 30-day intervals for KM curves |
| `models-logistic_coefs-{site}.csv` | Logistic regression coefficients (main analysis) |
| `models-rmst_coefs-{site}.csv` | RMST pseudo-regression coefficients |
| `subgroup-*-{site}.csv` | Subgroup analyses (CD19, BCMA) |
| `sensitivity-*-{site}.csv` | Sensitivity analyses (borderline reassignment, HFR collapse) |

### `upload_to_box/circadian/` (Shareable - for Circadian Paper)

| File Pattern | Description |
|--------------|-------------|
| `diagnostics-overall-{site}.csv` | Site-level summary statistics |
| `diagnostics-by_hour-{site}.csv` | N and outcomes by infusion hour |
| `table1-continuous-{site}.csv` | Continuous covariates by early/late infusion |
| `table1-categorical-{site}.csv` | Categorical covariates by early/late infusion |
| `table2-binary-{site}.csv` | Binary outcomes by early/late infusion |
| `survival-earlylate-{site}.csv` | Survival data for early vs late KM curves |
| `models-logistic_coefs-{site}.csv` | Logistic regression coefficients (time → outcomes) |
| `models-rmst_coefs-{site}.csv` | RMST coefficients (time → survival) |
| `mediation-logistic_coefs-{site}.csv` | Mediation analysis (total vs direct effects) |
| `mediation-rmst_coefs-{site}.csv` | RMST mediation analysis |
| `subgroup-*-{site}.csv` | Subgroup analyses (CD19, BCMA, male, female) |
| `sensitivity-*-{site}.csv` | Sensitivity analyses (cutpoint, HFR collapse) |

### `intermediate_output/` (Not Shareable - Local Only)

| File | Description |
|------|-------------|
| `cohort_analytic_{site}_{date}.csv` | Full patient-level analytic dataset |
| `model_data_{site}_{date}.csv` | Individual-level data for local modeling |
| `rmst_pseudo_365_{site}_{date}.csv` | Patient-level RMST pseudo-observations |

---

## Temperature Trajectories

Patients are assigned to one of four temperature trajectory groups based on hourly temperature data during the first 24 and 72 hours post-infusion:
| Trajectory | Description |
|------------|-------------|
| **NT** | `Normothermic`, Temperatures remain near normal throughout |
| **HT** | `Hypothermic`, Temperatures trend below normal |
| **HFR** | `Hyperthermic Fast-Resolver`, Early fever that resolves quickly |
| **HSR** | `Hyperthermic Slow-Resolver`, Prolonged fever or delayed resolution |

Assignment is performed in 02_traj.R by computing the sum of squared errors (SSE) between each patient's standardized hourly temperatures and four pre-specified polynomial reference curves (one per trajectory). Each patient is assigned to the trajectory with the lowest SSE. The second-best trajectory and a relative margin (margin_rel = (sse_second − sse_best) / sse_best) are retained; patients with margin_rel < 0.10 are flagged as borderline for sensitivity analyses.

02_traj.R also performs a 24h-vs-72h concordance analysis, repeating the SSE-based assignment using only the first 24 hours of data and comparing to the 72 h assignment. Pooled outputs include the cross-tabulation, the discordant-case detail table, and a site-level concordance summary.

Note: you may have no patients assigned to HFR group and that is not unexpected

---

## Key Outcomes

### Primary
- 90-day mortality
- 365-day mortality (RMST at τ=365 days)

### Secondary
- ICU admission within 30 days
- Invasive mechanical ventilation within 30 days
- Vasoactive support within 30 days
- Index hospital length of stay
- Tocilizumab administration (CRS proxy)
- Anakinra administration (ICANS proxy)
- Corticosteroid administration (toxicity proxy)
- Severe hematotoxicity (ICAHT grade ≥3)
- Hyperferritinemia
- IEC-HLH composite
- Infections within 30 days

---

## Analysis Notes

### Federated Analysis Approach
- Each site runs the pipeline locally
- Only aggregated summary statistics are shared (no patient-level data)
- Coefficient files include `estimate`, `std.error`, `conf.low`, `conf.high` for meta-analysis
- Meta-analysis uses inverse-variance weighted pooling via `metafor::rma()`

### Sensitivity Analyses
1. **Borderline Reassignment**: Patients with borderline trajectory fit reassigned to second-best trajectory
2. **HFR Collapse**: HFR trajectory collapsed into NT (due to typically small N)

### Subgroup Analyses
- CAR-T target: CD19 vs BCMA
- Sex: Male vs Female (circadian paper only)

---

### Package Dependencies
00_setup.R will attempt to install missing packages automatically. Key dependencies:

Data manipulation: `tidytable`, `collapse`, `data.table`, `tidyverse`, `arrow`, `janitor`
Survival / pseudo-observations: `pseudo`
Comorbidity scoring: `comorbidity`
Model output: `broom`
Configuration / paths: `yaml`, `rprojroot`, `here`
I/O: `readr`, `readxl`, `fst`, `arrow`
System / threading: `ps`, `lubridate`, `stringr`, `rlang`

The setup script auto-tunes thread counts based on detected cores and available RAM.

---

## Troubleshooting

### "Cannot find config directory"
- Ensure you opened the `.Rproj` file before running scripts
- Check that `config/config.yaml` exists in your project folder

### "site_lowercase is empty"
- Add `institution: "your_site"` to your config.yaml

### Models fail with convergence errors
- Common with small subgroups (e.g., BCMA)
- The pipeline uses Wald confidence intervals as fallback
- Results will still be generated

### Missing CLIF tables
- Check that all required tables exist in your `tables_location`
- Verify table names match CLIF 2.1 conventions
- Run `00_setup.R` to see validation output

### "Validation errors found"
- The setup script will list, per table, which required variables or category values are missing. Fix the underlying CLIF table or mapping and re-run 00_setup.R.

---

## Contact

For questions about this analysis, contact the study coordinating center.

For questions about CLIF, see: https://clif-consortium.github.io/
