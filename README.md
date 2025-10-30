#CLIF: CAR-T Temperature Trajectory

## CLIF VERSION 2.1

## Objective

This project provides an R-based pipeline for analyzing temperature trajectories in patients undergoing chimeric antigen receptor T-cell (CAR-T) therapy. The goal is to test the extent to which these trajectories predict CAR-T complications, including cytokine release syndrome, infections and ICU care.

## Required CLIF tables and fields

This project uses standardized CLIF tables:

Tables included:
1. **patient:** patient_id, race_category, ethnicity_category, sex_category, death_dttm
2. **hospitalization:** patient_id, hospitalization_id, age_at_admission, admission_dttm, discharge_dttm 
3. **adt:** hospitalization_id, location_category, in_dttm, out_dttm
   location_category = icu, ward
4. **vitals:** hospitalization_id, vital_category, vital_value, recorded_dttm
   vital_category = heart_rate, respiratory_rate, map, spo2, temp_c, height_cm, weight_kg
5. **labs:** hospitalization_id, lab_category, lab_value, lab_result_dttm
   lab_category = crp, wbc, creatinine, platelet_count, lymphocytes_absolute, ldh, ferritin
6. **medication_admin_continuous:** hospitalization_id, med_group, med_category
   med_group = vasoactives
   med_category = norepinephrine, vasopressin, epinephrine, phenylephrine
7. **respiratory_support:** hospitalization_id, device_category, recorded_dttm
   device_category = IMV, NIPPV, High Flow NC
8. **medication_admin_intermittent:** hospitalization_id, med_group
   med_group = steroid, analgesia, car_t
14. **microbiology_culture:** hospitalization_id, result_dttm, fluid_category, organism_category

## Cohort identification
Inclusion criteria:
Adult patients admitted to the hospital between 2018–2024 for CD19 or BCMA targeted CAR-T therapy infusion
Documented normal temperature (≥35°C and <38°C) within 4 hours prior to CAR-T infusion
At least 24 hours of recorded temperatures following infusion for trajectory assignment

Exclusion criteria: 
Patients discharged within 24 hours of CAR-T infusion
Patients who died within 24 hours of CAR-T infusion (due to inadequate temperature data)

## Expected Results
The initial results will be a database of CLIF tables: 
1. counts for exclusion flow diagram (proj_output/excl_summary)
2. summary statistics for vitals (proj_output/vitals_72h_summary)
3. hourly temperature summary by trajectory (proj_output/temp_hourly_summary)
4. summary statistics for continuous outcome and demographic variables (proj_output/continuous_all; proj_output/continuous_by_traj)
5. summary statistics for categorical outcome and demographic variables (proj_output/categorical_all; proj_output/categorical_by_traj)

Summary tables will be available in the output/final directory.

## Instructions for running the project

Quick start checklist

1. Populate config/config_clif_oncrisk.yaml.

2. Confirm clif_sites.csv contains your site_lowercase.

3. Place CLIF extracts in clif_data_location with required table names.

4. Run 00 → 01 → 02 → 03 in order.


