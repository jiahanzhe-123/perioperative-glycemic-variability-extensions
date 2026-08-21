-- 10_build_final_features.sql
-- Final grain: one row per subject_id + hadm_id + stay_id, with the first
-- eligible cardiac surgery ICU stay per subject. All upstream module tables
-- are one-row-per-stay before this join to prevent many-to-many inflation.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.cardiac_surgery_shr_gv_features_v2;

CREATE TABLE mimic_custom.cardiac_surgery_shr_gv_features_v2 AS
WITH icu_span AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           min(i.intime) AS all_icu_first_intime,
           max(i.outtime) AS all_icu_last_outtime,
           (EXTRACT(EPOCH FROM (max(i.outtime) - min(i.intime))) / 86400.0)::double precision AS total_icu_span_days
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_icu.icustays i
      ON i.subject_id = c.subject_id AND i.hadm_id = c.hadm_id
    GROUP BY c.subject_id, c.hadm_id, c.stay_id
)
SELECT
    -- identifiers and index-time fields
    c.subject_id,
    c.hadm_id,
    c.stay_id,
    c.admittime_ts,
    c.dischtime_ts,
    c.surgery_time,
    c.first_icu_intime,
    s.all_icu_last_outtime AS last_icu_outtime,

    -- exposure variables and candidates
    sg.shr AS SHR,
    sg.gv AS GV,
    sg.hba1c_pct,
    sg.glucose_postop_first,
    sg.glucose_adm_first,
    sg.glucose_max_postop_24h,
    sg.glucose_min_postop_24h,
    sg.glucose_median_postop_24h,
    sg.glucose_variability_delta_postop_24h,
    sg.glucose_mean_postop_24h,
    sg.glucose_sd_postop_24h,
    sg.shr_first_postop_24h,
    sg.shr_mean_postop_24h,
    sg.shr_median_postop_24h,
    sg.shr_max_postop_24h,
    sg.shr_mean_postop_48h,
    sg.shr_mean_postop_72h,
    sg.shr_mean_postop_7d,
    sg.shr_adm_first,
    sg.glucose_cv_postop_24h,
    sg.glucose_range_postop_24h,
    sg.glucose_masd_postop_24h,
    sg.eAG_mg_dl,
    sg.glucose_source_mix,
    sg.glucose_unit_inference_flag,
    sg.glucose_outlier_count,
    sg.n_glucose_postop_24h,
    sg.n_glucose_postop_48h,
    sg.n_glucose_postop_72h,
    sg.n_glucose_postop_7d,
    sg.hba1c_time,
    sg.hba1c_days_from_surgery,
    sg.hba1c_source,
    sg.hba1c_unit_inference_flag,
    sg.hba1c_window,
    sg.shr_formula_main,
    sg.shr_missing_reason,
    sg.gv_formula_main,
    sg.gv_missing_reason,

    -- demographics
    c.age_at_admission,
    c.age_ge_75,
    c.gender,
    cb.bmi,
    cb.height_cm,
    cb.weight_kg,
    cb.bmi_source,
    cb.bmi_missing_reason,

    -- postoperative vital signs
    v.hr_mean_postop_24h,
    v.rr_mean_postop_24h,
    v.sbp_mean_postop_24h,
    v.dbp_mean_postop_24h,
    v.temp_c_mean_postop_24h,
    v.hr_max_postop_24h,
    v.rr_max_postop_24h,
    v.sbp_max_postop_24h,
    v.dbp_max_postop_24h,
    v.temp_c_max_postop_24h,
    v.spo2_max_postop_24h,
    v.hr_min_postop_24h,
    v.rr_min_postop_24h,
    v.sbp_min_postop_24h,
    v.dbp_min_postop_24h,
    v.temp_c_min_postop_24h,
    v.spo2_min_postop_24h,
    v.fio2_mean_postop_24h,
    v.fio2_max_postop_24h,
    v.fio2_min_postop_24h,
    v.vitals_source,
    v.vitals_missing_reason,

    -- postoperative first labs
    lp.creat_postop_first,
    lp.bun_postop_first,
    lp.lactate_postop_first,
    lp.ph_postop_first,
    lp.wbc_postop_first,
    lp.hgb_postop_first,
    lp.rbc_postop_first,
    lp.hct_postop_first,
    lp.mcv_postop_first,
    lp.k_postop_first,
    lp.na_postop_first,
    lp.cl_postop_first,
    lp.ca_postop_first,
    lp.anion_gap_postop_first,
    lp.bicarb_postop_first,
    lp.pt_postop_first,
    lp.ptt_postop_first,
    lp.fibrinogen_postop_first,
    lp.platelets_postop_first,

    -- admission first/fallback labs
    la.creat_adm_first,
    la.bun_adm_first,
    la.lactate_adm_first,
    la.ph_adm_first,
    la.wbc_adm_first,
    la.hgb_adm_first,
    la.rbc_adm_first,
    la.hct_adm_first,
    la.mcv_adm_first,
    la.bili_adm_first,
    la.alt_adm_first,
    la.alp_adm_first,
    la.ast_adm_first,
    la.albumin_adm_first,
    la.k_adm_first,
    la.na_adm_first,
    la.cl_adm_first,
    la.ca_adm_first,
    la.anion_gap_adm_first,
    la.bicarb_adm_first,
    la.pt_adm_first,
    la.ptt_adm_first,
    la.fibrinogen_adm_first,
    la.platelets_adm_first,
    la.creat_adm_definition,
    la.albumin_adm_definition,

    -- comorbidities and cardiac disease
    cb.myocardial_infarct,
    cb.congestive_heart_failure,
    cb.peripheral_vascular_disease,
    cb.cerebrovascular_disease,
    cb.diabetes_without_cc,
    cb.diabetes_with_cc,
    cb.diabetes,
    cb.renal_disease,
    cb.malignant_cancer,
    cb.charlson_comorbidity_index,
    cb.hypertension,
    cb.chronic_pulmonary_disease,
    cb.atrial_fibrillation_flag,
    cb.heart_failure_flag,
    cb.valvular_disease_flag,
    cb.ischemic_heart_disease_flag,
    cb.cardiogenic_shock_flag,
    cb.cardiac_arrest_flag,
    cb.pulmonary_hypertension_flag,
    cb.endocarditis_flag,
    cb.prior_pci_or_cabg_flag,
    cb.acute_mi_current_admission_flag,
    cb.obesity_flag,
    cb.dyslipidemia_flag,
    cb.smoking_flag,
    cb.liver_disease_flag,
    cb.severe_liver_disease_flag,
    cb.renal_disease AS ckd_flag,
    cb.esrd_flag,
    cb.aki_flag,
    cb.stroke_flag,
    cb.sepsis_flag,
    cb.infection_flag,

    -- procedure phenotype
    c.procedure_group_main,
    c.procedure_group_all,
    c.procedure_icd_code_main,
    c.procedure_icd_version_main,
    c.procedure_title_main,
    c.n_cardiac_procedures,
    c.cabg_flag,
    c.pci_flag,
    c.open_valve_flag,
    c.transcatheter_valve_flag,
    c.aortic_surgery_flag,
    c.transplant_vad_flag,
    c.congenital_cardiac_flag,
    c.support_only_flag,
    c.open_surgery_flag,
    c.catheter_based_flag,
    c.mixed_open_and_catheter_flag,
    c.exclude_pci_sensitivity_flag,
    c.surgery_time_precision,
    c.surgery_to_icu_hours,
    c.icu_intime_to_surgery_hours,
    c.selected_stay_reason,

    -- treatment/support and severity variables
    t.vasopressor_any_postop_24h AS postop_24hr_vaso_flag,
    t.mechanical_ventilation_postop_24h_flag AS postop_24hr_invasive_vent_flag,
    t.rrt_flag,
    t.total_vent_duration_hours,
    sv.apsiii,
    (EXTRACT(EPOCH FROM (c.last_icu_outtime - c.first_icu_intime)) / 86400.0)::double precision AS icu_los_days,
    (EXTRACT(EPOCH FROM (c.dischtime_ts - c.admittime_ts)) / 86400.0)::double precision AS hosp_los_days,
    (EXTRACT(EPOCH FROM (c.last_icu_outtime - c.first_icu_intime)) / 86400.0)::double precision AS icu_los_days_derived,
    (EXTRACT(EPOCH FROM (c.dischtime_ts - c.admittime_ts)) / 86400.0)::double precision AS hosp_los_days_derived,
    c.adm_to_surgery_days,
    s.total_icu_span_days AS total_icu_span_days_derived,
    s.total_icu_span_days AS total_icu_span,
    t.insulin_any_postop_24h,
    t.insulin_any_postop_48h,
    t.insulin_any_postop_7d,
    t.medication_event_n,
    t.vent_event_n,
    t.rrt_event_n,
    t.insulin_infusion_postop_24h,
    t.insulin_subq_postop_24h,
    t.steroid_any_postop_48h,
    t.steroid_any_postop_7d,
    t.vasopressor_any_postop_24h,
    t.vasopressor_any_postop_48h,
    t.norepinephrine_any_postop_24h,
    t.epinephrine_any_postop_24h,
    t.vasopressin_any_postop_24h,
    t.phenylephrine_any_postop_24h,
    t.dobutamine_any_postop_24h,
    t.milrinone_any_postop_24h,
    t.inotrope_any_postop_24h,
    t.beta_blocker_preop_or_adm_flag,
    t.statin_preop_or_adm_flag,
    t.aspirin_preop_or_adm_flag,
    t.anticoagulant_adm_flag,
    t.antiplatelet_adm_flag,
    cb.aki_postop_48h_flag,
    t.rrt_postop_7d_flag,
    t.mechanical_ventilation_postop_24h_flag,
    t.mechanical_ventilation_postop_48h_flag,
    t.vent_duration_hours,
    t.vent_duration_reliable,
    t.vent_duration_unreliable,
    t.vaso_duration_hours,
    t.high_lactate_postop_24h_flag,
    t.severe_acidosis_postop_24h_flag,
    t.hypoxemia_postop_24h_flag,
    t.hypoglycemia_postop_24h_flag,
    t.hyperglycemia_180_postop_24h_flag,
    t.hyperglycemia_200_postop_24h_flag,
    sv.severity_proxy_score,
    sv.severity_proxy_components,
    sv.apsiii_available,
    sv.apsiii_missing_reason,

    -- outcomes; MIMIC dod is date-only for many post-discharge deaths.
    CASE WHEN COALESCE(c.deathtime, c.dod::timestamp) >= c.surgery_time
              AND COALESCE(c.deathtime, c.dod::timestamp) < c.surgery_time + interval '30 days'
         THEN 1 ELSE 0 END AS postop_30d_death_flag,
    CASE WHEN COALESCE(c.deathtime, c.dod::timestamp) >= c.surgery_time
              AND COALESCE(c.deathtime, c.dod::timestamp) < c.surgery_time + interval '365 days'
         THEN 1 ELSE 0 END AS one_year_death_flag,
    CASE WHEN COALESCE(c.deathtime, c.dod::timestamp) >= c.surgery_time
              AND COALESCE(c.deathtime, c.dod::timestamp) < c.surgery_time + interval '365 days'
         THEN 1 ELSE 0 END AS one_year_mortality,
    CASE WHEN COALESCE(c.deathtime, c.dod::timestamp) IS NULL THEN 365.0
         ELSE LEAST(365.0, GREATEST(0.0, EXTRACT(EPOCH FROM (COALESCE(c.deathtime, c.dod::timestamp) - c.surgery_time)) / 86400.0)) END AS survival_time_days,
    CASE WHEN c.deathtime IS NOT NULL THEN 'mimiciv_hosp.admissions.deathtime'
         WHEN c.dod IS NOT NULL THEN 'mimiciv_hosp.patients.dod'
         ELSE 'no_recorded_death_date' END AS death_time_source,

    -- module-level completeness/QC fields
    (sg.shr IS NULL) AS missing_shr_flag,
    (sg.gv IS NULL) AS missing_gv_flag,
    (sg.hba1c_pct IS NULL) AS missing_hba1c_flag,
    (sg.n_glucose_postop_24h IS NULL OR sg.n_glucose_postop_24h = 0) AS missing_glucose_postop_24h_flag,
    (v.hr_n IS NULL AND v.rr_n IS NULL AND v.sbp_n IS NULL AND v.dbp_n IS NULL) AS missing_vitals_flag,
    (lp.n_valid_postop_first_labs = 0) AS missing_postop_labs_flag,
    (la.n_valid_adm_labs = 0) AS missing_adm_labs_flag,
    (t.treatment_missing_reason IS NOT NULL) AS missing_treatment_flag,
    (cb.bmi IS NULL) AS missing_bmi_flag,
    (sv.apsiii_available = false) AS missing_apsiii_flag,
    jsonb_build_object(
      'shr', sg.shr_missing_reason,
      'gv', sg.gv_missing_reason,
      'hba1c', CASE WHEN sg.hba1c_pct IS NULL THEN 'unavailable_in_selected_HbA1c_window' ELSE NULL END,
      'glucose', CASE WHEN sg.n_glucose_postop_24h IS NULL THEN 'no_valid_24h_glucose' ELSE NULL END,
      'vitals', v.vitals_missing_reason,
      'postop_labs', lp.labs_missing_reason,
      'admission_labs', la.labs_missing_reason,
      'bmi', cb.bmi_missing_reason,
      'treatment', t.treatment_missing_reason,
      'apsiii', sv.apsiii_missing_reason
    ) AS missing_reason_jsonb,

    -- provenance
    c.cohort_source,
    'mimiciv_hosp.procedures_icd + mimiciv_hosp.d_icd_procedures + reviewed codebook'::text AS procedure_source,
    lp.labs_source AS labs_source,
    t.treatment_source,
    cb.comorbidity_source,
    'PostgreSQL local MIMIC-IV; all raw schemas read-only in this project'::text AS project_provenance
FROM mimic_custom.cardiac_surgery_cohort_v2 c
LEFT JOIN mimic_custom.shr_gv_candidates_v2 sg ON sg.stay_id = c.stay_id
LEFT JOIN mimic_custom.vitals_postop_24h_v2 v ON v.stay_id = c.stay_id
LEFT JOIN mimic_custom.labs_postop_first_v2 lp ON lp.stay_id = c.stay_id
LEFT JOIN mimic_custom.labs_adm_first_v2 la ON la.stay_id = c.stay_id
LEFT JOIN mimic_custom.comorbidities_bmi_v2 cb ON cb.stay_id = c.stay_id
LEFT JOIN mimic_custom.treatments_postop_v2 t ON t.stay_id = c.stay_id
LEFT JOIN mimic_custom.severity_v2 sv ON sv.stay_id = c.stay_id
LEFT JOIN icu_span s ON s.stay_id = c.stay_id;

CREATE UNIQUE INDEX cardiac_surgery_features_v2_pk
    ON mimic_custom.cardiac_surgery_shr_gv_features_v2 (subject_id, hadm_id, stay_id);

COMMIT;

\copy (SELECT * FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 ORDER BY subject_id) TO 'outputs/results/cardiac_surgery_shr_gv_features_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
