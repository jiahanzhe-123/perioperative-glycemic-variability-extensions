-- Step 8: 最终分析表 eicu_cardio_validation.cardiac_surgery_glucose_features_v1
BEGIN;
-- 40 条跨源(同分钟同值)重复: 优先保留中心实验室 'glucose'
DROP TABLE IF EXISTS eicu_cardio_validation.glucose_final;
CREATE TABLE eicu_cardio_validation.glucose_final AS
SELECT DISTINCT ON (patientunitstayid, offset_min, value_mgdl)
       patientunitstayid, offset_min, value_mgdl, unit_orig, source_name, source_table
FROM eicu_cardio_validation.glucose_clean
ORDER BY patientunitstayid, offset_min, value_mgdl,
         CASE source_name WHEN 'glucose' THEN 0 ELSE 1 END;

DROP TABLE IF EXISTS eicu_cardio_validation.glucose_features_final;
CREATE TABLE eicu_cardio_validation.glucose_features_final AS
SELECT patientunitstayid,
       count(*) AS glucose_n,
       round(avg(value_mgdl)::numeric, 2) AS glucose_mean_24h,
       round(stddev_samp(value_mgdl)::numeric, 2) AS glucose_sd_24h,
       round(min(value_mgdl)::numeric, 2) AS glucose_min_24h,
       round(max(value_mgdl)::numeric, 2) AS glucose_max_24h,
       round((max(value_mgdl) - min(value_mgdl))::numeric, 2) AS glucose_range_24h,
       round((stddev_samp(value_mgdl) / avg(value_mgdl))::numeric, 4) AS glucose_cv_24h,
       (array_agg(value_mgdl ORDER BY offset_min))[1] AS first_glucose,
       (array_agg(value_mgdl ORDER BY offset_min DESC))[1] AS last_glucose,
       (max(offset_min) - min(offset_min)) AS glucose_measurement_span,
       round((count(*) / nullif((max(offset_min) - min(offset_min))/60.0, 0))::numeric, 3) AS glucose_measurement_density,
       count(*) FILTER (WHERE source_name='bedside glucose') AS n_bedside,
       count(*) FILTER (WHERE source_name='glucose') AS n_lab,
       count(DISTINCT source_name) AS glucose_source_count
FROM eicu_cardio_validation.glucose_final
GROUP BY patientunitstayid;

DROP TABLE IF EXISTS eicu_cardio_validation.cardiac_surgery_glucose_features_v1;
CREATE TABLE eicu_cardio_validation.cardiac_surgery_glucose_features_v1 AS
SELECT
  c.patientunitstayid,
  c.patienthealthsystemstayid,
  c.hospitalid,
  -- 队列
  c.high_specificity_flag,
  c.broad_cohort_flag,
  c.open_extended_flag,
  c.open_core_flag,
  c.surgery_category,
  c.surgery_confidence,
  c.multi_procedure_flag,
  array_to_string(c.surgery_categories, ';') AS surgery_categories_all,
  c.evidence_summary,
  -- 人口学
  CASE WHEN c.age ~ '^[0-9]+$' THEN c.age::int WHEN c.age='> 89' THEN 90 ELSE NULL END AS age,
  (c.age = '> 89') AS age_censored_89plus,
  NULLIF(c.gender,'') AS sex,
  NULLIF(c.ethnicity,'') AS ethnicity,
  CASE WHEN c.admissionheight BETWEEN 100 AND 250 AND c.admissionweight BETWEEN 25 AND 400
       THEN round((c.admissionweight / ((c.admissionheight/100.0)^2))::numeric, 1) END AS bmi,
  -- 严重程度与入科
  v.apachescore, v.aps, v.apacheversion, v.predicted_hosp_mortality,
  v.elective_surgery,
  NULLIF(c.unitadmitsource,'') AS admission_source,
  c.unittype AS icu_unit_type,
  -- 血糖
  g.glucose_n, g.glucose_mean_24h, g.glucose_sd_24h, g.glucose_cv_24h,
  g.glucose_min_24h, g.glucose_max_24h, g.glucose_range_24h,
  g.first_glucose, g.last_glucose,
  g.glucose_measurement_span, g.glucose_measurement_density,
  g.n_bedside, g.n_lab, g.glucose_source_count,
  -- HbA1c(全库仅 9 条且均不在队列内 -> 不可得)
  NULL::numeric AS hba1c,
  NULL::int AS hba1c_offset,
  NULL::numeric AS eag,
  NULL::numeric AS shr,
  false AS hba1c_available,
  -- 实验室(0-24h 首次)
  v.lactate, v.creatinine, v.hemoglobin, v.wbc, v.platelets, v.albumin,
  -- 结局
  o.hospital_mortality, o.icu_mortality,
  round(o.icu_los_days::numeric, 2) AS icu_los_days,
  round(o.hospital_los_days::numeric, 2) AS hospital_los_days,
  -- 基础疾病与治疗
  (v.ph_diabetes OR v.apache_diabetes) AS diabetes,
  v.ph_renal AS renal_disease,
  v.ph_chf AS heart_failure,
  v.ph_pulmonary AS chronic_pulmonary_disease,
  v.ph_liver AS liver_disease,
  v.ph_hypertension AS hypertension,
  v.ph_afib AS atrial_fibrillation,
  v.mechanical_ventilation, v.vasopressor, v.insulin_infusion
FROM eicu_cardio_validation.cohort_base c
LEFT JOIN eicu_cardio_validation.glucose_features_final g USING (patientunitstayid)
LEFT JOIN eicu_cardio_validation.covariates v USING (patientunitstayid)
LEFT JOIN eicu_cardio_validation.outcomes o USING (patientunitstayid)
WHERE c.broad_cohort_flag AND c.rn_patient=1;
CREATE UNIQUE INDEX ON eicu_cardio_validation.cardiac_surgery_glucose_features_v1(patientunitstayid);
COMMIT;

SELECT count(*) AS n_final,
       count(*) FILTER (WHERE glucose_n>=2) AS gv_analysis_n,
       count(*) FILTER (WHERE high_specificity_flag AND glucose_n>=2) AS highspec_gv_n,
       count(*) FILTER (WHERE open_core_flag AND glucose_n>=2) AS open_core_gv_n,
       count(*) FILTER (WHERE hospital_mortality) AS hosp_deaths,
       count(*) FILTER (WHERE icu_mortality) AS icu_deaths
FROM eicu_cardio_validation.cardiac_surgery_glucose_features_v1;
\o /tmp/final_table.csv
SELECT * FROM eicu_cardio_validation.cardiac_surgery_glucose_features_v1 ORDER BY patientunitstayid;
\o
