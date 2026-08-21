-- 11_build_qc_report.sql
-- Produces machine-readable QC tables. scripts/render_qc_report.py converts
-- these tables into a Markdown report without inventing unrun results.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.qc_cohort_flow_v2;
DROP TABLE IF EXISTS mimic_custom.qc_table_counts_v2;
DROP TABLE IF EXISTS mimic_custom.qc_feature_missingness_v2;
DROP TABLE IF EXISTS mimic_custom.qc_procedure_distribution_v2;
DROP TABLE IF EXISTS mimic_custom.qc_outcome_summary_v2;
DROP TABLE IF EXISTS mimic_custom.qc_shr_gv_distribution_v2;
DROP TABLE IF EXISTS mimic_custom.qc_glucose_density_v2;
DROP TABLE IF EXISTS mimic_custom.qc_duplicate_keys_v2;
DROP TABLE IF EXISTS mimic_custom.qc_review_counts_v2;
DROP TABLE IF EXISTS mimic_custom.qc_treatment_reliability_v2;

CREATE TABLE mimic_custom.qc_cohort_flow_v2 AS
SELECT * FROM (VALUES
  ('procedure_codebook_rows', (SELECT count(*)::bigint FROM mimic_custom.cardiac_surgery_codebook_v2)),
  ('all_candidate_procedure_rows', (SELECT count(*)::bigint FROM mimic_custom.cardiac_surgery_cohort_all_candidates_v2)),
  ('all_candidate_subjects', (SELECT count(DISTINCT subject_id)::bigint FROM mimic_custom.cardiac_surgery_cohort_all_candidates_v2)),
  ('adult_candidate_rows', (SELECT count(*)::bigint FROM mimic_custom.cardiac_surgery_cohort_all_candidates_v2 WHERE adult_flag)),
  ('candidate_rows_with_icu_match', (SELECT count(*)::bigint FROM mimic_custom.cardiac_surgery_cohort_all_candidates_v2 WHERE has_icu_match)),
  ('adult_candidate_rows_with_icu', (SELECT count(*)::bigint FROM mimic_custom.cardiac_surgery_cohort_all_candidates_v2 WHERE adult_flag AND has_icu_match)),
  ('final_first_stay_rows', (SELECT count(*)::bigint FROM mimic_custom.cardiac_surgery_cohort_v2)),
  ('final_first_stay_subjects', (SELECT count(DISTINCT subject_id)::bigint FROM mimic_custom.cardiac_surgery_cohort_v2))
) AS x(step_name, n);

CREATE TABLE mimic_custom.qc_table_counts_v2 AS
SELECT table_name,
       count(*)::bigint AS n
FROM (
  SELECT 'cardiac_surgery_cohort_all_candidates_v2'::text AS table_name, subject_id FROM mimic_custom.cardiac_surgery_cohort_all_candidates_v2
  UNION ALL SELECT 'cardiac_surgery_cohort_v2', subject_id FROM mimic_custom.cardiac_surgery_cohort_v2
  UNION ALL SELECT 'vitals_postop_24h_v2', subject_id FROM mimic_custom.vitals_postop_24h_v2
  UNION ALL SELECT 'labs_postop_first_v2', subject_id FROM mimic_custom.labs_postop_first_v2
  UNION ALL SELECT 'labs_adm_first_v2', subject_id FROM mimic_custom.labs_adm_first_v2
  UNION ALL SELECT 'glucose_summary_v2', subject_id FROM mimic_custom.glucose_summary_v2
  UNION ALL SELECT 'hba1c_baseline_v2', subject_id FROM mimic_custom.hba1c_baseline_v2
  UNION ALL SELECT 'comorbidities_bmi_v2', subject_id FROM mimic_custom.comorbidities_bmi_v2
  UNION ALL SELECT 'treatments_postop_v2', subject_id FROM mimic_custom.treatments_postop_v2
  UNION ALL SELECT 'severity_v2', subject_id FROM mimic_custom.severity_v2
  UNION ALL SELECT 'shr_gv_candidates_v2', subject_id FROM mimic_custom.shr_gv_candidates_v2
  UNION ALL SELECT 'cardiac_surgery_shr_gv_features_v2', subject_id FROM mimic_custom.cardiac_surgery_shr_gv_features_v2
) t
GROUP BY table_name
ORDER BY table_name;

CREATE TABLE mimic_custom.qc_feature_missingness_v2 AS
WITH totals AS (SELECT count(*)::bigint AS n FROM mimic_custom.cardiac_surgery_shr_gv_features_v2), metrics AS (
    SELECT * FROM (VALUES
      ('shr','shr IS NULL'),('gv','gv IS NULL'),('hba1c_pct','hba1c_pct IS NULL'),
      ('glucose_mean_postop_24h','glucose_mean_postop_24h IS NULL'),
      ('hr_mean_postop_24h','hr_mean_postop_24h IS NULL'),('creat_postop_first','creat_postop_first IS NULL'),
      ('creat_adm_first','creat_adm_first IS NULL'),('bmi','bmi IS NULL'),
      ('charlson_comorbidity_index','charlson_comorbidity_index IS NULL'),
      ('severity_proxy_score','severity_proxy_score IS NULL'),('apsiii','apsiii IS NULL'),
      ('insulin_any_postop_48h','insulin_any_postop_48h IS NULL'),
      ('vasopressor_any_postop_24h','vasopressor_any_postop_24h IS NULL'),
      ('mechanical_ventilation_postop_24h_flag','mechanical_ventilation_postop_24h_flag IS NULL'),
      ('rrt_postop_7d_flag','rrt_postop_7d_flag IS NULL')
    ) AS m(variable_name, null_predicate)
)
SELECT m.variable_name,
       t.n AS total_n,
       CASE m.variable_name
         WHEN 'shr' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE shr IS NULL)
         WHEN 'gv' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE gv IS NULL)
         WHEN 'hba1c_pct' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE hba1c_pct IS NULL)
         WHEN 'glucose_mean_postop_24h' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE glucose_mean_postop_24h IS NULL)
         WHEN 'hr_mean_postop_24h' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE hr_mean_postop_24h IS NULL)
         WHEN 'creat_postop_first' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE creat_postop_first IS NULL)
         WHEN 'creat_adm_first' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE creat_adm_first IS NULL)
         WHEN 'bmi' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE bmi IS NULL)
         WHEN 'charlson_comorbidity_index' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE charlson_comorbidity_index IS NULL)
         WHEN 'severity_proxy_score' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE severity_proxy_score IS NULL)
         WHEN 'apsiii' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE apsiii IS NULL)
         WHEN 'insulin_any_postop_48h' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE insulin_any_postop_48h IS NULL)
         WHEN 'vasopressor_any_postop_24h' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE vasopressor_any_postop_24h IS NULL)
         WHEN 'mechanical_ventilation_postop_24h_flag' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE mechanical_ventilation_postop_24h_flag IS NULL)
         WHEN 'rrt_postop_7d_flag' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE rrt_postop_7d_flag IS NULL)
       END::bigint AS missing_n,
       CASE WHEN t.n = 0 THEN NULL ELSE (CASE m.variable_name
         WHEN 'shr' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE shr IS NULL)
         WHEN 'gv' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE gv IS NULL)
         WHEN 'hba1c_pct' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE hba1c_pct IS NULL)
         WHEN 'glucose_mean_postop_24h' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE glucose_mean_postop_24h IS NULL)
         WHEN 'hr_mean_postop_24h' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE hr_mean_postop_24h IS NULL)
         WHEN 'creat_postop_first' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE creat_postop_first IS NULL)
         WHEN 'creat_adm_first' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE creat_adm_first IS NULL)
         WHEN 'bmi' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE bmi IS NULL)
         WHEN 'charlson_comorbidity_index' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE charlson_comorbidity_index IS NULL)
         WHEN 'severity_proxy_score' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE severity_proxy_score IS NULL)
         WHEN 'apsiii' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE apsiii IS NULL)
         WHEN 'insulin_any_postop_48h' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE insulin_any_postop_48h IS NULL)
         WHEN 'vasopressor_any_postop_24h' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE vasopressor_any_postop_24h IS NULL)
         WHEN 'mechanical_ventilation_postop_24h_flag' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE mechanical_ventilation_postop_24h_flag IS NULL)
         WHEN 'rrt_postop_7d_flag' THEN (SELECT count(*) FROM mimic_custom.cardiac_surgery_shr_gv_features_v2 WHERE rrt_postop_7d_flag IS NULL)
       END)::numeric / t.n END AS missing_rate
FROM metrics m CROSS JOIN totals t;

CREATE TABLE mimic_custom.qc_procedure_distribution_v2 AS
SELECT procedure_group_main,
       count(*)::bigint AS n,
       avg(CASE WHEN open_surgery_flag THEN 1.0 ELSE 0.0 END)::numeric AS open_share,
       avg(CASE WHEN catheter_based_flag THEN 1.0 ELSE 0.0 END)::numeric AS catheter_share,
       avg(CASE WHEN exclude_pci_sensitivity_flag THEN 1.0 ELSE 0.0 END)::numeric AS pure_pci_share
FROM mimic_custom.cardiac_surgery_shr_gv_features_v2
GROUP BY procedure_group_main
ORDER BY n DESC;

CREATE TABLE mimic_custom.qc_outcome_summary_v2 AS
SELECT count(*)::bigint AS n,
       sum(postop_30d_death_flag)::bigint AS postop_30d_deaths,
       avg(postop_30d_death_flag::numeric) AS postop_30d_mortality,
       sum(one_year_death_flag)::bigint AS one_year_deaths,
       avg(one_year_death_flag::numeric) AS one_year_mortality
FROM mimic_custom.cardiac_surgery_shr_gv_features_v2;

CREATE TABLE mimic_custom.qc_shr_gv_distribution_v2 AS
SELECT 'shr'::text AS variable_name, count(shr)::bigint AS nonmissing_n,
       count(*) FILTER (WHERE shr IS NULL)::bigint AS missing_n,
       min(shr) AS min_value, percentile_cont(0.01) WITHIN GROUP (ORDER BY shr) AS p01,
       percentile_cont(0.25) WITHIN GROUP (ORDER BY shr) AS p25, percentile_cont(0.5) WITHIN GROUP (ORDER BY shr) AS median,
       percentile_cont(0.75) WITHIN GROUP (ORDER BY shr) AS p75, percentile_cont(0.99) WITHIN GROUP (ORDER BY shr) AS p99, max(shr) AS max_value
FROM mimic_custom.cardiac_surgery_shr_gv_features_v2
UNION ALL
SELECT 'gv', count(gv), count(*) FILTER (WHERE gv IS NULL),
       min(gv), percentile_cont(0.01) WITHIN GROUP (ORDER BY gv),
       percentile_cont(0.25) WITHIN GROUP (ORDER BY gv), percentile_cont(0.5) WITHIN GROUP (ORDER BY gv),
       percentile_cont(0.75) WITHIN GROUP (ORDER BY gv), percentile_cont(0.99) WITHIN GROUP (ORDER BY gv), max(gv)
FROM mimic_custom.cardiac_surgery_shr_gv_features_v2;

CREATE TABLE mimic_custom.qc_glucose_density_v2 AS
SELECT count(*)::bigint AS n,
       avg(n_glucose_postop_24h)::numeric AS mean_n_24h,
       percentile_cont(0.5) WITHIN GROUP (ORDER BY n_glucose_postop_24h) AS median_n_24h,
       percentile_cont(0.9) WITHIN GROUP (ORDER BY n_glucose_postop_24h) AS p90_n_24h,
       avg(CASE WHEN n_glucose_postop_24h >= 4 THEN 1.0 ELSE 0.0 END)::numeric AS share_at_least_4,
       avg(CASE WHEN n_glucose_postop_24h >= 6 THEN 1.0 ELSE 0.0 END)::numeric AS share_at_least_6
FROM mimic_custom.cardiac_surgery_shr_gv_features_v2;

CREATE TABLE mimic_custom.qc_duplicate_keys_v2 AS
SELECT count(*)::bigint AS n_rows,
       count(DISTINCT subject_id)::bigint AS n_subjects,
       count(DISTINCT hadm_id)::bigint AS n_hadm,
       count(DISTINCT stay_id)::bigint AS n_stay,
       (count(*) - count(DISTINCT (subject_id,hadm_id,stay_id)))::bigint AS duplicate_composite_key_rows
FROM mimic_custom.cardiac_surgery_shr_gv_features_v2;

CREATE TABLE mimic_custom.qc_review_counts_v2 AS
SELECT 'procedure_codebook' AS review_artifact, count(*)::bigint AS n_rows,
       count(*) FILTER (WHERE review_status = 'manual_review')::bigint AS manual_review_n,
       count(*) FILTER (WHERE include_primary = false)::bigint AS excluded_n
FROM mimic_custom.cardiac_surgery_codebook_v2
UNION ALL
SELECT 'itemid_review', count(*), count(*) FILTER (WHERE needs_manual_review), count(*) FILTER (WHERE include_candidate = false)
FROM mimic_custom.itemid_review_v2
UNION ALL
SELECT 'comorbidity_codebook', count(*), count(*) FILTER (WHERE needs_manual_review), count(*) FILTER (WHERE include_candidate = false)
FROM mimic_custom.comorbidity_codebook_v2
UNION ALL
SELECT 'drug_keyword_review', count(*), count(*) FILTER (WHERE needs_manual_review), count(*) FILTER (WHERE include_candidate = false)
FROM mimic_custom.drug_keyword_review_v2;

CREATE TABLE mimic_custom.qc_treatment_reliability_v2 AS
SELECT count(*)::bigint AS n,
       sum(CASE WHEN vent_duration_reliable THEN 1 ELSE 0 END)::bigint AS vent_duration_reliable_n,
       sum(CASE WHEN vent_duration_reliable = false THEN 1 ELSE 0 END)::bigint AS vent_duration_unreliable_n,
       sum(CASE WHEN apsiii_available THEN 1 ELSE 0 END)::bigint AS apsiii_available_n,
       sum(CASE WHEN apsiii_available = false THEN 1 ELSE 0 END)::bigint AS apsiii_unavailable_n
FROM mimic_custom.cardiac_surgery_shr_gv_features_v2;

COMMIT;

\copy (SELECT * FROM mimic_custom.qc_cohort_flow_v2 ORDER BY step_name) TO 'outputs/qc/qc_cohort_flow_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_table_counts_v2 ORDER BY table_name) TO 'outputs/qc/qc_table_counts_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_feature_missingness_v2 ORDER BY missing_rate DESC) TO 'outputs/qc/qc_feature_missingness_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_procedure_distribution_v2 ORDER BY n DESC) TO 'outputs/qc/qc_procedure_distribution_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_outcome_summary_v2) TO 'outputs/qc/qc_outcome_summary_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_shr_gv_distribution_v2) TO 'outputs/qc/qc_shr_gv_distribution_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_glucose_density_v2) TO 'outputs/qc/qc_glucose_density_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_duplicate_keys_v2) TO 'outputs/qc/qc_duplicate_keys_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_review_counts_v2 ORDER BY review_artifact) TO 'outputs/qc/qc_review_counts_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.qc_treatment_reliability_v2) TO 'outputs/qc/qc_treatment_reliability_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
