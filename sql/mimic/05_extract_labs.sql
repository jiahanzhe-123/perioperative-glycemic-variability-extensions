-- 05_extract_labs.sql
-- Laboratory events are restricted by subject_id/hadm_id and the cohort's
-- surgery-centered window. The long tables preserve the deduplication and
-- unit/outlier provenance used by the wide tables.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.hba1c_baseline_v2;
DROP TABLE IF EXISTS mimic_custom.glucose_summary_v2;
DROP TABLE IF EXISTS mimic_custom.glucose_long_v2;
DROP TABLE IF EXISTS mimic_custom.glucose_lab_values_v2;
DROP TABLE IF EXISTS mimic_custom.labs_adm_first_v2;
DROP TABLE IF EXISTS mimic_custom.labs_postop_first_v2;
DROP TABLE IF EXISTS mimic_custom.labs_long_v2;

CREATE TABLE mimic_custom.labs_long_v2 AS
WITH raw AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           c.admittime_ts,
           c.surgery_time,
           ir.variable_name,
           le.charttime,
           le.labevent_id,
           le.itemid,
           le.valuenum::double precision AS raw_value,
           le.valueuom,
           CASE WHEN le.valueuom IS NULL OR btrim(le.valueuom) = '' THEN true ELSE false END AS unit_inference_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_hosp.labevents le
      ON le.subject_id = c.subject_id
     AND le.hadm_id = c.hadm_id
     AND le.charttime >= c.surgery_time - interval '365 days'
     AND le.charttime < c.surgery_time + interval '24 hours'
    JOIN mimic_custom.itemid_review_v2 ir
      ON ir.source_table = 'mimiciv_hosp.labevents'
     AND ir.itemid = le.itemid
     AND ir.include_candidate
    WHERE le.valuenum IS NOT NULL
), normalized AS (
    SELECT r.*,
           CASE r.variable_name
             WHEN 'hba1c_pct' THEN CASE WHEN r.raw_value BETWEEN 2 AND 20 THEN r.raw_value END
             WHEN 'glucose' THEN CASE WHEN r.raw_value BETWEEN 20 AND 1000 THEN r.raw_value END
             WHEN 'creatinine' THEN CASE WHEN r.raw_value BETWEEN 0 AND 30 THEN r.raw_value END
             WHEN 'bun' THEN CASE WHEN r.raw_value BETWEEN 0 AND 200 THEN r.raw_value END
             WHEN 'lactate' THEN CASE WHEN r.raw_value BETWEEN 0 AND 50 THEN r.raw_value END
             WHEN 'ph' THEN CASE WHEN r.raw_value BETWEEN 6.5 AND 8.0 THEN r.raw_value END
             WHEN 'wbc' THEN CASE WHEN r.raw_value BETWEEN 0 AND 300 THEN r.raw_value END
             WHEN 'hemoglobin' THEN CASE WHEN r.raw_value BETWEEN 0 AND 30 THEN r.raw_value END
             WHEN 'rbc' THEN CASE WHEN r.raw_value BETWEEN 0 AND 15 THEN r.raw_value END
             WHEN 'hematocrit' THEN CASE WHEN r.raw_value BETWEEN 0 AND 100 THEN r.raw_value END
             WHEN 'mcv' THEN CASE WHEN r.raw_value BETWEEN 0 AND 200 THEN r.raw_value END
             WHEN 'bilirubin' THEN CASE WHEN r.raw_value BETWEEN 0 AND 100 THEN r.raw_value END
             WHEN 'alt' THEN CASE WHEN r.raw_value BETWEEN 0 AND 10000 THEN r.raw_value END
             WHEN 'alp' THEN CASE WHEN r.raw_value BETWEEN 0 AND 10000 THEN r.raw_value END
             WHEN 'ast' THEN CASE WHEN r.raw_value BETWEEN 0 AND 10000 THEN r.raw_value END
             WHEN 'albumin' THEN CASE WHEN r.raw_value BETWEEN 0 AND 10 THEN r.raw_value END
             WHEN 'potassium' THEN CASE WHEN r.raw_value BETWEEN 1 AND 15 THEN r.raw_value END
             WHEN 'sodium' THEN CASE WHEN r.raw_value BETWEEN 80 AND 200 THEN r.raw_value END
             WHEN 'chloride' THEN CASE WHEN r.raw_value BETWEEN 50 AND 150 THEN r.raw_value END
             WHEN 'calcium' THEN CASE WHEN r.raw_value BETWEEN 2 AND 20 THEN r.raw_value END
             WHEN 'anion_gap' THEN CASE WHEN r.raw_value BETWEEN 0 AND 80 THEN r.raw_value END
             WHEN 'bicarbonate' THEN CASE WHEN r.raw_value BETWEEN 0 AND 80 THEN r.raw_value END
             WHEN 'pt' THEN CASE WHEN r.raw_value BETWEEN 0 AND 200 THEN r.raw_value END
             WHEN 'ptt' THEN CASE WHEN r.raw_value BETWEEN 0 AND 300 THEN r.raw_value END
             WHEN 'fibrinogen' THEN CASE WHEN r.raw_value BETWEEN 0 AND 3000 THEN r.raw_value END
             WHEN 'platelets' THEN CASE WHEN r.raw_value BETWEEN 0 AND 3000 THEN r.raw_value END
           END AS value_valid
    FROM raw r
), dedup AS (
    SELECT n.*,
           row_number() OVER (
             PARTITION BY n.stay_id, n.variable_name, n.charttime
             ORDER BY n.labevent_id
           ) AS dedup_rank
    FROM normalized n
)
SELECT subject_id,
       hadm_id,
       stay_id,
       admittime_ts,
       surgery_time,
       variable_name,
       charttime,
       labevent_id,
       itemid,
       raw_value,
       value_valid,
       valueuom,
       unit_inference_flag,
       (value_valid IS NULL) AS outlier_flag,
       'mimiciv_hosp.labevents; duplicate rule = first labevent_id per stay/variable/charttime'::text AS source_provenance
FROM dedup
WHERE dedup_rank = 1;

CREATE INDEX labs_long_v2_stay_time_idx
    ON mimic_custom.labs_long_v2 (stay_id, charttime);
CREATE INDEX labs_long_v2_var_idx
    ON mimic_custom.labs_long_v2 (variable_name, stay_id);

-- Commit the expensive but reusable lab scan before the ICU chartevents
-- glucose scan. A separate indexed cohort-limited table makes the duplicate
-- check a targeted lookup rather than a repeated CTE re-scan.
COMMIT;
BEGIN;

CREATE TABLE mimic_custom.glucose_lab_values_v2 AS
SELECT subject_id,
       hadm_id,
       stay_id,
       charttime,
       labevent_id::bigint AS source_id,
       value_valid AS glucose_mg_dl,
       raw_value,
       valueuom,
       unit_inference_flag,
       outlier_flag
FROM mimic_custom.labs_long_v2
WHERE variable_name = 'glucose';

CREATE INDEX glucose_lab_values_v2_stay_time_idx
    ON mimic_custom.glucose_lab_values_v2 (stay_id, charttime, glucose_mg_dl);

-- First valid postoperative lab in the 24-hour index window.
CREATE TABLE mimic_custom.labs_postop_first_v2 AS
WITH firsts AS (
    SELECT DISTINCT ON (stay_id, variable_name)
           stay_id, variable_name, charttime, value_valid, valueuom, unit_inference_flag
    FROM mimic_custom.labs_long_v2
    WHERE value_valid IS NOT NULL
      AND charttime >= surgery_time
      AND charttime < surgery_time + interval '24 hours'
    ORDER BY stay_id, variable_name, charttime, labevent_id
)
SELECT c.subject_id,
       c.hadm_id,
       c.stay_id,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'creatinine') AS creat_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'bun') AS bun_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'lactate') AS lactate_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'ph') AS ph_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'wbc') AS wbc_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'hemoglobin') AS hgb_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'rbc') AS rbc_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'hematocrit') AS hct_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'mcv') AS mcv_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'potassium') AS k_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'sodium') AS na_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'chloride') AS cl_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'calcium') AS ca_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'anion_gap') AS anion_gap_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'bicarbonate') AS bicarb_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'pt') AS pt_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'ptt') AS ptt_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'fibrinogen') AS fibrinogen_postop_first,
       max(f.value_valid) FILTER (WHERE f.variable_name = 'platelets') AS platelets_postop_first,
       max(f.charttime) FILTER (WHERE f.variable_name = 'creatinine') AS creat_postop_first_time,
       max(f.charttime) FILTER (WHERE f.variable_name = 'lactate') AS lactate_postop_first_time,
       count(f.value_valid)::integer AS n_valid_postop_first_labs,
       CASE WHEN count(f.value_valid) = 0 THEN 'no valid lab in postoperative 24h window' ELSE NULL END AS labs_missing_reason,
       'mimiciv_hosp.labevents; first valid value by variable in surgery_time to +24h'::text AS labs_source
FROM mimic_custom.cardiac_surgery_cohort_v2 c
LEFT JOIN firsts f ON f.stay_id = c.stay_id
GROUP BY c.subject_id, c.hadm_id, c.stay_id;

-- Admission laboratory table keeps both admission-first and the latest
-- pre-procedure value; final *_adm_first uses admission-first when available,
-- otherwise the latest admission-period pre-procedure result.
CREATE TABLE mimic_custom.labs_adm_first_v2 AS
WITH admission_first AS (
    SELECT DISTINCT ON (stay_id, variable_name)
           stay_id, variable_name, charttime, value_valid, valueuom, unit_inference_flag
    FROM mimic_custom.labs_long_v2
    WHERE value_valid IS NOT NULL
      AND charttime >= admittime_ts
      AND charttime < admittime_ts + interval '24 hours'
    ORDER BY stay_id, variable_name, charttime, labevent_id
), preop_last AS (
    SELECT DISTINCT ON (stay_id, variable_name)
           stay_id, variable_name, charttime, value_valid, valueuom, unit_inference_flag
    FROM mimic_custom.labs_long_v2
    WHERE value_valid IS NOT NULL
      AND charttime >= admittime_ts
      AND charttime < surgery_time
    ORDER BY stay_id, variable_name, charttime DESC, labevent_id DESC
), combined AS (
    SELECT c.subject_id, c.hadm_id, c.stay_id,
           r.variable_name,
           COALESCE(a.value_valid, p.value_valid) AS value_valid,
           COALESCE(a.charttime, p.charttime) AS charttime,
           COALESCE(a.valueuom, p.valueuom) AS valueuom,
           COALESCE(a.unit_inference_flag, p.unit_inference_flag) AS unit_inference_flag,
           CASE WHEN a.value_valid IS NOT NULL THEN 'admission_24h_first'
                WHEN p.value_valid IS NOT NULL THEN 'admission_period_preop_last'
                ELSE 'missing' END AS definition
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    CROSS JOIN (SELECT DISTINCT variable_name FROM mimic_custom.labs_long_v2) r
    LEFT JOIN admission_first a ON a.stay_id = c.stay_id AND a.variable_name = r.variable_name
    LEFT JOIN preop_last p ON p.stay_id = c.stay_id AND p.variable_name = r.variable_name
)
SELECT c.subject_id,
       c.hadm_id,
       c.stay_id,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'creatinine') AS creat_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'bun') AS bun_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'lactate') AS lactate_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'ph') AS ph_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'wbc') AS wbc_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'hemoglobin') AS hgb_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'rbc') AS rbc_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'hematocrit') AS hct_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'mcv') AS mcv_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'bilirubin') AS bili_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'alt') AS alt_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'alp') AS alp_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'ast') AS ast_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'albumin') AS albumin_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'potassium') AS k_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'sodium') AS na_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'chloride') AS cl_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'calcium') AS ca_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'anion_gap') AS anion_gap_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'bicarbonate') AS bicarb_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'pt') AS pt_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'ptt') AS ptt_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'fibrinogen') AS fibrinogen_adm_first,
       max(co.value_valid) FILTER (WHERE co.variable_name = 'platelets') AS platelets_adm_first,
       max(co.charttime) FILTER (WHERE co.variable_name = 'creatinine') AS creat_adm_first_time,
       max(co.definition) FILTER (WHERE co.variable_name = 'creatinine') AS creat_adm_definition,
       max(co.definition) FILTER (WHERE co.variable_name = 'albumin') AS albumin_adm_definition,
       count(co.value_valid)::integer AS n_valid_adm_labs,
       CASE WHEN count(co.value_valid) = 0 THEN 'no valid lab in admission/period-preop definition' ELSE NULL END AS labs_missing_reason,
       'mimiciv_hosp.labevents; admission 24h first preferred, preop last fallback'::text AS labs_source
FROM mimic_custom.cardiac_surgery_cohort_v2 c
LEFT JOIN combined co ON co.stay_id = c.stay_id
GROUP BY c.subject_id, c.hadm_id, c.stay_id;

-- Glucose combines labevents and chartevents. Labevents are preferred when a
-- charted ICU value is within 2 minutes and within 1 mg/dL of a lab value.
CREATE TABLE mimic_custom.glucose_long_v2 AS
WITH lab_values AS (
    SELECT l.subject_id, l.hadm_id, l.stay_id, l.charttime, l.source_id,
           l.glucose_mg_dl,
           l.raw_value,
           l.valueuom,
           l.unit_inference_flag,
           'labevents'::text AS source_table,
           l.outlier_flag
    FROM mimic_custom.glucose_lab_values_v2 l
), chart_values AS (
    SELECT c.subject_id, c.hadm_id, c.stay_id, ce.charttime, ce.itemid::bigint AS source_id,
           CASE WHEN ce.valuenum BETWEEN 20 AND 1000 THEN ce.valuenum::double precision END AS glucose_mg_dl,
           ce.valuenum::double precision AS raw_value,
           ce.valueuom,
           (ce.valueuom IS NULL OR btrim(ce.valueuom) = '') AS unit_inference_flag,
           'chartevents'::text AS source_table,
           (ce.valuenum IS NULL OR ce.valuenum NOT BETWEEN 20 AND 1000) AS outlier_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_icu.chartevents ce
      ON ce.subject_id = c.subject_id
     AND ce.hadm_id = c.hadm_id
     AND ce.stay_id = c.stay_id
     AND ce.charttime >= c.surgery_time
     AND ce.charttime < c.surgery_time + interval '7 days'
    JOIN mimic_custom.itemid_review_v2 ir
      ON ir.source_table = 'mimiciv_icu.chartevents'
     AND ir.itemid = ce.itemid
     AND ir.include_candidate
     AND ir.variable_name = 'glucose'
    WHERE ce.valuenum IS NOT NULL
), all_values AS (
    SELECT * FROM lab_values
    UNION ALL
    SELECT * FROM chart_values
), dedup AS (
    SELECT v.*,
           CASE WHEN v.source_table = 'chartevents'
                     AND EXISTS (
                       SELECT 1 FROM mimic_custom.glucose_lab_values_v2 l
                       WHERE l.stay_id = v.stay_id
                         AND l.charttime BETWEEN v.charttime - interval '2 minutes' AND v.charttime + interval '2 minutes'
                         AND l.glucose_mg_dl IS NOT NULL
                         AND v.glucose_mg_dl IS NOT NULL
                         AND abs(l.glucose_mg_dl - v.glucose_mg_dl) <= 1
                     ) THEN 1 ELSE 0 END AS duplicate_of_lab
    FROM all_values v
)
SELECT subject_id,
       hadm_id,
       stay_id,
       charttime,
       source_id,
       glucose_mg_dl,
       raw_value,
       valueuom,
       unit_inference_flag,
       source_table,
       outlier_flag,
       (duplicate_of_lab = 1) AS duplicate_of_labevent,
       CASE WHEN duplicate_of_lab = 1 THEN 'chartevents removed when lab within 2 min and <=1 mg/dL'
            ELSE 'retained' END AS dedup_reason
FROM dedup
WHERE duplicate_of_lab = 0;

CREATE INDEX glucose_long_v2_stay_time_idx
    ON mimic_custom.glucose_long_v2 (stay_id, charttime);

-- Keep the expensive ICU glucose scan reusable if a later summary expression
-- fails validation.
COMMIT;
BEGIN;

CREATE TABLE mimic_custom.glucose_summary_v2 AS
WITH windows(window_name, window_hours) AS (
    VALUES ('24h',24),('48h',48),('72h',72),('7d',168)
), values_with_window AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           w.window_name,
           g.charttime,
           g.glucose_mg_dl,
           g.source_table,
           g.unit_inference_flag,
           g.outlier_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    CROSS JOIN windows w
    LEFT JOIN mimic_custom.glucose_long_v2 g
      ON g.stay_id = c.stay_id
     AND g.charttime >= c.surgery_time
     AND g.charttime < c.surgery_time + (w.window_hours || ' hours')::interval
), valid AS (
    SELECT * FROM values_with_window WHERE glucose_mg_dl IS NOT NULL
), ordered AS (
    SELECT v.*,
           lag(v.glucose_mg_dl) OVER (PARTITION BY v.stay_id, v.window_name ORDER BY v.charttime) AS previous_glucose
    FROM valid v
), stats AS (
    SELECT subject_id, hadm_id, stay_id, window_name,
           count(glucose_mg_dl)::integer AS n_glucose,
           avg(glucose_mg_dl) AS mean_glucose,
           stddev_samp(glucose_mg_dl) AS sd_glucose,
           min(glucose_mg_dl) AS min_glucose,
           max(glucose_mg_dl) AS max_glucose,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY glucose_mg_dl) AS median_glucose,
           avg(abs(glucose_mg_dl - previous_glucose)) FILTER (WHERE previous_glucose IS NOT NULL) AS masd_glucose,
           string_agg(DISTINCT source_table, '|' ORDER BY source_table) AS source_mix,
           bool_or(unit_inference_flag) AS unit_inference_flag
    FROM ordered
    GROUP BY subject_id, hadm_id, stay_id, window_name
), outliers AS (
    SELECT subject_id, hadm_id, stay_id, window_name,
           count(*) FILTER (WHERE outlier_flag)::integer AS outlier_count
    FROM values_with_window
    GROUP BY subject_id, hadm_id, stay_id, window_name
)
SELECT c.subject_id,
       c.hadm_id,
       c.stay_id,
       max(s.n_glucose) FILTER (WHERE s.window_name = '24h') AS n_glucose_postop_24h,
       max(s.n_glucose) FILTER (WHERE s.window_name = '48h') AS n_glucose_postop_48h,
       max(s.n_glucose) FILTER (WHERE s.window_name = '72h') AS n_glucose_postop_72h,
       max(s.n_glucose) FILTER (WHERE s.window_name = '7d') AS n_glucose_postop_7d,
       max(s.mean_glucose) FILTER (WHERE s.window_name = '24h') AS glucose_mean_postop_24h,
       max(s.sd_glucose) FILTER (WHERE s.window_name = '24h') AS glucose_sd_postop_24h,
       max(s.min_glucose) FILTER (WHERE s.window_name = '24h') AS glucose_min_postop_24h,
       max(s.max_glucose) FILTER (WHERE s.window_name = '24h') AS glucose_max_postop_24h,
       max(s.median_glucose) FILTER (WHERE s.window_name = '24h') AS glucose_median_postop_24h,
       max(s.sd_glucose / NULLIF(s.mean_glucose,0)) FILTER (WHERE s.window_name = '24h') AS glucose_cv_postop_24h,
       max(s.max_glucose - s.min_glucose) FILTER (WHERE s.window_name = '24h') AS glucose_range_postop_24h,
       max(s.max_glucose - s.min_glucose) FILTER (WHERE s.window_name = '24h') AS glucose_variability_delta_postop_24h,
       max(s.masd_glucose) FILTER (WHERE s.window_name = '24h') AS glucose_masd_postop_24h,
       max(s.mean_glucose) FILTER (WHERE s.window_name = '48h') AS glucose_mean_postop_48h,
       max(s.mean_glucose) FILTER (WHERE s.window_name = '72h') AS glucose_mean_postop_72h,
       max(s.mean_glucose) FILTER (WHERE s.window_name = '7d') AS glucose_mean_postop_7d,
       max(s.source_mix) FILTER (WHERE s.window_name = '24h') AS glucose_source_mix,
       bool_or(s.unit_inference_flag) FILTER (WHERE s.window_name = '24h') AS glucose_unit_inference_flag,
       COALESCE(max(o.outlier_count) FILTER (WHERE o.window_name = '24h'),0)::integer AS glucose_outlier_count,
       CASE WHEN max(s.n_glucose) FILTER (WHERE s.window_name = '24h') IS NULL THEN 'no valid postoperative 24h glucose'
            ELSE NULL END AS glucose_missing_reason,
       'mimiciv_hosp.labevents + mimiciv_icu.chartevents; labevent-preferred de-duplication'::text AS glucose_source
FROM mimic_custom.cardiac_surgery_cohort_v2 c
LEFT JOIN stats s ON s.stay_id = c.stay_id
LEFT JOIN outliers o ON o.stay_id = c.stay_id AND o.window_name = s.window_name
GROUP BY c.subject_id, c.hadm_id, c.stay_id;

CREATE TABLE mimic_custom.hba1c_baseline_v2 AS
WITH candidates AS (
    SELECT l.subject_id,
           l.hadm_id,
           l.stay_id,
           l.charttime AS hba1c_time,
           l.value_valid AS hba1c_pct,
           l.unit_inference_flag,
           CASE WHEN l.charttime < c.surgery_time
                     AND l.charttime >= c.surgery_time - interval '90 days' THEN 'pre_90d'
                WHEN l.charttime < c.surgery_time
                     AND l.charttime >= c.surgery_time - interval '365 days' THEN 'pre_365d'
                WHEN l.charttime >= c.surgery_time
                     AND l.charttime < c.surgery_time + interval '24 hours' THEN 'post_24h'
           END AS hba1c_window,
           CASE WHEN l.charttime < c.surgery_time
                     AND l.charttime >= c.surgery_time - interval '90 days' THEN 1
                WHEN l.charttime < c.surgery_time
                     AND l.charttime >= c.surgery_time - interval '365 days' THEN 2
                WHEN l.charttime >= c.surgery_time
                     AND l.charttime < c.surgery_time + interval '24 hours' THEN 3
           END AS window_priority
    FROM mimic_custom.labs_long_v2 l
    JOIN mimic_custom.cardiac_surgery_cohort_v2 c ON c.stay_id = l.stay_id
    WHERE l.variable_name = 'hba1c_pct'
      AND l.value_valid IS NOT NULL
), ranked AS (
    SELECT *,
           row_number() OVER (
             PARTITION BY stay_id
             ORDER BY window_priority,
                      CASE WHEN window_priority = 3 THEN hba1c_time END ASC,
                      CASE WHEN window_priority < 3 THEN hba1c_time END DESC
           ) AS rn
    FROM candidates
    WHERE window_priority IS NOT NULL
)
SELECT c.subject_id,
       c.hadm_id,
       c.stay_id,
       r.hba1c_pct,
       r.hba1c_time,
       EXTRACT(EPOCH FROM (r.hba1c_time - c.surgery_time)) / 86400.0 AS hba1c_days_from_surgery,
       r.hba1c_window,
       CASE WHEN r.hba1c_pct IS NULL THEN NULL ELSE 'mimiciv_hosp.labevents' END AS hba1c_source,
       r.unit_inference_flag AS hba1c_unit_inference_flag,
       CASE WHEN r.hba1c_pct IS NULL THEN 'no valid HbA1c from pre_90d/pre_365d/post_24h windows' ELSE NULL END AS hba1c_missing_reason
FROM mimic_custom.cardiac_surgery_cohort_v2 c
LEFT JOIN ranked r ON r.stay_id = c.stay_id AND r.rn = 1;

CREATE INDEX labs_postop_first_v2_stay_idx ON mimic_custom.labs_postop_first_v2 (stay_id);
CREATE INDEX labs_adm_first_v2_stay_idx ON mimic_custom.labs_adm_first_v2 (stay_id);
CREATE INDEX glucose_summary_v2_stay_idx ON mimic_custom.glucose_summary_v2 (stay_id);
CREATE INDEX hba1c_baseline_v2_stay_idx ON mimic_custom.hba1c_baseline_v2 (stay_id);

COMMIT;
