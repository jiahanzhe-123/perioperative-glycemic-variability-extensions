-- 04_extract_vitals_24h.sql
-- Every chartevents access is constrained by the final cohort's stay_id and
-- 24-hour index-time window before aggregation.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.vitals_postop_24h_v2;

CREATE TABLE mimic_custom.vitals_postop_24h_v2 AS
WITH raw AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           ir.variable_name,
           ce.charttime,
           ce.valuenum::double precision AS raw_value,
           ce.valueuom,
           ce.itemid
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_icu.chartevents ce
      ON ce.subject_id = c.subject_id
     AND ce.hadm_id = c.hadm_id
     AND ce.stay_id = c.stay_id
     AND ce.charttime >= c.surgery_time
     AND ce.charttime < c.surgery_time + interval '24 hours'
    JOIN mimic_custom.itemid_review_v2 ir
      ON ir.source_table = 'mimiciv_icu.chartevents'
     AND ir.itemid = ce.itemid
     AND ir.include_candidate
    WHERE ir.variable_name IN ('heart_rate','respiratory_rate','sbp','dbp',
                               'temperature_c','temperature_f','spo2','fio2')
      AND ce.valuenum IS NOT NULL
), normalized AS (
    SELECT r.*,
           CASE r.variable_name
             WHEN 'heart_rate' THEN CASE WHEN r.raw_value BETWEEN 20 AND 250 THEN r.raw_value END
             WHEN 'respiratory_rate' THEN CASE WHEN r.raw_value BETWEEN 4 AND 80 THEN r.raw_value END
             WHEN 'sbp' THEN CASE WHEN r.raw_value BETWEEN 30 AND 300 THEN r.raw_value END
             WHEN 'dbp' THEN CASE WHEN r.raw_value BETWEEN 10 AND 200 THEN r.raw_value END
             WHEN 'temperature_c' THEN CASE WHEN r.raw_value BETWEEN 25 AND 45 THEN r.raw_value END
             WHEN 'temperature_f' THEN CASE WHEN (r.raw_value - 32.0) * 5.0 / 9.0 BETWEEN 25 AND 45
                                             THEN (r.raw_value - 32.0) * 5.0 / 9.0 END
             WHEN 'spo2' THEN CASE WHEN r.raw_value BETWEEN 50 AND 100 THEN r.raw_value END
             WHEN 'fio2' THEN CASE WHEN r.raw_value BETWEEN 0 AND 1 THEN r.raw_value * 100.0
                                   WHEN r.raw_value > 1 AND r.raw_value <= 100 THEN r.raw_value END
           END AS value_normalized,
           (r.variable_name = 'fio2' AND (r.raw_value <= 1 OR r.valueuom IS NULL)) AS unit_inference_flag
    FROM raw r
)
SELECT c.subject_id,
       c.hadm_id,
       c.stay_id,
       avg(n.value_normalized) FILTER (WHERE n.variable_name = 'heart_rate') AS hr_mean_postop_24h,
       max(n.value_normalized) FILTER (WHERE n.variable_name = 'heart_rate') AS hr_max_postop_24h,
       min(n.value_normalized) FILTER (WHERE n.variable_name = 'heart_rate') AS hr_min_postop_24h,
       count(n.value_normalized) FILTER (WHERE n.variable_name = 'heart_rate')::integer AS hr_n,
       count(*) FILTER (WHERE n.variable_name = 'heart_rate' AND n.value_normalized IS NULL)::integer AS hr_outlier_count,
       avg(n.value_normalized) FILTER (WHERE n.variable_name = 'respiratory_rate') AS rr_mean_postop_24h,
       max(n.value_normalized) FILTER (WHERE n.variable_name = 'respiratory_rate') AS rr_max_postop_24h,
       min(n.value_normalized) FILTER (WHERE n.variable_name = 'respiratory_rate') AS rr_min_postop_24h,
       count(n.value_normalized) FILTER (WHERE n.variable_name = 'respiratory_rate')::integer AS rr_n,
       count(*) FILTER (WHERE n.variable_name = 'respiratory_rate' AND n.value_normalized IS NULL)::integer AS rr_outlier_count,
       avg(n.value_normalized) FILTER (WHERE n.variable_name = 'sbp') AS sbp_mean_postop_24h,
       max(n.value_normalized) FILTER (WHERE n.variable_name = 'sbp') AS sbp_max_postop_24h,
       min(n.value_normalized) FILTER (WHERE n.variable_name = 'sbp') AS sbp_min_postop_24h,
       count(n.value_normalized) FILTER (WHERE n.variable_name = 'sbp')::integer AS sbp_n,
       count(*) FILTER (WHERE n.variable_name = 'sbp' AND n.value_normalized IS NULL)::integer AS sbp_outlier_count,
       avg(n.value_normalized) FILTER (WHERE n.variable_name = 'dbp') AS dbp_mean_postop_24h,
       max(n.value_normalized) FILTER (WHERE n.variable_name = 'dbp') AS dbp_max_postop_24h,
       min(n.value_normalized) FILTER (WHERE n.variable_name = 'dbp') AS dbp_min_postop_24h,
       count(n.value_normalized) FILTER (WHERE n.variable_name = 'dbp')::integer AS dbp_n,
       count(*) FILTER (WHERE n.variable_name = 'dbp' AND n.value_normalized IS NULL)::integer AS dbp_outlier_count,
       avg(n.value_normalized) FILTER (WHERE n.variable_name IN ('temperature_c','temperature_f')) AS temp_c_mean_postop_24h,
       max(n.value_normalized) FILTER (WHERE n.variable_name IN ('temperature_c','temperature_f')) AS temp_c_max_postop_24h,
       min(n.value_normalized) FILTER (WHERE n.variable_name IN ('temperature_c','temperature_f')) AS temp_c_min_postop_24h,
       count(n.value_normalized) FILTER (WHERE n.variable_name IN ('temperature_c','temperature_f'))::integer AS temp_c_n,
       count(*) FILTER (WHERE n.variable_name IN ('temperature_c','temperature_f') AND n.value_normalized IS NULL)::integer AS temp_c_outlier_count,
       max(n.value_normalized) FILTER (WHERE n.variable_name = 'spo2') AS spo2_max_postop_24h,
       min(n.value_normalized) FILTER (WHERE n.variable_name = 'spo2') AS spo2_min_postop_24h,
       count(n.value_normalized) FILTER (WHERE n.variable_name = 'spo2')::integer AS spo2_n,
       count(*) FILTER (WHERE n.variable_name = 'spo2' AND n.value_normalized IS NULL)::integer AS spo2_outlier_count,
       avg(n.value_normalized) FILTER (WHERE n.variable_name = 'fio2') AS fio2_mean_postop_24h,
       max(n.value_normalized) FILTER (WHERE n.variable_name = 'fio2') AS fio2_max_postop_24h,
       min(n.value_normalized) FILTER (WHERE n.variable_name = 'fio2') AS fio2_min_postop_24h,
       count(n.value_normalized) FILTER (WHERE n.variable_name = 'fio2')::integer AS fio2_n,
       count(*) FILTER (WHERE n.variable_name = 'fio2' AND n.value_normalized IS NULL)::integer AS fio2_outlier_count,
       bool_or(n.unit_inference_flag) FILTER (WHERE n.variable_name = 'fio2') AS fio2_unit_inference_flag,
       'mimiciv_icu.chartevents; cohort-stay restricted; curated itemid_review_v2'::text AS vitals_source,
       CASE WHEN count(n.value_normalized) = 0 THEN 'no valid curated vital/FiO2 observations in postoperative 24h window'
            ELSE NULL END AS vitals_missing_reason
FROM mimic_custom.cardiac_surgery_cohort_v2 c
LEFT JOIN normalized n
  ON n.subject_id = c.subject_id
 AND n.hadm_id = c.hadm_id
 AND n.stay_id = c.stay_id
GROUP BY c.subject_id, c.hadm_id, c.stay_id;

CREATE INDEX vitals_postop_24h_v2_stay_idx
    ON mimic_custom.vitals_postop_24h_v2 (stay_id);

COMMIT;
