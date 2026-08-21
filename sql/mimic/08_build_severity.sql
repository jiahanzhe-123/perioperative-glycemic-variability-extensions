-- 08_build_severity.sql
-- APSIII is deliberately not fabricated. This table exposes a transparent
-- severity proxy and records that a full APSIII implementation is unavailable
-- without importing/rewriting a validated derived concept.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.severity_v2;

CREATE TABLE mimic_custom.severity_v2 AS
WITH components AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           c.age_ge_75,
           COALESCE(t.high_lactate_postop_24h_flag,false) AS high_lactate_postop_24h,
           (COALESCE(lp.creat_postop_first,0) >= 2 OR
            (lp.creat_postop_first IS NOT NULL AND la.creat_adm_first IS NOT NULL
             AND lp.creat_postop_first >= 1.5 * NULLIF(la.creat_adm_first,0))) AS creatinine_elevation,
           (lp.wbc_postop_first > 12 OR lp.wbc_postop_first < 4) AS wbc_elevation,
           (la.albumin_adm_first IS NOT NULL AND la.albumin_adm_first < 3) AS low_albumin,
           COALESCE(t.vasopressor_any_postop_24h,false) AS vaso_flag,
           COALESCE(t.mechanical_ventilation_postop_24h_flag,false) AS invasive_ventilation_flag,
           COALESCE(t.rrt_postop_7d_flag,false) AS rrt_flag,
           COALESCE(t.severe_acidosis_postop_24h_flag,false) AS acidosis_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN mimic_custom.labs_postop_first_v2 lp ON lp.stay_id = c.stay_id
    LEFT JOIN mimic_custom.labs_adm_first_v2 la ON la.stay_id = c.stay_id
    LEFT JOIN mimic_custom.treatments_postop_v2 t ON t.stay_id = c.stay_id
)
SELECT subject_id,
       hadm_id,
       stay_id,
       NULL::double precision AS apsiii,
       false AS apsiii_available,
       (CASE WHEN age_ge_75 THEN 1 ELSE 0 END)
         + (CASE WHEN high_lactate_postop_24h THEN 1 ELSE 0 END)
         + (CASE WHEN creatinine_elevation THEN 1 ELSE 0 END)
         + (CASE WHEN wbc_elevation THEN 1 ELSE 0 END)
         + (CASE WHEN low_albumin THEN 1 ELSE 0 END)
         + (CASE WHEN vaso_flag THEN 1 ELSE 0 END)
         + (CASE WHEN invasive_ventilation_flag THEN 1 ELSE 0 END)
         + (CASE WHEN rrt_flag THEN 1 ELSE 0 END)
         + (CASE WHEN acidosis_flag THEN 1 ELSE 0 END) AS severity_proxy_score,
       concat_ws('|',
           CASE WHEN age_ge_75 THEN 'age_ge_75' END,
           CASE WHEN high_lactate_postop_24h THEN 'lactate_ge_2' END,
           CASE WHEN creatinine_elevation THEN 'creatinine_elevation' END,
           CASE WHEN wbc_elevation THEN 'wbc_outside_4_12' END,
           CASE WHEN low_albumin THEN 'albumin_lt_3' END,
           CASE WHEN vaso_flag THEN 'vasopressor' END,
           CASE WHEN invasive_ventilation_flag THEN 'invasive_ventilation' END,
           CASE WHEN rrt_flag THEN 'rrt' END,
           CASE WHEN acidosis_flag THEN 'ph_lt_7_2' END) AS severity_proxy_components,
       'full APSIII not implemented; proxy is not interchangeable with APSIII'::text AS apsiii_missing_reason,
       'mimic_custom labs/vitals/treatments; transparent proxy only'::text AS severity_source
FROM components;

CREATE INDEX severity_v2_stay_idx ON mimic_custom.severity_v2 (stay_id);

COMMIT;
