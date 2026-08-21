-- Step 6a: 血糖提取与 GV 指标
-- 规则:
--  来源: lab 表 labname IN ('glucose','bedside glucose')(排除 'glucose - CSF');单位均为 mg/dL(Step1 已核实)
--  窗口: labresultoffset 0..1440 分钟(相对 ICU 入科)
--  合理性: 保留 20..1500 mg/dL;范围外视为录入错误并剔除(在 QC 中计数)
--  去重: 同 stay+同 offset+同值+同来源 的完全重复行只保留一条
--  nursecharting 'Bedside Glucose' 作为平行源在 QC 中做重叠分析,不混入主提取
BEGIN;
DROP TABLE IF EXISTS eicu_cardio_validation.glucose_raw;
CREATE TABLE eicu_cardio_validation.glucose_raw AS
WITH c AS (SELECT patientunitstayid FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND rn_patient=1)
SELECT l.patientunitstayid, l.labresultoffset AS offset_min,
       l.labresult AS value_mgdl,           -- original value
       l.labmeasurenamesystem AS unit_orig, -- original unit
       l.labname AS source_name,            -- 'glucose'(中心实验室) or 'bedside glucose'(POCT)
       'lab'::text AS source_table
FROM lab l
JOIN c USING (patientunitstayid)
WHERE l.labname IN ('glucose','bedside glucose')
  AND l.labresultoffset BETWEEN 0 AND 1440
  AND l.labresult IS NOT NULL;
CREATE INDEX ON eicu_cardio_validation.glucose_raw(patientunitstayid);

-- 剔除不合理值
DROP TABLE IF EXISTS eicu_cardio_validation.glucose_clean;
CREATE TABLE eicu_cardio_validation.glucose_clean AS
SELECT DISTINCT ON (patientunitstayid, offset_min, value_mgdl, source_name)
       *
FROM eicu_cardio_validation.glucose_raw
WHERE value_mgdl BETWEEN 20 AND 1500
ORDER BY patientunitstayid, offset_min, value_mgdl, source_name;
CREATE INDEX ON eicu_cardio_validation.glucose_clean(patientunitstayid);

-- 每 stay 汇总
DROP TABLE IF EXISTS eicu_cardio_validation.glucose_features;
CREATE TABLE eicu_cardio_validation.glucose_features AS
SELECT patientunitstayid,
       count(*) AS glucose_n,
       avg(value_mgdl) AS glucose_mean_24h,
       stddev_samp(value_mgdl) AS glucose_sd_24h,
       min(value_mgdl) AS glucose_min_24h,
       max(value_mgdl) AS glucose_max_24h,
       (max(value_mgdl) - min(value_mgdl)) AS glucose_range_24h,
       (stddev_samp(value_mgdl) / avg(value_mgdl)) AS glucose_cv_24h,
       (array_agg(value_mgdl ORDER BY offset_min))[1] AS first_glucose,
       (array_agg(value_mgdl ORDER BY offset_min DESC))[1] AS last_glucose,
       (max(offset_min) - min(offset_min)) AS measurement_span_minutes,
       count(*) / nullif((max(offset_min) - min(offset_min))/60.0, 0) AS glucose_measurement_density,
       count(*) FILTER (WHERE source_name='bedside glucose') AS n_bedside,
       count(*) FILTER (WHERE source_name='glucose') AS n_lab,
       count(DISTINCT source_name) AS glucose_source_count
FROM eicu_cardio_validation.glucose_clean
GROUP BY patientunitstayid;
CREATE UNIQUE INDEX ON eicu_cardio_validation.glucose_features(patientunitstayid);
COMMIT;

-- 汇总
SELECT count(*) AS stays_with_any_glucose,
       count(*) FILTER (WHERE glucose_n>=2) AS stays_ge2,
       percentile_cont(0.5) WITHIN GROUP (ORDER BY glucose_n) AS med_n,
       percentile_cont(0.25) WITHIN GROUP (ORDER BY glucose_n) AS q1_n,
       percentile_cont(0.75) WITHIN GROUP (ORDER BY glucose_n) AS q3_n,
       min(glucose_n), max(glucose_n)
FROM eicu_cardio_validation.glucose_features;
SELECT CASE WHEN glucose_n=2 THEN 'n=2' WHEN glucose_n BETWEEN 3 AND 5 THEN 'n=3-5' WHEN glucose_n BETWEEN 6 AND 10 THEN 'n=6-10' ELSE 'n>10' END AS grp, count(*)
FROM eicu_cardio_validation.glucose_features WHERE glucose_n>=2 GROUP BY 1 ORDER BY 1;
-- 剔除计数
SELECT (SELECT count(*) FROM eicu_cardio_validation.glucose_raw) AS raw_n,
       (SELECT count(*) FROM eicu_cardio_validation.glucose_clean) AS clean_n;
