-- eICU harmonized replication cohort
-- schema: multidatabase_glucose_validation (只新建,不动 eicu_cardio_validation)
BEGIN;
CREATE SCHEMA IF NOT EXISTS multidatabase_glucose_validation;
SET search_path TO multidatabase_glucose_validation, public, eicu_cardio_validation;

-- ========== 1. 队列帧(含敏感性命名的证据层) ==========
DROP TABLE IF EXISTS multidatabase_glucose_validation.eicu_cohort;
CREATE TABLE multidatabase_glucose_validation.eicu_cohort AS
WITH base AS (
  SELECT c.patientunitstayid, c.patienthealthsystemstayid, c.hospitalid, c.uniquepid,
         c.gender, c.age, c.unittype, c.unitadmitsource,
         c.hospitaladmitoffset, c.unitdischargeoffset, c.hospitaldischargeoffset,
         c.unitdischargestatus, c.hospitaldischargestatus,
         c.admissionheight, c.admissionweight,
         c.surgery_categories, c.surgery_category, c.surgery_confidence,
         c.high_specificity_flag, c.broad_cohort_flag
  FROM eicu_cardio_validation.cohort_base c
  WHERE c.rn_patient = 1
)
SELECT *,
  -- 主验证队列: CABG/open_valve 且 DEFINITE
  (surgery_category IN ('CABG','OPEN_VALVE') AND surgery_confidence='DEFINITE') AS in_main,
  -- Sensitivity A: + OPEN_AORTIC (DEFINITE)
  (surgery_category IN ('CABG','OPEN_VALVE','OPEN_AORTIC') AND surgery_confidence='DEFINITE') AS in_sensA,
  -- Sensitivity B: all high-specificity open cardiac
  (surgery_confidence='DEFINITE') AS in_sensB,
  -- Sensitivity C: broad phenotype CABG/open_valve (DEFINITE+PROBABLE)
  (surgery_category IN ('CABG','OPEN_VALVE') AND surgery_confidence IN ('DEFINITE','PROBABLE')) AS in_sensC,
  -- 手术类别: combined = CABG+VALVE 同时命中
  CASE
    WHEN 'CABG' = ANY(surgery_categories) AND 'OPEN_VALVE' = ANY(surgery_categories) THEN 'combined'
    WHEN 'CABG' = ANY(surgery_categories) THEN 'cabg'
    WHEN 'OPEN_VALVE' = ANY(surgery_categories) THEN 'valve'
    ELSE lower(surgery_category)
  END AS procedure_category,
  CASE WHEN age ~ '^[0-9]+$' THEN age::int WHEN age='> 89' THEN 90 ELSE NULL END AS age_years,
  (age = '> 89') AS age_censored
FROM base;
CREATE UNIQUE INDEX ON multidatabase_glucose_validation.eicu_cohort(patientunitstayid);

-- ========== 2. 血糖原始行(0-1440min, 两来源, 含 creatinine) ==========
DROP TABLE IF EXISTS multidatabase_glucose_validation.eicu_glucose_rows;
CREATE TABLE multidatabase_glucose_validation.eicu_glucose_rows AS
SELECT l.patientunitstayid, l.labresultoffset AS offset_min, l.labresult AS value_mgdl,
       CASE WHEN l.labname='bedside glucose' THEN 'poct' ELSE 'central_lab' END AS source_cat
FROM lab l
JOIN eicu_cardio_validation.cohort_base c ON c.patientunitstayid=l.patientunitstayid AND c.rn_patient=1
WHERE l.labname IN ('glucose','bedside glucose')
  AND l.labresultoffset BETWEEN 0 AND 1440
  AND l.labresult IS NOT NULL
  AND l.labresult BETWEEN 20 AND 1500;
CREATE INDEX ON multidatabase_glucose_validation.eicu_glucose_rows(patientunitstayid);

-- 去重: 完全重复(患者,时间,值,来源) -> 1条; 同分钟多值 -> 中位数
DROP TABLE IF EXISTS multidatabase_glucose_validation.eicu_glucose_minute;
CREATE TABLE multidatabase_glucose_validation.eicu_glucose_minute AS
WITH dedup AS (
  SELECT DISTINCT patientunitstayid, offset_min, value_mgdl, source_cat
  FROM multidatabase_glucose_validation.eicu_glucose_rows
)
SELECT patientunitstayid, offset_min,
       percentile_cont(0.5) WITHIN GROUP (ORDER BY value_mgdl) AS value_mgdl,
       string_agg(DISTINCT source_cat, '+') AS sources,
       count(*) AS n_records_in_minute
FROM dedup
GROUP BY patientunitstayid, offset_min;

-- 每 stay 特征
DROP TABLE IF EXISTS multidatabase_glucose_validation.eicu_glucose_features;
CREATE TABLE multidatabase_glucose_validation.eicu_glucose_features AS
SELECT m.patientunitstayid,
       (SELECT count(*) FROM multidatabase_glucose_validation.eicu_glucose_rows r WHERE r.patientunitstayid=m.patientunitstayid) AS glucose_n_raw,
       count(*) AS glucose_n,
       avg(value_mgdl) AS glucose_mean_24h,
       stddev_samp(value_mgdl) AS glucose_sd_24h,
       stddev_samp(value_mgdl)/avg(value_mgdl) AS glucose_cv_24h,
       min(value_mgdl) AS glucose_min_24h,
       max(value_mgdl) AS glucose_max_24h,
       max(value_mgdl)-min(value_mgdl) AS glucose_range_24h,
       (array_agg(value_mgdl ORDER BY offset_min))[1] AS glucose_first,
       (array_agg(value_mgdl ORDER BY offset_min DESC))[1] AS glucose_last,
       max(offset_min)-min(offset_min) AS measurement_span_minutes,
       sum((sources LIKE '%poct%')::int)::float/count(*) AS poct_fraction,
       sum((sources LIKE '%central_lab%')::int)::float/count(*) AS central_lab_fraction
FROM multidatabase_glucose_validation.eicu_glucose_minute m
GROUP BY m.patientunitstayid;

-- ========== 3. 协变量(creatinine 0-24h 首次, APACHE, diabetes, BMI) ==========
DROP TABLE IF EXISTS multidatabase_glucose_validation.eicu_cov;
CREATE TABLE multidatabase_glucose_validation.eicu_cov AS
WITH main_ids AS (SELECT patientunitstayid FROM multidatabase_glucose_validation.eicu_cohort WHERE in_sensC),  -- 覆盖所有敏感性
  creat AS (
    SELECT DISTINCT ON (l.patientunitstayid) l.patientunitstayid, l.labresult AS creatinine
    FROM lab l JOIN main_ids m ON m.patientunitstayid=l.patientunitstayid
    WHERE l.labname='creatinine' AND l.labresultoffset BETWEEN 0 AND 1440 AND l.labresult IS NOT NULL
    ORDER BY l.patientunitstayid, l.labresultoffset
  ),
  apr AS (
    SELECT DISTINCT ON (patientunitstayid) patientunitstayid, apachescore, acutephysiologyscore AS aps
    FROM apachepatientresult ORDER BY patientunitstayid, apachepatientresultsid
  ),
  ph AS (
    SELECT c.patientunitstayid,
      bool_or(lower(p.pasthistorypath) LIKE '%insulin dependent diabetes%' OR lower(p.pasthistorypath) LIKE '%non-insulin dependent diabetes%') AS ph_diabetes
    FROM main_ids c LEFT JOIN pasthistory p USING (patientunitstayid) GROUP BY c.patientunitstayid
  ),
  apv AS (SELECT DISTINCT ON (patientunitstayid) patientunitstayid, diabetes FROM apachepredvar ORDER BY patientunitstayid, apachepredvarid)
SELECT m.patientunitstayid, cr.creatinine, apr.apachescore, apr.aps,
       COALESCE(ph.ph_diabetes, false) OR (apv.diabetes=1) AS diabetes,
       CASE WHEN c.admissionheight BETWEEN 100 AND 250 AND c.admissionweight BETWEEN 25 AND 400
            THEN c.admissionweight / ((c.admissionheight/100.0)^2) END AS bmi
FROM main_ids m
JOIN multidatabase_glucose_validation.eicu_cohort c USING (patientunitstayid)
LEFT JOIN creat cr USING (patientunitstayid)
LEFT JOIN apr USING (patientunitstayid)
LEFT JOIN ph USING (patientunitstayid)
LEFT JOIN apv USING (patientunitstayid);

-- ========== 4. landmark 与结局 ==========
DROP TABLE IF EXISTS multidatabase_glucose_validation.eicu_outcomes;
CREATE TABLE multidatabase_glucose_validation.eicu_outcomes AS
SELECT patientunitstayid,
  (hospitaldischargestatus='Expired') AS hosp_mortality,
  (unitdischargestatus='Expired') AS icu_mortality,
  -- landmark: 24h 内死亡 or 24h 内出院(任何状态离开医院)
  (hospitaldischargeoffset < 1440 AND hospitaldischargestatus='Expired') AS died_within_24h,
  (hospitaldischargeoffset < 1440 AND hospitaldischargestatus<>'Expired') AS discharged_within_24h,
  (hospitaldischargeoffset >= 1440) AS in_landmark,
  (hospitaldischargeoffset >= 1440 AND hospitaldischargestatus='Expired') AS post_landmark_hosp_mortality,
  (unitdischargeoffset < 1440) AS icu_los_lt_24h,
  unitdischargeoffset/1440.0 AS icu_los_days,
  (hospitaldischargeoffset-hospitaladmitoffset)/1440.0 AS hosp_los_days
FROM multidatabase_glucose_validation.eicu_cohort;

COMMIT;

-- 概览
SELECT 'main' grp, count(*) FROM multidatabase_glucose_validation.eicu_cohort WHERE in_main
UNION ALL SELECT 'sensA', count(*) FROM multidatabase_glucose_validation.eicu_cohort WHERE in_sensA
UNION ALL SELECT 'sensB', count(*) FROM multidatabase_glucose_validation.eicu_cohort WHERE in_sensB
UNION ALL SELECT 'sensC', count(*) FROM multidatabase_glucose_validation.eicu_cohort WHERE in_sensC;
SELECT procedure_category, count(*) FROM multidatabase_glucose_validation.eicu_cohort WHERE in_main GROUP BY 1;
SELECT count(*) FILTER (WHERE in_main) main_n,
  count(*) FILTER (WHERE in_main AND died_within_24h) died24,
  count(*) FILTER (WHERE in_main AND discharged_within_24h) disch24,
  count(*) FILTER (WHERE in_main AND in_landmark) landmark_n,
  count(*) FILTER (WHERE in_main AND in_landmark AND post_landmark_hosp_mortality) landmark_deaths
FROM multidatabase_glucose_validation.eicu_cohort c
JOIN multidatabase_glucose_validation.eicu_outcomes o USING (patientunitstayid);
SELECT count(*) apache_n, count(*) FILTER (WHERE apachescore IS NOT NULL) apache_avail,
  count(*) FILTER (WHERE aps IS NOT NULL) aps_avail
FROM multidatabase_glucose_validation.eicu_cov;
