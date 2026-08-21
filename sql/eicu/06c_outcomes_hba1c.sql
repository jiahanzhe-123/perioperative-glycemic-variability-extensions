-- Step 6c: 结局 + HbA1c 覆盖核查
BEGIN;
DROP TABLE IF EXISTS eicu_cardio_validation.outcomes;
CREATE TABLE eicu_cardio_validation.outcomes AS
SELECT patientunitstayid,
       (hospitaldischargestatus = 'Expired') AS hospital_mortality,
       (unitdischargestatus = 'Expired') AS icu_mortality,
       CASE WHEN hospitaldischargestatus IN ('Alive','Expired') THEN 0 ELSE 1 END AS hosp_status_missing,
       unitdischargeoffset / 1440.0 AS icu_los_days,
       (hospitaldischargeoffset - hospitaladmitoffset) / 1440.0 AS hospital_los_days,
       hospitaldischargeoffset AS hosp_discharge_offset_min,
       unitdischargeoffset AS icu_discharge_offset_min
FROM eicu_cardio_validation.cohort_base
WHERE broad_cohort_flag AND rn_patient=1;
CREATE UNIQUE INDEX ON eicu_cardio_validation.outcomes(patientunitstayid);

-- HbA1c: lab 表无; 仅 customlab 有个位数记录
DROP TABLE IF EXISTS eicu_cardio_validation.hba1c_records;
CREATE TABLE eicu_cardio_validation.hba1c_records AS
SELECT patientunitstayid, labotheroffset AS offset_min, labothername, labotherresult, labothervaluetext
FROM customlab
WHERE lower(labothername) ~ 'a1c|glyc';
COMMIT;

SELECT count(*) FILTER (WHERE hospital_mortality) AS hosp_deaths,
       count(*) FILTER (WHERE icu_mortality) AS icu_deaths,
       count(*) AS n,
       sum(hosp_status_missing) AS hosp_status_missing,
       round(avg(icu_los_days)::numeric,2) AS mean_icu_los,
       round(avg(hospital_los_days)::numeric,2) AS mean_hosp_los
FROM eicu_cardio_validation.outcomes;
SELECT * FROM eicu_cardio_validation.hba1c_records;
SELECT count(*) AS hba1c_in_cohort FROM eicu_cardio_validation.hba1c_records h
JOIN eicu_cardio_validation.cohort_base c USING (patientunitstayid) WHERE c.broad_cohort_flag AND c.rn_patient=1;
