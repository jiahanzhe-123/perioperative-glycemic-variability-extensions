-- Step 4: 基础队列 + 去重 + cohort flow
BEGIN;
DROP TABLE IF EXISTS eicu_cardio_validation.cohort_base;
CREATE TABLE eicu_cardio_validation.cohort_base AS
SELECT p.patientunitstayid, p.patienthealthsystemstayid, p.uniquepid, p.hospitalid,
       p.gender, p.age, p.ethnicity, p.unittype, p.unitadmitsource, p.unitstaytype, p.unitvisitnumber,
       p.hospitaladmitoffset, p.unitdischargeoffset, p.hospitaldischargeoffset,
       p.unitdischargestatus, p.hospitaldischargestatus, p.admissionheight, p.admissionweight, p.dischargeweight,
       e.surgery_categories, e.surgery_category, e.multi_procedure_flag,
       e.evidence_source_count, e.evidence_summary, e.support_count, e.surgery_confidence,
       (e.surgery_confidence='DEFINITE' AND e.support_count >= 1) AS high_specificity_flag,
       (e.surgery_confidence IN ('DEFINITE','PROBABLE')) AS broad_cohort_flag,
       (e.surgery_category IN ('CABG','OPEN_VALVE','OPEN_AORTIC')) AS open_core_flag,
       (e.surgery_category IN ('CABG','OPEN_VALVE','OPEN_AORTIC','TRANSPLANT_VAD','CONGENITAL','OTHER_OPEN')) AS open_extended_flag,
       (CASE WHEN p.age ~ '^[0-9]+$' THEN p.age::int WHEN p.age = '> 89' THEN 90 ELSE NULL END) >= 18 AS adult_flag,
       row_number() OVER (PARTITION BY p.uniquepid ORDER BY p.unitvisitnumber, p.patientunitstayid) AS rn_patient
FROM eicu_cardio_validation.surgery_evidence e
JOIN patient p USING (patientunitstayid);
CREATE UNIQUE INDEX ON eicu_cardio_validation.cohort_base(patientunitstayid);
COMMIT;

-- flow
WITH f AS (
  SELECT
    (SELECT count(*) FROM eicu_cardio_validation.cohort_base) AS c1_any_evidence,
    (SELECT count(*) FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag) AS c2_definite_probable,
    (SELECT count(*) FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND adult_flag) AS c3_adult,
    (SELECT count(*) FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND adult_flag AND rn_patient=1) AS c4_first_stay,
    (SELECT count(*) FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND adult_flag AND rn_patient=1 AND open_core_flag) AS c5_open_core,
    (SELECT count(*) FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND adult_flag AND rn_patient=1 AND open_extended_flag) AS c6_open_extended,
    (SELECT count(*) FROM eicu_cardio_validation.cohort_base WHERE high_specificity_flag AND adult_flag AND rn_patient=1) AS c7_high_spec_first_stay,
    (SELECT count(DISTINCT hospitalid) FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND adult_flag AND rn_patient=1) AS hospitals
)
SELECT * FROM f;
-- 按信心级别看 flow
SELECT surgery_confidence, adult_flag, count(*) AS stays, count(*) FILTER (WHERE rn_patient=1) AS first_stays
FROM eicu_cardio_validation.cohort_base GROUP BY 1,2 ORDER BY 1,2;
