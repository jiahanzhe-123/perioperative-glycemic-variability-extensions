-- Step 6b: 协变量
BEGIN;
DROP TABLE IF EXISTS eicu_cardio_validation.covariates;
CREATE TABLE eicu_cardio_validation.covariates AS
WITH c AS (
  SELECT patientunitstayid FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND rn_patient=1
),
lab_first AS (
  SELECT patientunitstayid, labname, labresult
  FROM (
    SELECT l.patientunitstayid, l.labname, l.labresult,
           row_number() OVER (PARTITION BY l.patientunitstayid, l.labname ORDER BY l.labresultoffset) rn
    FROM lab l JOIN c USING (patientunitstayid)
    WHERE l.labname IN ('lactate','creatinine','Hgb','WBC x 1000','platelets x 1000','albumin')
      AND l.labresultoffset BETWEEN 0 AND 1440 AND l.labresult IS NOT NULL
  ) x WHERE rn=1
),
apr AS (  -- apachepatientresult 去重: 每 stay 一行
  SELECT DISTINCT ON (patientunitstayid) patientunitstayid, apachescore, acutephysiologyscore, apacheversion, predictedhospitalmortality
  FROM apachepatientresult ORDER BY patientunitstayid, apachepatientresultsid
),
apv AS (  -- apachepredvar 去重
  SELECT DISTINCT ON (patientunitstayid) patientunitstayid, electivesurgery, diabetes
  FROM apachepredvar ORDER BY patientunitstayid, apachepredvarid
),
aps AS (  -- apacheapsvar 去重
  SELECT DISTINCT ON (patientunitstayid) patientunitstayid, intubated, vent
  FROM apacheapsvar ORDER BY patientunitstayid, apacheapsvarid
),
lab_pivot AS (
  SELECT patientunitstayid,
         max(labresult) FILTER (WHERE labname='lactate') AS lactate,
         max(labresult) FILTER (WHERE labname='creatinine') AS creatinine,
         max(labresult) FILTER (WHERE labname='Hgb') AS hemoglobin,
         max(labresult) FILTER (WHERE labname='WBC x 1000') AS wbc,
         max(labresult) FILTER (WHERE labname='platelets x 1000') AS platelets,
         max(labresult) FILTER (WHERE labname='albumin') AS albumin
  FROM lab_first GROUP BY patientunitstayid
),
ph AS (
  SELECT c.patientunitstayid,
    bool_or(lower(pasthistorypath) LIKE '%insulin dependent diabetes%' OR lower(pasthistorypath) LIKE '%non-insulin dependent diabetes%') AS ph_diabetes,
    bool_or(lower(pasthistorypath) LIKE '%renal failure%' OR lower(pasthistorypath) LIKE '%renal insufficiency%') AS ph_renal,
    bool_or(lower(pasthistorypath) LIKE '%congestive heart failure%') AS ph_chf,
    bool_or(lower(pasthistorypath) LIKE '%copd%' OR lower(pasthistorypath) LIKE '%home oxygen%') AS ph_pulmonary,
    bool_or(lower(pasthistorypath) LIKE '%cirrhosis%') AS ph_liver,
    bool_or(lower(pasthistorypath) LIKE '%hypertension requiring treatment%') AS ph_hypertension,
    bool_or(lower(pasthistorypath) LIKE '%atrial fibrillation%') AS ph_afib
  FROM c LEFT JOIN pasthistory p USING (patientunitstayid)
  GROUP BY c.patientunitstayid
),
inf AS (
  SELECT c.patientunitstayid,
    bool_or(lower(drugname) ~ 'norepinephrine|epinephrine|dopamine|vasopressin|phenylephrine') AS vasopressor,
    bool_or(lower(drugname) LIKE 'insulin%') AS insulin_infusion
  FROM c LEFT JOIN infusiondrug i USING (patientunitstayid)
  GROUP BY c.patientunitstayid
)
SELECT c.patientunitstayid,
  r.apachescore, r.acutephysiologyscore AS aps, r.apacheversion,
  CASE WHEN r.predictedhospitalmortality ~ '^[0-9.]+$' THEN r.predictedhospitalmortality::numeric END AS predicted_hosp_mortality,
  (v.electivesurgery=1) AS elective_surgery,
  COALESCE((a.intubated=1) OR (a.vent=1), false) OR EXISTS(SELECT 1 FROM respiratorycare rc WHERE rc.patientunitstayid=c.patientunitstayid) AS mechanical_ventilation,
  (v.diabetes=1) AS apache_diabetes,
  lp.lactate, lp.creatinine, lp.hemoglobin, lp.wbc, lp.platelets, lp.albumin,
  ph.ph_diabetes, ph.ph_renal, ph.ph_chf, ph.ph_pulmonary, ph.ph_liver, ph.ph_hypertension, ph.ph_afib,
  inf.vasopressor, inf.insulin_infusion
FROM c
LEFT JOIN apr r USING (patientunitstayid)
LEFT JOIN apv v USING (patientunitstayid)
LEFT JOIN aps a USING (patientunitstayid)
LEFT JOIN lab_pivot lp USING (patientunitstayid)
LEFT JOIN ph USING (patientunitstayid)
LEFT JOIN inf USING (patientunitstayid);
CREATE UNIQUE INDEX ON eicu_cardio_validation.covariates(patientunitstayid);
COMMIT;

SELECT count(*) n, count(*) FILTER (WHERE lactate IS NOT NULL) lactate, count(*) FILTER (WHERE creatinine IS NOT NULL) creat,
 count(*) FILTER (WHERE albumin IS NOT NULL) albumin, count(*) FILTER (WHERE apachescore IS NOT NULL) apache,
 count(*) FILTER (WHERE mechanical_ventilation) vent, count(*) FILTER (WHERE vasopressor) vaso,
 count(*) FILTER (WHERE insulin_infusion) insulin, count(*) FILTER (WHERE ph_diabetes OR apache_diabetes) dm
FROM eicu_cardio_validation.covariates;
