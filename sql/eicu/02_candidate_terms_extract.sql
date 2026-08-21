-- Step 2: 提取心脏相关候选字符串(4个来源) -> /tmp/cand_terms_raw.csv
\pset footer off
\o /tmp/cand_terms_raw.csv
WITH src AS (
  SELECT 'patient' AS source_table, 'apacheadmissiondx' AS source_column,
         apacheadmissiondx AS original_value, count(DISTINCT patientunitstayid) AS frequency
  FROM patient
  WHERE lower(apacheadmissiondx) ~ 'cabg|bypass|valve|aort|coronary|cardiac|heart|transplant|assist device|septal|congenital|aneurysm|dissect|stent|cath|angio|tavr|ablat|pacemaker|defibrill|endarterect|sternot|thorac'
  GROUP BY 3
  UNION ALL
  SELECT 'admissiondx', 'admitdxpath', admitdxpath, count(DISTINCT patientunitstayid)
  FROM admissiondx
  WHERE lower(admitdxpath) ~ 'cabg|bypass|valve|aort|coronary|cardiac|heart|transplant|assist device|septal|congenital|aneurysm|dissect|stent|cath|angio|tavr|ablat|pacemaker|defibrill|endarterect|sternot|thorac'
  GROUP BY 3
  UNION ALL
  SELECT 'admissiondx', 'admitdxname', admitdxname, count(DISTINCT patientunitstayid)
  FROM admissiondx
  WHERE lower(admitdxname) ~ 'cabg|bypass|valve|aort|coronary|cardiac|heart|transplant|assist device|septal|congenital|aneurysm|dissect|stent|cath|angio|tavr|ablat|pacemaker|defibrill|endarterect|sternot|thorac'
  GROUP BY 3
  UNION ALL
  SELECT 'diagnosis', 'diagnosisstring', diagnosisstring, count(DISTINCT patientunitstayid)
  FROM diagnosis
  WHERE lower(diagnosisstring) ~ 'cabg|bypass|valve|aort|coronary|cardiac|heart|transplant|assist device|septal|congenital|aneurysm|dissect|stent|cath|angio|tavr|ablat|pacemaker|defibrill|endarterect|sternot|thorac'
  GROUP BY 3
  UNION ALL
  SELECT 'treatment', 'treatmentstring', treatmentstring, count(DISTINCT patientunitstayid)
  FROM treatment
  WHERE lower(treatmentstring) ~ 'cabg|bypass|valve|aort|coronary|cardiac|heart|transplant|assist device|septal|congenital|aneurysm|dissect|stent|cath|angio|tavr|ablat|pacemaker|defibrill|endarterect|sternot|thorac'
  GROUP BY 3
)
SELECT * FROM src ORDER BY source_table, frequency DESC;
\o
