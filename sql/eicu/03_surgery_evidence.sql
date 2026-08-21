-- Step 3: stay 级心脏手术证据表
-- 证据层级:
--   apache_dx  : patient.apacheadmissiondx 命中 INCLUDE 术语(本次住院的 APACHE 手术诊断)
--   admitdx    : admissiondx 中 Operative|Cardiovascular 路径命中 INCLUDE 术语
--   dx_text    : diagnosis 的 "cardiovascular|cardiac surgery|" 子树(限近期:<7天/未标时间的术后状态;>=7天仅作弱支持)
--   tx_text    : treatment 的 "cardiovascular|cardiac surgery|" 子树(本 ICU 期间实施的手术)
-- 支持条件: unitadmitsource in (OR/RR/PACU) | electivesurgery=1 | 心脏相关 unittype | apachepredvar.admitsource
-- 主类别优先级: TRANSPLANT_VAD > OPEN_AORTIC > CONGENITAL > OPEN_VALVE > CABG > OTHER_OPEN

BEGIN;
DROP TABLE IF EXISTS eicu_cardio_validation.surgery_term_hits;
CREATE TABLE eicu_cardio_validation.surgery_term_hits AS
-- 1) apache_dx
SELECT p.patientunitstayid, 'apache_dx'::text AS evidence_source, t.proposed_category,
       p.apacheadmissiondx AS matched_term, 1 AS term_strength  -- 1=强(本次住院手术诊断)
FROM patient p
JOIN eicu_cardio_validation.cand_terms_classified t
  ON t.source_table='patient' AND t.source_column='apacheadmissiondx'
 AND t.original_value = p.apacheadmissiondx AND t.proposed_category IN ('CABG','OPEN_VALVE','OPEN_AORTIC','TRANSPLANT_VAD','CONGENITAL','OTHER_OPEN')
UNION ALL
-- 2) admitdx(仅 Operative|Cardiovascular 路径)
SELECT a.patientunitstayid, 'admitdx', t.proposed_category, a.admitdxname, 1
FROM admissiondx a
JOIN eicu_cardio_validation.cand_terms_classified t
  ON t.source_table='admissiondx' AND t.source_column='admitdxname'
 AND t.original_value = a.admitdxname AND t.proposed_category IN ('CABG','OPEN_VALVE','OPEN_AORTIC','TRANSPLANT_VAD','CONGENITAL','OTHER_OPEN')
WHERE a.admitdxpath LIKE 'admission diagnosis|All Diagnosis|Operative|Diagnosis|Cardiovascular|%'
UNION ALL
-- 3) dx_text(cardiac surgery 子树; ">= 7 days" 为弱证据 strength=2)
SELECT d.patientunitstayid, 'dx_text', t.proposed_category, d.diagnosisstring,
       CASE WHEN lower(d.diagnosisstring) LIKE '%>= 7 days%' THEN 2 ELSE 1 END
FROM diagnosis d
JOIN eicu_cardio_validation.cand_terms_classified t
  ON t.source_table='diagnosis' AND t.source_column='diagnosisstring'
 AND t.original_value = d.diagnosisstring AND t.proposed_category IN ('CABG','OPEN_VALVE','OPEN_AORTIC','TRANSPLANT_VAD','CONGENITAL','OTHER_OPEN')
WHERE lower(d.diagnosisstring) LIKE 'cardiovascular|cardiac surgery|%'
UNION ALL
-- 4) tx_text(cardiac surgery 子树)
SELECT tr.patientunitstayid, 'tx_text', t.proposed_category, tr.treatmentstring, 1
FROM treatment tr
JOIN eicu_cardio_validation.cand_terms_classified t
  ON t.source_table='treatment' AND t.source_column='treatmentstring'
 AND t.original_value = tr.treatmentstring AND t.proposed_category IN ('CABG','OPEN_VALVE','OPEN_AORTIC','TRANSPLANT_VAD','CONGENITAL','OTHER_OPEN')
WHERE lower(tr.treatmentstring) LIKE 'cardiovascular|cardiac surgery|%';

CREATE INDEX ON eicu_cardio_validation.surgery_term_hits(patientunitstayid);

-- 支持条件
DROP TABLE IF EXISTS eicu_cardio_validation.surgery_support;
CREATE TABLE eicu_cardio_validation.surgery_support AS
SELECT p.patientunitstayid,
       (p.unitadmitsource IN ('Operating Room','Recovery Room','PACU')) AS sup_or_rr,
       (v.electivesurgery = 1) AS sup_elective,
       (p.unittype IN ('CSICU','CTICU','CCU-CTICU','Cardiac ICU')) AS sup_cardiac_unit
FROM patient p
LEFT JOIN apachepredvar v USING (patientunitstayid);

-- stay 级汇总
DROP TABLE IF EXISTS eicu_cardio_validation.surgery_evidence;
CREATE TABLE eicu_cardio_validation.surgery_evidence AS
WITH hits AS (
  SELECT patientunitstayid,
         array_agg(DISTINCT proposed_category) AS cats,
         array_agg(DISTINCT evidence_source) AS sources,
         string_agg(DISTINCT evidence_source || ':' || proposed_category || ':' || left(matched_term,80), ' || ') AS evidence_detail,
         bool_or(evidence_source IN ('apache_dx','admitdx') AND term_strength=1) AS has_primary_dx,
         bool_or(evidence_source='tx_text' AND term_strength=1) AS has_tx,
         bool_or(evidence_source='dx_text' AND term_strength=1) AS has_dx_recent,
         count(DISTINCT evidence_source) AS n_sources
  FROM eicu_cardio_validation.surgery_term_hits
  GROUP BY patientunitstayid
)
SELECT h.patientunitstayid,
       h.cats AS surgery_categories,
       -- 主类别(按预定优先级)
       CASE
         WHEN 'TRANSPLANT_VAD' = ANY(h.cats) THEN 'TRANSPLANT_VAD'
         WHEN 'OPEN_AORTIC'   = ANY(h.cats) THEN 'OPEN_AORTIC'
         WHEN 'CONGENITAL'    = ANY(h.cats) THEN 'CONGENITAL'
         WHEN 'OPEN_VALVE'    = ANY(h.cats) THEN 'OPEN_VALVE'
         WHEN 'CABG'          = ANY(h.cats) THEN 'CABG'
         ELSE 'OTHER_OPEN'
       END AS surgery_category,
       (array_length(h.cats,1) > 1) AS multi_procedure_flag,
       h.n_sources AS evidence_source_count,
       h.evidence_detail AS evidence_summary,
       s.sup_or_rr, s.sup_elective, s.sup_cardiac_unit,
       ((COALESCE(s.sup_or_rr,false))::int + (COALESCE(s.sup_elective,false))::int
        + (COALESCE(s.sup_cardiac_unit,false))::int
        + (h.has_tx)::int + (h.has_dx_recent)::int) AS support_count,
       CASE
         WHEN h.has_primary_dx THEN 'DEFINITE'
         WHEN (h.has_tx OR h.has_dx_recent) AND (COALESCE(s.sup_or_rr,false) OR COALESCE(s.sup_elective,false) OR COALESCE(s.sup_cardiac_unit,false)) THEN 'PROBABLE'
         WHEN (h.has_tx OR h.has_dx_recent) THEN 'POSSIBLE'
         ELSE 'POSSIBLE'
       END AS surgery_confidence
FROM hits h
LEFT JOIN eicu_cardio_validation.surgery_support s USING (patientunitstayid);

CREATE UNIQUE INDEX ON eicu_cardio_validation.surgery_evidence(patientunitstayid);
COMMIT;

SELECT surgery_confidence, count(*) FROM eicu_cardio_validation.surgery_evidence GROUP BY 1 ORDER BY 2 DESC;
SELECT surgery_category, surgery_confidence, count(*) FROM eicu_cardio_validation.surgery_evidence GROUP BY 1,2 ORDER BY 1,2;
