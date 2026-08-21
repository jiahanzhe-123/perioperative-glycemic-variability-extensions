-- 07_extract_treatments.sql
-- Medication and support events are cohort-restricted. Text medication
-- classes are retained in a keyword review table and are not causal-model
-- defaults; they are descriptive/sensitivity variables.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.treatments_postop_v2;
DROP TABLE IF EXISTS mimic_custom.drug_keyword_review_v2;

CREATE TABLE mimic_custom.drug_keyword_review_v2 (
    medication_domain text NOT NULL,
    keyword text NOT NULL,
    source_table text NOT NULL,
    include_candidate boolean NOT NULL,
    needs_manual_review boolean NOT NULL,
    notes text,
    review_version text NOT NULL DEFAULT 'v2_drug_keyword_review_2026-07-18'
);

INSERT INTO mimic_custom.drug_keyword_review_v2
    (medication_domain, keyword, source_table, include_candidate, needs_manual_review, notes)
VALUES
  ('insulin','insulin','mimiciv_hosp.prescriptions|mimiciv_hosp.emar|mimiciv_hosp.pharmacy|mimiciv_icu.inputevents',true,false,'includes named insulin products and ICU itemids'),
  ('steroid','hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol','mimiciv_hosp.prescriptions|mimiciv_hosp.emar|mimiciv_hosp.pharmacy|mimiciv_icu.inputevents',true,true,'text class; route/topical specificity requires review'),
  ('norepinephrine','norepinephrine|levophed','all medication sources',true,false,'vasopressor'),
  ('epinephrine','epinephrine','all medication sources',true,false,'vasopressor'),
  ('vasopressin','vasopressin','all medication sources',true,false,'vasopressor'),
  ('phenylephrine','phenylephrine','all medication sources',true,false,'vasopressor'),
  ('dopamine','dopamine','all medication sources',true,false,'vasopressor/inotrope'),
  ('dobutamine','dobutamine','all medication sources',true,false,'inotrope'),
  ('milrinone','milrinone','all medication sources',true,false,'inotrope'),
  ('beta_blocker','metoprolol|labetalol|esmolol|carvedilol|atenolol|propranolol','hospital medication sources',true,true,'preop/admission medication exposure'),
  ('statin','atorvastatin|simvastatin|rosuvastatin|pravastatin|fluvastatin|lovastatin','hospital medication sources',true,true,'preop/admission medication exposure'),
  ('aspirin','aspirin|acetylsalicylic','hospital medication sources',true,true,'preop/admission medication exposure'),
  ('anticoagulant','heparin|enoxaparin|warfarin|apixaban|rivaroxaban|dabigatran','hospital medication sources',true,true,'class text may include prophylaxis'),
  ('antiplatelet','clopidogrel|ticagrelor|prasugrel|abciximab|eptifibatide|tirofiban','hospital medication sources',true,true,'antiplatelet class text')
;

CREATE TABLE mimic_custom.treatments_postop_v2 AS
WITH input_drugs AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           ie.starttime,
           ie.endtime,
           ie.itemid,
           d.label AS medication_text,
           ie.rate,
           ie.rateuom,
           ie.ordercategoryname,
           NULL::text AS route_placeholder,
           CASE
             WHEN ie.itemid IN (223257,223258,223259,223260,223261,223262,229299,229619) THEN 'insulin'
             WHEN ie.itemid IN (221906) THEN 'norepinephrine'
             WHEN ie.itemid IN (221289,229617) THEN 'epinephrine'
             WHEN ie.itemid = 222315 THEN 'vasopressin'
             WHEN ie.itemid IN (221749,229630,229631,229632) THEN 'phenylephrine'
             WHEN ie.itemid = 221662 THEN 'dopamine'
             WHEN ie.itemid = 221653 THEN 'dobutamine'
             WHEN ie.itemid = 221986 THEN 'milrinone'
             WHEN lower(d.label) ~ '(hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol)' THEN 'steroid'
           END AS medication_domain,
           'mimiciv_icu.inputevents'::text AS source_table,
           (ie.rate IS NOT NULL OR lower(COALESCE(ie.ordercategoryname,'')) ~ '(continuous|infusion|drip)') AS infusion_proxy
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_icu.inputevents ie
      ON ie.subject_id = c.subject_id
     AND ie.hadm_id = c.hadm_id
     AND ie.stay_id = c.stay_id
     AND ie.starttime < c.surgery_time + interval '7 days'
     AND COALESCE(ie.endtime, ie.starttime) >= c.surgery_time
    JOIN mimiciv_icu.d_items d ON d.itemid = ie.itemid
    WHERE ie.itemid IN (223257,223258,223259,223260,223261,223262,229299,229619,
                        221906,221289,229617,222315,221749,229630,229631,229632,
                        221662,221653,221986)
       OR lower(d.label) ~ '(hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol)'
), prescription_drugs AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           p.starttime,
           p.stoptime AS endtime,
           NULL::integer AS itemid,
           p.drug::text AS medication_text,
           NULL::double precision AS rate,
           NULL::text AS rateuom,
           NULL::text AS ordercategoryname,
           p.route::text AS route_placeholder,
           CASE
             WHEN lower(p.drug) ~ 'insulin|novolog|humalog|lispro|aspart|glargine|nph' THEN 'insulin'
             WHEN lower(p.drug) ~ 'hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol' THEN 'steroid'
             WHEN lower(p.drug) ~ 'norepinephrine|levophed' THEN 'norepinephrine'
             WHEN lower(p.drug) ~ 'epinephrine' THEN 'epinephrine'
             WHEN lower(p.drug) ~ 'vasopressin' THEN 'vasopressin'
             WHEN lower(p.drug) ~ 'phenylephrine' THEN 'phenylephrine'
             WHEN lower(p.drug) ~ 'dopamine' THEN 'dopamine'
             WHEN lower(p.drug) ~ 'dobutamine' THEN 'dobutamine'
             WHEN lower(p.drug) ~ 'milrinone' THEN 'milrinone'
             WHEN lower(p.drug) ~ 'metoprolol|labetalol|esmolol|carvedilol|atenolol|propranolol' THEN 'beta_blocker'
             WHEN lower(p.drug) ~ 'atorvastatin|simvastatin|rosuvastatin|pravastatin|fluvastatin|lovastatin' THEN 'statin'
             WHEN lower(p.drug) ~ 'aspirin|acetylsalicylic' THEN 'aspirin'
             WHEN lower(p.drug) ~ 'heparin|enoxaparin|warfarin|apixaban|rivaroxaban|dabigatran' THEN 'anticoagulant'
             WHEN lower(p.drug) ~ 'clopidogrel|ticagrelor|prasugrel|abciximab|eptifibatide|tirofiban' THEN 'antiplatelet'
           END AS medication_domain,
           'mimiciv_hosp.prescriptions'::text AS source_table,
           false AS infusion_proxy
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_hosp.prescriptions p
      ON p.subject_id = c.subject_id
     AND p.hadm_id = c.hadm_id
     AND p.starttime < c.surgery_time + interval '7 days'
     AND COALESCE(p.stoptime, p.starttime) >= c.admittime_ts
    WHERE lower(p.drug) ~ '(insulin|novolog|humalog|lispro|aspart|glargine|nph|hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol|norepinephrine|levophed|epinephrine|vasopressin|phenylephrine|dopamine|dobutamine|milrinone|metoprolol|labetalol|esmolol|carvedilol|atenolol|propranolol|atorvastatin|simvastatin|rosuvastatin|pravastatin|fluvastatin|lovastatin|aspirin|acetylsalicylic|heparin|enoxaparin|warfarin|apixaban|rivaroxaban|dabigatran|clopidogrel|ticagrelor|prasugrel|abciximab|eptifibatide|tirofiban)'
), pharmacy_drugs AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           p.starttime,
           p.stoptime AS endtime,
           NULL::integer AS itemid,
           p.medication::text AS medication_text,
           NULL::double precision AS rate,
           NULL::text AS rateuom,
           NULL::text AS ordercategoryname,
           p.route::text AS route_placeholder,
           CASE
             WHEN lower(p.medication) ~ 'insulin|novolog|humalog|lispro|aspart|glargine|nph' THEN 'insulin'
             WHEN lower(p.medication) ~ 'hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol' THEN 'steroid'
             WHEN lower(p.medication) ~ 'norepinephrine|levophed' THEN 'norepinephrine'
             WHEN lower(p.medication) ~ 'epinephrine' THEN 'epinephrine'
             WHEN lower(p.medication) ~ 'vasopressin' THEN 'vasopressin'
             WHEN lower(p.medication) ~ 'phenylephrine' THEN 'phenylephrine'
             WHEN lower(p.medication) ~ 'dopamine' THEN 'dopamine'
             WHEN lower(p.medication) ~ 'dobutamine' THEN 'dobutamine'
             WHEN lower(p.medication) ~ 'milrinone' THEN 'milrinone'
             WHEN lower(p.medication) ~ 'metoprolol|labetalol|esmolol|carvedilol|atenolol|propranolol' THEN 'beta_blocker'
             WHEN lower(p.medication) ~ 'atorvastatin|simvastatin|rosuvastatin|pravastatin|fluvastatin|lovastatin' THEN 'statin'
             WHEN lower(p.medication) ~ 'aspirin|acetylsalicylic' THEN 'aspirin'
             WHEN lower(p.medication) ~ 'heparin|enoxaparin|warfarin|apixaban|rivaroxaban|dabigatran' THEN 'anticoagulant'
             WHEN lower(p.medication) ~ 'clopidogrel|ticagrelor|prasugrel|abciximab|eptifibatide|tirofiban' THEN 'antiplatelet'
           END AS medication_domain,
           'mimiciv_hosp.pharmacy'::text AS source_table,
           false AS infusion_proxy
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_hosp.pharmacy p
      ON p.subject_id = c.subject_id
     AND p.hadm_id = c.hadm_id
     AND p.starttime < c.surgery_time + interval '7 days'
     AND COALESCE(p.stoptime, p.starttime) >= c.admittime_ts
    WHERE lower(p.medication) ~ '(insulin|novolog|humalog|lispro|aspart|glargine|nph|hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol|norepinephrine|levophed|epinephrine|vasopressin|phenylephrine|dopamine|dobutamine|milrinone|metoprolol|labetalol|esmolol|carvedilol|atenolol|propranolol|atorvastatin|simvastatin|rosuvastatin|pravastatin|fluvastatin|lovastatin|aspirin|acetylsalicylic|heparin|enoxaparin|warfarin|apixaban|rivaroxaban|dabigatran|clopidogrel|ticagrelor|prasugrel|abciximab|eptifibatide|tirofiban)'
), emar_drugs AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           e.charttime AS starttime,
           e.charttime AS endtime,
           NULL::integer AS itemid,
           e.medication::text AS medication_text,
           NULL::double precision AS rate,
           NULL::text AS rateuom,
           NULL::text AS ordercategoryname,
           NULL::text AS route_placeholder,
           CASE
             WHEN lower(e.medication) ~ 'insulin|novolog|humalog|lispro|aspart|glargine|nph' THEN 'insulin'
             WHEN lower(e.medication) ~ 'hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol' THEN 'steroid'
             WHEN lower(e.medication) ~ 'norepinephrine|levophed' THEN 'norepinephrine'
             WHEN lower(e.medication) ~ 'epinephrine' THEN 'epinephrine'
             WHEN lower(e.medication) ~ 'vasopressin' THEN 'vasopressin'
             WHEN lower(e.medication) ~ 'phenylephrine' THEN 'phenylephrine'
             WHEN lower(e.medication) ~ 'dopamine' THEN 'dopamine'
             WHEN lower(e.medication) ~ 'dobutamine' THEN 'dobutamine'
             WHEN lower(e.medication) ~ 'milrinone' THEN 'milrinone'
             WHEN lower(e.medication) ~ 'metoprolol|labetalol|esmolol|carvedilol|atenolol|propranolol' THEN 'beta_blocker'
             WHEN lower(e.medication) ~ 'atorvastatin|simvastatin|rosuvastatin|pravastatin|fluvastatin|lovastatin' THEN 'statin'
             WHEN lower(e.medication) ~ 'aspirin|acetylsalicylic' THEN 'aspirin'
             WHEN lower(e.medication) ~ 'heparin|enoxaparin|warfarin|apixaban|rivaroxaban|dabigatran' THEN 'anticoagulant'
             WHEN lower(e.medication) ~ 'clopidogrel|ticagrelor|prasugrel|abciximab|eptifibatide|tirofiban' THEN 'antiplatelet'
           END AS medication_domain,
           'mimiciv_hosp.emar'::text AS source_table,
           false AS infusion_proxy
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_hosp.emar e
      ON e.subject_id = c.subject_id
     AND e.hadm_id = c.hadm_id
     AND e.charttime >= c.admittime_ts
     AND e.charttime < c.surgery_time + interval '7 days'
    WHERE lower(e.medication) ~ '(insulin|novolog|humalog|lispro|aspart|glargine|nph|hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol|norepinephrine|levophed|epinephrine|vasopressin|phenylephrine|dopamine|dobutamine|milrinone|metoprolol|labetalol|esmolol|carvedilol|atenolol|propranolol|atorvastatin|simvastatin|rosuvastatin|pravastatin|fluvastatin|lovastatin|aspirin|acetylsalicylic|heparin|enoxaparin|warfarin|apixaban|rivaroxaban|dabigatran|clopidogrel|ticagrelor|prasugrel|abciximab|eptifibatide|tirofiban)'
), medication_events AS (
    SELECT * FROM input_drugs
    UNION ALL SELECT * FROM prescription_drugs
    UNION ALL SELECT * FROM pharmacy_drugs
    UNION ALL SELECT * FROM emar_drugs
), vent_events AS (
    SELECT c.subject_id, c.hadm_id, c.stay_id, pe.starttime, pe.endtime, pe.itemid
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_icu.procedureevents pe
      ON pe.subject_id = c.subject_id
     AND pe.hadm_id = c.hadm_id
     AND pe.stay_id = c.stay_id
     AND pe.starttime < c.surgery_time + interval '7 days'
     AND COALESCE(pe.endtime, pe.starttime) >= c.surgery_time
    WHERE pe.itemid IN (225792,225794,224385)
), rrt_events AS (
    SELECT c.subject_id, c.hadm_id, c.stay_id, pe.starttime, pe.endtime, pe.itemid
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_icu.procedureevents pe
      ON pe.subject_id = c.subject_id
     AND pe.hadm_id = c.hadm_id
     AND pe.stay_id = c.stay_id
     AND pe.starttime < c.surgery_time + interval '7 days'
     AND COALESCE(pe.endtime, pe.starttime) >= c.surgery_time
    WHERE pe.itemid IN (225802,225803,225809,225955,225805,225441)
), drug_flags AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           count(m.stay_id)::bigint AS medication_event_n,
           bool_or(m.medication_domain = 'insulin' AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS insulin_any_postop_24h,
           bool_or(m.medication_domain = 'insulin' AND m.starttime < c.surgery_time + interval '48 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS insulin_any_postop_48h,
           bool_or(m.medication_domain = 'insulin' AND m.starttime < c.surgery_time + interval '7 days' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS insulin_any_postop_7d,
           bool_or(m.medication_domain = 'insulin' AND m.infusion_proxy AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS insulin_infusion_postop_24h,
           bool_or(m.medication_domain = 'insulin' AND lower(COALESCE(m.route_placeholder,'')) ~ '(subcutaneous|subq|\bsc\b)' AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS insulin_subq_postop_24h,
           bool_or(m.medication_domain = 'steroid' AND m.starttime < c.surgery_time + interval '48 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS steroid_any_postop_48h,
           bool_or(m.medication_domain = 'steroid' AND m.starttime < c.surgery_time + interval '7 days' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS steroid_any_postop_7d,
           bool_or(m.medication_domain IN ('norepinephrine','epinephrine','vasopressin','phenylephrine','dopamine') AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS vasopressor_any_postop_24h,
           bool_or(m.medication_domain IN ('norepinephrine','epinephrine','vasopressin','phenylephrine','dopamine') AND m.starttime < c.surgery_time + interval '48 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS vasopressor_any_postop_48h,
           bool_or(m.medication_domain = 'norepinephrine' AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS norepinephrine_any_postop_24h,
           bool_or(m.medication_domain = 'epinephrine' AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS epinephrine_any_postop_24h,
           bool_or(m.medication_domain = 'vasopressin' AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS vasopressin_any_postop_24h,
           bool_or(m.medication_domain = 'phenylephrine' AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS phenylephrine_any_postop_24h,
           bool_or(m.medication_domain = 'dobutamine' AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS dobutamine_any_postop_24h,
           bool_or(m.medication_domain = 'milrinone' AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS milrinone_any_postop_24h,
           bool_or(m.medication_domain IN ('dobutamine','milrinone') AND m.starttime < c.surgery_time + interval '24 hours' AND COALESCE(m.endtime,m.starttime) >= c.surgery_time) AS inotrope_any_postop_24h,
           bool_or(m.medication_domain = 'beta_blocker' AND m.starttime >= c.admittime_ts AND m.starttime <= c.surgery_time) AS beta_blocker_preop_or_adm_flag,
           bool_or(m.medication_domain = 'statin' AND m.starttime >= c.admittime_ts AND m.starttime <= c.surgery_time) AS statin_preop_or_adm_flag,
           bool_or(m.medication_domain = 'aspirin' AND m.starttime >= c.admittime_ts AND m.starttime <= c.surgery_time) AS aspirin_preop_or_adm_flag,
           bool_or(m.medication_domain = 'anticoagulant' AND m.starttime >= c.admittime_ts AND m.starttime <= c.surgery_time) AS anticoagulant_adm_flag,
           bool_or(m.medication_domain = 'antiplatelet' AND m.starttime >= c.admittime_ts AND m.starttime <= c.surgery_time) AS antiplatelet_adm_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN medication_events m ON m.stay_id = c.stay_id
    GROUP BY c.subject_id, c.hadm_id, c.stay_id
), vent_flags AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           count(v.stay_id)::bigint AS vent_event_n,
           bool_or(v.itemid IN (225792,224385) AND v.starttime < c.surgery_time + interval '24 hours' AND COALESCE(v.endtime,v.starttime) >= c.surgery_time) AS mechanical_ventilation_postop_24h_flag,
           bool_or(v.itemid IN (225792,224385) AND v.starttime < c.surgery_time + interval '48 hours' AND COALESCE(v.endtime,v.starttime) >= c.surgery_time) AS mechanical_ventilation_postop_48h_flag,
           bool_or(v.itemid = 225792 AND v.endtime IS NULL) AS vent_duration_unreliable,
           CASE WHEN count(v.stay_id) = 0 OR bool_or(v.itemid = 225792 AND v.endtime IS NULL) THEN NULL
                ELSE sum(GREATEST(0, EXTRACT(EPOCH FROM (LEAST(v.endtime, c.surgery_time + interval '7 days') - GREATEST(v.starttime, c.surgery_time))) / 3600.0) ) FILTER (WHERE v.itemid = 225792) END AS vent_duration_hours
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN vent_events v ON v.stay_id = c.stay_id
    GROUP BY c.subject_id, c.hadm_id, c.stay_id
), rrt_flags AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           count(r.stay_id)::bigint AS rrt_event_n,
           bool_or(r.stay_id IS NOT NULL) AS rrt_flag,
           bool_or(r.stay_id IS NOT NULL AND r.starttime < c.surgery_time + interval '7 days' AND COALESCE(r.endtime,r.starttime) >= c.surgery_time) AS rrt_postop_7d_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN rrt_events r ON r.stay_id = c.stay_id
    GROUP BY c.subject_id, c.hadm_id, c.stay_id
), glucose_flags AS (
    SELECT c.stay_id,
           bool_or(g.glucose_mg_dl < 70) AS hypoglycemia_postop_24h_flag,
           bool_or(g.glucose_mg_dl >= 180) AS hyperglycemia_180_postop_24h_flag,
           bool_or(g.glucose_mg_dl >= 200) AS hyperglycemia_200_postop_24h_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN mimic_custom.glucose_long_v2 g
      ON g.stay_id = c.stay_id
     AND g.charttime >= c.surgery_time
     AND g.charttime < c.surgery_time + interval '24 hours'
    GROUP BY c.stay_id
), lab_flags AS (
    SELECT c.stay_id,
           bool_or(l.variable_name = 'lactate' AND l.value_valid >= 2) AS high_lactate_postop_24h_flag,
           bool_or(l.variable_name = 'ph' AND l.value_valid < 7.2) AS severe_acidosis_postop_24h_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN mimic_custom.labs_long_v2 l
      ON l.stay_id = c.stay_id
     AND l.charttime >= c.surgery_time
     AND l.charttime < c.surgery_time + interval '24 hours'
     AND l.variable_name IN ('lactate','ph')
     AND l.value_valid IS NOT NULL
    GROUP BY c.stay_id
), support_flags AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           gf.hypoglycemia_postop_24h_flag,
           gf.hyperglycemia_180_postop_24h_flag,
           gf.hyperglycemia_200_postop_24h_flag,
           lf.high_lactate_postop_24h_flag,
           lf.severe_acidosis_postop_24h_flag,
           (v.spo2_min_postop_24h < 90) AS hypoxemia_postop_24h_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN glucose_flags gf ON gf.stay_id = c.stay_id
    LEFT JOIN lab_flags lf ON lf.stay_id = c.stay_id
    LEFT JOIN mimic_custom.vitals_postop_24h_v2 v ON v.stay_id = c.stay_id
), vaso_duration AS (
    SELECT c.stay_id,
           CASE WHEN bool_or(ie.endtime IS NULL) THEN NULL
                ELSE sum(GREATEST(0, EXTRACT(EPOCH FROM (LEAST(ie.endtime, c.surgery_time + interval '7 days') - GREATEST(ie.starttime, c.surgery_time))) / 3600.0) ) END AS vaso_duration_hours
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_icu.inputevents ie
      ON ie.stay_id = c.stay_id
     AND ie.starttime < c.surgery_time + interval '7 days'
     AND COALESCE(ie.endtime, ie.starttime) >= c.surgery_time
     AND ie.itemid IN (221906,221289,229617,222315,221749,229630,229631,229632,221662)
    GROUP BY c.stay_id, c.surgery_time
)
SELECT c.subject_id,
       c.hadm_id,
       c.stay_id,
       COALESCE(df.medication_event_n,0) AS medication_event_n,
       COALESCE(vf.vent_event_n,0) AS vent_event_n,
       COALESCE(rf.rrt_event_n,0) AS rrt_event_n,
       COALESCE(df.insulin_any_postop_24h,false) AS insulin_any_postop_24h,
       COALESCE(df.insulin_any_postop_48h,false) AS insulin_any_postop_48h,
       COALESCE(df.insulin_any_postop_7d,false) AS insulin_any_postop_7d,
       COALESCE(df.insulin_infusion_postop_24h,false) AS insulin_infusion_postop_24h,
       COALESCE(df.insulin_subq_postop_24h,false) AS insulin_subq_postop_24h,
       COALESCE(df.steroid_any_postop_48h,false) AS steroid_any_postop_48h,
       COALESCE(df.steroid_any_postop_7d,false) AS steroid_any_postop_7d,
       COALESCE(df.vasopressor_any_postop_24h,false) AS vasopressor_any_postop_24h,
       COALESCE(df.vasopressor_any_postop_48h,false) AS vasopressor_any_postop_48h,
       COALESCE(df.norepinephrine_any_postop_24h,false) AS norepinephrine_any_postop_24h,
       COALESCE(df.epinephrine_any_postop_24h,false) AS epinephrine_any_postop_24h,
       COALESCE(df.vasopressin_any_postop_24h,false) AS vasopressin_any_postop_24h,
       COALESCE(df.phenylephrine_any_postop_24h,false) AS phenylephrine_any_postop_24h,
       COALESCE(df.dobutamine_any_postop_24h,false) AS dobutamine_any_postop_24h,
       COALESCE(df.milrinone_any_postop_24h,false) AS milrinone_any_postop_24h,
       COALESCE(df.inotrope_any_postop_24h,false) AS inotrope_any_postop_24h,
       COALESCE(df.beta_blocker_preop_or_adm_flag,false) AS beta_blocker_preop_or_adm_flag,
       COALESCE(df.statin_preop_or_adm_flag,false) AS statin_preop_or_adm_flag,
       COALESCE(df.aspirin_preop_or_adm_flag,false) AS aspirin_preop_or_adm_flag,
       COALESCE(df.anticoagulant_adm_flag,false) AS anticoagulant_adm_flag,
       COALESCE(df.antiplatelet_adm_flag,false) AS antiplatelet_adm_flag,
       COALESCE(vf.mechanical_ventilation_postop_24h_flag,false) AS mechanical_ventilation_postop_24h_flag,
       COALESCE(vf.mechanical_ventilation_postop_48h_flag,false) AS mechanical_ventilation_postop_48h_flag,
       vf.vent_duration_hours,
       COALESCE(vf.vent_duration_unreliable,false) AS vent_duration_unreliable,
       COALESCE(rf.rrt_flag,false) AS rrt_flag,
       COALESCE(rf.rrt_postop_7d_flag,false) AS rrt_postop_7d_flag,
       sf.high_lactate_postop_24h_flag,
       sf.severe_acidosis_postop_24h_flag,
       sf.hypoxemia_postop_24h_flag,
       sf.hypoglycemia_postop_24h_flag,
       sf.hyperglycemia_180_postop_24h_flag,
       sf.hyperglycemia_200_postop_24h_flag,
       vd.vaso_duration_hours,
       vf.vent_duration_hours AS total_vent_duration_hours,
       (vf.vent_duration_hours IS NOT NULL) AS vent_duration_reliable,
       'inputevents + prescriptions + pharmacy + emar + procedureevents; text/ID review table retained'::text AS treatment_source,
       CASE WHEN COALESCE(df.medication_event_n,0) = 0
                  AND COALESCE(vf.vent_event_n,0) = 0
                  AND COALESCE(rf.rrt_event_n,0) = 0
            THEN 'no matching treatment/support event rows' ELSE NULL END AS treatment_missing_reason
FROM mimic_custom.cardiac_surgery_cohort_v2 c
LEFT JOIN drug_flags df ON df.stay_id = c.stay_id
LEFT JOIN vent_flags vf ON vf.stay_id = c.stay_id
LEFT JOIN rrt_flags rf ON rf.stay_id = c.stay_id
LEFT JOIN support_flags sf ON sf.stay_id = c.stay_id
LEFT JOIN vaso_duration vd ON vd.stay_id = c.stay_id;

CREATE INDEX treatments_postop_v2_stay_idx
    ON mimic_custom.treatments_postop_v2 (stay_id);

COMMIT;

\copy (SELECT * FROM mimic_custom.drug_keyword_review_v2 ORDER BY medication_domain) TO 'outputs/qc/drug_keyword_review_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
