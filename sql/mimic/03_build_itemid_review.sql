-- 03_build_itemid_review.sql
-- The curated IDs below are explicit, reviewable defaults. Broad label matches
-- are retained as manual-review rows but are not used by extraction by default.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.itemid_review_v2;
CREATE TABLE mimic_custom.itemid_review_v2 (
    variable_name text NOT NULL,
    source_table text NOT NULL,
    itemid integer,
    label text,
    fluid_or_category_or_unit text,
    include_candidate boolean NOT NULL,
    exclude_reason text,
    match_keyword text,
    needs_manual_review boolean NOT NULL,
    review_status text NOT NULL,
    mapping_version text NOT NULL DEFAULT 'v2_explicit_itemid_map_2026-07-18'
);

-- ICU charted items, input items, and procedure items.
WITH candidates AS (
    SELECT CASE
             WHEN d.itemid = 220045 THEN 'heart_rate'
             WHEN d.itemid = 220210 THEN 'respiratory_rate'
             WHEN d.itemid IN (220050,220179,224167,227243,225309) THEN 'sbp'
             WHEN d.itemid IN (220051,220180,224643,227242,225310) THEN 'dbp'
             WHEN d.itemid IN (223762,226329) THEN 'temperature_c'
             WHEN d.itemid = 223761 THEN 'temperature_f'
             WHEN d.itemid = 220277 THEN 'spo2'
             WHEN d.itemid = 223835 THEN 'fio2'
             WHEN d.itemid IN (220621,226537,228388,225664) THEN 'glucose'
             WHEN d.itemid IN (226730,226707) THEN 'height'
             WHEN d.itemid IN (226512,226531,224639) THEN 'weight'
             WHEN d.itemid IN (225792,225794,224385) THEN 'ventilation_marker'
             WHEN d.itemid IN (225802,225803,225809,225955,225805,225441) THEN 'rrt_marker'
             WHEN d.itemid IN (221906,221289,222315,221749,221662,221653,221986) THEN 'vaso_inotrope'
             WHEN d.itemid IN (223257,223258,223259,223260,223261,223262,229299,229619) THEN 'insulin'
             WHEN lower(d.label) ~ '(hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol)' THEN 'steroid_candidate'
           END AS variable_name,
           CASE d.linksto WHEN 'chartevents' THEN 'mimiciv_icu.chartevents'
                          WHEN 'inputevents' THEN 'mimiciv_icu.inputevents'
                          WHEN 'procedureevents' THEN 'mimiciv_icu.procedureevents'
                          ELSE 'mimiciv_icu.' || d.linksto END AS source_table,
           d.itemid,
           d.label,
           concat_ws(' | ', d.category, d.unitname, d.param_type) AS descriptor,
           (d.itemid IN (220045,220210,220050,220179,224167,227243,225309,
                         220051,220180,224643,227242,225310,223762,226329,
                         223761,220277,223835,220621,226537,228388,225664,
                         226730,226707,226512,226531,224639,225792,225794,
                         224385,225802,225803,225809,225955,225805,225441,
                         221906,221289,222315,221749,221662,221653,221986,
                         223257,223258,223259,223260,223261,223262,229299,229619)) AS exact_candidate,
           CASE WHEN d.itemid IN (225792,225794,224385,225802,225803,225809,225955,225805,225441,
                                  221906,221289,222315,221749,221662,221653,221986,
                                  223257,223258,223259,223260,223261,223262,229299,229619)
                THEN 'treatment/procedure label or explicit itemid'
                ELSE 'explicit vital/oxygen/glucose/anthropometry itemid' END AS match_keyword
    FROM mimiciv_icu.d_items d
    WHERE d.itemid IN (220045,220210,220050,220179,224167,227243,225309,
                       220051,220180,224643,227242,225310,223762,226329,
                       223761,220277,223835,220621,226537,228388,225664,
                       226730,226707,226512,226531,224639,225792,225794,
                       224385,225802,225803,225809,225955,225805,225441,
                       221906,221289,222315,221749,221662,221653,221986,
                       223257,223258,223259,223260,223261,223262,229299,229619)
       OR lower(d.label) ~ '(hydrocortisone|methylpred|dexamethasone|prednisone|prednisolone|solumedrol)'
)
INSERT INTO mimic_custom.itemid_review_v2
    (variable_name, source_table, itemid, label, fluid_or_category_or_unit,
     include_candidate, exclude_reason, match_keyword, needs_manual_review, review_status)
SELECT variable_name,
       source_table,
       itemid,
       label,
       descriptor,
       exact_candidate,
       CASE WHEN exact_candidate THEN NULL ELSE 'broad label match only; not used by default' END,
       match_keyword,
       NOT exact_candidate,
       CASE WHEN exact_candidate THEN 'auto_include' ELSE 'manual_review' END
FROM candidates
WHERE variable_name IS NOT NULL;

-- Hospital laboratory items. The blood-fluid restriction avoids urine/body-fluid
-- items that share names with the requested analytes.
WITH lab_map AS (
    SELECT d.itemid,
           d.label,
           concat_ws(' | ', d.fluid, d.category) AS descriptor,
           CASE
             WHEN d.itemid IN (50852,51631) THEN 'hba1c_pct'
             WHEN d.itemid IN (50809,50931,52569,52027) THEN 'glucose'
             WHEN d.itemid IN (50912,52546) THEN 'creatinine'
             WHEN d.itemid IN (51006,52647) THEN 'bun'
             WHEN d.itemid IN (50813,52442,53154) THEN 'lactate'
             WHEN d.itemid IN (50820) THEN 'ph'
             WHEN d.itemid IN (51300,51301,51755,51756) THEN 'wbc'
             WHEN d.itemid IN (51222,51640,50811) THEN 'hemoglobin'
             WHEN d.itemid IN (52170) THEN 'rbc'
             WHEN d.itemid IN (51221,51638,51639,52028,50810) THEN 'hematocrit'
             WHEN d.itemid IN (51250,51691) THEN 'mcv'
             WHEN d.itemid IN (50885,53089) THEN 'bilirubin'
             WHEN d.itemid IN (50861) THEN 'alt'
             WHEN d.itemid IN (50863,53086) THEN 'alp'
             WHEN d.itemid IN (50878) THEN 'ast'
             WHEN d.itemid IN (50862,53085) THEN 'albumin'
             WHEN d.itemid IN (50971,52610,50822,52452) THEN 'potassium'
             WHEN d.itemid IN (50983,52623,50824,52455) THEN 'sodium'
             WHEN d.itemid IN (50902,52535,50806,52434) THEN 'chloride'
             WHEN d.itemid IN (50893,52034,52035) THEN 'calcium'
             WHEN d.itemid IN (50868,52500) THEN 'anion_gap'
             WHEN d.itemid IN (50882,51739,50803,50804,52037) THEN 'bicarbonate'
             WHEN d.itemid IN (51274,52921) THEN 'pt'
             WHEN d.itemid IN (51275,52923) THEN 'ptt'
             WHEN d.itemid IN (51623,52116,51214) THEN 'fibrinogen'
             WHEN d.itemid IN (51265,53189) THEN 'platelets'
           END AS variable_name,
           d.itemid IN (50852,51631,50809,50931,52569,52027,50912,52546,51006,52647,50813,52442,53154,
                        50820,51300,51301,51755,51756,51222,51640,50811,52170,
                        51221,51638,51639,52028,50810,51250,51691,50885,53089,
                        50861,50863,53086,50878,50862,53085,50971,52610,50822,
                        52452,50983,52623,50824,52455,50902,52535,50806,52434,
                        50893,52034,52035,50868,52500,50882,51739,50803,50804,
                        52037,51274,52921,51275,52923,51623,52116,51214,51265,53189) AS exact_candidate
    FROM mimiciv_hosp.d_labitems d
    WHERE d.fluid = 'Blood'
      AND d.itemid IN (50852,51631,50809,50931,52569,52027,50912,52546,51006,52647,50813,52442,53154,
                       50820,51300,51301,51755,51756,51222,51640,50811,52170,
                       51221,51638,51639,52028,50810,51250,51691,50885,53089,
                       50861,50863,53086,50878,50862,53085,50971,52610,50822,
                       52452,50983,52623,50824,52455,50902,52535,50806,52434,
                       50893,52034,52035,50868,52500,50882,51739,50803,50804,
                       52037,51274,52921,51275,52923,51623,52116,51214,51265,53189)
)
INSERT INTO mimic_custom.itemid_review_v2
    (variable_name, source_table, itemid, label, fluid_or_category_or_unit,
     include_candidate, exclude_reason, match_keyword, needs_manual_review, review_status)
SELECT variable_name,
       'mimiciv_hosp.labevents',
       itemid,
       label,
       descriptor,
       exact_candidate,
       NULL,
       'explicit blood lab itemid',
       false,
       'auto_include'
FROM lab_map
WHERE variable_name IS NOT NULL;

-- Make missing dictionary mappings visible rather than silently dropping a
-- requested variable from downstream tables.
WITH required(variable_name, source_table) AS (
    VALUES
      ('heart_rate','mimiciv_icu.chartevents'),('respiratory_rate','mimiciv_icu.chartevents'),
      ('sbp','mimiciv_icu.chartevents'),('dbp','mimiciv_icu.chartevents'),
      ('temperature_c','mimiciv_icu.chartevents'),('temperature_f','mimiciv_icu.chartevents'),
      ('spo2','mimiciv_icu.chartevents'),('fio2','mimiciv_icu.chartevents'),
      ('glucose','mimiciv_icu.chartevents'),('glucose','mimiciv_hosp.labevents'),
      ('height','mimiciv_icu.chartevents'),
      ('weight','mimiciv_icu.chartevents'),('ventilation_marker','mimiciv_icu.procedureevents'),
      ('rrt_marker','mimiciv_icu.procedureevents'),('vaso_inotrope','mimiciv_icu.inputevents'),
      ('insulin','mimiciv_icu.inputevents'),('hba1c_pct','mimiciv_hosp.labevents'),
      ('creatinine','mimiciv_hosp.labevents'),('bun','mimiciv_hosp.labevents'),
      ('lactate','mimiciv_hosp.labevents'),('ph','mimiciv_hosp.labevents'),
      ('wbc','mimiciv_hosp.labevents'),('hemoglobin','mimiciv_hosp.labevents'),
      ('rbc','mimiciv_hosp.labevents'),('hematocrit','mimiciv_hosp.labevents'),
      ('mcv','mimiciv_hosp.labevents'),('bilirubin','mimiciv_hosp.labevents'),
      ('alt','mimiciv_hosp.labevents'),('alp','mimiciv_hosp.labevents'),
      ('ast','mimiciv_hosp.labevents'),('albumin','mimiciv_hosp.labevents'),
      ('potassium','mimiciv_hosp.labevents'),('sodium','mimiciv_hosp.labevents'),
      ('chloride','mimiciv_hosp.labevents'),('calcium','mimiciv_hosp.labevents'),
      ('anion_gap','mimiciv_hosp.labevents'),('bicarbonate','mimiciv_hosp.labevents'),
      ('pt','mimiciv_hosp.labevents'),('ptt','mimiciv_hosp.labevents'),
      ('fibrinogen','mimiciv_hosp.labevents'),('platelets','mimiciv_hosp.labevents')
), found AS (
    SELECT DISTINCT variable_name, source_table
    FROM mimic_custom.itemid_review_v2
    WHERE include_candidate
)
INSERT INTO mimic_custom.itemid_review_v2
    (variable_name, source_table, itemid, label, fluid_or_category_or_unit,
     include_candidate, exclude_reason, match_keyword, needs_manual_review, review_status)
SELECT r.variable_name,
       r.source_table,
       NULL,
       NULL,
       NULL,
       false,
       'no curated itemid found; variable unavailable until manual review',
       'required_variable_sentinel',
       true,
       'missing_mapping'
FROM required r
LEFT JOIN found f USING (variable_name, source_table)
WHERE f.variable_name IS NULL;

CREATE INDEX itemid_review_v2_variable_idx
    ON mimic_custom.itemid_review_v2 (variable_name, include_candidate);
CREATE INDEX itemid_review_v2_source_item_idx
    ON mimic_custom.itemid_review_v2 (source_table, itemid);

COMMIT;

\copy (SELECT * FROM mimic_custom.itemid_review_v2 ORDER BY source_table, variable_name, itemid NULLS LAST) TO 'outputs/qc/itemid_review_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
