-- 06_extract_comorbidities_bmi.sql
-- ICD phenotype rules are intentionally stored as a review table. The
-- Charlson result is labelled provisional because this local implementation
-- is a transparent ICD-pattern mapping, not an imported derived concept.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.comorbidities_bmi_v2;
DROP TABLE IF EXISTS mimic_custom.comorbidity_codebook_v2;

CREATE TABLE mimic_custom.comorbidity_codebook_v2 (
    icd_version smallint NOT NULL,
    icd_code text NOT NULL,
    diagnosis_title text,
    component text NOT NULL,
    match_keyword text NOT NULL,
    include_candidate boolean NOT NULL,
    needs_manual_review boolean NOT NULL,
    confidence text NOT NULL,
    rule_version text NOT NULL DEFAULT 'v2_icd_pattern_review_2026-07-18'
);

WITH rules(component, icd_version, code_regex, match_keyword, confidence) AS (
    VALUES
      ('myocardial_infarct',9,'^(410|412)','ICD9 410/412', 'high'),
      ('myocardial_infarct',10,'^(I21|I22|I25\.2)','ICD10 I21/I22/I25.2', 'high'),
      ('congestive_heart_failure',9,'^(428)','ICD9 428', 'high'),
      ('congestive_heart_failure',10,'^(I50|I11\.0|I13\.(0|2))','ICD10 I50/I11.0/I13.0-I13.2', 'high'),
      ('peripheral_vascular_disease',9,'^(440|441|443\.1|443\.9|785\.4|V43\.4)','ICD9 vascular disease pattern', 'moderate'),
      ('peripheral_vascular_disease',10,'^(I70|I71|I72|I73\.1|I73\.9|I77\.6|Z95\.82)','ICD10 vascular disease pattern', 'moderate'),
      ('cerebrovascular_disease',9,'^(430|431|432|433|434|435|436|437|438)','ICD9 cerebrovascular pattern', 'high'),
      ('cerebrovascular_disease',10,'^(G45|G46|I6[0-9])','ICD10 cerebrovascular pattern', 'high'),
      ('dementia',9,'^(290|294\.1)','ICD9 dementia pattern', 'moderate'),
      ('dementia',10,'^(F00|F01|F02|F03|G30)','ICD10 dementia pattern', 'moderate'),
      ('chronic_pulmonary_disease',9,'^(490|491|492|493|494|495|496)','ICD9 chronic pulmonary pattern', 'high'),
      ('chronic_pulmonary_disease',10,'^(J4|J6[0-9]|J70\.3)','ICD10 chronic pulmonary pattern', 'high'),
      ('connective_tissue_disease',9,'^(710|714|725)','ICD9 connective tissue pattern', 'moderate'),
      ('connective_tissue_disease',10,'^(M05|M06|M32|M33|M34|M35)','ICD10 connective tissue pattern', 'moderate'),
      ('ulcer_disease',9,'^(531|532|533|534)','ICD9 peptic ulcer pattern', 'high'),
      ('ulcer_disease',10,'^(K25|K26|K27|K28)','ICD10 peptic ulcer pattern', 'high'),
      ('mild_liver_disease',9,'^(571|070)','ICD9 liver disease pattern', 'moderate'),
      ('mild_liver_disease',10,'^(B18|K70|K73|K74)','ICD10 liver disease pattern', 'moderate'),
      ('diabetes_without_cc',9,'^(250\.[0-3])','ICD9 diabetes without complication', 'high'),
      ('diabetes_without_cc',10,'^(E1[0-4]\.(0|1|9))','ICD10 diabetes without complication', 'high'),
      ('diabetes_with_cc',9,'^(250\.[4-9])','ICD9 diabetes with complication', 'high'),
      ('diabetes_with_cc',10,'^(E1[0-4]\.(2|3|4|5|6|7|8))','ICD10 diabetes with complication', 'high'),
      ('renal_disease',9,'^(585|586|V42\.0|V45\.1|V56)','ICD9 renal disease pattern', 'high'),
      ('renal_disease',10,'^(N18|N19|Z49|Z94\.0|Z99\.2)','ICD10 renal disease pattern', 'high'),
      ('esrd',9,'^(585\.[5-6]|V45\.1|V56)','ICD9 ESRD pattern', 'high'),
      ('esrd',10,'^(N18\.[5-6]|N18\.0|Z99\.2)','ICD10 ESRD pattern', 'high'),
      ('malignant_cancer',9,'^(140|141|142|143|144|145|146|147|148|149|150|151|152|153|154|155|156|157|158|159|160|161|162|163|164|165|170|171|172|174|175|176|179|180|181|182|183|184|185|186|187|188|189|190|191|192|193|194|195|196|197|198|199|200|201|202|203|204|205|206|207|208)','ICD9 malignancy pattern', 'moderate'),
      ('malignant_cancer',10,'^C[0-9]','ICD10 malignancy pattern', 'moderate'),
      ('leukemia',9,'^(204|205|206|207)','ICD9 leukemia pattern', 'high'),
      ('leukemia',10,'^(C91|C92|C93|C94|C95)','ICD10 leukemia pattern', 'high'),
      ('lymphoma',9,'^(200|201|202|203)','ICD9 lymphoma pattern', 'high'),
      ('lymphoma',10,'^(C81|C82|C83|C84|C85|C86|C88)','ICD10 lymphoma pattern', 'high'),
      ('severe_liver_disease',9,'^(456\.[0-2]|572\.[2-8]|571\.[2-6])','ICD9 severe liver pattern', 'moderate'),
      ('severe_liver_disease',10,'^(K72|K76\.6|K76\.7|K74\.[4-6]|I85)','ICD10 severe liver pattern', 'moderate'),
      ('metastatic_cancer',9,'^(196|197|198|199)','ICD9 metastatic malignancy pattern', 'moderate'),
      ('metastatic_cancer',10,'^(C77|C78|C79|C80)','ICD10 metastatic malignancy pattern', 'moderate'),
      ('aids',9,'^(042|043|044)','ICD9 HIV/AIDS pattern', 'high'),
      ('aids',10,'^(B20|B21|B22|B23|B24)','ICD10 HIV/AIDS pattern', 'high'),
      ('hypertension',9,'^(401|402|403|404|405)','ICD9 hypertension pattern', 'high'),
      ('hypertension',10,'^(I10|I11|I12|I13|I15|I16)','ICD10 hypertension pattern', 'high'),
      ('obesity',9,'^(278\.0)','ICD9 obesity', 'moderate'),
      ('obesity',10,'^(E66)','ICD10 obesity', 'moderate'),
      ('dyslipidemia',9,'^(272)','ICD9 dyslipidemia', 'moderate'),
      ('dyslipidemia',10,'^(E78)','ICD10 dyslipidemia', 'moderate'),
      ('smoking',9,'^(305\.1|V15\.82)','ICD9 tobacco history', 'low'),
      ('smoking',10,'^(F17|Z72\.0|Z87\.891)','ICD10 tobacco history', 'low'),
      ('acute_kidney_injury',9,'^(584)','ICD9 AKI', 'high'),
      ('acute_kidney_injury',10,'^(N17)','ICD10 AKI', 'high'),
      ('stroke',9,'^(430|431|432|433|434|435|436|437|438)','ICD9 stroke/cerebrovascular', 'high'),
      ('stroke',10,'^(I6[0-9]|G45|G46)','ICD10 stroke/cerebrovascular', 'high'),
      ('sepsis',9,'^(038|995\.9|785\.52)','ICD9 sepsis', 'moderate'),
      ('sepsis',10,'^(A40|A41|R65\.2|R57\.2)','ICD10 sepsis/shock', 'moderate'),
      ('infection',9,'^(001|002|003|004|005|006|007|008|009|010|011|012|013|014|015|016|017|018|020|021|022|023|024|025|026|027|030|031|032|033|034|035|036|037|038|039|040|041|042|043|044|045|046|047|048|049|050|051|052|053|054|055|056|057|058|059|060|061|062|063|064|065|066|070|071|072|073|074|075|076|077|078|079|080|081|082|083|084|085|086|087|088|089|090|091|092|093|094|095|096|097|098|099)','ICD9 infection chapter', 'low'),
      ('infection',10,'^[A-B]','ICD10 infection chapter', 'low'),
      ('atrial_fibrillation',9,'^(427\.31|427\.32)','ICD9 AF', 'high'),
      ('atrial_fibrillation',10,'^(I48)','ICD10 AF', 'high'),
      ('valvular_disease',9,'^(394|395|396|397)','ICD9 valvular disease', 'moderate'),
      ('valvular_disease',10,'^(I05|I06|I07|I08|I34|I35|I36|I37|I38|Q23)','ICD10 valvular disease', 'moderate'),
      ('ischemic_heart_disease',9,'^(410|411|412|413|414)','ICD9 ischemic heart disease', 'high'),
      ('ischemic_heart_disease',10,'^(I20|I21|I22|I23|I24|I25)','ICD10 ischemic heart disease', 'high'),
      ('cardiogenic_shock',9,'^(785\.51)','ICD9 cardiogenic shock', 'high'),
      ('cardiogenic_shock',10,'^(R57\.0)','ICD10 cardiogenic shock', 'high'),
      ('cardiac_arrest',9,'^(427\.5)','ICD9 cardiac arrest', 'high'),
      ('cardiac_arrest',10,'^(I46)','ICD10 cardiac arrest', 'high'),
      ('pulmonary_hypertension',9,'^(416\.[0-8])','ICD9 pulmonary hypertension', 'moderate'),
      ('pulmonary_hypertension',10,'^(I27)','ICD10 pulmonary hypertension', 'moderate'),
      ('endocarditis',9,'^(421)','ICD9 endocarditis', 'moderate'),
      ('endocarditis',10,'^(I33|I38)','ICD10 endocarditis', 'moderate')
), dx AS (
    SELECT d.icd_version,
           btrim(d.icd_code::text) AS icd_code,
           btrim(d.long_title) AS diagnosis_title
    FROM mimiciv_hosp.d_icd_diagnoses d
)
INSERT INTO mimic_custom.comorbidity_codebook_v2
    (icd_version, icd_code, diagnosis_title, component, match_keyword,
     include_candidate, needs_manual_review, confidence)
SELECT dx.icd_version,
       dx.icd_code,
       dx.diagnosis_title,
       r.component,
       r.match_keyword,
       true,
       (r.confidence = 'low'),
       r.confidence
FROM dx
JOIN rules r
  ON r.icd_version = dx.icd_version
 AND dx.icd_code ~ r.code_regex;

CREATE INDEX comorbidity_codebook_v2_code_idx
    ON mimic_custom.comorbidity_codebook_v2 (icd_version, icd_code);
CREATE INDEX comorbidity_codebook_v2_component_idx
    ON mimic_custom.comorbidity_codebook_v2 (component);

CREATE TABLE mimic_custom.comorbidities_bmi_v2 AS
WITH current_dx AS (
    SELECT c.subject_id,
           c.hadm_id,
           cb.component,
           true AS flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_hosp.diagnoses_icd d
      ON d.subject_id = c.subject_id AND d.hadm_id = c.hadm_id
    JOIN mimic_custom.comorbidity_codebook_v2 cb
      ON cb.icd_version = d.icd_version
     AND cb.icd_code = btrim(d.icd_code::text)
    GROUP BY c.subject_id, c.hadm_id, cb.component
), dx_flags AS (
    SELECT c.subject_id,
           c.hadm_id,
           bool_or(cd.component = 'myocardial_infarct') AS myocardial_infarct,
           bool_or(cd.component = 'congestive_heart_failure') AS congestive_heart_failure,
           bool_or(cd.component = 'peripheral_vascular_disease') AS peripheral_vascular_disease,
           bool_or(cd.component = 'cerebrovascular_disease') AS cerebrovascular_disease,
           bool_or(cd.component = 'dementia') AS dementia,
           bool_or(cd.component = 'chronic_pulmonary_disease') AS chronic_pulmonary_disease,
           bool_or(cd.component = 'connective_tissue_disease') AS connective_tissue_disease,
           bool_or(cd.component = 'ulcer_disease') AS ulcer_disease,
           bool_or(cd.component = 'mild_liver_disease') AS mild_liver_raw,
           bool_or(cd.component = 'severe_liver_disease') AS severe_liver_disease,
           bool_or(cd.component = 'diabetes_without_cc') AS diabetes_without_cc,
           bool_or(cd.component = 'diabetes_with_cc') AS diabetes_with_cc,
           bool_or(cd.component = 'renal_disease') AS renal_disease,
           bool_or(cd.component = 'esrd') AS esrd,
           bool_or(cd.component = 'malignant_cancer') AS malignant_cancer,
           bool_or(cd.component = 'leukemia') AS leukemia,
           bool_or(cd.component = 'lymphoma') AS lymphoma,
           bool_or(cd.component = 'metastatic_cancer') AS metastatic_cancer,
           bool_or(cd.component = 'aids') AS aids,
           bool_or(cd.component = 'hypertension') AS hypertension,
           bool_or(cd.component = 'obesity') AS obesity_flag,
           bool_or(cd.component = 'dyslipidemia') AS dyslipidemia_flag,
           bool_or(cd.component = 'smoking') AS smoking_flag,
           bool_or(cd.component = 'acute_kidney_injury') AS aki_flag,
           bool_or(cd.component = 'stroke') AS stroke_flag,
           bool_or(cd.component = 'sepsis') AS sepsis_flag,
           bool_or(cd.component = 'infection') AS infection_flag,
           bool_or(cd.component = 'atrial_fibrillation') AS atrial_fibrillation_flag,
           bool_or(cd.component = 'valvular_disease') AS valvular_disease_flag,
           bool_or(cd.component = 'ischemic_heart_disease') AS ischemic_heart_disease_flag,
           bool_or(cd.component = 'cardiogenic_shock') AS cardiogenic_shock_flag,
           bool_or(cd.component = 'cardiac_arrest') AS cardiac_arrest_flag,
           bool_or(cd.component = 'pulmonary_hypertension') AS pulmonary_hypertension_flag,
           bool_or(cd.component = 'endocarditis') AS endocarditis_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN current_dx cd ON cd.subject_id = c.subject_id AND cd.hadm_id = c.hadm_id
    GROUP BY c.subject_id, c.hadm_id
), prior_revasc AS (
    SELECT c.subject_id,
           c.hadm_id,
           EXISTS (
             SELECT 1
             FROM mimiciv_hosp.procedures_icd p
             JOIN mimic_custom.cardiac_surgery_codebook_v2 cb
               ON cb.icd_version = p.icd_version
              AND cb.icd_code = btrim(p.icd_code::text)
              AND cb.procedure_group IN ('cabg','pci')
              AND cb.include_primary
             WHERE p.subject_id = c.subject_id
               AND p.chartdate < c.surgery_time::date
           ) AS prior_pci_or_cabg_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
), post_creat AS (
    SELECT c.stay_id,
           max(l.value_valid) AS creat_postop_48h_max
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimic_custom.labs_long_v2 l ON l.stay_id = c.stay_id
    WHERE l.variable_name = 'creatinine'
      AND l.value_valid IS NOT NULL
      AND l.charttime >= c.surgery_time
      AND l.charttime < c.surgery_time + interval '48 hours'
    GROUP BY c.stay_id
),
-- OMR names are stable in the imported MIMIC-IV table; values are parsed only
-- after stripping non-numeric characters and remain source-labelled.
omr_values AS (
    SELECT c.stay_id,
           max(CASE WHEN lower(o.result_name) IN ('bmi','bmi (kg/m2)')
                    THEN NULLIF(regexp_replace(o.result_value, '[^0-9.\-]', '', 'g'),'')::double precision END) AS omr_bmi,
           max(CASE WHEN lower(o.result_name) = 'height (inches)'
                    THEN NULLIF(regexp_replace(o.result_value, '[^0-9.\-]', '', 'g'),'')::double precision * 2.54 END) AS omr_height_cm,
           max(CASE WHEN lower(o.result_name) = 'weight (lbs)'
                    THEN NULLIF(regexp_replace(o.result_value, '[^0-9.\-]', '', 'g'),'')::double precision * 0.45359237 END) AS omr_weight_kg,
           max(o.chartdate) AS omr_chartdate
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_hosp.omr o
      ON o.subject_id = c.subject_id
     AND o.chartdate <= c.surgery_time::date
     AND o.chartdate >= (c.admittime_ts::date - 365)
    WHERE lower(o.result_name) IN ('bmi','bmi (kg/m2)','height (inches)','weight (lbs)')
    GROUP BY c.stay_id
), chart_hw AS (
    SELECT c.stay_id,
           max(CASE WHEN ce.itemid IN (226730,226707)
                    THEN CASE WHEN ce.itemid = 226730 THEN ce.valuenum ELSE ce.valuenum * 2.54 END END) AS chart_height_cm,
           max(CASE WHEN ce.itemid IN (226512,226531,224639)
                    THEN CASE WHEN ce.itemid = 226531 THEN ce.valuenum * 0.45359237 ELSE ce.valuenum END END) AS chart_weight_kg
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    JOIN mimiciv_icu.chartevents ce
      ON ce.stay_id = c.stay_id
     AND ce.charttime >= c.admittime_ts
     AND ce.charttime < c.surgery_time + interval '24 hours'
     AND ce.itemid IN (226730,226707,226512,226531,224639)
     AND ce.valuenum IS NOT NULL
    GROUP BY c.stay_id
), current_admission_mi AS (
    SELECT c.stay_id,
           bool_or(cb.component = 'myocardial_infarct') AS acute_mi_current_admission_flag
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN current_dx cb ON cb.subject_id = c.subject_id AND cb.hadm_id = c.hadm_id
    GROUP BY c.stay_id
)
SELECT c.subject_id,
       c.hadm_id,
       c.stay_id,
       COALESCE(f.myocardial_infarct,false) AS myocardial_infarct,
       COALESCE(f.congestive_heart_failure,false) AS congestive_heart_failure,
       COALESCE(f.peripheral_vascular_disease,false) AS peripheral_vascular_disease,
       COALESCE(f.cerebrovascular_disease,false) AS cerebrovascular_disease,
       COALESCE(f.diabetes_without_cc,false) AS diabetes_without_cc,
       COALESCE(f.diabetes_with_cc,false) AS diabetes_with_cc,
       (COALESCE(f.diabetes_without_cc,false) OR COALESCE(f.diabetes_with_cc,false)) AS diabetes,
       COALESCE(f.renal_disease,false) AS renal_disease,
       COALESCE(f.malignant_cancer,false) AS malignant_cancer,
       COALESCE(f.hypertension,false) AS hypertension,
       COALESCE(f.chronic_pulmonary_disease,false) AS chronic_pulmonary_disease,
       COALESCE(f.obesity_flag,false) AS obesity_flag,
       COALESCE(f.dyslipidemia_flag,false) AS dyslipidemia_flag,
       COALESCE(f.smoking_flag,false) AS smoking_flag,
       COALESCE(f.mild_liver_raw,false) AND NOT COALESCE(f.severe_liver_disease,false) AS liver_disease_flag,
       COALESCE(f.severe_liver_disease,false) AS severe_liver_disease_flag,
       COALESCE(f.esrd,false) AS esrd_flag,
       COALESCE(f.aki_flag,false) AS aki_flag,
       COALESCE(f.stroke_flag,false) AS stroke_flag,
       COALESCE(f.sepsis_flag,false) AS sepsis_flag,
       COALESCE(f.infection_flag,false) AS infection_flag,
       COALESCE(f.atrial_fibrillation_flag,false) AS atrial_fibrillation_flag,
       COALESCE(f.congestive_heart_failure,false) AS heart_failure_flag,
       COALESCE(f.valvular_disease_flag,false) AS valvular_disease_flag,
       COALESCE(f.ischemic_heart_disease_flag,false) AS ischemic_heart_disease_flag,
       COALESCE(f.cardiogenic_shock_flag,false) AS cardiogenic_shock_flag,
       COALESCE(f.cardiac_arrest_flag,false) AS cardiac_arrest_flag,
       COALESCE(f.pulmonary_hypertension_flag,false) AS pulmonary_hypertension_flag,
       COALESCE(f.endocarditis_flag,false) AS endocarditis_flag,
       COALESCE(pr.prior_pci_or_cabg_flag,false) AS prior_pci_or_cabg_flag,
       COALESCE(mi.acute_mi_current_admission_flag,false) AS acute_mi_current_admission_flag,
       CASE WHEN COALESCE(f.diabetes_with_cc,false) THEN 2
            WHEN COALESCE(f.diabetes_without_cc,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.myocardial_infarct,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.congestive_heart_failure,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.peripheral_vascular_disease,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.cerebrovascular_disease,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.dementia,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.chronic_pulmonary_disease,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.connective_tissue_disease,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.ulcer_disease,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.mild_liver_raw,false) AND NOT COALESCE(f.severe_liver_disease,false) THEN 1 ELSE 0 END
       + CASE WHEN COALESCE(f.renal_disease,false) THEN 2 ELSE 0 END
       + CASE WHEN COALESCE(f.malignant_cancer,false) THEN 2 ELSE 0 END
       + CASE WHEN COALESCE(f.leukemia,false) THEN 2 ELSE 0 END
       + CASE WHEN COALESCE(f.lymphoma,false) THEN 2 ELSE 0 END
       + CASE WHEN COALESCE(f.severe_liver_disease,false) THEN 3 ELSE 0 END
       + CASE WHEN COALESCE(f.metastatic_cancer,false) THEN 6 ELSE 0 END
       + CASE WHEN COALESCE(f.aids,false) THEN 6 ELSE 0 END AS charlson_comorbidity_index,
       true AS charlson_provisional,
       COALESCE(pc.creat_postop_48h_max >= 1.5 * NULLIF(ap.creat_adm_first,0), false) AS aki_postop_48h_flag,
       COALESCE(om.omr_bmi, CASE WHEN om.omr_height_cm > 0 AND om.omr_weight_kg > 0 THEN om.omr_weight_kg / power(om.omr_height_cm/100.0,2) END,
                CASE WHEN ch.chart_height_cm > 0 AND ch.chart_weight_kg > 0 THEN ch.chart_weight_kg / power(ch.chart_height_cm/100.0,2) END) AS bmi,
       COALESCE(om.omr_height_cm, ch.chart_height_cm) AS height_cm,
       COALESCE(om.omr_weight_kg, ch.chart_weight_kg) AS weight_kg,
       CASE WHEN om.omr_bmi IS NOT NULL THEN 'mimiciv_hosp.omr_bmi'
            WHEN om.omr_height_cm IS NOT NULL AND om.omr_weight_kg IS NOT NULL THEN 'mimiciv_hosp.omr_height_weight'
            WHEN ch.chart_height_cm IS NOT NULL AND ch.chart_weight_kg IS NOT NULL THEN 'mimiciv_icu.chartevents_height_weight'
            ELSE NULL END AS bmi_source,
       CASE WHEN COALESCE(om.omr_bmi, CASE WHEN om.omr_height_cm > 0 AND om.omr_weight_kg > 0 THEN om.omr_weight_kg / power(om.omr_height_cm/100.0,2) END,
                          CASE WHEN ch.chart_height_cm > 0 AND ch.chart_weight_kg > 0 THEN ch.chart_weight_kg / power(ch.chart_height_cm/100.0,2) END) IS NULL
            THEN 'no valid BMI or paired height/weight within admission-year window' ELSE NULL END AS bmi_missing_reason,
       CASE WHEN COALESCE(f.smoking_flag,false) THEN 'low_confidence_diagnosis_only' ELSE NULL END AS smoking_confidence,
       'mimiciv_hosp.diagnoses_icd + reviewed ICD pattern codebook; BMI from OMR then ICU chart fallback'::text AS comorbidity_source
FROM mimic_custom.cardiac_surgery_cohort_v2 c
LEFT JOIN dx_flags f ON f.subject_id = c.subject_id AND f.hadm_id = c.hadm_id
LEFT JOIN prior_revasc pr ON pr.subject_id = c.subject_id AND pr.hadm_id = c.hadm_id
LEFT JOIN current_admission_mi mi ON mi.stay_id = c.stay_id
LEFT JOIN post_creat pc ON pc.stay_id = c.stay_id
LEFT JOIN mimic_custom.labs_adm_first_v2 ap ON ap.stay_id = c.stay_id
LEFT JOIN omr_values om ON om.stay_id = c.stay_id
LEFT JOIN chart_hw ch ON ch.stay_id = c.stay_id;

CREATE INDEX comorbidities_bmi_v2_stay_idx
    ON mimic_custom.comorbidities_bmi_v2 (stay_id);

COMMIT;

\copy (SELECT * FROM mimic_custom.comorbidity_codebook_v2 ORDER BY icd_version, component, icd_code) TO 'outputs/qc/comorbidity_codebook_review.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
