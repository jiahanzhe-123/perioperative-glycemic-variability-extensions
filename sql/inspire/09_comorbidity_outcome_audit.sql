-- 09_comorbidity_outcome_audit.sql — 协变量(糖尿病/Charlson-wo-DM,ICD-10-CM)与结局结构审计
-- 说明:本步仅做数据结构审计(可用率、事件数、随访完整性),不运行任何结局模型。
BEGIN;

-- Deyo/Quan Charlson ICD-10(不含糖尿病及糖尿病并发症)
DROP TABLE IF EXISTS derived.comorbidity;
CREATE TABLE derived.comorbidity AS
WITH dx AS (
    SELECT DISTINCT subject_id, icd10_cm FROM raw.diagnosis WHERE icd10_cm IS NOT NULL
),
flags AS (
    SELECT subject_id,
      max((icd10_cm ~ '^(I21|I22|I252)')::int) AS mi,
      max((icd10_cm ~ '^(I43|I50|I110|I255|I42)')::int) AS chf,
      max((icd10_cm ~ '^(I70|I71|I731|I738|I739|I771|I790|I792|K551|K558|K559|Z958|Z959)')::int) AS pvd,
      max((icd10_cm ~ '^(G45|G46|H34|I60|I61|I62|I63|I64|I65|I66|I67|I68|I69)')::int) AS cvd,
      max((icd10_cm ~ '^(F00|F01|F02|F03|F051|G30|G311)')::int) AS dementia,
      max((icd10_cm ~ '^(I278|I279|J40|J41|J42|J43|J44|J45|J46|J47|J60|J61|J62|J63|J64|J65|J66|J67)')::int) AS copd,
      max((icd10_cm ~ '^(M05|M06|M315|M32|M33|M34|M351|M353|M360)')::int) AS rheumatic,
      max((icd10_cm ~ '^(K25|K26|K27|K28)')::int) AS pud,
      max((icd10_cm ~ '^(K70|K74|K760)')::int) AS mild_liver,
      max((icd10_cm ~ '^(E10|E11|E12|E13|E14)')::int) AS diabetes_any,
      max((icd10_cm ~ '^(E102|E103|E104|E105|E106|E107|E108|E112|E113|E114|E115|E116|E117|E118|E122|E123|E124|E125|E126|E127|E128|E132|E133|E134|E135|E136|E137|E138|E142|E143|E144|E145|E146|E147|E148)')::int) AS diabetes_complicated,
      max((icd10_cm ~ '^(G81|G82|G041|G114|G801|G802|G830|G831|G832|G833|G834|G839)')::int) AS hemiplegia,
      max((icd10_cm ~ '^(N18|N19|N25|Z49|Z940|Z992)')::int) AS renal,
      max((icd10_cm ~ '^(C0|C1|C2|C30|C31|C32|C33|C34|C37|C38|C39|C40|C41|C43|C45|C46|C47|C48|C49|C50|C51|C52|C53|C54|C55|C56|C57|C58|C60|C61|C62|C63|C64|C65|C66|C67|C68|C69|C70|C71|C72|C73|C74|C75|C76|C80|C81|C82|C83|C84|C85|C88|C90|C91|C92|C93|C94|C95|C96|C97)')::int) AS malignancy,
      max((icd10_cm ~ '^(K702|K703|K704|K705|K706|K707|K713|K714|K715|K717|K721|K729|K73|K74|K760|K762|K763|K764|K765|K766|K767|K768|K769|Z944)')::int) AS severe_liver,
      max((icd10_cm ~ '^(C77|C78|C79|C80)')::int) AS metastatic,
      max((icd10_cm ~ '^(B20|B21|B22|B23|B24)')::int) AS aids
    FROM dx GROUP BY subject_id
)
SELECT o.subject_id,
       (mi + chf + pvd + cvd + dementia + copd + rheumatic + pud + mild_liver
        + diabetes_any + hemiplegia + renal + malignancy + severe_liver*3 + metastatic*6 + aids*6
        + CASE WHEN o.age >= 80 THEN 4 WHEN o.age >= 70 THEN 3 WHEN o.age >= 60 THEN 2 WHEN o.age >= 50 THEN 1 ELSE 0 END) AS charlson_index,
       (mi + chf + pvd + cvd + dementia + copd + rheumatic + pud + mild_liver
        + hemiplegia + renal + malignancy + severe_liver*3 + metastatic*6 + aids*6
        + CASE WHEN o.age >= 80 THEN 4 WHEN o.age >= 70 THEN 3 WHEN o.age >= 60 THEN 2 WHEN o.age >= 50 THEN 1 ELSE 0 END) AS charlson_without_diabetes,
       diabetes_any AS diabetes,
       diabetes_complicated AS diabetes_complicated
FROM derived.cohort_first_operation o
LEFT JOIN flags USING (subject_id);

UPDATE derived.comorbidity SET
  charlson_index = COALESCE(charlson_index, 0),
  charlson_without_diabetes = COALESCE(charlson_without_diabetes, 0),
  diabetes = COALESCE(diabetes, 0),
  diabetes_complicated = COALESCE(diabetes_complicated, 0);

-- 结构审计:结局可用性与随访完整性(仅结构,不建模)
DROP TABLE IF EXISTS audit.outcome_structure;
CREATE TABLE audit.outcome_structure AS
SELECT count(*) AS n_cohort,
  count(*) FILTER (WHERE inhosp_death_time IS NOT NULL) AS n_inhosp_death_total,
  count(*) FILTER (WHERE allcause_death_time IS NOT NULL) AS n_allcause_death_total,
  count(*) FILTER (WHERE inhosp_death_time IS NOT NULL AND inhosp_death_time <= opend_time + 30*1440) AS n_inhosp_death_le_30d,
  count(*) FILTER (WHERE inhosp_death_time IS NOT NULL AND inhosp_death_time <= opend_time + 365*1440) AS n_inhosp_death_le_365d,
  count(*) FILTER (WHERE allcause_death_time IS NOT NULL AND allcause_death_time <= opend_time + 30*1440) AS n_allcause_death_le_30d,
  count(*) FILTER (WHERE allcause_death_time IS NOT NULL AND allcause_death_time <= opend_time + 365*1440) AS n_allcause_death_le_365d,
  count(*) FILTER (WHERE inhosp_death_time IS NOT NULL AND inhosp_death_time < opend_time + 1440) AS n_death_before_landmark_inhosp,
  count(*) FILTER (WHERE discharge_time < opend_time + 1440 AND inhosp_death_time IS NULL) AS n_discharge_before_landmark_alive,
  round(percentile_cont(0.5) WITHIN GROUP (ORDER BY discharge_time - opend_time)::numeric,1) AS median_opend_to_discharge_min,
  count(*) FILTER (WHERE discharge_time - opend_time >= 30*1440) AS n_observed_ge_30d,
  count(*) FILTER (WHERE discharge_time - opend_time >= 365*1440) AS n_observed_ge_365d,
  count(*) FILTER (WHERE inhosp_death_time IS NOT NULL AND allcause_death_time IS NOT NULL AND allcause_death_time < inhosp_death_time) AS n_allcause_before_inhosp,
  count(*) FILTER (WHERE allcause_death_time IS NULL AND inhosp_death_time IS NOT NULL) AS n_inhosp_without_allcause
FROM derived.cohort_first_operation;

-- 暴露可用率(结构)
DROP TABLE IF EXISTS audit.exposure_structure;
CREATE TABLE audit.exposure_structure AS
SELECT count(*) AS n_cohort,
  count(*) FILTER (WHERE b.n_glucose_0_24h >= 2) AS n_gv_ge2,
  count(*) FILTER (WHERE b.n_glucose_0_24h >= 2 AND COALESCE(b.inhosp_death_time, b.discharge_time) >= b.opend_time + 1440) AS n_gv_ge2_landmark,
  count(*) FILTER (WHERE b.n_glucose_0_24h >= 3) AS n_gv_ge3,
  count(*) FILTER (WHERE b.n_glucose_0_24h >= 3 AND b.n_distinct_values >= 3) AS n_gv_ge3_ge3distinct,
  count(*) FILTER (WHERE b.shr IS NOT NULL) AS n_shr,
  round(avg(b.age)::numeric,1) AS mean_age,
  count(*) FILTER (WHERE b.weight IS NULL OR b.height IS NULL) AS n_bmi_missing,
  count(*) FILTER (WHERE b.asa IS NULL) AS n_asa_missing,
  count(*) FILTER (WHERE coalesce(cm.diabetes,0)=1) AS n_diabetes,
  count(*) FILTER (WHERE b.cpbon_time IS NOT NULL) AS n_cpb
FROM analysis.inspire_base b LEFT JOIN derived.comorbidity cm USING (subject_id);

COMMIT;
