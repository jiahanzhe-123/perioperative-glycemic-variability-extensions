-- 03_schema_audit.sql — 导入后结构审计:重复行、ID 缺失、候选主键、时间/数值可解析性
-- 结果写 audit schema,供 feasibility_report 引用。
BEGIN;

DROP TABLE IF EXISTS audit.table_counts;
CREATE TABLE audit.table_counts AS
SELECT 'operations' AS table_name, count(*) AS db_rows,
       (SELECT csv_rows FROM meta.load_log WHERE table_name='operations') AS csv_rows,
       count(*) = (SELECT csv_rows FROM meta.load_log WHERE table_name='operations') AS match_ok, now() AS audited_at FROM raw.operations
UNION ALL SELECT 'diagnosis', count(*), (SELECT csv_rows FROM meta.load_log WHERE table_name='diagnosis'),
       count(*) = (SELECT csv_rows FROM meta.load_log WHERE table_name='diagnosis'), now() FROM raw.diagnosis
UNION ALL SELECT 'labs', count(*), (SELECT csv_rows FROM meta.load_log WHERE table_name='labs'),
       count(*) = (SELECT csv_rows FROM meta.load_log WHERE table_name='labs'), now() FROM raw.labs
UNION ALL SELECT 'medications', count(*), (SELECT csv_rows FROM meta.load_log WHERE table_name='medications'),
       count(*) = (SELECT csv_rows FROM meta.load_log WHERE table_name='medications'), now() FROM raw.medications
UNION ALL SELECT 'vitals', count(*), (SELECT csv_rows FROM meta.load_log WHERE table_name='vitals'),
       count(*) = (SELECT csv_rows FROM meta.load_log WHERE table_name='vitals'), now() FROM raw.vitals
UNION ALL SELECT 'ward_vitals', count(*), (SELECT csv_rows FROM meta.load_log WHERE table_name='ward_vitals'),
       count(*) = (SELECT csv_rows FROM meta.load_log WHERE table_name='ward_vitals'), now() FROM raw.ward_vitals;

DROP TABLE IF EXISTS audit.pk_candidates;
CREATE TABLE audit.pk_candidates AS
SELECT 'operations' AS table_name, 'op_id' AS candidate,
       count(*) AS n, count(DISTINCT op_id) AS n_distinct,
       count(*) = count(DISTINCT op_id) AS is_unique,
       count(*) FILTER (WHERE op_id IS NULL) AS n_null FROM raw.operations GROUP BY 1,2
UNION ALL
SELECT 'operations', 'subject_id', count(*), count(DISTINCT subject_id),
       count(*) = count(DISTINCT subject_id), count(*) FILTER (WHERE subject_id IS NULL) FROM raw.operations GROUP BY 1,2
UNION ALL
SELECT 'operations', 'hadm_id', count(*), count(DISTINCT hadm_id),
       count(*) = count(DISTINCT hadm_id), count(*) FILTER (WHERE hadm_id IS NULL) FROM raw.operations GROUP BY 1,2
UNION ALL
SELECT 'diagnosis', 'subject_id+chart_time+icd10_cm', count(*),
       count(DISTINCT (subject_id, chart_time, icd10_cm)),
       count(*) = count(DISTINCT (subject_id, chart_time, icd10_cm)),
       count(*) FILTER (WHERE subject_id IS NULL) FROM raw.diagnosis GROUP BY 1,2
UNION ALL
SELECT 'labs', 'subject_id+chart_time+item_name', count(*),
       count(DISTINCT (subject_id, chart_time, item_name)),
       count(*) = count(DISTINCT (subject_id, chart_time, item_name)),
       count(*) FILTER (WHERE subject_id IS NULL) FROM raw.labs GROUP BY 1,2
UNION ALL
SELECT 'vitals', 'op_id+chart_time+item_name', count(*),
       count(DISTINCT (op_id, chart_time, item_name)),
       count(*) = count(DISTINCT (op_id, chart_time, item_name)),
       count(*) FILTER (WHERE op_id IS NULL) FROM raw.vitals GROUP BY 1,2
UNION ALL
SELECT 'ward_vitals', 'subject_id+chart_time+item_name', count(*),
       count(DISTINCT (subject_id, chart_time, item_name)),
       count(*) = count(DISTINCT (subject_id, chart_time, item_name)),
       count(*) FILTER (WHERE subject_id IS NULL) FROM raw.ward_vitals GROUP BY 1,2;

DROP TABLE IF EXISTS audit.exact_duplicates;
CREATE TABLE audit.exact_duplicates AS
SELECT 'operations' AS table_name, count(*) - count(DISTINCT (op_id, subject_id, hadm_id, opdate, age, sex, icd10_pcs, orin_time, opend_time)) AS n_exact_dup FROM raw.operations
UNION ALL SELECT 'diagnosis', count(*) - count(DISTINCT (subject_id, chart_time, icd10_cm)) FROM raw.diagnosis
UNION ALL SELECT 'labs', count(*) - count(DISTINCT (subject_id, chart_time, item_name, value)) FROM raw.labs
UNION ALL SELECT 'medications', count(*) - count(DISTINCT (subject_id, chart_time, drug_name, route)) FROM raw.medications
UNION ALL SELECT 'vitals', count(*) - count(DISTINCT (op_id, subject_id, chart_time, item_name, value)) FROM raw.vitals
UNION ALL SELECT 'ward_vitals', count(*) - count(DISTINCT (subject_id, chart_time, item_name, value)) FROM raw.ward_vitals;

DROP TABLE IF EXISTS audit.missingness_operations;
CREATE TABLE audit.missingness_operations AS
SELECT
  count(*) AS n,
  round(100*avg((age IS NULL)::int),2) AS age_pct_null,
  round(100*avg((sex IS NULL)::int),2) AS sex_pct_null,
  round(100*avg((weight IS NULL)::int),2) AS weight_pct_null,
  round(100*avg((height IS NULL)::int),2) AS height_pct_null,
  round(100*avg((asa IS NULL)::int),2) AS asa_pct_null,
  round(100*avg((emop IS NULL)::int),2) AS emop_pct_null,
  round(100*avg((department IS NULL)::int),2) AS department_pct_null,
  round(100*avg((antype IS NULL)::int),2) AS antype_pct_null,
  round(100*avg((icd10_pcs IS NULL)::int),2) AS icd10_pcs_pct_null,
  round(100*avg((orin_time IS NULL)::int),2) AS orin_pct_null,
  round(100*avg((orout_time IS NULL)::int),2) AS orout_pct_null,
  round(100*avg((opstart_time IS NULL)::int),2) AS opstart_pct_null,
  round(100*avg((opend_time IS NULL)::int),2) AS opend_pct_null,
  round(100*avg((admission_time IS NULL)::int),2) AS admission_pct_null,
  round(100*avg((discharge_time IS NULL)::int),2) AS discharge_pct_null,
  round(100*avg((anstart_time IS NULL)::int),2) AS anstart_pct_null,
  round(100*avg((anend_time IS NULL)::int),2) AS anend_pct_null,
  round(100*avg((cpbon_time IS NULL)::int),2) AS cpbon_pct_null,
  round(100*avg((cpboff_time IS NULL)::int),2) AS cpboff_pct_null,
  round(100*avg((icuin_time IS NULL)::int),2) AS icuin_pct_null,
  round(100*avg((icuout_time IS NULL)::int),2) AS icuout_pct_null,
  round(100*avg((inhosp_death_time IS NULL)::int),2) AS inhosp_death_pct_null,
  round(100*avg((allcause_death_time IS NULL)::int),2) AS allcause_death_pct_null
FROM raw.operations;

COMMIT;
