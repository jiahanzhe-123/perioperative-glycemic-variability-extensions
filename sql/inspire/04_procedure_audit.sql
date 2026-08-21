-- 04_procedure_audit.sql — 手术代码频数、科室分布、时间锚点审计(不接触结局)
BEGIN;

DROP TABLE IF EXISTS audit.cardiac_code_frequency;
CREATE TABLE audit.cardiac_code_frequency AS
SELECT icd10_pcs,
       length(icd10_pcs) AS code_length,
       count(*) AS n_operations,
       count(DISTINCT subject_id) AS n_patients,
       string_agg(DISTINCT department, '|' ORDER BY department) AS departments
FROM raw.operations
GROUP BY icd10_pcs
ORDER BY n_operations DESC;

DROP TABLE IF EXISTS audit.code_length_distribution;
CREATE TABLE audit.code_length_distribution AS
SELECT length(icd10_pcs) AS code_length, count(DISTINCT icd10_pcs) AS n_distinct_codes,
       count(*) AS n_operations
FROM raw.operations
GROUP BY length(icd10_pcs)
ORDER BY 1;

DROP TABLE IF EXISTS audit.department_distribution;
CREATE TABLE audit.department_distribution AS
SELECT department, count(*) AS n_operations, count(DISTINCT subject_id) AS n_patients,
       count(DISTINCT icd10_pcs) AS n_distinct_pcs
FROM raw.operations
GROUP BY department
ORDER BY n_operations DESC;

-- CTS(Cardio-Thoracic Surgery)科的代码频数(表型参考)
DROP TABLE IF EXISTS audit.cts_code_frequency;
CREATE TABLE audit.cts_code_frequency AS
SELECT icd10_pcs, count(*) AS n_operations, count(DISTINCT subject_id) AS n_patients
FROM raw.operations
WHERE department = 'CTS'
GROUP BY icd10_pcs
ORDER BY n_operations DESC;

-- ---------- 时间锚点审计 ----------
DROP TABLE IF EXISTS audit.operation_time;
CREATE TABLE audit.operation_time AS
SELECT
  count(*) AS n,
  count(*) FILTER (WHERE orin_time IS NOT NULL AND orout_time IS NOT NULL AND orout_time < orin_time) AS n_negative_or_duration,
  count(*) FILTER (WHERE opstart_time IS NOT NULL AND opend_time IS NOT NULL AND opend_time < opstart_time) AS n_negative_op_duration,
  count(*) FILTER (WHERE anstart_time IS NOT NULL AND anend_time IS NOT NULL AND anend_time < anstart_time) AS n_negative_an_duration,
  count(*) FILTER (WHERE icuin_time IS NOT NULL AND icuout_time IS NOT NULL AND icuout_time < icuin_time) AS n_negative_icu_duration,
  count(*) FILTER (WHERE discharge_time IS NOT NULL AND admission_time IS NOT NULL AND discharge_time < admission_time) AS n_negative_hosp_duration,
  -- 相互时间顺序
  count(*) FILTER (WHERE orin_time IS NOT NULL AND opstart_time IS NOT NULL AND opstart_time < orin_time) AS n_opstart_before_orin,
  count(*) FILTER (WHERE opend_time IS NOT NULL AND orout_time IS NOT NULL AND orout_time < opend_time) AS n_orout_before_opend,
  count(*) FILTER (WHERE anstart_time IS NOT NULL AND orin_time IS NOT NULL AND orin_time < anstart_time) AS n_orin_before_anstart,
  count(*) FILTER (WHERE opend_time IS NOT NULL AND icuin_time IS NOT NULL AND icuin_time < opend_time) AS n_icuin_before_opend,
  count(*) FILTER (WHERE orout_time IS NOT NULL AND icuin_time IS NOT NULL AND icuin_time < orout_time) AS n_icuin_before_orout,
  -- 异常持续时间(>24h)
  count(*) FILTER (WHERE orout_time - orin_time > 1440) AS n_or_duration_over_24h,
  count(*) FILTER (WHERE opend_time - opstart_time > 1440) AS n_op_duration_over_24h,
  count(*) FILTER (WHERE anend_time - anstart_time > 1440) AS n_an_duration_over_24h,
  -- 锚点可用率
  round(100*avg((opend_time IS NOT NULL)::int)::numeric,2) AS opend_available_pct,
  round(100*avg((orout_time IS NOT NULL)::int)::numeric,2) AS orout_available_pct,
  round(100*avg((icuin_time IS NOT NULL)::int)::numeric,2) AS icuin_available_pct
FROM raw.operations;

DROP TABLE IF EXISTS audit.time_anchor_by_dept;
CREATE TABLE audit.time_anchor_by_dept AS
SELECT department,
       count(*) AS n,
       round(100*avg((opend_time IS NOT NULL)::int)::numeric,2) AS opend_avail_pct,
       round(100*avg((orout_time IS NOT NULL)::int)::numeric,2) AS orout_avail_pct,
       round(100*avg((icuin_time IS NOT NULL)::int)::numeric,2) AS icuin_avail_pct,
       round(percentile_cont(0.5) WITHIN GROUP (ORDER BY opend_time - opstart_time)::numeric,1) AS op_duration_median_min,
       round(percentile_cont(0.5) WITHIN GROUP (ORDER BY orout_time - orin_time)::numeric,1) AS or_duration_median_min
FROM raw.operations
GROUP BY department
ORDER BY n DESC;

COMMIT;
