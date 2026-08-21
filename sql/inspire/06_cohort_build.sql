-- 06_cohort_build.sql — 队列构建(表型 v1:OPEN_CABG + OPEN_VALVE 为 INCLUDE;
-- REVIEW_AORTIC_CANDIDATE 挂起为单独 pending 队列;每人仅首次合格手术;多手术审计)
BEGIN;

-- 表型 codebook 入库
CREATE TABLE IF NOT EXISTS meta.cardiac_procedure_codebook (
    icd10_pcs text PRIMARY KEY,
    n_operations bigint, n_patients bigint, departments text,
    action text, proposed_group text, reason text, confidence text, codebook_version text
);
-- codebook 由脚本经 STDIN 装载

-- 合格手术候选(INCLUDE 组)
DROP TABLE IF EXISTS derived.cohort_candidates;
CREATE TABLE derived.cohort_candidates AS
SELECT o.*, cb.proposed_group, cb.confidence AS codebook_confidence
FROM raw.operations o
JOIN meta.cardiac_procedure_codebook cb ON cb.icd10_pcs = o.icd10_pcs
WHERE cb.action = 'INCLUDE';

-- 多手术审计:同一患者多次合格手术
DROP TABLE IF EXISTS audit.multi_operation;
CREATE TABLE audit.multi_operation AS
SELECT subject_id, count(*) AS n_qualifying_operations,
       string_agg(icd10_pcs, '|' ORDER BY opstart_time) AS pcs_sequence,
       string_agg(DISTINCT proposed_group, '|') AS groups
FROM derived.cohort_candidates
GROUP BY subject_id
HAVING count(*) > 1;

-- 每人首次合格手术(opstart_time 最早;并列取 op_id 小者)
DROP TABLE IF EXISTS derived.cohort_first_operation;
CREATE TABLE derived.cohort_first_operation AS
SELECT DISTINCT ON (subject_id) *
FROM derived.cohort_candidates
ORDER BY subject_id, opstart_time NULLS LAST, op_id;

-- pending 主动脉候选队列(独立保存,不进主要表型,待字典确认)
DROP TABLE IF EXISTS derived.cohort_aortic_candidates;
CREATE TABLE derived.cohort_aortic_candidates AS
SELECT DISTINCT ON (o.subject_id) o.*, cb.proposed_group, cb.confidence AS codebook_confidence
FROM raw.operations o
JOIN meta.cardiac_procedure_codebook cb ON cb.icd10_pcs = o.icd10_pcs
WHERE cb.action = 'REVIEW_AORTIC_CANDIDATE'
ORDER BY o.subject_id, o.opstart_time NULLS LAST, o.op_id;

-- 队列流程
DROP TABLE IF EXISTS audit.cohort_flow;
CREATE TABLE audit.cohort_flow AS
SELECT 'all operations (raw.operations)' AS step, count(*) AS n FROM raw.operations
UNION ALL
SELECT 'with cardiac codebook match (any action)',
       (SELECT count(*) FROM raw.operations o JOIN meta.cardiac_procedure_codebook cb ON cb.icd10_pcs=o.icd10_pcs)
UNION ALL
SELECT 'INCLUDE phenotype (OPEN_CABG + OPEN_VALVE), operations', (SELECT count(*) FROM derived.cohort_candidates)
UNION ALL
SELECT 'INCLUDE phenotype, first operation per patient', (SELECT count(*) FROM derived.cohort_first_operation)
UNION ALL
SELECT 'patients with >1 qualifying operation (excluded beyond first)', (SELECT count(*) FROM audit.multi_operation)
UNION ALL
SELECT 'REVIEW_AORTIC_CANDIDATE (pending), first per patient', (SELECT count(*) FROM derived.cohort_aortic_candidates);

COMMIT;
