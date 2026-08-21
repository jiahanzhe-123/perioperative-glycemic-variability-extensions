-- 12_bloodonly_glucose_patch.sql
-- =============================================================================
-- Blood-only glucose patch(final statistical freeze, v1.0.0)
--
-- 目的:从实验室血糖序列中剔除非血液来源记录。冻结审计发现原 glucose 序列
-- 混入了约 1.21% 的尿糖记录(labevents itemid 51478, fluid = 'Urine'),
-- 影响 283 名患者、约 297 条 24h 窗口内记录,其中 1 名患者剔除后测量数 <2。
--
-- 使用方法:在 05_extract_labs.sql 之后运行本脚本,并以
--   mimic_custom.glucose_lab_values_bloodonly
-- 替换后续所有依赖 mimic_custom.glucose_lab_values_v2 的输入
-- (glucose_long_v2 及其下游 summary/features)。
-- 然后按 docs/final_glucose_variable_dictionary.csv 与
-- python/build_bloodonly_features.py 重建 stay 级特征。
--
-- 验证目标(与冻结清单不一致时不得继续):
--   (a) 排除的 24h 窗口记录 ≈ 297 条,涉及 283 个 stay;
--   (b) 最终序列中 itemid 51478 且 fluid='Urine' 的记录数 = 0;
--   (c) 1 个 stay 在排除后 24h 窗口测量数 <2。
--
-- 说明:本文件按冻结报告记录的确定性规则重建,未随仓库附带数据库重跑;
-- 请以上述验证目标为准进行核对。
-- =============================================================================

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

-- 1. 血液来源实验室血糖:labevent_id 回连 labevents 取 fluid,仅保留 Blood。
--    glucose_lab_values_v2.source_id = labevents.labevent_id(见 05_extract_labs.sql)。
DROP TABLE IF EXISTS mimic_custom.glucose_lab_values_bloodonly;
CREATE TABLE mimic_custom.glucose_lab_values_bloodonly AS
SELECT g.subject_id,
       g.hadm_id,
       g.stay_id,
       g.charttime,
       g.source_id,
       g.glucose_mg_dl,
       g.raw_value,
       g.valueuom,
       g.unit_inference_flag,
       g.outlier_flag
FROM mimic_custom.glucose_lab_values_v2 g
JOIN mimiciv_hosp.labevents le
  ON le.labevent_id = g.source_id
WHERE le.fluid = 'Blood'
  AND le.itemid <> 51478;   -- 双保险:Glucose, Urine

CREATE INDEX glucose_lab_values_bloodonly_idx
    ON mimic_custom.glucose_lab_values_bloodonly (stay_id, charttime, glucose_mg_dl);

-- 2. 验证查询(不写入表,结果应与冻结清单一致)
--    (a) 被排除的原始记录总数(全时间范围,仅供参考):
-- SELECT count(*) AS excluded_lab_rows_all_time
-- FROM mimic_custom.glucose_lab_values_v2 g
-- JOIN mimiciv_hosp.labevents le ON le.labevent_id = g.source_id
-- WHERE le.fluid <> 'Blood' OR le.itemid = 51478;
--
--    (b) 24h 窗口内排除记录与受影响 stay(目标 ≈297 条 / 283 stays):
-- SELECT count(*) AS excluded_rows_24h, count(DISTINCT g.stay_id) AS affected_stays
-- FROM mimic_custom.glucose_lab_values_v2 g
-- JOIN mimiciv_hosp.labevents le ON le.labevent_id = g.source_id
-- JOIN mimic_custom.cardiac_surgery_cohort_v2 c ON c.stay_id = g.stay_id
-- WHERE (le.fluid <> 'Blood' OR le.itemid = 51478)
--   AND g.charttime >= c.surgery_time
--   AND g.charttime <  c.surgery_time + interval '24 hours';
--
--    (c) 最终序列中不得存在尿糖:
-- SELECT count(*) AS urine_rows_remaining   -- 期望 0
-- FROM mimic_custom.glucose_lab_values_bloodonly g
-- JOIN mimiciv_hosp.labevents le ON le.labevent_id = g.source_id
-- WHERE le.fluid <> 'Blood' OR le.itemid = 51478;

COMMIT;
