-- 08_analysis_base.sql — analysis_inspire_base:一手术一行(本步不运行任何结局回归)
BEGIN;

DROP TABLE IF EXISTS analysis.inspire_base;
CREATE TABLE analysis.inspire_base AS
SELECT
    c.subject_id,
    c.op_id,
    c.hadm_id,
    c.proposed_group AS surgery_group,
    c.codebook_confidence,
    c.icd10_pcs,
    -- 人口学与术前
    c.age, c.sex, c.weight, c.height, c.race, c.asa, c.emop, c.department, c.antype,
    -- 时间锚点(分钟,相对时间)
    c.opdate, c.opstart_time, c.opend_time, c.orin_time, c.orout_time,
    c.anstart_time, c.anend_time, c.cpbon_time, c.cpboff_time,
    c.icuin_time, c.icuout_time, c.admission_time, c.discharge_time,
    (c.opend_time - c.opstart_time) AS op_duration_min,
    (c.orout_time - c.orin_time)     AS or_duration_min,
    (c.anend_time - c.anstart_time)  AS an_duration_min,
    (c.discharge_time - c.admission_time) AS hosp_duration_min,
    -- 测量过程
    f.n_0_12h, f.n_0_24h, f.n_0_48h, f.n_icu_0_24h,
    f.n_gv AS n_glucose_0_24h, f.n_distinct_values, f.span_hours, f.median_interval_min,
    -- 暴露(术后 0–24h,opend 锚点)
    f.mean_glucose, f.gv_sd, f.cv_pct, f.median_glucose, f.iqr_glucose, f.mad_glucose,
    f.min_glucose, f.max_glucose, f.arv, f.tw_arv,
    -- HbA1c / SHR
    h.hba1c_pct, h.hba1c_time, h.hba1c_days_before_opstart, h.eag_mg_dl,
    s.shr,
    -- 结局与 landmark
    c.inhosp_death_time, c.allcause_death_time,
    (c.inhosp_death_time IS NOT NULL) AS inhosp_death_flag,
    (c.allcause_death_time IS NOT NULL) AS allcause_death_flag,
    -- landmark = opend_time + 1440(术后 24h):在 landmark 时仍存活且有随访
    -- inhosp 视角:出院或死亡时间晚于 landmark;院外随访用 allcause_death_time(若存在)
    (COALESCE(c.inhosp_death_time, c.discharge_time) >= c.opend_time + 1440) AS landmark_eligible_inhosp,
    (COALESCE(c.allcause_death_time, c.discharge_time) >= c.opend_time + 1440) AS landmark_eligible_allcause,
    -- 结局时间(自 opend 起算分钟,供后续 landmark 分析;不预先计算生存模型)
    (LEAST(COALESCE(c.inhosp_death_time, c.discharge_time), c.opend_time + 30*1440) - c.opend_time) AS t30_from_opend_min,
    (COALESCE(c.inhosp_death_time, c.discharge_time) <= c.opend_time + 30*1440 AND c.inhosp_death_time IS NOT NULL) AS event_30d_inhosp
FROM derived.cohort_first_operation c
LEFT JOIN derived.glucose_features f USING (subject_id)
LEFT JOIN derived.hba1c_baseline h USING (subject_id)
LEFT JOIN derived.shr_features s USING (subject_id);

-- checksum 列(行级)
ALTER TABLE analysis.inspire_base ADD COLUMN row_checksum text;
UPDATE analysis.inspire_base
SET row_checksum = md5(subject_id::text || '|' || op_id::text || '|' || coalesce(icd10_pcs,'') || '|' || coalesce(gv_sd::text,''));

-- ---------- QC 断言(结果写 audit.qc_tests;任何 fail 必须解释) ----------
DROP TABLE IF EXISTS audit.qc_tests;
CREATE TABLE audit.qc_tests AS
SELECT 'one_operation_per_patient' AS test_name,
       (SELECT count(*) = (SELECT count(DISTINCT subject_id) FROM analysis.inspire_base) FROM analysis.inspire_base) AS pass,
       (SELECT count(*) FROM analysis.inspire_base)::text AS detail
UNION ALL
SELECT 'all_exposure_within_opend_0_24h',
       (SELECT count(*) = 0 FROM derived.glucose_long WHERE min_from_opend < 0 OR min_from_opend >= 1440),
       (SELECT count(*) FROM derived.glucose_long WHERE min_from_opend < 0 OR min_from_opend >= 1440)::text
UNION ALL
SELECT 'all_exposure_strictly_before_landmark',
       (SELECT count(*) = 0 FROM derived.glucose_long WHERE chart_time >= (SELECT opend_time FROM derived.cohort_first_operation c WHERE c.subject_id=derived.glucose_long.subject_id LIMIT 1) + 1440),
       'chart_time < opend_time + 1440 enforced by window filter'
UNION ALL
SELECT 'gv_cases_have_ge2_distinct_timepoints',
       (SELECT count(*) FILTER (WHERE gv_sd IS NOT NULL AND n_glucose_0_24h < 2) = 0 FROM analysis.inspire_base),
       (SELECT count(*) FILTER (WHERE gv_sd IS NOT NULL AND n_glucose_0_24h < 2) FROM analysis.inspire_base)::text
UNION ALL
SELECT 'shr_cases_hba1c_strictly_preop_1_90d',
       (SELECT count(*) = 0 FROM analysis.inspire_base WHERE shr IS NOT NULL AND (hba1c_days_before_opstart < 0 OR hba1c_days_before_opstart > 90)),
       (SELECT count(*) FROM analysis.inspire_base WHERE shr IS NOT NULL AND (hba1c_days_before_opstart < 0 OR hba1c_days_before_opstart > 90))::text
UNION ALL
SELECT 'eag_all_positive',
       (SELECT count(*) FILTER (WHERE hba1c_pct IS NOT NULL AND eag_mg_dl <= 0) = 0 FROM analysis.inspire_base),
       (SELECT count(*) FILTER (WHERE hba1c_pct IS NOT NULL AND eag_mg_dl <= 0) FROM analysis.inspire_base)::text
UNION ALL
SELECT 'cohort_flow_reproduces_analysis_rows',
       (SELECT count(*) FROM analysis.inspire_base) = (SELECT n FROM audit.cohort_flow WHERE step='INCLUDE phenotype, first operation per patient'),
       (SELECT count(*) FROM analysis.inspire_base)::text
UNION ALL
SELECT 'raw_import_counts_match',
       bool_and(match_ok), (SELECT string_agg(table_name||':'||csv_rows,';') FROM audit.table_counts)
FROM audit.table_counts;

DROP TABLE IF EXISTS audit.table_checksums;
CREATE TABLE audit.table_checksums AS
SELECT 'analysis.inspire_base' AS object_name, md5(string_agg(row_checksum, '' ORDER BY subject_id)) AS checksum,
       count(*) AS n_rows, now() AS computed_at FROM analysis.inspire_base
UNION ALL
SELECT 'derived.glucose_minute', md5(string_agg(md5(subject_id::text||chart_time::text||value::text), '' ORDER BY subject_id, chart_time)),
       count(*), now() FROM derived.glucose_minute
UNION ALL
SELECT 'derived.cohort_first_operation', md5(string_agg(md5(subject_id::text||op_id::text), '' ORDER BY subject_id)),
       count(*), now() FROM derived.cohort_first_operation;

COMMIT;
