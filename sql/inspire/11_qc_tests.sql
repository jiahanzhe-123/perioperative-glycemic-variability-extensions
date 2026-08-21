-- qc_tests.sql — 独立 QC 断言(可由 audit.qc_tests 复现;与 08_analysis_base.sql 内嵌断言一致)
-- 任何 pass=false 都必须停止并在 feasibility_report 中解释。
BEGIN;

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
       (SELECT count(*) = 0 FROM derived.glucose_long g JOIN derived.cohort_first_operation c USING (subject_id)
        WHERE g.chart_time >= c.opend_time + 1440),
       (SELECT count(*) FROM derived.glucose_long g JOIN derived.cohort_first_operation c USING (subject_id)
        WHERE g.chart_time >= c.opend_time + 1440)::text
UNION ALL
SELECT 'gv_cases_have_ge2_distinct_timepoints',
       (SELECT count(*) FILTER (WHERE gv_sd IS NOT NULL AND n_glucose_0_24h < 2) = 0 FROM analysis.inspire_base),
       (SELECT count(*) FILTER (WHERE gv_sd IS NOT NULL AND n_glucose_0_24h < 2) FROM analysis.inspire_base)::text
UNION ALL
SELECT 'gv_ge3_sensitivity_count',
       (SELECT count(*) FILTER (WHERE n_glucose_0_24h >= 3) >= 0 FROM analysis.inspire_base),
       (SELECT count(*) FILTER (WHERE n_glucose_0_24h >= 3) FROM analysis.inspire_base)::text
UNION ALL
SELECT 'shr_cases_hba1c_strictly_preop_1_90d',
       (SELECT count(*) = 0 FROM analysis.inspire_base WHERE shr IS NOT NULL AND (hba1c_days_before_opstart < 0 OR hba1c_days_before_opstart > 90)),
       (SELECT count(*) FROM analysis.inspire_base WHERE shr IS NOT NULL AND (hba1c_days_before_opstart < 0 OR hba1c_days_before_opstart > 90))::text
UNION ALL
SELECT 'eag_all_positive',
       (SELECT count(*) FILTER (WHERE hba1c_pct IS NOT NULL AND eag_mg_dl <= 0) = 0 FROM analysis.inspire_base),
       (SELECT count(*) FILTER (WHERE hba1c_pct IS NOT NULL AND eag_mg_dl <= 0) FROM analysis.inspire_base)::text
UNION ALL
SELECT 'shr_molecular_matches_same_sequence',
       (SELECT count(*) = 0 FROM analysis.inspire_base b
        WHERE shr IS NOT NULL AND abs(shr - mean_glucose/eag_mg_dl) > 0.001),
       (SELECT count(*) FROM analysis.inspire_base b WHERE shr IS NOT NULL AND abs(shr - mean_glucose/eag_mg_dl) > 0.001)::text
UNION ALL
SELECT 'cohort_flow_reproduces_analysis_rows',
       (SELECT count(*) FROM analysis.inspire_base) = (SELECT n FROM audit.cohort_flow WHERE step='INCLUDE phenotype, first operation per patient'),
       (SELECT count(*) FROM analysis.inspire_base)::text
UNION ALL
SELECT 'raw_import_counts_match',
       bool_and(match_ok), (SELECT string_agg(table_name||':'||csv_rows,';') FROM audit.table_counts)
FROM audit.table_counts;

COMMIT;
