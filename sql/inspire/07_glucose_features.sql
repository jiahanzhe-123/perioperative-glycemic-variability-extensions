-- 07_glucose_features.sql — 血糖长表(术后 0–24h)与 GV/SHR 特征
-- 锚点:opend_time(主);orout_time(敏感性)。单位:parameters.csv glucose=mg/dL。
BEGIN;

-- ---------- 血糖长表 ----------
DROP TABLE IF EXISTS derived.glucose_long;
CREATE TABLE derived.glucose_long AS
WITH g AS (
    SELECT c.subject_id, c.op_id, c.opstart_time, c.opend_time, c.orout_time, c.icuin_time,
           l.chart_time, l.value,
           l.chart_time - c.opend_time  AS min_from_opend,
           l.chart_time - c.opstart_time AS min_from_opstart,
           CASE WHEN c.icuin_time IS NOT NULL THEN l.chart_time - c.icuin_time END AS min_from_icuin
    FROM derived.cohort_first_operation c
    JOIN raw.labs l ON l.subject_id = c.subject_id AND l.item_name = 'glucose'
    WHERE l.value IS NOT NULL AND l.value BETWEEN 20 AND 1500
      AND l.chart_time >= c.opend_time
      AND l.chart_time <  c.opend_time + 1440
)
SELECT *, (min_from_opend >= 0 AND min_from_opend < 720)  AS win_0_12h,
       (min_from_opend >= 0 AND min_from_opend < 1440) AS win_0_24h,
       (min_from_opend >= 0 AND min_from_opend < 2880) AS win_0_48h,
       (min_from_icuin IS NOT NULL AND min_from_icuin >= 0 AND min_from_icuin < 1440) AS win_icu_0_24h
FROM g;

-- ---------- 去重:完全重复剔除,同分钟取中位数 ----------
DROP TABLE IF EXISTS audit.glucose_dedup;
CREATE TABLE audit.glucose_dedup AS
SELECT
  (SELECT count(*) FROM derived.glucose_long) AS n_raw_records,
  (SELECT count(*) FROM (SELECT DISTINCT subject_id, chart_time, value FROM derived.glucose_long) t) AS n_distinct_records,
  (SELECT count(*) FROM derived.glucose_long) - (SELECT count(*) FROM (SELECT DISTINCT subject_id, chart_time, value FROM derived.glucose_long) t) AS n_exact_duplicates_removed,
  (SELECT count(*) FROM (SELECT subject_id, chart_time FROM derived.glucose_long GROUP BY 1,2 HAVING count(*)>1) t) AS n_same_minute_multi,
  (SELECT count(DISTINCT subject_id) FROM derived.glucose_long) AS n_patients;

DROP TABLE IF EXISTS audit.glucose_resolution;
CREATE TABLE audit.glucose_resolution AS
SELECT subject_id, chart_time, count(*) AS n_values,
       min(value) AS min_v, max(value) AS max_v, max(value)-min(value) AS abs_diff
FROM derived.glucose_long
GROUP BY subject_id, chart_time
HAVING count(*) > 1;

DROP TABLE IF EXISTS derived.glucose_minute;
CREATE TABLE derived.glucose_minute AS
SELECT subject_id, op_id, chart_time, min_from_opend, min_from_opstart, min_from_icuin,
       percentile_cont(0.5) WITHIN GROUP (ORDER BY value) AS value,
       count(*) AS n_records_in_minute
FROM derived.glucose_long
GROUP BY subject_id, op_id, chart_time, min_from_opend, min_from_opstart, min_from_icuin;

-- ---------- 窗口特征 ----------
DROP TABLE IF EXISTS derived.glucose_features;
CREATE TABLE derived.glucose_features AS
WITH w AS (
    SELECT subject_id,
           count(*) FILTER (WHERE min_from_opend < 720)  AS n_0_12h,
           count(*) FILTER (WHERE min_from_opend < 1440) AS n_0_24h,
           count(*) FILTER (WHERE min_from_opend < 2880) AS n_0_48h,
           count(*) FILTER (WHERE min_from_icuin IS NOT NULL AND min_from_icuin >= 0 AND min_from_icuin < 1440) AS n_icu_0_24h
    FROM derived.glucose_minute GROUP BY subject_id
),
f24 AS (
    SELECT subject_id,
           count(*) AS n,
           count(DISTINCT value) AS n_distinct_values,
           round((max(chart_time)-min(chart_time))/60.0, 2) AS span_hours,
           round(percentile_cont(0.5) WITHIN GROUP (ORDER BY d)::numeric, 1) AS median_interval_min,
           round(avg(value)::numeric, 2) AS mean_glucose,
           round(stddev_samp(value)::numeric, 3) AS gv_sd,
           round((stddev_samp(value)/avg(value)*100)::numeric, 2) AS cv_pct,
           round(percentile_cont(0.5) WITHIN GROUP (ORDER BY value)::numeric, 2) AS median_glucose,
           round((percentile_cont(0.75) WITHIN GROUP (ORDER BY value) - percentile_cont(0.25) WITHIN GROUP (ORDER BY value))::numeric, 2) AS iqr_glucose,
           min(value) AS min_glucose, max(value) AS max_glucose
    FROM (
        SELECT subject_id, value, chart_time,
               chart_time - lag(chart_time) OVER (PARTITION BY subject_id ORDER BY chart_time) AS d
        FROM derived.glucose_minute
        WHERE min_from_opend >= 0 AND min_from_opend < 1440
    ) s
    GROUP BY subject_id
),
mad AS (
    SELECT subject_id,
           round(percentile_cont(0.5) WITHIN GROUP (ORDER BY abs(x.value - x.med))::numeric, 2) AS mad_glucose
    FROM (
        SELECT m.subject_id, m.value, mm.med
        FROM derived.glucose_minute m
        JOIN (
            SELECT subject_id, percentile_cont(0.5) WITHIN GROUP (ORDER BY value) AS med
            FROM derived.glucose_minute
            WHERE min_from_opend >= 0 AND min_from_opend < 1440
            GROUP BY subject_id
        ) mm USING (subject_id)
        WHERE m.min_from_opend >= 0 AND m.min_from_opend < 1440
    ) x
    GROUP BY subject_id
),
arv AS (
    SELECT subject_id,
           round(avg(abs_diff)::numeric, 3) AS arv,
           round((sum(abs_diff * dt) / nullif(sum(dt),0))::numeric, 3) AS tw_arv
    FROM (
        SELECT subject_id,
               abs(value - lag(value) OVER (PARTITION BY subject_id ORDER BY chart_time)) AS abs_diff,
               (chart_time - lag(chart_time) OVER (PARTITION BY subject_id ORDER BY chart_time))/60.0 AS dt
        FROM derived.glucose_minute
        WHERE min_from_opend >= 0 AND min_from_opend < 1440
    ) x
    WHERE abs_diff IS NOT NULL
    GROUP BY subject_id
)
SELECT w.subject_id, w.n_0_12h, w.n_0_24h, w.n_0_48h, w.n_icu_0_24h,
       f24.n AS n_gv, f24.n_distinct_values, f24.span_hours, f24.median_interval_min,
       f24.mean_glucose, f24.gv_sd, f24.cv_pct, f24.median_glucose, f24.iqr_glucose,
       f24.min_glucose, f24.max_glucose, mad.mad_glucose, arv.arv, arv.tw_arv
FROM w JOIN f24 USING (subject_id) LEFT JOIN arv USING (subject_id) LEFT JOIN mad USING (subject_id);

-- ---------- HbA1c(严格术前 1–90 日,opstart_time 为索引时刻) ----------
DROP TABLE IF EXISTS derived.hba1c_baseline;
CREATE TABLE derived.hba1c_baseline AS
WITH h AS (
    SELECT c.subject_id, l.chart_time, l.value, c.opstart_time,
           (c.opstart_time - l.chart_time)/1440.0 AS days_before_opstart,
           row_number() OVER (PARTITION BY c.subject_id ORDER BY l.chart_time DESC) AS rn
    FROM derived.cohort_first_operation c
    JOIN raw.labs l ON l.subject_id = c.subject_id AND l.item_name = 'hba1c'
    WHERE l.value IS NOT NULL AND l.value BETWEEN 3.0 AND 15.0
      AND l.chart_time <  c.opstart_time
      AND l.chart_time >= c.opstart_time - 90*1440
)
SELECT subject_id, value AS hba1c_pct, chart_time AS hba1c_time,
       round(days_before_opstart::numeric, 2) AS hba1c_days_before_opstart,
       round((28.7*value - 46.7)::numeric, 2) AS eag_mg_dl
FROM h WHERE rn = 1;

-- SHR = mean_glucose(同一术后序列)/ eAG;要求 eAG>0
DROP TABLE IF EXISTS derived.shr_features;
CREATE TABLE derived.shr_features AS
SELECT f.subject_id, h.hba1c_pct, h.hba1c_time, h.hba1c_days_before_opstart, h.eag_mg_dl,
       f.mean_glucose, f.gv_sd, f.n_gv,
       round((f.mean_glucose / h.eag_mg_dl)::numeric, 4) AS shr
FROM derived.glucose_features f
JOIN derived.hba1c_baseline h USING (subject_id)
WHERE h.eag_mg_dl > 0;

-- ---------- HbA1c 可用性审计 ----------
DROP TABLE IF EXISTS audit.hba1c_availability;
CREATE TABLE audit.hba1c_availability AS
SELECT count(*) AS n_cohort,
       count(*) FILTER (WHERE h.subject_id IS NOT NULL) AS n_with_strict_preop_hba1c,
       round(100*count(*) FILTER (WHERE h.subject_id IS NOT NULL)::numeric/count(*), 2) AS pct_with_hba1c,
       count(*) FILTER (WHERE h.eag_mg_dl > 0) AS n_eag_positive
FROM derived.cohort_first_operation c
LEFT JOIN derived.hba1c_baseline h USING (subject_id);

DROP TABLE IF EXISTS audit.hba1c_timing_distribution;
CREATE TABLE audit.hba1c_timing_distribution AS
SELECT 'strict pre-op 1-90d' AS window, count(*) AS n FROM derived.hba1c_baseline
UNION ALL
SELECT 'hba1c any time (cohort)', count(DISTINCT l.subject_id)
FROM derived.cohort_first_operation c JOIN raw.labs l ON l.subject_id=c.subject_id AND l.item_name='hba1c' AND l.value IS NOT NULL;

COMMIT;
