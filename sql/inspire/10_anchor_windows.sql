-- 10_anchor_windows.sql — 多锚点/多窗血糖序列与特征(协议 §7 冻结定义)
-- 窗口:opend 0-24h(已有,主)、opstart 0-24h、orout 0-24h、icuin 0-24h、opend 0-12h、opend 0-48h。
BEGIN;

DROP TABLE IF EXISTS derived.glucose_long_anchor;
CREATE TABLE derived.glucose_long_anchor AS
SELECT c.subject_id, w.anchor, w.t0, l.chart_time, l.value,
       l.chart_time - w.t0 AS min_offset
FROM derived.cohort_first_operation c
CROSS JOIN LATERAL (VALUES
    ('opend_0_24h',  c.opend_time,  1440),
    ('opstart_0_24h',c.opstart_time,1440),
    ('orout_0_24h',  c.orout_time,  1440),
    ('icuin_0_24h',  c.icuin_time,  1440),
    ('opend_0_12h',  c.opend_time,   720),
    ('opend_0_48h',  c.opend_time,  2880)
) AS w(anchor, t0, span)
JOIN raw.labs l ON l.subject_id = c.subject_id AND l.item_name = 'glucose'
WHERE w.t0 IS NOT NULL
  AND l.value IS NOT NULL AND l.value BETWEEN 20 AND 1500
  AND l.chart_time >= w.t0 AND l.chart_time < w.t0 + w.span;

-- 去重(同锚点同患者同分钟中位数)与同分钟审计
DROP TABLE IF EXISTS derived.glucose_minute_anchor;
CREATE TABLE derived.glucose_minute_anchor AS
SELECT subject_id, anchor, chart_time, min_offset,
       percentile_cont(0.5) WITHIN GROUP (ORDER BY value) AS value,
       count(*) AS n_records_in_minute
FROM derived.glucose_long_anchor
GROUP BY subject_id, anchor, chart_time, min_offset;

DROP TABLE IF EXISTS derived.glucose_features_anchor;
CREATE TABLE derived.glucose_features_anchor AS
WITH arv AS (
    SELECT subject_id, anchor,
           round(avg(abs_diff)::numeric, 3) AS arv
    FROM (
        SELECT subject_id, anchor,
               abs(value - lag(value) OVER (PARTITION BY subject_id, anchor ORDER BY chart_time)) AS abs_diff
        FROM derived.glucose_minute_anchor
    ) x WHERE abs_diff IS NOT NULL
    GROUP BY subject_id, anchor
),
mad AS (
    SELECT subject_id, anchor,
           round(percentile_cont(0.5) WITHIN GROUP (ORDER BY abs(x.value - x.med))::numeric, 2) AS mad_glucose
    FROM (
        SELECT m.subject_id, m.anchor, m.value, mm.med
        FROM derived.glucose_minute_anchor m
        JOIN (SELECT subject_id, anchor, percentile_cont(0.5) WITHIN GROUP (ORDER BY value) AS med
              FROM derived.glucose_minute_anchor GROUP BY subject_id, anchor) mm
        USING (subject_id, anchor)
    ) x GROUP BY subject_id, anchor
)
SELECT m.subject_id, m.anchor,
       count(*) AS n,
       count(DISTINCT m.value) AS n_distinct_values,
       round(avg(m.value)::numeric, 2) AS mean_glucose,
       round(stddev_samp(m.value)::numeric, 3) AS gv_sd,
       round((percentile_cont(0.75) WITHIN GROUP (ORDER BY m.value) - percentile_cont(0.25) WITHIN GROUP (ORDER BY m.value))::numeric, 2) AS iqr_glucose,
       round((max(m.chart_time)-min(m.chart_time))/60.0, 2) AS span_hours,
       arv.arv, mad.mad_glucose
FROM derived.glucose_minute_anchor m
LEFT JOIN arv USING (subject_id, anchor)
LEFT JOIN mad USING (subject_id, anchor)
GROUP BY m.subject_id, m.anchor, arv.arv, mad.mad_glucose;

-- 窗口重叠审计(opend 主窗 vs 其他窗:相同患者与相同记录比例)
DROP TABLE IF EXISTS audit.anchor_overlap;
CREATE TABLE audit.anchor_overlap AS
SELECT a.anchor,
       count(DISTINCT a.subject_id) AS n_patients,
       count(*) AS n_records,
       count(*) FILTER (WHERE a.subject_id IS NOT NULL AND o.subject_id IS NOT NULL) AS n_records_in_opend_main,
       count(DISTINCT a.subject_id) FILTER (WHERE o.subject_id IS NOT NULL) AS n_patients_in_opend_main,
       round(100*count(*) FILTER (WHERE a.subject_id IS NOT NULL AND o.subject_id IS NOT NULL)::numeric/nullif(count(*),0), 1) AS pct_records_overlap_opend
FROM derived.glucose_long_anchor a
LEFT JOIN derived.glucose_long o
  ON o.subject_id = a.subject_id AND o.chart_time = a.chart_time AND o.value = a.value
GROUP BY a.anchor;

COMMIT;
