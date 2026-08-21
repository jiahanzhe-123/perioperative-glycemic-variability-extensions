-- 05_lab_audit.sql — 血糖/HbA1c 候选项目字典(待人工确认,不自动认定)
BEGIN;

DROP TABLE IF EXISTS audit.lab_all_items;
CREATE TABLE audit.lab_all_items AS
SELECT item_name,
       count(*) AS n_records,
       count(DISTINCT subject_id) AS n_patients,
       count(DISTINCT value) AS n_distinct_values,
       round(min(value)::numeric, 3) AS min_value,
       round(max(value)::numeric, 3) AS max_value,
       round(avg(value)::numeric, 3) AS mean_value,
       round(100*avg((value IS NULL)::int)::numeric, 2) AS pct_null_value,
       min(chart_time) AS min_chart_time,
       max(chart_time) AS max_chart_time
FROM raw.labs
GROUP BY item_name
ORDER BY n_records DESC;

DROP TABLE IF EXISTS audit.lab_glucose_candidates;
CREATE TABLE audit.lab_glucose_candidates AS
SELECT item_name,
       count(*) AS n_records,
       count(DISTINCT subject_id) AS n_patients,
       count(DISTINCT value) AS n_distinct_values,
       round(min(value)::numeric, 3) AS min_value,
       round(percentile_cont(0.5) WITHIN GROUP (ORDER BY value)::numeric, 2) AS median_value,
       round(max(value)::numeric, 3) AS max_value,
       min(chart_time) AS min_chart_time,
       max(chart_time) AS max_chart_time,
       'pending manual confirmation' AS confirmation_status
FROM raw.labs
WHERE lower(item_name) ~ 'gluco|sugar|glyc|a1c|hba1c|hgb ?a1c|hemoglobin a1c|glycated'
GROUP BY item_name
ORDER BY n_records DESC;

-- 与 parameters.csv 对照:parameters 中声明的 glucose/hba1c 是否在数据里
DROP TABLE IF EXISTS audit.lab_parameter_crosscheck;
CREATE TABLE audit.lab_parameter_crosscheck AS
SELECT p.label AS parameter_label, p.unit AS parameter_unit,
       l.item_name IS NOT NULL AS present_in_data,
       count(l.*) AS n_records
FROM meta.parameter_dictionary p
LEFT JOIN raw.labs l ON l.item_name = p.label
WHERE p.table_name = 'labs'
GROUP BY p.label, p.unit, l.item_name IS NOT NULL
ORDER BY present_in_data, p.label;

COMMIT;
