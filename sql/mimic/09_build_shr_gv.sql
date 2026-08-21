-- 09_build_shr_gv.sql
-- Main SHR = mean postoperative 24h glucose / eAG derived from HbA1c.
-- Main GV = postoperative 24h glucose SD. Alternative candidates are retained
-- for sensitivity analysis and transparent phenotype review.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.qc_old_vs_new_shr_gv_v2;
DROP TABLE IF EXISTS mimic_custom.shr_gv_candidates_v2;

CREATE TABLE mimic_custom.shr_gv_candidates_v2 AS
WITH first_postop AS (
    SELECT DISTINCT ON (g.stay_id)
           g.stay_id,
           g.glucose_mg_dl AS glucose_postop_first,
           g.charttime AS glucose_postop_first_time
    FROM mimic_custom.glucose_long_v2 g
    JOIN mimic_custom.cardiac_surgery_cohort_v2 c ON c.stay_id = g.stay_id
    WHERE g.glucose_mg_dl IS NOT NULL
      AND g.charttime >= c.surgery_time
      AND g.charttime < c.surgery_time + interval '24 hours'
    ORDER BY g.stay_id, g.charttime, g.source_id
), first_adm AS (
    SELECT DISTINCT ON (l.stay_id)
           l.stay_id,
           l.value_valid AS glucose_adm_first,
           l.charttime AS glucose_adm_first_time
    FROM mimic_custom.labs_long_v2 l
    JOIN mimic_custom.cardiac_surgery_cohort_v2 c ON c.stay_id = l.stay_id
    WHERE l.variable_name = 'glucose'
      AND l.value_valid IS NOT NULL
      AND l.charttime >= c.admittime_ts
      AND l.charttime < c.admittime_ts + interval '24 hours'
    ORDER BY l.stay_id, l.charttime, l.labevent_id
), base AS (
    SELECT c.subject_id,
           c.hadm_id,
           c.stay_id,
           h.hba1c_pct,
           h.hba1c_time,
           h.hba1c_days_from_surgery,
           h.hba1c_source,
           h.hba1c_unit_inference_flag,
           h.hba1c_window,
           g.n_glucose_postop_24h,
           g.n_glucose_postop_48h,
           g.n_glucose_postop_72h,
           g.n_glucose_postop_7d,
           g.glucose_mean_postop_24h,
           g.glucose_sd_postop_24h,
           g.glucose_min_postop_24h,
           g.glucose_max_postop_24h,
           g.glucose_median_postop_24h,
           g.glucose_cv_postop_24h,
           g.glucose_range_postop_24h,
           g.glucose_variability_delta_postop_24h,
           g.glucose_masd_postop_24h,
           g.glucose_mean_postop_48h,
           g.glucose_mean_postop_72h,
           g.glucose_mean_postop_7d,
           g.glucose_source_mix,
           g.glucose_unit_inference_flag,
           g.glucose_outlier_count,
           fp.glucose_postop_first,
           fp.glucose_postop_first_time,
           fa.glucose_adm_first,
           fa.glucose_adm_first_time,
           (28.7 * h.hba1c_pct - 46.7) AS eAG_mg_dl
    FROM mimic_custom.cardiac_surgery_cohort_v2 c
    LEFT JOIN mimic_custom.hba1c_baseline_v2 h ON h.stay_id = c.stay_id
    LEFT JOIN mimic_custom.glucose_summary_v2 g ON g.stay_id = c.stay_id
    LEFT JOIN first_postop fp ON fp.stay_id = c.stay_id
    LEFT JOIN first_adm fa ON fa.stay_id = c.stay_id
)
SELECT subject_id,
       hadm_id,
       stay_id,
       eAG_mg_dl,
       glucose_postop_first,
       glucose_adm_first,
       glucose_mean_postop_24h,
       glucose_max_postop_24h,
       glucose_min_postop_24h,
       glucose_median_postop_24h,
       glucose_variability_delta_postop_24h,
       glucose_sd_postop_24h,
       glucose_cv_postop_24h,
       glucose_range_postop_24h,
       glucose_masd_postop_24h,
       glucose_mean_postop_48h,
       glucose_mean_postop_72h,
       glucose_mean_postop_7d,
       n_glucose_postop_24h,
       n_glucose_postop_48h,
       n_glucose_postop_72h,
       n_glucose_postop_7d,
       glucose_source_mix,
       glucose_unit_inference_flag,
       glucose_outlier_count,
       hba1c_pct,
       hba1c_time,
       hba1c_days_from_surgery,
       hba1c_source,
       hba1c_unit_inference_flag,
       hba1c_window,
       CASE WHEN eAG_mg_dl IS NULL THEN NULL
            WHEN eAG_mg_dl <= 0 THEN NULL
            ELSE glucose_postop_first / eAG_mg_dl END AS shr_first_postop_24h,
       CASE WHEN eAG_mg_dl IS NULL OR eAG_mg_dl <= 0 THEN NULL
            ELSE glucose_mean_postop_24h / eAG_mg_dl END AS shr_mean_postop_24h,
       CASE WHEN eAG_mg_dl IS NULL OR eAG_mg_dl <= 0 THEN NULL
            ELSE glucose_median_postop_24h / eAG_mg_dl END AS shr_median_postop_24h,
       CASE WHEN eAG_mg_dl IS NULL OR eAG_mg_dl <= 0 THEN NULL
            ELSE glucose_max_postop_24h / eAG_mg_dl END AS shr_max_postop_24h,
       CASE WHEN eAG_mg_dl IS NULL OR eAG_mg_dl <= 0 THEN NULL
            ELSE glucose_mean_postop_48h / eAG_mg_dl END AS shr_mean_postop_48h,
       CASE WHEN eAG_mg_dl IS NULL OR eAG_mg_dl <= 0 THEN NULL
            ELSE glucose_mean_postop_72h / eAG_mg_dl END AS shr_mean_postop_72h,
       CASE WHEN eAG_mg_dl IS NULL OR eAG_mg_dl <= 0 THEN NULL
            ELSE glucose_mean_postop_7d / eAG_mg_dl END AS shr_mean_postop_7d,
       CASE WHEN eAG_mg_dl IS NULL OR eAG_mg_dl <= 0 THEN NULL
            ELSE glucose_adm_first / eAG_mg_dl END AS shr_adm_first,
       CASE WHEN hba1c_pct IS NULL THEN 'missing_hba1c'
            WHEN eAG_mg_dl <= 0 THEN 'invalid_eAG_nonpositive'
            WHEN glucose_mean_postop_24h IS NULL THEN 'missing_glucose_mean_postop_24h'
            ELSE NULL END AS shr_missing_reason,
       CASE WHEN glucose_sd_postop_24h IS NULL THEN 'missing_or_single_postop_24h_glucose'
            ELSE NULL END AS gv_missing_reason,
       glucose_mean_postop_24h / NULLIF(eAG_mg_dl,0) AS shr,
       glucose_sd_postop_24h AS gv,
       'SHR = glucose_mean_postop_24h / (28.7*hba1c_pct - 46.7)'::text AS shr_formula_main,
       'GV = glucose_sd_postop_24h; CV/range/MASD retained as alternatives'::text AS gv_formula_main,
       glucose_source_mix AS glucose_source,
       'mimiciv_hosp.labevents + mimiciv_icu.chartevents; HbA1c pre_90d/pre_365d/post_24h hierarchy'::text AS exposure_source
FROM base;

CREATE INDEX shr_gv_candidates_v2_stay_idx
    ON mimic_custom.shr_gv_candidates_v2 (stay_id);

CREATE TABLE mimic_custom.qc_old_vs_new_shr_gv_v2 (
    old_table_name text,
    old_shr_column text,
    old_gv_column text,
    matched_n bigint,
    pearson_shr double precision,
    spearman_shr double precision,
    pearson_gv double precision,
    spearman_gv double precision,
    status text,
    notes text
);

-- Dynamic comparison is restricted to mimic_custom tables that contain all
-- three identifiers and columns whose names include SHR and GV. If none exist,
-- the table remains empty and the QC renderer states that no old table was found.
DO $$
DECLARE
    r record;
    q text;
BEGIN
    FOR r IN
        WITH cols AS (
            SELECT table_schema,
                   table_name,
                   max(column_name) FILTER (WHERE lower(column_name) ~ 'shr') AS shr_column,
                   max(column_name) FILTER (WHERE lower(column_name) ~ '(^|_)gv($|_)|variability') AS gv_column,
                   bool_or(column_name = 'subject_id') AS has_subject,
                   bool_or(column_name = 'hadm_id') AS has_hadm,
                   bool_or(column_name = 'stay_id') AS has_stay
            FROM information_schema.columns
            WHERE table_schema = 'mimic_custom'
            GROUP BY table_schema, table_name
        )
        SELECT * FROM cols
        WHERE shr_column IS NOT NULL
          AND gv_column IS NOT NULL
          AND has_subject AND has_hadm AND has_stay
          AND table_name NOT IN ('shr_gv_candidates_v2','qc_old_vs_new_shr_gv_v2')
    LOOP
        BEGIN
            q := format($fmt$
                WITH m AS (
                    SELECT n.shr::double precision AS new_shr,
                           n.gv::double precision AS new_gv,
                           o.%1$I::double precision AS old_shr,
                           o.%2$I::double precision AS old_gv
                    FROM mimic_custom.shr_gv_candidates_v2 n
                    JOIN mimic_custom.%3$I o
                      ON o.subject_id = n.subject_id
                     AND o.hadm_id = n.hadm_id
                     AND o.stay_id = n.stay_id
                    WHERE n.shr IS NOT NULL AND n.gv IS NOT NULL
                      AND o.%1$I IS NOT NULL AND o.%2$I IS NOT NULL
                ), rnk AS (
                    SELECT *,
                           rank() OVER (ORDER BY old_shr) AS old_shr_rank,
                           rank() OVER (ORDER BY new_shr) AS new_shr_rank,
                           rank() OVER (ORDER BY old_gv) AS old_gv_rank,
                           rank() OVER (ORDER BY new_gv) AS new_gv_rank
                    FROM m
                )
                INSERT INTO mimic_custom.qc_old_vs_new_shr_gv_v2
                    (old_table_name, old_shr_column, old_gv_column, matched_n,
                     pearson_shr, spearman_shr, pearson_gv, spearman_gv, status, notes)
                SELECT %4$L, %1$L, %2$L, count(*)::bigint,
                       corr(old_shr,new_shr), corr(old_shr_rank::double precision,new_shr_rank::double precision),
                       corr(old_gv,new_gv), corr(old_gv_rank::double precision,new_gv_rank::double precision),
                       'ok', 'ID join subject_id+hadm_id+stay_id'
                FROM rnk
            $fmt$, r.shr_column, r.gv_column, r.table_name, r.table_schema || '.' || r.table_name);
            EXECUTE q;
        EXCEPTION WHEN OTHERS THEN
            INSERT INTO mimic_custom.qc_old_vs_new_shr_gv_v2
                (old_table_name, old_shr_column, old_gv_column, status, notes)
            VALUES (r.table_schema || '.' || r.table_name, r.shr_column, r.gv_column,
                    'error', SQLERRM);
        END;
    END LOOP;
END $$;

COMMIT;

\copy (SELECT * FROM mimic_custom.qc_old_vs_new_shr_gv_v2 ORDER BY old_table_name) TO 'outputs/qc/qc_old_vs_new_shr_gv_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
