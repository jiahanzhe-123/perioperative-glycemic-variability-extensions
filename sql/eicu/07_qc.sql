-- Step 7: QC 检查
\pset footer off
\echo '=== QC1 每名患者仅一个 stay ==='
SELECT count(*) AS total_stays, count(DISTINCT uniquepid) AS unique_patients
FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND rn_patient=1;
\echo '=== QC2/3 暴露窗口 0-1440 且无入科前血糖 ==='
SELECT min(offset_min) AS min_off, max(offset_min) AS max_off FROM eicu_cardio_validation.glucose_clean;
\echo '=== QC4 血糖单位 ==='
SELECT unit_orig, count(*) FROM eicu_cardio_validation.glucose_clean GROUP BY 1;
\echo '=== QC6/7 同 stay 同 offset 同值的跨源重复(bedside+lab 同一分钟同值) ==='
SELECT count(*) AS cross_source_same_minute_same_value FROM (
  SELECT patientunitstayid, offset_min, value_mgdl
  FROM eicu_cardio_validation.glucose_clean
  GROUP BY 1,2,3 HAVING count(DISTINCT source_name) > 1
) x;
\echo '=== QC8/9 POCT 与中心实验室比例 ==='
SELECT source_name, count(*) AS n, round(100.0*count(*)/sum(count(*)) OVER (),1) AS pct
FROM eicu_cardio_validation.glucose_clean GROUP BY 1;
\echo '=== nursecharting Bedside Glucose 与 lab 的重叠(平行源核查) ==='
WITH nc AS (
  SELECT n.patientunitstayid, n.nursingchartoffset AS off, n.nursingchartvalue AS val
  FROM nursecharting n
  JOIN eicu_cardio_validation.cohort_base c ON c.patientunitstayid=n.patientunitstayid AND c.broad_cohort_flag AND c.rn_patient=1
  WHERE n.nursingchartcelltypevalname='Bedside Glucose' AND n.nursingchartoffset BETWEEN 0 AND 1440
)
SELECT count(*) AS nc_total,
       count(*) FILTER (WHERE EXISTS (
         SELECT 1 FROM eicu_cardio_validation.glucose_clean g
         WHERE g.patientunitstayid=nc.patientunitstayid AND g.offset_min=nc.off AND g.value_mgdl::text=nc.val)) AS nc_exact_dup_in_lab
FROM nc;
\echo '=== QC10 高 GV 是否由极端值驱动: sd 与 range 的相关 + 最大单跳占比 ==='
SELECT round(corr(glucose_sd_24h, glucose_range_24h)::numeric,3) AS corr_sd_range,
       round(corr(glucose_sd_24h, glucose_n)::numeric,3) AS corr_sd_n,
       round(corr(glucose_n, glucose_mean_24h)::numeric,3) AS corr_n_mean
FROM eicu_cardio_validation.glucose_features WHERE glucose_n>=2;
\echo '=== QC12 glucose_n>=2 筛选对死亡率影响 ==='
SELECT (g.glucose_n>=2) AS ge2, count(*) n,
       round(avg(o.hospital_mortality::int)*100,2) AS hosp_mort_pct
FROM eicu_cardio_validation.outcomes o
JOIN eicu_cardio_validation.glucose_features g USING (patientunitstayid)
GROUP BY 1;
SELECT count(*) n, round(avg(o.hospital_mortality::int)*100,2) AS hosp_mort_pct
FROM eicu_cardio_validation.outcomes o
WHERE NOT EXISTS (SELECT 1 FROM eicu_cardio_validation.glucose_features g WHERE g.patientunitstayid=o.patientunitstayid);
\echo '=== QC13 high-spec vs broad 死亡率 ==='
SELECT c.high_specificity_flag, count(*) n,
       round(avg(o.hospital_mortality::int)*100,2) AS hosp_mort_pct
FROM eicu_cardio_validation.cohort_base c
JOIN eicu_cardio_validation.outcomes o USING (patientunitstayid)
WHERE c.broad_cohort_flag AND c.rn_patient=1
GROUP BY 1;
\echo '=== 医院分布: 病例数过少的医院 ==='
SELECT count(DISTINCT hospitalid) AS hospitals FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND rn_patient=1;
SELECT count(*) AS hospitals_lt10_cases FROM (
  SELECT hospitalid FROM eicu_cardio_validation.cohort_base WHERE broad_cohort_flag AND rn_patient=1
  GROUP BY 1 HAVING count(*) < 10
) x;
