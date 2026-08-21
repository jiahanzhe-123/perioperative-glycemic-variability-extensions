# Analysis Workflow

Execution order for the full pipeline. Modules are independent: each database
module can be run without the others. Every step fails loudly (non-zero exit);
no step modifies raw data.

Conventions:

- All paths resolve through `config/paths.yml` (copied from
  `config/config.example.yml`).
- All analytic constants (windows, plausibility ranges, seeds, scaling) live in
  `config/constants.yml` and are computed once per target cohort, then frozen.
- Model definitions are machine-readable in `config/model_specifications.yml`
  and must agree with `docs/RESULT_PROVENANCE.md` and the manuscript.

## 0. Setup and checks (no clinical data required)

```bash
make lint        # sensitive-info scan + config/YAML checks
make test        # pytest + testthat unit tests
make synthetic   # simulated end-to-end workflow
```

## 1. MIMIC-IV module (primary)

Prerequisite: MIMIC-IV v3.1 loaded in PostgreSQL.

1. `sql/mimic/00_schema_inventory.sql` — schema inventory
2. `sql/mimic/01_build_cardiac_surgery_codebook.sql` — versioned procedure
   codebook (six categories; high-specificity open cardiac surgery)
3. `sql/mimic/02_build_cohort.sql` — index admission/stay cohort
4. `sql/mimic/03_build_itemid_review.sql` — laboratory item review
5. `sql/mimic/04_extract_vitals_24h.sql` — vitals
6. `sql/mimic/05_extract_labs.sql` — laboratory long tables (glucose, HbA1c,
   lactate, creatinine, hemoglobin, platelets, WBC, albumin)
7. `sql/mimic/06_extract_comorbidities_bmi.sql` — comorbidities, BMI
8. `sql/mimic/07_extract_treatments.sql` — treatments
9. `sql/mimic/08_build_severity.sql` — severity scores
10. `sql/mimic/09_build_shr_gv.sql` — SHR/GV inputs
11. `sql/mimic/10_build_final_features.sql` — final feature table
12. `sql/mimic/11_build_qc_report.sql` — extraction QC
13. `sql/mimic/12_bloodonly_glucose_patch.sql` — blood-only glucose patch
14. `analyses/02_glucose_processing/01_rebuild_glucose_series.py` — glucose
    series cleaning: plausibility 20–1500 mg/dL, exact-duplicate removal,
    2-minute/1 mg/dL cross-table near-duplicate suppression, source-class
    assignment, stay-minute/source-class medians, priority hierarchy
    (central laboratory > blood gas > POCT > ICU charted)
15. `analyses/01_cohort_construction/02_build_analysis_dataset.py` — analysis
    frames (date-anchored 24 h window, day-1 landmark, ≥2 measurements,
    alive at landmark)
16. `analyses/03_primary_mimic/03_run_primary_models_mice.R` — MICE
    (m=50, 20 iterations, PMM, seed 20260726) + Models A/B/C, 30-day and
    365-day; **Model B at 30 days is the unique primary analysis**
17. `analyses/03_primary_mimic/04_ph_rcs_absoluterisk.R` — PH diagnostics,
    RCS, standardized absolute risk
18. `analyses/03_primary_mimic/05_sensitivity_models.R` — sensitivity models
19. `analyses/03_primary_mimic/06_model_performance.R` — performance metrics
20. `analyses/03_primary_mimic/supplement/` — frame-fixed supplement refits
    (alternative GV metrics, measurement process, diabetes interaction,
    extremes, fixed-scale display constants)

## 2. Time-dependent module (MIMIC, 365-day)

1. `analyses/04_time_dependent/01_rcs_spline_framefix.R` — spline refit
2. `analyses/04_time_dependent/02_figures_framefix.R` — main figures
3. `analyses/04_time_dependent/03_figures_fixedscale.R` — fixed-scale figures
4. `analyses/04_time_dependent/04_rcs_curves.R` — RCS curves
5. `analyses/04_time_dependent/05_mice_pooled_interval_models.R` — **365d 区间
   估计(d1–7/d8–30/d31–365)+ GV×log(time)+ 30d tt 校正,MICE m=50 池化
   (analysis of record;含 sanity 复现断言)**
6. `analyses/04_time_dependent/09_covariate_tt_sensitivity.R` — 365d 协变量
   log(time) 敏感性(D1/D2),MICE 池化

Time-axis labels must distinguish the index-day axis from the post-landmark
axis; the post-landmark day 31–365 interval estimate is not a protective
effect (see manuscript).

## 3. Source sensitivity (MIMIC)

1. `analyses/05_source_sensitivity/01_build_source_series.py`
2. `analyses/05_source_sensitivity/02_fit_source_models.R`

## 4. SHR–GV joint module (MIMIC subset)

1. `analyses/08_shr_component/01_shr_hba1c_ipw.R` — HbA1c subset + IPW audit
2. `analyses/08_shr_component/02_joint_prep.py` — joint frame preparation
3. `analyses/08_shr_component/03_joint_models.R` — J0/J1, D1/D2 models and
   pooled Wald tests
4. `analyses/08_shr_component/04_joint_risk_performance.R` — standardized
   risk + paired bootstrap performance
5. `analyses/08_shr_component/05_joint_source_sensitivity.R`
6. `analyses/08_shr_component/06_joint_qc_assembly.py`

## 5. INSPIRE module (secondary cohort analysis; Level 3 only)

Prerequisite: INSPIRE v1.4.2 in PostgreSQL (authorized environment).

1. `sql/inspire/00_create_database.sql`, `02_indexes.sql`,
   `03_schema_audit.sql`, `04_procedure_audit.sql`, `05_lab_audit.sql`,
   `06_cohort_build.sql`, `07_glucose_features.sql`, `08_analysis_base.sql`,
   `09_comorbidity_outcome_audit.sql`, `10_anchor_windows.sql`,
   `11_qc_tests.sql`
2. `analyses/06_inspire/01_import_tables.sh`
3. `analyses/06_inspire/02_outcome_time_repairs.R` — reconciled mortality
   outcome; common postoperative-day-30 administrative censoring
4. `analyses/06_inspire/03_sensitivity_joint_repairs.R`
5. `analyses/06_inspire/04_final_models_v3.R` — corrected outcome/model frame
   (opend +24 h landmark; correct 48 h exposure → 48 h landmark; no
   exposure–risk overlap)
6. `analyses/06_inspire/05_uniform_admin_censoring_v5.R` — analysis of record
   for the reported 30-day and corrected 48-hour results; common day-30
   administrative censoring, never discharge censoring
7. `analyses/06_inspire/05_qc_assembly.py`

## 6. eICU-CRD module (harmonized comparison)

Prerequisite: eICU-CRD v2.0 in PostgreSQL.

1. `sql/eicu/` — `01_schema_audit` … `09_harmonized` in numeric order
2. `analyses/07_eicu/01_fit_harmonized_models.R` — M1–M4 modified Poisson
   with hospital-clustered robust variance
3. `analyses/07_eicu/02_harmonized_mimic_features.py`,
   `03_harmonized_mimic_model_data.py`, `04_models.R` — MIMIC-side
   harmonization for descriptive comparison
4. `analyses/07_eicu/05_tables_figures.R`, `06_qc.R`

## 7. Assembly, QC and figures

1. `analyses/09_quality_control/01_assemble_results.py`
2. `analyses/09_quality_control/02_missing_data_report.R`(MICE 诊断图)
3. `analyses/09_quality_control/03_frame_registry.R`
4. `analyses/09_quality_control/figure_build/`(最终绘图,全部配置驱动):
   - `build_main_figures.py` — Main F1–F5 + Supplementary S7/S8;
   - `build_fig23_main_figures.py` — Main F2/F3(投稿包版本,覆盖同名输出);
   - `build_supplement_figures.py` — S6 分布、S2(旧内部编号 S12)模型序列;
   - `build_ph_audit_figure.py` — S5 PH 审计;
   - `merge_s1_mice_figure.py` — S1 MICE 诊断纵向合成。
   输出统一写入 `results/figures/`(gitignored)。
5. `make validate` — frozen-result verification against
   `tests/expected/frozen_results.yml`

## Frozen constants

- MICE: m=50, 20 iterations, PMM, seed 20260726
- Glucose plausibility: 20–1500 mg/dL
- Effect scale: per 10 mg/dL (= 0.555 mmol/L)
- MIMIC window: date-anchored first 24 h; day-1 landmark
- INSPIRE window: [opend, opend+24 h); landmark opend+24 h
- eICU window: ICU admission 0–1440 min; 24 h landmark
