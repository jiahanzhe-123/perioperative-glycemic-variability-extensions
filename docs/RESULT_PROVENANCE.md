# Result Provenance

**Analysis of record (2026-08-01 起生效):** MIMIC-IV 主分析应用 BMI-DQ-1 规则
(BMI<10 或 >80 kg/m² → 缺失 → 冻结 MICE)后的结果:**Model B 30 天 HR 0.979
(0.922–1.040),P=0.495;365 天 HR 0.985 (0.942–1.030),P=0.506**(N=10,561;296/745 事件)。
INSPIRE 主结果为 v5 统一行政删失分析:**I2 30 天 HR 0.905 (0.653–1.254),P=0.547**。

**不保留(not retained)的旧估计:** BMI 规则前的 0.977 (0.921–1.036) / 0.992
(0.949–1.036) 仅为历史存档,留存于 legacy 工作区,不得再引用为结果(与
Additional file 4 的 "not retained" 记录一致)。

Mapping from manuscript/supplement result blocks to the code that generates
them. All paths are stable repository-relative paths. 展项级(逐表/逐图)溯源见
`docs/FINAL_OUTPUT_PROVENANCE.csv`(Main Tables 1–5、Main Figures 1–5、
Tables S1–S25、Figures S1–S8,共 43 行)。

投稿附加文件采用当前编号:**AF1 = Supplementary Material;AF2 = Figure Source
Map;AF3 = Model Specification Matrix;AF4 = Analysis Provenance Summary**。
早期重建时代的 "Additional file 9/10" 编号已废弃,不得引用。

Machine-readable model specifications: `config/model_specifications.yml`
(14 model rows)。Manuscript-side matrix: AF3。

## Main text

| Result block | Database | Generating script(s) | Key inputs | N / events | Role |
|---|---|---|---|---|---|
| Cohort flow | MIMIC-IV 3.1 | `sql/mimic/02_build_cohort.sql`; `analyses/01_cohort_construction/02_build_analysis_dataset.py` | MIMIC extraction | 12,992 → 10,561 | primary cohort |
| **BMI 合理性规则(数据准备)** | MIMIC-IV | **`analyses/01_cohort_construction/03_apply_bmi_plausibility.py`** | `analysis_base.csv`(只读)→ `analysis_base_bmi_repaired.csv` | 35 行转缺失(28 例 >80、7 例 <10;12,992 行) | analysis-of-record 数据步骤;输出与冻结 checksum 2163168d… 逐字节一致 |
| Baseline table | MIMIC-IV | `analyses/09_quality_control/01_assemble_results.py` | BMI-repaired analysis frame | N=10,561 | descriptive |
| **Primary: Model B, 30-day** | MIMIC-IV | `analyses/03_primary_mimic/03_run_primary_models_mice.R` | MICE frames (m=50, 20 it, PMM, seed 20260726) | N=10,561, 296 deaths | **unique primary analysis** — HR 0.979 (0.922–1.040) per 10 mg/dL, P=0.495 |
| Key secondary: Model B, 365-day | MIMIC-IV | same | same | 745 deaths | key secondary — HR 0.985 (0.942–1.030), P=0.506(PH violated; prespecified descriptive average) |
| Models A/C | MIMIC-IV | same | same | — | context / measurement-process |
| PH + RCS + absolute risk | MIMIC-IV | `analyses/03_primary_mimic/04_ph_rcs_absoluterisk.R` | MICE frames | — | diagnostics / secondary;30d RD −0.04 pp |
| 365-day time-dependence(intervals 1–7 / 8–30 / 31–365;GV×log(time);30d/365d 协变量 log(time) 敏感性) | MIMIC-IV | `analyses/04_time_dependent/05_mice_pooled_interval_models.R`、`09_covariate_tt_sensitivity.R`(已迁入并重跑;与冻结值全等) | analysis frame | 162/139/444 interval deaths | secondary;intervals 1.063/0.974/0.938;interaction 0.966;31–365 HR<1 不得解读为保护 |
| Main figures | MIMIC-IV | `analyses/09_quality_control/figure_build/`(build_main_figures / build_fig23 / build_supplement / build_ph_audit / merge_s1) | model outputs | — | display |

## SHR–GV joint module (pre-specified secondary; MIMIC HbA1c subset)

| Result block | Generating script(s) | N / events | Key values (analysis of record) |
|---|---|---|---|
| Joint frame + IPW audit | `analyses/08_shr_component/01_shr_hba1c_ipw.R`, `02_joint_prep.py` | N=4,779 (nested) | HbA1c strictly 1–90 d pre-index |
| J1 vs J0 2-df joint tests (30 d / 365 d) | `03_joint_models.R` | 99 / 303 events | omnibus F=1.844, P=0.158 / F=2.338, P=0.097 |
| D2 vs D1 GV added-information (30 d / 365 d) | `03_joint_models.R` | — | P=0.948 / P=0.476 |
| Standardized risk + paired bootstrap performance | `04_joint_risk_performance.R` | — | internal performance only |
| Source sensitivity | `05_joint_source_sensitivity.R` | — | sensitivity(blood-gas 限制框 30d HR 1.212,P=0.032——预设敏感性,不改变主模块结论) |
| Joint QC | `06_joint_qc_assembly.py` | — | QC |
| SHR×GV interaction (exploratory) | `03_joint_models.R` | — | P=0.105 / 0.107 |

Interpretation guard: the SHR ratio coefficient is never interpreted as an
independent biological effect; the high-high corner is not labelled a
high-risk phenotype.

## INSPIRE module (secondary exact-timestamp cohort analysis; not formal external validation)

**Analysis of record = v5 统一行政删失分析**(`analyses/06_inspire/05_uniform_admin_censoring_v5.R`,
2026-08-01 由仓库代码重跑,输出与 v5 包逐字节一致):死亡者用复合死亡时间,
所有非事件在统一术后第 30 日行政终点删失,绝不用出院删失。

| Result block | Generating script(s) | N / events | Key values |
|---|---|---|---|
| Cohort + reconciled outcome | `sql/inspire/06_cohort_build.sql`, `analyses/06_inspire/02_outcome_time_repairs.R` | eligible N=1,355; landmark N=1,353 | reconciled 30-day mortality |
| **Final models (30-day, v5)** | **`05_uniform_admin_censoring_v5.R`** | N=1,353, 27 events | I2 HR 0.905 (0.653–1.254), P=0.547 |
| 48-h exposure → 48-h landmark (v5) | same | N=1,511, 31 events | HR 1.104 (0.846–1.441), P=0.464 |
| 30-day standardized risk (v5) | same | N=1,353 | RD −0.33 pp (−1.80 to +1.25) |
| 365-day | — | — | **WITHDRAWN(v5)**:无统一 365 日登记覆盖终点;不报告 HR/绝对风险/RMST/校准 |
| QC assembly | `05_qc_assembly.py` | — | QC |

## eICU-CRD module (multicenter harmonized comparison — NOT external validation)

| Result block | Generating script(s) | N / events | Key values |
|---|---|---|---|
| Phenotype + extraction | `sql/eicu/01…09` in order | — | high-specificity CABG/open-valve/combined phenotype |
| Harmonized M1–M4 | `analyses/07_eicu/01_fit_harmonized_models.R` | N=7,115, 130 events, 67 hospitals | M3 RR 1.158 (1.025–1.308), P=0.019; M4 RR 1.133 (0.990–1.296) |
| Random-intercept sensitivity | same | same | OR 1.167 (1.036–1.316); convergence warning documented |
| MIMIC-side harmonization | `02_harmonized_mimic_features.py`, `03_harmonized_mimic_model_data.py`, `04_models.R` | N=8,117, 128 events | M3 1.024 / M4 1.047 |
| Tables/figures + QC | `05_tables_figures.R`, `06_qc.R` | — | — |

## Supplement refits (frame-fixed series)

| Block | Script |
|---|---|
| Primary refits (SHR-GV frames) | `analyses/03_primary_mimic/supplement/01_refit_primary_models.R` |
| Alternative GV metrics | `supplement/03_alt_gv_framefix.R` |
| Measurement process | `supplement/03_measurement_framefix.R` |
| Diabetes interaction | `supplement/04_diabetes_ix_framefix.R` |
| Extreme-glucose burden | `supplement/05_extremes_framefix.R` |
| Fixed-scale display constants | `supplement/06_fixed_scale_constants.py`, `07_fixed_scale_rescale.py`, `08_fixed_scale_display.py` |

## Verification status

- `make validate` 对 32 个哨兵值全部 PASS(锚定 = analysis of record;
  见 `results/qc/frozen_validation_report.csv`)。
- MIMIC 主分析与 INSPIRE v5 已于 2026-08-01 由仓库代码在本机授权环境重跑,
  输出与 analysis-of-record 逐字节一致(见 `results/qc/final_repository_qc.md`)。
- 展项级溯源表:`docs/FINAL_OUTPUT_PROVENANCE.csv`;figure_build 绘图脚本与
  365 天 MICE 池化区间脚本的迁入为开放事项(`docs/OPEN_ISSUES.md`)。
