# CODE INVENTORY(2026-07-29,整理前审计)

范围:围术期血糖变异性研究(MIMIC-IV / INSPIRE / eICU-CRD)。
清点规则:不按文件名(final/fixed/vN)判断版本,按调用关系、输出内容与投稿冻结数字确认。
状态列:USED-FINAL(被最终结果使用)/ SUPERSEDED(已被替代,留档)/ ARCHIVE(无最终结果引用,归档)。

## A. MIMIC-IV 主分析链(数据库:MIMIC-IV v3.1)

### A1. 数据提取与队列(SQL,原始项目 `~/Documents/Codex/2026-07-18/volumes-mimic-data/scripts/`)

| 原路径 | 功能 | 状态 | 建议新路径 |
|---|---|---|---|
| scripts/00_schema_inventory.sql | schema 盘点 | USED-FINAL(血缘) | sql/mimic/00_schema_inventory.sql |
| scripts/01_build_cardiac_surgery_codebook.sql | ICD 代码本(6,523 码人工筛查) | USED-FINAL | sql/mimic/01_build_cardiac_surgery_codebook.sql |
| scripts/02_build_cohort.sql | 队列构建(codebook 优先级+stay 链接) | USED-FINAL | sql/mimic/02_build_cohort.sql |
| scripts/03_build_itemid_review.sql | lab itemid 审查 | USED-FINAL | sql/mimic/03_build_itemid_review.sql |
| scripts/04_extract_vitals_24h.sql | 生命体征 | USED-FINAL | sql/mimic/04_extract_vitals_24h.sql |
| scripts/05_extract_labs.sql | 实验室(glucose_long_v2 等) | USED-FINAL | sql/mimic/05_extract_labs.sql |
| scripts/06_extract_comorbidities_bmi.sql | 合并症/BMI | USED-FINAL | sql/mimic/06_extract_comorbidities_bmi.sql |
| scripts/07_extract_treatments.sql | 治疗 | USED-FINAL | sql/mimic/07_extract_treatments.sql |
| scripts/08_build_severity.sql | SOFA | USED-FINAL | sql/mimic/08_build_severity.sql |
| scripts/09_build_shr_gv.sql | SHR–GV 特征 | USED-FINAL | sql/mimic/09_build_shr_gv.sql |
| scripts/10_build_final_features.sql | 终版特征 | USED-FINAL | sql/mimic/10_build_final_features.sql |
| scripts/11_build_qc_report.sql | QC 报告 | USED-FINAL | sql/mimic/11_build_qc_report.sql |
| (github_release)sql/mimic/12_bloodonly_glucose_patch.sql | 尿糖排除补丁(283 患者) | USED-FINAL | sql/mimic/12_bloodonly_glucose_patch.sql |

### A2. 当前主分析(rebuild 项目 `~/cardiac_glucose_rebuild_20260728/scripts/`)—— 最终版

| 原路径 | 功能 | 状态 | 建议新路径 |
|---|---|---|---|
| scripts/01_glucose_series_rebuild.py | 血糖序列三套(priority/allmedian/poctfirst)+ 去重审计 | USED-FINAL | analyses/02_glucose_processing/01_rebuild_glucose_series.py |
| scripts/02_build_analysis_dataset.py | 分析数据集(队列/landmark/手术类别/charlson_wo_dm/HbA1c 层级) | USED-FINAL | analyses/01_cohort_construction/02_build_analysis_dataset.py |
| scripts/03_primary_models.R | 主模型 CC + MICE(初版,loggedEvents 诊断) | SUPERSEDED(由 03b 完成 MICE) | archive/not_used_in_final_analysis/ |
| scripts/03b_mice_rerun.R | MICE m=50 修复版(0 loggedEvents)→ 主要结果 | USED-FINAL | analyses/03_primary_mimic/03_run_primary_models_mice.R |
| scripts/04_ph_rcs_absrisk.R | PH/RCS/RP/绝对风险 bootstrap(1000) | USED-FINAL | analyses/03_primary_mimic/04_ph_rcs_absoluterisk.R |
| scripts/06_shr_hba1c.R | SHR 模块 + IPW | USED-FINAL | analyses/08_shr_component/01_shr_hba1c_ipw.R |
| scripts/07a_source_series.py | 来源序列(POCT/lab/bloodgas/common/同患者) | USED-FINAL | analyses/05_source_sensitivity/01_build_source_series.py |
| scripts/07b_source_models.R | 来源敏感性模型 | USED-FINAL | analyses/05_source_sensitivity/02_fit_source_models.R |
| scripts/08_harmonized.R | harmonized MIMIC–eICU modified Poisson | USED-FINAL | analyses/07_eicu/01_fit_harmonized_models.R |
| scripts/09_sensitivity.R | open-core/替代 GV/极端/winsorized/expanded/minimal | USED-FINAL | analyses/03_primary_mimic/05_sensitivity_models.R |
| scripts/10_performance.R | 表现与增量(成对 bootstrap) | USED-FINAL | analyses/03_primary_mimic/06_model_performance.R |
| scripts/11_assembly.py | spec curve/矩阵/QC/number map | USED-FINAL | analyses/09_quality_control/01_assemble_results.py |
| scripts/12_missing_data_report.R | 缺失诊断 | USED-FINAL | analyses/09_quality_control/02_missing_data_report.R |

### A3. SHR–GV 联合模块(rebuild scripts/20–24)—— 最终版

| 原路径 | 功能 | 状态 | 建议新路径 |
|---|---|---|---|
| scripts/20_joint_prep.py | 联合队列/常数/分布/共线性/支持域 | USED-FINAL | analyses/08_shr_component/02_joint_prep.py |
| scripts/21_joint_models.R | J0/J1/N/D/I 模型 + D1 池化检验 + PH + timing + IPW | USED-FINAL | analyses/08_shr_component/03_joint_models.R |
| scripts/22_joint_risk_perf.R | 联合绝对风险+性能+稳定性 | USED-FINAL | analyses/08_shr_component/04_joint_risk_performance.R |
| scripts/23_joint_source.R | 联合来源敏感性 | USED-FINAL | analyses/08_shr_component/05_joint_source_sensitivity.R |
| scripts/24_joint_qc_assembly.py | 联合 QC + 矩阵 + number map | USED-FINAL | analyses/08_shr_component/06_joint_qc_assembly.py |

### A4. 历史/审计(非最终,superceded)

- `~/cardiac_glucose_supplementary_analyses/R/01–09`(landmark/alt_gv/measurement/diabetes_ix/extremes/ipw/missing/absrisk/figures):**SUPERSEDED**(frame 缺陷由 framefix 版替代)→ archive
- `~/cardiac_glucose_supplementary_analyses/corrected_scripts/*_bloodonly.R`(非 framefix 版):**SUPERSEDED** → archive
- `corrected_scripts/*_framefix.R` 与 `06_rcs_spline_framefix.R`、`08_figures_framefix.R`、`08b_figures_fixedscale.R`、`09_rcs_curves_framefix.R`、`07_frame_registry.R`:frame 修复版(补充材料来源)→ USED-FINAL(补充链)→ analyses/04_time_dependent 与 analyses/03_primary_mimic 的 supplement 部分
- `reviewer_issue_verification_20260728/scripts/fixed_scale_*.py`:fixed-scale 统一(表内 per-SD 尺度)→ USED-FINAL(补充链)
- `github_release/`(早期公开草稿):**SUPERSEDED**(由本次新仓库替代)

## B. INSPIRE(数据库:INSPIRE v1.4.2)

| 原路径 | 功能 | 状态 | 建议新路径 |
|---|---|---|---|
| inspire_cardiac_20260729/scripts/01_import_tables.sh | 导入(6 表,计数一致) | USED-FINAL | analyses/06_inspire/01_import_tables.sh |
| sql/00_create_database.sql | schemas+DDL | USED-FINAL | sql/inspire/00_create_database.sql |
| sql/02_indexes.sql | 索引 | USED-FINAL | sql/inspire/02_indexes.sql |
| sql/03_schema_audit.sql | 结构审计 | USED-FINAL | sql/inspire/03_schema_audit.sql |
| sql/04_procedure_audit.sql | 手术/时间审计 | USED-FINAL | sql/inspire/04_procedure_audit.sql |
| sql/05_lab_audit.sql | 实验室字典 | USED-FINAL | sql/inspire/05_lab_audit.sql |
| sql/06_cohort_build.sql | 表型+队列 | USED-FINAL | sql/inspire/06_cohort_build.sql |
| sql/07_glucose_features.sql | 血糖长表+GV/SHR | USED-FINAL | sql/inspire/07_glucose_features.sql |
| sql/08_analysis_base.sql | analysis_base | USED-FINAL | sql/inspire/08_analysis_base.sql |
| sql/09_comorbidity_outcome_audit.sql | 合并症+结局结构审计 | USED-FINAL | sql/inspire/09_comorbidity_outcome_audit.sql |
| sql/10_anchor_windows.sql | 多锚点血糖窗 | USED-FINAL | sql/inspire/10_anchor_windows.sql |
| sql/qc_tests.sql | QC | USED-FINAL | sql/inspire/11_qc_tests.sql |
| scripts/11_inspire_analysis.R | 主分析 I1–I3+MICE(旧结局,SUPERSEDED) | SUPERSEDED(结局修复) | archive |
| scripts/12_inspire_sensitivity.R | 锚点/分辨率(旧结局,SUPERSEDED) | SUPERSEDED | archive |
| scripts/13_inspire_joint_perf.R | 联合+性能(旧结局,SUPERSEDED) | SUPERSEDED | archive |
| scripts/14_inspire_qc_assembly.py | QC+矩阵(旧) | USED-FINAL(结构沿用) | analyses/06_inspire/05_qc_assembly.py |
| scripts/15_outcome_time_repairs_models.R | 结局复合+48h 正确 landmark | USED-FINAL | analyses/06_inspire/02_outcome_time_repairs.R |
| scripts/16_repairs_remaining_reruns.R | 修正结局下敏感性/联合 | USED-FINAL | analyses/06_inspire/03_sensitivity_joint_repairs.R |
| scripts/17_ccmask_fix_rerun.R | CC mask 修复(v3,最终 INSPIRE 结果) | USED-FINAL | analyses/06_inspire/04_final_models_v3.R |

## C. eICU-CRD(数据库:eICU-CRD v2.0)

| 原路径 | 功能 | 状态 | 建议新路径 |
|---|---|---|---|
| eicu_cardiac_glucose_validation/sql/01–08 | schema 审计/术语提取/分类/证据/队列/血糖/协变量/结局/QC/终表 | USED-FINAL | sql/eicu/01–08 |
| cardiac_glucose_multidatabase_replication_v1/sql/03_eicu_harmonized.sql | harmonized 数据 | USED-FINAL | sql/eicu/09_harmonized.sql |
| replication_v1/python/02_mimic_features.py、03_mimic_model_data.py | MIMIC 侧 harmonized 特征 | USED-FINAL | analyses/07_eicu/02_harmonized_mimic_data.py |
| replication_v1/R/05_models.R、06_tables_figures.R、04_qc.R | harmonized 模型/表图/QC | USED-FINAL | analyses/07_eicu/03_models_figures_qc.R |
| replication_v1/python/07_make_manuscript.py、08_make_supplement_docx.py、09_audit.py | 稿件生成(内部工具,含本地路径) | ARCHIVE(内部,不公开) | archive/not_used_in_final_analysis/internal_only |

## D. 无法追溯/待确认项

- 主文/补充 DOCX 的最终整合脚本(reviewer_fix/fix_*.py):内部 Word 工具,含本机路径 → 不公开(archive/internal_only);
- 早期 `cardiac_glucose_final_*` 与 `volumes-eicu-data` 等中间项目:中间态审计材料 → archive;
- 所有患者级数据文件(data/、outputs/ 大部分、00_glucose_series/、各 *dataset*.csv):**禁止公开**(.gitignore);
- 主文/Supplement DOCX 本体:投稿包,不放代码仓库(按作者决定)。

## E. 硬编码与敏感内容(整理时必须清除)

- 全部脚本含本地用户目录绝对路径(已统一净化为 `~/...` 记录) → 改读 config;
- docker 容器名/端口/密码(inspire-pg、mimic-pg、eicu-pg)→ config,密码不入库;
- 患者级 stay_id/subject_id 列(仅存在于禁止公开的数据文件,代码中不出现字面 ID);
- 数据库连接字符串(仅 inspire 脚本内嵌 docker exec;改为 config 驱动)。
