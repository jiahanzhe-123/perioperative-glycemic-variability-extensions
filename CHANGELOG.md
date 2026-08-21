# Changelog

All notable changes to this repository are documented here. Scientific results
are frozen; entries below concern code organization and release preparation only.

## [1.1.0] — 2026-08-22 (JAHA v5 Phase 3B extension release)

- Added the final-lineage Phase 3B aggregate evidence bundle, including source
  agreement, log-scale sensitivity, source-count closure, source-defined
  coefficient summaries, and the 18-row SD-GV analytic-context landscape.
- Added a public aggregate Figure 3 builder that reads only reviewed aggregate
  rows and does not refit models or access patient/stay-level values.
- Added aggregate-safe Figure 1 and Figure 3 files, panel-level source maps,
  and figure provenance. The controlled Figure 2 patient-level scatter and
  Bland–Altman source values remain excluded from the public repository.
- Updated release metadata to repository version `v1.1.0` and extension package
  citation version `0.2.0`. The GitHub release is published; no new archival DOI
  is claimed for v1.1.0 until a versioned archive record is independently verified.

## [1.0.1] — 2026-08-04 (Zenodo archival release)

- Added the Zenodo concept DOI [`10.5281/zenodo.21791846`](https://doi.org/10.5281/zenodo.21791846)
  and version DOI [`10.5281/zenodo.21791847`](https://doi.org/10.5281/zenodo.21791847)
  to the citation trail for the public code release.

## [1.0.0] — 2026-08-04 (public code release)

- Released the verified analysis code under the MIT License with the confirmed
  Department of Cardiology, The Second Qilu Hospital, Cheeloo College of
  Medicine, Shandong University copyright holder.
- Frozen validation remains **35/35 PASS**; the repository contains no
  patient-level data or credentials. The local INSPIRE PostgreSQL credential
  was rotated and is stored only in ignored local configuration.
- Published the versioned `v1.0.0` release independently of journal submission.

## [1.0.0-rc3] — 2026-08-01 (AI 收尾;未公开)

### Added
- 迁入最终绘图管线:`analyses/09_quality_control/figure_build/`(5 个脚本,
  覆盖 Main F1–F5 与 Supplementary S1–S8;剥离 daimon_runtime 依赖,配置驱动)。
- 迁入 365 天 MICE 池化分析:`analyses/04_time_dependent/05_mice_pooled_interval_models.R`
  (区间 1.063/0.974/0.938、GV×log(time) 0.966、30d tt 校正;sanity 复现主分析)
  与 `09_covariate_tt_sensitivity.R`(协变量 log(time) D1/D2)。
- 迁入 eICU 随机截距重跑:`analyses/07_eicu/07_random_intercept_rerun.R`
  (OR 1.167435;字节一致)。
- 新增发布一致性测试 `tests/pytest/test_release_consistency.py`
  (展项溯源 43/43、AF2 图源、AF 编号、数字一致性、AoR 集成;23/23 通过)。
- git 初始化(main 分支,12 个 commits；其中 10 个为项目分析/QC 主题提交，
  另 2 个为私有评审状态同步);提交前安全审查
  `results/qc/precommit_security_review.md` PASS。

### Changed
- 全部可运行模块在授权本机重跑:BMI 帧、MIMIC 主 MICE、PH/RCS/绝对风险、
  来源敏感性、SHR 联合模型、365d 区间/协变量 tt、INSPIRE v5、eICU 协调+RI——
  输出与 analysis of record 逐字节一致(flexsurv 诊断除外:本环境可执行
  30d RP 对照,不改变任何报告值)。
- 冻结锚定增至 35 项并全部指向仓库可复现位置(mimic_record_work/
  inspire_record_work);validator 35/35 PASS。
- AF2 重写为 13 个图展项(仓库相对脚本路径 + 验证状态);AF3 展开 12 处
  含糊规格;AF4 新增 INSPIRE v3 not-retained 行并标注当前值为 v5。
- 投稿包技术 QC 完成(PGV_SUBMISSION_RC3):主文修复 1 处杂散 `}`,
  元数据日期重置;Supplement 采用分页修复版并验证 S1–S25/S1–S8 连续。

### Notes
- GitHub 私有仓库已创建并推送；当前 review commit `7f5c811` 的
  documentation、lint、tests + synthetic 三个 Actions workflow 均通过。
  尚未创建版本标签、发布公开仓库或 Zenodo DOI。发布门槛见
  `docs/RELEASE_CHECKLIST.md`。

## [1.0.0-rc2] — 2026-08-01 (release-blocker fixes; not yet public)

### Analysis-of-record alignment (numbers re-anchored, no science changed)
- **BMI-DQ-1 规则接入正式代码路径**:新增
  `analyses/01_cohort_construction/03_apply_bmi_plausibility.py`
  (BMI<10 或 >80 → 缺失 → 冻结 MICE;35 行,输出与冻结帧 sha256 逐字节一致)。
- 主分析 `03_run_primary_models_mice.R` 及全部 MIMIC 下游脚本改为读取
  `analysis_base_bmi_repaired.csv`,路径全部改由 `config/paths.yml` 驱动。
- **重跑验证(授权本机)**:MIMIC 主分析输出与 analysis of record 逐字节一致
  (Model B 30d HR 0.9794 (0.9225–1.0398) P=0.4952;365d 0.9850
  (0.9420–1.0299) P=0.5062;N=10,561;296/745)。
- **INSPIRE v5 迁移**:新增 `analyses/06_inspire/05_uniform_admin_censoring_v5.R`
  (统一术后 30 日行政删失,绝不用于出院删失);重跑输出与 v5 包逐字节一致
  (I2 30d HR 0.9046 P=0.5471;48h 1.1044 P=0.4643)。
- 冻结锚定更新:`tests/expected/frozen_results.yml` 现锚定 analysis of record
  (0.9794/0.9850;SHR omnibus F=1.844 P=0.158;绝对风险 RD −0.00044;
  INSPIRE v5)。**0.977/0.992 为 BMI 规则前旧估计,不保留,仅存档。**
- INSPIRE 365d 哨兵随稿件撤回(v5)从冻结验证移除;旧值仅存于 legacy 工作区。
- validator 现 35/35 PASS。

### Infrastructure
- 新增 `config/paths.example.yml`;本地 `config/paths.yml`(gitignored)。
  `src/common/paths.R|py` 收紧:外部数据键必须显式配置,不再回退仓库内目录。
- 删除正式代码中全部本机工作区硬编码路径(约 25 个文件)。
- 新增管道级集成测试 `tests/pytest/test_bmi_rule_integration.py`
  (验证主分析输入帧确实执行 BMI 规则;17/17 passed)。
- 新增 `docs/FINAL_OUTPUT_PROVENANCE.csv`(43 展项溯源)与 AF1–AF4 编号统一。
- 仓库内 8 个指向旧工作区的 data/ symlink 已移除(配置驱动替代)。

### Known remaining gaps (see docs/OPEN_ISSUES.md)
- figure_build 绘图脚本与 365 天 MICE 池化区间脚本已迁入并完成重跑核对。
- GitHub 私有推送/CI、作者与治理元数据、数据库密码轮换仍需作者执行。

## [1.0.0-rc] — 2026-08-01 (release candidate, not yet public)

### Added
- Public repository structure assembled from the internal analysis workspaces:
  `sql/`, `analyses/`, `src/common/`, `config/`, `tests/`, `scripts/`, `docs/`.
- Central configuration system (`config/paths.yml` + `src/common/paths.R|py`);
  all hard-coded local paths and credentials removed from analysis code.
- Frozen-result validation: `tests/expected/frozen_results.yml` +
  `scripts/validate_results.py` (`make validate`).
- Synthetic dataset generator and end-to-end synthetic workflow
  (`data/synthetic/`, `analyses/09_quality_control/99_synthetic_workflow.R`).
- Unit tests (pytest + testthat) for core data rules.
- Machine-readable model specifications (`config/model_specifications.yml`).
- Governance files: LICENSE (MIT with confirmed copyright holder), CITATION.cff,
  CONTRIBUTING, CODE_OF_CONDUCT, SECURITY, DATA_AVAILABILITY, CHANGELOG.
- GitHub Actions: lint, tests + synthetic workflow, documentation checks.

### Changed
- File names normalized to stable, execution-ordered names
  (see `docs/FILE_MIGRATION_MAP.csv`).

### Archived (not deleted)
- Superseded internal script versions moved to
  `archive/not_used_in_final_analysis/` (see `archive/README.md` and
  `docs/ARCHIVE_AND_REMOVAL_LOG.md`).

### Scientific content
- No analysis logic, cohort definition, model, or result was modified.
