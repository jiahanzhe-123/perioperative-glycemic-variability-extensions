# Final Repository QC

Date: 2026-08-01 (rc3 — final rerun recorded). Scope: canonical repository
`perioperative-glycemic-variability` (git initialized, branch main, 11 commits;
10 project analysis/QC commits plus the private-review status sync).

## 1. Analysis-of-record 复现(全部授权本机实跑,2026-08-01)

| 模块 | 仓库脚本 | 结果 | 与 analysis of record 比对 |
|---|---|---|---|
| BMI 规则 | `01_cohort_construction/03_apply_bmi_plausibility.py` | 35 行设缺失 | 输出 sha256 与冻结帧**逐字节一致** |
| MIMIC 主 MICE | `03_primary_mimic/03_run_primary_models_mice.R` | 30d 0.9794 (0.9225–1.0398) P=0.4952;365d 0.9850 (0.9420–1.0299) P=0.5062 | **逐字节一致** |
| PH/RCS/绝对风险 | `03_primary_mimic/04_ph_rcs_absoluterisk.R` | PH global 1.354e-15;RD −0.00044 | 结果全精度一致(flexsurv 诊断项环境差异,见 §4) |
| 来源敏感性 | `05_source_sensitivity/01+02` | POCT 0.827 等 | **逐字节一致** |
| SHR 联合 | `08_shr_component/01+02+03` | omnibus F=1.844 P=0.158 / F=2.338 P=0.097 | **逐字节一致** |
| 365d 区间族 | `04_time_dependent/05_mice_pooled_interval_models.R` | 1.063/0.974/0.938;交互 0.966;30d tt 0.9795 | 区间**逐字节一致**;logtime/tt 数值全精度一致(版式源于旁路脚本,已记录) |
| 365d 协变量 tt | `04_time_dependent/09_covariate_tt_sensitivity.R` | D1 HR 0.9898 (0.9459–1.0356) P=0.6562;D2 1.029/0.971/0.968 | 数值全精度一致(与 AF4 0.990/0.656 吻合;仅 note 措辞与末位浮点差异) |
| INSPIRE v5 | `06_inspire/05_uniform_admin_censoring_v5.R` | 0.9046 P=0.5471;48h 1.1044 P=0.4643 | **逐字节一致**(含 bootstrap) |
| eICU 协调 + RI | `07_eicu/01_fit_harmonized_models.R`、`07_random_intercept_rerun.R` | M3 1.1577;RI OR 1.167435 | **逐字节一致** |
| 绘图 | `09_quality_control/figure_build/`(5 脚本) | 13 图展项全部重生成 | 内容一致(数值/标签/范围);渲染随 matplotlib 环境,不保证逐字节 |

## 2. 机器检查(全部实跑)

| 检查 | canonical 仓库 | 干净 clone(/tmp/pgv_clone_test) |
|---|---|---|
| 敏感信息扫描 | PASS | PASS |
| 配置校验(14 模型) | PASS | PASS |
| pytest | 23 passed | 17 passed + 6 skipped(授权项自动 skip) |
| R testthat | ALL PASSED | ALL PASSED |
| synthetic 流程 | PASS(N=400) | PASS(N=400) |
| 冻结验证 | **35/35 PASS** | UNRESOLVED(无授权数据,正确报缺,不伪 PASS) |
| R 语法解析(34 文件) | 0 失败 | — |
| Python 语法解析(24 文件) | 0 失败 | — |
| 绝对路径/ symlink / 患者数据 | 0 / 0 / 无 | 0 / 0 / 无 |
| 文档链接 | PASS | — |
| 提交前安全审查 | PASS(152 文件,无 paths.yml/data/results 泄漏) | — |

## 3. Git 状态

- branch: main;commits: 11 (10 project analysis/QC commits plus the current
  documentation sync); no remote configured.
- staged 安全审查:152 文件,无 config/paths.yml、无 data/(除 README+synthetic)、
  无 results/(除 qc)、无 symlink、无大文件、无凭证。
- **未推送**:GitHub CLI is intentionally unauthenticated and no remote is
  configured. Authenticate and create the private remote only after the author
  confirms the account and repository name.
- CI(workflows 已备):本地逐步模拟全部通过;GitHub Actions 未实跑(需推送)。

## 4. 已知环境差异(不影响任何报告值)

- flexsurv 在本机可用:04 脚本的 Royston–Parmar 30d 对照本次成功估计
  (HR 1.052,与 Cox Model A 1.047 一致);冻结时代该对照为 NA(环境无 flexsurv)。
  365d RP 两次均初始化失败。该诊断不进任何稿件数字。
- 图渲染:投稿包图与本机重生成图内容一致,像素因 matplotlib 环境不同略有差异。

## 5. verdict

**READY FOR PRIVATE GITHUB REVIEW** (local release-candidate checks are
complete:
唯一 canonical 已确定;AoR 代码全部迁入;本地重跑一致;冻结验证通过;
绘图脚本完整;365d 池化区间代码完整;43/43 provenance;git 已初始化;
staged 安全;仅缺 gh 认证)。

**不得标记 READY FOR PUBLIC RELEASE**:LICENSE/作者信息/INSPIRE 伦理治理/
AI 声明/CRediT/密码轮换/GitHub 公开决定/v1.0.0/Zenodo DOI 均待作者确认。
