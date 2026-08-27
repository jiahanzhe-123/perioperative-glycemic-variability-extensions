# Pre-commit Security Review

日期: 2026-08-01; 仓库: public repository (git init, branch main)

## 暂存概况

- staged files: 152
- config/paths.yml staged: 0 (必须=0)
- data/ staged (除 README+synthetic): 0 (必须=0)
- results/ staged (除 results/qc): 0 (必须=0)
- symlink staged: 0 (必须=0)
- 大文件 (>1MB) staged:

## 扫描

- The absolute-path regex appears once in its own test definition
  (`tests/pytest/test_release_consistency.py`); this is a scanner self-match,
  not a local path or credential in the repository.
- credential-literal scan on staged files: PASS (0 hits)
- patient-like columns occur only in the explicitly simulated files under
  `data/synthetic/`; no clinical records are present.

## 裁定

- sensitive_scan(含 staged 与全树): PASS
- 患者 ID 列命中仅限 data/synthetic/(模拟数据,README 明确标注,允许入库)。
- 无大文件、无 symlink、无 paths.yml、无 credentials、无患者级数据。
- 结论: 暂存内容可提交。
