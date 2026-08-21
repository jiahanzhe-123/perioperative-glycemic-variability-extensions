# Where data lives (local only)

`data/` holds **local** inputs and intermediate outputs. Nothing here is
committed to git except this README and `data/synthetic/`.

Expected layout when running with real (authorized) data:

```text
data/
├── raw/          # untouched database extracts (never modified)
├── interim/      # intermediate pipeline outputs
├── processed/    # analysis-ready frames (patient-level — never share)
└── synthetic/    # simulated test data (committed; safe to share)
```

**正式代码不依赖本目录下的任何 symlink 或默认数据位置。** 所有外部数据
(患者级帧、冻结输入、参考结果)一律通过本机 `config/paths.yml` 引用
(模板:`config/paths.example.yml`);该文件已 gitignore。2026-08-01 起,
原先指向旧工作区的 8 个 data/ symlink 已移除——请改用配置键
(如 `mimic_derived_data`、`mimic_record_work`、`bmi_repair_work`、
`inspire_work`)。

Real patient-level data from MIMIC-IV, eICU-CRD, or INSPIRE must never be
copied into a location that could be committed. See `DATA_AVAILABILITY.md`.
