#!/usr/bin/env bash
# MIMIC 主分析入口(需要授权数据;缺数据时明确报错退出非零)
set -euo pipefail
cd "$(dirname "$0")/.."
python3 - <<'PY'
import os, sys
sys.path.insert(0, "src/common")
from paths import PGV
need = [
    ("mimic_freeze", "final_dataset_for_v3_pipeline.csv"),
    ("rebuild_work", "data/analysis_base.csv"),
]
missing = []
for key, rel in need:
    p = os.path.join(PGV(key), rel)
    if not os.path.exists(p):
        missing.append(p)
if missing:
    print("[pgv] 缺少授权 MIMIC 输入文件(Level 2):")
    for m in missing: print("  -", m)
    print("请在 config/paths.yml 配置授权数据位置后重试;或运行 make synthetic。")
    sys.exit(2)
print("[pgv] MIMIC 输入齐备,依次执行主分析…")
PY
python3 analyses/02_glucose_processing/01_rebuild_glucose_series.py
python3 analyses/01_cohort_construction/02_build_analysis_dataset.py
Rscript analyses/03_primary_mimic/03_run_primary_models_mice.R
Rscript analyses/03_primary_mimic/04_ph_rcs_absoluterisk.R
python3 analyses/09_quality_control/01_assemble_results.py
