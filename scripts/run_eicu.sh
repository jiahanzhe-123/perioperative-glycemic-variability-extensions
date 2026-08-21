#!/usr/bin/env bash
# eICU harmonized 模块入口(需要授权数据)
set -euo pipefail
cd "$(dirname "$0")/.."
python3 - <<'PY'
import os, sys
sys.path.insert(0, "src/common")
from paths import PGV
cands = [
    os.path.join(PGV("replication_work"), "03_eicu_harmonized", "eicu_model_data.csv"),
]
missing = [p for p in cands if not os.path.exists(p)]
if missing:
    print("[pgv] 缺少授权 eICU 输入(Level 2):")
    for m in missing: print("  -", m)
    print("eICU-CRD 需 PhysioNet 授权;配置 config/paths.yml 后重试。")
    sys.exit(2)
print("[pgv] eICU 输入齐备,执行 harmonized 模块…")
PY
Rscript analyses/07_eicu/01_fit_harmonized_models.R
