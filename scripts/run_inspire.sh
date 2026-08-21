#!/usr/bin/env bash
# INSPIRE 分析入口(需要机构授权数据与数据库容器;缺时明确报错退出非零)
set -euo pipefail
cd "$(dirname "$0")/.."
python3 - <<'PY'
import os, sys
sys.path.insert(0, "src/common")
from paths import PGV
p = os.path.join(PGV("inspire_work"), "data", "inspire_base.csv")
if not os.path.exists(p):
    print("[pgv] 缺少授权 INSPIRE 输入(Level 3):", p)
    print("INSPIRE 数据访问取决于数据库管理方与机构治理;请在 config/paths.yml 配置后重试。")
    sys.exit(2)
print("[pgv] INSPIRE 输入齐备,依次执行(v3 frame repairs followed by v5 analysis-of-record)…")
PY
Rscript analyses/06_inspire/02_outcome_time_repairs.R
Rscript analyses/06_inspire/03_sensitivity_joint_repairs.R
# v3 repairs establish the corrected frame; v5 is the analysis-of-record
# script that produces the manuscript's reported INSPIRE estimates.
Rscript analyses/06_inspire/04_final_models_v3.R
Rscript analyses/06_inspire/05_uniform_admin_censoring_v5.R
python3 analyses/06_inspire/05_qc_assembly.py
