#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_bmi_rule_integration.py — BMI-DQ-1 规则的管道级集成测试。

与孤立函数测试不同,本文件验证:
1. 规则逻辑本身(fixture 帧);
2. 正式数据准备脚本 analyses/01_cohort_construction/03_apply_bmi_plausibility.py
   在真实授权输入上产生的帧确实应用了规则(计数/checksum/逐列一致);
3. 主分析脚本 03_run_primary_models_mice.R 读取的是修复后帧而非原始帧;
4. 主分析输出(已重跑)与 analysis of record 一致。

需要授权数据的用例在无 config/paths.yml 配置时自动 skip(CI 安全)。
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src", "common"))
from paths import PGV  # noqa: E402

FROZEN_SHA256 = "2163168d1c2cf07828dd10d3b79392ca94545fbf36e7fd68dcd9d2675eeefba4"
FROZEN_CHANGED = 35
EXPECTED_ROWS = 12992
BMI_SCRIPT = os.path.join(ROOT, "analyses", "01_cohort_construction", "03_apply_bmi_plausibility.py")
PRIMARY_SCRIPT = os.path.join(ROOT, "analyses", "03_primary_mimic", "03_run_primary_models_mice.R")


def apply_rule(df, lo=10.0, hi=80.0):
    """规则 BMI-DQ-1 的参考实现(与正式脚本同一逻辑)。"""
    out = df.copy()
    bmi = out["bmi"]
    bad = bmi.notna() & ((bmi < lo) | (bmi > hi) | ~np.isfinite(bmi))
    out.loc[bad, "bmi"] = np.nan
    return out, int(bad.sum())


# ---------- 1. 规则逻辑(fixture,任何环境可跑) ----------
def test_rule_fixture_boundaries():
    fx = pd.DataFrame({
        "stay_id": [1, 2, 3, 4, 5, 6, 7, 8],
        "bmi": [9.99, 10.0, 80.0, 80.01, np.inf, 0.4, 5565.2, 27.3],
        "age": [50, 60, 70, 55, 65, 75, 45, 58],
    })
    out, n_bad = apply_rule(fx)
    assert n_bad == 5                       # 9.99, 80.01, inf, 0.4, 5565.2
    assert out.loc[out["stay_id"] == 1, "bmi"].isna().all()
    assert out.loc[out["stay_id"] == 2, "bmi"].iloc[0] == 10.0   # 边界保留
    assert out.loc[out["stay_id"] == 3, "bmi"].iloc[0] == 80.0   # 边界保留
    assert out.loc[out["stay_id"] == 8, "bmi"].iloc[0] == 27.3
    assert out["age"].equals(fx["age"])     # 非 bmi 列不受影响


def test_rule_keeps_original_missing():
    fx = pd.DataFrame({"stay_id": [1, 2], "bmi": [np.nan, 200.0]})
    out, n_bad = apply_rule(fx)
    assert n_bad == 1
    assert out["bmi"].isna().sum() == 2


# ---------- 2. 真实管道集成(需要授权数据) ----------
def _authorized():
    try:
        src = os.path.join(PGV("mimic_derived_data"), "data", "analysis_base.csv")
        return os.path.exists(src)
    except KeyError:
        return False


requires_data = pytest.mark.skipif(not _authorized(), reason="需要授权 MIMIC 派生数据(config/paths.yml)")


@requires_data
def test_pipeline_frame_applies_rule(tmp_path=None):
    """运行正式脚本并验证输出帧:零不可能 BMI、35 行转缺失、checksum 命中冻结值。"""
    r = subprocess.run([sys.executable, BMI_SCRIPT], capture_output=True, text=True)
    assert r.returncode == 0, f"BMI script failed:\n{r.stdout}\n{r.stderr}"
    import hashlib
    def sha256(p):
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    src = os.path.join(PGV("mimic_derived_data"), "data", "analysis_base.csv")
    dst = os.path.join(PGV("mimic_record_work"), "data", "analysis_base_bmi_repaired.csv")
    assert os.path.exists(dst)
    a = pd.read_csv(src, low_memory=False)
    b = pd.read_csv(dst, low_memory=False)
    assert len(b) == EXPECTED_ROWS and b["stay_id"].is_unique
    plausible = b["bmi"].dropna()
    assert ((plausible >= 10) & (plausible <= 80)).all(), "输出帧仍含不可能 BMI"
    assert int(b["bmi"].isna().sum()) == int(a["bmi"].isna().sum()) + FROZEN_CHANGED
    assert a.drop(columns=["bmi"]).equals(b.drop(columns=["bmi"]))
    assert sha256(dst) == FROZEN_SHA256, "输出帧与 analysis-of-record 冻结帧不一致"


@requires_data
def test_primary_script_reads_repaired_frame():
    """主分析脚本必须读取 analysis_base_bmi_repaired.csv,不得读取原始帧。"""
    txt = open(PRIMARY_SCRIPT, encoding="utf-8").read()
    assert "analysis_base_bmi_repaired.csv" in txt
    assert '"analysis_base.csv"' not in txt.replace("analysis_base_bmi_repaired.csv", "")


@requires_data
def test_primary_output_is_analysis_of_record():
    """已重跑的主分析输出(Model B 30d/365d)与 analysis of record 一致。"""
    res = os.path.join(PGV("mimic_record_work"), "results", "04_primary_results.csv")
    if not os.path.exists(res):
        pytest.skip("主分析尚未在本机重跑(make mimic)")
    df = pd.read_csv(res)
    r30, r365 = df.iloc[0], df.iloc[2]
    assert int(r30["N"]) == 10561 and int(r30["events_30d"]) == 296 and int(r30["events_365d"]) == 745
    assert abs(r30["HR"] - 0.9794) < 0.005 and abs(r30["P"] - 0.4952) < 0.01
    assert abs(r365["HR"] - 0.9850) < 0.005 and abs(r365["P"] - 0.5062) < 0.01
