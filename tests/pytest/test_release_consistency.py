#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_release_consistency.py — 发布一致性测试(阶段四增强验证)。

覆盖:
1. 展项溯源完整性(43/43,唯一);
2. AF2 图源测试(脚本存在、仓库相对路径、无旧工作区路径、声明输出);
3. AF 编号测试(禁止旧 Additional file 9/10 出现在发布文件);
4. 主文-补充-CSV 一致性(冻结锚定 vs provenance 关键估计);
5. analysis-of-record 集成(区间模型、INSPIRE v5 与冻结值一致)。

需要授权数据/投稿包的用例在无配置时自动 skip(CI 安全)。
"""
import os
import re

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path.insert(0, os.path.join(ROOT, "src", "common"))
from paths import PGV  # noqa: E402

PROV = os.path.join(ROOT, "docs", "FINAL_OUTPUT_PROVENANCE.csv")

EXPECTED_EXHIBITS = (
    [f"Main Table {i}" for i in range(1, 6)]
    + [f"Main Figure {i}" for i in range(1, 6)]
    + [f"Table S{i}" for i in range(1, 26)]
    + [f"Figure S{i}" for i in range(1, 9)]
)

# ---------- 1. 展项溯源完整性 ----------
def test_exhibit_provenance_complete_and_unique():
    df = pd.read_csv(PROV)
    ids = df["output_id"].tolist()
    assert len(ids) == len(set(ids)) == 43, "展项不唯一或数量不为 43"
    missing = [e for e in EXPECTED_EXHIBITS if e not in ids]
    assert not missing, f"缺失展项: {missing}"
    for col in ["database", "generating_script", "source_result", "verification_status"]:
        assert df[col].notna().all(), f"列 {col} 存在空值"


def _pkg():
    try:
        p = PGV("manuscript_pkg")
        return p if os.path.isdir(p) else None
    except KeyError:
        return None


requires_pkg = pytest.mark.skipif(_pkg() is None, reason="未配置 manuscript_pkg(投稿包路径)")


# ---------- 2. AF2 图源测试 ----------
@requires_pkg
def test_figure_source_map_repo_relative():
    af2 = pd.read_csv(os.path.join(_pkg(), "Final_Additional_file_2_Figure_Source_Map.csv"))
    figs = af2["figure"].tolist()
    assert figs == [f"Main Figure {i}" for i in range(1, 6)] + \
                   [f"Supplementary Figure S{i}" for i in range(1, 9)], "AF2 必须恰为 13 个图展项"
    legacy = re.compile(r"/Users/|/home/|/Volumes/|cardiac_glucose|inspire_cardiac|Codex/|kimi/")
    for _, r in af2.iterrows():
        for script in re.split(r"\s*\+\s*", r["generating_script_repo"]):
            p = os.path.join(ROOT, script)
            assert os.path.exists(p), f"AF2 脚本不存在: {script}"
            txt = open(p, encoding="utf-8", errors="ignore").read()
            assert not legacy.search(txt), f"AF2 脚本含旧工作区路径: {script}"
        assert isinstance(r["output_repo"], str) and r["output_repo"], "AF2 缺少 output_repo"
        assert not r["output_repo"].startswith("/"), "output_repo 必须为相对路径或配置键形式"


# ---------- 3. AF 编号测试 ----------
OLD_AF = re.compile(r"Additional[_ ]file[_ ]?(9|10)\b", re.I)

def test_no_legacy_af_numbering_in_release_files():
    scan_roots = ["docs", "config", "tests", "README.md", "DATA_AVAILABILITY.md",
                  "CHANGELOG.md", os.path.join("results", "qc")]
    hits = []
    for sr in scan_roots:
        p = os.path.join(ROOT, sr)
        files = [p] if os.path.isfile(p) else [
            os.path.join(dp, f) for dp, _, fns in os.walk(p) for f in fns
            if f.endswith((".md", ".yml", ".csv", ".txt"))]
        for f in files:
            for i, ln in enumerate(open(f, encoding="utf-8", errors="ignore"), 1):
                if OLD_AF.search(ln) and not re.search(r"legacy|废弃|不得引用|已清除|~~", ln, re.I):
                    hits.append((f, i, ln.strip()[:80]))
    pkg = _pkg()
    if pkg:
        for fn in os.listdir(pkg):
            if fn.endswith(".csv"):
                for i, ln in enumerate(open(os.path.join(pkg, fn), encoding="utf-8", errors="ignore"), 1):
                    if OLD_AF.search(ln):
                        hits.append((fn, i, ln.strip()[:80]))
    assert not hits, f"旧 AF9/AF10 编号残留: {hits}"


# ---------- 4. 主文-补充-CSV 一致性 ----------
def _parse_hr_ci_p(s):
    m = re.search(r"HR\s*([0-9.]+)\s*\(([0-9.]+)[–-]([0-9.]+)\)\s*P\s*=\s*([0-9.]+)", str(s))
    return tuple(float(g) for g in m.groups()) if m else None


def test_provenance_key_estimates_match_frozen():
    prov = pd.read_csv(PROV).set_index("output_id")
    checks = {
        "Main Table 2": [(0.979, 0.922, 1.040, 0.495), (0.985, 0.942, 1.030, 0.506)],
        "Main Table 3": [(0.905, 0.653, 1.254, 0.547), (1.104, 0.846, 1.441, 0.464)],
    }
    for exhibit, expected_list in checks.items():
        ke = prov.loc[exhibit, "key_estimate"]
        found = [_parse_hr_ci_p(part) for part in re.split(r";", ke)]
        found = [f for f in found if f]
        assert len(found) == len(expected_list), f"{exhibit}: 解析到 {len(found)} 组 HR"
        for got, exp in zip(found, expected_list):
            for g, e in zip(got, exp):
                assert abs(g - e) < 0.0015, f"{exhibit}: {got} != {exp}"
    shr = prov.loc["Main Table 5", "key_estimate"]
    assert "1.844" in shr and "0.158" in shr and "2.338" in shr and "0.097" in shr
    t4 = prov.loc["Main Table 4", "key_estimate"]
    assert "1.158" in t4 and "1.025" in t4 and "0.019" in t4


def _authorized():
    try:
        return os.path.isdir(PGV("mimic_record_work"))
    except KeyError:
        return False


requires_data = pytest.mark.skipif(not _authorized(), reason="需要授权分析工作目录")


# ---------- 5. analysis-of-record 集成 ----------
@requires_data
def test_interval_and_logtime_are_aor():
    p = os.path.join(PGV("mimic_record_work"), "results", "PH365_INTERVAL_MICE_POOLED.csv")
    if not os.path.exists(p):
        pytest.skip("区间模型尚未重跑")
    df = pd.read_csv(p).set_index("interval")
    exp = {"d1-7": (1.063, 162, 10561), "d8-30": (0.974, 139, 10399), "d31-365": (0.938, 444, 10260)}
    for iv, (hr, ev, nr) in exp.items():
        r = df.loc[iv]
        assert abs(r["HR_per10"] - hr) < 0.005 and int(r["interval_events"]) == ev and int(r["n_at_risk_start"]) == nr
    lt = pd.read_csv(os.path.join(PGV("mimic_record_work"), "results", "PH365_GV_LOGTIME_SENSITIVITY_MICE_POOLED.csv"))
    tt = lt[lt["term"] == "tt(gv10)"].iloc[0]
    assert abs(tt["HR"] - 0.966) < 0.005


@requires_data
def test_inspire_v5_is_aor():
    try:
        d = PGV("inspire_record_work")
    except KeyError:
        pytest.skip("未配置 inspire_record_work")
    p30 = os.path.join(d, "30d_primary_results_v5.csv")
    if not os.path.exists(p30):
        pytest.skip("INSPIRE v5 尚未重跑")
    df = pd.read_csv(p30)
    r = df[df["model_id"] == "ADMINV5_I2_30d"].iloc[0]
    assert int(r["N"]) == 1353 and int(r["events"]) == 27
    assert abs(r["HR_per10"] - 0.9046) < 0.005 and abs(r["P"] - 0.5471) < 0.01
    w = pd.read_csv(os.path.join(d, "365d_withdrawal_v5.csv"))
    assert "WITHDRAWN" in w["status"].iloc[0]
