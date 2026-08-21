#!/usr/bin/env python3
# validate_results.py — 冻结结果验证:对照 tests/expected/frozen_results.yml 与授权数据中的冻结 CSV。
# 规则:N/事件数必须完全相等;HR/RR/OR/CI 容差 ±0.005;P 容差 ±0.01(绝对)。
# 失败不覆盖旧结果,输出 unresolved discrepancy 清单,退出码非零。
import os, sys, math
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "common"))
from paths import PGV, PGV_OUT

def load_expected(path):
    # 极简 YAML 解析(本仓库受控格式:一级 key 为节,二级为字段)
    import re
    exp, cur = {}, None
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"): continue
        m1 = re.match(r"^([A-Za-z0-9_]+):\s*$", ln)
        m2 = re.match(r"^  ([A-Za-z0-9_]+):\s*(.*)$", ln)
        if m1:
            cur = m1.group(1); exp[cur] = {}
        elif m2 and cur:
            k, v = m2.group(1), m2.group(2).strip().strip('"').strip("'")
            try: v = float(v) if ("." in v or "e" in v.lower()) else int(v)
            except ValueError: pass
            exp[cur][k] = v
    return exp

def read_section(key):
    src = EXPECTED[key].get("source_csv")
    if not src: return None
    root_key, rel = src.split(":", 1)
    try:
        p = os.path.join(PGV(root_key), rel)
    except KeyError:
        return None  # 键未在 config/paths.yml 配置(无授权数据环境)
    return pd.read_csv(p) if os.path.exists(p) else None

def close(a, b, tol):
    if a is None or b is None or (isinstance(a, float) and math.isnan(a)): return None
    return abs(float(a) - float(b)) <= tol

def eq(a, b):
    return int(a) == int(b)

report = []
def check(section, item, expected, actual, tol, exact=False, note=""):
    if actual is None:
        report.append((section, item, expected, "MISSING", "UNRESOLVED", note)); return False
    okc = eq(actual, expected) if exact else close(actual, expected, tol)
    report.append((section, item, expected, actual, "PASS" if okc else "MISMATCH", note))
    return bool(okc)

EXPECTED = load_expected(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "expected", "frozen_results.yml"))
allok = True

# --- MIMIC primary ---
df = read_section("mimic_primary")
if df is None:
    report.append(("mimic_primary", "source_csv", "mimic_record_work:results/04_primary_results.csv", "MISSING FILE", "UNRESOLVED", "需要授权数据路径(config/paths.yml)"))
    allok = False
else:
    r = df.iloc[0]; r365 = df.iloc[2]
    allok &= check("mimic_primary", "n", EXPECTED["mimic_primary"]["n"], r["N"], 0, exact=True)
    allok &= check("mimic_primary", "deaths_30d", EXPECTED["mimic_primary"]["deaths_30d"], r["events_30d"], 0, exact=True)
    allok &= check("mimic_primary", "deaths_365d", EXPECTED["mimic_primary"]["deaths_365d"], r["events_365d"], 0, exact=True)
    allok &= check("mimic_primary", "model_b_30d_hr", EXPECTED["mimic_primary"]["model_b_30d_hr"], r["HR"], 0.005)
    allok &= check("mimic_primary", "model_b_30d_ci_lower", EXPECTED["mimic_primary"]["model_b_30d_ci_lower"], r["lo"], 0.005)
    allok &= check("mimic_primary", "model_b_30d_ci_upper", EXPECTED["mimic_primary"]["model_b_30d_ci_upper"], r["hi"], 0.005)
    allok &= check("mimic_primary", "model_b_30d_p", EXPECTED["mimic_primary"]["model_b_30d_p"], r["P"], 0.01)
    allok &= check("mimic_primary", "model_b_365d_hr", EXPECTED["mimic_primary"]["model_b_365d_hr"], r365["HR"], 0.005)
    allok &= check("mimic_primary", "model_b_365d_p", EXPECTED["mimic_primary"]["model_b_365d_p"], r365["P"], 0.01)
    allok &= check("mimic_primary", "gv_sd_constant", EXPECTED["mimic_primary"]["gv_sd_constant"], r["gv_sd_constant"], 0.005)

# --- MIMIC absolute risk ---
df = read_section("mimic_absolute_risk_30d")
if df is not None:
    r = df[df["horizon_days"] == 30].iloc[0] if "horizon_days" in df.columns else df.iloc[0]
    allok &= check("mimic_absolute_risk_30d", "rd", EXPECTED["mimic_absolute_risk_30d"]["rd"], r["rd"], 0.005)
    allok &= check("mimic_absolute_risk_30d", "risk_q25", EXPECTED["mimic_absolute_risk_30d"]["risk_q25"], r["risk_q25"], 0.005)
    allok &= check("mimic_absolute_risk_30d", "risk_q75", EXPECTED["mimic_absolute_risk_30d"]["risk_q75"], r["risk_q75"], 0.005)

# --- MIMIC time-dependent ---
df = read_section("mimic_time_dependent_365d")
if df is not None:
    gv = df[(df["exposure"] == "gv_z") & (df["cohort"].isin(["A", "open_extended"]))]
    def iv(name):
        row = gv[gv["time_interval"] == name]
        return row["HR"].iloc[0] if len(row) else None
    allok &= check("mimic_time_dependent_365d", "interval_0_1_gv_hr",
                   EXPECTED["mimic_time_dependent_365d"]["interval_0_1_gv_hr"], iv("0-1 days"), 0.005)
    allok &= check("mimic_time_dependent_365d", "interval_2_7_gv_hr",
                   EXPECTED["mimic_time_dependent_365d"]["interval_2_7_gv_hr"], iv("2-7 days"), 0.005)
    allok &= check("mimic_time_dependent_365d", "interval_8_30_gv_hr",
                   EXPECTED["mimic_time_dependent_365d"]["interval_8_30_gv_hr"], iv("8-30 days"), 0.005)

# --- MIMIC SHR-GV joint ---
df = read_section("mimic_shr_gv_joint")
if df is not None:
    r30 = df[df["model_id"] == "OMNI_J1_vs_J0_30d"].iloc[0]
    r365 = df[df["model_id"] == "OMNI_J1_vs_J0_365d"].iloc[0]
    allok &= check("mimic_shr_gv_joint", "omnibus_30d_stat", EXPECTED["mimic_shr_gv_joint"]["omnibus_30d_stat"], r30["stat"], 0.005)
    allok &= check("mimic_shr_gv_joint", "omnibus_30d_p", EXPECTED["mimic_shr_gv_joint"]["omnibus_30d_p"], r30["P"], 0.01)
    allok &= check("mimic_shr_gv_joint", "omnibus_365d_stat", EXPECTED["mimic_shr_gv_joint"]["omnibus_365d_stat"], r365["stat"], 0.005)
    allok &= check("mimic_shr_gv_joint", "omnibus_365d_p", EXPECTED["mimic_shr_gv_joint"]["omnibus_365d_p"], r365["P"], 0.01)

# --- INSPIRE v5 ---
df = read_section("inspire_primary_v5")
if df is not None:
    r = df[df["model_id"] == "ADMINV5_I2_30d"].iloc[0]
    allok &= check("inspire_primary_v5", "landmark_n", EXPECTED["inspire_primary_v5"]["landmark_n"], r["N"], 0, exact=True)
    allok &= check("inspire_primary_v5", "deaths_30d", EXPECTED["inspire_primary_v5"]["deaths_30d"], r["events"], 0, exact=True)
    allok &= check("inspire_primary_v5", "i2_hr_30d", EXPECTED["inspire_primary_v5"]["i2_hr_30d"], r["HR_per10"], 0.005)
    allok &= check("inspire_primary_v5", "i2_p", EXPECTED["inspire_primary_v5"]["i2_p"], r["P"], 0.01)

df = read_section("inspire_48h")
if df is not None:
    r = df.iloc[0]
    allok &= check("inspire_48h", "hr", EXPECTED["inspire_48h"]["hr"], r["HR_per10"], 0.005)
    allok &= check("inspire_48h", "p", EXPECTED["inspire_48h"]["p"], r["P"], 0.01)

# --- eICU ---
df = read_section("eicu_harmonized_m3")
if df is not None:
    r = df[df["model_id"] == "HARM_eICU-CRD_M3"].iloc[0]
    allok &= check("eicu_harmonized_m3", "n", EXPECTED["eicu_harmonized_m3"]["n"], r["N"], 0, exact=True)
    allok &= check("eicu_harmonized_m3", "events", EXPECTED["eicu_harmonized_m3"]["events"], r["events"], 0, exact=True)
    allok &= check("eicu_harmonized_m3", "rr_per10", EXPECTED["eicu_harmonized_m3"]["rr_per10"], r["RR_per10"], 0.005)
    allok &= check("eicu_harmonized_m3", "p", EXPECTED["eicu_harmonized_m3"]["p"], r["P"], 0.01)

df = read_section("eicu_random_intercept")
if df is not None:
    r = df.iloc[0]
    allok &= check("eicu_random_intercept", "n", EXPECTED["eicu_random_intercept"]["n"], r["N"], 0, exact=True)
    allok &= check("eicu_random_intercept", "events", EXPECTED["eicu_random_intercept"]["events"], r["events"], 0, exact=True)
    allok &= check("eicu_random_intercept", "hospitals", EXPECTED["eicu_random_intercept"]["hospitals_expected"], r["hospitals"], 0, exact=True)
    allok &= check("eicu_random_intercept", "or_per10", EXPECTED["eicu_random_intercept"]["or_per10"], r["OR_per10"], 0.005)
    allok &= check("eicu_random_intercept", "p", EXPECTED["eicu_random_intercept"]["p"], r["P"], 0.01)

out_dir = PGV_OUT("results/qc")
os.makedirs(out_dir, exist_ok=True)
rep = pd.DataFrame(report, columns=["section", "item", "expected", "actual", "status", "note"])
rep.to_csv(os.path.join(out_dir, "frozen_validation_report.csv"), index=False)
n_pass = (rep["status"] == "PASS").sum()
n_mis = (rep["status"] == "MISMATCH").sum()
n_unres = (rep["status"] == "UNRESOLVED").sum()
print(rep.to_string(index=False))
print(f"\nfrozen validation: PASS={n_pass}, MISMATCH={n_mis}, UNRESOLVED={n_unres}")
sys.exit(0 if (n_mis == 0 and n_unres == 0 and n_pass > 0) else 1)
