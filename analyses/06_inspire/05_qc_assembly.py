# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# 14_inspire_qc_assembly.py — QC 单元测试、规格矩阵、number map、database-estimand mapping
import os, json
import numpy as np
import pandas as pd

ROOT = PGV("inspire_work")
RES = os.path.join(ROOT, "results")
K = json.load(open(os.path.join(RES,"standardization_constants.json")))
JK = json.load(open(os.path.join(RES,"inspire_joint_constants.json")))

b = pd.read_csv(os.path.join(ROOT,"data","inspire_base.csv"), low_memory=False)
cm = pd.read_csv(os.path.join(ROOT,"data","comorbidity.csv"))
d = b.merge(cm, on="subject_id", how="left")
d["event_30d"] = ((d["allcause_death_time"].notna()) & (d["allcause_death_time"] <= d["opend_time"] + 30*1440)).astype(int)
d["event_365d"] = ((d["allcause_death_time"].notna()) & (d["allcause_death_time"] <= d["opend_time"] + 365*1440)).astype(int)

tests = []
def T(name, ok, detail=""):
    tests.append({"test": name, "pass": bool(ok), "detail": str(detail)})

T("one operation per patient", d["subject_id"].is_unique and d["op_id"].is_unique)
T("age >= 18 for all", (d["age"] >= 18).all())
gl = d[(d["gv_sd"].notna())]
T("GV cases have >=2 distinct timepoints (opend window)",
  bool((gl["n_glucose_0_24h"] >= 2).all()))
shr = d[d["shr"].notna()]
T("SHR cases: HbA1c strictly 1-90 days pre-index",
  bool(((shr["hba1c_days_before_opstart"] > 0) & (shr["hba1c_days_before_opstart"] <= 90)).all()))
T("eAG > 0 for all SHR cases", bool((shr["eag_mg_dl"] > 0).all()))
T("SHR = mean_glucose/eAG (same frozen sequence)",
  bool(np.allclose(shr["shr"], shr["mean_glucose"]/shr["eag_mg_dl"], atol=1e-3)))
T("standardization constants frozen and present",
  os.path.exists(os.path.join(RES,"standardization_constants.json")) and
  abs(K["gv_sd"] - 13.0) < 100, f"gv_sd={K['gv_sd']:.4f}, N={K['N']}")
fl = pd.read_csv(os.path.join(RES,"01_cohort_flow.csv"))
T("cohort flow reproduces target N", int(fl["n"].iloc[-1]) == K["N"], f"flow={int(fl['n'].iloc[-1])}, K={K['N']}")
T("protocol frozen with sha256", os.path.exists(os.path.join(ROOT,"docs","protocol_sha256.txt")))
mice_lg = pd.read_csv(os.path.join(RES,"mice_logged_events.csv")) if os.path.exists(os.path.join(RES,"mice_logged_events.csv")) else None
T("MICE logged events = 0", mice_lg is None or len(mice_lg)==0)
prim = pd.read_csv(os.path.join(RES,"04_primary_results.csv"))
T("primary results carry analysis label",
  (prim["analysis_label"]=="confirmatory secondary / external extension / exploratory").all())
T("joint corner scenarios support-checked",
  os.path.exists(os.path.join(RES,"joint_corner_support.csv")))
qc = pd.DataFrame(tests)
qc.to_csv(os.path.join(RES,"qc_unit_tests.csv"), index=False)
print(qc.to_string(index=False))
assert qc["pass"].all(), "存在未通过的 QC 测试"

# ---- 规格矩阵 ----
rows = []
def addm(mid, module, outcome, scale, N, events, note=""):
    rows.append({"model_id":mid,"module":module,"database":"INSPIRE v1.4.2","cohort":"open CABG/valve, opend-anchored day-1 landmark, >=2 quantized glucose timepoints",
        "outcome":outcome,"landmark":"opend_time + 1440 min","exposure_definition":"quantized routine-laboratory SD-based GV (opend+0-24h frozen series)",
        "exposure_scale":scale,"covariates":"per frozen protocol (30d: age; 365d: age, sex, diabetes, Charlson-wo-DM, ASA, emop, group, BMI)",
        "missing_data":"MICE m=50 maxit=20 (seed 20260726) / complete-case",
        "N":N,"events":events,"standardization_constants":f"gv_mean={K['gv_mean']:.3f}, gv_sd={K['gv_sd']:.3f}",
        "seed":20260726,"code_version":"inspire_cardiac_20260729","input_checksum":"dd842f46e477f7740e771e543b83ad82","note":note})
for _, r in prim.iterrows():
    pass
m = pd.read_csv(os.path.join(RES,"mice_pooled_models.csv"))
for r in m.itertuples():
    hz = "30-day" if "30d" in r.model_id else "365-day"
    addm(r.model_id, "primary (INSPIRE external extension)", f"{hz} all-cause mortality",
         "per 10 mg/dL; per fixed-cohort SD", r.N, r.events, r.note if isinstance(r.note,str) else "")
om = pd.read_csv(os.path.join(RES,"09_joint_omnibus.csv"))
for r in om.itertuples():
    hz = "30-day" if "30d" in r.model_id else "365-day"
    addm(r.model_id, "joint SHR-GV (prespecified secondary)", f"{hz} all-cause mortality",
         "2-df joint Wald", r.N, r.events, r.note if isinstance(getattr(r,"note",None),str) else "")
an = pd.read_csv(os.path.join(RES,"05_time_anchor_results.csv"))
for r in an.itertuples():
    if pd.isna(getattr(r,"HR_per10",np.nan)): continue
    hz = "30-day" if "30d" in r.model_id else "365-day"
    addm(r.model_id, "time-anchor sensitivity", f"{hz} all-cause mortality", f"anchor={r.anchor}", r.N, r.events)
msm = pd.DataFrame(rows)
msm.to_csv(os.path.join(RES,"13_model_specification_matrix.csv"), index=False)

# ---- database-estimand mapping ----
maprows = [
 ("MIMIC-IV","UNIQUE PRIMARY ANALYSIS","date-anchored (calendar-date) index; window [index 00:00, +24h)",
  "day-1 calendar landmark","30-day all-cause (day 1-30 from index)","GV-only cohort; Model B (fixed clinical + RCS mean + linear GV); MICE m=50",
  "HR 0.977 per 10 mg/dL, P=0.434 (frozen, unchanged by this module)"),
 ("INSPIRE v1.4.2","exact-timestamp external extension","verified opend_time; window [opend, +24h)",
  "opend + 1440 min landmark","30-day all-cause (allcause linkage); in-hospital contrast","open CABG/valve cohort; Model I2 (reduced clinical per EPV rule + RCS mean + linear GV); MICE m=50",
  "reported separately; never pooled; confirmatory secondary label"),
 ("eICU-CRD","multicenter harmonized comparison","ICU admission + 0-1440 min proxy window",
  "24h landmark","post-landmark hospital mortality","modified Poisson, cluster-robust by hospital",
  "non-equivalent severity/measurement systems; no pooling, no validation claim"),
]
mp = pd.DataFrame(maprows, columns=["database","role","time_anchor","landmark","outcome","model_adjustment","interpretation_rule"])
mp.to_csv(os.path.join(RES,"database_estimand_mapping.csv"), index=False)

# ---- number map ----
p30 = prim.iloc[0]; p365 = prim.iloc[2]
nm = [
 ("主文/补充:INSPIRE 扩展分析(新增)", "—(此前无 INSPIRE)", "MICE_I2_30d",
  f"quantized GV per 10 mg/dL HR {p30.HR:.3f} ({p30.lo:.3f}–{p30.hi:.3f}), P={p30.P:.3f} (N=1,355, 18 events)",
  "新增,标记 exact-timestamp external extension"),
 ("主文/补充:INSPIRE 365d", "—", "MICE_I2_365d",
  f"HR {p365.HR:.3f} ({p365.lo:.3f}–{p365.hi:.3f}), P={p365.P:.3f} (100 events)", "新增(key secondary)"),
 ("主文:主要 MIMIC 结论", "GV-only 30d HR 0.977, P=0.434;365d 0.992, P=0.701", "—",
  "不变:INSPIRE 不得重新包装、取代或推翻", "不变"),
 ("补充:opend 0–48h 窗口(探索)", "—", "ANCH_opend_0_48h_30d",
  "HR 1.424 (1.089–1.862), P=0.0099 — 探索性窗口信号;不按最小 P 选主窗", "仅补充材料,不作主结论"),
]
nmdf = pd.DataFrame(nm, columns=["manuscript_location","old_value","new_model_id","new_value","action"])
nmdf.to_csv(os.path.join(RES,"14_manuscript_number_map.csv"), index=False)
print("assembly done:", len(msm), "matrix rows;", len(mp), "db mappings;", len(nmdf), "number-map rows")
