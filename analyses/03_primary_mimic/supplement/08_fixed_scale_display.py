# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# fixed_scale_display.py — 生成整合版展示 CSV(与 UPDATED_SUPPLEMENTARY_TABLE_DATA 现有格式一致)
import pandas as pd, numpy as np, os

RIV = PGV("mimic_riv")
OUT = os.path.join(RIV, "08_fixed_scale")
TBL = os.path.join(PGV("manuscript_work"), "UPDATED_SUPPLEMENTARY_TABLE_DATA")

def f2(x): return f"{x:.2f}"
def hr_ci(hr, lo, hi): return f"{f2(hr)} ({f2(lo)}\u2013{f2(hi)})"
def pv(p):
    if pd.isna(p): return ""
    return "<0.001" if p < 0.001 else f"{p:.3f}"
CO = {"SHR-GV":"SHR\u2013GV","open-core":"Open-core","GV-only":"GV-only"}
def outc(h): return "30-day mortality" if int(h)==30 else "1-year mortality"

# ---------- S26 ----------
ag = pd.read_csv(os.path.join(OUT, "s26_alt_gv_fixedscale_machine.csv"))
def role(r):
    if r["metric"]=="SD (mg/dL)" and r["min_n"]==2: return "Primary, full frame"
    if r["metric"]=="CV (%)": return "Alternative, full frame"
    if r["paired_role"]=="paired_SD_comparator":
        return "Paired comparator (ARV subset)" if "ARV (mg/dL)" in r["metric"] and "time-weighted" not in r["metric"] else "Paired comparator (time-weighted ARV subset)"
    if r["metric"]=="ARV (mg/dL)": return "Alternative, ARV subset"
    return "Alternative, time-weighted ARV subset"
def metric_disp(r):
    return "SD (mg/dL)" if r["paired_role"]=="paired_SD_comparator" else r["metric"]
rows = []
for (co,h), g in ag.groupby(["cohort","horizon_days"], sort=False):
    paired = g["paired_role"]=="paired_SD_comparator"
    order = [
        g[(g["metric"]=="SD (mg/dL)")&(g["min_n"]==2)],
        g[g["metric"]=="CV (%)"],
        g[paired & g["metric"].str.contains("ARV \(mg/dL\)", regex=True) & ~g["metric"].str.contains("time-weighted")],
        g[g["metric"]=="ARV (mg/dL)"],
        g[paired & g["metric"].str.contains("time-weighted")],
        g[g["metric"]=="time-weighted ARV (mg/dL)"],
    ]
    for part in order:
        for r in part.itertuples():
            fixed = (r.metric=="SD (mg/dL)" and r.min_n==2)
            hr  = r.HR_clinical_fixed if fixed else r.HR_clinical
            lo  = r.lo1_fixed if fixed else r.lo1
            hi  = r.hi1_fixed if fixed else r.hi1
            hr2 = r.HR_plus_mean_fixed if fixed else r.HR_plus_mean
            lo2 = r.lo2_fixed if fixed else r.lo2
            hi2 = r.hi2_fixed if fixed else r.hi2
            rd = r._asdict()
            rows.append({"Cohort":CO[r.cohort],"Outcome":outc(h),"GV metric":metric_disp(rd),"Comparator role":role(rd),
                "N":int(r.N),"Events":int(r.events),
                "HR (95% CI)":hr_ci(hr,lo,hi),"P value":pv(r.P_clinical),
                "Mean-glucose-adjusted HR (95% CI); P value":f"{hr_ci(hr2,lo2,hi2)}; {pv(r.P_plus_mean)}"})
s26 = pd.DataFrame(rows)
assert len(s26)==36, f"S26 行数异常: {len(s26)}"
s26.to_csv(os.path.join(TBL,"Supplementary_Table_S26_alternative_gv_metrics_fixedscale.csv"), index=False)

# ---------- S27A ----------
ms = pd.read_csv(os.path.join(OUT, "s27a_sequential_fixedscale_machine.csv"))
adj = {"none":"None","count":"+ measurement count","span":"+ measurement span",
       "count_span":"+ count and span","count_span_source":"+ count, span, and source mix"}
rows = [{"Cohort":CO[r.cohort],"Outcome":outc(r.horizon_days),"Adjustment":adj[r.adjustment],
         "N":int(r.N),"Events":int(r.events),"HR (95% CI)":hr_ci(r.HR_gv_fixed,r.lo_fixed,r.hi_fixed),
         "P value":pv(r.P)} for r in ms.itertuples()]
s27a = pd.DataFrame(rows)
with open(os.path.join(TBL,"Supplementary_Table_S27A_sequential_adjustment_fixedscale.csv"), "w", encoding="utf-8") as f:
    f.write(",".join(["Panel A. Sequential adjustment for measurement-process summaries"]*len(s27a.columns))+"\n")
    s27a.to_csv(f, index=False)

# ---------- S27B(数值不变,仅确认与现表一致) ----------
mr = pd.read_csv(os.path.join(RIV,"01_frame_repair/measurement/measurement_restriction_analyses_framefix.csv"))
rest = {">=3 glucose":"\u22653 glucose measurements",">=4 glucose":"\u22654 glucose measurements",
        "span>=12h":"Measurement span \u226512 h","POCT-dominant (>=50% POCT)":"POCT-dominant (\u226550% POCT)",
        "POCT-only GV (>=2 POCT)":"POCT-only GV (\u22652 POCT)","central-lab-only GV (>=2 lab)":"Central-laboratory-only GV (\u22652 lab)"}
rows = []
for r in mr.itertuples():
    if pd.isna(getattr(r,"HR_gv",np.nan)):
        cell, p = ("Not estimable (<100 N or <20 events)" if "NOT ESTIMABLE" in str(getattr(r,"note","")) else ""), ""
    else:
        cell, p = hr_ci(r.HR_gv, r.lo, r.hi), pv(r.P)
    rows.append({"Cohort":CO[r.cohort],"Outcome":outc(r.horizon_days),"Adjustment":rest[r.restriction],
                 "N":int(r.N),"Events":int(r.events),"HR (95% CI)":cell,"P value":p})
s27b = pd.DataFrame(rows)
with open(os.path.join(TBL,"Supplementary_Table_S27B_measurement_restrictions_fixedscale.csv"), "w", encoding="utf-8") as f:
    f.write(",".join(["Panel B. Sample restrictions by measurement availability and source"]*len(s27b.columns))+"\n")
    s27b.to_csv(f, index=False)

# ---------- S28 ----------
ex = pd.read_csv(os.path.join(OUT, "s28_extremes_fixedscale_machine.csv"))
adjx = {"none":"None","hypo":"+ hypoglycemia burden","hyper":"+ hyperglycemia burden",
        "minmax":"+ minimum and maximum glucose","burden":"+ extreme-glucose burden"}
rows = []
for r in ex.itertuples():
    if pd.isna(getattr(r,"HR_gv",np.nan)):
        cell, p = "", ""
    else:
        cell, p = hr_ci(r.HR_gv_fixed, r.lo_fixed, r.hi_fixed), pv(r.P_gv)
    rows.append({"Cohort":CO[r.cohort],"Outcome":outc(r.horizon_days),"Adjustment":adjx[r.adjustment],
                 "N":int(r.N),"Events":int(r.events),"HR (95% CI)":cell,"P value":p})
s28 = pd.DataFrame(rows)
s28.to_csv(os.path.join(TBL,"Supplementary_Table_S28_extreme_glucose_adjustments_fixedscale.csv"), index=False)

# ---------- S33 ----------
hb = pd.read_csv(os.path.join(OUT, "s33_hba1c_fixedscale_machine.csv"))
rows = []
for r in hb.itertuples():
    rule = "Main hierarchical rule" if r.analysis.startswith("all") else "Strictly pre-index"
    for ex_, hr, lo, hi, p in [("SHR",r.HR_shr_fixed,r.lo_shr_fixed,r.hi_shr_fixed,r.P_shr),
                               ("GV",r.HR_gv_fixed,r.lo_gv_fixed,r.hi_gv_fixed,r.P_gv)]:
        rows.append({"Outcome":outc(r.horizon_days),"HbA1c eligibility rule":rule,"N":int(r.N),"Events":int(r.events),
                     "Exposure":f"{ex_} per 1-SD","HR (95% CI)":hr_ci(hr,lo,hi),"P value":pv(p)})
s33 = pd.DataFrame(rows)
s33.to_csv(os.path.join(TBL,"Supplementary_Table_S33_hba1c_preindex_sensitivity_fixedscale.csv"), index=False)

# ---------- 对照:新旧展示值 ----------
print("=== S28 None 行(新 vs 旧) ===")
old28 = pd.read_csv(os.path.join(TBL,"Supplementary_Table_S28_extreme_glucose_adjustments.csv"), keep_default_na=False)
new28 = s28[s28["Adjustment"]=="None"]
for r in new28.itertuples():
    o = old28[(old28["Cohort"]==r.Cohort)&(old28["Outcome"]==r.Outcome)&(old28["Adjustment"]=="None")].iloc[0]
    print(f"  {r.Cohort:10s} {r.Outcome:18s}: 旧 {o['HR (95% CI)']}  →  新 {r._6}")
print("=== S33 全表(新) ===")
print(s33.to_string(index=False))
print("=== S26 SD 行(新) ===")
print(s26[s26["Comparator role"]=="Primary, full frame"][["Cohort","Outcome","HR (95% CI)","P value","Mean-glucose-adjusted HR (95% CI); P value"]].to_string(index=False))

# S33 主文句(3 位小数风格)
def f3(x): return f"{x:.3f}"
m = hb[hb["analysis"].str.startswith("strict")].iloc[0]
y = hb[hb["analysis"].str.startswith("strict")].iloc[1]
print("\n=== 主文 Results 新句(S33,3 位小数) ===")
print(f"30-day: N=4,613, 92 events, SHR HR {f3(m.HR_shr_fixed)} ({f3(m.lo_shr_fixed)}\u2013{f3(m.hi_shr_fixed)}), P=0.056; "
      f"GV HR {f3(m.HR_gv_fixed)} ({f3(m.lo_gv_fixed)}\u2013{f3(m.hi_gv_fixed)}), P=0.150; "
      f"1-year GV HR {f3(y.HR_gv_fixed)} ({f3(y.lo_gv_fixed)}\u2013{f3(y.hi_gv_fixed)}), P=0.014")
