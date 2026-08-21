# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# fixed_scale_rescale.py — 按主文(v3)固定标准化常数重缩放 S26/S27A/S28/S33 的 per-1-SD HR 与 CI
# 规则(作者 2026-07-28 裁定):
#   - 与主模型同一 frame 的 SD 模型 → 使用主文固定常数(队列池级 mean/sd);
#   - ARV/tw-ARV 配对子集与测量限制子集 → 保留子集内标准化,标注 "within the subset";
#   - 只缩放 log-HR 与 CI(k = s_main/s_fit),z 统计量与 P 不变。
import pandas as pd, numpy as np, os, json

PROJ = PGV("mimic_work")
RIV  = PGV("mimic_riv")
FR   = os.path.join(RIV, "01_frame_repair")
OUT  = os.path.join(RIV, "08_fixed_scale")
RF   = PGV("manuscript_work")
TBL  = os.path.join(RF, "UPDATED_SUPPLEMENTARY_TABLE_DATA")
os.makedirs(OUT, exist_ok=True)

C = json.load(open(os.path.join(OUT, "fixed_scale_constants.json")))
pool  = C["pool_constants"]
fsd   = C["frame_sd"]
sub33 = C["strictly_preindex_subset"]

def k_gv(co, h):  return pool[co]["gv"]["sd"]  / fsd[f"{co}|{h}"]["gv"]
def k_shr(co, h): return pool[co]["shr"]["sd"] / fsd[f"{co}|{h}"]["shr"]
def k_mean(co, h):return pool[co]["glucose_mean_postop_24h"]["sd"] / fsd[f"{co}|{h}"]["glucose_mean_postop_24h"]
def resc(hr, k):  return hr ** k

FIXED_LBL = lambda co, ex: f"per 1 SD (fixed primary-analysis constant; {ex} SD={pool[co][ex.lower()]['sd']:.4f}" + (" mg/dL)" if ex=="GV" else ")")
SUBSET_LBL = "per 1 SD within the analytic subset"

# ============ S26 alt_gv ============
ag = pd.read_csv(os.path.join(FR, "alt_gv/alt_gv_models_framefix.csv"))
ag["scale_rule"] = SUBSET_LBL; ag["k_exposure"] = 1.0; ag["k_mean"] = 1.0
mask = (ag["metric"]=="SD (mg/dL)") & (ag["min_n"]==2)
for i in ag[mask].index:
    co, h = ag.at[i,"cohort"], int(ag.at[i,"horizon_days"])
    ag.at[i,"k_exposure"] = k_gv(co,h); ag.at[i,"k_mean"] = k_mean(co,h)
    ag.at[i,"scale_rule"] = FIXED_LBL(co,"GV")
for col in ["HR_clinical","lo1","hi1","HR_plus_mean","lo2","hi2"]:
    ag[col+"_fixed"] = resc(ag[col], ag["k_exposure"])
ag["mean_HR_fixed"] = resc(ag["mean_HR"], ag["k_mean"])
ag.to_csv(os.path.join(OUT, "s26_alt_gv_fixedscale_machine.csv"), index=False)

# ============ S27A measurement sequential ============
ms = pd.read_csv(os.path.join(FR, "measurement/measurement_sequential_adjustment_framefix.csv"))
ms["k_exposure"] = [k_gv(r.cohort, int(r.horizon_days)) for r in ms.itertuples()]
for col in ["HR_gv","lo","hi"]:
    ms[col+"_fixed"] = resc(ms[col], ms["k_exposure"])
ms["scale_rule"] = [FIXED_LBL(r.cohort,"GV") for r in ms.itertuples()]
ms.to_csv(os.path.join(OUT, "s27a_sequential_fixedscale_machine.csv"), index=False)

# ============ S28 extremes ============
ex = pd.read_csv(os.path.join(FR, "extremes/extremes_extreme_adjustments_framefix.csv"))
ex["k_exposure"] = [k_gv(r.cohort, int(r.horizon_days)) if not pd.isna(r.HR_gv) else 1.0 for r in ex.itertuples()]
for col in ["HR_gv","lo","hi"]:
    ex[col+"_fixed"] = resc(ex[col], ex["k_exposure"])
ex["scale_rule"] = [FIXED_LBL(r.cohort,"GV") if not pd.isna(r.HR_gv) else "" for r in ex.itertuples()]
ex.to_csv(os.path.join(OUT, "s28_extremes_fixedscale_machine.csv"), index=False)

# ============ S33 HbA1c ============
hb = pd.read_csv(os.path.join(RIV, "02_hba1c_window_sensitivity/hba1c_window_sensitivity_models.csv"))
def k_hb(r):
    h = int(r.horizon_days)
    if r.analysis.startswith("all_windows"):
        return {"shr": k_shr("SHR-GV",h), "gv": k_gv("SHR-GV",h)}
    s = sub33[str(h)]
    return {"shr": pool["SHR-GV"]["shr"]["sd"]/s["shr"], "gv": pool["SHR-GV"]["gv"]["sd"]/s["gv"]}
hb["k_shr"] = [k_hb(r)["shr"] for r in hb.itertuples()]
hb["k_gv"]  = [k_hb(r)["gv"]  for r in hb.itertuples()]
for col in ["HR_shr","lo_shr","hi_shr"]:
    hb[col+"_fixed"] = resc(hb[col], hb["k_shr"])
for col in ["HR_gv","lo_gv","hi_gv"]:
    hb[col+"_fixed"] = resc(hb[col], hb["k_gv"])
hb["scale_rule"] = ["per 1 SD (fixed primary-analysis constant; both rules share the same constants)"
                    if r.analysis.startswith("all") else
                    "per 1 SD (fixed primary-analysis constant; rescaled from subset fit)" for r in hb.itertuples()]
hb.to_csv(os.path.join(OUT, "s33_hba1c_fixedscale_machine.csv"), index=False)

# ============ 验证锚点 ============
anchors = []
v3ref = pd.read_csv(os.path.join(OUT, "fixed_scale_validation_vs_v3.csv"))
def get_v3(co,h,ex_):
    r = v3ref[(v3ref["cohort"]==co)&(v3ref["horizon"]==h)&(v3ref["exposure"]==ex_)].iloc[0]
    return r["HR_v3"], r["lo_v3"], r["hi_v3"]
print("=== S28/S27A None 行 vs v3(应全等) ===")
for t, name in [(ms,"S27A"), (ex,"S28")]:
    none = t[t["adjustment"]=="none"]
    for r in none.itertuples():
        v = get_v3(r.cohort, int(r.horizon_days), "gv")
        ok = abs(r.HR_gv_fixed - v[0]) < 1e-6
        anchors.append((name, r.cohort, r.horizon_days, r.HR_gv_fixed, v[0], ok))
        print(f"  {name} {r.cohort:10s} {int(r.horizon_days):3d}d: {r.HR_gv_fixed:.6f} vs v3 {v[0]:.6f}  {'OK' if ok else 'FAIL'}")
print("=== S26 Primary SD 行 vs v3 ===")
for r in ag[mask].itertuples():
    v = get_v3(r.cohort, int(r.horizon_days), "gv")
    ok = abs(r.HR_clinical_fixed - v[0]) < 1e-6
    anchors.append(("S26", r.cohort, r.horizon_days, r.HR_clinical_fixed, v[0], ok))
    print(f"  S26 {r.cohort:10s} {int(r.horizon_days):3d}d: {r.HR_clinical_fixed:.6f} vs v3 {v[0]:.6f}  {'OK' if ok else 'FAIL'}")
print("=== S33 主规则行 vs v3 ===")
for r in hb[hb["analysis"].str.startswith("all")].itertuples():
    for ex_ in (["gv","shr"]):
        v = get_v3("SHR-GV", int(r.horizon_days), ex_)
        got = r.HR_gv_fixed if ex_=="gv" else r.HR_shr_fixed
        ok = abs(got - v[0]) < 1e-6
        anchors.append(("S33", "SHR-GV", r.horizon_days, got, v[0], ok))
        print(f"  S33 {int(r.horizon_days):3d}d {ex_}: {got:.6f} vs v3 {v[0]:.6f}  {'OK' if ok else 'FAIL'}")
assert all(a[5] for a in anchors), "锚点验证失败"
print("\n全部锚点通过。机器版 CSV 已写入", OUT)
