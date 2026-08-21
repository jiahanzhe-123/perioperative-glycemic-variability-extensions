# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# fixed_scale_constants.py — 计算主文(v3)固定标准化常数,验证后输出常数表
import pandas as pd, numpy as np, os, json

PROJ = PGV("mimic_work")
RIV  = PGV("mimic_riv")
FZ   = os.path.join(PROJ, "final_statistical_freeze_20260727")
OUT  = os.path.join(RIV, "08_fixed_scale")
os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(os.path.join(FZ, "01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"), low_memory=False)
for c in ["gv","shr","glucose_mean_postop_24h"]:
    d[c] = pd.to_numeric(d[c], errors="coerce")
v3 = pd.read_csv(os.path.join(FZ, "03_dependent_analyses_refit/v3_pipeline/tables/final_v3_primary_models.csv"))
prim = pd.read_csv(os.path.join(FZ, "02_primary_models_refit/final_primary_models.csv"))

pools = {"SHR-GV":"final_v2_A","open-core":"final_v2_open_core","GV-only":"final_v2_C"}
v3co  = {"SHR-GV":"A","open-core":"open_core","GV-only":"C"}

# ---- 1. 池级(主文/v3 固定)常数 ----
const = {}
for co, flag in pools.items():
    pool = d[d[flag]==1]
    const[co] = {v: {"mean": float(pool[v].mean()), "sd": float(pool[v].std()), "N": int(pool[v].notna().sum())}
                 for v in (["gv","shr","glucose_mean_postop_24h"] if co!="GV-only" else ["gv","glucose_mean_postop_24h"])}

# ---- 2. frame 级常数(用 frame 注册表的 stay 清单,精确) ----
reg = pd.read_csv(os.path.join(RIV, "04_frame_registry/frame_registry.csv"))
frame_sd = {}
for r in reg.itertuples():
    ids = set(int(x) for x in open(os.path.join(RIV, f"04_frame_registry/{r.frame_id}_stay_ids.txt")).read().split())
    sub = d[d["stay_id"].isin(ids)]
    frame_sd[(r.cohort, r.horizon_days)] = {v: float(sub[v].std()) for v in ["gv","shr","glucose_mean_postop_24h"]}

# ---- 3. strictly pre-index 子集常数(SHR-GV) ----
sub33 = {}
for h in (30, 365):
    ids = set(int(x) for x in open(os.path.join(RIV, f"04_frame_registry/FRAME_SHRGV_{h}d_stay_ids.txt")).read().split())
    sub = d[d["stay_id"].isin(ids) & (d["hba1c_window"]=="pre_90d")]
    sub33[h] = {"N": len(sub), "gv": float(sub["gv"].std()), "shr": float(sub["shr"].std())}

# ---- 4. 验证:HR_frame^(pool_sd/frame_sd) 必须等于 v3 ----
rows = []
ok_all = True
for co in ["SHR-GV","open-core","GV-only"]:
    for h in (30, 365):
        for ex in (["gv","shr"] if co!="GV-only" else ["gv"]):
            v3r = v3[(v3["cohort"]==v3co[co])&(v3["horizon_days"]==h)&(v3["exposure"]==("SHR per SD" if ex=="shr" else "GV per SD"))]
            pr  = prim[(prim["cohort"]==co)&(prim["horizon_days"]==h)]
            hr_v3 = float(v3r["HR"].iloc[0]); lo_v3 = float(v3r["lower95"].iloc[0]); hi_v3 = float(v3r["upper95"].iloc[0])
            col = "HR_shr_z" if ex=="shr" else "HR_gv_z"
            hr_fr = float(pr[col].iloc[0])
            k = const[co][ex]["sd"] / frame_sd[(co,h)][ex]
            hr_pred = hr_fr ** k
            ok = abs(hr_pred - hr_v3) < 1e-6
            ok_all &= ok
            rows.append({"cohort":co,"horizon":h,"exposure":ex,"k":k,"HR_frame":hr_fr,"HR_pred":hr_pred,
                         "HR_v3":hr_v3,"match":ok,"lo_v3":lo_v3,"hi_v3":hi_v3})
            print(f"{co:10s} {h:3d}d {ex}: k={k:.6f}  pred={hr_pred:.6f}  v3={hr_v3:.6f}  {'OK' if ok else 'MISMATCH'}")

assert ok_all, "v3 常数验证失败"
print("\n全部 10 行验证通过:主文固定常数 = 队列池级 mean/sd")

out = {"pool_constants": const,
       "frame_sd": {f"{co}|{h}": v for (co,h), v in frame_sd.items()},
       "strictly_preindex_subset": {str(h): v for h, v in sub33.items()}}
with open(os.path.join(OUT, "fixed_scale_constants.json"), "w") as f:
    json.dump(out, f, indent=2)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "fixed_scale_validation_vs_v3.csv"), index=False)
print("常数表与验证已写入", OUT)
print("\n=== strictly pre-index 子集 ===")
for h, v in sub33.items():
    print(f"  {h}d: N={v['N']}, sd_gv={v['gv']:.6f}, sd_shr={v['shr']:.6f}")
