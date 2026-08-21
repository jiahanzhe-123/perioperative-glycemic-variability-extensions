# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# 07a_source_series.py — 单一来源与共同来源血糖序列及 stay 级特征
# 输入:00_glucose_series per-record 抽取(与 Phase 1 相同清洗规则)
# 输出:data/features_source_{poct,centrallab,bloodgas,common}.csv
#       data/samepatient_poct_lab.csv(同时 >=2 POCT 且 >=2 central lab 患者的双 GV)
import os
import numpy as np
import pandas as pd

PROJ = PGV("mimic_work")
OUT  = PGV("mimic_record_work")

ce = pd.read_csv(os.path.join(PROJ,"00_glucose_series/chartevents_glucose.csv"), parse_dates=["charttime"])
lab = pd.read_csv(os.path.join(PROJ,"00_glucose_series/labevents_glucose.csv"), parse_dates=["charttime"])
g = pd.concat([ce, lab], ignore_index=True)
coh = pd.read_csv(os.path.join(PROJ,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"),
                  usecols=["stay_id","surgery_time"], low_memory=False)
coh["t0"] = pd.to_datetime(coh["surgery_time"]); coh = coh.drop_duplicates("stay_id")
g = g.merge(coh[["stay_id","t0"]], on="stay_id", how="inner")
g = g[(g["charttime"] >= g["t0"]) & (g["charttime"] < g["t0"] + pd.Timedelta(hours=24))]
g = g[(g["valuenum"] >= 20) & (g["valuenum"] <= 1500)]
g = g.drop_duplicates(subset=["stay_id","charttime","valuenum","source_cat"])
# 与 Phase 1 一致的跨表近似重复抑制(lab 优先,±2min ≤1)
labv = g[g["source_cat"].isin(["central_lab","blood_gas"])][["stay_id","charttime","valuenum"]].rename(
    columns={"charttime":"lab_time","valuenum":"lab_val"})
ces = g[g["source_cat"].isin(["poct","icu_charted"])][["stay_id","charttime","valuenum"]].copy()
ces["rid"] = ces.index
m = ces.merge(labv, on="stay_id")
m = m[(m["lab_time"] >= m["charttime"] - pd.Timedelta(minutes=2)) & (m["lab_time"] <= m["charttime"] + pd.Timedelta(minutes=2))]
m["adiff"] = (m["lab_val"] - m["valuenum"]).abs()
drop_ids = set(m.loc[m["adiff"] <= 1.0, "rid"])
g = g.drop(index=list(drop_ids))
g["minute"] = g["charttime"].dt.floor("min")

def build(mask):
    gg = g[mask].copy()
    s = gg.groupby(["stay_id","minute"])["valuenum"].median().reset_index().rename(columns={"valuenum":"value"})
    def one(x):
        v = np.sort(x["value"].to_numpy(float)); n = len(v)
        if n < 2: return pd.Series({"gv": np.nan, "mean_glu": v.mean(), "n_src": n})
        return pd.Series({"gv": v.std(ddof=1), "mean_glu": v.mean(), "n_src": n})
    return s.groupby("stay_id").apply(one, include_groups=False).reset_index()

specs = {
    "poct": g["source_cat"]=="poct",
    "centrallab": g["source_cat"]=="central_lab",
    "bloodgas": g["source_cat"]=="blood_gas",
    "common": g["source_cat"].isin(["poct","central_lab"]),
}
feats = {}
for name, mk in specs.items():
    f = build(mk); feats[name] = f
    f.to_csv(os.path.join(OUT,"data",f"features_source_{name}.csv"), index=False)
    print(name, "stays >=2:", int((f["n_src"]>=2).sum()), "; GV 非缺失:", int(f["gv"].notna().sum()))

sp = feats["poct"].merge(feats["centrallab"], on="stay_id", suffixes=("_poct","_lab"))
both = sp[(sp["n_src_poct"]>=2) & (sp["n_src_lab"]>=2)].copy()
both.to_csv(os.path.join(OUT,"data","samepatient_poct_lab.csv"), index=False)
print("同患者(>=2 POCT 且 >=2 central lab):", len(both))
print("PHASE7A_DONE")
