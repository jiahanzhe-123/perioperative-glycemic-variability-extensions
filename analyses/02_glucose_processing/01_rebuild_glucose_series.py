# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# 01_glucose_series_rebuild.py
# 血糖序列彻底重构:互斥来源分类、同分钟优先级、两套敏感性序列、分层去重审计。
# 输入(患者级,本机):00_glucose_series/{chartevents,labevents}_glucose.csv
#   + final_dataset_for_v3_pipeline.csv(surgery_time 锚点)
# 输出:data/series_{priority,allmedian,poctfirst}.csv
#       data/features_{priority,allmedian,poctfirst}.csv
#       results/glucose_dedup_audit.csv
#       data/excluded_glucose_records_audit.csv
# 验证:allmedian 序列特征必须复现冻结 final_glucose_features.csv(连续性检查)。
import os, json
import numpy as np
import pandas as pd

PROJ = PGV("mimic_work")
OUT  = PGV("mimic_record_work")
os.makedirs(os.path.join(OUT,"data"), exist_ok=True)
os.makedirs(os.path.join(OUT,"results"), exist_ok=True)

PRIO = {"central_lab":1, "blood_gas":2, "poct":3, "icu_charted":4}
LO, HI = 20.0, 1500.0

# ---------- 1. 读取 ----------
ce = pd.read_csv(os.path.join(PROJ,"00_glucose_series/chartevents_glucose.csv"), parse_dates=["charttime"])
lab = pd.read_csv(os.path.join(PROJ,"00_glucose_series/labevents_glucose.csv"), parse_dates=["charttime"])
g = pd.concat([ce, lab], ignore_index=True)
g["source_cat"] = g["source_cat"].astype(str)
assert set(g["source_cat"].unique()) <= set(PRIO), g["source_cat"].unique()
audit = {"raw_records_total": len(g),
         "raw_by_source": g["source_cat"].value_counts().to_dict()}

coh = pd.read_csv(os.path.join(PROJ,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"),
                  usecols=["stay_id","surgery_time"], low_memory=False)
coh["t0"] = pd.to_datetime(coh["surgery_time"])
coh = coh.drop_duplicates("stay_id")
assert coh["stay_id"].is_unique
g = g.merge(coh[["stay_id","t0"]], on="stay_id", how="inner")
audit["records_after_cohort_join"] = len(g)

# ---------- 2. 窗口 [t0, t0+24h) ----------
inwin = (g["charttime"] >= g["t0"]) & (g["charttime"] < g["t0"] + pd.Timedelta(hours=24))
excl_window = g[~inwin].copy(); excl_window["reason"] = "outside [surgical_date_index, +24h)"
g = g[inwin].copy()
audit["records_in_24h_window"] = len(g)
audit["excluded_outside_window"] = int((~inwin).sum())

# ---------- 3. 值域 20–1500 ----------
bad = (g["valuenum"] < LO) | (g["valuenum"] > HI)
excl_range = g[bad].copy(); excl_range["reason"] = "value outside 20-1500 mg/dL"
g = g[~bad].copy()
audit["excluded_out_of_range"] = int(bad.sum())

# ---------- 4. 完全重复 ----------
dup_exact = g.duplicated(subset=["stay_id","charttime","valuenum","source_cat"], keep="first")
excl_dup = g[dup_exact].copy(); excl_dup["reason"] = "exact duplicate (stay,time,value,source)"
g = g[~dup_exact].copy()
audit["excluded_exact_duplicates"] = int(dup_exact.sum())

# ---------- 5. 同 stay+time+value 跨表(仅审计,随后由优先级规则自然解决) ----------
xtab = g.duplicated(subset=["stay_id","charttime","valuenum"], keep=False)
audit["crosstable_same_time_value_record_rows"] = int(xtab.sum())

# ---------- 6. 跨表近似重复抑制(与冻结管线一致:chartevents 若在 ±2min 内有 ≤1 mg/dL 的 lab,剔除) ----------
labv = g[g["source_cat"].isin(["central_lab","blood_gas"])][["stay_id","charttime","valuenum"]].rename(
    columns={"charttime":"lab_time","valuenum":"lab_val"})
ce_idx = g["source_cat"].isin(["poct","icu_charted"])
ces = g[ce_idx][["stay_id","charttime","valuenum"]].copy()
ces["rid"] = ces.index
m = ces.merge(labv, on="stay_id")
m = m[(m["lab_time"] >= m["charttime"] - pd.Timedelta(minutes=2)) & (m["lab_time"] <= m["charttime"] + pd.Timedelta(minutes=2))]
m["adiff"] = (m["lab_val"] - m["valuenum"]).abs()
near_dup_ids = set(m.loc[m["adiff"] <= 1.0, "rid"])
audit["excluded_crosstable_near_duplicates_2min_1mgdl"] = len(near_dup_ids)
excl_near = g.loc[sorted(near_dup_ids)].copy(); excl_near["reason"] = "cross-table near-duplicate (lab within 2 min, |diff|<=1)"
g = g.drop(index=list(near_dup_ids))

# ±1min / ±5min 同值跨表(只报告,不删除)
rep = {}
for w in (1, 5):
    mm = ces.merge(labv, on="stay_id")
    mm = mm[(mm["lab_time"] >= mm["charttime"] - pd.Timedelta(minutes=w)) & (mm["lab_time"] <= mm["charttime"] + pd.Timedelta(minutes=w)) &
            (mm["lab_val"] == mm["valuenum"])]
    rep[f"same_value_cross_table_within_pm{w}min_rows"] = int(mm["rid"].nunique())
audit.update(rep)

# ---------- 7. 同分钟多记录统计 ----------
g["minute"] = g["charttime"].dt.floor("min")
per_min = g.groupby(["stay_id","minute"])
audit["stay_minutes_with_multi_records"] = int((per_min.size() > 1).sum())
audit["stay_minutes_total"] = int(per_min.size().sum())
diffs = per_min["valuenum"].agg(lambda v: v.max()-v.min() if len(v)>1 else 0.0)
audit["same_min_absdiff_median"] = float(diffs.median())
audit["same_min_absdiff_p95"] = float(diffs.quantile(.95))
audit["same_min_absdiff_max"] = float(diffs.max())

# ---------- 8. 三套序列 ----------
def series_priority(df):
    df = df.copy(); df["pr"] = df["source_cat"].map(PRIO)
    best = df.groupby(["stay_id","minute"])["pr"].transform("min")
    df = df[df["pr"] == best]
    s = df.groupby(["stay_id","minute"]).agg(value=("valuenum","median"), final_source_class=("source_cat","first")).reset_index()
    return s
def series_allmedian(df):
    s = df.groupby(["stay_id","minute"]).agg(value=("valuenum","median"),
        final_source_class=("source_cat", lambda x: "+".join(sorted(set(x))))).reset_index()
    return s
def series_poctfirst(df):
    out = []
    for (st, mi), gg in df.groupby(["stay_id","minute"]):
        p = gg[gg["source_cat"]=="poct"]
        use = p if len(p) else gg
        cls = "poct" if len(p) else "+".join(sorted(set(gg["source_cat"])))
        out.append((st, mi, use["valuenum"].median(), cls))
    return pd.DataFrame(out, columns=["stay_id","minute","value","final_source_class"])

ser = {}
ser["priority"]  = series_priority(g)
ser["allmedian"] = series_allmedian(g)
ser["poctfirst"] = series_poctfirst(g)

# ---------- 9. 特征 ----------
def feats(s):
    def one(gg):
        gg = gg.sort_values("minute"); v = gg["value"].to_numpy(float); n = len(v)
        mins = gg["minute"].values.astype("datetime64[s]").astype(np.int64)/60.0
        o = {"glucose_count": n, "mean_glucose": v.mean()}
        if n >= 2:
            sd = v.std(ddof=1)
            o.update(gv_sd=sd, cv=sd/v.mean(), min_glucose=v.min(), max_glucose=v.max(),
                     range_glucose=v.max()-v.min(), iqr_glucose=np.percentile(v,75)-np.percentile(v,25),
                     mad_glucose=np.median(np.abs(v-np.median(v))),
                     arv=np.abs(np.diff(v)).mean(),
                     tw_arv=(np.abs(np.diff(v))*np.diff(mins)).sum()/(mins[-1]-mins[0]) if mins[-1]>mins[0] else np.nan,
                     span_hours=(mins[-1]-mins[0])/60.0,
                     median_interval_min=np.median(np.diff(mins)) if n>=2 else np.nan)
        else:
            o.update(gv_sd=np.nan, cv=np.nan, min_glucose=v.min(), max_glucose=v.max(), range_glucose=0.0,
                     iqr_glucose=0.0, mad_glucose=0.0, arv=np.nan, tw_arv=np.nan, span_hours=np.nan, median_interval_min=np.nan)
        o.update(any_lt54=bool((v<54).any()), any_lt70=bool((v<70).any()),
                 any_gt180=bool((v>180).any()), any_gt250=bool((v>250).any()),
                 prop_lt70=(v<70).mean(), prop_gt180=(v>180).mean(), prop_70_180=((v>=70)&(v<=180)).mean())
        return pd.Series(o)
    f = s.groupby("stay_id").apply(one, include_groups=False).reset_index()
    # 互斥来源比例(每分钟归属)
    src = s.groupby("stay_id")["final_source_class"].apply(lambda x: x.value_counts(normalize=True))
    fr = src.unstack(fill_value=0.0)
    for c in ["central_lab","blood_gas","poct","icu_charted"]:
        f["frac_"+c] = f["stay_id"].map(fr[c] if c in fr.columns else 0.0).fillna(0.0)
    return f

feat = {k: feats(s) for k, s in ser.items()}
for k, s in ser.items():
    s.to_csv(os.path.join(OUT,"data",f"series_{k}.csv"), index=False)
    feat[k].to_csv(os.path.join(OUT,"data",f"features_{k}.csv"), index=False)

# ---------- 10. 排除记录审计 ----------
excl = pd.concat([excl_window, excl_range, excl_dup, excl_near], ignore_index=True)
excl.to_csv(os.path.join(OUT,"data","excluded_glucose_records_audit.csv"), index=False)
audit["excluded_total"] = len(excl)
audit["final_records_after_dedup"] = len(g)
audit["final_minutes_priority"] = len(ser["priority"])
audit["final_minutes_allmedian"] = len(ser["allmedian"])
audit["final_minutes_poctfirst"] = len(ser["poctfirst"])
audit["stays_with_any_record"] = int(g["stay_id"].nunique())

# 患者级测量数/跨度/间隔
cnt = ser["priority"].groupby("stay_id").size()
audit["per_patient_count_median"] = float(cnt.median()); audit["per_patient_count_p95"] = float(cnt.quantile(.95))
audit["per_patient_count_max"] = int(cnt.max())
sp = feat["priority"]["span_hours"]
audit["span_hours_median"] = float(sp.median()); audit["span_hours_p10"] = float(sp.quantile(.10))
audit["median_interval_min_median"] = float(feat["priority"]["median_interval_min"].median())
audit["priority_source_minutes_frac"] = ser["priority"]["final_source_class"].value_counts(normalize=True).to_dict()

pd.DataFrame([audit]).T.rename(columns={0:"value"}).to_csv(os.path.join(OUT,"results","glucose_dedup_audit.csv"))

# ---------- 11. 连续性验证:allmedian vs 冻结特征 ----------
ref = pd.read_csv(os.path.join(PROJ,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_glucose_features.csv"))
mm = feat["allmedian"].merge(ref, on="stay_id", suffixes=("","_ref"))
checks = {}
for c_new, c_ref in [("gv_sd","gv_sd"),("mean_glucose","mean_glucose"),("glucose_count","glucose_count"),
                     ("arv","arv"),("tw_arv","tw_arv"),("span_hours","span_hours")]:
    dlt = (mm[c_new] - mm[c_ref]).abs()
    checks[c_new] = {"n_mismatch_gt_1e-6": int((dlt > 1e-6).sum()), "max_abs_diff": float(dlt.max())}
ok = all(v["n_mismatch_gt_1e-6"]==0 for v in checks.values())
with open(os.path.join(OUT,"results","allmedian_vs_frozen_verification.json"),"w") as f:
    json.dump({"checks":checks, "all_pass":ok, "n_stays_merged":len(mm)}, f, indent=2)
print(json.dumps(checks, indent=1))
print("ALLMEDIAN_MATCH_FROZEN:", ok)
print("PHASE1_DONE")
