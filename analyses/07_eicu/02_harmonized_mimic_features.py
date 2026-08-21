# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
"""MIMIC-IV harmonized replication: glucose features + landmark outcomes.
输入: 02_mimic_harmonized/{frame_main,frame_sensA,frame_sensB}.csv,
      chartevents_glucose_raw.csv, labevents_glucreat_raw.csv,
      icustays.csv, admissions.csv (只读引用)
输出: 02_mimic_harmonized/mimic_glucose_features.csv, mimic_outcomes.csv, mimic_creatinine.csv,
      qc_mimic_extraction.csv
规则与 eICU 侧完全一致(protocol §4-5)。
"""
import pandas as pd, numpy as np, json, sys

BASE = _os.path.join(PGV("replication_work"), "02_mimic_harmonized")
MIMIC = PGV("mimic_raw_csv")

# ---------- frames ----------
frames = {}
for name in ["main","sensA","sensB"]:
    frames[name] = pd.read_csv(f"{BASE}/frame_{name}.csv", low_memory=False)
cohort = frames["sensB"][["subject_id","hadm_id","stay_id","first_icu_intime","admittime_ts","dischtime_ts"]].copy()
cohort["hadm_id"] = cohort["hadm_id"].astype(int); cohort["stay_id"] = cohort["stay_id"].astype(int)

# ---------- icustays intime 校验 ----------
icu = pd.read_csv(f"{MIMIC}/icu/icustays.csv", usecols=["stay_id","intime","outtime"], parse_dates=["intime","outtime"])
cohort = cohort.merge(icu, on="stay_id", how="left")
cohort["first_icu_intime"] = pd.to_datetime(cohort["first_icu_intime"])
intime_match = (cohort["first_icu_intime"] == cohort["intime"]).mean()
print(f"QC: first_icu_intime == icustays.intime 比例 = {intime_match:.4f}")
cohort["intime"] = cohort["intime"].fillna(cohort["first_icu_intime"])

# ---------- glucose 行 ----------
src_map_ce = {220621:"icu_charted",225664:"poct",226537:"icu_charted",228388:"icu_charted"}
src_map_le_glu = {50931:"central_lab",52569:"central_lab",50809:"blood_gas",52027:"blood_gas"}
creat_ids = {50912:"central_lab",52546:"central_lab",52024:"blood_gas"}

ce = pd.read_csv(f"{BASE}/chartevents_glucose_raw.csv")
ce["charttime"] = pd.to_datetime(ce["charttime"])
ce["source_cat"] = ce["itemid"].map(src_map_ce)
ce = ce[["stay_id","charttime","itemid","valuenum","valueuom","source_cat"]].dropna(subset=["valuenum"])

le = pd.read_csv(f"{BASE}/labevents_glucreat_raw.csv")
le["charttime"] = pd.to_datetime(le["charttime"])
le_glu = le[le["itemid"].isin(src_map_le_glu)].copy()
le_glu["source_cat"] = le_glu["itemid"].map(src_map_le_glu)
le_glu = le_glu[["hadm_id","charttime","itemid","valuenum","valueuom","source_cat"]].dropna(subset=["valuenum"])
le_creat = le[le["itemid"].isin(creat_ids)].copy()

# 单位核查
units = pd.concat([ce[["valueuom"]], le_glu[["valueuom"]]])["valueuom"].value_counts(dropna=False)
print("QC: glucose 单位分布:", units.to_dict())

stay_intime = cohort.set_index("stay_id")["intime"]
hadm_intime = cohort.drop_duplicates("hadm_id").set_index("hadm_id")["intime"]

ce["intime"] = ce["stay_id"].map(stay_intime)
le_glu["intime"] = le_glu["hadm_id"].map(hadm_intime)

def to_offset(df):
    df["offset_min"] = (df["charttime"] - df["intime"]).dt.total_seconds()/60.0
    return df
ce = to_offset(ce); le_glu = to_offset(le_glu)

ce_w = ce[(ce.offset_min>=0)&(ce.offset_min<=1440)&(ce.valuenum>=20)&(ce.valuenum<=1500)]
le_w = le_glu[(le_glu.offset_min>=0)&(le_glu.offset_min<=1440)&(le_glu.valuenum>=20)&(le_glu.valuenum<=1500)]
print(f"QC: chartevents 原始 {len(ce)}, 窗口+合理值后 {len(ce_w)}; labevents 原始 {len(le_glu)}, 窗口+合理值后 {len(le_w)}")

# labevents 经 hadm 映射到 cohort stay(每 hadm 一个 cohort stay)
hadm2stay = cohort.drop_duplicates("hadm_id").set_index("hadm_id")["stay_id"]
le_w["stay_id"] = le_w["hadm_id"].map(hadm2stay)
allg = pd.concat([
    ce_w.assign(row_src="chartevents")[["stay_id","offset_min","valuenum","source_cat","row_src"]].rename(columns={"valuenum":"value_mgdl"}),
    le_w.assign(row_src="labevents")[["stay_id","offset_min","valuenum","source_cat","row_src"]].rename(columns={"valuenum":"value_mgdl"})
], ignore_index=True).dropna(subset=["stay_id"])
allg["stay_id"] = allg["stay_id"].astype(int)
allg["offset_min"] = allg["offset_min"].round().astype(int)  # 分钟粒度
n_raw = len(allg)

# 去重: 完全重复 -> 同分钟中位数
dedup = allg.drop_duplicates(subset=["stay_id","offset_min","value_mgdl","source_cat"])
n_exact_dup = n_raw - len(dedup)
minute = (dedup.groupby(["stay_id","offset_min"], as_index=False)
          .agg(value_mgdl=("value_mgdl","median"),
               sources=("source_cat", lambda s: "+".join(sorted(set(s)))),
               n_records_in_minute=("value_mgdl","size")))
multi_min = int((minute.n_records_in_minute>1).sum())
print(f"QC: raw={n_raw}, 完全重复删除 {n_exact_dup}, 分钟时间点 {len(minute)}, 多记录分钟 {multi_min}")

def feats(g):
    v = g["value_mgdl"]
    return pd.Series({
        "glucose_n": len(v),
        "glucose_mean_24h": v.mean(),
        "glucose_sd_24h": v.std(ddof=1),
        "glucose_cv_24h": v.std(ddof=1)/v.mean() if len(v)>1 else np.nan,
        "glucose_min_24h": v.min(), "glucose_max_24h": v.max(),
        "glucose_range_24h": v.max()-v.min(),
        "glucose_first": g.loc[g.offset_min.idxmin(),"value_mgdl"],
        "glucose_last": g.loc[g.offset_min.idxmax(),"value_mgdl"],
        "measurement_span_minutes": g.offset_min.max()-g.offset_min.min(),
    })
gf = minute.groupby("stay_id").apply(feats, include_groups=False).reset_index()
raw_counts = allg.groupby("stay_id").size().rename("glucose_n_raw")
gf = gf.merge(raw_counts, on="stay_id", how="left")
src_frac = (minute.assign(poct=minute.sources.str.contains("poct").astype(int),
                          cl=minute.sources.str.contains("central_lab").astype(int),
                          bg=minute.sources.str.contains("blood_gas").astype(int),
                          icu=minute.sources.str.contains("icu_charted").astype(int))
            .groupby("stay_id")[["poct","cl","bg","icu"]].mean()
            .rename(columns={"poct":"poct_fraction","cl":"central_lab_fraction","bg":"blood_gas_fraction","icu":"icu_charted_fraction"}))
gf = gf.merge(src_frac, on="stay_id", how="left")
gf.to_csv(f"{BASE}/mimic_glucose_features.csv", index=False)
minute.to_csv(f"{BASE}/mimic_glucose_minute.csv", index=False)
print("glucose features stays:", len(gf), "ge2:", (gf.glucose_n>=2).sum())

# ---------- creatinine first 0-24h ----------
le_creat["intime"] = le_creat["hadm_id"].map(hadm_intime)
le_creat["offset_min"] = (pd.to_datetime(le_creat["charttime"]) - le_creat["intime"]).dt.total_seconds()/60.0
cr = le_creat[(le_creat.offset_min>=0)&(le_creat.offset_min<=1440)].dropna(subset=["valuenum"])
cr = cr.sort_values("offset_min").drop_duplicates("hadm_id", keep="first")[["hadm_id","valuenum"]].rename(columns={"valuenum":"creatinine"})
cr.to_csv(f"{BASE}/mimic_creatinine.csv", index=False)
print("creatinine avail:", len(cr))

# ---------- landmark & outcomes ----------
adm = pd.read_csv(f"{MIMIC}/hosp/admissions.csv", usecols=["hadm_id","dischtime","deathtime","hospital_expire_flag"], parse_dates=["dischtime","deathtime"])
o = cohort[["subject_id","hadm_id","stay_id","intime"]].merge(adm, on="hadm_id", how="left")
o["hosp_mortality"] = o.hospital_expire_flag==1
o["died_within_24h"] = o.hosp_mortality & (o.deathtime <= o.intime + pd.Timedelta(minutes=1440))
o["discharged_within_24h"] = (~o.hosp_mortality) & (o.dischtime <= o.intime + pd.Timedelta(minutes=1440))
o["in_landmark"] = ~(o.died_within_24h | o.discharged_within_24h)
o["post_landmark_hosp_mortality"] = o.hosp_mortality & o.in_landmark
# ICU mortality: MIMIC 无直接 icu expire 字段 -> 用 icu outtime 与 deathtime 对齐
o = o.merge(icu[["stay_id","outtime"]], on="stay_id", how="left")
o["icu_mortality"] = o.hosp_mortality & o.deathtime.notna() & (o.deathtime <= o.outtime + pd.Timedelta(minutes=60))
o["icu_los_lt_24h"] = (o.outtime - o.intime) < pd.Timedelta(minutes=1440)
o["icu_los_days"] = (o.outtime - o.intime).dt.total_seconds()/86400
o["hosp_los_days"] = (o.dischtime - o.intime).dt.total_seconds()/86400
o.to_csv(f"{BASE}/mimic_outcomes.csv", index=False)
for nm in ["main","sensA","sensB"]:
    ids = set(frames[nm].stay_id.astype(int))
    sub = o[o.stay_id.isin(ids)]
    print(nm, "n=",len(sub), "died24=",int(sub.died_within_24h.sum()), "disch24=",int(sub.discharged_within_24h.sum()),
          "landmark=",int(sub.in_landmark.sum()), "landmark_deaths=",int(sub.post_landmark_hosp_mortality.sum()))
with open(f"{BASE}/qc_extraction_summary.json","w") as f:
    json.dump({"intime_match":float(intime_match),"glucose_units":units.to_dict(),"raw_rows":n_raw,
               "exact_dups_removed":int(n_exact_dup),"minute_points":len(minute),"multi_record_minutes":multi_min,
               "units_ce": ce.valueuom.value_counts(dropna=False).to_dict(),
               "units_le": le_glu.valueuom.value_counts(dropna=False).to_dict()}, f, indent=2, default=str)
print("DONE")
