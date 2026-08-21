# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# 02_build_analysis_dataset.py
# 组装重构版分析数据集:互斥手术类别、diabetes + Charlson-excluding-diabetes、
# landmark 变量、固定协变量集、队列流程。
# 输入:final_dataset_for_v3_pipeline.csv + harmonized charlson_wo_diabetes +
#       data/features_priority.csv(Phase 1 输出)
import os, json, hashlib
import numpy as np
import pandas as pd

PROJ = PGV("mimic_work")
OUT  = PGV("mimic_record_work")

d = pd.read_csv(os.path.join(PROJ,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"), low_memory=False)
for c in ["gv","shr","glucose_mean_postop_24h","age_at_admission","bmi","charlson_comorbidity_index",
          "lactate_postop_first","creat_postop_first","wbc_postop_first","albumin_adm_first",
          "hgb_postop_first","platelets_postop_first","sofa_24h","survival_time_days",
          "postop_30d_death_flag","one_year_death_flag","diabetes","hba1c_pct","hba1c_days_from_surgery"]:
    if c in d.columns: d[c] = pd.to_numeric(d[c], errors="coerce")
d["t0"] = pd.to_datetime(d["surgery_time"])

# ---- 1. 互斥手术类别 ----
for f in ["cabg_flag","open_valve_flag","aortic_surgery_flag","transplant_vad_flag","congenital_cardiac_flag"]:
    d[f] = d[f].astype(str).str.lower().isin(["true","t","1"])
def proc_cat(r):
    if r["aortic_surgery_flag"]: return "open aortic surgery (+/- other)"
    if r["transplant_vad_flag"]: return "transplant/VAD"
    if r["cabg_flag"] and r["open_valve_flag"]: return "combined CABG + open valve"
    if r["cabg_flag"]: return "isolated CABG"
    if r["open_valve_flag"]: return "isolated open valve"
    return "congenital/other open cardiac"
d["procedure_cat6"] = d.apply(proc_cat, axis=1)
xtab = pd.crosstab(d["procedure_cat6"], d["procedure_group_main"], margins=True)
xtab.to_csv(os.path.join(OUT,"results","procedure_cat6_vs_legacy_crosstab.csv"))
assert (d.groupby("procedure_cat6").size().sum() == len(d)), "手术互斥类别未覆盖全队列"

# ---- 2. Charlson excluding diabetes(全队列可审计推导) ----
# 来源:fix_data_layer 的行级 ICD 修复(覆盖全部 12,992 stays)。
# 权重:complicated=+2,uncomplicated(any 且非 complicated)=+1,无糖尿病=0。
# charlson_without_diabetes = charlson_comorbidity_index - diabetes_weight。
# 另用 harmonized 文件(10,904 stays)做逐行交叉验证。
fx = pd.read_csv(os.path.join(PGV("mimic_sql_project"),"outputs","fix_data_layer","tables","fix_diabetes_icd_row_level.csv"))
fx = fx.drop_duplicates("stay_id")
for c in ["diabetes_icd_any_fixed","diabetes_icd_without_complication_fixed","diabetes_icd_with_complication_fixed"]:
    fx[c] = fx[c].astype(str).str.lower().isin(["true","t","1"])
d = d.merge(fx[["stay_id","diabetes_icd_any_fixed","diabetes_icd_without_complication_fixed","diabetes_icd_with_complication_fixed"]],
            on="stay_id", how="left", validate="many_to_one")
n_miss_fx = int(d["diabetes_icd_any_fixed"].isna().sum())
assert n_miss_fx == 0, f"ICD 修复表缺失 {n_miss_fx} stays"
# 当前 diabetes 列应与 ICD any 修复值一致(血缘核对)
cur_dm = pd.to_numeric(d["diabetes"], errors="coerce").fillna(0).astype(int)
assert (cur_dm == d["diabetes_icd_any_fixed"].astype(int)).all(), "当前 diabetes 列与 ICD 修复值不一致"
d["diabetes_weight"] = np.where(d["diabetes_icd_with_complication_fixed"], 2,
                        np.where(d["diabetes_icd_any_fixed"], 1, 0))
d["charlson_without_diabetes"] = d["charlson_comorbidity_index"] - d["diabetes_weight"]
assert (d["charlson_without_diabetes"] >= 0).all(), "charlson_wo_dm 出现负值"
# 交叉验证:与 harmonized 文件逐行一致
cc = pd.read_csv(os.path.join(PGV("method_audit_work"),"00_audit","mimic_model_data_charlson_wo_diabetes.csv"),
                 usecols=["stay_id","charlson_without_diabetes"])
cc = cc.drop_duplicates("stay_id")
mm = d.merge(cc, on="stay_id", how="inner", suffixes=("","_harm"))
n_x = len(mm); n_mismatch = int((mm["charlson_without_diabetes"] != mm["charlson_without_diabetes_harm"]).sum())
xval = {"harmonized_overlap_stays": n_x, "mismatch_rows": n_mismatch, "match": n_mismatch == 0}
assert n_mismatch == 0, f"charlson_wo_dm 与 harmonized 不一致 {n_mismatch} 行"
with open(os.path.join(OUT,"results","charlson_wo_dm_crossvalidation.json"),"w") as f:
    json.dump(xval, f, indent=2)

# ---- 3. landmark 变量 ----
d["landmark_time"] = d["t0"] + pd.Timedelta(hours=24)
d["t_lm_30"]  = np.minimum(d["survival_time_days"], 30.0) - 1.0
d["t_lm_365"] = np.minimum(d["survival_time_days"], 365.0) - 1.0
d["event_lm_30"]  = (pd.to_numeric(d["postop_30d_death_flag"], errors="coerce") == 1).astype(int)
d["event_lm_365"] = (pd.to_numeric(d["one_year_death_flag"], errors="coerce") == 1).astype(int)
d["landmark_eligible"] = d["survival_time_days"] >= 1.0   # landmark 时存活且有随访
assert (d.loc[d["landmark_eligible"], "t_lm_30"] >= 0).all(), "landmark 后随访时间<0"
assert (d.loc[~d["landmark_eligible"], "survival_time_days"] < 1).all()

# ---- 4. HbA1c 时间层级(互斥) ----
def a1c_tier(r):
    if pd.isna(r["hba1c_pct"]): return "no HbA1c"
    dd = r["hba1c_days_from_surgery"]
    if dd < 0 and dd >= -90: return "1-90 days pre-index"
    if dd < -90 and dd >= -365: return "91-365 days pre-index"
    if dd >= 0 and dd < 1: return "index calendar day"
    return "outside windows"
d["hba1c_tier"] = d.apply(a1c_tier, axis=1)

# ---- 5. 队列流程 ----
flow = []
def step(name, mask):
    flow.append({"step": name, "n": int(mask.sum()), "removed": int((~mask).sum() - (flow[-1]["n"] - len(d[mask]) if flow else 0))})
    return mask
allm = pd.Series(True, index=d.index)
step("open cardiac surgery cohort (v2 codebook, first qualifying stay, adult)", allm)
m2 = d["gv"].notna()
step(">=2 eligible glucose measurements in date-anchored 24h window (priority series pending substitution)", m2)
m3 = d["landmark_eligible"]
step("alive with follow-up at day-1 landmark", m2 & m3)
flow = pd.DataFrame(flow)
flow["removed_prev_step"] = flow["n"].diff().fillna(flow["n"].iloc[0]).astype(int) - flow["n"].iloc[0] + flow["n"].iloc[0]
flow = flow[["step","n"]]
flow["removed_from_previous"] = -flow["n"].diff().fillna(0).astype(int)
flow.loc[0,"removed_from_previous"] = 0
flow.to_csv(os.path.join(OUT,"results","cohort_flow_all.csv"), index=False)

# ---- 6. 保存 ----
keep = ["subject_id","hadm_id","stay_id","t0","landmark_time","gender","age_at_admission","bmi",
        "diabetes","charlson_comorbidity_index","charlson_without_diabetes","diabetes_weight","diabetes_icd_with_complication_fixed",
        "procedure_group_main","procedure_cat6","cabg_flag","open_valve_flag","aortic_surgery_flag",
        "transplant_vad_flag","congenital_cardiac_flag",
        "lactate_postop_first","creat_postop_first","wbc_postop_first","albumin_adm_first",
        "hgb_postop_first","platelets_postop_first","sofa_24h",
        "survival_time_days","t_lm_30","t_lm_365","event_lm_30","event_lm_365","landmark_eligible",
        "hba1c_pct","hba1c_days_from_surgery","hba1c_tier","hba1c_window","hba1c_source",
        "gv","shr","glucose_mean_postop_24h","n_glucose_postop_24h","span_hours","density_per_hour",
        "frac_poct","frac_central_lab","frac_blood_gas","frac_icu_charted",
        "dischtime_ts","hosp_los_days","icu_los_days","final_v2_A","final_v2_open_core","final_v2_C",
        "apachepatientresult_flag","hospitalid"]
keep = [c for c in keep if c in d.columns]
d[keep].to_csv(os.path.join(OUT,"data","analysis_base.csv"), index=False)

def sha(p):
    h = hashlib.sha256()
    with open(p,"rb") as f:
        for ch in iter(lambda: f.read(1<<20), b""): h.update(ch)
    return h.hexdigest()
cks = []
for p in [os.path.join(OUT,"data","analysis_base.csv")]:
    cks.append({"file": p, "rows": len(d), "cols": len(keep), "sha256": sha(p)})
pd.DataFrame(cks).to_csv(os.path.join(OUT,"results","dataset_checksums.csv"), index=False)
print("flow:\n", flow.to_string(index=False))
print("charlson_wo_dm 推导完成,交叉验证:", xval)
print("hba1c_tier 分布:\n", d["hba1c_tier"].value_counts().to_string())
print("procedure_cat6 分布:\n", d["procedure_cat6"].value_counts().to_string())
print("landmark_eligible:", int(d["landmark_eligible"].sum()), "/", len(d))
print("PHASE2_DONE")
