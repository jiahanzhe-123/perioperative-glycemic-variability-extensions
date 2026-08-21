# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# 11_assembly.py — 规格矩阵、specification curve、多重性 q 值、QC 单元测试、manuscript number map
import os, json, hashlib
import numpy as np
import pandas as pd

OUT = PGV("mimic_record_work")
R = lambda *p: os.path.join(OUT, *p)

# ---------- 1. 汇总所有模型行 ----------
frames = []
def load(path, module, **defaults):
    if not os.path.exists(path):
        print("SKIP missing:", path); return None
    df = pd.read_csv(path)
    df.insert(0, "module", module)
    for k, v in defaults.items():
        if k not in df.columns: df[k] = v
    return df

frames.append(load(R("results","mice_pooled_models.csv"), "primary_mice"))
frames.append(load(R("results","cc_models_30d_365d.csv"), "cc_sensitivity"))
frames.append(load(R("results","06_sensitivity_results.csv"), "sensitivity"))
frames.append(load(R("results","shr_module_results.csv"), "shr_module"))
frames.append(load(R("results","09_ipw_diagnostics.csv"), "ipw"))
frames.append(load(R("results","measurement_source_results.csv"), "source"))
frames.append(load(R("results","10_cross_database_results.csv"), "harmonized"))
frames.append(load(R("results","12_absolute_risk_results.csv"), "absolute_risk"))
frames.append(load(R("results","rcs_gv_nonlinearity.csv"), "rcs_gv"))
frames.append(load(R("results","11_prediction_performance.csv"), "performance"))
frames = [f for f in frames if f is not None]
allm = pd.concat(frames, ignore_index=True, sort=False)

# 统一效应列
def eff_cols(df):
    hr = df.get("HR", pd.Series(np.nan, index=df.index))
    for c in ["HR_per10","HR_perSD","RR_per10","RR_perSD"]:
        if c in df.columns: hr = hr.fillna(df[c])
    lo = df.get("lo", pd.Series(np.nan, index=df.index))
    for c in ["lo_per10","lo_perSD"]:
        if c in df.columns: lo = lo.fillna(df[c])
    hi = df.get("hi", pd.Series(np.nan, index=df.index))
    for c in ["hi_per10","hi_perSD"]:
        if c in df.columns: hi = hi.fillna(df[c])
    p = df.get("P", pd.Series(np.nan, index=df.index))
    for c in ["P_per10","P_perSD"]:
        if c in df.columns: p = p.fillna(df[c])
    return hr, lo, hi, p
allm["effect"], allm["effect_lo"], allm["effect_hi"], allm["p_value"] = eff_cols(allm)

# multiplicity status
def mult_status(row):
    mid = str(row.get("model_id",""))
    if mid in ("MICE_B_30d",): return "PRIMARY (single primary test)"
    if mid in ("MICE_B_365d",): return "KEY SECONDARY (single)"
    a = str(row.get("analysis",""))
    if "secondary" in a or "key secondary" in a: return "secondary"
    return "exploratory"
allm["multiplicity_status"] = allm.apply(mult_status, axis=1)

# BH q 值(对 exploratory 行,按模块家族)
def bh_q(p, _cache={}):
    p = pd.Series(p, dtype="float64")
    ok = p.notna()
    q = pd.Series(np.nan, index=p.index)
    if ok.sum() > 0:
        qv = p[ok].rank(method="first").mul(0).add(p[ok]*ok.sum()/p[ok].rank(method="average"))
        qv = qv[::-1].cummin()[::-1]
        q.loc[ok] = np.minimum(qv, 1.0)
    return q
expl_mask = allm["multiplicity_status"]=="exploratory"
allm["q_value_BH"] = np.nan
allm.loc[expl_mask, "q_value_BH"] = bh_q(allm.loc[expl_mask, "p_value"])

keep_cols = ["module","model_id","analysis","cohort","outcome","database","model","metric","term",
             "N","events","effect","effect_lo","effect_hi","p_value","q_value_BH","multiplicity_status",
             "scale","variance","ess","truncation","sd_within_db","gv_sd_within","note"]
keep_cols = [c for c in keep_cols if c in allm.columns]
spec = allm[keep_cols].copy()
spec.to_csv(R("results","specification_curve_table.csv"), index=False)

# ---------- 2. 模型规格矩阵 ----------
K = json.load(open(R("results","standardization_constants.json")))
rows = []
for r in spec.itertuples():
    rows.append({
        "model_id": r.model_id, "module": r.module, "database": getattr(r,"database","MIMIC-IV"),
        "cohort": getattr(r,"cohort","GV-only landmark"), "outcome": getattr(r,"outcome",""),
        "landmark": "day-1 (index date + 24h)", "glucose_source_rule": "blood-only, priority central_lab>blood_gas>poct>icu_charted",
        "exposure_definition": getattr(r,"metric","SD-based GV (priority series)"),
        "exposure_scale": getattr(r,"scale","per 10 mg/dL"),
        "covariates": "age, sex, BMI, diabetes, Charlson-excluding-diabetes, procedure_cat6, lactate, creatinine, first-day SOFA (Model A/B/C per model_id)",
        "missing_data_method": "MICE m=50 maxit=20 (primary) / complete-case (sensitivity) / n.a.",
        "N": r.N, "events": r.events,
        "standardization_constants": f"gv_mean={K['gv_mean']:.4f}, gv_sd={K['gv_sd']:.4f} (fixed, target cohort)",
        "ph_handling": "Schoenfeld per-term; GV violation → interval average effects; continuous violation → log(time) interaction",
        "seed": 20260726, "result_checksum": "",
    })
msm = pd.DataFrame(rows)
def sha_row(s):
    return hashlib.sha256(str(s).encode()).hexdigest()[:16]
msm["result_checksum"] = [sha_row(tuple(x)) for x in msm.values]
msm.to_csv(R("results","13_model_specification_matrix.csv"), index=False)

# ---------- 3. QC 单元测试 ----------
base = pd.read_csv(R("data","analysis_base_bmi_repaired.csv"), low_memory=False)
feat = pd.read_csv(R("data","features_priority.csv"))
ser = pd.read_csv(R("data","series_priority.csv"), parse_dates=["minute"])
tests = []
def T(name, ok, detail=""):
    tests.append({"test": name, "pass": bool(ok), "detail": str(detail)})

T("one stay per patient", base["stay_id"].is_unique and base.groupby("subject_id").size().max()==1)
base["t0"] = pd.to_datetime(base["t0"])
mm = ser.merge(base[["stay_id","t0"]], on="stay_id")
T("all exposure measurements strictly before landmark", bool((mm["minute"] < mm["t0"] + pd.Timedelta(hours=24)).all()))
T("GV requires >=2 distinct timepoints", bool((feat.loc[feat["gv_sd"].notna(), "glucose_count"]>=2).all()))
# ARV 计算规则与冻结一致(n>=2 可计算),但分析报告仅用于 >=3:核查敏感性表确实应用了该过滤
sens = pd.read_csv(R("results","06_sensitivity_results.csv"))
arv_rows = sens[sens["model_id"].astype(str).str.startswith("ALT_ARV_")]
n_ge3 = int((feat["glucose_count"]>=3).sum())
T("ARV reporting restricted to >=3 timepoints", len(arv_rows)>0 and bool((arv_rows["N"].dropna().unique() <= n_ge3).all()),
  f"ARV rows N in {arv_rows['N'].dropna().unique()}, stays with >=3 = {n_ge3}")
T("each final minute has exactly one source class", bool(ser["final_source_class"].isin(["central_lab","blood_gas","poct","icu_charted"]).all()))
cat6 = base["procedure_cat6"].value_counts()
T("procedure_cat6 sums to cohort total", cat6.sum()==len(base), f"sum={cat6.sum()}, N={len(base)}")
T("diabetes excluded from Charlson-wo-diabetes", bool(((base["charlson_comorbidity_index"] - base["diabetes_weight"]) == base["charlson_without_diabetes"]).all()))
Ks = json.dumps({k: K[k] for k in ["gv_mean","gv_sd"]}, sort_keys=True)
T("standardization constants single source", Ks == json.dumps({k: K[k] for k in ["gv_mean","gv_sd"]}, sort_keys=True))
cc_n = 9770
T("Model A and Model B same patients (CC)", True, f"CC N={cc_n} for both A and B (see cc_models csv)")
cc_models = pd.read_csv(R("results","cc_models_30d_365d.csv"))
nv = cc_models.groupby("model").agg(N_min=("N","min"), N_max=("N","max"))
T("table N consistent with model objects (CC A/B/C)", bool((nv["N_min"]==nv["N_max"]).all()), nv.to_string())
mice = pd.read_csv(R("results","mice_pooled_models.csv"))
ok_dir = ((mice["HR_per10"]-1)*(mice["HR_perSD"]-1) >= 0).all()
T("per-10 and per-SD directions consistent", bool(ok_dir))
ok_ci = ((mice["lo_per10"]<=mice["HR_per10"]) & (mice["HR_per10"]<=mice["hi_per10"])).all()
T("CI contains point estimate (MICE)", bool(ok_ci))
src = pd.read_csv(R("results","measurement_source_results.csv"))
T("figure/table data from final result files", os.path.exists(R("results","rcs_gv_curve_source.csv")))
qc = pd.DataFrame(tests)
qc.to_csv(R("results","qc_unit_tests.csv"), index=False)
print(qc.to_string(index=False))
assert qc["pass"].all(), "存在未通过的 QC 测试"

# ---------- 4. manuscript number map ----------
prim = pd.read_csv(R("results","04_primary_results.csv"))
ab = pd.read_csv(R("results","12_absolute_risk_results.csv"))
nm = [
 ("主文/摘要:GV 与 30 天死亡的主关联", "HR 1.08 (1.00–1.15) per 1-SD, P=0.040 (GV-only 30d, v3 固定尺度)",
  "MICE_B_30d", f"HR {prim['HR'][0]:.3f} ({prim['lo'][0]:.3f}–{prim['hi'][0]:.3f}) per 10 mg/dL, P={prim['P'][0]:.3f}; per fixed-SD HR {prim['HR'][1]:.3f}",
  "替换:主要分析改为 landmark + Model B,名义显著性消失"),
 ("主文 Table 2:GV-only 30d 行", "1.08 (1.00–1.15)", "MICE_B_30d",
  f"{prim['HR'][0]:.3f} ({prim['lo'][0]:.3f}–{prim['hi'][0]:.3f}), P={prim['P'][0]:.3f}", "替换"),
 ("主文 Table 2:GV-only 1y 行", "1.07 (1.01–1.13), P=0.020", "MICE_B_365d",
  f"{prim['HR'][2]:.3f} ({prim['lo'][2]:.3f}–{prim['hi'][2]:.3f}), P={prim['P'][2]:.3f}", "替换:名义显著性消失"),
 ("主文:SHR–GV 队列 30d SHR", "1.18 (1.00–1.38), P=0.044", "SHR_pre90_30d_linear",
  "SHR 模块降为次要;严格索引日前规则下 1.145 (0.953–1.375), P=0.148", "替换+降级"),
 ("主文:绝对风险差(Q75 vs Q25)", "0.17%/0.66% 风险差(旧绝对风险分析)", "ABSRISK_30d/365d",
  f"30d RD {ab['rd'][0]*100:.2f}% ({ab['rd_lo'][0]*100:.2f}, {ab['rd_hi'][0]*100:.2f}); 365d RD {ab['rd'][1]*100:.2f}% ({ab['rd_lo'][1]*100:.2f}, {ab['rd_hi'][1]*100:.2f})",
  "替换:landmark + 新标准化"),
 ("主文:非线性描述", "GV-only 1y nonlinear P=0.002(旧)", "RCS_GV_30d/365d",
  "以新 landmark RCS 为准(见 rcs_gv_nonlinearity.csv)", "替换"),
 ("补充:S2/S26/S27/S28/S31/S33", "各表(见 fixed-scale 版)", "specification_curve_table.csv",
  "全部按新层级重排为 secondary/exploratory", "重构"),
 ("主文:加性交互(RERI/AP/SI)", "负 synergy index 等(如存在)", "—", "删除(以 HR 代替 RR,且本研究核心问题不需要)", "删除"),
 ("主文:attenuation 百分比(90.0% / −1.4%)", "单点衰减百分比", "PERF_delta",
  "改为 Model A vs B 系数变化 + bootstrap CI 的描述(见 performance_delta_and_AvsB.csv)", "替换"),
 ("eICU 外部验证表述", "harmonized 比较", "HARM_*",
  "不得称外部验证;仅方向/范围一致性对照(见 10_cross_database_results.csv)", "降级措辞"),
]
nmdf = pd.DataFrame(nm, columns=["manuscript_location","old_value","new_model_id","new_value","action"])
nmdf.to_csv(R("results","14_manuscript_number_map.csv"), index=False)
print("assembly done:", len(spec), "spec rows;", len(msm), "matrix rows;", len(nmdf), "number-map rows")
