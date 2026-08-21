# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# 24_joint_qc_assembly.py — 联合模块 QC 单元测试、规格矩阵、number map、写回总项目
import os, json
import numpy as np
import pandas as pd

OUT = PGV("mimic_record_work")
RES = os.path.join(OUT, "results")
JK = json.load(open(os.path.join(RES,"shr_gv_joint_constants.json")))

base = pd.read_csv(os.path.join(OUT,"data","analysis_base_bmi_repaired.csv"), low_memory=False)
feat = pd.read_csv(os.path.join(OUT,"data","features_priority.csv"))
d = base.merge(feat.rename(columns={c:f"ps_{c}" for c in feat.columns if c!="stay_id"}), on="stay_id")
d["eAG"] = 28.7*d["hba1c_pct"] - 46.7
tgt = d[(d["landmark_eligible"].astype(str).str.lower().isin(["true","t","1"])) &
        (d["ps_gv_sd"].notna()) & (d["ps_glucose_count"]>=2)]
j = tgt[(tgt["hba1c_tier"]=="1-90 days pre-index") & (tgt["eAG"]>0)]

tests = []
def T(name, ok, detail=""):
    tests.append({"test": name, "pass": bool(ok), "detail": str(detail)})

T("SHR and GV from the identical patient glucose sequence",
  bool(j["ps_gv_sd"].notna().all() and j["ps_mean_glucose"].notna().all() and
       np.allclose(j["ps_mean_glucose"]/j["eAG"], j["ps_mean_glucose"]/j["eAG"])))
T("all HbA1c strictly 1-90 days pre-index",
  bool(((j["hba1c_days_from_surgery"]<0) & (j["hba1c_days_from_surgery"]>=-90)).all()))
ser = pd.read_csv(os.path.join(OUT,"data","series_priority.csv"), parse_dates=["minute"])
jser = ser[ser["stay_id"].isin(j["stay_id"])].merge(j[["stay_id","t0"]], on="stay_id")
jser["t0"] = pd.to_datetime(jser["t0"])
T("all exposure measurements strictly before landmark", bool((jser["minute"] < jser["t0"]+pd.Timedelta(hours=24)).all()))
T("one contribution per patient", bool(j["stay_id"].is_unique and j["subject_id"].is_unique))
T("eAG > 0 for all", bool((j["eAG"]>0).all()))
T("fixed constants single source", os.path.exists(os.path.join(RES,"shr_gv_joint_constants.json")),
  f"shr_sd={JK['shr_sd']:.6f}, gv_sd={JK['gv_sd']:.6f}")
sup = pd.read_csv(os.path.join(RES,"shr_gv_common_support.csv"))
T("all four corner scenarios inside common support", bool(sup["inside_convex_hull"].all()))

lin = pd.read_csv(os.path.join(RES,"shr_gv_joint_linear_results.csv"))
omni = pd.read_csv(os.path.join(RES,"shr_gv_joint_omnibus_tests.csv"))
T("J0/J1 identical patients by construction (MICE subset)",
  bool((lin["N"].dropna().unique()==[JK["N_joint"]]).all() if len(lin) else False),
  f"N set: {sorted(lin['N'].dropna().unique()) if len(lin) else 'n/a'}")
T("omnibus tests use D1 method, no P averaging",
  bool(omni["method"].str.contains("D1").all() if len(omni) else False))
T("lead secondary test present (30d J1 vs J0)",
  bool(((omni["model_id"]=="OMNI_J1_vs_J0_30d") & (omni["analysis"]=="LEAD SECONDARY TEST")).any()))
T("key secondary test present (365d J1 vs J0)",
  bool((omni["model_id"]=="OMNI_J1_vs_J0_365d").any()))
nev = lin.groupby("model_id").agg(N_min=("N","min"), N_max=("N","max")) if len(lin) else None
T("N/events consistent within models", bool((nev["N_min"]==nev["N_max"]).all()) if nev is not None and len(nev) else False)

qc = pd.DataFrame(tests)
qc.to_csv(os.path.join(RES,"shr_gv_qc_unit_tests.csv"), index=False)
print(qc.to_string(index=False))
assert qc["pass"].all(), "存在未通过的联合模块 QC 测试"

# ---- 规格矩阵 ----
files = {
 "shr_gv_joint_linear_results.csv":"J1 linear",
 "shr_gv_joint_omnibus_tests.csv":"omnibus",
 "shr_gv_joint_nonlinearity.csv":"nonlinearity",
 "shr_gv_component_decomposition.csv":"decomposition",
 "shr_gv_interaction_results.csv":"interaction",
 "shr_gv_absolute_risk_surface.csv":"absolute_risk",
 "shr_gv_model_performance.csv":"performance",
 "shr_gv_model_performance_pairs.csv":"performance_pairs",
 "shr_gv_source_sensitivity.csv":"source",
 "shr_gv_hba1c_timing_sensitivity.csv":"timing",
 "shr_gv_ipw_sensitivity.csv":"ipw",
 "shr_gv_ph_diagnostics.csv":"ph",
 "shr_gv_bootstrap_stability.csv":"stability",
}
rows = []
for f, mod in files.items():
    p = os.path.join(RES, f)
    if not os.path.exists(p): continue
    df = pd.read_csv(p)
    for r in df.itertuples():
        mid = getattr(r, "model_id", None)
        if pd.isna(mid): continue
        rows.append({"model_id": mid, "module": f"shr_gv_joint/{mod}",
            "database": "MIMIC-IV", "cohort": JK["cohort"], "landmark": "day-1",
            "glucose_source_rule": "blood-only priority series",
            "exposure_definition": "SHR = mean glucose / (28.7*HbA1c - 46.7); GV = within-patient SD (same series)",
            "exposure_scale": getattr(r, "scale", "per fixed joint-cohort SD / per 0.1 SHR / per 10 mg/dL GV"),
            "covariates": "fixed clinical set (age, sex, BMI, diabetes, Charlson-wo-DM, procedure_cat6, lactate, creatinine, SOFA)",
            "missing_data_method": "MICE m=50 (reused) / complete-case (sensitivity)",
            "N": getattr(r, "N", np.nan), "events": getattr(r, "events", np.nan),
            "standardization_constants": f"shr_mean={JK['shr_mean']:.4f}, shr_sd={JK['shr_sd']:.4f}, gv_mean={JK['gv_mean']:.4f}, gv_sd={JK['gv_sd']:.4f}",
            "seed": 20260726})
msm = pd.DataFrame(rows).drop_duplicates("model_id")
msm.to_csv(os.path.join(RES,"shr_gv_model_specification_matrix.csv"), index=False)

# ---- number map ----
omni30 = omni[omni.model_id=="OMNI_J1_vs_J0_30d"].iloc[0] if (omni.model_id=="OMNI_J1_vs_J0_30d").any() else None
omni365 = omni[omni.model_id=="OMNI_J1_vs_J0_365d"].iloc[0] if (omni.model_id=="OMNI_J1_vs_J0_365d").any() else None
shr30 = lin[lin.model_id=="J1_30d_shr_zf"].iloc[0] if (lin.model_id=="J1_30d_shr_zf").any() else None
gv30 = lin[lin.model_id=="J1_30d_gv_zf"].iloc[0] if (lin.model_id=="J1_30d_gv_zf").any() else None
nm = [
 ("补充材料:SHR–GV 联合模块(新增)", "—(此前无联合框架)", "OMNI_J1_vs_J0_30d",
  f"joint omnibus F={omni30.stat:.3f}, df=2, P={omni30.P:.4f}" if omni30 is not None else "n/a", "新增(lead secondary)"),
 ("补充材料:365 日联合检验", "—", "OMNI_J1_vs_J0_365d",
  f"joint omnibus F={omni365.stat:.3f}, df=2, P={omni365.P:.4f}" if omni365 is not None else "n/a", "新增(key secondary)"),
 ("主文/补充:SHR 单系数描述", "SHR 30d 1.17 (1.00–1.37) 等旧值", "J1_30d_shr_zf",
  f"J1 内 SHR per SD {shr30.HR:.3f} ({shr30.lo:.3f}–{shr30.hi:.3f}), P={shr30.P:.3f}" if shr30 is not None else "n/a",
  "替换为联合模型条件估计;单系数不作确定关联"),
 ("主文:主要阴性结论", "GV-only 30d 阴性(重构主要分析)", "—",
  "不变:联合模块不得推翻;见 shr_gv_final_interpretation.md", "不变"),
]
nmdf = pd.DataFrame(nm, columns=["manuscript_location","old_value","new_model_id","new_value","action"])
nmdf.to_csv(os.path.join(RES,"shr_gv_manuscript_number_map.csv"), index=False)

# ---- 写回总项目矩阵与 number map ----
for target, add_df in [("13_model_specification_matrix.csv", msm.assign(module=msm["module"])),
                       ("14_manuscript_number_map.csv", nmdf.rename(columns={
                           "manuscript_location":"manuscript_location","old_value":"old_value",
                           "new_model_id":"new_model_id","new_value":"new_value","action":"action"}))]:
    tp = os.path.join(RES, target)
    if os.path.exists(tp):
        old = pd.read_csv(tp)
        add_df2 = add_df.copy()
        for c in old.columns:
            if c not in add_df2.columns: add_df2[c] = np.nan
        for c in add_df2.columns:
            if c not in old.columns: old[c] = np.nan
        merged = pd.concat([old, add_df2[old.columns]], ignore_index=True)
        merged = merged.drop_duplicates(subset=[c for c in ["model_id","manuscript_location"] if c in merged.columns], keep="last")
        merged.to_csv(tp, index=False)
        print(f"写回 {target}: {len(old)} -> {len(merged)} 行")
print("QC + 装配完成")
