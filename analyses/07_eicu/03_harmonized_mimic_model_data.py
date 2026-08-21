# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
"""组装 MIMIC harmonized model dataset(与 eICU 同构)。"""
import pandas as pd, numpy as np, os
BASE = os.path.join(PGV("replication_work"), "02_mimic_harmonized")

frames = {n: pd.read_csv(f"{BASE}/frame_{n}.csv", low_memory=False) for n in ["main","sensA","sensB"]}
gf = pd.read_csv(f"{BASE}/mimic_glucose_features.csv")
cr = pd.read_csv(f"{BASE}/mimic_creatinine.csv")
oc = pd.read_csv(f"{BASE}/mimic_outcomes.csv")

base = frames["sensB"].copy()  # 最宽帧, 再用标记列区分
base = base.drop(columns=[c for c in ["icu_los_days","hosp_los_days"] if c in base.columns])
for name in ["main","sensA","sensB"]:
    base[f"in_{name}"] = base.stay_id.isin(set(frames[name].stay_id))

base["procedure_category"] = np.select(
    [(base.cabg_flag==1)&(base.open_valve_flag==1), base.cabg_flag==1, base.open_valve_flag==1],
    ["combined","cabg","valve"], default=base.procedure_group_main)

m = (base.merge(gf, on="stay_id", how="left")
         .merge(cr, on="hadm_id", how="left")
         .merge(oc[["stay_id","hosp_mortality","icu_mortality","died_within_24h","discharged_within_24h",
                    "in_landmark","post_landmark_hosp_mortality","icu_los_lt_24h","icu_los_days","hosp_los_days"]],
                on="stay_id", how="left"))
m = m.rename(columns={"age_at_admission":"age","gender":"sex","charlson_comorbidity_index":"charlson","sofa_24h":"sofa"})
keep = ["subject_id","hadm_id","stay_id","procedure_category","age","sex","in_main","in_sensA","in_sensB",
        "diabetes","bmi","creatinine","charlson","sofa",
        "glucose_n_raw","glucose_n","glucose_mean_24h","glucose_sd_24h","glucose_cv_24h",
        "glucose_min_24h","glucose_max_24h","glucose_range_24h","glucose_first","glucose_last",
        "measurement_span_minutes","poct_fraction","central_lab_fraction","blood_gas_fraction","icu_charted_fraction",
        "hosp_mortality","icu_mortality","died_within_24h","discharged_within_24h","in_landmark",
        "post_landmark_hosp_mortality","icu_los_lt_24h","icu_los_days","hosp_los_days"]
m = m[keep]
m.to_csv(os.path.join(BASE, "mimic_model_data.csv"), index=False)
for nm in ["main","sensA","sensB"]:
    s = m[m[f"in_{nm}"]]
    lm = s[s.in_landmark]
    print(nm, "n:", len(s), "ge2:", (s.glucose_n>=2).sum(), "landmark:", len(lm),
          "landmark_ge2:", (lm.glucose_n>=2).sum(),
          "landmark_deaths:", int(lm.post_landmark_hosp_mortality.sum()),
          "bmi_avail:", round(s.bmi.notna().mean(),3), "creat_avail:", round(s.creatinine.notna().mean(),3),
          "sofa_avail:", round(s.sofa.notna().mean(),3))
print("done")
