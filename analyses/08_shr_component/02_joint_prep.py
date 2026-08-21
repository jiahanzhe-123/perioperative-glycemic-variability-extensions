# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# 20_joint_prep.py — SHR-GV 联合模块:队列、常数、分布、共线性、二维共同支持
import os, json
import numpy as np
import pandas as pd
from scipy import stats

OUT = PGV("mimic_record_work")
RES = os.path.join(OUT, "results")
FIG = os.path.join(OUT, "figures")
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

base = pd.read_csv(os.path.join(OUT,"data","analysis_base_bmi_repaired.csv"), low_memory=False)
feat = pd.read_csv(os.path.join(OUT,"data","features_priority.csv"))
d = base.merge(feat.rename(columns={c: f"ps_{c}" for c in feat.columns if c!="stay_id"}), on="stay_id")
d["gv"] = d["ps_gv_sd"]; d["mean_glu"] = d["ps_mean_glucose"]; d["glucose_count"] = d["ps_glucose_count"]
d["eAG"] = 28.7*d["hba1c_pct"] - 46.7
d["shr"] = d["mean_glu"] / d["eAG"]

# ---- 联合队列(严格 1–90 日术前 HbA1c) ----
tgt = d[(d["landmark_eligible"].astype(str).str.lower().isin(["true","t","1"])) &
        (d["gv"].notna()) & (d["glucose_count"]>=2)]
flow = [
    ("open cardiac surgery cohort (rebuild)", len(d)),
    ("landmark-eligible + >=2 eligible glucose (GV-only target)", len(tgt)),
    ("HbA1c strictly 1-90 days pre-index", None),
    ("eAG > 0 (HbA1c QC range)", None),
]
j = tgt[tgt["hba1c_tier"]=="1-90 days pre-index"].copy()
flow[2] = ("HbA1c strictly 1-90 days pre-index", len(j))
j = j[j["eAG"] > 0]
flow[3] = ("eAG > 0 (HbA1c QC range)", len(j))
assert j["stay_id"].is_unique and j["subject_id"].is_unique
assert (j["hba1c_days_from_surgery"] < 0).all() and (j["hba1c_days_from_surgery"] >= -90).all()
assert (j["eAG"] > 0).all()
flow = pd.DataFrame(flow, columns=["step","n"])
flow["removed_from_previous"] = -flow["n"].diff().fillna(0).astype(int)
flow.to_csv(os.path.join(RES,"shr_gv_joint_cohort_flow.csv"), index=False)

# ---- 固定常数(联合队列一次计算) ----
K = json.load(open(os.path.join(RES,"standardization_constants.json")))
JK = dict(
  shr_mean=float(j["shr"].mean()), shr_sd=float(j["shr"].std()), shr_median=float(j["shr"].median()),
  shr_q25=float(j["shr"].quantile(.25)), shr_q75=float(j["shr"].quantile(.75)),
  gv_mean=float(j["gv"].mean()), gv_sd=float(j["gv"].std()), gv_median=float(j["gv"].median()),
  gv_q25=float(j["gv"].quantile(.25)), gv_q75=float(j["gv"].quantile(.75)),
  shr_center=float(j["shr"].median()), gv_center=float(j["gv"].median()),
  shr_knots3=[float(x) for x in j["shr"].quantile([.10,.50,.90])],
  gv_knots3=[float(x) for x in j["gv"].quantile([.10,.50,.90])],
  mean_knots3=[float(x) for x in j["mean_glu"].quantile([.10,.50,.90])],
  mean_knots4=[float(x) for x in j["mean_glu"].quantile([.05,.35,.65,.95])],
  a1c_knots3=[float(x) for x in j["hba1c_pct"].quantile([.10,.50,.90])],
  N_joint=int(len(j)), events_30d=int(j["event_lm_30"].sum()), events_365d=int(j["event_lm_365"].sum()),
  cohort="strictly pre-index 1-90d HbA1c, day-1 landmark survivors, >=2 eligible measurements",
  seed=20260726)
json.dump(JK, open(os.path.join(RES,"shr_gv_joint_constants.json"),"w"), indent=2)

# ---- 分布表 ----
def dist(v, name):
    qs = [0.01,0.025,0.05,0.25,0.50,0.75,0.95,0.975,0.99]
    o = {"variable":name, "n":int(v.notna().sum()), "mean":v.mean(), "sd":v.std(),
         "median":v.median(), "iqr":v.quantile(.75)-v.quantile(.25),
         "min":v.min(), "max":v.max(), "skew":stats.skew(v.dropna())}
    for q in qs: o[f"p{int(q*1000)/10}".replace(".0","").replace(".","_")] = v.quantile(q)
    return o
dist_rows = [dist(j[v], n) for v,n in [("shr","SHR"),("gv","GV (SD)"),("mean_glu","mean glucose"),
    ("hba1c_pct","HbA1c"),("eAG","eAG"),("glucose_count","glucose measurement count"),("ps_span_hours","measurement span (h)")]]
dist_tab = pd.DataFrame(dist_rows)
dist_tab.to_csv(os.path.join(RES,"shr_gv_exposure_distribution.csv"), index=False)
outliers = j.nlargest(15, "shr")[["stay_id","shr","gv","mean_glu","hba1c_pct","eAG","event_lm_30","event_lm_365"]]
outliers["reason"] = "top-15 SHR"
outliers.to_csv(os.path.join(RES,"shr_gv_outlier_listing.csv"), index=False)

# ---- 相关与共线性 ----
jj = j.copy()
jj["shr_z"] = (jj["shr"]-JK["shr_mean"])/JK["shr_sd"]; jj["gv_z"] = (jj["gv"]-JK["gv_mean"])/JK["gv_sd"]
jj["mean_z"] = (jj["mean_glu"]-jj["mean_glu"].mean())/jj["mean_glu"].std()
jj["a1c_z"] = (jj["hba1c_pct"]-jj["hba1c_pct"].mean())/jj["hba1c_pct"].std()
pairs = [("shr","gv"),("shr","mean_glu"),("shr","hba1c_pct"),("gv","mean_glu"),("gv","glucose_count"),("mean_glu","hba1c_pct")]
cor_rows = []
for a,b in pairs:
    cor_rows.append({"pair":f"{a} vs {b}",
        "pearson":jj[a].corr(jj[b]), "spearman":jj[a].corr(jj[b], method="spearman")})
# partial correlation(SHR~GV | covariates)用线性残差
import numpy.linalg as la
def partial(a, b, controls):
    X = np.column_stack([np.ones(len(jj))] + [jj[c].fillna(jj[c].mean()).values for c in controls])
    ra = jj[a].values - X @ la.lstsq(X, jj[a].values, rcond=None)[0]
    rb = jj[b].values - X @ la.lstsq(X, jj[b].values, rcond=None)[0]
    return np.corrcoef(ra, rb)[0,1]
ctrls = ["age_at_admission","bmi","charlson_without_diabetes","sofa_24h","lactate_postop_first","creat_postop_first"]
for a,b in [("shr","gv"),("shr","mean_glu"),("gv","mean_glu")]:
    cor_rows.append({"pair":f"{a} vs {b} (partial | clinical)", "pearson":partial(a,b,ctrls), "spearman":np.nan})
cor_tab = pd.DataFrame(cor_rows)

def vif_of(y, xs, dd):
    dd2 = dd[[y]+xs].dropna()
    X = np.column_stack([np.ones(len(dd2))] + [dd2[x].values for x in xs])
    b = la.lstsq(X, dd2[y].values, rcond=None)[0]
    r2 = 1 - ((dd2[y].values - X@b)**2).sum()/((dd2[y].values - dd2[y].mean())**2).sum()
    return 1/(1-r2)
numvars = ["shr_z","gv_z","mean_z","a1c_z","age_at_admission","bmi","charlson_without_diabetes","sofa_24h","lactate_postop_first","creat_postop_first"]
vif_rows = []
for v in ["shr_z","gv_z","mean_z","a1c_z"]:
    xs = [x for x in numvars if x != v]
    vif_rows.append({"term": v, "vif": vif_of(v, xs, jj)})
X = jj[numvars].dropna().values
X = (X - X.mean(0))/X.std(0)
sv = la.svd(X, compute_uv=False)
ci = sv.max()/sv
vdp = (ci**2)[:,None]**-1 / ((ci**2)**-1).sum()
vif_tab = pd.DataFrame(vif_rows)
vif_tab["condition_index_max"] = ci.max()
vif_tab["condition_index_min"] = ci.min()
cor_tab = pd.concat([cor_tab, pd.DataFrame([{"pair":"__VIF__","pearson":np.nan,"spearman":np.nan}])], ignore_index=True)
cor_tab.to_csv(os.path.join(RES,"shr_gv_correlation_collinearity.csv"), index=False)
vif_tab.to_csv(os.path.join(RES,"shr_gv_vif_condition.csv"), index=False)

# ---- 二维共同支持 ----
from scipy.spatial import ConvexHull
pts = jj[["shr","gv"]].dropna().values
hull = ConvexHull(pts)
inside = np.zeros(len(pts), dtype=bool)
# 逐点判 hull 包含(含边界):用 matplotlib path
from matplotlib.path import Path
hp = Path(pts[hull.vertices])
inside = hp.contains_points(pts, radius=1e-9)
jj_valid = jj.loc[jj[["shr","gv"]].notna().all(axis=1)].copy()
jj_valid["in_hull"] = inside
frac_in = float(np.mean(inside))
# 四角情景支持核查
corners = [("SHR P25 / GV P25", JK["shr_q25"], JK["gv_q25"]), ("SHR P75 / GV P25", JK["shr_q75"], JK["gv_q25"]),
           ("SHR P25 / GV P75", JK["shr_q25"], JK["gv_q75"]), ("SHR P75 / GV P75", JK["shr_q75"], JK["gv_q75"])]
sup_rows = []
for name, s_, g_ in corners:
    ok = hp.contains_points(np.array([[s_, g_]]), radius=1e-9)[0]
    # 近邻密度:最近 5% 距离内患者比例
    dd2 = np.sqrt(((pts - np.array([s_,g_])) / pts.std(0))**2).sum(1)
    sup_rows.append({"scenario":name, "shr":s_, "gv":g_, "inside_convex_hull":bool(ok),
                     "patients_within_0.5SD_radius":int((dd2<0.5).sum()), "pct":float((dd2<0.5).mean())})
sup_tab = pd.DataFrame(sup_rows)
sup_tab.to_csv(os.path.join(RES,"shr_gv_common_support.csv"), index=False)
# 网格支持(事件计数)
jj_valid["shr_bin"] = pd.qcut(jj_valid["shr"], 8, duplicates="drop")
jj_valid["gv_bin"] = pd.qcut(jj_valid["gv"], 8, duplicates="drop")
grid = jj_valid.groupby(["shr_bin","gv_bin"], observed=True).agg(
    n=("event_lm_30","size"), ev30=("event_lm_30","sum"), ev365=("event_lm_365","sum")).reset_index()
grid.to_csv(os.path.join(RES,"shr_gv_support_grid_counts.csv"), index=False)

# ---- 图 ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
for v, name in [("shr","SHR"),("gv","GV (SD, mg/dL)"),("mean_glu","mean glucose (mg/dL)"),("hba1c_pct","HbA1c (%)")]:
    fig, ax = plt.subplots(1, 2, figsize=(9,3.2))
    ax[0].hist(jj[v].dropna(), bins=50, color="#0072B2", alpha=0.8); ax[0].set_title(f"{name} histogram")
    jj[v].dropna().plot.kde(ax=ax[1], color="#0072B2"); ax[1].set_title(f"{name} density")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, f"joint_dist_{v}.png"), dpi=600)
    fig.savefig(os.path.join(FIG, f"joint_dist_{v}.pdf")); plt.close(fig)
fig, ax = plt.subplots(figsize=(6.5,5.5))
hb = ax.hexbin(jj_valid["shr"], jj_valid["gv"], gridsize=40, cmap="Blues", mincnt=1)
ax.plot(pts[hull.vertices,0], pts[hull.vertices,1], "r-", lw=1)
for name, s_, g_ in corners:
    ax.scatter([s_],[g_], marker="D", s=40, zorder=5, label=name)
ax.set_xlabel("SHR"); ax.set_ylabel("GV (SD, mg/dL)")
ax.legend(fontsize=7, loc="upper right")
ax.set_title("SHR–GV 2D density, convex hull, prespecified corners")
fig.colorbar(hb, label="patients")
fig.tight_layout(); fig.savefig(os.path.join(FIG,"joint_2d_support.png"), dpi=600)
fig.savefig(os.path.join(FIG,"joint_2d_support.pdf")); plt.close(fig)

print("joint N =", len(j), "; events 30d =", int(j["event_lm_30"].sum()), "; 365d =", int(j["event_lm_365"].sum()))
print("SHR: mean %.3f sd %.3f median %.3f q25 %.3f q75 %.3f" % (JK["shr_mean"],JK["shr_sd"],JK["shr_median"],JK["shr_q25"],JK["shr_q75"]))
print("GV: mean %.2f sd %.2f median %.2f q25 %.2f q75 %.2f" % (JK["gv_mean"],JK["gv_sd"],JK["gv_median"],JK["gv_q25"],JK["gv_q75"]))
print("corr(SHR,GV) pearson %.3f spearman %.3f" % (cor_tab.loc[0,"pearson"], cor_tab.loc[0,"spearman"]))
print("hull coverage: %.4f" % frac_in)
print(sup_tab.to_string(index=False))
print("PHASE20_DONE")
