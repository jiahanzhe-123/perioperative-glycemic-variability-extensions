# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV, PGV_OUT  # noqa: E402
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FR35 = Path(PGV("bmi_repair_work")) / "rerun_root" / "results"
REC  = Path(PGV("mimic_record_work")) / "results"
ph = pd.read_csv(FR35 / "07_ph_diagnostics.csv")
_iv_path = REC / "PH365_INTERVAL_MICE_POOLED.csv"
if not _iv_path.exists():  # 回退:phase1 冻结参考(只读)
    _iv_path = Path(PGV("stats_fix_work")) / "results" / "PH365_INTERVAL_MICE_POOLED.csv"
iv = pd.read_csv(_iv_path)

terms = [
    ("gv10", "GV per 0.555 mmol/L"),
    ("rms::rcs(mean_glu, c(113, 129.7273, 143.2308, 187))", "Mean glucose (RCS)"),
    ("age_at_admission", "Age"),
    ("gender", "Recorded sex"),
    ("bmi", "BMI"),
    ("diabetes", "Diabetes"),
    ("charlson_without_diabetes", "Charlson excl. diabetes"),
    ("procedure_cat6", "Procedure category"),
    ("lactate_postop_first", "Lactate"),
    ("creat_postop_first", "Creatinine"),
    ("sofa_24h", "First-day SOFA"),
    ("GLOBAL", "GLOBAL"),
]
p30 = ph[ph.model_id == "ModelB_30d"].set_index("term")["p"]
p365 = ph[ph.model_id == "ModelB_365d"].set_index("term")["p"]

C30, C365 = "#1f77b4", "#45a3ad"
fig, (axa, axb) = plt.subplots(1, 2, figsize=(13.5, 5.9), gridspec_kw={"width_ratios": [1.15, 1]})

y = np.arange(len(terms))[::-1]
labels = [t[1] for t in terms]
vals30 = [p30[t[0]] for t in terms]
vals365 = [p365[t[0]] for t in terms]
axa.scatter(vals30, y, s=90, color=C30, zorder=3, label="30 days")
axa.scatter(vals365, y, s=90, color=C365, zorder=3, label="365 days")
axa.axvline(0.05, color="#c44e52", ls="--", lw=1.2)
axa.set_yticks(y); axa.set_yticklabels(labels)
axa.set_xlim(-0.05, 1.05)
axa.set_xlabel("Schoenfeld residual P value (complete-case Model B)")
axa.legend(loc="upper left", frameon=False)
axa.set_title("a", loc="left", fontweight="bold", fontsize=14)

iv = iv.set_index("interval")
rows = [("d1-7", "365 d: post-landmark days 1–7 (162 events)"),
        ("d8-30", "365 d: post-landmark days 8–30 (139 events)"),
        ("d31-365", "365 d: post-landmark days 31–365 (444 events)")]
yb = np.arange(len(rows))[::-1]
plab = {"d1-7": "0.057", "d8-30": "0.569", "d31-365": "0.037"}
for yi, (key, lab) in zip(yb, rows):
    r = iv.loc[key]
    axb.errorbar(r.HR_per10, yi, xerr=[[r.HR_per10 - r.lo], [r.hi - r.HR_per10]],
                 fmt="o", color=C365, ecolor=C365, elinewidth=1.6, capsize=4, ms=9, zorder=3)
    axb.text(1.245, yi, f"HR {r.HR_per10:.3f} ({r.lo:.3f}–{r.hi:.3f}); P = {plab[key]}",
             va="center", fontsize=10)
axb.axvline(1.0, color="grey", ls="--", lw=1.2)
axb.set_yticks(yb); axb.set_yticklabels([r[1] for r in rows])
axb.set_xlim(0.85, 1.15)
axb.set_xlabel("GV interval-specific HR per 0.555 mmol/L (MICE-pooled, 95% CI)")
axb.set_title("b", loc="left", fontweight="bold", fontsize=14)

fig.tight_layout()
out = Path(PGV_OUT("results/figures"))
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "figS5_ph_audit.png", dpi=300, bbox_inches="tight")
print("saved", out / "figS5_ph_audit.png")
