# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV  # noqa: E402
#!/usr/bin/env python3
# build_bloodonly_features.py
# 由 blood-only 分钟级血糖序列计算全部 stay 级特征。
# 规则(final statistical freeze 冻结口径):
#   输入序列为已完成的 blood-only 分钟序列(上游 SQL/抽取负责:仅血液来源、
#   [surgical_date_index, +24h) 右端点不含、20-1500 mg/dL、完全重复剔除、同分钟中位数)。
#   本脚本不重做这些上游规则,只做 stay 级聚合,且可对照 final_glucose_features.csv 逐值验证。
# 用法: python build_bloodonly_features.py <minute_series.csv> <out_features.csv>
import sys
import numpy as np
import pandas as pd

THRESH = {"lt70": lambda v: v < 70, "lt54": lambda v: v < 54,
          "gt180": lambda v: v > 180, "gt250": lambda v: v > 250}

def per_stay(g: pd.DataFrame) -> dict:
    g = g.sort_values("minute")
    v = g["value"].to_numpy(dtype=float)
    n = len(v)
    minutes = g["minute"].values.astype("datetime64[s]").astype(np.int64) / 60.0
    out = {"glucose_count": n, "mean_glucose": v.mean()}
    if n >= 2:
        sd = v.std(ddof=1)
        out.update(gv_sd=sd, min_glucose=v.min(), max_glucose=v.max(),
                   range_glucose=v.max() - v.min(), cv=sd / v.mean(),
                   arv=np.abs(np.diff(v)).mean(),
                   tw_arv=(np.abs(np.diff(v)) * np.diff(minutes)).sum() / (minutes[-1] - minutes[0]) if minutes[-1] > minutes[0] else np.nan,
                   span_hours=(minutes[-1] - minutes[0]) / 60.0,
                   density_per_hour=n / ((minutes[-1] - minutes[0]) / 60.0) if minutes[-1] > minutes[0] else np.nan)
    else:
        out.update(gv_sd=np.nan, min_glucose=v.min(), max_glucose=v.max(),
                   range_glucose=0.0, cv=np.nan, arv=np.nan, tw_arv=np.nan,
                   span_hours=np.nan, density_per_hour=np.nan)
    src = g["sources"].str.split("+")
    for name in ("poct", "central_lab", "blood_gas", "icu_charted"):
        out[f"frac_{name}"] = src.apply(lambda s: name in s).mean()
    for key, fn in THRESH.items():
        out[f"any_{key}"] = bool(fn(v).any())
    out["prop_lt70"] = (v < 70).mean()
    out["prop_gt180"] = (v > 180).mean()
    out["prop_70_180"] = ((v >= 70) & (v <= 180)).mean()
    out["first_glucose_offset_min"] = np.nan  # 由上游窗口锚点决定,聚合层不可恢复
    return out

def main():
    src, dst = sys.argv[1], sys.argv[2]
    s = pd.read_csv(src, parse_dates=["minute"])
    feat = s.groupby("stay_id", sort=True).apply(lambda g: pd.Series(per_stay(g)), include_groups=False).reset_index()
    feat.to_csv(dst, index=False)
    print(f"stays={len(feat)} rows written -> {dst}")

if __name__ == "__main__":
    main()
