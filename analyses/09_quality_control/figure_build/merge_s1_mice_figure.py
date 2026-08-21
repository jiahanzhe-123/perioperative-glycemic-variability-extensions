#!/usr/bin/env python3
# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
"""merge_s1_mice_figure.py — Supplementary Figure S1(MICE 诊断)纵向合成。

将 02_missing_data_report.R 生成的 mice_trace.png(收敛轨迹)与
mice_density.png(observed vs imputed 密度)纵向合并为单一图版,
与投稿版 figS1_mice_merged_vertical.png 的版面一致(上 trace、下 density,
左对齐,白底,600dpi 源图缩放)。

输入:PGV("mimic_record_work")/figures/mice_trace.png, mice_density.png
输出:PGV_OUT("results/figures")/figS1_mice_merged_vertical.png
"""
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV, PGV_OUT  # noqa: E402
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

FIG = Path(PGV("mimic_record_work")) / "figures"
OUT = Path(PGV_OUT("results/figures"))
OUT.mkdir(parents=True, exist_ok=True)

trace = mpimg.imread(FIG / "mice_trace.png")
dens = mpimg.imread(FIG / "mice_density.png")

fig, axes = plt.subplots(2, 1, figsize=(7.0, 9.5))
for ax, img, lab in zip(axes, (trace, dens), ("a", "b")):
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(lab, loc="left", fontweight="bold", fontsize=12)
fig.suptitle("MICE m=50: convergence traces and observed-vs-imputed densities",
             fontsize=10, fontweight="bold", x=0.02, ha="left")
fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.01, hspace=0.03)
dst = OUT / "figS1_mice_merged_vertical.png"
fig.savefig(dst, dpi=600, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("saved", dst)
