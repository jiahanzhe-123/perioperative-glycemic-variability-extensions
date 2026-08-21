#!/usr/bin/env python3
# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV_OUT  # noqa: E402
import sys
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','axes.edgecolor':'black'})
import matplotlib.pyplot as plt

OUT = PGV_OUT("results/figures")
P4 = PGV_OUT("results/figures")
C_GRAY, C_BLUE, C_TEAL = '#7f7f7f', '#1f5f9f', '#2a9d8f'

def style_ax(ax):
    for s in ['top','right']: ax.spines[s].set_visible(False)
    ax.grid(axis='x', color='#e0e0e0', lw=1)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)

# ---------------- Figure 2 ----------------
fig, axes = plt.subplots(1, 2, figsize=(12.7, 4.75), dpi=300)
fig.subplots_adjust(left=0.16, right=0.9, top=0.82, bottom=0.18, wspace=0.75)

ax = axes[0]
rows = [('Model A\n(no mean glucose)', 1.047, 0.994, 1.103, 'P = 0.083', C_GRAY),
        ('Model B\n(primary)',        0.979, 0.922, 1.040, 'P = 0.495', C_BLUE),
        ('Model C\n(+ count/span)',    1.029, 0.967, 1.094, 'P = 0.370', C_TEAL)]
for i,(lab,hr,lo,hi,p,c) in enumerate(rows):
    y = len(rows)-1-i
    ax.errorbar(hr, y, xerr=[[hr-lo],[hi-hr]], fmt='o', color=c, ms=9, elinewidth=3.5, capsize=0)
    ax.annotate(p, (1.115, y), va='center', fontsize=10, color='#333333')
ax.axvline(1.0, ls='--', color='#999999', lw=1.2)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=10.5)
ax.set_xlim(0.87, 1.15); ax.set_ylim(-0.6, len(rows)-0.4)
ax.set_xticks([0.9,1.0,1.1])
ax.set_xlabel('Hazard ratio per 0.555 mmol/L', fontsize=11)
ax.set_title('30 day', fontsize=13, fontweight='bold', loc='center', pad=10)
ax.text(-0.42, 1.13, 'a', transform=ax.transAxes, fontsize=16, fontweight='bold')
style_ax(ax)

ax = axes[1]
rows = [('post-landmark\ndays 1-7\n(162 events)', 1.063, 0.998, 1.132, 'P = 0.057'),
        ('post-landmark\ndays 8-30\n(139 events)', 0.974, 0.889, 1.067, 'P = 0.569'),
        ('post-landmark\ndays 31-365\n(444 events)', 0.938, 0.883, 0.996, 'P = 0.037')]
for i,(lab,hr,lo,hi,p) in enumerate(rows):
    y = len(rows)-1-i
    ax.errorbar(hr, y, xerr=[[hr-lo],[hi-hr]], fmt='o', color=C_BLUE, ms=9, elinewidth=3.5, capsize=0)
    ax.annotate(p, (1.13, y), va='center', fontsize=10, color='#333333')
ax.axvline(1.0, ls='--', color='#999999', lw=1.2)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=9.5)
ax.set_xlim(0.85, 1.15); ax.set_ylim(-0.6, len(rows)-0.4)
ax.set_xticks([0.9,1.0,1.1])
ax.set_xlabel('Hazard ratio per 0.555 mmol/L (Model B)', fontsize=11)
ax.set_title('365 day: interval-specific estimates', fontsize=13, fontweight='bold', loc='center', pad=10)
ax.text(-0.5, 1.13, 'b', transform=ax.transAxes, fontsize=16, fontweight='bold')
style_ax(ax)

fig.savefig(f'{OUT}/fig2_primary_time.png', bbox_inches='tight')
for ext in ['png','pdf','svg']:
    fig.savefig(f'{P4}/Figure_2_MIMIC_Primary_GV.{ext}', bbox_inches='tight')
plt.close(fig)

# ---------------- Figure 3 ----------------
fig, axes = plt.subplots(1, 2, figsize=(12.85, 5.17), dpi=300)
fig.subplots_adjust(left=0.15, right=0.9, top=0.82, bottom=0.2, wspace=0.6)

ax = axes[0]
rows = [('30 day: Model A', 1.047, 0.994, 1.103, C_GRAY),
        ('30 day: Model B', 0.979, 0.922, 1.040, C_BLUE),
        ('365 day: Model A', 1.037, 0.997, 1.079, C_GRAY),
        ('365 day: Model B', 0.985, 0.942, 1.030, C_BLUE)]
for i,(lab,hr,lo,hi,c) in enumerate(rows):
    y = len(rows)-1-i
    ax.errorbar(hr, y, xerr=[[hr-lo],[hi-hr]], fmt='o', color=c, ms=8, elinewidth=3, capsize=0)
ax.axvline(1.0, ls='--', color='#999999', lw=1.2)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=11)
ax.set_xlim(0.88, 1.17); ax.set_ylim(-0.6, len(rows)-0.4)
ax.set_xticks([0.90,0.95,1.00,1.05,1.10,1.15])
ax.set_xlabel('Hazard ratio per 0.555 mmol/L', fontsize=11)
ax.set_title('Mean-glucose conditioning', fontsize=13, fontweight='bold', loc='center', pad=10)
ax.text(-0.38, 1.12, 'a', transform=ax.transAxes, fontsize=16, fontweight='bold')
style_ax(ax)

ax = axes[1]
ax.errorbar(-0.04, 0, xerr=[[0.24],[0.19]], fmt='o', color=C_BLUE, ms=9, elinewidth=3.5, capsize=0)
ax.annotate('\u22120.04 (\u22120.28 to +0.15)', (0.02, 0.28), fontsize=11, color='#333333')
ax.axvline(0.0, ls='--', color='#999999', lw=1.2)
ax.set_yticks([0]); ax.set_yticklabels(['30 day'], fontsize=11)
ax.set_xlim(-0.35, 0.3); ax.set_ylim(-0.7, 0.9)
ax.set_xticks([-0.25,0.00,0.25])
ax.set_xlabel('Risk difference, percentage points\n(Q75 vs Q25 GV)', fontsize=11)
ax.set_title('Absolute-risk contrast', fontsize=13, fontweight='bold', loc='center', pad=10)
ax.text(-0.28, 1.12, 'b', transform=ax.transAxes, fontsize=16, fontweight='bold')
style_ax(ax)

fig.text(0.155, 0.02, 'Model A\u2192B is a conditioning comparison, not mediation', fontsize=10, color='#888888')
fig.savefig(f'{OUT}/fig3_conditioning_risk.png', bbox_inches='tight')
for ext in ['png','pdf','svg']:
    fig.savefig(f'{P4}/Figure_3_Mean_Glucose_and_Absolute_Risk.{ext}', bbox_inches='tight')
plt.close(fig)
print('saved')
