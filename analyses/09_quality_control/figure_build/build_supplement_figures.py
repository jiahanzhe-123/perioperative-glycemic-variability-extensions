#!/usr/bin/env python3
# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV, PGV_OUT  # noqa: E402
"""Supplement figures S6/S7(old)/S11(old)/S12(old): drawn from the BMI-repaired
analysis frame and result CSVs. Journal-style, matplotlib Agg."""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FRAME = Path(PGV("mimic_record_work")) / "data" / "analysis_base_bmi_repaired.csv"
RES = Path(PGV("bmi_repair_work")) / "rerun_root" / "results"
OUT = Path(PGV_OUT("results/figures"))
OUT.mkdir(parents=True, exist_ok=True)

BLUE = '#0F4D92'; BLUE2 = '#3775BA'; TEAL = '#42949E'; RED = '#B64342'
GREY = '#767676'; DARK = '#272727'
plt.rcParams.update({'font.size': 7.0, 'axes.spines.right': False,
                     'axes.spines.top': False, 'axes.linewidth': 0.8,
                     'legend.frameon': False, 'figure.facecolor': 'white'})

df = pd.read_csv(FRAME)
prim = df[(df['gv'].notna()) & (df['survival_time_days'] >= 1)]
nest = prim[prim['hba1c_window'] == 'pre_90d'].copy()
assert len(prim) == 10561 and len(nest) == 4779, (len(prim), len(nest))
cc = nest[nest[['shr', 'gv']].notna().all(axis=1)]

def panel(ax, label):
    ax.text(-0.10, 1.04, label, transform=ax.transAxes, fontsize=9,
            weight='bold', va='bottom')

# ---------------- S6: nested SHR and GV distributions ----------------
fig, axes = plt.subplots(1, 4, figsize=(6.69, 1.95))
specs = [('shr', 'SHR', TEAL), ('gv', 'GV SD, mg/dL', BLUE),
         ('glucose_mean_postop_24h', 'Mean glucose, mg/dL', BLUE2),
         ('hba1c_pct', 'HbA1c, %', GREY)]
for ax, (col, lab, color) in zip(axes, specs):
    v = nest[col].dropna()
    ax.hist(v, bins=40, color=color, alpha=0.85, lw=0)
    ax.axvline(v.median(), color=DARK, ls='--', lw=0.8)
    ax.set_xlabel(lab, fontsize=6.5)
    ax.tick_params(labelsize=5.8)
axes[0].set_ylabel('Patients', fontsize=6.5)
fig.suptitle(f'Nested HbA1c cohort (N = {len(nest):,})'.replace(',', ' '),
             fontsize=7.5, weight='bold', x=0.02, ha='left')
for ax, lab in zip(axes, 'abcd'):
    panel(ax, lab)
fig.subplots_adjust(wspace=0.42, left=0.07, right=0.99, top=0.80, bottom=0.20)
fig.savefig(OUT / 'Supplementary_Figure_S6_Nested_Distributions.png', dpi=600,
            bbox_inches='tight', facecolor='white')
plt.close(fig)

# ---------------- S7: SHR-GV empirical support + corners ----------------
fig, ax = plt.subplots(figsize=(3.35, 3.0))
ax.scatter(cc['gv'], cc['shr'], s=1.5, color=BLUE, alpha=0.10, lw=0)
gv25, gv75 = cc['gv'].quantile([0.25, 0.75])
shr25, shr75 = cc['shr'].quantile([0.25, 0.75])
corner_style = [  # (x, y, label, dx, dy, ha)
    (gv25, shr25, 'SHR P25 / GV P25', -8, -14, 'right'),
    (gv25, shr75, 'SHR P75 / GV P25', -8, 8, 'right'),
    (gv75, shr25, 'SHR P25 / GV P75', 8, -14, 'left'),
    (gv75, shr75, 'SHR P75 / GV P75', 8, 8, 'left')]
for x, y, lab, dx, dy, ha in corner_style:
    ax.scatter([x], [y], s=42, color=RED, zorder=5, edgecolors='white', linewidths=0.6)
    ax.annotate(lab, (x, y), textcoords='offset points', xytext=(dx, dy), fontsize=5.6,
                color=RED, ha=ha)
ax.set_xlim(0, 80)
ax.set_ylim(0.4, 2.2)
ax.set_xlabel('GV SD, mg/dL', fontsize=6.5)
ax.set_ylabel('SHR', fontsize=6.5)
ax.set_title(f'Empirical SHR\u2013GV support (joint complete case, N = {len(cc):,})'.replace(',', ' '),
             fontsize=7.0, weight='bold', loc='left')
ax.tick_params(labelsize=5.8)
fig.subplots_adjust(left=0.14, right=0.97, top=0.90, bottom=0.16)
fig.savefig(OUT / 'Supplementary_Figure_S7_SHR_GV_Support.png', dpi=600,
            bbox_inches='tight', facecolor='white')
plt.close(fig)
print('S7 corners: gv25 %.3f gv75 %.3f shr25 %.3f shr75 %.3f' % (gv25, gv75, shr25, shr75))

# ---------------- S11: PH audit ----------------
ph = pd.read_csv(f'{RES}/07_ph_diagnostics.csv')
terms30 = ph[(ph['model_id'] == 'ModelB_30d')]
terms365 = ph[(ph['model_id'] == 'ModelB_365d')]
iv = ph[ph['model_id'].str.contains('GV_intervals')].copy()
fig, axes = plt.subplots(1, 2, figsize=(6.69, 2.9),
                         gridspec_kw={'width_ratios': [1.15, 1.0]})
ax = axes[0]
term_labels = {
    'gv10': 'GV per 0.555 mmol/L',
    'rms::rcs(mean_glu, c(113, 129.7273, 143.2308, 187))': 'Mean glucose (RCS)',
    'age_at_admission': 'Age', 'gender': 'Recorded sex', 'bmi': 'BMI',
    'diabetes': 'Diabetes', 'charlson_without_diabetes': 'Charlson excl. diabetes',
    'procedure_cat6': 'Procedure category', 'lactate_postop_first': 'Lactate',
    'creat_postop_first': 'Creatinine', 'sofa_24h': 'First-day SOFA',
    'GLOBAL': 'GLOBAL'}
terms = [t for t in terms30['term']]
y = np.arange(len(terms))[::-1]
p30 = [float(terms30[terms30['term'] == t]['p'].iloc[0]) for t in terms]
p365 = [float(terms365[terms365['term'] == t]['p'].iloc[0]) for t in terms]
ax.scatter(p30, y + 0.16, s=16, color=BLUE, label='30 days', zorder=3)
ax.scatter(p365, y - 0.16, s=16, color=TEAL, label='365 days', zorder=3)
ax.axvline(0.05, color=RED, ls='--', lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels([term_labels.get(t, t) for t in terms], fontsize=5.6)
ax.set_xlabel('Schoenfeld residual P value (complete-case Model B)', fontsize=6.3)
ax.set_xlim(-0.03, 1.03)
ax.legend(fontsize=5.8, loc='upper right', bbox_to_anchor=(1.0, 1.02))
ax.tick_params(axis='x', labelsize=5.8)
panel(ax, 'a')
ax = axes[1]
iv['horizon'] = iv['model_id'].str.extract(r'(30d|365d)')
order = [('30d', 'd1-7'), ('30d', 'd8-30'), ('365d', 'd1-7'),
         ('365d', 'd8-30'), ('365d', 'd31-365')]
labels = ['30 d: days 1\u20137', '30 d: days 8\u201330', '365 d: days 1\u20137',
          '365 d: days 8\u201330', '365 d: days 31\u2013365']
hrs = [float(iv[(iv['horizon'] == h) & (iv['interval'] == i)]['HR_per10'].iloc[0])
       for h, i in order]
evs = [int(iv[(iv['horizon'] == h) & (iv['interval'] == i)]['interval_events'].iloc[0])
       for h, i in order]
ps = [float(iv[(iv['horizon'] == h) & (iv['interval'] == i)]['p'].iloc[0])
      for h, i in order]
y = np.arange(len(order))[::-1]
colors = [BLUE if h == '30d' else TEAL for h, _ in order]
ax.scatter(hrs, y, s=18, color=colors, zorder=3)
for yi, hr, p in zip(y, hrs, ps):
    ax.text(1.24, yi, f'HR {hr:.3f}; P = {p:.3f}' if p >= 0.001 else f'HR {hr:.3f}; P < 0.001',
            fontsize=5.4, va='center', color=DARK)
ax.axvline(1, color=GREY, ls='--', lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels([f'{l} (n = {e})' for l, e in zip(labels, evs)], fontsize=5.6)
ax.set_xlabel('GV interval-specific HR per 0.555 mmol/L', fontsize=6.3)
ax.set_xlim(0.85, 1.22)
ax.tick_params(axis='x', labelsize=5.8)
panel(ax, 'b')
fig.subplots_adjust(wspace=0.72, left=0.14, right=0.86, top=0.93, bottom=0.18)
fig.savefig(OUT / 'Supplementary_Figure_S11_PH_Audit.png', dpi=600,
            bbox_inches='tight', facecolor='white')
plt.close(fig)

# ---------------- S12: prespecified model sequence (MICE) ----------------
mice = pd.read_csv(f'{RES}/mice_pooled_models.csv')
fig, axes = plt.subplots(1, 2, figsize=(6.69, 2.35))
for ax, horizon, ttl in [(axes[0], '30d', 'Mortality by index day 30 (Model B primary)'),
                         (axes[1], '365d', 'Mortality by index day 365 (Model B key secondary)')]:
    sub = mice[mice['model_id'].str.endswith(horizon)]
    labels = ['Model A', 'Model B', 'Model C']
    est = sub['HR_per10'].astype(float).values
    lo = sub['lo_per10'].astype(float).values
    hi = sub['hi_per10'].astype(float).values
    ps = sub['P_per10'].astype(float).values
    y = np.arange(3)[::-1]
    ax.scatter(est, y, s=22, color=[GREY, BLUE, BLUE2], zorder=3)
    for yi, e, l, h in zip(y, est, lo, hi):
        ax.plot([l, h], [yi, yi], color=DARK, lw=1.0, zorder=2)
    for yi, p in zip(y, ps):
        ax.text(1.16, yi, f'P = {p:.3f}' if p >= 0.001 else 'P < 0.001',
                fontsize=5.6, va='center', color=DARK)
    ax.axvline(1, color=GREY, ls='--', lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.0)
    ax.set_xlabel('HR per 0.555 mmol/L', fontsize=6.3)
    ax.set_title(ttl, fontsize=7.0, weight='bold', loc='left')
    ax.set_xlim(0.88, 1.15)
    ax.tick_params(labelsize=5.8)
panel(axes[0], 'a'); panel(axes[1], 'b')
fig.subplots_adjust(wspace=0.52, left=0.14, right=0.90, top=0.86, bottom=0.22)
fig.savefig(OUT / 'Supplementary_Figure_S12_Model_Sequence.png', dpi=600,
            bbox_inches='tight', facecolor='white')
plt.close(fig)
print('S6/S7/S11/S12 drawn')
