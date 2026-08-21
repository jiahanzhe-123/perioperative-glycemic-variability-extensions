#!/usr/bin/env python3
"""Generate GV-centred main and supplementary figures from current aggregate results."""
from __future__ import annotations

# [pgv] 路径经配置驱动(src/common/paths.py);不得在此写入本机绝对路径。
import os as _os, sys as _sys
for _cand in ("src/common", "../src/common", "../../src/common", "../../../src/common"):
    if _os.path.exists(_os.path.join(_cand, "paths.py")):
        _sys.path.insert(0, _cand); break
else:
    raise RuntimeError("src/common/paths.py not found; run from repository root")
from paths import PGV_OUT  # noqa: E402

import csv
import os
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np


OUT = Path(PGV_OUT("results/figures"))
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7.0,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "figure.facecolor": "white",
})

BLUE = "#0F4D92"
BLUE2 = "#3775BA"
TEAL = "#42949E"
RED = "#B64342"
GREY = "#767676"
DARK = "#272727"
LIGHT = "#EAF1F8"
PALE = "#F5F6F7"


def save(fig: plt.Figure, stem: str) -> None:
    for ext, kwargs in {
        ".png": {"dpi": 600},
        ".pdf": {},
        ".svg": {},
    }.items():
        fig.savefig(OUT / f"{stem}{ext}", bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def panel(ax, label: str) -> None:
    ax.text(-0.10, 1.04, label, transform=ax.transAxes, fontsize=9, weight="bold", va="bottom")


def forest(ax, labels, est, lo, hi, xlabel, xlim, colors=None, annotation=None):
    y = np.arange(len(labels))[::-1]
    colors = colors or [BLUE] * len(labels)
    for yi, e, l, h, color in zip(y, est, lo, hi, colors):
        ax.plot([l, h], [yi, yi], color=color, lw=1.4, solid_capstyle="round")
        ax.scatter(e, yi, s=20, color=color, zorder=3, edgecolors="white", linewidths=0.4)
    ax.axvline(1, color=GREY, ls="--", lw=0.8, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    ax.set_ylim(-0.8, len(labels) - 0.2)
    ax.grid(axis="x", color="#E3E5E7", lw=0.6)
    ax.tick_params(axis="y", length=0)
    if annotation:
        for yi, text in zip(y, annotation):
            ax.text(xlim[1] * 0.995, yi, text, ha="right", va="center", fontsize=6.2, color=DARK)


def draw_box(ax, xy, width, height, header, lines, color=BLUE, fontsize=6.3):
    x, y = xy
    patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.015",
                           ec=color, fc="white", lw=1.1)
    ax.add_patch(patch)
    ax.text(x + 0.018, y + height - 0.04, header, va="top", ha="left", weight="bold", color=color, fontsize=7)
    ax.text(x + 0.018, y + height - 0.085, "\n".join(lines), va="top", ha="left", color=DARK,
            fontsize=fontsize, linespacing=1.3)
    return (x + width / 2, y + height / 2)


def arrow(ax, start, end, color=GREY):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9,
                                 color=color, lw=1.0, connectionstyle="arc3,rad=0"))


def figure_1():
    fig, ax = plt.subplots(figsize=(6.69, 4.30))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    mimic = draw_box(ax, (0.05, 0.56), 0.40, 0.32, "MIMIC-IV | primary", [
        "Date-anchored 24-h window",
        "Day-1 landmark",
        "10 561 patients; 296 deaths by day 30",
        "Model B: GV + flexible mean glucose",
        "Primary test: HR 0.979 (0.922–1.040)",
    ], BLUE)
    inspire = draw_box(ax, (0.55, 0.56), 0.40, 0.32, "INSPIRE | exact-timestamp extension", [
        "High-specificity open CABG / valve cohort",
        "opend + 0–24 h; 24-h landmark",
        "1 353 post-landmark patients; 27 deaths",
        "Model I2 24-h: HR 0.905 (0.653–1.254)",
        "Precision limited; not external validation",
    ], TEAL)
    eicu = draw_box(ax, (0.05, 0.15), 0.40, 0.26, "eICU-CRD | harmonized comparison", [
        "ICU-admission anchor; hospital mortality",
        "M3: RR 1.158 (1.025–1.308)",
        "M4 + count/span: RR 1.133 (0.990–1.296)",
        "No pooling; database-specific estimands",
    ], BLUE2)
    shr = draw_box(ax, (0.55, 0.15), 0.40, 0.26, "Nested MIMIC HbA1c component analysis", [
        "4 779 patients; SHR secondary only",
        "Joint SHR–GV: P = 0.158 (30 day)",
        "No incremental value beyond components",
        "Not repeated as a three-database exposure",
    ], GREY)
    arrow(ax, (0.45, 0.72), (0.54, 0.72), GREY)
    ax.text(0.5, 0.46, "Robustness and transportability checks", ha="center", va="center", fontsize=7.0,
            color=GREY)
    save(fig, "Figure_1_Study_Design")


def figure_2():
    fig, axes = plt.subplots(1, 2, figsize=(6.69, 2.72), sharex=True)
    rows = {
        "30 day": (
            ["Model A\n(no mean glucose)", "Model B\n(primary)", "Model C\n(+ count/span)"],
            [1.0472, 0.9794, 1.0286], [0.9940, 0.9225, 0.9671], [1.1032, 1.0398, 1.0940],
            ["P = 0.083", "P = 0.495", "P = 0.370"],
        ),
        "365 day": (
            ["Model A\n(no mean glucose)", "Model B\n(key secondary)", "Model C\n(+ count/span)"],
            [1.0368, 0.9850, 1.0243], [0.9966, 0.9420, 0.9794], [1.0787, 1.0299, 1.0713],
            ["P = 0.073", "P = 0.506", "P = 0.294"],
        ),
    }
    for i, (title, data) in enumerate(rows.items()):
        labels, e, lo, hi, ann = data
        forest(axes[i], labels, e, lo, hi, "Hazard ratio per 0.555 mmol/L", (0.82, 1.18),
               [GREY, BLUE, TEAL], ann)
        axes[i].set_title(title, fontsize=8, weight="bold", loc="left")
        panel(axes[i], chr(ord("a") + i))
    fig.subplots_adjust(wspace=0.62, left=0.19, right=0.98, top=0.85, bottom=0.24)
    save(fig, "Figure_2_MIMIC_Primary_GV")


def figure_3():
    fig = plt.figure(figsize=(6.69, 2.72))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.62)
    ax1 = fig.add_subplot(gs[0, 0])
    forest(ax1,
           ["30 day: Model A", "30 day: Model B", "365 day: Model A", "365 day: Model B"],
           [1.0472, 0.9794, 1.0368, 0.9850],
           [0.9940, 0.9225, 0.9966, 0.9420],
           [1.1032, 1.0398, 1.0787, 1.0299],
           "Hazard ratio per 0.555 mmol/L", (0.88, 1.16), [GREY, BLUE, GREY, BLUE])
    ax1.set_title("Mean-glucose conditioning", fontsize=8, weight="bold", loc="left")
    ax1.text(0.5, -0.35, "Model A→B is a conditioning comparison, not mediation", transform=ax1.transAxes,
             ha="center", va="top", fontsize=5.9, color=GREY)
    panel(ax1, "a")
    ax2 = fig.add_subplot(gs[0, 1])
    horizons = ["30 day", "365 day"]
    rd = np.array([-0.0442, -0.0268])
    lo = np.array([-0.2818, -0.5548])
    hi = np.array([0.1455, 0.4103])
    y = np.array([1, 0])
    ax2.axvline(0, color=GREY, ls="--", lw=0.8)
    ax2.errorbar(rd, y, xerr=[rd - lo, hi - rd], fmt="o", color=BLUE, markersize=5,
                 lw=1.4, capsize=2)
    ax2.set_yticks(y)
    ax2.set_yticklabels(horizons)
    ax2.set_xlabel("Risk difference, percentage points\n(Q75 vs Q25 GV)")
    ax2.set_xlim(-0.65, 0.65)
    ax2.set_ylim(-0.6, 1.6)
    ax2.grid(axis="x", color="#E3E5E7", lw=0.6)
    ax2.tick_params(axis="y", length=0)
    ax2.set_title("Absolute-risk contrasts", fontsize=8, weight="bold", loc="left")
    for yy, val in zip(y, ["−0.04 (−0.28 to +0.15)", "−0.03 (−0.55 to +0.41)"]):
        ax2.text(0.62, yy + 0.20, val, ha="right", va="center", fontsize=5.9)
    panel(ax2, "b")
    fig.subplots_adjust(left=0.17, right=0.98, top=0.84, bottom=0.26)
    save(fig, "Figure_3_Mean_Glucose_and_Absolute_Risk")


def figure_4():
    fig = plt.figure(figsize=(6.69, 3.26))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.42)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(-2, 74)
    ax1.set_ylim(-0.7, 2.2)
    ax1.set_yticks([1.55, 0.45])
    ax1.set_yticklabels(["Current 24-h analysis", "Superseded old 48-h analysis"])
    ax1.set_xticks([0, 24, 48, 72])
    ax1.set_xlabel("Hours after verified surgical completion")
    for spine in ["top", "right", "left"]:
        ax1.spines[spine].set_visible(False)
    ax1.tick_params(axis="y", length=0)
    ax1.hlines(1.55, 0, 24, color=TEAL, lw=7, alpha=0.9)
    ax1.hlines(1.55, 24, 72, color=LIGHT, lw=7)
    ax1.text(12, 1.78, "exposure", ha="center", fontsize=6.2)
    ax1.text(48, 1.78, "risk starts", ha="center", fontsize=6.2)
    ax1.vlines(24, 1.26, 1.84, color=DARK, lw=0.8)
    ax1.hlines(0.45, 0, 48, color=RED, lw=7, alpha=0.65)
    ax1.hlines(0.45, 24, 72, color=LIGHT, lw=7)
    ax1.axvspan(24, 48, ymin=0.08, ymax=0.50, color="#FDE7E5", zorder=-1)
    ax1.text(36, 0.78, "exposure–risk overlap\ninvalidates prior result", ha="center", va="bottom",
             color=RED, fontsize=6.1, weight="bold")
    ax1.vlines(24, 0.16, 0.74, color=DARK, lw=0.8)
    ax1.set_title("Landmark audit", fontsize=8, weight="bold", loc="left")
    panel(ax1, "a")
    ax2 = fig.add_subplot(gs[0, 1])
    forest(ax2,
           ["24-h landmark\ncurrent", "48-h landmark\ncorrected"],
           [0.9045542, 1.1043840], [0.6525616, 0.8464951], [1.2538560, 1.4408410],
           "Hazard ratio per 0.555 mmol/L", (0.45, 1.95), [TEAL, BLUE],
           ["N = 1 353; 27 deaths", "N = 1 511; 31 deaths"])
    ax2.set_title("Current estimates", fontsize=8, weight="bold", loc="left")
    panel(ax2, "b")
    fig.subplots_adjust(left=0.16, right=0.98, top=0.87, bottom=0.25)
    save(fig, "Figure_4_INSPIRE_Landmark_Audit")


def figure_5():
    fig, axes = plt.subplots(1, 2, figsize=(6.69, 3.10), gridspec_kw={"width_ratios": [1.0, 1.08]})
    labels = ["POCT only", "Central laboratory", "Blood gas", "Common source", "Same patients: POCT", "Same patients: laboratory"]
    est = [0.8268, 0.9785, 1.0572, 0.9775, 0.7589, 0.9709]
    lo = [0.7401, 0.9333, 0.9681, 0.9230, 0.6388, 0.9103]
    hi = [0.9238, 1.0260, 1.1544, 1.0353, 0.9015, 1.0355]
    forest(axes[0], labels, est, lo, hi, "Hazard ratio per 0.555 mmol/L", (0.55, 1.25),
           [TEAL, BLUE, BLUE2, GREY, TEAL, BLUE])
    axes[0].set_title("MIMIC-IV measurement-source estimates", fontsize=8, weight="bold", loc="left")
    axes[0].text(0.50, -0.20, "30-day Model B; source analyses are dependence diagnostics", transform=axes[0].transAxes,
                 ha="center", va="top", fontsize=5.7, color=GREY)
    panel(axes[0], "a")
    labels2 = ["MIMIC M3\n(+ mean glucose)", "MIMIC M4\n(+ count/span)",
               "eICU M3\n(+ mean glucose)", "eICU M4\n(+ count/span)"]
    est2 = [1.0241, 1.0471, 1.1577, 1.1328]
    lo2 = [0.9460, 0.9663, 1.0246, 0.9902]
    hi2 = [1.1087, 1.1347, 1.3082, 1.2960]
    forest(axes[1], labels2, est2, lo2, hi2, "Risk ratio per 0.555 mmol/L", (0.85, 1.40),
           [BLUE, BLUE, TEAL, TEAL])
    axes[1].set_title("Non-pooled harmonized comparison", fontsize=8, weight="bold", loc="left")
    axes[1].text(0.52, -0.20, "Post-landmark hospital mortality; different estimands", transform=axes[1].transAxes,
                 ha="center", va="top", fontsize=5.7, color=GREY)
    panel(axes[1], "b")
    fig.subplots_adjust(wspace=0.78, left=0.22, right=0.98, top=0.85, bottom=0.27)
    save(fig, "Figure_5_Measurement_Context_and_Databases")


def supplementary_figure_s1():
    fig = plt.figure(figsize=(6.69, 3.20))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.24, 1.0], wspace=0.52)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis("off")
    panel(ax1, "a")
    draw_box(ax1, (0.05, 0.69), 0.40, 0.19, "Eligible cohort", [
        "1 355 patients",
        "27 reconciled 30-day deaths",
    ], TEAL)
    draw_box(ax1, (0.55, 0.69), 0.40, 0.19, "Post-landmark", [
        "1 353 patients; 27 deaths",
        "2 excluded before landmark",
    ], BLUE)
    arrow(ax1, (0.45, 0.785), (0.54, 0.785))
    draw_box(ax1, (0.55, 0.28), 0.40, 0.23, "Final model frame", [
        "I1/I2/I3: 1 353 / 27",
        "Formula variables complete",
        "Common day-30 administrative end",
    ], BLUE)
    arrow(ax1, (0.75, 0.69), (0.75, 0.52), BLUE)
    draw_box(ax1, (0.05, 0.22), 0.40, 0.29, "Earlier mask (not retained)", [
        "1 343 / 23 earlier",
        "Excluded 10 patients",
        "and 4 events (ASA/BMI)",
        "not in 30-day formula",
    ], RED, fontsize=5.8)
    ax1.add_patch(FancyArrowPatch((0.43, 0.66), (0.43, 0.52), arrowstyle="-|>",
                                  mutation_scale=9, color=RED, lw=1.0,
                                  linestyle="--", connectionstyle="arc3,rad=0.12"))
    ax1.text(0.50, 0.58, "earlier-frame audit", ha="center", va="center", fontsize=5.7, color=RED)

    ax2 = fig.add_subplot(gs[0, 1])
    forest(ax2,
           ["I1", "I2 (primary)", "I3"],
           [0.8459585, 0.9045542, 0.9284988],
           [0.6340848, 0.6525616, 0.6723939],
           [1.1286280, 1.2538560, 1.2821500],
           "Hazard ratio per 0.555 mmol/L", (0.45, 1.50),
           [GREY, TEAL, BLUE],
           ["P = 0.255", "P = 0.547", "P = 0.652"])
    ax2.set_title("Final 30-day models", fontsize=8, weight="bold", loc="left")
    ax2.text(0.50, -0.23, "All models: N = 1 353; 27 deaths", transform=ax2.transAxes,
             ha="center", va="top", fontsize=5.8, color=GREY)
    panel(ax2, "b")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.26)
    save(fig, "Supplementary_Figure_S1_INSPIRE_Model_Frame_Audit")


def supplementary_figure_s2():
    fig = plt.figure(figsize=(6.69, 3.18))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.20, 1.0], wspace=0.53)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis("off")
    panel(ax1, "a")
    ax1.text(0.02, 0.91, "365-day follow-up audit", fontsize=8, weight="bold", ha="left")
    labels = [
        ("100 deaths", 100, RED),
        ("151 contact ≥365 d", 151, TEAL),
        ("1 102 no confirmed 365 d", 1102, GREY),
    ]
    yvals = [0.69, 0.48, 0.27]
    max_count = 1353
    for (label, count, color), y in zip(labels, yvals):
        width = 0.52 * count / max_count
        ax1.add_patch(FancyBboxPatch((0.05, y - 0.05), width, 0.10,
                                     boxstyle="round,pad=0.008,rounding_size=0.012",
                                     ec=color, fc=color, alpha=0.82, lw=0.8))
        ax1.text(0.60, y, label, ha="left", va="center", fontsize=5.8, color=DARK)
    ax1.text(0.05, 0.10, "The 1 102 include 11 discharged and 1 091 with last\nclinical contact before 365 days.",
             ha="left", va="bottom", fontsize=5.6, color=GREY, linespacing=1.25)
    ax1.text(0.05, 0.02, "Audit categories are not the Cox censoring mechanism.",
             ha="left", va="bottom", fontsize=5.6, color=RED)

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis("off")
    ax2.text(0.03, 0.91, "365-day analysis withdrawn", fontsize=8, weight="bold", ha="left")
    ax2.add_patch(FancyBboxPatch((0.05, 0.21), 0.87, 0.56,
                                 boxstyle="round,pad=0.02,rounding_size=0.02",
                                 ec=RED, fc="#FFF7F6", lw=1.0))
    ax2.text(0.10, 0.68, "No current INSPIRE 365-day Cox estimate", fontsize=6.6,
             weight="bold", color=RED, va="top")
    ax2.text(0.10, 0.56,
             "The earlier model allowed post-discharge deaths to\ncontribute time but censored non-events at discharge.\nNo common 365-day registry-coverage end is documented.\n\nNot carried forward: HR, absolute risk, RMST,\ncalibration, or performance.",
             fontsize=5.8, color=DARK, va="top", linespacing=1.30)
    panel(ax2, "b")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.18)
    save(fig, "Supplementary_Figure_S2_INSPIRE_365d_Analytic_Boundary")


def write_source_data():
    rows = [
        ["Figure", "Panel_or_element", "Measure", "Estimate", "Lower_CI", "Upper_CI", "P_value", "N", "Events", "Unit", "Source", "Notes"],
        ["2", "a", "30-day Model A", "1.047166", "0.993999", "1.103175", "0.083", "10 561", "296", "HR per 0.555 mmol/L", "MIMIC primary results", "No mean-glucose conditioning"],
        ["2", "a", "30-day Model B", "0.979381", "0.922463", "1.039810", "0.495", "10 561", "296", "HR per 0.555 mmol/L", "MIMIC primary results", "Primary model"],
        ["2", "a", "30-day Model C", "1.028573", "0.967085", "1.093971", "0.370", "10 561", "296", "HR per 0.555 mmol/L", "MIMIC primary results", "Adds measurement count and span"],
        ["2", "b", "365-day Model A", "1.036818", "0.996598", "1.078661", "0.073", "10 561", "745", "HR per 0.555 mmol/L", "MIMIC primary results", "No mean-glucose conditioning"],
        ["2", "b", "365-day Model B", "0.984981", "0.941996", "1.029927", "0.506", "10 561", "745", "HR per 0.555 mmol/L", "MIMIC primary results", "Key secondary model"],
        ["2", "b", "365-day Model C", "1.024297", "0.979398", "1.071253", "0.294", "10 561", "745", "HR per 0.555 mmol/L", "MIMIC primary results", "Adds measurement count and span"],
        ["3", "a", "30-day Model A", "1.047166", "0.993999", "1.103175", "0.083", "10 561", "296", "HR per 0.555 mmol/L", "MIMIC primary results", "Conditioning comparison"],
        ["3", "a", "30-day Model B", "0.979381", "0.922463", "1.039810", "0.495", "10 561", "296", "HR per 0.555 mmol/L", "MIMIC primary results", "Conditioning comparison"],
        ["3", "a", "365-day Model A", "1.036818", "0.996598", "1.078661", "0.073", "10 561", "745", "HR per 0.555 mmol/L", "MIMIC primary results", "Conditioning comparison"],
        ["3", "a", "365-day Model B", "0.984981", "0.941996", "1.029927", "0.506", "10 561", "745", "HR per 0.555 mmol/L", "MIMIC primary results", "Conditioning comparison"],
        ["3", "b", "30-day Q75-versus-Q25 risk difference", "-0.044153", "-0.281825", "0.145538", "", "9 751", "240", "percentage-point risk difference", "MIMIC absolute risk", "Model-standardized; 1000 bootstrap refits"],
        ["3", "b", "365-day Q75-versus-Q25 risk difference", "-0.026772", "-0.554800", "0.410297", "", "9 751", "633", "percentage-point risk difference", "MIMIC absolute risk", "Model-standardized; 1000 bootstrap refits"],
        ["4", "a", "Current exposure window", "0-24", "", "", "", "", "", "hours after verified surgical completion", "INSPIRE design", "Risk begins at 24 h"],
        ["4", "a", "Superseded 48-h exposure window", "0-48", "", "", "", "", "", "hours after verified surgical completion", "Historical audit", "Risk had incorrectly begun at 24 h"],
        ["4", "b", "Current 24-h Model I2", "0.904554", "0.652562", "1.253856", "0.547", "1 353", "27", "HR per 0.555 mmol/L", "INSPIRE administrative-censoring result", "Formula-specific full frame; common day-30 administrative end"],
        ["4", "b", "Corrected 48-h landmark", "1.104384", "0.846495", "1.440841", "0.464", "1 511", "31", "HR per 0.555 mmol/L", "INSPIRE administrative-censoring result", "Risk begins at 48 h; common day-30 administrative end"],
        ["5", "a", "POCT only", "0.826838", "0.740076", "0.923772", "0.000774", "7 602", "107", "HR per 0.555 mmol/L", "MIMIC source analysis", "30-day Model B"],
        ["5", "a", "Central laboratory", "0.978542", "0.933313", "1.025962", "0.369", "826", "129", "HR per 0.555 mmol/L", "MIMIC source analysis", "30-day Model B"],
        ["5", "a", "Blood gas", "1.057164", "0.968110", "1.154410", "0.216", "9 431", "156", "HR per 0.555 mmol/L", "MIMIC source analysis", "30-day Model B"],
        ["5", "a", "Common source", "0.977530", "0.922999", "1.035283", "0.438", "8 355", "202", "HR per 0.555 mmol/L", "MIMIC source analysis", "30-day Model B"],
        ["5", "a", "Same patients: POCT", "0.758857", "0.638762", "0.901532", "0.002", "409", "49", "HR per 0.555 mmol/L", "MIMIC source analysis", "30-day Model B"],
        ["5", "a", "Same patients: laboratory", "0.970897", "0.910334", "1.035489", "0.369", "409", "49", "HR per 0.555 mmol/L", "MIMIC source analysis", "30-day Model B"],
        ["5", "b", "MIMIC M3", "1.024100", "0.946000", "1.108700", "0.556", "8 117", "128", "RR per 0.555 mmol/L", "Harmonized comparison", "Mean-glucose conditioned"],
        ["5", "b", "MIMIC M4", "1.047100", "0.966300", "1.134700", "0.261", "8 117", "128", "RR per 0.555 mmol/L", "Harmonized comparison", "Adds measurement count and span"],
        ["5", "b", "eICU M3", "1.157748", "1.024578", "1.308226", "0.019", "7 115", "130", "RR per 0.555 mmol/L", "Harmonized comparison", "Mean-glucose conditioned"],
        ["5", "b", "eICU M4", "1.132796", "0.990176", "1.295958", "0.069", "7 115", "130", "RR per 0.555 mmol/L", "Harmonized comparison", "Adds measurement count and span"],
        ["S1", "a", "Eligible cohort", "1 355", "", "", "", "1 355", "27", "patients / deaths", "INSPIRE frame audit", "Eligible 24-h landmark cohort"],
        ["S1", "a", "Pre-landmark exclusions", "2", "", "", "", "", "", "patients", "INSPIRE frame audit", "Deaths before the 24-h landmark"],
        ["S1", "a", "Post-landmark analytic frame", "1 353", "", "", "", "1 353", "27", "patients / deaths", "INSPIRE frame audit", "All 30-day formula variables complete; common day-30 administrative end"],
        ["S1", "a", "Superseded wide complete-case mask", "1 343", "", "", "", "1 343", "23", "patients / events", "Historical audit", "Excluded 10 patients and 4 events because ASA/BMI were retained"],
        ["S1", "b", "I1 24-h", "0.845959", "0.634085", "1.128628", "0.255", "1 353", "27", "HR per 0.555 mmol/L", "INSPIRE administrative-censoring result", "Formula-specific full frame"],
        ["S1", "b", "I2 24-h primary", "0.904554", "0.652562", "1.253856", "0.547", "1 353", "27", "HR per 0.555 mmol/L", "INSPIRE administrative-censoring result", "Formula-specific full frame"],
        ["S1", "b", "I3 24-h", "0.928499", "0.672394", "1.282150", "0.652", "1 353", "27", "HR per 0.555 mmol/L", "INSPIRE administrative-censoring result", "Formula-specific full frame"],
        ["S2", "a", "Reconciled 365-day deaths", "100", "", "", "", "", "100", "deaths", "INSPIRE follow-up audit", "Audit category; not the model censoring mechanism"],
        ["S2", "a", "Documented clinical contact at least 365 days", "151", "", "", "", "", "", "patients", "INSPIRE follow-up audit", "Audit category; not the model censoring mechanism"],
        ["S2", "a", "No confirmed 365-day survival", "1 102", "", "", "", "", "", "patients", "INSPIRE follow-up audit", "11 discharged and 1 091 last clinical contact before 365 days"],
        ["S2", "a", "Discharged within no-confirmation category", "11", "", "", "", "", "", "patients", "INSPIRE follow-up audit", "Subset of 1 102"],
        ["S2", "a", "Last clinical contact before 365 days", "1 091", "", "", "", "", "", "patients", "INSPIRE follow-up audit", "Subset of 1 102"],
        ["S2", "b", "365-day Cox results", "", "", "", "", "", "", "withdrawn", "INSPIRE analytic boundary", "earlier event-dependent discharge censoring; no documented common 365-day registry-coverage end"],
        ["S2", "b", "365-day derived measures", "", "", "", "", "", "", "withdrawn", "INSPIRE analytic boundary", "No HR, absolute risk, RMST, calibration, or performance result carried forward"],
    ]
    with (OUT / "Figure_Source_Data.csv").open("w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(rows)


if __name__ == "__main__":
    figure_1()
    figure_2()
    figure_3()
    figure_4()
    figure_5()
    supplementary_figure_s1()
    supplementary_figure_s2()
    write_source_data()
