#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild Figure 3 from the associated retained plotting-source TSV.

This file performs figure generation only. It does not rerun GWAS, LDSC, or
GenomicSEM. The original data table requires separate redistribution approval.
Run: python draw_figure3.py  (from this script directory)
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "FIGURE3_PLOT_SOURCE.tsv")

MINUS = "\u2212"   # typographic minus, used in panel d interval strings


def fmt_range(lo, hi, nd=3, negative_possible=False):
    """Format a confidence interval without ever producing '--'."""
    if negative_possible:
        return ("(" + ("%." + str(nd) + "f") % lo + " to "
                + ("%." + str(nd) + "f") % hi + ")").replace("-", MINUS)
    return "(" + ("%." + str(nd) + "f") % lo + "\u2013" + ("%." + str(nd) + "f") % hi + ")"


def load_source():
    with open(SRC, encoding="utf-8-sig", newline="") as fh:
        recs = [r for r in csv.DictReader(fh, delimiter="\t") if r["panel_figure"] == "Fig3"]
    if not recs:
        raise SystemExit("no Fig3 rows in " + SRC)

    rows = []          # panel a/b, ordered right-factor, right-residual, left-factor, left-residual
    rg_d = {}
    cv = {}
    for side in ("right", "left"):
        for construct in ("Regional FFM factor", "Residual FFM"):
            a = [r for r in recs if r["panel"] == "a" and r["side"] == side and r["construct"] == construct][0]
            b = [r for r in recs if r["panel"] == "b" and r["side"] == side and r["construct"] == construct][0]
            rows.append({
                "side": side, "construct": construct,
                "rg": float(a["value"]), "rg_se": float(a["se"]),
                "rg_lo": float(a["ci_low"]), "rg_hi": float(a["ci_high"]),
                "cov": float(b["value"]), "cov_se": float(b["se"]),
                "cov_lo": float(b["ci_low"]), "cov_hi": float(b["ci_high"]),
            })
    for side in ("right", "left"):
        c = [r for r in recs if r["panel"] == "c" and r["side"] == side][0]
        d = [r for r in recs if r["panel"] == "d" and r["side"] == side][0]
        rg_d[side] = {"estimate": float(c["value"]), "CI_low": float(c["ci_low"]),
                      "CI_high": float(c["ci_high"]), "p": float(c["p_two_sided"])}
        cv[side] = {"estimate": float(d["value"]), "CI_low": float(d["ci_low"]),
                    "CI_high": float(d["ci_high"]), "p": float(d["p_two_sided"])}
    return rows, rg_d, cv


rows, rg_d, cv = load_source()

plt.rcParams.update({"font.family": ["Arimo", "Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 8.6, "svg.fonttype": "none",
                     "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.6,
                     "savefig.facecolor": "white"})
fig = plt.figure(figsize=(180 / 25.4, 179 / 25.4), facecolor="white")
NAVY = "#1E4E83"
LIGHT = "#dbe7f1"
BLACK = "#111111"


def txt(x, y, s, size=8.6, bold=False, ha="left", color=BLACK):
    fig.text(x, y, s, fontsize=size, fontweight="bold" if bold else "normal",
             va="center", ha=ha, color=color)


def line(y):
    fig.add_artist(Line2D([.035, .973], [y, y], transform=fig.transFigure,
                          color="#1a1a1a", lw=.6))


def panel_letter(y, k, title):
    txt(.038, y, k, 16, False)
    txt(.08, y, title, 10.1, True)


def ax_forest(x, y, w, h, xmin, xmax, xticks, ys):
    ax = fig.add_axes((x, y, w, h))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0.35, len(ys) + 0.62)
    ax.set_yticks([])
    ax.set_xticks(xticks)
    ax.tick_params(axis="x", labelsize=8, length=3, pad=3, width=.55)
    for t in xticks:
        ax.axvline(t, color=LIGHT, lw=.52, zorder=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#777777")
    ax.spines["bottom"].set_linewidth(.55)
    return ax


# ---------------------------------------------------------------- panel a
panel_letter(.968, "a", "Genetic correlation with adjusted ALM")
txt(.054, .923, "Model", 8.5, True)
txt(.176, .923, "Construct", 8.5, True)
txt(.400, .923, "Genetic correlation", 8.5, True)
txt(.680, .923, "r_g (95% CI)", 8.5, True)
txt(.932, .923, "SE", 8.5, True, ha="right")
line(.907)
ys = [4, 3, 2, 1]
ax = ax_forest(.388, .737, .222, .163, .69, 1.04, [.70, .80, .90, 1.00], ys)
for i, (r, y) in enumerate(zip(rows, ys)):
    yy = [.867, .826, .787, .748][i]
    if i % 2 == 0:
        txt(.054, yy, r["side"].title(), 8.8, True)
    txt(.176, yy, r["construct"], 8.7)
    ax.errorbar(r["rg"], y,
                xerr=[[r["rg"] - r["rg_lo"]], [r["rg_hi"] - r["rg"]]],
                fmt="D" if i % 2 == 0 else "o",
                mfc="white" if i % 2 == 0 else NAVY, mec=NAVY, ecolor=NAVY,
                ms=5.3, capsize=3, lw=1.4, zorder=4)
    txt(.68, yy, "%.3f " % r["rg"] + fmt_range(r["rg_lo"], r["rg_hi"]), 8.3)
    txt(.948, yy, "%.4f" % r["rg_se"], 8.3, ha="right")

# ---------------------------------------------------------------- panel b
panel_letter(.700, "b", "Unstandardized genetic covariance with adjusted ALM")
txt(.054, .658, "Model", 8.5, True)
txt(.176, .658, "Construct", 8.5, True)
txt(.397, .658, "Genetic covariance", 8.5, True)
txt(.680, .658, "Covariance (95% CI)", 8.5, True)
txt(.948, .658, "SE", 8.5, True, ha="right")
line(.643)
ax = ax_forest(.395, .473, .215, .160, .17, .285, [.18, .21, .24, .27], ys)
for i, (r, y) in enumerate(zip(rows, ys)):
    yy = [.611, .572, .532, .493][i]
    if i % 2 == 0:
        txt(.054, yy, r["side"].title(), 8.8, True)
    txt(.176, yy, r["construct"], 8.7)
    ax.errorbar(r["cov"], y,
                xerr=[[r["cov"] - r["cov_lo"]], [r["cov_hi"] - r["cov"]]],
                fmt="D" if i % 2 == 0 else "o",
                mfc="white" if i % 2 == 0 else NAVY, mec=NAVY, ecolor=NAVY,
                ms=5.1, capsize=3, lw=1.4, zorder=4)
    txt(.680, yy, "%.3f " % r["cov"] + fmt_range(r["cov_lo"], r["cov_hi"]), 8.3)
    txt(.948, yy, "%.4f" % r["cov_se"], 8.3, ha="right")

# ---------------------------------------------------------------- panel c
panel_letter(.429, "c", "Paired residual-minus-factor genetic correlation")
txt(.054, .391, "Model", 8.4, True)
txt(.37, .391, "Paired \u0394r_g", 8.4, True)
txt(.650, .391, "\u0394r_g (95% CI)", 8.4, True)
txt(.945, .391, "P", 8.4, True, ha="right")
line(.378)
ax = ax_forest(.365, .287, .250, .077, 0, .29, [0, .1, .2], (2, 1))
ax.axvline(0, color="#666666", ls="--", lw=.7, zorder=1)
for i, side in enumerate(("right", "left")):
    yy = [.344, .306][i]
    y = 2 - i
    x = rg_d[side]
    txt(.054, yy, "Right primary" if i == 0 else "Left-side model", 8.7, True)
    ax.errorbar(x["estimate"], y,
                xerr=[[x["estimate"] - x["CI_low"]], [x["CI_high"] - x["estimate"]]],
                fmt="o", mfc=NAVY, mec=NAVY, ecolor=NAVY, ms=4.9, capsize=3, lw=1.4)
    txt(.65, yy, "%+.3f " % x["estimate"] + fmt_range(x["CI_low"], x["CI_high"], 3), 8.2)
    txt(.946, yy, "<0.001" if x["p"] < 1e-3 else "%.3f" % x["p"], 8.3, ha="right")

# ---------------------------------------------------------------- panel d
panel_letter(.204, "d", "Paired residual-minus-factor genetic covariance")
txt(.054, .164, "Model", 8.4, True)
txt(.37, .164, "Paired \u0394Cov_g", 8.4, True)
txt(.650, .164, "\u0394Cov_g (95% CI)", 8.4, True)
txt(.945, .164, "P", 8.4, True, ha="right")
line(.152)
ax = ax_forest(.365, .060, .250, .078, -.043, .005, [-.04, -.02, 0], (2, 1))
ax.axvline(0, color="#666666", ls="--", lw=.7, zorder=1)
for i, side in enumerate(("right", "left")):
    yy = [.116, .078][i]
    y = 2 - i
    h = cv[side]
    e, lo, hi = h["estimate"], h["CI_low"], h["CI_high"]
    txt(.054, yy, "Right primary" if i == 0 else "Left-side model", 8.7, True)
    ax.errorbar(e, y, xerr=[[e - lo], [hi - e]], fmt="o", mfc=NAVY, mec=NAVY,
                ecolor=NAVY, ms=4.9, capsize=3, lw=1.4)
    txt(.650, yy, ("%.4f " % e).replace("-", MINUS) + fmt_range(lo, hi, 4, negative_possible=True), 8.15)
    txt(.946, yy, "<0.001" if h["p"] < 1e-3 else "%.3f" % h["p"], 8.3, ha="right")

# ---------------------------------------------------------------- output
fig.savefig(os.path.join(HERE, "Fig3.svg"), format="svg", facecolor="white")
fig.savefig(os.path.join(HERE, "Fig3.pdf"), format="pdf", facecolor="white")
fig.savefig(os.path.join(HERE, "Fig3.eps"), format="eps", facecolor="white")
fig.savefig(os.path.join(HERE, "Fig3_600dpi.png"), format="png", dpi=600, facecolor="white")
plt.close(fig)
print("Figure 3 written to", HERE)
