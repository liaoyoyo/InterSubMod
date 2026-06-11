"""
C3_AF_4bin_sweep.py
====================
Phase 2 completeness audit C3 — AF (allele frequency) 4-bin full sweep
for baseline vs V6 (longphase NG≥3 marker filter).

Motivation
----------
The upstream 6-bin counterexample TSV (`step1_counterexample_per_af_bin.tsv`)
shows AF<0.1 is a V6 reversal stratum (Δ marker rate = -0.0304).
This script collapses to the canonical 4 bins and reports per-bin:

  • n_regions / n_TP / n_FP (all regions, regardless of marker status)
  • baseline marker count + marker_rate (= TP / (TP+FP) among baseline markers)
  • V6 marker count + marker_rate (same for V6)
  • Δ marker count / Δ marker rate (pp)
  • ΔF1 proxy = V6_marker_rate − baseline_marker_rate (per bin)

"Marker" definition (matches upstream convention)
-------------------------------------------------
A region is a marker if `{bam}_off_NG >= 3` (HP-off NumGroups ≥ 3), i.e.
the region passes longphase's NG≥3 subclone-marker filter on the HP-off
side. Matches `baseline_vs_v6_counterexample.py` line 100.

AF bins
-------
  Low      : caller_af < 0.10
  Mid-Low  : 0.10 ≤ caller_af < 0.30
  Mid-High : 0.30 ≤ caller_af < 0.50
  High     : caller_af ≥ 0.50

Inputs
------
  research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/
    step1_master_four_way.tsv  (HCC1395 paired-pileup master, 4-way joined)

Outputs
-------
  phase2_completeness_audit/C3_AF_4bin_summary.tsv
  phase2_completeness_audit/figures/C3_AF_4bin_sweep.png  (1920×1080)
  assets/figures/target2/C3_AF_4bin_sweep.png             (copy)
"""

from __future__ import annotations

import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
sys.path.insert(0, str(REPO_ROOT))
from scripts.lib.plot_setup import setup_plot_style  # noqa: E402

setup_plot_style()

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------
MASTER_PATH = (
    REPO_ROOT
    / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_four_way.tsv"
)
OUT_DIR = REPO_ROOT / "phase2_completeness_audit"
FIG_DIR = OUT_DIR / "figures"
ASSETS_DIR = REPO_ROOT / "assets/figures/target2"
OUT_TSV = OUT_DIR / "C3_AF_4bin_summary.tsv"
OUT_FIG = FIG_DIR / "C3_AF_4bin_sweep.png"
ASSETS_FIG = ASSETS_DIR / "C3_AF_4bin_sweep.png"

for d in (OUT_DIR, FIG_DIR, ASSETS_DIR):
    d.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------------
# AF bin definition (4 bin)
# -------------------------------------------------------------------------
BIN_ORDER = ["Low", "Mid-Low", "Mid-High", "High"]
BIN_LABELS = {
    "Low":      "AF<0.10",
    "Mid-Low":  "0.10≤AF<0.30",
    "Mid-High": "0.30≤AF<0.50",
    "High":     "AF≥0.50",
}


def af_bin4(af: float | None) -> str | None:
    if af is None or not np.isfinite(af):
        return None
    if af < 0.10:
        return "Low"
    if af < 0.30:
        return "Mid-Low"
    if af < 0.50:
        return "Mid-High"
    return "High"


def safe_float(s: str) -> float | None:
    try:
        v = float(s)
        if not np.isfinite(v):
            return None
        return v
    except (TypeError, ValueError):
        return None


def safe_int(s: str) -> int:
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return 0


# -------------------------------------------------------------------------
# Load + bin
# -------------------------------------------------------------------------
print("=" * 78)
print("Phase 2 audit C3 — AF 4-bin sweep (baseline vs V6 markers)")
print("=" * 78)
print(f"  source : {MASTER_PATH}")

# regions[bin] = {"TP": n, "FP": n}
regions_count: dict[str, dict[str, int]] = {b: {"TP": 0, "FP": 0} for b in BIN_ORDER}
# markers[bin][bam] = {"TP": n, "FP": n}
markers: dict[str, dict[str, dict[str, int]]] = {
    b: {"baseline": {"TP": 0, "FP": 0}, "V6": {"TP": 0, "FP": 0}} for b in BIN_ORDER
}
# Per-bin AF distribution for histogram (TP/FP × baseline-marker/V6-marker/all)
af_by_bin_label_mark: dict[str, list[float]] = defaultdict(list)
n_total = 0
n_skipped_af = 0
n_skipped_bin = 0

with MASTER_PATH.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        n_total += 1
        af = safe_float(row.get("caller_af", ""))
        if af is None:
            n_skipped_af += 1
            continue
        ab = af_bin4(af)
        if ab is None:
            n_skipped_bin += 1
            continue
        label = row["label"]  # "TP" or "FP"
        if label not in ("TP", "FP"):
            continue

        regions_count[ab][label] += 1

        ng_off_base = safe_int(row.get("baseline_off_NG", "0"))
        ng_off_v6 = safe_int(row.get("V6_off_NG", "0"))
        is_marker_base = ng_off_base >= 3
        is_marker_v6 = ng_off_v6 >= 3

        if is_marker_base:
            markers[ab]["baseline"][label] += 1
            af_by_bin_label_mark[f"baseline|{label}"].append(af)
        if is_marker_v6:
            markers[ab]["V6"][label] += 1
            af_by_bin_label_mark[f"V6|{label}"].append(af)
        af_by_bin_label_mark[f"ALL|{label}"].append(af)

print(f"  scanned: {n_total:,} rows  (af-na skipped: {n_skipped_af}, bin-skip: {n_skipped_bin})")

# -------------------------------------------------------------------------
# Build summary TSV
# -------------------------------------------------------------------------
rows_out: list[dict] = []
for ab in BIN_ORDER:
    n_tp = regions_count[ab]["TP"]
    n_fp = regions_count[ab]["FP"]
    n_reg = n_tp + n_fp

    b_tp = markers[ab]["baseline"]["TP"]
    b_fp = markers[ab]["baseline"]["FP"]
    b_n = b_tp + b_fp
    b_rate = (b_tp / b_n) if b_n > 0 else float("nan")

    v_tp = markers[ab]["V6"]["TP"]
    v_fp = markers[ab]["V6"]["FP"]
    v_n = v_tp + v_fp
    v_rate = (v_tp / v_n) if v_n > 0 else float("nan")

    d_count = v_n - b_n
    d_count_pct = (d_count / b_n * 100.0) if b_n > 0 else float("nan")
    d_rate_pp = (v_rate - b_rate) * 100.0 if (np.isfinite(b_rate) and np.isfinite(v_rate)) else float("nan")
    df1_proxy = v_rate - b_rate if (np.isfinite(b_rate) and np.isfinite(v_rate)) else float("nan")

    rows_out.append({
        "AF_bin": ab,
        "AF_range": BIN_LABELS[ab],
        "n_regions": n_reg,
        "n_TP": n_tp,
        "n_FP": n_fp,
        "baseline_marker_count": b_n,
        "baseline_marker_TP": b_tp,
        "baseline_marker_FP": b_fp,
        "baseline_marker_rate": round(b_rate, 6),
        "V6_marker_count": v_n,
        "V6_marker_TP": v_tp,
        "V6_marker_FP": v_fp,
        "V6_marker_rate": round(v_rate, 6),
        "delta_marker_count": d_count,
        "delta_marker_count_pct": round(d_count_pct, 3),
        "delta_marker_rate_pp": round(d_rate_pp, 4),
        "delta_F1_proxy": round(df1_proxy, 6),
    })

df_out = pd.DataFrame(rows_out)
df_out.to_csv(OUT_TSV, sep="\t", index=False)
print(f"\nWrote: {OUT_TSV}")

print("\n=== 4-bin summary ===")
print(df_out[[
    "AF_bin", "n_regions", "baseline_marker_count", "baseline_marker_rate",
    "V6_marker_count", "V6_marker_rate", "delta_marker_count_pct",
    "delta_marker_rate_pp", "delta_F1_proxy",
]].to_string(index=False))

# -------------------------------------------------------------------------
# Figure (1920×1080)
# -------------------------------------------------------------------------
fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
gs = GridSpec(
    nrows=2, ncols=2, figure=fig,
    height_ratios=[1.0, 1.0], width_ratios=[1.05, 1.0],
    wspace=0.22, hspace=0.32,
    left=0.055, right=0.985, top=0.92, bottom=0.075,
)

ax_a = fig.add_subplot(gs[0, 0])   # Panel A grouped bar
ax_b = fig.add_subplot(gs[0, 1])   # Panel B AF histogram TP/FP
ax_c = fig.add_subplot(gs[1, :])   # Panel C ΔF1 proxy bar (full width)

bins = BIN_ORDER
bin_disp = [BIN_LABELS[b] for b in bins]
x = np.arange(len(bins))
width = 0.36

base_rates = df_out["baseline_marker_rate"].to_numpy()
v6_rates = df_out["V6_marker_rate"].to_numpy()
delta_pp = df_out["delta_marker_rate_pp"].to_numpy()
df1_proxy_arr = df_out["delta_F1_proxy"].to_numpy()
base_counts = df_out["baseline_marker_count"].to_numpy()
v6_counts = df_out["V6_marker_count"].to_numpy()
delta_count_pct = df_out["delta_marker_count_pct"].to_numpy()

# ---------- Panel A: grouped bar (marker rate baseline vs V6) ----------
color_base = "#6B7280"   # slate-gray (baseline)
color_v6 = "#2563EB"     # blue-600 (V6)

bars_b = ax_a.bar(x - width / 2, base_rates, width,
                  label="baseline marker rate", color=color_base,
                  edgecolor="#374151", linewidth=0.8)
bars_v = ax_a.bar(x + width / 2, v6_rates, width,
                  label="V6 marker rate", color=color_v6,
                  edgecolor="#1E3A8A", linewidth=0.8)

# Bar labels: rate + count
for rect, n in zip(bars_b, base_counts):
    h = rect.get_height()
    ax_a.text(rect.get_x() + rect.get_width() / 2, h + 0.012,
              f"{h:.3f}\n(n={n:,})", ha="center", va="bottom",
              fontsize=9, color="#1F2937")
for rect, n in zip(bars_v, v6_counts):
    h = rect.get_height()
    ax_a.text(rect.get_x() + rect.get_width() / 2, h + 0.012,
              f"{h:.3f}\n(n={n:,})", ha="center", va="bottom",
              fontsize=9, color="#1E3A8A")

# Δ annotation above each pair
y_top = max(np.nanmax(base_rates), np.nanmax(v6_rates))
ann_y = min(y_top + 0.16, 1.16)
for i, (bp, vp, dp, dcp) in enumerate(zip(base_rates, v6_rates, delta_pp, delta_count_pct)):
    sign = "+" if dp >= 0 else ""
    color_ann = "#16A34A" if dp > 0 else ("#DC2626" if dp < 0 else "#6B7280")
    ax_a.text(i, ann_y, f"Δrate = {sign}{dp:.2f} pp",
              ha="center", va="bottom", fontsize=10.5, fontweight="bold",
              color=color_ann)
    ax_a.text(i, ann_y - 0.045,
              f"Δcount = {dcp:+.1f}%",
              ha="center", va="bottom", fontsize=9, color="#374151")

ax_a.set_xticks(x)
ax_a.set_xticklabels(bin_disp, fontsize=10.5)
ax_a.set_ylabel("Marker precision (TP / (TP+FP) within marker set)", fontsize=11)
ax_a.set_ylim(0, 1.30)
ax_a.set_yticks(np.arange(0, 1.01, 0.2))
ax_a.grid(axis="y", linestyle=":", color="#9CA3AF", alpha=0.55)
ax_a.set_axisbelow(True)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
ax_a.set_title("Panel A — Baseline vs V6 marker precision per AF bin (HCC1395)",
               fontsize=12, fontweight="bold", pad=8)
ax_a.legend(loc="upper left", frameon=False, fontsize=10)

# ---------- Panel B: AF histogram TP vs FP, baseline-marker vs V6-marker ----------
af_tp_base = np.array(af_by_bin_label_mark["baseline|TP"])
af_fp_base = np.array(af_by_bin_label_mark["baseline|FP"])
af_tp_v6 = np.array(af_by_bin_label_mark["V6|TP"])
af_fp_v6 = np.array(af_by_bin_label_mark["V6|FP"])

af_bins_edges = np.linspace(0.0, 1.0, 41)  # 0.025 width

# 2 stacked histograms (step-style) to compare baseline vs V6
ax_b.hist(af_tp_base, bins=af_bins_edges, histtype="step", linewidth=1.8,
          color="#2563EB", linestyle="--",
          label=f"baseline marker TP (n={len(af_tp_base):,})")
ax_b.hist(af_tp_v6, bins=af_bins_edges, histtype="step", linewidth=2.0,
          color="#1D4ED8", linestyle="-",
          label=f"V6 marker TP (n={len(af_tp_v6):,})")
ax_b.hist(af_fp_base, bins=af_bins_edges, histtype="step", linewidth=1.8,
          color="#DC2626", linestyle="--",
          label=f"baseline marker FP (n={len(af_fp_base):,})")
ax_b.hist(af_fp_v6, bins=af_bins_edges, histtype="step", linewidth=2.0,
          color="#7F1D1D", linestyle="-",
          label=f"V6 marker FP (n={len(af_fp_v6):,})")

# AF bin separators (4 bin)
for bd in (0.10, 0.30, 0.50):
    ax_b.axvline(bd, color="#9CA3AF", linestyle=":", linewidth=1.0, alpha=0.7)
ax_b.text(0.05, ax_b.get_ylim()[1] * 0.95 if False else 0, "", )

ax_b.set_xlabel("caller_af", fontsize=11)
ax_b.set_ylabel("Marker count", fontsize=11)
ax_b.set_xlim(0.0, 1.0)
ax_b.grid(axis="y", linestyle=":", color="#9CA3AF", alpha=0.55)
ax_b.set_axisbelow(True)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)
ax_b.set_title("Panel B — Marker AF distribution: baseline vs V6 (TP/FP step histogram)",
               fontsize=12, fontweight="bold", pad=8)
ax_b.legend(loc="upper right", frameon=False, fontsize=9)

# Annotate bin separator labels along the top
ymax_b = ax_b.get_ylim()[1]
for xc, lab in zip([0.05, 0.20, 0.40, 0.75], bin_disp):
    ax_b.text(xc, ymax_b * 0.99, lab, ha="center", va="top",
              fontsize=9, color="#374151",
              bbox=dict(boxstyle="round,pad=0.2", facecolor="#F3F4F6",
                        edgecolor="#D1D5DB", linewidth=0.6))

# ---------- Panel C: ΔF1 proxy bar (4 bin) ----------
colors_c = ["#DC2626" if v < 0 else "#16A34A" for v in df1_proxy_arr]
bars_c = ax_c.bar(x, df1_proxy_arr, width=0.55,
                  color=colors_c, edgecolor="#1F2937", linewidth=0.9)

for rect, v, dp in zip(bars_c, df1_proxy_arr, delta_pp):
    h = rect.get_height()
    sign = "+" if v >= 0 else ""
    voffset = 0.0035 if h >= 0 else -0.0035
    va = "bottom" if h >= 0 else "top"
    ax_c.text(rect.get_x() + rect.get_width() / 2, h + voffset,
              f"ΔF1_proxy = {sign}{v:.4f}\n(Δrate = {sign}{dp:.2f} pp)",
              ha="center", va=va, fontsize=10.5, fontweight="bold",
              color="#1F2937")

ax_c.axhline(0.0, color="#374151", linewidth=1.0)

# Highlight reversal bins with explicit text
for i, v in enumerate(df1_proxy_arr):
    if v < 0:
        ax_c.text(i, ax_c.get_ylim()[0] * 0.85 if ax_c.get_ylim()[0] < 0 else -0.005,
                  "REVERSAL", ha="center", va="top", fontsize=9.5,
                  color="#7F1D1D", fontweight="bold")

ax_c.set_xticks(x)
ax_c.set_xticklabels(bin_disp, fontsize=11)
ax_c.set_ylabel("ΔF1 proxy = V6_marker_rate − baseline_marker_rate", fontsize=11)
ymin = min(np.nanmin(df1_proxy_arr) - 0.012, -0.015)
ymax = max(np.nanmax(df1_proxy_arr) + 0.012, 0.020)
ax_c.set_ylim(ymin, ymax)
ax_c.grid(axis="y", linestyle=":", color="#9CA3AF", alpha=0.55)
ax_c.set_axisbelow(True)
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.set_title("Panel C — ΔF1 proxy per AF bin (red = V6 reversal, green = V6 improvement)",
               fontsize=12, fontweight="bold", pad=8)

# Master title + footer
fig.suptitle(
    "C3 — AF 4-bin sweep · baseline vs V6 marker precision · HCC1395 paired-pileup",
    fontsize=14, fontweight="bold", y=0.975,
)
fig.text(
    0.5, 0.012,
    "Source: research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_four_way.tsv | "
    "marker = {bam}_off_NG ≥ 3 | rate = TP/(TP+FP) among markers",
    ha="center", va="bottom", fontsize=9, color="#4B5563",
)

fig.savefig(OUT_FIG, dpi=100, bbox_inches="tight")
plt.close(fig)
print(f"Wrote: {OUT_FIG}")

# Copy to assets/figures/target2/
shutil.copy2(OUT_FIG, ASSETS_FIG)
print(f"Wrote: {ASSETS_FIG}")

print("\nDone.")
