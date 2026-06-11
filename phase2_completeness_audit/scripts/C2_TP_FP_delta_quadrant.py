"""
C2_TP_FP_delta_quadrant.py
===========================
Figure agent C2 — V6 vs baseline TP/FP delta 拆分 (TP gained / TP lost / FP removed / FP added)

目的：補 D3 dimension — 不只看 aggregate marker rate，要看 V6 對每類 region 的「實際 delta」。
驗證宣告：V6 net TP:FP trade ratio ~14:1。

Definitions
-----------
對每 region (35,332):
  baseline_is_marker = (baseline_off_NG ≥ 3)
  V6_is_marker       = (V6_off_NG ≥ 3)

4 mutually-exclusive transition classes (cross-tab on label ∈ {TP, FP}):
  TP gained  = (not baseline_marker) & (V6_marker)     & (label == TP)
  TP lost    = (baseline_marker)     & (not V6_marker) & (label == TP)
  FP removed = (baseline_marker)     & (not V6_marker) & (label == FP)
  FP added   = (not baseline_marker) & (V6_marker)     & (label == FP)

Net deltas:
  Net TP = TP_gained - TP_lost     (V6 對 TP 的淨貢獻)
  Net FP = FP_added  - FP_removed  (V6 對 FP 的淨代價)
  Trade ratio = Net TP / Net FP     (越大代表 trade-off 越划算)

Expected (counterexample_summary.json q1_marker_membership):
  TP gained = 8043 / TP lost = 350 / FP removed = 11 / FP added = 560
  Net TP = 7693 / Net FP = 549 / Trade ratio ≈ 14.0 : 1

Inputs
------
 - master TSV: research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/
               step1_master_four_way.tsv  (35,332 regions)

Outputs
-------
 - phase2_completeness_audit/scripts/C2_TP_FP_delta_quadrant.py     (this script)
 - phase2_completeness_audit/figures/C2_TP_FP_delta_quadrant.png    (1920×1080)
 - phase2_completeness_audit/C2_TP_FP_delta_4class.tsv              (4-row summary)
 - assets/figures/target2/C2_TP_FP_delta_quadrant.png               (copy)

Panels (1920×1080)
-------------------
 A (top-left, wide):  4 個 vertical bar (TP gained / TP lost / FP removed / FP added)
                       色彩編碼: 綠 = V6 對 ground truth 有利、紅 = 不利
                       (TP gained 綠 / TP lost 紅 / FP removed 綠 / FP added 紅)
 B (top-right):       net TP vs net FP grouped bar (含 trade ratio annotation)
 C (bottom-left):     TP gained vs TP lost paired comparison
 D (bottom-right):    FP removed vs FP added paired comparison
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# CJK font setup
# ---------------------------------------------------------------------------
REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
sys.path.insert(0, str(REPO_ROOT))
from scripts.lib.plot_setup import setup_plot_style  # noqa: E402

setup_plot_style(base_size=11, dpi=120)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
MASTER_TSV = (
    REPO_ROOT
    / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_four_way.tsv"
)
OUT_DIR = REPO_ROOT / "phase2_completeness_audit"
FIG_DIR = OUT_DIR / "figures"
FIG_OUT = FIG_DIR / "C2_TP_FP_delta_quadrant.png"
TSV_OUT = OUT_DIR / "C2_TP_FP_delta_4class.tsv"
ASSETS_COPY = REPO_ROOT / "assets/figures/target2/C2_TP_FP_delta_quadrant.png"

FIG_DIR.mkdir(parents=True, exist_ok=True)
ASSETS_COPY.parent.mkdir(parents=True, exist_ok=True)

# Color scheme — 綠/紅 follow `favorable / unfavorable to ground truth`
COLOR_FAV = "#2E7D32"   # green — TP gained, FP removed
COLOR_UNFAV = "#C62828"  # red   — TP lost, FP added
EDGE_FAV = "#1B5E20"
EDGE_UNFAV = "#7F1D1D"

# ---------------------------------------------------------------------------
# Load + compute 4-class cross-tab
# ---------------------------------------------------------------------------
print("=" * 72)
print("C2 — V6 vs baseline TP/FP delta quadrant (4-class cross-tab)")
print("=" * 72)
print(f"Loading master TSV: {MASTER_TSV.name} ...")

usecols = ["region_id", "chr", "pos", "label", "baseline_off_NG", "V6_off_NG"]
df = pd.read_csv(MASTER_TSV, sep="\t", usecols=usecols, low_memory=False)
print(f"  rows = {len(df):,}")

df["baseline_off_NG"] = pd.to_numeric(df["baseline_off_NG"], errors="coerce").fillna(0).astype(int)
df["V6_off_NG"] = pd.to_numeric(df["V6_off_NG"], errors="coerce").fillna(0).astype(int)

baseline_is_marker = df["baseline_off_NG"] >= 3
v6_is_marker = df["V6_off_NG"] >= 3
is_tp = df["label"] == "TP"
is_fp = df["label"] == "FP"

mask_tp_gained = (~baseline_is_marker) & v6_is_marker & is_tp
mask_tp_lost = baseline_is_marker & (~v6_is_marker) & is_tp
mask_fp_removed = baseline_is_marker & (~v6_is_marker) & is_fp
mask_fp_added = (~baseline_is_marker) & v6_is_marker & is_fp

n_tp_gained = int(mask_tp_gained.sum())
n_tp_lost = int(mask_tp_lost.sum())
n_fp_removed = int(mask_fp_removed.sum())
n_fp_added = int(mask_fp_added.sum())

# Denominators for rates (label-conditioned)
n_tp_total = int(is_tp.sum())
n_fp_total = int(is_fp.sum())
n_baseline_marker_tp = int((baseline_is_marker & is_tp).sum())
n_baseline_marker_fp = int((baseline_is_marker & is_fp).sum())
n_v6_marker_tp = int((v6_is_marker & is_tp).sum())
n_v6_marker_fp = int((v6_is_marker & is_fp).sum())
n_both_marker_tp = int((baseline_is_marker & v6_is_marker & is_tp).sum())
n_both_marker_fp = int((baseline_is_marker & v6_is_marker & is_fp).sum())

net_tp = n_tp_gained - n_tp_lost
net_fp = n_fp_added - n_fp_removed
trade_ratio = (net_tp / net_fp) if net_fp != 0 else float("inf")

print()
print("Counts:")
print(f"  TP gained  (baseline-, V6+, TP) = {n_tp_gained:>6,d}")
print(f"  TP lost    (baseline+, V6-, TP) = {n_tp_lost:>6,d}")
print(f"  FP removed (baseline+, V6-, FP) = {n_fp_removed:>6,d}")
print(f"  FP added   (baseline-, V6+, FP) = {n_fp_added:>6,d}")
print(f"  Net TP     = {net_tp:>6,d}")
print(f"  Net FP     = {net_fp:>6,d}")
print(f"  Trade ratio (Net TP : Net FP) = {trade_ratio:.2f} : 1")
print(f"  TP total = {n_tp_total:,}  FP total = {n_fp_total:,}")

# ---------------------------------------------------------------------------
# Persist TSV
# ---------------------------------------------------------------------------
tsv_rows = [
    {
        "class": "TP_gained",
        "definition": "baseline_NG<3 & V6_NG>=3 & label=TP",
        "count": n_tp_gained,
        "label_total": n_tp_total,
        "rate_of_label_total_pct": round(100.0 * n_tp_gained / n_tp_total, 4) if n_tp_total else np.nan,
        "favorable_to_v6": True,
    },
    {
        "class": "TP_lost",
        "definition": "baseline_NG>=3 & V6_NG<3 & label=TP",
        "count": n_tp_lost,
        "label_total": n_tp_total,
        "rate_of_label_total_pct": round(100.0 * n_tp_lost / n_tp_total, 4) if n_tp_total else np.nan,
        "favorable_to_v6": False,
    },
    {
        "class": "FP_removed",
        "definition": "baseline_NG>=3 & V6_NG<3 & label=FP",
        "count": n_fp_removed,
        "label_total": n_fp_total,
        "rate_of_label_total_pct": round(100.0 * n_fp_removed / n_fp_total, 4) if n_fp_total else np.nan,
        "favorable_to_v6": True,
    },
    {
        "class": "FP_added",
        "definition": "baseline_NG<3 & V6_NG>=3 & label=FP",
        "count": n_fp_added,
        "label_total": n_fp_total,
        "rate_of_label_total_pct": round(100.0 * n_fp_added / n_fp_total, 4) if n_fp_total else np.nan,
        "favorable_to_v6": False,
    },
    {
        "class": "Net_TP",
        "definition": "TP_gained - TP_lost",
        "count": net_tp,
        "label_total": n_tp_total,
        "rate_of_label_total_pct": round(100.0 * net_tp / n_tp_total, 4) if n_tp_total else np.nan,
        "favorable_to_v6": net_tp > 0,
    },
    {
        "class": "Net_FP",
        "definition": "FP_added - FP_removed",
        "count": net_fp,
        "label_total": n_fp_total,
        "rate_of_label_total_pct": round(100.0 * net_fp / n_fp_total, 4) if n_fp_total else np.nan,
        "favorable_to_v6": net_fp < 0,
    },
    {
        "class": "Trade_Ratio",
        "definition": "Net_TP / Net_FP",
        "count": round(trade_ratio, 4),
        "label_total": None,
        "rate_of_label_total_pct": None,
        "favorable_to_v6": trade_ratio > 1.0,
    },
]
tsv_df = pd.DataFrame(tsv_rows)
tsv_df.to_csv(TSV_OUT, sep="\t", index=False)
print(f"\nWrote TSV: {TSV_OUT}")

# ---------------------------------------------------------------------------
# Plot 1920×1080
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(16, 9), dpi=120)
gs = GridSpec(
    nrows=2,
    ncols=2,
    figure=fig,
    width_ratios=[1.25, 1.0],
    height_ratios=[1.0, 1.0],
    hspace=0.42,
    wspace=0.28,
    left=0.07,
    right=0.97,
    top=0.90,
    bottom=0.09,
)

fig.suptitle(
    "C2 — V6 vs baseline: TP/FP transition delta (4-class cross-tab)\n"
    f"Net TP = {net_tp:,}  ·  Net FP = {net_fp:,}  ·  Trade ratio = {trade_ratio:.2f} : 1",
    fontsize=15,
    fontweight="bold",
    y=0.975,
)

# ------------- Panel A (top-left): 4-class vertical bars -------------------
axA = fig.add_subplot(gs[0, 0])
labels_a = ["TP gained\n(救回)", "TP lost\n(漏失)", "FP removed\n(去除)", "FP added\n(引入)"]
counts_a = [n_tp_gained, n_tp_lost, n_fp_removed, n_fp_added]
colors_a = [COLOR_FAV, COLOR_UNFAV, COLOR_FAV, COLOR_UNFAV]
edges_a = [EDGE_FAV, EDGE_UNFAV, EDGE_FAV, EDGE_UNFAV]

x_a = np.arange(len(labels_a))
bars_a = axA.bar(
    x_a, counts_a, width=0.6, color=colors_a, edgecolor=edges_a, linewidth=1.0
)
axA.set_xticks(x_a)
axA.set_xticklabels(labels_a, fontsize=10.5)
axA.set_ylabel("Region count", fontsize=11)
axA.set_title("A · 4-class transition counts (35,332 regions)", fontsize=12, pad=10)
ymax_a = max(counts_a) * 1.18
axA.set_ylim(0, ymax_a)
axA.grid(axis="y", linestyle=":", color="#9CA3AF", alpha=0.5)
axA.set_axisbelow(True)
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)

# Bar labels: count + pct of label total
for rect, c, lbl_total in zip(
    bars_a,
    counts_a,
    [n_tp_total, n_tp_total, n_fp_total, n_fp_total],
):
    pct = 100.0 * c / lbl_total if lbl_total else 0.0
    h = rect.get_height()
    axA.text(
        rect.get_x() + rect.get_width() / 2,
        h + ymax_a * 0.015,
        f"{c:,}\n({pct:.2f}%)",
        ha="center",
        va="bottom",
        fontsize=10,
        color="#111827",
        fontweight="bold",
    )

# Legend strip
from matplotlib.patches import Patch
legend_a = [
    Patch(facecolor=COLOR_FAV, edgecolor=EDGE_FAV, label="V6 對 ground truth 有利"),
    Patch(facecolor=COLOR_UNFAV, edgecolor=EDGE_UNFAV, label="V6 對 ground truth 不利"),
]
axA.legend(handles=legend_a, loc="upper right", fontsize=9.5, frameon=False)

# ------------- Panel B (top-right): Net TP vs Net FP -----------------------
axB = fig.add_subplot(gs[0, 1])
labels_b = ["Net TP\n(V6 淨貢獻)", "Net FP\n(V6 淨代價)"]
counts_b = [net_tp, net_fp]
colors_b = [COLOR_FAV if net_tp > 0 else COLOR_UNFAV,
            COLOR_UNFAV if net_fp > 0 else COLOR_FAV]
edges_b = [EDGE_FAV if net_tp > 0 else EDGE_UNFAV,
           EDGE_UNFAV if net_fp > 0 else EDGE_FAV]
x_b = np.arange(len(labels_b))
bars_b = axB.bar(x_b, counts_b, width=0.5, color=colors_b, edgecolor=edges_b, linewidth=1.0)
axB.set_xticks(x_b)
axB.set_xticklabels(labels_b, fontsize=10.5)
axB.set_ylabel("Net region count", fontsize=11)
axB.set_title("B · Net delta + trade ratio", fontsize=12, pad=10)
ymax_b = max(counts_b) * 1.25 if max(counts_b) > 0 else 100
axB.set_ylim(0, ymax_b)
axB.grid(axis="y", linestyle=":", color="#9CA3AF", alpha=0.5)
axB.set_axisbelow(True)
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
for rect, c in zip(bars_b, counts_b):
    h = rect.get_height()
    axB.text(
        rect.get_x() + rect.get_width() / 2,
        h + ymax_b * 0.015,
        f"{c:,}",
        ha="center",
        va="bottom",
        fontsize=11,
        color="#111827",
        fontweight="bold",
    )

# Trade ratio annotation centered between bars
axB.annotate(
    f"Net TP : Net FP\n= {net_tp:,} : {net_fp:,}\n≈ {trade_ratio:.2f} : 1",
    xy=(0.5, 0.65),
    xycoords="axes fraction",
    ha="center",
    va="center",
    fontsize=13,
    fontweight="bold",
    color="#0F172A",
    bbox=dict(
        boxstyle="round,pad=0.55",
        facecolor="#FEF3C7",
        edgecolor="#D97706",
        linewidth=1.2,
    ),
)

# ------------- Panel C (bottom-left): TP gained vs TP lost -----------------
axC = fig.add_subplot(gs[1, 0])
labels_c = ["TP gained", "TP lost"]
counts_c = [n_tp_gained, n_tp_lost]
colors_c = [COLOR_FAV, COLOR_UNFAV]
edges_c = [EDGE_FAV, EDGE_UNFAV]
x_c = np.arange(len(labels_c))
bars_c = axC.bar(x_c, counts_c, width=0.5, color=colors_c, edgecolor=edges_c, linewidth=1.0)
axC.set_xticks(x_c)
axC.set_xticklabels(labels_c, fontsize=10.5)
axC.set_ylabel("Region count (TP only)", fontsize=11)
ratio_tp = (n_tp_gained / n_tp_lost) if n_tp_lost else float("inf")
axC.set_title(
    f"C · TP gained vs TP lost (ratio = {ratio_tp:.1f} : 1)", fontsize=12, pad=10
)
ymax_c = max(counts_c) * 1.18
axC.set_ylim(0, ymax_c)
axC.grid(axis="y", linestyle=":", color="#9CA3AF", alpha=0.5)
axC.set_axisbelow(True)
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)
for rect, c in zip(bars_c, counts_c):
    pct = 100.0 * c / n_tp_total
    h = rect.get_height()
    axC.text(
        rect.get_x() + rect.get_width() / 2,
        h + ymax_c * 0.015,
        f"{c:,}\n({pct:.2f}% of {n_tp_total:,} TP)",
        ha="center",
        va="bottom",
        fontsize=9.5,
        color="#111827",
        fontweight="bold",
    )

# ------------- Panel D (bottom-right): FP removed vs FP added --------------
axD = fig.add_subplot(gs[1, 1])
labels_d = ["FP removed", "FP added"]
counts_d = [n_fp_removed, n_fp_added]
colors_d = [COLOR_FAV, COLOR_UNFAV]
edges_d = [EDGE_FAV, EDGE_UNFAV]
x_d = np.arange(len(labels_d))
bars_d = axD.bar(x_d, counts_d, width=0.5, color=colors_d, edgecolor=edges_d, linewidth=1.0)
axD.set_xticks(x_d)
axD.set_xticklabels(labels_d, fontsize=10.5)
axD.set_ylabel("Region count (FP only)", fontsize=11)
ratio_fp = (n_fp_added / n_fp_removed) if n_fp_removed else float("inf")
axD.set_title(
    f"D · FP removed vs FP added (added : removed = {ratio_fp:.1f} : 1)",
    fontsize=12,
    pad=10,
)
ymax_d = max(counts_d) * 1.22
axD.set_ylim(0, ymax_d)
axD.grid(axis="y", linestyle=":", color="#9CA3AF", alpha=0.5)
axD.set_axisbelow(True)
axD.spines["top"].set_visible(False)
axD.spines["right"].set_visible(False)
for rect, c in zip(bars_d, counts_d):
    pct = 100.0 * c / n_fp_total
    h = rect.get_height()
    axD.text(
        rect.get_x() + rect.get_width() / 2,
        h + ymax_d * 0.015,
        f"{c:,}\n({pct:.2f}% of {n_fp_total:,} FP)",
        ha="center",
        va="bottom",
        fontsize=9.5,
        color="#111827",
        fontweight="bold",
    )

# Footer: provenance / definition reminder
fig.text(
    0.5,
    0.025,
    "Source: step1_master_four_way.tsv (35,332 regions)  ·  Marker rule: NG_off ≥ 3  ·  "
    "Trade ratio = (TP_gained − TP_lost) / (FP_added − FP_removed)",
    ha="center",
    va="bottom",
    fontsize=9.5,
    color="#374151",
)

fig.savefig(FIG_OUT, dpi=120, bbox_inches="tight")
plt.close(fig)
print(f"Wrote figure: {FIG_OUT}")

# ---------------------------------------------------------------------------
# Copy to assets/figures/target2/
# ---------------------------------------------------------------------------
shutil.copy2(FIG_OUT, ASSETS_COPY)
print(f"Copied to: {ASSETS_COPY}")

print("\nDone.")
