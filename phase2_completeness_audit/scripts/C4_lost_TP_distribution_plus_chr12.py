"""
C4_lost_TP_distribution_plus_chr12.py
====================================
Figure agent C4 — 350 lost TP 染色體+位置分佈 + chr12 single-chr deep dive

目的：counterexample (chr12 ΔF1=-0.0158 / 350 lost TP) 視覺化，讓 PI 知道
      V6 的 trade-off 範圍與機制。

Inputs
------
 - master TSV: research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/
               step1_master_four_way.tsv   (35,332 regions)
 - per-chr derived TSV: step1_baseline_vs_v6_per_chr.tsv

Outputs
-------
 - phase2_completeness_audit/figures/C4_lost_TP_distribution_plus_chr12.png  (1920×1080)
 - phase2_completeness_audit/C4_lost_TP_350_metadata.tsv  (350 rows × metadata)
 - assets/figures/target2/C4_lost_TP_distribution_plus_chr12.png (copy)

Definition of 350 lost TP
-------------------------
    baseline NG (off) ≥ 3      (baseline 是 marker)
  & V6 NG (off)        < 3      (V6 不是 marker)
  & label              = TP
  → 350 regions confirmed (matches step1_counterexample_summary.json q5)

Panels (1920×1080, 2×2 grid)
----------------------------
 A (top-left)  : 350 lost TP 染色體分佈 bar (24 chr)
 B (top-right) : 350 lost TP LOH zone + AF + baseline NG 三聯小圖
 C (bot-left)  : chr12 single-chr metric (marker count baseline vs V6
                  + rate + ratio + lost TP count)
 D (bot-right) : chr12 vs chr8 vs chr19 對比 (chr12=V6 反向, chr8/19=hotspot 改善)
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# CJK font setup via shared helper (feedback_matplotlib_cjk_font_rule.md)
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
PER_CHR_TSV = (
    REPO_ROOT
    / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_counterexample_per_chr.tsv"
)
OUT_DIR = REPO_ROOT / "phase2_completeness_audit"
FIG_DIR = OUT_DIR / "figures"
FIG_OUT = FIG_DIR / "C4_lost_TP_distribution_plus_chr12.png"
META_OUT = OUT_DIR / "C4_lost_TP_350_metadata.tsv"
ASSETS_COPY = REPO_ROOT / "assets/figures/target2/C4_lost_TP_distribution_plus_chr12.png"

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def derive_loh_zone(row: pd.Series) -> str:
    """Reproduce LOH_zone label from LOH_Bed_Annotation + loh_side.

    Reproduces semantics used by step1_counterexample_lost_tp_regions.tsv
    (NA when LOH_Bed_Annotation missing/NA; otherwise inner/outer/boundary
    derived from loh_side).
    """
    ann = row.get("LOH_Bed_Annotation")
    if ann is None or (isinstance(ann, float) and np.isnan(ann)) or str(ann).strip() in ("", "NA", "nan"):
        return "NA"
    side = str(row.get("loh_side", "")).strip().lower()
    if side in ("inner", "outer", "boundary"):
        return side
    return "NA"


def af_bin(af: float) -> str:
    if af is None or np.isnan(af):
        return "NA"
    if af < 0.1:
        return "0.0-0.1"
    if af < 0.2:
        return "0.1-0.2"
    if af < 0.3:
        return "0.2-0.3"
    if af < 0.5:
        return "0.3-0.5"
    if af < 0.8:
        return "0.5-0.8"
    return "0.8-1.0"


# ---------------------------------------------------------------------------
# 1. Load master TSV and filter 350 lost TP
# ---------------------------------------------------------------------------
print(f"[load] master TSV: {MASTER_TSV}")
df = pd.read_csv(MASTER_TSV, sep="\t", low_memory=False)
print(f"[load] master rows={len(df):,}  cols={len(df.columns)}")

# Filter: baseline NG (off) ≥3, V6 NG (off) <3, label=TP
mask = (
    (df["label"] == "TP")
    & (df["baseline_off_NG"] >= 3)
    & (df["V6_off_NG"] < 3)
)
lost = df.loc[mask].copy()
print(f"[filter] 350 lost TP rows actual count={len(lost)}")

# Derive LOH zone + AF bin
lost["LOH_zone"] = lost.apply(derive_loh_zone, axis=1)
lost["AF_bin"] = lost["caller_af"].apply(af_bin)

# Save metadata TSV (compact: 350 rows × key cols)
META_COLS = [
    "region_id",
    "chr",
    "pos",
    "label",
    "caller_af",
    "AF_bin",
    "LOH_Bed_Annotation",
    "LOH_Subtype",
    "loh_side",
    "LOH_zone",
    "Coverage_Category",
    "baseline_off_NG",
    "V3F_off_NG",
    "V5_off_NG",
    "V6_off_NG",
]
META_COLS = [c for c in META_COLS if c in lost.columns]
lost[META_COLS].to_csv(META_OUT, sep="\t", index=False)
print(f"[write] metadata TSV: {META_OUT}  rows={len(lost)}")


# ---------------------------------------------------------------------------
# 2. Aggregations for panels
# ---------------------------------------------------------------------------
# A: lost TP per chr count
chr_counts = (
    lost["chr"].value_counts().reindex(CHR_ORDER, fill_value=0)
)
top3_chr = chr_counts.sort_values(ascending=False).head(3)
print(f"[agg] Panel A top-3 chr: {top3_chr.to_dict()}")

# B: LOH zone, AF bin, baseline NG breakdown
loh_zone_counts = lost["LOH_zone"].value_counts()
af_bin_order = ["0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.5", "0.5-0.8", "0.8-1.0", "NA"]
af_bin_counts = lost["AF_bin"].value_counts().reindex(af_bin_order, fill_value=0)
base_ng_counts = (
    lost["baseline_off_NG"].astype(int).value_counts().sort_index()
)
print(f"[agg] LOH zone={loh_zone_counts.to_dict()}")
print(f"[agg] baseline NG={base_ng_counts.to_dict()}")

# C: chr12 single-chr deep dive (recompute from master TSV)
chr12 = df[df["chr"] == "chr12"]
chr12_base_marker = (chr12["baseline_off_NG"] >= 3).sum()
chr12_base_tp = ((chr12["baseline_off_NG"] >= 3) & (chr12["label"] == "TP")).sum()
chr12_base_fp = ((chr12["baseline_off_NG"] >= 3) & (chr12["label"] == "FP")).sum()
chr12_v6_marker = (chr12["V6_off_NG"] >= 3).sum()
chr12_v6_tp = ((chr12["V6_off_NG"] >= 3) & (chr12["label"] == "TP")).sum()
chr12_v6_fp = ((chr12["V6_off_NG"] >= 3) & (chr12["label"] == "FP")).sum()
chr12_base_rate = chr12_base_tp / max(chr12_base_marker, 1)
chr12_v6_rate = chr12_v6_tp / max(chr12_v6_marker, 1)
chr12_lost = len(lost[lost["chr"] == "chr12"])
chr12_delta = chr12_v6_rate - chr12_base_rate
print(
    f"[chr12] base marker={chr12_base_marker} TP/FP={chr12_base_tp}/{chr12_base_fp} "
    f"rate={chr12_base_rate:.4f}"
)
print(
    f"[chr12] V6   marker={chr12_v6_marker} TP/FP={chr12_v6_tp}/{chr12_v6_fp} "
    f"rate={chr12_v6_rate:.4f}  Δ={chr12_delta:+.4f}  lost TP={chr12_lost}"
)

# D: chr12 vs chr8 vs chr19 (load pre-computed per_chr rates)
per_chr = pd.read_csv(PER_CHR_TSV, sep="\t")
# coerce signed-string columns (e.g. "+0.0119") to numeric
for col in ["baseline_rate", "V6_rate", "delta_V6_minus_baseline_rate"]:
    per_chr[col] = pd.to_numeric(per_chr[col].astype(str).str.replace("+", "", regex=False), errors="coerce")
focus = per_chr[per_chr["chr"].isin(["chr8", "chr12", "chr19"])].copy()
focus = focus.set_index("chr").reindex(["chr8", "chr12", "chr19"])
print(f"[focus] chr8/12/19 per-chr table:\n{focus[['baseline_rate','V6_rate','delta_V6_minus_baseline_rate']]}")


# ---------------------------------------------------------------------------
# 3. Plot — 1920×1080
# ---------------------------------------------------------------------------
FIG_W_IN = 19.2  # 1920 px @ 100 dpi
FIG_H_IN = 10.8  # 1080 px @ 100 dpi
fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN), dpi=100)
gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.22, left=0.05, right=0.97, top=0.93, bottom=0.07)

fig.suptitle(
    "C4 — 350 lost TP 染色體分佈 + chr12 反向 deep dive  "
    "(V6 marker (NG≥3) 比 baseline 少, label=TP, n=350)",
    fontsize=14,
    fontweight="bold",
    y=0.985,
)

# ---- Panel A: lost TP per chr ----
ax_a = fig.add_subplot(gs[0, 0])
bars = ax_a.bar(
    range(len(CHR_ORDER)),
    chr_counts.values,
    color=["#d62728" if c == "chr12" else "#1f77b4" for c in CHR_ORDER],
    edgecolor="black",
    linewidth=0.4,
)
ax_a.set_xticks(range(len(CHR_ORDER)))
ax_a.set_xticklabels([c.replace("chr", "") for c in CHR_ORDER], fontsize=9)
ax_a.set_ylabel("lost TP count", fontsize=11)
ax_a.set_title("Panel A — 350 lost TP 染色體分佈 (chr12 紅色 = ΔF1 -0.0158 反向)", fontsize=12, fontweight="bold")
ax_a.grid(axis="y", linestyle="--", alpha=0.4)
# annotate top-3
top3_idx = top3_chr.index.tolist()
for i, c in enumerate(CHR_ORDER):
    v = chr_counts.iloc[i]
    if c in top3_idx and v > 0:
        rank = top3_idx.index(c) + 1
        ax_a.text(
            i,
            v + 1,
            f"#{rank}\n{int(v)}",
            ha="center",
            va="bottom",
            fontsize=8.5,
            fontweight="bold",
            color="#d62728" if c == "chr12" else "#333",
        )
    elif v > 0:
        ax_a.text(i, v + 0.5, f"{int(v)}", ha="center", va="bottom", fontsize=7.5, color="#555")
ax_a.set_ylim(0, chr_counts.max() * 1.25)

# ---- Panel B: LOH zone + AF + baseline NG ----
ax_b = fig.add_subplot(gs[0, 1])
ax_b.axis("off")
ax_b.set_title("Panel B — 350 lost TP 屬性分解 (LOH / AF / baseline NG)", fontsize=12, fontweight="bold")

# 3 sub-axes inside Panel B region
gs_b = gs[0, 1].subgridspec(1, 3, wspace=0.55)
ax_b1 = fig.add_subplot(gs_b[0, 0])
ax_b2 = fig.add_subplot(gs_b[0, 1])
ax_b3 = fig.add_subplot(gs_b[0, 2])

# B1: LOH zone bar
zone_order = ["NA", "inner", "outer", "boundary"]
zone_vals = [int(loh_zone_counts.get(z, 0)) for z in zone_order]
colors_zone = ["#9467bd", "#2ca02c", "#ff7f0e", "#8c564b"]
b1 = ax_b1.bar(zone_order, zone_vals, color=colors_zone, edgecolor="black", linewidth=0.4)
ax_b1.set_title("LOH zone", fontsize=10.5)
ax_b1.set_ylabel("count", fontsize=9.5)
ax_b1.tick_params(axis="x", labelsize=9)
for rect, v in zip(b1, zone_vals):
    pct = 100.0 * v / max(sum(zone_vals), 1)
    ax_b1.text(
        rect.get_x() + rect.get_width() / 2,
        v + 2,
        f"{v}\n({pct:.0f}%)",
        ha="center",
        va="bottom",
        fontsize=8.5,
        fontweight="bold" if zone_order[zone_vals.index(v)] == "NA" else "normal",
    )
ax_b1.set_ylim(0, max(zone_vals) * 1.28)
ax_b1.grid(axis="y", linestyle="--", alpha=0.35)

# B2: AF bin
b2 = ax_b2.bar(af_bin_counts.index, af_bin_counts.values, color="#17becf", edgecolor="black", linewidth=0.4)
ax_b2.set_title("caller_af 分布", fontsize=10.5)
ax_b2.set_ylabel("count", fontsize=9.5)
ax_b2.tick_params(axis="x", labelsize=8, rotation=35)
for rect, v in zip(b2, af_bin_counts.values):
    if v > 0:
        ax_b2.text(rect.get_x() + rect.get_width() / 2, v + 1, f"{int(v)}", ha="center", va="bottom", fontsize=8)
ax_b2.set_ylim(0, max(af_bin_counts.values) * 1.18 if max(af_bin_counts.values) > 0 else 1)
ax_b2.grid(axis="y", linestyle="--", alpha=0.35)

# B3: baseline NG distribution
ng_order = sorted(base_ng_counts.index.tolist())
ng_vals = [int(base_ng_counts.get(n, 0)) for n in ng_order]
b3 = ax_b3.bar([str(n) for n in ng_order], ng_vals, color="#bcbd22", edgecolor="black", linewidth=0.4)
ax_b3.set_title("baseline NG (這 350 個 TP)", fontsize=10.5)
ax_b3.set_xlabel("baseline NG (off)", fontsize=9.5)
ax_b3.set_ylabel("count", fontsize=9.5)
for rect, v in zip(b3, ng_vals):
    if v > 0:
        ax_b3.text(rect.get_x() + rect.get_width() / 2, v + 2, f"{int(v)}", ha="center", va="bottom", fontsize=8.5)
ax_b3.set_ylim(0, max(ng_vals) * 1.18)
ax_b3.grid(axis="y", linestyle="--", alpha=0.35)

# ---- Panel C: chr12 single-chr metric ----
ax_c = fig.add_subplot(gs[1, 0])
ax_c.set_title("Panel C — chr12 single-chr metric (baseline vs V6)", fontsize=12, fontweight="bold")

# Grouped bar: marker count (TP+FP) split TP/FP
labels = ["baseline\n(NG≥3)", "V6\n(NG≥3)"]
tp_vals = [chr12_base_tp, chr12_v6_tp]
fp_vals = [chr12_base_fp, chr12_v6_fp]
x = np.arange(len(labels))
w = 0.55
ax_c.bar(x, tp_vals, w, label="TP", color="#2ca02c", edgecolor="black", linewidth=0.4)
ax_c.bar(x, fp_vals, w, bottom=tp_vals, label="FP", color="#d62728", edgecolor="black", linewidth=0.4)
ax_c.set_xticks(x)
ax_c.set_xticklabels(labels, fontsize=10)
ax_c.set_ylabel("region count", fontsize=11)
ax_c.legend(loc="upper left", fontsize=10)
ax_c.grid(axis="y", linestyle="--", alpha=0.4)

# annotate counts + rate
for i, (tp, fp) in enumerate(zip(tp_vals, fp_vals)):
    total = tp + fp
    rate = tp / max(total, 1)
    ax_c.text(i, tp / 2, f"TP={tp}", ha="center", va="center", fontsize=9.5, color="white", fontweight="bold")
    ax_c.text(i, tp + fp / 2, f"FP={fp}", ha="center", va="center", fontsize=9.5, color="white", fontweight="bold")
    ax_c.text(
        i,
        total + max(tp_vals + fp_vals) * 0.04,
        f"total={total}\nTP rate={rate:.4f}",
        ha="center",
        va="bottom",
        fontsize=9.5,
        fontweight="bold",
    )
ax_c.set_ylim(0, max([t + f for t, f in zip(tp_vals, fp_vals)]) * 1.25)

# delta annotation
ax_c.text(
    0.5,
    0.97,
    f"ΔTP rate (V6 − baseline) = {chr12_delta:+.4f}     lost TP (baseline-only) = {chr12_lost}",
    transform=ax_c.transAxes,
    ha="center",
    va="top",
    fontsize=11,
    fontweight="bold",
    color="#d62728",
    bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff0f0", edgecolor="#d62728", linewidth=1.2),
)

# ---- Panel D: chr12 vs chr8 vs chr19 ----
ax_d = fig.add_subplot(gs[1, 1])
ax_d.set_title("Panel D — chr12 反向 vs chr8 / chr19 V6 改善 (TP rate)", fontsize=12, fontweight="bold")
focus_chrs = ["chr8", "chr12", "chr19"]
base_rates = focus["baseline_rate"].values
v6_rates = focus["V6_rate"].values
deltas = focus["delta_V6_minus_baseline_rate"].values

x = np.arange(len(focus_chrs))
w = 0.36
b_base = ax_d.bar(x - w / 2, base_rates, w, label="baseline TP rate", color="#7f7f7f", edgecolor="black", linewidth=0.4)
b_v6 = ax_d.bar(x + w / 2, v6_rates, w, label="V6 TP rate", color="#1f77b4", edgecolor="black", linewidth=0.4)

ax_d.set_xticks(x)
ax_d.set_xticklabels(focus_chrs, fontsize=11, fontweight="bold")
ax_d.set_ylabel("TP rate (marker NG≥3)", fontsize=11)
ax_d.set_ylim(0, 1.08)
ax_d.legend(loc="lower right", fontsize=10)
ax_d.grid(axis="y", linestyle="--", alpha=0.4)

for i, (b, v, d) in enumerate(zip(base_rates, v6_rates, deltas)):
    ax_d.text(i - w / 2, b + 0.012, f"{b:.4f}", ha="center", va="bottom", fontsize=8.5)
    ax_d.text(i + w / 2, v + 0.012, f"{v:.4f}", ha="center", va="bottom", fontsize=8.5)
    color = "#d62728" if d < 0 else "#2ca02c"
    arrow = "↓" if d < 0 else "↑"
    ax_d.text(
        i,
        max(b, v) + 0.07,
        f"Δ={d:+.4f} {arrow}",
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
        color=color,
    )

# context footnote
fig.text(
    0.5,
    0.015,
    f"Source: step1_master_four_way.tsv ({len(df):,} regions)  ·  "
    f"350 lost TP filter: baseline_off_NG≥3 & V6_off_NG<3 & label=TP  ·  "
    f"chr12 ΔF1=-0.0158 (TP rate 0.9445→0.9287)  ·  "
    f"chr8 ΔF1=+0.2111 (V6 hotspot 改善)  ·  chr19 ΔF1=-0.0026 (微負)",
    ha="center",
    va="bottom",
    fontsize=9,
    color="#555",
)

# ---------------------------------------------------------------------------
# 4. Save
# ---------------------------------------------------------------------------
FIG_DIR.mkdir(parents=True, exist_ok=True)
ASSETS_COPY.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(FIG_OUT, dpi=100, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"[write] figure: {FIG_OUT}")

shutil.copy2(FIG_OUT, ASSETS_COPY)
print(f"[copy ] {ASSETS_COPY}")

# ---------------------------------------------------------------------------
# 5. Status summary print
# ---------------------------------------------------------------------------
print("\n=== STATUS SUMMARY ===")
print(f"350 lost TP top-3 chr: {[(c, int(v)) for c, v in top3_chr.items()]}")
print(
    "chr12 vs hotspot: chr12 Δ={:+.4f} (lost {} TP), chr8 Δ={:+.4f} (hotspot), chr19 Δ={:+.4f}".format(
        chr12_delta, chr12_lost, focus.loc["chr8", "delta_V6_minus_baseline_rate"], focus.loc["chr19", "delta_V6_minus_baseline_rate"]
    )
)
print(f"LOH zone breakdown: {loh_zone_counts.to_dict()}")
