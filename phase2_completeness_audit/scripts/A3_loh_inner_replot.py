"""
A3_loh_inner_replot.py
======================
Re-visualize `loh_inner_flag` (binary {0, 1}) as a grouped bar chart of
TP-inner% vs FP-inner% across 5 samples. The earlier A3 boxplot grid is
misleading for this feature because boxplot median = 0 for both classes
(a single flat line) even though the signal exists in *prevalence*, not
distribution shape. Effect size is preserved via AUC + Cohen's d (which
the boxplot grid couldn't show visually).

Inputs
------
* HCC1395 master  : research/methyl_augmented_filter_phase2/cycle2/data/per_bam_master_V6.tsv
* HCC1937 master  : research/.../HCC1937_master_augmented.tsv
* HCC1954 master  : research/.../HCC1954_master_augmented.tsv
* H1437   master  : research/.../H1437_master_augmented.tsv
* H2009   master  : research/.../H2009_master_augmented.tsv
* AUC/Cohen's d   : phase2_completeness_audit/A3_per_feature_AUC_cohend_5sample.tsv

Outputs
-------
* phase2_completeness_audit/A3_loh_inner_proportion_5sample.tsv
* phase2_completeness_audit/figures/A3_loh_inner_flag_grouped_bar_5sample.png

Rationale
---------
* Binary 0/1 -> proportion (% positive in each class) is the natural summary.
* Gap (pp) = TP% - FP% is the visual analogue of Cohen's d sign/magnitude
  for this binary feature.
* 5/5 samples sign-consistent (TP > FP) -> ⭐3 anchor (per Phase 2 A3 review).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
sys.path.insert(0, str(REPO_ROOT))
from scripts.lib.plot_setup import setup_plot_style  # noqa: E402

setup_plot_style()

OUT_TSV_DIR = REPO_ROOT / "phase2_completeness_audit"
OUT_FIG_DIR = OUT_TSV_DIR / "figures"
OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = REPO_ROOT / "research/methyl_augmented_filter_phase2/cycle2/data"
A3_AUC_TSV = OUT_TSV_DIR / "A3_per_feature_AUC_cohend_5sample.tsv"

SAMPLES = [
    ("HCC1395", DATA_DIR / "per_bam_master_V6.tsv"),
    ("HCC1937", DATA_DIR / "HCC1937_master_augmented.tsv"),
    ("HCC1954", DATA_DIR / "HCC1954_master_augmented.tsv"),
    ("H1437",   DATA_DIR / "H1437_master_augmented.tsv"),
    ("H2009",   DATA_DIR / "H2009_master_augmented.tsv"),
]

# -------------------------------------------------------------------------
# Load each sample, count TP/FP and loh_inner_flag==1 within each class
# -------------------------------------------------------------------------
print("=" * 72)
print("Phase A3 replot — loh_inner_flag grouped bar (5 samples)")
print("=" * 72)

# Load AUC/Cohen's d reference from existing A3 TSV
auc_df = pd.read_csv(A3_AUC_TSV, sep="\t")
loh_stats = auc_df[auc_df["feature"] == "loh_inner_flag"].set_index("sample")[["AUC", "Cohens_d"]]
print(f"\nAUC/Cohen's d source rows for loh_inner_flag: {len(loh_stats)}")

rows = []
for sample, path in SAMPLES:
    print(f"  loading {sample} from {path.name} ...")
    if not path.exists():
        raise FileNotFoundError(f"Missing source TSV: {path}")
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in ("label", "loh_inner_flag"),
                     low_memory=False)
    if "loh_inner_flag" not in df.columns:
        raise KeyError(f"{sample}: loh_inner_flag column not found in {path}")
    if "label" not in df.columns:
        raise KeyError(f"{sample}: label column not found in {path}")

    df["loh_inner_flag"] = pd.to_numeric(df["loh_inner_flag"], errors="coerce").fillna(0).astype(int)

    tp_mask = df["label"] == "TP"
    fp_mask = df["label"] == "FP"
    n_tp = int(tp_mask.sum())
    n_fp = int(fp_mask.sum())
    n_tp_inner = int((tp_mask & (df["loh_inner_flag"] == 1)).sum())
    n_fp_inner = int((fp_mask & (df["loh_inner_flag"] == 1)).sum())

    tp_pct = (n_tp_inner / n_tp * 100.0) if n_tp > 0 else np.nan
    fp_pct = (n_fp_inner / n_fp * 100.0) if n_fp > 0 else np.nan
    gap_pp = tp_pct - fp_pct

    auc_val = float(loh_stats.loc[sample, "AUC"]) if sample in loh_stats.index else np.nan
    d_val = float(loh_stats.loc[sample, "Cohens_d"]) if sample in loh_stats.index else np.nan

    rows.append({
        "sample": sample,
        "n_TP": n_tp,
        "n_FP": n_fp,
        "n_TP_inner": n_tp_inner,
        "n_FP_inner": n_fp_inner,
        "TP_inner_pct": round(tp_pct, 3),
        "FP_inner_pct": round(fp_pct, 3),
        "Gap_pp": round(gap_pp, 3),
        "AUC": round(auc_val, 4),
        "Cohens_d": round(d_val, 4),
    })
    print(f"    -> TP={n_tp:>6d} ({n_tp_inner:>5d} inner = {tp_pct:5.2f}%)  "
          f"FP={n_fp:>5d} ({n_fp_inner:>4d} inner = {fp_pct:5.2f}%)  Gap={gap_pp:+5.2f}pp")

out_df = pd.DataFrame(rows)
out_tsv = OUT_TSV_DIR / "A3_loh_inner_proportion_5sample.tsv"
out_df.to_csv(out_tsv, sep="\t", index=False)
print(f"\nWrote: {out_tsv}")

# -------------------------------------------------------------------------
# Plot grouped bar: TP_inner_pct (blue) vs FP_inner_pct (red)
# -------------------------------------------------------------------------
samples = out_df["sample"].tolist()
tp_pcts = out_df["TP_inner_pct"].to_numpy()
fp_pcts = out_df["FP_inner_pct"].to_numpy()
gaps    = out_df["Gap_pp"].to_numpy()
aucs    = out_df["AUC"].to_numpy()
ds      = out_df["Cohens_d"].to_numpy()

x = np.arange(len(samples))
width = 0.36

fig, ax = plt.subplots(figsize=(11.5, 6.0), dpi=150)

bars_tp = ax.bar(x - width / 2, tp_pcts, width,
                 label="TP — inner%", color="#2563EB", edgecolor="#1E3A8A", linewidth=0.8)
bars_fp = ax.bar(x + width / 2, fp_pcts, width,
                 label="FP — inner%", color="#DC2626", edgecolor="#7F1D1D", linewidth=0.8)

# Value labels on top of each bar
for rect in list(bars_tp) + list(bars_fp):
    h = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, h + 0.6, f"{h:.1f}%",
            ha="center", va="bottom", fontsize=9, color="#111827")

# Gap annotation + AUC + Cohen's d sub-line above each sample group
y_top = max(tp_pcts.max(), fp_pcts.max())
gap_y = y_top + 6.5
for i, (g, a, d) in enumerate(zip(gaps, aucs, ds)):
    ax.text(i, gap_y, f"Gap = {g:+.1f} pp",
            ha="center", va="bottom", fontsize=10, fontweight="bold", color="#111827")
    ax.text(i, gap_y - 2.6, f"AUC={a:.2f}  d={d:+.2f}",
            ha="center", va="bottom", fontsize=8.5, color="#374151")

ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("% of class located in LOH inner zone", fontsize=11)
y_max = max(65.0, gap_y + 8.0)
ax.set_ylim(0, y_max)
ax.set_yticks(np.arange(0, int(y_max) + 1, 10))
ax.grid(axis="y", linestyle=":", color="#9CA3AF", alpha=0.55)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Place legend below title to avoid overlap with right-most Gap annotation
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2,
          frameon=False, fontsize=10.5)

fig.suptitle("loh_inner_flag — TP vs FP proportion in LOH inner zone (binary feature)",
             fontsize=13, fontweight="bold", y=0.985)
ax.set_title("5/5 sign-consistent (TP > FP all samples); AUC 0.54–0.72; boxplot 不適用 binary feature",
             fontsize=10.5, color="#374151", pad=12)

fig.tight_layout(rect=(0, 0.04, 1, 0.97))

out_png = OUT_FIG_DIR / "A3_loh_inner_flag_grouped_bar_5sample.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Wrote: {out_png}")

print("\nDone.")
