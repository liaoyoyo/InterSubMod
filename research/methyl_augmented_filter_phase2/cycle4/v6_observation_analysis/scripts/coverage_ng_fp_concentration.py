#!/usr/bin/env python3
"""User 2026-05-20 質問：Coverage_Multiple + V6_NG 是否能用於 FP 移除？

Pre-reg hypothesis:
  H_COV: Coverage_Multiple AUC≈0.5 + KS 0.43-0.48 偏高 → 形狀不同非單調
         若 FP 集中在 X 軸中間（如 cov 0.8-1.2 plateau），則 bin-gated filter 仍可移除部分 FP
         即使 global AUC 0.5
  H_NG1: V6_NG=1 是否 FP-enriched？(NG=1 是 marker filter NG≥3 切掉的最低 NG bucket)

Method: 不用 density (mass-loss for FP small samples); 用 absolute count + FP fraction
  per bin = FP / (TP + FP)。

  - Coverage_Multiple_imp: 30 bin in [0.5, 2.5]
  - V6_off_NG: 7 bin {0, 1, 2, 3, 4, 5, 6+}

Output: figures/09_coverage_fp_concentration.png + figures/10_v6_ng_fp_concentration.png
        + data/coverage_ng_per_bin_fp_fraction.tsv
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
CYCLE1_MASTER = REPO / "research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv"
CYCLE2_DATA = REPO / "research/methyl_augmented_filter_phase2/cycle2/data"
OUT = REPO / "research/methyl_augmented_filter_phase2/cycle4/v6_observation_analysis"
FIG = OUT / "figures"
DATA = OUT / "data"

SAMPLE_ORDER = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]
SAMPLE_COLORS = {"HCC1395": "#D97757", "HCC1937": "#15803D", "HCC1954": "#A16207",
                 "H1437": "#1E3A8A", "H2009": "#6B7280"}


def load_sample(sample: str) -> pd.DataFrame:
    if sample == "HCC1395":
        df = pd.read_csv(CYCLE1_MASTER, sep="\t", low_memory=False)
    else:
        df = pd.read_csv(CYCLE2_DATA / f"{sample}_master_augmented.tsv", sep="\t", low_memory=False)
    if "y" not in df.columns and "label" in df.columns:
        df["y"] = (df["label"] == "TP").astype(int)
    if "Coverage_Multiple_imp" not in df.columns and "Coverage_Multiple" in df.columns:
        df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(df["Coverage_Multiple"].median())
    return df


def coverage_fp_concentration(sample_dfs):
    """Figure 09: Coverage_Multiple TP/FP absolute count + FP fraction per bin."""
    fig, axes = plt.subplots(2, 5, figsize=(22, 9))

    bin_edges = np.linspace(0.5, 2.5, 31)  # 30 bins
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    all_rows = []
    for col_idx, sample in enumerate(SAMPLE_ORDER):
        df = sample_dfs[sample]
        col = "Coverage_Multiple_imp"
        if col not in df.columns:
            continue
        sub = df[[col, "y"]].dropna()
        tp_vals = sub.loc[sub["y"]==1, col].values
        fp_vals = sub.loc[sub["y"]==0, col].values

        tp_hist, _ = np.histogram(tp_vals, bins=bin_edges)
        fp_hist, _ = np.histogram(fp_vals, bins=bin_edges)

        # Row 1: Absolute count, same y-scale per sample (TP+FP combined max)
        ax = axes[0, col_idx]
        width = bin_edges[1] - bin_edges[0]
        ax.bar(bin_centers, tp_hist, width=width * 0.85, alpha=0.55,
               color="steelblue", label=f"TP n={len(tp_vals):,}", edgecolor="black", linewidth=0.3)
        ax.bar(bin_centers, fp_hist, width=width * 0.85, alpha=0.85,
               color="#C2410C", label=f"FP n={len(fp_vals):,}", edgecolor="black", linewidth=0.3)
        ax.set_title(f"{sample} — Coverage TP vs FP (count)", fontsize=10, fontweight="bold")
        ax.set_xlabel("Coverage_Multiple")
        ax.set_ylabel("count")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(axis="y", alpha=0.3)

        # Row 2: FP fraction per bin = FP / (TP + FP)
        ax2 = axes[1, col_idx]
        total = tp_hist + fp_hist
        with np.errstate(divide="ignore", invalid="ignore"):
            fp_frac = np.where(total > 0, fp_hist / total, np.nan)
        baseline_fp_frac = len(fp_vals) / (len(tp_vals) + len(fp_vals))

        # Color by enrichment (above/below baseline FP rate)
        bar_colors = ["#C2410C" if (not np.isnan(f) and f > baseline_fp_frac) else "#15803D"
                      for f in fp_frac]
        ax2.bar(bin_centers, fp_frac, width=width * 0.85, color=bar_colors,
                edgecolor="black", linewidth=0.3, alpha=0.85)
        ax2.axhline(baseline_fp_frac, color="black", linestyle="--", linewidth=1.5,
                    label=f"sample-wide FP rate {baseline_fp_frac:.3f}")
        ax2.set_xlabel("Coverage_Multiple")
        ax2.set_ylabel("FP fraction (FP/(TP+FP))")
        ax2.set_ylim(0, min(1.0, max(0.3, np.nanmax(fp_frac) * 1.1)))
        ax2.set_title(f"{sample} — FP fraction (紅=enriched, 綠=低於 baseline)",
                      fontsize=10, fontweight="bold")
        ax2.legend(fontsize=7)
        ax2.grid(axis="y", alpha=0.3)

        # Record per-bin TSV
        for i, c in enumerate(bin_centers):
            all_rows.append({
                "sample": sample, "feature": "Coverage_Multiple_imp",
                "bin_center": c, "TP_count": int(tp_hist[i]), "FP_count": int(fp_hist[i]),
                "FP_fraction": float(fp_frac[i]) if not np.isnan(fp_frac[i]) else None,
                "baseline_FP_fraction": baseline_fp_frac,
                "enrichment": float(fp_frac[i] / baseline_fp_frac) if not np.isnan(fp_frac[i]) and baseline_fp_frac > 0 else None,
            })

    fig.suptitle("Figure 9 — Coverage_Multiple TP/FP per-bin count + FP fraction\n"
                 "上排：absolute count（紅=FP, 藍=TP）／下排：FP/(TP+FP) per bin（紅=enriched, 綠=低於 baseline FP rate）",
                 fontsize=12, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(FIG / "09_coverage_fp_concentration.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ 09_coverage_fp_concentration.png")
    return all_rows


def ng_fp_concentration(sample_dfs):
    """Figure 10: V6_off_NG TP/FP per NG-bucket count + FP fraction."""
    fig, axes = plt.subplots(2, 5, figsize=(22, 9))

    ng_buckets = [0, 1, 2, 3, 4, 5, 6, 7, 8]  # 9 buckets

    all_rows = []
    for col_idx, sample in enumerate(SAMPLE_ORDER):
        df = sample_dfs[sample]
        col = "V6_off_NG"
        if col not in df.columns:
            continue
        sub = df[[col, "y"]].dropna()
        tp_vals = sub.loc[sub["y"]==1, col].astype(int).values
        fp_vals = sub.loc[sub["y"]==0, col].astype(int).values

        tp_counts = np.array([np.sum(tp_vals == k) for k in ng_buckets])
        fp_counts = np.array([np.sum(fp_vals == k) for k in ng_buckets])

        # Row 1: Absolute count
        ax = axes[0, col_idx]
        x = np.arange(len(ng_buckets))
        width = 0.38
        ax.bar(x - width/2, tp_counts, width, color="steelblue",
               label=f"TP n={len(tp_vals):,}", edgecolor="black", linewidth=0.3)
        ax.bar(x + width/2, fp_counts, width, color="#C2410C",
               label=f"FP n={len(fp_vals):,}", edgecolor="black", linewidth=0.3)
        ax.set_xticks(x)
        ax.set_xticklabels(ng_buckets)
        ax.set_xlabel("V6_off_NG")
        ax.set_ylabel("count")
        ax.set_title(f"{sample} — NG TP vs FP (count)", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7)
        ax.grid(axis="y", alpha=0.3)

        # Row 2: FP fraction per NG bucket
        ax2 = axes[1, col_idx]
        total = tp_counts + fp_counts
        with np.errstate(divide="ignore", invalid="ignore"):
            fp_frac = np.where(total > 0, fp_counts / total, np.nan)
        baseline_fp_frac = len(fp_vals) / (len(tp_vals) + len(fp_vals))

        bar_colors = ["#C2410C" if (not np.isnan(f) and f > baseline_fp_frac) else "#15803D"
                      for f in fp_frac]
        bars = ax2.bar(x, fp_frac, width=0.75, color=bar_colors, edgecolor="black", linewidth=0.3, alpha=0.85)
        ax2.axhline(baseline_fp_frac, color="black", linestyle="--", linewidth=1.5,
                    label=f"sample-wide FP rate {baseline_fp_frac:.3f}")
        for i, (bar, f) in enumerate(zip(bars, fp_frac)):
            if not np.isnan(f) and total[i] > 0:
                ax2.text(bar.get_x() + bar.get_width()/2, f + 0.01,
                         f"{f:.3f}\nn={int(total[i])}", ha="center", fontsize=7.5)
        ax2.set_xticks(x)
        ax2.set_xticklabels(ng_buckets)
        ax2.set_xlabel("V6_off_NG")
        ax2.set_ylabel("FP fraction (FP/(TP+FP))")
        ax2.set_ylim(0, min(1.0, max(0.5, np.nanmax(fp_frac) * 1.15)))
        ax2.set_title(f"{sample} — NG bucket FP fraction", fontsize=10, fontweight="bold")
        ax2.legend(fontsize=7)
        ax2.grid(axis="y", alpha=0.3)

        # Record
        for i, k in enumerate(ng_buckets):
            all_rows.append({
                "sample": sample, "feature": "V6_off_NG", "NG_value": k,
                "TP_count": int(tp_counts[i]), "FP_count": int(fp_counts[i]),
                "FP_fraction": float(fp_frac[i]) if not np.isnan(fp_frac[i]) else None,
                "baseline_FP_fraction": baseline_fp_frac,
                "enrichment": float(fp_frac[i] / baseline_fp_frac) if not np.isnan(fp_frac[i]) and baseline_fp_frac > 0 else None,
            })

    fig.suptitle("Figure 10 — V6_off_NG per-bucket count + FP fraction\n"
                 "上排：absolute count／下排：FP/(TP+FP) per NG bucket（紅=enriched, 綠=低 baseline FP rate）",
                 fontsize=12, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(FIG / "10_v6_ng_fp_concentration.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ 10_v6_ng_fp_concentration.png")
    return all_rows


if __name__ == "__main__":
    print("Loading 5 sample master TSVs ...")
    sample_dfs = {s: load_sample(s) for s in SAMPLE_ORDER}
    for s, df in sample_dfs.items():
        print(f"  {s}: rows={len(df):,} TP={int((df['y']==1).sum()):,} FP={int((df['y']==0).sum()):,}")

    print("\n=== Figure 09 — Coverage_Multiple FP concentration ===")
    cov_rows = coverage_fp_concentration(sample_dfs)

    print("\n=== Figure 10 — V6_off_NG FP concentration ===")
    ng_rows = ng_fp_concentration(sample_dfs)

    # Save TSV
    all_rows = cov_rows + ng_rows
    out_tsv = DATA / "coverage_ng_per_bin_fp_fraction.tsv"
    pd.DataFrame(all_rows).to_csv(out_tsv, sep="\t", index=False)
    print(f"\n  ✓ {out_tsv}")

    # Quick verdict text
    print("\n=== QUICK VERDICTS ===")
    df_all = pd.DataFrame(all_rows)
    print("\n[Coverage_Multiple] FP enrichment 最高 bin (per sample):")
    for s in SAMPLE_ORDER:
        sub = df_all[(df_all["sample"]==s) & (df_all["feature"]=="Coverage_Multiple_imp") & df_all["enrichment"].notna()]
        if len(sub) > 0:
            top = sub.nlargest(3, "enrichment")
            for _, r in top.iterrows():
                print(f"  {s} bin={r['bin_center']:.2f}: TP={r['TP_count']} FP={r['FP_count']} "
                      f"FP_frac={r['FP_fraction']:.3f} ({r['enrichment']:.2f}× baseline)")

    print("\n[V6_off_NG] FP fraction per NG bucket:")
    for s in SAMPLE_ORDER:
        sub = df_all[(df_all["sample"]==s) & (df_all["feature"]=="V6_off_NG")]
        for _, r in sub.iterrows():
            if r["TP_count"] + r["FP_count"] > 0:
                marker = " ⚠" if (r["enrichment"] is not None and r["enrichment"] > 1.5) else ""
                print(f"  {s} NG={int(r['NG_value'])}: TP={r['TP_count']:>6} FP={r['FP_count']:>5} "
                      f"FP_frac={r['FP_fraction'] if r['FP_fraction'] is not None else 'NA':.3f} "
                      f"enrich={r['enrichment'] if r['enrichment'] is not None else 'NA':.2f}×{marker}")
