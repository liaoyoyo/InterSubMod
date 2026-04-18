#!/usr/bin/env python3
"""G2: Germline Signature Profiling — PON membership, GT, H flag analysis.

Reads G1 output (full features for PASS variants + lightweight all-variants detail).
Produces 10 figures and data TSVs characterizing germline-like patterns in TO FP.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    OUTPUT_ROOT, SAMPLE_ORDER, setup_plot_style, save_figure, add_caption,
    compute_auc, compute_effect_size, ensure_dir,
    COLOR_TP, COLOR_FP, TRUTH_PALETTE,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

WORKSPACE = OUTPUT_ROOT / "20260401_germline_fp_identification"
G1_DIR = WORKSPACE / "G1_vcf_features"
G2_DIR = ensure_dir(WORKSPACE / "G2_germline_signatures")
FIG_DIR = ensure_dir(G2_DIR / "figures")
DATA_DIR = ensure_dir(G2_DIR / "data")


def _print(msg: str):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_pass_features() -> pd.DataFrame:
    """Load G1 merged PASS features (full feature set)."""
    path = G1_DIR / "all_samples_merged.tsv.gz"
    _print(f"[G2] Loading PASS features from {path.name}")
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    _print(f"  -> {len(df):,} rows x {len(df.columns)} cols")
    return df


def load_all_variants_lite(sample: str) -> pd.DataFrame:
    """Load lightweight all-variants detail for one sample."""
    path = G1_DIR / f"{sample}_all_variants_lite.tsv.gz"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)


def load_pon_summaries() -> pd.DataFrame:
    """Load and merge all per-sample PoN summaries."""
    frames = []
    for sample in SAMPLE_ORDER:
        path = G1_DIR / f"{sample}_pon_summary.tsv"
        if path.exists():
            frames.append(pd.read_csv(path, sep="\t"))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Analysis functions
# ---------------------------------------------------------------------------


def analyze_pon_membership(df_pass: pd.DataFrame, df_lite_all: pd.DataFrame):
    """Analyze PoN hit count distribution across groups."""
    _print("\n[G2] === PoN Membership Analysis ===")

    # --- Fig01: PoN count distribution (PASS FP vs TP) ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: PASS variants
    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()
    for ax_idx, (label, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        sub = pass_data[pass_data["truth_label"] == label]
        counts = sub["pon_count"].value_counts().sort_index()
        axes[0].bar(counts.index + (ax_idx * 0.35 - 0.175), counts.values,
                     width=0.35, color=color, label=label, alpha=0.8)
    axes[0].set_xlabel("PoN Hit Count")
    axes[0].set_ylabel("Number of Variants")
    axes[0].set_title("PASS Variants: PoN Count Distribution")
    axes[0].legend()
    axes[0].set_xticks(range(5))

    # Right: all variants (from lite data)
    if not df_lite_all.empty:
        for ax_idx, (label, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
            sub = df_lite_all[df_lite_all["truth_label"] == label]
            counts = sub["pon_count"].value_counts().sort_index()
            total = len(sub)
            axes[1].bar(counts.index + (ax_idx * 0.35 - 0.175),
                        counts.values / total * 100,
                        width=0.35, color=color, label=label, alpha=0.8)
        axes[1].set_xlabel("PoN Hit Count")
        axes[1].set_ylabel("Fraction (%)")
        axes[1].set_title("ALL Variants: PoN Count Fraction")
        axes[1].legend()
        axes[1].set_xticks(range(5))

    fig.suptitle("G2-Fig01: PoN Hit Count Distribution", fontsize=14, fontweight="bold")
    add_caption(fig, "Source: G1 VCF features | 7 samples pooled")
    save_figure(fig, FIG_DIR / "G2_fig01_pon_count_distribution.png")


def analyze_gt_distribution(df_pass: pd.DataFrame):
    """Analyze genotype distribution."""
    _print("\n[G2] === GT Distribution Analysis ===")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()

    # --- Fig03: GT distribution per sample ---
    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    axes = axes.flatten()

    summary_rows = []
    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        sub = pass_data[pass_data["sample"] == sample]
        if sub.empty:
            ax.set_visible(False)
            continue

        for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
            grp = sub[sub["truth_label"] == label]
            n_het = (grp["is_het"] == 1).sum()
            n_hom = (grp["is_hom_alt"] == 1).sum()
            n_total = len(grp)
            het_frac = n_het / n_total * 100 if n_total > 0 else 0
            hom_frac = n_hom / n_total * 100 if n_total > 0 else 0
            summary_rows.append({
                "sample": sample, "truth_label": label,
                "n_het": n_het, "n_hom_alt": n_hom, "n_total": n_total,
                "het_frac": het_frac, "hom_frac": hom_frac,
            })
            ax.bar([f"{label}\nhet", f"{label}\nhom"],
                   [het_frac, hom_frac], color=color, alpha=0.7, width=0.6)

        ax.set_title(sample, fontsize=10)
        ax.set_ylabel("% of variants")
        ax.set_ylim(0, 100)

    if len(SAMPLE_ORDER) < 8:
        axes[-1].set_visible(False)

    fig.suptitle("G2-Fig03: GT Distribution (het vs hom) — PASS TP vs FP", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G2_fig03_gt_distribution_per_sample.png")

    # Save summary
    pd.DataFrame(summary_rows).to_csv(DATA_DIR / "gt_distribution.tsv", sep="\t", index=False)
    _print(f"  [G2] GT distribution summary saved")

    return pd.DataFrame(summary_rows)


def analyze_h_flag(df_pass: pd.DataFrame):
    """Analyze H flag prevalence."""
    _print("\n[G2] === H Flag Analysis ===")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()

    # --- Fig05: H flag prevalence per sample ---
    fig, ax = plt.subplots(figsize=(12, 5))
    x_positions = np.arange(len(SAMPLE_ORDER))
    width = 0.35

    summary_rows = []
    for offset, (label, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        h_rates = []
        for sample in SAMPLE_ORDER:
            sub = pass_data[(pass_data["sample"] == sample) & (pass_data["truth_label"] == label)]
            if len(sub) == 0:
                h_rates.append(0)
                continue
            rate = (sub["h_flag"] == 1).mean() * 100
            h_rates.append(rate)
            summary_rows.append({
                "sample": sample, "truth_label": label,
                "h_flag_count": (sub["h_flag"] == 1).sum(),
                "total": len(sub), "h_flag_rate": rate,
            })
        ax.bar(x_positions + offset * width - width / 2, h_rates,
               width=width, color=color, label=label, alpha=0.8)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right")
    ax.set_ylabel("H Flag Prevalence (%)")
    ax.set_title("G2-Fig05: H Flag Prevalence (Single Haplotype Support)", fontsize=13, fontweight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    add_caption(fig, "H flag = variant found only in one haplotype | Higher H → more likely somatic")
    save_figure(fig, FIG_DIR / "G2_fig05_h_flag_prevalence.png")

    pd.DataFrame(summary_rows).to_csv(DATA_DIR / "h_flag_prevalence.tsv", sep="\t", index=False)


def analyze_af_by_gt(df_pass: pd.DataFrame):
    """Analyze AF density by GT type."""
    _print("\n[G2] === AF × GT Analysis ===")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()
    pass_data["af"] = pd.to_numeric(pass_data["af"], errors="coerce")

    # --- Fig04: AF density by GT type ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax_idx, (label, title) in enumerate([("TP", "True Positives"), ("FP", "False Positives")]):
        ax = axes[ax_idx]
        sub = pass_data[pass_data["truth_label"] == label]
        for gt_label, gt_col, color, ls in [
            ("het (0/1)", "is_het", "#2196F3", "-"),
            ("hom-alt (1/1)", "is_hom_alt", "#F44336", "--"),
        ]:
            gt_sub = sub[sub[gt_col] == 1]["af"].dropna()
            if len(gt_sub) > 10 and gt_sub.std() > 1e-6:
                gt_sub.plot.kde(ax=ax, label=f"{gt_label} (n={len(gt_sub):,})",
                                color=color, linestyle=ls, bw_method=0.05)
            elif len(gt_sub) > 0:
                ax.axvline(gt_sub.mean(), color=color, linestyle=ls, alpha=0.7,
                          label=f"{gt_label} (n={len(gt_sub):,}, const)")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Allele Fraction")
        ax.set_ylabel("Density")
        ax.set_title(title)
        ax.legend(fontsize=8)

    fig.suptitle("G2-Fig04: AF Distribution by Genotype Type", fontsize=14, fontweight="bold")
    add_caption(fig, "PASS variants only | Germline het expected at AF~0.5, hom at AF~1.0")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G2_fig04_af_density_by_gt.png")


def analyze_pon_per_sample_heatmap(df_lite_all: pd.DataFrame):
    """Heatmap of PoN membership per sample."""
    _print("\n[G2] === Per-Sample PoN Heatmap ===")

    # --- Fig02: Per-sample PoN heatmap ---
    rows = []
    for sample in SAMPLE_ORDER:
        sub = df_lite_all[df_lite_all["sample"] == sample] if "sample" in df_lite_all.columns else pd.DataFrame()
        if sub.empty:
            continue
        for label in ["TP", "FP"]:
            grp = sub[sub["truth_label"] == label]
            n = len(grp)
            if n == 0:
                continue
            for pc in range(5):
                count = (grp["pon_count"] == pc).sum()
                rows.append({
                    "sample": sample, "truth_label": label,
                    "pon_count": pc, "fraction": count / n * 100,
                })

    if not rows:
        _print("  [G2] No data for PoN heatmap")
        return

    df_heat = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax_idx, label in enumerate(["TP", "FP"]):
        ax = axes[ax_idx]
        sub = df_heat[df_heat["truth_label"] == label]
        pivot = sub.pivot_table(index="sample", columns="pon_count", values="fraction", fill_value=0)
        pivot = pivot.reindex(SAMPLE_ORDER).fillna(0)
        sns.heatmap(pivot, annot=True, fmt=".1f", cmap="YlOrRd", ax=ax, vmin=0, vmax=100)
        ax.set_title(f"{label} — PoN Count Fraction (%)")
        ax.set_ylabel("")

    fig.suptitle("G2-Fig02: Per-Sample PoN Membership Heatmap", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G2_fig02_pon_per_sample_heatmap.png")


def analyze_combined_germline_roc(df_pass: pd.DataFrame):
    """ROC analysis for combined germline indicators."""
    _print("\n[G2] === Combined Germline ROC ===")
    from sklearn.metrics import roc_curve, roc_auc_score

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()
    pass_data["af"] = pd.to_numeric(pass_data["af"], errors="coerce")
    pass_data = pass_data.dropna(subset=["af"])

    y_true = (pass_data["truth_label"] == "FP").astype(int).values

    # Individual feature scores (higher = more FP-like / germline-like)
    scores = {}

    # 1. H_flag absent (germline has both haplotypes → H=0)
    scores["H_flag_absent"] = (1 - pass_data["h_flag"].astype(float)).values

    # 2. hom_alt (GT=1/1 → germline hom)
    scores["is_hom_alt"] = pass_data["is_hom_alt"].astype(float).values

    # 3. AF proximity to 0.5 or 1.0 (germline modes)
    af = pass_data["af"].values
    scores["AF_germline_mode"] = np.maximum(
        1 - np.abs(af - 0.5) * 4,      # peak at 0.5
        1 - np.abs(af - 1.0) * 4,      # peak at 1.0
    ).clip(0, 1)

    # 4. Combined score
    scores["Combined"] = (
        scores["H_flag_absent"] * 0.3 +
        scores["is_hom_alt"] * 0.3 +
        scores["AF_germline_mode"] * 0.4
    )

    # --- Fig07: ROC for PoN_count (not really applicable for PASS variants) ---
    # --- Fig09: Combined germline indicator ROC ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Individual features
    ax = axes[0]
    for name, score in scores.items():
        if name == "Combined":
            continue
        auc = compute_auc(y_true, score)
        fpr, tpr, _ = roc_curve(y_true, score)
        ax.plot(fpr, tpr, label=f"{name} (AUC={auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Individual Features")
    ax.legend(fontsize=8)

    # Right: Combined
    ax = axes[1]
    auc_combined = compute_auc(y_true, scores["Combined"])
    fpr, tpr, _ = roc_curve(y_true, scores["Combined"])
    ax.plot(fpr, tpr, color="#E65100", linewidth=2,
            label=f"Combined (AUC={auc_combined:.3f})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Combined Germline Score")
    ax.legend()

    fig.suptitle("G2-Fig09: Germline Indicator ROC — PASS Variants",
                 fontsize=14, fontweight="bold")
    add_caption(fig, "Target: predict FP (germline) among PASS variants | 7 samples pooled")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G2_fig09_combined_germline_roc.png")

    # Save AUC summary
    auc_rows = []
    for name, score in scores.items():
        auc_rows.append({"feature": name, "auc": compute_auc(y_true, score)})
    pd.DataFrame(auc_rows).to_csv(DATA_DIR / "germline_indicator_auc.tsv", sep="\t", index=False)
    _print(f"  [G2] Combined AUC = {auc_combined:.4f}")


def analyze_per_sample_consistency(df_pass: pd.DataFrame):
    """Per-sample consistency of germline indicators."""
    _print("\n[G2] === Per-Sample Consistency ===")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()
    pass_data["af"] = pd.to_numeric(pass_data["af"], errors="coerce")

    metrics = {
        "H_flag_rate": lambda g: (g["h_flag"] == 1).mean(),
        "hom_alt_rate": lambda g: (g["is_hom_alt"] == 1).mean(),
        "het_rate": lambda g: (g["is_het"] == 1).mean(),
        "mean_AF": lambda g: g["af"].mean(),
    }

    rows = []
    for sample in SAMPLE_ORDER:
        for label in ["TP", "FP"]:
            sub = pass_data[(pass_data["sample"] == sample) & (pass_data["truth_label"] == label)]
            if sub.empty:
                continue
            row = {"sample": sample, "truth_label": label, "n": len(sub)}
            for mname, fn in metrics.items():
                try:
                    row[mname] = fn(sub)
                except Exception:
                    row[mname] = float("nan")
            rows.append(row)

    df_metrics = pd.DataFrame(rows)

    # --- Fig10: Forest plot of FP - TP delta per sample ---
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 6))
    if len(metrics) == 1:
        axes = [axes]

    for ax_idx, mname in enumerate(metrics):
        ax = axes[ax_idx]
        deltas = []
        for sample in SAMPLE_ORDER:
            tp_val = df_metrics.loc[(df_metrics["sample"] == sample) & (df_metrics["truth_label"] == "TP"), mname]
            fp_val = df_metrics.loc[(df_metrics["sample"] == sample) & (df_metrics["truth_label"] == "FP"), mname]
            if not tp_val.empty and not fp_val.empty:
                delta = fp_val.values[0] - tp_val.values[0]
                deltas.append(delta)
                ax.barh(sample, delta, color=COLOR_FP if delta > 0 else COLOR_TP, alpha=0.7)
            else:
                deltas.append(0)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_title(mname, fontsize=10)
        ax.set_xlabel("FP - TP delta")

    fig.suptitle("G2-Fig10: Per-Sample FP-TP Delta (Consistency Check)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G2_fig10_per_sample_consistency.png")

    df_metrics.to_csv(DATA_DIR / "per_sample_metrics.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    setup_plot_style()
    _print("=" * 70)
    _print("G2: Germline Signature Profiling")
    _print("=" * 70)

    # Load data
    df_pass = load_pass_features()

    # Load all-variants lite data (for PoN analysis)
    lite_frames = []
    for sample in SAMPLE_ORDER:
        df_lite = load_all_variants_lite(sample)
        if not df_lite.empty:
            df_lite["sample"] = sample
            lite_frames.append(df_lite)
    df_lite_all = pd.concat(lite_frames, ignore_index=True) if lite_frames else pd.DataFrame()
    _print(f"  Loaded {len(df_lite_all):,} all-variants-lite rows")

    # Run analyses
    analyze_pon_membership(df_pass, df_lite_all)
    analyze_pon_per_sample_heatmap(df_lite_all)
    analyze_gt_distribution(df_pass)
    analyze_af_by_gt(df_pass)
    analyze_h_flag(df_pass)
    analyze_combined_germline_roc(df_pass)
    analyze_per_sample_consistency(df_pass)

    _print(f"\n[G2] DONE. Figures: {FIG_DIR}")
    _print(f"[G2] Data: {DATA_DIR}")


if __name__ == "__main__":
    main()
