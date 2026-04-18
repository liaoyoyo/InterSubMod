#!/usr/bin/env python3
"""G3: Allele Fraction Architecture — AF distribution, GMM fit, AF×GT interaction.

Reads G1 output to characterize AF patterns that distinguish germline FP from somatic TP.
Produces 8 figures.
"""

from __future__ import annotations

import sys
from pathlib import Path

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
    compute_auc, ensure_dir,
    COLOR_TP, COLOR_FP,
)

WORKSPACE = OUTPUT_ROOT / "20260401_germline_fp_identification"
G1_DIR = WORKSPACE / "G1_vcf_features"
G3_DIR = ensure_dir(WORKSPACE / "G3_af_architecture")
FIG_DIR = ensure_dir(G3_DIR / "figures")
DATA_DIR = ensure_dir(G3_DIR / "data")


def _print(msg: str):
    print(msg, flush=True)


def load_pass_features() -> pd.DataFrame:
    path = G1_DIR / "all_samples_merged.tsv.gz"
    _print(f"[G3] Loading from {path.name}")
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    df["af"] = pd.to_numeric(df["af"], errors="coerce")
    _print(f"  -> {len(df):,} rows")
    return df


def load_all_lite() -> pd.DataFrame:
    """Load all-variants-lite for AF of PoN-filtered variants.

    Only loads the columns we need and samples down to keep memory manageable.
    """
    needed_cols = ["variant_key", "truth_label", "is_pass", "af", "pon_count", "gt", "h_flag"]
    frames = []
    for sample in SAMPLE_ORDER:
        path = G1_DIR / f"{sample}_all_variants_lite.tsv.gz"
        if path.exists():
            _print(f"  Loading lite: {sample}")
            df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False,
                            usecols=lambda c: c in needed_cols)
            df["sample"] = sample
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df["af"] = pd.to_numeric(df["af"], errors="coerce")
    _print(f"  Loaded {len(df):,} all-variants-lite rows")
    # Sample down for KDE (keep all FP and TP, sample others)
    if len(df) > 500_000:
        labeled = df[df["truth_label"].isin(["TP", "FP"])]
        unlabeled = df[~df["truth_label"].isin(["TP", "FP"])]
        if len(unlabeled) > 500_000:
            unlabeled = unlabeled.sample(500_000, random_state=42)
        df = pd.concat([labeled, unlabeled], ignore_index=True)
        _print(f"  After sampling: {len(df):,} rows")
    return df


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def fig01_af_density_overlay(df_pass: pd.DataFrame, df_lite: pd.DataFrame):
    """AF density: PON-filtered FP vs residual (PASS) FP vs TP."""
    _print("\n[G3] Fig01: AF density overlay")

    fig, ax = plt.subplots(figsize=(10, 5))

    # 1. PASS TP
    tp = df_pass[(df_pass["filter_is_pass"] == 1) & (df_pass["truth_label"] == "TP")]["af"].dropna()
    tp.plot.kde(ax=ax, label=f"PASS TP (n={len(tp):,})", color=COLOR_TP, linewidth=2, bw_method=0.03)

    # 2. PASS FP (residual)
    fp_pass = df_pass[(df_pass["filter_is_pass"] == 1) & (df_pass["truth_label"] == "FP")]["af"].dropna()
    fp_pass.plot.kde(ax=ax, label=f"PASS FP / residual (n={len(fp_pass):,})", color=COLOR_FP, linewidth=2, bw_method=0.03)

    # 3. PoN-filtered FP (from lite)
    if not df_lite.empty:
        pon_fp = df_lite[(df_lite["truth_label"] == "FP") & (df_lite["is_pass"] == 0)]["af"].dropna()
        n_pon_fp = len(pon_fp)
        if n_pon_fp > 100:
            # Sample for KDE to avoid slow computation
            if n_pon_fp > 50_000:
                pon_fp = pon_fp.sample(50_000, random_state=42)
            pon_fp.plot.kde(ax=ax, label=f"PoN-filtered FP (n={n_pon_fp:,})",
                           color="#FF9800", linewidth=1.5, linestyle="--", bw_method=0.03)

    ax.set_xlim(0, 1)
    ax.set_xlabel("Allele Fraction")
    ax.set_ylabel("Density")
    ax.legend(fontsize=9)
    ax.set_title("G3-Fig01: AF Distribution — PASS TP vs PASS FP vs PoN-Filtered FP",
                 fontsize=13, fontweight="bold")
    add_caption(fig, "All 7 samples pooled | Germline het→AF~0.5, hom→AF~1.0")
    save_figure(fig, FIG_DIR / "G3_fig01_af_density_overlay.png")


def fig02_af_gmm_fit(df_pass: pd.DataFrame):
    """Gaussian mixture model fit on residual FP AF."""
    _print("\n[G3] Fig02: AF GMM fit")
    from sklearn.mixture import GaussianMixture

    fp_pass = df_pass[(df_pass["filter_is_pass"] == 1) & (df_pass["truth_label"] == "FP")]["af"].dropna()
    tp_pass = df_pass[(df_pass["filter_is_pass"] == 1) & (df_pass["truth_label"] == "TP")]["af"].dropna()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax_idx, (data, label, color) in enumerate([
        (fp_pass, "PASS FP (residual)", COLOR_FP),
        (tp_pass, "PASS TP", COLOR_TP),
    ]):
        ax = axes[ax_idx]
        X = data.values.reshape(-1, 1)

        # Fit 2-component GMM
        gmm = GaussianMixture(n_components=2, random_state=42).fit(X)
        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
        weights = gmm.weights_

        # Histogram
        ax.hist(data, bins=50, density=True, alpha=0.5, color=color, edgecolor="white")

        # GMM components
        x_range = np.linspace(0, 1, 200).reshape(-1, 1)
        for i in range(2):
            comp = weights[i] * sp_stats.norm.pdf(x_range.flatten(), means[i], stds[i])
            ax.plot(x_range, comp, "--", linewidth=1.5,
                    label=f"Component {i+1}: μ={means[i]:.3f}, σ={stds[i]:.3f}, w={weights[i]:.2f}")

        # Total density
        total = np.exp(gmm.score_samples(x_range))
        ax.plot(x_range, total, "k-", linewidth=2, label="GMM total")

        ax.set_xlim(0, 1)
        ax.set_xlabel("Allele Fraction")
        ax.set_ylabel("Density")
        ax.set_title(f"{label} (n={len(data):,})")
        ax.legend(fontsize=7)

    fig.suptitle("G3-Fig02: Gaussian Mixture Model Fit on AF Distribution",
                 fontsize=14, fontweight="bold")
    add_caption(fig, "2-component GMM | Germline FP expected: mode at 0.5 (het) + 1.0 (hom)")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G3_fig02_af_gmm_fit.png")


def fig03_af_gt_joint(df_pass: pd.DataFrame):
    """AF × GT joint distribution."""
    _print("\n[G3] Fig03: AF × GT joint plot")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()
    pass_data["gt_type"] = "other"
    pass_data.loc[pass_data["is_het"] == 1, "gt_type"] = "het (0/1)"
    pass_data.loc[pass_data["is_hom_alt"] == 1, "gt_type"] = "hom (1/1)"

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax_idx, (label, color_main) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        ax = axes[ax_idx]
        sub = pass_data[pass_data["truth_label"] == label]
        for gt, marker, c in [("het (0/1)", "o", "#42A5F5"), ("hom (1/1)", "s", "#EF5350")]:
            gt_sub = sub[sub["gt_type"] == gt]
            ax.scatter(gt_sub["af"], np.random.normal(0, 0.1, len(gt_sub)),
                      alpha=0.05, s=2, c=c, label=f"{gt} (n={len(gt_sub):,})")
        # AF marginal
        ax2 = ax.twinx()
        sub["af"].dropna().plot.kde(ax=ax2, color=color_main, linewidth=2, bw_method=0.03)
        ax2.set_ylabel("")
        ax2.set_yticks([])

        ax.set_xlim(0, 1)
        ax.set_xlabel("Allele Fraction")
        ax.set_ylabel("GT type (jittered)")
        ax.set_title(f"{label} — AF × GT")
        ax.legend(fontsize=8, loc="upper left")

    fig.suptitle("G3-Fig03: AF × GT 2D Distribution", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G3_fig03_af_gt_joint.png")


def fig04_fp_rate_by_af_bin(df_pass: pd.DataFrame):
    """FP rate by AF bin, per sample."""
    _print("\n[G3] Fig04: FP rate by AF bin")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()
    pass_data["af_bin"] = pd.cut(pass_data["af"], bins=20, labels=False)
    bin_edges = np.linspace(0, 1, 21)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    fig, ax = plt.subplots(figsize=(12, 6))

    for sample in SAMPLE_ORDER:
        sub = pass_data[pass_data["sample"] == sample]
        if sub.empty:
            continue
        rates = []
        for b in range(20):
            bin_sub = sub[sub["af_bin"] == b]
            if len(bin_sub) < 5:
                rates.append(float("nan"))
            else:
                rates.append((bin_sub["truth_label"] == "FP").mean() * 100)
        ax.plot(bin_centers, rates, "o-", label=sample, alpha=0.7, markersize=3)

    ax.set_xlabel("Allele Fraction")
    ax.set_ylabel("FP Rate (%)")
    ax.set_title("G3-Fig04: FP Rate by AF Bin (PASS Variants)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 1)
    add_caption(fig, "Higher FP rate at AF~0.5 and AF>0.8 → germline het and hom modes")
    save_figure(fig, FIG_DIR / "G3_fig04_fp_rate_by_af_bin.png")


def fig05_af_threshold_sweep(df_pass: pd.DataFrame):
    """Precision-recall at various AF thresholds for germline FP identification."""
    _print("\n[G3] Fig05: AF threshold sweep PR")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()
    y_true = (pass_data["truth_label"] == "FP").astype(int).values
    af = pass_data["af"].values

    # Score: closer to 0.5 or 1.0 → more germline-like
    germline_score = np.maximum(
        1 - np.abs(af - 0.5) * 4,
        1 - np.abs(af - 1.0) * 4,
    ).clip(0, 1)

    from sklearn.metrics import precision_recall_curve, roc_curve
    precision, recall, thresholds_pr = precision_recall_curve(y_true, germline_score)
    fpr, tpr, thresholds_roc = roc_curve(y_true, germline_score)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # PR curve
    ax = axes[0]
    ax.plot(recall, precision, color=COLOR_FP, linewidth=2)
    ax.set_xlabel("Recall (FP found)")
    ax.set_ylabel("Precision (of flagged)")
    ax.set_title("Precision-Recall")
    ax.grid(alpha=0.3)

    # ROC curve
    ax = axes[1]
    auc = compute_auc(y_true, germline_score)
    ax.plot(fpr, tpr, color=COLOR_FP, linewidth=2, label=f"AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC")
    ax.legend()

    fig.suptitle("G3-Fig05: AF-Based Germline Score — PR & ROC",
                 fontsize=14, fontweight="bold")
    add_caption(fig, "Score = max(proximity to AF=0.5, proximity to AF=1.0)")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G3_fig05_af_threshold_sweep.png")

    _print(f"  [G3] AF germline score AUC = {auc:.4f}")


def fig06_per_sample_af(df_pass: pd.DataFrame):
    """Per-sample AF distribution comparison."""
    _print("\n[G3] Fig06: Per-sample AF distributions")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()

    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    axes = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        sub = pass_data[pass_data["sample"] == sample]
        for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
            grp = sub[sub["truth_label"] == label]["af"].dropna()
            if len(grp) > 10:
                grp.plot.kde(ax=ax, label=f"{label} (n={len(grp):,})",
                            color=color, bw_method=0.03, linewidth=1.5)
        ax.set_xlim(0, 1)
        ax.set_title(sample, fontsize=10)
        ax.set_xlabel("AF")
        ax.legend(fontsize=7)

    if len(SAMPLE_ORDER) < 8:
        axes[-1].set_visible(False)

    fig.suptitle("G3-Fig06: Per-Sample AF Distribution (PASS TP vs FP)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G3_fig06_per_sample_af.png")


def fig07_af_roc(df_pass: pd.DataFrame):
    """ROC: various AF-derived features."""
    _print("\n[G3] Fig07: AF-based ROC")
    from sklearn.metrics import roc_curve

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()
    y_true = (pass_data["truth_label"] == "FP").astype(int).values
    af = pass_data["af"].values

    scores = {
        "AF (raw)": af,
        "1-AF (inverted)": 1 - af,
        "|AF-0.5|": np.abs(af - 0.5),
        "Germline mode score": np.maximum(
            1 - np.abs(af - 0.5) * 4,
            1 - np.abs(af - 1.0) * 4,
        ).clip(0, 1),
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    for name, score in scores.items():
        mask = np.isfinite(score)
        auc = compute_auc(y_true[mask], score[mask])
        fpr, tpr, _ = roc_curve(y_true[mask], score[mask])
        ax.plot(fpr, tpr, label=f"{name} (AUC={auc:.3f})", linewidth=1.5)

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("G3-Fig07: AF-Based Feature ROC", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    add_caption(fig, "PASS variants, 7 samples pooled")
    save_figure(fig, FIG_DIR / "G3_fig07_af_roc.png")


def fig08_high_af_characterization(df_pass: pd.DataFrame):
    """Characterize variants with AF > 0.8."""
    _print("\n[G3] Fig08: High-AF characterization")

    pass_data = df_pass[df_pass["filter_is_pass"] == 1].copy()

    # AF bins
    bins = [(0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)]
    bin_labels = ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]

    rows = []
    for lo, hi in bins:
        sub = pass_data[(pass_data["af"] >= lo) & (pass_data["af"] < hi)]
        n = len(sub)
        n_tp = (sub["truth_label"] == "TP").sum()
        n_fp = (sub["truth_label"] == "FP").sum()
        rows.append({
            "af_bin": f"{lo:.1f}-{hi:.1f}",
            "n_total": n,
            "n_tp": n_tp,
            "n_fp": n_fp,
            "fp_rate": n_fp / n * 100 if n > 0 else 0,
            "tp_rate": n_tp / n * 100 if n > 0 else 0,
        })

    df_bins = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Stacked bar
    ax = axes[0]
    x = np.arange(len(bin_labels))
    ax.bar(x, df_bins["n_tp"], color=COLOR_TP, label="TP")
    ax.bar(x, df_bins["n_fp"], bottom=df_bins["n_tp"], color=COLOR_FP, label="FP")
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels)
    ax.set_xlabel("AF Bin")
    ax.set_ylabel("Count")
    ax.set_title("Variant Count by AF Bin")
    ax.legend()

    # Right: FP rate
    ax = axes[1]
    ax.bar(x, df_bins["fp_rate"], color=COLOR_FP, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels)
    ax.set_xlabel("AF Bin")
    ax.set_ylabel("FP Rate (%)")
    ax.set_title("FP Rate by AF Bin")
    ax.axhline(50, color="gray", linestyle="--", alpha=0.5, label="50%")
    ax.legend()

    fig.suptitle("G3-Fig08: High-AF Variant Characterization",
                 fontsize=14, fontweight="bold")
    add_caption(fig, "PASS variants, all 7 samples pooled")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G3_fig08_high_af_characterization.png")

    df_bins.to_csv(DATA_DIR / "af_bin_summary.tsv", sep="\t", index=False)
    _print(f"  [G3] AF bin summary:\n{df_bins.to_string(index=False)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    setup_plot_style()
    _print("=" * 70)
    _print("G3: Allele Fraction Architecture")
    _print("=" * 70)

    df_pass = load_pass_features()
    df_lite = load_all_lite()

    fig01_af_density_overlay(df_pass, df_lite)
    fig02_af_gmm_fit(df_pass)
    fig03_af_gt_joint(df_pass)
    fig04_fp_rate_by_af_bin(df_pass)
    fig05_af_threshold_sweep(df_pass)
    fig06_per_sample_af(df_pass)
    fig07_af_roc(df_pass)
    fig08_high_af_characterization(df_pass)

    _print(f"\n[G3] DONE. Output: {G3_DIR}")


if __name__ == "__main__":
    main()
