#!/usr/bin/env python3
"""O11: Methylation Heterogeneity — Germline vs Somatic Discrimination.

Tests whether within-group (same haplotype) methylation heterogeneity
differs between germline variants (FP) and somatic variants (TP).

Hypothesis:
  - Germline ASM (mQTL-driven): stable, concordant within each haplotype
  - Cancer ASM (stochastic): variable, discordant within each haplotype
  => Within-group heterogeneity should be LOWER for germline (FP) than somatic (TP)

Data source:
  Phase 1 shard manifest (637 paired regions) -> raw methylation matrices,
  reads metadata, and BERNOULLI distance matrices per region.

Outputs 8 figures + region_heterogeneity_features.tsv + observation_report.md
"""

from __future__ import annotations

import json
import sys
import traceback
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

# -- shared observation utilities --
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")
from observation_common import (
    setup_plot_style, save_figure, add_caption,
    compute_auc, compute_effect_size, mannwhitney_test,
    bootstrap_ci, format_p, format_ci,
    ensure_dir, safe_div, write_round_context, write_feature_statistics,
    OUTPUT_ROOT,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO,
    SAMPLE_ORDER,
)

warnings.filterwarnings("ignore", category=FutureWarning)

# ===== Configuration =====
OBS_ID = "O11"
RUN_DATE = datetime.now().strftime("%Y%m%d")
WORKSPACE_NAME = f"{RUN_DATE}_O11_methylation_heterogeneity"
OUT_DIR = OUTPUT_ROOT / WORKSPACE_NAME

MANIFEST_PATH = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_manifest_shard_export_paired_multibio_sample637/"
    "phase1_shard_manifest.tsv"
)

# Minimum reads per HP group to compute within-group features
MIN_GROUP_READS = 5
# Methylation binarization threshold
METHYL_THRESHOLD = 0.5
# Minimum CpG sites with data to compute epipolymorphism
MIN_CPGS_FOR_EPIPOLY = 5

# New feature names
HETERO_FEATURES = [
    "within_dist_var_mean",    # A1: mean of within-group pairwise distance variance
    "within_dist_mean",        # A1b: mean of within-group mean distance
    "per_cpg_concordance",     # A2: mean per-CpG concordance within groups
    "per_cpg_entropy",         # A3: mean per-CpG binary entropy within groups
    "epipolymorphism",         # A4: mean epipolymorphism score within groups
    "between_within_ratio",    # B1: within_mean / (within_mean + between_mean)
]


# ===== Feature Computation =====

def _upper_tri_values(matrix: np.ndarray) -> np.ndarray:
    """Extract upper triangle values (excluding diagonal) from square matrix."""
    n = matrix.shape[0]
    if n < 2:
        return np.array([])
    idx = np.triu_indices(n, k=1)
    return matrix[idx]


def _binary_entropy(p: float) -> float:
    """Binary entropy H(p) = -p*log2(p) - (1-p)*log2(1-p)."""
    if p <= 0 or p >= 1:
        return 0.0
    return -p * np.log2(p) - (1.0 - p) * np.log2(1.0 - p)


def compute_region_features(
    methyl_path: Path,
    reads_path: Path,
    dist_path: Path,
) -> Optional[Dict[str, Any]]:
    """Compute within-group heterogeneity features for one region.

    Returns dict of feature values, or None if data insufficient.
    """
    # Load methylation matrix (reads x CpGs, raw probabilities)
    try:
        methyl_df = pd.read_csv(methyl_path)
    except Exception:
        return None

    if methyl_df.empty or len(methyl_df) < 2 * MIN_GROUP_READS:
        return None

    # Load reads metadata
    try:
        reads_df = pd.read_csv(reads_path, sep="\t")
    except Exception:
        return None

    if "hp" not in reads_df.columns or "read_id" not in reads_df.columns:
        return None

    # Load BERNOULLI distance matrix
    try:
        dist_df = pd.read_csv(dist_path, index_col=0)
        dist_matrix = dist_df.values.astype(float)
    except Exception:
        dist_matrix = None

    # Map HP tags to groups (HP1-family vs HP2-family)
    hp_col = reads_df["hp"].astype(str)
    hp1_mask = hp_col.isin(["1", "1-1"])
    hp2_mask = hp_col.isin(["2", "2-1"])

    n_hp1 = int(hp1_mask.sum())
    n_hp2 = int(hp2_mask.sum())

    if n_hp1 < MIN_GROUP_READS or n_hp2 < MIN_GROUP_READS:
        return None

    # Get read indices
    hp1_idx = np.where(hp1_mask.values)[0]
    hp2_idx = np.where(hp2_mask.values)[0]

    # CpG columns (skip read_id column)
    cpg_cols = [c for c in methyl_df.columns if c != "read_id"]
    num_cpgs = len(cpg_cols)
    if num_cpgs < MIN_CPGS_FOR_EPIPOLY:
        return None

    # Raw methylation values for each group
    methyl_vals = methyl_df[cpg_cols].values  # (n_reads, n_cpgs)

    result = {
        "n_reads": len(methyl_df),
        "n_hp1": n_hp1,
        "n_hp2": n_hp2,
        "num_cpgs": num_cpgs,
    }

    # --- A1: Within-Group Pairwise Distance Variance ---
    if dist_matrix is not None and dist_matrix.shape[0] == len(methyl_df):
        within_vars = []
        within_means = []
        between_vals_list = []
        for idx in [hp1_idx, hp2_idx]:
            if len(idx) < 2:
                continue
            sub = dist_matrix[np.ix_(idx, idx)]
            tri = _upper_tri_values(sub)
            tri = tri[np.isfinite(tri)]
            if len(tri) > 0:
                within_vars.append(float(np.var(tri)))
                within_means.append(float(np.mean(tri)))

        if len(within_vars) == 2:
            result["within_dist_var_mean"] = float(np.mean(within_vars))
            result["within_dist_mean"] = float(np.mean(within_means))

            # Between-group distances
            between_sub = dist_matrix[np.ix_(hp1_idx, hp2_idx)]
            bvals = between_sub.ravel()
            bvals = bvals[np.isfinite(bvals)]
            if len(bvals) > 0:
                between_mean = float(np.mean(bvals))
                w_mean = result["within_dist_mean"]
                denom = w_mean + between_mean
                result["between_within_ratio"] = w_mean / denom if denom > 0 else np.nan
            else:
                result["between_within_ratio"] = np.nan
        else:
            result["within_dist_var_mean"] = np.nan
            result["within_dist_mean"] = np.nan
            result["between_within_ratio"] = np.nan
    else:
        result["within_dist_var_mean"] = np.nan
        result["within_dist_mean"] = np.nan
        result["between_within_ratio"] = np.nan

    # --- A2: Per-CpG Concordance ---
    concordances = []
    entropies = []
    for grp_idx in [hp1_idx, hp2_idx]:
        grp_methyl = methyl_vals[grp_idx]  # (n_group, n_cpgs)
        for j in range(num_cpgs):
            col_vals = grp_methyl[:, j]
            # Filter out NA
            valid = col_vals[~np.isnan(col_vals)] if hasattr(col_vals, '__len__') else np.array([])
            # Handle string "NA"
            try:
                valid = valid.astype(float)
                valid = valid[np.isfinite(valid)]
            except (ValueError, TypeError):
                continue
            if len(valid) < MIN_GROUP_READS:
                continue
            # Binarize
            frac_meth = float(np.mean(valid > METHYL_THRESHOLD))
            concordance = max(frac_meth, 1.0 - frac_meth)
            concordances.append(concordance)
            entropies.append(_binary_entropy(frac_meth))

    result["per_cpg_concordance"] = float(np.mean(concordances)) if concordances else np.nan
    result["per_cpg_entropy"] = float(np.mean(entropies)) if entropies else np.nan

    # --- A4: Epipolymorphism ---
    epipoly_scores = []
    for grp_idx in [hp1_idx, hp2_idx]:
        grp_methyl = methyl_vals[grp_idx]  # (n_group, n_cpgs)
        # Binarize the matrix, treat NA as separate category
        binary_rows = []
        for i in range(len(grp_idx)):
            row = grp_methyl[i]
            # Convert to binary string pattern
            pattern = []
            for v in row:
                try:
                    v_f = float(v)
                    if np.isnan(v_f):
                        pattern.append("N")
                    elif v_f > METHYL_THRESHOLD:
                        pattern.append("1")
                    else:
                        pattern.append("0")
                except (ValueError, TypeError):
                    pattern.append("N")
            binary_rows.append("".join(pattern))

        if len(binary_rows) < MIN_GROUP_READS:
            continue

        # Count pattern frequencies
        from collections import Counter
        counts = Counter(binary_rows)
        n_total = len(binary_rows)
        freqs = np.array([c / n_total for c in counts.values()])
        epipoly = 1.0 - float(np.sum(freqs ** 2))
        epipoly_scores.append(epipoly)

    result["epipolymorphism"] = float(np.mean(epipoly_scores)) if epipoly_scores else np.nan

    return result


# ===== Data Loading =====

def load_manifest() -> pd.DataFrame:
    """Load the Phase 1 shard manifest."""
    print(f"[{OBS_ID}] Loading manifest from {MANIFEST_PATH} ...")
    df = pd.read_csv(MANIFEST_PATH, sep="\t")
    print(f"[{OBS_ID}]   -> {len(df)} regions")
    print(f"[{OBS_ID}]   TP: {(df['truth_status']=='TP').sum()}, "
          f"FP: {(df['truth_status']=='FP').sum()}")
    return df


def extract_sample_name(dataset_label: str) -> str:
    """Extract sample name from dataset_label like 'HCC1395 5kHz paired'."""
    parts = dataset_label.split()
    if len(parts) >= 1:
        name = parts[0]
        if len(parts) >= 2 and parts[1].upper() == "DORADO":
            name = f"{name}_DORADO"
        return name
    return dataset_label


def compute_all_region_features(manifest: pd.DataFrame) -> pd.DataFrame:
    """Compute heterogeneity features for all regions in the manifest."""
    rows = []
    n_total = len(manifest)
    n_success = 0
    n_skip = 0

    for i, mrow in manifest.iterrows():
        region_dir = Path(mrow["region_dir"])
        region_key = mrow["region_key"]
        truth = mrow["truth_status"]
        sample = extract_sample_name(mrow["dataset_label"])

        methyl_path = region_dir / "methylation" / "methylation.csv"
        reads_path = region_dir / "reads" / "reads.tsv"
        dist_path = region_dir / "distance" / "BERNOULLI" / "matrix.csv"

        if not methyl_path.exists() or not reads_path.exists():
            n_skip += 1
            continue

        try:
            features = compute_region_features(methyl_path, reads_path, dist_path)
        except Exception as e:
            n_skip += 1
            if n_skip <= 5:
                print(f"[{OBS_ID}]   WARN: {region_key}: {e}")
            continue

        if features is None:
            n_skip += 1
            continue

        features["region_key"] = region_key
        features["truth_status"] = truth
        features["sample"] = sample
        features["quality_score"] = mrow.get("Quality_Score", np.nan)
        rows.append(features)
        n_success += 1

        if (n_success + n_skip) % 50 == 0:
            print(f"[{OBS_ID}]   Progress: {n_success + n_skip}/{n_total} "
                  f"(success={n_success}, skip={n_skip})")

    print(f"[{OBS_ID}] Feature extraction complete: {n_success} regions "
          f"({n_skip} skipped)")
    return pd.DataFrame(rows)


# ===== Figures =====

def fig01_within_dist_variance(df: pd.DataFrame) -> Path:
    """Violin: within-group distance variance by truth label."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [3, 1]})

    ax = axes[0]
    plot_df = df.dropna(subset=["within_dist_var_mean"])
    sns.violinplot(
        data=plot_df, x="sample", y="within_dist_var_mean",
        hue="truth_status", split=True,
        palette={"TP": COLOR_TP, "FP": COLOR_FP},
        inner="quartile", order=SAMPLE_ORDER, ax=ax, cut=0,
    )
    ax.set_title("Fig 01: Within-Group Distance Variance by Truth Label")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Mean within-group pairwise distance variance")
    ax.tick_params(axis="x", rotation=30)

    # Stats panel
    ax2 = axes[1]
    ax2.axis("off")
    tp_v = plot_df.loc[plot_df["truth_status"] == "TP", "within_dist_var_mean"].values
    fp_v = plot_df.loc[plot_df["truth_status"] == "FP", "within_dist_var_mean"].values
    mw = mannwhitney_test(tp_v, fp_v)
    es = compute_effect_size(tp_v, fp_v)
    auc = compute_auc(
        (plot_df["truth_status"] == "TP").astype(int).values,
        plot_df["within_dist_var_mean"].values,
    )
    stats_text = (
        f"TP mean: {np.nanmean(tp_v):.4f} (n={len(tp_v)})\n"
        f"FP mean: {np.nanmean(fp_v):.4f} (n={len(fp_v)})\n"
        f"Delta: {np.nanmean(tp_v) - np.nanmean(fp_v):.4f}\n"
        f"{'─' * 25}\n"
        f"Mann-Whitney p: {format_p(mw['p_value'])}\n"
        f"Cohen's d: {es['d']:.4f}\n"
        f"  95% CI: {format_ci(es['ci_lo'], es['ci_hi'])}\n"
        f"AUC(TP): {auc:.4f}\n"
        f"AUC(best): {max(auc, 1-auc):.4f}"
    )
    ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes,
             va="top", fontsize=9, fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#f5f5f5"))

    add_caption(fig, f"Source: Phase 1 manifest (paired) | {len(plot_df)} regions | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig01_within_dist_variance.png")


def fig02_per_cpg_concordance(df: pd.DataFrame) -> Path:
    """Violin: per-CpG concordance by truth label, per sample."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [3, 1]})

    ax = axes[0]
    plot_df = df.dropna(subset=["per_cpg_concordance"])
    sns.violinplot(
        data=plot_df, x="sample", y="per_cpg_concordance",
        hue="truth_status", split=True,
        palette={"TP": COLOR_TP, "FP": COLOR_FP},
        inner="quartile", order=SAMPLE_ORDER, ax=ax, cut=0,
    )
    ax.set_title("Fig 02: Per-CpG Concordance by Truth Label")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Mean per-CpG within-group concordance")
    ax.tick_params(axis="x", rotation=30)

    ax2 = axes[1]
    ax2.axis("off")
    tp_v = plot_df.loc[plot_df["truth_status"] == "TP", "per_cpg_concordance"].values
    fp_v = plot_df.loc[plot_df["truth_status"] == "FP", "per_cpg_concordance"].values
    mw = mannwhitney_test(tp_v, fp_v)
    es = compute_effect_size(tp_v, fp_v)
    auc = compute_auc(
        (plot_df["truth_status"] == "TP").astype(int).values,
        plot_df["per_cpg_concordance"].values,
    )
    stats_text = (
        f"TP mean: {np.nanmean(tp_v):.4f} (n={len(tp_v)})\n"
        f"FP mean: {np.nanmean(fp_v):.4f} (n={len(fp_v)})\n"
        f"{'─' * 25}\n"
        f"Mann-Whitney p: {format_p(mw['p_value'])}\n"
        f"Cohen's d: {es['d']:.4f}\n"
        f"  95% CI: {format_ci(es['ci_lo'], es['ci_hi'])}\n"
        f"AUC(TP): {auc:.4f}\n"
        f"AUC(best): {max(auc, 1-auc):.4f}"
    )
    ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes,
             va="top", fontsize=9, fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#f5f5f5"))

    add_caption(fig, f"Source: Phase 1 manifest (paired) | {len(plot_df)} regions | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig02_per_cpg_concordance.png")


def fig03_epipolymorphism(df: pd.DataFrame) -> Path:
    """Box plot: epipolymorphism score by truth label with effect sizes."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [3, 1]})

    ax = axes[0]
    plot_df = df.dropna(subset=["epipolymorphism"])
    sns.boxplot(
        data=plot_df, x="sample", y="epipolymorphism",
        hue="truth_status",
        palette={"TP": COLOR_TP, "FP": COLOR_FP},
        order=SAMPLE_ORDER, ax=ax, fliersize=1,
    )
    ax.set_title("Fig 03: Epipolymorphism Score by Truth Label")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Epipolymorphism (1 - Σp²)")
    ax.tick_params(axis="x", rotation=30)

    # Per-sample effect sizes
    ax2 = axes[1]
    ax2.axis("off")
    lines = []
    for s in SAMPLE_ORDER:
        sub = plot_df[plot_df["sample"] == s]
        tp_v = sub.loc[sub["truth_status"] == "TP", "epipolymorphism"].values
        fp_v = sub.loc[sub["truth_status"] == "FP", "epipolymorphism"].values
        es = compute_effect_size(tp_v, fp_v)
        a = compute_auc(
            (sub["truth_status"] == "TP").astype(int).values,
            sub["epipolymorphism"].values,
        )
        lines.append(f"{s[:8]:8s} d={es['d']:+.3f} AUC={max(a, 1-a):.3f}")
    # Pooled
    tp_all = plot_df.loc[plot_df["truth_status"] == "TP", "epipolymorphism"].values
    fp_all = plot_df.loc[plot_df["truth_status"] == "FP", "epipolymorphism"].values
    es_all = compute_effect_size(tp_all, fp_all)
    auc_all = compute_auc(
        (plot_df["truth_status"] == "TP").astype(int).values,
        plot_df["epipolymorphism"].values,
    )
    lines.append(f"{'─' * 30}")
    lines.append(f"{'POOLED':8s} d={es_all['d']:+.3f} AUC={max(auc_all, 1-auc_all):.3f}")
    ax2.text(0.05, 0.95, "\n".join(lines), transform=ax2.transAxes,
             va="top", fontsize=8, fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#f5f5f5"))

    add_caption(fig, f"Source: Phase 1 manifest (paired) | {len(plot_df)} regions | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig03_epipolymorphism.png")


def fig04_feature_auc_ranking(df: pd.DataFrame) -> Path:
    """Horizontal bar chart: AUC ranking for all heterogeneity features."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    results = []
    is_tp = (df["truth_status"] == "TP").astype(int).values
    for feat in HETERO_FEATURES:
        if feat not in df.columns:
            continue
        vals = df[feat].values
        auc_raw = compute_auc(is_tp, vals)
        auc_best = max(auc_raw, 1 - auc_raw) if not np.isnan(auc_raw) else np.nan
        es = compute_effect_size(
            df.loc[df["truth_status"] == "TP", feat].dropna().values,
            df.loc[df["truth_status"] == "FP", feat].dropna().values,
        )
        results.append({
            "feature": feat,
            "auc_raw": auc_raw,
            "auc_best": auc_best,
            "cohens_d": es["d"],
            "d_ci_lo": es["ci_lo"],
            "d_ci_hi": es["ci_hi"],
        })

    rdf = pd.DataFrame(results).sort_values("auc_best", ascending=True)
    colors = [COLOR_TP if a >= 0.55 else COLOR_FP if a < 0.52 else "#9E9E9E"
              for a in rdf["auc_best"]]
    ax.barh(rdf["feature"], rdf["auc_best"], color=colors, alpha=0.8)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5, label="Random")
    ax.axvline(0.55, color="orange", linestyle=":", alpha=0.5, label="Threshold (0.55)")
    ax.set_xlim(0.45, max(0.65, rdf["auc_best"].max() + 0.03))
    ax.set_title("Fig 04: Heterogeneity Feature AUC Ranking (Paired TP/FP)")
    ax.set_xlabel("AUC (best direction)")
    ax.legend(fontsize=8)

    for idx_i, (_, row) in enumerate(rdf.iterrows()):
        ax.text(row["auc_best"] + 0.003, idx_i,
                f"{row['auc_best']:.3f} (d={row['cohens_d']:+.3f})",
                va="center", fontsize=8)

    add_caption(fig, f"Source: Phase 1 manifest (paired) | {len(df)} regions | {RUN_DATE}")
    fig.tight_layout()

    # Save table
    rdf.sort_values("auc_best", ascending=False).to_csv(
        OUT_DIR / "fig04_auc_ranking.tsv", sep="\t", index=False)
    return save_figure(fig, OUT_DIR / "fig04_feature_auc_ranking.png")


def fig05_cpg_density_stratified(df: pd.DataFrame) -> Path:
    """3-panel AUC by CpG density bin."""
    bins = [("Low (<20)", df["num_cpgs"] < 20),
            ("Medium (20-50)", (df["num_cpgs"] >= 20) & (df["num_cpgs"] < 50)),
            ("High (≥50)", df["num_cpgs"] >= 50)]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for ax, (bin_name, mask) in zip(axes, bins):
        sub = df[mask].dropna(subset=HETERO_FEATURES[:1])  # at least one feature
        n_tp = (sub["truth_status"] == "TP").sum()
        n_fp = (sub["truth_status"] == "FP").sum()

        if n_tp < 10 or n_fp < 10:
            ax.set_title(f"{bin_name}\n(n_tp={n_tp}, n_fp={n_fp})\nInsufficient data")
            ax.set_xlim(0.4, 0.7)
            continue

        is_tp = (sub["truth_status"] == "TP").astype(int).values
        results = []
        for feat in HETERO_FEATURES:
            if feat not in sub.columns:
                continue
            auc_raw = compute_auc(is_tp, sub[feat].values)
            auc_best = max(auc_raw, 1 - auc_raw) if not np.isnan(auc_raw) else np.nan
            results.append({"feature": feat, "auc_best": auc_best})

        rdf = pd.DataFrame(results).sort_values("auc_best", ascending=True)
        colors = [COLOR_TP if a >= 0.55 else "#9E9E9E" for a in rdf["auc_best"]]
        ax.barh(rdf["feature"], rdf["auc_best"], color=colors, alpha=0.8)
        ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5)
        ax.set_xlim(0.4, max(0.7, rdf["auc_best"].max() + 0.03))
        ax.set_title(f"{bin_name}\n(n={len(sub)}, TP={n_tp}, FP={n_fp})")
        ax.set_xlabel("AUC (best)")

        for idx_i, (_, row) in enumerate(rdf.iterrows()):
            ax.text(row["auc_best"] + 0.003, idx_i, f"{row['auc_best']:.3f}",
                    va="center", fontsize=7)

    fig.suptitle("Fig 05: CpG Density Stratified AUC", fontsize=14, y=1.02)
    add_caption(fig, f"Source: Phase 1 manifest (paired) | {len(df)} regions | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig05_cpg_density_stratified.png")


def fig06_correlation_with_existing(df: pd.DataFrame) -> Path:
    """Heatmap: correlation between new and existing ISM features."""
    existing_features = ["quality_score"]
    all_feats = HETERO_FEATURES + existing_features
    available = [f for f in all_feats if f in df.columns and df[f].notna().sum() > 10]

    if len(available) < 2:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "Insufficient features for correlation",
                ha="center", va="center", fontsize=12)
        ax.set_title("Fig 06: Feature Correlation Matrix")
        fig.tight_layout()
        return save_figure(fig, OUT_DIR / "fig06_feature_correlation.png")

    corr = df[available].corr(method="spearman")

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(
        corr, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
        vmin=-1, vmax=1, ax=ax, square=True,
        linewidths=0.5, cbar_kws={"shrink": 0.8},
    )
    ax.set_title("Fig 06: Spearman Correlation — New vs Existing Features")

    add_caption(fig, f"Source: Phase 1 manifest (paired) | {len(df)} regions | {RUN_DATE}")
    fig.tight_layout()

    corr.to_csv(OUT_DIR / "fig06_correlation_matrix.tsv", sep="\t")
    return save_figure(fig, OUT_DIR / "fig06_feature_correlation.png")


def fig07_cross_feature_scatter(df: pd.DataFrame) -> Path:
    """Scatter: within_dist_var_mean vs per_cpg_entropy, colored by truth."""
    feat_x = "within_dist_var_mean"
    feat_y = "per_cpg_entropy"

    if feat_x not in df.columns or feat_y not in df.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "Features not available", ha="center", va="center")
        return save_figure(fig, OUT_DIR / "fig07_cross_feature_scatter.png")

    plot_df = df.dropna(subset=[feat_x, feat_y])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Scatter by truth
    ax = axes[0]
    for label, color in [("FP", COLOR_FP), ("TP", COLOR_TP)]:
        sub = plot_df[plot_df["truth_status"] == label]
        ax.scatter(sub[feat_x], sub[feat_y], c=color, alpha=0.5,
                   s=20, label=f"{label} (n={len(sub)})", edgecolors="none")
    ax.set_xlabel("Within-Group Distance Variance")
    ax.set_ylabel("Per-CpG Entropy")
    ax.set_title("Heterogeneity Features by Truth Label")
    ax.legend()

    # Scatter by sample
    ax = axes[1]
    sample_colors = plt.cm.tab10(np.linspace(0, 1, len(SAMPLE_ORDER)))
    for s, c in zip(SAMPLE_ORDER, sample_colors):
        sub = plot_df[plot_df["sample"] == s]
        ax.scatter(sub[feat_x], sub[feat_y], c=[c], alpha=0.5,
                   s=20, label=f"{s} ({len(sub)})", edgecolors="none")
    ax.set_xlabel("Within-Group Distance Variance")
    ax.set_ylabel("Per-CpG Entropy")
    ax.set_title("Heterogeneity Features by Sample")
    ax.legend(fontsize=7)

    fig.suptitle("Fig 07: 2D Feature Space Exploration", fontsize=14, y=1.02)
    add_caption(fig, f"Source: Phase 1 manifest (paired) | {len(plot_df)} regions | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig07_cross_feature_scatter.png")


def fig08_combined_roc(df: pd.DataFrame) -> Path:
    """ROC curve for individual and combined features."""
    from sklearn.metrics import roc_curve, roc_auc_score
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    fig, ax = plt.subplots(figsize=(8, 8))

    is_tp = (df["truth_status"] == "TP").astype(int).values
    colors = plt.cm.Set1(np.linspace(0, 1, len(HETERO_FEATURES) + 1))

    # Individual features
    for i, feat in enumerate(HETERO_FEATURES):
        if feat not in df.columns:
            continue
        vals = df[feat].values
        mask = np.isfinite(vals) & np.isfinite(is_tp.astype(float))
        if mask.sum() < 20:
            continue
        y = is_tp[mask]
        x = vals[mask]
        if len(np.unique(y)) < 2:
            continue
        auc = roc_auc_score(y, x)
        # Use best direction
        if auc < 0.5:
            x = -x
            auc = 1 - auc
        fpr, tpr, _ = roc_curve(y, x)
        ax.plot(fpr, tpr, color=colors[i], alpha=0.7,
                label=f"{feat} (AUC={auc:.3f})")

    # Combined logistic regression
    avail = [f for f in HETERO_FEATURES if f in df.columns]
    combo_df = df[avail + ["truth_status"]].dropna()
    if len(combo_df) >= 30 and len(combo_df["truth_status"].unique()) == 2:
        X = combo_df[avail].values
        y = (combo_df["truth_status"] == "TP").astype(int).values
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        try:
            lr = LogisticRegression(max_iter=1000, random_state=42)
            lr.fit(X_scaled, y)
            y_prob = lr.predict_proba(X_scaled)[:, 1]
            auc_combo = roc_auc_score(y, y_prob)
            fpr, tpr, _ = roc_curve(y, y_prob)
            ax.plot(fpr, tpr, color="black", linewidth=2,
                    label=f"Combined LR (AUC={auc_combo:.3f})")
        except Exception as e:
            print(f"[{OBS_ID}] Combined LR failed: {e}")

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="Random (0.500)")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Fig 08: ROC — Individual & Combined Heterogeneity Features")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)

    add_caption(fig, f"Source: Phase 1 manifest (paired) | {len(df)} regions | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig08_combined_roc.png")


# ===== Report =====

def write_observation_report(df: pd.DataFrame, auc_table: pd.DataFrame) -> Path:
    """Write a concise observation report in markdown."""
    path = OUT_DIR / "observation_report.md"

    n_tp = (df["truth_status"] == "TP").sum()
    n_fp = (df["truth_status"] == "FP").sum()

    lines = [
        f"# O11: Methylation Heterogeneity — Observation Report",
        f"",
        f"**Date**: {RUN_DATE}",
        f"**Regions**: {len(df)} (TP={n_tp}, FP={n_fp})",
        f"**Hypothesis**: Germline ASM → lower within-group heterogeneity than somatic ASM",
        f"",
        f"## Feature AUC Summary (Paired, TP/FP)",
        f"",
    ]

    if auc_table is not None and len(auc_table) > 0:
        lines.append("| Feature | AUC (best) | Cohen's d | Direction |")
        lines.append("|---------|-----------|-----------|-----------|")
        for _, row in auc_table.iterrows():
            direction = "TP>FP" if row["auc_raw"] >= 0.5 else "FP>TP"
            lines.append(
                f"| {row['feature']} | {row['auc_best']:.3f} | "
                f"{row['cohens_d']:+.3f} | {direction} |"
            )
    lines.append("")

    # Go/No-Go assessment
    if auc_table is not None and len(auc_table) > 0:
        best_auc = auc_table["auc_best"].max()
        if best_auc >= 0.65:
            verdict = "**Outcome A: Strong signal** — proceed to pipeline integration"
        elif best_auc >= 0.55:
            verdict = "**Outcome B: Moderate signal** — consider as ML ensemble auxiliary feature"
        else:
            verdict = "**Outcome C: No signal** — within-group heterogeneity does not distinguish germline/somatic"
        lines.append(f"## Go/No-Go Decision")
        lines.append(f"")
        lines.append(f"Best AUC: {best_auc:.3f}")
        lines.append(f"")
        lines.append(verdict)
        lines.append("")

    # Consistency checks
    lines.append("## CONSISTENCY Checks")
    lines.append("")

    # Check: FP AlleleDelta > TP baseline (via QS as proxy)
    tp_qs = df.loc[df["truth_status"] == "TP", "quality_score"].dropna()
    fp_qs = df.loc[df["truth_status"] == "FP", "quality_score"].dropna()
    if len(tp_qs) > 0 and len(fp_qs) > 0:
        lines.append(f"- QS baseline: TP median={tp_qs.median():.1f}, "
                      f"FP median={fp_qs.median():.1f} "
                      f"{'✓ TP>FP as expected' if tp_qs.median() > fp_qs.median() else '⚠ unexpected'}")

    # Check: feature direction consistency
    if auc_table is not None:
        n_flip = (auc_table["auc_raw"] < 0.5).sum()
        lines.append(f"- Feature direction: {n_flip}/{len(auc_table)} features have FP>TP direction")

    lines.append("")
    lines.append(f"---")
    lines.append(f"*Generated by build_observation_O11_methylation_heterogeneity.py*")

    with path.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[{OBS_ID}] Report saved: {path.name}")
    return path


# ===== Main =====

def main():
    print("=" * 70)
    print(f"O11: Methylation Heterogeneity — Germline vs Somatic Discrimination")
    print(f"Run date: {RUN_DATE}")
    print(f"Output: {OUT_DIR}")
    print("=" * 70)

    setup_plot_style()
    ensure_dir(OUT_DIR)

    # Step 1: Load manifest and compute features
    manifest = load_manifest()
    print(f"\n[{OBS_ID}] Step 1: Computing region-level heterogeneity features...")
    df = compute_all_region_features(manifest)

    if len(df) < 20:
        print(f"[{OBS_ID}] ERROR: Only {len(df)} regions with valid features. Aborting.")
        return

    # Save features
    feat_path = OUT_DIR / "region_heterogeneity_features.tsv"
    df.to_csv(feat_path, sep="\t", index=False)
    print(f"[{OBS_ID}] Saved {len(df)} regions to {feat_path.name}")

    # Feature statistics
    write_feature_statistics(df, HETERO_FEATURES, OUT_DIR / "feature_statistics.tsv")

    # Round context
    write_round_context(
        out_dir=OUT_DIR,
        observation_id=OBS_ID,
        title="Methylation Heterogeneity — Germline vs Somatic Discrimination",
        description=(
            "Tests whether within-group (same haplotype) methylation heterogeneity "
            "differs between germline FP and somatic TP variants, using raw per-region "
            "methylation matrices from paired-mode data."
        ),
        script_path="scripts/analysis/build_observation_O11_methylation_heterogeneity.py",
        data_sources=[{
            "name": "Phase 1 shard manifest (paired multibio 637)",
            "path": str(MANIFEST_PATH),
            "rows": len(manifest),
        }],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "n_tp": int((df["truth_status"] == "TP").sum()),
            "n_fp": int((df["truth_status"] == "FP").sum()),
            "features_computed": HETERO_FEATURES,
            "min_group_reads": MIN_GROUP_READS,
            "methyl_threshold": METHYL_THRESHOLD,
        },
    )

    # Step 2: Generate figures
    print(f"\n[{OBS_ID}] Step 2: Generating figures...")
    fig01_within_dist_variance(df)
    fig02_per_cpg_concordance(df)
    fig03_epipolymorphism(df)
    auc_path = fig04_feature_auc_ranking(df)
    fig05_cpg_density_stratified(df)
    fig06_correlation_with_existing(df)
    fig07_cross_feature_scatter(df)
    fig08_combined_roc(df)

    # Load AUC table for report
    auc_table = pd.read_csv(OUT_DIR / "fig04_auc_ranking.tsv", sep="\t")

    # Step 3: Write report
    print(f"\n[{OBS_ID}] Step 3: Writing observation report...")
    write_observation_report(df, auc_table)

    # Summary
    print(f"\n[{OBS_ID}] === Summary ===")
    print(f"Regions processed: {len(df)}")
    print(f"  TP: {(df['truth_status']=='TP').sum()}, FP: {(df['truth_status']=='FP').sum()}")
    if auc_table is not None and len(auc_table) > 0:
        best = auc_table.loc[auc_table["auc_best"].idxmax()]
        print(f"Best feature: {best['feature']} (AUC={best['auc_best']:.3f}, d={best['cohens_d']:+.3f})")
        best_auc = best["auc_best"]
        if best_auc >= 0.65:
            print("→ Outcome A: STRONG SIGNAL — proceed to integration")
        elif best_auc >= 0.55:
            print("→ Outcome B: MODERATE SIGNAL — consider as auxiliary ML feature")
        else:
            print("→ Outcome C: NO SIGNAL — confirm Phase 2 normal ref is the only path")

    print(f"\n[{OBS_ID}] All outputs saved to: {OUT_DIR}")
    print(f"[{OBS_ID}] Done.")


if __name__ == "__main__":
    main()
