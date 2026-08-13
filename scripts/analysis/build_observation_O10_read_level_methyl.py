#!/usr/bin/env python3
"""O10: Read-Level Methylation Feature Observation.

Analyses read-level methylation features from Phase 1A read training tables
to understand ALT/REF read methylation distributions, HP-tag stratification,
truth-label differences, and feature importance for read classification.

Data source:
  Phase 1A paired-mode read training table (86,521 reads across 7 samples).

Outputs 8 figures + feature_statistics.tsv + round_context.json
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# -- shared observation utilities --
from scripts.analysis.observation_common import (  # noqa: E402
    setup_plot_style, save_figure, add_caption,
    compute_auc, compute_effect_size, mannwhitney_test, chi2_test,
    bootstrap_ci, format_p, format_ci,
    ensure_dir, safe_div, write_round_context, write_feature_statistics,
    markdown_table, OUTPUT_ROOT,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO,
)

from scripts.lib.verification_schema_contract import (  # noqa: E402
    LEGACY_CLASSES,
    SchemaContractError,
    select_legacy_view,
)

warnings.filterwarnings("ignore", category=FutureWarning)

# ===== Configuration =====
OBS_ID = "O10"
RUN_DATE = "20260401"
WORKSPACE_NAME = f"{RUN_DATE}_O10_read_level_methyl"
OUT_DIR = OUTPUT_ROOT / WORKSPACE_NAME

READ_TABLE_PATH = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_manifest_shard_export_paired_multibio_sample637/"
    "phase1_shard_read_training_table.tsv"
)

# Colors
COLOR_ALT = "#4CAF50"
COLOR_REF = "#9C27B0"
ALT_REF_PALETTE = {"ALT": COLOR_ALT, "REF": COLOR_REF}

METHYL_FEATURES = [
    "methyl_mean", "methyl_std", "methyl_median",
    "methyl_high_fraction", "methyl_low_fraction",
    "methyl_na_fraction", "num_cpg_observed",
]

SAMPLE_ORDER = [
    "HCC1395", "HCC1395_DORADO", "COLO829",
    "H1437", "H2009", "HCC1937", "HCC1954",
]


# ===== Data Loading =====
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read-table", type=Path, default=READ_TABLE_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize an historical unversioned VerificationClass column; "
            "unknown values are excluded from the legacy four-state figure and counted."
        ),
    )
    return parser.parse_args()


def attach_legacy_verification_view(
    df: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> tuple[pd.DataFrame, Dict[str, object]]:
    """Require L4 provenance for the historical read-level panel."""
    view = select_legacy_view(
        df,
        allow_unversioned_v1=allow_unversioned_v1,
        unversioned_unknown_policy="exclude",
    )
    selected = df.copy()
    selected["VerificationClass_SourceValue"] = selected[view.field]
    selected["VerificationClass"] = pd.Categorical(
        view.values,
        categories=LEGACY_CLASSES,
        ordered=True,
    )
    selected["VerificationSourceField"] = view.field
    selected["VerificationSchemaStatus"] = view.schema_status
    metadata = view.metadata()
    metadata.update(
        {
            "requested_view": "legacy",
            "input_rows": int(len(selected)),
            "excluded_unknown_rows": int(selected["VerificationClass"].isna().sum()),
        }
    )
    return selected, metadata


def load_read_training_table(
    path: Path = READ_TABLE_PATH,
    allow_unversioned_v1: bool = False,
) -> tuple[pd.DataFrame, Dict[str, object]]:
    """Load the Phase 1A read training table."""
    print(f"[O10] Loading read training table from {path} ...")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    try:
        df, verification_metadata = attach_legacy_verification_view(
            df,
            allow_unversioned_v1=allow_unversioned_v1,
        )
    except SchemaContractError as exc:
        raise SystemExit(f"[O10][schema-contract] {path}: {exc}") from exc
    print(f"[O10]   -> {len(df):,} reads x {len(df.columns)} cols")
    print(f"[O10]   Samples: {df['sample'].nunique()}, Regions: {df['region_key'].nunique()}")
    print(f"[O10]   TP reads: {(df['truth_status']=='TP').sum():,}, FP reads: {(df['truth_status']=='FP').sum():,}")
    print(f"[O10]   ALT reads: {(df['alt_support']=='ALT').sum():,}, REF reads: {(df['alt_support']=='REF').sum():,}")
    print(f"[O10]   Verification contract: {verification_metadata}")
    return df, verification_metadata


# ===== Derived Columns =====
def add_derived_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add useful derived columns."""
    # Binary labels
    df["is_alt"] = (df["alt_support"] == "ALT").astype(int)
    df["is_tp"] = (df["truth_status"] == "TP").astype(int)

    # Simplified HP tag: 1, 2, 0 (unphased), other
    hp_map = {"1": "HP1", "2": "HP2", "0": "HP0_unphased",
              "1-1": "HP1(ALT)", "2-1": "HP2(ALT)", "3": "HP_other"}
    df["hp_simple"] = df["hp"].astype(str).map(hp_map).fillna("HP_other")

    # CpG observation rate
    df["cpg_obs_rate"] = df["num_cpg_observed"] / df["num_cpg_total"].clip(lower=1)

    # Methylation bimodality: high_fraction + low_fraction (close to 1 = bimodal)
    df["methyl_bimodal_score"] = df["methyl_high_fraction"] + df["methyl_low_fraction"]

    # Intermediate fraction = 1 - high - low (methylation in middle range)
    df["methyl_intermediate_fraction"] = 1.0 - df["methyl_bimodal_score"]

    return df


# ===== Figure 1: methyl_mean by ALT/REF violin =====
def fig01_methyl_mean_by_alt_support(df: pd.DataFrame) -> Path:
    """Violin plot of methyl_mean split by ALT vs REF support."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [3, 1]})

    # Main violin
    ax = axes[0]
    sns.violinplot(
        data=df, x="sample", y="methyl_mean", hue="alt_support",
        split=True, palette=ALT_REF_PALETTE, inner="quartile",
        order=SAMPLE_ORDER, ax=ax, cut=0,
    )
    ax.set_title("Fig 01: Read-Level Methylation Mean by ALT/REF Support")
    ax.set_xlabel("Sample")
    ax.set_ylabel("methyl_mean (per read)")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(title="Support", loc="upper right")

    # Summary statistics panel
    ax2 = axes[1]
    ax2.axis("off")
    alt_vals = df.loc[df["alt_support"] == "ALT", "methyl_mean"]
    ref_vals = df.loc[df["alt_support"] == "REF", "methyl_mean"]
    mw = mannwhitney_test(alt_vals.values, ref_vals.values)
    es = compute_effect_size(alt_vals.values, ref_vals.values)
    auc = compute_auc(df["is_alt"].values, df["methyl_mean"].values)

    stats_text = (
        f"ALT mean: {alt_vals.mean():.4f} (n={len(alt_vals):,})\n"
        f"REF mean: {ref_vals.mean():.4f} (n={len(ref_vals):,})\n"
        f"Delta: {alt_vals.mean() - ref_vals.mean():.4f}\n"
        f"───────────────────\n"
        f"Mann-Whitney p: {format_p(mw['p_value'])}\n"
        f"Effect r: {mw['effect_size_r']:.4f}\n"
        f"Cohen's d: {es['d']:.4f}\n"
        f"  95% CI: {format_ci(es['ci_lo'], es['ci_hi'])}\n"
        f"AUC(ALT): {auc:.4f}"
    )
    ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes,
             verticalalignment="top", fontsize=9, fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#f5f5f5"))

    add_caption(fig, f"Source: Phase 1A read training table | {len(df):,} reads | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig01_read_methyl_mean_by_alt_support.png")


# ===== Figure 2: methyl_na_fraction histogram =====
def fig02_methyl_na_fraction(df: pd.DataFrame) -> Path:
    """Histogram of NA fraction by ALT/REF and TP/FP."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # By ALT/REF
    ax = axes[0]
    for label, color in ALT_REF_PALETTE.items():
        vals = df.loc[df["alt_support"] == label, "methyl_na_fraction"]
        ax.hist(vals, bins=50, alpha=0.5, color=color, label=label, density=True)
    ax.set_title("NA Fraction by ALT/REF")
    ax.set_xlabel("methyl_na_fraction")
    ax.set_ylabel("Density")
    ax.legend()
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5, label="50%")

    # By TP/FP
    ax = axes[1]
    tp_palette = {"TP": COLOR_TP, "FP": COLOR_FP}
    for label, color in tp_palette.items():
        vals = df.loc[df["truth_status"] == label, "methyl_na_fraction"]
        ax.hist(vals, bins=50, alpha=0.5, color=color, label=label, density=True)
    ax.set_title("NA Fraction by TP/FP")
    ax.set_xlabel("methyl_na_fraction")
    ax.set_ylabel("Density")
    ax.legend()
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5)

    # Stats annotation
    alt_na = df.loc[df["alt_support"] == "ALT", "methyl_na_fraction"]
    ref_na = df.loc[df["alt_support"] == "REF", "methyl_na_fraction"]
    mw = mannwhitney_test(alt_na.values, ref_na.values)
    axes[0].text(0.98, 0.95, f"MW p={format_p(mw['p_value'])}\nr={mw['effect_size_r']:.3f}",
                 transform=axes[0].transAxes, ha="right", va="top", fontsize=8,
                 bbox=dict(facecolor="white", alpha=0.8))

    tp_na = df.loc[df["truth_status"] == "TP", "methyl_na_fraction"]
    fp_na = df.loc[df["truth_status"] == "FP", "methyl_na_fraction"]
    mw2 = mannwhitney_test(tp_na.values, fp_na.values)
    axes[1].text(0.98, 0.95, f"MW p={format_p(mw2['p_value'])}\nr={mw2['effect_size_r']:.3f}",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=8,
                 bbox=dict(facecolor="white", alpha=0.8))

    fig.suptitle("Fig 02: Methylation NA Fraction Distribution", fontsize=13, y=1.02)
    add_caption(fig, f"Source: Phase 1A read training table | {len(df):,} reads | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig02_read_methyl_na_fraction.png")


# ===== Figure 3: CpG count distribution =====
def fig03_cpg_count_distribution(df: pd.DataFrame) -> Path:
    """Distribution of CpG observed per read, by ALT/REF and TP/FP."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Overall distribution
    ax = axes[0]
    ax.hist(df["num_cpg_observed"], bins=80, color="#607D8B", alpha=0.7, edgecolor="white")
    ax.axvline(df["num_cpg_observed"].median(), color="red", linestyle="--",
               label=f"Median={df['num_cpg_observed'].median():.0f}")
    ax.set_title("Overall CpG Observed/Read")
    ax.set_xlabel("num_cpg_observed")
    ax.set_ylabel("Count")
    ax.legend(fontsize=8)

    # By ALT/REF
    ax = axes[1]
    for label, color in ALT_REF_PALETTE.items():
        vals = df.loc[df["alt_support"] == label, "num_cpg_observed"]
        ax.hist(vals, bins=60, alpha=0.5, color=color, label=label, density=True)
    ax.set_title("By ALT/REF Support")
    ax.set_xlabel("num_cpg_observed")
    ax.set_ylabel("Density")
    ax.legend()

    # CpG obs rate by sample
    ax = axes[2]
    sns.boxplot(data=df, x="sample", y="cpg_obs_rate", order=SAMPLE_ORDER,
                color="#B0BEC5", ax=ax, fliersize=1)
    ax.set_title("CpG Observation Rate by Sample")
    ax.set_xlabel("Sample")
    ax.set_ylabel("cpg_obs_rate (observed/total)")
    ax.tick_params(axis="x", rotation=30)

    # Stats
    alt_cpg = df.loc[df["alt_support"] == "ALT", "num_cpg_observed"]
    ref_cpg = df.loc[df["alt_support"] == "REF", "num_cpg_observed"]
    mw = mannwhitney_test(alt_cpg.values, ref_cpg.values)
    axes[1].text(0.98, 0.95, f"MW p={format_p(mw['p_value'])}\nr={mw['effect_size_r']:.3f}",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=8,
                 bbox=dict(facecolor="white", alpha=0.8))

    fig.suptitle("Fig 03: CpG Count Distribution per Read", fontsize=13, y=1.02)
    add_caption(fig, f"Source: Phase 1A read training table | {len(df):,} reads | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig03_read_cpg_count_distribution.png")


# ===== Figure 4: Methylation by HP tag =====
def fig04_methyl_by_hp_tag(df: pd.DataFrame) -> Path:
    """Methylation features stratified by HP tag."""
    hp_order = ["HP1", "HP2", "HP1(ALT)", "HP2(ALT)", "HP0_unphased", "HP_other"]
    hp_colors = {
        "HP1": "#1976D2", "HP2": "#D32F2F",
        "HP1(ALT)": "#64B5F6", "HP2(ALT)": "#EF9A9A",
        "HP0_unphased": "#9E9E9E", "HP_other": "#BDBDBD",
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    features = ["methyl_mean", "methyl_std", "methyl_high_fraction", "methyl_low_fraction"]
    titles = ["Methylation Mean", "Methylation Std", "High Fraction (>0.5)", "Low Fraction (<0.5)"]

    for ax, feat, title in zip(axes.flat, features, titles):
        sns.boxplot(
            data=df, x="hp_simple", y=feat, order=hp_order,
            palette=hp_colors, ax=ax, fliersize=0.5,
        )
        ax.set_title(title)
        ax.set_xlabel("HP Tag")
        ax.set_ylabel(feat)
        ax.tick_params(axis="x", rotation=30)

    fig.suptitle("Fig 04: Read-Level Methylation by HP Tag", fontsize=14, y=1.01)
    add_caption(fig, f"Source: Phase 1A read training table | {len(df):,} reads | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig04_read_methyl_by_hp_tag.png")


# ===== Figure 5: Methylation by truth label (TP/FP region) =====
def fig05_methyl_by_truth_label(df: pd.DataFrame) -> Path:
    """Read methylation in TP vs FP regions."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    features = ["methyl_mean", "methyl_std", "methyl_median",
                "methyl_high_fraction", "methyl_low_fraction", "methyl_na_fraction"]
    tp_palette = {"TP": COLOR_TP, "FP": COLOR_FP}

    stats_rows = []
    for ax, feat in zip(axes.flat, features):
        sns.violinplot(
            data=df, x="truth_status", y=feat, hue="truth_status",
            palette=tp_palette, inner="quartile", ax=ax, cut=0,
        )
        ax.set_title(feat)
        ax.set_xlabel("")

        # Stats
        tp_v = df.loc[df["truth_status"] == "TP", feat].dropna().values
        fp_v = df.loc[df["truth_status"] == "FP", feat].dropna().values
        mw = mannwhitney_test(tp_v, fp_v)
        es = compute_effect_size(tp_v, fp_v)
        auc_val = compute_auc((df["truth_status"] == "TP").astype(int).values,
                              df[feat].values)

        ax.text(0.98, 0.95, f"p={format_p(mw['p_value'])}\nd={es['d']:.3f}\nAUC={auc_val:.3f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7,
                bbox=dict(facecolor="white", alpha=0.8))

        stats_rows.append({
            "feature": feat, "mw_p": mw["p_value"], "effect_r": mw["effect_size_r"],
            "cohens_d": es["d"], "auc_tp": auc_val,
        })

    fig.suptitle("Fig 05: Read-Level Methylation in TP vs FP Regions", fontsize=14, y=1.01)
    add_caption(fig, f"Source: Phase 1A read training table | {len(df):,} reads | {RUN_DATE}")
    fig.tight_layout()

    # Save stats
    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(OUT_DIR / "fig05_tp_fp_stats.tsv", sep="\t", index=False)

    return save_figure(fig, OUT_DIR / "fig05_read_methyl_by_truth_label.png")


# ===== Figure 6: Methylation by VerificationClass =====
def fig06_methyl_by_verification_class(df: pd.DataFrame) -> Path:
    """Read methylation across four VerificationClass categories."""
    vc_order = list(LEGACY_CLASSES)
    vc_colors = {"Strong": "#1B5E20", "Weak": "#FFC107", "Subclone": "#FF9800", "Noise": "#D32F2F"}

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    features = ["methyl_mean", "methyl_std", "methyl_median",
                "methyl_high_fraction", "methyl_low_fraction", "num_cpg_observed"]

    for ax, feat in zip(axes.flat, features):
        present_order = [v for v in vc_order if v in df["VerificationClass"].unique()]
        sns.boxplot(
            data=df, x="VerificationClass", y=feat, order=present_order,
            palette=vc_colors, ax=ax, fliersize=0.3,
        )
        ax.set_title(feat)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=20)

    fig.suptitle("Fig 06: Read-Level Methylation by VerificationClass", fontsize=14, y=1.01)
    add_caption(fig, f"Source: Phase 1A read training table | {len(df):,} reads | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig06_read_methyl_by_verification_class.png")


# ===== Figure 7: Paired vs TO comparison (this dataset is paired-only) =====
def fig07_methyl_paired_vs_to_or_sample(df: pd.DataFrame) -> Path:
    """Since this dataset is paired-only, compare across samples instead.

    Show ALT-REF delta per sample to reveal sample-specific effects.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    features = ["methyl_mean", "methyl_std", "methyl_high_fraction", "methyl_low_fraction"]

    for ax, feat in zip(axes.flat, features):
        deltas = []
        for sample in SAMPLE_ORDER:
            sdf = df[df["sample"] == sample]
            alt_v = sdf.loc[sdf["alt_support"] == "ALT", feat].dropna()
            ref_v = sdf.loc[sdf["alt_support"] == "REF", feat].dropna()
            if len(alt_v) > 10 and len(ref_v) > 10:
                delta = alt_v.mean() - ref_v.mean()
                _, ci_lo, ci_hi = bootstrap_ci(np.mean, alt_v.values - np.mean(ref_v.values))
                deltas.append({"sample": sample, "delta": delta,
                               "ci_lo": ci_lo, "ci_hi": ci_hi})

        if deltas:
            ddf = pd.DataFrame(deltas)
            colors = [COLOR_TP if d >= 0 else COLOR_FP for d in ddf["delta"]]
            ax.barh(ddf["sample"], ddf["delta"], color=colors, alpha=0.7)
            ax.errorbar(ddf["delta"], ddf["sample"],
                        xerr=[ddf["delta"] - ddf["ci_lo"], ddf["ci_hi"] - ddf["delta"]],
                        fmt="none", color="black", capsize=3)
            ax.axvline(0, color="gray", linestyle="--", alpha=0.5)
        ax.set_title(f"ALT-REF Delta: {feat}")
        ax.set_xlabel("Delta (ALT - REF)")

    fig.suptitle("Fig 07: Sample-Specific ALT-REF Methylation Delta\n"
                 "(Paired-mode only; no TO data in read training table)",
                 fontsize=13, y=1.02)
    add_caption(fig, f"Source: Phase 1A read training table | paired-only | {RUN_DATE}")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "fig07_read_methyl_paired_vs_to.png")


# ===== Figure 8: Feature importance AUC ranking =====
def fig08_feature_auc_ranking(df: pd.DataFrame) -> Path:
    """AUC ranking of read-level features for ALT/REF classification and TP/FP classification."""
    all_features = [
        "methyl_mean", "methyl_std", "methyl_median",
        "methyl_high_fraction", "methyl_low_fraction",
        "methyl_na_fraction", "num_cpg_observed", "num_cpg_total",
        "cpg_obs_rate", "methyl_bimodal_score", "methyl_intermediate_fraction",
        "mapq",
    ]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # --- Panel A: ALT/REF classification ---
    ax = axes[0]
    results_alt = []
    for feat in all_features:
        if feat not in df.columns:
            continue
        vals = df[feat].values
        labels = df["is_alt"].values
        mask = np.isfinite(vals) & np.isfinite(labels)
        auc_val = compute_auc(labels[mask], vals[mask])
        # Also compute 1-AUC to handle reversed direction
        auc_directed = max(auc_val, 1 - auc_val) if not np.isnan(auc_val) else np.nan
        results_alt.append({"feature": feat, "auc_raw": auc_val, "auc_best": auc_directed})

    rdf_alt = pd.DataFrame(results_alt).sort_values("auc_best", ascending=True)
    colors_alt = [COLOR_ALT if a >= 0.5 else COLOR_REF for a in rdf_alt["auc_raw"]]
    ax.barh(rdf_alt["feature"], rdf_alt["auc_best"], color=colors_alt, alpha=0.8)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlim(0.45, max(0.65, rdf_alt["auc_best"].max() + 0.02))
    ax.set_title("AUC for ALT/REF Read Classification")
    ax.set_xlabel("AUC (best direction)")
    for i, row in rdf_alt.iterrows():
        ax.text(row["auc_best"] + 0.002, list(rdf_alt["feature"]).index(row["feature"]),
                f"{row['auc_best']:.3f}", va="center", fontsize=8)

    # --- Panel B: TP/FP region classification ---
    ax = axes[1]
    results_tp = []
    for feat in all_features:
        if feat not in df.columns:
            continue
        vals = df[feat].values
        labels = df["is_tp"].values
        mask = np.isfinite(vals) & np.isfinite(labels)
        auc_val = compute_auc(labels[mask], vals[mask])
        auc_directed = max(auc_val, 1 - auc_val) if not np.isnan(auc_val) else np.nan
        results_tp.append({"feature": feat, "auc_raw": auc_val, "auc_best": auc_directed})

    rdf_tp = pd.DataFrame(results_tp).sort_values("auc_best", ascending=True)
    colors_tp = [COLOR_TP if a >= 0.5 else COLOR_FP for a in rdf_tp["auc_raw"]]
    ax.barh(rdf_tp["feature"], rdf_tp["auc_best"], color=colors_tp, alpha=0.8)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlim(0.45, max(0.65, rdf_tp["auc_best"].max() + 0.02))
    ax.set_title("AUC for TP/FP Region Classification")
    ax.set_xlabel("AUC (best direction)")
    for i, row in rdf_tp.iterrows():
        ax.text(row["auc_best"] + 0.002, list(rdf_tp["feature"]).index(row["feature"]),
                f"{row['auc_best']:.3f}", va="center", fontsize=8)

    fig.suptitle("Fig 08: Read-Level Feature AUC Ranking", fontsize=14, y=1.02)
    add_caption(fig, f"Source: Phase 1A read training table | {len(df):,} reads | {RUN_DATE}")
    fig.tight_layout()

    # Save AUC tables
    pd.DataFrame(results_alt).sort_values("auc_best", ascending=False).to_csv(
        OUT_DIR / "fig08_auc_alt_ref.tsv", sep="\t", index=False)
    pd.DataFrame(results_tp).sort_values("auc_best", ascending=False).to_csv(
        OUT_DIR / "fig08_auc_tp_fp.tsv", sep="\t", index=False)

    return save_figure(fig, OUT_DIR / "fig08_read_feature_importance_for_alt_support.png")


# ===== Comprehensive Statistics =====
def compute_comprehensive_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Compute comprehensive statistics for all methylation features."""
    all_features = METHYL_FEATURES + ["cpg_obs_rate", "methyl_bimodal_score",
                                       "methyl_intermediate_fraction", "num_cpg_total", "mapq"]

    rows = []
    for feat in all_features:
        if feat not in df.columns:
            continue

        # ALT vs REF
        alt_v = df.loc[df["alt_support"] == "ALT", feat].dropna().values
        ref_v = df.loc[df["alt_support"] == "REF", feat].dropna().values
        mw_ar = mannwhitney_test(alt_v, ref_v)
        es_ar = compute_effect_size(alt_v, ref_v)
        auc_ar = compute_auc(df["is_alt"].values, df[feat].values)

        # TP vs FP
        tp_v = df.loc[df["truth_status"] == "TP", feat].dropna().values
        fp_v = df.loc[df["truth_status"] == "FP", feat].dropna().values
        mw_tf = mannwhitney_test(tp_v, fp_v)
        es_tf = compute_effect_size(tp_v, fp_v)
        auc_tf = compute_auc(df["is_tp"].values, df[feat].values)

        rows.append({
            "feature": feat,
            "alt_mean": float(np.mean(alt_v)),
            "ref_mean": float(np.mean(ref_v)),
            "alt_ref_delta": float(np.mean(alt_v) - np.mean(ref_v)),
            "alt_ref_mw_p": mw_ar["p_value"],
            "alt_ref_effect_r": mw_ar["effect_size_r"],
            "alt_ref_cohens_d": es_ar["d"],
            "alt_ref_d_ci_lo": es_ar["ci_lo"],
            "alt_ref_d_ci_hi": es_ar["ci_hi"],
            "alt_ref_auc": auc_ar,
            "tp_mean": float(np.mean(tp_v)),
            "fp_mean": float(np.mean(fp_v)),
            "tp_fp_delta": float(np.mean(tp_v) - np.mean(fp_v)),
            "tp_fp_mw_p": mw_tf["p_value"],
            "tp_fp_effect_r": mw_tf["effect_size_r"],
            "tp_fp_cohens_d": es_tf["d"],
            "tp_fp_d_ci_lo": es_tf["ci_lo"],
            "tp_fp_d_ci_hi": es_tf["ci_hi"],
            "tp_fp_auc": auc_tf,
        })

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(OUT_DIR / "comprehensive_stats.tsv", sep="\t", index=False)
    print(f"[O10] Saved comprehensive_stats.tsv ({len(rows)} features)")
    return stats_df


# ===== Main =====
def main():
    global READ_TABLE_PATH, OUT_DIR
    args = parse_args()
    READ_TABLE_PATH = args.read_table.resolve()
    OUT_DIR = args.output_dir.resolve()
    print("=" * 70)
    print(f"O10: Read-Level Methylation Feature Observation")
    print(f"Run date: {RUN_DATE}")
    print(f"Output: {OUT_DIR}")
    print("=" * 70)

    ensure_dir(OUT_DIR)

    # Load data
    df, verification_metadata = load_read_training_table(
        READ_TABLE_PATH,
        allow_unversioned_v1=args.allow_unversioned_v1,
    )
    df = add_derived_columns(df)

    # Feature statistics
    all_cols = METHYL_FEATURES + [
        "num_cpg_total", "cpg_obs_rate", "methyl_bimodal_score",
        "methyl_intermediate_fraction", "mapq",
    ]
    write_feature_statistics(df, all_cols, OUT_DIR / "feature_statistics.tsv")

    # Round context
    write_round_context(
        out_dir=OUT_DIR,
        observation_id=OBS_ID,
        title="Read-Level Methylation Feature Observation",
        description=(
            "Analyses read-level methylation features from Phase 1A "
            "read training tables to understand ALT/REF, TP/FP, HP-tag, "
            "and VerificationClass distributions at individual read granularity."
        ),
        script_path="scripts/analysis/build_observation_O10_read_level_methyl.py",
        data_sources=[{
            "name": "Phase 1A read training table (paired multibio 637)",
            "path": str(READ_TABLE_PATH),
            "rows": len(df),
            "cols": len(df.columns),
        }],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "verification_contract": verification_metadata,
            "samples": list(df["sample"].unique()),
            "modes": list(df["mode"].unique()),
            "n_alt": int((df["alt_support"] == "ALT").sum()),
            "n_ref": int((df["alt_support"] == "REF").sum()),
            "n_tp": int((df["truth_status"] == "TP").sum()),
            "n_fp": int((df["truth_status"] == "FP").sum()),
        },
    )

    # Generate all figures
    print("\n[O10] Generating figures...")
    fig01_methyl_mean_by_alt_support(df)
    fig02_methyl_na_fraction(df)
    fig03_cpg_count_distribution(df)
    fig04_methyl_by_hp_tag(df)
    fig05_methyl_by_truth_label(df)
    fig06_methyl_by_verification_class(df)
    fig07_methyl_paired_vs_to_or_sample(df)
    fig08_feature_auc_ranking(df)

    # Comprehensive stats
    print("\n[O10] Computing comprehensive statistics...")
    comp_stats = compute_comprehensive_stats(df)

    print("\n[O10] === Summary ===")
    print(comp_stats[["feature", "alt_ref_auc", "alt_ref_cohens_d",
                       "tp_fp_auc", "tp_fp_cohens_d"]].to_string(index=False))

    print(f"\n[O10] All outputs saved to: {OUT_DIR}")
    print("[O10] Done.")


if __name__ == "__main__":
    main()
