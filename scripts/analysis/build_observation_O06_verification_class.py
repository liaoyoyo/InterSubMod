#!/usr/bin/env python3
"""O6: VerificationClass Structure Analysis — Systematic Observation.

Analyses the four ISM VerificationClass categories (Strong, Subclone, Weak, Noise)
across samples, modes (paired/TO), truth labels, LOH subtypes, and feature boundaries.

Generates 6 figures, 5 data tables, and round_context.json.

Usage:
    python build_observation_O06_verification_class.py
"""

from __future__ import annotations

import sys
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")
from observation_common import *

from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats
from sklearn.metrics import cohen_kappa_score

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OBS_ID = "O06"
OBS_TITLE = "VerificationClass Structure Analysis"
OUT_DIR = OUTPUT_ROOT / "20260401_O06_verification_class"
DATA_DIR = OUT_DIR / "data"

VC_ORDER = ["Strong", "Subclone", "Weak", "Noise"]
VC_COLORS = {
    "Strong": "#1565C0",   # deep blue
    "Subclone": "#2E7D32", # green
    "Weak": "#F9A825",     # yellow
    "Noise": "#9E9E9E",    # grey
}
VC_PALETTE = [VC_COLORS[c] for c in VC_ORDER]

DL_ORDER = ["hp", "allele", "none"]
DL_COLORS = {"hp": "#1565C0", "allele": "#E65100", "none": "#9E9E9E"}

LOH_ORDER = ["None", "LOH_Strong", "LOH_Subclone", "LOH_Weak", "LOH_Noise"]

# Features for boundary analysis (fig06)
BOUNDARY_FEATURES = ["PairwiseMedianDist", "CramersV", "AlleleDelta", "HPMergedDelta"]

# Columns to load
USECOLS = [
    "VerificationClass", "DominantLabel", "Stability", "verification_class",
    "LOH_Subtype", "Quality_Tier", "Quality_Score",
    "truth_label", "variant_key", "sample", "mode", "dataset_id", "dataset_label",
    "PairwiseMedianDist", "CramersV", "AlleleDelta", "HPMergedDelta",
    "HPMergedSig", "AlleleSig", "HPMergedP", "AlleleP",
    "NumReads", "NumCpGs", "Potential_LOH", "Coverage_Multiple",
    "caller_af",
]

CAPTION = f"Source: all_region_rows.tsv.gz (748,391 rows) | {OBS_ID} | InterSubMod"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_data() -> pd.DataFrame:
    """Load master dataset with required columns."""
    df = load_master_dataset(usecols=USECOLS)

    # Ensure VerificationClass is categorical with correct order
    df["VerificationClass"] = pd.Categorical(df["VerificationClass"], categories=VC_ORDER, ordered=True)

    # Ensure LOH_Subtype has consistent values
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("None").replace("", "None")
    df["LOH_Subtype"] = pd.Categorical(df["LOH_Subtype"], categories=LOH_ORDER, ordered=True)

    # Ensure DominantLabel
    df["DominantLabel"] = df["DominantLabel"].fillna("none")

    print(f"[O06] Loaded {len(df):,} rows")
    print(f"[O06] VerificationClass distribution:")
    print(df["VerificationClass"].value_counts().to_string())
    print(f"[O06] truth_label distribution:")
    print(df["truth_label"].value_counts().to_string())
    return df


# ---------------------------------------------------------------------------
# Figure 1: VerificationClass composition by sample x mode
# ---------------------------------------------------------------------------

def fig01_verification_class_composition(df: pd.DataFrame) -> plt.Figure:
    """Stacked bar chart: VC proportions per dataset_id (14 datasets)."""
    print("\n[O06] === fig01: VerificationClass composition ===")

    # Compute proportions per dataset
    ct = pd.crosstab(df["dataset_id"], df["VerificationClass"], normalize="index")
    ct = ct.reindex(columns=VC_ORDER, fill_value=0)

    # Sort datasets by sample order and mode
    dataset_ids = sorted(ct.index, key=lambda x: (
        next((i for i, s in enumerate(SAMPLE_ORDER) if s.lower().replace("_", "") in x.replace("_", "")), 99),
        0 if "paired" in x else 1
    ))
    ct = ct.loc[dataset_ids]

    fig, ax = plt.subplots(figsize=(14, 6))
    bottom = np.zeros(len(ct))

    for vc in VC_ORDER:
        vals = ct[vc].values
        bars = ax.bar(range(len(ct)), vals, bottom=bottom, label=vc,
                      color=VC_COLORS[vc], edgecolor="white", linewidth=0.5)
        # Add percentage labels for large enough segments
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 0.05:
                ax.text(i, b + v / 2, f"{v:.0%}", ha="center", va="center",
                        fontsize=7, color="white" if vc in ("Strong", "Subclone") else "black")
        bottom += vals

    ax.set_xticks(range(len(ct)))
    ax.set_xticklabels([d.replace("_", "\n") for d in ct.index],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Proportion")
    ax.set_title(f"[{OBS_ID}] VerificationClass Composition by Dataset")
    ax.legend(loc="upper right", framealpha=0.9)
    ax.set_ylim(0, 1.02)

    # Add total counts on top
    for i, did in enumerate(ct.index):
        n = (df["dataset_id"] == did).sum()
        ax.text(i, 1.01, f"n={n:,}", ha="center", va="bottom", fontsize=6, color="#666")

    add_caption(fig, CAPTION)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 2: TP/FP rate by VerificationClass
# ---------------------------------------------------------------------------

def fig02_verification_class_tp_fp_rate(df: pd.DataFrame) -> Dict:
    """Bar chart: TP/FP counts and precision per VC, faceted by mode."""
    print("\n[O06] === fig02: VerificationClass TP/FP rates ===")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
    stats_rows = []
    chi2_results = {}

    for idx, mode_val in enumerate(["paired", "to"]):
        ax = axes[idx]
        sub = df[df["mode"] == mode_val].copy()

        # Count TP/FP per VC
        ct_counts = pd.crosstab(sub["VerificationClass"], sub["truth_label"])
        for col in ["TP", "FP"]:
            if col not in ct_counts.columns:
                ct_counts[col] = 0
        ct_counts = ct_counts.reindex(VC_ORDER, fill_value=0)

        tp_vals = ct_counts["TP"].values
        fp_vals = ct_counts["FP"].values
        total_vals = tp_vals + fp_vals
        precision_vals = np.where(total_vals > 0, tp_vals / total_vals, 0)
        fp_rate_vals = np.where(total_vals > 0, fp_vals / total_vals, 0)

        x = np.arange(len(VC_ORDER))
        w = 0.35
        bars_tp = ax.bar(x - w / 2, tp_vals, w, label="TP", color=COLOR_TP, edgecolor="white")
        bars_fp = ax.bar(x + w / 2, fp_vals, w, label="FP", color=COLOR_FP, edgecolor="white")

        # Add precision annotation
        for i in range(len(VC_ORDER)):
            y_max = max(tp_vals[i], fp_vals[i])
            ax.text(i, y_max + max(total_vals) * 0.02,
                    f"prec={precision_vals[i]:.2%}",
                    ha="center", va="bottom", fontsize=7, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels(VC_ORDER)
        ax.set_ylabel("Count")
        ax.set_title(f"{mode_val.upper()}")
        ax.legend(loc="upper right")

        # Chi-squared test: VC x truth_label
        contingency = ct_counts[["TP", "FP"]].values
        chi2_res = chi2_test(contingency)
        chi2_results[mode_val] = chi2_res
        ax.text(0.02, 0.97,
                f"Chi2={chi2_res['chi2']:.1f}, p={format_p(chi2_res['p_value'])}\n"
                f"Cramer's V={chi2_res['cramers_v']:.3f}",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

        # Collect stats
        for i, vc in enumerate(VC_ORDER):
            stats_rows.append({
                "mode": mode_val, "VerificationClass": vc,
                "TP": int(tp_vals[i]), "FP": int(fp_vals[i]),
                "Total": int(total_vals[i]),
                "Precision": round(float(precision_vals[i]), 4),
                "FP_Rate": round(float(fp_rate_vals[i]), 4),
            })

    fig.suptitle(f"[{OBS_ID}] VerificationClass — TP/FP Counts & Precision", fontsize=13, y=1.02)
    add_caption(fig, CAPTION)
    fig.tight_layout()

    return {"fig": fig, "stats": stats_rows, "chi2": chi2_results}


# ---------------------------------------------------------------------------
# Figure 3: DominantLabel by VerificationClass
# ---------------------------------------------------------------------------

def fig03_dominant_label_by_verification_class(df: pd.DataFrame) -> plt.Figure:
    """Stacked bar: DominantLabel distribution within each VC, split paired/TO."""
    print("\n[O06] === fig03: DominantLabel by VerificationClass ===")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)

    for idx, mode_val in enumerate(["paired", "to"]):
        ax = axes[idx]
        sub = df[df["mode"] == mode_val]

        ct = pd.crosstab(sub["VerificationClass"], sub["DominantLabel"], normalize="index")
        for dl in DL_ORDER:
            if dl not in ct.columns:
                ct[dl] = 0
        ct = ct.reindex(index=VC_ORDER, columns=DL_ORDER, fill_value=0)

        bottom = np.zeros(len(VC_ORDER))
        x = np.arange(len(VC_ORDER))

        for dl in DL_ORDER:
            vals = ct[dl].values
            ax.bar(x, vals, bottom=bottom, label=dl, color=DL_COLORS[dl],
                   edgecolor="white", linewidth=0.5)
            for i, (v, b) in enumerate(zip(vals, bottom)):
                if v > 0.05:
                    ax.text(i, b + v / 2, f"{v:.0%}", ha="center", va="center",
                            fontsize=7, color="white" if dl != "none" else "black")
            bottom += vals

        ax.set_xticks(x)
        ax.set_xticklabels(VC_ORDER)
        ax.set_ylabel("Proportion" if idx == 0 else "")
        ax.set_title(f"{mode_val.upper()}")
        ax.legend(loc="upper right", framealpha=0.9)
        ax.set_ylim(0, 1.05)

        # Add counts
        for i, vc in enumerate(VC_ORDER):
            n = (sub["VerificationClass"] == vc).sum()
            ax.text(i, 1.01, f"n={n:,}", ha="center", va="bottom", fontsize=6, color="#666")

    fig.suptitle(f"[{OBS_ID}] DominantLabel Distribution by VerificationClass", fontsize=13, y=1.02)
    add_caption(fig, CAPTION)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 4: Class transition matrix (paired -> TO)
# ---------------------------------------------------------------------------

def fig04_class_transition_paired_to_to(df: pd.DataFrame) -> Dict:
    """Heatmap: paired VC -> TO VC transition matrix via variant_key join."""
    print("\n[O06] === fig04: Class transition paired -> TO ===")

    paired = df[df["mode"] == "paired"][["variant_key", "VerificationClass"]].copy()
    paired = paired.rename(columns={"VerificationClass": "VC_paired"})
    to = df[df["mode"] == "to"][["variant_key", "VerificationClass"]].copy()
    to = to.rename(columns={"VerificationClass": "VC_to"})

    merged = paired.merge(to, on="variant_key", how="inner")
    print(f"[O06] Paired-TO matched variants: {len(merged):,}")

    # Build transition matrix
    trans = pd.crosstab(merged["VC_paired"], merged["VC_to"])
    trans = trans.reindex(index=VC_ORDER, columns=VC_ORDER, fill_value=0)

    # Also compute row-normalized (proportion)
    trans_norm = trans.div(trans.sum(axis=1), axis=0).fillna(0)

    # Cohen's Kappa
    # Encode to numeric for kappa
    vc_map = {vc: i for i, vc in enumerate(VC_ORDER)}
    paired_enc = merged["VC_paired"].map(vc_map).dropna().astype(int)
    to_enc = merged["VC_to"].map(vc_map).dropna().astype(int)
    common_idx = paired_enc.index.intersection(to_enc.index)
    kappa = cohen_kappa_score(paired_enc.loc[common_idx], to_enc.loc[common_idx])
    print(f"[O06] Cohen's Kappa (paired vs TO): {kappa:.4f}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Left: raw counts
    sns.heatmap(trans, annot=True, fmt="d", cmap="Blues", ax=axes[0],
                xticklabels=VC_ORDER, yticklabels=VC_ORDER)
    axes[0].set_xlabel("TO VerificationClass")
    axes[0].set_ylabel("Paired VerificationClass")
    axes[0].set_title("Counts")

    # Right: row-normalized proportions
    sns.heatmap(trans_norm, annot=True, fmt=".2f", cmap="YlOrRd", ax=axes[1],
                xticklabels=VC_ORDER, yticklabels=VC_ORDER, vmin=0, vmax=1)
    axes[1].set_xlabel("TO VerificationClass")
    axes[1].set_ylabel("Paired VerificationClass")
    axes[1].set_title("Row-normalized (proportion)")

    fig.suptitle(
        f"[{OBS_ID}] Paired → TO VerificationClass Transition\n"
        f"(n={len(merged):,} matched variants, Cohen's κ = {kappa:.3f})",
        fontsize=12, y=1.04
    )
    add_caption(fig, CAPTION)
    fig.tight_layout()

    return {"fig": fig, "trans": trans, "trans_norm": trans_norm, "kappa": kappa, "n_matched": len(merged)}


# ---------------------------------------------------------------------------
# Figure 5: VerificationClass x LOH_Subtype
# ---------------------------------------------------------------------------

def fig05_verification_class_by_loh_subtype(df: pd.DataFrame) -> plt.Figure:
    """Heatmap: VC x LOH_Subtype counts, split by paired/TO."""
    print("\n[O06] === fig05: VerificationClass by LOH_Subtype ===")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for idx, mode_val in enumerate(["paired", "to"]):
        ax = axes[idx]
        sub = df[df["mode"] == mode_val]

        ct = pd.crosstab(sub["VerificationClass"], sub["LOH_Subtype"])
        ct = ct.reindex(index=VC_ORDER, columns=LOH_ORDER, fill_value=0)

        sns.heatmap(ct, annot=True, fmt="d", cmap="YlGnBu", ax=ax,
                    xticklabels=[s.replace("LOH_", "") for s in LOH_ORDER],
                    yticklabels=VC_ORDER)
        ax.set_xlabel("LOH Subtype")
        ax.set_ylabel("VerificationClass" if idx == 0 else "")
        ax.set_title(f"{mode_val.upper()} (n={len(sub):,})")

    fig.suptitle(f"[{OBS_ID}] VerificationClass × LOH Subtype", fontsize=13, y=1.02)
    add_caption(fig, CAPTION)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 6: Class boundary features (violin plots)
# ---------------------------------------------------------------------------

def fig06_class_boundary_features(df: pd.DataFrame) -> Dict:
    """Violin plots: 4 features x 4 VCs, with Kruskal-Wallis tests."""
    print("\n[O06] === fig06: Class boundary features ===")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    kw_results = {}

    for ax_idx, (feat, ax) in enumerate(zip(BOUNDARY_FEATURES, axes.flat)):
        # Get data for each VC
        plot_data = df[["VerificationClass", feat]].dropna(subset=[feat]).copy()
        plot_data[feat] = pd.to_numeric(plot_data[feat], errors="coerce")
        plot_data = plot_data.dropna(subset=[feat])

        # Violin plot
        sns.violinplot(data=plot_data, x="VerificationClass", y=feat,
                       order=VC_ORDER, palette=VC_COLORS, ax=ax,
                       inner="quartile", cut=0, linewidth=0.8)

        # Kruskal-Wallis test across 4 groups
        groups = [plot_data.loc[plot_data["VerificationClass"] == vc, feat].values
                  for vc in VC_ORDER]
        groups_valid = [g for g in groups if len(g) >= 2]
        if len(groups_valid) >= 2:
            kw_stat, kw_p = sp_stats.kruskal(*groups_valid)
            # Effect size: eta-squared = (H - k + 1) / (n - k)
            n_total = sum(len(g) for g in groups_valid)
            k = len(groups_valid)
            eta_sq = (kw_stat - k + 1) / (n_total - k) if n_total > k else 0
            eta_sq = max(0, eta_sq)
        else:
            kw_stat, kw_p, eta_sq = float("nan"), float("nan"), float("nan")

        kw_results[feat] = {
            "H_statistic": float(kw_stat),
            "p_value": float(kw_p),
            "eta_squared": float(eta_sq),
        }

        ax.set_title(f"{feat}\nKW H={kw_stat:.1f}, p={format_p(kw_p)}, η²={eta_sq:.4f}",
                      fontsize=10)
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelsize=8)

    fig.suptitle(f"[{OBS_ID}] Feature Distributions by VerificationClass (Strong vs Weak boundary)",
                 fontsize=12, y=1.01)
    add_caption(fig, CAPTION)
    fig.tight_layout()

    return {"fig": fig, "kw_results": kw_results}


# ---------------------------------------------------------------------------
# Data output helpers
# ---------------------------------------------------------------------------

def write_source_summary(df: pd.DataFrame) -> None:
    """Write source_summary.tsv."""
    rows = []
    for mode_val in ["paired", "to"]:
        sub = df[df["mode"] == mode_val]
        for vc in VC_ORDER:
            vc_sub = sub[sub["VerificationClass"] == vc]
            rows.append({
                "mode": mode_val,
                "VerificationClass": vc,
                "count": len(vc_sub),
                "pct_of_mode": round(safe_div(len(vc_sub), len(sub)) * 100, 2),
                "TP": (vc_sub["truth_label"] == "TP").sum(),
                "FP": (vc_sub["truth_label"] == "FP").sum(),
            })

    summary_df = pd.DataFrame(rows)
    out_path = DATA_DIR / "source_summary.tsv"
    summary_df.to_csv(out_path, sep="\t", index=False)
    print(f"[O06] Written: {out_path.name}")


def write_verification_class_by_sample_mode(df: pd.DataFrame) -> None:
    """Write per-sample/mode VC counts."""
    ct = pd.crosstab(
        [df["sample"], df["mode"]],
        df["VerificationClass"],
    ).reset_index()
    ct.columns.name = None
    out_path = DATA_DIR / "verification_class_by_sample_mode.tsv"
    ct.to_csv(out_path, sep="\t", index=False)
    print(f"[O06] Written: {out_path.name}")


def write_class_tp_fp_rates(stats_rows: List[Dict]) -> None:
    """Write TP/FP rates TSV."""
    out_path = DATA_DIR / "class_tp_fp_rates.tsv"
    pd.DataFrame(stats_rows).to_csv(out_path, sep="\t", index=False)
    print(f"[O06] Written: {out_path.name}")


def write_class_transition_matrix(trans: pd.DataFrame, trans_norm: pd.DataFrame,
                                   kappa: float, n_matched: int) -> None:
    """Write transition matrix TSV."""
    out_path = DATA_DIR / "class_transition_matrix.tsv"
    with open(out_path, "w") as f:
        f.write(f"# Cohen's Kappa = {kappa:.4f}\n")
        f.write(f"# n_matched = {n_matched}\n")
        f.write(f"# Rows = Paired VC, Columns = TO VC\n")
        f.write("# --- Raw Counts ---\n")
    trans.to_csv(out_path, sep="\t", mode="a")
    with open(out_path, "a") as f:
        f.write("\n# --- Row-Normalized Proportions ---\n")
    trans_norm.to_csv(out_path, sep="\t", mode="a")
    print(f"[O06] Written: {out_path.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"{'=' * 70}")
    print(f"  O6: {OBS_TITLE}")
    print(f"  Output: {OUT_DIR}")
    print(f"{'=' * 70}")

    setup_plot_style()
    ensure_dir(OUT_DIR)
    ensure_dir(DATA_DIR)

    # Load data
    df = load_data()

    # ----- Data outputs -----
    write_source_summary(df)
    write_verification_class_by_sample_mode(df)

    # Feature statistics
    stat_cols = [
        "VerificationClass", "DominantLabel", "Stability", "LOH_Subtype",
        "Quality_Tier", "Quality_Score",
        "PairwiseMedianDist", "CramersV", "AlleleDelta", "HPMergedDelta",
        "NumReads", "NumCpGs", "caller_af",
    ]
    write_feature_statistics(df, stat_cols, DATA_DIR / "feature_statistics.tsv")

    # ----- Figure 1 -----
    fig1 = fig01_verification_class_composition(df)
    save_figure(fig1, OUT_DIR / "fig01_verification_class_composition.png")

    # ----- Figure 2 -----
    res2 = fig02_verification_class_tp_fp_rate(df)
    save_figure(res2["fig"], OUT_DIR / "fig02_verification_class_tp_fp_rate.png")
    write_class_tp_fp_rates(res2["stats"])

    # ----- Figure 3 -----
    fig3 = fig03_dominant_label_by_verification_class(df)
    save_figure(fig3, OUT_DIR / "fig03_dominant_label_by_verification_class.png")

    # ----- Figure 4 -----
    res4 = fig04_class_transition_paired_to_to(df)
    save_figure(res4["fig"], OUT_DIR / "fig04_class_transition_paired_to_to.png")
    write_class_transition_matrix(res4["trans"], res4["trans_norm"],
                                   res4["kappa"], res4["n_matched"])

    # ----- Figure 5 -----
    fig5 = fig05_verification_class_by_loh_subtype(df)
    save_figure(fig5, OUT_DIR / "fig05_verification_class_by_loh_subtype.png")

    # ----- Figure 6 -----
    res6 = fig06_class_boundary_features(df)
    save_figure(res6["fig"], OUT_DIR / "fig06_class_boundary_features.png")

    # ----- Round context -----
    write_round_context(
        out_dir=OUT_DIR,
        observation_id=OBS_ID,
        title=OBS_TITLE,
        description=(
            "Systematic analysis of the four ISM VerificationClass categories "
            "(Strong, Subclone, Weak, Noise) across samples, modes, truth labels, "
            "LOH subtypes, and feature boundaries. Includes paired-TO transition "
            "analysis and class boundary feature distributions."
        ),
        script_path="scripts/analysis/build_observation_O06_verification_class.py",
        data_sources=[{
            "name": "all_region_rows.tsv.gz",
            "path": str(MASTER_DATASET_PATH),
            "rows": len(df),
        }],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "figures": [
                "fig01_verification_class_composition.png",
                "fig02_verification_class_tp_fp_rate.png",
                "fig03_dominant_label_by_verification_class.png",
                "fig04_class_transition_paired_to_to.png",
                "fig05_verification_class_by_loh_subtype.png",
                "fig06_class_boundary_features.png",
            ],
            "statistical_tests": {
                "fig02_chi2_paired": res2["chi2"].get("paired", {}),
                "fig02_chi2_to": res2["chi2"].get("to", {}),
                "fig04_cohens_kappa": res4["kappa"],
                "fig04_n_matched_variants": res4["n_matched"],
                "fig06_kruskal_wallis": res6["kw_results"],
            },
        },
    )

    print(f"\n{'=' * 70}")
    print(f"  O6 complete. Output: {OUT_DIR}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
