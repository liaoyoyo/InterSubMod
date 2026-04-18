#!/usr/bin/env python3
"""O1: Master Dataset Global Distribution Observation.

Generates 12 figures and 4 data files that characterize the full 748K-row
master dataset: numeric distributions, categorical compositions, correlations,
missing-data patterns, discriminative power (AUC), and paired-vs-TO shifts.

Usage:
    python3 build_observation_O01_master_distribution.py
"""

from __future__ import annotations

import sys
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")
from observation_common import *

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

# ── Column lists ──────────────────────────────────────────────────────────
NUMERIC_COLS_20 = [
    "NumReads", "NumCpGs", "GlobalP", "CramersV",
    "HeuristicScore", "PairwiseMeanDist", "PairwiseMedianDist",
    "HP_Ratio", "HP1FamilyN", "HP2FamilyN",
    "NHP0", "NHP3", "Quality_Score", "Coverage_Multiple",
    "caller_af", "caller_gq", "AlleleDelta", "HPMergedDelta",
    "effective_hp_reads", "hp_assign_rate",
]

CATEGORICAL_COLS = [
    "VerificationClass", "Quality_Tier", "Coverage_Category",
    "LOH_Subtype", "DominantLabel", "truth_label", "mode",
]

CAPTION_SRC = (
    "Source: all_region_rows.tsv.gz "
    "(20260327_loh_round1_cross_sample_audit) | "
)


# ── Helper: AUC via scipy (fallback-safe) ─────────────────────────────────
def _auc_scipy(y_true, y_score):
    """ROC AUC using Mann-Whitney U statistic (no sklearn dependency)."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_score)
    y_true, y_score = y_true[mask], y_score[mask]
    classes = np.unique(y_true)
    if len(classes) < 2 or len(y_true) < 10:
        return float("nan")
    try:
        return compute_auc(y_true, y_score)
    except Exception:
        # fallback
        pos = y_score[y_true == 1]
        neg = y_score[y_true == 0]
        u, _ = sp_stats.mannwhitneyu(pos, neg, alternative="two-sided")
        return u / (len(pos) * len(neg))


# ══════════════════════════════════════════════════════════════════════════
# Figure generators
# ══════════════════════════════════════════════════════════════════════════


def fig01_numeric_distributions(df, out):
    """20 core numeric feature histograms (4x5 subplot panel)."""
    fig, axes = plt.subplots(4, 5, figsize=(22, 16))
    axes = axes.flatten()
    for i, col in enumerate(NUMERIC_COLS_20):
        ax = axes[i]
        if col not in df.columns:
            ax.set_title(f"{col}\n(missing)", fontsize=9)
            ax.axis("off")
            continue
        vals = df[col].dropna()
        if len(vals) == 0:
            ax.set_title(f"{col}\n(all NaN)", fontsize=9)
            ax.axis("off")
            continue
        # Clip extreme tails for readability (1st/99th percentile)
        lo, hi = np.nanpercentile(vals, [1, 99])
        clipped = vals.clip(lo, hi)
        ax.hist(clipped, bins=60, color="#546E7A", edgecolor="none", alpha=0.8)
        ax.set_title(col, fontsize=10, fontweight="bold")
        ax.set_xlabel("")
        # Annotate stats
        ax.text(
            0.97, 0.95,
            f"n={len(vals):,}\nmed={vals.median():.3g}\nmean={vals.mean():.3g}",
            transform=ax.transAxes, fontsize=7, va="top", ha="right",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
        )
    fig.suptitle("Fig 01 — Numeric Feature Distributions (20 core columns)", fontsize=14, y=1.01)
    add_caption(fig, CAPTION_SRC + f"{len(df):,} rows", y=-0.01)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig01_numeric_feature_distributions.png")


def fig02_categorical_barplots(df, out):
    """Categorical feature count bar plots."""
    n_cats = len(CATEGORICAL_COLS)
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()
    for i, col in enumerate(CATEGORICAL_COLS):
        ax = axes[i]
        if col not in df.columns:
            ax.set_title(f"{col} (missing)")
            ax.axis("off")
            continue
        vc = df[col].value_counts().head(15)
        bars = ax.barh(range(len(vc)), vc.values, color="#78909C", edgecolor="none")
        ax.set_yticks(range(len(vc)))
        ax.set_yticklabels(vc.index, fontsize=8)
        ax.invert_yaxis()
        ax.set_title(col, fontsize=11, fontweight="bold")
        ax.set_xlabel("Count")
        for bar, v in zip(bars, vc.values):
            ax.text(bar.get_width() + max(vc.values) * 0.01, bar.get_y() + bar.get_height() / 2,
                    f"{v:,}", va="center", fontsize=7)
    # hide unused axes
    for j in range(len(CATEGORICAL_COLS), len(axes)):
        axes[j].axis("off")
    fig.suptitle("Fig 02 — Categorical Feature Bar Plots", fontsize=14, y=1.01)
    add_caption(fig, CAPTION_SRC + f"{len(df):,} rows", y=-0.01)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig02_categorical_feature_barplots.png")


def fig03_sample_mode_heatmap(df, out):
    """sample x mode region counts heatmap."""
    ct = df.groupby(["sample", "mode"]).size().unstack(fill_value=0)
    # reindex for consistent order
    samples = [s for s in SAMPLE_ORDER if s in ct.index]
    modes = [m for m in MODE_ORDER if m in ct.columns]
    ct = ct.reindex(index=samples, columns=modes, fill_value=0)

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(ct, annot=True, fmt=",", cmap="YlOrRd", ax=ax, linewidths=0.5)
    ax.set_title("Fig 03 — Region Counts per Sample × Mode", fontsize=13)
    ax.set_ylabel("Sample")
    ax.set_xlabel("Mode")
    add_caption(fig, CAPTION_SRC + f"{len(df):,} rows", y=-0.02)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig03_sample_mode_region_counts.png")


def fig04_truth_composition(df, out):
    """TP/FP stacked bar per sample/mode."""
    sub = df[df["truth_label"].isin(["TP", "FP"])].copy()
    ct = sub.groupby(["sample", "mode", "truth_label"]).size().unstack(fill_value=0)
    ct["total"] = ct.sum(axis=1)
    for lbl in ["TP", "FP"]:
        if lbl not in ct.columns:
            ct[lbl] = 0
        ct[f"{lbl}_pct"] = ct[lbl] / ct["total"] * 100

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=False)
    for idx_m, mode in enumerate(MODE_ORDER):
        ax = axes[idx_m]
        if mode not in ct.index.get_level_values("mode"):
            ax.set_title(f"{mode} (no data)")
            continue
        sub_ct = ct.xs(mode, level="mode")
        samples = [s for s in SAMPLE_ORDER if s in sub_ct.index]
        sub_ct = sub_ct.reindex(samples)
        tp_pct = sub_ct["TP_pct"].values
        fp_pct = sub_ct["FP_pct"].values
        x = np.arange(len(samples))
        ax.bar(x, tp_pct, color=COLOR_TP, label="TP", edgecolor="none")
        ax.bar(x, fp_pct, bottom=tp_pct, color=COLOR_FP, label="FP", edgecolor="none")
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Percentage (%)")
        ax.set_title(f"{mode.upper()} — TP/FP Composition", fontsize=12)
        ax.legend(loc="upper right")
        ax.set_ylim(0, 105)
        # annotate counts
        for i, s in enumerate(samples):
            tp_n = int(sub_ct.loc[s, "TP"]) if s in sub_ct.index else 0
            fp_n = int(sub_ct.loc[s, "FP"]) if s in sub_ct.index else 0
            ax.text(i, 102, f"{tp_n}/{fp_n}", ha="center", fontsize=7, color="#555")
    fig.suptitle("Fig 04 — TP/FP Composition by Sample and Mode", fontsize=14, y=1.02)
    add_caption(fig, CAPTION_SRC + f"TP+FP rows = {len(sub):,}", y=-0.02)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig04_truth_label_composition.png")


def fig05_correlation_heatmap(df, out):
    """Spearman correlation matrix for 20 numeric columns."""
    cols_present = [c for c in NUMERIC_COLS_20 if c in df.columns]
    corr = df[cols_present].corr(method="spearman")

    fig, ax = plt.subplots(figsize=(14, 12))
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(
        corr, mask=mask, annot=True, fmt=".2f", cmap="RdBu_r",
        center=0, vmin=-1, vmax=1, ax=ax, linewidths=0.3,
        annot_kws={"fontsize": 7},
    )
    ax.set_title("Fig 05 — Spearman Correlation Matrix (20 core numeric features)", fontsize=13)
    add_caption(fig, CAPTION_SRC + f"{len(df):,} rows", y=-0.02)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig05_feature_correlation_heatmap.png")


def fig06_missing_pattern(df, out):
    """Missing / NaN / zero proportion per column."""
    all_cols = NUMERIC_COLS_20 + CATEGORICAL_COLS
    rows_list = []
    for col in all_cols:
        if col not in df.columns:
            continue
        s = df[col]
        n = len(s)
        na_count = int(s.isna().sum())
        zero_count = int((s == 0).sum()) if pd.api.types.is_numeric_dtype(s) else 0
        rows_list.append({
            "feature": col,
            "missing_pct": na_count / n * 100,
            "zero_pct": zero_count / n * 100,
        })
    miss_df = pd.DataFrame(rows_list).set_index("feature")

    fig, ax = plt.subplots(figsize=(14, 8))
    x = np.arange(len(miss_df))
    w = 0.35
    ax.bar(x - w / 2, miss_df["missing_pct"], w, label="Missing/NaN %", color="#EF5350")
    ax.bar(x + w / 2, miss_df["zero_pct"], w, label="Zero %", color="#FFA726")
    ax.set_xticks(x)
    ax.set_xticklabels(miss_df.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Percentage (%)")
    ax.set_title("Fig 06 — Missing / NaN / Zero Proportion per Feature", fontsize=13)
    ax.legend()
    add_caption(fig, CAPTION_SRC + f"{len(df):,} rows", y=-0.05)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig06_missing_nan_pattern.png")


def fig07_qs_violin_by_sample(df, out):
    """Quality Score violin plots by sample, split by mode."""
    sub = df[df["Quality_Score"].notna()].copy()
    fig, ax = plt.subplots(figsize=(16, 7))
    samples = [s for s in SAMPLE_ORDER if s in sub["sample"].unique()]
    sub_plot = sub[sub["sample"].isin(samples)]
    sns.violinplot(
        data=sub_plot, x="sample", y="Quality_Score", hue="mode",
        order=samples, hue_order=MODE_ORDER,
        palette=MODE_PALETTE, split=True, inner="quartile",
        cut=0, ax=ax, density_norm="width",
    )
    ax.set_title("Fig 07 — Quality Score Distribution by Sample (paired vs TO)", fontsize=13)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Quality Score")
    ax.legend(title="Mode", loc="upper right")
    add_caption(fig, CAPTION_SRC + f"non-NaN QS rows = {len(sub_plot):,}", y=-0.03)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig07_quality_score_distribution_by_sample.png")


def fig08_numreads_numcpgs_scatter(df, out):
    """NumReads vs NumCpGs scatter, coloured by truth."""
    sub = df[df["truth_label"].isin(["TP", "FP"]) & df["NumReads"].notna() & df["NumCpGs"].notna()].copy()
    # Downsample if too large (> 50K points per class)
    MAX_PER_CLASS = 50000
    parts = []
    for lbl in ["TP", "FP"]:
        part = sub[sub["truth_label"] == lbl]
        if len(part) > MAX_PER_CLASS:
            part = part.sample(MAX_PER_CLASS, random_state=42)
        parts.append(part)
    sub_s = pd.concat(parts)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        ms = sub_s[sub_s["mode"] == mode]
        for lbl, color in TRUTH_PALETTE.items():
            d = ms[ms["truth_label"] == lbl]
            ax.scatter(d["NumReads"], d["NumCpGs"], c=color, alpha=0.08,
                       s=6, label=f"{lbl} (n={len(d):,})", edgecolors="none")
        ax.set_xlabel("NumReads")
        ax.set_ylabel("NumCpGs")
        ax.set_title(f"{mode.upper()}", fontsize=12)
        ax.legend(fontsize=8, markerscale=4)
    fig.suptitle("Fig 08 — NumReads vs NumCpGs (by truth label)", fontsize=14, y=1.02)
    add_caption(fig, CAPTION_SRC + f"sampled {len(sub_s):,} / {len(sub):,} TP+FP rows", y=-0.02)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig08_numreads_numcpgs_scatter.png")


def fig09_coverage_category(df, out):
    """Coverage category stacked bar per sample."""
    sub = df[df["Coverage_Category"].notna()].copy()
    ct = sub.groupby(["sample", "Coverage_Category"]).size().unstack(fill_value=0)
    # Normalise to pct
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
    samples = [s for s in SAMPLE_ORDER if s in ct_pct.index]
    ct_pct = ct_pct.reindex(samples)
    cat_order = sorted(ct_pct.columns)

    fig, ax = plt.subplots(figsize=(14, 7))
    bottom = np.zeros(len(samples))
    cmap = plt.cm.Set2(np.linspace(0, 1, len(cat_order)))
    for i, cat in enumerate(cat_order):
        vals = ct_pct[cat].values if cat in ct_pct.columns else np.zeros(len(samples))
        ax.bar(samples, vals, bottom=bottom, label=cat, color=cmap[i], edgecolor="none")
        bottom += vals
    ax.set_ylabel("Percentage (%)")
    ax.set_title("Fig 09 — Coverage Category Composition by Sample", fontsize=13)
    ax.legend(title="Coverage Category", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    ax.set_xticklabels(samples, rotation=30, ha="right")
    add_caption(fig, CAPTION_SRC + f"{len(sub):,} rows with Coverage_Category", y=-0.04)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig09_coverage_category_composition.png")


def fig10_discriminative_auc(df, out):
    """Per-feature AUC for TP vs FP, split by mode. Bar chart."""
    sub = df[df["truth_label"].isin(["TP", "FP"])].copy()
    sub["y"] = encode_truth_binary(sub)

    results = []
    for mode in MODE_ORDER:
        ms = sub[sub["mode"] == mode]
        for col in NUMERIC_COLS_20:
            if col not in ms.columns:
                continue
            auc_val = _auc_scipy(ms["y"].values, ms[col].values)
            results.append({"feature": col, "mode": mode, "auc": auc_val})
    auc_df = pd.DataFrame(results)

    # Save data
    # (will be saved in main())

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)
    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        sub_a = auc_df[auc_df["mode"] == mode].copy()
        sub_a["auc_abs"] = (sub_a["auc"] - 0.5).abs() + 0.5  # distance from 0.5
        sub_a = sub_a.sort_values("auc_abs", ascending=True)
        colors = [COLOR_TP if a >= 0.5 else COLOR_FP for a in sub_a["auc"]]
        ax.barh(range(len(sub_a)), sub_a["auc"].values, color=colors, edgecolor="none")
        ax.set_yticks(range(len(sub_a)))
        ax.set_yticklabels(sub_a["feature"].values, fontsize=8)
        ax.axvline(0.5, color="grey", linestyle="--", alpha=0.5)
        ax.set_xlabel("AUC")
        ax.set_title(f"{mode.upper()} — Feature AUC (TP vs FP)", fontsize=12)
        ax.set_xlim(0, 1)
        # annotate values
        for j, (feat, auc_v) in enumerate(zip(sub_a["feature"], sub_a["auc"])):
            if not np.isnan(auc_v):
                ax.text(auc_v + 0.01, j, f"{auc_v:.3f}", va="center", fontsize=7)
    fig.suptitle("Fig 10 — Top Discriminative Features (AUC for TP vs FP)", fontsize=14, y=1.02)
    add_caption(fig, CAPTION_SRC + f"TP+FP rows = {len(sub):,}", y=-0.02)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig10_top_discriminative_features.png")
    return auc_df


def fig11_paired_vs_to_shifts(df, out):
    """Paired vs TO feature distribution comparison for 6 key features."""
    # Pick 6 features that typically differ most between paired and TO
    key_features = ["Quality_Score", "NumReads", "CramersV",
                    "hp_assign_rate", "caller_af", "HPMergedDelta"]
    key_features = [c for c in key_features if c in df.columns]

    n_feat = len(key_features)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    for i, col in enumerate(key_features):
        ax = axes[i]
        sub = df[df["mode"].isin(MODE_ORDER) & df[col].notna()].copy()
        sns.violinplot(
            data=sub, x="mode", y=col, order=MODE_ORDER,
            palette=MODE_PALETTE, inner="quartile", cut=0,
            density_norm="width", ax=ax,
        )
        # Mann-Whitney test
        paired_vals = sub.loc[sub["mode"] == "paired", col].dropna().values
        to_vals = sub.loc[sub["mode"] == "to", col].dropna().values
        mw = mannwhitney_test(paired_vals, to_vals)
        ax.set_title(f"{col}\n(MW p={format_p(mw['p_value'])}, r={mw['effect_size_r']:.3f})",
                     fontsize=10)
        ax.set_xlabel("")
    for j in range(len(key_features), len(axes)):
        axes[j].axis("off")
    fig.suptitle("Fig 11 — Paired vs TO Feature Distribution Shifts", fontsize=14, y=1.02)
    add_caption(fig, CAPTION_SRC + f"{len(df):,} rows", y=-0.02)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig11_paired_vs_to_feature_shifts.png")


def fig12_dataset_summary_heatmap(df, out):
    """Dataset-level TP/FP/TP-rate/F1 summary heatmap."""
    sub = df[df["truth_label"].isin(["TP", "FP"])].copy()
    grp = sub.groupby("dataset_id")["truth_label"].value_counts().unstack(fill_value=0)
    if "TP" not in grp.columns:
        grp["TP"] = 0
    if "FP" not in grp.columns:
        grp["FP"] = 0
    grp["total"] = grp["TP"] + grp["FP"]
    grp["TP_rate"] = (grp["TP"] / grp["total"] * 100).round(1)
    # F1: precision = TP/(TP+FP), recall is unknown (no FN), so show precision as proxy
    grp["precision"] = (grp["TP"] / grp["total"] * 100).round(1)

    display_cols = ["TP", "FP", "total", "TP_rate"]
    fig, ax = plt.subplots(figsize=(10, max(6, len(grp) * 0.5)))
    sns.heatmap(
        grp[display_cols].astype(float), annot=True, fmt=".0f",
        cmap="YlGnBu", ax=ax, linewidths=0.5,
    )
    ax.set_title("Fig 12 — Dataset Summary (TP / FP / TP Rate %)", fontsize=13)
    ax.set_ylabel("dataset_id")
    add_caption(fig, CAPTION_SRC + f"TP+FP rows = {len(sub):,}", y=-0.03)
    fig.tight_layout()
    save_figure(fig, out / "figures" / "fig12_dataset_id_summary_table.png")
    return grp


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    setup_plot_style()
    OUT = OUTPUT_ROOT / "20260401_O01_master_distribution"
    ensure_dir(OUT / "figures")
    ensure_dir(OUT / "data")

    print("=" * 70)
    print("O1: Master Dataset Global Distribution Observation")
    print("=" * 70)

    # ── Load data ─────────────────────────────────────────────────────
    df = load_master_dataset()
    n_rows, n_cols = df.shape
    print(f"Dataset shape: {n_rows:,} rows × {n_cols} columns\n")

    # ── Fig 01: Numeric distributions ─────────────────────────────────
    print("[1/12] Numeric feature distributions ...")
    fig01_numeric_distributions(df, OUT)

    # ── Fig 02: Categorical bar plots ─────────────────────────────────
    print("[2/12] Categorical feature bar plots ...")
    fig02_categorical_barplots(df, OUT)

    # ── Fig 03: Sample × Mode region counts ───────────────────────────
    print("[3/12] Sample × Mode region counts heatmap ...")
    fig03_sample_mode_heatmap(df, OUT)

    # ── Fig 04: Truth label composition ───────────────────────────────
    print("[4/12] Truth label composition ...")
    fig04_truth_composition(df, OUT)

    # ── Fig 05: Correlation heatmap ───────────────────────────────────
    print("[5/12] Feature correlation heatmap ...")
    fig05_correlation_heatmap(df, OUT)

    # ── Fig 06: Missing/NaN pattern ───────────────────────────────────
    print("[6/12] Missing/NaN pattern ...")
    fig06_missing_pattern(df, OUT)

    # ── Fig 07: QS violin by sample ───────────────────────────────────
    print("[7/12] Quality Score violin by sample ...")
    fig07_qs_violin_by_sample(df, OUT)

    # ── Fig 08: NumReads vs NumCpGs scatter ───────────────────────────
    print("[8/12] NumReads vs NumCpGs scatter ...")
    fig08_numreads_numcpgs_scatter(df, OUT)

    # ── Fig 09: Coverage category composition ─────────────────────────
    print("[9/12] Coverage category composition ...")
    fig09_coverage_category(df, OUT)

    # ── Fig 10: Discriminative AUC ────────────────────────────────────
    print("[10/12] Discriminative AUC ...")
    auc_df = fig10_discriminative_auc(df, OUT)

    # ── Fig 11: Paired vs TO shifts ───────────────────────────────────
    print("[11/12] Paired vs TO feature shifts ...")
    fig11_paired_vs_to_shifts(df, OUT)

    # ── Fig 12: Dataset summary ───────────────────────────────────────
    print("[12/12] Dataset summary heatmap ...")
    dataset_grp = fig12_dataset_summary_heatmap(df, OUT)

    # ══════════════════════════════════════════════════════════════════
    # Data outputs
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Data outputs ---")

    # source_summary.tsv
    sub_tl = df[df["truth_label"].isin(["TP", "FP"])]
    src_summary = sub_tl.groupby("dataset_id")["truth_label"].value_counts().unstack(fill_value=0)
    for c in ["TP", "FP"]:
        if c not in src_summary.columns:
            src_summary[c] = 0
    src_summary["total"] = src_summary.sum(axis=1)
    src_summary.to_csv(OUT / "data" / "source_summary.tsv", sep="\t")
    print(f"  -> source_summary.tsv ({len(src_summary)} datasets)")

    # feature_statistics.tsv
    write_feature_statistics(df, NUMERIC_COLS_20, OUT / "data" / "feature_statistics.tsv")

    # feature_auc_by_mode.tsv
    auc_df.to_csv(OUT / "data" / "feature_auc_by_mode.tsv", sep="\t", index=False)
    print(f"  -> feature_auc_by_mode.tsv ({len(auc_df)} rows)")

    # round_context.json
    write_round_context(
        out_dir=OUT,
        observation_id="O01",
        title="Master Dataset Global Distribution Observation",
        description=(
            "Characterise the full 748K-row master dataset: numeric and "
            "categorical distributions, correlations, missing data patterns, "
            "per-feature discriminative power (AUC), and paired-vs-TO shifts."
        ),
        script_path="scripts/analysis/build_observation_O01_master_distribution.py",
        data_sources=[{
            "name": "all_region_rows.tsv.gz",
            "path": str(MASTER_DATASET_PATH),
            "rows": n_rows,
            "cols": n_cols,
        }],
        row_count=n_rows,
        col_count=n_cols,
        extra={"figures_count": 12, "data_files_count": 4},
    )
    print(f"  -> round_context.json")

    print("\n" + "=" * 70)
    print(f"O1 complete. All outputs in: {OUT}")
    print("=" * 70)


if __name__ == "__main__":
    main()
