#!/usr/bin/env python3
"""M2: Non-LOH region TP/FP discrimination analysis (Wave 3 Script D).

Quantifies TP/FP discrimination power in non-LOH (Q4) regions vs LOH (Q1+Q2),
focusing on TO mode. Paired mode is observational only (too few FP for judgment).

Output: {OUTPUT_ROOT}/20260406_non_loh_discrimination/
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

sys.path.insert(0, str(Path(__file__).parent))
from observation_common import (
    OUTPUT_ROOT, SAMPLE_ORDER, MODE_ORDER,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO,
    load_master_dataset, setup_plot_style, save_figure, add_caption,
    compute_auc, encode_truth_binary, ensure_dir,
)

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260406_non_loh_discrimination")

# All features to evaluate
FEATURES = [
    "CramersV", "PairwiseMedianDist", "PairwiseMeanDist",
    "HPMergedDelta", "AlleleDelta",
    "LabelHPPermanovaF", "LabelAllelePermanovaF",
    "ClusterPermanovaF",
    "caller_af", "NumReads", "NumCpGs",
    "Coverage_Multiple",
    "HeuristicScore", "Quality_Score",
    "HP1FamilyN", "HP2FamilyN",
    "HPFineF", "HPFineNGroups",
    "GlobalP",
    "UnassignedAffinity",
    "Stability",
    "NHP3", "NHP0",
]


def require_explicit_truth_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Require real binary truth; VerificationClass is annotation only."""
    if "truth_label" not in df.columns:
        raise ValueError(
            "Non-LOH discrimination requires an explicit truth_label column with TP/FP values; "
            "VerificationClass cannot be used to manufacture truth."
        )
    truth = df["truth_label"].astype("string")
    invalid_mask = truth.isna() | ~truth.isin(["TP", "FP"])
    if invalid_mask.any():
        invalid = truth[invalid_mask].fillna("<MISSING>").value_counts().to_dict()
        raise ValueError(f"truth_label contains missing or non-binary values: {invalid}")
    validated = df.copy()
    validated["truth_label"] = truth
    validated["truth_selection_field"] = "truth_label"
    validated["truth_schema_status"] = "EXPLICIT_BINARY_TRUTH"
    return validated


def write_truth_schema_provenance(df: pd.DataFrame) -> None:
    counts = df["truth_label"].value_counts()
    pd.DataFrame(
        [
            {
                "truth_selection_field": "truth_label",
                "truth_schema_status": "EXPLICIT_BINARY_TRUTH",
                "verification_class_role": "annotation_only",
                "row_count": len(df),
                "tp_count": int(counts.get("TP", 0)),
                "fp_count": int(counts.get("FP", 0)),
            }
        ]
    ).to_csv(OUT_DIR / "truth_schema_provenance.tsv", sep="\t", index=False)


def classify_loh(df: pd.DataFrame) -> pd.DataFrame:
    """Add is_loh column based on Potential_LOH."""
    df = df.copy()
    df["is_loh"] = df["Potential_LOH"].astype(str).str.lower().isin(["true", "1", "yes"])
    return df


def d1_auc_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """D1: Compute AUC for all features, LOH vs non-LOH, TO mode."""
    print("\n" + "=" * 60)
    print("D1: Non-LOH vs LOH AUC Matrix")
    print("=" * 60)

    df["truth_binary"] = encode_truth_binary(df)
    dv = df[df["truth_binary"].notna()].copy()

    results = []
    for mode in MODE_ORDER:
        dm = dv[dv["mode"] == mode]
        for loh_status in [True, False]:
            dl = dm[dm["is_loh"] == loh_status]
            loh_label = "LOH" if loh_status else "Non-LOH"

            y = dl["truth_binary"].values
            tp_count = int(y.sum())
            fp_count = int(len(y) - y.sum())

            for feat in FEATURES:
                if feat not in dl.columns:
                    continue
                score = pd.to_numeric(dl[feat], errors="coerce").values
                auc = compute_auc(y, score)

                results.append({
                    "mode": mode,
                    "loh_status": loh_label,
                    "feature": feat,
                    "auc": round(auc, 4) if not np.isnan(auc) else np.nan,
                    "n_total": len(dl),
                    "n_tp": tp_count,
                    "n_fp": fp_count,
                    "direction": "TP-favored" if auc > 0.52 else (
                        "FP-favored" if auc < 0.48 else "Neutral"),
                })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "d1_auc_loh_vs_nonloh.tsv", sep="\t", index=False)
    print(f"  Saved {len(results_df)} rows to d1_auc_loh_vs_nonloh.tsv")

    # Print summary for TO
    to_results = results_df[results_df["mode"] == "to"]
    for loh_label in ["LOH", "Non-LOH"]:
        sub = to_results[to_results["loh_status"] == loh_label]
        above_058 = sub[sub["auc"] > 0.58]
        print(f"\n  TO {loh_label}: {len(above_058)}/{len(sub)} features AUC > 0.58")
        if len(above_058) > 0:
            for _, r in above_058.sort_values("auc", ascending=False).iterrows():
                print(f"    {r['feature']}: AUC={r['auc']:.4f} ({r['direction']})")

    return results_df


def d2_solvable_regions(df: pd.DataFrame) -> pd.DataFrame:
    """D2: Quantify solvable vs unsolvable regions."""
    print("\n" + "=" * 60)
    print("D2: Solvable Region Quantification")
    print("=" * 60)

    df["truth_binary"] = encode_truth_binary(df)
    dv = df[df["truth_binary"].notna()].copy()
    dv_to = dv[dv["mode"] == "to"]

    # Per-region best AUC is not meaningful at single-site level
    # Instead, count TP/FP ratios in LOH vs non-LOH
    results = []
    for sample in SAMPLE_ORDER:
        ds = dv_to[dv_to["sample"] == sample]
        for is_loh in [True, False]:
            dl = ds[ds["is_loh"] == is_loh]
            loh_label = "LOH" if is_loh else "Non-LOH"
            n_total = len(dl)
            y = dl["truth_binary"].values
            n_tp = int(y.sum())
            n_fp = n_total - n_tp
            tp_rate = n_tp / n_total if n_total > 0 else 0
            fp_rate = n_fp / n_total if n_total > 0 else 0

            results.append({
                "sample": sample,
                "loh_status": loh_label,
                "n_total": n_total,
                "n_tp": n_tp,
                "n_fp": n_fp,
                "tp_rate": round(tp_rate, 4),
                "fp_rate": round(fp_rate, 4),
            })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "d2_region_tp_fp_counts.tsv", sep="\t", index=False)
    print(f"  Saved to d2_region_tp_fp_counts.tsv")

    # Print summary
    for loh_label in ["LOH", "Non-LOH"]:
        sub = results_df[results_df["loh_status"] == loh_label]
        avg_fp_rate = sub["fp_rate"].mean()
        total_n = sub["n_total"].sum()
        print(f"  {loh_label}: avg FP rate = {avg_fp_rate:.4f}, total N = {total_n:,}")

    return results_df


def d3_per_sample_nonloh(df: pd.DataFrame) -> pd.DataFrame:
    """D3: Per-sample non-LOH AUC for top features (TO mode only)."""
    print("\n" + "=" * 60)
    print("D3: Per-Sample Non-LOH AUC")
    print("=" * 60)

    df["truth_binary"] = encode_truth_binary(df)
    dv = df[df["truth_binary"].notna()].copy()
    dv_to = dv[dv["mode"] == "to"]

    top_features = [
        "CramersV", "PairwiseMedianDist", "AlleleDelta", "HPMergedDelta",
        "caller_af", "NumReads", "HeuristicScore", "Quality_Score",
        "ClusterPermanovaF", "LabelAllelePermanovaF",
    ]

    results = []
    for sample in SAMPLE_ORDER:
        ds = dv_to[dv_to["sample"] == sample]
        for is_loh in [True, False]:
            dl = ds[ds["is_loh"] == is_loh]
            loh_label = "LOH" if is_loh else "Non-LOH"
            y = dl["truth_binary"].values

            for feat in top_features:
                if feat not in dl.columns:
                    continue
                score = pd.to_numeric(dl[feat], errors="coerce").values
                auc = compute_auc(y, score)
                results.append({
                    "sample": sample,
                    "loh_status": loh_label,
                    "feature": feat,
                    "auc": round(auc, 4) if not np.isnan(auc) else np.nan,
                    "n": len(dl),
                })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "d3_per_sample_nonloh_auc.tsv", sep="\t", index=False)

    # Check 7/7 direction consistency for non-LOH
    print("\n  Non-LOH 7/7 direction consistency (TO):")
    nonloh = results_df[results_df["loh_status"] == "Non-LOH"]
    for feat in top_features:
        feat_data = nonloh[nonloh["feature"] == feat]
        n_above = (feat_data["auc"] > 0.5).sum()
        n_below = (feat_data["auc"] < 0.5).sum()
        consistency = max(n_above, n_below)
        print(f"    {feat}: {consistency}/7 consistent "
              f"({'TP' if n_above > n_below else 'FP'}-favored, "
              f"median AUC={feat_data['auc'].median():.4f})")

    return results_df


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig01_auc_comparison_heatmap(auc_df: pd.DataFrame):
    """F1: AUC comparison heatmap — LOH vs non-LOH × features (TO mode)."""
    setup_plot_style()

    to_data = auc_df[auc_df["mode"] == "to"].copy()

    # Pivot: rows=features, columns=LOH status
    pivot = to_data.pivot_table(
        index="feature", columns="loh_status", values="auc", aggfunc="first"
    )
    pivot = pivot.reindex(columns=["LOH", "Non-LOH"])
    pivot = pivot.dropna(how="all")
    pivot = pivot.sort_values("Non-LOH", ascending=False)

    fig, ax = plt.subplots(figsize=(8, max(8, len(pivot) * 0.35)))
    sns.heatmap(
        pivot, annot=True, fmt=".3f", cmap="RdYlGn", center=0.5,
        vmin=0.35, vmax=0.65, linewidths=0.5, ax=ax,
        cbar_kws={"label": "AUC (TP=1)"},
    )
    ax.set_title("TO Mode: Feature AUC — LOH vs Non-LOH\n"
                 "Green=TP-favored (>0.5), Red=FP-favored (<0.5)")
    ax.set_ylabel("Feature")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "01_auc_loh_vs_nonloh_heatmap.png")


def fig02_nonloh_top10_bar(auc_df: pd.DataFrame):
    """F2: Non-LOH top-10 features ranked bar (TO mode)."""
    setup_plot_style()

    to_nonloh = auc_df[
        (auc_df["mode"] == "to") & (auc_df["loh_status"] == "Non-LOH")
    ].copy()
    to_nonloh["auc_abs"] = (to_nonloh["auc"] - 0.5).abs()
    top10 = to_nonloh.nlargest(10, "auc_abs")

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#4CAF50" if a > 0.5 else "#F44336" for a in top10["auc"]]
    bars = ax.barh(range(len(top10)), top10["auc"].values, color=colors)
    ax.set_yticks(range(len(top10)))
    ax.set_yticklabels(top10["feature"].values)
    ax.axvline(0.5, color="black", linestyle="--", alpha=0.5)
    ax.axvline(0.58, color="green", linestyle=":", alpha=0.5, label="0.58 threshold")
    ax.set_xlabel("AUC")
    ax.set_title("Top-10 Features by Discrimination Power in Non-LOH Regions (TO)")
    ax.legend()

    for i, (_, r) in enumerate(top10.iterrows()):
        ax.text(r["auc"] + 0.005, i, f"{r['auc']:.3f}", va="center", fontsize=9)

    ax.invert_yaxis()
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "02_nonloh_top10_features.png")


def fig03_per_sample_heatmap(sample_df: pd.DataFrame):
    """F3: Per-sample non-LOH AUC heatmap."""
    setup_plot_style()

    nonloh = sample_df[sample_df["loh_status"] == "Non-LOH"].copy()
    pivot = nonloh.pivot_table(
        index="feature", columns="sample", values="auc", aggfunc="first"
    )
    pivot = pivot.reindex(columns=SAMPLE_ORDER)

    # Sort by median AUC
    pivot["median"] = pivot.median(axis=1)
    pivot = pivot.sort_values("median", ascending=False)
    pivot = pivot.drop(columns="median")

    fig, ax = plt.subplots(figsize=(12, max(6, len(pivot) * 0.4)))
    sns.heatmap(
        pivot, annot=True, fmt=".3f", cmap="RdYlGn", center=0.5,
        vmin=0.35, vmax=0.65, linewidths=0.5, ax=ax,
    )
    ax.set_title("Per-Sample Feature AUC in Non-LOH Regions (TO mode)")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "03_per_sample_nonloh_auc.png")


def fig04_tp_fp_ratio_comparison(region_df: pd.DataFrame):
    """F4: TP/FP ratio comparison LOH vs non-LOH per sample."""
    setup_plot_style()

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(SAMPLE_ORDER))
    width = 0.35

    loh_fp = []
    nonloh_fp = []
    for sample in SAMPLE_ORDER:
        for is_loh, store in [(True, loh_fp), (False, nonloh_fp)]:
            label = "LOH" if is_loh else "Non-LOH"
            row = region_df[
                (region_df["sample"] == sample) & (region_df["loh_status"] == label)
            ]
            store.append(row["fp_rate"].values[0] if len(row) > 0 else 0)

    ax.bar(x - width / 2, loh_fp, width, label="LOH FP rate", color="#EF5350")
    ax.bar(x + width / 2, nonloh_fp, width, label="Non-LOH FP rate", color="#42A5F5")

    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLE_ORDER, rotation=15, ha="right")
    ax.set_ylabel("FP Rate")
    ax.set_title("FP Rate: LOH vs Non-LOH (TO mode)\n"
                 "Lower FP rate in LOH confirms TP-enrichment (J9)")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "04_fp_rate_loh_vs_nonloh.png")


def fig05_auc_delta_distribution(auc_df: pd.DataFrame):
    """F5: AUC delta (LOH minus Non-LOH) distribution."""
    setup_plot_style()

    to_data = auc_df[auc_df["mode"] == "to"].copy()
    loh_auc = to_data[to_data["loh_status"] == "LOH"].set_index("feature")["auc"]
    nonloh_auc = to_data[to_data["loh_status"] == "Non-LOH"].set_index("feature")["auc"]

    common = loh_auc.index.intersection(nonloh_auc.index)
    delta = (loh_auc.loc[common] - nonloh_auc.loc[common]).dropna()
    delta = delta.sort_values()

    fig, ax = plt.subplots(figsize=(10, max(5, len(delta) * 0.3)))
    colors = ["#4CAF50" if d > 0 else "#F44336" for d in delta.values]
    ax.barh(range(len(delta)), delta.values, color=colors)
    ax.set_yticks(range(len(delta)))
    ax.set_yticklabels(delta.index)
    ax.axvline(0, color="black", linestyle="-", alpha=0.5)
    ax.set_xlabel("AUC Delta (LOH - Non-LOH)")
    ax.set_title("Feature AUC: LOH minus Non-LOH (TO mode)\n"
                 "Green=feature works better in LOH, Red=better in non-LOH")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "05_auc_delta_loh_minus_nonloh.png")


def fig06_stacked_category(region_df: pd.DataFrame):
    """F6: TP vs FP stacked bar per sample, split by LOH status."""
    setup_plot_style()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for idx, loh_label in enumerate(["LOH", "Non-LOH"]):
        ax = axes[idx]
        sub = region_df[region_df["loh_status"] == loh_label]

        tp_vals = []
        fp_vals = []
        for sample in SAMPLE_ORDER:
            row = sub[sub["sample"] == sample]
            tp_vals.append(row["n_tp"].values[0] if len(row) > 0 else 0)
            fp_vals.append(row["n_fp"].values[0] if len(row) > 0 else 0)

        x = np.arange(len(SAMPLE_ORDER))
        ax.bar(x, tp_vals, label="TP", color=COLOR_TP, alpha=0.8)
        ax.bar(x, fp_vals, bottom=tp_vals, label="FP", color=COLOR_FP, alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right", fontsize=8)
        ax.set_title(f"{loh_label} Regions (TO)")
        ax.set_ylabel("Count")
        ax.legend()

    fig.suptitle("TP / FP Distribution: LOH vs Non-LOH (TO mode)")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "06_tp_fp_stacked_by_loh.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Wave 3 M2: Non-LOH Region TP/FP Discrimination Analysis")
    print("=" * 70)

    df = require_explicit_truth_labels(load_master_dataset())
    df = classify_loh(df)
    write_truth_schema_provenance(df)

    # Add mode column if missing
    if "mode" not in df.columns:
        df["mode"] = df["sample_mode"] if "sample_mode" in df.columns else "unknown"

    print("Truth source: truth_label (EXPLICIT_BINARY_TRUTH)")

    print(f"\nDataset: {len(df):,} rows")
    print(f"LOH: {df['is_loh'].sum():,} ({df['is_loh'].mean()*100:.1f}%)")
    print(f"Non-LOH: {(~df['is_loh']).sum():,} ({(~df['is_loh']).mean()*100:.1f}%)")

    # D1: AUC matrix
    auc_df = d1_auc_matrix(df)

    # D2: Region counts
    region_df = d2_solvable_regions(df)

    # D3: Per-sample
    sample_df = d3_per_sample_nonloh(df)

    # Figures
    print("\n--- Generating figures ---")
    fig01_auc_comparison_heatmap(auc_df)
    fig02_nonloh_top10_bar(auc_df)
    fig03_per_sample_heatmap(sample_df)
    fig04_tp_fp_ratio_comparison(region_df)
    fig05_auc_delta_distribution(auc_df)
    fig06_stacked_category(region_df)

    print(f"\nAll outputs saved to: {OUT_DIR}")
    print("M2 complete.")


if __name__ == "__main__":
    main()
