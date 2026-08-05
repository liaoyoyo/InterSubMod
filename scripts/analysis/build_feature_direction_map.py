#!/usr/bin/env python3
"""M3: Feature direction map + correctability assessment (Wave 3 Script E).

Creates 4-way AUC matrix (Paired/TO × LOH/Non-LOH) for 31+ features,
identifies mode-flip features, evaluates AF-bin cross-validation,
and simulates corrections.

Output: {OUTPUT_ROOT}/20260406_feature_direction_map/
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

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260406_feature_direction_map")

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

AF_BINS = [(0, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 1.01)]
AF_LABELS = ["0-0.2", "0.2-0.5", "0.5-0.8", "0.8-1.0"]


def require_explicit_truth_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Require real binary truth; VerificationClass is annotation only."""
    if "truth_label" not in df.columns:
        raise ValueError(
            "Feature direction analysis requires an explicit truth_label column with TP/FP values; "
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
    df = df.copy()
    df["is_loh"] = df["Potential_LOH"].astype(str).str.lower().isin(["true", "1", "yes"])
    return df


def classify_direction(auc: float) -> str:
    if np.isnan(auc):
        return "N/A"
    if auc > 0.52:
        return "Positive"
    if auc < 0.48:
        return "Inverted"
    return "Neutral"


# ---------------------------------------------------------------------------
# E1: 4-way AUC matrix
# ---------------------------------------------------------------------------

def e1_four_way_auc(df: pd.DataFrame) -> pd.DataFrame:
    """Compute AUC for all features × (Paired/TO) × (LOH/Non-LOH)."""
    print("\n" + "=" * 60)
    print("E1: 4-Way AUC Matrix")
    print("=" * 60)

    df["truth_binary"] = encode_truth_binary(df)
    dv = df[df["truth_binary"].notna()].copy()

    results = []
    for mode in MODE_ORDER:
        dm = dv[dv["mode"] == mode]
        for is_loh in [True, False]:
            dl = dm[dm["is_loh"] == is_loh]
            loh_label = "LOH" if is_loh else "Non-LOH"
            y = dl["truth_binary"].values

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
                    "direction": classify_direction(auc),
                    "n": len(dl),
                    "condition": f"{mode.upper()}_{loh_label}",
                })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "e1_four_way_auc.tsv", sep="\t", index=False)

    # Identify mode-flip features (Paired positive, TO inverted)
    print("\nMode-flip features (Paired positive → TO inverted):")
    for feat in FEATURES:
        paired_rows = results_df[
            (results_df["feature"] == feat) & (results_df["mode"] == "paired")
        ]
        to_rows = results_df[
            (results_df["feature"] == feat) & (results_df["mode"] == "to")
        ]
        if len(paired_rows) == 0 or len(to_rows) == 0:
            continue

        paired_auc_mean = paired_rows["auc"].mean()
        to_auc_mean = to_rows["auc"].mean()

        if paired_auc_mean > 0.52 and to_auc_mean < 0.48:
            print(f"  FLIP: {feat} — Paired={paired_auc_mean:.4f}, TO={to_auc_mean:.4f}")
        elif paired_auc_mean < 0.48 and to_auc_mean > 0.52:
            print(f"  FLIP (reverse): {feat} — Paired={paired_auc_mean:.4f}, TO={to_auc_mean:.4f}")

    return results_df


# ---------------------------------------------------------------------------
# E3: AF-bin cross-validation
# ---------------------------------------------------------------------------

def e3_af_bin_validation(df: pd.DataFrame) -> pd.DataFrame:
    """Cross-validate AUC by AF bins to detect confound."""
    print("\n" + "=" * 60)
    print("E3: AF-Bin Cross-Validation")
    print("=" * 60)

    df["truth_binary"] = encode_truth_binary(df)
    dv = df[df["truth_binary"].notna()].copy()
    dv["caller_af_numeric"] = pd.to_numeric(dv["caller_af"], errors="coerce")

    # Focus on TO mode (where inversions happen)
    dv_to = dv[dv["mode"] == "to"]

    # Key features to check (especially inverted ones)
    check_features = [
        "caller_af", "HPMergedDelta", "AlleleDelta", "CramersV",
        "PairwiseMedianDist", "Quality_Score", "NumReads",
    ]

    results = []
    for feat in check_features:
        if feat not in dv_to.columns:
            continue
        for i, (lo, hi) in enumerate(AF_BINS):
            af_bin = dv_to[
                (dv_to["caller_af_numeric"] >= lo) &
                (dv_to["caller_af_numeric"] < hi)
            ]
            y = af_bin["truth_binary"].values
            score = pd.to_numeric(af_bin[feat], errors="coerce").values
            auc = compute_auc(y, score)

            results.append({
                "feature": feat,
                "af_bin": AF_LABELS[i],
                "af_lo": lo,
                "af_hi": hi,
                "auc": round(auc, 4) if not np.isnan(auc) else np.nan,
                "direction": classify_direction(auc),
                "n": len(af_bin),
            })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "e3_af_bin_auc.tsv", sep="\t", index=False)

    # Check for direction flips across AF bins
    print("\nAF-bin direction flips (TO mode):")
    for feat in check_features:
        feat_data = results_df[results_df["feature"] == feat]
        directions = feat_data["direction"].values
        if len(set(d for d in directions if d != "N/A")) > 1:
            print(f"  FLIP: {feat} — {list(zip(feat_data['af_bin'], feat_data['auc'], directions))}")
        else:
            d = directions[0] if len(directions) > 0 else "N/A"
            aucs = feat_data["auc"].dropna()
            print(f"  STABLE: {feat} — {d}, AUC range [{aucs.min():.3f}, {aucs.max():.3f}]")

    return results_df


# ---------------------------------------------------------------------------
# E4: Correction simulation
# ---------------------------------------------------------------------------

def e4_correction_simulation(df: pd.DataFrame) -> pd.DataFrame:
    """Simulate corrections for mode-flip features."""
    print("\n" + "=" * 60)
    print("E4: Correction Simulation")
    print("=" * 60)

    df["truth_binary"] = encode_truth_binary(df)
    dv = df[df["truth_binary"].notna()].copy()
    dv_to = dv[dv["mode"] == "to"]
    dv_to = dv_to.copy()

    results = []

    # Correction 1: HPMergedDelta → use AlleleDelta in LOH regions
    feat = "HPMergedDelta"
    y = dv_to["truth_binary"].values

    # Original
    score_orig = pd.to_numeric(dv_to[feat], errors="coerce").values
    auc_orig = compute_auc(y, score_orig)

    # Corrected: LOH uses AlleleDelta, non-LOH uses HPMergedDelta
    score_corrected = np.where(
        dv_to["is_loh"].values,
        pd.to_numeric(dv_to["AlleleDelta"], errors="coerce").values,
        score_orig,
    )
    auc_corrected = compute_auc(y, score_corrected)

    results.append({
        "feature": feat,
        "correction": "LOH→AlleleDelta",
        "auc_original": round(auc_orig, 4),
        "auc_corrected": round(auc_corrected, 4),
        "auc_delta": round(auc_corrected - auc_orig, 4),
    })
    print(f"  {feat}: {auc_orig:.4f} → {auc_corrected:.4f} (Δ={auc_corrected-auc_orig:+.4f})")

    # Correction 2: caller_af → abs(af - 0.5) symmetrized
    feat = "caller_af"
    score_orig = pd.to_numeric(dv_to[feat], errors="coerce").values
    auc_orig = compute_auc(y, score_orig)

    # Symmetrized: higher distance from 0.5 → more separable
    score_sym = np.abs(score_orig - 0.5)
    auc_sym = compute_auc(y, score_sym)

    # Inverted symmetrized: lower distance from 0.5 → more like heterozygous
    score_inv_sym = -score_sym
    auc_inv_sym = compute_auc(y, score_inv_sym)

    results.append({
        "feature": feat,
        "correction": "abs(af-0.5)",
        "auc_original": round(auc_orig, 4),
        "auc_corrected": round(auc_sym, 4),
        "auc_delta": round(auc_sym - auc_orig, 4),
    })
    results.append({
        "feature": feat,
        "correction": "-abs(af-0.5)",
        "auc_original": round(auc_orig, 4),
        "auc_corrected": round(auc_inv_sym, 4),
        "auc_delta": round(auc_inv_sym - auc_orig, 4),
    })
    print(f"  {feat}: {auc_orig:.4f} → abs(af-0.5)={auc_sym:.4f}, -abs(af-0.5)={auc_inv_sym:.4f}")

    # Correction 3: LOH-aware QS (remove LOH penalty, use Allele in LOH)
    if "Quality_Score" in dv_to.columns:
        qs_orig = pd.to_numeric(dv_to["Quality_Score"], errors="coerce").values
        auc_qs_orig = compute_auc(y, qs_orig)

        # Simple correction: in LOH regions, add Allele bonus
        allele_delta = pd.to_numeric(dv_to["AlleleDelta"], errors="coerce").fillna(0).values
        qs_corrected = qs_orig.copy()
        loh_mask = dv_to["is_loh"].values
        qs_corrected[loh_mask] = qs_corrected[loh_mask] + allele_delta[loh_mask] * 20
        auc_qs_corrected = compute_auc(y, qs_corrected)

        results.append({
            "feature": "Quality_Score",
            "correction": "LOH+AlleleDelta*20",
            "auc_original": round(auc_qs_orig, 4),
            "auc_corrected": round(auc_qs_corrected, 4),
            "auc_delta": round(auc_qs_corrected - auc_qs_orig, 4),
        })
        print(f"  Quality_Score: {auc_qs_orig:.4f} → {auc_qs_corrected:.4f} "
              f"(Δ={auc_qs_corrected-auc_qs_orig:+.4f})")

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "e4_correction_simulation.tsv", sep="\t", index=False)
    return results_df


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig01_four_way_heatmap(auc_df: pd.DataFrame):
    """F1: 4-way AUC heatmap."""
    setup_plot_style()

    pivot = auc_df.pivot_table(
        index="feature", columns="condition", values="auc", aggfunc="first"
    )
    col_order = ["PAIRED_LOH", "PAIRED_Non-LOH", "TO_LOH", "TO_Non-LOH"]
    pivot = pivot.reindex(columns=[c for c in col_order if c in pivot.columns])

    # Sort by TO overall AUC
    to_cols = [c for c in pivot.columns if c.startswith("TO")]
    pivot["sort_key"] = pivot[to_cols].mean(axis=1) if to_cols else 0
    pivot = pivot.sort_values("sort_key", ascending=False).drop(columns="sort_key")

    fig, ax = plt.subplots(figsize=(10, max(8, len(pivot) * 0.35)))
    sns.heatmap(
        pivot, annot=True, fmt=".3f", cmap="RdYlGn", center=0.5,
        vmin=0.35, vmax=0.65, linewidths=0.5, ax=ax,
    )
    ax.set_title("4-Way AUC: (Paired/TO) × (LOH/Non-LOH)\n"
                 "Green=TP-favored, Red=FP-favored")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "01_four_way_auc_heatmap.png")


def fig02_mode_flip_scatter(auc_df: pd.DataFrame):
    """F2: Paired AUC vs TO AUC scatter, annotating mode-flip features."""
    setup_plot_style()

    # Average across LOH states for each mode
    avg = auc_df.groupby(["feature", "mode"])["auc"].mean().reset_index()
    paired = avg[avg["mode"] == "paired"].set_index("feature")["auc"]
    to = avg[avg["mode"] == "to"].set_index("feature")["auc"]
    common = paired.index.intersection(to.index)

    fig, ax = plt.subplots(figsize=(8, 8))
    for feat in common:
        p_auc = paired[feat]
        t_auc = to[feat]
        if np.isnan(p_auc) or np.isnan(t_auc):
            continue

        # Color by flip type
        if p_auc > 0.52 and t_auc < 0.48:
            color = "#F44336"
            marker = "s"
        elif p_auc < 0.48 and t_auc > 0.52:
            color = "#FF9800"
            marker = "D"
        else:
            color = "#757575"
            marker = "o"

        ax.scatter(p_auc, t_auc, c=color, marker=marker, s=60, zorder=5)
        ax.annotate(feat, (p_auc, t_auc), fontsize=7, alpha=0.8,
                    xytext=(3, 3), textcoords="offset points")

    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.3)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.3)
    ax.plot([0.3, 0.7], [0.3, 0.7], "k--", alpha=0.2)

    # Legend
    ax.scatter([], [], c="#F44336", marker="s", label="Paired+ → TO- (mode-flip)")
    ax.scatter([], [], c="#FF9800", marker="D", label="Paired- → TO+ (reverse flip)")
    ax.scatter([], [], c="#757575", marker="o", label="Consistent direction")
    ax.legend(loc="lower right")

    ax.set_xlabel("Paired Mode AUC (avg LOH+non-LOH)")
    ax.set_ylabel("TO Mode AUC (avg LOH+non-LOH)")
    ax.set_title("Feature Direction: Paired vs TO Mode")
    ax.set_xlim(0.3, 0.7)
    ax.set_ylim(0.3, 0.7)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "02_mode_flip_scatter.png")


def fig03_af_bin_lines(af_df: pd.DataFrame):
    """F3: AF-bin AUC line chart for key features (TO mode)."""
    setup_plot_style()

    features_to_plot = af_df["feature"].unique()
    n_feat = len(features_to_plot)

    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.cm.Set2(np.linspace(0, 1, n_feat))

    for i, feat in enumerate(features_to_plot):
        feat_data = af_df[af_df["feature"] == feat].sort_values("af_lo")
        ax.plot(
            range(len(AF_LABELS)), feat_data["auc"].values,
            marker="o", label=feat, color=cmap[i], linewidth=2,
        )

    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.axhline(0.58, color="green", linestyle=":", alpha=0.5)
    ax.set_xticks(range(len(AF_LABELS)))
    ax.set_xticklabels(AF_LABELS)
    ax.set_xlabel("Allele Frequency Bin")
    ax.set_ylabel("AUC")
    ax.set_title("Feature AUC by AF Bin (TO mode)\nL3 Cross-Validation for Confound Detection")
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "03_af_bin_auc_lines.png")


def fig04_correction_comparison(corr_df: pd.DataFrame):
    """F4: Before vs after correction AUC comparison."""
    setup_plot_style()

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(corr_df))
    width = 0.35

    ax.bar(x - width / 2, corr_df["auc_original"], width,
           label="Original", color="#90A4AE")
    ax.bar(x + width / 2, corr_df["auc_corrected"], width,
           label="Corrected", color="#4CAF50")

    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    labels = [f"{r['feature']}\n({r['correction']})" for _, r in corr_df.iterrows()]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("AUC")
    ax.set_title("Feature Correction Simulation (TO mode)\n"
                 "Gray=original, Green=corrected")
    ax.legend()

    for i, (_, r) in enumerate(corr_df.iterrows()):
        delta = r["auc_delta"]
        ax.text(i + width / 2, r["auc_corrected"] + 0.005,
                f"Δ={delta:+.4f}", ha="center", fontsize=8)

    fig.tight_layout()
    save_figure(fig, OUT_DIR / "04_correction_comparison.png")


def fig05_per_sample_consistency(df: pd.DataFrame):
    """F5: Per-sample direction consistency heatmap (TO mode, non-LOH)."""
    setup_plot_style()

    df["truth_binary"] = encode_truth_binary(df)
    dv = df[df["truth_binary"].notna()].copy()
    dv_to = dv[(dv["mode"] == "to") & (~dv["is_loh"])]

    key_features = [
        "CramersV", "PairwiseMedianDist", "AlleleDelta", "HPMergedDelta",
        "caller_af", "NumReads", "HeuristicScore", "Quality_Score",
    ]

    results = []
    for sample in SAMPLE_ORDER:
        ds = dv_to[dv_to["sample"] == sample]
        y = ds["truth_binary"].values
        for feat in key_features:
            if feat not in ds.columns:
                continue
            score = pd.to_numeric(ds[feat], errors="coerce").values
            auc = compute_auc(y, score)
            results.append({"sample": sample, "feature": feat, "auc": auc})

    res_df = pd.DataFrame(results)
    pivot = res_df.pivot_table(index="feature", columns="sample", values="auc")
    pivot = pivot.reindex(columns=SAMPLE_ORDER)

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(
        pivot, annot=True, fmt=".3f", cmap="RdYlGn", center=0.5,
        vmin=0.35, vmax=0.65, linewidths=0.5, ax=ax,
    )
    ax.set_title("Per-Sample AUC in Non-LOH Regions (TO mode)\n"
                 "7/7 same direction = confirmed")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "05_per_sample_consistency.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Wave 3 M3: Feature Direction Map + Correctability Assessment")
    print("=" * 70)

    df = require_explicit_truth_labels(load_master_dataset())
    df = classify_loh(df)
    write_truth_schema_provenance(df)

    if "mode" not in df.columns:
        df["mode"] = df["sample_mode"] if "sample_mode" in df.columns else "unknown"
    print("Truth source: truth_label (EXPLICIT_BINARY_TRUTH)")

    # E1
    auc_df = e1_four_way_auc(df)

    # E3
    af_df = e3_af_bin_validation(df)

    # E4
    corr_df = e4_correction_simulation(df)

    # Figures
    print("\n--- Generating figures ---")
    fig01_four_way_heatmap(auc_df)
    fig02_mode_flip_scatter(auc_df)
    fig03_af_bin_lines(af_df)
    fig04_correction_comparison(corr_df)
    fig05_per_sample_consistency(df)

    print(f"\nAll outputs saved to: {OUT_DIR}")
    print("M3 complete.")


if __name__ == "__main__":
    main()
