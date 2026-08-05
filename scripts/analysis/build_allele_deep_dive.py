#!/usr/bin/env python3
"""M4: AlleleDelta deep dive + multi-feature combination (Wave 3 Script F).

Investigates whether multi-feature combinations can improve discrimination
in LOH regions beyond the single AlleleDelta AUC=0.556.

Output: {OUTPUT_ROOT}/20260406_allele_deep_dive/
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

sys.path.insert(0, str(Path(__file__).parent))
from observation_common import (
    OUTPUT_ROOT, SAMPLE_ORDER,
    load_master_dataset, setup_plot_style, save_figure, ensure_dir,
    compute_auc, encode_truth_binary, compute_effect_size,
)

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260406_allele_deep_dive")


def require_explicit_truth_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Require real binary truth; VerificationClass is never a truth proxy."""
    if "truth_label" not in df.columns:
        raise ValueError(
            "Allele deep dive requires an explicit truth_label column with TP/FP values; "
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


# ---------------------------------------------------------------------------
# F1: Identify hopeless vs separable subgroups in LOH
# ---------------------------------------------------------------------------

def f1_hopeless_vs_separable(df: pd.DataFrame) -> pd.DataFrame:
    """Identify subgroups where AlleleDelta has zero discrimination."""
    print("\n" + "=" * 60)
    print("F1: Hopeless vs Separable Subgroup Identification")
    print("=" * 60)

    dv = df.copy()
    dv["truth_binary"] = encode_truth_binary(dv)
    dv = dv[dv["truth_binary"].notna()].copy()
    loh_to = dv[(dv["mode"] == "to") & (dv["is_loh"])].copy()

    ad = pd.to_numeric(loh_to["AlleleDelta"], errors="coerce")

    # Define hopeless: AlleleDelta very close to 0 (effect size indistinguishable)
    # Use absolute value <= 0.01 as "near zero"
    loh_to["ad_group"] = "separable"
    loh_to.loc[ad.abs() <= 0.01, "ad_group"] = "hopeless"

    results = []
    for group in ["hopeless", "separable"]:
        sub = loh_to[loh_to["ad_group"] == group]
        y = sub["truth_binary"].values
        n_tp = int(y.sum())
        n_fp = int(len(y) - y.sum())

        for feat in ["AlleleDelta", "CramersV", "PairwiseMedianDist",
                      "caller_af", "NumReads"]:
            if feat not in sub.columns:
                continue
            score = pd.to_numeric(sub[feat], errors="coerce").values
            auc = compute_auc(y, score)
            results.append({
                "group": group,
                "feature": feat,
                "auc": round(auc, 4) if not np.isnan(auc) else np.nan,
                "n_total": len(sub),
                "n_tp": n_tp,
                "n_fp": n_fp,
            })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "f1_hopeless_vs_separable.tsv", sep="\t", index=False)

    print(f"\n  Hopeless (|AlleleDelta| <= 0.01): "
          f"{(loh_to['ad_group'] == 'hopeless').sum():,}")
    print(f"  Separable (|AlleleDelta| > 0.01): "
          f"{(loh_to['ad_group'] == 'separable').sum():,}")
    print(f"\n  AUC comparison:")
    for _, r in results_df.iterrows():
        print(f"    [{r['group']}] {r['feature']}: AUC={r['auc']}")

    # Also try percentile-based definition
    ad_q25 = ad.quantile(0.25)
    print(f"\n  Alternative: AlleleDelta Q25={ad_q25:.4f}")
    below_q25 = loh_to[ad <= ad_q25]
    above_q25 = loh_to[ad > ad_q25]
    for label, sub in [("below_Q25", below_q25), ("above_Q25", above_q25)]:
        y = sub["truth_binary"].values
        score = pd.to_numeric(sub["AlleleDelta"], errors="coerce").values
        auc = compute_auc(y, score)
        print(f"    AlleleDelta AUC ({label}): {auc:.4f}, N={len(sub):,}")

    return results_df


# ---------------------------------------------------------------------------
# F3: Simple combination models
# ---------------------------------------------------------------------------

def f3_combination_models(df: pd.DataFrame) -> pd.DataFrame:
    """Test simple multi-feature combinations in LOH TO regions."""
    print("\n" + "=" * 60)
    print("F3: Multi-Feature Combination Models")
    print("=" * 60)

    dv = df.copy()
    dv["truth_binary"] = encode_truth_binary(dv)
    dv = dv[dv["truth_binary"].notna()].copy()
    loh_to = dv[(dv["mode"] == "to") & (dv["is_loh"])].copy()

    y = loh_to["truth_binary"].values

    # Normalize features to [0, 1] range
    features = {}
    for feat in ["AlleleDelta", "CramersV", "PairwiseMedianDist"]:
        vals = pd.to_numeric(loh_to[feat], errors="coerce").values
        vmin, vmax = np.nanmin(vals), np.nanmax(vals)
        if vmax > vmin:
            features[feat] = (vals - vmin) / (vmax - vmin)
        else:
            features[feat] = np.zeros_like(vals)

    # Grid search for linear combination weights
    results = []
    best_auc = 0
    best_weights = (0, 0, 0)

    weights_range = np.arange(0, 1.1, 0.2)  # Coarse grid: 0, 0.2, 0.4, 0.6, 0.8, 1.0

    for w1 in weights_range:
        for w2 in weights_range:
            for w3 in weights_range:
                if w1 + w2 + w3 == 0:
                    continue
                # Normalize weights
                wsum = w1 + w2 + w3
                w1n, w2n, w3n = w1 / wsum, w2 / wsum, w3 / wsum

                combo = (w1n * features["AlleleDelta"] +
                         w2n * features["CramersV"] +
                         w3n * features["PairwiseMedianDist"])
                auc = compute_auc(y, combo)

                results.append({
                    "w_AlleleDelta": round(w1n, 2),
                    "w_CramersV": round(w2n, 2),
                    "w_PairwiseMedianDist": round(w3n, 2),
                    "auc": round(auc, 4) if not np.isnan(auc) else np.nan,
                })

                if not np.isnan(auc) and auc > best_auc:
                    best_auc = auc
                    best_weights = (w1n, w2n, w3n)

    results_df = pd.DataFrame(results).dropna(subset=["auc"])
    results_df.to_csv(OUT_DIR / "f3_combo_grid_search.tsv", sep="\t", index=False)

    print(f"\n  Grid search: {len(results_df)} combinations tested")
    print(f"  Best AUC: {best_auc:.4f} (w={best_weights})")
    print(f"  Single feature baselines:")
    for feat in features:
        auc = compute_auc(y, features[feat])
        print(f"    {feat}: AUC={auc:.4f}")

    # Rank-based combination
    ranks = {}
    for feat in features:
        r = pd.Series(features[feat]).rank(method="average")
        ranks[feat] = r.values

    rank_combo = sum(ranks.values()) / len(ranks)
    rank_auc = compute_auc(y, rank_combo)
    print(f"  Rank-based AUC: {rank_auc:.4f}")

    # Voting (2/3 positive)
    # For each feature, above median = positive
    votes = np.zeros(len(y))
    for feat in features:
        med = np.nanmedian(features[feat])
        votes += (features[feat] > med).astype(float)
    vote_auc = compute_auc(y, votes)
    print(f"  Voting (2/3) AUC: {vote_auc:.4f}")

    # Save summary
    summary = pd.DataFrame([
        {"method": "Best linear combo", "auc": best_auc,
         "weights": str(best_weights)},
        {"method": "Rank-based", "auc": round(rank_auc, 4), "weights": "equal"},
        {"method": "Voting 2/3", "auc": round(vote_auc, 4), "weights": "equal"},
    ])
    summary.to_csv(OUT_DIR / "f3_combo_summary.tsv", sep="\t", index=False)

    return results_df


# ---------------------------------------------------------------------------
# F4: Methylation features in LOH subregions
# ---------------------------------------------------------------------------

def f4_methylation_loh_subregions(df: pd.DataFrame) -> pd.DataFrame:
    """Check if methylation features work in cnLOH vs deletion LOH."""
    print("\n" + "=" * 60)
    print("F4: Methylation Features in LOH Subregions")
    print("=" * 60)

    dv = df.copy()
    dv["truth_binary"] = encode_truth_binary(dv)
    dv = dv[dv["truth_binary"].notna()].copy()
    loh_to = dv[(dv["mode"] == "to") & (dv["is_loh"])].copy()

    cov = pd.to_numeric(loh_to["Coverage_Multiple"], errors="coerce")

    # Split by coverage: cnLOH (0.8-1.2) vs deletion LOH (<0.8)
    loh_to["loh_type"] = "unknown"
    loh_to.loc[(cov >= 0.8) & (cov <= 1.2), "loh_type"] = "cnLOH"
    loh_to.loc[cov < 0.8, "loh_type"] = "deletion_LOH"
    loh_to.loc[cov > 1.2, "loh_type"] = "gain_LOH"

    methyl_features = ["CramersV", "PairwiseMedianDist", "PairwiseMeanDist",
                        "AlleleDelta", "HPMergedDelta"]

    results = []
    for loh_type in ["cnLOH", "deletion_LOH", "gain_LOH"]:
        sub = loh_to[loh_to["loh_type"] == loh_type]
        y = sub["truth_binary"].values
        n_tp = int(y.sum())
        n_fp = int(len(y) - y.sum())

        for feat in methyl_features:
            if feat not in sub.columns:
                continue
            score = pd.to_numeric(sub[feat], errors="coerce").values
            auc = compute_auc(y, score)
            results.append({
                "loh_type": loh_type,
                "feature": feat,
                "auc": round(auc, 4) if not np.isnan(auc) else np.nan,
                "n_total": len(sub),
                "n_tp": n_tp,
                "n_fp": n_fp,
            })
            auc_str = f"{auc:.4f}" if not np.isnan(auc) else "NaN"
            print(f"  [{loh_type}] {feat}: AUC={auc_str}, "
                  f"N={len(sub):,} (TP={n_tp}, FP={n_fp})")

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "f4_methylation_loh_subregions.tsv", sep="\t", index=False)
    return results_df


# ---------------------------------------------------------------------------
# F5: Confound check for combinations
# ---------------------------------------------------------------------------

def f5_confound_check(df: pd.DataFrame) -> pd.DataFrame:
    """Residualize combo scores against NumReads and Coverage."""
    print("\n" + "=" * 60)
    print("F5: Confound Check for Combinations")
    print("=" * 60)

    dv = df.copy()
    dv["truth_binary"] = encode_truth_binary(dv)
    dv = dv[dv["truth_binary"].notna()].copy()
    loh_to = dv[(dv["mode"] == "to") & (dv["is_loh"])].copy()

    y = loh_to["truth_binary"].values

    features_to_check = ["AlleleDelta", "CramersV", "PairwiseMedianDist"]
    results = []

    for feat in features_to_check:
        score = pd.to_numeric(loh_to[feat], errors="coerce").values
        auc_raw = compute_auc(y, score)

        # Residualize against NumReads
        nr = pd.to_numeric(loh_to["NumReads"], errors="coerce").values
        mask = np.isfinite(score) & np.isfinite(nr) & np.isfinite(y)
        if mask.sum() > 10:
            slope, intercept = np.polyfit(nr[mask], score[mask], 1)
            resid_nr = score - (slope * nr + intercept)
            auc_resid_nr = compute_auc(y, resid_nr)
        else:
            auc_resid_nr = np.nan

        # Residualize against Coverage
        cov = pd.to_numeric(loh_to["Coverage_Multiple"], errors="coerce").values
        mask2 = np.isfinite(score) & np.isfinite(cov) & np.isfinite(y)
        if mask2.sum() > 10:
            slope2, intercept2 = np.polyfit(cov[mask2], score[mask2], 1)
            resid_cov = score - (slope2 * cov + intercept2)
            auc_resid_cov = compute_auc(y, resid_cov)
        else:
            auc_resid_cov = np.nan

        drop_nr = abs(auc_raw - auc_resid_nr) if not np.isnan(auc_resid_nr) else np.nan
        drop_cov = abs(auc_raw - auc_resid_cov) if not np.isnan(auc_resid_cov) else np.nan

        results.append({
            "feature": feat,
            "auc_raw": round(auc_raw, 4),
            "auc_resid_numreads": round(auc_resid_nr, 4) if not np.isnan(auc_resid_nr) else np.nan,
            "auc_resid_coverage": round(auc_resid_cov, 4) if not np.isnan(auc_resid_cov) else np.nan,
            "drop_numreads": round(drop_nr, 4) if not np.isnan(drop_nr) else np.nan,
            "drop_coverage": round(drop_cov, 4) if not np.isnan(drop_cov) else np.nan,
            "confound_pass": (drop_nr < 0.05 if not np.isnan(drop_nr) else True) and
                             (drop_cov < 0.05 if not np.isnan(drop_cov) else True),
        })

        nr_str = f"{auc_resid_nr:.4f}" if not np.isnan(auc_resid_nr) else "NaN"
        cov_str = f"{auc_resid_cov:.4f}" if not np.isnan(auc_resid_cov) else "NaN"
        print(f"  {feat}: raw={auc_raw:.4f}, resid_NR={nr_str}, resid_Cov={cov_str}")

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "f5_confound_check.tsv", sep="\t", index=False)
    return results_df


# ---------------------------------------------------------------------------
# Per-sample consistency
# ---------------------------------------------------------------------------

def per_sample_consistency(df: pd.DataFrame) -> pd.DataFrame:
    """Check 7/7 sample consistency for key features in LOH TO."""
    print("\n" + "=" * 60)
    print("Per-Sample Direction Consistency (LOH TO)")
    print("=" * 60)

    dv = df.copy()
    dv["truth_binary"] = encode_truth_binary(dv)
    dv = dv[dv["truth_binary"].notna()].copy()

    results = []
    for sample in SAMPLE_ORDER:
        ds = dv[(dv["mode"] == "to") & (dv["sample"] == sample) & (dv["is_loh"])]
        y = ds["truth_binary"].values

        for feat in ["AlleleDelta", "CramersV", "PairwiseMedianDist",
                      "caller_af", "NumReads"]:
            if feat not in ds.columns:
                continue
            score = pd.to_numeric(ds[feat], errors="coerce").values
            auc = compute_auc(y, score)
            results.append({
                "sample": sample,
                "feature": feat,
                "auc": round(auc, 4) if not np.isnan(auc) else np.nan,
                "n": len(ds),
            })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "per_sample_consistency.tsv", sep="\t", index=False)

    # Print consistency
    for feat in ["AlleleDelta", "CramersV", "PairwiseMedianDist"]:
        feat_data = results_df[results_df["feature"] == feat]
        n_above = (feat_data["auc"] > 0.5).sum()
        n_total = len(feat_data.dropna(subset=["auc"]))
        print(f"  {feat}: {n_above}/{n_total} > 0.5, "
              f"median={feat_data['auc'].median():.4f}")

    return results_df


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig01_hopeless_separable(results_df: pd.DataFrame):
    setup_plot_style()
    pivot = results_df.pivot_table(index="feature", columns="group", values="auc")
    fig, ax = plt.subplots(figsize=(8, 5))
    pivot.plot(kind="bar", ax=ax, color={"hopeless": "#EF5350", "separable": "#66BB6A"})
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.axhline(0.58, color="green", linestyle=":", alpha=0.5, label="0.58 threshold")
    ax.set_ylabel("AUC")
    ax.set_title("Feature AUC: Hopeless vs Separable LOH Subgroups (TO)")
    ax.legend()
    ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha="right")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "01_hopeless_vs_separable.png")


def fig02_combo_heatmap(combo_df: pd.DataFrame):
    setup_plot_style()
    # Show best combinations as heatmap
    top20 = combo_df.nlargest(20, "auc")
    fig, ax = plt.subplots(figsize=(8, 6))

    data = top20[["w_AlleleDelta", "w_CramersV", "w_PairwiseMedianDist", "auc"]].copy()
    labels = [f"w=({r['w_AlleleDelta']:.1f},{r['w_CramersV']:.1f},{r['w_PairwiseMedianDist']:.1f})"
              for _, r in data.iterrows()]

    ax.barh(range(len(top20)), top20["auc"].values, color="#42A5F5")
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(0.58, color="green", linestyle=":", alpha=0.5)
    ax.set_xlabel("AUC")
    ax.set_title("Top-20 Linear Combinations (LOH TO)")
    ax.invert_yaxis()
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "02_combo_top20.png")


def fig03_per_sample_bar(consistency_df: pd.DataFrame):
    setup_plot_style()
    # Per-sample AlleleDelta AUC in LOH
    ad_data = consistency_df[consistency_df["feature"] == "AlleleDelta"]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(SAMPLE_ORDER))
    aucs = []
    for sample in SAMPLE_ORDER:
        row = ad_data[ad_data["sample"] == sample]
        aucs.append(row["auc"].values[0] if len(row) > 0 else np.nan)

    colors = ["#4CAF50" if a > 0.58 else ("#FFC107" if a > 0.5 else "#F44336")
              for a in aucs]
    ax.bar(x, aucs, color=colors)
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.axhline(0.58, color="green", linestyle=":", alpha=0.5, label="0.58 threshold")
    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLE_ORDER, rotation=15, ha="right")
    ax.set_ylabel("AUC")
    ax.set_title("AlleleDelta AUC per Sample (LOH TO)")
    ax.legend()
    for i, a in enumerate(aucs):
        if not np.isnan(a):
            ax.text(i, a + 0.005, f"{a:.3f}", ha="center", fontsize=8)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "03_per_sample_allele_auc.png")


def fig04_loh_subregion(subregion_df: pd.DataFrame):
    setup_plot_style()
    pivot = subregion_df.pivot_table(index="feature", columns="loh_type", values="auc")
    fig, ax = plt.subplots(figsize=(10, 5))
    pivot.plot(kind="bar", ax=ax)
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_ylabel("AUC")
    ax.set_title("Feature AUC by LOH Subtype (cnLOH vs Deletion LOH)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha="right")
    ax.legend(title="LOH Type")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "04_loh_subregion_auc.png")


def fig05_confound(confound_df: pd.DataFrame):
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(confound_df))
    width = 0.25

    ax.bar(x - width, confound_df["auc_raw"], width, label="Raw", color="#90A4AE")
    ax.bar(x, confound_df["auc_resid_numreads"], width,
           label="Resid(NumReads)", color="#42A5F5")
    ax.bar(x + width, confound_df["auc_resid_coverage"], width,
           label="Resid(Coverage)", color="#66BB6A")

    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(confound_df["feature"], rotation=15, ha="right")
    ax.set_ylabel("AUC")
    ax.set_title("Confound Check: Raw vs Residualized AUC (LOH TO)")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "05_confound_check.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Wave 3 M4: AlleleDelta Deep Dive + Multi-Feature Combination")
    print("=" * 70)

    df = require_explicit_truth_labels(load_master_dataset())
    df = classify_loh(df)
    write_truth_schema_provenance(df)

    if "mode" not in df.columns:
        df["mode"] = df["sample_mode"] if "sample_mode" in df.columns else "unknown"
    print("Truth source: truth_label (EXPLICIT_BINARY_TRUTH)")

    # F1: Hopeless vs separable
    hopeless_df = f1_hopeless_vs_separable(df)

    # F3: Combination models
    combo_df = f3_combination_models(df)

    # F4: Methylation in LOH subregions
    subregion_df = f4_methylation_loh_subregions(df)

    # F5: Confound check
    confound_df = f5_confound_check(df)

    # Per-sample consistency
    consistency_df = per_sample_consistency(df)

    # Figures
    print("\n--- Generating figures ---")
    fig01_hopeless_separable(hopeless_df)
    fig02_combo_heatmap(combo_df)
    fig03_per_sample_bar(consistency_df)
    fig04_loh_subregion(subregion_df)
    fig05_confound(confound_df)

    # Final judgment
    best_combo_auc = combo_df["auc"].max()
    print(f"\n{'=' * 60}")
    print(f"JUDGMENT:")
    if best_combo_auc > 0.60:
        print(f"  Best combo AUC = {best_combo_auc:.4f} > 0.60 → FEASIBLE")
        print(f"  Recommend Step 3.4 Allele scoring weight design")
    elif best_combo_auc > 0.58:
        print(f"  Best combo AUC = {best_combo_auc:.4f} in [0.58, 0.60] → MARGINAL")
        print(f"  Need more sample validation before proceeding")
    else:
        print(f"  Best combo AUC = {best_combo_auc:.4f} < 0.58 → NOT FEASIBLE")
        print(f"  LOH methylation/Allele dimension definitively closed")
    print(f"{'=' * 60}")

    print(f"\nAll outputs saved to: {OUT_DIR}")
    print("M4 complete.")


if __name__ == "__main__":
    main()
