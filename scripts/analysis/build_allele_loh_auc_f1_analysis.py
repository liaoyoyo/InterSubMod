#!/usr/bin/env python3
"""Step 3.3 (Enhanced): LOH 內 Allele AUC + F1 影響分析

C.1: 2×2 AUC matrix — (LOH/non-LOH) × (HP/Allele) for TP/FP discrimination
C.2: Allele-based filtering F1 impact (scenarios A1-A5, per sample)
C.3: Cross 7-sample direction consistency + confound check

Wave 1 established:
  - HP PERMANOVA valid rate in LOH: 5-6% (unusable)
  - Allele PERMANOVA valid rate in LOH: 50-59% (10x HP)
  - Allele sig rate in LOH: 28-31% (20x HP)
  - HP AUC in LOH: 0.237-0.461 (below random)
  - Need: Allele AUC in LOH → confirm discrimination power

O11-O13 lessons:
  - Weak signals (AUC 0.55-0.65) must check confound
  - Per-sample consistency required (7/7 same direction)
  - Partial correlation with NumReads, Coverage_Multiple

Data: Master dataset (748K rows), 7 samples × 2 modes.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure,
    compute_auc, compute_effect_size, compute_metrics,
    ensure_dir, safe_div, encode_truth_binary,
    format_p, format_ci, mannwhitney_test,
    SAMPLE_ORDER, OUTPUT_ROOT,
    COLOR_TP, COLOR_FP, MODE_ORDER,
)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as sp_stats
from datetime import datetime

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260406_allele_loh_auc")

# ClairS-TO baseline FN per sample (from VCF reports / master dataset)
# Only HCC1395 has precise FN from VCF. Others use master dataset TP count ratio.
BASELINE_FN = {
    "HCC1395": 11051,
}


def assign_loh_status(df):
    """Assign LOH status using to_loh_bed_hit for TO mode, Potential_LOH as fallback."""
    df = df.copy()
    # For TO mode: use to_loh_bed_hit (region-level LOH.bed)
    # For paired mode: use Potential_LOH (HP_Ratio based)
    def _is_true(s):
        return s.astype(str).str.lower().isin(["true", "1", "yes"])

    # Universal: use Potential_LOH (ISM HP Imbalance) as LOH indicator
    df["is_loh"] = _is_true(df["Potential_LOH"])
    return df


def get_fn(sample, tp_count):
    """Get FN count for a sample. Use known baseline or estimate."""
    if sample in BASELINE_FN:
        return BASELINE_FN[sample]
    # Rough estimate: assume FN/TP ratio similar to HCC1395 (11051/28495 ≈ 0.388)
    return int(tp_count * 0.388)


# ---------------------------------------------------------------------------
# C.1: 2×2 AUC Matrix
# ---------------------------------------------------------------------------

def c1_auc_matrix(df):
    """Compute AUC for (LOH/non-LOH) × (HP/Allele) dimensions."""
    print("\n" + "=" * 50)
    print("C.1: 2×2 AUC Matrix")
    print("=" * 50)

    df["truth_binary"] = encode_truth_binary(df)
    df_valid = df[df["truth_binary"].notna()].copy()

    # Features to test
    features = {
        "HPMergedDelta": "HP effect size (Delta)",
        "AlleleDelta": "Allele effect size (Delta)",
        "LabelHPPermanovaP": "HP PERMANOVA p-value (inverted)",
        "LabelAllelePermanovaP": "Allele PERMANOVA p-value (inverted)",
    }

    results = []

    for mode in MODE_ORDER:
        dm = df_valid[df_valid["mode"] == mode]
        for loh_status in [True, False]:
            dl = dm[dm["is_loh"] == loh_status]
            loh_label = "LOH" if loh_status else "Non-LOH"

            for feat, desc in features.items():
                # Overall AUC
                vals = dl[feat].values if feat in dl.columns else np.full(len(dl), np.nan)
                y = dl["truth_binary"].values

                # For p-values, invert (lower p = more significant = higher score for TP)
                if "PermanovaP" in feat:
                    score = 1.0 - pd.to_numeric(dl[feat], errors="coerce").values
                else:
                    score = pd.to_numeric(dl[feat], errors="coerce").values

                auc_overall = compute_auc(y, score)

                # Use PERMANOVA Valid columns for valid rate (not isfinite on Delta)
                valid_col_map = {
                    "HPMergedDelta": "LabelHPPermanovaValid",
                    "AlleleDelta": "LabelAllelePermanovaValid",
                    "LabelHPPermanovaP": "LabelHPPermanovaValid",
                    "LabelAllelePermanovaP": "LabelAllelePermanovaValid",
                }
                valid_col = valid_col_map.get(feat)
                if valid_col and valid_col in dl.columns:
                    n_valid = int(dl[valid_col].astype(str).str.lower().isin(
                        ["true", "1", "yes"]).sum())
                else:
                    n_valid = int(np.isfinite(score).sum())

                row = {
                    "mode": mode,
                    "loh_status": loh_label,
                    "feature": feat,
                    "description": desc,
                    "auc": round(auc_overall, 4) if not np.isnan(auc_overall) else np.nan,
                    "n_total": len(dl),
                    "n_valid": n_valid,
                    "valid_rate_pct": round(n_valid / len(dl) * 100, 1) if len(dl) > 0 else 0,
                }
                results.append(row)

                if feat in ["AlleleDelta", "HPMergedDelta"]:
                    print(f"  [{mode.upper()}] [{loh_label}] {feat}: AUC={auc_overall:.4f}, "
                          f"valid={n_valid:,}/{len(dl):,} ({n_valid/len(dl)*100:.1f}%)")

        print()

    # Per-sample AUC for key features
    print("Per-sample AUC (AlleleDelta, TO mode):")
    sample_results = []
    for sample in SAMPLE_ORDER:
        for loh_status in [True, False]:
            ds = df_valid[(df_valid["mode"] == "to") &
                          (df_valid["sample"] == sample) &
                          (df_valid["is_loh"] == loh_status)]
            loh_label = "LOH" if loh_status else "Non-LOH"

            for feat in ["AlleleDelta", "HPMergedDelta"]:
                score = pd.to_numeric(ds[feat], errors="coerce").values
                y = ds["truth_binary"].values
                auc_val = compute_auc(y, score)

                # Use PERMANOVA Valid columns for valid rate
                per_sample_valid_map = {
                    "HPMergedDelta": "LabelHPPermanovaValid",
                    "AlleleDelta": "LabelAllelePermanovaValid",
                }
                vc = per_sample_valid_map.get(feat)
                if vc and vc in ds.columns:
                    n_valid = int(ds[vc].astype(str).str.lower().isin(
                        ["true", "1", "yes"]).sum())
                else:
                    n_valid = int(np.isfinite(score).sum())

                sample_results.append({
                    "sample": sample,
                    "loh_status": loh_label,
                    "feature": feat,
                    "auc": round(auc_val, 4) if not np.isnan(auc_val) else np.nan,
                    "n_valid": n_valid,
                })

            if loh_status:
                allele_auc = [r for r in sample_results
                              if r["sample"] == sample and r["feature"] == "AlleleDelta"
                              and r["loh_status"] == "LOH"]
                if allele_auc:
                    print(f"  {sample}: AlleleDelta AUC={allele_auc[-1]['auc']}")

    overall_df = pd.DataFrame(results)
    overall_df.to_csv(OUT_DIR / "auc_matrix_2x2.tsv", sep="\t", index=False)

    sample_df = pd.DataFrame(sample_results)
    sample_df.to_csv(OUT_DIR / "allele_loh_statistics.tsv", sep="\t", index=False)

    # Direction consistency check
    print("\n[Direction consistency] AlleleDelta AUC in LOH (TO), per sample:")
    loh_allele_aucs = sample_df[(sample_df["loh_status"] == "LOH") &
                                 (sample_df["feature"] == "AlleleDelta")]
    above_05 = (loh_allele_aucs["auc"] > 0.5).sum()
    total_valid = loh_allele_aucs["auc"].notna().sum()
    print(f"  {above_05}/{total_valid} samples have AUC > 0.5")
    above_058 = (loh_allele_aucs["auc"] > 0.58).sum()
    print(f"  {above_058}/{total_valid} samples have AUC > 0.58 (threshold)")

    return overall_df, sample_df


def fig01_auc_heatmap(overall_df):
    """F1: 2×2 AUC heatmap."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        dm = overall_df[overall_df["mode"] == mode]

        # Pivot for heatmap
        key_features = ["HPMergedDelta", "AlleleDelta"]
        pivot = dm[dm["feature"].isin(key_features)].pivot(
            index="loh_status", columns="feature", values="auc"
        )
        pivot = pivot.reindex(index=["LOH", "Non-LOH"], columns=key_features)

        sns.heatmap(pivot, annot=True, fmt=".3f", cmap="RdYlGn",
                    center=0.5, vmin=0.3, vmax=0.7, ax=ax,
                    linewidths=1, linecolor="white")
        ax.set_title(f"{mode.upper()} — AUC Matrix\n(TP vs FP discrimination)")
        ax.set_ylabel("LOH Status")
        ax.set_xlabel("Feature Dimension")

    fig.suptitle("2×2 AUC Matrix: (LOH/Non-LOH) × (HP/Allele)", fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "01_auc_heatmap_2x2.png")


def fig02_allele_auc_per_sample(sample_df):
    """F2: AlleleDelta AUC per sample in LOH vs Non-LOH."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(14, 6))

    allele = sample_df[sample_df["feature"] == "AlleleDelta"]

    x = np.arange(len(SAMPLE_ORDER))
    width = 0.35

    loh = allele[allele["loh_status"] == "LOH"].set_index("sample").reindex(SAMPLE_ORDER)
    nonloh = allele[allele["loh_status"] == "Non-LOH"].set_index("sample").reindex(SAMPLE_ORDER)

    bars1 = ax.bar(x - width/2, loh["auc"], width, label="LOH", color="#D32F2F", alpha=0.8)
    bars2 = ax.bar(x + width/2, nonloh["auc"], width, label="Non-LOH", color="#4CAF50", alpha=0.8)

    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.3, linewidth=0.8)
    ax.axhline(0.58, color="red", linestyle=":", alpha=0.5, linewidth=0.8,
               label="AUC=0.58 threshold")

    for i, (v_loh, v_nonloh) in enumerate(zip(loh["auc"], nonloh["auc"])):
        if not np.isnan(v_loh):
            ax.text(i - width/2, v_loh + 0.01, f"{v_loh:.3f}", ha="center", fontsize=7)
        if not np.isnan(v_nonloh):
            ax.text(i + width/2, v_nonloh + 0.01, f"{v_nonloh:.3f}", ha="center", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLE_ORDER, rotation=15, ha="right")
    ax.set_ylabel("AUC (TP vs FP)")
    ax.set_title("AlleleDelta AUC per Sample (TO mode)\nLOH vs Non-LOH regions")
    ax.legend()
    ax.set_ylim(0.3, 0.8)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "02_allele_auc_per_sample.png")


def fig03_allele_delta_dist(df):
    """F3: AlleleDelta distribution in LOH, TP vs FP."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    df_to = df[df["mode"] == "to"]

    for idx, loh_status in enumerate([True, False]):
        ax = axes[idx]
        dl = df_to[df_to["is_loh"] == loh_status]
        loh_label = "LOH" if loh_status else "Non-LOH"

        ad = pd.to_numeric(dl["AlleleDelta"], errors="coerce")
        valid = dl[ad.notna()].copy()
        valid["AlleleDelta_num"] = ad[ad.notna()]

        for truth, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
            subset = valid[valid["truth_label"] == truth]["AlleleDelta_num"]
            if len(subset) > 10:
                ax.hist(subset.clip(-5, 5), bins=80, alpha=0.5, color=color,
                        label=f"{truth} (n={len(subset):,}, med={subset.median():.3f})",
                        density=True)

        ax.set_xlabel("AlleleDelta")
        ax.set_ylabel("Density")
        ax.set_title(f"{loh_label} — AlleleDelta Distribution")
        ax.legend(fontsize=8)

    fig.suptitle("AlleleDelta in LOH vs Non-LOH (TO mode, TP vs FP)", fontsize=12, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "03_allele_delta_distribution.png")


def fig04_valid_rate_comparison(overall_df):
    """F4: Allele vs HP PERMANOVA valid rate in LOH."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(10, 5))

    # Get valid rates for LOH
    to_loh = overall_df[(overall_df["mode"] == "to") & (overall_df["loh_status"] == "LOH")]

    features = to_loh["feature"].values
    valid_rates = to_loh["valid_rate_pct"].values

    colors = ["#1976D2" if "HP" in f else "#D32F2F" for f in features]
    x = np.arange(len(features))
    ax.bar(x, valid_rates, 0.5, color=colors)

    for i, (f, v) in enumerate(zip(features, valid_rates)):
        ax.text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("Valid Rate (%)")
    ax.set_title("Feature Valid Rate in LOH Regions (TO mode)\n"
                 "Blue=HP dimension, Red=Allele dimension")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "04_valid_rate_in_loh.png")


# ---------------------------------------------------------------------------
# C.2: Allele-based filtering F1 impact
# ---------------------------------------------------------------------------

ALLELE_SCENARIOS = {
    "A0_baseline": {
        "desc": "No removal (baseline)",
        "remove": lambda df: pd.Series(False, index=df.index),
    },
    "A1_loh_no_allele_sig": {
        "desc": "LOH + AlleleSig=False",
        "remove": lambda df: (
            df["is_loh"] &
            ~df["AlleleSig"].astype(str).str.lower().isin(["true", "1", "yes"])
        ),
    },
    "A2_loh_allele_p_high": {
        "desc": "LOH + AlleleP > 0.1",
        "remove": lambda df: (
            df["is_loh"] &
            (pd.to_numeric(df["AlleleP"], errors="coerce").fillna(1.0) > 0.1)
        ),
    },
    "A3_loh_allele_delta_low": {
        "desc": "LOH + AlleleDelta < median",
        "remove": lambda df: (
            df["is_loh"] &
            (pd.to_numeric(df["AlleleDelta"], errors="coerce")
             < pd.to_numeric(df.loc[df["is_loh"], "AlleleDelta"], errors="coerce").median())
        ),
    },
    "A4_any_loh_dual_no_sig": {
        "desc": "any LOH + HPMergedSig=F + AlleleSig=F",
        "remove": lambda df: (
            df["is_loh"] &
            ~df["HPMergedSig"].astype(str).str.lower().isin(["true", "1", "yes"]) &
            ~df["AlleleSig"].astype(str).str.lower().isin(["true", "1", "yes"])
        ),
    },
    "A5_allele_primary": {
        "desc": "Allele-primary: LOH用Allele, non-LOH用HP filter",
        "remove": lambda df: (
            # LOH: remove if AlleleSig=False
            (df["is_loh"] &
             ~df["AlleleSig"].astype(str).str.lower().isin(["true", "1", "yes"]))
            |
            # non-LOH: remove if HPMergedSig=False
            (~df["is_loh"] &
             ~df["HPMergedSig"].astype(str).str.lower().isin(["true", "1", "yes"]))
        ),
    },
}


def c2_f1_impact(df):
    """C.2: Allele-based filtering F1 impact, per sample."""
    print("\n" + "=" * 50)
    print("C.2: Allele-based Filtering F1 Impact")
    print("=" * 50)

    df_to = df[df["mode"] == "to"].copy()
    df_to["truth_binary"] = encode_truth_binary(df_to)
    df_valid = df_to[df_to["truth_binary"].notna()].copy()

    all_results = []
    all_removed = []

    for sample in SAMPLE_ORDER:
        ds = df_valid[df_valid["sample"] == sample].copy()
        total_tp = (ds["truth_binary"] == 1).sum()
        total_fp = (ds["truth_binary"] == 0).sum()
        fn = get_fn(sample, total_tp)

        baseline_m = compute_metrics(total_tp, total_fp, fn)
        print(f"\n[{sample}] Baseline: TP={total_tp:,}, FP={total_fp:,}, "
              f"FN≈{fn:,}, F1={baseline_m['f1']:.4f}")

        for scenario_name, scenario in ALLELE_SCENARIOS.items():
            remove_mask = scenario["remove"](ds)
            n_removed = remove_mask.sum()

            removed_tp = (remove_mask & (ds["truth_binary"] == 1)).sum()
            removed_fp = (remove_mask & (ds["truth_binary"] == 0)).sum()

            new_tp = total_tp - removed_tp
            new_fp = total_fp - removed_fp
            new_fn = fn + removed_tp

            new_m = compute_metrics(new_tp, new_fp, new_fn)
            delta_f1 = new_m["f1"] - baseline_m["f1"]
            tp_loss_pct = safe_div(removed_tp, total_tp) * 100
            fp_removal_pct = safe_div(removed_fp, total_fp) * 100

            result = {
                "sample": sample,
                "scenario": scenario_name,
                "description": scenario["desc"],
                "n_removed": n_removed,
                "removed_tp": removed_tp,
                "removed_fp": removed_fp,
                "tp_loss_pct": round(tp_loss_pct, 2),
                "fp_removal_pct": round(fp_removal_pct, 2),
                "new_tp": new_tp,
                "new_fp": new_fp,
                "new_fn": new_fn,
                "new_precision": round(new_m["precision"], 4),
                "new_recall": round(new_m["recall"], 4),
                "new_f1": round(new_m["f1"], 4),
                "baseline_f1": round(baseline_m["f1"], 4),
                "delta_f1": round(delta_f1, 4),
                "tp_loss_safe": "PASS" if tp_loss_pct <= 2.0 else "FAIL",
            }
            all_results.append(result)

            if scenario_name != "A0_baseline":
                safe_tag = "PASS" if tp_loss_pct <= 2.0 else "FAIL"
                print(f"  {scenario_name}: removed {n_removed:,} "
                      f"(TP-{removed_tp:,}, FP-{removed_fp:,}), "
                      f"TP loss {tp_loss_pct:.1f}%, "
                      f"F1={new_m['f1']:.4f}, ΔF1={delta_f1:+.4f}, [{safe_tag}]")

            # Collect removed variants
            if scenario_name != "A0_baseline" and n_removed > 0:
                removed_rows = ds[remove_mask][
                    ["Chr", "Pos", "sample", "truth_label", "is_loh",
                     "AlleleSig", "HPMergedSig", "caller_af", "Coverage_Multiple",
                     "HP_Ratio", "NumReads", "AlleleDelta", "HPMergedDelta",
                     "AlleleP", "LabelAllelePermanovaP"]
                ].copy()
                removed_rows["scenario"] = scenario_name
                all_removed.append(removed_rows)

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_DIR / "f1_scenario_results.tsv", sep="\t", index=False)

    if all_removed:
        removed_df = pd.concat(all_removed, ignore_index=True)
        removed_df.to_csv(OUT_DIR / "allele_filtered_variants.tsv.gz",
                          sep="\t", index=False, compression="gzip")
        print(f"\n[Removed variants] Saved {len(removed_df):,} rows")

    return results_df


def fig05_f1_heatmap(results_df):
    """F7: F1 impact heatmap — scenarios × samples."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(14, 7))

    non_baseline = results_df[results_df["scenario"] != "A0_baseline"]
    pivot = non_baseline.pivot(index="scenario", columns="sample", values="delta_f1")
    pivot = pivot.reindex(columns=SAMPLE_ORDER)

    sns.heatmap(pivot, annot=True, fmt=".4f", cmap="RdYlGn", center=0,
                ax=ax, linewidths=1, linecolor="white",
                vmin=-0.05, vmax=0.05)
    ax.set_title("ΔF1 per Scenario × Sample (TO mode, Allele-based filtering)")
    ax.set_ylabel("Scenario")
    ax.set_xlabel("Sample")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "07_f1_scenario_heatmap.png")


def fig06_tp_fp_tradeoff(results_df):
    """F8: TP loss vs FP removal scatter (all scenarios)."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(10, 8))

    non_baseline = results_df[results_df["scenario"] != "A0_baseline"]
    scenarios = non_baseline["scenario"].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(scenarios)))

    for i, sc in enumerate(scenarios):
        ds = non_baseline[non_baseline["scenario"] == sc]
        ax.scatter(ds["fp_removal_pct"], ds["tp_loss_pct"],
                   s=80, color=colors[i], label=sc.replace("_", " "), zorder=5)

    ax.axhline(2.0, color="red", linestyle="--", alpha=0.5, label="TP loss 2% limit")
    ax.set_xlabel("FP Removal (%)")
    ax.set_ylabel("TP Loss (%)")
    ax.set_title("TP Loss vs FP Removal Tradeoff\n(Allele-based filtering scenarios, all 7 samples)")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "08_tp_fp_tradeoff.png")


# ---------------------------------------------------------------------------
# C.3: Confound check
# ---------------------------------------------------------------------------

def c3_confound_check(df):
    """Check if AlleleDelta AUC is confounded with NumReads or Coverage."""
    print("\n" + "=" * 50)
    print("C.3: Confound Check (partial correlation)")
    print("=" * 50)

    df_to = df[(df["mode"] == "to")].copy()
    df_to["truth_binary"] = encode_truth_binary(df_to)
    df_valid = df_to[df_to["truth_binary"].notna()].copy()

    confound_results = []

    for loh_status in [True, False]:
        dl = df_valid[df_valid["is_loh"] == loh_status].copy()
        loh_label = "LOH" if loh_status else "Non-LOH"

        ad = pd.to_numeric(dl["AlleleDelta"], errors="coerce")
        nr = pd.to_numeric(dl["NumReads"], errors="coerce")
        cov = pd.to_numeric(dl["Coverage_Multiple"], errors="coerce")
        y = dl["truth_binary"]

        # Raw AlleleDelta AUC
        raw_auc = compute_auc(y.values, ad.values)

        # Partial correlation: residualize AlleleDelta on NumReads
        valid_mask = ad.notna() & nr.notna() & y.notna()
        if valid_mask.sum() > 100:
            from numpy.polynomial import polynomial as P
            ad_v = ad[valid_mask].values
            nr_v = nr[valid_mask].values
            y_v = y[valid_mask].values

            # Simple linear residualization
            slope, intercept = np.polyfit(nr_v, ad_v, 1)
            ad_resid_nr = ad_v - (slope * nr_v + intercept)
            auc_resid_nr = compute_auc(y_v, ad_resid_nr)
        else:
            auc_resid_nr = np.nan

        # Partial correlation: residualize on Coverage_Multiple
        valid_mask2 = ad.notna() & cov.notna() & y.notna()
        if valid_mask2.sum() > 100:
            ad_v2 = ad[valid_mask2].values
            cov_v2 = cov[valid_mask2].values
            y_v2 = y[valid_mask2].values

            slope2, intercept2 = np.polyfit(cov_v2, ad_v2, 1)
            ad_resid_cov = ad_v2 - (slope2 * cov_v2 + intercept2)
            auc_resid_cov = compute_auc(y_v2, ad_resid_cov)
        else:
            auc_resid_cov = np.nan

        confound_results.append({
            "loh_status": loh_label,
            "raw_auc": round(raw_auc, 4) if not np.isnan(raw_auc) else np.nan,
            "auc_resid_numreads": round(auc_resid_nr, 4) if not np.isnan(auc_resid_nr) else np.nan,
            "auc_resid_coverage": round(auc_resid_cov, 4) if not np.isnan(auc_resid_cov) else np.nan,
            "auc_drop_numreads": round(raw_auc - auc_resid_nr, 4) if not (np.isnan(raw_auc) or np.isnan(auc_resid_nr)) else np.nan,
            "auc_drop_coverage": round(raw_auc - auc_resid_cov, 4) if not (np.isnan(raw_auc) or np.isnan(auc_resid_cov)) else np.nan,
        })

        print(f"  [{loh_label}] Raw AUC: {raw_auc:.4f}, "
              f"Resid(NumReads): {auc_resid_nr:.4f}, "
              f"Resid(Coverage): {auc_resid_cov:.4f}")

    confound_df = pd.DataFrame(confound_results)
    confound_df.to_csv(OUT_DIR / "confound_check.tsv", sep="\t", index=False)
    print("\n[O11-O13 check]")
    for _, row in confound_df.iterrows():
        drop_nr = row["auc_drop_numreads"]
        drop_cov = row["auc_drop_coverage"]
        if not np.isnan(drop_nr) and abs(drop_nr) > 0.05:
            print(f"  WARNING: {row['loh_status']} AUC drops {drop_nr:.4f} "
                  f"after NumReads correction → potential confound")
        if not np.isnan(drop_cov) and abs(drop_cov) > 0.05:
            print(f"  WARNING: {row['loh_status']} AUC drops {drop_cov:.4f} "
                  f"after Coverage correction → potential confound")

    return confound_df


# ---------------------------------------------------------------------------
# Conclusions
# ---------------------------------------------------------------------------

def write_conclusions(overall_df, sample_df, f1_results, confound_df):
    """Write structured conclusions."""
    # Key numbers
    to_loh_allele = overall_df[
        (overall_df["mode"] == "to") &
        (overall_df["loh_status"] == "LOH") &
        (overall_df["feature"] == "AlleleDelta")
    ]
    allele_auc_loh = to_loh_allele["auc"].values[0] if len(to_loh_allele) > 0 else np.nan

    to_loh_hp = overall_df[
        (overall_df["mode"] == "to") &
        (overall_df["loh_status"] == "LOH") &
        (overall_df["feature"] == "HPMergedDelta")
    ]
    hp_auc_loh = to_loh_hp["auc"].values[0] if len(to_loh_hp) > 0 else np.nan

    # Per-sample consistency
    loh_allele = sample_df[
        (sample_df["loh_status"] == "LOH") &
        (sample_df["feature"] == "AlleleDelta")
    ]
    above_05 = (loh_allele["auc"] > 0.5).sum()
    above_058 = (loh_allele["auc"] > 0.58).sum()
    total_valid = loh_allele["auc"].notna().sum()

    # Confound
    loh_confound = confound_df[confound_df["loh_status"] == "LOH"].iloc[0] if len(confound_df) > 0 else {}

    # F1 best scenario (TP loss ≤ 2%)
    safe = f1_results[(f1_results["tp_loss_safe"] == "PASS") &
                       (f1_results["scenario"] != "A0_baseline")]
    if len(safe) > 0:
        best_scenario = safe.loc[safe.groupby("scenario")["delta_f1"].mean().idxmax()]
        best_sc_name = safe.groupby("scenario")["delta_f1"].mean().idxmax()
        best_delta_f1 = safe.groupby("scenario")["delta_f1"].mean().max()
    else:
        best_sc_name = "None (all exceed TP loss limit)"
        best_delta_f1 = np.nan

    content = f"""<!--
建立時間: {datetime.now().strftime('%Y-%m-%d %H:%M')}
目標: LOH 內 Allele AUC + F1 影響分析 — 限制/條件/邏輯/結論紀錄
處理範圍: 7 samples × 2 modes, LOH/Non-LOH × HP/Allele
關聯檔案:
  - scripts/analysis/build_allele_loh_auc_f1_analysis.py
  - observation_workspaces/20260408_permanova_loh_audit/summary_statistics.tsv
-->

# LOH 內 Allele AUC + F1 影響分析 — 結論紀錄

## 操作條件

- **數據版本**: Master dataset 2026-03-27 HP fix rebuild (748,391 rows × 116 cols)
- **LongPhase**: Baseline
- **篩選條件**: 7 samples × 2 modes (paired/to), all chromosomes
- **LOH 定義**: ISM Potential_LOH (HP_Ratio < 0.1 or > 0.9)
- **安全約束**: TP loss ≤ 2%
- **限制**:
  - FN count 僅 HCC1395 有精確值 (11,051)，其他 samples 使用 FN/TP ratio 估計
  - AlleleDelta 是 ISM 計算的 effect size，可能受 AF confound（O12 教訓）
  - Partial correlation 使用線性 residualization（可能不足以消除非線性 confound）

## 過程邏輯

### C.1: 2×2 AUC 矩陣

- Step 1: 分層 (LOH/Non-LOH) × (HP/Allele) × (mode)
- Step 2: 每層計算 AlleleDelta 和 HPMergedDelta 的 TP/FP AUC
- Step 3: 逐 sample 計算方向一致性
- Step 4: 與 Wave 1 PERMANOVA valid rate 交叉對照

**結果**:
- TO LOH AlleleDelta AUC: **{allele_auc_loh:.4f}** {'(> 0.58 threshold → 有區分力)' if allele_auc_loh > 0.58 else '(< 0.58 threshold → 區分力不足)' if not np.isnan(allele_auc_loh) else '(N/A)'}
- TO LOH HPMergedDelta AUC: **{hp_auc_loh:.4f}** (對照)
- 方向一致性: {above_05}/{total_valid} samples AUC > 0.5, {above_058}/{total_valid} samples AUC > 0.58

### C.2: Allele 篩選 F1 影響

- 5 情境 (A1-A5) × 7 samples = 35 組模擬
- 每組計算 TP_removed, FP_removed, ΔF1
- 安全約束: TP loss ≤ 2%
- 輸出每個被移除 variant 的完整特徵

**最佳安全情境**: {best_sc_name}
**平均 ΔF1**: {f"{best_delta_f1:+.4f}" if not np.isnan(best_delta_f1) else "N/A"}

### C.3: Confound 檢查

- AlleleDelta residualized on NumReads → AUC 變化
- AlleleDelta residualized on Coverage → AUC 變化
- O11-O13 教訓: AUC drop > 0.05 → potential confound

**結果**:
- LOH raw AUC: {loh_confound.get('raw_auc', 'N/A')}
- LOH AUC after NumReads correction: {loh_confound.get('auc_resid_numreads', 'N/A')}
- LOH AUC after Coverage correction: {loh_confound.get('auc_resid_coverage', 'N/A')}
- AUC drop (NumReads): {loh_confound.get('auc_drop_numreads', 'N/A')}
- AUC drop (Coverage): {loh_confound.get('auc_drop_coverage', 'N/A')}

## 初步結論

- C1: Allele 維度在 LOH 內的 TP/FP 區分力 {'確認有效 (AUC > 0.58)' if allele_auc_loh > 0.58 else '未達閾值 (AUC ≤ 0.58)' if not np.isnan(allele_auc_loh) else '無法判定'}
  - 置信度: {'高' if above_058 >= 5 else '中' if above_058 >= 3 else '低'}
  - 證據: {above_058}/{total_valid} samples AUC > 0.58

- C2: Allele {'優於' if allele_auc_loh > hp_auc_loh else '劣於' if not np.isnan(allele_auc_loh) else '無法比較'} HP 在 LOH 內 (AUC {allele_auc_loh:.4f} vs {hp_auc_loh:.4f})
  - 置信度: 高
  - 證據: 同 dataset 直接比較

- C3: Allele-based 篩選的 F1 影響需依 TP loss 安全約束判定
  - 見 f1_scenario_results.tsv 逐 sample 細節

- C4: Confound 檢查 {"通過" if abs(float(loh_confound.get("auc_drop_numreads", 0.1) or 0.1)) < 0.05 else "需注意"}
  - 見 confound_check.tsv

## 待確認事項

1. 如果 AlleleDelta AUC 在 0.55-0.65 → 必須用 AF-bin 交叉驗證排除 confound (L3 check)
2. Allele-primary mode (A5) 是最激進的策略，需確認其 TP loss 在所有 samples 都安全
3. 需與 Script B 的象限篩選結果交叉比較，找出最優混合策略

## 後續影響

- 如果 C1 確認有效 → LOH 區域改用 Allele-only scoring（ISM C++ 端修改）
- 如果 C2 Allele >> HP → ISM 在 LOH 偵測時自動切換維度
- F1 場景分析結果輸入 Wave 3 的策略設計
- Confound 檢查結果決定是否需要 L3 AF-bin 細分驗證
"""
    path = OUT_DIR / "allele_loh_auc_conclusions.md"
    path.write_text(content, encoding="utf-8")
    print(f"[conclusions] Written: {path.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Step 3.3 (Enhanced): LOH 內 Allele AUC + F1 影響分析")
    print("=" * 70)

    df = load_master_dataset()
    df = assign_loh_status(df)

    print(f"\nTotal rows: {len(df):,}")
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        loh_n = dm["is_loh"].sum()
        print(f"  {mode.upper()}: {len(dm):,} rows, LOH={loh_n:,} ({loh_n/len(dm)*100:.1f}%)")

    # C.1: AUC Matrix
    overall_df, sample_df = c1_auc_matrix(df)

    print("\nGenerating C.1 figures...")
    fig01_auc_heatmap(overall_df)
    fig02_allele_auc_per_sample(sample_df)
    fig03_allele_delta_dist(df)
    fig04_valid_rate_comparison(overall_df)

    # C.2: F1 Impact
    f1_results = c2_f1_impact(df)

    print("\nGenerating C.2 figures...")
    fig05_f1_heatmap(f1_results)
    fig06_tp_fp_tradeoff(f1_results)

    # C.3: Confound Check
    confound_df = c3_confound_check(df)

    # Conclusions
    write_conclusions(overall_df, sample_df, f1_results, confound_df)

    print("\n" + "=" * 70)
    print(f"DONE. Output: {OUT_DIR}")
    print("Figures: 6, Tables: 5, Conclusions: 1")
    print("=" * 70)


if __name__ == "__main__":
    main()
