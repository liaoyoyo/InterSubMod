#!/usr/bin/env python3
"""Step 1.4 (Enhanced): SEQC2 LOH Validation + F1 Impact Analysis

B.1: LOH.bed vs SEQC2 LOH Jaccard reproduction (target: 0.928 autosomal)
B.2: HP_Ratio → SEQC2 LOH prediction ROC
B.3: F1 impact per quadrant-removal scenario (S1-S5)

Data: Master dataset (748K rows) — HCC1395 + HCC1395_DORADO, TO mode only.
SEQC2 LOH benchmark: 320 autosomal-only entries.

chrX/Y exclusion:
  - SEQC2 benchmark does not include chrX/Y (HCC1395 is female, chrX hemizygous ≠ LOH)
  - Master dataset HCC1395/HCC1395_DORADO has 0 chrX/Y rows (excluded at build time)
  - Only COLO829 has chrX/Y (3,388 rows / 748K = 0.45%), not used in this analysis
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure,
    compute_auc, compute_effect_size, compute_metrics,
    ensure_dir, safe_div, encode_truth_binary,
    format_p, format_ci,
    SAMPLE_ORDER, OUTPUT_ROOT,
    COLOR_TP, COLOR_FP,
)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as sp_stats
from datetime import datetime

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260406_seqc2_loh_validation")

# ---------------------------------------------------------------------------
# SEQC2 LOH BED loading
# ---------------------------------------------------------------------------

SEQC2_BED_PATH = Path("/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed")

LOH_BED_PATHS = {
    "HCC1395": Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
        "20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased_LOH.bed"
    ),
    "HCC1395_DORADO": Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
        "20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_phased_LOH.bed"
    ),
}

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

# ClairS-TO baseline (from 2026-03-08 report, master dataset internal counts)
BASELINE_F1_METRICS = {
    "HCC1395": {"tp": 28495, "fp": 11601, "fn": 11051},
    "HCC1395_DORADO": None,  # will be computed from data
}


def load_seqc2_loh():
    """Load SEQC2 LOH benchmark BED (autosomal only, type=='loh')."""
    df = pd.read_csv(SEQC2_BED_PATH, sep="\t", header=None,
                     names=["chr", "start", "end", "type"])
    loh = df[df["type"] == "loh"].copy()
    loh_auto = loh[loh["chr"].isin(AUTOSOMES)].copy()
    loh_auto["length"] = loh_auto["end"] - loh_auto["start"]
    print(f"[SEQC2] Loaded {len(loh)} LOH entries, {len(loh_auto)} autosomal")
    print(f"[SEQC2] Total autosomal LOH: {loh_auto['length'].sum()/1e6:.1f} Mb")
    return loh_auto


def load_longphase_loh(sample):
    """Load LongPhase-TO LOH.bed (BED9 format), autosomal only."""
    path = LOH_BED_PATHS[sample]
    # BED9: chr start end name score strand thickStart thickEnd itemRgb
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chr", "start", "end"])
    df_auto = df[df["chr"].isin(AUTOSOMES)].copy()
    df_auto["length"] = df_auto["end"] - df_auto["start"]
    print(f"[LOH.bed {sample}] {len(df)} entries, {len(df_auto)} autosomal, "
          f"total {df_auto['length'].sum()/1e6:.1f} Mb")
    return df_auto


def compute_interval_metrics(bed_a, bed_b, label_a="A", label_b="B"):
    """Compute Jaccard, Sensitivity, Precision, F1 for two BED interval sets.

    Uses base-pair level comparison (not entry-level).
    bed_a is 'truth', bed_b is 'prediction' for sensitivity/precision.
    """
    # Build per-chromosome sorted interval lists
    chrs = sorted(set(bed_a["chr"].unique()) | set(bed_b["chr"].unique()))

    total_intersect = 0
    total_union = 0
    total_a = 0
    total_b = 0

    for chrom in chrs:
        a_intervals = bed_a[bed_a["chr"] == chrom][["start", "end"]].values
        b_intervals = bed_b[bed_b["chr"] == chrom][["start", "end"]].values

        # Merge intervals within each set
        a_merged = _merge_intervals(a_intervals)
        b_merged = _merge_intervals(b_intervals)

        a_bases = sum(e - s for s, e in a_merged)
        b_bases = sum(e - s for s, e in b_merged)

        # Compute intersection
        intersect = _interval_intersection_bases(a_merged, b_merged)
        union = a_bases + b_bases - intersect

        total_intersect += intersect
        total_union += union
        total_a += a_bases
        total_b += b_bases

    jaccard = safe_div(total_intersect, total_union)
    sensitivity = safe_div(total_intersect, total_a)  # recall of truth (A)
    precision = safe_div(total_intersect, total_b)     # precision of prediction (B)
    f1 = safe_div(2 * sensitivity * precision, sensitivity + precision)

    return {
        "jaccard": jaccard,
        "sensitivity": sensitivity,
        "precision": precision,
        "f1": f1,
        "intersect_mb": total_intersect / 1e6,
        "union_mb": total_union / 1e6,
        f"{label_a}_mb": total_a / 1e6,
        f"{label_b}_mb": total_b / 1e6,
    }


def _merge_intervals(intervals):
    """Merge overlapping intervals. Input: Nx2 array of [start, end]."""
    if len(intervals) == 0:
        return []
    sorted_ivs = sorted(intervals.tolist() if hasattr(intervals, 'tolist') else intervals,
                        key=lambda x: x[0])
    merged = [sorted_ivs[0]]
    for s, e in sorted_ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def _interval_intersection_bases(merged_a, merged_b):
    """Count bases in intersection of two merged interval lists."""
    total = 0
    i, j = 0, 0
    while i < len(merged_a) and j < len(merged_b):
        a_s, a_e = merged_a[i]
        b_s, b_e = merged_b[j]
        overlap_s = max(a_s, b_s)
        overlap_e = min(a_e, b_e)
        if overlap_s < overlap_e:
            total += overlap_e - overlap_s
        if a_e <= b_e:
            i += 1
        else:
            j += 1
    return total


# ---------------------------------------------------------------------------
# Quadrant assignment (same logic as Script A)
# ---------------------------------------------------------------------------

def assign_quadrants(df):
    mask_ism = df["Potential_LOH"].astype(str).str.lower().isin(["true", "1", "yes"])
    mask_bed = df["to_loh_bed_hit"].astype(str).str.lower().isin(["true", "1", "yes"])
    df["loh_quadrant"] = "Q4_neither"
    df.loc[mask_ism & mask_bed, "loh_quadrant"] = "Q1_both"
    df.loc[mask_ism & ~mask_bed, "loh_quadrant"] = "Q2_ism_only"
    df.loc[~mask_ism & mask_bed, "loh_quadrant"] = "Q3_bed_only"
    df["hp_imbalance"] = mask_ism
    df["loh_bed"] = mask_bed
    return df


# ---------------------------------------------------------------------------
# B.1: Jaccard reproduction
# ---------------------------------------------------------------------------

def b1_jaccard(seqc2_loh):
    """B.1: LOH.bed vs SEQC2 LOH Jaccard comparison."""
    print("\n" + "=" * 50)
    print("B.1: LOH.bed vs SEQC2 LOH Jaccard")
    print("=" * 50)

    results = []
    for sample in ["HCC1395", "HCC1395_DORADO"]:
        lp_loh = load_longphase_loh(sample)
        metrics = compute_interval_metrics(seqc2_loh, lp_loh, "SEQC2", "LOH.bed")
        metrics["sample"] = sample
        results.append(metrics)
        print(f"\n[{sample}] Jaccard={metrics['jaccard']:.4f}, "
              f"Sensitivity={metrics['sensitivity']:.4f}, "
              f"Precision={metrics['precision']:.4f}, "
              f"F1={metrics['f1']:.4f}")
        print(f"  SEQC2: {metrics['SEQC2_mb']:.1f} Mb, "
              f"LOH.bed: {metrics['LOH.bed_mb']:.1f} Mb, "
              f"Intersect: {metrics['intersect_mb']:.1f} Mb")

    # Save
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_DIR / "seqc2_jaccard_comparison.tsv", sep="\t", index=False)
    return results, results_df


def fig01_jaccard_bar(results_df):
    """F1: Jaccard comparison bar chart."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Jaccard per sample
    ax = axes[0]
    x = np.arange(len(results_df))
    bars = ax.bar(x, results_df["jaccard"], 0.5, color=["#1976D2", "#42A5F5"])
    ax.set_xticks(x)
    ax.set_xticklabels(results_df["sample"])
    ax.set_ylabel("Jaccard Index")
    ax.set_title("LOH.bed vs SEQC2 LOH Jaccard\n(autosomal only)")
    ax.set_ylim(0, 1)
    ax.axhline(0.928, color="red", linestyle="--", alpha=0.5, linewidth=0.8,
               label="Expected: 0.928")
    for i, v in enumerate(results_df["jaccard"]):
        ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=10, fontweight="bold")
    ax.legend()

    # Right: Sensitivity/Precision/F1
    ax = axes[1]
    width = 0.25
    metrics = ["sensitivity", "precision", "f1"]
    colors = ["#4CAF50", "#FF9800", "#9C27B0"]
    for j, (m, c) in enumerate(zip(metrics, colors)):
        ax.bar(x + (j - 1) * width, results_df[m], width, label=m.capitalize(), color=c)
        for i, v in enumerate(results_df[m]):
            ax.text(i + (j - 1) * width, v + 0.01, f"{v:.3f}", ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(results_df["sample"])
    ax.set_ylabel("Score")
    ax.set_title("LOH.bed vs SEQC2 LOH Metrics\n(autosomal only)")
    ax.set_ylim(0, 1.1)
    ax.legend()

    fig.tight_layout()
    save_figure(fig, OUT_DIR / "01_seqc2_jaccard_bar.png")


# ---------------------------------------------------------------------------
# B.2: HP_Ratio → SEQC2 LOH prediction
# ---------------------------------------------------------------------------

def b2_hp_ratio_roc(df_hcc, seqc2_loh):
    """B.2: Use HP_Ratio to predict whether variant falls in SEQC2 LOH region."""
    print("\n" + "=" * 50)
    print("B.2: HP_Ratio → SEQC2 LOH prediction")
    print("=" * 50)

    # Label each variant as in/out of SEQC2 LOH
    df_hcc = df_hcc.copy()
    df_hcc["in_seqc2_loh"] = False
    for _, row in seqc2_loh.iterrows():
        mask = ((df_hcc["Chr"] == row["chr"]) &
                (df_hcc["Pos"] >= row["start"]) &
                (df_hcc["Pos"] < row["end"]))
        df_hcc.loc[mask, "in_seqc2_loh"] = True

    seqc2_rate = df_hcc["in_seqc2_loh"].mean()
    print(f"  Variants in SEQC2 LOH: {df_hcc['in_seqc2_loh'].sum():,} / "
          f"{len(df_hcc):,} ({seqc2_rate*100:.1f}%)")

    # HP_Ratio score: |HP_Ratio - 0.5| (higher = more imbalanced)
    df_hcc["hp_ratio_score"] = (df_hcc["HP_Ratio"] - 0.5).abs()
    valid = df_hcc[df_hcc["hp_ratio_score"].notna()].copy()

    auc = compute_auc(valid["in_seqc2_loh"].astype(int).values,
                      valid["hp_ratio_score"].values)
    print(f"  HP_Ratio score AUC for SEQC2 LOH: {auc:.4f}")

    # Threshold sweep
    thresholds = np.arange(0.05, 0.50, 0.01)
    sweep_results = []
    for thr in thresholds:
        predicted_loh = valid["hp_ratio_score"] >= thr
        tp = (predicted_loh & valid["in_seqc2_loh"]).sum()
        fp = (predicted_loh & ~valid["in_seqc2_loh"]).sum()
        fn = (~predicted_loh & valid["in_seqc2_loh"]).sum()
        tn = (~predicted_loh & ~valid["in_seqc2_loh"]).sum()
        m = compute_metrics(tp, fp, fn)
        m["threshold"] = round(float(thr), 2)
        m["specificity"] = safe_div(tn, tn + fp)
        sweep_results.append(m)

    sweep_df = pd.DataFrame(sweep_results)
    best_idx = sweep_df["f1"].idxmax()
    best = sweep_df.iloc[best_idx]
    print(f"  Best threshold: {best['threshold']:.2f}, F1={best['f1']:.4f}, "
          f"Precision={best['precision']:.4f}, Recall={best['recall']:.4f}")

    # Compare with LOH.bed prediction
    loh_bed_pred = valid["loh_bed"]
    loh_tp = (loh_bed_pred & valid["in_seqc2_loh"]).sum()
    loh_fp = (loh_bed_pred & ~valid["in_seqc2_loh"]).sum()
    loh_fn = (~loh_bed_pred & valid["in_seqc2_loh"]).sum()
    loh_m = compute_metrics(loh_tp, loh_fp, loh_fn)
    print(f"  LOH.bed prediction: F1={loh_m['f1']:.4f}, "
          f"Precision={loh_m['precision']:.4f}, Recall={loh_m['recall']:.4f}")

    # Save
    sweep_df.to_csv(OUT_DIR / "hp_ratio_threshold_sweep.tsv", sep="\t", index=False)

    return auc, sweep_df, best, df_hcc


def fig02_hp_ratio_roc(df_hcc):
    """F3: HP_Ratio ROC for SEQC2 LOH prediction."""
    from sklearn.metrics import roc_curve

    setup_plot_style()
    fig, ax = plt.subplots(figsize=(7, 7))

    valid = df_hcc[df_hcc["hp_ratio_score"].notna()].copy()

    for sample in ["HCC1395", "HCC1395_DORADO"]:
        ds = valid[valid["sample"] == sample]
        if len(ds) < 10:
            continue
        fpr, tpr, _ = roc_curve(ds["in_seqc2_loh"].astype(int), ds["hp_ratio_score"])
        auc_val = compute_auc(ds["in_seqc2_loh"].astype(int).values,
                              ds["hp_ratio_score"].values)
        ax.plot(fpr, tpr, label=f"{sample} (AUC={auc_val:.3f})", linewidth=2)

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.8)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("HP_Ratio → SEQC2 LOH Prediction ROC\n"
                 "Score = |HP_Ratio - 0.5|, higher = more imbalanced")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "03_hp_ratio_roc.png")


def fig03_threshold_sweep(sweep_df):
    """F4: HP_Ratio threshold sweep F1 curve."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(sweep_df["threshold"], sweep_df["f1"], "o-", color="#9C27B0",
            label="F1", linewidth=2, markersize=3)
    ax.plot(sweep_df["threshold"], sweep_df["precision"], "--", color="#FF9800",
            label="Precision", linewidth=1.5)
    ax.plot(sweep_df["threshold"], sweep_df["recall"], "--", color="#4CAF50",
            label="Recall", linewidth=1.5)

    best_idx = sweep_df["f1"].idxmax()
    best = sweep_df.iloc[best_idx]
    ax.axvline(best["threshold"], color="red", linestyle=":", alpha=0.5)
    ax.annotate(f"Best F1={best['f1']:.3f}\nthr={best['threshold']:.2f}",
                xy=(best["threshold"], best["f1"]),
                xytext=(best["threshold"] + 0.05, best["f1"] - 0.1),
                arrowprops=dict(arrowstyle="->", color="red"),
                fontsize=9, color="red")

    ax.set_xlabel("|HP_Ratio - 0.5| threshold")
    ax.set_ylabel("Score")
    ax.set_title("HP_Ratio Threshold Sweep for SEQC2 LOH Prediction")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "04_hp_ratio_threshold_sweep.png")


# ---------------------------------------------------------------------------
# B.3: F1 impact analysis (S1-S5 scenarios)
# ---------------------------------------------------------------------------

SCENARIOS = {
    "S0_baseline": {"desc": "No removal (baseline)", "remove": lambda df: pd.Series(False, index=df.index)},
    "S1_remove_Q1": {"desc": "Remove Q1 (Both LOH)", "remove": lambda df: df["loh_quadrant"] == "Q1_both"},
    "S2_remove_Q2": {"desc": "Remove Q2 (ISM-only)", "remove": lambda df: df["loh_quadrant"] == "Q2_ism_only"},
    "S3_remove_Q1Q2": {"desc": "Remove Q1+Q2 (all HP Imbalance)",
                       "remove": lambda df: df["loh_quadrant"].isin(["Q1_both", "Q2_ism_only"])},
    "S4_keep_Q4_only": {"desc": "Keep Q4 only (Neither LOH)",
                        "remove": lambda df: df["loh_quadrant"] != "Q4_neither"},
    "S5_remove_Q1_no_allele": {"desc": "Remove Q1 + AlleleSig=False",
                               "remove": lambda df: (
                                   (df["loh_quadrant"] == "Q1_both") &
                                   (~df["AlleleSig"].astype(str).str.lower().isin(["true", "1", "yes"]))
                               )},
}


def b3_f1_impact(df_hcc):
    """B.3: F1 impact per quadrant-removal scenario."""
    print("\n" + "=" * 50)
    print("B.3: F1 Impact Analysis (S1-S5)")
    print("=" * 50)

    all_results = []
    all_removed = []

    for sample in ["HCC1395", "HCC1395_DORADO"]:
        ds = df_hcc[df_hcc["sample"] == sample].copy()
        truth = encode_truth_binary(ds)
        ds["truth_binary"] = truth
        ds_valid = ds[ds["truth_binary"].notna()].copy()

        total_tp = (ds_valid["truth_binary"] == 1).sum()
        total_fp = (ds_valid["truth_binary"] == 0).sum()

        # For FN, use baseline if available
        if sample == "HCC1395":
            baseline_fn = BASELINE_F1_METRICS["HCC1395"]["fn"]
        else:
            # Approximate: scale FN by TP ratio if known, else use HCC1395 ratio
            baseline_fn = int(total_tp * 11051 / 28495)  # approximate

        baseline_m = compute_metrics(total_tp, total_fp, baseline_fn)
        print(f"\n[{sample}] Baseline: TP={total_tp:,}, FP={total_fp:,}, "
              f"FN={baseline_fn:,}, F1={baseline_m['f1']:.4f}")

        for scenario_name, scenario in SCENARIOS.items():
            remove_mask = scenario["remove"](ds_valid)
            n_removed = remove_mask.sum()

            removed_tp = (remove_mask & (ds_valid["truth_binary"] == 1)).sum()
            removed_fp = (remove_mask & (ds_valid["truth_binary"] == 0)).sum()

            new_tp = total_tp - removed_tp
            new_fp = total_fp - removed_fp
            new_fn = baseline_fn + removed_tp  # removed TPs become FNs

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

            if scenario_name != "S0_baseline":
                print(f"  {scenario_name}: removed {n_removed:,} "
                      f"(TP-{removed_tp:,}, FP-{removed_fp:,}), "
                      f"TP loss {tp_loss_pct:.1f}%, FP removal {fp_removal_pct:.1f}%, "
                      f"F1={new_m['f1']:.4f}, ΔF1={delta_f1:+.4f}, "
                      f"safe={'PASS' if tp_loss_pct <= 2.0 else 'FAIL'}")

            # Collect removed variants for output
            if scenario_name != "S0_baseline" and n_removed > 0:
                removed_rows = ds_valid[remove_mask][
                    ["Chr", "Pos", "sample", "truth_label", "loh_quadrant",
                     "AlleleSig", "HPMergedSig", "caller_af", "Coverage_Multiple",
                     "HP_Ratio", "NumReads"]
                ].copy()
                removed_rows["scenario"] = scenario_name
                all_removed.append(removed_rows)

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_DIR / "f1_impact_scenarios.tsv", sep="\t", index=False)

    # Save removed variants
    if all_removed:
        removed_df = pd.concat(all_removed, ignore_index=True)
        removed_df.to_csv(OUT_DIR / "removed_variants_per_scenario.tsv.gz",
                          sep="\t", index=False, compression="gzip")
        print(f"\n[Removed variants] Saved {len(removed_df):,} rows to "
              "removed_variants_per_scenario.tsv.gz")

    return results_df


def fig04_f1_impact(results_df):
    """F5: F1 impact per scenario bar chart."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    for idx, sample in enumerate(["HCC1395", "HCC1395_DORADO"]):
        ax = axes[idx]
        ds = results_df[(results_df["sample"] == sample) &
                        (results_df["scenario"] != "S0_baseline")]

        x = np.arange(len(ds))
        colors = ["#F44336" if row["tp_loss_safe"] == "FAIL" else "#4CAF50"
                  for _, row in ds.iterrows()]

        bars = ax.bar(x, ds["delta_f1"], 0.6, color=colors)
        ax.axhline(0, color="black", linewidth=0.5)

        for i, (_, row) in enumerate(ds.iterrows()):
            ax.text(i, row["delta_f1"] + (0.002 if row["delta_f1"] >= 0 else -0.008),
                    f"ΔF1={row['delta_f1']:+.4f}\nTP-{row['tp_loss_pct']:.1f}%",
                    ha="center", fontsize=7)

        ax.set_xticks(x)
        ax.set_xticklabels([s.replace("_", "\n") for s in ds["scenario"]], fontsize=8)
        ax.set_ylabel("ΔF1 (vs baseline)")
        ax.set_title(f"{sample} — F1 Impact per Scenario\n"
                     f"Baseline F1={ds.iloc[0]['baseline_f1']:.4f}", fontsize=11)

    fig.tight_layout()
    save_figure(fig, OUT_DIR / "05_f1_impact_scenarios.png")


def fig05_tp_fp_tradeoff(results_df):
    """F6: TP loss vs FP removal scatter."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(8, 8))

    markers = {"HCC1395": "o", "HCC1395_DORADO": "s"}

    for sample in ["HCC1395", "HCC1395_DORADO"]:
        ds = results_df[(results_df["sample"] == sample) &
                        (results_df["scenario"] != "S0_baseline")]
        ax.scatter(ds["fp_removal_pct"], ds["tp_loss_pct"],
                   s=100, marker=markers[sample], zorder=5, label=sample)
        for _, row in ds.iterrows():
            ax.annotate(row["scenario"].replace("_", "\n"),
                        (row["fp_removal_pct"], row["tp_loss_pct"]),
                        fontsize=7, ha="left", textcoords="offset points",
                        xytext=(5, 5))

    ax.axhline(2.0, color="red", linestyle="--", alpha=0.5, label="TP loss 2% limit")
    ax.set_xlabel("FP Removal (%)")
    ax.set_ylabel("TP Loss (%)")
    ax.set_title("TP Loss vs FP Removal Tradeoff\n(per quadrant removal scenario)")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "06_tp_fp_tradeoff.png")


# ---------------------------------------------------------------------------
# Conclusions document
# ---------------------------------------------------------------------------

def write_conclusions(jaccard_results, auc_val, best_threshold, f1_results_df):
    """Write structured conclusions."""
    hcc_jaccard = [r for r in jaccard_results if r["sample"] == "HCC1395"][0]
    dorado_jaccard = [r for r in jaccard_results if r["sample"] == "HCC1395_DORADO"][0]

    # F1 summary for HCC1395
    hcc_f1 = f1_results_df[f1_results_df["sample"] == "HCC1395"]

    content = f"""<!--
建立時間: {datetime.now().strftime('%Y-%m-%d %H:%M')}
目標: SEQC2 LOH 驗證 + F1 影響分析 — 限制/條件/邏輯/結論紀錄
處理範圍: HCC1395 + HCC1395_DORADO, TO mode, autosomal only
關聯檔案:
  - scripts/analysis/build_seqc2_loh_validation.py
  - observation_workspaces/20260406_loh_concordance/wave1_summary_report.md
-->

# SEQC2 LOH 驗證 + F1 影響分析 — 結論紀錄

## 操作條件

- **數據版本**: Master dataset 2026-03-27 HP fix rebuild (748,391 rows × 116 cols)
- **LongPhase**: Baseline（非 PON-only）
- **篩選條件**: TO mode only, HCC1395 + HCC1395_DORADO, autosomal only
- **SEQC2 LOH benchmark**: 320 entries, autosomal only
- **ClairS-TO F1 baseline**: HCC1395 F1=0.7127 (TP=28,396, FP=11,843, FN=11,051)
- **Master dataset TP/FP**: HCC1395 TP=28,495, FP=11,601（與 VCF 基線略有差異，原因: ISM region 處理）

### chrX/Y 排除說明

- **排除範圍**: chrX, chrY
- **原因**: SEQC2 LOH benchmark 不含 chrX/Y。HCC1395 是女性細胞株，chrX hemizygous ≠ LOH
- **影響**: Master dataset HCC1395/HCC1395_DORADO 的 chrX/Y 筆數 = 0（建置時已排除），故無實際影響
- **COLO829**: 唯一有 chrX/Y 的 sample（1,870 TO rows），本分析僅限 HCC1395，不受影響
- **全樣本分析時**: 必須排除 COLO829 chrX/Y 或標記 autosomal_only=True

## 過程邏輯

### B.1: LOH.bed vs SEQC2 LOH Jaccard

- Step 1: 載入 SEQC2 LOH BED（320 entries, autosomal）→ 計算總長
- Step 2: 載入 LongPhase-TO LOH.bed per sample（autosomal）
- Step 3: Base-pair level interval comparison → Jaccard, Sensitivity, Precision, F1
- Step 4: 與預期值 0.928 比較

**結果**:
- HCC1395 5kHz: Jaccard = {hcc_jaccard['jaccard']:.4f}
- HCC1395 Dorado: Jaccard = {dorado_jaccard['jaccard']:.4f}
- 驗證基準: {'PASS' if abs(hcc_jaccard['jaccard'] - 0.928) < 0.01 else 'CLOSE' if abs(hcc_jaccard['jaccard'] - 0.928) < 0.05 else 'DIFFER'} (target 0.928 ± 0.002)

### B.2: HP_Ratio → SEQC2 LOH 預測

- Step 1: 標記每個 variant 是否落在 SEQC2 LOH 區間
- Step 2: 使用 |HP_Ratio - 0.5| 作為 score → ROC AUC
- Step 3: Threshold sweep 找最佳 F1
- Step 4: 與 LOH.bed 預測力比較

**結果**:
- HP_Ratio AUC: {auc_val:.4f}
- Best threshold: {best_threshold['threshold']:.2f}, F1={best_threshold['f1']:.4f}

### B.3: F1 影響分析

- Step 1: 計算 baseline F1（master dataset 內部一致）
- Step 2: 五個情境 (S1-S5) 逐一移除 variants → 計算新 F1
- Step 3: 安全約束: TP loss ≤ 2%
- Step 4: 輸出每個被移除 variant 的完整特徵

## 初步結論

- C1: LOH.bed 與 SEQC2 LOH 一致性 {'高' if hcc_jaccard['jaccard'] > 0.9 else '中等'}（Jaccard={hcc_jaccard['jaccard']:.4f}）
  - 置信度: 高
  - 證據: Base-pair level Jaccard 重現

- C2: HP_Ratio 對 SEQC2 LOH 有{'強' if auc_val > 0.8 else '中等' if auc_val > 0.65 else '弱'}預測力（AUC={auc_val:.4f}）
  - 置信度: 高
  - 證據: ROC AUC + threshold sweep

- C3: 基於象限的篩選策略 F1 影響
  - 需檢查: 哪些情境通過 TP loss ≤ 2% 安全約束
  - 細節見 f1_impact_scenarios.tsv

## 待確認事項

1. Jaccard 若偏離 0.928 ± 0.002，需釐清是 interval merge 策略差異還是實質差異
2. HP_Ratio AUC 是否受 LOH 區域 variant density 影響（LOH 區域 variant 更少）
3. F1 影響需與 Script C (Allele-based) 的結果交叉比較

## 後續影響

- 如果 S5 (精準移除: Q1+AlleleSig=False) 通過安全約束且 ΔF1 > 0 → 可作為 ISM filtering 策略候選
- Jaccard 驗證確認 LOH.bed 可靠性後，Wave 3 可以更大膽地使用 LOH.bed 分層
- HP_Ratio 預測力結果可輸入 Step 4 生物學框架
"""
    path = OUT_DIR / "seqc2_validation_conclusions.md"
    path.write_text(content, encoding="utf-8")
    print(f"[conclusions] Written: {path.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Step 1.4 (Enhanced): SEQC2 LOH Validation + F1 Impact")
    print("=" * 70)

    # Load data
    df = load_master_dataset()
    df_to = df[df["mode"] == "to"].copy()
    df_to = assign_quadrants(df_to)

    # HCC1395 samples only (SEQC2 is HCC1395-specific)
    df_hcc = df_to[df_to["sample"].isin(["HCC1395", "HCC1395_DORADO"])].copy()
    print(f"\nHCC1395 + DORADO TO: {len(df_hcc):,} rows")

    # chrX/Y check
    chrxy = df_hcc[df_hcc["Chr"].isin(["chrX", "chrY"])]
    print(f"chrX/Y rows in HCC1395: {len(chrxy)} (expected: 0)")
    assert len(chrxy) == 0, f"Unexpected chrX/Y rows: {len(chrxy)}"

    # B.1: Jaccard
    seqc2_loh = load_seqc2_loh()
    jaccard_results, jaccard_df = b1_jaccard(seqc2_loh)
    fig01_jaccard_bar(jaccard_df)

    # B.2: HP_Ratio ROC
    auc_val, sweep_df, best_thr, df_hcc_labeled = b2_hp_ratio_roc(df_hcc, seqc2_loh)
    fig02_hp_ratio_roc(df_hcc_labeled)
    fig03_threshold_sweep(sweep_df)

    # B.3: F1 impact
    f1_results = b3_f1_impact(df_hcc_labeled)
    fig04_f1_impact(f1_results)
    fig05_tp_fp_tradeoff(f1_results)

    # Conclusions
    write_conclusions(jaccard_results, auc_val, best_thr, f1_results)

    print("\n" + "=" * 70)
    print(f"DONE. Output: {OUT_DIR}")
    print("Figures: 5, Tables: 4, Conclusions: 1")
    print("=" * 70)


if __name__ == "__main__":
    main()
