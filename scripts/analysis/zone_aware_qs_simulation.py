#!/usr/bin/env python3
"""
Zone-Aware QS Adjustment Simulation

Simulates the F1 impact of zone-aware QS adjustment on the TO mode.
Approach: For each region, compute QS_adjusted = QS + zone_delta.
Then sweep a QS threshold and compare F1 curves.

Methodology:
  - Original F1: all caller results (TP+FP), no ISM filter
  - QS-filtered F1: only keep variants with QS >= threshold
  - Zone-aware QS F1: adjust QS by zone, then apply same threshold

Key insight: ISM QS is used as a post-filter on caller results.
Zone-aware adjustment changes which variants survive the filter.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/" \
         "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/zone_aware_validation"
FIG_DIR = os.path.join(OUT_DIR, "figures")

SAMPLE_ORDER = ["H2009", "H1437", "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "COLO829"]

# ── Load ───────────────────────────────────────────────────────────
print("Loading...")
cols = ["sample", "mode", "truth_label", "Quality_Score",
        "core_loh_like", "to_loh_bed_hit",
        "caller_af", "HPFineNGroups", "NumReads",
        "Coverage_Multiple", "context_truth_total"]
df = pd.read_csv(MASTER, sep='\t', usecols=cols, low_memory=False)
print(f"Loaded {len(df):,} rows")

# Type conversions
for col in ["Quality_Score", "caller_af", "HPFineNGroups", "NumReads",
            "Coverage_Multiple", "context_truth_total"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
for col in ["core_loh_like", "to_loh_bed_hit"]:
    df[col] = df[col].map({"True": True, "False": False, True: True, False: False}).fillna(False)
df["is_tp"] = df["truth_label"] == "TP"
df["loh_flag"] = np.where(df["mode"] == "paired", df["core_loh_like"], df["to_loh_bed_hit"])

# AF classification
def is_intermediate_af(af):
    if pd.isna(af): return False
    return (0.1 <= af <= 0.4) or (0.6 <= af <= 0.9)

def is_extreme_af(af):
    if pd.isna(af): return False
    return (af < 0.1) or (af > 0.9)

df["af_intermediate"] = df["caller_af"].apply(is_intermediate_af)
df["af_extreme"] = df["caller_af"].apply(is_extreme_af)

# ── Zone Assignment (updated Z1 definition) ────────────────────────
# Z1: LOH + Intermediate AF + NGroups >= 2 (relaxed from >=3)
df["zone_z1"] = df["loh_flag"] & df["af_intermediate"] & (df["HPFineNGroups"] >= 2)
df["zone_z2"] = (df["HPFineNGroups"] >= 4) & (df["NumReads"] >= 80)
df["zone_z3"] = df["loh_flag"] & df["af_extreme"] & (df["HPFineNGroups"] <= 1)
df["zone_z5"] = (df["Coverage_Multiple"] > 1.5) & (df["HPFineNGroups"] <= 2) & (~df["loh_flag"])
df["zone_z4"] = (~df["zone_z1"]) & (~df["zone_z2"]) & (~df["zone_z3"]) & \
                (~df["zone_z5"]) & (df["Coverage_Multiple"].between(0.7, 1.3))

def assign_zone(row):
    if row["zone_z1"]: return "Z1"
    if row["zone_z2"]: return "Z2"
    if row["zone_z3"]: return "Z3"
    if row["zone_z5"]: return "Z5"
    if row["zone_z4"]: return "Z4"
    return "Unassigned"

df["zone"] = df.apply(assign_zone, axis=1)

# ── F1 Calculation Helper ──────────────────────────────────────────
def compute_f1(tp_count, fp_count, fn_count):
    """Compute F1 from TP, FP, FN counts."""
    if tp_count == 0:
        return 0.0
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


# ── QS Threshold Sweep ─────────────────────────────────────────────
# Zone delta values to test (TO mode only)
# Based on zone TP rates: Z1=0.965, Z2=0.891, Z3=0.608, Z4=0.694, Z5=0.667
# Strategy: boost Z1/Z2 (TP-rich), penalize Z3/Z5 (TP-poor)
delta_configs = {
    "No Adjustment": {"Z1": 0, "Z2": 0, "Z3": 0, "Z4": 0, "Z5": 0, "Unassigned": 0},
    "Conservative (+5/-5)": {"Z1": 5, "Z2": 5, "Z3": -5, "Z4": 0, "Z5": -5, "Unassigned": 0},
    "Moderate (+10/-10)": {"Z1": 10, "Z2": 10, "Z3": -10, "Z4": 0, "Z5": -10, "Unassigned": 0},
    "Aggressive (+15/-15)": {"Z1": 15, "Z2": 15, "Z3": -15, "Z4": 0, "Z5": -15, "Unassigned": 0},
    "Asymmetric (+10/-15)": {"Z1": 10, "Z2": 10, "Z3": -15, "Z4": 0, "Z5": -10, "Unassigned": 0},
}

# Focus on TO mode
to_df = df[df["mode"] == "to"].copy()
print(f"\nTO mode: {len(to_df):,} regions")

# Sweep thresholds
thresholds = list(range(0, 105, 5))

results_all = []
for sample in SAMPLE_ORDER:
    sm = to_df[to_df["sample"] == sample]
    if len(sm) == 0:
        continue

    truth_total = sm["context_truth_total"].iloc[0]
    original_tp = sm["is_tp"].sum()
    original_fp = len(sm) - original_tp
    original_fn = truth_total - original_tp
    original_f1 = compute_f1(original_tp, original_fp, original_fn)

    print(f"\n{sample}: N={len(sm):,} TP={original_tp:,} FP={original_fp:,} FN={original_fn:,} "
          f"F1(no filter)={original_f1:.4f} Truth={truth_total:,}")

    for config_name, deltas in delta_configs.items():
        # Compute adjusted QS
        qs_adj = sm["Quality_Score"].copy()
        for zone_label, delta in deltas.items():
            mask = sm["zone"] == zone_label
            qs_adj.loc[mask] += delta

        for t in thresholds:
            # Apply threshold: keep variants with QS >= t
            pass_mask = qs_adj >= t
            passed = sm[pass_mask]
            tp_pass = passed["is_tp"].sum()
            fp_pass = len(passed) - tp_pass
            fn = truth_total - tp_pass  # FN = truth - TP that survived filter
            f1 = compute_f1(tp_pass, fp_pass, fn)

            results_all.append({
                "sample": sample, "config": config_name,
                "threshold": t,
                "tp": tp_pass, "fp": fp_pass, "fn": fn,
                "f1": f1, "original_f1": original_f1,
                "delta_f1": f1 - original_f1,
                "n_passed": len(passed), "n_total": len(sm)
            })

rdf = pd.DataFrame(results_all)
rdf.to_csv(os.path.join(OUT_DIR, "zone_aware_qs_simulation.tsv"), sep='\t', index=False)

# ── Find Optimal Threshold per Config per Sample ──────────────────
print("\n" + "="*80)
print("OPTIMAL THRESHOLD ANALYSIS")
print("="*80)

optimal_results = []
for sample in SAMPLE_ORDER:
    sdf = rdf[rdf["sample"] == sample]
    if len(sdf) == 0: continue
    print(f"\n{sample}:")
    for config in delta_configs.keys():
        cdf = sdf[sdf["config"] == config]
        best = cdf.loc[cdf["f1"].idxmax()]
        print(f"  {config:<25s}: best_T={best['threshold']:3.0f} F1={best['f1']:.4f} "
              f"(delta={best['delta_f1']:+.4f}) TP={best['tp']:,} FP={best['fp']:,}")
        optimal_results.append({
            "sample": sample, "config": config,
            "best_threshold": best["threshold"],
            "best_f1": best["f1"],
            "delta_f1_vs_nofilter": best["delta_f1"],
            "tp": best["tp"], "fp": best["fp"], "fn": best["fn"]
        })

opt_df = pd.DataFrame(optimal_results)
opt_df.to_csv(os.path.join(OUT_DIR, "zone_aware_optimal_thresholds.tsv"), sep='\t', index=False)

# ── Figure: F1 vs Threshold Curves per Sample ─────────────────────
fig, axes = plt.subplots(2, 4, figsize=(24, 12))
axes = axes.flatten()

config_colors = {
    "No Adjustment": 'gray',
    "Conservative (+5/-5)": '#3498db',
    "Moderate (+10/-10)": '#2ecc71',
    "Aggressive (+15/-15)": '#e74c3c',
    "Asymmetric (+10/-15)": '#f39c12'
}

for idx, sample in enumerate(SAMPLE_ORDER):
    ax = axes[idx]
    sdf = rdf[rdf["sample"] == sample]
    if len(sdf) == 0:
        ax.set_visible(False)
        continue

    orig_f1 = sdf["original_f1"].iloc[0]

    for config, color in config_colors.items():
        cdf = sdf[sdf["config"] == config]
        ax.plot(cdf["threshold"], cdf["f1"], color=color, label=config,
               linewidth=2 if "No" not in config else 1.5,
               linestyle='--' if "No" in config else '-')

    ax.axhline(y=orig_f1, color='black', linestyle=':', alpha=0.5,
               label=f'No filter ({orig_f1:.4f})')
    ax.set_xlabel("QS Threshold")
    ax.set_ylabel("F1")
    ax.set_title(sample, fontweight='bold')
    ax.grid(alpha=0.3)
    if idx == 0:
        ax.legend(fontsize=6, loc='lower left')

# Hide unused subplot
axes[7].set_visible(False)

plt.suptitle("Zone-Aware QS Adjustment: F1 vs Threshold — TO Mode",
             fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "zone_aware_qs_f1_curves.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"\n→ Figure: {FIG_DIR}/zone_aware_qs_f1_curves.png")

# ── Summary: Best Delta F1 Across Configs ─────────────────────────
print("\n" + "="*80)
print("SUMMARY: Best F1 Improvement per Config (across samples)")
print("="*80)

for config in delta_configs.keys():
    copt = opt_df[opt_df["config"] == config]
    mean_f1 = copt["best_f1"].mean()
    mean_delta = copt["delta_f1_vs_nofilter"].mean()
    n_positive = (copt["delta_f1_vs_nofilter"] > 0).sum()
    print(f"  {config:<25s}: Mean F1={mean_f1:.4f} Mean Delta={mean_delta:+.5f} "
          f"({n_positive}/{len(copt)} positive)")

# ── Compare Zone-Aware to non-Zone at fixed threshold ─────────────
print("\n" + "="*80)
print("FIXED THRESHOLD COMPARISON (T=50, T=70)")
print("="*80)

for t in [50, 70]:
    print(f"\n--- Threshold = {t} ---")
    for sample in SAMPLE_ORDER:
        sdf = rdf[(rdf["sample"] == sample) & (rdf["threshold"] == t)]
        if len(sdf) == 0: continue
        no_adj = sdf[sdf["config"] == "No Adjustment"].iloc[0]
        print(f"  {sample:16s} No-adj F1={no_adj['f1']:.4f}", end="")
        for config in ["Moderate (+10/-10)", "Asymmetric (+10/-15)"]:
            adj = sdf[sdf["config"] == config].iloc[0]
            print(f" | {config}: {adj['f1']:.4f} (Δ={adj['f1'] - no_adj['f1']:+.5f})", end="")
        print()

print("\nDone.")
