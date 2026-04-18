#!/usr/bin/env python3
"""
LOH + AF Filtering Validation, Tier Cutpoint Verification, Quality Score Analysis
(2026-03-31)

Part 1: Validate AF-based TP/FP separation within LOH-like regions
Part 2: effective_hp_reads distribution histograms for tier cutpoint verification
Part 3: Quality Score component analysis and potential redesign directions
"""

import os
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")

ROUND1_WS = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260327_loh_round1_cross_sample_audit"
)
OUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260331_loh_af_tier_qs_analysis"
)
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
SAMPLE_LABELS = {
    "HCC1395": "HCC1395_5kHz", "HCC1395_DORADO": "HCC1395_DORADO",
    "COLO829": "COLO829", "H1437": "H1437", "H2009": "H2009",
    "HCC1937": "HCC1937", "HCC1954": "HCC1954",
}
TIER_ORDER = ["C0", "C", "B", "A", "A+"]
TRUTH_COLORS = {"TP": "#4daf4a", "FP": "#e41a1c"}

def assign_tier(n):
    if pd.isna(n): return "C0"
    n = int(n)
    if n == 0: return "C0"
    elif n < 10: return "C"
    elif n < 30: return "B"
    elif n < 50: return "A"
    else: return "A+"

def write_tsv(df, path, desc=""):
    df.to_csv(path, sep="\t", index=False)
    print(f"  {path.name} ({len(df)} rows) — {desc}")


# ============================================================
# Load Data
# ============================================================
print("[0] Loading master dataset ...")
df = pd.read_csv(ROUND1_WS / "all_region_rows.tsv.gz", sep="\t", low_memory=False)
df["tier"] = df["effective_hp_reads"].apply(assign_tier)
df["tier"] = pd.Categorical(df["tier"], categories=TIER_ORDER, ordered=True)
df["core_loh_like"] = df["core_loh_like"].astype(bool)
df["HPMergedSig"] = df["HPMergedSig"].map({"True": True, "False": False, True: True, False: False}).fillna(False).astype(bool)
print(f"  Loaded {len(df):,} rows")


# ============================================================
# Part 1: LOH + AF Filtering Validation
# ============================================================
print("\n[1/3] Part 1 — LOH + AF filtering validation ...")

AF_BINS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.01]
AF_LABELS = ["0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6",
             "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.95", "0.95-1.0"]

# --- P1a: AF distribution within LOH-like, by mode × truth ---
rows_p1a = []
for mode in ["paired", "to"]:
    dm = df[(df["mode"] == mode) & (df["core_loh_like"] == True)]
    for truth in ["TP", "FP"]:
        sub = dm[dm["truth_label"] == truth]
        af_vals = sub["caller_af"].dropna()
        binned = pd.cut(af_vals, bins=AF_BINS, labels=AF_LABELS, right=False)
        counts = binned.value_counts().reindex(AF_LABELS, fill_value=0)
        for af_label, count in counts.items():
            rows_p1a.append({
                "mode": mode, "truth": truth, "af_bin": af_label,
                "n": int(count), "frac": round(count / len(sub), 4) if len(sub) > 0 else 0,
            })
df_p1a = pd.DataFrame(rows_p1a)
write_tsv(df_p1a, OUT_DIR / "p1a_loh_af_distribution.tsv", "AF distribution in LOH-like by mode × truth")

# --- P1b: AF threshold sweep — TP precision/recall ---
rows_p1b = []
for mode in ["paired", "to"]:
    dm = df[(df["mode"] == mode) & (df["core_loh_like"] == True)]
    for tier in ["A", "A+", "All"]:
        if tier == "All":
            dt = dm
        else:
            dt = dm[dm["tier"] == tier]
        n_tp_total = (dt["truth_label"] == "TP").sum()
        n_fp_total = (dt["truth_label"] == "FP").sum()
        n_total = len(dt)
        for af_thresh in [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95]:
            # LOH-like AND AF < threshold → predict TP
            selected = dt[dt["caller_af"] < af_thresh]
            tp_in = (selected["truth_label"] == "TP").sum()
            fp_in = (selected["truth_label"] == "FP").sum()
            n_selected = len(selected)
            precision = tp_in / n_selected if n_selected > 0 else 0
            recall = tp_in / n_tp_total if n_tp_total > 0 else 0
            # LOH-like AND AF > threshold → predict FP
            fp_selected = dt[dt["caller_af"] >= af_thresh]
            fp_caught = (fp_selected["truth_label"] == "FP").sum()
            tp_lost = (fp_selected["truth_label"] == "TP").sum()
            fp_precision = fp_caught / len(fp_selected) if len(fp_selected) > 0 else 0
            fp_recall = fp_caught / n_fp_total if n_fp_total > 0 else 0
            rows_p1b.append({
                "mode": mode, "tier": tier, "af_threshold": af_thresh,
                "n_loh_total": n_total, "n_tp_total": n_tp_total, "n_fp_total": n_fp_total,
                "tp_predict_n": n_selected, "tp_predict_tp": tp_in, "tp_predict_fp": fp_in,
                "tp_precision": round(precision, 4), "tp_recall": round(recall, 4),
                "fp_predict_n": len(fp_selected), "fp_caught": fp_caught, "tp_lost": tp_lost,
                "fp_precision": round(fp_precision, 4), "fp_recall": round(fp_recall, 4),
            })
df_p1b = pd.DataFrame(rows_p1b)
write_tsv(df_p1b, OUT_DIR / "p1b_af_threshold_sweep.tsv", "AF threshold sweep in LOH-like")

# --- P1c: Per-sample verification ---
rows_p1c = []
for sample in SAMPLES:
    for mode in ["paired", "to"]:
        dm = df[(df["sample"] == sample) & (df["mode"] == mode) & (df["core_loh_like"] == True)]
        if len(dm) < 10:
            continue
        for af_thresh in [0.8, 0.95]:
            low_af = dm[dm["caller_af"] < af_thresh]
            high_af = dm[dm["caller_af"] >= af_thresh]
            rows_p1c.append({
                "sample": sample, "mode": mode, "af_threshold": af_thresh,
                "n_loh": len(dm),
                "low_af_n": len(low_af),
                "low_af_tp": (low_af["truth_label"] == "TP").sum(),
                "low_af_fp": (low_af["truth_label"] == "FP").sum(),
                "low_af_tp_purity": round((low_af["truth_label"] == "TP").sum() / len(low_af), 4) if len(low_af) > 0 else 0,
                "high_af_n": len(high_af),
                "high_af_tp": (high_af["truth_label"] == "TP").sum(),
                "high_af_fp": (high_af["truth_label"] == "FP").sum(),
                "high_af_fp_rate": round((high_af["truth_label"] == "FP").sum() / len(high_af), 4) if len(high_af) > 0 else 0,
            })
df_p1c = pd.DataFrame(rows_p1c)
write_tsv(df_p1c, OUT_DIR / "p1c_per_sample_af_verification.tsv", "Per-sample AF threshold verification")

# --- Fig01: AF histogram in LOH-like (paired vs TO, TP vs FP) ---
fig, axes = plt.subplots(2, 2, figsize=(16, 10))
for i, mode in enumerate(["paired", "to"]):
    for j, truth in enumerate(["TP", "FP"]):
        ax = axes[i, j]
        sub = df[(df["mode"] == mode) & (df["core_loh_like"] == True) & (df["truth_label"] == truth)]
        af_vals = sub["caller_af"].dropna()
        ax.hist(af_vals, bins=50, color=TRUTH_COLORS[truth], alpha=0.7, edgecolor="black", linewidth=0.3)
        ax.axvline(x=0.8, color="blue", linestyle="--", linewidth=1.5, label="AF=0.8")
        ax.axvline(x=0.95, color="red", linestyle="--", linewidth=1.5, label="AF=0.95")
        n = len(sub)
        med = af_vals.median() if len(af_vals) > 0 else 0
        ax.set_title(f"{mode.upper()} {truth} LOH-like (n={n:,}, median AF={med:.3f})", fontsize=11, fontweight="bold")
        ax.set_xlabel("Caller AF")
        ax.set_ylabel("Count")
        ax.legend(fontsize=8)
fig.suptitle("Part 1: AF Distribution in LOH-like Regions", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig01_loh_af_histogram.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig01_loh_af_histogram.png")

# --- Fig02: AF threshold sweep — precision/recall curves (TO) ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for i, tier in enumerate(["A", "A+", "All"]):
    ax = axes[i]
    sub = df_p1b[(df_p1b["mode"] == "to") & (df_p1b["tier"] == tier)]
    thresholds = sub["af_threshold"].values
    ax.plot(thresholds, sub["tp_precision"].values, "g-o", label="TP precision (AF<thresh)", markersize=5)
    ax.plot(thresholds, sub["tp_recall"].values, "g--s", label="TP recall (AF<thresh)", markersize=5)
    ax.plot(thresholds, sub["fp_precision"].values, "r-o", label="FP precision (AF≥thresh)", markersize=5)
    ax.plot(thresholds, sub["fp_recall"].values, "r--s", label="FP recall (AF≥thresh)", markersize=5)
    ax.axvline(x=0.8, color="blue", linestyle=":", alpha=0.5)
    ax.axvline(x=0.95, color="red", linestyle=":", alpha=0.5)
    n_total = sub["n_loh_total"].values[0] if len(sub) > 0 else 0
    ax.set_title(f"TO Tier {tier} (n={n_total:,})", fontsize=12, fontweight="bold")
    ax.set_xlabel("AF Threshold")
    ax.set_ylabel("Precision / Recall")
    ax.legend(fontsize=7, loc="center left")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
fig.suptitle("Part 1: AF Threshold Sweep in TO LOH-like — TP/FP Precision & Recall", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig02_af_threshold_sweep_to.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig02_af_threshold_sweep_to.png")


# ============================================================
# Part 2: effective_hp_reads Distribution & Tier Cutpoint Verification
# ============================================================
print("\n[2/3] Part 2 — effective_hp_reads distribution ...")

# --- P2a: effective_hp_reads statistics by mode × truth ---
rows_p2a = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for truth in ["TP", "FP"]:
        sub = dm[dm["truth_label"] == truth]
        ehr = sub["effective_hp_reads"]
        rows_p2a.append({
            "mode": mode, "truth": truth, "n": len(sub),
            "mean": round(ehr.mean(), 2), "median": round(ehr.median(), 2),
            "q10": round(ehr.quantile(0.10), 2), "q25": round(ehr.quantile(0.25), 2),
            "q75": round(ehr.quantile(0.75), 2), "q90": round(ehr.quantile(0.90), 2),
            "pct_eq0": round((ehr == 0).mean() * 100, 2),
            "pct_lt10": round((ehr < 10).mean() * 100, 2),
            "pct_lt30": round((ehr < 30).mean() * 100, 2),
            "pct_ge50": round((ehr >= 50).mean() * 100, 2),
            "pct_ge80": round((ehr >= 80).mean() * 100, 2),
        })
df_p2a = pd.DataFrame(rows_p2a)
write_tsv(df_p2a, OUT_DIR / "p2a_ehr_statistics.tsv", "effective_hp_reads statistics")

# --- P2b: Per-sample tier distribution ---
rows_p2b = []
for sample in SAMPLES:
    for mode in ["paired", "to"]:
        dm = df[(df["sample"] == sample) & (df["mode"] == mode)]
        if len(dm) == 0:
            continue
        for truth in ["TP", "FP"]:
            sub = dm[dm["truth_label"] == truth]
            n = len(sub)
            for tier in TIER_ORDER:
                tier_n = (sub["tier"] == tier).sum()
                rows_p2b.append({
                    "sample": sample, "mode": mode, "truth": truth,
                    "tier": tier, "n": int(tier_n),
                    "frac": round(tier_n / n, 4) if n > 0 else 0,
                })
df_p2b = pd.DataFrame(rows_p2b)
write_tsv(df_p2b, OUT_DIR / "p2b_per_sample_tier_dist.tsv", "Per-sample tier distribution")

# --- Fig03: effective_hp_reads histogram (TP vs FP, paired vs TO) ---
fig, axes = plt.subplots(2, 2, figsize=(16, 10))
tier_boundaries = [0, 10, 30, 50]
tier_labels_at = [5, 20, 40, 75]
tier_names = ["C0/C", "B", "A", "A+"]

for i, mode in enumerate(["paired", "to"]):
    for j, truth in enumerate(["TP", "FP"]):
        ax = axes[i, j]
        sub = df[(df["mode"] == mode) & (df["truth_label"] == truth)]
        ehr = sub["effective_hp_reads"].clip(upper=150)
        ax.hist(ehr, bins=100, color=TRUTH_COLORS[truth], alpha=0.7, edgecolor="none")
        for b in tier_boundaries:
            ax.axvline(x=b, color="black", linestyle="--", linewidth=1, alpha=0.5)
        for pos, name in zip(tier_labels_at, tier_names):
            ax.text(pos, ax.get_ylim()[1] * 0.95, name, ha="center", fontsize=9,
                    fontweight="bold", bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
        n = len(sub)
        med = sub["effective_hp_reads"].median()
        ax.set_title(f"{mode.upper()} {truth} (n={n:,}, median={med:.0f})", fontsize=11, fontweight="bold")
        ax.set_xlabel("effective_hp_reads")
        ax.set_ylabel("Count")
        ax.set_xlim(-2, 155)
fig.suptitle("Part 2: effective_hp_reads Distribution with Tier Boundaries", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig03_ehr_histogram_with_tiers.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig03_ehr_histogram_with_tiers.png")

# --- Fig04: Per-sample tier stacked bar (TO mode, TP vs FP side by side) ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
tier_colors = {"C0": "#d9d9d9", "C": "#fc8d59", "B": "#fee090", "A": "#4575b4", "A+": "#313695"}
for j, truth in enumerate(["TP", "FP"]):
    ax = axes[j]
    bottom = np.zeros(len(SAMPLES))
    for tier in TIER_ORDER:
        vals = []
        for sample in SAMPLES:
            row = df_p2b[(df_p2b["sample"] == sample) & (df_p2b["mode"] == "to") &
                         (df_p2b["truth"] == truth) & (df_p2b["tier"] == tier)]
            vals.append(row["frac"].values[0] if len(row) > 0 else 0)
        vals = np.array(vals)
        ax.bar(range(len(SAMPLES)), vals, bottom=bottom, color=tier_colors[tier], label=tier, edgecolor="white", linewidth=0.3)
        bottom += vals
    ax.set_xticks(range(len(SAMPLES)))
    ax.set_xticklabels([SAMPLE_LABELS.get(s, s) for s in SAMPLES], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Fraction")
    ax.set_title(f"TO {truth} — Tier Distribution", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(0, 1.05)
fig.suptitle("Part 2: TO Per-Sample Tier Distribution (TP vs FP)", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig04_to_per_sample_tier_stacked.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig04_to_per_sample_tier_stacked.png")

# --- Fig05: FP rate by effective_hp_reads bins (fine-grained, TO mode) ---
fig, ax = plt.subplots(figsize=(14, 5))
dm = df[df["mode"] == "to"]
bins_fine = list(range(0, 105, 5)) + [150, 200, 500]
dm_copy = dm.copy()
dm_copy["ehr_bin"] = pd.cut(dm_copy["effective_hp_reads"], bins=[-1] + bins_fine, labels=[f"{bins_fine[i-1] if i>0 else 0}-{bins_fine[i]}" for i in range(len(bins_fine))])
bin_stats = dm_copy.groupby("ehr_bin", observed=False).agg(
    n=("truth_label", "size"),
    fp_n=("truth_label", lambda s: (s == "FP").sum()),
).reset_index()
bin_stats["fp_rate"] = bin_stats["fp_n"] / bin_stats["n"]
valid = bin_stats[bin_stats["n"] >= 50]
ax.bar(range(len(valid)), valid["fp_rate"].values, color="#e41a1c", alpha=0.7, edgecolor="black", linewidth=0.3)
ax.set_xticks(range(len(valid)))
ax.set_xticklabels(valid["ehr_bin"].values, rotation=60, ha="right", fontsize=7)
for b_val, b_name in zip([0, 2, 6, 10], ["C0", "B→", "A→", "A+→"]):
    if b_val < len(valid):
        ax.axvline(x=b_val - 0.5, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.set_ylabel("FP Rate")
ax.set_title("TO Mode: FP Rate by effective_hp_reads (fine-grained bins)", fontsize=13, fontweight="bold")
# Add n labels
for k, (_, row) in enumerate(valid.iterrows()):
    if row["n"] >= 100:
        ax.text(k, row["fp_rate"] + 0.01, f'n={row["n"]:,}', ha="center", fontsize=5, rotation=90)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig05_to_fp_rate_by_ehr_fine.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig05_to_fp_rate_by_ehr_fine.png")


# ============================================================
# Part 3: Quality Score Analysis
# ============================================================
print("\n[3/3] Part 3 — Quality Score analysis ...")

# --- P3a: Quality Score vs truth (by mode) ---
rows_p3a = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for truth in ["TP", "FP"]:
        sub = dm[dm["truth_label"] == truth]
        qs = sub["quality_score"].dropna()
        rows_p3a.append({
            "mode": mode, "truth": truth, "n": len(sub),
            "qs_mean": round(qs.mean(), 2), "qs_median": round(qs.median(), 2),
            "qs_q25": round(qs.quantile(0.25), 2), "qs_q75": round(qs.quantile(0.75), 2),
        })
df_p3a = pd.DataFrame(rows_p3a)
write_tsv(df_p3a, OUT_DIR / "p3a_quality_score_by_truth.tsv", "Quality Score by truth")

# --- P3b: Quality Score components decomposition ---
# Reconstruct QS components
dm = df.copy()
dm["qs_read_penalty"] = np.where(dm["num_reads"] < 30, -20, np.where(dm["num_reads"] < 50, -10, 0))
dm["qs_cpg_penalty"] = np.where(dm["num_cpgs"] < 20, -15, np.where(dm["num_cpgs"] < 30, -10, 0))
cov_mult = dm["num_reads"] / 75.0  # approximate coverage multiple
dm["qs_cov_penalty"] = np.where(cov_mult < 0.5, -30, np.where(cov_mult < 0.8, -15, np.where(cov_mult > 2.0, -20, 0)))
dm["qs_loh_penalty"] = np.where(dm["core_loh_like"], -25, 0)
dm["qs_verify_bonus"] = np.where(dm["HPMergedSig"] & dm["AlleleSig"].map({"True": True, "False": False, True: True, False: False}).fillna(False), 10, 0)
dm["qs_effect_bonus"] = np.where(dm["cramers_v"] >= 0.3, 15, np.where(dm["cramers_v"] >= 0.2, 5, 0))

rows_p3b = []
for mode in ["paired", "to"]:
    mm = dm[dm["mode"] == mode]
    for truth in ["TP", "FP"]:
        sub = mm[mm["truth_label"] == truth]
        row = {"mode": mode, "truth": truth, "n": len(sub)}
        for comp in ["qs_read_penalty", "qs_cpg_penalty", "qs_cov_penalty", "qs_loh_penalty", "qs_verify_bonus", "qs_effect_bonus"]:
            row[f"{comp}_mean"] = round(sub[comp].mean(), 3)
            row[f"{comp}_nonzero_pct"] = round((sub[comp] != 0).mean() * 100, 2)
        rows_p3b.append(row)
df_p3b = pd.DataFrame(rows_p3b)
write_tsv(df_p3b, OUT_DIR / "p3b_qs_components.tsv", "QS component decomposition")

# --- P3c: LOH penalty impact analysis ---
# What if we remove the LOH penalty?
rows_p3c = []
for mode in ["paired", "to"]:
    mm = dm[dm["mode"] == mode]
    for truth in ["TP", "FP"]:
        sub = mm[mm["truth_label"] == truth]
        current_qs = sub["quality_score"]
        no_loh_qs = (current_qs - sub["qs_loh_penalty"]).clip(0, 100)
        rows_p3c.append({
            "mode": mode, "truth": truth,
            "current_qs_median": round(current_qs.median(), 2),
            "no_loh_qs_median": round(no_loh_qs.median(), 2),
            "delta_median": round(no_loh_qs.median() - current_qs.median(), 2),
            "loh_penalty_affected_pct": round((sub["qs_loh_penalty"] != 0).mean() * 100, 2),
        })
df_p3c = pd.DataFrame(rows_p3c)
write_tsv(df_p3c, OUT_DIR / "p3c_loh_penalty_impact.tsv", "Impact of removing LOH penalty")

# --- P3d: Methylation features as QS predictors ---
# Can methylation features predict TP/FP better than current QS?
METHYL_FEATURES = ["quality_score", "pairwise_median_dist", "cramers_v", "allele_delta", "caller_af", "hp0_ratio", "hp_assign_rate"]
rows_p3d = []
for mode in ["paired", "to"]:
    mm = df[df["mode"] == mode]
    for feat in METHYL_FEATURES:
        vals_tp = mm[mm["truth_label"] == "TP"][feat].dropna()
        vals_fp = mm[mm["truth_label"] == "FP"][feat].dropna()
        if len(vals_tp) < 10 or len(vals_fp) < 10:
            continue
        # Mann-Whitney U test
        stat, pval = stats.mannwhitneyu(vals_tp, vals_fp, alternative="two-sided")
        # Effect size (rank-biserial correlation)
        n1, n2 = len(vals_tp), len(vals_fp)
        r = 1 - (2 * stat) / (n1 * n2)
        # AUC approximation
        auc = stat / (n1 * n2)
        rows_p3d.append({
            "mode": mode, "feature": feat,
            "tp_median": round(vals_tp.median(), 4),
            "fp_median": round(vals_fp.median(), 4),
            "delta_median": round(vals_tp.median() - vals_fp.median(), 4),
            "mannwhitney_p": f"{pval:.2e}",
            "rank_biserial_r": round(r, 4),
            "auc": round(auc, 4),
        })
df_p3d = pd.DataFrame(rows_p3d)
write_tsv(df_p3d, OUT_DIR / "p3d_feature_discriminability.tsv", "Feature TP/FP discriminability")

# --- Fig06: Quality Score distribution (TP vs FP) ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    dm_mode = df[df["mode"] == mode]
    for truth in ["TP", "FP"]:
        sub = dm_mode[dm_mode["truth_label"] == truth]
        ax.hist(sub["quality_score"].dropna(), bins=50, alpha=0.6, color=TRUTH_COLORS[truth],
                label=f"{truth} (n={len(sub):,})", density=True, edgecolor="none")
    ax.set_xlabel("Quality Score")
    ax.set_ylabel("Density")
    ax.set_title(f"{mode.upper()} — Quality Score Distribution", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
fig.suptitle("Part 3: Quality Score — TP vs FP", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig06_quality_score_distribution.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig06_quality_score_distribution.png")

# --- Fig07: QS component impact (TO mode) ---
fig, ax = plt.subplots(figsize=(12, 5))
comps = ["qs_read_penalty", "qs_cpg_penalty", "qs_cov_penalty", "qs_loh_penalty", "qs_verify_bonus", "qs_effect_bonus"]
comp_labels = ["Read Count\nPenalty", "CpG Count\nPenalty", "Coverage\nPenalty", "LOH\nPenalty", "Verify\nBonus", "Effect Size\nBonus"]
to_tp = df_p3b[(df_p3b["mode"] == "to") & (df_p3b["truth"] == "TP")]
to_fp = df_p3b[(df_p3b["mode"] == "to") & (df_p3b["truth"] == "FP")]
x = np.arange(len(comps))
w = 0.35
tp_vals = [to_tp[f"{c}_mean"].values[0] for c in comps]
fp_vals = [to_fp[f"{c}_mean"].values[0] for c in comps]
ax.bar(x - w/2, tp_vals, w, color=TRUTH_COLORS["TP"], label="TP mean", edgecolor="black", linewidth=0.3)
ax.bar(x + w/2, fp_vals, w, color=TRUTH_COLORS["FP"], label="FP mean", edgecolor="black", linewidth=0.3)
for k in range(len(comps)):
    ax.text(x[k] - w/2, tp_vals[k] + (0.3 if tp_vals[k] >= 0 else -0.8), f'{tp_vals[k]:.1f}', ha="center", fontsize=7)
    ax.text(x[k] + w/2, fp_vals[k] + (0.3 if fp_vals[k] >= 0 else -0.8), f'{fp_vals[k]:.1f}', ha="center", fontsize=7)
ax.set_xticks(x)
ax.set_xticklabels(comp_labels, fontsize=9)
ax.set_ylabel("Mean Score Impact")
ax.set_title("TO Mode: Quality Score Component Impact — TP vs FP", fontsize=13, fontweight="bold")
ax.axhline(y=0, color="black", linewidth=0.5)
ax.legend()
plt.tight_layout()
fig.savefig(FIG_DIR / "fig07_qs_component_impact_to.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig07_qs_component_impact_to.png")

# --- Fig08: Feature AUC comparison ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_p3d[df_p3d["mode"] == mode].sort_values("auc", ascending=True)
    colors = ["#e41a1c" if auc > 0.5 else "#4daf4a" for auc in sub["auc"]]
    ax.barh(range(len(sub)), sub["auc"].values, color=colors, edgecolor="black", linewidth=0.3)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(sub["feature"].values, fontsize=9)
    ax.axvline(x=0.5, color="black", linestyle="--", linewidth=1)
    for k, (_, row) in enumerate(sub.iterrows()):
        ax.text(row["auc"] + 0.005, k, f'{row["auc"]:.3f}', va="center", fontsize=8)
    ax.set_xlabel("AUC (>0.5 = FP higher)")
    ax.set_title(f"{mode.upper()} — Feature Discriminability", fontsize=12, fontweight="bold")
fig.suptitle("Part 3: Feature AUC for TP/FP Discrimination", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig08_feature_auc_comparison.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig08_feature_auc_comparison.png")

print(f"\n{'='*60}")
print(f"Analysis complete!")
print(f"  Output: {OUT_DIR}")
print(f"  TSV: {len(list(OUT_DIR.glob('*.tsv')))}")
print(f"  Figures: {len(list(FIG_DIR.glob('*.png')))}")
print(f"{'='*60}")
