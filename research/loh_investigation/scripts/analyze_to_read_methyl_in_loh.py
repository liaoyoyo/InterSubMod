#!/usr/bin/env python3
"""
Analyze TO mode read-level methylation inside LOH zones.

Extends the paired-only LOH methyl_mean analysis to TO mode data.
Uses sample400 shard export (31,317 TO reads from HCC1395 + HCC1395_DORADO)
plus the paired-only export (86,521 reads) for comparison.

Key question: Does methyl_mean still separate TP/FP inside LOH zones
in TO mode, where ISM methylation analysis is known to fail?

Output:
  - research/loh_investigation/figures/to_read_methyl_fig{01..12}*.png
  - research/loh_investigation/data/to_read_methyl_stats.tsv
  - research/loh_investigation/reports/20260402_to_read_methyl_in_loh.md
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ──────────────────────────────────────────────────────────────────

PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
RESEARCH_DIR = PROJECT_ROOT / "research" / "loh_investigation"
FIG_DIR = RESEARCH_DIR / "figures"
DATA_DIR = RESEARCH_DIR / "data"
REPORT_DIR = RESEARCH_DIR / "reports"

# Read training tables
TO_TABLE = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
                "20260325_phase1_manifest_shard_export_sample400/"
                "phase1_shard_read_training_table.tsv")

PAIRED_TABLE = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
                    "20260325_phase1_manifest_shard_export_paired_multibio_sample637/"
                    "phase1_shard_read_training_table.tsv")

# LOH.bed paths (same as O15)
LOH_BED_PATHS = {
    "HCC1395": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1395_DORADO": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
}

# ── LOH zone assignment ────────────────────────────────────────────────────

def load_loh_intervals(bed_path):
    """Load LOH.bed as list of (chr, start, end)."""
    intervals = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            intervals.setdefault(chrom, []).append((start, end))
    # Sort intervals for binary search
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def is_in_loh(chrom, pos, intervals):
    """Check if a position falls inside any LOH interval."""
    if chrom not in intervals:
        return False
    for start, end in intervals[chrom]:
        if start <= pos < end:
            return True
        if pos < start:
            break
    return False


def assign_loh_zone(df, sample_col="sample"):
    """Assign LOH zone to each read based on region position."""
    # Extract region position from region_key (chr:pos:ref:alt)
    def parse_region_pos(rk):
        parts = rk.split(":")
        return parts[0], int(parts[1])

    loh_cache = {}
    for sample in df[sample_col].unique():
        bed_path = LOH_BED_PATHS.get(sample)
        if bed_path and os.path.exists(bed_path):
            loh_cache[sample] = load_loh_intervals(bed_path)
        else:
            loh_cache[sample] = {}

    zones = []
    for _, row in df.iterrows():
        sample = row[sample_col]
        chrom, pos = parse_region_pos(row["region_key"])
        intervals = loh_cache.get(sample, {})
        zones.append("inside" if is_in_loh(chrom, pos, intervals) else "outside")
    df["loh_zone"] = zones
    return df


# ── Load data ──────────────────────────────────────────────────────────────

def load_data():
    """Load both TO and paired read training tables."""
    print("Loading TO table (sample400)...")
    to_df = pd.read_csv(TO_TABLE, sep="\t")
    to_df = to_df[to_df["mode"] == "to-pure"].copy()
    to_df["mode_label"] = "TO"
    print(f"  TO reads: {len(to_df)}")

    # Also load paired from sample400 for same-sample comparison
    paired_from_400 = pd.read_csv(TO_TABLE, sep="\t")
    paired_from_400 = paired_from_400[paired_from_400["mode"] == "paired-pure"].copy()
    paired_from_400["mode_label"] = "Paired"
    print(f"  Paired reads (from sample400): {len(paired_from_400)}")

    # Also load full paired for broader comparison
    print("Loading full paired table...")
    paired_full = pd.read_csv(PAIRED_TABLE, sep="\t")
    paired_full["mode_label"] = "Paired"
    print(f"  Paired reads (full 637): {len(paired_full)}")

    # Map sample names
    sample_map = {
        "hcc1395_5khz_paired": "HCC1395",
        "hcc1395_5khz_to": "HCC1395",
        "hcc1395_dorado_paired": "HCC1395_DORADO",
        "hcc1395_dorado_to": "HCC1395_DORADO",
    }
    to_df["sample"] = to_df["dataset_id"].map(sample_map)
    paired_from_400["sample"] = paired_from_400["dataset_id"].map(sample_map)

    # Full paired has different dataset_id format — check
    paired_ds = paired_full["dataset_id"].unique()
    print(f"  Paired full dataset_ids: {paired_ds[:5]}")

    return to_df, paired_from_400, paired_full


def compute_auc_safe(y_true, y_score):
    """Compute AUC safely, returning NaN if impossible."""
    if len(set(y_true)) < 2:
        return np.nan
    try:
        return roc_auc_score(y_true, y_score)
    except ValueError:
        return np.nan


def cohens_d(x, y):
    """Compute Cohen's d with pooled standard deviation."""
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    mx, my = np.mean(x), np.mean(y)
    sx, sy = np.std(x, ddof=1), np.std(y, ddof=1)
    sp = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
    if sp < 1e-12:
        return np.nan
    return (mx - my) / sp


# ── Analysis ───────────────────────────────────────────────────────────────

def analyze_methyl_features(df, mode_label, report_lines):
    """Analyze methylation features by LOH zone and truth status."""
    features = ["methyl_mean", "methyl_std", "methyl_median",
                 "methyl_high_fraction", "methyl_low_fraction", "methyl_na_fraction"]

    results = []
    for zone in ["inside", "outside"]:
        zone_df = df[df["loh_zone"] == zone]
        tp = zone_df[zone_df["truth_status"] == "TP"]
        fp = zone_df[zone_df["truth_status"] == "FP"]

        for feat in features:
            tp_vals = tp[feat].dropna()
            fp_vals = fp[feat].dropna()
            if len(tp_vals) < 5 or len(fp_vals) < 5:
                continue

            auc = compute_auc_safe(
                np.concatenate([np.ones(len(tp_vals)), np.zeros(len(fp_vals))]),
                np.concatenate([tp_vals.values, fp_vals.values])
            )
            d = cohens_d(tp_vals.values, fp_vals.values)

            results.append({
                "mode": mode_label,
                "zone": zone,
                "feature": feat,
                "n_tp": len(tp_vals),
                "n_fp": len(fp_vals),
                "tp_median": np.median(tp_vals),
                "fp_median": np.median(fp_vals),
                "auc": auc,
                "cohens_d": d,
            })

    # Region-level aggregation
    for zone in ["inside", "outside"]:
        zone_df = df[df["loh_zone"] == zone]
        region_agg = zone_df.groupby(["region_key", "truth_status"]).agg(
            methyl_mean_region=("methyl_mean", "mean"),
            n_reads=("methyl_mean", "count")
        ).reset_index()

        tp_regions = region_agg[region_agg["truth_status"] == "TP"]
        fp_regions = region_agg[region_agg["truth_status"] == "FP"]

        if len(tp_regions) >= 3 and len(fp_regions) >= 3:
            auc_region = compute_auc_safe(
                np.concatenate([np.ones(len(tp_regions)), np.zeros(len(fp_regions))]),
                np.concatenate([tp_regions["methyl_mean_region"].values,
                                fp_regions["methyl_mean_region"].values])
            )
            d_region = cohens_d(tp_regions["methyl_mean_region"].values,
                                fp_regions["methyl_mean_region"].values)

            results.append({
                "mode": mode_label,
                "zone": zone,
                "feature": "region_methyl_mean",
                "n_tp": len(tp_regions),
                "n_fp": len(fp_regions),
                "tp_median": np.median(tp_regions["methyl_mean_region"]),
                "fp_median": np.median(fp_regions["methyl_mean_region"]),
                "auc": auc_region,
                "cohens_d": d_region,
            })

    return pd.DataFrame(results)


# ── Figures ────────────────────────────────────────────────────────────────

def fig01_methyl_mean_violin(to_df, paired_df, out):
    """Violin: methyl_mean by TP/FP × LOH zone × mode."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for col, mode_label in enumerate(["TO", "Paired"]):
        df = to_df if mode_label == "TO" else paired_df
        for row, zone in enumerate(["inside", "outside"]):
            ax = axes[row, col]
            zone_df = df[df["loh_zone"] == zone].copy()
            if len(zone_df) == 0:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                ax.set_title(f"{mode_label} — LOH {zone}")
                continue

            sns.violinplot(data=zone_df, x="truth_status", y="methyl_mean",
                           order=["TP", "FP"], palette={"TP": "#2196F3", "FP": "#F44336"},
                           ax=ax, inner="quartile", cut=0)

            # AUC annotation
            tp_vals = zone_df[zone_df["truth_status"] == "TP"]["methyl_mean"].dropna()
            fp_vals = zone_df[zone_df["truth_status"] == "FP"]["methyl_mean"].dropna()
            if len(tp_vals) >= 5 and len(fp_vals) >= 5:
                auc = compute_auc_safe(
                    np.concatenate([np.ones(len(tp_vals)), np.zeros(len(fp_vals))]),
                    np.concatenate([tp_vals.values, fp_vals.values])
                )
                ax.text(0.98, 0.95, f"AUC={auc:.3f}\nTP={len(tp_vals)}, FP={len(fp_vals)}",
                        transform=ax.transAxes, ha="right", va="top", fontsize=9,
                        bbox=dict(boxstyle="round", fc="wheat", alpha=0.8))

            ax.set_title(f"{mode_label} — LOH {zone}")
            ax.set_ylabel("methyl_mean")
            ax.set_ylim(-0.05, 1.05)

    fig.suptitle("Read-Level methyl_mean by Truth × LOH Zone × Mode", fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig02_auc_comparison_bar(stats_df, out):
    """Bar: AUC comparison across modes × zones × features."""
    plot_df = stats_df[stats_df["feature"].isin(["methyl_mean", "methyl_high_fraction",
                                                   "methyl_low_fraction", "region_methyl_mean"])]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for i, zone in enumerate(["inside", "outside"]):
        ax = axes[i]
        zone_data = plot_df[plot_df["zone"] == zone]
        if len(zone_data) == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            continue

        x = np.arange(len(zone_data["feature"].unique()))
        width = 0.35
        features = zone_data["feature"].unique()

        for j, mode in enumerate(["TO", "Paired"]):
            mode_data = zone_data[zone_data["mode"] == mode]
            aucs = [mode_data[mode_data["feature"] == f]["auc"].values[0]
                    if len(mode_data[mode_data["feature"] == f]) > 0 else 0
                    for f in features]
            offset = -width/2 + j * width
            bars = ax.bar(x + offset, aucs, width, label=mode,
                          color=["#FF9800", "#4CAF50"][j], alpha=0.8)
            for bar, auc_val in zip(bars, aucs):
                if not np.isnan(auc_val):
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                            f"{auc_val:.3f}", ha="center", va="bottom", fontsize=8)

        ax.set_xticks(x)
        ax.set_xticklabels([f.replace("_", "\n") for f in features], fontsize=8)
        ax.set_ylabel("AUC")
        ax.set_title(f"LOH {zone}")
        ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
        ax.set_ylim(0.3, 1.0)
        ax.legend()

    fig.suptitle("Read-Level Methylation AUC: TO vs Paired × LOH Zone", fontsize=13)
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig03_roc_inside_loh(to_df, paired_df, out):
    """ROC curves for methyl_mean inside LOH, both modes."""
    fig, ax = plt.subplots(figsize=(8, 8))

    for mode_label, df, color in [("TO", to_df, "#FF9800"), ("Paired", paired_df, "#4CAF50")]:
        zone_df = df[df["loh_zone"] == "inside"]
        tp = zone_df[zone_df["truth_status"] == "TP"]["methyl_mean"].dropna()
        fp = zone_df[zone_df["truth_status"] == "FP"]["methyl_mean"].dropna()
        if len(tp) < 5 or len(fp) < 5:
            continue

        y_true = np.concatenate([np.ones(len(tp)), np.zeros(len(fp))])
        y_score = np.concatenate([tp.values, fp.values])
        # methyl_mean: TP < FP, so invert for ROC
        fpr, tpr, _ = roc_curve(y_true, -y_score)
        auc = roc_auc_score(y_true, -y_score)
        ax.plot(fpr, tpr, color=color, lw=2,
                label=f"{mode_label} (AUC={auc:.3f}, n_TP={len(tp)}, n_FP={len(fp)})")

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.set_title("ROC: methyl_mean Inside LOH (TP<FP → inverted)")
    ax.legend(loc="lower right")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig04_methyl_mean_histogram(to_df, paired_df, out):
    """Histogram: methyl_mean distribution TP vs FP, inside LOH, both modes."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for i, (mode_label, df) in enumerate([("TO", to_df), ("Paired", paired_df)]):
        ax = axes[i]
        zone_df = df[df["loh_zone"] == "inside"]
        tp = zone_df[zone_df["truth_status"] == "TP"]["methyl_mean"].dropna()
        fp = zone_df[zone_df["truth_status"] == "FP"]["methyl_mean"].dropna()

        bins = np.linspace(0, 1, 51)
        ax.hist(tp, bins=bins, alpha=0.6, color="#2196F3", label=f"TP (n={len(tp)})", density=True)
        ax.hist(fp, bins=bins, alpha=0.6, color="#F44336", label=f"FP (n={len(fp)})", density=True)
        ax.set_xlabel("methyl_mean")
        ax.set_ylabel("Density")
        ax.set_title(f"{mode_label} — LOH inside")
        ax.legend()

    fig.suptitle("methyl_mean Distribution: TP vs FP Inside LOH", fontsize=13)
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig05_region_level_boxplot(to_df, paired_df, out):
    """Boxplot: region-level methyl_mean inside LOH."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for i, (mode_label, df) in enumerate([("TO", to_df), ("Paired", paired_df)]):
        ax = axes[i]
        zone_df = df[df["loh_zone"] == "inside"]
        region_agg = zone_df.groupby(["region_key", "truth_status"]).agg(
            methyl_mean_region=("methyl_mean", "mean"),
            n_reads=("methyl_mean", "count")
        ).reset_index()

        tp_regions = region_agg[region_agg["truth_status"] == "TP"]
        fp_regions = region_agg[region_agg["truth_status"] == "FP"]

        data = []
        for _, r in tp_regions.iterrows():
            data.append({"truth": "TP", "methyl_mean": r["methyl_mean_region"]})
        for _, r in fp_regions.iterrows():
            data.append({"truth": "FP", "methyl_mean": r["methyl_mean_region"]})
        plot_df = pd.DataFrame(data)

        if len(plot_df) > 0:
            sns.boxplot(data=plot_df, x="truth", y="methyl_mean",
                        palette={"TP": "#2196F3", "FP": "#F44336"}, ax=ax)
            sns.stripplot(data=plot_df, x="truth", y="methyl_mean",
                          color="black", alpha=0.4, size=3, ax=ax)

        auc_region = np.nan
        if len(tp_regions) >= 3 and len(fp_regions) >= 3:
            auc_region = compute_auc_safe(
                np.concatenate([np.ones(len(tp_regions)), np.zeros(len(fp_regions))]),
                np.concatenate([tp_regions["methyl_mean_region"].values,
                                fp_regions["methyl_mean_region"].values])
            )
        ax.text(0.98, 0.95,
                f"AUC={auc_region:.3f}\nTP={len(tp_regions)} regions\nFP={len(fp_regions)} regions",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round", fc="wheat", alpha=0.8))
        ax.set_title(f"{mode_label} — Region-Level Inside LOH")
        ax.set_ylabel("Region mean(methyl_mean)")

    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig06_per_sample_auc(to_df, paired_df, out):
    """Per-sample AUC for methyl_mean inside LOH."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for i, (mode_label, df) in enumerate([("TO", to_df), ("Paired", paired_df)]):
        ax = axes[i]
        zone_df = df[df["loh_zone"] == "inside"]
        sample_aucs = []
        for sample in sorted(zone_df["sample"].unique()):
            s_df = zone_df[zone_df["sample"] == sample]
            tp = s_df[s_df["truth_status"] == "TP"]["methyl_mean"].dropna()
            fp = s_df[s_df["truth_status"] == "FP"]["methyl_mean"].dropna()
            if len(tp) >= 5 and len(fp) >= 5:
                auc = compute_auc_safe(
                    np.concatenate([np.ones(len(tp)), np.zeros(len(fp))]),
                    np.concatenate([tp.values, fp.values])
                )
                sample_aucs.append({"sample": sample, "auc": auc,
                                     "n_tp": len(tp), "n_fp": len(fp)})

        if sample_aucs:
            sa_df = pd.DataFrame(sample_aucs)
            bars = ax.bar(range(len(sa_df)), sa_df["auc"], color="#FF9800" if mode_label == "TO" else "#4CAF50")
            ax.set_xticks(range(len(sa_df)))
            ax.set_xticklabels(sa_df["sample"], rotation=45, ha="right", fontsize=9)
            for bar, (_, r) in zip(bars, sa_df.iterrows()):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f"{r['auc']:.3f}\n({r['n_tp']}/{r['n_fp']})",
                        ha="center", va="bottom", fontsize=8)
            ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
            ax.set_ylim(0.3, 1.0)

        ax.set_title(f"{mode_label} — Per-Sample AUC Inside LOH")
        ax.set_ylabel("AUC (methyl_mean)")

    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig07_to_vs_paired_scatter(to_df, paired_df, out):
    """Scatter: region methyl_mean TO vs Paired for shared regions."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Get region-level means for both modes
    for i, zone in enumerate(["inside", "outside"]):
        ax = axes[i]
        to_zone = to_df[to_df["loh_zone"] == zone]
        pa_zone = paired_df[paired_df["loh_zone"] == zone]

        to_reg = to_zone.groupby(["region_key", "truth_status"])["methyl_mean"].mean().reset_index()
        to_reg.columns = ["region_key", "truth_status", "to_methyl"]
        pa_reg = pa_zone.groupby(["region_key", "truth_status"])["methyl_mean"].mean().reset_index()
        pa_reg.columns = ["region_key", "truth_status", "pa_methyl"]

        merged = to_reg.merge(pa_reg, on=["region_key", "truth_status"])
        if len(merged) == 0:
            ax.text(0.5, 0.5, "No shared regions", transform=ax.transAxes, ha="center")
            ax.set_title(f"LOH {zone}")
            continue

        for truth, color, marker in [("TP", "#2196F3", "o"), ("FP", "#F44336", "x")]:
            sub = merged[merged["truth_status"] == truth]
            ax.scatter(sub["to_methyl"], sub["pa_methyl"], c=color, marker=marker,
                       alpha=0.5, s=20, label=f"{truth} (n={len(sub)})")

        ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
        ax.set_xlabel("TO methyl_mean")
        ax.set_ylabel("Paired methyl_mean")
        ax.set_title(f"LOH {zone} — Region methyl_mean")
        ax.legend()
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)

        if len(merged) >= 3:
            r, p = stats.pearsonr(merged["to_methyl"], merged["pa_methyl"])
            ax.text(0.02, 0.95, f"r={r:.3f}, p={p:.1e}", transform=ax.transAxes,
                    va="top", fontsize=9)

    fig.suptitle("Region methyl_mean: TO vs Paired (Shared Regions)", fontsize=13)
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig08_loh_zone_distribution(to_df, paired_df, out):
    """Stacked bar: LOH zone × truth_status distribution."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for i, (mode_label, df) in enumerate([("TO", to_df), ("Paired", paired_df)]):
        ax = axes[i]
        counts = df.groupby(["loh_zone", "truth_status"]).size().unstack(fill_value=0)
        counts.plot(kind="bar", stacked=True, ax=ax,
                    color={"TP": "#2196F3", "FP": "#F44336"})
        ax.set_title(f"{mode_label} — Read Count by Zone × Truth")
        ax.set_ylabel("n_reads")
        ax.set_xlabel("")
        for container in ax.containers:
            ax.bar_label(container, fmt="%d", fontsize=8)

    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig09_potential_loh_cross_check(to_df, out):
    """Cross-check: ISM Potential_LOH vs LOH.bed zone in TO mode."""
    fig, ax = plt.subplots(figsize=(8, 5))
    zone_df = to_df.copy()
    zone_df["Potential_LOH_str"] = zone_df["Potential_LOH"].astype(str)
    ct = pd.crosstab(zone_df["loh_zone"], zone_df["Potential_LOH_str"])
    ct.plot(kind="bar", ax=ax, color={"True": "#FF5722", "False": "#8BC34A"})
    ax.set_title("TO Mode: ISM Potential_LOH vs LOH.bed Zone")
    ax.set_ylabel("n_reads")
    ax.set_xlabel("LOH.bed Zone")
    for container in ax.containers:
        ax.bar_label(container, fmt="%d", fontsize=8)
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig10_methyl_mean_by_potential_loh(to_df, out):
    """Violin: methyl_mean by Potential_LOH × truth inside LOH.bed zone."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    inside_df = to_df[to_df["loh_zone"] == "inside"].copy()
    inside_df["Potential_LOH_str"] = inside_df["Potential_LOH"].astype(str)

    for i, ploh in enumerate(["True", "False"]):
        ax = axes[i]
        sub = inside_df[inside_df["Potential_LOH_str"] == ploh]
        if len(sub) < 10:
            ax.text(0.5, 0.5, f"n={len(sub)}", transform=ax.transAxes, ha="center")
            ax.set_title(f"Potential_LOH={ploh}")
            continue

        sns.violinplot(data=sub, x="truth_status", y="methyl_mean",
                       order=["TP", "FP"], palette={"TP": "#2196F3", "FP": "#F44336"},
                       ax=ax, inner="quartile", cut=0)
        tp = sub[sub["truth_status"] == "TP"]["methyl_mean"].dropna()
        fp = sub[sub["truth_status"] == "FP"]["methyl_mean"].dropna()
        if len(tp) >= 5 and len(fp) >= 5:
            auc = compute_auc_safe(
                np.concatenate([np.ones(len(tp)), np.zeros(len(fp))]),
                np.concatenate([tp.values, fp.values])
            )
            ax.text(0.98, 0.95, f"AUC={auc:.3f}\nTP={len(tp)}, FP={len(fp)}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=9,
                    bbox=dict(boxstyle="round", fc="wheat", alpha=0.8))

        ax.set_title(f"LOH inside, Potential_LOH={ploh}")
        ax.set_ylabel("methyl_mean")
        ax.set_ylim(-0.05, 1.05)

    fig.suptitle("TO: methyl_mean by Truth × Potential_LOH (Inside LOH.bed)", fontsize=13)
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig11_effect_size_heatmap(stats_df, out):
    """Heatmap: Cohen's d for all features × mode × zone."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    for i, zone in enumerate(["inside", "outside"]):
        ax = axes[i]
        zone_data = stats_df[stats_df["zone"] == zone]
        if len(zone_data) == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            continue

        pivot = zone_data.pivot(index="feature", columns="mode", values="cohens_d")
        if len(pivot) == 0:
            continue

        sns.heatmap(pivot, annot=True, fmt=".2f", cmap="RdYlGn_r", center=0,
                    ax=ax, vmin=-1.5, vmax=1.5, linewidths=0.5)
        ax.set_title(f"LOH {zone} — Cohen's d")

    fig.suptitle("Effect Size (Cohen's d): TP - FP by Feature", fontsize=13)
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def fig12_sample_consistency(to_df, out):
    """Heatmap: per-sample per-feature AUC inside LOH (TO mode)."""
    features = ["methyl_mean", "methyl_high_fraction", "methyl_low_fraction",
                 "methyl_std", "methyl_na_fraction"]
    zone_df = to_df[to_df["loh_zone"] == "inside"]

    results = []
    for sample in sorted(zone_df["sample"].unique()):
        s_df = zone_df[zone_df["sample"] == sample]
        for feat in features:
            tp = s_df[s_df["truth_status"] == "TP"][feat].dropna()
            fp = s_df[s_df["truth_status"] == "FP"][feat].dropna()
            auc = compute_auc_safe(
                np.concatenate([np.ones(len(tp)), np.zeros(len(fp))]),
                np.concatenate([tp.values, fp.values])
            ) if len(tp) >= 5 and len(fp) >= 5 else np.nan
            results.append({"sample": sample, "feature": feat, "auc": auc})

    res_df = pd.DataFrame(results)
    pivot = res_df.pivot(index="sample", columns="feature", values="auc")

    fig, ax = plt.subplots(figsize=(10, 4))
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap="RdYlGn", center=0.5,
                ax=ax, vmin=0.3, vmax=0.9, linewidths=0.5)
    ax.set_title("TO Inside LOH: Per-Sample Per-Feature AUC")
    plt.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ── Report ─────────────────────────────────────────────────────────────────

def generate_report(stats_df, to_df, paired_df, out):
    """Generate markdown report."""
    lines = []
    lines.append("<!--")
    lines.append("建立時間: 2026-04-02 22:30")
    lines.append("目標: 擴展 methyl_mean LOH 分析到 TO mode read-level 數據")
    lines.append("處理範圍: sample400 shard (HCC1395 + HCC1395_DORADO) × TO + Paired")
    lines.append("關聯檔案:")
    lines.append("  - research/loh_investigation/scripts/analyze_to_read_methyl_in_loh.py")
    lines.append("  - research/loh_investigation/figures/to_read_methyl_fig*.png")
    lines.append("  - research/loh_investigation/data/to_read_methyl_stats.tsv")
    lines.append("-->")
    lines.append("")
    lines.append("# TO Mode Read-Level Methylation Inside LOH")
    lines.append("")
    lines.append(f"**Date**: 2026-04-02")
    lines.append(f"**Data**: sample400 shard export (HCC1395 × 2 platforms × 2 modes)")
    lines.append(f"**LOH zones**: per-sample LongPhase-TO LOH.bed (binary inside/outside)")
    lines.append("")

    # Data summary
    lines.append("## Data Summary")
    lines.append("")
    for mode_label, df in [("TO", to_df), ("Paired", paired_df)]:
        n_total = len(df)
        n_inside = len(df[df["loh_zone"] == "inside"])
        n_outside = len(df[df["loh_zone"] == "outside"])
        for zone in ["inside", "outside"]:
            z = df[df["loh_zone"] == zone]
            n_tp = len(z[z["truth_status"] == "TP"])
            n_fp = len(z[z["truth_status"] == "FP"])
            n_regions_tp = z[z["truth_status"] == "TP"]["region_key"].nunique()
            n_regions_fp = z[z["truth_status"] == "FP"]["region_key"].nunique()
            lines.append(f"- **{mode_label} {zone}**: {n_tp} TP reads ({n_regions_tp} regions) + "
                         f"{n_fp} FP reads ({n_regions_fp} regions)")
    lines.append("")

    # Main results
    lines.append("## Key Results: methyl_mean Inside LOH")
    lines.append("")
    lines.append("| Mode | Zone | Feature | n_TP | n_FP | TP median | FP median | AUC | Cohen's d |")
    lines.append("|------|------|---------|------|------|-----------|-----------|-----|-----------|")
    key_features = ["methyl_mean", "methyl_high_fraction", "methyl_low_fraction", "region_methyl_mean"]
    for _, row in stats_df[stats_df["feature"].isin(key_features)].iterrows():
        lines.append(f"| {row['mode']} | {row['zone']} | {row['feature']} | "
                     f"{row['n_tp']} | {row['n_fp']} | {row['tp_median']:.4f} | "
                     f"{row['fp_median']:.4f} | {row['auc']:.4f} | {row['cohens_d']:.4f} |")
    lines.append("")

    # All features
    lines.append("## Full Feature Table")
    lines.append("")
    lines.append("| Mode | Zone | Feature | n_TP | n_FP | TP median | FP median | AUC | Cohen's d |")
    lines.append("|------|------|---------|------|------|-----------|-----------|-----|-----------|")
    for _, row in stats_df.iterrows():
        lines.append(f"| {row['mode']} | {row['zone']} | {row['feature']} | "
                     f"{row['n_tp']} | {row['n_fp']} | {row['tp_median']:.4f} | "
                     f"{row['fp_median']:.4f} | {row['auc']:.4f} | {row['cohens_d']:.4f} |")
    lines.append("")

    # Conclusions
    lines.append("## Conclusions")
    lines.append("")

    # Extract key AUCs
    to_inside = stats_df[(stats_df["mode"] == "TO") & (stats_df["zone"] == "inside") &
                          (stats_df["feature"] == "methyl_mean")]
    pa_inside = stats_df[(stats_df["mode"] == "Paired") & (stats_df["zone"] == "inside") &
                          (stats_df["feature"] == "methyl_mean")]
    to_region = stats_df[(stats_df["mode"] == "TO") & (stats_df["zone"] == "inside") &
                          (stats_df["feature"] == "region_methyl_mean")]
    pa_region = stats_df[(stats_df["mode"] == "Paired") & (stats_df["zone"] == "inside") &
                          (stats_df["feature"] == "region_methyl_mean")]

    if len(to_inside) > 0:
        r = to_inside.iloc[0]
        lines.append(f"1. **TO methyl_mean inside LOH**: AUC={r['auc']:.4f}, d={r['cohens_d']:.4f} "
                     f"(n_TP={r['n_tp']}, n_FP={r['n_fp']})")
    if len(pa_inside) > 0:
        r = pa_inside.iloc[0]
        lines.append(f"2. **Paired methyl_mean inside LOH**: AUC={r['auc']:.4f}, d={r['cohens_d']:.4f} "
                     f"(n_TP={r['n_tp']}, n_FP={r['n_fp']})")
    if len(to_region) > 0:
        r = to_region.iloc[0]
        lines.append(f"3. **TO region-level inside LOH**: AUC={r['auc']:.4f}, d={r['cohens_d']:.4f} "
                     f"(n_TP={r['n_tp']}, n_FP={r['n_fp']})")
    if len(pa_region) > 0:
        r = pa_region.iloc[0]
        lines.append(f"4. **Paired region-level inside LOH**: AUC={r['auc']:.4f}, d={r['cohens_d']:.4f} "
                     f"(n_TP={r['n_tp']}, n_FP={r['n_fp']})")
    lines.append("")

    with open(out, "w") as f:
        f.write("\n".join(lines))
    print(f"  Saved: {out}")


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("TO Read-Level Methylation Inside LOH Analysis")
    print("=" * 60)

    # Load data
    to_df, paired_from_400, paired_full = load_data()

    # Use paired_from_400 for same-sample comparison (HCC1395 only)
    paired_df = paired_from_400

    # Assign LOH zones
    print("\nAssigning LOH zones...")
    to_df = assign_loh_zone(to_df)
    paired_df = assign_loh_zone(paired_df)

    print(f"\n  TO zone distribution:")
    print(to_df.groupby(["loh_zone", "truth_status"]).size())
    print(f"\n  Paired zone distribution:")
    print(paired_df.groupby(["loh_zone", "truth_status"]).size())

    # Analyze
    print("\nComputing statistics...")
    to_stats = analyze_methyl_features(to_df, "TO", [])
    paired_stats = analyze_methyl_features(paired_df, "Paired", [])
    stats_df = pd.concat([to_stats, paired_stats], ignore_index=True)

    # Save stats
    stats_path = DATA_DIR / "to_read_methyl_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False)
    print(f"  Stats saved: {stats_path}")

    # Print key results
    print("\n" + "=" * 60)
    print("KEY RESULTS: methyl_mean Inside LOH")
    print("=" * 60)
    key = stats_df[(stats_df["feature"].isin(["methyl_mean", "region_methyl_mean"])) &
                    (stats_df["zone"] == "inside")]
    for _, row in key.iterrows():
        print(f"  {row['mode']:8s} {row['feature']:20s}  AUC={row['auc']:.4f}  "
              f"d={row['cohens_d']:.4f}  TP_med={row['tp_median']:.3f}  FP_med={row['fp_median']:.3f}  "
              f"n_TP={row['n_tp']}  n_FP={row['n_fp']}")

    # Generate figures
    print("\nGenerating figures...")
    fig01_methyl_mean_violin(to_df, paired_df, FIG_DIR / "to_read_methyl_fig01_violin.png")
    fig02_auc_comparison_bar(stats_df, FIG_DIR / "to_read_methyl_fig02_auc_bar.png")
    fig03_roc_inside_loh(to_df, paired_df, FIG_DIR / "to_read_methyl_fig03_roc.png")
    fig04_methyl_mean_histogram(to_df, paired_df, FIG_DIR / "to_read_methyl_fig04_histogram.png")
    fig05_region_level_boxplot(to_df, paired_df, FIG_DIR / "to_read_methyl_fig05_region_boxplot.png")
    fig06_per_sample_auc(to_df, paired_df, FIG_DIR / "to_read_methyl_fig06_per_sample_auc.png")
    fig07_to_vs_paired_scatter(to_df, paired_df, FIG_DIR / "to_read_methyl_fig07_scatter.png")
    fig08_loh_zone_distribution(to_df, paired_df, FIG_DIR / "to_read_methyl_fig08_zone_dist.png")
    fig09_potential_loh_cross_check(to_df, FIG_DIR / "to_read_methyl_fig09_potential_loh.png")
    fig10_methyl_mean_by_potential_loh(to_df, FIG_DIR / "to_read_methyl_fig10_potential_loh_violin.png")
    fig11_effect_size_heatmap(stats_df, FIG_DIR / "to_read_methyl_fig11_effect_size.png")
    fig12_sample_consistency(to_df, FIG_DIR / "to_read_methyl_fig12_sample_consistency.png")

    # Generate report
    print("\nGenerating report...")
    generate_report(stats_df, to_df, paired_df,
                    REPORT_DIR / "20260402_to_read_methyl_in_loh.md")

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
