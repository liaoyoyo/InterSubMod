#!/usr/bin/env python3
"""
Analyze whether TP/FP methylation level difference persists INSIDE LOH zones.

Background (O10): TP methyl_mean=0.463 vs FP=0.679, AUC=0.733.
Core question: Does this separation survive inside LOH regions?

Output:
  - 10 figures  -> research/loh_investigation/figures/loh_methyl_level_fig*.png
  - Stats TSV   -> research/loh_investigation/data/loh_methyl_level_stats.tsv
  - Report MD   -> research/loh_investigation/reports/20260402_loh_methyl_level_inside_loh.md
"""

import sys
import warnings
from pathlib import Path
from textwrap import dedent

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------- style ----------
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 180,
    "savefig.bbox": "tight",
})

COL_TP = "#2196F3"
COL_FP = "#F44336"
COL_INSIDE = "#E53935"
COL_OUTSIDE = "#1E88E5"
TRUTH_PAL = {"TP": COL_TP, "FP": COL_FP}
ZONE_PAL = {"inside": COL_INSIDE, "outside": COL_OUTSIDE}

# ---------- paths ----------
ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation")
FIG_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"
RPT_DIR = ROOT / "reports"
for d in (FIG_DIR, DATA_DIR, RPT_DIR):
    d.mkdir(parents=True, exist_ok=True)

TRAIN_TABLE = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_manifest_shard_export_paired_multibio_sample637/"
    "phase1_shard_read_training_table.tsv"
)

LOH_BED_PATHS = {
    "HCC1395": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1395_DORADO": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "COLO829": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260317_colo829_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed",
    "H1437": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_h1437_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "H2009": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1937": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_hcc1937_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1954": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_hcc1954_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
}

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]

METHYL_FEATURES = [
    "methyl_mean", "methyl_std", "methyl_median",
    "methyl_high_fraction", "methyl_low_fraction", "methyl_na_fraction",
]

# ============================================================
# Helper functions
# ============================================================

def load_bed_intervals(bed_path):
    """Load BED intervals into dict {chrom: [(start, end), ...]} sorted."""
    intervals = {}
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            intervals.setdefault(chrom, []).append((s, e))
    for ch in intervals:
        intervals[ch].sort()
    return intervals


def assign_loh_zone(df, intervals):
    """Assign 'inside' or 'outside' LOH zone for each read via binary search."""
    zones = pd.Series("outside", index=df.index)
    for chrom, grp in df.groupby("chr"):
        if chrom not in intervals:
            continue
        ivs = intervals[chrom]
        starts = np.array([s for s, e in ivs])
        ends = np.array([e for s, e in ivs])
        positions = grp["start"].values
        insert_idx = np.searchsorted(starts, positions, side="right") - 1
        in_range = (insert_idx >= 0) & (insert_idx < len(ivs))
        mask = np.zeros(len(positions), dtype=bool)
        valid = np.where(in_range)[0]
        if len(valid) > 0:
            vi = insert_idx[valid]
            mask[valid] = (positions[valid] >= starts[vi]) & (positions[valid] < ends[vi])
        zones.loc[grp.index[mask]] = "inside"
    return zones


def safe_auc(y_true_binary, scores):
    """Return AUC or NaN if degenerate."""
    if len(np.unique(y_true_binary)) < 2 or len(scores) < 2:
        return np.nan
    return roc_auc_score(y_true_binary, scores)


def cohens_d(x, y):
    """Cohen's d with pooled std."""
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    pooled = np.sqrt(((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2))
    if pooled == 0:
        return 0.0
    return (np.mean(x) - np.mean(y)) / pooled


def bootstrap_d_ci(x, y, n_boot=1000, seed=42, alpha=0.05):
    """Bootstrap 95% CI for Cohen's d."""
    rng = np.random.RandomState(seed)
    ds = []
    for _ in range(n_boot):
        xi = rng.choice(x, size=len(x), replace=True)
        yi = rng.choice(y, size=len(y), replace=True)
        ds.append(cohens_d(xi, yi))
    ds = np.array(ds)
    lo = np.nanpercentile(ds, 100 * alpha / 2)
    hi = np.nanpercentile(ds, 100 * (1 - alpha / 2))
    return lo, hi


def compute_metrics(tp_vals, fp_vals, label=""):
    """Compute AUC, Cohen's d, CI, medians for TP vs FP values."""
    # AUC: TP=1 means LOWER methyl_mean is predictive of TP
    # We use 1-value as score so that lower = higher score
    y = np.concatenate([np.ones(len(tp_vals)), np.zeros(len(fp_vals))])
    scores = np.concatenate([tp_vals, fp_vals])
    # For methyl_mean: FP has HIGHER values, so score = value predicts FP (AUC with FP=1)
    # But we want AUC with TP=1 and a score that is high for TP.
    # Flip: use 1-methyl_mean as score for TP=1.
    auc_val = safe_auc(y, 1.0 - scores)
    d = cohens_d(tp_vals, fp_vals)
    d_lo, d_hi = bootstrap_d_ci(tp_vals, fp_vals) if len(tp_vals) >= 5 and len(fp_vals) >= 5 else (np.nan, np.nan)
    return {
        "label": label,
        "n_tp": len(tp_vals),
        "n_fp": len(fp_vals),
        "tp_mean": np.mean(tp_vals),
        "tp_median": np.median(tp_vals),
        "tp_iqr": np.subtract(*np.percentile(tp_vals, [75, 25])),
        "fp_mean": np.mean(fp_vals),
        "fp_median": np.median(fp_vals),
        "fp_iqr": np.subtract(*np.percentile(fp_vals, [75, 25])),
        "auc": auc_val,
        "cohens_d": d,
        "d_ci_lo": d_lo,
        "d_ci_hi": d_hi,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("LOH Methylation Level Analysis")
    print("=" * 70)

    # ---- Load data ----
    print("\n[1/7] Loading read training table ...")
    cols_needed = [
        "sample", "chr", "start", "truth_status", "region_key", "region_id",
        "alt_support", "hp", "num_cpg_observed", "Potential_LOH",
        "VerificationClass", "Quality_Score",
    ] + METHYL_FEATURES
    df = pd.read_csv(TRAIN_TABLE, sep="\t", usecols=cols_needed)
    df = df[df["truth_status"].isin(["TP", "FP"])].copy()
    print(f"  Loaded {len(df):,} reads ({df['truth_status'].value_counts().to_dict()})")

    # ---- Assign LOH zones per sample ----
    print("\n[2/7] Assigning LOH zones ...")
    df["loh_zone"] = "outside"
    for sample_name in SAMPLE_ORDER:
        bed_path = LOH_BED_PATHS.get(sample_name)
        if bed_path is None or not Path(bed_path).exists():
            print(f"  WARNING: No LOH.bed for {sample_name}")
            continue
        intervals = load_bed_intervals(bed_path)
        mask = df["sample"] == sample_name
        df.loc[mask, "loh_zone"] = assign_loh_zone(df.loc[mask], intervals).values
        n_in = (df.loc[mask, "loh_zone"] == "inside").sum()
        n_out = (df.loc[mask, "loh_zone"] == "outside").sum()
        print(f"  {sample_name:20s}: inside={n_in:>6,}  outside={n_out:>6,}")

    zone_counts = df["loh_zone"].value_counts()
    print(f"\n  Total: inside={zone_counts.get('inside',0):,}  outside={zone_counts.get('outside',0):,}")

    # ---- Part 1: Zone Distribution ----
    print("\n[3/7] Part 1 — Zone distribution ...")
    dist_reads = df.groupby(["sample", "loh_zone", "truth_status"]).size().reset_index(name="n_reads")
    dist_regions = df.groupby(["sample", "loh_zone", "truth_status"])["region_key"].nunique().reset_index(name="n_regions")
    dist = dist_reads.merge(dist_regions, on=["sample", "loh_zone", "truth_status"])
    print(dist.to_string(index=False))

    # Fig 01 — Zone distribution stacked bar
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ax, zone in zip(axes, ["inside", "outside"]):
        sub = dist[dist["loh_zone"] == zone].copy()
        sub["sample"] = pd.Categorical(sub["sample"], categories=SAMPLE_ORDER, ordered=True)
        sub = sub.sort_values("sample")
        pivot = sub.pivot(index="sample", columns="truth_status", values="n_reads").fillna(0)
        pivot = pivot.reindex(columns=["TP", "FP"])
        pivot.plot(kind="bar", stacked=True, color=[COL_TP, COL_FP], ax=ax, edgecolor="white")
        ax.set_title(f"LOH {zone}", fontsize=14, fontweight="bold")
        ax.set_ylabel("Read count")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=45)
        for container in ax.containers:
            ax.bar_label(container, fmt="%d", fontsize=7, label_type="center")
    fig.suptitle("Fig01: Read Counts by LOH Zone x Truth Status x Sample", fontsize=15, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG_DIR / "loh_methyl_level_fig01_zone_distribution.png")
    plt.close(fig)
    print("  -> fig01 saved")

    # ---- Part 2: Methylation Level Comparison ----
    print("\n[4/7] Part 2 — Core methylation comparison ...")
    all_metrics = []
    for zone in ["inside", "outside"]:
        sub = df[df["loh_zone"] == zone]
        tp_vals = sub.loc[sub["truth_status"] == "TP", "methyl_mean"].dropna().values
        fp_vals = sub.loc[sub["truth_status"] == "FP", "methyl_mean"].dropna().values
        m = compute_metrics(tp_vals, fp_vals, label=f"methyl_mean|{zone}")
        all_metrics.append(m)
        print(f"  {zone:8s} methyl_mean: AUC={m['auc']:.4f}  d={m['cohens_d']:.3f} [{m['d_ci_lo']:.3f}, {m['d_ci_hi']:.3f}]  "
              f"TP_med={m['tp_median']:.3f}  FP_med={m['fp_median']:.3f}  n_TP={m['n_tp']:,}  n_FP={m['n_fp']:,}")

    # Fig 02 — THE KEY FIGURE: Violin 2x2
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, zone in zip(axes, ["inside", "outside"]):
        sub = df[df["loh_zone"] == zone].copy()
        sns.violinplot(data=sub, x="truth_status", y="methyl_mean", hue="truth_status",
                       palette=TRUTH_PAL, order=["TP", "FP"],
                       inner="quartile", cut=0, ax=ax, density_norm="width", legend=False)
        # Overlay strip
        sns.stripplot(data=sub.sample(min(2000, len(sub)), random_state=42),
                      x="truth_status", y="methyl_mean", hue="truth_status",
                      palette=TRUTH_PAL, order=["TP", "FP"],
                      size=1.5, alpha=0.15, jitter=0.25, ax=ax, legend=False)
        # Annotate AUC
        m = [x for x in all_metrics if x["label"] == f"methyl_mean|{zone}"][0]
        ax.set_title(f"LOH {zone}\nAUC={m['auc']:.3f}  d={m['cohens_d']:.3f}\nn_TP={m['n_tp']:,}  n_FP={m['n_fp']:,}",
                     fontsize=12, fontweight="bold")
        ax.set_ylabel("methyl_mean" if zone == "inside" else "")
        ax.set_xlabel("")
    fig.suptitle("Fig02: methyl_mean Distribution by LOH Zone x Truth Status (Core Question)", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "loh_methyl_level_fig02_methyl_mean_by_zone_truth.png")
    plt.close(fig)
    print("  -> fig02 saved")

    # Compute all features x zones
    feature_zone_metrics = []
    for feat in METHYL_FEATURES:
        for zone in ["inside", "outside"]:
            sub = df[df["loh_zone"] == zone]
            tp_vals = sub.loc[sub["truth_status"] == "TP", feat].dropna().values
            fp_vals = sub.loc[sub["truth_status"] == "FP", feat].dropna().values
            m = compute_metrics(tp_vals, fp_vals, label=f"{feat}|{zone}")
            feature_zone_metrics.append(m)
    all_metrics.extend(feature_zone_metrics)

    # Fig 03 — 6-panel violin inside LOH only
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes_flat = axes.flatten()
    for i, feat in enumerate(METHYL_FEATURES):
        ax = axes_flat[i]
        sub = df[df["loh_zone"] == "inside"].copy()
        sns.violinplot(data=sub, x="truth_status", y=feat, hue="truth_status",
                       palette=TRUTH_PAL, order=["TP", "FP"],
                       inner="quartile", cut=0, ax=ax, density_norm="width", legend=False)
        m = [x for x in feature_zone_metrics if x["label"] == f"{feat}|inside"][0]
        ax.set_title(f"{feat}\nAUC={m['auc']:.3f}  d={m['cohens_d']:.3f}", fontsize=11, fontweight="bold")
        ax.set_xlabel("")
    fig.suptitle("Fig03: Methylation Features TP vs FP — LOH Inside Only", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "loh_methyl_level_fig03_methyl_features_inside_loh.png")
    plt.close(fig)
    print("  -> fig03 saved")

    # Fig 04 — AUC comparison bar chart
    auc_rows = []
    for m in feature_zone_metrics:
        feat, zone = m["label"].split("|")
        auc_rows.append({"feature": feat, "zone": zone, "auc": m["auc"]})
    auc_df = pd.DataFrame(auc_rows)
    fig, ax = plt.subplots(figsize=(14, 6))
    pivot = auc_df.pivot(index="feature", columns="zone", values="auc")
    pivot = pivot.reindex(METHYL_FEATURES)
    x = np.arange(len(pivot))
    w = 0.35
    bars1 = ax.bar(x - w / 2, pivot["inside"], w, label="inside LOH", color=COL_INSIDE, alpha=0.85)
    bars2 = ax.bar(x + w / 2, pivot["outside"], w, label="outside LOH", color=COL_OUTSIDE, alpha=0.85)
    ax.axhline(0.733, color="gray", ls="--", lw=1.5, label="O10 baseline (0.733)")
    ax.axhline(0.5, color="lightgray", ls=":", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("methyl_", "m_") for f in METHYL_FEATURES], rotation=30, ha="right")
    ax.set_ylabel("AUC (TP=1, higher = better separation)")
    ax.set_title("Fig04: AUC Inside vs Outside LOH by Methyl Feature", fontsize=14, fontweight="bold")
    ax.legend(loc="lower right")
    for bar_group in [bars1, bars2]:
        for bar in bar_group:
            h = bar.get_height()
            if not np.isnan(h):
                ax.annotate(f"{h:.3f}", xy=(bar.get_x() + bar.get_width() / 2, h),
                            xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "loh_methyl_level_fig04_auc_comparison.png")
    plt.close(fig)
    print("  -> fig04 saved")

    # ---- Part 3: Per-Sample Breakdown ----
    print("\n[5/7] Part 3 — Per-sample methyl_mean AUC inside LOH ...")
    per_sample_metrics = []
    for sample in SAMPLE_ORDER:
        sub = df[(df["sample"] == sample) & (df["loh_zone"] == "inside")]
        tp_vals = sub.loc[sub["truth_status"] == "TP", "methyl_mean"].dropna().values
        fp_vals = sub.loc[sub["truth_status"] == "FP", "methyl_mean"].dropna().values
        if len(tp_vals) < 5 or len(fp_vals) < 5:
            print(f"  {sample:20s}: SKIPPED (n_TP={len(tp_vals)}, n_FP={len(fp_vals)})")
            per_sample_metrics.append({"sample": sample, "auc": np.nan, "cohens_d": np.nan,
                                       "n_tp": len(tp_vals), "n_fp": len(fp_vals),
                                       "tp_median": np.nan, "fp_median": np.nan})
            continue
        m = compute_metrics(tp_vals, fp_vals, label=f"methyl_mean|inside|{sample}")
        per_sample_metrics.append({"sample": sample, "auc": m["auc"], "cohens_d": m["cohens_d"],
                                   "n_tp": m["n_tp"], "n_fp": m["n_fp"],
                                   "tp_median": m["tp_median"], "fp_median": m["fp_median"]})
        print(f"  {sample:20s}: AUC={m['auc']:.4f}  d={m['cohens_d']:.3f}  n_TP={m['n_tp']:,}  n_FP={m['n_fp']:,}")
    ps_df = pd.DataFrame(per_sample_metrics)

    # Fig 05 — Faceted violin per sample inside LOH
    inside_df = df[df["loh_zone"] == "inside"].copy()
    inside_df["sample"] = pd.Categorical(inside_df["sample"], categories=SAMPLE_ORDER, ordered=True)
    valid_samples = inside_df.groupby("sample")["truth_status"].apply(lambda x: x.nunique() == 2)
    valid_samples = valid_samples[valid_samples].index.tolist()
    plot_df = inside_df[inside_df["sample"].isin(valid_samples)]

    n_valid = len(valid_samples)
    ncols = min(4, n_valid)
    nrows = (n_valid + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), sharey=True)
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]
    for i, sample in enumerate(valid_samples):
        ax = axes_flat[i]
        s_df = plot_df[plot_df["sample"] == sample]
        sns.violinplot(data=s_df, x="truth_status", y="methyl_mean", hue="truth_status",
                       palette=TRUTH_PAL, order=["TP", "FP"],
                       inner="quartile", cut=0, ax=ax, density_norm="width", legend=False)
        row = ps_df[ps_df["sample"] == sample].iloc[0]
        ax.set_title(f"{sample}\nAUC={row['auc']:.3f}  n={row['n_tp']+row['n_fp']:,}", fontsize=10, fontweight="bold")
        ax.set_xlabel("")
        if i % ncols != 0:
            ax.set_ylabel("")
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)
    fig.suptitle("Fig05: methyl_mean TP vs FP per Sample — LOH Inside Only", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "loh_methyl_level_fig05_per_sample_methyl_mean.png")
    plt.close(fig)
    print("  -> fig05 saved")

    # ---- Part 4: Region-Level Aggregation ----
    print("\n[6/7] Part 4 — Region-level aggregation ...")
    region_agg_metrics = {}
    for zone in ["inside", "outside"]:
        sub = df[df["loh_zone"] == zone]
        region_means = sub.groupby(["region_key", "truth_status"])["methyl_mean"].mean().reset_index()
        tp_rm = region_means.loc[region_means["truth_status"] == "TP", "methyl_mean"].dropna().values
        fp_rm = region_means.loc[region_means["truth_status"] == "FP", "methyl_mean"].dropna().values
        m = compute_metrics(tp_rm, fp_rm, label=f"region_methyl_mean|{zone}")
        region_agg_metrics[zone] = m
        n_tp_reg = len(tp_rm)
        n_fp_reg = len(fp_rm)
        print(f"  {zone:8s} region-level: AUC={m['auc']:.4f}  d={m['cohens_d']:.3f}  "
              f"n_TP_regions={n_tp_reg}  n_FP_regions={n_fp_reg}")
        all_metrics.append(m)

    # Fig 06 — Region-level vs read-level comparison
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    # Panel A: inside LOH region-level distributions
    ax = axes[0]
    sub = df[df["loh_zone"] == "inside"]
    region_means = sub.groupby(["region_key", "truth_status"])["methyl_mean"].mean().reset_index()
    sns.violinplot(data=region_means, x="truth_status", y="methyl_mean", hue="truth_status",
                   palette=TRUTH_PAL, order=["TP", "FP"],
                   inner="quartile", cut=0, ax=ax, density_norm="width", legend=False)
    m_in = region_agg_metrics["inside"]
    ax.set_title(f"Region-Level Inside LOH\nAUC={m_in['auc']:.3f}  d={m_in['cohens_d']:.3f}", fontsize=11, fontweight="bold")
    ax.set_xlabel("")

    # Panel B: outside LOH region-level
    ax = axes[1]
    sub = df[df["loh_zone"] == "outside"]
    region_means = sub.groupby(["region_key", "truth_status"])["methyl_mean"].mean().reset_index()
    sns.violinplot(data=region_means, x="truth_status", y="methyl_mean", hue="truth_status",
                   palette=TRUTH_PAL, order=["TP", "FP"],
                   inner="quartile", cut=0, ax=ax, density_norm="width", legend=False)
    m_out = region_agg_metrics["outside"]
    ax.set_title(f"Region-Level Outside LOH\nAUC={m_out['auc']:.3f}  d={m_out['cohens_d']:.3f}", fontsize=11, fontweight="bold")
    ax.set_xlabel("")

    # Panel C: comparison table
    ax = axes[2]
    ax.axis("off")
    # Build comparison data
    read_inside = [x for x in all_metrics if x["label"] == "methyl_mean|inside"][0]
    read_outside = [x for x in all_metrics if x["label"] == "methyl_mean|outside"][0]
    table_data = [
        ["", "Inside LOH", "Outside LOH"],
        ["Read-level AUC", f"{read_inside['auc']:.3f}", f"{read_outside['auc']:.3f}"],
        ["Region-level AUC", f"{m_in['auc']:.3f}", f"{m_out['auc']:.3f}"],
        ["Read-level d", f"{read_inside['cohens_d']:.3f}", f"{read_outside['cohens_d']:.3f}"],
        ["Region-level d", f"{m_in['cohens_d']:.3f}", f"{m_out['cohens_d']:.3f}"],
        ["n (read TP/FP)", f"{read_inside['n_tp']}/{read_inside['n_fp']}", f"{read_outside['n_tp']}/{read_outside['n_fp']}"],
        ["n (region TP/FP)", f"{m_in['n_tp']}/{m_in['n_fp']}", f"{m_out['n_tp']}/{m_out['n_fp']}"],
    ]
    tbl = ax.table(cellText=table_data, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.2, 1.6)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#EEEEEE")
            cell.set_text_props(fontweight="bold")
        if col == 0:
            cell.set_text_props(fontweight="bold")
    ax.set_title("Read vs Region Level Comparison", fontsize=11, fontweight="bold")

    fig.suptitle("Fig06: Region-Level vs Read-Level methyl_mean — LOH Stratified", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "loh_methyl_level_fig06_region_level_vs_read_level.png")
    plt.close(fig)
    print("  -> fig06 saved")

    # ---- Fig 07: methyl_mean by alt_support inside LOH ----
    print("\n  Generating fig07-09 ...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    inside = df[df["loh_zone"] == "inside"].copy()
    for ax, truth in zip(axes, ["TP", "FP"]):
        sub = inside[inside["truth_status"] == truth]
        sns.violinplot(data=sub, x="alt_support", y="methyl_mean",
                       hue="alt_support", palette={"ALT": "#FF9800", "REF": "#4CAF50"},
                       inner="quartile", cut=0, ax=ax, density_norm="width", legend=True)
        ax.set_title(f"{truth} — Inside LOH", fontsize=12, fontweight="bold")
        ax.set_xlabel("alt_support")
        # Add counts
        counts = sub["alt_support"].value_counts()
        ax.set_xlabel(f"alt_support  (ALT={counts.get('ALT',0):,} / REF={counts.get('REF',0):,})")
    fig.suptitle("Fig07: methyl_mean by ALT/REF Support — LOH Inside", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "loh_methyl_level_fig07_methyl_by_alt_support_loh.png")
    plt.close(fig)
    print("  -> fig07 saved")

    # ---- Fig 08: methyl_mean by HP inside LOH ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    # HP values: 0=unphased, 1=hp1, 2=hp2, 1-1=hp1(alt), 2-1=hp2(alt), 3=ambiguous
    hp_pal = {"0": "#BDBDBD", "1": "#9C27B0", "2": "#009688",
              "1-1": "#CE93D8", "2-1": "#80CBC4", "3": "#757575"}
    for ax, truth in zip(axes, ["TP", "FP"]):
        sub = inside[inside["truth_status"] == truth]
        hp_vals = sorted(sub["hp"].dropna().unique(), key=str)
        hp_order = [h for h in ["0", "1", "1-1", "2", "2-1", "3"] if h in hp_vals]
        pal_sub = {k: v for k, v in hp_pal.items() if k in hp_order}
        sns.violinplot(data=sub, x="hp", y="methyl_mean",
                       hue="hp", palette=pal_sub, order=hp_order,
                       inner="quartile", cut=0, ax=ax, density_norm="width", legend=True)
        ax.set_title(f"{truth} — Inside LOH", fontsize=12, fontweight="bold")
        ax.set_xlabel("HP tag")
    fig.suptitle("Fig08: methyl_mean by HP Tag — LOH Inside", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "loh_methyl_level_fig08_methyl_by_hp_loh.png")
    plt.close(fig)
    print("  -> fig08 saved")

    # ---- Fig 09: Confound check ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    # (a) num_cpg_observed
    ax = axes[0]
    sns.violinplot(data=inside, x="truth_status", y="num_cpg_observed", hue="truth_status",
                   palette=TRUTH_PAL, order=["TP", "FP"],
                   inner="quartile", cut=0, ax=ax, density_norm="width", legend=False)
    tp_cpg = inside.loc[inside["truth_status"] == "TP", "num_cpg_observed"].dropna().values
    fp_cpg = inside.loc[inside["truth_status"] == "FP", "num_cpg_observed"].dropna().values
    cpg_d = cohens_d(tp_cpg, fp_cpg)
    ax.set_title(f"num_cpg_observed (d={cpg_d:.3f})", fontsize=12, fontweight="bold")
    ax.set_xlabel("")

    # (b) reads per region
    ax = axes[1]
    reads_per_region = inside.groupby(["region_key", "truth_status"]).size().reset_index(name="reads_per_region")
    sns.violinplot(data=reads_per_region, x="truth_status", y="reads_per_region", hue="truth_status",
                   palette=TRUTH_PAL, order=["TP", "FP"],
                   inner="quartile", cut=0, ax=ax, density_norm="width", legend=False)
    tp_rpr = reads_per_region.loc[reads_per_region["truth_status"] == "TP", "reads_per_region"].values
    fp_rpr = reads_per_region.loc[reads_per_region["truth_status"] == "FP", "reads_per_region"].values
    rpr_d = cohens_d(tp_rpr, fp_rpr)
    ax.set_title(f"Reads per Region (d={rpr_d:.3f})", fontsize=12, fontweight="bold")
    ax.set_xlabel("")

    fig.suptitle("Fig09: Confound Check — LOH Inside", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "loh_methyl_level_fig09_confound_check.png")
    plt.close(fig)
    print("  -> fig09 saved")

    # ---- Part 5: Potential_LOH distribution check ----
    print("\n  Potential_LOH distribution inside LOH by truth_status:")
    pot_loh = inside.groupby(["truth_status", "Potential_LOH"]).size().unstack(fill_value=0)
    print(pot_loh)

    # ---- Fig 10: Summary Table ----
    print("\n[7/7] Part 5 — Summary table ...")
    summary_rows = []
    for m in all_metrics:
        if "|" not in m["label"]:
            continue
        parts = m["label"].split("|")
        feat = parts[0]
        zone = parts[1]
        summary_rows.append({
            "feature": feat,
            "zone": zone,
            "n_tp": m["n_tp"],
            "n_fp": m["n_fp"],
            "tp_median": m["tp_median"],
            "fp_median": m["fp_median"],
            "auc": m["auc"],
            "cohens_d": m["cohens_d"],
            "d_ci_lo": m["d_ci_lo"],
            "d_ci_hi": m["d_ci_hi"],
        })
    summary_df = pd.DataFrame(summary_rows)

    # Save stats TSV
    summary_df.to_csv(DATA_DIR / "loh_methyl_level_stats.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"  Stats saved to {DATA_DIR / 'loh_methyl_level_stats.tsv'}")

    # Render as figure
    fig, ax = plt.subplots(figsize=(20, max(6, 0.5 * len(summary_df) + 2)))
    ax.axis("off")
    # Format for display
    disp = summary_df.copy()
    for c in ["tp_median", "fp_median", "auc", "cohens_d", "d_ci_lo", "d_ci_hi"]:
        disp[c] = disp[c].apply(lambda v: f"{v:.3f}" if pd.notnull(v) else "N/A")
    disp["d_ci"] = disp["d_ci_lo"] + " .. " + disp["d_ci_hi"]
    disp = disp.drop(columns=["d_ci_lo", "d_ci_hi"])
    header = disp.columns.tolist()
    rows = disp.values.tolist()
    tbl = ax.table(cellText=rows, colLabels=header, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.4)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#BBDEFB")
            cell.set_text_props(fontweight="bold")
        # Highlight inside LOH rows
        if row > 0 and len(rows[row-1]) > 1 and rows[row-1][1] == "inside":
            cell.set_facecolor("#FFF9C4")
    ax.set_title("Fig10: Complete Metrics Summary (yellow = inside LOH)", fontsize=14, fontweight="bold", pad=20)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "loh_methyl_level_fig10_summary_table.png")
    plt.close(fig)
    print("  -> fig10 saved")

    # ---- Generate report ----
    read_in = [x for x in all_metrics if x["label"] == "methyl_mean|inside"][0]
    read_out = [x for x in all_metrics if x["label"] == "methyl_mean|outside"][0]
    reg_in = region_agg_metrics["inside"]
    reg_out = region_agg_metrics["outside"]

    ps_lines = []
    for _, row in ps_df.iterrows():
        auc_str = f"{row['auc']:.3f}" if pd.notnull(row['auc']) else "N/A"
        d_str = f"{row['cohens_d']:.3f}" if pd.notnull(row['cohens_d']) else "N/A"
        ps_lines.append(f"| {row['sample']:20s} | {auc_str:>8s} | {d_str:>8s} | {int(row['n_tp']):>6d} | {int(row['n_fp']):>6d} |")

    report = dedent(f"""\
    <!--
    建立時間: 2026-04-02
    目標: 分析 TP/FP 甲基化水平差異在 LOH 區域內是否持續存在
    處理範圍: 86,521 reads x 7 samples, paired-only, LOH zone overlay from TO LOH.bed
    關聯檔案:
      - research/loh_investigation/figures/loh_methyl_level_fig*.png
      - research/loh_investigation/data/loh_methyl_level_stats.tsv
    -->

    # LOH 區域內甲基化水平 TP/FP 差異分析

    ## 背景

    O10 觀察發現 read-level methyl_mean 在 TP (0.463) vs FP (0.679) 之間有顯著差異 (AUC=0.733)。
    本分析的核心問題：**這個差異在 LOH 區域內是否仍然存在？**

    ## 數據

    - 來源: phase1_shard_read_training_table.tsv (paired-only, 86,521 reads)
    - LOH 區域: 7 samples 的 TO-mode LOH.bed
    - LOH 判定: 以 read 的 chr + start 座標是否落在 LOH.bed 區間內

    ## 核心結果

    ### LOH Zone Distribution

    | Zone    | n_reads | % of total |
    |---------|---------|------------|
    | inside  | {zone_counts.get("inside", 0):>7,} | {zone_counts.get("inside", 0)/len(df)*100:.1f}% |
    | outside | {zone_counts.get("outside", 0):>7,} | {zone_counts.get("outside", 0)/len(df)*100:.1f}% |

    ### methyl_mean AUC (TP=1)

    | Level        | Inside LOH | Outside LOH | O10 Baseline |
    |--------------|------------|-------------|--------------|
    | Read-level   | **{read_in['auc']:.3f}**    | {read_out['auc']:.3f}       | 0.733        |
    | Region-level | **{reg_in['auc']:.3f}**    | {reg_out['auc']:.3f}       | —            |

    ### Cohen's d (TP - FP)

    | Level        | Inside LOH | Outside LOH |
    |--------------|------------|-------------|
    | Read-level   | {read_in['cohens_d']:.3f} [{read_in['d_ci_lo']:.3f}, {read_in['d_ci_hi']:.3f}] | {read_out['cohens_d']:.3f} [{read_out['d_ci_lo']:.3f}, {read_out['d_ci_hi']:.3f}] |
    | Region-level | {reg_in['cohens_d']:.3f} [{reg_in['d_ci_lo']:.3f}, {reg_in['d_ci_hi']:.3f}] | {reg_out['cohens_d']:.3f} [{reg_out['d_ci_lo']:.3f}, {reg_out['d_ci_hi']:.3f}] |

    ### Per-Sample Inside LOH

    | Sample               |      AUC | Cohen's d | n_TP   | n_FP   |
    |----------------------|----------|-----------|--------|--------|
    {chr(10).join(ps_lines)}

    ## 圖表清單

    1. fig01 — Zone distribution by sample
    2. fig02 — **Core figure**: methyl_mean violin by zone x truth_status
    3. fig03 — All methyl features inside LOH
    4. fig04 — AUC inside vs outside comparison
    5. fig05 — Per-sample violin inside LOH
    6. fig06 — Region-level vs read-level comparison
    7. fig07 — ALT/REF support inside LOH
    8. fig08 — HP tag inside LOH
    9. fig09 — Confound check (num_cpg, reads per region)
    10. fig10 — Complete summary table

    ## 結論

    **Inside LOH read-level AUC = {read_in['auc']:.3f}** (region-level = {reg_in['auc']:.3f})
    **Outside LOH read-level AUC = {read_out['auc']:.3f}** (region-level = {reg_out['auc']:.3f})

    """)

    # Add interpretation
    if read_in['auc'] >= 0.70:
        report += "**判定: methyl_mean TP/FP 差異在 LOH 區域內 *持續存在*，AUC 維持在高水平。** 甲基化水平是一個穩健的 discriminator，即使在 LOH 區域內也有效。\n"
    elif read_in['auc'] >= 0.60:
        report += "**判定: methyl_mean TP/FP 差異在 LOH 區域內 *部分衰減*，AUC 下降但仍有中等 discriminability。** 甲基化水平在 LOH 區域內仍有一定價值，但效力降低。\n"
    else:
        report += "**判定: methyl_mean TP/FP 差異在 LOH 區域內 *消失或大幅衰減*。** LOH 區域的特殊性質可能使甲基化水平不再具有區分 TP/FP 的能力。\n"

    with open(RPT_DIR / "20260402_loh_methyl_level_inside_loh.md", "w") as f:
        f.write(report)
    print(f"\n  Report saved to {RPT_DIR / '20260402_loh_methyl_level_inside_loh.md'}")

    # ---- Final Summary ----
    print("\n" + "=" * 70)
    print("FINAL RESULTS SUMMARY")
    print("=" * 70)
    print(f"\nTotal reads: {len(df):,}")
    print(f"  Inside LOH:  {zone_counts.get('inside',0):,} ({zone_counts.get('inside',0)/len(df)*100:.1f}%)")
    print(f"  Outside LOH: {zone_counts.get('outside',0):,} ({zone_counts.get('outside',0)/len(df)*100:.1f}%)")
    print(f"\nmethyl_mean AUC (TP=1):")
    print(f"  Inside LOH  — Read-level: {read_in['auc']:.4f}  Region-level: {reg_in['auc']:.4f}")
    print(f"  Outside LOH — Read-level: {read_out['auc']:.4f}  Region-level: {reg_out['auc']:.4f}")
    print(f"  O10 baseline: 0.733")
    print(f"\nCohen's d (TP - FP methyl_mean):")
    print(f"  Inside LOH:  {read_in['cohens_d']:.3f}")
    print(f"  Outside LOH: {read_out['cohens_d']:.3f}")
    print(f"\nPer-sample AUC inside LOH:")
    for _, row in ps_df.iterrows():
        auc_str = f"{row['auc']:.3f}" if pd.notnull(row['auc']) else "N/A"
        print(f"  {row['sample']:20s}: {auc_str}")
    print("\nDone.")


if __name__ == "__main__":
    main()
