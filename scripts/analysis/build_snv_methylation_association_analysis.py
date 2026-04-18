#!/usr/bin/env python3
"""
SNV-Methylation Association Analysis: Multi-Method ASM Detection

This script quantifies what proportion of SNV sites show significant
allele-specific methylation (ASM) using 5 independent methods,
compares ISM results with independent methods, and analyzes paired vs TO modes.

Methods:
  M1: ISM PERMANOVA (HPMergedSig, AlleleSig) — already computed
  M2: ALT vs REF mean methylation Wilcoxon test
  M3: Per-CpG Fisher's exact test (binarized methylation × allele)
  M4: Point-biserial correlation (allele × mean methylation)
  M5: Bootstrap delta confidence interval

Output:
  - per_region_multi_method.tsv
  - concordance_matrix.tsv
  - proportion_summary.tsv
  - 8+ figures
  - 20260401_SNV甲基化關聯性定量分析報告_01.md
"""

import sys
import os
import warnings
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from functools import partial

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu, fisher_exact, pointbiserialr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MASTER_DATASET_PATH = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260327_loh_round1_cross_sample_audit/"
    "all_region_rows.tsv.gz"
)

CANONICAL_ROOT = Path("/big7_disk/liaoyoyo2001/big7_disk_output/canonical")

OUTPUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260401_snv_methylation_association"
)

SAMPLE_ORDER = [
    "HCC1395", "HCC1395_DORADO", "COLO829",
    "H1437", "H2009", "HCC1937", "HCC1954",
]

# Run config for paired mode (sample -> run path)
PAIRED_RUN_PATHS = {
    "HCC1395": CANONICAL_ROOT / "HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix",
    "HCC1395_DORADO": CANONICAL_ROOT / "HCC1395_DORADO/paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix",
    "COLO829": CANONICAL_ROOT / "COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix",
    "H1437": CANONICAL_ROOT / "H1437/paired_full/20260315_H1437_paired_full_full_complete_matrix",
    "H2009": CANONICAL_ROOT / "H2009/paired_full/20260315_H2009_paired_full_full_complete_matrix",
    "HCC1937": CANONICAL_ROOT / "HCC1937/paired_full/20260315_HCC1937_paired_full_full_complete_matrix",
    "HCC1954": CANONICAL_ROOT / "HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix",
}

# TO mode: no complete_matrix runs with read-level data available
# TO analysis relies on master dataset M1 statistics only
TO_RUN_PATHS = {}

N_BOOTSTRAP = 1000
ALPHA = 0.05
MIN_READS_PER_ALLELE = 5
METHYL_BINARIZE_THRESHOLD = 0.5
MAX_REGIONS_PER_SAMPLE = 500  # subsample for speed; set 0 for unlimited
N_WORKERS = max(1, cpu_count() - 2)  # use all but 2 cores

np.random.seed(42)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def load_region_data(region_dir: Path):
    """Load reads.tsv and methylation.csv for a region."""
    reads_path = region_dir / "reads" / "reads.tsv"
    methyl_path = region_dir / "methylation" / "methylation.csv"

    if not reads_path.exists() or not methyl_path.exists():
        return None, None

    reads = pd.read_csv(reads_path, sep="\t")
    methyl = pd.read_csv(methyl_path, index_col=0)

    return reads, methyl


def compute_read_mean_methylation(methyl_df):
    """Compute per-read mean methylation from the CpG matrix."""
    # Convert to numeric, coerce errors
    numeric = methyl_df.apply(pd.to_numeric, errors="coerce")
    return numeric.mean(axis=1)


def method_m2_wilcoxon(alt_methyl, ref_methyl):
    """M2: Wilcoxon rank-sum test on mean methylation ALT vs REF."""
    if len(alt_methyl) < MIN_READS_PER_ALLELE or len(ref_methyl) < MIN_READS_PER_ALLELE:
        return np.nan, np.nan, False
    try:
        stat, pval = mannwhitneyu(alt_methyl, ref_methyl, alternative="two-sided")
    except ValueError:
        return np.nan, np.nan, False
    delta = np.median(alt_methyl) - np.median(ref_methyl)
    return delta, pval, pval < ALPHA


def method_m3_per_cpg_fisher(methyl_df, allele_labels):
    """M3: Per-CpG Fisher's exact test (binarized methylation × allele).

    Returns: n_significant_cpgs, n_tested_cpgs, max_odds_ratio, min_pval
    """
    n_sig = 0
    n_tested = 0
    min_p = 1.0
    max_or = 1.0

    alt_mask = allele_labels == "ALT"
    ref_mask = allele_labels == "REF"

    for col in methyl_df.columns:
        vals = pd.to_numeric(methyl_df[col], errors="coerce")
        valid = vals.dropna()
        if len(valid) < 10:
            continue

        # Binarize
        binary = (vals > METHYL_BINARIZE_THRESHOLD).astype(int)
        valid_idx = binary.dropna().index

        a_idx = valid_idx[valid_idx.isin(alt_mask[alt_mask].index)]
        r_idx = valid_idx[valid_idx.isin(ref_mask[ref_mask].index)]

        if len(a_idx) < 3 or len(r_idx) < 3:
            continue

        # 2x2 table: allele (ALT/REF) × methylation (high/low)
        alt_high = binary.loc[a_idx].sum()
        alt_low = len(a_idx) - alt_high
        ref_high = binary.loc[r_idx].sum()
        ref_low = len(r_idx) - ref_high

        table = [[alt_high, alt_low], [ref_high, ref_low]]
        try:
            odds_ratio, pval = fisher_exact(table)
        except ValueError:
            continue

        n_tested += 1
        if pval < ALPHA:
            n_sig += 1
        if pval < min_p:
            min_p = pval
            max_or = odds_ratio

    frac_sig = n_sig / n_tested if n_tested > 0 else 0.0
    return n_sig, n_tested, frac_sig, max_or, min_p


def method_m4_point_biserial(mean_methyl, allele_binary):
    """M4: Point-biserial correlation between allele (0/1) and mean methylation."""
    valid = ~(np.isnan(mean_methyl) | np.isnan(allele_binary))
    if valid.sum() < 10:
        return np.nan, np.nan, False
    try:
        r, pval = pointbiserialr(allele_binary[valid], mean_methyl[valid])
    except Exception:
        return np.nan, np.nan, False
    return r, pval, pval < ALPHA


def method_m5_bootstrap_delta(alt_methyl, ref_methyl, n_boot=N_BOOTSTRAP):
    """M5: Bootstrap CI for mean methylation delta (ALT - REF)."""
    if len(alt_methyl) < MIN_READS_PER_ALLELE or len(ref_methyl) < MIN_READS_PER_ALLELE:
        return np.nan, np.nan, np.nan, False

    observed_delta = np.mean(alt_methyl) - np.mean(ref_methyl)
    boot_deltas = np.empty(n_boot)
    for i in range(n_boot):
        a_boot = np.random.choice(alt_methyl, size=len(alt_methyl), replace=True)
        r_boot = np.random.choice(ref_methyl, size=len(ref_methyl), replace=True)
        boot_deltas[i] = np.mean(a_boot) - np.mean(r_boot)

    ci_lo = np.percentile(boot_deltas, 2.5)
    ci_hi = np.percentile(boot_deltas, 97.5)
    sig = not (ci_lo <= 0 <= ci_hi)  # CI doesn't contain 0
    return observed_delta, ci_lo, ci_hi, sig


def analyze_single_region(region_dir: Path):
    """Run all independent methods (M2-M5) on a single region.

    Returns a dict of results or None if insufficient data.
    """
    reads, methyl = load_region_data(region_dir)
    if reads is None or methyl is None:
        return None

    # Align reads and methylation
    if "alt_support" not in reads.columns:
        return None

    # Filter to ALT/REF reads only
    reads_filtered = reads[reads["alt_support"].isin(["ALT", "REF"])].copy()
    if len(reads_filtered) < 2 * MIN_READS_PER_ALLELE:
        return None

    # Compute per-read mean methylation
    common_ids = reads_filtered["read_id"].values
    methyl_sub = methyl.loc[methyl.index.isin(common_ids)]
    if len(methyl_sub) < 2 * MIN_READS_PER_ALLELE:
        return None

    mean_methyl = compute_read_mean_methylation(methyl_sub)

    # Align
    reads_aligned = reads_filtered.set_index("read_id").loc[methyl_sub.index]
    allele_labels = reads_aligned["alt_support"]

    alt_mean = mean_methyl[allele_labels == "ALT"].dropna().values
    ref_mean = mean_methyl[allele_labels == "REF"].dropna().values

    if len(alt_mean) < MIN_READS_PER_ALLELE or len(ref_mean) < MIN_READS_PER_ALLELE:
        return None

    # M2: Wilcoxon
    m2_delta, m2_p, m2_sig = method_m2_wilcoxon(alt_mean, ref_mean)

    # M3: Per-CpG Fisher
    m3_n_sig, m3_n_tested, m3_frac, m3_or, m3_min_p = method_m3_per_cpg_fisher(
        methyl_sub, allele_labels
    )

    # M4: Point-biserial
    allele_binary = (allele_labels == "ALT").astype(float).values
    m4_r, m4_p, m4_sig = method_m4_point_biserial(mean_methyl.values, allele_binary)

    # M5: Bootstrap delta
    m5_delta, m5_ci_lo, m5_ci_hi, m5_sig = method_m5_bootstrap_delta(alt_mean, ref_mean)

    return {
        "n_alt": len(alt_mean),
        "n_ref": len(ref_mean),
        "n_total": len(alt_mean) + len(ref_mean),
        "alt_mean_methyl": np.mean(alt_mean),
        "ref_mean_methyl": np.mean(ref_mean),
        "n_cpgs": methyl_sub.shape[1],
        # M2
        "m2_delta": m2_delta,
        "m2_pval": m2_p,
        "m2_sig": m2_sig,
        # M3
        "m3_n_sig_cpgs": m3_n_sig,
        "m3_n_tested_cpgs": m3_n_tested,
        "m3_frac_sig": m3_frac,
        "m3_max_or": m3_or,
        "m3_min_pval": m3_min_p,
        "m3_sig": m3_frac > 0.1,  # >10% of CpGs significant
        # M4
        "m4_r": m4_r,
        "m4_pval": m4_p,
        "m4_sig": m4_sig,
        # M5
        "m5_delta": m5_delta,
        "m5_ci_lo": m5_ci_lo,
        "m5_ci_hi": m5_ci_hi,
        "m5_sig": m5_sig,
    }


def find_region_dirs(run_path: Path, label: str, max_n: int = MAX_REGIONS_PER_SAMPLE):
    """Find region directories under a run path for given label (tp/fp)."""
    base = run_path / f"intersubmod_{label}" / f"filtered_snv_{label}"
    if not base.exists():
        return []

    dirs = []
    for chr_dir in sorted(base.iterdir()):
        if not chr_dir.is_dir():
            continue
        for snv_dir in sorted(chr_dir.iterdir()):
            if not snv_dir.is_dir():
                continue
            for region_dir in sorted(snv_dir.iterdir()):
                if region_dir.is_dir() and (region_dir / "reads").exists():
                    dirs.append(region_dir)

    # Subsample if too many
    if len(dirs) > max_n:
        rng = np.random.RandomState(42)
        dirs = list(rng.choice(dirs, size=max_n, replace=False))

    return dirs


def extract_region_key(region_dir: Path):
    """Extract chr_pos identifier from region directory path."""
    # path like: .../chr1/chr1_877772/chr1_875772_879772
    parts = region_dir.parts
    snv_dir = parts[-2]  # chr1_877772
    return snv_dir


def _worker_analyze(args):
    """Worker function for multiprocessing: analyze one region."""
    rdir_str, sample, mode_name, label = args
    rdir = Path(rdir_str)
    try:
        result = analyze_single_region(rdir)
        if result is None:
            return None
        result["sample"] = sample
        result["mode"] = mode_name
        result["truth_label"] = label.upper()
        result["region_key"] = extract_region_key(rdir)
        return result
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Main Analysis
# ---------------------------------------------------------------------------


def run_multi_method_analysis():
    """Run the full multi-method ASM analysis."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Phase 1: ISM data from master dataset ---
    print("=" * 60)
    print("Phase 1: Loading ISM master dataset for M1 statistics")
    print("=" * 60)

    master = pd.read_csv(
        MASTER_DATASET_PATH, sep="\t", compression="gzip",
        usecols=[
            "RegionID", "Chr", "Pos", "NumReads", "NumCpGs",
            "HPMergedDelta", "HPMergedP", "HPMergedSig",
            "AlleleDelta", "AlleleP", "AlleleSig",
            "HPFineSig", "HPFineP", "GlobalP",
            "ClusterPermanovaP", "LabelAllelePermanovaP",
            "Significant", "PassedGating",
            "sample", "mode", "truth_label",
            "Potential_LOH", "VerificationClass",
        ],
        low_memory=False,
    )
    print(f"  Master: {len(master):,} rows")

    # M1 proportion table
    m1_results = []
    for mode in ["paired", "to"]:
        for label in ["TP", "FP"]:
            sub = master[(master["mode"] == mode) & (master["truth_label"] == label)]
            n = len(sub)
            if n == 0:
                continue

            gated = sub[sub["PassedGating"] == True]
            ng = len(gated)

            row = {
                "mode": mode,
                "truth_label": label,
                "n_total": n,
                "HPMergedSig_rate": sub["HPMergedSig"].mean(),
                "AlleleSig_rate": sub["AlleleSig"].mean(),
                "HPFineSig_rate": sub["HPFineSig"].mean(),
                "Significant_rate": sub["Significant"].mean(),
                "GlobalP_lt005_rate": (sub["GlobalP"] < 0.05).mean(),
                "stringent_ASM_rate": ((sub["HPMergedSig"] == True) & (sub["HPMergedDelta"].abs() > 0.1)).mean(),
                "strong_ASM_rate": ((sub["HPMergedSig"] == True) & (sub["HPMergedDelta"].abs() > 0.2)).mean(),
                "n_gated": ng,
                "gated_HPMergedSig_rate": gated["HPMergedSig"].mean() if ng > 0 else 0,
            }
            m1_results.append(row)

    m1_df = pd.DataFrame(m1_results)
    m1_df.to_csv(OUTPUT_DIR / "m1_ism_proportion_summary.tsv", sep="\t", index=False)
    print("\n  M1 ISM proportions saved.")
    print(m1_df.to_string(index=False))

    # Per-sample M1
    m1_per_sample = []
    for mode in ["paired", "to"]:
        for sample in SAMPLE_ORDER:
            for label in ["TP", "FP"]:
                sub = master[
                    (master["mode"] == mode)
                    & (master["sample"] == sample)
                    & (master["truth_label"] == label)
                ]
                if len(sub) < 10:
                    continue
                m1_per_sample.append({
                    "sample": sample,
                    "mode": mode,
                    "truth_label": label,
                    "n": len(sub),
                    "HPMergedSig_rate": sub["HPMergedSig"].mean(),
                    "AlleleSig_rate": sub["AlleleSig"].mean(),
                    "stringent_ASM": ((sub["HPMergedSig"] == True) & (sub["HPMergedDelta"].abs() > 0.1)).mean(),
                    "median_HPMergedDelta": sub["HPMergedDelta"].median(),
                    "median_AlleleDelta": sub["AlleleDelta"].median(),
                })
    m1_sample_df = pd.DataFrame(m1_per_sample)
    m1_sample_df.to_csv(OUTPUT_DIR / "m1_per_sample.tsv", sep="\t", index=False)

    # --- Phase 2: Independent methods (M2-M5) on read-level data ---
    print("\n" + "=" * 60)
    print("Phase 2: Independent methods (M2-M5) on read-level data")
    print("=" * 60)

    # Collect all region tasks
    all_tasks = []
    for mode_name, run_paths in [("paired", PAIRED_RUN_PATHS), ("to", TO_RUN_PATHS)]:
        for sample, run_path in run_paths.items():
            if not run_path.exists():
                print(f"  [SKIP] {sample} {mode_name}: run path not found")
                continue

            for label in ["tp", "fp"]:
                region_dirs = find_region_dirs(run_path, label)
                print(f"  [{mode_name}] {sample} {label.upper()}: {len(region_dirs)} regions")
                for rdir in region_dirs:
                    all_tasks.append((str(rdir), sample, mode_name, label))

    print(f"\n  Total tasks: {len(all_tasks)}, using {N_WORKERS} workers")

    # Run in parallel
    all_results = []
    with Pool(N_WORKERS) as pool:
        for i, result in enumerate(pool.imap_unordered(_worker_analyze, all_tasks, chunksize=20)):
            if result is not None:
                all_results.append(result)
            if (i + 1) % 500 == 0:
                print(f"    -> {i+1}/{len(all_tasks)} tasks done, {len(all_results)} valid results")

    print(f"    -> Final: {len(all_results)} valid results from {len(all_tasks)} tasks")

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUTPUT_DIR / "per_region_multi_method.tsv", sep="\t", index=False)
    print(f"\n  Total regions analyzed: {len(results_df):,}")

    # --- Phase 3: Concordance analysis ---
    print("\n" + "=" * 60)
    print("Phase 3: Concordance analysis (ISM vs Independent methods)")
    print("=" * 60)

    if len(results_df) == 0:
        print("  No results to analyze!")
        return

    # Compute proportions by method
    method_cols = ["m2_sig", "m3_sig", "m4_sig", "m5_sig"]
    method_names = ["M2:Wilcoxon", "M3:Fisher", "M4:PtBiserial", "M5:Bootstrap"]

    proportion_rows = []
    for mode in ["paired", "to"]:
        for label in ["TP", "FP"]:
            sub = results_df[(results_df["mode"] == mode) & (results_df["truth_label"] == label)]
            if len(sub) < 10:
                continue
            row = {"mode": mode, "truth_label": label, "n": len(sub)}
            for mc, mn in zip(method_cols, method_names):
                row[mn] = sub[mc].mean()
            # Add combined: significant in >=2 methods
            combined = (sub[method_cols].sum(axis=1) >= 2)
            row["Combined_2of4"] = combined.mean()
            combined3 = (sub[method_cols].sum(axis=1) >= 3)
            row["Combined_3of4"] = combined3.mean()
            proportion_rows.append(row)

    prop_df = pd.DataFrame(proportion_rows)
    prop_df.to_csv(OUTPUT_DIR / "proportion_summary_independent.tsv", sep="\t", index=False)
    print("\n  Independent method proportions:")
    print(prop_df.to_string(index=False))

    # Per-sample proportions for independent methods
    per_sample_indep = []
    for sample in SAMPLE_ORDER:
        for label in ["TP", "FP"]:
            sub = results_df[(results_df["sample"] == sample) &
                             (results_df["truth_label"] == label) &
                             (results_df["mode"] == "paired")]
            if len(sub) < 5:
                continue
            row = {"sample": sample, "truth_label": label, "n": len(sub)}
            for mc, mn in zip(method_cols, method_names):
                row[mn] = sub[mc].mean()
            delta = sub["m2_delta"].dropna().abs()
            row["median_abs_delta"] = delta.median() if len(delta) > 0 else np.nan
            row["stringent_m2"] = ((sub["m2_sig"] == True) & (sub["m2_delta"].abs() > 0.1)).mean()
            per_sample_indep.append(row)
    if per_sample_indep:
        ps_df = pd.DataFrame(per_sample_indep)
        ps_df.to_csv(OUTPUT_DIR / "per_sample_independent_proportions.tsv", sep="\t", index=False)
        print("\n  Per-sample independent proportions:")
        print(ps_df.to_string(index=False))

    # Method-method concordance matrix
    print("\n  Method concordance matrix:")
    for mode in ["paired", "to"]:
        sub = results_df[results_df["mode"] == mode]
        if len(sub) < 20:
            continue
        print(f"\n  --- {mode} mode (n={len(sub)}) ---")
        for i, (mc1, mn1) in enumerate(zip(method_cols, method_names)):
            for j, (mc2, mn2) in enumerate(zip(method_cols, method_names)):
                if j <= i:
                    continue
                both_sig = ((sub[mc1] == True) & (sub[mc2] == True)).sum()
                either_sig = ((sub[mc1] == True) | (sub[mc2] == True)).sum()
                jaccard = both_sig / either_sig if either_sig > 0 else 0
                print(f"    {mn1} vs {mn2}: Jaccard={jaccard:.3f}, both={both_sig}, either={either_sig}")

    # --- Phase 4: Effect size analysis ---
    print("\n" + "=" * 60)
    print("Phase 4: Effect size distribution")
    print("=" * 60)

    for mode in ["paired", "to"]:
        for label in ["TP", "FP"]:
            sub = results_df[(results_df["mode"] == mode) & (results_df["truth_label"] == label)]
            if len(sub) < 10:
                continue
            delta = sub["m2_delta"].dropna()
            abs_delta = delta.abs()
            print(f"\n  [{mode}] {label} (n={len(sub)}): |delta| median={abs_delta.median():.4f}, "
                  f">0.05: {(abs_delta>0.05).mean()*100:.1f}%, "
                  f">0.1: {(abs_delta>0.1).mean()*100:.1f}%, "
                  f">0.2: {(abs_delta>0.2).mean()*100:.1f}%")

    # --- Phase 5: Figures ---
    print("\n" + "=" * 60)
    print("Phase 5: Generating figures")
    print("=" * 60)

    generate_figures(results_df, m1_df, m1_sample_df, OUTPUT_DIR)

    print("\n  Done! All outputs in:", OUTPUT_DIR)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def generate_figures(results_df, m1_df, m1_sample_df, out_dir):
    """Generate all figures."""
    sns.set_style("whitegrid")

    # Fig 1: ISM ASM proportions by mode and truth label
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax_idx, mode in enumerate(["paired", "to"]):
        ax = axes[ax_idx]
        sub = m1_df[m1_df["mode"] == mode]
        metrics = ["HPMergedSig_rate", "AlleleSig_rate", "stringent_ASM_rate", "strong_ASM_rate"]
        labels_m = ["HPMergedSig", "AlleleSig", "Stringent\n(Sig+|Δ|>0.1)", "Strong\n(Sig+|Δ|>0.2)"]
        x = np.arange(len(metrics))
        w = 0.35
        for i, label in enumerate(["TP", "FP"]):
            row = sub[sub["truth_label"] == label]
            if len(row) == 0:
                continue
            vals = [row.iloc[0][m] * 100 for m in metrics]
            color = "#2196F3" if label == "TP" else "#F44336"
            bars = ax.bar(x + (i - 0.5) * w, vals, w, label=label, color=color, alpha=0.8)
            for bar, val in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                        f"{val:.1f}%", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels_m, fontsize=9)
        ax.set_ylabel("Proportion (%)")
        ax.set_title(f"ISM ASM Indicators — {mode.upper()} mode")
        ax.legend()
        ax.set_ylim(0, max(60, ax.get_ylim()[1] * 1.15))
    plt.tight_layout()
    plt.savefig(out_dir / "fig01_ism_asm_proportions.png", dpi=150)
    plt.close()
    print("  fig01 saved")

    # Fig 2: Independent methods ASM proportions
    if len(results_df) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        method_cols = ["m2_sig", "m3_sig", "m4_sig", "m5_sig"]
        method_labels = ["M2:Wilcoxon", "M3:Fisher\n(>10% CpGs)", "M4:PtBiserial", "M5:Bootstrap"]

        for ax_idx, mode in enumerate(["paired", "to"]):
            ax = axes[ax_idx]
            x = np.arange(len(method_cols))
            w = 0.35
            for i, label in enumerate(["TP", "FP"]):
                sub = results_df[(results_df["mode"] == mode) & (results_df["truth_label"] == label)]
                if len(sub) < 5:
                    continue
                vals = [sub[mc].mean() * 100 for mc in method_cols]
                color = "#2196F3" if label == "TP" else "#F44336"
                bars = ax.bar(x + (i - 0.5) * w, vals, w, label=f"{label} (n={len(sub)})",
                              color=color, alpha=0.8)
                for bar, val in zip(bars, vals):
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                            f"{val:.1f}%", ha="center", va="bottom", fontsize=8)
            ax.set_xticks(x)
            ax.set_xticklabels(method_labels, fontsize=9)
            ax.set_ylabel("Proportion with significant ASM (%)")
            ax.set_title(f"Independent Methods — {mode.upper()} mode")
            ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "fig02_independent_methods_proportions.png", dpi=150)
        plt.close()
        print("  fig02 saved")

    # Fig 3: Effect size (delta) distribution
    if len(results_df) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        for ax_idx, mode in enumerate(["paired", "to"]):
            ax = axes[ax_idx]
            for label, color in [("TP", "#2196F3"), ("FP", "#F44336")]:
                sub = results_df[(results_df["mode"] == mode) & (results_df["truth_label"] == label)]
                delta = sub["m2_delta"].dropna()
                if len(delta) > 0:
                    ax.hist(delta, bins=50, alpha=0.5, color=color, label=f"{label} (n={len(delta)})",
                            density=True)
            ax.axvline(0, color="black", ls="--", alpha=0.5)
            ax.axvline(-0.1, color="gray", ls=":", alpha=0.5)
            ax.axvline(0.1, color="gray", ls=":", alpha=0.5)
            ax.set_xlabel("ALT - REF mean methylation delta")
            ax.set_ylabel("Density")
            ax.set_title(f"Effect Size Distribution — {mode.upper()}")
            ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "fig03_delta_distribution.png", dpi=150)
        plt.close()
        print("  fig03 saved")

    # Fig 4: Per-sample ASM rates (ISM)
    if len(m1_sample_df) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        for ax_idx, mode in enumerate(["paired", "to"]):
            ax = axes[ax_idx]
            sub = m1_sample_df[m1_sample_df["mode"] == mode]
            for label, color in [("TP", "#2196F3"), ("FP", "#F44336")]:
                s = sub[sub["truth_label"] == label]
                if len(s) == 0:
                    continue
                ax.barh(
                    np.arange(len(s)) + (0.2 if label == "FP" else -0.2),
                    s["HPMergedSig_rate"] * 100,
                    height=0.35, color=color, alpha=0.8, label=label,
                )
            samples_in = sub[sub["truth_label"] == "TP"]["sample"].values
            ax.set_yticks(np.arange(len(samples_in)))
            ax.set_yticklabels(samples_in, fontsize=9)
            ax.set_xlabel("HPMergedSig Rate (%)")
            ax.set_title(f"Per-Sample ISM ASM Rate — {mode.upper()}")
            ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "fig04_per_sample_ism_asm.png", dpi=150)
        plt.close()
        print("  fig04 saved")

    # Fig 5: Method concordance heatmap
    if len(results_df) > 20:
        method_cols = ["m2_sig", "m3_sig", "m4_sig", "m5_sig"]
        method_names = ["M2:Wilcoxon", "M3:Fisher", "M4:PtBiserial", "M5:Bootstrap"]

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax_idx, mode in enumerate(["paired", "to"]):
            ax = axes[ax_idx]
            sub = results_df[results_df["mode"] == mode]
            if len(sub) < 10:
                continue

            n = len(method_cols)
            concordance = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    both = ((sub[method_cols[i]] == True) & (sub[method_cols[j]] == True)).sum()
                    either = ((sub[method_cols[i]] == True) | (sub[method_cols[j]] == True)).sum()
                    concordance[i, j] = both / either if either > 0 else 0

            sns.heatmap(concordance, annot=True, fmt=".2f", cmap="YlOrRd",
                        xticklabels=method_names, yticklabels=method_names,
                        vmin=0, vmax=1, ax=ax)
            ax.set_title(f"Method Concordance (Jaccard) — {mode.upper()}")
        plt.tight_layout()
        plt.savefig(out_dir / "fig05_method_concordance.png", dpi=150)
        plt.close()
        print("  fig05 saved")

    # Fig 6: Stringent ASM by truth label and LOH status (ISM)
    # (from master data, not subsampled)

    # Fig 7: CpG Fisher sig fraction vs Wilcoxon delta
    if len(results_df) > 20:
        fig, ax = plt.subplots(figsize=(8, 6))
        for label, color, marker in [("TP", "#2196F3", "o"), ("FP", "#F44336", "s")]:
            sub = results_df[results_df["truth_label"] == label]
            ax.scatter(sub["m2_delta"], sub["m3_frac_sig"],
                       c=color, alpha=0.3, s=15, marker=marker, label=label)
        ax.axhline(0.1, color="gray", ls=":", alpha=0.5, label="M3 threshold")
        ax.axvline(-0.1, color="gray", ls="--", alpha=0.3)
        ax.axvline(0.1, color="gray", ls="--", alpha=0.3)
        ax.set_xlabel("M2: ALT-REF mean methylation delta")
        ax.set_ylabel("M3: Fraction of CpGs with significant Fisher test")
        ax.set_title("Method Cross-Validation: M2 effect size vs M3 CpG-level significance")
        ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "fig06_m2_vs_m3_scatter.png", dpi=150)
        plt.close()
        print("  fig06 saved")

    # Fig 7: Number-of-methods agreement distribution
    if len(results_df) > 20:
        method_cols = ["m2_sig", "m3_sig", "m4_sig", "m5_sig"]
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax_idx, mode in enumerate(["paired", "to"]):
            ax = axes[ax_idx]
            sub = results_df[results_df["mode"] == mode]
            if len(sub) < 10:
                continue
            n_methods_sig = sub[method_cols].sum(axis=1)
            for label, color in [("TP", "#2196F3"), ("FP", "#F44336")]:
                s = n_methods_sig[sub["truth_label"] == label]
                counts = s.value_counts().sort_index()
                total = len(s)
                ax.bar(counts.index + (-0.15 if label == "TP" else 0.15),
                       counts.values / total * 100, width=0.3,
                       color=color, alpha=0.8, label=label)
            ax.set_xlabel("Number of methods agreeing (out of 4)")
            ax.set_ylabel("Proportion (%)")
            ax.set_title(f"Multi-Method Agreement — {mode.upper()}")
            ax.set_xticks([0, 1, 2, 3, 4])
            ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "fig07_method_agreement_distribution.png", dpi=150)
        plt.close()
        print("  fig07 saved")

    # Fig 8: Per-sample independent method ASM rates (M2 Wilcoxon)
    if len(results_df) > 20:
        samples_in_data = [s for s in SAMPLE_ORDER if s in results_df["sample"].values]
        if len(samples_in_data) > 1:
            fig, ax = plt.subplots(figsize=(14, 7))
            x = np.arange(len(samples_in_data))
            w = 0.35
            for i, label in enumerate(["TP", "FP"]):
                rates = []
                ns = []
                for s in samples_in_data:
                    sub = results_df[(results_df["sample"] == s) &
                                     (results_df["truth_label"] == label) &
                                     (results_df["mode"] == "paired")]
                    rates.append(sub["m2_sig"].mean() * 100 if len(sub) >= 5 else 0)
                    ns.append(len(sub))
                color = "#2196F3" if label == "TP" else "#F44336"
                bars = ax.bar(x + (i - 0.5) * w, rates, w, label=label, color=color, alpha=0.8)
                for bar, val, n in zip(bars, rates, ns):
                    if n >= 5:
                        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                                f"{val:.1f}%\n(n={n})", ha="center", va="bottom", fontsize=7)
            ax.set_xticks(x)
            ax.set_xticklabels(samples_in_data, fontsize=9, rotation=30)
            ax.set_ylabel("M2 Wilcoxon ASM Rate (%)")
            ax.set_title("Per-Sample Independent ASM Rate (M2 Wilcoxon, Paired mode)")
            ax.legend()
            plt.tight_layout()
            plt.savefig(out_dir / "fig08_per_sample_m2_asm.png", dpi=150)
            plt.close()
            print("  fig08 saved")

    # Fig 9: Per-sample median |delta| comparison
    if len(results_df) > 20:
        samples_in_data = [s for s in SAMPLE_ORDER if s in results_df["sample"].values]
        if len(samples_in_data) > 1:
            fig, ax = plt.subplots(figsize=(14, 7))
            x = np.arange(len(samples_in_data))
            w = 0.35
            for i, label in enumerate(["TP", "FP"]):
                deltas = []
                for s in samples_in_data:
                    sub = results_df[(results_df["sample"] == s) &
                                     (results_df["truth_label"] == label) &
                                     (results_df["mode"] == "paired")]
                    d = sub["m2_delta"].dropna().abs()
                    deltas.append(d.median() if len(d) >= 5 else 0)
                color = "#2196F3" if label == "TP" else "#F44336"
                bars = ax.bar(x + (i - 0.5) * w, deltas, w, label=label, color=color, alpha=0.8)
                for bar, val in zip(bars, deltas):
                    if val > 0:
                        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.001,
                                f"{val:.3f}", ha="center", va="bottom", fontsize=8)
            ax.set_xticks(x)
            ax.set_xticklabels(samples_in_data, fontsize=9, rotation=30)
            ax.set_ylabel("Median |ALT-REF methylation delta|")
            ax.set_title("Per-Sample Effect Size (Paired mode)")
            ax.legend()
            plt.tight_layout()
            plt.savefig(out_dir / "fig09_per_sample_median_delta.png", dpi=150)
            plt.close()
            print("  fig09 saved")


# ---------------------------------------------------------------------------
# Entry
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    run_multi_method_analysis()
