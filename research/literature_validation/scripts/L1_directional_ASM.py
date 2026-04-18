#!/usr/bin/env python3
"""L1_directional_ASM.py — Test the epiTRACERx directional ASM hypothesis.

Literature prediction (Basilicata et al., Gaiti et al., epiTRACERx):
  At somatic SNV sites, the ALT allele should show systematically DIFFERENT
  methylation than the REF allele. Specifically, at promoter/CpG-island loci,
  ALT-carrying reads are predicted to be more methylated (gene silencing).

This script tests:
  1. Whether AlleleDelta (PERMANOVA) differs between TP and FP (Part A)
  2. Whether directional methylation (ALT_mean - REF_mean) is systematic at TP (Part B)
  3. Whether CpG-level variability (fCpG concept from EVOFLUx) differs TP vs FP (Part C)

Data sources:
  - significance_summary.csv: per-region ISM features (all regions, fast)
  - reads.tsv + methylation.csv: per-read raw data (sampled regions, slower)

Usage:
  python3 research/literature_validation/scripts/L1_directional_ASM.py
"""

import os
import sys
import glob
import random
import warnings
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================================
# Configuration
# ============================================================================
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
CANONICAL_ROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "research", "literature_validation")
FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures")
DATA_DIR = os.path.join(OUTPUT_DIR, "data")

SAMPLES = [
    "HCC1395", "COLO829", "H1437", "H2009",
    "HCC1937", "HCC1954", "HCC1395_DORADO",
]

# Run ID patterns (complete_matrix runs)
RUN_PATTERNS = {
    "HCC1395": "20260314_HCC1395_paired_pileup_pileup_complete_matrix",
    "COLO829": "20260315_COLO829_paired_pileup_pileup_complete_matrix",
    "H1437": "20260315_H1437_paired_pileup_pileup_complete_matrix",
    "H2009": "20260315_H2009_paired_pileup_pileup_complete_matrix",
    "HCC1937": "20260315_HCC1937_paired_pileup_pileup_complete_matrix",
    "HCC1954": "20260315_HCC1954_paired_pileup_pileup_complete_matrix",
    "HCC1395_DORADO": "20260315_HCC1395_DORADO_paired_pileup_pileup_complete_matrix",
}

N_SAMPLE_REGIONS = 500  # regions to sample for Part B per TP/FP per sample
RANDOM_SEED = 42


# ============================================================================
# Part A: Aggregate analysis from significance_summary.csv
# ============================================================================
def load_significance_summary(sample, category):
    """Load significance_summary.csv for a sample and category (tp/fp)."""
    run_id = RUN_PATTERNS.get(sample)
    if not run_id:
        return None

    path = os.path.join(
        CANONICAL_ROOT, sample, "paired_pileup", run_id,
        f"intersubmod_{category}", "significance_summary.csv"
    )
    if not os.path.isfile(path):
        print(f"  [WARN] Missing: {path}")
        return None

    try:
        df = pd.read_csv(path)
        df["sample"] = sample
        df["category"] = category.upper()
        return df
    except Exception as e:
        print(f"  [ERROR] Failed to load {path}: {e}")
        return None


def part_a_aggregate_analysis():
    """Analyze AlleleDelta and allele significance across TP/FP for all samples."""
    print("=" * 70)
    print("PART A: Aggregate Allele Analysis from significance_summary.csv")
    print("=" * 70)

    all_data = []
    for sample in SAMPLES:
        for cat in ["tp", "fp"]:
            df = load_significance_summary(sample, cat)
            if df is not None:
                all_data.append(df)
                print(f"  {sample}/{cat}: {len(df)} regions", flush=True)

    if not all_data:
        print("[ERROR] No data loaded")
        return None

    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal regions: {len(combined)}", flush=True)

    # --- A1: AlleleDelta distribution TP vs FP ---
    print("\n--- A1: AlleleDelta Distribution ---")
    tp_delta = combined.loc[combined["category"] == "TP", "AlleleDelta"].dropna()
    fp_delta = combined.loc[combined["category"] == "FP", "AlleleDelta"].dropna()

    print(f"  TP AlleleDelta: n={len(tp_delta)}, "
          f"mean={tp_delta.mean():.4f}, median={tp_delta.median():.4f}, "
          f"std={tp_delta.std():.4f}")
    print(f"  FP AlleleDelta: n={len(fp_delta)}, "
          f"mean={fp_delta.mean():.4f}, median={fp_delta.median():.4f}, "
          f"std={fp_delta.std():.4f}")

    if len(tp_delta) > 0 and len(fp_delta) > 0:
        u_stat, u_p = stats.mannwhitneyu(tp_delta, fp_delta, alternative="two-sided")
        cohens_d = (tp_delta.mean() - fp_delta.mean()) / np.sqrt(
            (tp_delta.std()**2 + fp_delta.std()**2) / 2
        ) if (tp_delta.std() > 0 or fp_delta.std() > 0) else 0
        print(f"  Mann-Whitney U: stat={u_stat:.0f}, p={u_p:.2e}")
        print(f"  Cohen's d: {cohens_d:.4f}")

    # --- A2: Allele significance rates ---
    print("\n--- A2: Allele Significance Rates ---")
    for cat in ["TP", "FP"]:
        subset = combined[combined["category"] == cat]
        if "AlleleSig" in subset.columns:
            sig_rate = subset["AlleleSig"].mean()
            n_sig = subset["AlleleSig"].sum()
            print(f"  {cat}: {n_sig:.0f}/{len(subset)} = {sig_rate:.4f} "
                  f"({sig_rate*100:.1f}%) allele-significant")

    # --- A3: Per-sample breakdown ---
    print("\n--- A3: Per-Sample AlleleDelta ---")
    print(f"  {'Sample':<18} {'TP_mean':>8} {'FP_mean':>8} {'TP_sig%':>8} {'FP_sig%':>8} {'n_TP':>6} {'n_FP':>6}")
    print("  " + "-" * 65)

    per_sample_stats = []
    for sample in SAMPLES:
        tp_s = combined[(combined["sample"] == sample) & (combined["category"] == "TP")]
        fp_s = combined[(combined["sample"] == sample) & (combined["category"] == "FP")]

        tp_mean = tp_s["AlleleDelta"].mean() if len(tp_s) > 0 else float("nan")
        fp_mean = fp_s["AlleleDelta"].mean() if len(fp_s) > 0 else float("nan")
        tp_sig = tp_s["AlleleSig"].mean() * 100 if "AlleleSig" in tp_s.columns and len(tp_s) > 0 else float("nan")
        fp_sig = fp_s["AlleleSig"].mean() * 100 if "AlleleSig" in fp_s.columns and len(fp_s) > 0 else float("nan")

        print(f"  {sample:<18} {tp_mean:8.4f} {fp_mean:8.4f} {tp_sig:7.1f}% {fp_sig:7.1f}% {len(tp_s):6d} {len(fp_s):6d}")
        per_sample_stats.append({
            "sample": sample, "tp_mean_delta": tp_mean, "fp_mean_delta": fp_mean,
            "tp_sig_pct": tp_sig, "fp_sig_pct": fp_sig,
            "n_tp": len(tp_s), "n_fp": len(fp_s),
        })

    # --- A4: AlleleDelta by VerificationClass ---
    print("\n--- A4: AlleleDelta by VerificationClass ---")
    if "VerificationClass" in combined.columns:
        for vc in sorted(combined["VerificationClass"].dropna().unique()):
            subset = combined[combined["VerificationClass"] == vc]
            tp_v = subset[subset["category"] == "TP"]["AlleleDelta"].dropna()
            fp_v = subset[subset["category"] == "FP"]["AlleleDelta"].dropna()
            if len(tp_v) > 0 and len(fp_v) > 0:
                print(f"  {vc:20s}: TP={tp_v.mean():.4f}(n={len(tp_v)}) "
                      f"FP={fp_v.mean():.4f}(n={len(fp_v)}) "
                      f"delta={tp_v.mean()-fp_v.mean():.4f}")

    # --- A5: AlleleDelta by LOH status ---
    print("\n--- A5: AlleleDelta by LOH Status ---")
    if "Potential_LOH" in combined.columns:
        for loh in sorted(combined["Potential_LOH"].dropna().unique()):
            subset = combined[combined["Potential_LOH"] == loh]
            tp_v = subset[subset["category"] == "TP"]["AlleleDelta"].dropna()
            fp_v = subset[subset["category"] == "FP"]["AlleleDelta"].dropna()
            if len(tp_v) > 10 and len(fp_v) > 10:
                print(f"  LOH={loh}: TP={tp_v.mean():.4f}(n={len(tp_v)}) "
                      f"FP={fp_v.mean():.4f}(n={len(fp_v)})")

    # Save aggregate data
    agg_path = os.path.join(DATA_DIR, "part_a_aggregate.csv")
    pd.DataFrame(per_sample_stats).to_csv(agg_path, index=False)
    print(f"\n  Saved: {agg_path}")

    return combined


# ============================================================================
# Part B: Directional ASM from raw reads
# ============================================================================
def build_region_paths_from_csv(sample, category, sig_df=None):
    """Build region directory paths directly from significance_summary.csv.

    Avoids slow glob scanning by constructing paths from Chr/Pos columns.
    Directory structure: filtered_snv_{cat}/{Chr}/{Chr}_{Pos}/{Chr}_{Pos-5000}_{Pos+5000}/
    """
    run_id = RUN_PATTERNS.get(sample)
    if not run_id:
        return []

    base = os.path.join(
        CANONICAL_ROOT, sample, "paired_pileup", run_id,
        f"intersubmod_{category}", f"filtered_snv_{category}"
    )
    if not os.path.isdir(base):
        return []

    # Load significance_summary.csv if not provided
    if sig_df is None:
        csv_path = os.path.join(
            CANONICAL_ROOT, sample, "paired_pileup", run_id,
            f"intersubmod_{category}", "significance_summary.csv"
        )
        if not os.path.isfile(csv_path):
            return []
        sig_df = pd.read_csv(csv_path)

    regions = []
    window_half = 5000
    for _, row in sig_df.iterrows():
        chrom = str(row.get("Chr", ""))
        pos = row.get("Pos", 0)
        if not chrom or pd.isna(pos):
            continue
        pos = int(pos)
        start = pos - window_half
        end = pos + window_half

        region_dir = os.path.join(
            base, chrom, f"{chrom}_{pos}", f"{chrom}_{start}_{end}"
        )
        reads_path = os.path.join(region_dir, "reads", "reads.tsv")
        methyl_path = os.path.join(region_dir, "methylation", "methylation.csv")

        if os.path.isfile(reads_path) and os.path.isfile(methyl_path):
            regions.append({
                "region_id": f"{chrom}_{pos}",
                "reads_path": reads_path,
                "methyl_path": methyl_path,
                "chr": chrom,
                "pos": pos,
            })

    return regions


def compute_directional_asm(reads_path, methyl_path):
    """Compute directional ASM for a single region.

    Returns dict with:
      alt_mean_methyl: mean methylation of ALT reads
      ref_mean_methyl: mean methylation of REF reads
      signed_delta: ALT_mean - REF_mean (positive = ALT more methylated)
      n_alt: number of ALT reads
      n_ref: number of REF reads
      n_cpg: number of CpGs
      cpg_variance_mean: mean per-CpG variance (fCpG proxy)
      cpg_variance_high_frac: fraction of CpGs with variance > 0.05 (fCpG candidates)
    """
    try:
        reads = pd.read_csv(reads_path, sep="\t")
        methyl = pd.read_csv(methyl_path, index_col=0)
    except Exception:
        return None

    if len(reads) < 4 or methyl.shape[1] < 2:
        return None

    # Match reads by read_id
    alt_ids = set(reads.loc[reads["alt_support"] == "ALT", "read_id"].values)
    ref_ids = set(reads.loc[reads["alt_support"] == "REF", "read_id"].values)

    if len(alt_ids) < 2 or len(ref_ids) < 2:
        return None

    # Compute per-read mean methylation (averaging across CpGs, ignoring NA)
    methyl_numeric = methyl.apply(pd.to_numeric, errors="coerce")
    read_means = methyl_numeric.mean(axis=1)

    alt_means = read_means[read_means.index.isin(alt_ids)]
    ref_means = read_means[read_means.index.isin(ref_ids)]

    if len(alt_means) < 2 or len(ref_means) < 2:
        return None

    alt_mean = alt_means.mean()
    ref_mean = ref_means.mean()

    # CpG-level variance (fCpG proxy)
    cpg_var = methyl_numeric.var(axis=0).dropna()
    cpg_var_mean = cpg_var.mean() if len(cpg_var) > 0 else float("nan")
    cpg_high_var_frac = (cpg_var > 0.05).mean() if len(cpg_var) > 0 else float("nan")

    # Per-allele CpG methylation for allele-specific CpG variance
    alt_methyl = methyl_numeric.loc[methyl_numeric.index.isin(alt_ids)]
    ref_methyl = methyl_numeric.loc[methyl_numeric.index.isin(ref_ids)]

    # Per-CpG allele-specific difference
    alt_cpg_means = alt_methyl.mean(axis=0)
    ref_cpg_means = ref_methyl.mean(axis=0)
    cpg_diff = (alt_cpg_means - ref_cpg_means).dropna()
    cpg_diff_mean = cpg_diff.mean() if len(cpg_diff) > 0 else float("nan")
    cpg_diff_std = cpg_diff.std() if len(cpg_diff) > 0 else float("nan")

    return {
        "alt_mean_methyl": alt_mean,
        "ref_mean_methyl": ref_mean,
        "signed_delta": alt_mean - ref_mean,
        "n_alt": len(alt_means),
        "n_ref": len(ref_means),
        "n_cpg": methyl.shape[1],
        "cpg_variance_mean": cpg_var_mean,
        "cpg_high_var_frac": cpg_high_var_frac,
        "cpg_diff_mean": cpg_diff_mean,
        "cpg_diff_std": cpg_diff_std,
    }


def part_b_directional_asm(n_sample=N_SAMPLE_REGIONS):
    """Compute directional ASM for sampled regions across all samples."""
    print("\n" + "=" * 70)
    print("PART B: Directional ASM Analysis (ALT_mean_methyl - REF_mean_methyl)")
    print("=" * 70)

    random.seed(RANDOM_SEED)
    all_results = []

    for sample in SAMPLES:
        for cat in ["tp", "fp"]:
            print(f"\n  Processing {sample}/{cat}...", flush=True)
            regions = build_region_paths_from_csv(sample, cat)
            print(f"    Found {len(regions)} regions with reads+methylation data", flush=True)

            if not regions:
                continue

            # Sample regions
            if len(regions) > n_sample:
                sampled = random.sample(regions, n_sample)
            else:
                sampled = regions

            processed = 0
            skipped = 0
            for reg in sampled:
                result = compute_directional_asm(reg["reads_path"], reg["methyl_path"])
                if result is not None:
                    result["sample"] = sample
                    result["category"] = cat.upper()
                    result["region_id"] = reg["region_id"]
                    all_results.append(result)
                    processed += 1
                else:
                    skipped += 1

            print(f"    Processed: {processed}, Skipped: {skipped}")

    if not all_results:
        print("[ERROR] No directional ASM data computed")
        return None

    df = pd.DataFrame(all_results)

    # --- B1: Overall directional ASM ---
    print("\n--- B1: Directional ASM (ALT - REF) Distribution ---")
    tp_delta = df.loc[df["category"] == "TP", "signed_delta"].dropna()
    fp_delta = df.loc[df["category"] == "FP", "signed_delta"].dropna()

    print(f"  TP signed_delta: n={len(tp_delta)}, "
          f"mean={tp_delta.mean():.4f}, median={tp_delta.median():.4f}, "
          f"std={tp_delta.std():.4f}")
    print(f"  FP signed_delta: n={len(fp_delta)}, "
          f"mean={fp_delta.mean():.4f}, median={fp_delta.median():.4f}, "
          f"std={fp_delta.std():.4f}")

    # Is the mean significantly different from 0?
    if len(tp_delta) > 10:
        t_stat, t_p = stats.ttest_1samp(tp_delta, 0)
        print(f"  TP one-sample t-test (H0: delta=0): t={t_stat:.2f}, p={t_p:.2e}")

    if len(fp_delta) > 10:
        t_stat, t_p = stats.ttest_1samp(fp_delta, 0)
        print(f"  FP one-sample t-test (H0: delta=0): t={t_stat:.2f}, p={t_p:.2e}")

    # TP vs FP comparison
    if len(tp_delta) > 10 and len(fp_delta) > 10:
        u_stat, u_p = stats.mannwhitneyu(tp_delta, fp_delta, alternative="two-sided")
        print(f"  TP vs FP Mann-Whitney: p={u_p:.2e}")

    # Direction counts
    print(f"\n  TP: ALT>REF = {(tp_delta > 0).sum()} ({(tp_delta > 0).mean()*100:.1f}%), "
          f"ALT<REF = {(tp_delta < 0).sum()} ({(tp_delta < 0).mean()*100:.1f}%)")
    print(f"  FP: ALT>REF = {(fp_delta > 0).sum()} ({(fp_delta > 0).mean()*100:.1f}%), "
          f"ALT<REF = {(fp_delta < 0).sum()} ({(fp_delta < 0).mean()*100:.1f}%)")

    # --- B2: Per-sample directional ASM ---
    print("\n--- B2: Per-Sample Directional ASM ---")
    print(f"  {'Sample':<18} {'TP_mean':>8} {'FP_mean':>8} {'TP_ALT>REF%':>12} {'FP_ALT>REF%':>12}")
    print("  " + "-" * 62)
    for sample in SAMPLES:
        tp_s = df[(df["sample"] == sample) & (df["category"] == "TP")]["signed_delta"]
        fp_s = df[(df["sample"] == sample) & (df["category"] == "FP")]["signed_delta"]
        tp_pct = (tp_s > 0).mean() * 100 if len(tp_s) > 0 else float("nan")
        fp_pct = (fp_s > 0).mean() * 100 if len(fp_s) > 0 else float("nan")
        print(f"  {sample:<18} {tp_s.mean():8.4f} {fp_s.mean():8.4f} "
              f"{tp_pct:11.1f}% {fp_pct:11.1f}%")

    # --- B3: CpG-level variability (fCpG proxy) ---
    print("\n--- B3: CpG Variability (fCpG proxy from EVOFLUx) ---")
    tp_var = df.loc[df["category"] == "TP", "cpg_variance_mean"].dropna()
    fp_var = df.loc[df["category"] == "FP", "cpg_variance_mean"].dropna()
    print(f"  TP cpg_variance_mean: {tp_var.mean():.4f} (std={tp_var.std():.4f})")
    print(f"  FP cpg_variance_mean: {fp_var.mean():.4f} (std={fp_var.std():.4f})")

    tp_hv = df.loc[df["category"] == "TP", "cpg_high_var_frac"].dropna()
    fp_hv = df.loc[df["category"] == "FP", "cpg_high_var_frac"].dropna()
    print(f"  TP high-var CpG fraction: {tp_hv.mean():.4f}")
    print(f"  FP high-var CpG fraction: {fp_hv.mean():.4f}")

    if len(tp_var) > 10 and len(fp_var) > 10:
        u_stat, u_p = stats.mannwhitneyu(tp_var, fp_var, alternative="two-sided")
        print(f"  CpG variance TP vs FP Mann-Whitney: p={u_p:.2e}")

    # --- B4: Per-CpG allele-specific difference (ASM magnitude) ---
    print("\n--- B4: Per-CpG ASM Magnitude ---")
    tp_cpg_diff = df.loc[df["category"] == "TP", "cpg_diff_mean"].dropna()
    fp_cpg_diff = df.loc[df["category"] == "FP", "cpg_diff_mean"].dropna()
    print(f"  TP cpg_diff_mean (ALT-REF per CpG): mean={tp_cpg_diff.mean():.4f}, "
          f"std={tp_cpg_diff.std():.4f}")
    print(f"  FP cpg_diff_mean (ALT-REF per CpG): mean={fp_cpg_diff.mean():.4f}, "
          f"std={fp_cpg_diff.std():.4f}")

    # Absolute ASM magnitude
    tp_abs = tp_cpg_diff.abs()
    fp_abs = fp_cpg_diff.abs()
    print(f"  TP |ASM|: mean={tp_abs.mean():.4f}")
    print(f"  FP |ASM|: mean={fp_abs.mean():.4f}")

    if len(tp_abs) > 10 and len(fp_abs) > 10:
        u_stat, u_p = stats.mannwhitneyu(tp_abs, fp_abs, alternative="two-sided")
        print(f"  |ASM| TP vs FP Mann-Whitney: p={u_p:.2e}")

    # Save detailed results
    detail_path = os.path.join(DATA_DIR, "part_b_directional_asm.csv")
    df.to_csv(detail_path, index=False)
    print(f"\n  Saved: {detail_path}")

    return df


# ============================================================================
# Part C: Generate summary plots
# ============================================================================
def part_c_plots(agg_df, dir_df):
    """Generate summary visualizations."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print("\n[WARN] matplotlib/seaborn not available, skipping plots")
        return

    print("\n" + "=" * 70)
    print("PART C: Generating Plots")
    print("=" * 70)

    os.makedirs(FIGURES_DIR, exist_ok=True)

    # --- Plot 1: AlleleDelta TP vs FP violin ---
    if agg_df is not None:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left: AlleleDelta distribution
        ax = axes[0]
        tp_data = agg_df.loc[agg_df["category"] == "TP", "AlleleDelta"].dropna()
        fp_data = agg_df.loc[agg_df["category"] == "FP", "AlleleDelta"].dropna()
        # Clip for visibility
        tp_clip = tp_data.clip(-0.5, 0.5)
        fp_clip = fp_data.clip(-0.5, 0.5)

        parts = ax.violinplot([tp_clip, fp_clip], positions=[1, 2], showmeans=True, showmedians=True)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["TP", "FP"])
        ax.set_ylabel("AlleleDelta (PERMANOVA)")
        ax.set_title("A1: AlleleDelta Distribution")
        ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)

        # Right: Per-sample AlleleDelta comparison
        ax = axes[1]
        sample_data = []
        for sample in SAMPLES:
            for cat in ["TP", "FP"]:
                subset = agg_df[(agg_df["sample"] == sample) & (agg_df["category"] == cat)]
                mean_d = subset["AlleleDelta"].mean()
                sample_data.append({"sample": sample, "category": cat, "mean_delta": mean_d})
        sdf = pd.DataFrame(sample_data)
        for i, cat in enumerate(["TP", "FP"]):
            sub = sdf[sdf["category"] == cat]
            x = np.arange(len(SAMPLES))
            ax.bar(x + i * 0.35, sub["mean_delta"].values, 0.3, label=cat,
                   color=["#2196F3", "#F44336"][i], alpha=0.7)
        ax.set_xticks(np.arange(len(SAMPLES)) + 0.175)
        ax.set_xticklabels([s[:8] for s in SAMPLES], rotation=45, ha="right")
        ax.set_ylabel("Mean AlleleDelta")
        ax.set_title("A3: Per-Sample AlleleDelta")
        ax.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "01_allele_delta_aggregate.png"), dpi=150)
        plt.close()
        print("  Saved: figures/01_allele_delta_aggregate.png")

    # --- Plot 2: Directional ASM ---
    if dir_df is not None:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Top-left: signed delta histogram
        ax = axes[0, 0]
        tp_d = dir_df.loc[dir_df["category"] == "TP", "signed_delta"].dropna()
        fp_d = dir_df.loc[dir_df["category"] == "FP", "signed_delta"].dropna()
        bins = np.linspace(-0.3, 0.3, 61)
        ax.hist(tp_d.clip(-0.3, 0.3), bins=bins, alpha=0.6, label=f"TP (n={len(tp_d)})", color="#2196F3")
        ax.hist(fp_d.clip(-0.3, 0.3), bins=bins, alpha=0.6, label=f"FP (n={len(fp_d)})", color="#F44336")
        ax.axvline(x=0, color="black", linestyle="--")
        ax.axvline(x=tp_d.mean(), color="#1565C0", linestyle="-", label=f"TP mean={tp_d.mean():.4f}")
        ax.axvline(x=fp_d.mean(), color="#C62828", linestyle="-", label=f"FP mean={fp_d.mean():.4f}")
        ax.set_xlabel("Signed Delta (ALT mean methyl - REF mean methyl)")
        ax.set_ylabel("Count")
        ax.set_title("B1: Directional ASM Distribution")
        ax.legend(fontsize=8)

        # Top-right: per-sample directional
        ax = axes[0, 1]
        for i, cat in enumerate(["TP", "FP"]):
            means = []
            for sample in SAMPLES:
                sub = dir_df[(dir_df["sample"] == sample) & (dir_df["category"] == cat)]
                means.append(sub["signed_delta"].mean())
            x = np.arange(len(SAMPLES))
            ax.bar(x + i * 0.35, means, 0.3, label=cat,
                   color=["#2196F3", "#F44336"][i], alpha=0.7)
        ax.set_xticks(np.arange(len(SAMPLES)) + 0.175)
        ax.set_xticklabels([s[:8] for s in SAMPLES], rotation=45, ha="right")
        ax.set_ylabel("Mean Signed Delta")
        ax.set_title("B2: Per-Sample Directional ASM")
        ax.axhline(y=0, color="gray", linestyle="--")
        ax.legend()

        # Bottom-left: CpG variance
        ax = axes[1, 0]
        tp_v = dir_df.loc[dir_df["category"] == "TP", "cpg_variance_mean"].dropna()
        fp_v = dir_df.loc[dir_df["category"] == "FP", "cpg_variance_mean"].dropna()
        bins = np.linspace(0, 0.15, 61)
        ax.hist(tp_v.clip(0, 0.15), bins=bins, alpha=0.6, label=f"TP (n={len(tp_v)})", color="#2196F3")
        ax.hist(fp_v.clip(0, 0.15), bins=bins, alpha=0.6, label=f"FP (n={len(fp_v)})", color="#F44336")
        ax.set_xlabel("Mean CpG Variance (fCpG proxy)")
        ax.set_ylabel("Count")
        ax.set_title("B3: CpG Variability (EVOFLUx fCpG concept)")
        ax.legend(fontsize=8)

        # Bottom-right: |ASM| magnitude
        ax = axes[1, 1]
        tp_abs = dir_df.loc[dir_df["category"] == "TP", "cpg_diff_mean"].abs().dropna()
        fp_abs = dir_df.loc[dir_df["category"] == "FP", "cpg_diff_mean"].abs().dropna()
        bins = np.linspace(0, 0.2, 61)
        ax.hist(tp_abs.clip(0, 0.2), bins=bins, alpha=0.6, label=f"TP (mean={tp_abs.mean():.4f})", color="#2196F3")
        ax.hist(fp_abs.clip(0, 0.2), bins=bins, alpha=0.6, label=f"FP (mean={fp_abs.mean():.4f})", color="#F44336")
        ax.set_xlabel("|ASM| per CpG (ALT-REF absolute)")
        ax.set_ylabel("Count")
        ax.set_title("B4: ASM Magnitude per CpG")
        ax.legend(fontsize=8)

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "02_directional_asm_detail.png"), dpi=150)
        plt.close()
        print("  Saved: figures/02_directional_asm_detail.png")

    # --- Plot 3: ALT/REF mean methylation scatter ---
    if dir_df is not None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for i, cat in enumerate(["TP", "FP"]):
            ax = axes[i]
            sub = dir_df[dir_df["category"] == cat]
            ax.scatter(sub["ref_mean_methyl"], sub["alt_mean_methyl"],
                       alpha=0.15, s=5, color=["#2196F3", "#F44336"][i])
            ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="ALT=REF")
            ax.set_xlabel("REF Mean Methylation")
            ax.set_ylabel("ALT Mean Methylation")
            ax.set_title(f"{cat}: ALT vs REF Methylation (n={len(sub)})")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_aspect("equal")
            ax.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "03_alt_ref_scatter.png"), dpi=150)
        plt.close()
        print("  Saved: figures/03_alt_ref_scatter.png")


# ============================================================================
# Main
# ============================================================================
def main():
    print("=" * 70)
    print("L1 Directional ASM Analysis — Literature Validation")
    print("Testing epiTRACERx/Gaiti/Basilicata predictions against ISM data")
    print("=" * 70)

    # Part A: fast aggregate analysis
    agg_df = part_a_aggregate_analysis()

    # Part B: directional ASM from raw reads (sampled)
    dir_df = part_b_directional_asm()

    # Part C: plots
    part_c_plots(agg_df, dir_df)

    # --- Final Summary ---
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print("""
    L1 Hypothesis (epiTRACERx): ALT allele shows systematic methylation
    directionality at TP somatic SNV sites.

    Key questions answered:
    1. Does AlleleDelta (PERMANOVA) differ TP vs FP?
    2. Is there a consistent directional ASM (ALT>REF or ALT<REF)?
    3. Does CpG-level variability (fCpG proxy) differ TP vs FP?
    4. Is |ASM| (absolute allele-specific methylation difference) discriminative?

    See figures/ for visualizations.
    See data/ for detailed per-region results.
    """)


if __name__ == "__main__":
    main()
