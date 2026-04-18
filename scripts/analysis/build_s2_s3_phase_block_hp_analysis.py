#!/usr/bin/env python3
"""
S2+S3 Phase Block / HP Tag Analysis
====================================
Agent B: LongPhase-TO vs LongPhase-S phase block structure and HP concordance.

Sections:
  S2.1   TO phase block size distribution
  S2.1b  Paired phase block size distribution + comparison plot
  S2.2   LOH.bed genome coverage
  S3.1   Same-locus HP concordance (from master dataset)
  S3.2   HP assignment rate stratified analysis
  S3.3   HP imbalance at TP sites
"""

import gzip
import os
import sys
import warnings
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

OUT_DIR = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260402_phasing_causal_chain"
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# hg38 autosomal + chrX lengths (GRCh38 reference)
# ============================================================
HG38_CHR_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895,
}
HG38_AUTOSOMAL_X_TOTAL = sum(HG38_CHR_LENGTHS.values())

SAMPLE_ROUND_IDS = {
    "HCC1395": "20260315_hcc1395_to_pilot",
    "HCC1395_DORADO": "20260315_hcc1395_dorado_to_pilot",
    "COLO829": "20260317_colo829_to_pilot_1",
    "H1437": "20260318_h1437_to_pilot_fastresume",
    "H2009": "20260318_h2009_to_pilot_fastresume",
    "HCC1937": "20260318_hcc1937_to_pilot_fastresume",
    "HCC1954": "20260318_hcc1954_to_pilot_fastresume",
}

MASTER_DATASET = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"

# ============================================================
# Utility: parse VCF for phase blocks
# ============================================================
def parse_vcf_phase_blocks(vcf_path, is_gzipped=False):
    """
    Parse a VCF file and extract phase block information.
    Returns dict: {(chr, PS): [(pos1, pos2, ...), ...]}
    """
    blocks = defaultdict(list)  # key=(chrom, ps_value) -> list of positions
    opener = gzip.open if is_gzipped else open
    mode = "rt" if is_gzipped else "r"

    with opener(vcf_path, mode) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            fmt_keys = fields[8].split(":")
            fmt_vals = fields[9].split(":")

            # find PS index
            ps_val = None
            gt_val = None
            for i, k in enumerate(fmt_keys):
                if k == "PS" and i < len(fmt_vals):
                    ps_val = fmt_vals[i]
                if k == "GT" and i < len(fmt_vals):
                    gt_val = fmt_vals[i]

            # Only count phased variants (contain |)
            if gt_val is None or "|" not in gt_val:
                continue
            if ps_val is None or ps_val == ".":
                continue

            blocks[(chrom, ps_val)].append(pos)

    return blocks


def compute_block_stats(blocks):
    """
    Given blocks dict, compute per-block stats.
    Returns list of dicts with chr, ps, start, end, size_bp, n_variants.
    """
    rows = []
    for (chrom, ps), positions in blocks.items():
        positions_sorted = sorted(positions)
        start = positions_sorted[0]
        end = positions_sorted[-1]
        size_bp = end - start + 1 if len(positions_sorted) > 1 else 1
        rows.append({
            "chr": chrom,
            "ps": ps,
            "start": start,
            "end": end,
            "size_bp": size_bp,
            "n_variants": len(positions_sorted),
        })
    return rows


def compute_per_chr_summary(block_rows):
    """Compute per-chromosome summary: N50, median, #blocks, max."""
    df = pd.DataFrame(block_rows)
    if df.empty:
        return pd.DataFrame()

    results = []
    for chrom, grp in df.groupby("chr"):
        sizes = grp["size_bp"].values
        sizes_sorted = np.sort(sizes)[::-1]
        total = sizes_sorted.sum()
        cumsum = np.cumsum(sizes_sorted)
        n50_idx = np.searchsorted(cumsum, total / 2)
        n50 = sizes_sorted[min(n50_idx, len(sizes_sorted) - 1)]
        results.append({
            "chr": chrom,
            "n_blocks": len(grp),
            "total_variants": grp["n_variants"].sum(),
            "median_block_bp": int(np.median(sizes)),
            "mean_block_bp": int(np.mean(sizes)),
            "n50_block_bp": int(n50),
            "max_block_bp": int(sizes.max()),
            "min_block_bp": int(sizes.min()),
        })
    return pd.DataFrame(results).sort_values("chr", key=lambda x: x.map(
        {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24}
    ))


# ============================================================
# S2.1: TO Phase Block Stats
# ============================================================
def run_s2_1():
    print("=" * 70)
    print("S2.1: TO Phase Block Size Distribution (HCC1395)")
    print("=" * 70)

    vcf_path = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased.vcf"
    print(f"  Parsing: {vcf_path}")
    blocks = parse_vcf_phase_blocks(vcf_path, is_gzipped=False)
    block_rows = compute_block_stats(blocks)
    print(f"  Total phase blocks: {len(block_rows)}")

    df_blocks = pd.DataFrame(block_rows)
    if df_blocks.empty:
        print("  WARNING: No phase blocks found!")
        return df_blocks, pd.DataFrame()

    phased_variants = df_blocks["n_variants"].sum()
    print(f"  Total phased variants in blocks: {phased_variants}")
    print(f"  Overall median block size: {df_blocks['size_bp'].median():.0f} bp")
    print(f"  Overall mean block size: {df_blocks['size_bp'].mean():.0f} bp")
    print(f"  Overall max block size: {df_blocks['size_bp'].max():.0f} bp")

    chr_summary = compute_per_chr_summary(block_rows)
    out_path = os.path.join(OUT_DIR, "s2_1_to_phase_block_stats.tsv")
    chr_summary.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")
    print(chr_summary.to_string(index=False))
    print()

    return df_blocks, chr_summary


# ============================================================
# S2.1b: Paired Phase Block Stats + Comparison
# ============================================================
def run_s2_1b(to_block_df):
    print("=" * 70)
    print("S2.1b: Paired Phase Block Size Distribution (HCC1395)")
    print("=" * 70)

    vcf_path = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
    print(f"  Parsing: {vcf_path}")
    blocks = parse_vcf_phase_blocks(vcf_path, is_gzipped=True)
    block_rows = compute_block_stats(blocks)
    print(f"  Total phase blocks: {len(block_rows)}")

    df_blocks = pd.DataFrame(block_rows)
    if df_blocks.empty:
        print("  WARNING: No phase blocks found!")
        return

    phased_variants = df_blocks["n_variants"].sum()
    print(f"  Total phased variants in blocks: {phased_variants}")
    print(f"  Overall median block size: {df_blocks['size_bp'].median():.0f} bp")
    print(f"  Overall mean block size: {df_blocks['size_bp'].mean():.0f} bp")
    print(f"  Overall max block size: {df_blocks['size_bp'].max():.0f} bp")

    chr_summary = compute_per_chr_summary(block_rows)
    out_path = os.path.join(OUT_DIR, "s2_1b_paired_phase_block_stats.tsv")
    chr_summary.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")
    print(chr_summary.to_string(index=False))
    print()

    # --- Comparison plot ---
    print("  Generating comparison plot...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram comparison
    ax = axes[0]
    to_sizes = to_block_df["size_bp"].values
    paired_sizes = df_blocks["size_bp"].values

    # Use log-scale bins
    all_sizes = np.concatenate([to_sizes[to_sizes > 0], paired_sizes[paired_sizes > 0]])
    if len(all_sizes) > 0:
        bins = np.logspace(np.log10(max(1, all_sizes.min())), np.log10(all_sizes.max()), 60)
    else:
        bins = np.logspace(0, 8, 60)

    ax.hist(to_sizes[to_sizes > 0], bins=bins, alpha=0.6, label=f"TO (n={len(to_sizes)})", color="tab:orange", edgecolor="none")
    ax.hist(paired_sizes[paired_sizes > 0], bins=bins, alpha=0.6, label=f"Paired (n={len(paired_sizes)})", color="tab:blue", edgecolor="none")
    ax.set_xscale("log")
    ax.set_xlabel("Phase Block Size (bp)")
    ax.set_ylabel("Count")
    ax.set_title("Phase Block Size Distribution\n(HCC1395, TO vs Paired)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Cumulative distribution
    ax2 = axes[1]
    for sizes, label, color in [
        (to_sizes, "TO", "tab:orange"),
        (paired_sizes, "Paired", "tab:blue"),
    ]:
        sorted_s = np.sort(sizes[sizes > 0])
        cdf = np.arange(1, len(sorted_s) + 1) / len(sorted_s)
        ax2.plot(sorted_s, cdf, label=label, color=color, linewidth=1.5)
    ax2.set_xscale("log")
    ax2.set_xlabel("Phase Block Size (bp)")
    ax2.set_ylabel("Cumulative Fraction")
    ax2.set_title("CDF of Phase Block Sizes")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(OUT_DIR, "s2_1_phase_block_size_comparison.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fig_path}")
    print()

    # Summary comparison table
    print("  --- Phase Block Summary Comparison ---")
    for label, sizes in [("TO", to_sizes), ("Paired", paired_sizes)]:
        s = sizes[sizes > 0]
        ss = np.sort(s)[::-1]
        total = ss.sum()
        cumsum = np.cumsum(ss)
        n50_idx = np.searchsorted(cumsum, total / 2)
        n50 = ss[min(n50_idx, len(ss) - 1)]
        print(f"  {label:8s}: N={len(s):>8,}  median={np.median(s):>12,.0f} bp  "
              f"mean={np.mean(s):>12,.0f} bp  N50={n50:>12,} bp  max={s.max():>12,} bp")
    print()


# ============================================================
# S2.2: LOH.bed Genome Coverage
# ============================================================
def run_s2_2():
    print("=" * 70)
    print("S2.2: LOH.bed Genome Coverage")
    print("=" * 70)

    results = []
    for sample, round_id in SAMPLE_ROUND_IDS.items():
        bed_path = (
            f"/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            f"{round_id}/step03_longphase_to/tumor_phased_LOH.bed"
        )
        if not os.path.exists(bed_path):
            print(f"  {sample}: LOH.bed not found at {bed_path}")
            continue

        # Show first 10 lines for first sample
        if sample == "HCC1395":
            print(f"\n  First 10 lines of {sample} LOH.bed:")
            with open(bed_path) as fh:
                for i, line in enumerate(fh):
                    if i >= 10:
                        break
                    print(f"    {line.rstrip()}")
            print()

        total_bp = 0
        n_regions = 0
        per_chr_bp = defaultdict(int)
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#") or line.strip() == "":
                    continue
                parts = line.strip().split("\t")
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                region_bp = end - start
                total_bp += region_bp
                per_chr_bp[chrom] += region_bp
                n_regions += 1

        frac = total_bp / HG38_AUTOSOMAL_X_TOTAL
        results.append({
            "sample": sample,
            "round_id": round_id,
            "n_loh_regions": n_regions,
            "total_loh_bp": total_bp,
            "genome_total_bp": HG38_AUTOSOMAL_X_TOTAL,
            "loh_coverage_fraction": frac,
            "loh_coverage_pct": frac * 100,
        })
        print(f"  {sample:18s}: {n_regions:>6,} regions  {total_bp:>14,} bp  "
              f"coverage={frac*100:.2f}%")

    df = pd.DataFrame(results)
    out_path = os.path.join(OUT_DIR, "s2_2_loh_bed_genome_coverage.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n  Saved: {out_path}")
    print()
    return df


# ============================================================
# S3.1: Same-Locus HP Concordance
# ============================================================
def run_s3_1(master):
    print("=" * 70)
    print("S3.1: Same-Locus HP Concordance")
    print("=" * 70)

    # Separate by mode
    paired = master[master["mode"] == "paired"].copy()
    to = master[master["mode"] == "to"].copy()

    print(f"  Paired rows: {len(paired)}")
    print(f"  TO rows: {len(to)}")

    # Create same-locus pairs
    paired_key = paired.set_index(["variant_key", "sample"])
    to_key = to.set_index(["variant_key", "sample"])
    common_idx = paired_key.index.intersection(to_key.index)
    print(f"  Same-locus pairs: {len(common_idx)}")

    if len(common_idx) == 0:
        print("  No same-locus pairs found!")
        return

    p_df = paired_key.loc[common_idx].reset_index()
    t_df = to_key.loc[common_idx].reset_index()

    # Merge on variant_key + sample
    merged = p_df[["variant_key", "sample", "HP_Ratio", "core_loh_like", "truth_label",
                    "effective_hp_reads", "HP1FamilyN", "HP2FamilyN"]].merge(
        t_df[["variant_key", "sample", "HP_Ratio", "core_loh_like", "truth_label",
               "effective_hp_reads", "HP1FamilyN", "HP2FamilyN"]],
        on=["variant_key", "sample"],
        suffixes=("_paired", "_to"),
    )
    print(f"  Merged pairs: {len(merged)}")

    # --- HP_Ratio Spearman correlation (per sample), with swap-adjustment ---
    corr_results = []
    samples = sorted(merged["sample"].unique())
    for s in samples:
        sub = merged[merged["sample"] == s].dropna(subset=["HP_Ratio_paired", "HP_Ratio_to"])
        if len(sub) < 10:
            continue
        x = sub["HP_Ratio_paired"].values
        y = sub["HP_Ratio_to"].values

        # Raw correlation
        rho_raw, p_raw = stats.spearmanr(x, y)
        # Swap-adjusted: try 1-y
        rho_swap, p_swap = stats.spearmanr(x, 1 - y)
        # Pick higher absolute
        if abs(rho_swap) > abs(rho_raw):
            rho_best, p_best = rho_swap, p_swap
            best_type = "swapped"
        else:
            rho_best, p_best = rho_raw, p_raw
            best_type = "raw"

        corr_results.append({
            "sample": s,
            "n_pairs": len(sub),
            "spearman_raw": rho_raw,
            "p_raw": p_raw,
            "spearman_swapped": rho_swap,
            "p_swapped": p_swap,
            "best_rho": rho_best,
            "best_type": best_type,
        })
        print(f"  {s:18s}  n={len(sub):>5,}  rho_raw={rho_raw:+.4f}  "
              f"rho_swap={rho_swap:+.4f}  best={rho_best:+.4f} ({best_type})")

    df_corr = pd.DataFrame(corr_results)
    out_path = os.path.join(OUT_DIR, "s3_1_hp_concordance_by_sample.tsv")
    df_corr.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")

    # --- LOH concordance matrix ---
    print("\n  LOH binary concordance:")
    merged["loh_paired"] = merged["core_loh_like_paired"].astype(bool)
    merged["loh_to"] = merged["core_loh_like_to"].astype(bool)

    loh_matrix_rows = []
    for s in samples:
        sub = merged[merged["sample"] == s]
        both_loh = ((sub["loh_paired"]) & (sub["loh_to"])).sum()
        both_non = ((~sub["loh_paired"]) & (~sub["loh_to"])).sum()
        flip_to_loh = ((~sub["loh_paired"]) & (sub["loh_to"])).sum()
        flip_paired_loh = ((sub["loh_paired"]) & (~sub["loh_to"])).sum()
        total = len(sub)
        concordance = (both_loh + both_non) / total if total > 0 else 0
        loh_matrix_rows.append({
            "sample": s,
            "n_pairs": total,
            "both_LOH": both_loh,
            "both_nonLOH": both_non,
            "flip_TO_LOH_only": flip_to_loh,
            "flip_Paired_LOH_only": flip_paired_loh,
            "concordance": concordance,
        })
        print(f"  {s:18s}  concordance={concordance:.3f}  "
              f"both_LOH={both_loh:>5}  both_non={both_non:>5}  "
              f"TO_only_LOH={flip_to_loh:>5}  Paired_only_LOH={flip_paired_loh:>5}")

    df_loh = pd.DataFrame(loh_matrix_rows)
    out_path2 = os.path.join(OUT_DIR, "s3_1_loh_concordance_matrix.tsv")
    df_loh.to_csv(out_path2, sep="\t", index=False)
    print(f"  Saved: {out_path2}")

    # --- Stratified analysis: HP concordance in non-LOH subset ---
    print("\n  Stratified: HP correlation at non-LOH sites only:")
    for s in samples:
        sub = merged[(merged["sample"] == s) &
                     (~merged["loh_paired"]) &
                     (~merged["loh_to"])].dropna(subset=["HP_Ratio_paired", "HP_Ratio_to"])
        if len(sub) < 10:
            print(f"  {s:18s}  n={len(sub)}  (too few)")
            continue
        x = sub["HP_Ratio_paired"].values
        y = sub["HP_Ratio_to"].values
        rho_raw, _ = stats.spearmanr(x, y)
        rho_swap, _ = stats.spearmanr(x, 1 - y)
        rho_best = rho_raw if abs(rho_raw) >= abs(rho_swap) else rho_swap
        print(f"  {s:18s}  n={len(sub):>5,}  best_rho={rho_best:+.4f}")

    # --- Scatter plot ---
    print("\n  Generating scatter plot...")
    # Determine grid: one facet per sample
    n_samples = len(samples)
    ncols = min(4, n_samples)
    nrows = (n_samples + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows),
                             squeeze=False)

    truth_colors = {"TP": "tab:green", "FP": "tab:red", "FN": "tab:purple"}

    for idx, s in enumerate(samples):
        ax = axes[idx // ncols][idx % ncols]
        sub = merged[merged["sample"] == s].dropna(subset=["HP_Ratio_paired", "HP_Ratio_to"])
        for truth in ["FP", "TP", "FN"]:
            t_sub = sub[sub["truth_label_paired"] == truth]
            if len(t_sub) == 0:
                continue
            ax.scatter(
                t_sub["HP_Ratio_paired"], t_sub["HP_Ratio_to"],
                c=truth_colors.get(truth, "gray"), alpha=0.3, s=8, label=truth,
                edgecolors="none",
            )
        ax.set_xlabel("Paired HP_Ratio")
        ax.set_ylabel("TO HP_Ratio")
        ax.set_title(s, fontsize=10)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.5)
        ax.plot([0, 1], [1, 0], "k:", alpha=0.3, linewidth=0.5)  # swap line
        if idx == 0:
            ax.legend(fontsize=7, markerscale=2)
        ax.set_aspect("equal")

    # Hide unused axes
    for idx in range(n_samples, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle("S3.1: HP_Ratio Scatter (Paired vs TO, per sample)", fontsize=12, y=1.02)
    plt.tight_layout()
    fig_path = os.path.join(OUT_DIR, "s3_1_hp_ratio_scatter_by_sample.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fig_path}")
    print()

    return merged


# ============================================================
# S3.2: HP Assignment Rate Stratified Analysis
# ============================================================
def run_s3_2(master):
    print("=" * 70)
    print("S3.2: HP Assignment Rate by Mode x Sample")
    print("=" * 70)

    cols_needed = ["sample", "mode", "hp_assign_rate", "hp0_ratio", "hp3_ratio"]
    sub = master[cols_needed].dropna(subset=["hp_assign_rate"])

    agg_rows = []
    for (sample, mode), grp in sub.groupby(["sample", "mode"]):
        agg_rows.append({
            "sample": sample,
            "mode": mode,
            "n": len(grp),
            "mean_hp_assign_rate": grp["hp_assign_rate"].mean(),
            "median_hp_assign_rate": grp["hp_assign_rate"].median(),
            "mean_hp0_ratio": grp["hp0_ratio"].mean(),
            "median_hp0_ratio": grp["hp0_ratio"].median(),
            "mean_hp3_ratio": grp["hp3_ratio"].mean(),
            "median_hp3_ratio": grp["hp3_ratio"].median(),
        })

    df_agg = pd.DataFrame(agg_rows)
    out_path = os.path.join(OUT_DIR, "s3_2_hp_assignment_rate_by_mode_sample.tsv")
    df_agg.to_csv(out_path, sep="\t", index=False)
    print(df_agg.to_string(index=False))
    print(f"\n  Saved: {out_path}")

    # Summary: NHP0 (unphased) and NHP3 (conflict) proportions
    print("\n  NHP0 (unphased) and NHP3 (conflict) comparison:")
    for mode in ["paired", "to"]:
        mode_sub = sub[sub["mode"] == mode]
        print(f"  {mode:8s}: mean hp0_ratio={mode_sub['hp0_ratio'].mean():.4f}  "
              f"mean hp3_ratio={mode_sub['hp3_ratio'].mean():.4f}  "
              f"mean hp_assign_rate={mode_sub['hp_assign_rate'].mean():.4f}")
    print()


# ============================================================
# S3.3: HP Imbalance at TP Sites
# ============================================================
def run_s3_3(master):
    print("=" * 70)
    print("S3.3: HP Imbalance at TP Sites")
    print("=" * 70)

    master = master.copy()
    master["hp_imbalance"] = (master["HP_Ratio"] - 0.5).abs()

    # tier by effective_hp_reads
    bins = [0, 10, 30, 50, float("inf")]
    labels = ["<10", "10-29", "30-49", ">=50"]
    master["ehp_tier"] = pd.cut(master["effective_hp_reads"], bins=bins, labels=labels, right=False)

    # Stratified analysis: mode x truth_label x sample
    cols = ["sample", "mode", "truth_label", "ehp_tier"]
    rows = []
    for keys, grp in master.groupby(cols, observed=True):
        sample, mode, truth, tier = keys
        rows.append({
            "sample": sample,
            "mode": mode,
            "truth_label": truth,
            "ehp_tier": tier,
            "n": len(grp),
            "mean_hp_imbalance": grp["hp_imbalance"].mean(),
            "median_hp_imbalance": grp["hp_imbalance"].median(),
            "mean_hp_ratio": grp["HP_Ratio"].mean(),
        })
    df_strat = pd.DataFrame(rows)
    out_path = os.path.join(OUT_DIR, "s3_3_hp_imbalance_stratified.tsv")
    df_strat.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")

    # Print summary: overall mode x truth
    print("\n  Overall HP Imbalance by mode x truth:")
    for (mode, truth), grp in master.groupby(["mode", "truth_label"]):
        print(f"  {mode:8s} {truth:4s}  n={len(grp):>6,}  "
              f"mean |HP-0.5|={grp['hp_imbalance'].mean():.4f}  "
              f"median={grp['hp_imbalance'].median():.4f}")

    # Verify: TO TP >> Paired TP
    print("\n  Verification: TO TP vs Paired TP imbalance")
    to_tp = master[(master["mode"] == "to") & (master["truth_label"] == "TP")]["hp_imbalance"]
    paired_tp = master[(master["mode"] == "paired") & (master["truth_label"] == "TP")]["hp_imbalance"]
    if len(to_tp) > 0 and len(paired_tp) > 0:
        u_stat, u_p = stats.mannwhitneyu(to_tp, paired_tp, alternative="greater")
        print(f"  TO TP mean={to_tp.mean():.4f} vs Paired TP mean={paired_tp.mean():.4f}")
        print(f"  Mann-Whitney U (TO > Paired): U={u_stat:.0f}, p={u_p:.2e}")
        print(f"  Ratio TO/Paired: {to_tp.mean()/paired_tp.mean():.2f}x")

    # By tier
    print("\n  By effective_hp_reads tier:")
    for tier in labels:
        to_t = master[(master["mode"] == "to") & (master["truth_label"] == "TP") &
                       (master["ehp_tier"] == tier)]["hp_imbalance"]
        pa_t = master[(master["mode"] == "paired") & (master["truth_label"] == "TP") &
                       (master["ehp_tier"] == tier)]["hp_imbalance"]
        if len(to_t) > 5 and len(pa_t) > 5:
            ratio = to_t.mean() / pa_t.mean() if pa_t.mean() > 0 else float("inf")
            print(f"  tier={tier:>5s}  TO_TP n={len(to_t):>5,} mean={to_t.mean():.4f}  "
                  f"Paired_TP n={len(pa_t):>5,} mean={pa_t.mean():.4f}  ratio={ratio:.2f}x")

    # --- Violin plot ---
    print("\n  Generating violin plot...")
    # Create mode_truth label
    plot_df = master[master["truth_label"].isin(["TP", "FP"])].copy()
    plot_df["mode_truth"] = plot_df["mode"] + " " + plot_df["truth_label"]

    samples = sorted(plot_df["sample"].unique())
    n_samples = len(samples)
    ncols = min(4, n_samples)
    nrows = (n_samples + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows),
                             squeeze=False)

    order = ["paired TP", "paired FP", "to TP", "to FP"]
    palette = {"paired TP": "tab:green", "paired FP": "tab:red",
               "to TP": "tab:cyan", "to FP": "tab:orange"}

    for idx, s in enumerate(samples):
        ax = axes[idx // ncols][idx % ncols]
        s_df = plot_df[plot_df["sample"] == s]
        # Filter to categories present
        present_order = [o for o in order if o in s_df["mode_truth"].values]
        if len(present_order) == 0:
            ax.set_visible(False)
            continue
        sns.violinplot(
            data=s_df, x="mode_truth", y="hp_imbalance", order=present_order,
            palette=palette, ax=ax, cut=0, inner="quartile", linewidth=0.5,
            density_norm="width",
        )
        ax.set_title(s, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("|HP_Ratio - 0.5|")
        ax.set_ylim(-0.02, 0.55)
        ax.tick_params(axis="x", rotation=30, labelsize=8)

    for idx in range(n_samples, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle("S3.3: HP Imbalance by Mode x Truth (per sample)", fontsize=12, y=1.02)
    plt.tight_layout()
    fig_path = os.path.join(OUT_DIR, "s3_3_hp_imbalance_by_mode_truth.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fig_path}")
    print()


# ============================================================
# MAIN
# ============================================================
def main():
    print("#" * 70)
    print("# S2+S3 Phase Block / HP Tag Analysis")
    print(f"# Output: {OUT_DIR}")
    print("#" * 70)
    print()

    # --- S2.1 ---
    to_block_df, to_chr_summary = run_s2_1()

    # --- S2.1b ---
    if not to_block_df.empty:
        run_s2_1b(to_block_df)

    # --- S2.2 ---
    loh_cov_df = run_s2_2()

    # --- Load master dataset ---
    print("=" * 70)
    print("Loading master dataset...")
    print("=" * 70)
    master = pd.read_csv(MASTER_DATASET, sep="\t")
    print(f"  Rows: {len(master):,}")
    print(f"  Samples: {sorted(master['sample'].unique())}")
    print(f"  Modes: {sorted(master['mode'].unique())}")
    print()

    # --- S3.1 ---
    merged = run_s3_1(master)

    # --- S3.2 ---
    run_s3_2(master)

    # --- S3.3 ---
    run_s3_3(master)

    # ============================================================
    # SUMMARY
    # ============================================================
    print()
    print("#" * 70)
    print("# SUMMARY")
    print("#" * 70)
    print()

    if not to_block_df.empty:
        to_sizes = to_block_df["size_bp"].values
        print(f"S2.1  TO phase blocks: {len(to_block_df):,}, "
              f"median={np.median(to_sizes):,.0f} bp, "
              f"mean={np.mean(to_sizes):,.0f} bp")

    if loh_cov_df is not None and len(loh_cov_df) > 0:
        for _, row in loh_cov_df.iterrows():
            print(f"S2.2  {row['sample']:18s} LOH coverage: "
                  f"{row['loh_coverage_pct']:.2f}% ({row['n_loh_regions']:,} regions)")

    # Key S3 numbers
    if merged is not None and len(merged) > 0:
        print(f"\nS3.1  Same-locus pairs: {len(merged):,}")

    # HP imbalance key number
    to_tp_imb = master[(master["mode"] == "to") & (master["truth_label"] == "TP")]["HP_Ratio"].apply(lambda x: abs(x - 0.5))
    pa_tp_imb = master[(master["mode"] == "paired") & (master["truth_label"] == "TP")]["HP_Ratio"].apply(lambda x: abs(x - 0.5))
    if len(to_tp_imb) > 0 and len(pa_tp_imb) > 0:
        print(f"S3.3  TO TP |HP-0.5|={to_tp_imb.mean():.4f} vs "
              f"Paired TP={pa_tp_imb.mean():.4f} "
              f"(ratio={to_tp_imb.mean()/pa_tp_imb.mean():.2f}x)")

    print(f"\nAll outputs saved to: {OUT_DIR}")
    print("DONE.")


if __name__ == "__main__":
    main()
