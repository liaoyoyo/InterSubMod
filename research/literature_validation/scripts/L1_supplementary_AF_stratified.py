#!/usr/bin/env python3
"""L1 supplementary: AF-stratified ASM analysis.

Check if the |ASM| TP<FP finding is independent of AF confound.
Also test whether directional ASM changes with AF bin.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
CANONICAL_ROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
DATA_DIR = os.path.join(PROJECT_ROOT, "research", "literature_validation", "data")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "research", "literature_validation", "figures")

SAMPLES = ["HCC1395", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954", "HCC1395_DORADO"]
RUN_PATTERNS = {
    "HCC1395": "20260314_HCC1395_paired_pileup_pileup_complete_matrix",
    "COLO829": "20260315_COLO829_paired_pileup_pileup_complete_matrix",
    "H1437": "20260315_H1437_paired_pileup_pileup_complete_matrix",
    "H2009": "20260315_H2009_paired_pileup_pileup_complete_matrix",
    "HCC1937": "20260315_HCC1937_paired_pileup_pileup_complete_matrix",
    "HCC1954": "20260315_HCC1954_paired_pileup_pileup_complete_matrix",
    "HCC1395_DORADO": "20260315_HCC1395_DORADO_paired_pileup_pileup_complete_matrix",
}


def load_all_sig_summaries():
    """Load significance_summary.csv for all samples, TP and FP."""
    frames = []
    for sample in SAMPLES:
        run_id = RUN_PATTERNS[sample]
        for cat in ["tp", "fp"]:
            path = os.path.join(
                CANONICAL_ROOT, sample, "paired_pileup", run_id,
                f"intersubmod_{cat}", "significance_summary.csv"
            )
            if not os.path.isfile(path):
                continue
            df = pd.read_csv(path)
            df["sample"] = sample
            df["category"] = cat.upper()
            frames.append(df)
    return pd.concat(frames, ignore_index=True)


def compute_af_proxy(df):
    """Compute AF proxy from read counts.

    AF ≈ n_ALT / (n_ALT + n_REF). We can approximate from HP family counts
    or use a simpler heuristic.
    """
    # HP1FamilyN and HP2FamilyN represent the two haplotype family sizes
    # In paired mode, ALT tends to be on one haplotype
    # A rough AF proxy: min(HP1, HP2) / (HP1 + HP2) for heterozygous sites
    hp1 = pd.to_numeric(df.get("HP1FamilyN", 0), errors="coerce").fillna(0)
    hp2 = pd.to_numeric(df.get("HP2FamilyN", 0), errors="coerce").fillna(0)
    total = hp1 + hp2
    af_proxy = np.where(total > 0, np.minimum(hp1, hp2) / total, np.nan)
    return af_proxy


def main():
    print("=" * 70)
    print("L1 Supplementary: AF-Stratified ASM Analysis")
    print("=" * 70, flush=True)

    df = load_all_sig_summaries()
    print(f"Total regions: {len(df)}", flush=True)

    # Compute AF proxy
    df["af_proxy"] = compute_af_proxy(df)

    # Define AF bins
    af_bins = [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)]

    # --- S1: AlleleDelta stratified by AF ---
    print("\n--- S1: AlleleDelta by AF Bin ---")
    print(f"  {'AF Bin':>12} | {'TP_mean':>8} {'TP_n':>7} | {'FP_mean':>8} {'FP_n':>7} | {'delta':>8} {'p':>10}")
    print("  " + "-" * 75)

    for lo, hi in af_bins:
        mask = (df["af_proxy"] >= lo) & (df["af_proxy"] < hi)
        tp_vals = df.loc[mask & (df["category"] == "TP"), "AlleleDelta"].dropna()
        fp_vals = df.loc[mask & (df["category"] == "FP"), "AlleleDelta"].dropna()
        if len(tp_vals) > 10 and len(fp_vals) > 10:
            _, p = stats.mannwhitneyu(tp_vals, fp_vals, alternative="two-sided")
            diff = tp_vals.mean() - fp_vals.mean()
            print(f"  [{lo:.1f}, {hi:.1f}) | {tp_vals.mean():8.4f} {len(tp_vals):7d} | "
                  f"{fp_vals.mean():8.4f} {len(fp_vals):7d} | {diff:8.4f} {p:10.2e}")
        else:
            print(f"  [{lo:.1f}, {hi:.1f}) | n_TP={len(tp_vals)}, n_FP={len(fp_vals)} (insufficient)")

    # --- S2: AlleleSig rate by AF ---
    print("\n--- S2: AlleleSig Rate by AF Bin ---")
    print(f"  {'AF Bin':>12} | {'TP_sig%':>8} {'TP_n':>7} | {'FP_sig%':>8} {'FP_n':>7}")
    print("  " + "-" * 50)

    for lo, hi in af_bins:
        mask = (df["af_proxy"] >= lo) & (df["af_proxy"] < hi)
        tp_sub = df.loc[mask & (df["category"] == "TP")]
        fp_sub = df.loc[mask & (df["category"] == "FP")]
        if len(tp_sub) > 10 and len(fp_sub) > 10:
            tp_sig = tp_sub["AlleleSig"].mean() * 100
            fp_sig = fp_sub["AlleleSig"].mean() * 100
            print(f"  [{lo:.1f}, {hi:.1f}) | {tp_sig:7.1f}% {len(tp_sub):7d} | "
                  f"{fp_sig:7.1f}% {len(fp_sub):7d}")

    # --- S3: HP_Ratio relationship with ASM ---
    print("\n--- S3: HP_Ratio and AlleleDelta Correlation ---")
    hp_ratio = pd.to_numeric(df.get("HP_Ratio", np.nan), errors="coerce")
    allele_delta = pd.to_numeric(df.get("AlleleDelta", np.nan), errors="coerce")

    valid = hp_ratio.notna() & allele_delta.notna() & (hp_ratio > 0)
    if valid.sum() > 100:
        for cat in ["TP", "FP"]:
            cat_mask = valid & (df["category"] == cat)
            r, p = stats.spearmanr(hp_ratio[cat_mask], allele_delta[cat_mask])
            print(f"  {cat}: Spearman(HP_Ratio, AlleleDelta) r={r:.4f}, p={p:.2e}, n={cat_mask.sum()}")

    # --- S4: Directional ASM from Part B data, stratified ---
    print("\n--- S4: Directional ASM from Part B Data ---")
    dir_path = os.path.join(DATA_DIR, "part_b_directional_asm.csv")
    if os.path.isfile(dir_path):
        dir_df = pd.read_csv(dir_path)

        # Stratify by mean methylation level
        for cat in ["TP", "FP"]:
            sub = dir_df[dir_df["category"] == cat]
            overall_methyl = (sub["alt_mean_methyl"] + sub["ref_mean_methyl"]) / 2

            # Low methylation regions (< 0.3)
            low = sub[overall_methyl < 0.3]
            # Medium (0.3-0.7)
            mid = sub[(overall_methyl >= 0.3) & (overall_methyl < 0.7)]
            # High (>= 0.7)
            high = sub[overall_methyl >= 0.7]

            print(f"\n  {cat} by methylation level:")
            for name, subset in [("Low (<0.3)", low), ("Mid (0.3-0.7)", mid), ("High (>=0.7)", high)]:
                if len(subset) > 10:
                    sd = subset["signed_delta"]
                    alt_pct = (sd > 0).mean() * 100
                    print(f"    {name:15s}: n={len(subset):5d}, mean_delta={sd.mean():+.4f}, "
                          f"ALT>REF={alt_pct:.1f}%")

        # Stratify by n_reads
        print(f"\n  Directional ASM by read count:")
        for cat in ["TP", "FP"]:
            sub = dir_df[dir_df["category"] == cat]
            n_reads = sub["n_alt"] + sub["n_ref"]
            for lo, hi, label in [(4, 10, "4-9"), (10, 20, "10-19"), (20, 50, "20-49"), (50, 999, "50+")]:
                subset = sub[(n_reads >= lo) & (n_reads < hi)]
                if len(subset) > 20:
                    sd = subset["signed_delta"]
                    print(f"    {cat} reads={label:5s}: n={len(subset):5d}, "
                          f"mean_delta={sd.mean():+.4f}, |ASM|={sd.abs().mean():.4f}")

    # --- S5: LOH-stratified directional analysis ---
    print("\n--- S5: Non-LOH Only Analysis ---")
    non_loh = df[df["Potential_LOH"] == False]
    loh = df[df["Potential_LOH"] == True]

    for name, subset in [("Non-LOH", non_loh), ("LOH", loh)]:
        tp_d = subset.loc[subset["category"] == "TP", "AlleleDelta"].dropna()
        fp_d = subset.loc[subset["category"] == "FP", "AlleleDelta"].dropna()
        if len(tp_d) > 10 and len(fp_d) > 10:
            u, p = stats.mannwhitneyu(tp_d, fp_d, alternative="two-sided")
            print(f"  {name:8s}: TP_delta={tp_d.mean():.4f}(n={len(tp_d)}), "
                  f"FP_delta={fp_d.mean():.4f}(n={len(fp_d)}), p={p:.2e}")

    # --- S6: AUC of |AlleleDelta| for TP/FP discrimination ---
    print("\n--- S6: Discrimination Power (AUC) ---")
    for subset_name, subset in [("All", df), ("Non-LOH", non_loh)]:
        tp_abs = subset.loc[subset["category"] == "TP", "AlleleDelta"].abs().dropna()
        fp_abs = subset.loc[subset["category"] == "FP", "AlleleDelta"].abs().dropna()
        if len(tp_abs) > 100 and len(fp_abs) > 100:
            # Simple AUC via Mann-Whitney
            u, _ = stats.mannwhitneyu(fp_abs, tp_abs, alternative="two-sided")
            auc = u / (len(fp_abs) * len(tp_abs))
            print(f"  {subset_name:8s} |AlleleDelta| AUC(FP>TP): {auc:.4f} "
                  f"(n_TP={len(tp_abs)}, n_FP={len(fp_abs)})")

    # --- S7: Generate AF-stratified plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        # Plot 1: AlleleDelta by AF bin
        ax = axes[0]
        af_bin_labels = []
        tp_means = []
        fp_means = []
        for lo, hi in af_bins:
            mask = (df["af_proxy"] >= lo) & (df["af_proxy"] < hi)
            tp_v = df.loc[mask & (df["category"] == "TP"), "AlleleDelta"].dropna()
            fp_v = df.loc[mask & (df["category"] == "FP"), "AlleleDelta"].dropna()
            if len(tp_v) > 10 and len(fp_v) > 10:
                af_bin_labels.append(f"[{lo:.1f},{hi:.1f})")
                tp_means.append(tp_v.mean())
                fp_means.append(fp_v.mean())

        x = np.arange(len(af_bin_labels))
        ax.bar(x - 0.175, tp_means, 0.3, label="TP", color="#2196F3", alpha=0.7)
        ax.bar(x + 0.175, fp_means, 0.3, label="FP", color="#F44336", alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(af_bin_labels, rotation=30)
        ax.set_ylabel("Mean AlleleDelta")
        ax.set_title("S1: AlleleDelta by AF Bin")
        ax.legend()

        # Plot 2: AlleleSig rate by AF
        ax = axes[1]
        tp_sigs = []
        fp_sigs = []
        af_labels2 = []
        for lo, hi in af_bins:
            mask = (df["af_proxy"] >= lo) & (df["af_proxy"] < hi)
            tp_sub = df.loc[mask & (df["category"] == "TP")]
            fp_sub = df.loc[mask & (df["category"] == "FP")]
            if len(tp_sub) > 10 and len(fp_sub) > 10:
                af_labels2.append(f"[{lo:.1f},{hi:.1f})")
                tp_sigs.append(tp_sub["AlleleSig"].mean() * 100)
                fp_sigs.append(fp_sub["AlleleSig"].mean() * 100)

        x = np.arange(len(af_labels2))
        ax.bar(x - 0.175, tp_sigs, 0.3, label="TP", color="#2196F3", alpha=0.7)
        ax.bar(x + 0.175, fp_sigs, 0.3, label="FP", color="#F44336", alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(af_labels2, rotation=30)
        ax.set_ylabel("Allele Sig Rate (%)")
        ax.set_title("S2: AlleleSig Rate by AF Bin")
        ax.legend()

        # Plot 3: Directional ASM by methylation level
        ax = axes[2]
        if os.path.isfile(dir_path):
            dir_df = pd.read_csv(dir_path)
            bins_meth = [(0, 0.3, "Low"), (0.3, 0.7, "Mid"), (0.7, 1.01, "High")]
            tp_bars = []
            fp_bars = []
            labels3 = []
            for lo, hi, label in bins_meth:
                for cat, bars in [("TP", tp_bars), ("FP", fp_bars)]:
                    sub = dir_df[dir_df["category"] == cat]
                    overall = (sub["alt_mean_methyl"] + sub["ref_mean_methyl"]) / 2
                    subset = sub[(overall >= lo) & (overall < hi)]
                    bars.append(subset["signed_delta"].mean() if len(subset) > 5 else 0)
                labels3.append(label)

            x = np.arange(len(labels3))
            ax.bar(x - 0.175, tp_bars, 0.3, label="TP", color="#2196F3", alpha=0.7)
            ax.bar(x + 0.175, fp_bars, 0.3, label="FP", color="#F44336", alpha=0.7)
            ax.set_xticks(x)
            ax.set_xticklabels(labels3)
            ax.set_ylabel("Mean Signed Delta (ALT-REF)")
            ax.set_title("S4: Directional ASM by Methylation Level")
            ax.axhline(y=0, color="gray", linestyle="--")
            ax.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "04_supplementary_af_stratified.png"), dpi=150)
        plt.close()
        print(f"\n  Saved: figures/04_supplementary_af_stratified.png")

    except ImportError:
        print("[WARN] matplotlib not available")

    print("\n" + "=" * 70)
    print("SUPPLEMENTARY ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
