#!/usr/bin/env python3
"""L2/L3 feasibility check — assess data availability for PMD and Normal baseline.

L2: PMD (Partially Methylated Domains) stratification
  - Hypothesis: Regions within PMDs show higher subclonal methylation heterogeneity
  - Requires: PMD annotation BED file or computation from methylation data
  - Alternative: Use mean methylation level as PMD proxy (PMDs have lower mean)

L3: Normal BAM methylation baseline
  - Hypothesis: Subtracting normal BAM methylation reduces germline ASM confound
  - Requires: Normal BAM with methylation calls for each sample
  - Check: Do we have normal BAMs available?
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
CANONICAL_ROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
FIGURES_DIR = os.path.join(PROJECT_ROOT, "research", "literature_validation", "figures")
DATA_DIR = os.path.join(PROJECT_ROOT, "research", "literature_validation", "data")

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


def load_sig_summary(sample, category):
    run_id = RUN_PATTERNS[sample]
    path = os.path.join(
        CANONICAL_ROOT, sample, "paired_pileup", run_id,
        f"intersubmod_{category}", "significance_summary.csv"
    )
    if os.path.isfile(path):
        df = pd.read_csv(path)
        df["sample"] = sample
        df["category"] = category.upper()
        return df
    return None


def check_normal_bam_availability():
    """Check if normal BAM files exist for each sample."""
    print("\n--- L3 Check: Normal BAM Availability ---")
    # Known paths for normal BAMs
    normal_bam_patterns = [
        "/big8_disk/liaoyoyo2001/{sample}/normal/*.bam",
        "/big8_disk/liaoyoyo2001/{sample}/normal.bam",
        "/big8_disk/liaoyoyo2001/Data/{sample}/normal/*.bam",
        "/big7_disk/liaoyoyo2001/Data/{sample}/*normal*.bam",
    ]

    # Also check config.sh for sample paths
    config_path = os.path.join(PROJECT_ROOT, "scripts", "pipeline", "config.sh")
    if os.path.isfile(config_path):
        print(f"  Pipeline config: {config_path}")
        with open(config_path) as f:
            content = f.read()
        # Look for normal BAM paths
        import re
        normal_matches = re.findall(r'normal.*bam|NORMAL.*BAM', content, re.IGNORECASE)
        if normal_matches:
            print(f"  Found {len(normal_matches)} normal BAM references in config.sh")
            for m in normal_matches[:5]:
                print(f"    {m[:100]}")

    # Check actual filesystem
    import glob
    for sample in SAMPLES:
        found = []
        for pattern in normal_bam_patterns:
            matches = glob.glob(pattern.format(sample=sample))
            found.extend(matches)
        if found:
            print(f"  {sample}: FOUND normal BAM(s) — {found[0]}")
        else:
            # Try broader search
            search_dirs = [
                f"/big8_disk/liaoyoyo2001/{sample}",
                f"/big7_disk/liaoyoyo2001/Data/{sample}",
            ]
            for d in search_dirs:
                if os.path.isdir(d):
                    bams = glob.glob(os.path.join(d, "**/*normal*.bam"), recursive=True)
                    if bams:
                        print(f"  {sample}: FOUND — {bams[0]}")
                        found.extend(bams)
                        break
            if not found:
                print(f"  {sample}: NOT FOUND")


def pmd_proxy_analysis():
    """Use mean methylation level as PMD proxy.

    PMDs (Partially Methylated Domains) are characterized by:
    - Lower mean methylation (0.3-0.7 vs >0.8 for fully methylated)
    - Higher methylation variability (entropy)
    - Often in gene-poor, late-replicating regions

    We use mean methylation per region as a PMD proxy:
    - "PMD-like": mean_methyl 0.3-0.7
    - "Non-PMD-like": mean_methyl > 0.7 or < 0.3
    """
    print("\n" + "=" * 70)
    print("L2: PMD Proxy Analysis (mean methylation as PMD indicator)")
    print("=" * 70, flush=True)

    # Load Part B directional ASM data (has mean methylation per region)
    dir_path = os.path.join(DATA_DIR, "part_b_directional_asm.csv")
    if not os.path.isfile(dir_path):
        print("[ERROR] Part B data not available. Run L1_directional_ASM.py first.")
        return

    df = pd.read_csv(dir_path)

    # Compute overall methylation per region
    df["overall_methyl"] = (df["alt_mean_methyl"] + df["ref_mean_methyl"]) / 2

    # PMD classification
    df["methyl_zone"] = pd.cut(
        df["overall_methyl"],
        bins=[0, 0.3, 0.7, 1.01],
        labels=["Hypomethylated", "PMD-like", "Fully_methylated"]
    )

    # --- L2-1: Methylation variability by zone ---
    print("\n--- L2-1: CpG Variance by Methylation Zone ---")
    print(f"  {'Zone':>20} | {'TP var':>8} {'TP n':>6} | {'FP var':>8} {'FP n':>6} | {'TP |ASM|':>9} {'FP |ASM|':>9}")
    print("  " + "-" * 80)

    for zone in ["Hypomethylated", "PMD-like", "Fully_methylated"]:
        tp = df[(df["category"] == "TP") & (df["methyl_zone"] == zone)]
        fp = df[(df["category"] == "FP") & (df["methyl_zone"] == zone)]

        tp_var = tp["cpg_variance_mean"].mean() if len(tp) > 0 else float("nan")
        fp_var = fp["cpg_variance_mean"].mean() if len(fp) > 0 else float("nan")
        tp_asm = tp["cpg_diff_mean"].abs().mean() if len(tp) > 0 else float("nan")
        fp_asm = fp["cpg_diff_mean"].abs().mean() if len(fp) > 0 else float("nan")

        print(f"  {zone:>20} | {tp_var:8.4f} {len(tp):6d} | {fp_var:8.4f} {len(fp):6d} | "
              f"{tp_asm:9.4f} {fp_asm:9.4f}")

    # --- L2-2: AlleleDelta by methylation zone (from aggregate data) ---
    print("\n--- L2-2: Subclonal Markers by Methylation Zone ---")
    for zone in ["Hypomethylated", "PMD-like", "Fully_methylated"]:
        tp = df[(df["category"] == "TP") & (df["methyl_zone"] == zone)]
        fp = df[(df["category"] == "FP") & (df["methyl_zone"] == zone)]

        if len(tp) > 10 and len(fp) > 10:
            # Signed delta comparison
            tp_sd = tp["signed_delta"]
            fp_sd = fp["signed_delta"]
            u, p = stats.mannwhitneyu(tp_sd, fp_sd, alternative="two-sided")
            print(f"  {zone:>20}: TP_signed={tp_sd.mean():+.4f}(n={len(tp)}), "
                  f"FP_signed={fp_sd.mean():+.4f}(n={len(fp)}), MWU p={p:.2e}")

    # --- L2-3: Literature prediction test ---
    print("\n--- L2-3: PMD Hypothesis Test ---")
    print("  Literature prediction: PMD regions should show HIGHER methylation")
    print("  heterogeneity and stronger subclonal signals.")

    pmd_tp = df[(df["category"] == "TP") & (df["methyl_zone"] == "PMD-like")]
    non_pmd_tp = df[(df["category"] == "TP") & (df["methyl_zone"] == "Fully_methylated")]

    if len(pmd_tp) > 10 and len(non_pmd_tp) > 10:
        # CpG variance comparison
        u, p = stats.mannwhitneyu(
            pmd_tp["cpg_variance_mean"], non_pmd_tp["cpg_variance_mean"],
            alternative="greater"
        )
        print(f"  CpG variance: PMD={pmd_tp['cpg_variance_mean'].mean():.4f} vs "
              f"Non-PMD={non_pmd_tp['cpg_variance_mean'].mean():.4f}, "
              f"MWU(PMD>Non-PMD) p={p:.2e}")

        # High-var CpG fraction
        u2, p2 = stats.mannwhitneyu(
            pmd_tp["cpg_high_var_frac"], non_pmd_tp["cpg_high_var_frac"],
            alternative="greater"
        )
        print(f"  High-var frac: PMD={pmd_tp['cpg_high_var_frac'].mean():.4f} vs "
              f"Non-PMD={non_pmd_tp['cpg_high_var_frac'].mean():.4f}, "
              f"MWU(PMD>Non-PMD) p={p2:.2e}")

    # --- L2-4: Can PMD zone improve TP/FP discrimination? ---
    print("\n--- L2-4: Discrimination Power by Zone ---")
    for zone in ["Hypomethylated", "PMD-like", "Fully_methylated"]:
        tp_abs_asm = df.loc[(df["category"] == "TP") & (df["methyl_zone"] == zone), "cpg_diff_mean"].abs()
        fp_abs_asm = df.loc[(df["category"] == "FP") & (df["methyl_zone"] == zone), "cpg_diff_mean"].abs()
        if len(tp_abs_asm) > 20 and len(fp_abs_asm) > 20:
            u, _ = stats.mannwhitneyu(fp_abs_asm, tp_abs_asm)
            auc = u / (len(fp_abs_asm) * len(tp_abs_asm))
            print(f"  {zone:>20}: |ASM| AUC(FP>TP)={auc:.4f} "
                  f"(n_TP={len(tp_abs_asm)}, n_FP={len(fp_abs_asm)})")

    # --- Generate plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        # Plot 1: CpG variance by methylation zone
        ax = axes[0]
        zones = ["Hypomethylated", "PMD-like", "Fully_methylated"]
        zone_labels = ["Hypo\n(<0.3)", "PMD-like\n(0.3-0.7)", "Fully meth\n(>=0.7)"]
        for i, cat in enumerate(["TP", "FP"]):
            means = []
            for z in zones:
                sub = df[(df["category"] == cat) & (df["methyl_zone"] == z)]
                means.append(sub["cpg_variance_mean"].mean() if len(sub) > 0 else 0)
            x = np.arange(len(zones))
            ax.bar(x + i * 0.35, means, 0.3, label=cat,
                   color=["#2196F3", "#F44336"][i], alpha=0.7)
        ax.set_xticks(np.arange(len(zones)) + 0.175)
        ax.set_xticklabels(zone_labels)
        ax.set_ylabel("Mean CpG Variance")
        ax.set_title("L2-1: CpG Variance by Methylation Zone")
        ax.legend()

        # Plot 2: |ASM| by zone
        ax = axes[1]
        for i, cat in enumerate(["TP", "FP"]):
            means = []
            for z in zones:
                sub = df[(df["category"] == cat) & (df["methyl_zone"] == z)]
                means.append(sub["cpg_diff_mean"].abs().mean() if len(sub) > 0 else 0)
            x = np.arange(len(zones))
            ax.bar(x + i * 0.35, means, 0.3, label=cat,
                   color=["#2196F3", "#F44336"][i], alpha=0.7)
        ax.set_xticks(np.arange(len(zones)) + 0.175)
        ax.set_xticklabels(zone_labels)
        ax.set_ylabel("|ASM| (ALT-REF per CpG)")
        ax.set_title("L2-3: ASM Magnitude by Zone")
        ax.legend()

        # Plot 3: Overall methylation distribution
        ax = axes[2]
        tp_meth = df.loc[df["category"] == "TP", "overall_methyl"].dropna()
        fp_meth = df.loc[df["category"] == "FP", "overall_methyl"].dropna()
        bins = np.linspace(0, 1, 51)
        ax.hist(tp_meth, bins=bins, alpha=0.6, label=f"TP (n={len(tp_meth)})", color="#2196F3", density=True)
        ax.hist(fp_meth, bins=bins, alpha=0.6, label=f"FP (n={len(fp_meth)})", color="#F44336", density=True)
        ax.set_xlabel("Overall Methylation Level")
        ax.set_ylabel("Density")
        ax.set_title("L2: Methylation Level Distribution")
        ax.axvspan(0.3, 0.7, alpha=0.1, color="orange", label="PMD-like zone")
        ax.legend(fontsize=8)

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "05_L2_pmd_proxy.png"), dpi=150)
        plt.close()
        print(f"\n  Saved: figures/05_L2_pmd_proxy.png")

    except ImportError:
        print("[WARN] matplotlib not available")


def main():
    print("=" * 70)
    print("L2/L3 Feasibility Check")
    print("=" * 70, flush=True)

    # L3: Check normal BAM availability
    check_normal_bam_availability()

    # L2: PMD proxy analysis
    pmd_proxy_analysis()

    print("\n" + "=" * 70)
    print("L2/L3 FEASIBILITY CHECK COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
