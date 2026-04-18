#!/usr/bin/env python3
"""
Step 1: LOH 區域 AF 分布描述
==============================
研究計劃：LOH 區域 Intermediate AF 與亞克隆偵測
目的：建立 baseline — LOH 區域的 variants AF 實際分布如何？
      是否有 intermediate AF 群集？與 non-LOH 比較。

輸入：all_region_rows.tsv.gz (master dataset)
輸出：
  - 圖表：LOH vs non-LOH AF 分布 (per-sample)
  - 圖表：Intermediate AF 比例跨樣本統計
  - 統計摘要 TSV

備註：
  Cell line purity ≈ 1.0 (no normal cell contamination)
  Deletion LOH (CN=1): expected somatic AF = 1.0 or 0.0
  cnLOH (CN=2): expected somatic AF = 0, 0.5, or 1.0
  Intermediate AF → possible subclonal events
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ── Paths ──
MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUTDIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_subclone_af")
FIGDIR = OUTDIR / "figures"
DATADIR = OUTDIR / "data"
FIGDIR.mkdir(parents=True, exist_ok=True)
DATADIR.mkdir(parents=True, exist_ok=True)

# ── Sample order (fixed) ──
SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
SAMPLE_LABELS = {
    "HCC1395": "HCC1395\n(breast)",
    "HCC1395_DORADO": "HCC1395\n(Dorado)",
    "COLO829": "COLO829\n(melanoma)",
    "H1437": "H1437\n(lung)",
    "H2009": "H2009\n(lung)",
    "HCC1937": "HCC1937\n(BRCA1)",
    "HCC1954": "HCC1954\n(breast)",
}

# Cell line purity = 1.0 (pure tumor culture)
PURITY = 1.0

# AF bins for intermediate detection
# "Extreme": AF < 0.1 or AF > 0.9 (expected for clonal LOH at purity=1.0)
# "Near-half": 0.4 <= AF <= 0.6 (expected for cnLOH post-duplication or diploid het)
# "Intermediate": 0.1 <= AF < 0.4 or 0.6 < AF <= 0.9 (unexpected in clonal LOH)
def classify_af(af):
    if af < 0.1 or af > 0.9:
        return "Extreme (0-0.1 / 0.9-1.0)"
    elif 0.4 <= af <= 0.6:
        return "Near-half (0.4-0.6)"
    else:
        return "Intermediate (0.1-0.4 / 0.6-0.9)"


def load_data():
    """Load master dataset, filter TO mode."""
    print("Loading master dataset...")
    cols = [
        "sample", "mode", "truth_label", "caller_af",
        "to_loh_bed_hit", "HPFineNGroups", "Coverage_Multiple",
        "NumReads", "HP_Ratio", "AlleleDelta",
        "Potential_LOH", "LOH_Subtype", "Quality_Score",
    ]
    df = pd.read_csv(MASTER, sep="\t", usecols=cols, low_memory=False)
    print(f"  Total rows: {len(df):,}")

    # Focus on TO mode (LOH.bed comes from LongPhase-TO)
    df_to = df[df["mode"] == "to"].copy()
    print(f"  TO mode rows: {len(df_to):,}")

    # Convert LOH column
    df_to["loh"] = df_to["to_loh_bed_hit"].astype(str).str.lower() == "true"
    df_to["af_class"] = df_to["caller_af"].apply(classify_af)

    return df_to


def fig1_af_histograms(df):
    """Figure 1: Per-sample AF histograms, LOH vs non-LOH, TP vs FP."""
    print("\n=== Figure 1: AF distributions by LOH status ===")
    fig, axes = plt.subplots(len(SAMPLE_ORDER), 2, figsize=(16, 4 * len(SAMPLE_ORDER)))
    fig.suptitle("Step 1: Caller AF Distribution — LOH vs Non-LOH (TO mode)",
                 fontsize=16, fontweight="bold", y=0.995)

    bins = np.linspace(0, 1, 51)

    for i, sample in enumerate(SAMPLE_ORDER):
        ds = df[df["sample"] == sample]
        if len(ds) == 0:
            continue

        for j, (loh_val, loh_label) in enumerate([(True, "LOH"), (False, "Non-LOH")]):
            ax = axes[i, j]
            sub = ds[ds["loh"] == loh_val]

            tp = sub[sub["truth_label"] == "TP"]["caller_af"]
            fp = sub[sub["truth_label"] == "FP"]["caller_af"]

            if len(tp) > 0:
                ax.hist(tp, bins=bins, alpha=0.6, label=f"TP (n={len(tp):,})",
                        color="#2196F3", density=True, edgecolor="none")
            if len(fp) > 0:
                ax.hist(fp, bins=bins, alpha=0.6, label=f"FP (n={len(fp):,})",
                        color="#F44336", density=True, edgecolor="none")

            # Annotate expected AF values for purity=1.0
            if loh_val:
                # LOH: expected AF = 0 or 1 (deletion) or 0/0.5/1 (cnLOH)
                ax.axvline(0.5, color="green", ls="--", alpha=0.5, lw=1.5,
                           label="cnLOH post-dup (0.5)")
                # Shade intermediate region
                ax.axvspan(0.1, 0.4, alpha=0.08, color="orange")
                ax.axvspan(0.6, 0.9, alpha=0.08, color="orange")

            ax.set_title(f"{SAMPLE_LABELS[sample].replace(chr(10), ' ')} — {loh_label}",
                         fontsize=11)
            ax.legend(fontsize=8, loc="upper right")
            ax.set_xlabel("Caller AF")
            ax.set_ylabel("Density")
            ax.set_xlim(0, 1)

    plt.tight_layout(rect=[0, 0, 1, 0.99])
    out = FIGDIR / "01_af_distribution_loh_vs_nonloh.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig2_intermediate_af_proportion(df):
    """Figure 2: Proportion of intermediate AF variants in LOH vs non-LOH."""
    print("\n=== Figure 2: Intermediate AF proportion ===")

    records = []
    for sample in SAMPLE_ORDER:
        ds = df[df["sample"] == sample]
        for loh_val, loh_label in [(True, "LOH"), (False, "Non-LOH")]:
            for truth in ["TP", "FP"]:
                sub = ds[(ds["loh"] == loh_val) & (ds["truth_label"] == truth)]
                n_total = len(sub)
                if n_total == 0:
                    continue
                n_inter = len(sub[sub["af_class"].str.startswith("Intermediate")])
                n_extreme = len(sub[sub["af_class"].str.startswith("Extreme")])
                n_near_half = len(sub[sub["af_class"].str.startswith("Near-half")])

                records.append({
                    "sample": sample,
                    "loh": loh_label,
                    "truth": truth,
                    "n_total": n_total,
                    "n_intermediate": n_inter,
                    "n_extreme": n_extreme,
                    "n_near_half": n_near_half,
                    "pct_intermediate": 100 * n_inter / n_total,
                    "pct_extreme": 100 * n_extreme / n_total,
                    "pct_near_half": 100 * n_near_half / n_total,
                    "median_af": sub["caller_af"].median(),
                    "mean_af": sub["caller_af"].mean(),
                    "std_af": sub["caller_af"].std(),
                })

    stats = pd.DataFrame(records)
    stats.to_csv(DATADIR / "step1_af_class_statistics.tsv", sep="\t", index=False)
    print(f"  Stats saved: {DATADIR / 'step1_af_class_statistics.tsv'}")

    # Plot: grouped bar chart
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    for ax_idx, truth in enumerate(["TP", "FP"]):
        ax = axes[ax_idx]
        sub_stats = stats[stats["truth"] == truth]

        x = np.arange(len(SAMPLE_ORDER))
        width = 0.35

        loh_pct = []
        nonloh_pct = []
        for sample in SAMPLE_ORDER:
            row_loh = sub_stats[(sub_stats["sample"] == sample) & (sub_stats["loh"] == "LOH")]
            row_nonloh = sub_stats[(sub_stats["sample"] == sample) & (sub_stats["loh"] == "Non-LOH")]
            loh_pct.append(row_loh["pct_intermediate"].values[0] if len(row_loh) > 0 else 0)
            nonloh_pct.append(row_nonloh["pct_intermediate"].values[0] if len(row_nonloh) > 0 else 0)

        bars1 = ax.bar(x - width/2, loh_pct, width, label="LOH", color="#FF9800", alpha=0.8)
        bars2 = ax.bar(x + width/2, nonloh_pct, width, label="Non-LOH", color="#607D8B", alpha=0.8)

        # Add value labels
        for bar in bars1:
            if bar.get_height() > 0:
                ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.5,
                        f'{bar.get_height():.1f}%', ha='center', va='bottom', fontsize=8)
        for bar in bars2:
            if bar.get_height() > 0:
                ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.5,
                        f'{bar.get_height():.1f}%', ha='center', va='bottom', fontsize=8)

        ax.set_xlabel("Sample")
        ax.set_ylabel("% Intermediate AF (0.1-0.4 / 0.6-0.9)")
        ax.set_title(f"Intermediate AF Proportion — {truth}", fontsize=14, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([s.replace("_", "\n") for s in SAMPLE_ORDER], fontsize=9)
        ax.legend()
        ax.set_ylim(0, max(max(loh_pct), max(nonloh_pct)) * 1.25 + 1)

    plt.tight_layout()
    out = FIGDIR / "02_intermediate_af_proportion.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    return stats


def fig3_af_fine_grain_loh(df):
    """Figure 3: Fine-grain AF distribution in LOH regions with AF class annotations."""
    print("\n=== Figure 3: Fine-grain LOH AF distribution ===")

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()

    bins = np.linspace(0, 1, 101)  # 0.01 resolution

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        ds = df[(df["sample"] == sample) & (df["loh"] == True)]

        tp = ds[ds["truth_label"] == "TP"]["caller_af"]
        fp = ds[ds["truth_label"] == "FP"]["caller_af"]

        if len(tp) > 0:
            ax.hist(tp, bins=bins, alpha=0.6, label=f"TP (n={len(tp):,})",
                    color="#2196F3", edgecolor="none")
        if len(fp) > 0:
            ax.hist(fp, bins=bins, alpha=0.6, label=f"FP (n={len(fp):,})",
                    color="#F44336", edgecolor="none")

        # Expected values for purity=1.0
        ax.axvline(0.5, color="green", ls="--", alpha=0.7, lw=1.5)

        # Shade intermediate zones
        ax.axvspan(0.1, 0.4, alpha=0.1, color="orange", label="Intermediate zone")
        ax.axvspan(0.6, 0.9, alpha=0.1, color="orange")

        # Count intermediates
        n_tp_inter = len(tp[(tp >= 0.1) & (tp <= 0.4) | (tp >= 0.6) & (tp <= 0.9)]) if len(tp) > 0 else 0
        n_fp_inter = len(fp[(fp >= 0.1) & (fp <= 0.4) | (fp >= 0.6) & (fp <= 0.9)]) if len(fp) > 0 else 0

        ax.set_title(f"{SAMPLE_LABELS[sample].replace(chr(10), ' ')}\n"
                     f"Intermediate: TP={n_tp_inter}, FP={n_fp_inter}",
                     fontsize=10)
        ax.legend(fontsize=7, loc="upper left")
        ax.set_xlabel("Caller AF")
        ax.set_ylabel("Count")
        ax.set_xlim(0, 1)

    # Hide unused subplot
    axes[-1].axis("off")

    fig.suptitle("Step 1: LOH Region AF Distribution (Fine-grain, TO mode)\n"
                 "Orange = intermediate AF zone (0.1-0.4 / 0.6-0.9) — potential subclone signal",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "03_loh_af_fine_grain.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig4_af_vs_coverage_multiple(df):
    """Figure 4: AF vs Coverage_Multiple in LOH regions — check CN relationship."""
    print("\n=== Figure 4: AF vs Coverage_Multiple in LOH ===")

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        ds = df[(df["sample"] == sample) & (df["loh"] == True)]
        ds = ds.dropna(subset=["Coverage_Multiple", "caller_af"])

        if len(ds) == 0:
            ax.text(0.5, 0.5, "No LOH data", transform=ax.transAxes, ha="center")
            continue

        tp = ds[ds["truth_label"] == "TP"]
        fp = ds[ds["truth_label"] == "FP"]

        if len(tp) > 0:
            ax.scatter(tp["caller_af"], tp["Coverage_Multiple"], alpha=0.15, s=5,
                       c="#2196F3", label=f"TP (n={len(tp):,})")
        if len(fp) > 0:
            ax.scatter(fp["caller_af"], fp["Coverage_Multiple"], alpha=0.15, s=5,
                       c="#F44336", label=f"FP (n={len(fp):,})")

        # Reference lines
        ax.axhline(0.5, color="gray", ls=":", alpha=0.5, label="CN=1 (deletion LOH)")
        ax.axhline(1.0, color="gray", ls="--", alpha=0.5, label="CN=2 (cnLOH/diploid)")

        ax.set_title(SAMPLE_LABELS[sample].replace(chr(10), " "), fontsize=10)
        ax.legend(fontsize=7, loc="upper left", markerscale=3)
        ax.set_xlabel("Caller AF")
        ax.set_ylabel("Coverage_Multiple")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, min(ds["Coverage_Multiple"].quantile(0.99) * 1.2, 4))

    axes[-1].axis("off")

    fig.suptitle("Step 1: LOH Regions — AF vs Coverage_Multiple (CN proxy)\n"
                 "Deletion LOH (CM≈0.5) vs cnLOH (CM≈1.0) — different expected AF",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "04_loh_af_vs_coverage_multiple.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig5_af_class_by_cn_tier(df):
    """Figure 5: AF class distribution stratified by CN tier in LOH regions."""
    print("\n=== Figure 5: AF class by CN tier in LOH ===")

    # Define CN tiers based on Coverage_Multiple
    def cn_tier(cm):
        if pd.isna(cm):
            return "Unknown"
        if cm < 0.75:
            return "CN≈1 (deletion)"
        elif cm < 1.25:
            return "CN≈2 (cnLOH/diploid)"
        elif cm < 1.75:
            return "CN≈3 (gain)"
        else:
            return "CN≥4 (high gain)"

    df_loh = df[df["loh"] == True].copy()
    df_loh["cn_tier"] = df_loh["Coverage_Multiple"].apply(cn_tier)

    cn_order = ["CN≈1 (deletion)", "CN≈2 (cnLOH/diploid)", "CN≈3 (gain)", "CN≥4 (high gain)"]

    # Summary table
    records = []
    for sample in SAMPLE_ORDER:
        ds = df_loh[df_loh["sample"] == sample]
        for cn in cn_order:
            sub = ds[ds["cn_tier"] == cn]
            for truth in ["TP", "FP"]:
                sub_t = sub[sub["truth_label"] == truth]
                n = len(sub_t)
                if n == 0:
                    continue
                n_inter = len(sub_t[sub_t["af_class"].str.startswith("Intermediate")])
                records.append({
                    "sample": sample,
                    "cn_tier": cn,
                    "truth": truth,
                    "n_total": n,
                    "n_intermediate": n_inter,
                    "pct_intermediate": 100 * n_inter / n if n > 0 else 0,
                    "median_af": sub_t["caller_af"].median(),
                })

    cn_stats = pd.DataFrame(records)
    cn_stats.to_csv(DATADIR / "step1_loh_af_by_cn_tier.tsv", sep="\t", index=False)
    print(f"  CN-tier stats saved: {DATADIR / 'step1_loh_af_by_cn_tier.tsv'}")

    # Plot: AF histogram per CN tier (aggregated across samples)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    bins = np.linspace(0, 1, 51)

    for ax_idx, cn in enumerate(cn_order):
        ax = axes.flatten()[ax_idx]
        sub = df_loh[df_loh["cn_tier"] == cn]

        tp = sub[sub["truth_label"] == "TP"]["caller_af"]
        fp = sub[sub["truth_label"] == "FP"]["caller_af"]

        if len(tp) > 0:
            ax.hist(tp, bins=bins, alpha=0.6, label=f"TP (n={len(tp):,})",
                    color="#2196F3", density=True, edgecolor="none")
        if len(fp) > 0:
            ax.hist(fp, bins=bins, alpha=0.6, label=f"FP (n={len(fp):,})",
                    color="#F44336", density=True, edgecolor="none")

        # Expected AF annotation
        if "deletion" in cn:
            ax.axvline(1.0, color="purple", ls="--", alpha=0.7, label="Expected (AF=1.0)")
            ax.axvline(0.0, color="purple", ls="--", alpha=0.7)
        elif "cnLOH" in cn:
            ax.axvline(0.5, color="green", ls="--", alpha=0.7, label="Post-dup (AF=0.5)")
            ax.axvline(1.0, color="purple", ls="--", alpha=0.5)

        ax.axvspan(0.1, 0.4, alpha=0.08, color="orange")
        ax.axvspan(0.6, 0.9, alpha=0.08, color="orange")

        n_inter_tp = len(tp[(tp >= 0.1) & (tp <= 0.4) | (tp >= 0.6) & (tp <= 0.9)]) if len(tp) > 0 else 0
        n_inter_fp = len(fp[(fp >= 0.1) & (fp <= 0.4) | (fp >= 0.6) & (fp <= 0.9)]) if len(fp) > 0 else 0

        ax.set_title(f"{cn}\n(Intermediate: TP={n_inter_tp:,}, FP={n_inter_fp:,})",
                     fontsize=12, fontweight="bold")
        ax.legend(fontsize=9)
        ax.set_xlabel("Caller AF")
        ax.set_ylabel("Density")
        ax.set_xlim(0, 1)

    fig.suptitle("Step 1: LOH AF Distribution by CN Tier (Coverage_Multiple proxy)\n"
                 "All 7 samples aggregated — orange = intermediate AF zones",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "05_loh_af_by_cn_tier.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    return cn_stats


def print_summary(stats, cn_stats):
    """Print key statistics summary."""
    print("\n" + "=" * 80)
    print("STEP 1 SUMMARY: LOH AF Distribution Baseline")
    print("=" * 80)

    print("\n--- Overall AF Class Distribution (TO mode) ---")
    for loh in ["LOH", "Non-LOH"]:
        sub = stats[stats["loh"] == loh]
        if len(sub) == 0:
            continue
        print(f"\n  {loh}:")
        for truth in ["TP", "FP"]:
            sub_t = sub[sub["truth"] == truth]
            if len(sub_t) == 0:
                continue
            total = sub_t["n_total"].sum()
            inter = sub_t["n_intermediate"].sum()
            extreme = sub_t["n_extreme"].sum()
            near_half = sub_t["n_near_half"].sum()
            print(f"    {truth}: N={total:,}")
            print(f"      Extreme (0-0.1/0.9-1.0): {extreme:,} ({100*extreme/total:.1f}%)")
            print(f"      Near-half (0.4-0.6):      {near_half:,} ({100*near_half/total:.1f}%)")
            print(f"      Intermediate (0.1-0.4/0.6-0.9): {inter:,} ({100*inter/total:.1f}%)")

    print("\n--- LOH Intermediate AF by CN Tier ---")
    if cn_stats is not None and len(cn_stats) > 0:
        for cn in cn_stats["cn_tier"].unique():
            sub = cn_stats[cn_stats["cn_tier"] == cn]
            for truth in ["TP", "FP"]:
                sub_t = sub[sub["truth"] == truth]
                if len(sub_t) == 0:
                    continue
                total = sub_t["n_total"].sum()
                inter = sub_t["n_intermediate"].sum()
                pct = 100 * inter / total if total > 0 else 0
                print(f"  {cn} × {truth}: {inter:,}/{total:,} ({pct:.1f}%) intermediate")

    print("\n--- Per-Sample LOH Intermediate AF (TP) ---")
    tp_loh = stats[(stats["truth"] == "TP") & (stats["loh"] == "LOH")]
    for _, row in tp_loh.iterrows():
        print(f"  {row['sample']}: {row['pct_intermediate']:.1f}% intermediate "
              f"({row['n_intermediate']:,}/{row['n_total']:,}), "
              f"median AF={row['median_af']:.3f}")


def main():
    df = load_data()

    # Generate figures
    fig1_af_histograms(df)
    stats = fig2_intermediate_af_proportion(df)
    fig3_af_fine_grain_loh(df)
    fig4_af_vs_coverage_multiple(df)
    cn_stats = fig5_af_class_by_cn_tier(df)

    # Summary
    print_summary(stats, cn_stats)

    print(f"\n✓ All figures saved to: {FIGDIR}")
    print(f"✓ All data saved to: {DATADIR}")


if __name__ == "__main__":
    main()
