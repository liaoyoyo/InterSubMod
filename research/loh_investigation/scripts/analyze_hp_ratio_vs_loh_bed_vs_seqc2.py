#!/usr/bin/env python3
"""
HP_Ratio (ISM site-level) vs LOH.bed (region-level) vs SEQC2 (ground truth) comparison.

Analyses:
  A. Three-way Venn for HCC1395 (ISM vs LOH.bed vs SEQC2)
  B. Two-way comparison for all 7 TO samples (ISM vs LOH.bed)
  C. TP/FP rates by LOH classification category
  D. HP_Ratio distribution by zone and truth label
  E. Chromosomal position plots showing clustering patterns
  F. Discordant site characteristics

Output: figures + TSV + report  →  research/loh_investigation/
"""

import sys, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pathlib import Path

warnings.filterwarnings("ignore")
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
    "font.size": 10, "axes.titlesize": 12, "figure.facecolor": "white"
})

# === Paths ===
MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
SEQC2_VARIANT = "/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/data/seqc2_vs_to_validation/hcc1395_variant_loh_zone.tsv"
SEQC2_TP_BED = "/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/data/seqc2_vs_to_validation/loh_classified_TP.bed"
SEQC2_FP_BED = "/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/data/seqc2_vs_to_validation/loh_classified_FP.bed"
SEQC2_FN_BED = "/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/data/seqc2_vs_to_validation/loh_classified_FN.bed"
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation")
FIG_DIR = OUT_DIR / "figures"
DATA_DIR = OUT_DIR / "data"
REPORT_DIR = OUT_DIR / "reports"
FIG_DIR.mkdir(parents=True, exist_ok=True)

PREFIX = "loh3way"

# === Human-readable chromosome order ===
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHR_SIZES_HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]


def load_master_to():
    """Load master dataset, TO mode only."""
    print("[1/8] Loading master dataset (TO only)...")
    cols_needed = [
        "RegionID", "Chr", "Pos", "NumReads", "HP_Ratio", "Potential_LOH",
        "LOH_Subtype", "VerificationClass", "Quality_Score", "truth_label",
        "sample", "mode", "effective_hp_reads", "core_loh_like",
        "to_loh_bed_hit", "caller_af", "caller_gq", "CramersV",
        "HPMergedDelta", "AlleleDelta", "hp_ratio_core"
    ]
    df = pd.read_csv(MASTER, sep="\t", usecols=cols_needed)
    df = df[df["mode"] == "to"].copy()
    df["HP_Ratio"] = pd.to_numeric(df["HP_Ratio"], errors="coerce")
    df["Potential_LOH"] = df["Potential_LOH"].astype(str).str.lower().isin(["true", "1"])
    df["to_loh_bed_hit"] = df["to_loh_bed_hit"].astype(str).str.lower().isin(["true", "1"])
    df["core_loh_like"] = df["core_loh_like"].astype(str).str.lower().isin(["true", "1"])
    df["is_tp"] = df["truth_label"] == "TP"
    df["is_fp"] = df["truth_label"] == "FP"
    print(f"   Loaded {len(df):,} TO rows, {df['sample'].nunique()} samples")
    return df


def load_seqc2_variant():
    """Load SEQC2 variant-level LOH zone data for HCC1395."""
    print("[2/8] Loading SEQC2 variant-level data...")
    df = pd.read_csv(SEQC2_VARIANT, sep="\t")
    df["Potential_LOH"] = df["Potential_LOH"].astype(str).str.lower().isin(["true", "1"])
    # loh_zone categories: both_LOH, SEQC2_only, TO_only_in_gain_loss, TO_only_novel, neither
    df["seqc2_loh"] = df["loh_zone"].isin(["both_LOH", "SEQC2_only"])
    df["to_loh_bed"] = df["loh_zone"].isin(["both_LOH", "TO_only_in_gain_loss", "TO_only_novel"])
    df["ism_loh"] = df["Potential_LOH"]
    df["is_tp"] = df["truth_label"] == "TP"
    df["is_fp"] = df["truth_label"] == "FP"
    print(f"   Loaded {len(df):,} HCC1395 variant rows")
    return df


def load_seqc2_beds():
    """Load SEQC2 classified BED regions."""
    beds = {}
    for label, path in [("TP", SEQC2_TP_BED), ("FP", SEQC2_FP_BED), ("FN", SEQC2_FN_BED)]:
        if os.path.exists(path):
            bed = pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end"])
            bed["size_mb"] = (bed["end"] - bed["start"]) / 1e6
            beds[label] = bed
    return beds


# ==============================================================================
# Analysis A: Three-Way Venn for HCC1395
# ==============================================================================
def analysis_a_threeway_venn(seqc2_df):
    """Three-way classification: ISM Potential_LOH vs LOH.bed vs SEQC2."""
    print("[3/8] Analysis A: Three-way Venn (HCC1395)...")

    # Build 3-bit category
    seqc2_df = seqc2_df.copy()
    seqc2_df["cat3"] = (
        seqc2_df["ism_loh"].astype(int) * 4 +
        seqc2_df["to_loh_bed"].astype(int) * 2 +
        seqc2_df["seqc2_loh"].astype(int)
    )
    labels = {
        0: "None",
        1: "SEQC2 only",
        2: "LOH.bed only",
        3: "LOH.bed ∩ SEQC2",
        4: "ISM only",
        5: "ISM ∩ SEQC2",
        6: "ISM ∩ LOH.bed",
        7: "All three"
    }
    seqc2_df["cat3_label"] = seqc2_df["cat3"].map(labels)

    # Count and TP/FP rates
    summary = []
    for cat_val in sorted(labels.keys()):
        sub = seqc2_df[seqc2_df["cat3"] == cat_val]
        n = len(sub)
        n_tp = sub["is_tp"].sum()
        n_fp = sub["is_fp"].sum()
        tp_rate = n_tp / n if n > 0 else 0
        fp_rate = n_fp / n if n > 0 else 0
        summary.append({
            "category": labels[cat_val],
            "cat_code": cat_val,
            "n_sites": n,
            "n_TP": n_tp,
            "n_FP": n_fp,
            "TP_rate": tp_rate,
            "FP_rate": fp_rate,
            "pct_of_total": n / len(seqc2_df) * 100
        })
    sum_df = pd.DataFrame(summary)
    sum_df.to_csv(DATA_DIR / f"{PREFIX}_threeway_venn_hcc1395.tsv", sep="\t", index=False)

    # --- Figure A1: Stacked bar of category counts ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Left: counts
    ax = axes[0]
    ordered = sum_df.sort_values("n_sites", ascending=True)
    colors = plt.cm.Set2(np.linspace(0, 1, 8))
    bars = ax.barh(range(len(ordered)), ordered["n_sites"], color=colors)
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels(ordered["category"], fontsize=9)
    for i, (_, row) in enumerate(ordered.iterrows()):
        ax.text(row["n_sites"] + 50, i, f'{row["n_sites"]:,} ({row["pct_of_total"]:.1f}%)',
                va="center", fontsize=8)
    ax.set_xlabel("Number of SNV sites")
    ax.set_title("A1: Three-Way LOH Classification\n(HCC1395 TO, ISM vs LOH.bed vs SEQC2)")

    # Right: TP/FP rates
    ax = axes[1]
    x = np.arange(len(sum_df))
    w = 0.35
    ordered2 = sum_df[sum_df["n_sites"] > 10].sort_values("cat_code")
    x2 = np.arange(len(ordered2))
    ax.bar(x2 - w/2, ordered2["TP_rate"], w, label="TP rate", color="#2196F3", alpha=0.8)
    ax.bar(x2 + w/2, ordered2["FP_rate"], w, label="FP rate", color="#F44336", alpha=0.8)
    ax.set_xticks(x2)
    ax.set_xticklabels(ordered2["category"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Rate")
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.set_title("A2: TP/FP Rate by LOH Category\n(excluding categories with ≤10 sites)")
    for i, (_, row) in enumerate(ordered2.iterrows()):
        ax.text(i - w/2, row["TP_rate"] + 0.02, f'{row["TP_rate"]:.2f}', ha="center", fontsize=7)
        ax.text(i + w/2, row["FP_rate"] + 0.02, f'{row["FP_rate"]:.2f}', ha="center", fontsize=7)

    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig01_threeway_venn.png")
    plt.close(fig)
    print(f"   Saved fig01 + threeway_venn TSV")
    return seqc2_df, sum_df


# ==============================================================================
# Analysis B: Two-Way Comparison Across All Samples
# ==============================================================================
def analysis_b_twoway_all_samples(master_df):
    """ISM Potential_LOH vs LOH.bed hit for all TO samples."""
    print("[4/8] Analysis B: Two-way comparison (all samples)...")

    results = []
    for sample in SAMPLE_ORDER:
        sub = master_df[master_df["sample"] == sample]
        if len(sub) == 0:
            continue
        # 4 categories
        for ism_val, bed_val, cat_label in [
            (True, True, "Both LOH"),
            (True, False, "ISM only"),
            (False, True, "LOH.bed only"),
            (False, False, "Neither"),
        ]:
            mask = (sub["Potential_LOH"] == ism_val) & (sub["to_loh_bed_hit"] == bed_val)
            ss = sub[mask]
            n = len(ss)
            n_tp = ss["is_tp"].sum()
            n_fp = ss["is_fp"].sum()
            results.append({
                "sample": sample,
                "category": cat_label,
                "n_sites": n,
                "n_TP": n_tp,
                "n_FP": n_fp,
                "TP_rate": n_tp / n if n > 0 else 0,
                "FP_rate": n_fp / n if n > 0 else 0,
                "pct_of_sample": n / len(sub) * 100
            })
    res_df = pd.DataFrame(results)
    res_df.to_csv(DATA_DIR / f"{PREFIX}_twoway_all_samples.tsv", sep="\t", index=False)

    # --- Figure B1: Heatmap of category proportions by sample ---
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Left: proportion heatmap
    pivot_pct = res_df.pivot(index="sample", columns="category", values="pct_of_sample")
    pivot_pct = pivot_pct.reindex(SAMPLE_ORDER).reindex(
        columns=["Both LOH", "ISM only", "LOH.bed only", "Neither"])
    ax = axes[0]
    sns.heatmap(pivot_pct, annot=True, fmt=".1f", cmap="YlOrRd", ax=ax,
                cbar_kws={"label": "% of sample"})
    ax.set_title("B1: LOH Category Distribution (%)\n(ISM Potential_LOH vs LOH.bed hit)")
    ax.set_ylabel("")

    # Middle: TP rate heatmap
    pivot_tp = res_df.pivot(index="sample", columns="category", values="TP_rate")
    pivot_tp = pivot_tp.reindex(SAMPLE_ORDER).reindex(
        columns=["Both LOH", "ISM only", "LOH.bed only", "Neither"])
    ax = axes[1]
    sns.heatmap(pivot_tp, annot=True, fmt=".3f", cmap="Blues", ax=ax, vmin=0, vmax=1,
                cbar_kws={"label": "TP rate"})
    ax.set_title("B2: TP Rate by Category × Sample")
    ax.set_ylabel("")

    # Right: FP rate heatmap
    pivot_fp = res_df.pivot(index="sample", columns="category", values="FP_rate")
    pivot_fp = pivot_fp.reindex(SAMPLE_ORDER).reindex(
        columns=["Both LOH", "ISM only", "LOH.bed only", "Neither"])
    ax = axes[2]
    sns.heatmap(pivot_fp, annot=True, fmt=".3f", cmap="Reds", ax=ax, vmin=0, vmax=1,
                cbar_kws={"label": "FP rate"})
    ax.set_title("B3: FP Rate by Category × Sample")
    ax.set_ylabel("")

    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig02_twoway_heatmaps.png")
    plt.close(fig)
    print(f"   Saved fig02 + twoway TSV")
    return res_df


# ==============================================================================
# Analysis C: HP_Ratio Distribution by Zone and Truth
# ==============================================================================
def analysis_c_hp_ratio_distribution(master_df):
    """HP_Ratio distribution by LOH zone and truth label."""
    print("[5/8] Analysis C: HP_Ratio distributions...")

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()

    summary_rows = []
    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        sub = master_df[master_df["sample"] == sample].copy()
        if len(sub) == 0:
            ax.set_visible(False)
            continue

        # 4 groups: TP+LOH.bed, TP+no_LOH.bed, FP+LOH.bed, FP+no_LOH.bed
        groups = [
            ("TP in LOH.bed", sub[(sub["is_tp"]) & (sub["to_loh_bed_hit"])], "#2196F3"),
            ("TP outside", sub[(sub["is_tp"]) & (~sub["to_loh_bed_hit"])], "#64B5F6"),
            ("FP in LOH.bed", sub[(sub["is_fp"]) & (sub["to_loh_bed_hit"])], "#F44336"),
            ("FP outside", sub[(sub["is_fp"]) & (~sub["to_loh_bed_hit"])], "#EF9A9A"),
        ]

        positions = []
        for j, (label, grp, color) in enumerate(groups):
            if len(grp) > 0:
                hp = grp["HP_Ratio"].dropna()
                # Violin
                parts = ax.violinplot([hp.values], positions=[j], showmedians=True, widths=0.7)
                for pc in parts["bodies"]:
                    pc.set_facecolor(color)
                    pc.set_alpha(0.6)
                parts["cmedians"].set_color("black")

                # Stats
                extreme = ((hp < 0.1) | (hp > 0.9)).mean()
                summary_rows.append({
                    "sample": sample, "group": label,
                    "n": len(hp), "median": hp.median(), "mean": hp.mean(),
                    "extreme_rate": extreme
                })

        ax.set_xticks(range(4))
        ax.set_xticklabels(["TP\nin LOH", "TP\noutside", "FP\nin LOH", "FP\noutside"],
                          fontsize=7)
        ax.axhline(0.1, color="gray", ls="--", lw=0.5, alpha=0.5)
        ax.axhline(0.9, color="gray", ls="--", lw=0.5, alpha=0.5)
        ax.set_ylim(-0.05, 1.05)
        ax.set_ylabel("HP_Ratio" if i % 4 == 0 else "")
        n_total = len(sub)
        ax.set_title(f"{sample}\n(n={n_total:,})", fontsize=10)

    # Last panel: legend
    ax = axes[7]
    ax.axis("off")
    legend_elements = [
        mpatches.Patch(facecolor="#2196F3", alpha=0.6, label="TP in LOH.bed zone"),
        mpatches.Patch(facecolor="#64B5F6", alpha=0.6, label="TP outside LOH.bed"),
        mpatches.Patch(facecolor="#F44336", alpha=0.6, label="FP in LOH.bed zone"),
        mpatches.Patch(facecolor="#EF9A9A", alpha=0.6, label="FP outside LOH.bed"),
    ]
    ax.legend(handles=legend_elements, loc="center", fontsize=11)
    ax.text(0.5, 0.15, "Dashed lines = ISM LOH threshold (0.1/0.9)",
            ha="center", va="center", fontsize=9, transform=ax.transAxes, color="gray")

    plt.suptitle("C: HP_Ratio Distribution by LOH.bed Zone × Truth Label (TO mode)",
                 fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig03_hp_ratio_by_zone_truth.png")
    plt.close(fig)

    sum_df = pd.DataFrame(summary_rows)
    sum_df.to_csv(DATA_DIR / f"{PREFIX}_hp_ratio_by_zone_truth.tsv", sep="\t", index=False)
    print(f"   Saved fig03 + hp_ratio TSV")
    return sum_df


# ==============================================================================
# Analysis D: Chromosomal Position Plot (HCC1395)
# ==============================================================================
def analysis_d_chromosome_position(seqc2_df, master_df):
    """Chromosome ideogram showing positions of discordant sites."""
    print("[6/8] Analysis D: Chromosomal position plots...")

    # === D1: HCC1395 three-way (using seqc2_df) — vectorized ===
    fig, ax = plt.subplots(figsize=(18, 10))

    chrs_present = [c for c in CHR_ORDER if c in seqc2_df["Chr"].values]
    y_positions = {c: i for i, c in enumerate(reversed(chrs_present))}

    # Draw chromosome outlines with correct per-chr length
    max_chr_mb = 0
    for chrom, y in y_positions.items():
        size_bp = CHR_SIZES_HG38.get(chrom, 200_000_000)
        size_mb = size_bp / 1e6
        if size_mb > max_chr_mb:
            max_chr_mb = size_mb
        ax.barh(y, size_mb, height=0.6, color="#F0F0F0", edgecolor="#888888", linewidth=0.8)
        # Right-end tick for chromosome boundary
        ax.plot([size_mb, size_mb], [y - 0.3, y + 0.3], color="#888888", linewidth=0.8)
    ax.set_xlim(0, max_chr_mb * 1.04)

    # Build display_cat column
    seqc2_df = seqc2_df.copy()
    seqc2_df["display_cat"] = seqc2_df["loh_zone"].copy()
    ism_only_mask = (seqc2_df["ism_loh"]) & (seqc2_df["loh_zone"] == "neither")
    seqc2_df.loc[ism_only_mask, "display_cat"] = "ISM_only"

    cat_colors = {
        "both_LOH": ("#7B1FA2", 0.3, "Both LOH (ISM+LOH.bed+SEQC2)"),
        "neither": ("#9E9E9E", 0.05, "Neither"),
        "SEQC2_only": ("#FF9800", 0.8, "SEQC2 only (FN of LOH detection)"),
        "TO_only_novel": ("#4CAF50", 0.9, "TO LOH.bed only (novel)"),
        "TO_only_in_gain_loss": ("#8BC34A", 0.7, "TO LOH.bed in gain/loss"),
        "ISM_only": ("#E91E63", 0.6, "ISM Potential_LOH only"),
    }

    # Vectorized scatter: plot each category with all points at once
    plot_order = ["neither", "both_LOH", "ISM_only", "SEQC2_only",
                  "TO_only_novel", "TO_only_in_gain_loss"]

    np.random.seed(42)
    for cat in plot_order:
        sub = seqc2_df[seqc2_df["display_cat"] == cat].copy()
        if len(sub) == 0 or cat not in cat_colors:
            continue
        color, alpha, label = cat_colors[cat]
        sub = sub[sub["Chr"].isin(y_positions)]
        if len(sub) == 0:
            continue
        xs = sub["Pos"].values / 1e6
        ys = sub["Chr"].map(y_positions).values + np.random.uniform(-0.15, 0.15, len(sub))
        size = 8 if cat != "neither" else 2

        # TP and FP separate markers: TP=circle, FP=square
        tp_mask = sub["is_tp"].values
        fp_mask = sub["is_fp"].values
        if tp_mask.sum() > 0:
            ax.scatter(xs[tp_mask], ys[tp_mask], c=color, marker="o",
                      s=size, alpha=alpha, linewidths=0, rasterized=True)
        if fp_mask.sum() > 0:
            ax.scatter(xs[fp_mask], ys[fp_mask], c=color, marker="s",
                      s=size, alpha=alpha, linewidths=0, rasterized=True)

    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()), fontsize=8)
    ax.set_xlabel("Position (Mb)")
    ax.set_title("D1: HCC1395 TO — SNV Positions by LOH Classification\n"
                 "(○=TP, □=FP; colored by LOH category)")

    handles = []
    for cat in ["both_LOH", "ISM_only", "SEQC2_only", "TO_only_novel", "TO_only_in_gain_loss", "neither"]:
        if cat in cat_colors:
            c, a, l = cat_colors[cat]
            handles.append(mpatches.Patch(facecolor=c, alpha=a, label=l))
    handles.append(plt.Line2D([0], [0], marker="o", color="k", linestyle="None",
                              markersize=6, label="TP (somatic)"))
    handles.append(plt.Line2D([0], [0], marker="s", color="k", linestyle="None",
                              markersize=6, label="FP (germline)"))
    ax.legend(handles=handles, loc="upper right", fontsize=7, ncol=2)

    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig04_chr_position_hcc1395.png")
    plt.close(fig)

    # === D2: Discordant site density (ISM-only LOH) per chromosome ===
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # ISM-only LOH density
    ism_only = seqc2_df[seqc2_df["display_cat"] == "ISM_only"]
    chr_counts = ism_only.groupby("Chr").agg(
        n_sites=("Pos", "count"),
        n_tp=("is_tp", "sum"),
        n_fp=("is_fp", "sum")
    ).reindex([c for c in CHR_ORDER if c in ism_only["Chr"].values]).fillna(0)

    ax = axes[0]
    x = np.arange(len(chr_counts))
    w = 0.4
    ax.bar(x - w/2, chr_counts["n_tp"], w, label="TP", color="#2196F3", alpha=0.8)
    ax.bar(x + w/2, chr_counts["n_fp"], w, label="FP", color="#F44336", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(chr_counts.index, rotation=90, fontsize=7)
    ax.set_ylabel("Count")
    ax.set_title("D2a: ISM-only LOH Sites by Chromosome\n(HCC1395 TO, not in LOH.bed or SEQC2)")
    ax.legend()

    # LOH.bed vs SEQC2 only (FP/FN of LOH detection)
    seqc2_only = seqc2_df[seqc2_df["display_cat"] == "SEQC2_only"]
    to_only = seqc2_df[seqc2_df["display_cat"].isin(["TO_only_novel", "TO_only_in_gain_loss"])]
    chr_seqc2 = seqc2_only.groupby("Chr").size()
    chr_to = to_only.groupby("Chr").size()
    all_chrs = sorted(set(chr_seqc2.index) | set(chr_to.index),
                      key=lambda c: CHR_ORDER.index(c) if c in CHR_ORDER else 99)

    ax = axes[1]
    x2 = np.arange(len(all_chrs))
    ax.bar(x2 - w/2, [chr_seqc2.get(c, 0) for c in all_chrs], w,
           label="SEQC2-only (FN)", color="#FF9800", alpha=0.8)
    ax.bar(x2 + w/2, [chr_to.get(c, 0) for c in all_chrs], w,
           label="TO LOH.bed-only (FP)", color="#4CAF50", alpha=0.8)
    ax.set_xticks(x2)
    ax.set_xticklabels(all_chrs, rotation=90, fontsize=7)
    ax.set_ylabel("Count")
    ax.set_title("D2b: LOH Detection Discordance by Chromosome\n(SEQC2 vs LongPhase-TO)")
    ax.legend()

    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig05_chr_density_hcc1395.png")
    plt.close(fig)
    print(f"   Saved fig04 + fig05")
    return ism_only


# ==============================================================================
# Analysis E: ISM-only LOH Characteristics
# ==============================================================================
def analysis_e_discordant_characteristics(master_df):
    """Characteristics of discordant sites (ISM-only vs LOH.bed-only) across samples."""
    print("[7/8] Analysis E: Discordant site characteristics...")

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))

    features_to_compare = [
        ("effective_hp_reads", "Effective HP Reads"),
        ("HP_Ratio", "HP_Ratio"),
        ("caller_af", "Caller AF"),
        ("caller_gq", "Caller GQ"),
    ]

    summary_rows = []
    for fi, (feat, feat_label) in enumerate(features_to_compare):
        ax_row = fi // 4
        ax_col = fi % 4

        for si, sample in enumerate(SAMPLE_ORDER):
            sub = master_df[master_df["sample"] == sample]
            if len(sub) == 0:
                continue

            both = sub[(sub["Potential_LOH"]) & (sub["to_loh_bed_hit"])][feat].dropna()
            ism_only = sub[(sub["Potential_LOH"]) & (~sub["to_loh_bed_hit"])][feat].dropna()
            bed_only = sub[(~sub["Potential_LOH"]) & (sub["to_loh_bed_hit"])][feat].dropna()
            neither = sub[(~sub["Potential_LOH"]) & (~sub["to_loh_bed_hit"])][feat].dropna()

            for cat, vals in [("Both", both), ("ISM only", ism_only),
                              ("LOH.bed only", bed_only), ("Neither", neither)]:
                if len(vals) > 0:
                    summary_rows.append({
                        "sample": sample, "feature": feat,
                        "category": cat, "n": len(vals),
                        "median": vals.median(), "mean": vals.mean(),
                        "std": vals.std(), "q25": vals.quantile(0.25),
                        "q75": vals.quantile(0.75)
                    })

    char_df = pd.DataFrame(summary_rows)
    char_df.to_csv(DATA_DIR / f"{PREFIX}_discordant_characteristics.tsv", sep="\t", index=False)

    # --- Figure E1: Box plot comparison for key features ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    cat_order = ["Both", "ISM only", "LOH.bed only", "Neither"]
    cat_colors = {"Both": "#7B1FA2", "ISM only": "#E91E63",
                  "LOH.bed only": "#4CAF50", "Neither": "#9E9E9E"}

    for fi, (feat, feat_label) in enumerate(features_to_compare):
        ax = axes[fi]
        feat_data = char_df[char_df["feature"] == feat]

        # Group by sample x category
        x_labels = []
        data_to_plot = []
        colors = []
        for si, sample in enumerate(SAMPLE_ORDER):
            for ci, cat in enumerate(cat_order):
                row = feat_data[(feat_data["sample"] == sample) & (feat_data["category"] == cat)]
                if len(row) > 0:
                    x_labels.append(f"{sample[:6]}\n{cat[:4]}")
                    data_to_plot.append(row["median"].values[0])
                    colors.append(cat_colors[cat])

        # Grouped bar - samples as groups, categories as bars within
        n_samples = len(SAMPLE_ORDER)
        n_cats = len(cat_order)
        x = np.arange(n_samples)
        width = 0.2

        for ci, cat in enumerate(cat_order):
            vals = []
            for sample in SAMPLE_ORDER:
                row = feat_data[(feat_data["sample"] == sample) & (feat_data["category"] == cat)]
                vals.append(row["median"].values[0] if len(row) > 0 else 0)
            offset = (ci - 1.5) * width
            bars = ax.bar(x + offset, vals, width, label=cat if fi == 0 else "",
                         color=cat_colors[cat], alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels([s[:8] for s in SAMPLE_ORDER], rotation=45, ha="right", fontsize=7)
        ax.set_ylabel(f"Median {feat_label}")
        ax.set_title(f"E{fi+1}: {feat_label} by LOH Category")
        if fi == 0:
            ax.legend(fontsize=7, loc="upper right")

    plt.suptitle("E: Feature Characteristics of Discordant LOH Sites\n"
                 "(ISM Potential_LOH vs LOH.bed hit, TO mode)", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig06_discordant_features.png")
    plt.close(fig)
    print(f"   Saved fig06 + discordant characteristics TSV")
    return char_df


# ==============================================================================
# Analysis F: Inclusion Relationship and Transition Zone
# ==============================================================================
def analysis_f_inclusion_and_transition(master_df, seqc2_df):
    """Analyze inclusion relationships and boundary/transition zones."""
    print("[8/8] Analysis F: Inclusion relationships...")

    # --- F1: Concordance matrix per sample ---
    concordance_rows = []
    for sample in SAMPLE_ORDER:
        sub = master_df[master_df["sample"] == sample]
        if len(sub) == 0:
            continue
        n = len(sub)
        both = ((sub["Potential_LOH"]) & (sub["to_loh_bed_hit"])).sum()
        ism_only = ((sub["Potential_LOH"]) & (~sub["to_loh_bed_hit"])).sum()
        bed_only = ((~sub["Potential_LOH"]) & (sub["to_loh_bed_hit"])).sum()
        neither = ((~sub["Potential_LOH"]) & (~sub["to_loh_bed_hit"])).sum()

        # Jaccard-like overlap
        union = both + ism_only + bed_only
        jaccard = both / union if union > 0 else 0

        # ISM LOH rate vs LOH.bed rate
        ism_rate = (both + ism_only) / n
        bed_rate = (both + bed_only) / n

        # ISM→LOH.bed precision/recall
        ism_in_bed_recall = both / (both + bed_only) if (both + bed_only) > 0 else 0
        ism_in_bed_precision = both / (both + ism_only) if (both + ism_only) > 0 else 0

        concordance_rows.append({
            "sample": sample, "n_total": n,
            "both_loh": both, "ism_only": ism_only,
            "bed_only": bed_only, "neither": neither,
            "ism_loh_rate": ism_rate, "bed_loh_rate": bed_rate,
            "excess_ism": ism_rate - bed_rate,
            "jaccard": jaccard,
            "ism_precision_vs_bed": ism_in_bed_precision,
            "ism_recall_vs_bed": ism_in_bed_recall
        })
    conc_df = pd.DataFrame(concordance_rows)
    conc_df.to_csv(DATA_DIR / f"{PREFIX}_concordance_matrix.tsv", sep="\t", index=False)

    # --- Figure F1: Concordance summary ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Left: LOH rates comparison
    ax = axes[0]
    x = np.arange(len(conc_df))
    w = 0.3
    ax.bar(x - w/2, conc_df["ism_loh_rate"], w, label="ISM Potential_LOH", color="#E91E63", alpha=0.8)
    ax.bar(x + w/2, conc_df["bed_loh_rate"], w, label="LOH.bed hit", color="#4CAF50", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(conc_df["sample"].values, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("LOH Rate")
    ax.legend(fontsize=8)
    ax.set_title("F1a: LOH Rate by System\n(ISM site-level vs LOH.bed region-level)")
    for i, row in conc_df.iterrows():
        ax.text(i - w/2, row["ism_loh_rate"] + 0.01,
                f'{row["ism_loh_rate"]:.2f}', ha="center", fontsize=7)
        ax.text(i + w/2, row["bed_loh_rate"] + 0.01,
                f'{row["bed_loh_rate"]:.2f}', ha="center", fontsize=7)

    # Middle: Excess ISM LOH
    ax = axes[1]
    colors = ["#E91E63" if v > 0 else "#4CAF50" for v in conc_df["excess_ism"]]
    ax.bar(x, conc_df["excess_ism"], color=colors, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(conc_df["sample"].values, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("ISM Rate − LOH.bed Rate")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_title("F1b: ISM Excess LOH Rate\n(positive = ISM reports more LOH)")
    for i, row in conc_df.iterrows():
        ax.text(i, row["excess_ism"] + 0.005 * (1 if row["excess_ism"] > 0 else -1),
                f'{row["excess_ism"]:.3f}', ha="center", fontsize=7)

    # Right: Jaccard overlap
    ax = axes[2]
    ax.bar(x, conc_df["jaccard"], color="#7B1FA2", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(conc_df["sample"].values, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Jaccard Index")
    ax.set_ylim(0, 1)
    ax.set_title("F1c: ISM ∩ LOH.bed Jaccard Overlap")
    for i, row in conc_df.iterrows():
        ax.text(i, row["jaccard"] + 0.02, f'{row["jaccard"]:.3f}', ha="center", fontsize=7)

    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig07_concordance_summary.png")
    plt.close(fig)

    # --- Figure F2: HP_Ratio transition zone (ISM-only LOH) ---
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        sub = master_df[master_df["sample"] == sample]
        if len(sub) == 0:
            ax.set_visible(False)
            continue

        # ISM-only LOH: HP_Ratio is extreme but NOT in LOH.bed
        ism_only = sub[(sub["Potential_LOH"]) & (~sub["to_loh_bed_hit"])]
        both = sub[(sub["Potential_LOH"]) & (sub["to_loh_bed_hit"])]

        if len(ism_only) > 0:
            hp = ism_only["HP_Ratio"].dropna()
            ehp = ism_only["effective_hp_reads"].dropna()
            ax.scatter(ehp, hp, c=ism_only["is_tp"].map({True: "#2196F3", False: "#F44336"}),
                      alpha=0.3, s=5, label="ISM-only LOH")

        ax.axhline(0.1, color="gray", ls="--", lw=0.5)
        ax.axhline(0.9, color="gray", ls="--", lw=0.5)
        ax.set_xlabel("Effective HP Reads" if i >= 4 else "")
        ax.set_ylabel("HP_Ratio" if i % 4 == 0 else "")
        n_ism_only = len(ism_only) if len(sub) > 0 else 0
        ax.set_title(f"{sample}\nISM-only LOH: {n_ism_only:,}", fontsize=9)

    ax = axes[7]
    ax.axis("off")
    handles = [
        plt.Line2D([0], [0], marker="o", color="#2196F3", linestyle="None", markersize=5, label="TP"),
        plt.Line2D([0], [0], marker="o", color="#F44336", linestyle="None", markersize=5, label="FP"),
    ]
    ax.legend(handles=handles, loc="center", fontsize=11)
    ax.text(0.5, 0.2, "ISM-only LOH = HP_Ratio < 0.1 or > 0.9\nbut NOT in LOH.bed region",
            ha="center", va="center", fontsize=9, transform=ax.transAxes)

    plt.suptitle("F2: ISM-only LOH Sites — HP_Ratio vs Effective HP Reads\n"
                 "(sites where ISM calls LOH but LOH.bed does not)", fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig08_ism_only_loh_scatter.png")
    plt.close(fig)
    print(f"   Saved fig07 + fig08 + concordance TSV")
    return conc_df


# ==============================================================================
# Analysis G: HCC1395 SEQC2 Deep Dive — TP/FP in each LOH zone
# ==============================================================================
def analysis_g_seqc2_deep_dive(seqc2_df):
    """Deep dive into HCC1395 SEQC2 zones with TP/FP analysis."""
    print("[+] Analysis G: HCC1395 SEQC2 deep dive...")

    # Crosstab: loh_zone × truth_label
    ct = pd.crosstab(seqc2_df["loh_zone"], seqc2_df["truth_label"], margins=True)

    # HP_Ratio stats by loh_zone × truth
    hp_stats = seqc2_df.groupby(["loh_zone", "truth_label"]).agg(
        n=("HP_Ratio", "count"),
        hp_median=("HP_Ratio", "median"),
        hp_mean=("HP_Ratio", "mean"),
        extreme_rate=("Potential_LOH", "mean"),
        af_median=("caller_af", "median"),
        ehp_median=("effective_hp_reads", "median")
    ).reset_index()
    hp_stats.to_csv(DATA_DIR / f"{PREFIX}_seqc2_zone_truth_stats.tsv", sep="\t", index=False)

    # --- Figure G1: Comprehensive SEQC2 zone comparison ---
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    # G1a: Site counts by zone × truth
    ax = axes[0, 0]
    zones = ["both_LOH", "neither", "SEQC2_only", "TO_only_novel", "TO_only_in_gain_loss"]
    zone_labels = ["Both LOH", "Neither", "SEQC2\nonly", "TO novel", "TO gain/loss"]
    for zi, zone in enumerate(zones):
        sub = seqc2_df[seqc2_df["loh_zone"] == zone]
        n_tp = sub["is_tp"].sum()
        n_fp = sub["is_fp"].sum()
        ax.bar(zi - 0.2, n_tp, 0.35, color="#2196F3", alpha=0.8)
        ax.bar(zi + 0.2, n_fp, 0.35, color="#F44336", alpha=0.8)
        ax.text(zi - 0.2, n_tp + 20, str(n_tp), ha="center", fontsize=7)
        ax.text(zi + 0.2, n_fp + 20, str(n_fp), ha="center", fontsize=7)
    ax.set_xticks(range(len(zones)))
    ax.set_xticklabels(zone_labels, fontsize=8)
    ax.set_ylabel("Count")
    ax.set_title("G1a: TP/FP Counts by SEQC2 Zone")
    ax.legend(["TP", "FP"], fontsize=8)

    # G1b: HP_Ratio violin by zone
    ax = axes[0, 1]
    zone_data = []
    zone_names = []
    for zone in zones:
        sub = seqc2_df[seqc2_df["loh_zone"] == zone]["HP_Ratio"].dropna()
        if len(sub) > 0:
            zone_data.append(sub.values)
            zone_names.append(zone.replace("_", "\n")[:12])
    if zone_data:
        parts = ax.violinplot(zone_data, showmedians=True)
        ax.set_xticks(range(1, len(zone_names) + 1))
        ax.set_xticklabels(zone_names, fontsize=7)
    ax.axhline(0.1, color="red", ls="--", lw=0.5, alpha=0.5)
    ax.axhline(0.9, color="red", ls="--", lw=0.5, alpha=0.5)
    ax.set_ylabel("HP_Ratio")
    ax.set_title("G1b: HP_Ratio Distribution by SEQC2 Zone")

    # G1c: ISM Potential_LOH rate by zone
    ax = axes[0, 2]
    loh_rates = []
    for zone in zones:
        sub = seqc2_df[seqc2_df["loh_zone"] == zone]
        loh_rates.append(sub["ism_loh"].mean() if len(sub) > 0 else 0)
    ax.bar(range(len(zones)), loh_rates, color="#E91E63", alpha=0.8)
    ax.set_xticks(range(len(zones)))
    ax.set_xticklabels(zone_labels, fontsize=8)
    ax.set_ylabel("ISM Potential_LOH Rate")
    ax.set_title("G1c: ISM LOH Rate by SEQC2 Zone")
    for i, r in enumerate(loh_rates):
        ax.text(i, r + 0.02, f"{r:.3f}", ha="center", fontsize=8)

    # G1d: Caller AF by zone × truth (box)
    ax = axes[1, 0]
    plot_data = []
    labels_combined = []
    colors_combined = []
    for zone in ["both_LOH", "neither"]:
        for truth, color in [("TP", "#2196F3"), ("FP", "#F44336")]:
            sub = seqc2_df[(seqc2_df["loh_zone"] == zone) & (seqc2_df["truth_label"] == truth)]
            af = sub["caller_af"].dropna()
            if len(af) > 0:
                plot_data.append(af.values)
                labels_combined.append(f"{zone[:8]}\n{truth}")
                colors_combined.append(color)
    if plot_data:
        bp = ax.boxplot(plot_data, patch_artist=True, widths=0.6)
        for patch, color in zip(bp["boxes"], colors_combined):
            patch.set_facecolor(color)
            patch.set_alpha(0.5)
        ax.set_xticklabels(labels_combined, fontsize=7)
    ax.set_ylabel("Caller AF")
    ax.set_title("G1d: Caller AF by Zone × Truth")

    # G1e: Effective HP reads by zone
    ax = axes[1, 1]
    ehp_data = []
    ehp_labels = []
    for zone in zones:
        sub = seqc2_df[seqc2_df["loh_zone"] == zone]["effective_hp_reads"]
        sub = pd.to_numeric(sub, errors="coerce").dropna()
        if len(sub) > 0:
            ehp_data.append(sub.values)
            ehp_labels.append(zone.replace("_", "\n")[:12])
    if ehp_data:
        bp = ax.boxplot(ehp_data, patch_artist=True, widths=0.6)
        for patch in bp["boxes"]:
            patch.set_facecolor("#7B1FA2")
            patch.set_alpha(0.5)
        ax.set_xticklabels(ehp_labels, fontsize=7)
    ax.set_ylabel("Effective HP Reads")
    ax.set_title("G1e: Effective HP Reads by SEQC2 Zone")

    # G1f: VerificationClass distribution by zone
    ax = axes[1, 2]
    vc_cats = ["Strong", "Subclone", "Weak", "Noise"]
    zone_vc = {}
    for zone in ["both_LOH", "neither"]:
        sub = seqc2_df[seqc2_df["loh_zone"] == zone]
        vc_counts = sub["VerificationClass"].value_counts()
        zone_vc[zone] = [vc_counts.get(vc, 0) / len(sub) * 100 for vc in vc_cats]

    x = np.arange(len(vc_cats))
    w = 0.35
    if "both_LOH" in zone_vc:
        ax.bar(x - w/2, zone_vc["both_LOH"], w, label="Both LOH", color="#7B1FA2", alpha=0.8)
    if "neither" in zone_vc:
        ax.bar(x + w/2, zone_vc["neither"], w, label="Neither", color="#9E9E9E", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(vc_cats, fontsize=8)
    ax.set_ylabel("Percentage (%)")
    ax.set_title("G1f: VerificationClass by Zone")
    ax.legend(fontsize=8)

    plt.suptitle("G: HCC1395 SEQC2 LOH Zone Deep Dive (TO mode)", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{PREFIX}_fig09_seqc2_deep_dive.png")
    plt.close(fig)
    print(f"   Saved fig09 + seqc2 zone stats TSV")
    return hp_stats


# ==============================================================================
# Main
# ==============================================================================
def main():
    print("=" * 70)
    print("HP_Ratio vs LOH.bed vs SEQC2 — Three-Way LOH Comparison")
    print("=" * 70)

    master_df = load_master_to()
    seqc2_df = load_seqc2_variant()

    # Run all analyses
    seqc2_df_annotated, venn_df = analysis_a_threeway_venn(seqc2_df)
    twoway_df = analysis_b_twoway_all_samples(master_df)
    hp_dist_df = analysis_c_hp_ratio_distribution(master_df)
    ism_only_df = analysis_d_chromosome_position(seqc2_df_annotated, master_df)
    char_df = analysis_e_discordant_characteristics(master_df)
    conc_df = analysis_f_inclusion_and_transition(master_df, seqc2_df_annotated)
    seqc2_stats = analysis_g_seqc2_deep_dive(seqc2_df_annotated)

    # === Print summary ===
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("\n--- A. Three-Way Venn (HCC1395) ---")
    for _, row in venn_df.iterrows():
        print(f"  {row['category']:25s}: {row['n_sites']:6,} sites "
              f"(TP={row['n_TP']:5,}, FP={row['n_FP']:5,}, "
              f"TP rate={row['TP_rate']:.3f}, FP rate={row['FP_rate']:.3f})")

    print("\n--- F. Concordance (ISM vs LOH.bed) ---")
    for _, row in conc_df.iterrows():
        print(f"  {row['sample']:16s}: ISM={row['ism_loh_rate']:.3f} "
              f"LOH.bed={row['bed_loh_rate']:.3f} "
              f"excess={row['excess_ism']:+.3f} "
              f"Jaccard={row['jaccard']:.3f}")

    print(f"\n--- Output ---")
    print(f"  Figures: {FIG_DIR}/{PREFIX}_fig01-09_*.png")
    print(f"  Data:    {DATA_DIR}/{PREFIX}_*.tsv")

    print("\nDone.")


if __name__ == "__main__":
    main()
