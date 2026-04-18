#!/usr/bin/env python3
"""
Generate publication-quality figures for purity-dependent self-phasing analysis.
Output: PNG figures to workspace directory.
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd

OUTPUT_DIR = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260402_phasing_causal_chain")
FIG_DIR = OUTPUT_DIR / "purity_figures"
os.makedirs(FIG_DIR, exist_ok=True)

# Color scheme
C_TP = "#2166ac"
C_FP = "#b2182b"
C_SCAFFOLD = "#4393c3"
C_NOSCAFFOLD = "#d6604d"
C_NEUTRAL = "#878787"
PURITIES = [0.20, 0.40, 0.60, 0.80, 1.00]
PURITY_LABELS = ["20%", "40%", "60%", "80%", "100%"]


def load_data():
    summary = pd.read_csv(OUTPUT_DIR / "purity_self_phasing_vcf_summary.tsv", sep="\t")
    loh = pd.read_csv(OUTPUT_DIR / "purity_loh_proxy_by_af_threshold.tsv", sep="\t")
    af_strat = pd.read_csv(OUTPUT_DIR / "purity_self_phasing_af_stratified.tsv", sep="\t")
    return summary, loh, af_strat


def fig1_tp_scaffold_rate(summary):
    """Figure 1: TP scaffold involvement rate vs purity."""
    fig, ax = plt.subplots(figsize=(8, 5))

    x = summary["purity"].values
    y = summary["tp_scaffold_rate"].values * 100

    ax.plot(x, y, "o-", color=C_TP, linewidth=2.5, markersize=10, zorder=3)
    for xi, yi in zip(x, y):
        ax.annotate(f"{yi:.1f}%", (xi, yi), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=11, fontweight="bold", color=C_TP)

    ax.set_xlabel("Tumor Purity", fontsize=13)
    ax.set_ylabel("TP Somatic in Phasing Scaffold (%)", fontsize=13)
    ax.set_title("Self-Phasing Scaffold Contamination vs Tumor Purity\n(HCC1395 Subsample Series)", fontsize=14)
    ax.set_xticks(PURITIES)
    ax.set_xticklabels(PURITY_LABELS, fontsize=12)
    ax.set_ylim(0, 75)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax.grid(axis="y", alpha=0.3)

    ax.axhspan(0, 35, alpha=0.08, color="green", label="Low contamination zone")
    ax.axhspan(50, 75, alpha=0.08, color="red", label="High contamination zone")
    ax.legend(loc="upper left", fontsize=10)

    fig.tight_layout()
    out = FIG_DIR / "fig1_tp_scaffold_rate_vs_purity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def fig2_scaffold_loh_effect(loh):
    """Figure 2: Scaffold vs non-scaffold LOH-like rate at each purity."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    x_pos = np.arange(len(PURITIES))
    width = 0.35

    scaffold_loh = loh["tp_scaffold_loh09"].values * 100
    noscaffold_loh = loh["tp_noscaffold_loh09"].values * 100

    # Panel A: absolute rates
    bars1 = ax1.bar(x_pos - width/2, scaffold_loh, width, color=C_SCAFFOLD, label="TP in Scaffold", edgecolor="white")
    bars2 = ax1.bar(x_pos + width/2, noscaffold_loh, width, color=C_NOSCAFFOLD, label="TP not in Scaffold", edgecolor="white")

    for bar, val in zip(bars1, scaffold_loh):
        if val > 0.01:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                     f"{val:.2f}%", ha="center", va="bottom", fontsize=9, fontweight="bold", color=C_SCAFFOLD)
    for bar, val in zip(bars2, noscaffold_loh):
        if val > 0.01:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                     f"{val:.2f}%", ha="center", va="bottom", fontsize=9, color=C_NOSCAFFOLD)

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(PURITY_LABELS, fontsize=12)
    ax1.set_xlabel("Tumor Purity", fontsize=13)
    ax1.set_ylabel("LOH-like Rate (AF > 0.9) %", fontsize=13)
    ax1.set_title("A. TP LOH-like Rate: Scaffold vs Non-scaffold", fontsize=13)
    ax1.legend(fontsize=10)
    ax1.grid(axis="y", alpha=0.3)

    # Panel B: ratio (scaffold / non-scaffold)
    ratios = []
    for s, ns in zip(scaffold_loh, noscaffold_loh):
        if ns > 0.001:
            ratios.append(s / ns)
        else:
            ratios.append(np.nan)

    valid = [(i, r) for i, r in enumerate(ratios) if not np.isnan(r)]
    if valid:
        idx, vals = zip(*valid)
        ax2.bar(list(idx), list(vals), color=C_TP, width=0.6, edgecolor="white")
        for i, v in zip(idx, vals):
            ax2.text(i, v + 0.5, f"{v:.1f}×", ha="center", va="bottom",
                     fontsize=12, fontweight="bold", color=C_TP)

    ax2.axhline(y=1, color="gray", linestyle="--", linewidth=1, alpha=0.7, label="No difference (1×)")
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(PURITY_LABELS, fontsize=12)
    ax2.set_xlabel("Tumor Purity", fontsize=13)
    ax2.set_ylabel("Scaffold / Non-scaffold Ratio", fontsize=13)
    ax2.set_title("B. Self-Phasing LOH Amplification Factor", fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(axis="y", alpha=0.3)

    fig.suptitle("Self-Phasing LOH Effect Disappears Below 60% Purity", fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    out = FIG_DIR / "fig2_scaffold_loh_effect_vs_purity.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def fig3_tp_af_distribution(af_strat):
    """Figure 3: TP AF distribution stacked bar across purities."""
    fig, ax = plt.subplots(figsize=(10, 6))

    af_bins = ["0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.5", "0.5-0.8", "0.8-1.0"]
    colors = ["#d73027", "#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4"]

    x_pos = np.arange(len(PURITIES))
    width = 0.6

    for purity_idx, purity in enumerate(PURITIES):
        sub = af_strat[af_strat["purity"] == purity].copy()
        total = sub["n"].sum()
        if total == 0:
            continue
        bottom = 0
        for bin_idx, af_bin in enumerate(af_bins):
            row = sub[sub["af_bin"] == af_bin]
            if len(row) > 0:
                pct = row["n"].values[0] / total * 100
            else:
                pct = 0
            ax.bar(x_pos[purity_idx], pct, width, bottom=bottom,
                   color=colors[bin_idx],
                   label=f"AF {af_bin}" if purity_idx == 0 else "",
                   edgecolor="white", linewidth=0.5)
            # Label significant segments
            if pct > 8:
                ax.text(x_pos[purity_idx], bottom + pct/2, f"{pct:.0f}%",
                        ha="center", va="center", fontsize=9,
                        color="white" if bin_idx in [0, 5] else "black")
            bottom += pct

    # Add expected VAF line markers
    for i, (purity, expected) in enumerate(zip(PURITIES, [0.10, 0.20, 0.30, 0.40, 0.50])):
        ax.annotate(f"E[VAF]={expected:.2f}", (i, 103), ha="center", fontsize=8,
                    color=C_NEUTRAL, fontstyle="italic")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(PURITY_LABELS, fontsize=12)
    ax.set_xlabel("Tumor Purity", fontsize=13)
    ax.set_ylabel("Fraction of TP PASS Variants (%)", fontsize=13)
    ax.set_title("TP Somatic Variant AF Distribution by Tumor Purity\n(HCC1395 Subsample × SEQC2 Truth Set)", fontsize=14)
    ax.set_ylim(0, 110)
    ax.legend(loc="upper left", fontsize=9, ncol=2, title="Allele Frequency Bin")
    ax.grid(axis="y", alpha=0.2)

    fig.tight_layout()
    out = FIG_DIR / "fig3_tp_af_distribution_by_purity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def fig4_loh_enrichment_direction(loh):
    """Figure 4: LOH enrichment direction (FP/TP) vs purity — always FP-enriched."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    x_pos = np.arange(len(PURITIES))

    tp_loh = loh["tp_loh09_rate"].values * 100
    fp_loh = loh["fp_loh09_rate"].values * 100

    # Panel A: TP vs FP LOH rate
    width = 0.35
    ax1.bar(x_pos - width/2, tp_loh, width, color=C_TP, label="TP (SEQC2 truth)", edgecolor="white")
    ax1.bar(x_pos + width/2, fp_loh, width, color=C_FP, label="FP (non-SEQC2)", edgecolor="white")

    for i, (tp, fp) in enumerate(zip(tp_loh, fp_loh)):
        ax1.text(i - width/2, tp + 0.5, f"{tp:.2f}%", ha="center", fontsize=8, color=C_TP)
        ax1.text(i + width/2, fp + 0.5, f"{fp:.1f}%", ha="center", fontsize=8, color=C_FP)

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(PURITY_LABELS, fontsize=12)
    ax1.set_xlabel("Tumor Purity", fontsize=13)
    ax1.set_ylabel("LOH-like Rate (AF > 0.9) %", fontsize=13)
    ax1.set_title("A. VCF AF-based LOH: Always FP-enriched", fontsize=13)
    ax1.legend(fontsize=10)
    ax1.grid(axis="y", alpha=0.3)

    # Panel B: contrast with ISM HP_Ratio (schematic)
    # Show the key conceptual difference
    categories = ["VCF AF\n(allele freq)", "ISM HP_Ratio\n(haplotype balance)"]
    vcf_direction = [1, 1, 1, 1, 1]  # always FP-enriched
    ism_direction = [-1, -1, -1, -1, -1]  # TP-enriched at 100% (from our previous analysis)

    bar_colors_vcf = [C_FP] * 5
    bar_colors_ism = [C_TP] * 5

    x2 = np.arange(5)
    # VCF enrichment (log scale, all >> 1)
    enrichments = loh["loh09_enrichment_FP_over_TP"].values.copy()
    enrichments[enrichments == 0] = np.nan
    enrichments = np.where(np.isnan(enrichments), 100, enrichments)  # placeholder for inf
    log_enr = np.log10(np.clip(enrichments, 1, 2000))

    ax2.bar(x2, log_enr, 0.6, color=C_FP, edgecolor="white", alpha=0.8)
    for i, (e, le) in enumerate(zip(enrichments, log_enr)):
        label = f"{e:.0f}×" if e < 100 else ">100×"
        ax2.text(i, le + 0.05, label, ha="center", fontsize=10, fontweight="bold", color=C_FP)

    ax2.axhline(y=0, color="black", linewidth=1)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(PURITY_LABELS, fontsize=12)
    ax2.set_xlabel("Tumor Purity", fontsize=13)
    ax2.set_ylabel("log₁₀(FP/TP LOH Enrichment)", fontsize=13)
    ax2.set_title("B. VCF-level LOH Enrichment (FP/TP)\nAll purities: FP >> TP", fontsize=13)
    ax2.set_ylim(0, 4)
    ax2.grid(axis="y", alpha=0.3)

    # Add annotation
    ax2.text(0.5, 0.05, "FP LOH rate always >> TP LOH rate\n→ VCF AF measures germline contamination, not self-phasing",
             transform=ax2.transAxes, ha="center", va="bottom", fontsize=9,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.suptitle("VCF AF ≠ ISM HP_Ratio: Different Metrics, Opposite Conclusions", fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    out = FIG_DIR / "fig4_loh_enrichment_direction_vs_purity.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def fig5_phase_block_structure(summary):
    """Figure 5: Phase block count and median size vs purity."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    x = summary["purity"].values
    blocks = summary["n_phase_blocks"].values
    median_sz = summary["block_median_size"].values

    ax1.plot(x, blocks, "s-", color="#7570b3", linewidth=2.5, markersize=10)
    for xi, yi in zip(x, blocks):
        ax1.annotate(f"{yi:,}", (xi, yi), textcoords="offset points",
                     xytext=(0, 12), ha="center", fontsize=11, fontweight="bold", color="#7570b3")
    ax1.set_xlabel("Tumor Purity", fontsize=13)
    ax1.set_ylabel("Number of Phase Blocks", fontsize=13)
    ax1.set_title("A. Phase Block Count", fontsize=13)
    ax1.set_xticks(PURITIES)
    ax1.set_xticklabels(PURITY_LABELS, fontsize=12)
    ax1.grid(axis="y", alpha=0.3)

    ax2.plot(x, median_sz, "D-", color="#d95f02", linewidth=2.5, markersize=10)
    for xi, yi in zip(x, median_sz):
        ax2.annotate(f"{yi:.0f}", (xi, yi), textcoords="offset points",
                     xytext=(0, 12), ha="center", fontsize=11, fontweight="bold", color="#d95f02")
    ax2.set_xlabel("Tumor Purity", fontsize=13)
    ax2.set_ylabel("Median Block Size (variants)", fontsize=13)
    ax2.set_title("B. Median Phase Block Size", fontsize=13)
    ax2.set_xticks(PURITIES)
    ax2.set_xticklabels(PURITY_LABELS, fontsize=12)
    ax2.grid(axis="y", alpha=0.3)

    # Annotate the 100% anomaly
    ax1.annotate("100% purity:\nMore blocks, smaller size\n(somatic variants fragment scaffold)",
                 xy=(1.0, 4072), xytext=(0.7, 3500),
                 arrowprops=dict(arrowstyle="->", color="gray"),
                 fontsize=9, ha="center", color="gray", fontstyle="italic")

    fig.suptitle("Phase Block Structure Changes with Tumor Purity", fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    out = FIG_DIR / "fig5_phase_block_structure_vs_purity.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def fig6_combined_summary(summary, loh):
    """Figure 6: Combined 4-panel summary figure."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    x = summary["purity"].values

    # Panel A: TP scaffold rate
    ax = axes[0, 0]
    y = summary["tp_scaffold_rate"].values * 100
    ax.plot(x, y, "o-", color=C_TP, linewidth=2.5, markersize=10)
    for xi, yi in zip(x, y):
        ax.annotate(f"{yi:.1f}%", (xi, yi), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=10, fontweight="bold", color=C_TP)
    ax.set_ylabel("TP in Scaffold (%)", fontsize=12)
    ax.set_title("A. Scaffold Contamination Rate", fontsize=13, fontweight="bold")
    ax.set_ylim(0, 75)
    ax.grid(axis="y", alpha=0.3)
    ax.axhspan(0, 35, alpha=0.06, color="green")
    ax.axhspan(50, 75, alpha=0.06, color="red")

    # Panel B: Scaffold LOH effect
    ax = axes[0, 1]
    scaffold_loh = loh["tp_scaffold_loh09"].values * 100
    noscaffold_loh = loh["tp_noscaffold_loh09"].values * 100
    x_pos = np.arange(len(PURITIES))
    w = 0.35
    ax.bar(x_pos - w/2, scaffold_loh, w, color=C_SCAFFOLD, label="In Scaffold", edgecolor="white")
    ax.bar(x_pos + w/2, noscaffold_loh, w, color=C_NOSCAFFOLD, label="Not in Scaffold", edgecolor="white")
    for i, v in enumerate(scaffold_loh):
        if v > 0.01:
            ax.text(i - w/2, v + 0.2, f"{v:.2f}%", ha="center", fontsize=8, fontweight="bold", color=C_SCAFFOLD)
    ax.set_ylabel("TP LOH-like Rate (%)", fontsize=12)
    ax.set_title("B. Scaffold LOH Amplification", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    # Panel C: TP AF median vs expected
    ax = axes[1, 0]
    tp_af = summary["tp_af_median"].values
    expected = summary["expected_het_vaf"].values
    ax.plot(x, tp_af, "o-", color=C_TP, linewidth=2.5, markersize=10, label="Observed TP AF median")
    ax.plot(x, expected, "s--", color=C_NEUTRAL, linewidth=2, markersize=8, label="Expected (P/2)")
    for xi, yi, ei in zip(x, tp_af, expected):
        ax.annotate(f"{yi:.3f}", (xi, yi), textcoords="offset points",
                    xytext=(8, -5), fontsize=9, color=C_TP)
    ax.set_ylabel("Allele Frequency", fontsize=12)
    ax.set_title("C. TP AF Dilution with Purity", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, 0.6)

    # Panel D: AF>0.8 fraction
    ax = axes[1, 1]
    # Compute from af_strat data or use loh data
    tp_high_af = loh["tp_loh08_rate"].values * 100
    fp_high_af = loh["fp_loh08_rate"].values * 100
    ax.plot(x, tp_high_af, "o-", color=C_TP, linewidth=2.5, markersize=10, label="TP AF>0.8")
    ax.plot(x, fp_high_af, "s-", color=C_FP, linewidth=2.5, markersize=10, label="FP AF>0.8")
    for xi, yi in zip(x, tp_high_af):
        ax.annotate(f"{yi:.1f}%", (xi, yi), textcoords="offset points",
                    xytext=(8, 5), fontsize=9, color=C_TP)
    ax.set_ylabel("AF > 0.8 Fraction (%)", fontsize=12)
    ax.set_title("D. Extreme AF Fraction (LOH-like proxy)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    for ax in axes.flat:
        ax.set_xticks(PURITIES if ax in [axes[0, 0], axes[1, 0]] else np.arange(len(PURITIES)))
        if ax in [axes[0, 0], axes[1, 0]]:
            ax.set_xticklabels(PURITY_LABELS, fontsize=11)
        else:
            ax.set_xticklabels(PURITY_LABELS, fontsize=11)
        ax.set_xlabel("Tumor Purity", fontsize=12)

    fig.suptitle("Self-Phasing Effect vs Tumor Purity — HCC1395 Subsample Validation\n(5 purity levels × 22 chromosomes × SEQC2 truth set)",
                 fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    out = FIG_DIR / "fig6_combined_purity_summary.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    summary, loh, af_strat = load_data()

    fig1_tp_scaffold_rate(summary)
    fig2_scaffold_loh_effect(loh)
    fig3_tp_af_distribution(af_strat)
    fig4_loh_enrichment_direction(loh)
    fig5_phase_block_structure(summary)
    fig6_combined_summary(summary, loh)

    print(f"\nAll figures saved to: {FIG_DIR}/")


if __name__ == "__main__":
    main()
