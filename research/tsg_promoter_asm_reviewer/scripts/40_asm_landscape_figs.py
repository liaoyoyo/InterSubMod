#!/usr/bin/env python3
"""
40 — A6 figures part 1: Fig1 genome-wide ASM landscape (Manhattan) + Fig4 coverage vs |Δβ|
(regression-to-extreme). Depends ONLY on asm_dualaxis_tp.tsv (no BAM, no regime-A gates) =
zero rework risk. Fig2 (gate funnel) + Fig3 (IGV triptych) are built separately after the
residual-controls job returns.

Single-sample HCC1395, HP-axis. Output: figures/fig1_asm_landscape.png, figures/fig4_regression_to_extreme.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts")
from lib.plot_setup import setup_plot_style
setup_plot_style()  # per-glyph CJK fallback chain (DejaVu->Noto CJK JP->Droid); fixes tofu

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
TPP = WS / "genome_survey_v2/asm_dualaxis_tp.tsv"
FIG1 = WS / "figures/fig1_asm_landscape.png"
FIG4 = WS / "figures/fig4_regression_to_extreme.png"
HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
BLUE, RED, GREEN, ORANGE, GREY = "#1E3A8A", "#C2410C", "#15803D", "#D97757", "#9CA3AF"


def load():
    d = pd.read_csv(TPP, sep="\t").dropna(subset=["wilcoxon_p", "mean_delta"])
    d = d[d.axis.isin(HP_AXES)].copy()
    d["abs_delta"] = d.mean_delta.abs()
    d["extremity"] = (d.mean_germ_beta - 0.5).abs()
    d["sig"] = d.wilcoxon_p < 0.05
    d["neglogp"] = -np.log10(d.wilcoxon_p.clip(lower=1e-300))
    d = d[d.chrom.isin(CHROM_ORDER)].copy()
    return d


def regime(row):
    hi = row.n_paired_cpg >= 100; ex = row.extremity > 0.3; loh = row.loh_status == "LOH"
    if hi and not ex and not loh:
        return "A_credible"
    if (not hi) and ex:
        return "B_regression"
    if loh:
        return "C_LOH"
    return "D_inter"


def fig1_manhattan(d):
    sig = d[d.sig].copy()
    sig["regime"] = sig.apply(regime, axis=1)
    # cumulative x
    offset = {}; cum = 0; ticks = []; ticklabels = []
    for c in CHROM_ORDER:
        sub = sig[sig.chrom == c]
        offset[c] = cum
        n = sub.somatic_pos.max() if len(sub) else 1e6
        ticks.append(cum + n / 2); ticklabels.append(c.replace("chr", ""))
        cum += (n if n > 0 else 1e6) * 1.02
    sig["x"] = sig.apply(lambda r: offset[r.chrom] + r.somatic_pos, axis=1)

    fig, ax = plt.subplots(figsize=(16, 6), dpi=140)
    # background all-sig in alternating grey
    for i, c in enumerate(CHROM_ORDER):
        sub = sig[sig.chrom == c]
        ax.scatter(sub.x, sub.neglogp, s=6, c=(GREY if i % 2 == 0 else "#C9CDD4"), alpha=0.5, edgecolors="none")
    # regime-A credible overlaid green
    A = sig[sig.regime == "A_credible"]
    ax.scatter(A.x, A.neglogp, s=22, c=GREEN, alpha=0.9, edgecolors="white", linewidths=0.3,
               label=f"regime-A credible (n={len(A)})", zorder=5)
    # BRCA2 annotate
    b = sig[(sig.chrom == "chr13") & (sig.somatic_pos == 32315128)]
    if len(b):
        bx, by = float(b.x.iloc[0]), float(b.neglogp.iloc[0])
        ax.scatter([bx], [by], s=90, marker="*", c=ORANGE, edgecolors="black", linewidths=0.5, zorder=6)
        ax.annotate("BRCA2/ZAR1L\n(rank 25)", (bx, by), textcoords="offset points", xytext=(8, 8),
                    fontsize=10, color=ORANGE, fontweight="bold")
    ax.set_xticks(ticks); ax.set_xticklabels(ticklabels, fontsize=8)
    ax.set_ylabel("−log10 Wilcoxon p", fontsize=12)
    ax.set_xlabel("染色體 (genome position)", fontsize=12)
    ax.set_title("Fig1 · 全基因組 HP-axis ASM landscape (HCC1395 單樣本, sig p<0.05)\n"
                 "灰=全顯著位點 (artifact-dominated); 綠=regime-A credible 子集 (高cov+中baseline+nonLOH)",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="upper right", fontsize=10)
    ax.set_ylim(bottom=-np.log10(0.05) - 0.3)
    ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.6)
    plt.tight_layout(); plt.savefig(FIG1, facecolor="white", bbox_inches="tight")
    print(f"Fig1 -> {FIG1} ({FIG1.stat().st_size//1024} KB), n_sig={len(sig)}, regimeA={len(A)}")


def fig4_regression(d):
    sig = d[d.sig].copy()
    fig, ax = plt.subplots(figsize=(9, 7), dpi=140)
    nonloh = sig[sig.loh_status != "LOH"]; loh = sig[sig.loh_status == "LOH"]
    ax.scatter(nonloh.n_paired_cpg, nonloh.abs_delta, s=10, c=BLUE, alpha=0.35, label=f"nonLOH (n={len(nonloh)})", edgecolors="none")
    ax.scatter(loh.n_paired_cpg, loh.abs_delta, s=14, c=RED, alpha=0.45, label=f"LOH (n={len(loh)})", edgecolors="none")
    ax.axvline(100, color=GREEN, ls="--", lw=1.2, alpha=0.8)
    ax.text(105, ax.get_ylim()[1] * 0.93, "n_cpg=100\n(credible 門檻)", fontsize=9, color=GREEN)
    # shade regression-to-extreme zone (low cov + high |Δβ|)
    ax.axvspan(5, 100, ymin=0.0, ymax=1.0, color=RED, alpha=0.04)
    ax.set_xscale("log"); ax.set_xlim(5, 600)
    ax.set_xlabel("n_paired_cpg (coverage, log)", fontsize=12)
    ax.set_ylabel("|Δβ| (somatic − germline, HP-axis)", fontsize=12)
    ax.set_title("Fig4 · regression-to-extreme: 低 coverage 放大 |Δβ|\n"
                 "LOH(紅) 集中在低 cov 高 |Δβ| 區 = baseline 壓到極端→假性大效應",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10, loc="upper right")
    plt.tight_layout(); plt.savefig(FIG4, facecolor="white", bbox_inches="tight")
    print(f"Fig4 -> {FIG4} ({FIG4.stat().st_size//1024} KB)")


def main():
    d = load()
    print(f"loaded HP-axis TP loci: n={len(d)}, sig={int(d.sig.sum())}")
    fig1_manhattan(d)
    fig4_regression(d)


if __name__ == "__main__":
    main()
