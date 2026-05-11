#!/usr/bin/env python3
"""
Step 5: TP/FP discrimination figures.

Figures produced:
  F1: TO mode TP rate heatmap (LOH_Subtype x AF_class) across CN tiers
  F2: TO mode FP rate heatmap (same axes) -- artifact regions
  F3: Filter scheme performance bar (TO mode): TP_rate, TP_recall, FP_reduction
  F4: Per-sample S1-S6 comparison (paired_full + TO)
  F5: TP vs FP recall ROC-like (scheme operating points)
  F6: Biology module interpretation summary (text-based heatmap)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import DATA_DIR, FIG_DIR, PALETTE, SAMPLE_ORDER


FIG_DIR.mkdir(parents=True, exist_ok=True)

LOH_ORDER = ["None", "LOH_Noise", "LOH_Weak", "LOH_Subclone", "LOH_Strong"]
AF_ORDER = ["Extreme", "Intermediate", "Near-half"]
CN_ORDER = ["T0", "T1", "T2", "T3", "T4"]
CN_LABEL = {
    "T0": "CN<1\n(del/LOH)",
    "T1": "CN~2\n(diploid)",
    "T2": "CN~3\n(gain)",
    "T3": "CN~4",
    "T4": "CN>=5\n(amp)",
}


def set_style() -> None:
    plt.rcParams.update({
        "figure.facecolor": PALETTE["bg"],
        "axes.facecolor": PALETTE["bg"],
        "axes.edgecolor": PALETTE["dark"],
        "axes.labelcolor": PALETTE["dark"],
        "xtick.color": PALETTE["dark"],
        "ytick.color": PALETTE["dark"],
        "text.color": PALETTE["dark"],
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
    })


def fig1_to_tp_heatmap() -> None:
    df = pd.read_csv(DATA_DIR / "tpfp_cube_coarse_TO.tsv", sep="\t")
    fig, axes = plt.subplots(1, 5, figsize=(18, 4.5),
                             gridspec_kw={"wspace": 0.18, "hspace": 0.1})

    for ax, cn_t in zip(axes, CN_ORDER):
        sub = df[df["cn_tier_F"] == cn_t]
        mat = np.full((len(LOH_ORDER), len(AF_ORDER)), np.nan)
        ncount = np.zeros_like(mat, dtype=int)
        for i, loh in enumerate(LOH_ORDER):
            for j, af in enumerate(AF_ORDER):
                cell = sub[(sub["LOH_Subtype"] == loh) & (sub["AF_class"] == af)]
                if not cell.empty and cell["n"].iloc[0] >= 20:
                    mat[i, j] = cell["tp_rate"].iloc[0]
                    ncount[i, j] = int(cell["n"].iloc[0])
        im = ax.imshow(mat, cmap="RdYlGn", vmin=0.3, vmax=1.0, aspect="auto")
        ax.set_xticks(range(len(AF_ORDER)))
        ax.set_xticklabels(AF_ORDER, rotation=30, ha="right")
        ax.set_yticks(range(len(LOH_ORDER)))
        ax.set_yticklabels(LOH_ORDER if cn_t == "T0" else [])
        ax.set_title(CN_LABEL[cn_t], fontsize=11)
        for i in range(len(LOH_ORDER)):
            for j in range(len(AF_ORDER)):
                v = mat[i, j]
                if not np.isnan(v):
                    col = "white" if v < 0.55 else "black"
                    ax.text(j, i, f"{v*100:.0f}%\nn={ncount[i,j]}",
                            ha="center", va="center", fontsize=8, color=col)
    fig.suptitle("F1. HCC1395 TO mode: TP rate by LOH_Subtype x AF_class across CN tiers\n"
                 f"(baseline TP rate = 71.1%)", fontsize=13, fontweight="bold")
    cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="TP rate")
    plt.savefig(FIG_DIR / "fig_v2_1_to_tp_heatmap.png", dpi=150, bbox_inches="tight",
                facecolor=PALETTE["bg"])
    plt.close()
    print(f"[viz] wrote fig_v2_1_to_tp_heatmap.png")


def fig2_to_fp_heatmap() -> None:
    df = pd.read_csv(DATA_DIR / "tpfp_cube_coarse_TO.tsv", sep="\t")
    df["fp_rate"] = df["n_fp"] / df["n"]
    fig, axes = plt.subplots(1, 5, figsize=(18, 4.5),
                             gridspec_kw={"wspace": 0.18})

    for ax, cn_t in zip(axes, CN_ORDER):
        sub = df[df["cn_tier_F"] == cn_t]
        mat = np.full((len(LOH_ORDER), len(AF_ORDER)), np.nan)
        ncount = np.zeros_like(mat, dtype=int)
        for i, loh in enumerate(LOH_ORDER):
            for j, af in enumerate(AF_ORDER):
                cell = sub[(sub["LOH_Subtype"] == loh) & (sub["AF_class"] == af)]
                if not cell.empty and cell["n"].iloc[0] >= 20:
                    mat[i, j] = cell["fp_rate"].iloc[0]
                    ncount[i, j] = int(cell["n"].iloc[0])
        im = ax.imshow(mat, cmap="Reds", vmin=0.0, vmax=0.55, aspect="auto")
        ax.set_xticks(range(len(AF_ORDER)))
        ax.set_xticklabels(AF_ORDER, rotation=30, ha="right")
        ax.set_yticks(range(len(LOH_ORDER)))
        ax.set_yticklabels(LOH_ORDER if cn_t == "T0" else [])
        ax.set_title(CN_LABEL[cn_t], fontsize=11)
        for i in range(len(LOH_ORDER)):
            for j in range(len(AF_ORDER)):
                v = mat[i, j]
                if not np.isnan(v):
                    col = "black" if v < 0.3 else "white"
                    ax.text(j, i, f"{v*100:.0f}%\nn={ncount[i,j]}",
                            ha="center", va="center", fontsize=8, color=col)
    fig.suptitle("F2. HCC1395 TO mode: FP rate (artifact/germline-leak) by biology module\n"
                 "(baseline FP rate = 28.9%)", fontsize=13, fontweight="bold")
    cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="FP rate")
    plt.savefig(FIG_DIR / "fig_v2_2_to_fp_heatmap.png", dpi=150, bbox_inches="tight",
                facecolor=PALETTE["bg"])
    plt.close()
    print(f"[viz] wrote fig_v2_2_to_fp_heatmap.png")


def fig3_filter_scheme_bar() -> None:
    df = pd.read_csv(DATA_DIR / "tpfp_stratified_filter_schemes_TO.tsv", sep="\t")
    order = ["S1_LOH_Strong_Extreme", "S2_Subclonal_LOH_Intermediate",
             "S3_Diploid_Het", "S4_NonLOH_Extreme", "S5_Combo",
             "S6_S1_NG>=3", "S7_Combo_NG>=3"]
    df = df[df["scheme"].isin(order)].set_index("scheme").reindex(order).reset_index()

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8),
                             gridspec_kw={"wspace": 0.3})
    x = np.arange(len(df))
    colors = [PALETTE["green"], PALETTE["blue"], PALETTE["accent"],
              PALETTE["red"], PALETTE["dark"], PALETTE["green"], PALETTE["dark"]]

    ax = axes[0]
    bars = ax.bar(x, df["tp_rate"], color=colors, edgecolor=PALETTE["dark"])
    ax.axhline(0.7107, ls="--", color=PALETTE["grey"], label="baseline 71.1%")
    for b, v, n in zip(bars, df["tp_rate"], df["n_included"]):
        ax.text(b.get_x()+b.get_width()/2, v+0.005, f"{v*100:.1f}%\nn={n}",
                ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("TP rate (TP / (TP+FP))")
    ax.set_xticks(x); ax.set_xticklabels(df["scheme"], rotation=35, ha="right", fontsize=9)
    ax.set_title("TP purity of each module")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    ax.bar(x, df["tp_recall"]*100, color=colors, edgecolor=PALETTE["dark"])
    for xi, v in zip(x, df["tp_recall"]*100):
        ax.text(xi, v+0.3, f"{v:.1f}%", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("TP recall (% of all TP retained)")
    ax.set_xticks(x); ax.set_xticklabels(df["scheme"], rotation=35, ha="right", fontsize=9)
    ax.set_title("Coverage of module (TP recall)")

    ax = axes[2]
    ax.bar(x, df["fp_reduction"]*100, color=colors, edgecolor=PALETTE["dark"])
    for xi, v in zip(x, df["fp_reduction"]*100):
        ax.text(xi, v+0.3, f"{v:.1f}%", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("FP reduction (%)")
    ax.set_xticks(x); ax.set_xticklabels(df["scheme"], rotation=35, ha="right", fontsize=9)
    ax.set_title("FP removal efficiency")

    fig.suptitle("F3. TO mode filter scheme metrics (each row = biology-informed module)",
                 fontsize=13, fontweight="bold")
    plt.savefig(FIG_DIR / "fig_v2_3_filter_scheme_bar.png", dpi=150,
                bbox_inches="tight", facecolor=PALETTE["bg"])
    plt.close()
    print(f"[viz] wrote fig_v2_3_filter_scheme_bar.png")


def fig4_per_sample_schemes() -> None:
    df = pd.read_csv(DATA_DIR / "tpfp_per_sample_schemes.tsv", sep="\t")
    scheme_order = ["S1_LOH_Strong_Extreme", "S2_Subclonal_LOH_Inter",
                    "S3_Diploid_Het", "S4_NonLOH_Extreme", "S5_Combo_S1_S2_S3",
                    "S6_S1_NG>=3"]
    df = df[df["scheme"].isin(scheme_order)]

    samples_present = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    fig, axes = plt.subplots(1, len(samples_present), figsize=(3.2*len(samples_present), 4.2),
                             sharey=True, gridspec_kw={"wspace": 0.08})
    if len(samples_present) == 1:
        axes = [axes]

    for ax, samp in zip(axes, samples_present):
        sub = df[df["sample"] == samp].set_index("scheme").reindex(scheme_order).reset_index()
        x = np.arange(len(scheme_order))
        colors = [PALETTE["green"], PALETTE["blue"], PALETTE["accent"],
                  PALETTE["red"], PALETTE["dark"], PALETTE["green"]]
        bars = ax.bar(x, sub["tp_rate"], color=colors, edgecolor=PALETTE["dark"])
        for b, v, n in zip(bars, sub["tp_rate"], sub["n"]):
            if n > 0:
                ax.text(b.get_x()+b.get_width()/2, v+0.005, f"n={n}",
                        ha="center", va="bottom", fontsize=7)
        ax.set_ylim(0, 1.05)
        ax.set_title(samp, fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels([s.replace("_", "\n") for s in scheme_order],
                           rotation=0, fontsize=7)
        ax.axhline(0.95, ls=":", color=PALETTE["grey"], alpha=0.5)

    axes[0].set_ylabel("TP rate")
    fig.suptitle("F4. Per-sample filter scheme TP rate (paired_full; dotted=95% reference)",
                 fontsize=12, fontweight="bold")
    plt.savefig(FIG_DIR / "fig_v2_4_per_sample_schemes.png", dpi=150,
                bbox_inches="tight", facecolor=PALETTE["bg"])
    plt.close()
    print(f"[viz] wrote fig_v2_4_per_sample_schemes.png")


def fig5_operating_point() -> None:
    """Scheme operating points: x=TP recall, y=TP rate, size=n."""
    to = pd.read_csv(DATA_DIR / "tpfp_stratified_filter_schemes_TO.tsv", sep="\t")
    pf = pd.read_csv(DATA_DIR / "tpfp_stratified_filter_schemes.tsv", sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5),
                             gridspec_kw={"wspace": 0.3})

    for ax, df, title, baseline in [
        (axes[0], to, "HCC1395 TO mode (baseline TP=71.1%)", 0.7107),
        (axes[1], pf, "All paired_full pooled (baseline TP=95.9%)", 0.9593),
    ]:
        df = df[df["n_included"] >= 50].copy()  # skip tiny
        # Size proportional to log n
        sizes = 50 + 150 * np.log10(df["n_included"] + 1)
        colors = [PALETTE["green"] if ("S1" in s or "S3" in s or "S5" in s or "S6" in s or "S7" in s) and "S4" not in s
                  else (PALETTE["red"] if "S4" in s else PALETTE["accent"])
                  for s in df["scheme"]]
        ax.scatter(df["tp_recall"]*100, df["tp_rate"]*100, s=sizes, c=colors,
                   alpha=0.75, edgecolors=PALETTE["dark"], linewidths=1.2)
        for _, r in df.iterrows():
            ax.annotate(r["scheme"].replace("S", "S"), (r["tp_recall"]*100, r["tp_rate"]*100),
                        fontsize=8, xytext=(6, -2), textcoords="offset points")
        ax.axhline(baseline*100, ls="--", color=PALETTE["grey"], label=f"baseline {baseline*100:.1f}%")
        ax.set_xlabel("TP recall (% of all TP retained)")
        ax.set_ylabel("TP rate (%)")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(loc="lower left", fontsize=9)

    fig.suptitle("F5. Filter scheme operating points\n"
                 "(top-left = high purity, low coverage; top-right = ideal)",
                 fontsize=13, fontweight="bold")
    plt.savefig(FIG_DIR / "fig_v2_5_operating_points.png", dpi=150,
                bbox_inches="tight", facecolor=PALETTE["bg"])
    plt.close()
    print(f"[viz] wrote fig_v2_5_operating_points.png")


def fig6_biology_module_interpretation() -> None:
    """Text summary: biology module -> mechanism -> expected TP rate."""
    modules = [
        ("S1 LOH_Strong + Extreme", "Deletion LOH / cnLOH -- one allele purged;\nremaining haplotype at AF~0 or ~1", 90.1, 0.92, 99.75, PALETTE["green"]),
        ("S2 Subclonal LOH + Inter", "Subclonal/weak LOH; admixed tumor populations\nwith partial allele loss", 87.4, 0.66, 99.77, PALETTE["blue"]),
        ("S3 None + Near-half CN~2/3", "Canonical heterozygous somatic in diploid CN;\nbalanced read support", 95.5, 1.27, 99.85, PALETTE["accent"]),
        ("S4 None + Extreme AF", "Non-LOH extreme AF; high germline-leak /\nmapping bias / low-freq artifact risk", 71.1, 75.95, 24.35, PALETTE["red"]),
        ("S5 Combo (S1 or S2 or S3) - S4", "Union of high-confidence biology;\nexcludes S4 high-risk bucket", 91.8, 2.85, 99.37, PALETTE["dark"]),
    ]

    fig, ax = plt.subplots(figsize=(16, 6.5))
    ax.set_xlim(0, 10); ax.set_ylim(0, len(modules)+1)
    ax.axis("off")

    ax.text(0.1, len(modules)+0.6, "Module", fontsize=11, fontweight="bold")
    ax.text(2.2, len(modules)+0.6, "Biological mechanism", fontsize=11, fontweight="bold")
    ax.text(6.2, len(modules)+0.6, "TP rate", fontsize=11, fontweight="bold")
    ax.text(7.3, len(modules)+0.6, "TP recall", fontsize=11, fontweight="bold")
    ax.text(8.5, len(modules)+0.6, "FP reduc.", fontsize=11, fontweight="bold")

    for i, (name, mech, tp, rec, fpred, col) in enumerate(modules):
        y = len(modules) - i - 0.3
        rect = plt.Rectangle((0.05, y-0.15), 9.9, 0.95, facecolor=col, alpha=0.12,
                             edgecolor=col, linewidth=2)
        ax.add_patch(rect)
        ax.text(0.15, y+0.35, name, fontsize=10, fontweight="bold", color=col)
        ax.text(2.2, y+0.35, mech, fontsize=9, verticalalignment="center")
        ax.text(6.2, y+0.35, f"{tp:.1f}%", fontsize=11, fontweight="bold",
                color=PALETTE["green"] if tp >= 90 else (PALETTE["red"] if tp < 75 else PALETTE["dark"]))
        ax.text(7.3, y+0.35, f"{rec:.2f}%", fontsize=10)
        ax.text(8.5, y+0.35, f"{fpred:.2f}%", fontsize=10,
                color=PALETTE["green"] if fpred >= 99 else (PALETTE["red"] if fpred < 50 else PALETTE["dark"]))

    ax.set_title("F6. Biology module interpretation (HCC1395 TO mode)\n"
                 "Each module defines a parameter space with distinct biology and TP/FP expectation",
                 fontsize=13, fontweight="bold", pad=14)
    plt.savefig(FIG_DIR / "fig_v2_6_biology_module_interpretation.png", dpi=150,
                bbox_inches="tight", facecolor=PALETTE["bg"])
    plt.close()
    print(f"[viz] wrote fig_v2_6_biology_module_interpretation.png")


def main() -> None:
    set_style()
    fig1_to_tp_heatmap()
    fig2_to_fp_heatmap()
    fig3_filter_scheme_bar()
    fig4_per_sample_schemes()
    fig5_operating_point()
    fig6_biology_module_interpretation()
    print(f"\n[viz] all figures written to {FIG_DIR}")


if __name__ == "__main__":
    main()
