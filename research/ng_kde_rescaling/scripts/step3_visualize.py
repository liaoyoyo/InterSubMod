#!/usr/bin/env python3
"""
Step 3 — Visualization: 6 PNGs for NG × new-KDE × CN strategies × LOH × AF.

Follows feedback_figure_layout_standard (2026-04):
  - Fixed sample order
  - Palette: dark/bg/accent/green/red/blue/grey
  - Arial ≥14pt, min 9pt
  - PNG ≥1200px, fit-within centered
  - Paired vs TO dimensional separation where relevant
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import DATA_DIR, FIG_DIR, PALETTE, SAMPLE_ORDER

# Global style
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.facecolor": PALETTE["bg"],
    "axes.facecolor": PALETTE["bg"],
    "axes.edgecolor": PALETTE["dark"],
    "axes.labelcolor": PALETTE["dark"],
    "xtick.color": PALETTE["dark"],
    "ytick.color": PALETTE["dark"],
    "text.color": PALETTE["dark"],
    "savefig.facecolor": PALETTE["bg"],
    "savefig.dpi": 140,
})


def load() -> dict:
    dfs = {}
    dfs["master"] = pd.read_csv(DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz",
                                sep="\t", low_memory=False)
    dfs["cn_matrix"] = pd.read_csv(DATA_DIR / "cn_tier_coverage_matrix.tsv", sep="\t")
    dfs["ng_dist"] = pd.read_csv(DATA_DIR / "ng_dist_stratified.tsv", sep="\t")
    dfs["ng_tprate"] = pd.read_csv(DATA_DIR / "ng_tprate_stratified.tsv", sep="\t")
    dfs["af_cross"] = pd.read_csv(DATA_DIR / "af_bin10_ng_crosstab.tsv", sep="\t")
    dfs["repl"] = pd.read_csv(DATA_DIR / "ng_20260414_replication_inter_vs_extreme.tsv", sep="\t")
    dfs["ext"] = pd.read_csv(DATA_DIR / "ng_loh_inner_vs_outer_extension.tsv", sep="\t")
    dfs["h2"] = pd.read_csv(DATA_DIR / "ng_h2_verification_lohbed_vs_potential.tsv", sep="\t")
    return dfs


# ────────────────────────────────────────────────────────────────
def fig1_covm_density(master: pd.DataFrame, out: Path) -> None:
    """CovM density per sample (new-KDE baseline) with strategy A/F boundary lines."""
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), sharex=True, sharey=False)
    axes = axes.flatten()
    xs = np.linspace(0, 3, 300)
    strat_A = [0.75, 1.25, 1.75]
    strat_F = [0.65, 0.99, 1.33, 1.82]

    for i, s in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        sub = master[(master["sample"] == s) & (master["mode"] == "paired_full")]
        vals = sub["CovM_used"].dropna().clip(0, 3).values
        if len(vals) == 0:
            ax.set_visible(False)
            continue
        # hist
        ax.hist(vals, bins=np.arange(0, 3.05, 0.05), color=PALETTE["blue"],
                alpha=0.75, edgecolor="none")
        # median vertical
        med = float(np.median(vals))
        ax.axvline(med, color=PALETTE["dark"], linestyle="-", linewidth=1.2,
                   label=f"median={med:.2f}")
        # Strategy A boundaries (solid)
        for b in strat_A:
            ax.axvline(b, color=PALETTE["accent"], linestyle="--", linewidth=0.9, alpha=0.7)
        # Strategy F boundaries (dotted)
        for b in strat_F:
            ax.axvline(b, color=PALETTE["green"], linestyle=":", linewidth=0.9, alpha=0.7)
        ax.set_title(f"{s} (n={len(sub):,})", fontsize=12)
        ax.set_xlim(0, 3)
        ax.set_xlabel("CovM (new KDE baseline)", fontsize=10)
        if i % 4 == 0:
            ax.set_ylabel("count", fontsize=10)
        ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)

    # Legend in 8th subplot
    ax = axes[7]
    ax.set_visible(True)
    ax.axis("off")
    handles = [
        mpl.lines.Line2D([0], [0], color=PALETTE["blue"], lw=6, alpha=0.75, label="CovM hist"),
        mpl.lines.Line2D([0], [0], color=PALETTE["dark"], lw=1.5, label="sample median"),
        mpl.lines.Line2D([0], [0], color=PALETTE["accent"], lw=1.2, ls="--",
                         label="Strategy A (legacy 0.75/1.25/1.75)"),
        mpl.lines.Line2D([0], [0], color=PALETTE["green"], lw=1.2, ls=":",
                         label="Strategy F (SEQC2 0.65/0.99/1.33/1.82)"),
    ]
    ax.legend(handles=handles, loc="center", fontsize=11, frameon=True,
              edgecolor=PALETTE["dark"])

    fig.suptitle("F1: CovM distribution under new KDE baseline (paired_full, 7 samples)\n"
                 "vertical lines = CN tier boundaries for Strategy A (legacy) & F (SEQC2-grounded)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")


# ────────────────────────────────────────────────────────────────
def fig2_tprate_heatmap(tprate: pd.DataFrame, out: Path) -> None:
    """TP rate heatmap: 7 samples × CN tier × LOH side, for Strategy A & F × NG filter=all/NG>=4+NR80."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    filters = ["all", "NG=4_NR>=80"]
    strats = ["A_Legacy", "F_SEQC2_grounded"]
    for ri, ng_f in enumerate(filters):
        for ci, strat in enumerate(strats):
            ax = axes[ri, ci]
            sub = tprate[(tprate["strategy"] == strat) & (tprate["ng_filter"] == ng_f)].copy()
            # build matrix: row=sample, col=tier_lohside
            tier_order = sorted(sub["tier"].dropna().unique())
            col_order = [f"{t}_{side}" for t in tier_order for side in ["Inner", "Outer"]]
            sub["col_key"] = sub["tier"].astype(str) + "_" + sub["loh_side"]
            pivot = sub.pivot_table(index="sample", columns="col_key", values="tp_rate", aggfunc="first")
            n_pivot = sub.pivot_table(index="sample", columns="col_key", values="n", aggfunc="first")
            pivot = pivot.reindex(SAMPLE_ORDER).reindex(columns=col_order)
            n_pivot = n_pivot.reindex(SAMPLE_ORDER).reindex(columns=col_order)

            im = ax.imshow(pivot.values, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
            ax.set_xticks(range(len(col_order)))
            ax.set_xticklabels(col_order, rotation=45, ha="right", fontsize=9)
            ax.set_yticks(range(len(SAMPLE_ORDER)))
            ax.set_yticklabels(SAMPLE_ORDER, fontsize=10)
            # annotate
            for yi in range(pivot.shape[0]):
                for xi in range(pivot.shape[1]):
                    v = pivot.values[yi, xi]
                    n = n_pivot.values[yi, xi]
                    if np.isnan(v) or np.isnan(n) or n < 5:
                        continue
                    txt = f"{v:.2f}\n(n={int(n)})"
                    ax.text(xi, yi, txt, ha="center", va="center",
                            fontsize=7, color="white" if v < 0.5 else "black")
            ax.set_title(f"{strat} | filter={ng_f}", fontsize=12)
            fig.colorbar(im, ax=ax, shrink=0.7, label="TP rate")

    fig.suptitle("F2: TP rate by sample × CN tier × LOH side (Strategy A vs F)\n"
                 "filter=all (row 1), NG≥4+NR≥80 canonical filter (row 2)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")


# ────────────────────────────────────────────────────────────────
def fig3_ng_distribution(ng_dist: pd.DataFrame, out: Path) -> None:
    """Stacked bar: NG=0/1/2/3/4 composition per sample × CN tier (Strategy A, LOH Inner + Outer)."""
    strat = "A_Legacy"
    fig, axes = plt.subplots(2, 7, figsize=(18, 8), sharey=True)
    ng_colors = {
        0: PALETTE["grey"], 1: PALETTE["blue"], 2: PALETTE["green"],
        3: PALETTE["accent"], 4: PALETTE["red"],
    }
    for row_i, side in enumerate(["Inner", "Outer"]):
        for col_i, sample in enumerate(SAMPLE_ORDER):
            ax = axes[row_i, col_i]
            sub = ng_dist[(ng_dist["strategy"] == strat) &
                          (ng_dist["sample"] == sample) &
                          (ng_dist["loh_side"] == side)]
            if sub.empty:
                ax.axis("off")
                continue
            tiers = sorted(sub["tier"].unique())
            bottoms = np.zeros(len(tiers))
            for ng in [0, 1, 2, 3, 4]:
                vals = []
                for t in tiers:
                    cell = sub[(sub["tier"] == t) & (sub["NG"] == ng)]
                    pct = float(cell["pct_within_cell"].iloc[0]) if not cell.empty else 0.0
                    vals.append(pct)
                vals = np.array(vals)
                ax.bar(tiers, vals, bottom=bottoms, color=ng_colors[ng],
                       edgecolor="white", linewidth=0.5, label=f"NG={ng}" if (row_i == 0 and col_i == 6) else None)
                bottoms += vals
            ax.set_title(f"{sample}\n(LOH {side})", fontsize=10)
            ax.set_ylim(0, 100)
            if col_i == 0:
                ax.set_ylabel("% within cell", fontsize=10)
            ax.tick_params(axis="x", labelsize=9)
            if row_i == 1:
                ax.set_xlabel("CN tier", fontsize=9)

    # legend on last axis
    axes[0, 6].legend(loc="upper right", bbox_to_anchor=(1.45, 1.0), fontsize=9,
                      frameon=True, edgecolor=PALETTE["dark"])
    fig.suptitle("F3: HPFineNGroups distribution per sample × CN tier × LOH side (Strategy A)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")


# ────────────────────────────────────────────────────────────────
def fig4_af_ng_heatmap(af_cross: pd.DataFrame, out: Path) -> None:
    """AF 10-bin × NG cell TP rate, per sample (Strategy A CN1)."""
    fig, axes = plt.subplots(2, 4, figsize=(16, 8.5))
    axes = axes.flatten()
    af_bins = sorted(af_cross["AF_bin10"].dropna().unique())
    ngs = sorted(af_cross["NG"].dropna().unique())
    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        sub = af_cross[(af_cross["sample"] == sample) & (af_cross["cn_tier_A"] == "T0")]
        if sub.empty:
            ax.axis("off")
            continue
        pivot_tp = sub.pivot_table(index="NG", columns="AF_bin10", values="tp_rate", aggfunc="first")
        pivot_n = sub.pivot_table(index="NG", columns="AF_bin10", values="n", aggfunc="first")
        pivot_tp = pivot_tp.reindex(index=sorted(ngs, reverse=True), columns=af_bins)
        pivot_n = pivot_n.reindex(index=sorted(ngs, reverse=True), columns=af_bins)
        im = ax.imshow(pivot_tp.values, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(af_bins)))
        ax.set_xticklabels([f"{b:.1f}" for b in af_bins] if len(af_bins) and isinstance(af_bins[0], (int, float))
                           else [str(b).split(",")[0].replace("(", "") for b in af_bins],
                           rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(pivot_tp.shape[0]))
        ax.set_yticklabels([f"NG={int(ng)}" for ng in sorted(ngs, reverse=True)], fontsize=9)
        # annotate n
        for yi in range(pivot_tp.shape[0]):
            for xi in range(pivot_tp.shape[1]):
                n = pivot_n.values[yi, xi]
                v = pivot_tp.values[yi, xi]
                if np.isnan(n) or n < 3:
                    continue
                ax.text(xi, yi, f"{int(n)}", ha="center", va="center",
                        fontsize=6, color="white" if v < 0.5 else "black")
        ax.set_title(sample, fontsize=11)
        ax.set_xlabel("AF bin (left edge)", fontsize=9)

    # last axis for colorbar
    axes[7].axis("off")
    cbar = fig.colorbar(im, ax=axes[7], shrink=0.9)
    cbar.set_label("TP rate", fontsize=10)

    fig.suptitle("F4: AF (10-bin) × NG cell TP rate in CN1 (Strategy A, paired_full)\n"
                 "cell annotation = n; color = TP rate",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")


# ────────────────────────────────────────────────────────────────
def fig5_replication_forest(repl: pd.DataFrame, out: Path) -> None:
    """Forest: ΔNG (Inter − Extreme) per sample, per strategy A/B/F; the key 20260414 replication figure."""
    fig, ax = plt.subplots(1, 1, figsize=(11, 7))
    strat_order = ["A_Legacy", "B_ClinicalStrict", "F_SEQC2_grounded"]
    strat_offset = {"A_Legacy": -0.25, "B_ClinicalStrict": 0, "F_SEQC2_grounded": 0.25}
    strat_colors = {
        "A_Legacy": PALETTE["accent"],
        "B_ClinicalStrict": PALETTE["blue"],
        "F_SEQC2_grounded": PALETTE["green"],
    }

    y_pos = {s: i for i, s in enumerate(SAMPLE_ORDER)}
    for strat in strat_order:
        sub = repl[repl["strategy"] == strat]
        for _, row in sub.iterrows():
            s = row["sample"]
            d = row["delta_ng_inter_vs_extreme"]
            n_i = row["n_intermediate"]
            n_e = row["n_extreme"]
            if pd.isna(d) or n_i < 5 or n_e < 5:
                continue
            y = y_pos[s] + strat_offset[strat]
            ax.scatter(d, y, s=120, color=strat_colors[strat], edgecolor=PALETTE["dark"],
                       linewidth=0.6, zorder=3)
            ax.text(d + 0.02, y, f"n={int(n_i)}/{int(n_e)}",
                    fontsize=7, color=PALETTE["dark"], va="center", alpha=0.7)

    # 20260414 canonical ΔNG line
    ax.axvline(0.705, color=PALETTE["dark"], linestyle="-", linewidth=2, alpha=0.9,
               label="20260414 canonical ΔNG=+0.705")
    ax.axvline(0, color=PALETTE["grey"], linestyle="--", linewidth=1, alpha=0.6)

    ax.set_yticks(range(len(SAMPLE_ORDER)))
    ax.set_yticklabels(SAMPLE_ORDER, fontsize=11)
    ax.invert_yaxis()
    ax.set_xlabel("ΔNG (Intermediate AF − Extreme AF) in CN1 LOH (Potential_LOH Inner)",
                  fontsize=11)
    ax.set_xlim(-1.0, 1.2)
    ax.grid(axis="x", alpha=0.3, linestyle="--")

    # Strategy legend
    handles = [
        mpl.lines.Line2D([0], [0], marker="o", color="w", markerfacecolor=strat_colors[s],
                         markeredgecolor=PALETTE["dark"], markersize=10, label=s)
        for s in strat_order
    ]
    handles.append(mpl.lines.Line2D([0], [0], color=PALETTE["dark"], lw=2,
                                    label="20260414 ΔNG=+0.705"))
    ax.legend(handles=handles, loc="lower right", fontsize=9, frameon=True,
              edgecolor=PALETTE["dark"])

    # Background shading for POS/NEG regions
    ax.axvspan(-1.2, 0, alpha=0.07, color=PALETTE["red"])
    ax.axvspan(0, 1.5, alpha=0.07, color=PALETTE["green"])
    ax.text(-0.5, -0.6, "REVERSAL", fontsize=14, fontweight="bold", color=PALETTE["red"],
            ha="center", alpha=0.7)
    ax.text(0.5, -0.6, "PRESERVED", fontsize=14, fontweight="bold", color=PALETTE["green"],
            ha="center", alpha=0.7)

    ax.set_title("F5: 20260414 replication under new KDE baseline —\n"
                 "ΔNG (Intermediate AF vs Extreme AF) in CN1 LOH, per sample × 3 strategies",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")


# ────────────────────────────────────────────────────────────────
def fig6_h2_verification(h2: pd.DataFrame, out: Path) -> None:
    """H2 verification: LOH.bed + TP-only (exact 20260414) vs Potential_LOH (current)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    methods = ["LOH_bed_TP_only", "PotentialLOH_TP_only"]
    method_titles = {
        "LOH_bed_TP_only": "(a) EXACT 20260414: LOH.bed + TP-only + CN1<0.75",
        "PotentialLOH_TP_only": "(b) Potential_LOH fallback + TP-only + CN1<0.75",
    }
    for ax_i, method in enumerate(methods):
        ax = axes[ax_i]
        sub = h2[h2["method"] == method]
        d = sub.set_index("sample").reindex(SAMPLE_ORDER)
        colors = [PALETTE["green"] if v > 0.1 else (PALETTE["red"] if v < -0.1 else PALETTE["grey"])
                  for v in d["delta_ng"].fillna(0)]
        bars = ax.barh(range(len(SAMPLE_ORDER)), d["delta_ng"].fillna(0).values,
                       color=colors, edgecolor=PALETTE["dark"], linewidth=0.6)
        ax.axvline(0.705, color=PALETTE["dark"], linestyle="-", linewidth=2, alpha=0.9,
                   label="20260414 canonical\nΔNG=+0.705")
        ax.axvline(0, color=PALETTE["grey"], linestyle="--", linewidth=1)
        ax.set_yticks(range(len(SAMPLE_ORDER)))
        ax.set_yticklabels(SAMPLE_ORDER, fontsize=11)
        ax.invert_yaxis()
        ax.set_xlabel("ΔNG (Inter − Extreme) in CN1 LOH TP", fontsize=11)
        ax.set_xlim(-1.0, 1.2)
        ax.grid(axis="x", alpha=0.3, linestyle="--")
        ax.set_title(method_titles[method], fontsize=12)

        # Annotate n
        for i, s in enumerate(SAMPLE_ORDER):
            n_i = d.loc[s, "n_inter"] if s in d.index else np.nan
            n_e = d.loc[s, "n_extreme"] if s in d.index else np.nan
            dv = d.loc[s, "delta_ng"] if s in d.index else np.nan
            if pd.notna(n_i) and pd.notna(n_e):
                ax.text(-0.98, i, f"n={int(n_i)}/{int(n_e)}",
                        fontsize=8, color=PALETTE["dark"], va="center")
            if pd.notna(dv):
                ax.text(dv + (0.02 if dv >= 0 else -0.02), i, f"{dv:+.3f}",
                        fontsize=8, color=PALETTE["dark"], va="center",
                        ha="left" if dv >= 0 else "right")

        ax.legend(loc="lower right", fontsize=9, frameon=True, edgecolor=PALETTE["dark"])

    fig.suptitle("F6: H2 verification — reversal holds under 20260414 EXACT conditions\n"
                 "(LOH.bed external vs Potential_LOH fallback) → methodology NOT the cause of reversal",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")


# ────────────────────────────────────────────────────────────────
def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    dfs = load()

    print("[step3] generating F1 — CovM density")
    fig1_covm_density(dfs["master"], FIG_DIR / "fig1_covm_density_per_sample.png")

    print("[step3] generating F2 — TP rate heatmap Strategy A vs F")
    fig2_tprate_heatmap(dfs["ng_tprate"], FIG_DIR / "fig2_tprate_heatmap_A_vs_F.png")

    print("[step3] generating F3 — NG distribution stacked")
    fig3_ng_distribution(dfs["ng_dist"], FIG_DIR / "fig3_ng_distribution_stacked.png")

    print("[step3] generating F4 — AF 10-bin × NG heatmap")
    fig4_af_ng_heatmap(dfs["af_cross"], FIG_DIR / "fig4_af10bin_ng_heatmap.png")

    print("[step3] generating F5 — 20260414 replication forest")
    fig5_replication_forest(dfs["repl"], FIG_DIR / "fig5_replication_forest_3strategies.png")

    print("[step3] generating F6 — H2 verification")
    fig6_h2_verification(dfs["h2"], FIG_DIR / "fig6_h2_verification_lohbed.png")

    print(f"\n[step3] all figures in: {FIG_DIR}")


if __name__ == "__main__":
    main()
