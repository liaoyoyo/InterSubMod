#!/usr/bin/env python3
"""
Step 7: E1-E5 deep-validation figures (v2 report §11-14).

Figures:
  F7: Baseline vs scheme TP:FP ratio + fold-improvement (HCC1395 TO + COLO829 pf)
  F8: Per-sample x scheme TP_rate heatmap with n annotation (7 samples x 7 schemes)
  F9: S4 feature violin comparison (TP vs FP) for HCC1395 TO
  F10: Sub-scheme operating points (tp_recall vs tp_rate)
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

SCHEME_ORDER = [
    "S1_LOH_Strong_Extreme",
    "S2_Subclonal_LOH_Inter",
    "S3_Diploid_Het",
    "S4_NonLOH_Extreme",
    "S5_Combo",
    "S6_S1_NG>=3",
    "S7_Combo_NG>=3",
]
SCHEME_LABEL_SHORT = {
    "S1_LOH_Strong_Extreme": "S1\nLOH_Strong\n+Extreme",
    "S2_Subclonal_LOH_Inter": "S2\nSub/Weak\n+Inter",
    "S3_Diploid_Het": "S3\nDiploid\nHet",
    "S4_NonLOH_Extreme": "S4\nNonLOH\n+Extreme",
    "S5_Combo": "S5\nCombo\n(¬S4)",
    "S6_S1_NG>=3": "S6\nS1+NG>=3",
    "S7_Combo_NG>=3": "S7\nCombo+NG>=3",
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


def fig7_baseline_vs_scheme() -> None:
    df = pd.read_csv(DATA_DIR / "tpfp_per_sample_scheme_full.tsv", sep="\t")
    panels = [
        ("HCC1395", "to_pileup", "HCC1395 TO mode (baseline TP:FP = 2.47:1)"),
        ("COLO829", "paired_full", "COLO829 paired_full (baseline = 15.5:1)"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 9),
                             gridspec_kw={"hspace": 0.5, "wspace": 0.28})

    for row, (samp, mode, title) in enumerate(panels):
        sub = df[(df["sample"] == samp) & (df["mode"] == mode)].set_index("scheme")
        sub = sub.reindex(SCHEME_ORDER)

        # Left: TP:FP ratio
        ax = axes[row, 0]
        scheme_vals = sub["scheme_ratio"].replace([np.inf, -np.inf], np.nan).values
        baseline = sub["baseline_ratio"].iloc[0] if not sub["baseline_ratio"].isna().all() else np.nan
        scheme_vals_plot = np.where(np.isfinite(scheme_vals), scheme_vals, 0.0)
        colors = [
            PALETTE["red"] if (v is None or not np.isfinite(v) or v < baseline * 0.95)
            else PALETTE["green"] if v >= baseline * 1.5
            else PALETTE["blue"]
            for v in scheme_vals
        ]
        x = np.arange(len(SCHEME_ORDER))
        ax.bar(x, scheme_vals_plot, color=colors, alpha=0.85, edgecolor=PALETTE["dark"])
        ax.axhline(baseline, color=PALETTE["accent"], linestyle="--", linewidth=2,
                   label=f"baseline = {baseline:.1f}")
        for xi, (val, n) in enumerate(zip(scheme_vals, sub["n"].values)):
            if not np.isfinite(val):
                ax.text(xi, ax.get_ylim()[1] * 0.05, "n/a", ha="center",
                        fontsize=8, color=PALETTE["grey"])
            else:
                ax.text(xi, val * 1.02, f"{val:.1f}\n(n={int(n)})",
                        ha="center", va="bottom", fontsize=8, color=PALETTE["dark"])
        ax.set_xticks(x)
        ax.set_xticklabels([SCHEME_LABEL_SHORT[s] for s in SCHEME_ORDER], fontsize=8)
        ax.set_ylabel("TP:FP ratio", fontsize=11)
        ax.set_title(f"{title}\nabsolute TP:FP ratio", fontsize=11)
        ax.legend(loc="upper left", fontsize=9)

        # Right: fold-improvement
        ax = axes[row, 1]
        folds = sub["fold_improvement"].replace([np.inf, -np.inf], np.nan).values
        fold_plot = np.where(np.isfinite(folds), folds, 0.0)
        colors = [
            PALETTE["green"] if (np.isfinite(f) and f >= 3.0)
            else PALETTE["blue"] if (np.isfinite(f) and f >= 1.5)
            else PALETTE["grey"] if (np.isfinite(f) and f >= 0.95)
            else PALETTE["red"]
            for f in folds
        ]
        ax.bar(x, fold_plot, color=colors, alpha=0.85, edgecolor=PALETTE["dark"])
        ax.axhline(1.0, color=PALETTE["accent"], linestyle="--", linewidth=2,
                   label="baseline (1.0x)")
        for xi, f in enumerate(folds):
            if not np.isfinite(f):
                ax.text(xi, 0.05, "n/a", ha="center", fontsize=8, color=PALETTE["grey"])
            else:
                ax.text(xi, f * 1.02, f"{f:.2f}x", ha="center", va="bottom",
                        fontsize=8, color=PALETTE["dark"], fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([SCHEME_LABEL_SHORT[s] for s in SCHEME_ORDER], fontsize=8)
        ax.set_ylabel("fold-improvement (scheme_ratio / baseline_ratio)", fontsize=11)
        ax.set_title(f"{title}\nfold-improvement over baseline", fontsize=11)
        ax.legend(loc="upper right", fontsize=9)

    fig.suptitle("Figure v2-7 | Baseline vs Scheme TP:FP + Fold-Improvement",
                 fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(FIG_DIR / "fig_v2_7a_paired_full_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[step7] wrote fig_v2_7a_paired_full_comparison.png")


def fig8_per_sample_heatmap() -> None:
    df = pd.read_csv(DATA_DIR / "tpfp_per_sample_scheme_full.tsv", sep="\t")
    # Build sample_mode rows: HCC1395 TO first, then paired_full in SAMPLE_ORDER
    sample_modes = [("HCC1395", "to_pileup", "HCC1395_TO")]
    for s in SAMPLE_ORDER:
        sample_modes.append((s, "paired_full", f"{s}_pf"))

    n_rows = len(sample_modes); n_cols = len(SCHEME_ORDER)
    tp_mat = np.full((n_rows, n_cols), np.nan)
    n_mat = np.zeros((n_rows, n_cols), dtype=int)

    for i, (samp, mode, _) in enumerate(sample_modes):
        for j, scheme in enumerate(SCHEME_ORDER):
            cell = df[(df["sample"] == samp) & (df["mode"] == mode) & (df["scheme"] == scheme)]
            if not cell.empty and int(cell["n"].iloc[0]) >= 10:
                tp_mat[i, j] = float(cell["scheme_precision"].iloc[0])
                n_mat[i, j] = int(cell["n"].iloc[0])

    fig, ax = plt.subplots(figsize=(14, 7))
    im = ax.imshow(tp_mat, cmap="RdYlGn", vmin=0.3, vmax=1.0, aspect="auto")

    for i in range(n_rows):
        for j in range(n_cols):
            if np.isnan(tp_mat[i, j]):
                ax.text(j, i, "n<10" if n_mat[i, j] > 0 else "n=0",
                        ha="center", va="center", fontsize=8, color=PALETTE["grey"])
            else:
                val = tp_mat[i, j]
                txt = f"{val:.2f}\n(n={n_mat[i, j]})"
                color = "white" if val < 0.55 or val > 0.85 else PALETTE["dark"]
                ax.text(j, i, txt, ha="center", va="center", fontsize=9,
                        color=color, fontweight="bold" if val >= 0.9 else "normal")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([SCHEME_LABEL_SHORT[s] for s in SCHEME_ORDER], fontsize=9)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([lbl for _, _, lbl in sample_modes], fontsize=10)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("TP rate (precision)", fontsize=10)
    ax.set_title("Figure v2-8 | Per-sample × Scheme TP Rate\n"
                 "(cells with n≥10; top row HCC1395 TO is main testbed)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_v2_8_per_sample_scheme_heatmap.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[step7] wrote fig_v2_8_per_sample_scheme_heatmap.png")


def fig9_s4_feature_violin() -> None:
    """S4 TP vs FP feature distributions on HCC1395 TO."""
    master = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
    df = pd.read_csv(master, sep="\t", compression="gzip", low_memory=False)
    df = df[(df["sample"] == "HCC1395") & (df["mode"] == "to_pileup")].copy()
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("None")
    m_s4 = (df["LOH_Subtype"] == "None") & (df["AF_class"] == "Extreme")
    s4 = df[m_s4]
    tp = s4[s4["tp_label"] == 1]
    fp = s4[s4["tp_label"] == 0]
    print(f"[fig9] S4 TO: TP={len(tp)}, FP={len(fp)}")

    features = [
        ("NumReads", "Read count", "linear"),
        ("NumCpGs", "CpG count", "linear"),
        ("AlleleDelta", "Allele delta (~AF)", "linear"),
        ("CovM_used", "Coverage multiple", "linear"),
        ("HPFineNGroups", "HP fine N groups (1-4)", "linear"),
    ]

    e4 = pd.read_csv(DATA_DIR / "tpfp_feature_diffs.tsv", sep="\t")
    e4_s4 = e4[(e4["sample"] == "HCC1395") & (e4["mode"] == "to_pileup") &
               (e4["scheme"] == "S4_NonLOH_Extreme")].set_index("feature")

    fig, axes = plt.subplots(1, 5, figsize=(20, 5))
    for ax, (feat, label, scale) in zip(axes, features):
        tp_v = tp[feat].dropna().values
        fp_v = fp[feat].dropna().values
        parts = ax.violinplot([tp_v, fp_v], positions=[0, 1],
                               showmedians=True, widths=0.75)
        for pc, color in zip(parts["bodies"], [PALETTE["green"], PALETTE["red"]]):
            pc.set_facecolor(color); pc.set_alpha(0.7)
            pc.set_edgecolor(PALETTE["dark"])
        for key in ("cmaxes", "cmins", "cbars", "cmedians"):
            if key in parts:
                parts[key].set_edgecolor(PALETTE["dark"])
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"TP\n(n={len(tp_v)})", f"FP\n(n={len(fp_v)})"])
        ax.set_ylabel(label, fontsize=10)
        if scale == "log":
            ax.set_yscale("log")

        if feat in e4_s4.index:
            d = e4_s4.loc[feat, "cliffs_delta"]
            p = e4_s4.loc[feat, "mannwhitney_p"]
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            ax.set_title(f"{feat}\nCliff's δ={d:.2f} {sig}",
                         fontsize=11, fontweight="bold")
        else:
            ax.set_title(feat, fontsize=11)

    fig.suptitle("Figure v2-9 | S4 (NonLOH+Extreme) TP vs FP feature distributions (HCC1395 TO)\n"
                 "TP medians green, FP red; Cliff's δ < 0 means TP smaller than FP",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_v2_9_feature_violin_in_S4.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[step7] wrote fig_v2_9_feature_violin_in_S4.png")


def fig10_subscheme_operating() -> None:
    """Sub-scheme operating points for HCC1395 TO overlaid on S1-S7."""
    full = pd.read_csv(DATA_DIR / "tpfp_per_sample_scheme_full.tsv", sep="\t")
    subs = pd.read_csv(DATA_DIR / "tpfp_subschemes.tsv", sep="\t")

    base = full[(full["sample"] == "HCC1395") & (full["mode"] == "to_pileup")].copy()
    base["type"] = "S1-S7"
    subs_to = subs[(subs["sample"] == "HCC1395") & (subs["mode"] == "to_pileup")].copy()
    subs_to = subs_to.rename(columns={"scheme_id": "scheme"})
    subs_to["type"] = "sub-scheme"
    subs_to["scheme_precision"] = subs_to["tp_rate"]

    fig, ax = plt.subplots(figsize=(12, 7))

    # baseline horizontal line
    baseline_prec = float(base["baseline_precision"].iloc[0])
    ax.axhline(baseline_prec, color=PALETTE["accent"], linestyle="--", linewidth=2,
               label=f"baseline precision = {baseline_prec:.3f}")

    def plot_group(sub_df, marker, label, color_factor):
        xs = sub_df["tp_recall"].values
        ys = sub_df["scheme_precision"].values if "scheme_precision" in sub_df.columns else sub_df["tp_rate"].values
        ns = sub_df["n"].values
        sizes = np.clip(np.sqrt(ns) * 3, 30, 350)
        ax.scatter(xs, ys, s=sizes, marker=marker, alpha=0.75,
                   edgecolors=PALETTE["dark"], linewidths=1.2,
                   c=color_factor, label=label)
        for x, y, name, n in zip(xs, ys, sub_df["scheme"].values, ns):
            if n < 5:
                continue
            ax.annotate(name.split("_")[0], (x, y),
                        textcoords="offset points", xytext=(6, 6), fontsize=8,
                        color=PALETTE["dark"])

    plot_group(base[base["scheme"].isin(SCHEME_ORDER)], "o", "S1-S7 (primary)",
               PALETTE["blue"])
    plot_group(subs_to, "^", "sub-schemes (E5)", PALETTE["green"])

    ax.set_xlabel("TP recall (fraction of caller TPs retained)", fontsize=11)
    ax.set_ylabel("Precision (TP rate within scheme)", fontsize=11)
    ax.set_xlim(-0.02, 1.0)
    ax.set_ylim(0.4, 1.03)
    ax.set_title("Figure v2-10 | Sub-scheme operating points vs S1-S7 (HCC1395 TO)\n"
                 "marker size ∝ √n; sub-schemes (▲) mostly fail to beat S5 (●)",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="lower left", fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_v2_10_subscheme_operating_points.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[step7] wrote fig_v2_10_subscheme_operating_points.png")


def fig7b_to_cross_sample() -> None:
    """2026-04-22: 6-sample TO cross-sample baseline-vs-scheme + fold-improvement.

    After archive TO rerun, we have 6 TO samples (HCC1395, HCC1395_DORADO,
    HCC1937, HCC1954, H2009, H1437). Layout: 3 rows x 4 cols (2 samples per row,
    each sample gets TP:FP ratio + fold-improvement subplots).
    """
    df = pd.read_csv(DATA_DIR / "tpfp_per_sample_scheme_full.tsv", sep="\t")

    to_samples = [
        ("HCC1395", "HCC1395 TO"),
        ("HCC1395_DORADO", "HCC1395_DORADO TO"),
        ("HCC1937", "HCC1937 TO"),
        ("HCC1954", "HCC1954 TO"),
        ("H2009", "H2009 TO"),
        ("H1437", "H1437 TO"),
    ]

    n_samples = len(to_samples)
    ncols = 4  # 2 samples per row * 2 subplots per sample
    nrows = (n_samples + 1) // 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 4.5 * nrows),
                             gridspec_kw={"hspace": 0.55, "wspace": 0.32})

    x = np.arange(len(SCHEME_ORDER))

    for idx, (samp, label) in enumerate(to_samples):
        row = idx // 2
        col_pair = idx % 2
        ax_ratio = axes[row, col_pair * 2]
        ax_fold = axes[row, col_pair * 2 + 1]

        sub = df[(df["sample"] == samp) & (df["mode"] == "to_pileup")].set_index("scheme")
        sub = sub.reindex(SCHEME_ORDER)

        if sub["baseline_ratio"].isna().all():
            for ax, txt in [(ax_ratio, "NO DATA"), (ax_fold, "NO DATA")]:
                ax.text(0.5, 0.5, f"{label}\n{txt}", ha="center", va="center",
                        fontsize=14, color=PALETTE["grey"], transform=ax.transAxes)
                ax.set_xticks([]); ax.set_yticks([])
            continue

        baseline = float(sub["baseline_ratio"].iloc[0])
        baseline_tp = int(sub["baseline_tp"].iloc[0])
        baseline_fp = int(sub["baseline_fp"].iloc[0])

        # TP:FP ratio panel
        scheme_vals = sub["scheme_ratio"].replace([np.inf, -np.inf], np.nan).values
        scheme_vals_plot = np.where(np.isfinite(scheme_vals), scheme_vals, 0.0)
        colors = [
            PALETTE["red"] if (v is None or not np.isfinite(v) or v < baseline * 0.95)
            else PALETTE["green"] if v >= baseline * 1.5
            else PALETTE["blue"]
            for v in scheme_vals
        ]
        ax_ratio.bar(x, scheme_vals_plot, color=colors, alpha=0.85, edgecolor=PALETTE["dark"])
        ax_ratio.axhline(baseline, color=PALETTE["accent"], linestyle="--", linewidth=2,
                         label=f"baseline={baseline:.1f}")
        for xi, (val, n) in enumerate(zip(scheme_vals, sub["n"].values)):
            if np.isfinite(val):
                ax_ratio.text(xi, val * 1.02, f"{val:.1f}\nn={int(n)}",
                              ha="center", va="bottom", fontsize=7, color=PALETTE["dark"])
            elif n > 0:
                ax_ratio.text(xi, 0, f"n={int(n)}",
                              ha="center", va="bottom", fontsize=7, color=PALETTE["grey"])
        ax_ratio.set_xticks(x)
        ax_ratio.set_xticklabels([SCHEME_LABEL_SHORT[s] for s in SCHEME_ORDER], fontsize=7)
        ax_ratio.set_ylabel("TP:FP ratio", fontsize=10)
        ax_ratio.set_title(f"{label}\nTP:FP ratio (TP={baseline_tp}, FP={baseline_fp})",
                           fontsize=10)
        ax_ratio.legend(loc="upper left", fontsize=8)

        # Fold-improvement panel
        folds = sub["fold_improvement"].replace([np.inf, -np.inf], np.nan).values
        fold_plot = np.where(np.isfinite(folds), folds, 0.0)
        colors_f = [
            PALETTE["green"] if (np.isfinite(f) and f >= 3.0)
            else PALETTE["blue"] if (np.isfinite(f) and f >= 1.5)
            else PALETTE["grey"] if (np.isfinite(f) and f >= 0.95)
            else PALETTE["red"]
            for f in folds
        ]
        ax_fold.bar(x, fold_plot, color=colors_f, alpha=0.85, edgecolor=PALETTE["dark"])
        ax_fold.axhline(1.0, color=PALETTE["accent"], linestyle="--", linewidth=2,
                        label="baseline (1.0x)")
        for xi, f in enumerate(folds):
            if np.isfinite(f):
                ax_fold.text(xi, f * 1.02, f"{f:.2f}x",
                             ha="center", va="bottom", fontsize=7,
                             color=PALETTE["dark"], fontweight="bold")
        ax_fold.set_xticks(x)
        ax_fold.set_xticklabels([SCHEME_LABEL_SHORT[s] for s in SCHEME_ORDER], fontsize=7)
        ax_fold.set_ylabel("fold-improvement", fontsize=10)
        ax_fold.set_title(f"{label}\nfold-improvement", fontsize=10)
        ax_fold.legend(loc="upper right", fontsize=8)

    # Hide any unused subplots
    total_subplots = nrows * ncols
    used = n_samples * 2
    for k in range(used, total_subplots):
        r, c = divmod(k, ncols)
        axes[r, c].set_visible(False)

    fig.suptitle("Figure v2-7 | Cross-sample TO: Baseline vs Scheme TP:FP + Fold-Improvement\n"
                 "(6 samples post-archive ISM rerun 2026-04-22; S1-S7 schemes)",
                 fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    out_path = FIG_DIR / "fig_v2_7_baseline_tp_fp_fold.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[step7] wrote {out_path.name} (6 TO panels)")


def main() -> None:
    set_style()
    fig7_baseline_vs_scheme()
    fig7b_to_cross_sample()
    fig8_per_sample_heatmap()
    fig9_s4_feature_violin()
    fig10_subscheme_operating()
    print(f"[step7] all figures -> {FIG_DIR}")


if __name__ == "__main__":
    main()
