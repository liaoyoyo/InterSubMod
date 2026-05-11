#!/usr/bin/env python3
"""obs14 | Cross-sample TO observation using extended master (archive-enriched).

Dimensions available from archive: HPFineNGroups x AF_class x NumReads.
LOH_Subtype / CN tier are NaN in archive rows (pre-KDE pipeline).

Figures:
  - ``archive_to/fp_scale_comparison.png``  — FP counts per-sample + ratio TO/paired
  - ``archive_to/L1_af_class_to.png``       — AF_class x sample (TO only, 6 samples)
  - ``archive_to/L1_hpfinengroups_to.png``  — HPFineNGroups x sample (TO only)
  - ``archive_to/L2_af_x_ng_to.png``        — AF_class x HPFineNGroups heatmap per sample
  - ``archive_to/numreads_density_to.png``  — NumReads distribution (tp vs fp) per sample
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (
    AF_ORDER, DATA_DIR, FIG_NEW_DIR, PALETTE, apply_style,
)

EXT_MASTER = DATA_DIR / "master_extended.tsv.gz"
OUT_DIR = FIG_NEW_DIR / "archive_to"

TO_SAMPLES = ["HCC1395", "HCC1395_DORADO", "HCC1937",
              "HCC1954", "H2009", "H1437"]


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def load_extended() -> pd.DataFrame:
    df = pd.read_csv(EXT_MASTER, sep="\t", compression="gzip", low_memory=False)
    return df


def fig_fp_scale(df: pd.DataFrame) -> Path:
    summary = (df.groupby(["sample", "mode"])["tp_label"]
               .agg(["count", lambda x: (x == 0).sum(), lambda x: (x == 1).sum()])
               .rename(columns={"<lambda_0>": "n_fp", "<lambda_1>": "n_tp"})
               .reset_index())
    summary["fp_ratio"] = summary["n_fp"] / summary["count"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    colors = {"paired_full": PALETTE["blue"], "to_pileup": PALETTE["red"]}
    x = np.arange(len(TO_SAMPLES))
    width = 0.38
    for i, mode in enumerate(["paired_full", "to_pileup"]):
        fps = [summary[(summary["sample"] == s) & (summary["mode"] == mode)]["n_fp"].sum()
               for s in TO_SAMPLES]
        offset = (i - 0.5) * width
        bars = ax1.bar(x + offset, fps, width, label=mode, color=colors[mode])
        for b, v in zip(bars, fps):
            ax1.text(b.get_x() + b.get_width() / 2, v, f"{int(v):,}",
                     ha="center", va="bottom", fontsize=8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(TO_SAMPLES, rotation=20, ha="right", fontsize=9)
    ax1.set_ylabel("FP count", fontsize=10)
    ax1.set_title("FP count per sample: paired vs TO", fontsize=11, color=PALETTE["dark"])
    ax1.set_yscale("log")
    ax1.legend(fontsize=9)
    ax1.grid(axis="y", alpha=0.3)

    # Ratio panel
    ratios = []
    for s in TO_SAMPLES:
        fp_p = summary[(summary["sample"] == s) & (summary["mode"] == "paired_full")]["n_fp"].sum()
        fp_t = summary[(summary["sample"] == s) & (summary["mode"] == "to_pileup")]["n_fp"].sum()
        ratios.append(fp_t / max(fp_p, 1))
    bars = ax2.bar(x, ratios, color=PALETTE["accent"])
    for b, v in zip(bars, ratios):
        ax2.text(b.get_x() + b.get_width() / 2, v, f"{v:.0f}x",
                 ha="center", va="bottom", fontsize=9, fontweight="bold",
                 color=PALETTE["dark"])
    ax2.set_xticks(x)
    ax2.set_xticklabels(TO_SAMPLES, rotation=20, ha="right", fontsize=9)
    ax2.set_ylabel("TO/paired FP ratio", fontsize=10)
    ax2.set_title("TO-mode FP inflation ratio (TO FPs / paired FPs)",
                  fontsize=11, color=PALETTE["dark"])
    ax2.axhline(1, color=PALETTE["dark"], linestyle="--", alpha=0.5)
    ax2.grid(axis="y", alpha=0.3)

    fig.suptitle("TO-mode FP inflation across 6 samples (archive data, 2026-03-30)",
                 fontsize=13, fontweight="bold", color=PALETTE["dark"], y=1.02)
    fig.tight_layout()
    out = OUT_DIR / "fp_scale_comparison.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_l1_af_to(df: pd.DataFrame) -> Path:
    to = df[df["mode"] == "to_pileup"].copy()
    to = to[to["sample"].isin(TO_SAMPLES)]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=True)
    axes = axes.flatten()
    for idx, sample in enumerate(TO_SAMPLES):
        ax = axes[idx]
        sub = to[to["sample"] == sample]
        baseline = (sub["tp_label"] == 1).mean()
        rows = []
        for af_cls in AF_ORDER:
            s = sub[sub["AF_class"] == af_cls]
            n = len(s)
            n_tp = (s["tp_label"] == 1).sum()
            tp_rate = n_tp / n if n > 0 else float("nan")
            lo, hi = wilson_ci(n_tp, n)
            rows.append((af_cls, n, tp_rate, lo, hi))
        xs = np.arange(len(rows))
        ys = [r[2] for r in rows]
        lows = [r[2] - r[3] if not np.isnan(r[2]) else 0 for r in rows]
        highs = [r[4] - r[2] if not np.isnan(r[2]) else 0 for r in rows]
        bars = ax.bar(xs, ys, color=PALETTE["green"], alpha=0.8,
                      yerr=[lows, highs], capsize=4, error_kw={"alpha": 0.5})
        for b, (lab, n, r, lo, hi) in zip(bars, rows):
            if n > 0:
                ax.text(b.get_x() + b.get_width() / 2, min(r + 0.03, 0.95),
                        f"{r*100:.0f}%\n(n={n:,})", ha="center", va="bottom",
                        fontsize=7)
        ax.axhline(baseline, color=PALETTE["red"], linestyle="--",
                   linewidth=1, alpha=0.7, label=f"baseline {baseline*100:.0f}%")
        ax.set_xticks(xs)
        ax.set_xticklabels(AF_ORDER, fontsize=8, rotation=15)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"{sample} TO (n_FP={int((sub['tp_label']==0).sum()):,})",
                     fontsize=10, color=PALETTE["dark"])
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(axis="y", alpha=0.3)
        if idx % 3 == 0:
            ax.set_ylabel("TP rate", fontsize=9)

    fig.suptitle("L1 AF_class x TP-rate per TO sample (archive data)",
                 fontsize=13, fontweight="bold", color=PALETTE["dark"], y=1.00)
    fig.tight_layout()
    out = OUT_DIR / "L1_af_class_to.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_l1_ng_to(df: pd.DataFrame) -> Path:
    to = df[df["mode"] == "to_pileup"].copy()
    to = to[to["sample"].isin(TO_SAMPLES)]
    ng_values = sorted(to["HPFineNGroups"].dropna().unique())[:6]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=True)
    axes = axes.flatten()
    for idx, sample in enumerate(TO_SAMPLES):
        ax = axes[idx]
        sub = to[to["sample"] == sample]
        baseline = (sub["tp_label"] == 1).mean()
        rows = []
        for ng in ng_values:
            s = sub[sub["HPFineNGroups"] == ng]
            n = len(s)
            n_tp = (s["tp_label"] == 1).sum()
            tp_rate = n_tp / n if n > 0 else float("nan")
            lo, hi = wilson_ci(n_tp, n)
            rows.append((ng, n, tp_rate, lo, hi))
        xs = np.arange(len(rows))
        ys = [r[2] for r in rows]
        lows = [r[2] - r[3] if not np.isnan(r[2]) else 0 for r in rows]
        highs = [r[4] - r[2] if not np.isnan(r[2]) else 0 for r in rows]
        bars = ax.bar(xs, ys, color=PALETTE["blue"], alpha=0.8,
                      yerr=[lows, highs], capsize=4, error_kw={"alpha": 0.5})
        for b, (lab, n, r, lo, hi) in zip(bars, rows):
            if n > 0 and not np.isnan(r):
                ax.text(b.get_x() + b.get_width() / 2, min(r + 0.03, 0.95),
                        f"{r*100:.0f}%\n(n={n:,})", ha="center", va="bottom",
                        fontsize=7)
        ax.axhline(baseline, color=PALETTE["red"], linestyle="--",
                   linewidth=1, alpha=0.7, label=f"baseline {baseline*100:.0f}%")
        ax.set_xticks(xs)
        ax.set_xticklabels([f"NG={int(v)}" for v in ng_values], fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"{sample} TO (n_FP={int((sub['tp_label']==0).sum()):,})",
                     fontsize=10, color=PALETTE["dark"])
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(axis="y", alpha=0.3)
        if idx % 3 == 0:
            ax.set_ylabel("TP rate", fontsize=9)

    fig.suptitle("L1 HPFineNGroups x TP-rate per TO sample (archive data)",
                 fontsize=13, fontweight="bold", color=PALETTE["dark"], y=1.00)
    fig.tight_layout()
    out = OUT_DIR / "L1_hpfinengroups_to.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_l2_af_ng_to(df: pd.DataFrame) -> Path:
    to = df[df["mode"] == "to_pileup"].copy()
    to = to[to["sample"].isin(TO_SAMPLES)]
    ng_values = [int(v) for v in sorted(to["HPFineNGroups"].dropna().unique())[:5]]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
    axes = axes.flatten()
    norm = Normalize(vmin=0.3, vmax=1.0)
    cmap = plt.get_cmap("RdYlGn")
    cmap.set_bad(color=PALETTE["grey"], alpha=0.3)
    last_im = None

    for idx, sample in enumerate(TO_SAMPLES):
        ax = axes[idx]
        sub = to[to["sample"] == sample]
        mat_rate = np.full((len(AF_ORDER), len(ng_values)), np.nan)
        mat_n = np.zeros_like(mat_rate, dtype=int)
        for i, af_cls in enumerate(AF_ORDER):
            for j, ng in enumerate(ng_values):
                s = sub[(sub["AF_class"] == af_cls) & (sub["HPFineNGroups"] == ng)]
                n = len(s)
                if n < 20:
                    continue
                mat_rate[i, j] = (s["tp_label"] == 1).mean()
                mat_n[i, j] = n
        last_im = ax.imshow(np.ma.masked_invalid(mat_rate), cmap=cmap,
                            norm=norm, aspect="auto")
        ax.set_xticks(range(len(ng_values)))
        ax.set_yticks(range(len(AF_ORDER)))
        ax.set_xticklabels([f"NG={v}" for v in ng_values], fontsize=8)
        ax.set_yticklabels(AF_ORDER, fontsize=8)
        ax.set_title(f"{sample} TO", fontsize=10, color=PALETTE["dark"])
        for yi in range(len(AF_ORDER)):
            for xi in range(len(ng_values)):
                r = mat_rate[yi, xi]
                n = mat_n[yi, xi]
                if n == 0:
                    continue
                col = "black" if (not np.isnan(r)) and r > 0.55 else "white"
                ax.text(xi, yi, f"{r*100:.0f}\n{n}",
                        ha="center", va="center", fontsize=6.5, color=col)

    if last_im is not None:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.7])
        fig.colorbar(last_im, cax=cbar_ax, label="TP rate")

    fig.suptitle("L2 AF_class x HPFineNGroups heatmap per TO sample (archive)",
                 fontsize=13, fontweight="bold", color=PALETTE["dark"], y=0.98)
    fig.tight_layout(rect=[0, 0, 0.9, 0.96])
    out = OUT_DIR / "L2_af_x_ng_to.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_extended()
    print(f"[obs14] extended master: {len(df):,} rows")

    figs = [
        fig_fp_scale(df),
        fig_l1_af_to(df),
        fig_l1_ng_to(df),
        fig_l2_af_ng_to(df),
    ]
    for f in figs:
        print(f"[obs14] wrote {f}")


if __name__ == "__main__":
    main()
