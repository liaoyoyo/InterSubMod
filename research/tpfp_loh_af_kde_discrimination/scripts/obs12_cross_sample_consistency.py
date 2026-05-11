#!/usr/bin/env python3
"""
obs12 | Cross-sample consistency of 5D cells (2 figures + TSV).

Aggregate 5D-cell TP rate across all 8 sample-modes.  For each cell produce:
  - ``n_{sample}`` count in each sample
  - ``tp_rate_{sample}`` TP rate if n>=20 else NaN
  - ``n_samples_n20`` number of samples with n>=20
  - ``n_samples_high`` number of samples where tp_rate > baseline(sample) + 0.10
  - ``cross_sample_avg_tp`` mean of available tp_rate across samples with n>=20

Figures:
  - ``consistency/cell_consistency_5d.png`` — Top-30 cells sorted by
    ``n_samples_high`` descending; stacked bar of per-sample TP rate
  - ``consistency/cell_consistency_3d_locn.png`` — collapse to LOH x AF x CN
    (fig_v2_1 coordinate system), scatter sized by avg_n, coloured by
    ``n_samples_high`` (0-8)

Writes TSV: ``data/cell_consistency_matrix.tsv``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (
    DATA_DIR, DIM_ORDERS, DIM_SHORT_LABEL, FIG_NEW_DIR, PALETTE,
    SAMPLE_MODE_ORDER, SAMPLE_MODE_TO, apply_style, build_5d_key,
    compute_tp_rate_by_key, load_master, research_suptitle,
    sample_label, takeaway_caption,
)

CONS_DIR = FIG_NEW_DIR / "consistency"
KEY_DIMS = ["LOH_Subtype", "AF_class", "cn_tier_F", "HPFineNGroups", "nr_band"]
# "high" means: cell FP rate <= 0.5 * baseline FP rate
# equivalent to: (1 - tp_rate) <= 0.5 * (1 - baseline)
# rationale: +0.10 absolute is impossible when baseline > 0.90
HIGH_FP_HALVING = 0.5


def build_consistency_matrix(df5: pd.DataFrame, sm_list: list | None = None) -> tuple[pd.DataFrame, dict]:
    """Return (long-format tidy TP-rate df per cell, dict of per-sample baseline).

    ``sm_list`` allows restricting the computation to a subset of sample-modes
    (e.g., TO-only).  Default = all 13 SAMPLE_MODE_ORDER.
    """
    sm_list = sm_list if sm_list is not None else SAMPLE_MODE_ORDER
    baselines = {}
    rows = []
    for sample, mode in sm_list:
        sub = df5[(df5["sample"] == sample) & (df5["mode"] == mode)]
        if len(sub) == 0:
            continue
        baselines[(sample, mode)] = sub["tp_label"].mean()
        agg = compute_tp_rate_by_key(sub, KEY_DIMS, min_n=20)
        agg["sample_label"] = sample_label(sample, mode)
        rows.append(agg)
    long_df = pd.concat(rows, ignore_index=True)
    return long_df, baselines


def pivot_matrix(long_df: pd.DataFrame, baselines: dict) -> pd.DataFrame:
    """Convert long-format per-cell/per-sample counts into wide matrix + aggregates."""
    # pivot n
    n_wide = long_df.pivot_table(
        index=KEY_DIMS, columns="sample_label", values="n",
        fill_value=0, aggfunc="sum",
    )
    n_wide.columns = [f"n_{c}" for c in n_wide.columns]
    # pivot tp_rate (NaN for n<20)
    tp_wide = long_df.pivot_table(
        index=KEY_DIMS, columns="sample_label", values="tp_rate",
        aggfunc="first",
    )
    tp_wide.columns = [f"tp_{c}" for c in tp_wide.columns]
    merged = n_wide.join(tp_wide, how="outer").reset_index()

    # derived fields
    n_cols = [c for c in merged.columns if c.startswith("n_")]
    tp_cols = [c for c in merged.columns if c.startswith("tp_")]
    merged["n_samples_n20"] = (merged[n_cols] >= 20).sum(axis=1)
    merged["cross_sample_avg_tp"] = merged[tp_cols].mean(axis=1, skipna=True)
    merged["cross_sample_min_tp"] = merged[tp_cols].min(axis=1, skipna=True)

    # high-threshold counter: cell FP rate <= 0.5 x baseline FP rate
    def count_high(row):
        cnt = 0
        for (sample, mode), b in baselines.items():
            lbl = sample_label(sample, mode)
            tp_val = row.get(f"tp_{lbl}")
            if pd.isna(tp_val):
                continue
            baseline_fp = max(1e-6, 1.0 - b)
            cell_fp = max(0.0, 1.0 - tp_val)
            if cell_fp <= HIGH_FP_HALVING * baseline_fp:
                cnt += 1
        return cnt
    merged["n_samples_high"] = merged.apply(count_high, axis=1)

    # global weighted n (for scatter size)
    merged["total_n"] = merged[n_cols].sum(axis=1)
    return merged


def plot_top_cells_5d(matrix: pd.DataFrame, baselines: dict,
                      sm_list: list | None = None,
                      out_name: str = "cell_consistency_5d.png",
                      title_suffix: str = "") -> tuple[Path, int, int]:
    """Top-30 cells by n_samples_high then cross_sample_avg_tp (desc), stacked bar.

    Returns (output_path, count_high_ge_5, count_high_ge_3).
    """
    sm_list = sm_list if sm_list is not None else SAMPLE_MODE_ORDER
    eligible = matrix[matrix["n_samples_n20"] >= 2].copy()
    eligible = eligible.sort_values(
        ["n_samples_high", "cross_sample_avg_tp", "total_n"],
        ascending=[False, False, False],
    ).head(30).reset_index(drop=True)

    fig = plt.figure(figsize=(18, 11))
    gs = gridspec.GridSpec(1, 2, width_ratios=[6, 1], wspace=0.05, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    ax_info = fig.add_subplot(gs[0, 1])
    ax_info.axis("off")

    norm = Normalize(vmin=0.3, vmax=1.0)
    cmap = plt.get_cmap("RdYlGn")
    n_samples = len(sm_list)
    bar_width = 0.9 / n_samples
    ys = np.arange(len(eligible))

    for j, (s, m) in enumerate(sm_list):
        lbl = sample_label(s, m)
        col = f"tp_{lbl}"
        n_col = f"n_{lbl}"
        if col not in eligible.columns:
            continue
        values = eligible[col].values.astype(float)
        colors = [cmap(norm(v)) if not np.isnan(v) else PALETTE["grey"] for v in values]
        ax.barh(ys + (j - n_samples / 2) * bar_width, np.nan_to_num(values, nan=0.01),
                height=bar_width, color=colors,
                edgecolor=PALETTE["dark"], linewidth=0.2, alpha=0.95)

    def fmt_cell(row):
        return (
            f"{row['LOH_Subtype']}|{row['AF_class']}|{row['cn_tier_F']}|"
            f"NG={row['HPFineNGroups']}|NR={row['nr_band']}"
        )
    ax.set_yticks(ys)
    ax.set_yticklabels([fmt_cell(r) for _, r in eligible.iterrows()], fontsize=7)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("TP rate per sample (RdYlGn)", fontsize=10)
    ax.invert_yaxis()
    ax.grid(True, axis="x", alpha=0.25)
    ax.axvline(0.5, color=PALETTE["grey"], linestyle=":", linewidth=0.8)

    legend_txt = " | ".join([sample_label(s, m) for s, m in sm_list])
    ax.set_title(
        f"Bars within each row (top→bottom) = {legend_txt}",
        fontsize=8, color=PALETTE["grey"], loc="left",
    )

    lines = [f"cell aggregate info{title_suffix}:"]
    for i, row in eligible.iterrows():
        lines.append(
            f"#{i+1}: high={int(row['n_samples_high'])}/{n_samples}  "
            f"avg_tp={row['cross_sample_avg_tp']:.3f}  "
            f"n20={int(row['n_samples_n20'])}  N={int(row['total_n'])}"
        )
    ax_info.text(0.0, 0.995, "\n".join(lines[:31]), transform=ax_info.transAxes,
                 va="top", ha="left", fontsize=6.5, family="monospace",
                 color=PALETTE["dark"])

    # Research-question suptitle + takeaway
    high_ge_5 = int((matrix["n_samples_high"] >= 5).sum())
    high_ge_3 = int((matrix["n_samples_high"] >= 3).sum())
    if sm_list is SAMPLE_MODE_ORDER or (sm_list and len(sm_list) == 13):
        question = "哪些 5D cells 在 13 個 sample-mode 中跨樣本 FP 半減（TP 高）一致？"
        ctx_note = f"judgment: cell FP ≤ {HIGH_FP_HALVING:.1f}× baseline FP per sample"
    else:
        question = f"6 個 TO 樣本之間的 5D cells 共識：哪些 cells FP 半減跨 TO 樣本一致？"
        ctx_note = f"TO-only consensus · judgment: cell FP ≤ {HIGH_FP_HALVING:.1f}× baseline FP"
    research_suptitle(fig, question, context=ctx_note, y=0.975, fontsize=17)

    # takeaway with Hard Gate status
    top_cell = eligible.iloc[0] if len(eligible) else None
    if top_cell is None:
        takeaway_text = "○ Takeaway：無可用 cell。"
        colour = "grey"
    else:
        top_label = fmt_cell(top_cell)
        top_hits = int(top_cell["n_samples_high"])
        if n_samples == 13:
            hg1 = "✓ 觸發" if high_ge_5 >= 10 else "✗ 未觸發"
            hg2 = "✓ 觸發" if high_ge_3 >= 50 else "✗ 未觸發"
            takeaway_text = (f"Hard Gate | n_samples_high ≥ 5: {high_ge_5} cells ({hg1}, 閾值 10) · "
                             f"≥ 3: {high_ge_3} cells ({hg2}, 閾值 50) · "
                             f"Top-1: ({top_label}) 達 {top_hits}/{n_samples}")
            colour = "green" if high_ge_5 >= 10 else "blue"
        else:
            # TO-only (n_samples=6)
            # Scaled thresholds: ≥ 5 cells with high >= 4/6 → strong; ≥ 20 cells with high >= 3/6 → moderate
            strong_cnt = int((matrix["n_samples_high"] >= 4).sum())
            mod_cnt = int((matrix["n_samples_high"] >= 3).sum())
            takeaway_text = (f"TO-only | cells n_samples_high ≥ 4/6: {strong_cnt} · ≥ 3/6: {mod_cnt} · "
                             f"Top-1: ({top_label}) 達 {top_hits}/{n_samples}")
            colour = "green" if strong_cnt >= 5 else "blue"

    takeaway_caption(fig, takeaway_text, color=colour, y=0.015, fontsize=12)

    fig.tight_layout(rect=[0, 0.055, 1, 0.93])
    out = CONS_DIR / out_name
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out, high_ge_5, high_ge_3


def plot_3d_locn(matrix: pd.DataFrame) -> Path:
    """Collapse to LOH x AF x CN.  Scatter coloured by n_samples_high."""
    loh_order = DIM_ORDERS["LOH_Subtype"]
    af_order = DIM_ORDERS["AF_class"]
    cn_order = DIM_ORDERS["cn_tier_F"]

    # aggregate across NG/NR for each (LOH, AF, CN)
    grouped = (
        matrix.groupby(["LOH_Subtype", "AF_class", "cn_tier_F"], observed=True)
        .agg(
            n_samples_high=("n_samples_high", "max"),  # cell best-case
            n_samples_high_mean=("n_samples_high", "mean"),
            cross_avg_tp=("cross_sample_avg_tp", "mean"),
            total_n=("total_n", "sum"),
        )
        .reset_index()
    )

    # Each CN facet: ~3.8" wide x 6.5" tall so 5 LOH rows get proper vertical space.
    fig, axes = plt.subplots(1, len(cn_order), figsize=(4.0 * len(cn_order), 7.5),
                             sharey=True)
    n_samples_total = len(SAMPLE_MODE_ORDER)
    research_suptitle(
        fig,
        "LOH × AF × CN 空間中哪些格位跨 13 sample-mode 一致 FP 半減？",
        context=f"per CN facet: marker size = total n · color = max n_samples_high (0-{n_samples_total})",
        y=0.98, fontsize=17,
    )

    norm = Normalize(vmin=0, vmax=n_samples_total)
    cmap = plt.get_cmap("viridis")
    last_scatter = None

    for ax, cn in zip(axes, cn_order):
        sub = grouped[grouped["cn_tier_F"] == cn]
        if len(sub) == 0:
            ax.set_title(f"CN {cn} (empty)", fontsize=9)
            ax.set_xticks(range(len(af_order)))
            ax.set_xticklabels(af_order, rotation=30, ha="right", fontsize=7)
            ax.set_yticks(range(len(loh_order)))
            ax.set_yticklabels(loh_order, fontsize=7)
            continue
        for _, row in sub.iterrows():
            x = af_order.index(row["AF_class"]) if row["AF_class"] in af_order else None
            y = loh_order.index(row["LOH_Subtype"]) if row["LOH_Subtype"] in loh_order else None
            if x is None or y is None:
                continue
            size = max(30, min(1200, row["total_n"] / 4))
            sc = ax.scatter(x, y, s=size, c=[row["n_samples_high"]],
                            cmap=cmap, norm=norm,
                            edgecolor=PALETTE["dark"], linewidth=0.6, alpha=0.9)
            last_scatter = sc
            ax.text(x, y, f"{int(row['n_samples_high'])}",
                    ha="center", va="center", fontsize=7, color="white", fontweight="bold")
        ax.set_xticks(range(len(af_order)))
        ax.set_xticklabels(af_order, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(len(loh_order)))
        ax.set_yticklabels(loh_order, fontsize=8)
        ax.set_title(f"CN {cn}", fontsize=10, fontweight="bold", color=PALETTE["dark"])
        ax.grid(True, alpha=0.2)
        ax.set_xlim(-0.6, len(af_order) - 0.4)
        ax.set_ylim(-0.6, len(loh_order) - 0.4)

    if last_scatter is not None:
        cbar_ax = fig.add_axes([0.94, 0.18, 0.012, 0.68])
        cbar = fig.colorbar(last_scatter, cax=cbar_ax)
        cbar.set_label(f"n_samples_high (0-{n_samples_total})", fontsize=9)

    # Takeaway: best cell across 5 CN facets
    if grouped["n_samples_high"].max() > 0:
        best = grouped.loc[grouped["n_samples_high"].idxmax()]
        takeaway_text = (f"○ Takeaway：fig_v2_1 座標最一致格 = "
                         f"(LOH={best['LOH_Subtype']}, AF={best['AF_class']}, CN={best['cn_tier_F']}) "
                         f"達 {int(best['n_samples_high'])}/{n_samples_total} 樣本 FP 半減。")
        colour = "blue"
    else:
        takeaway_text = "○ Takeaway：無可用 cell 觸發 FP 半減。"
        colour = "grey"
    takeaway_caption(fig, takeaway_text, color=colour, y=0.015, fontsize=12)

    fig.tight_layout(rect=[0, 0.07, 0.92, 0.92])
    out = CONS_DIR / "cell_consistency_3d_locn.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    apply_style()
    CONS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    df = load_master()
    df5 = build_5d_key(df)
    df5 = df5[df5["HPFineNGroups"].isin(DIM_ORDERS["HPFineNGroups"])]

    # === All 13 sample-modes (main consistency) ===
    long_df, baselines = build_consistency_matrix(df5)
    matrix = pivot_matrix(long_df, baselines)

    tsv_path = DATA_DIR / "cell_consistency_matrix.tsv"
    matrix.sort_values(
        ["n_samples_high", "cross_sample_avg_tp", "total_n"],
        ascending=[False, False, False],
    ).to_csv(tsv_path, sep="\t", index=False)
    print(f"[obs12] wrote {tsv_path}  ({len(matrix)} cells)")

    high5 = (matrix["n_samples_high"] >= 5).sum()
    high3 = (matrix["n_samples_high"] >= 3).sum()
    print(f"[obs12] Hard Gate cells n_samples_high>=5: {high5}; >=3: {high3}")

    out1, _, _ = plot_top_cells_5d(matrix, baselines, sm_list=SAMPLE_MODE_ORDER,
                                   out_name="cell_consistency_5d.png")
    print(f"[obs12] wrote {out1}")
    out2 = plot_3d_locn(matrix)
    print(f"[obs12] wrote {out2}")

    # === TO-only consistency (6 samples; new 2026-04-22 V2) ===
    long_to, baselines_to = build_consistency_matrix(df5, sm_list=SAMPLE_MODE_TO)
    matrix_to = pivot_matrix(long_to, baselines_to)

    tsv_to = DATA_DIR / "cell_consistency_matrix_TO.tsv"
    matrix_to.sort_values(
        ["n_samples_high", "cross_sample_avg_tp", "total_n"],
        ascending=[False, False, False],
    ).to_csv(tsv_to, sep="\t", index=False)
    high5_to = (matrix_to["n_samples_high"] >= 5).sum()
    high4_to = (matrix_to["n_samples_high"] >= 4).sum()
    high3_to = (matrix_to["n_samples_high"] >= 3).sum()
    print(f"[obs12] TO-only: {len(matrix_to)} cells · n_samples_high ≥5:{high5_to} ≥4:{high4_to} ≥3:{high3_to}")

    out3, _, _ = plot_top_cells_5d(matrix_to, baselines_to, sm_list=SAMPLE_MODE_TO,
                                   out_name="cell_consistency_TO_only.png",
                                   title_suffix=" (TO-only, 6 samples)")
    print(f"[obs12] wrote {out3}")


if __name__ == "__main__":
    main()
