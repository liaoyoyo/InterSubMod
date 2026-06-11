#!/usr/bin/env python3
"""
C1: Per-chromosome HCC1395 baseline vs V6 ISM full-genome sweep.

Purpose:
  Fill the "step-by-step full-genome observation" gap in PI 4-goal report.
  Show which chromosomes V6 improves, which regress (counter-example
  candidates), and confirm chr8 / chr19 hotspot signals.

Data source (frozen 2026-05-22):
  /big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn/
    step1_v3f_v5_v6_three_way/step1_master_four_way.tsv
  35,332 regions x 112 cols, HCC1395 4-way: baseline / V3F / V5 / V6.

Per-chromosome metrics (24 expected; this dataset has 22: chr1-22; no chrX/chrY):
  - n_regions    : total regions per chr
  - n_TP / n_FP  : TP/FP counts per chr
  - baseline_marker_count : regions where baseline_off_NG >= 3
  - V6_marker_count       : regions where V6_off_NG >= 3
  - baseline_marker_rate  : marker / n_regions (NG>=3 fraction)
  - V6_marker_rate        : marker / n_regions
  - baseline_marker_TP_purity : TP / total markers (NG>=3) for baseline
  - V6_marker_TP_purity       : same for V6
  - baseline_HP1_HP2_ratio    : sum(baseline_off_1) / sum(baseline_off_2)
  - V6_HP1_HP2_ratio          : sum(V6_off_1)       / sum(V6_off_2)
  - delta_marker_rate         : V6 - baseline
  - delta_marker_TP_purity    : V6 - baseline (F1 proxy; caveat: NOT actual caller F1)
  - delta_HP1_HP2_ratio       : V6 - baseline (closer to 1.0 = better balance)

Figure layout (1920x1080 16:9):
  Panel A (top, ~55% height) :
    heatmap 22 chr (cols) x 4 metric rows
      row 0: baseline marker count (log10 scale)
      row 1: V6 marker count       (log10 scale)
      row 2: delta marker rate     (V6 - baseline, diverging)
      row 3: delta TP purity       (V6 - baseline, diverging) -- F1 proxy
  Panel B (bottom, ~40% height):
    bar chart chr1-22 x delta marker TP purity (F1 proxy)
    red = V6 regression (counter-example candidate)
    green = V6 improvement
    annotate chr8 (LOH hotspot), chr19 (priority-bug hotspot),
    chr12 (counter), and top 3 +/- regressors.

CJK font: Droid Sans Fallback (matplotlib rcParams injected).

Outputs:
  - figures/C1_per_chr_V6_vs_baseline_grid.png
  - C1_per_chr_summary.tsv
  Copy of PNG to docs/.../target2/.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib import gridspec

# ---------- CJK font (per CLAUDE.md feedback_matplotlib_cjk_font_rule) ----------
plt.rcParams["font.sans-serif"] = [
    "Droid Sans Fallback", "Noto Sans CJK TC", "Noto Sans CJK SC", "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42

# ---------- Paths ----------
ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
SRC_TSV = ROOT / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_four_way.tsv"
AUDIT_DIR = ROOT / "research/methyl_augmented_filter_phase2/phase2_completeness_audit"
OUT_TSV = AUDIT_DIR / "C1_per_chr_summary.tsv"
OUT_PNG = AUDIT_DIR / "figures/C1_per_chr_V6_vs_baseline_grid.png"
COPY_PNG = ROOT / (
    "docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/"
    "assets/figures/target2/C1_per_chr_V6_vs_baseline_grid.png"
)


def chr_sort_key(c: str) -> int:
    """chr1, chr2, ..., chr22, chrX, chrY canonical ordering."""
    s = c.replace("chr", "")
    if s == "X":
        return 23
    if s == "Y":
        return 24
    try:
        return int(s)
    except ValueError:
        return 999


def main() -> int:
    print(f"[C1] loading {SRC_TSV}", flush=True)
    df = pd.read_csv(SRC_TSV, sep="\t", low_memory=False)
    print(f"[C1] loaded shape={df.shape}", flush=True)

    # Sanity checks
    assert "chr" in df.columns and "label" in df.columns
    assert "baseline_off_NG" in df.columns and "V6_off_NG" in df.columns
    assert "baseline_off_1" in df.columns and "baseline_off_2" in df.columns
    assert "V6_off_1" in df.columns and "V6_off_2" in df.columns

    # Coerce numerics
    for col in [
        "baseline_off_NG", "V6_off_NG",
        "baseline_off_1", "baseline_off_2",
        "V6_off_1", "V6_off_2",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    df["is_TP"] = (df["label"].astype(str).str.upper() == "TP").astype(int)
    df["baseline_marker"] = (df["baseline_off_NG"] >= 3).astype(int)
    df["V6_marker"] = (df["V6_off_NG"] >= 3).astype(int)

    # Per-chr aggregation
    grouped = df.groupby("chr", sort=False)
    rows: list[dict] = []
    for chrom, g in grouped:
        n_regions = len(g)
        n_TP = int(g["is_TP"].sum())
        n_FP = n_regions - n_TP

        base_mark = int(g["baseline_marker"].sum())
        v6_mark = int(g["V6_marker"].sum())
        base_rate = base_mark / n_regions if n_regions else np.nan
        v6_rate = v6_mark / n_regions if n_regions else np.nan

        # marker TP purity = TP fraction among markers (NG>=3 regions)
        base_mark_rows = g[g["baseline_marker"] == 1]
        v6_mark_rows = g[g["V6_marker"] == 1]
        base_tp_purity = (
            base_mark_rows["is_TP"].sum() / len(base_mark_rows)
            if len(base_mark_rows) else np.nan
        )
        v6_tp_purity = (
            v6_mark_rows["is_TP"].sum() / len(v6_mark_rows)
            if len(v6_mark_rows) else np.nan
        )

        # HP1:HP2 ratio (off-region read totals)
        b1, b2 = g["baseline_off_1"].sum(), g["baseline_off_2"].sum()
        v1, v2 = g["V6_off_1"].sum(), g["V6_off_2"].sum()
        base_ratio = b1 / b2 if b2 > 0 else np.nan
        v6_ratio = v1 / v2 if v2 > 0 else np.nan

        rows.append({
            "chr": chrom,
            "n_regions": n_regions,
            "n_TP": n_TP,
            "n_FP": n_FP,
            "baseline_marker_count": base_mark,
            "V6_marker_count": v6_mark,
            "baseline_marker_rate": base_rate,
            "V6_marker_rate": v6_rate,
            "baseline_marker_TP_purity": base_tp_purity,
            "V6_marker_TP_purity": v6_tp_purity,
            "baseline_HP1_HP2_ratio": base_ratio,
            "V6_HP1_HP2_ratio": v6_ratio,
            "delta_marker_rate": v6_rate - base_rate,
            "delta_marker_TP_purity": v6_tp_purity - base_tp_purity,
            "delta_HP1_HP2_ratio": v6_ratio - base_ratio,
        })

    summary = pd.DataFrame(rows)
    summary["sort_key"] = summary["chr"].map(chr_sort_key)
    summary = summary.sort_values("sort_key").drop(columns=["sort_key"]).reset_index(drop=True)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6f")
    print(f"[C1] wrote {OUT_TSV} ({len(summary)} chr rows)", flush=True)

    # ---------- Figure ----------
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    gs = gridspec.GridSpec(
        2, 1,
        height_ratios=[1.0, 0.95],
        hspace=0.40,
        left=0.06, right=0.985,
        top=0.92, bottom=0.07,
    )

    chrs = summary["chr"].tolist()
    n_chr = len(chrs)
    x = np.arange(n_chr)

    # ---- Panel A: heatmap ----
    ax_a = fig.add_subplot(gs[0])

    # Build 4 x n_chr matrix; mixed scales -> normalize per row for color, annotate raw value.
    row_labels = [
        "baseline marker count (NG≥3, log10)",
        "V6 marker count (NG≥3, log10)",
        "Δ marker rate (V6 − baseline)",
        "Δ marker TP purity (V6 − baseline)  [F1 proxy]",
    ]
    mat_raw = np.zeros((4, n_chr), dtype=float)
    mat_raw[0] = summary["baseline_marker_count"].values
    mat_raw[1] = summary["V6_marker_count"].values
    mat_raw[2] = summary["delta_marker_rate"].values
    mat_raw[3] = summary["delta_marker_TP_purity"].values

    # Color matrix: rows 0/1 sequential viridis on log10, rows 2/3 diverging RdBu_r centered 0
    color_mat = np.zeros_like(mat_raw)
    color_mat[0] = np.log10(mat_raw[0] + 1.0)
    color_mat[1] = np.log10(mat_raw[1] + 1.0)
    color_mat[2] = mat_raw[2]
    color_mat[3] = mat_raw[3]

    # Draw 4 rows individually so colormaps differ
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(chrs, rotation=45, ha="right", fontsize=10)
    ax_a.set_yticks(np.arange(4))
    ax_a.set_yticklabels(row_labels, fontsize=11)
    ax_a.set_xlim(-0.5, n_chr - 0.5)
    ax_a.set_ylim(-0.5, 3.5)
    ax_a.invert_yaxis()
    ax_a.set_title(
        "HCC1395 per-chromosome V6 vs baseline ISM full genome sweep  (Panel A: metric heatmap)",
        fontsize=13, fontweight="bold", pad=10,
    )

    # row 0 + 1 viridis log10
    cmap_seq = plt.get_cmap("viridis")
    seq_vmax = max(color_mat[0].max(), color_mat[1].max())
    for r in (0, 1):
        for j in range(n_chr):
            v = color_mat[r, j]
            face = cmap_seq(v / seq_vmax if seq_vmax > 0 else 0.0)
            ax_a.add_patch(plt.Rectangle((j - 0.5, r - 0.5), 1, 1, facecolor=face, edgecolor="white", lw=0.4))
            ax_a.text(j, r, f"{int(mat_raw[r, j])}", ha="center", va="center",
                      color="white" if (v / max(seq_vmax, 1e-9)) > 0.45 else "black",
                      fontsize=8)

    # rows 2/3 diverging RdBu_r centered 0
    cmap_div = plt.get_cmap("RdBu_r")
    for r in (2, 3):
        vmax = float(np.nanmax(np.abs(color_mat[r])))
        if vmax == 0 or np.isnan(vmax):
            vmax = 1.0
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        for j in range(n_chr):
            v = color_mat[r, j]
            if np.isnan(v):
                face = (0.85, 0.85, 0.85, 1.0)
                label_txt = "NA"
            else:
                face = cmap_div(norm(v))
                label_txt = f"{v:+.3f}"
            ax_a.add_patch(plt.Rectangle((j - 0.5, r - 0.5), 1, 1, facecolor=face, edgecolor="white", lw=0.4))
            ax_a.text(j, r, label_txt, ha="center", va="center",
                      color="black", fontsize=7.5)

    # Highlight chr8 + chr19 + chr12 columns
    def _col_idx(name: str) -> int | None:
        try:
            return chrs.index(name)
        except ValueError:
            return None

    for name, color in [("chr8", "#d62728"), ("chr19", "#ff7f0e"), ("chr12", "#1f77b4")]:
        idx = _col_idx(name)
        if idx is not None:
            ax_a.add_patch(plt.Rectangle((idx - 0.5, -0.5), 1, 4,
                                         fill=False, edgecolor=color, lw=2.0, zorder=10))

    # ---- Panel B: bar chart delta marker TP purity ----
    ax_b = fig.add_subplot(gs[1])
    delta_tp = summary["delta_marker_TP_purity"].values
    colors = ["#2ca02c" if v > 0 else "#d62728" if v < 0 else "#7f7f7f" for v in delta_tp]
    bars = ax_b.bar(x, delta_tp, color=colors, edgecolor="black", lw=0.5)
    ax_b.axhline(0, color="black", lw=0.8)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(chrs, rotation=45, ha="right", fontsize=10)
    ax_b.set_ylabel("Δ marker TP purity\n(V6 − baseline)  [F1 proxy]", fontsize=11)
    ax_b.set_title(
        "Panel B: per-chr ΔTP-purity within NG≥3 markers  "
        "(green = V6 improvement, red = V6 regression / counter-example candidate)",
        fontsize=12, fontweight="bold", pad=8,
    )
    ax_b.grid(True, axis="y", alpha=0.3)

    # Annotate top 3 + / - regressors
    order_pos = np.argsort(-delta_tp)
    order_neg = np.argsort(delta_tp)
    annotated = set()
    for idx in list(order_pos[:3]) + list(order_neg[:3]):
        if idx in annotated or np.isnan(delta_tp[idx]):
            continue
        annotated.add(idx)
        sign = "+" if delta_tp[idx] >= 0 else ""
        ax_b.annotate(
            f"{sign}{delta_tp[idx]:.3f}",
            xy=(idx, delta_tp[idx]),
            xytext=(0, 6 if delta_tp[idx] >= 0 else -16),
            textcoords="offset points",
            ha="center", fontsize=8.5, fontweight="bold",
        )

    # Highlight chr8 / chr19 / chr12 labels in red box
    for name, color in [("chr8", "#d62728"), ("chr19", "#ff7f0e"), ("chr12", "#1f77b4")]:
        idx = _col_idx(name)
        if idx is not None:
            ax_b.get_xticklabels()[idx].set_color(color)
            ax_b.get_xticklabels()[idx].set_fontweight("bold")

    # Footer caveat
    fig.text(
        0.06, 0.012,
        "Source: research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_four_way.tsv  |  "
        "Caveat: marker-TP-purity is an F1 proxy (no per-chr caller F1 column); "
        "highlighted chr8 (LOH hotspot) / chr19 (priority-bug hotspot) / chr12 (counter).",
        fontsize=8, color="#555555",
    )

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"[C1] wrote {OUT_PNG}", flush=True)

    # Copy to PI presentation assets
    COPY_PNG.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(OUT_PNG, COPY_PNG)
    print(f"[C1] copied to {COPY_PNG}", flush=True)

    # Console summary: top 3 +/- chr
    s = summary.dropna(subset=["delta_marker_TP_purity"]).copy()
    top_pos = s.nlargest(3, "delta_marker_TP_purity")[["chr", "delta_marker_TP_purity"]]
    top_neg = s.nsmallest(3, "delta_marker_TP_purity")[["chr", "delta_marker_TP_purity"]]
    print("\n[C1] Top 3 V6 improvements (Δ TP purity):")
    for _, r in top_pos.iterrows():
        print(f"  {r['chr']:>6s}  Δ={r['delta_marker_TP_purity']:+.4f}")
    print("\n[C1] Top 3 V6 regressions (Δ TP purity):")
    for _, r in top_neg.iterrows():
        print(f"  {r['chr']:>6s}  Δ={r['delta_marker_TP_purity']:+.4f}")
    print(f"\n[C1] n_chr={len(summary)}  total_regions={int(summary['n_regions'].sum())}")

    # chr8 / chr19 specific
    for name in ("chr8", "chr19", "chr12"):
        row = summary[summary["chr"] == name]
        if row.empty:
            continue
        r = row.iloc[0]
        print(f"[C1] {name}: n={r['n_regions']:.0f} "
              f"base_marker={r['baseline_marker_count']:.0f} "
              f"V6_marker={r['V6_marker_count']:.0f} "
              f"Δrate={r['delta_marker_rate']:+.4f} "
              f"ΔTPpurity={r['delta_marker_TP_purity']:+.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
