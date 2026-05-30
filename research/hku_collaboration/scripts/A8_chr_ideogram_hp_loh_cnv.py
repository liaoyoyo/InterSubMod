#!/usr/bin/env python3
"""A8 — 24-panel chromosome ideogram + LOH/CNV shading + HP tag stacked bar per 10 Mb bin.

Inputs
------
- BAM: /big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam
- BED: /big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed (chrom, start, end, type)

Scope: 22 autosomes + chrX (per HG38_CHR_LENGTHS).

Outputs (under research/hku_collaboration/)
  figures/A8_chr_ideogram_hp_tag.png         — 24-panel grid
  figures/A8_per_chr_hp_summary_table.png    — normalized stacked HP composition
  figures/A8_chr_loh_cnv_coverage.png        — per-chr LOH & CNV coverage %
  data/A8_per_chr_hp_loh_cnv.tsv             — bin-level table
  data/A8_per_chr_summary.tsv                — per-chr summary
"""
from __future__ import annotations

import bisect
import os
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import pysam

# Local
sys.path.insert(0, str(Path(__file__).parent))
from _common import BAM_PATH, LOH_BED, DATA_DIR, FIG_DIR, setup_cjk_font, log_timing  # noqa: E402

# ============================================================
# Config
# ============================================================
HG38_CHR_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895,
}
CHROM_ORDER = list(HG38_CHR_LENGTHS.keys())  # chr1..chr22, chrX

BIN_SIZE = 10_000_000  # 10 Mb
HP_KEYS = ["1", "2", "1-1", "2-1", "3", "no_HP"]
HP_LABELS = {"1": "HP1", "2": "HP2", "1-1": "HP1-1", "2-1": "HP2-1", "3": "HP3", "no_HP": "no-HP"}
HP_COLORS = {
    "1": "#1f77b4",     # blue
    "2": "#d62728",     # red
    "1-1": "#9ecae1",   # light blue
    "2-1": "#fcae91",   # light red
    "3": "#7f7f7f",     # grey
    "no_HP": "#dddddd", # very light grey
}

# CNV/LOH shading colors
COL_LOH = "#ff9999"      # light red
COL_GAIN = "#cce5ff"     # light blue
COL_LOSS = "#ccffcc"     # light green
COL_CHR_BAR = "#bdbdbd"  # grey for chromosome bar


# ============================================================
# BED loader (full BED including all 3 types)
# ============================================================
def load_bed_typed(bed_path: str) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """Return dict[chrom][type] = sorted list of (start, end). Types: loh / gain / loss."""
    regs: dict[str, dict[str, list[tuple[int, int]]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            typ = parts[3].strip().lower()
            if typ not in ("loh", "gain", "loss"):
                continue
            regs.setdefault(chrom, {}).setdefault(typ, []).append((s, e))
    for c in regs:
        for t in regs[c]:
            regs[c][t].sort()
    return regs


def intervals_overlap_any(intervals: list[tuple[int, int]], qs: int, qe: int) -> bool:
    """Boolean overlap of [qs, qe) vs sorted intervals."""
    if not intervals:
        return False
    starts = [iv[0] for iv in intervals]
    idx = bisect.bisect_right(starts, qe) - 1
    while idx >= 0:
        s, e = intervals[idx]
        if e <= qs:
            return False
        if s < qe and e > qs:
            return True
        idx -= 1
    return False


def intervals_coverage(intervals: list[tuple[int, int]], chr_length: int) -> int:
    """Total covered bp (merged) within [0, chr_length)."""
    if not intervals:
        return 0
    merged: list[tuple[int, int]] = []
    for s, e in sorted(intervals):
        s = max(0, s)
        e = min(chr_length, e)
        if e <= s:
            continue
        if not merged or s > merged[-1][1]:
            merged.append((s, e))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
    return sum(e - s for s, e in merged)


# ============================================================
# Per-chr worker: stream BAM, bin reads, count HP tags
# ============================================================
def _normalize_hp(val) -> str:
    """Normalize HP tag value to one of HP_KEYS, else 'no_HP'."""
    if val is None:
        return "no_HP"
    s = str(val)
    if s in ("1", "2", "1-1", "2-1", "3"):
        return s
    return "no_HP"


def process_chrom(args) -> tuple[str, list[dict]]:
    """Stream one chrom, bin HP counts. Returns (chrom, list of bin dicts)."""
    chrom, chr_length, bam_path = args
    n_bins = (chr_length + BIN_SIZE - 1) // BIN_SIZE
    counters: list[Counter] = [Counter() for _ in range(n_bins)]

    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for r in bam.fetch(chrom, 0, chr_length):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                continue
            # Use read midpoint to assign to one bin (avoid double counting reads
            # spanning bin boundaries — typical ONT read length << 10 Mb).
            start = r.reference_start
            end = r.reference_end if r.reference_end is not None else start + 1
            mid = (start + end) // 2
            if mid < 0 or mid >= chr_length:
                continue
            bin_idx = mid // BIN_SIZE
            if bin_idx >= n_bins:
                bin_idx = n_bins - 1
            tag = None
            try:
                if r.has_tag("HP"):
                    tag = r.get_tag("HP")
            except KeyError:
                tag = None
            counters[bin_idx][_normalize_hp(tag)] += 1
    finally:
        bam.close()

    rows = []
    for i, c in enumerate(counters):
        bs = i * BIN_SIZE
        be = min((i + 1) * BIN_SIZE, chr_length)
        total = sum(c.values())
        rows.append({
            "chrom": chrom,
            "bin_start": bs,
            "bin_end": be,
            "hp1": c.get("1", 0),
            "hp2": c.get("2", 0),
            "hp1_1": c.get("1-1", 0),
            "hp2_1": c.get("2-1", 0),
            "hp3": c.get("3", 0),
            "no_hp": c.get("no_HP", 0),
            "total_reads": total,
        })
    return chrom, rows


# ============================================================
# Figure 1 — 24-panel ideogram
# ============================================================
def plot_ideogram_grid(df_bins: pd.DataFrame, bed_regs: dict, out_path: Path) -> None:
    """4 row x 6 col grid. Each panel: chr bar + LOH/CNV shading + HP stacked bars."""
    n_rows, n_cols = 4, 6
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 14))
    axes = axes.flatten()

    max_total_global = int(df_bins["total_reads"].max()) if len(df_bins) else 1

    for idx, chrom in enumerate(CHROM_ORDER):
        ax = axes[idx]
        chr_len = HG38_CHR_LENGTHS[chrom]
        sub = df_bins[df_bins["chrom"] == chrom].sort_values("bin_start").reset_index(drop=True)
        total_reads_chr = int(sub["total_reads"].sum())

        # X axis in Mb
        x_max_mb = chr_len / 1e6

        # Top zone (chromosome bar) at y = max_total_global * 1.10..1.18
        bar_top = max_total_global * 1.18
        bar_bot = max_total_global * 1.10
        ax.add_patch(mpatches.Rectangle(
            (0, bar_bot), x_max_mb, bar_top - bar_bot,
            facecolor=COL_CHR_BAR, edgecolor="black", linewidth=0.5,
        ))

        # Middle zone: CNV/LOH shading from y=0 to y=max_total_global * 1.05
        shade_top = max_total_global * 1.05
        # gain first (lowest z), then loss, then loh (most prominent)
        for typ, color in (("gain", COL_GAIN), ("loss", COL_LOSS), ("loh", COL_LOH)):
            for s, e in bed_regs.get(chrom, {}).get(typ, []):
                s_mb = max(0, s) / 1e6
                e_mb = min(chr_len, e) / 1e6
                if e_mb <= s_mb:
                    continue
                ax.add_patch(mpatches.Rectangle(
                    (s_mb, 0), e_mb - s_mb, shade_top,
                    facecolor=color, edgecolor="none", alpha=0.45, zorder=1,
                ))

        # Bottom zone: HP stacked bars per 10 Mb bin
        if len(sub):
            bin_widths_mb = (sub["bin_end"] - sub["bin_start"]) / 1e6
            x_left_mb = sub["bin_start"] / 1e6
            x_centers_mb = (sub["bin_start"] + sub["bin_end"]) / 2 / 1e6
            bottoms = np.zeros(len(sub))
            for hp in HP_KEYS:
                col = {"1": "hp1", "2": "hp2", "1-1": "hp1_1",
                       "2-1": "hp2_1", "3": "hp3", "no_HP": "no_hp"}[hp]
                vals = sub[col].values.astype(float)
                ax.bar(
                    x_centers_mb, vals, width=bin_widths_mb * 0.92,
                    bottom=bottoms, color=HP_COLORS[hp], edgecolor="none",
                    linewidth=0, zorder=3,
                )
                bottoms = bottoms + vals

        ax.set_xlim(0, x_max_mb)
        ax.set_ylim(0, max_total_global * 1.22)
        ax.set_title(f"{chrom}  (n={total_reads_chr:,})", fontsize=9)
        ax.set_xlabel("Position (Mb)", fontsize=7)
        if idx % n_cols == 0:
            ax.set_ylabel("Reads per 10 Mb", fontsize=7)
        ax.tick_params(axis="both", labelsize=6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Hide unused panels (24 chr fit exactly into 4x6)
    for j in range(len(CHROM_ORDER), n_rows * n_cols):
        axes[j].axis("off")

    # Legend (HP keys + CNV/LOH shading) in a separate figure-level region
    hp_handles = [mpatches.Patch(color=HP_COLORS[k], label=HP_LABELS[k]) for k in HP_KEYS]
    shade_handles = [
        mpatches.Patch(facecolor=COL_LOH, alpha=0.45, label="LOH"),
        mpatches.Patch(facecolor=COL_GAIN, alpha=0.45, label="CNV gain"),
        mpatches.Patch(facecolor=COL_LOSS, alpha=0.45, label="CNV loss"),
        mpatches.Patch(facecolor=COL_CHR_BAR, label="chr bar"),
    ]
    fig.legend(
        handles=hp_handles + shade_handles,
        loc="lower center", ncol=10, fontsize=9, frameon=False,
        bbox_to_anchor=(0.5, -0.005),
    )

    fig.suptitle(
        "A8 — HCC1395 T-mode chromosome ideogram: HP tag composition per 10 Mb bin "
        "with SEQC2 LOH / CNV shading (22 autosomes + chrX)",
        fontsize=13, y=0.995,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.98))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Figure 2 — per-chr HP composition normalized stacked bar
# ============================================================
def plot_per_chr_summary(df_summary: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, 6))
    chroms = df_summary["chrom"].tolist()
    x = np.arange(len(chroms))
    bottoms = np.zeros(len(chroms))
    for hp in HP_KEYS:
        col = f"{hp.replace('-', '_').replace('no_HP', 'no_hp').replace('1', '1').lower()}_pct"
        # Simpler mapping
        col = {"1": "hp1_pct", "2": "hp2_pct", "1-1": "hp1_1_pct",
               "2-1": "hp2_1_pct", "3": "hp3_pct", "no_HP": "no_hp_pct"}[hp]
        vals = df_summary[col].values.astype(float)
        ax.bar(x, vals, bottom=bottoms, color=HP_COLORS[hp], label=HP_LABELS[hp], edgecolor="white", linewidth=0.3)
        bottoms = bottoms + vals
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Fraction of reads (normalized to 1)")
    ax.set_ylim(0, 1.0)
    ax.set_title("A8 — HCC1395 T-mode HP tag composition per chromosome (normalized)")
    ax.legend(loc="upper right", fontsize=8, ncol=6, frameon=False, bbox_to_anchor=(1.0, 1.12))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Figure 3 — per-chr LOH + CNV coverage %
# ============================================================
def plot_loh_cnv_coverage(df_summary: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, 6))
    chroms = df_summary["chrom"].tolist()
    x = np.arange(len(chroms))
    w = 0.4
    ax.bar(x - w/2, df_summary["loh_coverage_pct"], width=w, color=COL_LOH, edgecolor="darkred",
           linewidth=0.4, label="LOH coverage %")
    ax.bar(x + w/2, df_summary["cnv_coverage_pct"], width=w, color="#5b9bd5", edgecolor="navy",
           linewidth=0.4, label="CNV (gain+loss) coverage %")
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("% of chromosome covered")
    ax.set_title("A8 — Per-chromosome LOH and CNV coverage (HCC1395, SEQC2 truth)")
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, max(105, df_summary[["loh_coverage_pct", "cnv_coverage_pct"]].values.max() + 5))
    # Annotate top coverage bars
    for i, c in enumerate(chroms):
        v = df_summary.iloc[i]["loh_coverage_pct"]
        if v >= 20:
            ax.text(i - w/2, v + 1, f"{v:.0f}", ha="center", va="bottom", fontsize=7, color="darkred")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Main
# ============================================================
def main():
    setup_cjk_font()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[A8] BAM={BAM_PATH}", flush=True)
    print(f"[A8] BED={LOH_BED}", flush=True)
    print(f"[A8] chroms: {len(CHROM_ORDER)} (22 autosomes + chrX), bin={BIN_SIZE/1e6:.0f} Mb", flush=True)

    # ---- Load BED (all types) ----
    t0 = time.time()
    bed_regs = load_bed_typed(LOH_BED)
    log_timing("load_bed_typed", t0, time.time())
    n_loh = sum(len(v.get("loh", [])) for v in bed_regs.values())
    n_gain = sum(len(v.get("gain", [])) for v in bed_regs.values())
    n_loss = sum(len(v.get("loss", [])) for v in bed_regs.values())
    print(f"[A8] BED loaded: loh={n_loh}, gain={n_gain}, loss={n_loss}", flush=True)

    # ---- Stream BAM in parallel by chrom ----
    t0 = time.time()
    tasks = [(c, HG38_CHR_LENGTHS[c], BAM_PATH) for c in CHROM_ORDER]
    n_workers = min(len(tasks), 12)  # 12 chrom workers, BAM IO bound
    print(f"[A8] streaming BAM with {n_workers} workers...", flush=True)
    with Pool(processes=n_workers) as pool:
        results = pool.map(process_chrom, tasks)
    log_timing("BAM streaming (all 24 chrom)", t0, time.time())

    # ---- Build bin-level df ----
    bin_rows: list[dict] = []
    for chrom, rows in sorted(results, key=lambda x: CHROM_ORDER.index(x[0])):
        for row in rows:
            bs, be = row["bin_start"], row["bin_end"]
            loh_iv = bed_regs.get(chrom, {}).get("loh", [])
            gain_iv = bed_regs.get(chrom, {}).get("gain", [])
            loss_iv = bed_regs.get(chrom, {}).get("loss", [])
            row["in_loh"] = int(intervals_overlap_any(loh_iv, bs, be))
            row["in_cnv_gain"] = int(intervals_overlap_any(gain_iv, bs, be))
            row["in_cnv_loss"] = int(intervals_overlap_any(loss_iv, bs, be))
            bin_rows.append(row)
    df_bins = pd.DataFrame(bin_rows)
    df_bins = df_bins[[
        "chrom", "bin_start", "bin_end",
        "hp1", "hp2", "hp1_1", "hp2_1", "hp3", "no_hp",
        "in_loh", "in_cnv_gain", "in_cnv_loss", "total_reads",
    ]]

    bin_tsv = DATA_DIR / "A8_per_chr_hp_loh_cnv.tsv"
    df_bins.to_csv(bin_tsv, sep="\t", index=False)
    print(f"[A8] wrote {bin_tsv} ({len(df_bins)} rows)", flush=True)

    # ---- Per-chr summary ----
    summary_rows = []
    for chrom in CHROM_ORDER:
        chr_len = HG38_CHR_LENGTHS[chrom]
        sub = df_bins[df_bins["chrom"] == chrom]
        total = int(sub["total_reads"].sum())
        if total == 0:
            pcts = {f"{k}_pct": 0.0 for k in ["hp1", "hp2", "hp1_1", "hp2_1", "hp3", "no_hp"]}
        else:
            pcts = {
                "hp1_pct": sub["hp1"].sum() / total,
                "hp2_pct": sub["hp2"].sum() / total,
                "hp1_1_pct": sub["hp1_1"].sum() / total,
                "hp2_1_pct": sub["hp2_1"].sum() / total,
                "hp3_pct": sub["hp3"].sum() / total,
                "no_hp_pct": sub["no_hp"].sum() / total,
            }
        loh_bp = intervals_coverage(bed_regs.get(chrom, {}).get("loh", []), chr_len)
        gain_bp = intervals_coverage(bed_regs.get(chrom, {}).get("gain", []), chr_len)
        loss_bp = intervals_coverage(bed_regs.get(chrom, {}).get("loss", []), chr_len)
        # CNV total coverage = gain U loss (merge)
        cnv_iv_all = bed_regs.get(chrom, {}).get("gain", []) + bed_regs.get(chrom, {}).get("loss", [])
        cnv_bp = intervals_coverage(cnv_iv_all, chr_len)
        summary_rows.append({
            "chrom": chrom,
            "chr_length": chr_len,
            "total_reads": total,
            **pcts,
            "loh_coverage_bp": loh_bp,
            "cnv_gain_coverage_bp": gain_bp,
            "cnv_loss_coverage_bp": loss_bp,
            "loh_coverage_pct": 100.0 * loh_bp / chr_len,
            "cnv_coverage_pct": 100.0 * cnv_bp / chr_len,
        })
    df_summary = pd.DataFrame(summary_rows)
    df_summary["chrom"] = pd.Categorical(df_summary["chrom"], categories=CHROM_ORDER, ordered=True)
    df_summary = df_summary.sort_values("chrom").reset_index(drop=True)

    summary_tsv = DATA_DIR / "A8_per_chr_summary.tsv"
    df_summary.to_csv(summary_tsv, sep="\t", index=False)
    print(f"[A8] wrote {summary_tsv} ({len(df_summary)} rows)", flush=True)

    # ---- Figures ----
    t0 = time.time()
    plot_ideogram_grid(df_bins, bed_regs, FIG_DIR / "A8_chr_ideogram_hp_tag.png")
    log_timing("Figure 1 (24-panel ideogram)", t0, time.time())

    t0 = time.time()
    plot_per_chr_summary(df_summary, FIG_DIR / "A8_per_chr_hp_summary_table.png")
    log_timing("Figure 2 (HP composition summary)", t0, time.time())

    t0 = time.time()
    plot_loh_cnv_coverage(df_summary, FIG_DIR / "A8_chr_loh_cnv_coverage.png")
    log_timing("Figure 3 (LOH/CNV coverage)", t0, time.time())

    # ---- Console quick stats for finding writeup ----
    print("\n[A8] === Quick stats for §10 narrative ===", flush=True)
    print(df_summary[["chrom", "total_reads", "hp1_pct", "hp2_pct", "no_hp_pct",
                      "loh_coverage_pct", "cnv_coverage_pct"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
