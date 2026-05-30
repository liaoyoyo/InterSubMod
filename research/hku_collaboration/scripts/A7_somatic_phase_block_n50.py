#!/usr/bin/env python3
"""
A7 — HKU collaboration somatic phase block (PS) N50 + size distribution
=========================================================================

Goal:
  Compute per-chromosome PS-tag-based phase block size distribution + N50 for
  the HCC1395 LongPhase --somaticMode tagged BAM, across 22 autosomes + chrX.

Inputs:
  BAM: HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam (indexed)
  Per-read PS (int) + HP (str e.g. "1", "2", "1-1", "2-1") tags.

Method (per chrom in parallel):
  1) Stream reads, filter MAPQ>=20, drop dup/supp/secondary.
  2) For each PS group, accumulate (min_pos, max_pos, n_reads,
     n_somatic_hp = count(HP in {'1-1','2-1'})).
  3) block_size_bp = max_pos - min_pos + 1
     somatic_density = n_somatic_hp / n_reads
  4) Per-chr stats: N50, median, IQR, max, n_blocks, coverage_ratio (sum(size) / chr_len).
  5) Global stats across all chr.

Outputs (under research/hku_collaboration/):
  data/A7_ps_block_sizes.tsv
  data/A7_per_chr_n50_summary.tsv
  data/A7_run.log
  figures/A7_phase_block_size_distribution.png (24-panel 4x6 hist)
  figures/A7_n50_per_chr_bar.png
  figures/A7_block_density_vs_chr_length.png

CJK font injection via _common.setup_cjk_font (Droid Sans Fallback).
"""
from __future__ import annotations

import csv
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

# Ensure local _common.py is importable
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from _common import BAM_PATH, FIG_DIR, DATA_DIR, setup_cjk_font  # noqa: E402

# --- HG38 lengths (autosomes + chrX), aligned with build_s2_s3_phase_block_hp_analysis.py ---
HG38_CHR_LENGTHS: Dict[str, int] = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895,
}
# Canonical chr order for plotting / TSV: chr1..chr22, chrX
CHROMS_ORDER: List[str] = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

SOMATIC_HP_VALUES = {"1-1", "2-1"}  # HP1_1 / HP2_1 -> somatic-side haplotype
MAPQ_MIN = 20


# ============================================================
# Per-chr worker — return list[dict] with one row per PS block
# ============================================================
def process_chrom(args: Tuple[str, str]) -> Tuple[str, List[dict], dict]:
    """Stream BAM for one chrom, aggregate per-PS block stats.

    Returns:
        chrom, list of block dicts, error/info dict.
    """
    chrom, bam_path = args
    t0 = time.time()
    blocks_min: Dict[int, int] = {}
    blocks_max: Dict[int, int] = {}
    blocks_nreads: Dict[int, int] = defaultdict(int)
    blocks_nsomatic: Dict[int, int] = defaultdict(int)

    n_reads_seen = 0
    n_reads_ps = 0

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return chrom, [], {"error": f"open BAM failed: {e}"}

    try:
        for rec in bam.fetch(chrom):
            n_reads_seen += 1
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.is_duplicate:
                continue
            if rec.mapping_quality < MAPQ_MIN:
                continue
            if not rec.has_tag("PS"):
                continue
            try:
                ps = int(rec.get_tag("PS"))
            except (ValueError, TypeError):
                continue
            if ps <= 0:
                continue

            pos = rec.reference_start  # 0-based
            end = rec.reference_end if rec.reference_end is not None else pos + 1

            cur_min = blocks_min.get(ps)
            if cur_min is None or pos < cur_min:
                blocks_min[ps] = pos
            cur_max = blocks_max.get(ps)
            if cur_max is None or end > cur_max:
                blocks_max[ps] = end

            blocks_nreads[ps] += 1
            n_reads_ps += 1

            if rec.has_tag("HP"):
                hp_val = rec.get_tag("HP")
                if isinstance(hp_val, (int, float)):
                    hp_str = str(int(hp_val))
                else:
                    hp_str = str(hp_val)
                if hp_str in SOMATIC_HP_VALUES:
                    blocks_nsomatic[ps] += 1
    finally:
        bam.close()

    rows: List[dict] = []
    for ps, mn in blocks_min.items():
        mx = blocks_max[ps]
        nr = blocks_nreads[ps]
        nsom = blocks_nsomatic.get(ps, 0)
        size = mx - mn  # half-open: end - start
        rows.append({
            "chrom": chrom,
            "ps_id": ps,
            "n_reads": nr,
            "block_size_bp": size,
            "somatic_density": (nsom / nr) if nr > 0 else 0.0,
        })

    info = {
        "elapsed_s": time.time() - t0,
        "n_reads_seen": n_reads_seen,
        "n_reads_ps": n_reads_ps,
        "n_ps_blocks": len(rows),
    }
    return chrom, rows, info


# ============================================================
# N50 + per-chr summary
# ============================================================
def calc_n50(sizes: List[int]) -> int:
    if not sizes:
        return 0
    s = sorted(sizes, reverse=True)
    total = sum(s)
    half = total / 2.0
    cum = 0
    for v in s:
        cum += v
        if cum >= half:
            return v
    return s[-1]


def per_chr_summary(chrom: str, rows: List[dict]) -> dict:
    sizes = [r["block_size_bp"] for r in rows if r["block_size_bp"] > 0]
    chr_len = HG38_CHR_LENGTHS.get(chrom, 0)
    if not sizes:
        return {
            "chrom": chrom,
            "chr_length": chr_len,
            "n_ps_blocks": 0,
            "median_block_size": 0,
            "q25_block_size": 0,
            "q75_block_size": 0,
            "n50": 0,
            "max_block": 0,
            "max_block_ps": "",
            "coverage_ratio": 0.0,
            "total_block_span_bp": 0,
        }
    arr = np.array(sizes, dtype=np.int64)
    q25, med, q75 = np.percentile(arr, [25, 50, 75])
    # find PS id of max block
    max_idx = int(np.argmax(arr))
    # walk rows to find matching ps_id (first match)
    max_block = int(arr[max_idx])
    max_block_ps = ""
    for r in rows:
        if r["block_size_bp"] == max_block:
            max_block_ps = str(r["ps_id"])
            break
    n50 = calc_n50(sizes)
    total_span = int(arr.sum())
    cov = (total_span / chr_len) if chr_len > 0 else 0.0
    return {
        "chrom": chrom,
        "chr_length": chr_len,
        "n_ps_blocks": len(sizes),
        "median_block_size": int(med),
        "q25_block_size": int(q25),
        "q75_block_size": int(q75),
        "n50": int(n50),
        "max_block": max_block,
        "max_block_ps": max_block_ps,
        "coverage_ratio": round(cov, 4),
        "total_block_span_bp": total_span,
    }


# ============================================================
# Plots
# ============================================================
def plot_size_distribution(chrom_rows: Dict[str, List[dict]],
                           chrom_summary: Dict[str, dict],
                           out_path: Path,
                           x_max_bp: int) -> None:
    """24-panel grid (4 row x 6 col) per-chr histogram of block sizes."""
    fig, axes = plt.subplots(4, 6, figsize=(20, 12), constrained_layout=True)
    axes = axes.flatten()

    # shared bins
    bins = np.linspace(0, x_max_bp, 41)

    for i, chrom in enumerate(CHROMS_ORDER):
        ax = axes[i]
        rows = chrom_rows.get(chrom, [])
        sizes = np.array([r["block_size_bp"] for r in rows if r["block_size_bp"] > 0],
                         dtype=np.int64)
        if len(sizes):
            ax.hist(sizes, bins=bins, color="#3b7dd8", edgecolor="black", linewidth=0.3)
        s = chrom_summary[chrom]
        ax.set_title(
            f"{chrom}\nN50={s['n50']/1000:.0f} kb  med={s['median_block_size']/1000:.0f} kb  n={s['n_ps_blocks']}",
            fontsize=10,
        )
        ax.set_xlim(0, x_max_bp)
        ax.tick_params(labelsize=8)
        if i % 6 == 0:
            ax.set_ylabel("PS block count", fontsize=9)
        if i >= 18:
            ax.set_xlabel("Block size (bp)", fontsize=9)
        # x ticks in Mb
        xt = np.linspace(0, x_max_bp, 5)
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{int(v/1e6)}M" if v >= 1e6 else f"{int(v/1e3)}k" for v in xt],
                            fontsize=7, rotation=0)

    # hide unused panels (24 panels, 23 chroms used -> 1 spare)
    for j in range(len(CHROMS_ORDER), len(axes)):
        axes[j].axis("off")

    fig.suptitle(
        "A7 — Somatic phase block (PS) size distribution per chromosome\n"
        "HCC1395 Tmode tagged ClairS pileup v040 (LongPhase --somaticMode), MAPQ>=20",
        fontsize=13,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_n50_per_chr_bar(chrom_summary: Dict[str, dict], out_path: Path) -> None:
    """Bar chart of N50 per chr (chr1..chr22, chrX), colour by chr length."""
    chroms = CHROMS_ORDER
    n50_vals = [chrom_summary[c]["n50"] for c in chroms]
    chr_lens = [HG38_CHR_LENGTHS[c] for c in chroms]
    n_blocks = [chrom_summary[c]["n_ps_blocks"] for c in chroms]

    # Global N50 across all blocks pooled
    all_sizes: List[int] = []
    for c in chroms:
        # Reconstruct from summary not possible -> store earlier; we pass via attr.
        pass
    # we will compute global later in main and pass directly; here use mean N50 as proxy if not set
    fig, ax = plt.subplots(1, 1, figsize=(14, 6), constrained_layout=True)
    norm = plt.Normalize(vmin=min(chr_lens), vmax=max(chr_lens))
    cmap = plt.cm.viridis
    colors = [cmap(norm(l)) for l in chr_lens]

    x = np.arange(len(chroms))
    bars = ax.bar(x, [v / 1000 for v in n50_vals], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("PS block N50 (kb)", fontsize=11)
    ax.set_title("A7 — Somatic phase block N50 per chromosome (HCC1395 Tmode tagged)",
                  fontsize=12)

    # annotate n_blocks on top of bars
    for xi, val, nb in zip(x, n50_vals, n_blocks):
        ax.text(xi, val / 1000, f"n={nb}", ha="center", va="bottom", fontsize=7)

    # colour bar for chr length
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7)
    cbar.set_label("Chromosome length (bp)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.grid(axis="y", alpha=0.3)

    # global N50 horizontal line is added in main via setattr trick; use attribute
    global_n50 = getattr(plot_n50_per_chr_bar, "_global_n50", None)
    if global_n50:
        ax.axhline(global_n50 / 1000, color="red", ls="--", lw=1.2,
                   label=f"Global N50 = {global_n50/1000:.0f} kb")
        ax.legend(loc="upper right", fontsize=10)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_block_density_vs_chr_length(chrom_summary: Dict[str, dict], out_path: Path) -> None:
    chroms = CHROMS_ORDER
    chr_lens = np.array([HG38_CHR_LENGTHS[c] for c in chroms], dtype=np.float64)
    cov = np.array([chrom_summary[c]["coverage_ratio"] for c in chroms], dtype=np.float64)

    fig, ax = plt.subplots(1, 1, figsize=(10, 7), constrained_layout=True)
    ax.scatter(chr_lens / 1e6, cov, s=80, c="#3b7dd8", edgecolor="black", linewidth=0.5)
    for c, x, y in zip(chroms, chr_lens / 1e6, cov):
        ax.annotate(c, (x, y), fontsize=8, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Chromosome length (Mb)", fontsize=11)
    ax.set_ylabel("PS coverage ratio (sum block span / chr length)", fontsize=11)
    ax.set_title("A7 — PS block coverage vs chromosome length (HCC1395 Tmode)", fontsize=12)
    ax.grid(alpha=0.3)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Main
# ============================================================
def main() -> int:
    setup_cjk_font()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    log_path = DATA_DIR / "A7_run.log"
    log_fh = open(log_path, "w")

    def log(msg: str) -> None:
        print(msg, file=sys.stderr)
        log_fh.write(msg + "\n")
        log_fh.flush()

    t_start = time.time()
    log(f"[A7] start  BAM={BAM_PATH}")
    log(f"[A7] chroms={len(CHROMS_ORDER)}  MAPQ>={MAPQ_MIN}  somatic_HP={sorted(SOMATIC_HP_VALUES)}")

    # parallel per chr
    n_workers = min(len(CHROMS_ORDER), max(1, mp.cpu_count() - 2))
    log(f"[A7] spawning {n_workers} workers")

    chrom_rows: Dict[str, List[dict]] = {}
    chrom_info: Dict[str, dict] = {}

    work = [(c, BAM_PATH) for c in CHROMS_ORDER]
    with mp.Pool(processes=n_workers) as pool:
        for chrom, rows, info in pool.imap_unordered(process_chrom, work):
            chrom_rows[chrom] = rows
            chrom_info[chrom] = info
            if info.get("error"):
                log(f"  [WARN] {chrom}: {info['error']}")
            else:
                log(f"  [done] {chrom}: {info['n_ps_blocks']} blocks, "
                    f"{info['n_reads_ps']}/{info['n_reads_seen']} reads w/ PS, "
                    f"t={info['elapsed_s']:.1f}s")

    # ----- per-chr summary -----
    chrom_summary: Dict[str, dict] = {}
    for c in CHROMS_ORDER:
        chrom_summary[c] = per_chr_summary(c, chrom_rows.get(c, []))

    # ----- write TSV 1: per-block sizes -----
    blocks_tsv = DATA_DIR / "A7_ps_block_sizes.tsv"
    with open(blocks_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["chrom", "ps_id", "n_reads",
                                         "block_size_bp", "somatic_density"],
                            delimiter="\t")
        w.writeheader()
        for c in CHROMS_ORDER:
            # sort within chrom by ps_id for stable output
            for r in sorted(chrom_rows.get(c, []), key=lambda x: x["ps_id"]):
                w.writerow({
                    "chrom": r["chrom"],
                    "ps_id": r["ps_id"],
                    "n_reads": r["n_reads"],
                    "block_size_bp": r["block_size_bp"],
                    "somatic_density": f"{r['somatic_density']:.4f}",
                })
    log(f"[A7] wrote {blocks_tsv}")

    # ----- write TSV 2: per-chr N50 summary -----
    n50_tsv = DATA_DIR / "A7_per_chr_n50_summary.tsv"
    with open(n50_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "chrom", "chr_length", "n_ps_blocks", "median_block_size",
            "q25_block_size", "q75_block_size", "n50",
            "max_block", "max_block_ps", "coverage_ratio", "total_block_span_bp",
        ], delimiter="\t")
        w.writeheader()
        for c in CHROMS_ORDER:
            w.writerow(chrom_summary[c])
    log(f"[A7] wrote {n50_tsv}")

    # ----- global stats -----
    all_sizes: List[int] = []
    for c in CHROMS_ORDER:
        all_sizes.extend([r["block_size_bp"] for r in chrom_rows.get(c, [])
                          if r["block_size_bp"] > 0])
    global_n50 = calc_n50(all_sizes)
    global_median = int(np.median(all_sizes)) if all_sizes else 0
    global_mean = int(np.mean(all_sizes)) if all_sizes else 0
    total_blocks = len(all_sizes)
    # max block + its chr
    max_chrom = ""
    max_block_size = 0
    max_block_ps = ""
    for c in CHROMS_ORDER:
        s = chrom_summary[c]
        if s["max_block"] > max_block_size:
            max_block_size = s["max_block"]
            max_chrom = c
            max_block_ps = s["max_block_ps"]

    log("")
    log("=" * 60)
    log("[A7] GLOBAL SUMMARY")
    log(f"  total PS blocks       : {total_blocks:,}")
    log(f"  global N50            : {global_n50:,} bp ({global_n50/1000:.1f} kb)")
    log(f"  median block size     : {global_median:,} bp ({global_median/1000:.1f} kb)")
    log(f"  mean block size       : {global_mean:,} bp ({global_mean/1000:.1f} kb)")
    log(f"  max block             : {max_block_size:,} bp ({max_block_size/1e6:.2f} Mb) "
        f"on {max_chrom} (PS={max_block_ps})")
    log("=" * 60)

    # ----- plots -----
    # x_max for histogram panels: use 95th percentile of sizes (cap visual)
    if all_sizes:
        x_max_bp = int(np.percentile(all_sizes, 99))
        # at least 1 Mb visual to avoid super thin histograms
        x_max_bp = max(x_max_bp, 1_000_000)
        # round up to nearest 0.5 Mb
        x_max_bp = int(np.ceil(x_max_bp / 500_000) * 500_000)
    else:
        x_max_bp = 1_000_000

    fig1 = FIG_DIR / "A7_phase_block_size_distribution.png"
    plot_size_distribution(chrom_rows, chrom_summary, fig1, x_max_bp=x_max_bp)
    log(f"[A7] wrote {fig1}  (x_max={x_max_bp/1e6:.1f} Mb)")

    # pass global N50 to bar plot via attribute
    plot_n50_per_chr_bar._global_n50 = global_n50  # type: ignore[attr-defined]
    fig2 = FIG_DIR / "A7_n50_per_chr_bar.png"
    plot_n50_per_chr_bar(chrom_summary, fig2)
    log(f"[A7] wrote {fig2}")

    fig3 = FIG_DIR / "A7_block_density_vs_chr_length.png"
    plot_block_density_vs_chr_length(chrom_summary, fig3)
    log(f"[A7] wrote {fig3}")

    log("")
    log(f"[A7] total wall: {time.time() - t_start:.1f}s")

    # one-liner stdout
    print(
        f"A7 done: global N50={global_n50:,} bp ({global_n50/1000:.1f} kb), "
        f"median={global_median:,} bp, max block on {max_chrom} "
        f"({max_block_size/1e6:.2f} Mb), total blocks={total_blocks:,}"
    )

    log_fh.close()
    # Return for embedding consumers
    return 0


if __name__ == "__main__":
    sys.exit(main())
