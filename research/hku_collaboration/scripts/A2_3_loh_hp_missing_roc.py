#!/usr/bin/env python3
"""Analysis 3: LOH HP-missing rate ROC vs SEQC2 truth.

For 1kb bins on chr1 + chr8:
  - hp_missing_rate = reads_without_HP / total_reads
  - label = bin overlaps SEQC2 LOH region
ROC + AUC.

Output:
  - data/A2_3_bin_hp_missing.tsv
  - figures/A2_3_loh_hp_missing_vs_seqc2.png
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from sklearn.metrics import roc_auc_score, roc_curve

sys.path.insert(0, str(Path(__file__).parent))
from _common import (
    BAM_PATH, DATA_DIR, FIG_DIR,
    load_loh_regions, interval_overlap, setup_cjk_font, log_timing,
)
import matplotlib.pyplot as plt

PILOT_CHROMS_LOH = ["chr1", "chr8"]
BIN_SIZE = 1000  # 1 kb
MIN_READS = 5     # skip low-coverage bins to avoid noise


def get_chrom_length(bam, chrom):
    for ref, length in zip(bam.references, bam.lengths):
        if ref == chrom:
            return length
    return None


def main():
    setup_cjk_font()
    t0 = time.time()

    print(f"[A2_3] Loading LOH regions for {PILOT_CHROMS_LOH} ...",
          file=sys.stderr, flush=True)
    loh = load_loh_regions(chroms=PILOT_CHROMS_LOH)
    for c in PILOT_CHROMS_LOH:
        n = len(loh.get(c, []))
        total = sum(e - s for s, e in loh.get(c, []))
        print(f"  {c}: {n} LOH regions, total {total/1e6:.2f} Mb",
              file=sys.stderr, flush=True)
    log_timing("load_loh", t0, time.time())

    t1 = time.time()
    bam = pysam.AlignmentFile(BAM_PATH, "rb")

    rows = []
    for chrom in PILOT_CHROMS_LOH:
        L = get_chrom_length(bam, chrom)
        if L is None:
            print(f"[A2_3] WARN: {chrom} not found in BAM", file=sys.stderr)
            continue
        print(f"[A2_3] {chrom}: length {L}, bins ~{L//BIN_SIZE}",
              file=sys.stderr, flush=True)
        loh_intervals = loh.get(chrom, [])

        # Use count_coverage? simpler: iterate fetch per bin, count reads
        # Faster: iterate reads once, bin them
        # Reset per-bin counters
        n_bins = (L + BIN_SIZE - 1) // BIN_SIZE
        n_total = np.zeros(n_bins, dtype=np.int32)
        n_hp_missing = np.zeros(n_bins, dtype=np.int32)

        c_t0 = time.time()
        n_reads = 0
        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate or read.mapping_quality < 20:
                continue
            # Assign to bin by read mid position (alt: by start)
            mid = (read.reference_start + read.reference_end) // 2 if read.reference_end else read.reference_start
            bidx = mid // BIN_SIZE
            if bidx < 0 or bidx >= n_bins:
                continue
            n_total[bidx] += 1
            if not read.has_tag("HP"):
                n_hp_missing[bidx] += 1
            n_reads += 1
            if n_reads % 1_000_000 == 0:
                print(f"  {chrom} {n_reads:,} reads scanned ({time.time()-c_t0:.1f}s)",
                      file=sys.stderr, flush=True)
        print(f"  {chrom} total reads scanned: {n_reads:,} ({time.time()-c_t0:.1f}s)",
              file=sys.stderr, flush=True)

        # Build rows
        for bidx in range(n_bins):
            if n_total[bidx] < MIN_READS:
                continue
            start = bidx * BIN_SIZE
            end = min(start + BIN_SIZE, L)
            in_loh = interval_overlap(loh_intervals, start, end)
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "n_reads": int(n_total[bidx]),
                "n_hp_missing": int(n_hp_missing[bidx]),
                "hp_missing_rate": float(n_hp_missing[bidx]) / int(n_total[bidx]),
                "in_loh": int(in_loh),
            })

    bam.close()
    log_timing("bam_scan", t1, time.time())

    df = pd.DataFrame(rows)
    out_tsv = DATA_DIR / "A2_3_bin_hp_missing.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[A2_3] Wrote {out_tsv} ({len(df)} bins, "
          f"{df['in_loh'].sum()} LOH bins, {(1-df['in_loh']).sum()} non-LOH bins)",
          file=sys.stderr, flush=True)

    if df["in_loh"].nunique() < 2:
        print("[A2_3] FATAL: only one class in label", file=sys.stderr)
        return

    # ROC + AUC overall and per-chrom
    y = df["in_loh"].values
    score = df["hp_missing_rate"].values
    auc_all = roc_auc_score(y, score)
    fpr_all, tpr_all, _ = roc_curve(y, score)

    per_chrom = {}
    for c in PILOT_CHROMS_LOH:
        m = df["chrom"] == c
        if m.sum() < 10 or df.loc[m, "in_loh"].nunique() < 2:
            continue
        per_chrom[c] = {
            "n_bins": int(m.sum()),
            "auc": float(roc_auc_score(df.loc[m, "in_loh"], df.loc[m, "hp_missing_rate"])),
            "loh_frac": float(df.loc[m, "in_loh"].mean()),
            "median_hp_missing_loh": float(df.loc[m & (df["in_loh"]==1), "hp_missing_rate"].median()),
            "median_hp_missing_nonloh": float(df.loc[m & (df["in_loh"]==0), "hp_missing_rate"].median()),
        }

    stats = {
        "n_bins_total": int(len(df)),
        "auc_overall": float(auc_all),
        "frac_loh": float(df["in_loh"].mean()),
        "median_hp_missing_loh": float(df.loc[df["in_loh"]==1, "hp_missing_rate"].median()),
        "median_hp_missing_nonloh": float(df.loc[df["in_loh"]==0, "hp_missing_rate"].median()),
        "per_chrom": per_chrom,
    }
    print("[A2_3] STATS:", stats, file=sys.stderr, flush=True)

    # Save summary (flatten per_chrom)
    flat = {k: v for k, v in stats.items() if k != "per_chrom"}
    for c, sub in per_chrom.items():
        for k, v in sub.items():
            flat[f"{c}_{k}"] = v
    pd.Series(flat).to_csv(DATA_DIR / "A2_3_summary_stats.tsv", sep="\t")

    # Plot: ROC + violin
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    ax.plot(fpr_all, tpr_all, color="#264653", linewidth=2,
            label=f"Overall AUC = {auc_all:.3f}")
    for c, sub in per_chrom.items():
        m = df["chrom"] == c
        fpr_c, tpr_c, _ = roc_curve(df.loc[m, "in_loh"], df.loc[m, "hp_missing_rate"])
        ax.plot(fpr_c, tpr_c, alpha=0.7,
                label=f"{c} AUC = {sub['auc']:.3f} (n={sub['n_bins']})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(f"ROC: HP-missing rate vs SEQC2 LOH\n"
                 f"1 kb bins, ≥{MIN_READS} reads")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Violin
    ax2 = axes[1]
    data = [
        df.loc[df["in_loh"] == 0, "hp_missing_rate"].values,
        df.loc[df["in_loh"] == 1, "hp_missing_rate"].values,
    ]
    parts = ax2.violinplot(data, positions=[0, 1], showmedians=True, widths=0.7)
    for pc, color in zip(parts["bodies"], ["#2a9d8f", "#e76f51"]):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels([
        f"Non-LOH\n(n={(df['in_loh']==0).sum()})",
        f"LOH\n(n={(df['in_loh']==1).sum()})",
    ])
    ax2.set_ylabel("HP-missing rate")
    ax2.set_title(f"HP-missing rate by LOH status\n"
                  f"median: LOH={stats['median_hp_missing_loh']:.3f} vs "
                  f"non-LOH={stats['median_hp_missing_nonloh']:.3f}")
    ax2.grid(True, alpha=0.3, axis="y")

    fig.suptitle(
        "Analysis 3 — LOH HP-missing rate ROC vs SEQC2 truth (chr1+chr8 pilot)\n"
        "BAM: HCC1395 ClairS-TO ssrs v0.4 → LongPhase v1.7.3 haplotag --somaticMode",
        fontsize=10,
    )
    fig.tight_layout()
    out_png = FIG_DIR / "A2_3_loh_hp_missing_vs_seqc2.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[A2_3] Wrote {out_png}", file=sys.stderr, flush=True)

    log_timing("A2_3 total", t0, time.time())


if __name__ == "__main__":
    main()
