#!/usr/bin/env python3
"""Analysis 1: PS set TP somatic位點 數量分佈.

對每個 PS set，count 跨該 PS 的 TP somatic 位點數。
目標：證明 PS set 內有足夠 TP 位點可做 read-level 分群。

Output:
  - data/A2_1_ps_tp_counts.tsv
  - figures/A2_1_ps_set_tp_count_dist.png
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).parent))
from _common import (
    BAM_PATH, DATA_DIR, FIG_DIR, PILOT_CHROMS,
    load_tp_positions, setup_cjk_font, log_timing,
)
import matplotlib.pyplot as plt


def main():
    setup_cjk_font()
    t0 = time.time()
    print(f"[A2_1] Loading TP positions for {PILOT_CHROMS} ...", file=sys.stderr, flush=True)
    tp = load_tp_positions(chroms=PILOT_CHROMS)
    for c, plist in tp.items():
        print(f"  {c}: {len(plist)} TP", file=sys.stderr, flush=True)
    log_timing("load_tp", t0, time.time())

    # For each PS: track set of TP positions covered + set of unique reads
    # ps_id is per-chrom unique (PS is just the first variant pos in set; combine chrom+ps)
    ps_tp: dict[tuple[str, int], set[int]] = defaultdict(set)
    ps_reads: dict[tuple[str, int], set[str]] = defaultdict(set)
    ps_chrom: dict[tuple[str, int], str] = {}

    t1 = time.time()
    bam = pysam.AlignmentFile(BAM_PATH, "rb")

    for chrom in PILOT_CHROMS:
        if chrom not in tp:
            continue
        positions = tp[chrom]
        print(f"[A2_1] Processing {chrom}: {len(positions)} TP positions ...", file=sys.stderr, flush=True)
        chrom_t0 = time.time()
        # Iterate TP positions; for each, find all reads spanning it and their PS tags.
        for i, pos in enumerate(positions):
            if i % 2000 == 0 and i > 0:
                print(f"  {chrom} {i}/{len(positions)} ({time.time()-chrom_t0:.1f}s)",
                      file=sys.stderr, flush=True)
            # pysam uses 0-based; VCF pos is 1-based
            zpos = pos - 1
            try:
                for read in bam.fetch(chrom, zpos, zpos + 1):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    if read.is_duplicate or read.mapping_quality < 20:
                        continue
                    if not read.has_tag("PS"):
                        continue
                    ps = read.get_tag("PS")
                    key = (chrom, int(ps))
                    ps_tp[key].add(pos)
                    ps_reads[key].add(read.query_name)
                    ps_chrom[key] = chrom
            except (ValueError, OSError) as e:
                print(f"  fetch error at {chrom}:{pos}: {e}", file=sys.stderr)
                continue
        print(f"  {chrom} done in {time.time()-chrom_t0:.1f}s",
              file=sys.stderr, flush=True)

    bam.close()
    log_timing("bam_scan", t1, time.time())

    # Build dataframe
    rows = []
    for key, tps in ps_tp.items():
        chrom, ps = key
        rows.append({
            "chrom": chrom,
            "ps_id": ps,
            "n_tp": len(tps),
            "n_reads": len(ps_reads[key]),
        })
    df = pd.DataFrame(rows).sort_values(["chrom", "ps_id"]).reset_index(drop=True)
    out_tsv = DATA_DIR / "A2_1_ps_tp_counts.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[A2_1] Wrote {out_tsv} ({len(df)} PS sets)", file=sys.stderr, flush=True)

    # Stats
    n_tp = df["n_tp"].values
    stats = {
        "n_ps_sets": int(len(df)),
        "median_tp_per_ps": float(np.median(n_tp)),
        "mean_tp_per_ps": float(np.mean(n_tp)),
        "p25": float(np.percentile(n_tp, 25)),
        "p75": float(np.percentile(n_tp, 75)),
        "max_tp_per_ps": int(np.max(n_tp)) if len(n_tp) else 0,
        "n_ps_with_ge2_tp": int((n_tp >= 2).sum()),
        "frac_ps_with_ge2_tp": float((n_tp >= 2).mean()) if len(n_tp) else 0,
    }
    print("[A2_1] STATS:", stats, file=sys.stderr, flush=True)
    pd.Series(stats).to_csv(DATA_DIR / "A2_1_summary_stats.tsv", sep="\t")

    # Plot: 2x2 dual-panel layout (all linear axes per 2026-05-23 review point #2)
    # Row 1: full-range histogram + full-range CDF
    # Row 2: zoom-in [0, 20] histogram + zoom-in [0, 20] CDF (focus on main peak)
    max_tp = int(np.max(n_tp)) if len(n_tp) else 10
    full_clip = max(max_tp, 10)  # full range up to actual max
    zoom_clip = 20  # zoom-in to main peak

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # --- Row 1, Col 1: full-range histogram (linear) ---
    ax = axes[0, 0]
    bins_full = np.arange(0, full_clip + 2) - 0.5
    ax.hist(n_tp, bins=bins_full, color="#3a7ca5",
            edgecolor="white", linewidth=0.4)
    ax.set_xlabel(f"PS set 內 TP somatic 位點數 (full range 0-{full_clip})")
    ax.set_ylabel("PS set 數 (linear)")
    ax.set_title(
        f"[Full range] PS set TP 位點分佈\n"
        f"N={stats['n_ps_sets']} PS sets, median={stats['median_tp_per_ps']:.0f}, "
        f"IQR=[{stats['p25']:.0f},{stats['p75']:.0f}], max={stats['max_tp_per_ps']}"
    )
    ax.axvline(stats["median_tp_per_ps"], color="red", linestyle="--",
               label=f"median={stats['median_tp_per_ps']:.0f}")
    ax.axvline(2, color="orange", linestyle=":",
               label=f"≥2 TP: {stats['frac_ps_with_ge2_tp']*100:.1f}% PS sets")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Row 1, Col 2: full-range CDF (linear) ---
    ax2 = axes[0, 1]
    sorted_n = np.sort(n_tp)
    cdf = np.arange(1, len(sorted_n) + 1) / len(sorted_n)
    ax2.plot(sorted_n, cdf, color="#2a9d8f", linewidth=1.8)
    ax2.set_xlabel(f"PS 內 TP 位點數 (linear, full 0-{full_clip})")
    ax2.set_ylabel("Cumulative fraction of PS sets")
    ax2.set_title("[Full range] PS set TP count CDF")
    ax2.axhline(0.5, color="red", linestyle="--", alpha=0.6, label="median")
    ax2.axvline(stats["median_tp_per_ps"], color="red", linestyle="--", alpha=0.6)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # --- Row 2, Col 1: zoom-in [0, 20] histogram (linear) ---
    ax3 = axes[1, 0]
    bins_zoom = np.arange(0, zoom_clip + 2) - 0.5
    n_tp_zoom = n_tp[n_tp <= zoom_clip]
    ax3.hist(n_tp_zoom, bins=bins_zoom, color="#3a7ca5",
             edgecolor="white", linewidth=0.4)
    ax3.set_xlim(-0.5, zoom_clip + 0.5)
    ax3.set_xlabel(f"PS set 內 TP somatic 位點數 (zoom [0, {zoom_clip}])")
    ax3.set_ylabel("PS set 數 (linear)")
    ax3.set_title(
        f"[Zoom main peak] PS set TP 位點分佈 0-{zoom_clip}\n"
        f"{(n_tp <= zoom_clip).sum()} / {len(n_tp)} PS sets in this range "
        f"({(n_tp <= zoom_clip).mean()*100:.1f}%)"
    )
    ax3.axvline(stats["median_tp_per_ps"], color="red", linestyle="--",
                label=f"median={stats['median_tp_per_ps']:.0f}")
    ax3.axvline(2, color="orange", linestyle=":",
                label=f"≥2 TP cut")
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # --- Row 2, Col 2: zoom-in [0, 20] CDF (linear) ---
    ax4 = axes[1, 1]
    ax4.plot(sorted_n, cdf, color="#2a9d8f", linewidth=1.8)
    ax4.set_xlim(0, zoom_clip)
    ax4.set_xlabel(f"PS 內 TP 位點數 (linear, zoom [0, {zoom_clip}])")
    ax4.set_ylabel("Cumulative fraction of PS sets")
    ax4.set_title(f"[Zoom main peak] PS set TP count CDF 0-{zoom_clip}")
    ax4.axhline(0.5, color="red", linestyle="--", alpha=0.6, label="median")
    ax4.axvline(stats["median_tp_per_ps"], color="red", linestyle="--", alpha=0.6)
    ax4.axvline(2, color="orange", linestyle=":", alpha=0.6, label="≥2 TP cut")
    ax4.grid(True, alpha=0.3)
    ax4.legend()

    fig.suptitle(
        "Analysis 1 — PS set 內 TP somatic 位點分佈 (linear axes, dual-panel)\n"
        "BAM: HCC1395 ClairS-TO ssrs v0.4 → LongPhase v1.7.3 haplotag --somaticMode",
        fontsize=10,
    )
    fig.tight_layout()
    out_png = FIG_DIR / "A2_1_ps_set_tp_count_dist.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[A2_1] Wrote {out_png}", file=sys.stderr, flush=True)

    log_timing("A2_1 total", t0, time.time())


if __name__ == "__main__":
    main()
