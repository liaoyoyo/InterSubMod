#!/usr/bin/env python3
"""Analysis 2: 相鄰兩 TP 位點 common reads (distance vs common_reads 散點).

對每個 PS set 內相鄰 TP somatic pair (i, j) — by chrom+pos sort：
  distance = pos_j - pos_i
  common_reads = 跨 pos_i 與 pos_j 兩個位點的 read 數
畫散點 + ±16 bp caller window 紅線參考.

Output:
  - data/A2_2_pair_common_reads.tsv
  - figures/A2_2_common_reads_xy.png
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

# ClairS-TO ssrs caller window from A1 report Q1: ±16 bp = 33 bp total
CALLER_WIN = 16


def get_reads_at_pos(bam, chrom, pos1based):
    """Return (set of read names with PS tag, dict[read_name]->ps) covering pos (1-based)."""
    zpos = pos1based - 1
    reads = {}
    try:
        for read in bam.fetch(chrom, zpos, zpos + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate or read.mapping_quality < 20:
                continue
            if not read.has_tag("PS"):
                continue
            reads[read.query_name] = int(read.get_tag("PS"))
    except (ValueError, OSError):
        pass
    return reads


def main():
    setup_cjk_font()
    t0 = time.time()

    # Load TP positions
    print(f"[A2_2] Loading TP positions for {PILOT_CHROMS} ...", file=sys.stderr, flush=True)
    tp = load_tp_positions(chroms=PILOT_CHROMS)
    for c, plist in tp.items():
        print(f"  {c}: {len(plist)} TP", file=sys.stderr, flush=True)
    log_timing("load_tp", t0, time.time())

    # First pass: build mapping of TP_pos -> dominant PS id (the PS shared by most reads).
    # Then for adjacent pairs sharing same PS, compute common reads.
    t1 = time.time()
    bam = pysam.AlignmentFile(BAM_PATH, "rb")

    rows = []
    for chrom in PILOT_CHROMS:
        if chrom not in tp:
            continue
        positions = tp[chrom]
        if len(positions) < 2:
            continue
        print(f"[A2_2] {chrom}: scanning {len(positions)} TP positions, "
              f"computing {len(positions)-1} adjacent pairs",
              file=sys.stderr, flush=True)

        # Cache reads@pos sliding window: only need previous + current
        prev_reads = None
        prev_pos = None
        c_t0 = time.time()
        for i, pos in enumerate(positions):
            if i % 2000 == 0 and i > 0:
                print(f"  {chrom} {i}/{len(positions)} ({time.time()-c_t0:.1f}s, "
                      f"{len(rows)} pairs so far)", file=sys.stderr, flush=True)
            cur_reads = get_reads_at_pos(bam, chrom, pos)
            if prev_reads is not None and prev_pos is not None:
                # Find PS sets that have reads at both pos_i (prev) and pos_j (cur)
                # Per-PS common reads
                # Group prev_reads by PS
                prev_by_ps = defaultdict(set)
                for rname, ps in prev_reads.items():
                    prev_by_ps[ps].add(rname)
                cur_by_ps = defaultdict(set)
                for rname, ps in cur_reads.items():
                    cur_by_ps[ps].add(rname)
                shared_ps = set(prev_by_ps) & set(cur_by_ps)
                for ps in shared_ps:
                    common = prev_by_ps[ps] & cur_by_ps[ps]
                    if len(common) == 0:
                        continue
                    rows.append({
                        "chrom": chrom,
                        "pos_i": prev_pos,
                        "pos_j": pos,
                        "distance": pos - prev_pos,
                        "common_reads": len(common),
                        "ps_id": ps,
                        "n_reads_i": len(prev_by_ps[ps]),
                        "n_reads_j": len(cur_by_ps[ps]),
                    })
            prev_reads = cur_reads
            prev_pos = pos
        print(f"  {chrom} done in {time.time()-c_t0:.1f}s, "
              f"total pairs: {len(rows)}", file=sys.stderr, flush=True)

    bam.close()
    log_timing("bam_scan", t1, time.time())

    df = pd.DataFrame(rows)
    out_tsv = DATA_DIR / "A2_2_pair_common_reads.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[A2_2] Wrote {out_tsv} ({len(df)} pairs)", file=sys.stderr, flush=True)

    if len(df) == 0:
        print("[A2_2] WARN: no pairs", file=sys.stderr)
        return

    # Stats
    in_win = df["distance"] <= CALLER_WIN
    stats = {
        "n_pairs_total": int(len(df)),
        "median_distance": float(df["distance"].median()),
        "median_common_reads": float(df["common_reads"].median()),
        "mean_common_reads": float(df["common_reads"].mean()),
        "n_pairs_in_caller_win": int(in_win.sum()),
        "frac_pairs_in_caller_win": float(in_win.mean()),
        "median_common_reads_in_win": float(df.loc[in_win, "common_reads"].median()) if in_win.any() else float("nan"),
        "median_common_reads_out_win": float(df.loc[~in_win, "common_reads"].median()) if (~in_win).any() else float("nan"),
    }
    print("[A2_2] STATS:", stats, file=sys.stderr, flush=True)
    pd.Series(stats).to_csv(DATA_DIR / "A2_2_summary_stats.tsv", sep="\t")

    # Plot: 2x2 dual-panel layout (all linear axes per 2026-05-23 review point #2)
    # Row 1: full-range hexbin (distance ≤ 500) + full-range boxplot (linear Y)
    # Row 2: zoom-in hexbin (distance ≤ 100 bp, common ≤ 100) + caller-win focus boxplot
    full_dist_clip = 500
    zoom_dist_clip = 100   # focus around ±16 bp caller window
    zoom_common_clip = 100  # focus low-common region (most pairs)

    fig, axes = plt.subplots(2, 2, figsize=(15, 11))

    # --- Row 1, Col 1: full-range hexbin (linear) ---
    ax = axes[0, 0]
    mask_full = df["distance"] <= full_dist_clip
    if mask_full.sum() > 0:
        x = df.loc[mask_full, "distance"].values
        y = df.loc[mask_full, "common_reads"].values
        hb = ax.hexbin(x, y, gridsize=50, cmap="viridis", mincnt=1)  # linear (no yscale)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("pair count")
    ax.axvline(CALLER_WIN, color="red", linestyle="--", linewidth=1.5,
               label=f"ClairS-TO ssrs caller win = ±{CALLER_WIN} bp")
    ax.set_xlabel(f"相鄰 TP pair 距離 (bp, full range ≤ {full_dist_clip})")
    ax.set_ylabel("Common reads 數 (linear)")
    ax.set_title(
        f"[Full range] distance × common reads (linear)\n"
        f"N={stats['n_pairs_total']} pairs, "
        f"median common={stats['median_common_reads']:.0f}, "
        f"in-window={stats['n_pairs_in_caller_win']} ({stats['frac_pairs_in_caller_win']*100:.1f}%)"
    )
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    # --- Row 1, Col 2: full-range boxplot (linear Y) ---
    ax2 = axes[0, 1]
    data_box = [
        df.loc[in_win, "common_reads"].values,
        df.loc[(~in_win) & (df["distance"] <= full_dist_clip), "common_reads"].values,
        df.loc[df["distance"] > full_dist_clip, "common_reads"].values,
    ]
    labels = [
        f"≤±{CALLER_WIN} bp\n(n={in_win.sum()})",
        f"({CALLER_WIN},{full_dist_clip}] bp\n(n={((~in_win) & (df['distance']<=full_dist_clip)).sum()})",
        f">{full_dist_clip} bp\n(n={(df['distance']>full_dist_clip).sum()})",
    ]
    bp = ax2.boxplot(data_box, labels=labels, showfliers=False, patch_artist=True,
                     medianprops=dict(color="red", linewidth=2))
    colors = ["#e76f51", "#f4a261", "#2a9d8f"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax2.set_ylabel("Common reads (linear)")
    ax2.set_title("[Full range] Common reads by distance bucket")
    ax2.grid(True, alpha=0.3, axis="y")

    # --- Row 2, Col 1: zoom-in hexbin (distance ≤ 100 bp, common ≤ 100) — caller-win focus ---
    ax3 = axes[1, 0]
    mask_zoom = (df["distance"] <= zoom_dist_clip) & (df["common_reads"] <= zoom_common_clip)
    if mask_zoom.sum() > 0:
        x = df.loc[mask_zoom, "distance"].values
        y = df.loc[mask_zoom, "common_reads"].values
        hb = ax3.hexbin(x, y, gridsize=40, cmap="viridis", mincnt=1)  # linear
        cb = fig.colorbar(hb, ax=ax3)
        cb.set_label("pair count")
    ax3.axvline(CALLER_WIN, color="red", linestyle="--", linewidth=1.5,
                label=f"caller win = ±{CALLER_WIN} bp")
    ax3.set_xlabel(f"相鄰 TP pair 距離 (bp, zoom ≤ {zoom_dist_clip})")
    ax3.set_ylabel(f"Common reads (linear, zoom ≤ {zoom_common_clip})")
    ax3.set_xlim(0, zoom_dist_clip)
    ax3.set_ylim(0, zoom_common_clip)
    ax3.set_title(
        f"[Zoom caller-win focus] distance ≤ {zoom_dist_clip} bp × common ≤ {zoom_common_clip}\n"
        f"{mask_zoom.sum()} / {len(df)} pairs in this range "
        f"({mask_zoom.mean()*100:.1f}%)"
    )
    ax3.legend(loc="upper right")
    ax3.grid(True, alpha=0.3)

    # --- Row 2, Col 2: zoom boxplot focused on caller window vs just-outside ---
    ax4 = axes[1, 1]
    # Finer-grained distance buckets in caller-window focus zone
    bucket_le_win = df["distance"] <= CALLER_WIN
    bucket_win_to_50 = (df["distance"] > CALLER_WIN) & (df["distance"] <= 50)
    bucket_50_to_100 = (df["distance"] > 50) & (df["distance"] <= 100)
    bucket_100_to_500 = (df["distance"] > 100) & (df["distance"] <= 500)
    data_zoom_box = [
        df.loc[bucket_le_win, "common_reads"].values,
        df.loc[bucket_win_to_50, "common_reads"].values,
        df.loc[bucket_50_to_100, "common_reads"].values,
        df.loc[bucket_100_to_500, "common_reads"].values,
    ]
    labels_zoom = [
        f"≤{CALLER_WIN} bp\n(n={bucket_le_win.sum()})",
        f"({CALLER_WIN},50] bp\n(n={bucket_win_to_50.sum()})",
        f"(50,100] bp\n(n={bucket_50_to_100.sum()})",
        f"(100,500] bp\n(n={bucket_100_to_500.sum()})",
    ]
    bp2 = ax4.boxplot(data_zoom_box, labels=labels_zoom, showfliers=False, patch_artist=True,
                      medianprops=dict(color="red", linewidth=2))
    colors2 = ["#e76f51", "#f4a261", "#e9c46a", "#2a9d8f"]
    for patch, color in zip(bp2["boxes"], colors2):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax4.set_ylabel("Common reads (linear)")
    ax4.set_ylim(0, zoom_common_clip)
    ax4.set_title(f"[Zoom common ≤ {zoom_common_clip}] Common reads, finer distance buckets")
    ax4.grid(True, alpha=0.3, axis="y")

    fig.suptitle(
        "Analysis 2 — 相鄰 TP somatic pair 的 common reads 分佈 (linear axes, dual-panel)\n"
        "BAM: HCC1395 ClairS-TO ssrs v0.4 → LongPhase v1.7.3 haplotag --somaticMode  |  "
        f"caller win ±{CALLER_WIN} bp from A1 Q1",
        fontsize=10,
    )
    fig.tight_layout()
    out_png = FIG_DIR / "A2_2_common_reads_xy.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[A2_2] Wrote {out_png}", file=sys.stderr, flush=True)

    log_timing("A2_2 total", t0, time.time())


if __name__ == "__main__":
    main()
