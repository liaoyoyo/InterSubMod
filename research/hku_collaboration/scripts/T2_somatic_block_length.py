#!/usr/bin/env python3
"""
T2 — Somatic phase-block length characterization (LongPhase-S Reviewer 1 Minor #3).

Scan an LongPhase-tagged BAM, group reads by PS tag, compute per-PS-block
extent (max_end - min_start) and read count, then summarize per chromosome
and globally. Output: per-block TSV, summary TSV, two-panel PNG.

Usage:
    python3 T2_somatic_block_length.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, rcParams
import numpy as np
import pandas as pd
import pysam

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BAM_PATH = Path(
    "/big7_disk/liaoyoyo2001/data/bam/"
    "HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam"
)
OUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration/findings_5_24"
)
FIG_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration/figures"
)
PER_BLOCK_TSV = OUT_DIR / "T2_somatic_block_per_block.tsv"
SUMMARY_TSV = OUT_DIR / "T2_somatic_block_summary.tsv"
FIG_PATH = FIG_DIR / "T2_somatic_block_length_dist.png"

CHROMS = ["chr1", "chr8", "chr19"]
MAPQ_MIN = 20
# Median TP pair distance reference (verified upstream artifact)
TP_PAIR_REF_BP = 17_400

# ---------------------------------------------------------------------------
# CJK font setup (Droid Sans Fallback verified available on host)
# ---------------------------------------------------------------------------
_CJK_FONT_PATH = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
if Path(_CJK_FONT_PATH).exists():
    font_manager.fontManager.addfont(_CJK_FONT_PATH)
rcParams["font.family"] = ["DejaVu Sans", "Droid Sans Fallback"]
rcParams["axes.unicode_minus"] = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def n50(lengths: np.ndarray) -> float:
    """Compute N50 over a vector of block lengths."""
    if lengths.size == 0:
        return float("nan")
    sorted_desc = np.sort(lengths)[::-1]
    cum = np.cumsum(sorted_desc)
    half = cum[-1] / 2.0
    idx = np.searchsorted(cum, half)
    return float(sorted_desc[idx])


def summarize(label: str, df_blocks: pd.DataFrame) -> dict:
    """Compute summary dict over a per-block dataframe slice."""
    if df_blocks.empty:
        return {
            "chromosome": label, "n_blocks": 0,
            "median_len": float("nan"), "q1_len": float("nan"),
            "q3_len": float("nan"), "n50": float("nan"),
            "max_len": float("nan"), "min_len": float("nan"),
            "mean_len": float("nan"),
            "median_reads_per_ps": float("nan"),
        }
    lens = df_blocks["block_length"].to_numpy()
    reads = df_blocks["n_reads"].to_numpy()
    return {
        "chromosome": label,
        "n_blocks": int(lens.size),
        "median_len": float(np.median(lens)),
        "q1_len": float(np.percentile(lens, 25)),
        "q3_len": float(np.percentile(lens, 75)),
        "n50": n50(lens),
        "max_len": int(np.max(lens)),
        "min_len": int(np.min(lens)),
        "mean_len": float(np.mean(lens)),
        "median_reads_per_ps": float(np.median(reads)),
    }


def scan_chrom(bam: pysam.AlignmentFile, chrom: str) -> pd.DataFrame:
    """
    Iterate reads on `chrom`, group by PS tag.

    Returns dataframe with columns:
        chromosome, PS_id, block_start, block_end, block_length, n_reads
    """
    # PS_id -> [min_start, max_end, n_reads]
    blocks: dict[str, list[int]] = {}
    n_total = 0
    n_kept = 0
    n_with_ps = 0
    t0 = time.time()

    for read in bam.fetch(chrom, multiple_iterators=False):
        n_total += 1
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_duplicate or read.is_unmapped:
            continue
        if read.mapping_quality < MAPQ_MIN:
            continue
        n_kept += 1
        if not read.has_tag("PS"):
            continue
        n_with_ps += 1
        # PS may be int (PS:i:) or str (PS:Z:); normalize to str.
        ps_id = str(read.get_tag("PS"))
        start = int(read.reference_start)
        end = int(read.reference_end) if read.reference_end is not None else start + 1
        rec = blocks.get(ps_id)
        if rec is None:
            blocks[ps_id] = [start, end, 1]
        else:
            if start < rec[0]:
                rec[0] = start
            if end > rec[1]:
                rec[1] = end
            rec[2] += 1

    dt = time.time() - t0
    print(
        f"[{chrom}] reads_total={n_total:,} kept_mapq{MAPQ_MIN}+_primary="
        f"{n_kept:,} with_PS={n_with_ps:,} blocks={len(blocks):,} "
        f"({dt:.1f}s)",
        flush=True,
    )

    if not blocks:
        return pd.DataFrame(
            columns=[
                "chromosome", "PS_id", "block_start", "block_end",
                "block_length", "n_reads",
            ]
        )

    rows = []
    for ps_id, (mn, mx, nr) in blocks.items():
        rows.append(
            (chrom, ps_id, mn, mx, mx - mn, nr)
        )
    df = pd.DataFrame(
        rows,
        columns=[
            "chromosome", "PS_id", "block_start", "block_end",
            "block_length", "n_reads",
        ],
    )
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"BAM:    {BAM_PATH}", flush=True)
    print(f"Chroms: {CHROMS}", flush=True)
    print(f"MAPQ:   >= {MAPQ_MIN}", flush=True)

    bam = pysam.AlignmentFile(str(BAM_PATH), "rb")
    per_chrom = []
    for ch in CHROMS:
        df_ch = scan_chrom(bam, ch)
        per_chrom.append(df_ch)
    bam.close()

    df_blocks = pd.concat(per_chrom, ignore_index=True)
    df_blocks = df_blocks.sort_values(
        ["chromosome", "block_start"], ignore_index=True
    )
    df_blocks.to_csv(PER_BLOCK_TSV, sep="\t", index=False)
    print(f"per-block TSV -> {PER_BLOCK_TSV} ({len(df_blocks):,} rows)", flush=True)

    # Summaries: per chromosome + ALL
    summaries = []
    for ch in CHROMS:
        summaries.append(summarize(ch, df_blocks[df_blocks["chromosome"] == ch]))
    summaries.append(summarize("ALL", df_blocks))

    # Histogram counts for ALL (1k / 10k / 100k / 1M / 10M edges)
    bin_edges = [0, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, np.inf]
    bin_labels = ["<1kb", "1-10kb", "10-100kb", "100kb-1Mb", "1-10Mb", ">10Mb"]
    if not df_blocks.empty:
        hist_counts, _ = np.histogram(
            df_blocks["block_length"].to_numpy(), bins=bin_edges
        )
        hist_summary = {f"hist_{lab}": int(c) for lab, c in zip(bin_labels, hist_counts)}
    else:
        hist_summary = {f"hist_{lab}": 0 for lab in bin_labels}

    df_summary = pd.DataFrame(summaries)
    # Append histogram columns only on the ALL row (others NaN)
    for lab in bin_labels:
        col = f"hist_{lab}"
        df_summary[col] = float("nan")
    if not df_blocks.empty:
        for lab in bin_labels:
            df_summary.loc[df_summary["chromosome"] == "ALL", f"hist_{lab}"] = (
                hist_summary[f"hist_{lab}"]
            )

    df_summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
    print(f"summary TSV   -> {SUMMARY_TSV}", flush=True)

    # ---------------- Figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax_h, ax_b = axes

    all_lens = df_blocks["block_length"].to_numpy()
    if all_lens.size:
        # Log-spaced bins from 100 bp to 100 Mb (covers all observed)
        log_bins = np.logspace(2, 8, 60)
        ax_h.hist(all_lens, bins=log_bins, color="steelblue", edgecolor="white")
        ax_h.set_xscale("log")
        med_all = np.median(all_lens)
        n50_all = n50(all_lens)
        ax_h.axvline(med_all, color="crimson", linestyle="--", lw=1.5,
                     label=f"median = {med_all:,.0f} bp")
        ax_h.axvline(n50_all, color="darkorange", linestyle="--", lw=1.5,
                     label=f"N50 = {n50_all:,.0f} bp")
        ax_h.axvline(TP_PAIR_REF_BP, color="green", linestyle=":", lw=1.5,
                     label=f"TP pair median ref = {TP_PAIR_REF_BP:,} bp")
        # Mark decade ticks for context
        for x in [1e3, 1e4, 1e5, 1e6, 1e7]:
            ax_h.axvline(x, color="lightgray", linestyle="-", lw=0.5, zorder=0)
        ax_h.legend(fontsize=9, loc="upper right")
    ax_h.set_xlabel("PS block length (bp, log scale)")
    ax_h.set_ylabel("Number of PS blocks")
    ax_h.set_title("Somatic phase-block length distribution (ALL chroms)")

    # Box plot per chromosome (log y)
    box_data = [
        df_blocks[df_blocks["chromosome"] == ch]["block_length"].to_numpy()
        for ch in CHROMS
    ]
    bp = ax_b.boxplot(
        box_data, labels=CHROMS, showfliers=False,
        patch_artist=True,
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("steelblue")
        patch.set_alpha(0.6)
    ax_b.set_yscale("log")
    ax_b.set_ylabel("PS block length (bp, log scale)")
    ax_b.set_title("Per-chromosome block length (box, fliers hidden)")
    ax_b.axhline(TP_PAIR_REF_BP, color="green", linestyle=":", lw=1.2,
                 label=f"TP pair median ref = {TP_PAIR_REF_BP:,} bp")
    ax_b.legend(fontsize=9, loc="lower right")

    fig.suptitle(
        "HCC1395 Tagged BAM — Somatic Haplotype Block Length "
        "(LongPhase-S, MAPQ>=20, primary)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_PATH, dpi=150)
    plt.close(fig)
    print(f"figure        -> {FIG_PATH}", flush=True)

    # Console summary print
    print("\n=== Summary ===")
    cols_to_show = [
        "chromosome", "n_blocks", "median_len", "q1_len", "q3_len",
        "n50", "max_len", "median_reads_per_ps",
    ]
    with pd.option_context(
        "display.float_format", lambda v: f"{v:,.0f}" if pd.notna(v) else "NA"
    ):
        print(df_summary[cols_to_show].to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
