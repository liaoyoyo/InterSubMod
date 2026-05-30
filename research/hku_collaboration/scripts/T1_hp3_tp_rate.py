#!/usr/bin/env python3
"""
T1: HP3 read fraction + HP3 TP rate per chromosome (chr1 / chr8 / chr19).

Background
----------
Reviewer 3 Major #5 (LongPhase-S paper) asks: 在 LongPhase --somaticMode
輸出 tagged BAM 中，HP:Z:3 reads 佔總 reads 比例為何？HP3 中落在 SEQC2 truth
set 的 somatic TP 比例又為何？這份腳本針對 HCC1395 Tmode tagged BAM 的
三條既有 deliverable 已掃描染色體 (chr1 / chr8 / chr19) 量化兩個數字。

Inputs
------
- Tagged BAM: /big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam
- SEQC2 truth VCF: /big7_disk/liaoyoyo2001/InterSubMod/data/answer/SEQC2/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz

Outputs
-------
- TSV: findings_5_24/T1_HP3_TP_rate_per_chr.tsv
- PNG: figures/T1_HP3_TP_rate_per_chr.png

Read filter
-----------
- mapq >= 20
- exclude secondary / supplementary / duplicate / unmapped

HP tag notes
------------
LongPhase --somaticMode 寫出 HP 為 string tag (HP:Z:1 / HP:Z:2 /
HP:Z:1-1 / HP:Z:2-1 / HP:Z:3)。HP3 嚴格匹配字串 "3"。
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam

# ----- CJK font -----
from matplotlib import font_manager as fm

_CJK_FONT_PATH = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
if Path(_CJK_FONT_PATH).exists():
    fm.fontManager.addfont(_CJK_FONT_PATH)
    plt.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


BAM_PATH = "/big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam"
VCF_PATH = "/big7_disk/liaoyoyo2001/InterSubMod/data/answer/SEQC2/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"

CHROMS = ["chr1", "chr8", "chr19"]
EXPECTED_PASS = {"chr1": 3440, "chr8": 2605, "chr19": 874}  # sanity check

MIN_MAPQ = 20


def load_truth_positions(vcf_path: str, chroms: list[str]) -> dict[str, set[int]]:
    """Load SEQC2 PASS-only positions per chromosome.

    Returns dict {chr: set of 1-based VCF positions}.
    """
    truth: dict[str, set[int]] = {c: set() for c in chroms}
    vcf = pysam.VariantFile(vcf_path)
    for chrom in chroms:
        for rec in vcf.fetch(chrom):
            # PASS-only: this VCF tags every record with multiple filter labels
            # such as ['PASS', 'HighConf'] (PASS is just one of several tags).
            # Accept any record whose filter list either is empty or contains 'PASS'.
            filt = list(rec.filter)
            if filt and "PASS" not in filt:
                continue
            truth[chrom].add(rec.pos)  # 1-based
    vcf.close()
    return truth


def count_hp3_and_tp(
    bam_path: str,
    chrom: str,
    truth_positions: set[int],
) -> tuple[int, int, int]:
    """Per chromosome counts.

    Returns (total_reads, hp3_reads, hp3_tp_reads):
      total_reads  = filtered reads on chrom
      hp3_reads    = filtered reads with HP:Z:3
      hp3_tp_reads = HP3 reads whose alignment covers at least one truth pos
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    total = 0
    hp3 = 0
    hp3_tp = 0

    # Sort truth_positions to enable fast linear sweep via bisect if needed.
    # For per-read covered-position lookup, set intersection over read.get_reference_positions
    # is acceptable (truth size <= 3,440).
    truth_set = truth_positions  # already set

    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_duplicate:
            continue
        if read.mapping_quality < MIN_MAPQ:
            continue

        total += 1

        # HP tag (string)
        try:
            hp = read.get_tag("HP")
        except KeyError:
            hp = None

        if hp is None:
            continue
        # HP:Z:3 exact match (cast to str just in case some encoder writes int 3)
        if str(hp) != "3":
            continue

        hp3 += 1

        # Check truth coverage: any reference position covered by this read
        # falls into truth_set?
        # get_reference_positions returns 0-based; truth_set is 1-based VCF pos.
        # Convert by adding 1.
        ref_positions = read.get_reference_positions(full_length=False)
        if not ref_positions:
            continue
        # Quick range bail-out
        lo = ref_positions[0] + 1
        hi = ref_positions[-1] + 1
        # If no truth pos in [lo, hi] range, skip.
        # (Use any() over a set for correctness rather than range scan.)
        # Simple membership check:
        for p0 in ref_positions:
            if (p0 + 1) in truth_set:
                hp3_tp += 1
                break

    bam.close()
    return total, hp3, hp3_tp


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", default=BAM_PATH)
    p.add_argument("--vcf", default=VCF_PATH)
    p.add_argument(
        "--out-tsv",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration/findings_5_24/T1_HP3_TP_rate_per_chr.tsv",
    )
    p.add_argument(
        "--out-png",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration/figures/T1_HP3_TP_rate_per_chr.png",
    )
    args = p.parse_args(argv)

    print(f"[T1] BAM: {args.bam}", file=sys.stderr)
    print(f"[T1] VCF: {args.vcf}", file=sys.stderr)

    # 1. Load truth
    t0 = time.time()
    truth = load_truth_positions(args.vcf, CHROMS)
    for c in CHROMS:
        n = len(truth[c])
        exp = EXPECTED_PASS[c]
        ok = "OK" if n == exp else f"WARN(expected {exp})"
        print(f"[T1] truth {c}: {n} PASS positions  [{ok}]", file=sys.stderr)
    print(f"[T1] truth load done in {time.time()-t0:.1f}s", file=sys.stderr)

    # 2. Scan BAM per chrom
    rows = []
    for chrom in CHROMS:
        t_chr = time.time()
        total, hp3, hp3_tp = count_hp3_and_tp(args.bam, chrom, truth[chrom])
        elapsed = time.time() - t_chr
        hp3_frac = hp3 / total if total else 0.0
        hp3_tp_rate = hp3_tp / hp3 if hp3 else 0.0
        print(
            f"[T1] {chrom}: total={total:,} hp3={hp3:,} ({hp3_frac:.6%}) "
            f"hp3_tp={hp3_tp:,} ({hp3_tp_rate:.4%}) in {elapsed:.1f}s",
            file=sys.stderr,
        )
        rows.append(
            {
                "chromosome": chrom,
                "total_reads": total,
                "hp3_reads": hp3,
                "hp3_fraction": hp3_frac,
                "hp3_tp_reads": hp3_tp,
                "hp3_tp_rate": hp3_tp_rate,
            }
        )

    # 3. TSV
    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w") as fh:
        fh.write(
            "chromosome\ttotal_reads\thp3_reads\thp3_fraction\thp3_tp_reads\thp3_tp_rate\n"
        )
        for r in rows:
            fh.write(
                f"{r['chromosome']}\t{r['total_reads']}\t{r['hp3_reads']}\t"
                f"{r['hp3_fraction']:.8f}\t{r['hp3_tp_reads']}\t{r['hp3_tp_rate']:.8f}\n"
            )
    print(f"[T1] TSV written: {out_tsv}", file=sys.stderr)

    # 4. Plot
    fig, ax1 = plt.subplots(figsize=(8, 5))
    chroms = [r["chromosome"] for r in rows]
    fracs = [r["hp3_fraction"] for r in rows]
    rates = [r["hp3_tp_rate"] for r in rows]

    x = list(range(len(chroms)))
    bar_w = 0.35
    x_left = [i - bar_w / 2 for i in x]
    x_right = [i + bar_w / 2 for i in x]

    color_l = "#1f77b4"
    color_r = "#d62728"

    b1 = ax1.bar(x_left, fracs, bar_w, color=color_l, label="HP3 fraction (log)")
    ax1.set_yscale("log")
    ax1.set_ylabel("HP3 read fraction (log scale)", color=color_l)
    ax1.tick_params(axis="y", labelcolor=color_l)
    # Log-scale floor for visualization
    nonzero_fracs = [f for f in fracs if f > 0]
    if nonzero_fracs:
        ax1.set_ylim(min(nonzero_fracs) / 5, max(max(fracs), 1e-3))

    ax2 = ax1.twinx()
    b2 = ax2.bar(x_right, rates, bar_w, color=color_r, label="HP3 TP rate (linear)")
    ax2.set_ylabel("HP3 TP read rate (linear)", color=color_r)
    ax2.tick_params(axis="y", labelcolor=color_r)
    if rates:
        ax2.set_ylim(0, max(max(rates) * 1.25, 0.01))

    ax1.set_xticks(x)
    ax1.set_xticklabels(chroms)
    ax1.set_xlabel("Chromosome")
    ax1.set_title(
        "HCC1395 Tmode tagged BAM — HP3 fraction vs HP3 TP rate per chromosome"
    )

    # Value labels
    for rect, v in zip(b1, fracs):
        ax1.text(
            rect.get_x() + rect.get_width() / 2,
            v,
            f"{v:.2%}" if v >= 1e-4 else f"{v:.2e}",
            ha="center",
            va="bottom",
            fontsize=8,
            color=color_l,
        )
    for rect, v in zip(b2, rates):
        ax2.text(
            rect.get_x() + rect.get_width() / 2,
            v,
            f"{v:.2%}",
            ha="center",
            va="bottom",
            fontsize=8,
            color=color_r,
        )

    lines = [b1, b2]
    labels = ["HP3 fraction (log)", "HP3 TP rate (linear)"]
    ax1.legend(lines, labels, loc="upper left", fontsize=9)

    fig.tight_layout()
    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"[T1] PNG written: {out_png}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
