#!/usr/bin/env python3
"""
T10: HP3 read fraction + HP3 TP rate per chromosome — full 24-chr extension.

Background
----------
Extends T1 (chr1 / chr8 / chr19 pilot) to all 24 chromosomes (chr1..chr22 + chrX;
chrY skipped because HCC1395 is a female-derived cell line and SEQC2 truth set
covers chr1-22 + chrX only). Reviewer 3 Major #5 (LongPhase-S paper) asks for
genome-wide quantification of HP:Z:3 fraction in tagged BAM and how many of
those HP3 reads cover SEQC2 PASS somatic positions (HP3 TP rate).

Bug fix vs T1
-------------
- SEQC2 VCF has multi-label FILTER (e.g. ('PASS', 'HighConf')).
  Correct check: 'PASS' in list(rec.filter), NOT filt == ['PASS'].
  (Template T1_hp3_tp_rate.py already applies this fix.)

Inputs
------
- Tagged BAM: /big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam
- SEQC2 truth VCF: /big7_disk/liaoyoyo2001/InterSubMod/data/answer/SEQC2/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz

Outputs
-------
- TSV: findings_5_24/T10_HP3_TP_rate_24chr.tsv
- PNG: figures/T10_HP3_TP_rate_24chr.png

Read filter
-----------
- mapq >= 20
- exclude secondary / supplementary / duplicate / unmapped (primary only)

HP tag notes
------------
LongPhase --somaticMode 寫出 HP 為 string tag (HP:Z:1 / HP:Z:2 /
HP:Z:1-1 / HP:Z:2-1 / HP:Z:3). HP3 嚴格匹配字串 "3".
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

# 24 chrs: chr1..22 autosomes + chrX. chrY skipped (female sample).
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

# Sanity baseline from T1 (3-chr pilot)
T1_BASELINE = {
    "chr1": {"total": 2236562, "hp3": 6419, "hp3_tp": 5658, "truth_pass": 3440},
    "chr8": {"total": 1653223, "hp3": 5539, "hp3_tp": 4454, "truth_pass": 2605},
    "chr19": {"total": 535588, "hp3": 1264, "hp3_tp": 1184, "truth_pass": 874},
}

MIN_MAPQ = 20

# Outlier / annotation thresholds (per task spec)
HP3_FRAC_OUTLIER_THRESHOLD = 0.05  # > 5% HP3 fraction = red highlight
CHR8_LOH_HOTSPOT = "chr8"  # light-blue frame annotation


def load_truth_positions(vcf_path: str, chroms: list[str]) -> dict[str, set[int]]:
    """Load SEQC2 PASS positions per chromosome.

    Bug fix: SEQC2 VCF tags each record with multiple FILTER labels
    (e.g. ['PASS', 'HighConf']). Accept records where 'PASS' is in the
    filter list (NOT filt == ['PASS']).
    """
    truth: dict[str, set[int]] = {c: set() for c in chroms}
    vcf = pysam.VariantFile(vcf_path)
    for chrom in chroms:
        try:
            for rec in vcf.fetch(chrom):
                filt = list(rec.filter)
                # Empty filter or PASS-containing accepted
                if filt and "PASS" not in filt:
                    continue
                truth[chrom].add(rec.pos)  # 1-based
        except ValueError as e:
            print(f"[T10] WARN: fetch failed for {chrom}: {e}", file=sys.stderr)
    vcf.close()
    return truth


def count_hp3_and_tp(
    bam_path: str,
    chrom: str,
    truth_positions: set[int],
) -> tuple[int, int, int]:
    """Per chromosome counts: (total_reads, hp3_reads, hp3_tp_reads)."""
    bam = pysam.AlignmentFile(bam_path, "rb")

    total = 0
    hp3 = 0
    hp3_tp = 0

    truth_set = truth_positions

    try:
        read_iter = bam.fetch(chrom)
    except ValueError as e:
        print(f"[T10] WARN: BAM fetch failed for {chrom}: {e}", file=sys.stderr)
        bam.close()
        return 0, 0, 0

    for read in read_iter:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_duplicate:
            continue
        if read.mapping_quality < MIN_MAPQ:
            continue

        total += 1

        try:
            hp = read.get_tag("HP")
        except KeyError:
            hp = None

        if hp is None:
            continue
        if str(hp) != "3":
            continue

        hp3 += 1

        # Skip truth coverage check if no truth on this chrom
        if not truth_set:
            continue

        ref_positions = read.get_reference_positions(full_length=False)
        if not ref_positions:
            continue
        for p0 in ref_positions:
            if (p0 + 1) in truth_set:
                hp3_tp += 1
                break

    bam.close()
    return total, hp3, hp3_tp


def make_plot(rows: list[dict], out_png: Path) -> None:
    """Dual bar chart over 24 chrs with outlier annotation."""
    fig, ax1 = plt.subplots(figsize=(16, 6.5))

    chroms = [r["chromosome"] for r in rows]
    fracs = [r["hp3_fraction"] * 100 for r in rows]  # percent
    rates = [r["hp3_tp_rate"] * 100 for r in rows]   # percent

    x = list(range(len(chroms)))
    bar_w = 0.4
    x_left = [i - bar_w / 2 for i in x]
    x_right = [i + bar_w / 2 for i in x]

    color_l = "#1f77b4"
    color_r = "#d62728"

    # Highlight outlier chrs (HP3 fraction > 5%) with red background span
    for i, r in enumerate(rows):
        if r["hp3_fraction"] > HP3_FRAC_OUTLIER_THRESHOLD:
            ax1.axvspan(i - 0.5, i + 0.5, facecolor="#ffcccc", alpha=0.5, zorder=0)
        if r["chromosome"] == CHR8_LOH_HOTSPOT:
            # light blue frame: draw rectangle border
            ax1.axvspan(i - 0.5, i + 0.5, facecolor="none",
                        edgecolor="#1f9bd6", linewidth=2.0, linestyle="--", zorder=1)

    b1 = ax1.bar(x_left, fracs, bar_w, color=color_l, label="HP3 fraction (%)", zorder=3)
    ax1.set_ylabel("HP3 read fraction (%)", color=color_l, fontsize=11)
    ax1.tick_params(axis="y", labelcolor=color_l)
    # Linear scale, dynamic upper bound (task spec says 0-2% but outliers can exceed)
    ymax_l = max(2.0, max(fracs) * 1.10)
    ax1.set_ylim(0, ymax_l)

    ax2 = ax1.twinx()
    b2 = ax2.bar(x_right, rates, bar_w, color=color_r, label="HP3 TP rate (%)", zorder=3)
    ax2.set_ylabel("HP3 TP read rate (%)", color=color_r, fontsize=11)
    ax2.tick_params(axis="y", labelcolor=color_r)
    ax2.set_ylim(0, 100)

    ax1.set_xticks(x)
    ax1.set_xticklabels(chroms, rotation=45, ha="right", fontsize=9)
    ax1.set_xlabel("Chromosome", fontsize=11)
    ax1.set_title(
        "HCC1395 Tmode tagged BAM — HP3 fraction vs HP3 TP rate (24-chr genome-wide)\n"
        "Red span = HP3 fraction outlier (>5%); dashed blue = chr8 LOH 99% hotspot",
        fontsize=12,
    )

    # Value labels — only for outliers / extreme cases to avoid clutter
    for rect, v, r in zip(b1, fracs, rows):
        if r["hp3_fraction"] > HP3_FRAC_OUTLIER_THRESHOLD or r["chromosome"] in {"chr1", "chr8", "chr19"}:
            ax1.text(
                rect.get_x() + rect.get_width() / 2,
                v + ymax_l * 0.01,
                f"{v:.2f}%",
                ha="center", va="bottom", fontsize=7, color=color_l,
                rotation=0,
            )
    for rect, v, r in zip(b2, rates, rows):
        if r["hp3_reads"] > 0:
            ax2.text(
                rect.get_x() + rect.get_width() / 2,
                v + 1.0,
                f"{v:.0f}%",
                ha="center", va="bottom", fontsize=7, color=color_r,
                rotation=0,
            )

    # Legend with outlier annotation entries
    outlier_patch = mpatches.Patch(facecolor="#ffcccc", alpha=0.5, label="HP3 fraction outlier (>5%)")
    hotspot_patch = mpatches.Patch(facecolor="none", edgecolor="#1f9bd6", linewidth=2.0,
                                   linestyle="--", label="chr8 LOH 99% hotspot")
    ax1.legend(handles=[b1, b2, outlier_patch, hotspot_patch],
               loc="upper left", fontsize=9, framealpha=0.95)

    ax1.grid(axis="y", linestyle=":", alpha=0.4, zorder=0)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", default=BAM_PATH)
    p.add_argument("--vcf", default=VCF_PATH)
    p.add_argument(
        "--out-tsv",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration/findings_5_24/T10_HP3_TP_rate_24chr.tsv",
    )
    p.add_argument(
        "--out-png",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration/figures/T10_HP3_TP_rate_24chr.png",
    )
    p.add_argument(
        "--chroms",
        default=",".join(CHROMS),
        help="Comma-separated chr list to scan (default: all 24). Use for resume.",
    )
    p.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip PNG generation (TSV-only run, for partial / resume runs).",
    )
    p.add_argument(
        "--append-tsv",
        action="store_true",
        help="Append to existing TSV instead of overwriting (resume mode).",
    )
    args = p.parse_args(argv)

    chroms_run = [c.strip() for c in args.chroms.split(",") if c.strip()]

    print(f"[T10] BAM: {args.bam}", file=sys.stderr)
    print(f"[T10] VCF: {args.vcf}", file=sys.stderr)
    print(f"[T10] CHROMS to scan: {len(chroms_run)} = {chroms_run}", file=sys.stderr)

    # 1. Load truth for requested chrs
    t0 = time.time()
    truth = load_truth_positions(args.vcf, chroms_run)
    truth_counts = {c: len(truth[c]) for c in chroms_run}
    print(f"[T10] truth load done in {time.time()-t0:.1f}s; total PASS = {sum(truth_counts.values()):,}",
          file=sys.stderr)
    for c in chroms_run:
        baseline = T1_BASELINE.get(c)
        if baseline:
            exp = baseline["truth_pass"]
            ok = "OK" if truth_counts[c] == exp else f"MISMATCH(expected {exp})"
            print(f"[T10] truth {c}: {truth_counts[c]:,}  [baseline {ok}]", file=sys.stderr)
        else:
            print(f"[T10] truth {c}: {truth_counts[c]:,}", file=sys.stderr)

    # 2. Scan BAM per chrom (stream once per chr; do not load full BAM)
    rows = []
    t_scan = time.time()
    for chrom in chroms_run:
        t_chr = time.time()
        total, hp3, hp3_tp = count_hp3_and_tp(args.bam, chrom, truth[chrom])
        elapsed = time.time() - t_chr
        hp3_frac = hp3 / total if total else 0.0
        hp3_tp_rate = hp3_tp / hp3 if hp3 else 0.0
        # Sanity vs T1
        baseline = T1_BASELINE.get(chrom)
        sanity = ""
        if baseline:
            checks = []
            if total != baseline["total"]:
                checks.append(f"total {total:,} != T1 {baseline['total']:,}")
            if hp3 != baseline["hp3"]:
                checks.append(f"hp3 {hp3:,} != T1 {baseline['hp3']:,}")
            if hp3_tp != baseline["hp3_tp"]:
                checks.append(f"hp3_tp {hp3_tp:,} != T1 {baseline['hp3_tp']:,}")
            sanity = " [T1 baseline " + ("OK" if not checks else "MISMATCH: " + "; ".join(checks)) + "]"
        print(
            f"[T10] {chrom}: total={total:,} hp3={hp3:,} ({hp3_frac:.4%}) "
            f"hp3_tp={hp3_tp:,} ({hp3_tp_rate:.4%}) in {elapsed:.1f}s{sanity}",
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
                "seqc2_truth_count": truth_counts[chrom],
            }
        )
    print(f"[T10] scan done in {(time.time()-t_scan)/60:.1f} min ({len(chroms_run)} chr)", file=sys.stderr)

    # 3. TSV
    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    header = ("chromosome\ttotal_reads\thp3_reads\thp3_fraction\thp3_tp_reads\t"
              "hp3_tp_rate\tseqc2_truth_count\n")
    if args.append_tsv and out_tsv.exists():
        mode = "a"
    else:
        mode = "w"
    with out_tsv.open(mode) as fh:
        if mode == "w":
            fh.write(header)
        for r in rows:
            fh.write(
                f"{r['chromosome']}\t{r['total_reads']}\t{r['hp3_reads']}\t"
                f"{r['hp3_fraction']:.8f}\t{r['hp3_tp_reads']}\t{r['hp3_tp_rate']:.8f}\t"
                f"{r['seqc2_truth_count']}\n"
            )
    print(f"[T10] TSV {mode}-mode write: {out_tsv}", file=sys.stderr)

    # 4. Plot (skip if requested or partial run)
    if args.no_plot:
        print(f"[T10] PNG skipped (--no-plot)", file=sys.stderr)
    else:
        out_png = Path(args.out_png)
        make_plot(rows, out_png)
        print(f"[T10] PNG written: {out_png}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
