"""Shared utilities for HKU collaboration analyses (5/23).

BAM provenance (from @PG):
  - SNVs from ClairS-TO ssrs v0.4.x
  - Tagging by LongPhase mainline v1.7.3 `haplotag --somaticMode`
  - Paired-with-normal: -s HCC1395BL_methyl_phase.vcf.gz -b HCC1395BL.bam
  - SEQC2 --highCon-snp guidance
  - HP tag values: 1/2/3/5/7 (H1/H2/H3/H1_1/H2_1)
"""
from __future__ import annotations

import gzip
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# --- CJK font injection (per memory feedback_matplotlib_cjk_font_rule.md) ---
def setup_cjk_font():
    rcParams["font.sans-serif"] = [
        "Droid Sans Fallback",
        "Noto Sans CJK TC",
        "Noto Sans CJK SC",
        "DejaVu Sans",
        "sans-serif",
    ]
    rcParams["axes.unicode_minus"] = False
    rcParams["figure.dpi"] = 150
    rcParams["savefig.dpi"] = 150


# --- Paths ---
BAM_PATH = "/big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam"
TP_VCF = "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
LOH_BED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
REF_FASTA = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"

OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration")
FIG_DIR = OUT_DIR / "figures"
DATA_DIR = OUT_DIR / "data"

PILOT_CHROMS = ["chr1", "chr8", "chr19"]  # chr1+chr8 main, chr19 reference


def load_tp_positions(vcf_path: str = TP_VCF, chroms: list[str] | None = None) -> dict[str, list[int]]:
    """Load TP somatic SNV positions from SEQC2 VCF. PASS;HighConf only.

    Returns: dict[chrom] -> sorted list of positions (1-based VCF coords).
    """
    positions: dict[str, list[int]] = {}
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, flt = parts[:7]
            if "PASS" not in flt:
                continue
            # Keep only SNV (single base ref & alt)
            if len(ref) != 1 or len(alt) != 1:
                continue
            if chroms and chrom not in chroms:
                continue
            positions.setdefault(chrom, []).append(int(pos))
    for c in positions:
        positions[c].sort()
    return positions


def load_loh_regions(bed_path: str = LOH_BED, chroms: list[str] | None = None) -> dict[str, list[tuple[int, int]]]:
    """Load LOH regions (type=='loh'). Returns dict[chrom] -> list of (start, end) (BED 0-based)."""
    regions: dict[str, list[tuple[int, int]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, typ = parts[0], int(parts[1]), int(parts[2]), parts[3]
            if typ != "loh":
                continue
            if chroms and chrom not in chroms:
                continue
            regions.setdefault(chrom, []).append((start, end))
    for c in regions:
        regions[c].sort()
    return regions


def interval_overlap(intervals: list[tuple[int, int]], qstart: int, qend: int) -> bool:
    """Check if [qstart, qend) overlaps any (start, end) in intervals (sorted)."""
    if not intervals:
        return False
    import bisect
    starts = [iv[0] for iv in intervals]
    idx = bisect.bisect_right(starts, qend) - 1
    while idx >= 0:
        s, e = intervals[idx]
        if e <= qstart:
            return False
        if s < qend and e > qstart:
            return True
        idx -= 1
    return False


def log_timing(label: str, t0: float, t1: float):
    print(f"[TIMING] {label}: {t1 - t0:.2f}s", file=sys.stderr, flush=True)
