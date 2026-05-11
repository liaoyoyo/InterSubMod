#!/usr/bin/env python3
"""
Per-region BAM read-level QC extraction for feature-layered observation (G4).

For each (sample, mode) combination in the master_extended dataset, open the
corresponding tumor_tagged.bam and compute per-region QC statistics centred on
each variant position with a fixed window (default 5000 bp).

Features extracted per region:
    n_reads            : number of primary reads overlapping the window
    MapQ_mean          : mean mapping quality
    MapQ_median        : median mapping quality
    MapQ_p10           : 10th percentile mapping quality
    MapQ_p90           : 90th percentile mapping quality
    LowMQ20_Frac       : fraction of reads with MapQ < 20
    NM_mean            : mean number of mismatches (NM tag)
    Softclip_Frac      : fraction of reads containing any soft-clipped base
    Strand_Bias        : |n_forward - n_reverse| / n_reads
    Read_Length_mean   : mean read length

Output: one TSV per (sample, mode):
    bam_readqc_{sample}_{mode}.tsv

The script is idempotent — if an output TSV already exists and is non-empty the
corresponding (sample, mode) pair is skipped. Use --force to overwrite.

Usage:
    bam_readqc_per_region.py --sample HCC1395 --mode paired_full [--limit N] [--chr chr19]
    bam_readqc_per_region.py --all                                  # run all combos sequentially
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pysam

# ---------------------------------------------------------------------------
# Paths & configuration
# ---------------------------------------------------------------------------

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER_EXT = ROOT / "research/tpfp_loh_af_kde_discrimination/data/master_extended.tsv.gz"
OUT_DIR = ROOT / "research/feature_layered_observation/data/bam_readqc"
LOG_DIR = ROOT / "research/feature_layered_observation/data"
DEFAULT_WINDOW = 5000  # bp total window centred on Pos

# BAM path table: (sample, mode) -> absolute path to tumor_tagged.bam
BAM_PATHS: dict[tuple[str, str], Path] = {
    # -------- paired_full (canonical complete_matrix directory) -----------
    ("HCC1395", "paired_full"): ROOT / "output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam",
    ("HCC1395_DORADO", "paired_full"): ROOT / "output/canonical/HCC1395_DORADO/paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix/longphase_s/HCC1395_DORADO_tagged.bam",
    ("H1437", "paired_full"): ROOT / "output/canonical/H1437/paired_full/20260315_H1437_paired_full_full_complete_matrix/longphase_s/H1437_tagged.bam",
    ("H2009", "paired_full"): ROOT / "output/canonical/H2009/paired_full/20260315_H2009_paired_full_full_complete_matrix/longphase_s/H2009_tagged.bam",
    ("HCC1937", "paired_full"): ROOT / "output/canonical/HCC1937/paired_full/20260315_HCC1937_paired_full_full_complete_matrix/longphase_s/HCC1937_tagged.bam",
    ("HCC1954", "paired_full"): ROOT / "output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix/longphase_s/HCC1954_tagged.bam",
    ("COLO829", "paired_full"): ROOT / "output/canonical/COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/longphase_s/COLO829_tagged.bam",
    # -------- to_pileup (TO pipeline fastresume / pilot) ------------------
    ("HCC1395", "to_pileup"): ROOT / "output/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam",
    ("HCC1395_DORADO", "to_pileup"): ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_tagged.bam",
    ("H1437", "to_pileup"): ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step03_longphase_to/tumor_tagged.bam",
    ("H2009", "to_pileup"): ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_tagged.bam",
    ("HCC1937", "to_pileup"): ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step03_longphase_to/tumor_tagged.bam",
    ("HCC1954", "to_pileup"): ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step03_longphase_to/tumor_tagged.bam",
    ("COLO829", "to_pileup"): ROOT / "output/synthesis/research_rounds/20260423_colo829_to_pilot/step03_longphase_to/tumor_tagged.bam",
}

OUTPUT_COLUMNS = [
    "RegionID", "sample", "mode", "Chr", "Pos",
    "win_start", "win_end",
    "n_reads",
    "MapQ_mean", "MapQ_median", "MapQ_p10", "MapQ_p90", "LowMQ20_Frac",
    "NM_mean",
    "Softclip_Frac",
    "Strand_Bias",
    "Read_Length_mean",
]


# ---------------------------------------------------------------------------
# Region iteration
# ---------------------------------------------------------------------------

@dataclass
class Region:
    region_id: str
    chrom: str
    pos: int

    def window(self, half: int) -> tuple[int, int]:
        start = max(0, self.pos - half)
        end = self.pos + half
        return start, end


def iter_regions(sample: str, mode: str, chr_filter: str | None = None, limit: int | None = None) -> Iterator[Region]:
    """Stream (RegionID, Chr, Pos) rows from the master_extended TSV for the
    target (sample, mode). Reads the gzipped file line-by-line for low memory
    footprint."""
    with gzip.open(MASTER_EXT, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx_rid = header.index("RegionID")
        idx_chr = header.index("Chr")
        idx_pos = header.index("Pos")
        idx_sample = header.index("sample")
        idx_mode = header.index("mode")
        n = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[idx_sample] != sample or parts[idx_mode] != mode:
                continue
            chrom = parts[idx_chr]
            if chr_filter is not None and chrom != chr_filter:
                continue
            yield Region(parts[idx_rid], chrom, int(parts[idx_pos]))
            n += 1
            if limit is not None and n >= limit:
                return


# ---------------------------------------------------------------------------
# Per-region BAM QC extraction
# ---------------------------------------------------------------------------

def compute_region_qc(bam: pysam.AlignmentFile, region: Region, window_half: int) -> dict:
    """Fetch primary reads in [pos - half, pos + half] and compute QC stats."""
    start, end = region.window(window_half)
    try:
        reads = bam.fetch(region.chrom, start, end, multiple_iterators=False)
    except (ValueError, KeyError):
        # chrom may not be present in BAM header; return empty stats
        return _empty_stats(region, start, end)

    mapqs: list[int] = []
    nms: list[int] = []
    softclip_flag: list[bool] = []
    forward = 0
    reverse = 0
    lengths: list[int] = []

    for r in reads:
        # Skip secondary, supplementary, duplicates, unmapped
        if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_unmapped:
            continue
        mapqs.append(r.mapping_quality)
        nm = r.get_tag("NM") if r.has_tag("NM") else None
        if nm is not None:
            nms.append(int(nm))
        cigar = r.cigartuples or []
        has_sc = any(op == 4 for op, _ in cigar)  # BAM_CSOFT_CLIP = 4
        softclip_flag.append(has_sc)
        if r.is_reverse:
            reverse += 1
        else:
            forward += 1
        if r.query_length:
            lengths.append(r.query_length)

    n = len(mapqs)
    if n == 0:
        return _empty_stats(region, start, end)

    mapqs_a = np.asarray(mapqs)
    stats = {
        "RegionID": region.region_id,
        "Chr": region.chrom,
        "Pos": region.pos,
        "win_start": start,
        "win_end": end,
        "n_reads": n,
        "MapQ_mean": float(mapqs_a.mean()),
        "MapQ_median": float(np.median(mapqs_a)),
        "MapQ_p10": float(np.percentile(mapqs_a, 10)),
        "MapQ_p90": float(np.percentile(mapqs_a, 90)),
        "LowMQ20_Frac": float((mapqs_a < 20).sum() / n),
        "NM_mean": float(np.mean(nms)) if nms else float("nan"),
        "Softclip_Frac": float(sum(softclip_flag) / n),
        "Strand_Bias": float(abs(forward - reverse) / n),
        "Read_Length_mean": float(np.mean(lengths)) if lengths else float("nan"),
    }
    return stats


def _empty_stats(region: Region, start: int, end: int) -> dict:
    return {
        "RegionID": region.region_id,
        "Chr": region.chrom,
        "Pos": region.pos,
        "win_start": start,
        "win_end": end,
        "n_reads": 0,
        "MapQ_mean": float("nan"),
        "MapQ_median": float("nan"),
        "MapQ_p10": float("nan"),
        "MapQ_p90": float("nan"),
        "LowMQ20_Frac": float("nan"),
        "NM_mean": float("nan"),
        "Softclip_Frac": float("nan"),
        "Strand_Bias": float("nan"),
        "Read_Length_mean": float("nan"),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def output_path(sample: str, mode: str, chr_filter: str | None = None) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix = f"_{chr_filter}" if chr_filter else ""
    return OUT_DIR / f"bam_readqc_{sample}_{mode}{suffix}.tsv"


def process_one(sample: str, mode: str, window_half: int, threads: int,
                chr_filter: str | None, limit: int | None, force: bool,
                log: callable) -> int:
    out = output_path(sample, mode, chr_filter)
    if out.exists() and out.stat().st_size > 0 and not force:
        log(f"[skip] {sample}/{mode}: output already exists ({out})")
        return 0

    bam_path = BAM_PATHS.get((sample, mode))
    if bam_path is None or not bam_path.exists():
        log(f"[ERROR] BAM not found for {sample}/{mode}: {bam_path}")
        return 1

    log(f"[start] {sample}/{mode} bam={bam_path}")
    t0 = time.time()
    bam = pysam.AlignmentFile(str(bam_path), "rb", threads=threads)

    tmp = out.with_suffix(out.suffix + ".tmp")
    written = 0
    with open(tmp, "w") as fh:
        fh.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for i, region in enumerate(iter_regions(sample, mode, chr_filter, limit)):
            stats = compute_region_qc(bam, region, window_half)
            row = [stats["RegionID"], sample, mode, stats["Chr"], str(stats["Pos"]),
                   str(stats["win_start"]), str(stats["win_end"]),
                   str(stats["n_reads"])]
            for c in ("MapQ_mean", "MapQ_median", "MapQ_p10", "MapQ_p90", "LowMQ20_Frac",
                      "NM_mean", "Softclip_Frac", "Strand_Bias", "Read_Length_mean"):
                v = stats[c]
                row.append("" if (v != v) else f"{v:.6g}")  # NaN check
            fh.write("\t".join(row) + "\n")
            written += 1
            if written % 2000 == 0:
                elapsed = time.time() - t0
                log(f"[progress] {sample}/{mode}: {written} regions, {elapsed:.1f}s, "
                    f"{written/elapsed:.1f} reg/s")
    bam.close()
    tmp.replace(out)
    elapsed = time.time() - t0
    log(f"[done] {sample}/{mode}: {written} regions in {elapsed:.1f}s -> {out}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", help="sample name (e.g. HCC1395)")
    ap.add_argument("--mode", help="mode (paired_full | to_pileup)")
    ap.add_argument("--all", action="store_true", help="run all (sample, mode) combos in BAM_PATHS")
    ap.add_argument("--window", type=int, default=DEFAULT_WINDOW, help="full window size in bp (default 5000)")
    ap.add_argument("--threads", type=int, default=4, help="pysam I/O threads per BAM (default 4)")
    ap.add_argument("--chr", dest="chr_filter", default=None, help="restrict to one chromosome (for smoke test)")
    ap.add_argument("--limit", type=int, default=None, help="cap number of regions (for smoke test)")
    ap.add_argument("--force", action="store_true", help="overwrite existing output")
    ap.add_argument("--log", default=None, help="append log to this file (stdout+file)")
    args = ap.parse_args()

    def log(msg: str) -> None:
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
        print(line, flush=True)
        if args.log:
            with open(args.log, "a") as fh:
                fh.write(line + "\n")

    half = args.window // 2
    status = 0
    if args.all:
        pairs = list(BAM_PATHS.keys())
    else:
        if not args.sample or not args.mode:
            ap.error("either --all, or both --sample and --mode are required")
        pairs = [(args.sample, args.mode)]

    for sample, mode in pairs:
        rc = process_one(sample, mode, half, args.threads, args.chr_filter,
                         args.limit, args.force, log)
        if rc != 0:
            status = rc
    return status


if __name__ == "__main__":
    sys.exit(main())
