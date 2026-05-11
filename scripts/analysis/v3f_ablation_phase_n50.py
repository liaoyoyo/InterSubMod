#!/usr/bin/env python3
"""V3-Fixed Ablation — Phase block N50 calculation from phased VCF.

Computes phase block size distribution (PS field) for 5 phased VCFs:
  A1 baseline / A2 V2b / A3 v3f_no_pononly / A4 V3F+PONonly / A5 V5

A5 uses V2b's phased VCF (V5 reuses V2b phase output for haplotag);
so A5 = A2 for phase metrics.

Output:
  InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/phase_block_stats.tsv
"""

from __future__ import annotations

import csv
import gzip
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation"
os.makedirs(OUT_DIR, exist_ok=True)

VCFS: Dict[str, str] = {
    "A1_baseline": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased.vcf",
    "A2_v2b": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_phased.vcf",
    "A3_v3f_no_pononly": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/tumor_phased.vcf",
    "A4_v3f_with_pononly": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_phased.vcf",  # V3F reuses V2b VCF
    "A5_v5": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_phased.vcf",  # V5 reuses V2b VCF
}


def opener(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_phased_vcf(path: str) -> Tuple[Dict[Tuple[str, str], List[int]], int, int, int]:
    """Parse VCF, return:
      blocks: {(chrom, ps_id): [pos1, pos2, ...]}
      total_variants, phased_variants, unphased_variants
    """
    blocks: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    total = phased = unphased = 0

    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            fmt_keys = fields[8].split(":")
            sample = fields[9].split(":")
            fmt = dict(zip(fmt_keys, sample))

            gt = fmt.get("GT", "")
            ps = fmt.get("PS", ".")

            total += 1

            # Phased = "|" in GT
            if "|" in gt and ps and ps != "." and ps != "0":
                phased += 1
                blocks[(chrom, ps)].append(pos)
            else:
                unphased += 1

    return blocks, total, phased, unphased


def calc_n50(sizes: List[int]) -> int:
    """Calculate N50 from list of block sizes (in bp)."""
    if not sizes:
        return 0
    sizes_sorted = sorted(sizes, reverse=True)
    total = sum(sizes_sorted)
    half = total / 2
    cumsum = 0
    for s in sizes_sorted:
        cumsum += s
        if cumsum >= half:
            return s
    return sizes_sorted[-1]


def block_stats(blocks: Dict[Tuple[str, str], List[int]]) -> Dict:
    """Compute block stats: count, sizes, N50."""
    sizes = []
    n_variants_per_block = []
    for (chrom, ps), positions in blocks.items():
        if len(positions) < 1:
            continue
        size = max(positions) - min(positions) + 1
        sizes.append(size)
        n_variants_per_block.append(len(positions))

    return {
        "n_blocks": len(blocks),
        "total_block_span_bp": sum(sizes),
        "median_block_size_bp": sorted(sizes)[len(sizes) // 2] if sizes else 0,
        "max_block_size_bp": max(sizes) if sizes else 0,
        "N50_bp": calc_n50(sizes),
        "median_variants_per_block": sorted(n_variants_per_block)[len(n_variants_per_block) // 2] if n_variants_per_block else 0,
    }


def main():
    rows = []

    for label, vcf_path in VCFS.items():
        if not os.path.exists(vcf_path):
            print(f"WARN: {label} VCF not found at {vcf_path}", file=sys.stderr)
            continue

        print(f"Processing {label}...", file=sys.stderr)
        blocks, total, phased, unphased = parse_phased_vcf(vcf_path)
        stats = block_stats(blocks)

        phased_rate = phased / total * 100 if total > 0 else 0

        row = {
            "version": label,
            "vcf_path": vcf_path,
            "total_variants": total,
            "phased_variants": phased,
            "unphased_variants": unphased,
            "phased_rate_pct": f"{phased_rate:.2f}",
            **stats,
        }
        rows.append(row)

    # Write output
    out_path = os.path.join(OUT_DIR, "phase_block_stats.tsv")
    if rows:
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {out_path}", file=sys.stderr)

        print("\n=== Phase block stats summary ===", file=sys.stderr)
        print(f"{'version':<25} {'phased%':>8} {'n_blocks':>10} {'N50_bp':>12} {'max_bp':>14}", file=sys.stderr)
        for r in rows:
            print(
                f"{r['version']:<25} {r['phased_rate_pct']:>8} {r['n_blocks']:>10} "
                f"{r['N50_bp']:>12} {r['max_block_size_bp']:>14}",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
