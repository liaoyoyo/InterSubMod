#!/usr/bin/env python3
"""V3-Fixed Ablation @ 0.6 — Phase block N50 from phased VCF (3 versions)."""

from __future__ import annotations
import csv, gzip, os, sys
from collections import defaultdict
from typing import Dict, List, Tuple

OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation"
os.makedirs(OUT_DIR, exist_ok=True)

VCFS = {
    "B1_baseline_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/baseline_06/tumor_phased.vcf",
    "B3_v3f_no_pononly_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v3f_no_pononly_06/tumor_phased.vcf",
    "B5_v5_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v5_06/tumor_phased.vcf",
}


def parse_phased_vcf(path):
    blocks = defaultdict(list)
    total = phased = unphased = 0
    f = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    with f:
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
            if "|" in gt and ps and ps != "." and ps != "0":
                phased += 1
                blocks[(chrom, ps)].append(pos)
            else:
                unphased += 1
    return blocks, total, phased, unphased


def calc_n50(sizes):
    if not sizes:
        return 0
    s = sorted(sizes, reverse=True)
    half = sum(s) / 2
    cum = 0
    for x in s:
        cum += x
        if cum >= half:
            return x
    return s[-1]


def main():
    rows = []
    for label, vcf_path in VCFS.items():
        if not os.path.exists(vcf_path):
            print(f"WARN: {label} VCF missing", file=sys.stderr)
            continue
        print(f"Processing {label}...", file=sys.stderr)
        blocks, total, phased, unphased = parse_phased_vcf(vcf_path)
        sizes = [max(p) - min(p) + 1 for p in blocks.values() if len(p) >= 1]
        rows.append({
            "version": label,
            "total_variants": total,
            "phased_variants": phased,
            "phased_rate_pct": f"{phased/total*100:.2f}" if total else "0",
            "n_blocks": len(blocks),
            "N50_bp": calc_n50(sizes),
            "max_block_bp": max(sizes) if sizes else 0,
        })

    out = os.path.join(OUT_DIR, "phase_block_stats_06.tsv")
    if rows:
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {out}", file=sys.stderr)
        print(f"\n{'version':<28} {'phased%':>8} {'n_blocks':>10} {'N50_bp':>14} {'max_bp':>14}", file=sys.stderr)
        for r in rows:
            print(f"{r['version']:<28} {r['phased_rate_pct']:>8} {r['n_blocks']:>10} {r['N50_bp']:>14} {r['max_block_bp']:>14}", file=sys.stderr)


if __name__ == "__main__":
    main()
