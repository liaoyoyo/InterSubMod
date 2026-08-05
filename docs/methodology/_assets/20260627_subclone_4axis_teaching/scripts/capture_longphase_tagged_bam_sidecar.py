#!/usr/bin/env python3
"""Capture an exact alignment-level HP/PS sidecar from a streamed tagged BAM."""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import pysam


EXPECTED_HP = {".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-bam", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    hp_counts = Counter()
    hp_ps_counts = Counter()
    alignment_classes = Counter()
    unknown_hp = Counter()
    rows = 0
    sorted_ok = True
    previous = None
    with pysam.AlignmentFile(str(args.input_bam), "rb", check_sq=False) as bam, \
            args.output.open("w", encoding="utf-8") as output:
        output.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
        rank = {name: index for index, name in enumerate(bam.references)}
        for alignment in bam.fetch(until_eof=True):
            alignment_classes["total"] += 1
            if alignment.is_unmapped:
                alignment_classes["unmapped"] += 1
                continue
            if alignment.is_secondary:
                alignment_classes["secondary"] += 1
            elif alignment.is_supplementary:
                alignment_classes["supplementary"] += 1
            else:
                alignment_classes["primary"] += 1
            chrom = alignment.reference_name
            start = int(alignment.reference_start)
            end = int(alignment.reference_end or start + 1)
            hp = str(alignment.get_tag("HP")) if alignment.has_tag("HP") else "."
            ps = str(alignment.get_tag("PS")) if alignment.has_tag("PS") else "."
            cigar_digest = hashlib.blake2b((alignment.cigarstring or "*").encode(), digest_size=8).hexdigest()
            hp_counts[hp] += 1
            if ps != ".":
                hp_ps_counts[hp] += 1
            if hp not in EXPECTED_HP:
                unknown_hp[hp] += 1
            order = (rank.get(chrom, 10**9), start)
            if previous is not None and order < previous:
                sorted_ok = False
            previous = order
            output.write(f"{chrom}\t{start}\t{end}\t{alignment.query_name}\t{alignment.flag}\t"
                         f"{alignment.mapping_quality}\t{cigar_digest}\t{hp}\t{ps}\n")
            rows += 1
    summary = {
        "schema_version": "1.0",
        "input_stream": str(args.input_bam),
        "output": str(args.output),
        "rows_mapped": rows,
        "alignment_classes": dict(alignment_classes),
        "HP_counts": dict(hp_counts),
        "HP_with_PS_counts": dict(hp_ps_counts),
        "unknown_HP_counts": dict(unknown_hp),
        "coordinate_sorted": sorted_ok,
        "pass": rows == alignment_classes["total"] - alignment_classes["unmapped"]
                and not unknown_hp and sorted_ok,
    }
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"CAPTURED tagged BAM stream: mapped={rows} HP={dict(hp_counts)} -> {args.output}", flush=True)
    raise SystemExit(0 if summary["pass"] else 1)


if __name__ == "__main__":
    main()
