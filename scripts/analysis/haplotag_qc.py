#!/usr/bin/env python3
"""Summarize HP-tag assignment from a BAM file into haplotag_qc.tsv."""

import argparse
import csv
import os
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, help="Tagged BAM path")
    parser.add_argument("--sample", default="unknown", help="Sample name")
    parser.add_argument("--output-tsv", required=True, help="Output TSV path")
    parser.add_argument("--samtools", default="samtools", help="samtools binary")
    return parser.parse_args()


def extract_hp_tag(fields):
    for field in fields[11:]:
        if field.startswith("HP:"):
            parts = field.split(":", 2)
            if len(parts) == 3:
                return parts[2]
    return ""


def main():
    args = parse_args()

    counts = {
        "reads_total": 0,
        "hp1": 0,
        "hp2": 0,
        "hp_other": 0,
        "hp_missing": 0,
    }

    proc = subprocess.Popen(
        [args.samtools, "view", args.bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert proc.stdout is not None
    for line in proc.stdout:
        if not line.strip():
            continue
        counts["reads_total"] += 1
        fields = line.rstrip("\n").split("\t")
        hp = extract_hp_tag(fields)
        if hp in {"1", "HP1"}:
            counts["hp1"] += 1
        elif hp in {"2", "HP2"}:
            counts["hp2"] += 1
        elif hp:
            counts["hp_other"] += 1
        else:
            counts["hp_missing"] += 1

    stderr = proc.communicate()[1]
    if proc.returncode != 0:
        sys.stderr.write(stderr)
        raise SystemExit(proc.returncode)

    assigned = counts["hp1"] + counts["hp2"] + counts["hp_other"]
    assign_rate = assigned / counts["reads_total"] if counts["reads_total"] else 0.0

    os.makedirs(os.path.dirname(args.output_tsv) or ".", exist_ok=True)
    with open(args.output_tsv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample",
                "bam",
                "reads_total",
                "hp1",
                "hp2",
                "hp_other",
                "hp_missing",
                "hp_assign_rate",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample": args.sample,
                "bam": args.bam,
                "reads_total": counts["reads_total"],
                "hp1": counts["hp1"],
                "hp2": counts["hp2"],
                "hp_other": counts["hp_other"],
                "hp_missing": counts["hp_missing"],
                "hp_assign_rate": f"{assign_rate:.6f}",
            }
        )

    print(f"[haplotag_qc] Wrote {args.output_tsv}")


if __name__ == "__main__":
    main()
