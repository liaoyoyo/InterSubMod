#!/usr/bin/env python3
"""Fast region-level NG aggregation using os.scandir + parallel processes.

Computes per-region NG_off/NG_on for each (BAM, label) combination, then aggregates
marker filter (NG≥3) TP rates for V3F vs V5 vs V6 head-to-head.
"""
from __future__ import annotations

import csv
import os
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
BASE = REPO / "research/paired_priority_bug_audit/phaseC_genome_three_way"

# Buckets that count toward NG (excluding "0" unphased and empty)
NG_VALID = {"1", "2", "1-1", "2-1", "3", "11", "21", "33"}


def get_region_ng(reads_tsv: str) -> tuple[str, int]:
    """Return (region_id, NG count)."""
    region_id = os.path.basename(os.path.dirname(os.path.dirname(reads_tsv)))
    buckets: set[str] = set()
    try:
        with open(reads_tsv, "rb") as fh:
            header = fh.readline().decode().rstrip().split("\t")
            try:
                hp_idx = header.index("hp")
            except ValueError:
                return (region_id, 0)
            for line in fh:
                cols = line.decode().rstrip().split("\t")
                if len(cols) > hp_idx:
                    val = cols[hp_idx]
                    if val in NG_VALID:
                        buckets.add(val)
    except Exception:
        pass
    return (region_id, len(buckets))


def find_reads_tsv(run_dir: Path) -> list[str]:
    """Find all reads.tsv in run_dir using os.walk (faster than rglob)."""
    out = []
    for root, _, files in os.walk(run_dir):
        if "reads.tsv" in files:
            out.append(os.path.join(root, "reads.tsv"))
    return out


def aggregate_run(run_dir: Path) -> dict[str, int]:
    files = find_reads_tsv(run_dir)
    region_ng: dict[str, int] = {}
    with ProcessPoolExecutor(max_workers=8) as ex:
        for region_id, ng in ex.map(get_region_ng, files, chunksize=200):
            region_ng[region_id] = ng
    return region_ng


def main() -> int:
    print("[Phase C region-level NG fast agg] starting...")
    bams = ["V3F", "V5", "V6"]
    runs = ["off_tp", "off_fp", "on_tp", "on_fp"]
    data: dict[str, dict[str, dict[str, int]]] = {}
    for bam in bams:
        data[bam] = {}
        for run in runs:
            run_dir = BASE / f"{bam}_{run}"
            if not run_dir.exists():
                print(f"  [SKIP] {bam}_{run} missing")
                continue
            print(f"  Processing {bam}_{run}...", flush=True)
            data[bam][run] = aggregate_run(run_dir)
            print(f"    {len(data[bam][run])} regions")

    # Build per-region NG_off / NG_on table per BAM × label
    summary_rows = []
    for bam in bams:
        if bam not in data:
            continue
        for label, off_run, on_run in [("TP", "off_tp", "on_tp"), ("FP", "off_fp", "on_fp")]:
            off = data[bam].get(off_run, {})
            on = data[bam].get(on_run, {})
            for region_id in sorted(set(off) | set(on)):
                summary_rows.append({
                    "BAM": bam,
                    "label": label,
                    "region_id": region_id,
                    "ng_off": off.get(region_id, -1),
                    "ng_on": on.get(region_id, -1),
                })

    # Output per-row TSV
    out_tsv = BASE / "v3f_vs_v5_vs_v6_region_ng.tsv"
    with out_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["BAM", "label", "region_id", "ng_off", "ng_on"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"\nPer-region table → {out_tsv} ({len(summary_rows)} rows)")

    # Aggregate marker (NG≥3) per BAM
    print("\n=== Marker filter (NG_off ≥ 3) — three-way (genome-wide) ===")
    print(f"{'BAM':<6} {'NG≥3 N':>10} {'TP':>7} {'FP':>5} {'TP rate':>9} | {'NG_on=2 TP':>11} {'FP':>5} {'TP rate':>9}")
    summary_agg = []
    for bam in bams:
        if bam not in data:
            continue
        rows = [r for r in summary_rows if r["BAM"] == bam]
        tp_off = sum(1 for r in rows if r["label"] == "TP" and r["ng_off"] >= 3)
        fp_off = sum(1 for r in rows if r["label"] == "FP" and r["ng_off"] >= 3)
        tp_on2 = sum(1 for r in rows if r["label"] == "TP" and r["ng_on"] == 2)
        fp_on2 = sum(1 for r in rows if r["label"] == "FP" and r["ng_on"] == 2)
        rate_off = tp_off / max(tp_off + fp_off, 1)
        rate_on = tp_on2 / max(tp_on2 + fp_on2, 1)
        print(
            f"{bam:<6} {tp_off+fp_off:>10} {tp_off:>7} {fp_off:>5} "
            f"{rate_off:>9.4f} | {tp_on2:>11} {fp_on2:>5} {rate_on:>9.4f}"
        )
        summary_agg.append({
            "BAM": bam,
            "marker_off_n": tp_off + fp_off,
            "marker_off_tp": tp_off,
            "marker_off_fp": fp_off,
            "marker_off_rate": rate_off,
            "ng_on2_tp": tp_on2,
            "ng_on2_fp": fp_on2,
            "ng_on2_rate": rate_on,
        })

    out_agg = BASE / "v3f_vs_v5_vs_v6_genome_summary.tsv"
    with out_agg.open("w", newline="") as fh:
        if summary_agg:
            writer = csv.DictWriter(fh, fieldnames=list(summary_agg[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_agg)
    print(f"\nMarker summary → {out_agg}")

    # Per-cell TP rate per BAM
    for bam in bams:
        if bam not in data:
            continue
        rows = [r for r in summary_rows if r["BAM"] == bam]
        cell_tp: dict[tuple[int, int], int] = defaultdict(int)
        cell_fp: dict[tuple[int, int], int] = defaultdict(int)
        for r in rows:
            key = (r["ng_off"], r["ng_on"])
            if r["label"] == "TP":
                cell_tp[key] += 1
            else:
                cell_fp[key] += 1
        print(f"\n=== {bam} BAM per-cell TP rate (genome-wide) ===")
        print(f"{'NG_off':>7} {'NG_on':>6} {'TP':>7} {'FP':>5} {'TP rate':>9}")
        for noff, non in sorted(set(cell_tp) | set(cell_fp)):
            tp = cell_tp[(noff, non)]
            fp = cell_fp[(noff, non)]
            rate = tp / max(tp + fp, 1)
            print(f"{noff:>7} {non:>6} {tp:>7} {fp:>5} {rate:>9.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
