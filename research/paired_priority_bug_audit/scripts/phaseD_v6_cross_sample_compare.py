#!/usr/bin/env python3
"""Phase D V6 cross-sample comparison.

Reads V6 ISM outputs from 4-5 samples, computes per-sample marker filter
metrics (NG≥3 coverage, TP rate, hp distribution) to validate V6 cross-sample
consistency.
"""
from __future__ import annotations

import csv
import os
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
BASE = REPO / "research/paired_priority_bug_audit/phaseD_v6_5sample"
SAMPLES = ["H1437", "H2009", "HCC1954", "HCC1937", "COLO829"]
RUNS = ["off_tp", "off_fp", "on_tp", "on_fp"]
NG_VALID = {"1", "2", "1-1", "2-1", "3", "11", "21", "33"}


def get_region_data(reads_tsv: str) -> tuple[str, int, dict[str, int]]:
    region_id = os.path.basename(os.path.dirname(os.path.dirname(reads_tsv)))
    buckets: set[str] = set()
    hp_counter: Counter[str] = Counter()
    try:
        with open(reads_tsv, "rb") as fh:
            header = fh.readline().decode().rstrip().split("\t")
            try:
                hp_idx = header.index("hp")
            except ValueError:
                return (region_id, 0, dict(hp_counter))
            for line in fh:
                cols = line.decode().rstrip().split("\t")
                if len(cols) > hp_idx:
                    val = cols[hp_idx]
                    hp_counter[val] += 1
                    if val in NG_VALID:
                        buckets.add(val)
    except Exception:
        pass
    return (region_id, len(buckets), dict(hp_counter))


def find_reads_tsv(run_dir: Path) -> list[str]:
    out = []
    for root, _, files in os.walk(run_dir):
        if "reads.tsv" in files:
            out.append(os.path.join(root, "reads.tsv"))
    return out


def aggregate_run(run_dir: Path) -> tuple[dict[str, int], dict[str, dict[str, int]]]:
    files = find_reads_tsv(run_dir)
    region_ng: dict[str, int] = {}
    region_hp: dict[str, dict[str, int]] = {}
    with ProcessPoolExecutor(max_workers=8) as ex:
        for region_id, ng, hp in ex.map(get_region_data, files, chunksize=200):
            region_ng[region_id] = ng
            region_hp[region_id] = hp
    return region_ng, region_hp


def main() -> int:
    print("[Phase D V6 cross-sample] starting...")
    summary_rows = []
    for sample in SAMPLES:
        sample_dir = BASE / sample
        if not sample_dir.exists():
            print(f"  [SKIP] {sample}: no output dir")
            continue
        sample_data = {}
        for run in RUNS:
            run_dir = sample_dir / run
            if not run_dir.exists():
                print(f"  [SKIP] {sample}_{run}: missing")
                continue
            ng, hp = aggregate_run(run_dir)
            sample_data[run] = (ng, hp)
            print(f"  {sample}_{run}: {len(ng)} regions")

        if not all(r in sample_data for r in RUNS):
            print(f"  [INCOMPLETE] {sample} skipped")
            continue

        # Marker filter (NG_off ≥ 3)
        tp_off = sum(1 for r in sample_data["off_tp"][0] if sample_data["off_tp"][0][r] >= 3)
        fp_off = sum(1 for r in sample_data["off_fp"][0] if sample_data["off_fp"][0][r] >= 3)
        rate_off = tp_off / max(tp_off + fp_off, 1)

        # NG_on=2 cell
        tp_on2 = sum(1 for r in sample_data["on_tp"][0] if sample_data["on_tp"][0][r] == 2)
        fp_on2 = sum(1 for r in sample_data["on_fp"][0] if sample_data["on_fp"][0][r] == 2)
        rate_on2 = tp_on2 / max(tp_on2 + fp_on2, 1)

        # hp distribution (off, TP+FP)
        hp_total: Counter[str] = Counter()
        for r, hp in sample_data["off_tp"][1].items():
            for k, v in hp.items():
                hp_total[k] += v
        for r, hp in sample_data["off_fp"][1].items():
            for k, v in hp.items():
                hp_total[k] += v

        h11 = hp_total.get("1-1", 0) + hp_total.get("11", 0)
        h21 = hp_total.get("2-1", 0) + hp_total.get("21", 0)
        h33 = hp_total.get("3", 0) + hp_total.get("33", 0)
        ratio = h11 / max(h21, 1)

        n_tp_total = len(sample_data["off_tp"][0])
        n_fp_total = len(sample_data["off_fp"][0])
        summary_rows.append({
            "sample": sample,
            "n_tp": n_tp_total,
            "n_fp": n_fp_total,
            "marker_off_n": tp_off + fp_off,
            "marker_off_tp": tp_off,
            "marker_off_fp": fp_off,
            "marker_off_rate": rate_off,
            "ng_on2_tp": tp_on2,
            "ng_on2_fp": fp_on2,
            "ng_on2_rate": rate_on2,
            "hp_1_1_count": h11,
            "hp_2_1_count": h21,
            "hp_3_count": h33,
            "hp_1_1_2_1_ratio": ratio,
        })

    # Output
    print("\n=== V6 cross-sample marker stats ===")
    print(f"{'sample':<10} {'TP':>6} {'FP':>5} {'NG≥3 N':>8} {'rate':>7} {'NG_on=2 rate':>13} {'h33':>10} {'1-1:2-1':>10}")
    for r in summary_rows:
        print(
            f"{r['sample']:<10} {r['n_tp']:>6} {r['n_fp']:>5} "
            f"{r['marker_off_n']:>8} {r['marker_off_rate']:>7.3f} "
            f"{r['ng_on2_rate']:>13.3f} {r['hp_3_count']:>10} {r['hp_1_1_2_1_ratio']:>10.3f}"
        )

    # Phase D 驗收
    print("\n=== Phase D 驗收（V6 跨樣本一致性）===")
    n_pass_marker = sum(1 for r in summary_rows if r["marker_off_rate"] >= 0.85)
    n_pass_ng_on2 = sum(1 for r in summary_rows if r["ng_on2_rate"] >= 0.85)
    print(f"  Marker rate ≥ 0.85: {n_pass_marker}/{len(summary_rows)}")
    print(f"  NG_on=2 rate ≥ 0.85: {n_pass_ng_on2}/{len(summary_rows)}")

    out_tsv = BASE / "v6_cross_sample_summary.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        if summary_rows:
            writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)
    print(f"\nSummary → {out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
