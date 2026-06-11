#!/usr/bin/env python3
"""V6 per-chr somatic ALT HP1:HP2 ratio analysis.

Tests the 3 priority-bug-disproof arguments for V6:
  ① Biology: 23 chr 應無系統偏移 → V6 per-chr 是否中性
  ② cross-chr consistency: baseline 23 chr 全部偏 HP1 (94.6%) → V6 是否消除
  ③ Paired control: paired mode HP1:HP2 ≈ 1:1 → V6 是否接近

Filter: alt_support='ALT' AND is_tumor=1
Group by chr. Count hp ∈ {1, 1-1, 11} (HP1 系列) vs {2, 2-1, 21} (HP2 系列).
"""
from __future__ import annotations

import csv
import os
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

REPO = "/big7_disk/liaoyoyo2001/InterSubMod"

SOURCES = {
    "V3F_HCC1395": REPO + "/research/paired_priority_bug_audit/phaseC_genome_three_way/V3F_off_tp",
    "V5_HCC1395": REPO + "/research/paired_priority_bug_audit/phaseC_genome_three_way/V5_off_tp",
    "V6_HCC1395": REPO + "/research/paired_priority_bug_audit/phaseC_genome_three_way/V6_off_tp",
    "V6_H1437": REPO + "/research/paired_priority_bug_audit/phaseD_v6_5sample/H1437/off_tp",
    "V6_H2009": REPO + "/research/paired_priority_bug_audit/phaseD_v6_5sample/H2009/off_tp",
    "V6_HCC1954": REPO + "/research/paired_priority_bug_audit/phaseD_v6_5sample/HCC1954/off_tp",
    "V6_HCC1937": REPO + "/research/paired_priority_bug_audit/phaseD_v6_5sample/HCC1937/off_tp",
}

HP1_SERIES = {"1", "1-1", "11"}
HP2_SERIES = {"2", "2-1", "21"}

CHRS_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def count_per_chr(reads_tsv: str) -> dict[str, tuple[int, int]]:
    """Returns {chr: (hp1_alt_count, hp2_alt_count)} for ALT+tumor reads only."""
    counts: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    try:
        with open(reads_tsv, "rb") as fh:
            header = fh.readline().decode().rstrip().split("\t")
            try:
                chr_idx = header.index("chr")
                hp_idx = header.index("hp")
                alt_idx = header.index("alt_support")
                tumor_idx = header.index("is_tumor")
            except ValueError:
                return {}
            for line in fh:
                cols = line.decode().rstrip().split("\t")
                if len(cols) <= max(chr_idx, hp_idx, alt_idx, tumor_idx):
                    continue
                if cols[alt_idx] != "ALT" or cols[tumor_idx] != "1":
                    continue
                hp = cols[hp_idx]
                if hp in HP1_SERIES:
                    counts[cols[chr_idx]][0] += 1
                elif hp in HP2_SERIES:
                    counts[cols[chr_idx]][1] += 1
    except Exception:
        pass
    return {k: tuple(v) for k, v in counts.items()}


def find_reads_tsv(base: str) -> list[str]:
    out = []
    for root, _, files in os.walk(base):
        if "reads.tsv" in files:
            out.append(os.path.join(root, "reads.tsv"))
    return out


def aggregate_source(name: str, base: str) -> dict[str, tuple[int, int]]:
    """Aggregate per-chr HP1/HP2 counts across all region reads.tsv."""
    files = find_reads_tsv(base)
    print(f"  {name}: {len(files)} reads.tsv files")
    chr_counts: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    with ProcessPoolExecutor(max_workers=8) as ex:
        for region_counts in ex.map(count_per_chr, files, chunksize=200):
            for chrom, (h1, h2) in region_counts.items():
                chr_counts[chrom][0] += h1
                chr_counts[chrom][1] += h2
    return {k: tuple(v) for k, v in chr_counts.items()}


def main() -> int:
    print("[V6 per-chr somatic ALT HP1:HP2 ratio]")
    print("Filter: alt_support='ALT' AND is_tumor=1")
    print("HP1 series = {1, 1-1, 11}; HP2 series = {2, 2-1, 21}")
    print()

    all_data: dict[str, dict[str, tuple[int, int]]] = {}
    for name, base in SOURCES.items():
        if not os.path.exists(base):
            print(f"  [SKIP] {name}: {base} missing")
            continue
        all_data[name] = aggregate_source(name, base)

    # Per-chr cross-source matrix (BAM × chr ratio)
    print("\n=== Per-chr ALT HP1:HP2 ratio (somatic ALT reads, tumor only) ===")
    print(f"{'chr':<8}", end="")
    for name in all_data:
        print(f" {name:>15}", end="")
    print()

    for chrom in CHRS_ORDER:
        print(f"{chrom:<8}", end="")
        for name, data in all_data.items():
            h1, h2 = data.get(chrom, (0, 0))
            if h1 + h2 == 0:
                print(f" {'n/a':>15}", end="")
            else:
                ratio = h1 / max(h2, 1)
                print(f" {h1:>6}:{h2:<6}={ratio:>4.2f}", end="")
        print()

    # Summary: count of chr biased > 5:1 toward HP1 (argument ② threshold)
    print("\n=== Argument ② (cross-chr 偏 HP1) Summary ===")
    print(f"{'Source':<20} {'genome HP1:HP2':>16} {'chr ratio>5':>12} {'chr ratio<0.2':>14} {'chr balanced':>14}")
    summary_rows = []
    for name, data in all_data.items():
        genome_h1 = sum(v[0] for v in data.values())
        genome_h2 = sum(v[1] for v in data.values())
        genome_ratio = genome_h1 / max(genome_h2, 1)
        chr_ratios = []
        for chrom in CHRS_ORDER:
            h1, h2 = data.get(chrom, (0, 0))
            if h1 + h2 > 100:  # min coverage threshold
                chr_ratios.append((chrom, h1 / max(h2, 1)))
        n_biased_hp1 = sum(1 for _, r in chr_ratios if r > 5.0)
        n_biased_hp2 = sum(1 for _, r in chr_ratios if r < 0.2)
        n_balanced = sum(1 for _, r in chr_ratios if 0.5 <= r <= 2.0)
        print(
            f"{name:<20} {f'{genome_h1}:{genome_h2}':>10}={genome_ratio:>5.2f} "
            f"{n_biased_hp1:>11}/{len(chr_ratios)} {n_biased_hp2:>13}/{len(chr_ratios)} "
            f"{n_balanced:>11}/{len(chr_ratios)}"
        )
        summary_rows.append({
            "source": name,
            "genome_hp1": genome_h1,
            "genome_hp2": genome_h2,
            "genome_ratio": genome_ratio,
            "n_chr_total": len(chr_ratios),
            "n_chr_biased_hp1_gt5": n_biased_hp1,
            "n_chr_biased_hp2_lt0_2": n_biased_hp2,
            "n_chr_balanced_05_to_2": n_balanced,
        })

    # Save TSV
    out_tsv = REPO + "/research/paired_priority_bug_audit/v6_per_chr_alt_ratio_summary.tsv"
    with open(out_tsv, "w", newline="") as fh:
        if summary_rows:
            writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)
    print(f"\nSummary → {out_tsv}")

    # Per-chr detail TSV
    detail_tsv = REPO + "/research/paired_priority_bug_audit/v6_per_chr_alt_ratio_detail.tsv"
    with open(detail_tsv, "w", newline="") as fh:
        fields = ["source", "chr", "hp1_count", "hp2_count", "ratio"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for name, data in all_data.items():
            for chrom in CHRS_ORDER:
                h1, h2 = data.get(chrom, (0, 0))
                if h1 + h2 == 0:
                    continue
                writer.writerow({
                    "source": name,
                    "chr": chrom,
                    "hp1_count": h1,
                    "hp2_count": h2,
                    "ratio": h1 / max(h2, 1),
                })
    print(f"Detail → {detail_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
