#!/usr/bin/env python3
"""V5 + force highPurity=false (path 2 only) HP metrics — 15-site sample.

Compares with existing 5 versions @ 0.93 to evaluate whether path 3 抵消 path 2 反轉.
"""
from __future__ import annotations
import csv, os, sys
from typing import Dict, List, Tuple
import pysam

OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260501_v5_force_path2only"
os.makedirs(OUT_DIR, exist_ok=True)

BAMS: Dict[str, str] = {
    "A1_baseline_old": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam",
    "A2_v2b_old": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_tagged.bam",
    "A5_v5_old": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam",
    "NEW_baseline_two_pass": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/baseline_09/tumor_tagged.bam",
    "NEW_v5_flag_three_pass": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_tagged.bam",
    "NEW_v5_flag_force_path2only": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v5_flag_force_path2only/tumor_tagged.bam",
}

SITES: List[Tuple[str, str, int, str]] = [
    ("A_TP01", "chr6", 145444893, "Phase4 TP_01"),
    ("A_TP02", "chr4", 70548355, "Phase4 TP_02"),
    ("A_TP03", "chr5", 153209947, "Phase4 TP_03"),
    ("A_TP04", "chr16", 35118902, "Phase4 TP_04"),
    ("A_TP05", "chr7", 109185781, "Phase4 TP_05"),
    ("B_FPA1", "chr8", 93565727, "Phase4 FPA1"),
    ("B_FPA2", "chr9", 137953060, "Phase4 FPA2"),
    ("B_FPB1", "chr7", 52087777, "Phase4 FPB1"),
    ("B_FPB2", "chr9", 75383880, "Phase4 FPB2"),
    ("C_FPSF1", "chr2", 16800066, "Self-FP_01"),
    ("C_FPSF2", "chr11", 65528013, "Self-FP_02"),
    ("D_GE1", "chr19", 11220247, "Phase3 LOH GE_01"),
    ("D_GE2", "chr5", 35850791, "Phase3 LOH GE_02"),
    ("E_AMB1", "chr1", 156785432, "Self-AMB_01"),
    ("E_AMB2", "chr12", 89213411, "Self-AMB_02"),
]


def get_hp_counts(bam: pysam.AlignmentFile, chrom: str, pos: int) -> Dict[int, int]:
    counts = {0: 0, 1: 0, 2: 0, 11: 0, 21: 0, 33: 0}
    for read in bam.fetch(chrom, pos - 1, pos):
        if read.is_secondary or read.is_unmapped:
            continue
        try:
            hp = int(read.get_tag("HP"))
        except KeyError:
            hp = 0
        counts[hp] = counts.get(hp, 0) + 1
    return counts


def main():
    summary_rows = []
    for version, bam_path in BAMS.items():
        if not os.path.exists(bam_path):
            print(f"WARN: {version} not found", file=sys.stderr)
            continue
        print(f"Processing {version}...", file=sys.stderr)
        agg = {0: 0, 1: 0, 2: 0, 11: 0, 21: 0, 33: 0}
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for _, chrom, pos, _ in SITES:
                counts = get_hp_counts(bam, chrom, pos)
                for k, v in counts.items():
                    agg[k] += v
        hp1f = agg[1] + agg[11]
        hp2f = agg[2] + agg[21]
        amb = agg[33]
        total = sum(agg.values())
        ratio = float('inf') if hp2f == 0 else hp1f / hp2f
        amb_pct = (amb / (hp1f + hp2f + amb) * 100) if (hp1f + hp2f + amb) > 0 else 0.0
        summary_rows.append({
            "version": version,
            "HP_0": agg[0], "HP_1": agg[1], "HP_2": agg[2],
            "HP_11": agg[11], "HP_21": agg[21], "HP_33": agg[33],
            "HP1_family": hp1f, "HP2_family": hp2f,
            "ratio_HP1:HP2": f"{ratio:.3f}",
            "AMB_pct": f"{amb_pct:.2f}",
            "total_reads": total,
        })

    out_path = os.path.join(OUT_DIR, "summary_ratio_path2only.tsv")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\nWrote {out_path}")
    print(f"\n{'version':<32} {'HP1_family':>10} {'HP2_family':>10} {'HP_33':>7} {'ratio':>8} {'AMB%':>6}")
    for r in summary_rows:
        print(f"{r['version']:<32} {r['HP1_family']:>10} {r['HP2_family']:>10} {r['HP_33']:>7} {r['ratio_HP1:HP2']:>8} {r['AMB_pct']:>6}")


if __name__ == "__main__":
    main()
