#!/usr/bin/env python3
"""V3-Fixed getVote-only Ablation Metrics — 5-BAM comparison.

Compares 5 versions (A1-A5) on 15 cherry-picked sites + summary stats:
- A1 baseline
- A2 V2b (PON-only ON, original getVote)
- A3 v3f_no_pononly (PON-only OFF, V3-Fixed getVote) ← target
- A4 V3-Fixed (PON-only ON, V3-Fixed getVote)
- A5 V5 (PON-only ON, V5 Layer 1.5)

Outputs:
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/per_site_hp_counts.tsv
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/summary_ratio.tsv
"""

from __future__ import annotations

import csv
import os
import sys
from typing import Dict, List, Tuple

try:
    import pysam
except ImportError:
    print("ERROR: pysam not found. Install via: pip install pysam", file=sys.stderr)
    sys.exit(1)

OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation"
os.makedirs(OUT_DIR, exist_ok=True)

BAMS: Dict[str, str] = {
    "A1_baseline": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam",
    "A2_v2b": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_tagged.bam",
    "A3_v3f_no_pononly": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/tumor_tagged.bam",
    "A4_v3f_with_pononly": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam",
    "A5_v5": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam",
}

# 15 cherry-picked sites (from v5_imbalance_improvement.py)
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
    ("C_V5max1", "chr19", 4639528, "V5 reassign 1"),
    ("C_V5max2", "chr19", 2235521, "V5 reassign 2"),
    ("C_V5max3", "chr19", 7405500, "V5 reassign 3"),
    ("D_SP1", "chr19", 17565944, "Self-phasing extreme 1"),
    ("D_SP2", "chr19", 12452332, "Self-phasing extreme 2"),
    ("D_SP3", "chr19", 12467180, "Self-phasing extreme 3"),
]


def get_hp_counts(bam: pysam.AlignmentFile, chrom: str, pos: int) -> Dict[int, int]:
    """Count HP tag distribution for reads overlapping site (1-based pos)."""
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
    per_site_rows = []
    summary_rows = []

    for version, bam_path in BAMS.items():
        if not os.path.exists(bam_path):
            print(f"WARN: {version} BAM not found at {bam_path}", file=sys.stderr)
            continue

        print(f"Processing {version}...", file=sys.stderr)

        # Aggregate counts across all 15 sites
        agg = {0: 0, 1: 0, 2: 0, 11: 0, 21: 0, 33: 0}

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for case, chrom, pos, _ in SITES:
                counts = get_hp_counts(bam, chrom, pos)
                hp1_family = counts[1] + counts[11]
                hp2_family = counts[2] + counts[21]
                amb = counts[33]
                untagged = counts[0]
                total = hp1_family + hp2_family + amb + untagged

                ratio_str = "inf" if hp2_family == 0 else f"{hp1_family/hp2_family:.2f}"
                amb_pct = (amb / (hp1_family + hp2_family + amb) * 100) if (hp1_family + hp2_family + amb) > 0 else 0.0

                per_site_rows.append({
                    "version": version,
                    "case": case,
                    "chrom": chrom,
                    "pos": pos,
                    "HP_0": counts[0],
                    "HP_1": counts[1],
                    "HP_2": counts[2],
                    "HP_11": counts[11],
                    "HP_21": counts[21],
                    "HP_33": counts[33],
                    "HP1_family": hp1_family,
                    "HP2_family": hp2_family,
                    "ratio_HP1:HP2": ratio_str,
                    "AMB_pct": f"{amb_pct:.2f}",
                    "total_reads": total,
                })

                for k, v in counts.items():
                    agg[k] += v

        # Summary per version (15-site aggregate)
        hp1f_agg = agg[1] + agg[11]
        hp2f_agg = agg[2] + agg[21]
        amb_agg = agg[33]
        total_agg = sum(agg.values())
        ratio_agg = "inf" if hp2f_agg == 0 else f"{hp1f_agg/hp2f_agg:.3f}"
        amb_pct_agg = (amb_agg / (hp1f_agg + hp2f_agg + amb_agg) * 100) if (hp1f_agg + hp2f_agg + amb_agg) > 0 else 0.0

        summary_rows.append({
            "version": version,
            "HP_0": agg[0],
            "HP_1": agg[1],
            "HP_2": agg[2],
            "HP_11": agg[11],
            "HP_21": agg[21],
            "HP_33": agg[33],
            "HP1_family": hp1f_agg,
            "HP2_family": hp2f_agg,
            "ratio_HP1:HP2": ratio_agg,
            "AMB_pct": f"{amb_pct_agg:.2f}",
            "total_reads": total_agg,
        })

    # Write outputs
    per_site_path = os.path.join(OUT_DIR, "per_site_hp_counts.tsv")
    summary_path = os.path.join(OUT_DIR, "summary_ratio.tsv")

    if per_site_rows:
        with open(per_site_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(per_site_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(per_site_rows)
        print(f"Wrote {per_site_path} ({len(per_site_rows)} rows)", file=sys.stderr)

    if summary_rows:
        with open(summary_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(summary_rows)
        print(f"Wrote {summary_path} ({len(summary_rows)} rows)", file=sys.stderr)

        print("\n=== 15-site aggregate summary ===", file=sys.stderr)
        print(f"{'version':<25} {'HP1_fam':>10} {'HP2_fam':>10} {'HP_33':>8} {'ratio':>8} {'AMB%':>6}", file=sys.stderr)
        for r in summary_rows:
            print(f"{r['version']:<25} {r['HP1_family']:>10} {r['HP2_family']:>10} {r['HP_33']:>8} {r['ratio_HP1:HP2']:>8} {r['AMB_pct']:>6}", file=sys.stderr)


if __name__ == "__main__":
    main()
