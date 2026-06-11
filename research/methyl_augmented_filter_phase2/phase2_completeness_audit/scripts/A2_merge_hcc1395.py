#!/usr/bin/env python3
"""
A2 merge: prepend HCC1395 baseline + V6 rows (from A2_scan.log first-run kill) to the
4-sample partial output, then regenerate final TSV + PNG.

HCC1395 baseline + V6 raw counts were lost when the background task was killed
mid-run; only the printed ratio summary survived. To preserve those two rows we
re-derive HP1/HP2 from the log lines and fill remaining HP bins as 'unknown'
(NaN placeholder; we still expose HP1 + HP2 + ratio which is what the report
cares about). This is acceptable because:

  - PI 4-29 cohort is defined as HP1:HP2 ratio (other HP bins are auxiliary)
  - Full HP6 breakdown for HCC1395 is captured by A1 separately

If full re-scan is needed later, run:
  python3 scripts/A2_ALT_only_HP_ratio.py --skip-samples "HCC1937/V6,HCC1954/V6,H1437/V6,H2009/V6" \
       --seed-json A2_ALT_only_HP_ratio_5sample.json \
       --out-tsv A2_ALT_only_HP_ratio_5sample.tsv \
       --out-png figures/A2_ALT_only_vs_all_reads_ratio.png \
       --out-json A2_ALT_only_HP_ratio_5sample.json
"""

import json
import sys
from pathlib import Path

# Source: A2_scan.log line capture (2026-05-21 02:36 - 03:32 CST)
# HCC1395/baseline -> n_SNV=3,575 ALT_HP1=39,732 ALT_HP2=9,014 ALT_ratio=4.4078 ALL_ratio=3.4359 elapsed=1294.4s
# HCC1395/V6       -> n_SNV=3,575 ALT_HP1=327    ALT_HP2=769   ALT_ratio=0.4252 ALL_ratio=3.3627 elapsed=1315.6s
HCC1395_ROWS = [
    {
        "sample": "HCC1395",
        "bam_version": "baseline",
        "scan_scope": "chr8+chr19",
        "n_SNV_positions": 3575,
        "n_SNV_with_coverage": None,  # not captured in log
        "ALT_reads_HP1": 39732,
        "ALT_reads_HP2": 9014,
        "ALT_reads_HP0": None,
        "ALT_reads_HP11": None,
        "ALT_reads_HP21": None,
        "ALT_reads_HP33": None,
        "ALT_reads_untagged": None,
        "ALT_reads_total": None,
        "ALL_reads_HP1": None,
        "ALL_reads_HP2": None,
        "ALL_reads_total": None,
        "ALT_only_ratio": 4.4078,
        "all_reads_ratio": 3.4359,
        "elapsed_sec": 1294.4,
        "bam_path": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam",
        "vcf_path": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased.vcf",
    },
    {
        "sample": "HCC1395",
        "bam_version": "V6",
        "scan_scope": "chr8+chr19",
        "n_SNV_positions": 3575,
        "n_SNV_with_coverage": None,
        "ALT_reads_HP1": 327,
        "ALT_reads_HP2": 769,
        "ALT_reads_HP0": None,
        "ALT_reads_HP11": None,
        "ALT_reads_HP21": None,
        "ALT_reads_HP33": None,
        "ALT_reads_untagged": None,
        "ALT_reads_total": None,
        "ALL_reads_HP1": None,
        "ALL_reads_HP2": None,
        "ALL_reads_total": None,
        "ALT_only_ratio": 0.4252,
        "all_reads_ratio": 3.3627,
        "elapsed_sec": 1315.6,
        "bam_path": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam",
        "vcf_path": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased.vcf",
    },
]

TSV_COLS = [
    "sample", "bam_version", "scan_scope", "n_SNV_positions", "n_SNV_with_coverage",
    "ALT_reads_HP1", "ALT_reads_HP2", "ALT_reads_HP0",
    "ALT_reads_HP11", "ALT_reads_HP21", "ALT_reads_HP33",
    "ALT_reads_untagged", "ALT_reads_total",
    "ALL_reads_HP1", "ALL_reads_HP2", "ALL_reads_total",
    "ALT_only_ratio", "all_reads_ratio",
    "elapsed_sec", "bam_path", "vcf_path",
]


def main():
    audit_dir = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit")
    part2_tsv = audit_dir / "A2_ALT_only_HP_ratio_5sample.tsv"
    part2_json = audit_dir / "A2_ALT_only_HP_ratio_5sample.json"
    final_tsv = audit_dir / "A2_ALT_only_HP_ratio_5sample.tsv"
    final_json = audit_dir / "A2_ALT_only_HP_ratio_5sample.json"

    if not part2_tsv.exists():
        print(f"[ERR] part2 TSV not found: {part2_tsv}", file=sys.stderr)
        return 1

    # Read part2 rows
    part2_rows = []
    with open(part2_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            vals = line.rstrip("\n").split("\t")
            d = dict(zip(header, vals))
            part2_rows.append(d)
    print(f"[merge] part2 rows = {len(part2_rows)}")

    # Final ordering: HCC1395 baseline, HCC1395 V6, HCC1937/V6, HCC1954/V6, H1437/V6, H2009/V6
    sample_order = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]
    version_order = ["baseline", "V6"]

    all_rows = HCC1395_ROWS + [
        {k: r.get(k, "") for k in TSV_COLS} for r in part2_rows
    ]
    # Sort by predefined order
    def _key(r):
        s = r["sample"]
        v = r["bam_version"]
        try:
            si = sample_order.index(s)
        except ValueError:
            si = 99
        try:
            vi = version_order.index(v)
        except ValueError:
            vi = 99
        return (si, vi)
    all_rows.sort(key=_key)

    # Write final TSV
    with open(final_tsv, "w") as f:
        f.write("\t".join(TSV_COLS) + "\n")
        for r in all_rows:
            f.write("\t".join("" if r.get(k) is None else str(r.get(k, "")) for k in TSV_COLS) + "\n")
    print(f"[merge] wrote final TSV -> {final_tsv} ({len(all_rows)} rows)")

    # Merge JSON
    if part2_json.exists():
        with open(part2_json) as f:
            part2_data = json.load(f)
        # We prepend the HCC1395 rows as 'log-derived' partial raw entries
        hcc_raw = [
            {
                "sample": r["sample"], "version": r["bam_version"],
                "bam_path": r["bam_path"], "vcf_path": r["vcf_path"],
                "scan_scope": r["scan_scope"], "n_SNV_positions": r["n_SNV_positions"],
                "alt_hp_counts": {"1": r["ALT_reads_HP1"], "2": r["ALT_reads_HP2"]},
                "ref_hp_counts": {},
                "all_hp_counts": {},
                "n_loci_processed": None,
                "n_loci_no_coverage": None,
                "n_reads_no_base": None,
                "n_reads_total": None,
                "elapsed_sec": r["elapsed_sec"],
                "log_derived": True,
                "note": "HP1+HP2 + ratio extracted from A2_scan.log line; other HP bins not captured (kill mid-run); see A1 for full HP6 breakdown",
            } for r in HCC1395_ROWS
        ]
        merged = {
            "scan_scope": part2_data.get("scan_scope", "chr8+chr19"),
            "hp_categories": part2_data.get("hp_categories"),
            "rows": hcc_raw + part2_data.get("rows", []),
            "total_elapsed_sec": part2_data.get("total_elapsed_sec"),
            "generated_at": part2_data.get("generated_at"),
            "partial": False,
            "merge_note": "HCC1395 baseline + V6 rows derived from A2_scan.log (first-run kill mid-run); 4 V6 extension rows from A2_scan_part2.log",
        }
        with open(final_json, "w") as f:
            json.dump(merged, f, indent=2)
        print(f"[merge] wrote final JSON -> {final_json} ({len(merged['rows'])} rows)")

    print("[merge] done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
