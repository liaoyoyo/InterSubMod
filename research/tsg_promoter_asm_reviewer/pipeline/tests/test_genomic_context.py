#!/usr/bin/env python3
"""
P2 genomic-context regression — LOH + CNV integration + axis-gate.
Needs SEQC2 beds on /big8 (SKIPs cleanly if absent).

Confirms the audit fix: BRCA2 = nonLOH + integer CN=5 + HP-axis valid;
LOH loci switch primary axis to ALLELE (HP-axis structurally invalid in LOH).

    python3 pipeline/tests/test_genomic_context.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from lib.genomic_context import GenomicContext, GAIN_CN  # noqa: E402

LOH_BED = os.path.join(HERE, "..", "..", "output", "seqc2_loh_only.bed")

# golden: (chrom, pos, total_cn, loh_status, primary_axis, hp_axis_valid)
GOLDEN = [
    ("BRCA2",      "chr13", 32315128, 5.0, "nonLOH", "HP",     True),
    ("TP_LOH_chr4", "chr4", 10739757, 2.0, "LOH",    "ALLELE", False),
    ("het_chr2",    "chr2", 58616453, 3.0, "nonLOH", "HP",     True),
]


def run():
    if not (os.path.exists(GAIN_CN) and os.path.exists(LOH_BED)):
        print("SKIP — SEQC2 CN/LOH beds not available")
        return 0
    gc = GenomicContext(loh_bed=LOH_BED)
    failed = 0
    print("P2 genomic-context regression (CN + LOH + axis-gate)")
    for tag, chrom, pos, cn, loh, axis, valid in GOLDEN:
        a = gc.annotate(chrom, pos)
        ok = (a["total_cn"] == cn and a["loh_status"] == loh
              and a["primary_axis"] == axis and a["hp_axis_valid"] is valid)
        print(f"  [{'PASS' if ok else 'FAIL'}] {tag:12s} cn={a['total_cn']} loh={a['loh_status']} "
              f"axis={a['primary_axis']} hp_valid={a['hp_axis_valid']}")
        if not ok:
            failed += 1
    print(f"\n{'== ALL PASS ==' if failed == 0 else '== %d FAILED ==' % failed}")
    return failed


def test_genomic_context():
    assert run() == 0


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
