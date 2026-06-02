#!/usr/bin/env python3
"""
P3a causation regression — mechanical-cis (CpG create/destroy). Needs ref FASTA on /big8
(SKIPs cleanly if absent). Golden values cross-checked vs feature_matrix_trial.json.

    python3 pipeline/tests/test_causation.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from lib.causation import RefContext, causation_summary  # noqa: E402

# golden: (tag, chrom, pos, ref, alt, expected_mechanical)
GOLDEN = [
    ("BRCA2",     "chr13", 32315128, "G", "A", "NEUTRAL"),       # regulatory cis, not mechanical
    ("het_DYTN",  "chr2",  206705390, "G", "A", "DESTROYS_CpG"),  # germline het destroying a CpG
    ("FP_chr9",   "chr9",  41788031, "A", "T", "CREATES_CpG"),    # FP creating a CpG
]


def run():
    rc = RefContext()
    if not rc.available:
        print("SKIP — reference FASTA not available")
        return 0
    failed = 0
    print("P3a causation regression (mechanical-cis CpG create/destroy)")
    for tag, chrom, pos, ref, alt, exp in GOLDEN:
        got = rc.mechanical_cis(chrom, pos, ref, alt)
        ok = (got == exp)
        print(f"  [{'PASS' if ok else 'FAIL'}] {tag:10s} {ref}>{alt}  got={got}  golden={exp}")
        if not ok:
            failed += 1
    # BRCA2 cis-candidate -> regulatory (not mechanical) causation summary
    cs = causation_summary("NEUTRAL", "T3")
    ok = "regulatory" in cs["causation_mechanism"]
    print(f"  [{'PASS' if ok else 'FAIL'}] BRCA2 T3 causation = {cs['causation_mechanism']}")
    if not ok:
        failed += 1
    print(f"\n{'== ALL PASS ==' if failed == 0 else '== %d FAILED ==' % failed}")
    return failed


def test_causation():
    assert run() == 0


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
