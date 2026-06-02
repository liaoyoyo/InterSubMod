#!/usr/bin/env python3
"""
BRCA2 canonical regression test — freezes the golden cis-ASM metrics so any future
change to cis_asm_core is caught by git-diff. Self-contained (no pytest needed):

    python3 pipeline/tests/test_brca2_regression.py    # -> PASS/FAIL + exit code

Also pytest-compatible (def test_*). Golden values are L1-verified (grep-able in
genome_survey_v2/feature_matrix_trial.json + control_cohesion_cistest.json), NOT
fabricated. Tolerances are non-zero only to absorb float rounding.

If a golden value legitimately changes (method improvement), update GOLDEN here AND
commit with message tagged [golden-update] + rationale.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from lib.cis_asm_core import load_level1, dbeta_axis, cis_test, cohesion_per_tag  # noqa: E402

FIXTURE = os.path.join(HERE, "..", "fixtures", "brca2_level1.tsv.gz")
SPOS = "32315128"  # BRCA2/ZAR1L chr13:32,315,128

# golden: (value, tolerance) — frozen 2026-06-03, source feature_matrix_trial.json
GOLDEN = {
    "dbeta_HP":   (-0.122, 0.005),   # HP1 vs HP1-1 paired any-mod Δβ (cross-extractor validated)
    "d_cis":      (-0.1418, 0.005),  # tumor HP1-1 vs normal HP1 (somatic vs baseline)
    "d_drift":    (-0.0217, 0.005),  # tumor HP1 vs normal HP1 (germline drift, small)
    "d_somatic":  (-0.1201, 0.005),  # tumor HP1-1 vs tumor HP1 (within-tumor ASM)
    "p_cis":      (0.0,     0.001),  # cis significant
    "sil_HP11":   (0.313,   0.010),  # somatic-subclone cohesion (most cohesive)
    "sil_HP1":    (0.119,   0.010),  # germline cohesion
}


def _compute():
    D = load_level1(FIXTURE)
    coh = cohesion_per_tag(D, SPOS)
    ct = cis_test(D, SPOS)
    return {
        "dbeta_HP": dbeta_axis(D, SPOS, {"1"}, {"1-1"}),
        "d_cis": ct.get("d_cis"), "d_drift": ct.get("d_drift"),
        "d_somatic": ct.get("d_somatic"), "p_cis": ct.get("p_cis"),
        "sil_HP11": coh.get("1-1"), "sil_HP1": coh.get("1"),
    }


def run():
    got = _compute()
    failed = 0
    print(f"BRCA2 regression ({SPOS}) — fixture {os.path.basename(FIXTURE)}")
    for name, (exp, tol) in GOLDEN.items():
        v = got.get(name)
        ok = v is not None and abs(v - exp) <= tol
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:11s} got={v}  golden={exp}  (tol {tol})")
        if not ok:
            failed += 1
    print(f"\n{'== ALL PASS ==' if failed == 0 else '== %d FAILED ==' % failed}  ({len(GOLDEN)} checks)")
    return failed


def test_brca2_golden():
    """pytest entry point."""
    assert run() == 0


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
