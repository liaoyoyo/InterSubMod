#!/usr/bin/env python3
"""
P1 BAM-direct extractor regression — BAM-dependent (SKIPs cleanly if BAM absent,
so the hermetic Level-1 test stays the portable golden).

Confirms:
  - 5mC / 5hmC separated at source (mod_code column), window-restricted (PoC-36 fix)
  - cross-extractor INVARIANT: dbeta_mod(any) on BAM == MSA Level-1 dbeta_HP (-0.122)
  - BRCA2 ASM is pure 5mC: dbeta_5mC=-0.121, dbeta_5hmC=-0.001
  - cache HIT on 2nd call (no BAM re-walk)

    python3 pipeline/tests/test_brca2_bam_extract.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from stages.stage_extract import extract_or_cache  # noqa: E402
from lib.cis_asm_core import load_level1_plus, dbeta_mod  # noqa: E402

# TODO(P2): move BAM paths to loci.yaml / a config module
TUMOR = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
         "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam")
NORMAL = "/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam"
CACHE = os.path.join(HERE, "..", "cache", "level1")  # gitignored (research/**/cache/**)

SPOS = 32315128  # BRCA2/ZAR1L
GOLDEN = {
    "dbeta_5mC":  (-0.121, 0.005),   # ASM is pure 5mC
    "dbeta_5hmC": (-0.001, 0.005),   # 5hmC negligible at BRCA2
    "dbeta_any":  (-0.122, 0.005),   # == MSA Level-1 dbeta_HP (cross-extractor invariant)
}


def run():
    if not (os.path.exists(TUMOR) and os.path.exists(NORMAL)):
        print("SKIP — BAM not available (hermetic Level-1 test in test_brca2_regression.py covers the core)")
        return 0
    path, st1 = extract_or_cache("chr13", SPOS, "G", "A", TUMOR, NORMAL, cache_dir=CACHE)
    Dp = load_level1_plus(path)
    s = str(SPOS)
    got = {
        "dbeta_5mC": dbeta_mod(Dp, s, "m"),
        "dbeta_5hmC": dbeta_mod(Dp, s, "h"),
        "dbeta_any": dbeta_mod(Dp, s, "any"),
    }
    _, st2 = extract_or_cache("chr13", SPOS, "G", "A", TUMOR, NORMAL, cache_dir=CACHE)

    failed = 0
    print(f"P1 BAM-extract regression (chr13:{SPOS})  extract#1={st1}  extract#2={st2}")
    for name, (exp, tol) in GOLDEN.items():
        v = got[name]
        ok = v is not None and abs(v - exp) <= tol
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:11s} got={v}  golden={exp}  (tol {tol})")
        if not ok:
            failed += 1
    cache_ok = (st2 == "HIT")
    print(f"  [{'PASS' if cache_ok else 'FAIL'}] cache HIT on 2nd call (st2={st2})")
    if not cache_ok:
        failed += 1
    print(f"\n{'== ALL PASS ==' if failed == 0 else '== %d FAILED ==' % failed}")
    return failed


def test_brca2_bam_extract():
    assert run() == 0


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
