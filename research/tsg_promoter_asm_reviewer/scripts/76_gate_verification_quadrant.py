#!/usr/bin/env python3
"""
76 - Gate verification: is each ISM "Significant" call a GOOD (real ASM) or a
threshold-passing artifact? And are there GOOD signals that MISSED the gate?

Uses ONLY independent signals already in significance_summary (post-hoc, no re-run,
no circular use of the Significant flag itself):
  - cross-test robustness: how many of 4 independent tests are p<=0.05
      (ClusterPermanovaP, LabelHPPermanovaP, LabelAllelePermanovaP, HPFineP)
  - dispersion-artifact flag: ClusterDispersionWarn / LabelHPDispersionWarn
      (PERMANOVA can be significant due to group-DISPERSION not location -> FP detector)
  - reproducibility: Stability (bootstrap)
  - validity: *Valid flags
  - effect size: |HPMergedDelta|

4 quadrants:
  A real-good      : Significant + robust(>=3/4) + no dispersion warn + stable
  B suspect-FP     : Significant but fragile(1/4) OR dispersion warn OR low stability
  C missed-good    : NOT Significant but strong (|HPMergedDelta|>=0.1 + >=2 tests + no warn
                     + CramersV in [0.05,0.1) -> just below the cutoff)
  D correct-reject : NOT Significant + weak

USAGE: python3 76_gate_verification_quadrant.py <sample> <class>   (default HCC1395 tp)
"""
import sys, csv
import numpy as np

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
sample = sys.argv[1] if len(sys.argv) > 1 else "HCC1395"
cls = sys.argv[2] if len(sys.argv) > 2 else "tp"
path = f"{EX}/{sample}_{cls}/significance_summary.csv"
rows = list(csv.DictReader(open(path)))


def f(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


TESTS = ["ClusterPermanovaP", "LabelHPPermanovaP", "LabelAllelePermanovaP", "HPFineP"]


def n_tests_sig(r):
    return sum(1 for t in TESTS if (f(r, t) is not None and f(r, t) <= 0.05))


def has_dispersion_warn(r):
    return tb(r, "ClusterDispersionWarn") or tb(r, "LabelHPDispersionWarn") or tb(r, "LabelAlleleDispersionWarn")


def stability(r):
    return f(r, "Stability")


def cramers(r):
    # the value the gate uses (best reliable across axes); approximate with CramersV column
    return f(r, "CramersV")


n = len(rows)
sig = [r for r in rows if tb(r, "Significant")]
nonsig = [r for r in rows if not tb(r, "Significant")]

# Stability not populated in this run? check
st_pop = sum(1 for r in rows if stability(r) is not None and stability(r) > 0)
USE_STABILITY = st_pop > 0.05 * n   # only use if actually populated

# robustness/clean among Significant (valid signals: cross-test + dispersion [+ stability if populated])
A = B = MID = 0
B_reasons = {"fragile_<=1test": 0, "dispersion_warn": 0, "low_stability": 0}
for r in sig:
    nt = n_tests_sig(r)
    warn = has_dispersion_warn(r)
    st = stability(r)
    lowstab = USE_STABILITY and (st is not None and st < 0.5)
    if nt <= 1 or warn or lowstab:
        B += 1
        if nt <= 1:
            B_reasons["fragile_<=1test"] += 1
        if warn:
            B_reasons["dispersion_warn"] += 1
        if lowstab:
            B_reasons["low_stability"] += 1
    elif nt >= 3:
        A += 1
    else:  # 2/4 tests, no warn
        MID += 1

# missed-good among non-Significant: strong by OTHER signals but missed gate (CramersV/coverage)
C = 0
C_examples = []
for r in nonsig:
    d = f(r, "HPMergedDelta")
    nt = n_tests_sig(r)
    warn = has_dispersion_warn(r)
    if d is not None and abs(d) >= 0.10 and nt >= 2 and not warn:
        C += 1
        if len(C_examples) < 8:
            C_examples.append((r.get("Chr"), r.get("Pos"), round(d, 3),
                               round(cramers(r), 3) if cramers(r) is not None else None, nt))
D = len(nonsig) - C
print(f"[note] Stability column populated: {st_pop}/{n} (USE_STABILITY={USE_STABILITY}); "
      f"若需 bootstrap stability 訊號要重跑 ISM 開 bootstrap")

print(f"=== Gate verification quadrant: {sample} {cls} (n={n}) ===")
print(f"  Significant={len(sig)}  non-Significant={len(nonsig)}")
print()
print(f"  [A] real-good   (Sig + robust>=3/4 tests + no dispersion warn): {A}  ({A/len(sig)*100:.1f}% of Sig)")
print(f"  [MID] moderate  (Sig + 2/4 tests + no warn)                  : {MID}  ({MID/len(sig)*100:.1f}% of Sig)")
print(f"  [B] suspect-FP  (Sig but <=1 test OR dispersion-warn)        : {B}  ({B/len(sig)*100:.1f}% of Sig)")
print(f"        reasons: {B_reasons}")
print(f"  [C] missed-good (NOT Sig but |Δβ|>=0.1 + >=2 tests + no warn) : {C}")
print(f"  [D] correct-reject (NOT Sig + weak)                          : {D}")
print()
# distributions for context
sig_ntests = [n_tests_sig(r) for r in sig]
sig_warn = sum(1 for r in sig if has_dispersion_warn(r))
print(f"  among Significant: dispersion-warn={sig_warn} ({sig_warn/len(sig)*100:.1f}%); "
      f"n_tests_sig dist {dict(zip(*np.unique(sig_ntests, return_counts=True)))}")
st_vals = [stability(r) for r in sig if stability(r) is not None]
if st_vals:
    print(f"  Significant Stability: median={np.median(st_vals):.3f} (n_with_stability={len(st_vals)})")
if C_examples:
    print(f"  missed-good 範例 (Chr,Pos,Δβ,CramersV,n_tests): {C_examples[:6]}")
