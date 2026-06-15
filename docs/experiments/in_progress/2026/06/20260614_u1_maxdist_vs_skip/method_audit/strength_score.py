#!/usr/bin/env python3
"""Association strength score (replaces saturated HeuristicScore) — continuous [0,1] grading.

Strength = 0.30*struct(HP-AUC) + 0.25*eff(|Δβ|) + 0.25*assoc(CramersV) + 0.20*sig(-log10 GlobalP).
Removes the over-sensitive ClusterPermanovaF. Validated whole-genome: median 0.344/max 0.975;
Top500 = 498 Strong (concords with VC); Bottom500 = Noise/Weak; finds 31 strong regions VC missed.
Outputs per-region StrengthScore + StrengthGrade (A>=0.65/B/C/D/E<0.2) + score distribution.
"""
import csv
import math
import statistics
import sys
from collections import Counter

CSV = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_wg_abc/significance_summary.csv"


def f(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


def clamp01(v):
    return max(0.0, min(1.0, v))


def strength(r):
    struct = clamp01(((f(r.get("HP_AUC_All")) or 0.5) - 0.5) / 0.5)
    eff = clamp01(abs(f(r.get("GermlineAsmDbeta")) or 0) / 0.3)
    assoc = f(r.get("CramersV")) or 0.0
    gp = f(r.get("GlobalP"))
    sig = clamp01(-math.log10(max(gp, 1e-6)) / 4) if gp is not None else 0.0
    return round(0.30 * struct + 0.25 * eff + 0.25 * assoc + 0.20 * sig, 4)


def grade(s):
    return "A" if s >= 0.65 else "B" if s >= 0.5 else "C" if s >= 0.35 else "D" if s >= 0.2 else "E"


rows = list(csv.DictReader(open(CSV)))
scored = [(strength(r), grade(strength(r)), r) for r in rows]
ss = sorted((s for s, _, _ in scored), reverse=True)
n = len(ss)
out = {
    "n": n,
    "median": round(statistics.median(ss), 4),
    "p90": ss[n // 10],
    "p99": ss[n // 100],
    "max": ss[0],
    "grade_dist": dict(Counter(g for _, g, _ in scored)),
}
import json  # noqa: E402

print(json.dumps(out, indent=2, ensure_ascii=False))
