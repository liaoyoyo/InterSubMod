#!/usr/bin/env python3
"""Quantify how much of the observed read-AF is explainable by copy number alone.

Motivation
----------
The frozen read-AF score orders candidate edges by "ancestor AF >= descendant
AF".  That monotonicity argument assumes a copy-neutral, single-copy mutation
model.  Under gain/LOH the *same* clone can produce a low AF purely from
mutation multiplicity, which would be misread as a smaller (later) subclone.

This script does NOT recompute CCF or re-rank any tree.  It only measures how
close the observed HP-family-specific read-AF sits to the discrete grid that a
pure copy-number explanation would predict, so the audit can state how much of
the AF signal is confounded rather than informative.

Model used for the expected grid
--------------------------------
read-AF here is family-specific: ALT/(REF+ALT) within one germline HP family.
Let c be the number of copies of that HP present, and m the number of those
copies carrying the mutation (1 <= m <= c).  Then

    E[read-AF] = m / c

SEQC2 truth supplies total CN and an LOH flag, but not allele-specific CN, so c
is derived under two explicitly labelled assumptions:

  * LOH segment      -> all copies sit on the retained haplotype:  c = total_CN
  * non-LOH segment  -> copies split as evenly as possible:        c = round(total_CN / 2)

Both are approximations and are reported as such.  A site is called
"CN-explainable" when the observed AF lies within `TOL` of some m/c grid point.
Because a coarse grid is easy to hit by chance, the script also reports the
grid-hit rate expected from AF values drawn uniformly at random, giving an
excess over chance rather than a bare percentage.
"""

from __future__ import annotations

import json
import os
from collections import Counter, defaultdict
from fractions import Fraction

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
SITES = os.path.join(DATA, "hcc1395_site_cn_annotation.jsonl")
OUT = os.path.join(DATA, "af_vs_cn_expectation.json")

TOL = 0.05  # absolute AF tolerance for calling a grid hit
MIN_DEPTH = 10  # minimum ALT+REF reads for the AF to be worth interpreting


def hp_copies(total_cn, is_loh):
    """Copies of the observed HP family, under the stated assumptions."""
    if total_cn is None:
        return None
    if is_loh:
        c = int(round(total_cn))
    else:
        c = int(round(total_cn / 2.0))
    return max(c, 1)


def grid_points(c):
    return [m / c for m in range(1, c + 1)]


def nearest_grid(af, c):
    pts = grid_points(c)
    best = min(pts, key=lambda p: abs(p - af))
    return best, abs(best - af)


def chance_hit_rate(c, tol=TOL):
    """P(a uniform-random AF in (0,1] lands within tol of some m/c point)."""
    covered = 0.0
    for m in range(1, c + 1):
        p = m / c
        lo, hi = max(0.0, p - tol), min(1.0, p + tol)
        covered += hi - lo
    # grid spacing is 1/c; for c large the windows merge, so cap at 1
    return min(covered, 1.0)


def main():
    per_class = defaultdict(list)
    per_cn = defaultdict(list)
    rows = 0
    usable = 0

    hits = Counter()
    total = Counter()
    chance_sum = defaultdict(float)
    dist_by_class = defaultdict(list)
    af_near_third = 0
    af_near_half = 0
    af_near_one = 0

    with open(SITES) as fh:
        for ln in fh:
            r = json.loads(ln)
            rows += 1
            af = r.get("read_af")
            a, ref = r.get("alt_reads"), r.get("ref_reads")
            if af is None or a is None or ref is None:
                continue
            depth = a + ref
            if depth < MIN_DEPTH:
                continue
            usable += 1

            cls = r["seqc2_class"]
            cn = r["seqc2_total_cn"]
            is_loh = bool(r["seqc2_loh"])
            c = hp_copies(cn, is_loh)

            per_class[cls].append(af)
            per_cn[cn].append(af)

            if c:
                pt, dist = nearest_grid(af, c)
                total[cls] += 1
                dist_by_class[cls].append(dist)
                chance_sum[cls] += chance_hit_rate(c)
                if dist <= TOL:
                    hits[cls] += 1

            if abs(af - 1 / 3) <= TOL:
                af_near_third += 1
            if abs(af - 0.5) <= TOL:
                af_near_half += 1
            if abs(af - 1.0) <= TOL:
                af_near_one += 1

    def stats(v):
        if not v:
            return None
        s = sorted(v)
        n = len(s)
        return {
            "n": n,
            "median": round(s[n // 2], 4),
            "q1": round(s[n // 4], 4),
            "q3": round(s[(3 * n) // 4], 4),
            "mean": round(sum(s) / n, 4),
        }

    explain = {}
    for cls in sorted(total):
        n = total[cls]
        h = hits[cls]
        exp_chance = chance_sum[cls] / n if n else 0.0
        d = sorted(dist_by_class[cls])
        explain[cls] = {
            "n_sites": n,
            "grid_hits": h,
            "grid_hit_rate_percent": round(100.0 * h / n, 2) if n else None,
            "chance_hit_rate_percent": round(100.0 * exp_chance, 2),
            "excess_over_chance_percent": round(100.0 * (h / n - exp_chance), 2)
            if n
            else None,
            "median_distance_to_grid": round(d[len(d) // 2], 4) if d else None,
        }

    allcls_n = sum(total.values())
    allcls_h = sum(hits.values())
    allcls_chance = sum(chance_sum.values()) / allcls_n if allcls_n else 0.0

    out = {
        "generated_by": os.path.basename(__file__),
        "model": {
            "read_af_definition": "family-specific ALT/(REF+ALT) within one germline HP family",
            "expected_grid": "E[read-AF] = m/c, m in 1..c",
            "c_assumption_loh": "c = total_CN (all copies on retained haplotype)",
            "c_assumption_non_loh": "c = round(total_CN / 2) (even split; allele-specific CN unavailable)",
            "tolerance": TOL,
            "min_depth": MIN_DEPTH,
            "caveat": "SEQC2 supplies total CN + LOH only; allele-specific CN is assumed, not measured. Grid-hit rate is model-conditional evidence of CN explainability, not proof that no subclone exists.",
        },
        "counts": {
            "site_rows_read": rows,
            "sites_with_af_and_depth_ge_min": usable,
        },
        "af_by_cn_class": {k: stats(v) for k, v in sorted(per_class.items())},
        "af_by_total_cn": {str(k): stats(v) for k, v in sorted(per_cn.items())},
        "cn_explainability": explain,
        "cn_explainability_all": {
            "n_sites": allcls_n,
            "grid_hits": allcls_h,
            "grid_hit_rate_percent": round(100.0 * allcls_h / allcls_n, 2)
            if allcls_n
            else None,
            "chance_hit_rate_percent": round(100.0 * allcls_chance, 2),
            "excess_over_chance_percent": round(
                100.0 * (allcls_h / allcls_n - allcls_chance), 2
            )
            if allcls_n
            else None,
        },
        "af_landmark_counts": {
            "within_tol_of_one_third": af_near_third,
            "within_tol_of_one_half": af_near_half,
            "within_tol_of_one": af_near_one,
            "denominator": usable,
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"rows={rows} usable(depth>={MIN_DEPTH})={usable}")
    print(f"wrote {OUT}")
    for cls, v in explain.items():
        print(
            f"  {cls:10s} n={v['n_sites']:6d} hit={v['grid_hit_rate_percent']:6.2f}% "
            f"chance={v['chance_hit_rate_percent']:6.2f}% "
            f"excess={v['excess_over_chance_percent']:+6.2f}%"
        )


if __name__ == "__main__":
    main()
