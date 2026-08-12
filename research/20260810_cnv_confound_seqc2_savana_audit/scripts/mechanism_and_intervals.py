#!/usr/bin/env python3
"""Mechanism check and interval estimates for the CN-inflated uniqueness finding.

The stratified test showed read-AF breaks ties far more often on CN-altered
ground (56.2% vs 25.5%, CMH OR ~3.5 after controlling candidate-tree count and
k).  That is an association.  This script tests the proposed *mechanism* and
puts intervals on the headline numbers.

Proposed mechanism
------------------
read-AF is family-specific ALT/(REF+ALT) within one HP family.  Where that
haplotype carries c > 1 copies, a mutation present on m of them reads out at
m/c rather than at 1.  Two mutations in the same unit can then differ in AF
purely from multiplicity, giving the monotonicity score a gradient to sort on
that has nothing to do with which mutation came first.

Testable prediction: the *within-unit spread* of read-AF should be larger on
CN-altered ground, and largest under gain+LOH, where c is largest.  If the
spread were flat across CN classes the mechanism would be refuted and the
association would need another explanation.

Also computes Wilson score intervals and an excess estimate: how many
"unique" trees would not be unique if CN-altered ground behaved like neutral
ground.
"""

from __future__ import annotations

import json
import math
import os
from collections import defaultdict

from scipy.stats import mannwhitneyu

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
ANNOT = os.path.join(DATA, "hcc1395_unit_cn_annotation.jsonl")
TOPO = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
    "samples/HCC1395/HCC1395.topology.jsonl"
)
OUT = os.path.join(DATA, "mechanism_and_intervals.json")

MIN_DEPTH = 10


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {
        "point_percent": round(100 * p, 2),
        "lo_percent": round(100 * max(0.0, centre - half), 2),
        "hi_percent": round(100 * min(1.0, centre + half), 2),
        "n": n,
    }


def main():
    cn_by_region = {}
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") == "ANNOTATED":
                cn_by_region[r["region_id"]] = r["seqc2_class"]

    spread_by_class = defaultdict(list)
    spread_by_altered = defaultdict(list)
    contested_flag = {}
    rows = []

    with open(TOPO) as fh:
        for ln in fh:
            u = json.loads(ln)
            rid = u.get("region_id")
            cls = cn_by_region.get(rid)
            if cls is None or u.get("unit_status") != "ranked":
                continue
            ttc = int(u.get("total_tree_count") or 0)
            # Restrict to active sites: inactive af_coverage rows are all-reference
            # and are not seen by the score, so including them would inflate spread.
            active = set(u.get("active_positions") or [])
            afs = []
            for c in u.get("af_coverage") or []:
                if c.get("position") not in active:
                    continue
                a, r_ = c.get("alt_reads"), c.get("ref_reads")
                if a is None or r_ is None or (a + r_) < MIN_DEPTH:
                    continue
                afs.append(a / (a + r_))
            rows.append(
                {
                    "cn_class": cls,
                    "altered": cls != "neutral",
                    "total_tree_count": ttc,
                    "best_tree_unique": bool(u.get("best_tree_unique")),
                    "af_broke_tie": ttc > 1 and bool(u.get("best_tree_unique")),
                    "n_af": len(afs),
                }
            )
            if len(afs) >= 2:
                spread = max(afs) - min(afs)
                spread_by_class[cls].append(spread)
                spread_by_altered["altered" if cls != "neutral" else "neutral"].append(
                    spread
                )

    def desc(v):
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
            "fraction_above_0.2": round(
                sum(1 for x in s if x > 0.2) / n, 4
            ),
        }

    spread_out = {k: desc(v) for k, v in sorted(spread_by_class.items())}

    mech_test = None
    if spread_by_altered.get("altered") and spread_by_altered.get("neutral"):
        stat, p = mannwhitneyu(
            spread_by_altered["altered"],
            spread_by_altered["neutral"],
            alternative="greater",
        )
        mech_test = {
            "test": "Mann-Whitney U, altered > neutral, one-sided",
            "n_altered": len(spread_by_altered["altered"]),
            "n_neutral": len(spread_by_altered["neutral"]),
            "median_altered": desc(spread_by_altered["altered"])["median"],
            "median_neutral": desc(spread_by_altered["neutral"])["median"],
            "u_statistic": float(stat),
            "p_value": float(p),
        }

    # ---- intervals on the contested-layer rates ----
    contested = [r for r in rows if r["total_tree_count"] > 1]
    neu = [r for r in contested if not r["altered"]]
    alt = [r for r in contested if r["altered"]]
    neu_k = sum(1 for r in neu if r["af_broke_tie"])
    alt_k = sum(1 for r in alt if r["af_broke_tie"])

    ci_neu = wilson(neu_k, len(neu))
    ci_alt = wilson(alt_k, len(alt))

    # ---- excess estimate: how many "unique" trees are CN-attributable ----
    excess = None
    if ci_neu and alt:
        point = alt_k - len(alt) * (ci_neu["point_percent"] / 100.0)
        lo = alt_k - len(alt) * (ci_neu["hi_percent"] / 100.0)
        hi = alt_k - len(alt) * (ci_neu["lo_percent"] / 100.0)
        total_unique = sum(1 for r in rows if r["best_tree_unique"])
        excess = {
            "logic": "expected af_broke_tie count on CN-altered ground if it behaved like CN-neutral ground, subtracted from observed",
            "observed_af_broke_tie_on_altered": alt_k,
            "expected_at_neutral_rate": round(
                len(alt) * (ci_neu["point_percent"] / 100.0), 1
            ),
            "excess_point": round(point, 1),
            "excess_lo": round(lo, 1),
            "excess_hi": round(hi, 1),
            "total_unique_best_tree_units": total_unique,
            "excess_as_percent_of_unique_units": {
                "point": round(100 * point / total_unique, 2),
                "lo": round(100 * lo / total_unique, 2),
                "hi": round(100 * hi / total_unique, 2),
            },
            "caveat": "Assumes the CN-neutral rate is the CN-free baseline. The neutral contested sample is small (n=%d), so the interval is wide; and neutral units are not guaranteed free of unmodelled subclonal CN."
            % len(neu),
        }

    out = {
        "generated_by": os.path.basename(__file__),
        "mechanism": {
            "prediction": "within-unit read-AF spread larger on CN-altered ground, largest under gain+LOH",
            "within_unit_af_spread_by_cn_class": spread_out,
            "test_altered_vs_neutral": mech_test,
        },
        "contested_layer_rates_with_intervals": {
            "neutral": ci_neu,
            "cn_altered": ci_alt,
        },
        "cn_attributable_uniqueness_excess": excess,
        "min_depth_for_af": MIN_DEPTH,
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print("within-unit read-AF spread by CN class:")
    for k, v in spread_out.items():
        if v:
            print(
                f"  {k:10s} n={v['n']:5d} median={v['median']:.4f} "
                f"q3={v['q3']:.4f} frac>0.2={v['fraction_above_0.2']:.3f}"
            )
    print(f"\nmechanism test: {mech_test}")
    print(f"\nneutral contested rate : {ci_neu}")
    print(f"altered contested rate : {ci_alt}")
    print(f"\nexcess: {json.dumps(excess, indent=1, ensure_ascii=False)}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
