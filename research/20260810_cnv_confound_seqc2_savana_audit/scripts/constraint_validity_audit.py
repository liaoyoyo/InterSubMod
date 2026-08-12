#!/usr/bin/env python3
"""Do the linkage constraints actually hold, and do they survive copy number?

Two constraints were claimed in the previous step.  They do not have the same
standing and must be audited separately.

ORDERING ("A before B", from patterns {RR, AR, AA})
    B is only ever seen on molecules that already carry A.  Under an
    infinite-sites model mutations only accumulate along a molecular lineage, so
    an A+B molecule must descend from an A molecule.  Nothing in that argument
    references copy number, ploidy or purity -- it is a statement about molecular
    descent.  The claim is therefore copy-number-free BY CONSTRUCTION, and the
    real threat is sampling: concluding "RA was never seen" when RA exists but
    was missed.  That is a power question, quantified below.

EXCLUSION ("A and B are on different branches", from patterns {RR, AR, RA})
    No molecule carries both.  On copy-neutral 1+1 ground each cell contributes
    one molecule per haplotype, so "no molecule carries both" does imply "no cell
    carries both" -- a genuine cellular-branch statement.
    Under gain the haplotype carries c > 1 copies per cell.  Mutation A can sit
    on copy 1 and B on copy 2 OF THE SAME CELL; no single molecule carries both,
    yet the cell carries both.  The molecular observation is unchanged but the
    cellular inference collapses.

So the honest position is that ordering transfers to CN-altered ground and
exclusion does not.  This script measures how much of each claim survives.
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
ANNOT = os.path.join(DATA, "hcc1395_unit_cn_annotation.jsonl")
TOPO = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
    "samples/HCC1395/HCC1395.topology.jsonl"
)
OUT = os.path.join(DATA, "constraint_validity_audit.json")

CONF = 0.95  # confidence for the "pattern truly absent" exclusion bound


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {
        "count": k,
        "n": n,
        "percent": round(100 * p, 2),
        "lo": round(100 * max(0.0, c - h), 2),
        "hi": round(100 * min(1.0, c + h), 2),
    }


def min_detectable_fraction(depth, conf=CONF):
    """Smallest true pattern fraction that would have been seen with prob >= conf.

    If a pattern has true fraction p, the chance of missing it in `depth`
    independent molecules is (1-p)^depth.  Requiring that to be <= 1-conf gives
    p >= 1 - (1-conf)^(1/depth).
    """
    if depth <= 0:
        return 1.0
    return 1.0 - (1.0 - conf) ** (1.0 / depth)


def main():
    cn = {}
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") == "ANNOTATED":
                cn[r["region_id"]] = r["seqc2_class"]

    rows = []
    with open(TOPO) as fh:
        for ln in fh:
            u = json.loads(ln)
            rid = u.get("region_id")
            cls = cn.get(rid)
            if cls is None or u.get("unit_status") != "ranked":
                continue
            k = int(u.get("active_bit_count") or 0)
            if k < 2:
                continue
            depths = []
            for c in u.get("af_coverage") or []:
                a, r_ = c.get("alt_reads"), c.get("ref_reads")
                if a is not None and r_ is not None:
                    depths.append(a + r_)
            if not depths:
                continue
            # linkage across the unit is bounded by its least-covered site
            depth = min(depths)
            morph = u.get("representative_best_morphology") or {}
            edges = u.get("representative_best_edges") or []
            deep = sum(1 for e in edges if e.get("parent_label") != "ROOT")
            rows.append(
                {
                    "cn_class": cls,
                    "altered": cls != "neutral",
                    "k": k,
                    "depth": depth,
                    "has_sister": bool(morph.get("has_sister")),
                    "deep_edges": deep,
                    "mdf": min_detectable_fraction(depth),
                }
            )

    # ---------- exclusion claims: where can they reach the cellular level? ----------
    excl = [r for r in rows if r["has_sister"]]
    excl_neu = [r for r in excl if not r["altered"]]
    excl_alt = [r for r in excl if r["altered"]]
    cellular_valid = wilson(len(excl_neu), len(excl)) if excl else None

    # ---------- ordering claims: all CN states, but bounded by sampling ----------
    order_total = sum(r["deep_edges"] for r in rows)
    order_neu = sum(r["deep_edges"] for r in rows if not r["altered"])
    order_alt = sum(r["deep_edges"] for r in rows if r["altered"])

    def depth_profile(sub):
        if not sub:
            return None
        d = sorted(r["depth"] for r in sub)
        m = sorted(r["mdf"] for r in sub)
        n = len(d)
        return {
            "n_units": n,
            "median_linking_depth": d[n // 2],
            "q1_depth": d[n // 4],
            "q3_depth": d[(3 * n) // 4],
            "median_min_detectable_pattern_fraction": round(m[n // 2], 4),
            "units_with_depth_ge_20": sum(1 for x in d if x >= 20),
            "units_with_depth_ge_50": sum(1 for x in d if x >= 50),
            "fraction_depth_ge_20": round(sum(1 for x in d if x >= 20) / n, 4),
        }

    # how strong is the "absent pattern" bound, banded by depth
    bands = {}
    for lo, hi, lab in [(3, 9, "3-9"), (10, 19, "10-19"), (20, 49, "20-49"), (50, 10**9, ">=50")]:
        sub = [r for r in rows if lo <= r["depth"] <= hi]
        if not sub:
            continue
        bands[lab] = {
            "units": len(sub),
            "percent_of_units": round(100 * len(sub) / len(rows), 2),
            "min_detectable_pattern_fraction_at_band_low": round(
                min_detectable_fraction(lo), 4
            ),
            "min_detectable_pattern_fraction_at_band_high": round(
                min_detectable_fraction(min(hi, 200)), 4
            ),
        }

    out = {
        "generated_by": os.path.basename(__file__),
        "scope": "HCC1395 chr1-22, canonical ranked units with k>=2 and per-site depth available",
        "units_audited": len(rows),
        "claim_1_ordering": {
            "statement": "A precedes B along the molecular lineage",
            "why_copy_number_free": "an A+B molecule must descend from an A molecule under infinite sites; the argument never invokes copy number, ploidy or purity",
            "valid_on_cn_altered_ground": True,
            "level_of_the_claim": "molecular lineage, not cellular lineage",
            "total_ordering_constraints": order_total,
            "on_cn_neutral": order_neu,
            "on_cn_altered": order_alt,
            "residual_threat": "sampling, not copy number: concluding a pattern is absent when it was merely missed",
            "sampling_power": depth_profile(rows),
            "power_by_depth_band": bands,
        },
        "claim_2_exclusion": {
            "statement": "A and B are on different branches",
            "molecular_level": "no molecule carries both - holds wherever it is observed",
            "cellular_level_requires": "one copy of the haplotype per cell, i.e. copy-neutral 1+1 ground",
            "why_it_breaks_under_gain": "with c>1 copies per cell, A on copy 1 and B on copy 2 of the SAME cell produce no A+B molecule; the molecular observation is unchanged but the cellular inference is invalid",
            "exclusion_units_total": len(excl),
            "exclusion_units_cn_neutral_cellular_valid": len(excl_neu),
            "exclusion_units_cn_altered_molecular_only": len(excl_alt),
            "share_reaching_cellular_level": cellular_valid,
            "verdict": "the exclusion constraint does NOT transfer to CN-altered ground at the cellular level; on CN-altered ground it may only be stated as a molecular-level observation",
        },
        "corrected_position": {
            "ordering": "usable genome-wide, including CN-altered ground, as a molecular-lineage constraint",
            "exclusion": "usable as a cellular-branch constraint ONLY on copy-neutral ground; elsewhere demote to molecular-level phrasing",
            "note": "this corrects the previous framing, which presented the exclusion constraint as the more valuable of the two without qualifying its cellular reach",
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"units audited: {len(rows)}")
    print("\n=== claim 1: ordering ===")
    print(f"  total constraints {order_total} (neutral {order_neu} / altered {order_alt})")
    p = out["claim_1_ordering"]["sampling_power"]
    print(f"  median linking depth {p['median_linking_depth']} (q1 {p['q1_depth']} q3 {p['q3_depth']})")
    print(f"  median min-detectable pattern fraction: {p['median_min_detectable_pattern_fraction']}")
    print(f"  depth>=20: {p['fraction_depth_ge_20']*100:.1f}% of units")
    print("  power by depth band:")
    for k, v in bands.items():
        print(
            f"    depth {k:6s} units={v['units']:5d} ({v['percent_of_units']:5.2f}%)  "
            f"可排除比例 > {v['min_detectable_pattern_fraction_at_band_low']*100:.1f}%–{v['min_detectable_pattern_fraction_at_band_high']*100:.1f}%"
        )
    print("\n=== claim 2: exclusion ===")
    e = out["claim_2_exclusion"]
    print(f"  exclusion units total       : {e['exclusion_units_total']}")
    print(f"  reach CELLULAR level (neutral): {e['exclusion_units_cn_neutral_cellular_valid']}")
    print(f"  molecular-only (CN-altered) : {e['exclusion_units_cn_altered_molecular_only']}")
    print(f"  share reaching cellular     : {cellular_valid}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
