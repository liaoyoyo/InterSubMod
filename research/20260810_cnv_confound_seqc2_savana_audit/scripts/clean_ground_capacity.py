#!/usr/bin/env python3
"""How much can the CN-clean ground actually carry, and what is special about CN ground?

Two questions the first pass left unanswered:

A. Is the CN-neutral subset informative enough to stand on its own?
   4.56% of units survive a copy-neutral filter.  "How many" is not the same as
   "enough": what matters is what those units can and cannot support -- their k
   and candidate-tree distribution (are they trivially small?), their genomic
   spread (do they represent the genome?), how many carry a genuinely subclonal
   read-AF, and what effect size the surviving sample can still detect.

B. Is there anything structurally distinctive about CN ground, beyond "it is
   confounded"?  Specifically whether copy number changes the *shape* of the
   reconstructed trees (coarse_class from the canonical census) rather than only
   the ability to pick one, and whether read-AF collapses onto the m/c grid in a
   way that tracks the copy number itself.

Descriptive throughout; no hypothesis is tested here that is not stated.
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
ANNOT = os.path.join(DATA, "hcc1395_unit_cn_annotation.jsonl")
SITES = os.path.join(DATA, "hcc1395_site_cn_annotation.jsonl")
TOPO = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
    "samples/HCC1395/HCC1395.topology.jsonl"
)
CENSUS = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1/HCC1395.census.jsonl"
)
OUT = os.path.join(DATA, "clean_ground_capacity.json")

CLONAL_AF = 0.95  # read-AF at or above this reads as clonal on a single-copy haplotype


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
        "point_percent": round(100 * p, 2),
        "lo_percent": round(100 * max(0.0, c - h), 2),
        "hi_percent": round(100 * min(1.0, c + h), 2),
    }


def detectable_effect(n1, n2, p_base, alpha=0.05, power=0.80):
    """Smallest second-arm proportion detectable at the given n, two-sided."""
    za, zb = 1.959963985, 0.841621234
    lo, hi = p_base, 0.999
    for _ in range(200):
        mid = (lo + hi) / 2
        pbar = (p_base * n1 + mid * n2) / (n1 + n2)
        se0 = math.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
        se1 = math.sqrt(p_base * (1 - p_base) / n1 + mid * (1 - mid) / n2)
        if se1 == 0:
            break
        z = (abs(mid - p_base) - za * se0) / se1
        if z >= zb:
            hi = mid
        else:
            lo = mid
    p2 = hi
    odds = (p2 / (1 - p2)) / (p_base / (1 - p_base)) if 0 < p_base < 1 and p2 < 1 else None
    return {
        "baseline_percent": round(100 * p_base, 2),
        "min_detectable_percent": round(100 * p2, 2),
        "min_detectable_odds_ratio": round(odds, 3) if odds else None,
        "alpha": alpha,
        "power": power,
    }


def main():
    # ---- census: tree shape per unit ----
    shape = {}
    with open(CENSUS) as fh:
        for ln in fh:
            r = json.loads(ln)
            classes = r.get("coarse_class_tree_counts") or {}
            shape[r["region_id"]] = {
                "resolution_class": r.get("resolution_class"),
                "coarse_classes": list(classes.keys()),
                "coarse_class_count": r.get("coarse_class_count"),
                "signatures": [
                    s.get("shape_signature") for s in (r.get("topology_signatures") or [])
                ],
                "signature_count": r.get("topology_signature_count"),
            }

    cn_by_region = {}
    chrom_by_region = {}
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") == "ANNOTATED":
                cn_by_region[r["region_id"]] = r["seqc2_class"]
                chrom_by_region[r["region_id"]] = r["chrom"]

    units = []
    with open(TOPO) as fh:
        for ln in fh:
            u = json.loads(ln)
            rid = u.get("region_id")
            cls = cn_by_region.get(rid)
            if cls is None or u.get("unit_status") != "ranked":
                continue
            ttc = int(u.get("total_tree_count") or 0)
            sh = shape.get(rid, {})
            units.append(
                {
                    "region_id": rid,
                    "chrom": chrom_by_region.get(rid),
                    "cn_class": cls,
                    "altered": cls != "neutral",
                    "k": int(u.get("active_bit_count") or 0),
                    "total_tree_count": ttc,
                    "best_tree_unique": bool(u.get("best_tree_unique")),
                    "recurrence_required": bool(u.get("recurrence_required")),
                    "min_vertex_sets": int(u.get("minimum_vertex_set_count") or 0),
                    "coarse_classes": sh.get("coarse_classes") or [],
                    "signatures": sh.get("signatures") or [],
                    "signature_count": sh.get("signature_count"),
                }
            )

    neu = [u for u in units if not u["altered"]]
    alt = [u for u in units if u["altered"]]

    def dist(vals):
        c = Counter(vals)
        tot = sum(c.values())
        return {
            "counts": {str(k): v for k, v in sorted(c.items(), key=lambda kv: str(kv[0]))},
            "percent": {
                str(k): round(100 * v / tot, 2)
                for k, v in sorted(c.items(), key=lambda kv: str(kv[0]))
            },
            "n": tot,
        }

    def kbin(k):
        return str(k) if k <= 4 else "5-8" if k <= 8 else ">8"

    def tbin(t):
        return "1" if t <= 1 else "2-3" if t <= 3 else "4-10" if t <= 10 else ">10"

    # ---------------- A. capacity of the clean ground ----------------
    neu_sites = []
    alt_sites_by_cn = defaultdict(list)
    site_chrom_neu = Counter()
    with open(SITES) as fh:
        for ln in fh:
            r = json.loads(ln)
            af = r.get("read_af")
            if af is None:
                continue
            if r["seqc2_class"] == "neutral":
                neu_sites.append(af)
                site_chrom_neu[r["chrom"]] += 1
            else:
                alt_sites_by_cn[r["seqc2_total_cn"]].append(af)

    subclonal = sum(1 for a in neu_sites if a < CLONAL_AF)

    capacity = {
        "units": {
            "neutral_ranked": len(neu),
            "altered_ranked": len(alt),
            "neutral_share_percent": round(100 * len(neu) / len(units), 2),
        },
        "k_distribution": {"neutral": dist([kbin(u["k"]) for u in neu]),
                           "altered": dist([kbin(u["k"]) for u in alt])},
        "candidate_tree_distribution": {
            "neutral": dist([tbin(u["total_tree_count"]) for u in neu]),
            "altered": dist([tbin(u["total_tree_count"]) for u in alt]),
        },
        "genomic_representativeness": {
            "chromosomes_with_any_neutral_unit": len({u["chrom"] for u in neu}),
            "chromosomes_total": len({u["chrom"] for u in units}),
            "neutral_units_per_chrom": dict(
                sorted(Counter(u["chrom"] for u in neu).items(), key=lambda kv: -kv[1])
            ),
        },
        "neutral_site_read_af": {
            "n_sites_with_af": len(neu_sites),
            "clonal_at_or_above_%.2f" % CLONAL_AF: wilson(
                len(neu_sites) - subclonal, len(neu_sites)
            ),
            "subclonal_below_%.2f" % CLONAL_AF: wilson(subclonal, len(neu_sites)),
            "note": "on copy-neutral non-LOH ground in a pure tumour line, a single-copy clonal mutation reads 1/1; a materially lower AF is the only place a subclonal fraction can be read without a copy-number explanation",
        },
        "statistical_power": {
            "contested_neutral_n": sum(1 for u in neu if u["total_tree_count"] > 1),
            "contested_altered_n": sum(1 for u in alt if u["total_tree_count"] > 1),
            "detectable_at_current_n": detectable_effect(
                sum(1 for u in neu if u["total_tree_count"] > 1),
                sum(1 for u in alt if u["total_tree_count"] > 1),
                0.2549,
            ),
        },
    }

    # ---------------- B. what is distinctive about CN ground ----------------
    # B1. tree shape composition
    def shape_profile(subset):
        cc = Counter()
        sig = Counter()
        multi = 0
        rec = 0
        for u in subset:
            for c in u["coarse_classes"]:
                cc[c] += 1
            for s in u["signatures"]:
                sig[s] += 1
            if (u["signature_count"] or 0) > 1:
                multi += 1
            if u["recurrence_required"]:
                rec += 1
        n = len(subset)
        return {
            "n_units": n,
            "coarse_class_share_percent": {
                k: round(100 * v / n, 2) for k, v in cc.most_common()
            },
            "top_shape_signatures_percent": {
                k: round(100 * v / n, 2) for k, v in sig.most_common(6)
            },
            "multi_signature_percent": round(100 * multi / n, 2) if n else None,
            "recurrence_required_percent": round(100 * rec / n, 2) if n else None,
        }

    shapes = {
        "neutral": shape_profile(neu),
        "altered": shape_profile(alt),
        "by_cn_class": {
            c: shape_profile([u for u in units if u["cn_class"] == c])
            for c in sorted({u["cn_class"] for u in units})
        },
    }

    # B2. does read-AF track the copy number itself?
    af_by_cn = {}
    for cnv, vals in sorted(alt_sites_by_cn.items()):
        if len(vals) < 30:
            continue
        v = sorted(vals)
        n = len(v)
        # distance to the nearest 1/cn grid point, the multiplicity prediction
        grid = [m / cnv for m in range(1, int(round(cnv)) + 1)] if cnv >= 1 else [1.0]
        near = [min(abs(a - g) for g in grid) for a in v]
        near.sort()
        af_by_cn[str(cnv)] = {
            "n": n,
            "median_af": round(v[n // 2], 4),
            "q1": round(v[n // 4], 4),
            "q3": round(v[(3 * n) // 4], 4),
            "predicted_lowest_grid_point": round(1.0 / cnv, 4) if cnv else None,
            "median_distance_to_nearest_grid_point": round(near[len(near) // 2], 4),
            "fraction_within_0.05_of_a_grid_point": round(
                sum(1 for d in near if d <= 0.05) / n, 4
            ),
        }
    neu_v = sorted(neu_sites)
    af_by_cn["2.0_neutral_reference"] = {
        "n": len(neu_v),
        "median_af": round(neu_v[len(neu_v) // 2], 4) if neu_v else None,
        "q1": round(neu_v[len(neu_v) // 4], 4) if neu_v else None,
        "q3": round(neu_v[(3 * len(neu_v)) // 4], 4) if neu_v else None,
    }

    # B3. resolution as a function of copy number
    res_by_cn = {}
    for c in sorted({u["cn_class"] for u in units}):
        sub = [u for u in units if u["cn_class"] == c]
        con = [u for u in sub if u["total_tree_count"] > 1]
        res_by_cn[c] = {
            "n_units": len(sub),
            "structure_unique_percent": round(
                100 * sum(1 for u in sub if u["total_tree_count"] == 1) / len(sub), 2
            )
            if sub
            else None,
            "median_candidate_trees": sorted(u["total_tree_count"] for u in sub)[
                len(sub) // 2
            ]
            if sub
            else None,
            "contested_n": len(con),
            "af_broke_tie_percent": round(
                100 * sum(1 for u in con if u["best_tree_unique"]) / len(con), 2
            )
            if con
            else None,
        }

    out = {
        "generated_by": os.path.basename(__file__),
        "scope": "HCC1395 chr1-22, canonical ranked units with SEQC2 CN annotation",
        "clonal_af_threshold": CLONAL_AF,
        "A_clean_ground_capacity": capacity,
        "B_what_is_distinctive_about_cn_ground": {
            "tree_shape_profiles": shapes,
            "read_af_by_total_copy_number": af_by_cn,
            "resolution_by_cn_class": res_by_cn,
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print("=== A. clean ground ===")
    print(json.dumps(capacity["units"], indent=1))
    print("k dist neutral :", capacity["k_distribution"]["neutral"]["percent"])
    print("k dist altered :", capacity["k_distribution"]["altered"]["percent"])
    print("trees neutral  :", capacity["candidate_tree_distribution"]["neutral"]["percent"])
    print("trees altered  :", capacity["candidate_tree_distribution"]["altered"]["percent"])
    print("neutral sites  :", capacity["neutral_site_read_af"])
    print("power          :", capacity["statistical_power"])
    print("\n=== B. shape ===")
    print("neutral:", shapes["neutral"]["coarse_class_share_percent"],
          "multi-sig", shapes["neutral"]["multi_signature_percent"])
    print("altered:", shapes["altered"]["coarse_class_share_percent"],
          "multi-sig", shapes["altered"]["multi_signature_percent"])
    print("\n=== B. af by CN ===")
    for k, v in af_by_cn.items():
        print(" ", k, v)
    print("\n=== B. resolution by CN ===")
    for k, v in res_by_cn.items():
        print(" ", k, v)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
