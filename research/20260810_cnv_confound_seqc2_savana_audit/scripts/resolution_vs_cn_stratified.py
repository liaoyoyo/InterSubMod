#!/usr/bin/env python3
"""Does copy number inflate the reported tree-uniqueness rates?

The canonical cohort reports 55.10% unique_best_tree and 88.26% one rooted
unlabeled topology.  Those rates are produced by the read-AF score breaking
ties between structurally equivalent candidates.  Under gain/LOH the AF of two
mutations can differ purely from copy number, so the score may be breaking ties
on copy-number artefacts rather than on lineage signal -- which would inflate
both headline numbers.

Uniqueness is therefore decomposed into two layers:

  structure_unique : total_tree_count == 1
      The enumeration itself returned one tree.  Read-AF played no role, so
      this layer is copy-number robust.

  af_broke_tie     : total_tree_count > 1 and best_tree_unique
      Multiple structurally valid trees existed and the read-AF score picked
      one.  This is exactly the layer copy number can corrupt.

  still_tied       : total_tree_count > 1 and not best_tree_unique

The comparison of interest is af_broke_tie rate in CN-neutral versus CN-altered
ground.  Local mutation density is a known confounder in this project, so the
test is stratified on both the number of candidate trees and the active bit
count k, and combined with Cochran-Mantel-Haenszel.  A marginal difference that
disappears under stratification is reported as such.
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict

from scipy.stats import fisher_exact

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
ANNOT = os.path.join(DATA, "hcc1395_unit_cn_annotation.jsonl")
TOPO = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
    "samples/HCC1395/HCC1395.topology.jsonl"
)
OUT = os.path.join(DATA, "resolution_vs_cn_stratified.json")


def tree_bin(n):
    if n <= 1:
        return "1"
    if n == 2:
        return "2"
    if n == 3:
        return "3"
    if n <= 5:
        return "4-5"
    if n <= 10:
        return "6-10"
    return ">10"


def k_bin(k):
    if k <= 2:
        return str(k)
    if k <= 4:
        return "3-4"
    if k <= 6:
        return "5-6"
    return ">6"


def cmh(tables):
    """Cochran-Mantel-Haenszel common odds ratio + chi-square test.

    tables: list of (a, b, c, d) = (exposed_yes, exposed_no, unexposed_yes, unexposed_no)
    """
    num = den = 0.0
    stat_num = 0.0
    stat_var = 0.0
    used = 0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        r1, r2 = a + b, c + d
        c1 = a + c
        if r1 == 0 or r2 == 0 or c1 == 0 or (n - c1) == 0:
            continue
        used += 1
        num += a * d / n
        den += b * c / n
        exp = r1 * c1 / n
        var = r1 * r2 * c1 * (n - c1) / (n * n * (n - 1)) if n > 1 else 0.0
        stat_num += a - exp
        stat_var += var
    if den == 0 or stat_var == 0:
        return {"strata_used": used, "common_odds_ratio": None, "p_value": None}
    or_mh = num / den
    chi = (abs(stat_num) - 0.5) ** 2 / stat_var
    p = math.erfc(math.sqrt(chi / 2.0))
    return {
        "strata_used": used,
        "common_odds_ratio": round(or_mh, 4),
        "chi_square": round(chi, 4),
        "p_value": p,
    }


def main():
    # ---- load CN annotation ----
    cn_by_region = {}
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") == "ANNOTATED":
                cn_by_region[r["region_id"]] = r["seqc2_class"]

    # ---- join with topology facts ----
    rows = []
    with open(TOPO) as fh:
        for ln in fh:
            u = json.loads(ln)
            rid = u.get("region_id")
            cls = cn_by_region.get(rid)
            if cls is None:
                continue
            if u.get("unit_status") != "ranked":
                continue
            ttc = u.get("total_tree_count")
            if ttc is None:
                continue
            ttc = int(ttc)
            rows.append(
                {
                    "region_id": rid,
                    "cn_class": cls,
                    "altered": cls != "neutral",
                    "k": int(u.get("active_bit_count") or 0),
                    "total_tree_count": ttc,
                    "best_tree_unique": bool(u.get("best_tree_unique")),
                    "structure_unique": ttc == 1,
                    "af_broke_tie": ttc > 1 and bool(u.get("best_tree_unique")),
                    "still_tied": ttc > 1 and not bool(u.get("best_tree_unique")),
                    "min_vertex_sets": u.get("minimum_vertex_set_count"),
                }
            )

    n = len(rows)

    # ---- layer decomposition per CN class ----
    layers = defaultdict(Counter)
    for r in rows:
        c = layers[r["cn_class"]]
        c["n"] += 1
        if r["structure_unique"]:
            c["structure_unique"] += 1
        elif r["af_broke_tie"]:
            c["af_broke_tie"] += 1
        else:
            c["still_tied"] += 1

    layer_out = {}
    for cls, c in sorted(layers.items()):
        tot = c["n"]
        layer_out[cls] = {
            "n": tot,
            "structure_unique": c["structure_unique"],
            "structure_unique_percent": round(100 * c["structure_unique"] / tot, 2),
            "af_broke_tie": c["af_broke_tie"],
            "af_broke_tie_percent": round(100 * c["af_broke_tie"] / tot, 2),
            "still_tied": c["still_tied"],
            "still_tied_percent": round(100 * c["still_tied"] / tot, 2),
            "overall_unique_percent": round(
                100 * (c["structure_unique"] + c["af_broke_tie"]) / tot, 2
            ),
        }

    # ---- the contested layer: among units where AF had a job to do ----
    contested = [r for r in rows if r["total_tree_count"] > 1]

    def rate(subset):
        if not subset:
            return None
        w = sum(1 for r in subset if r["af_broke_tie"])
        return {
            "n": len(subset),
            "af_broke_tie": w,
            "percent": round(100 * w / len(subset), 2),
        }

    contested_by_class = {
        cls: rate([r for r in contested if r["cn_class"] == cls])
        for cls in sorted({r["cn_class"] for r in contested})
    }

    neu = [r for r in contested if not r["altered"]]
    alt = [r for r in contested if r["altered"]]
    marginal = {
        "neutral": rate(neu),
        "cn_altered": rate(alt),
    }
    a = sum(1 for r in alt if r["af_broke_tie"])
    b = len(alt) - a
    c_ = sum(1 for r in neu if r["af_broke_tie"])
    d = len(neu) - c_
    odds, p_marginal = fisher_exact([[a, b], [c_, d]])
    marginal["fisher_odds_ratio_altered_vs_neutral"] = round(odds, 4)
    marginal["fisher_p_value"] = p_marginal
    marginal["absolute_difference_pp"] = round(
        (100 * a / len(alt)) - (100 * c_ / len(neu)), 2
    ) if alt and neu else None

    # ---- stratified: candidate-tree count, then k, then both ----
    def stratify(keyfn, label):
        strata = defaultdict(lambda: {"alt": [0, 0], "neu": [0, 0]})
        for r in contested:
            s = strata[keyfn(r)]
            side = "alt" if r["altered"] else "neu"
            s[side][0 if r["af_broke_tie"] else 1] += 1
        tables = []
        detail = {}
        for kk, s in sorted(strata.items(), key=lambda kv: str(kv[0])):
            aa, bb = s["alt"]
            cc, dd = s["neu"]
            tables.append((aa, bb, cc, dd))
            detail[str(kk)] = {
                "altered_n": aa + bb,
                "altered_af_broke_tie_percent": round(100 * aa / (aa + bb), 2)
                if (aa + bb)
                else None,
                "neutral_n": cc + dd,
                "neutral_af_broke_tie_percent": round(100 * cc / (cc + dd), 2)
                if (cc + dd)
                else None,
            }
        return {"label": label, "strata": detail, "cmh": cmh(tables)}

    strat_tree = stratify(lambda r: tree_bin(r["total_tree_count"]), "total_tree_count")
    strat_k = stratify(lambda r: k_bin(r["k"]), "active_bit_count")
    strat_both = stratify(
        lambda r: (tree_bin(r["total_tree_count"]), k_bin(r["k"])),
        "total_tree_count x active_bit_count",
    )

    # ---- what the headline rates become if restricted to CN-neutral ground ----
    def headline(subset):
        if not subset:
            return None
        tot = len(subset)
        uniq = sum(1 for r in subset if r["best_tree_unique"])
        struct = sum(1 for r in subset if r["structure_unique"])
        return {
            "n": tot,
            "unique_best_tree": uniq,
            "unique_best_tree_percent": round(100 * uniq / tot, 2),
            "structure_unique": struct,
            "structure_unique_percent": round(100 * struct / tot, 2),
        }

    headlines = {
        "all_ranked_units": headline(rows),
        "cn_neutral_only": headline([r for r in rows if not r["altered"]]),
        "cn_altered_only": headline([r for r in rows if r["altered"]]),
    }

    out = {
        "generated_by": os.path.basename(__file__),
        "scope": "HCC1395, chr1-22, canonical ranked units with SEQC2 CN annotation",
        "definitions": {
            "structure_unique": "total_tree_count == 1; enumeration alone gave one tree; read-AF played no role; copy-number robust",
            "af_broke_tie": "total_tree_count > 1 and best_tree_unique; read-AF selected the winner; copy-number vulnerable",
            "still_tied": "total_tree_count > 1 and read-AF did not produce a unique winner",
        },
        "counts": {"ranked_units_with_cn": n, "contested_units": len(contested)},
        "layer_decomposition_by_cn_class": layer_out,
        "af_broke_tie_by_cn_class_contested_only": contested_by_class,
        "marginal_test": marginal,
        "stratified_by_candidate_tree_count": strat_tree,
        "stratified_by_k": strat_k,
        "stratified_by_both": strat_both,
        "headline_rates": headlines,
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"ranked units with CN: {n}; contested (>1 candidate tree): {len(contested)}")
    print("\nlayer decomposition by CN class:")
    for cls, v in layer_out.items():
        print(
            f"  {cls:10s} n={v['n']:5d} struct_uniq={v['structure_unique_percent']:6.2f}% "
            f"af_broke={v['af_broke_tie_percent']:6.2f}% tied={v['still_tied_percent']:6.2f}%"
        )
    print("\ncontested-only af_broke_tie rate:")
    for cls, v in contested_by_class.items():
        if v:
            print(f"  {cls:10s} n={v['n']:5d} {v['percent']:6.2f}%")
    print(f"\nmarginal: {marginal}")
    print(f"\nCMH by tree count : {strat_tree['cmh']}")
    print(f"CMH by k          : {strat_k['cmh']}")
    print(f"CMH by both       : {strat_both['cmh']}")
    print(f"\nheadlines: {json.dumps(headlines, indent=1)}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
