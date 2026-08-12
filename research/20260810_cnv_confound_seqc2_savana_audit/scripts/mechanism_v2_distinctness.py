#!/usr/bin/env python3
"""Revised mechanism test: AF *distinctness*, not AF spread.

The first mechanism hypothesis -- that copy number widens the within-unit
read-AF spread and so gives the score more to sort on -- was refuted by its own
test: neutral units have a LARGER median spread (0.783) than CN-altered ones
(0.522), the opposite of the prediction (Mann-Whitney p = 0.9997).

Re-reading the numbers points at a different mechanism.  The read-AF score is
compared as exact rationals, so a tie happens when candidate trees produce
*equal* scores.  What matters is therefore not how far apart the AF values are
but whether they are distinct at all.

On CN-neutral ground in a pure tumour cell line, a clonal mutation on a
single-copy haplotype reads out at exactly 1/1.  Many sites in the same unit
then share the identical value, every parent assignment scores the same, and
the unit stays tied.  On gain/LOH ground the same clone yields m/c values at
different grid points (1/3, 2/3, 1/1 ...), the scores separate, and the tie
breaks -- on copy number, not on lineage.

Prediction: units whose read-AF values are all identical should be far more
common on CN-neutral ground, and such units should almost never break ties.

Uses the exact `fraction` strings emitted by the C++ solver, so distinctness is
tested in the same arithmetic the solver used.
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict
from fractions import Fraction

from scipy.stats import fisher_exact

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
ANNOT = os.path.join(DATA, "hcc1395_unit_cn_annotation.jsonl")
TOPO = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
    "samples/HCC1395/HCC1395.topology.jsonl"
)
OUT = os.path.join(DATA, "mechanism_v2_distinctness.json")


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {
        "point_percent": round(100 * p, 2),
        "lo_percent": round(100 * max(0.0, c - h), 2),
        "hi_percent": round(100 * min(1.0, c + h), 2),
        "n": n,
    }


def cmh(tables):
    num = den = 0.0
    sn = sv = 0.0
    used = 0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n <= 1:
            continue
        r1, r2, c1 = a + b, c + d, a + c
        if r1 == 0 or r2 == 0 or c1 == 0 or (n - c1) == 0:
            continue
        used += 1
        num += a * d / n
        den += b * c / n
        sn += a - r1 * c1 / n
        sv += r1 * r2 * c1 * (n - c1) / (n * n * (n - 1))
    if den == 0 or sv == 0:
        return {"strata_used": used, "common_odds_ratio": None, "p_value": None}
    chi = (abs(sn) - 0.5) ** 2 / sv
    return {
        "strata_used": used,
        "common_odds_ratio": round(num / den, 4),
        "chi_square": round(chi, 4),
        "p_value": math.erfc(math.sqrt(chi / 2.0)),
    }


def tree_bin(n):
    if n == 2:
        return "2"
    if n == 3:
        return "3"
    if n <= 5:
        return "4-5"
    if n <= 10:
        return "6-10"
    return ">10"


def main():
    cn_by_region = {}
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") == "ANNOTATED":
                cn_by_region[r["region_id"]] = r["seqc2_class"]

    rows = []
    with open(TOPO) as fh:
        for ln in fh:
            u = json.loads(ln)
            rid = u.get("region_id")
            cls = cn_by_region.get(rid)
            if cls is None or u.get("unit_status") != "ranked":
                continue
            # Only active sites enter the score s(p,c); af_coverage also carries
            # all-reference rows (fraction "0/1") for inactive positions, which
            # can never contribute a term and must not count toward distinctness.
            active = set(u.get("active_positions") or [])
            fracs = []
            for c in u.get("af_coverage") or []:
                if c.get("position") not in active:
                    continue
                fr = c.get("fraction")
                if not fr:
                    continue
                try:
                    fracs.append(Fraction(fr))
                except (ValueError, ZeroDivisionError):
                    continue
            if len(fracs) < 2:
                continue
            distinct = len(set(fracs))
            ttc = int(u.get("total_tree_count") or 0)
            rows.append(
                {
                    "cn_class": cls,
                    "altered": cls != "neutral",
                    "n_af": len(fracs),
                    "n_distinct_af": distinct,
                    "all_identical": distinct == 1,
                    "distinct_ratio": distinct / len(fracs),
                    "total_tree_count": ttc,
                    "best_tree_unique": bool(u.get("best_tree_unique")),
                    "af_broke_tie": ttc > 1 and bool(u.get("best_tree_unique")),
                    "contested": ttc > 1,
                }
            )

    # ---- prediction 1: all-identical AF is more common on neutral ground ----
    by_class = defaultdict(lambda: {"n": 0, "identical": 0, "ratios": []})
    for r in rows:
        b = by_class[r["cn_class"]]
        b["n"] += 1
        b["identical"] += 1 if r["all_identical"] else 0
        b["ratios"].append(r["distinct_ratio"])

    ident_out = {}
    for cls, b in sorted(by_class.items()):
        rs = sorted(b["ratios"])
        ident_out[cls] = {
            "n_units": b["n"],
            "all_af_identical": b["identical"],
            "all_af_identical_percent": round(100 * b["identical"] / b["n"], 2),
            "median_distinct_ratio": round(rs[len(rs) // 2], 4),
        }

    neu_rows = [r for r in rows if not r["altered"]]
    alt_rows = [r for r in rows if r["altered"]]
    a = sum(1 for r in neu_rows if r["all_identical"])
    b_ = len(neu_rows) - a
    c_ = sum(1 for r in alt_rows if r["all_identical"])
    d_ = len(alt_rows) - c_
    or1, p1 = fisher_exact([[a, b_], [c_, d_]])
    pred1 = {
        "statement": "all read-AF values identical within a unit",
        "neutral": wilson(a, len(neu_rows)),
        "cn_altered": wilson(c_, len(alt_rows)),
        "fisher_odds_ratio_neutral_vs_altered": round(or1, 4),
        "fisher_p_value": p1,
        "verdict": "SUPPORTED" if (a / max(len(neu_rows), 1)) > (c_ / max(len(alt_rows), 1)) and p1 < 0.05 else "NOT SUPPORTED",
    }

    # ---- prediction 2: identical-AF units almost never break ties ----
    contested = [r for r in rows if r["contested"]]
    ident_contested = [r for r in contested if r["all_identical"]]
    dist_contested = [r for r in contested if not r["all_identical"]]
    ia = sum(1 for r in ident_contested if r["af_broke_tie"])
    da = sum(1 for r in dist_contested if r["af_broke_tie"])
    or2, p2 = fisher_exact(
        [[da, len(dist_contested) - da], [ia, len(ident_contested) - ia]]
    )
    pred2 = {
        "statement": "among contested units, distinct-AF units break ties more often than identical-AF units",
        "identical_af_units": wilson(ia, len(ident_contested)),
        "distinct_af_units": wilson(da, len(dist_contested)),
        "fisher_odds_ratio_distinct_vs_identical": round(or2, 4),
        "fisher_p_value": p2,
        "verdict": "SUPPORTED"
        if len(ident_contested)
        and (da / max(len(dist_contested), 1)) > (ia / max(len(ident_contested), 1))
        and p2 < 0.05
        else "NOT SUPPORTED",
    }

    # ---- does distinctness explain away the CN association? ----
    # Stratify the CN -> af_broke_tie association on AF distinctness.
    strata = defaultdict(lambda: {"alt": [0, 0], "neu": [0, 0]})
    for r in contested:
        key = "identical" if r["all_identical"] else "distinct"
        s = strata[(key, tree_bin(r["total_tree_count"]))]
        side = "alt" if r["altered"] else "neu"
        s[side][0 if r["af_broke_tie"] else 1] += 1
    tables = []
    detail = {}
    for k, s in sorted(strata.items(), key=lambda kv: str(kv[0])):
        aa, bb = s["alt"]
        cc, dd = s["neu"]
        tables.append((aa, bb, cc, dd))
        detail[str(k)] = {
            "altered_n": aa + bb,
            "altered_rate": round(100 * aa / (aa + bb), 2) if (aa + bb) else None,
            "neutral_n": cc + dd,
            "neutral_rate": round(100 * cc / (cc + dd), 2) if (cc + dd) else None,
        }
    mediation = {
        "logic": "if AF distinctness is the mechanism, controlling for it should shrink the CN odds ratio toward 1",
        "strata": detail,
        "cmh_controlling_distinctness_and_tree_count": cmh(tables),
    }

    # ---- the contrast that remains once arithmetically-unresolvable units are removed ----
    nd = [r for r in contested if not r["all_identical"]]
    na = sum(1 for r in nd if r["altered"] and r["af_broke_tie"])
    nb = sum(1 for r in nd if r["altered"] and not r["af_broke_tie"])
    nc = sum(1 for r in nd if not r["altered"] and r["af_broke_tie"])
    ndd = sum(1 for r in nd if not r["altered"] and not r["af_broke_tie"])
    nor, np_ = fisher_exact([[na, nb], [nc, ndd]])
    non_degenerate = {
        "logic": "units whose active read-AF values are all identical cannot break a tie by construction; excluding them isolates whatever CN effect remains",
        "cn_altered": wilson(na, na + nb),
        "cn_neutral": wilson(nc, nc + ndd),
        "fisher_odds_ratio": round(nor, 4),
        "fisher_p_value": np_,
        "table_altered_broke_not": [na, nb],
        "table_neutral_broke_not": [nc, ndd],
        "reading": "no CN effect detectable once arithmetic resolvability is accounted for"
        if np_ > 0.05
        else "a residual CN effect survives",
    }

    # ---- how many contested units are arithmetically unresolvable, by arm ----
    alt_c = [r for r in contested if r["altered"]]
    neu_c = [r for r in contested if not r["altered"]]
    unresolvable = {
        "cn_altered": wilson(sum(1 for r in alt_c if r["all_identical"]), len(alt_c)),
        "cn_neutral": wilson(sum(1 for r in neu_c if r["all_identical"]), len(neu_c)),
    }

    out = {
        "generated_by": os.path.basename(__file__),
        "supersedes": "mechanism_and_intervals.json spread hypothesis (refuted)",
        "non_degenerate_contrast": non_degenerate,
        "arithmetically_unresolvable_share": unresolvable,
        "prediction_1_identical_af_more_common_on_neutral": pred1,
        "prediction_2_identical_af_units_stay_tied": pred2,
        "identical_af_by_cn_class": ident_out,
        "mediation_check": mediation,
        "counts": {
            "units_with_ge2_af_values": len(rows),
            "contested": len(contested),
            "contested_identical_af": len(ident_contested),
            "contested_distinct_af": len(dist_contested),
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print("all-AF-identical rate by CN class:")
    for k, v in ident_out.items():
        print(
            f"  {k:10s} n={v['n_units']:5d} identical={v['all_af_identical_percent']:6.2f}% "
            f"median_distinct_ratio={v['median_distinct_ratio']:.3f}"
        )
    print(f"\nprediction 1: {json.dumps(pred1, indent=1, ensure_ascii=False)}")
    print(f"\nprediction 2: {json.dumps(pred2, indent=1, ensure_ascii=False)}")
    print(f"\nmediation CMH: {mediation['cmh_controlling_distinctness_and_tree_count']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
