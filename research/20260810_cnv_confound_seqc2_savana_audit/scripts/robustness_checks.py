#!/usr/bin/env python3
"""Two robustness checks the main analysis owes the reader.

Check A -- is the CN effect confounded by chromosome?
    The CN-neutral comparison group is not a random sample of the genome: 87.7% of
    neutral units sit on chr21, chr2, chr6 and chr16, while six chromosomes have no
    neutral unit at all.  A neutral-vs-altered contrast could therefore be picking up
    chromosome identity rather than copy number.  The test is repeated stratified by
    chromosome, and separately restricted to the four chromosomes that contain both
    arms of the comparison, so the contrast is made within chromosome.

Check B -- is "identical AF units never break ties" a tautology?
    The frozen score is s(p,c) = sum_{i in p}(AF_i - AF_j(c)), summed over edges.  If every
    read-AF in a unit equals v, each term is (v - v) = 0, so every candidate tree scores
    exactly 0 and a tie is mathematically forced.  If so, the 0.00% observation is a
    consistency check on the implementation, not empirical evidence for the mechanism, and
    must not be presented as the latter.  This is verified directly against the recorded
    best_score_fraction values.
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
OUT = os.path.join(DATA, "robustness_checks.json")


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


def main():
    cn_by_region = {}
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") == "ANNOTATED":
                cn_by_region[r["region_id"]] = (r["seqc2_class"], r["chrom"])

    rows = []
    with open(TOPO) as fh:
        for ln in fh:
            u = json.loads(ln)
            rid = u.get("region_id")
            hit = cn_by_region.get(rid)
            if hit is None or u.get("unit_status") != "ranked":
                continue
            cls, chrom = hit
            ttc = int(u.get("total_tree_count") or 0)
            active = set(u.get("active_positions") or [])
            fracs = []
            for c in u.get("af_coverage") or []:
                if c.get("position") not in active:
                    continue
                fr = c.get("fraction")
                if fr:
                    try:
                        fracs.append(Fraction(fr))
                    except (ValueError, ZeroDivisionError):
                        pass
            rows.append(
                {
                    "chrom": chrom,
                    "cn_class": cls,
                    "altered": cls != "neutral",
                    "total_tree_count": ttc,
                    "contested": ttc > 1,
                    "af_broke_tie": ttc > 1 and bool(u.get("best_tree_unique")),
                    "n_af": len(fracs),
                    "all_identical": len(fracs) >= 2 and len(set(fracs)) == 1,
                    "best_score_fraction": u.get("best_score_fraction"),
                }
            )

    contested = [r for r in rows if r["contested"]]

    # ---- Check A: chromosome as confounder ----
    neutral_by_chrom = Counter(r["chrom"] for r in rows if not r["altered"])
    total_neutral = sum(neutral_by_chrom.values())
    top4 = neutral_by_chrom.most_common(4)
    concentration = {
        "neutral_units_total": total_neutral,
        "top4_chromosomes": [
            {"chrom": c, "n": n, "percent_of_neutral": round(100 * n / total_neutral, 2)}
            for c, n in top4
        ],
        "top4_combined_percent": round(
            100 * sum(n for _, n in top4) / total_neutral, 2
        ),
        "chromosomes_with_zero_neutral": sorted(
            {r["chrom"] for r in rows} - set(neutral_by_chrom)
        ),
    }

    # stratified by chromosome across all chromosomes that have both groups
    strata = defaultdict(lambda: {"alt": [0, 0], "neu": [0, 0]})
    for r in contested:
        s = strata[r["chrom"]]
        side = "alt" if r["altered"] else "neu"
        s[side][0 if r["af_broke_tie"] else 1] += 1

    usable_chroms = []
    tables = []
    detail = {}
    for chrom, s in sorted(strata.items()):
        aa, bb = s["alt"]
        cc, dd = s["neu"]
        if (aa + bb) > 0 and (cc + dd) > 0:
            usable_chroms.append(chrom)
            tables.append((aa, bb, cc, dd))
            detail[chrom] = {
                "altered_n": aa + bb,
                "altered_rate_percent": round(100 * aa / (aa + bb), 2),
                "neutral_n": cc + dd,
                "neutral_rate_percent": round(100 * cc / (cc + dd), 2),
                "difference_pp": round(
                    100 * aa / (aa + bb) - 100 * cc / (cc + dd), 2
                ),
            }

    within_chrom = {
        "logic": "restrict the contrast to chromosomes containing both CN-neutral and CN-altered contested units, so chromosome identity cannot drive the difference",
        "chromosomes_usable": usable_chroms,
        "per_chromosome": detail,
        "cmh_stratified_by_chromosome": cmh(tables),
    }

    # pooled comparison restricted to those chromosomes only
    sub = [r for r in contested if r["chrom"] in set(usable_chroms)]
    sa = sum(1 for r in sub if r["altered"] and r["af_broke_tie"])
    sb = sum(1 for r in sub if r["altered"] and not r["af_broke_tie"])
    sc = sum(1 for r in sub if not r["altered"] and r["af_broke_tie"])
    sd = sum(1 for r in sub if not r["altered"] and not r["af_broke_tie"])
    pooled_or, pooled_p = (
        fisher_exact([[sa, sb], [sc, sd]]) if (sa + sb) and (sc + sd) else (None, None)
    )
    within_chrom["pooled_on_usable_chromosomes"] = {
        "altered_n": sa + sb,
        "altered_rate_percent": round(100 * sa / (sa + sb), 2) if (sa + sb) else None,
        "neutral_n": sc + sd,
        "neutral_rate_percent": round(100 * sc / (sc + sd), 2) if (sc + sd) else None,
        "fisher_odds_ratio": round(pooled_or, 4) if pooled_or else None,
        "fisher_p_value": pooled_p,
    }

    # ---- Check B: tautology test ----
    ident = [r for r in rows if r["all_identical"]]
    scores = Counter(str(r["best_score_fraction"]) for r in ident)
    all_zero = all(
        r["best_score_fraction"] in ("0", "0/1", 0) or
        (isinstance(r["best_score_fraction"], str) and Fraction(r["best_score_fraction"]) == 0)
        for r in ident
        if r["best_score_fraction"] is not None
    )
    non_ident = [r for r in rows if r["n_af"] >= 2 and not r["all_identical"]]
    nz = sum(
        1
        for r in non_ident
        if r["best_score_fraction"] is not None
        and Fraction(str(r["best_score_fraction"])) != 0
    )

    tautology = {
        "argument": "s(p,c) = sum_{i in p}(AF_i - AF_j(c)); if all AF equal v every term is 0, so every candidate tree scores 0 and a tie is forced",
        "identical_af_units": len(ident),
        "their_best_score_values": dict(scores.most_common(10)),
        "all_identical_af_units_score_zero": all_zero,
        "verdict": (
            "TAUTOLOGY CONFIRMED - the 0.00% tie-breaking rate among identical-AF units is "
            "mathematically forced by the score definition, not empirical evidence. It is a "
            "valid implementation consistency check only."
            if all_zero
            else "NOT a pure tautology - some identical-AF units carry a non-zero score, so the implementation does something beyond the stated formula; investigate."
        ),
        "distinct_af_units_with_nonzero_score": nz,
        "distinct_af_units_total": len(non_ident),
        "implication": "Prediction 1 (CN-neutral ground more often yields identical AF: 20.80% vs 11.78%) remains a genuine empirical observation and is the load-bearing half of the mechanism; Prediction 2 must be demoted to a consistency check.",
    }

    out = {
        "generated_by": os.path.basename(__file__),
        "check_a_chromosome_confounding": {
            "neutral_group_concentration": concentration,
            "within_chromosome_contrast": within_chrom,
        },
        "check_b_tautology": tautology,
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print("=== Check A: neutral group concentration ===")
    print(json.dumps(concentration, indent=1, ensure_ascii=False))
    print("\n=== within-chromosome contrast ===")
    for c, v in detail.items():
        print(
            f"  {c:7s} alt n={v['altered_n']:5d} {v['altered_rate_percent']:6.2f}%  "
            f"neu n={v['neutral_n']:4d} {v['neutral_rate_percent']:6.2f}%  "
            f"diff={v['difference_pp']:+6.2f}pp"
        )
    print(f"\nCMH by chromosome: {within_chrom['cmh_stratified_by_chromosome']}")
    print(f"pooled on usable chroms: {within_chrom['pooled_on_usable_chromosomes']}")
    print(f"\n=== Check B: tautology ===")
    print(f"identical-AF units: {len(ident)}")
    print(f"their best_score values: {dict(scores.most_common(5))}")
    print(f"verdict: {tautology['verdict']}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
