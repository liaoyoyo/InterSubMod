#!/usr/bin/env python3
"""Can read-level co-occurrence constrain what VAF-based subclone methods infer?

VAF-based reconstruction (PyClone / SciClone / DPClust / CITUP and friends)
works from *marginal* variant allele fractions.  It never observes which
mutations sit on the same DNA molecule, so it must infer co-occurrence from the
fractions themselves, typically via the pigeonhole / sum rule: if f_A >= f_B the
pair is allowed to be nested with A ancestral, and clusters of similar f are
assumed to be one clone.

Long reads give the missing observable directly.  For a pair of mutations
sharing reads, the observed haplotype patterns settle the relationship without
any frequency argument at all:

    {RR, AR, AA}          nested, A before B
    {RR, RA, AA}          nested, B before A
    {RR, AR, RA}          mutually exclusive - different branches
    all four gametes      incompatible with a perfect phylogeny (recurrence or
                          loss required)

This script measures how often that structural verdict and the VAF-based
verdict would disagree, which is the concrete value the linkage layer adds:

  * mutually exclusive pairs carrying a large AF gap -- a VAF-only method sees
    an ancestor/descendant ordering that the molecules rule out
  * nested pairs violating AF monotonicity -- the ordering VAF would infer is
    the reverse of the observed one
  * four-gamete violations -- invisible to marginal VAF entirely

Everything is stratified by copy-number state, because the whole point is that
the structural layer keeps working where the frequency layer degrades.
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
OUT = os.path.join(DATA, "linkage_as_constraint_on_vaf.json")

MIN_DEPTH = 10
BIG_GAP = 0.20  # AF gap a VAF-only method would read as a real ordering


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
            if int(u.get("active_bit_count") or 0) != 2:
                continue  # k=2 isolates one pairwise relationship per unit

            af = {}
            for c in u.get("af_coverage") or []:
                p, a, r_ = c.get("position"), c.get("alt_reads"), c.get("ref_reads")
                if p is not None and a is not None and r_ is not None and (a + r_) >= MIN_DEPTH:
                    af[p] = a / (a + r_)
            pos = u.get("active_positions") or []
            if len(pos) != 2 or not all(p in af for p in pos):
                continue

            morph = u.get("representative_best_morphology") or {}
            has_sister = bool(morph.get("has_sister"))
            has_direct = bool(morph.get("has_direct"))
            recurrence = bool(u.get("recurrence_required"))

            if recurrence:
                rel = "four_gamete_violation"
            elif has_sister and not has_direct:
                rel = "mutually_exclusive"
            elif has_direct and not has_sister:
                rel = "nested"
            else:
                rel = "mixed_or_unresolved"

            a1, a2 = af[pos[0]], af[pos[1]]
            gap = abs(a1 - a2)

            # what the representative tree says about ordering, when nested
            order_ok = None
            if rel == "nested":
                # parent acquires first; find the deeper acquisition
                deep = [
                    e
                    for e in (u.get("representative_best_edges") or [])
                    if e.get("parent_label") != "ROOT"
                ]
                if deep:
                    child_pos = deep[0].get("acquired_position")
                    parent_pos = [p for p in pos if p != child_pos]
                    if child_pos in af and parent_pos and parent_pos[0] in af:
                        order_ok = af[parent_pos[0]] >= af[child_pos]

            rows.append(
                {
                    "cn_class": cls,
                    "altered": cls != "neutral",
                    "relation": rel,
                    "af_gap": gap,
                    "af_monotonic_ok": order_ok,
                    "determined": int(u.get("total_tree_count") or 0) == 1
                    or bool(u.get("best_tree_unique")),
                }
            )

    n = len(rows)

    def by(subset, label):
        if not subset:
            return {"label": label, "n": 0}
        rel = Counter(r["relation"] for r in subset)
        excl = [r for r in subset if r["relation"] == "mutually_exclusive"]
        nested = [r for r in subset if r["relation"] == "nested"]
        # a VAF-only method would call a large-gap pair an ordering; the
        # molecules say these are on different branches
        misread = [r for r in excl if r["af_gap"] >= BIG_GAP]
        nested_checked = [r for r in nested if r["af_monotonic_ok"] is not None]
        nested_violating = [r for r in nested_checked if r["af_monotonic_ok"] is False]
        gaps = sorted(r["af_gap"] for r in excl)
        return {
            "label": label,
            "n_pairs": len(subset),
            "relation_counts": dict(rel),
            "relation_percent": {
                k: round(100 * v / len(subset), 2) for k, v in rel.items()
            },
            "mutually_exclusive_with_large_af_gap": wilson(len(misread), len(excl))
            if excl
            else None,
            "median_af_gap_in_exclusive_pairs": round(gaps[len(gaps) // 2], 4)
            if gaps
            else None,
            "nested_pairs_checked": len(nested_checked),
            "nested_violating_af_monotonicity": wilson(
                len(nested_violating), len(nested_checked)
            )
            if nested_checked
            else None,
            "four_gamete_violations": rel.get("four_gamete_violation", 0),
        }

    out = {
        "generated_by": os.path.basename(__file__),
        "scope": "HCC1395 chr1-22, canonical ranked units with exactly k=2 active mutations and read depth >= %d at both sites" % MIN_DEPTH,
        "why_k2": "a two-mutation unit isolates exactly one pairwise relationship, so the structural verdict and the VAF verdict can be compared without confounding from tree size",
        "large_af_gap_threshold": BIG_GAP,
        "what_vaf_only_methods_see": "marginal AF per mutation; co-occurrence must be inferred from the fractions (pigeonhole / sum rule)",
        "what_long_reads_add": "the haplotype pattern set, which settles nested vs mutually exclusive vs four-gamete violation without any frequency argument",
        "all": by(rows, "all k=2 pairs"),
        "cn_neutral": by([r for r in rows if not r["altered"]], "CN-neutral"),
        "cn_altered": by([r for r in rows if r["altered"]], "CN-altered"),
        "by_cn_class": {
            c: by([r for r in rows if r["cn_class"] == c], c)
            for c in sorted({r["cn_class"] for r in rows})
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    for key in ["all", "cn_neutral", "cn_altered"]:
        p = out[key]
        if not p.get("n_pairs"):
            continue
        print(f"\n=== {p['label']} (n={p['n_pairs']}) ===")
        print("  relation:", p["relation_percent"])
        me = p["mutually_exclusive_with_large_af_gap"]
        if me:
            print(
                f"  互斥但 AF 差距 >= {BIG_GAP}: {me['count']}/{me['n']} = {me['percent']}% ({me['lo']}-{me['hi']})"
            )
            print(f"  互斥對的 AF 差距中位: {p['median_af_gap_in_exclusive_pairs']}")
        nv = p["nested_violating_af_monotonicity"]
        if nv:
            print(
                f"  巢狀但違反 AF 單調: {nv['count']}/{nv['n']} = {nv['percent']}% ({nv['lo']}-{nv['hi']})"
            )
        print(f"  四配子違反: {p['four_gamete_violations']}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
