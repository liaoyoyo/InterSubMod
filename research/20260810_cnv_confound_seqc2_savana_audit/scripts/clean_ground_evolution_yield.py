#!/usr/bin/env python3
"""If copy-number-altered ground is discarded, what evolutionary signal is left?

"177 CN-neutral units with k>=2" is not yet an answer, because units differ in
how much lineage they can express:

  k = 1  no tree is definable at all
  k = 2  one ordering question only (A before B, B before A, or a branch)
  k >= 3 multi-step lineage: chains, branches, and their combination

and because a unit only contributes if its tree is actually determined:

  structure-unique  the enumeration returned one tree; read-AF never used
  AF-resolved       several candidates, read-AF picked one.  On CN-neutral
                    ground this selection is NOT copy-number confounded, so it
                    is usable here even though the same step is unreliable on
                    CN-altered ground
  still tied        no single tree; contributes shape/candidate-set only

This script counts the yield along those two axes, and converts it into the
quantity that actually matters for a lineage claim: how many parent -> child
mutation-acquisition edges are determined on copy-number-clean ground.
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
CENSUS = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1/HCC1395.census.jsonl"
)
OUT = os.path.join(DATA, "clean_ground_evolution_yield.json")


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
    shape = {}
    with open(CENSUS) as fh:
        for ln in fh:
            r = json.loads(ln)
            shape[r["region_id"]] = list((r.get("coarse_class_tree_counts") or {}).keys())

    cn, chrom = {}, {}
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") == "ANNOTATED":
                cn[r["region_id"]] = r["seqc2_class"]
                chrom[r["region_id"]] = r["chrom"]

    units = []
    with open(TOPO) as fh:
        for ln in fh:
            u = json.loads(ln)
            rid = u.get("region_id")
            if rid not in cn or u.get("unit_status") != "ranked":
                continue
            k = int(u.get("active_bit_count") or 0)
            ttc = int(u.get("total_tree_count") or 0)
            uniq = bool(u.get("best_tree_unique"))
            edges = u.get("representative_best_edges") or []
            # edges from ROOT are the first acquisition; deeper edges express order
            deep_edges = [e for e in edges if e.get("parent_label") != "ROOT"]
            hidden = [
                e
                for e in edges
                if str(e.get("child_label", "")).startswith("H_")
                or str(e.get("parent_label", "")).startswith("H_")
            ]
            if ttc == 1:
                state = "structure_unique"
            elif uniq:
                state = "af_resolved"
            else:
                state = "still_tied"
            units.append(
                {
                    "rid": rid,
                    "chrom": chrom[rid],
                    "cls": cn[rid],
                    "altered": cn[rid] != "neutral",
                    "k": k,
                    "ttc": ttc,
                    "state": state,
                    "determined": state in ("structure_unique", "af_resolved"),
                    "n_edges": len(edges),
                    "n_deep_edges": len(deep_edges),
                    "n_hidden_edges": len(hidden),
                    "recurrence": bool(u.get("recurrence_required")),
                    "shapes": shape.get(rid, []),
                    "branching": any("Sister" in s for s in shape.get(rid, [])),
                }
            )

    def kband(k):
        return "k=1" if k <= 1 else "k=2" if k == 2 else "k=3" if k == 3 else "k>=4"

    def profile(sub, label):
        n = len(sub)
        if n == 0:
            return {"label": label, "n": 0}
        det = [u for u in sub if u["determined"]]
        su = [u for u in sub if u["state"] == "structure_unique"]
        af = [u for u in sub if u["state"] == "af_resolved"]
        tied = [u for u in sub if u["state"] == "still_tied"]
        return {
            "label": label,
            "n_units": n,
            "structure_unique": wilson(len(su), n),
            "af_resolved": wilson(len(af), n),
            "still_tied": wilson(len(tied), n),
            "determined_total": wilson(len(det), n),
            "determined_edges": sum(u["n_edges"] for u in det),
            "determined_deep_edges": sum(u["n_deep_edges"] for u in det),
            "determined_with_branching": wilson(
                sum(1 for u in det if u["branching"]), len(det)
            )
            if det
            else None,
            "determined_with_hidden_ancestor": sum(
                1 for u in det if u["n_hidden_edges"] > 0
            ),
            "chromosomes_covered": len({u["chrom"] for u in det}),
        }

    neu = [u for u in units if not u["altered"]]
    alt = [u for u in units if u["altered"]]

    out = {
        "generated_by": os.path.basename(__file__),
        "definitions": {
            "structure_unique": "total_tree_count == 1; read-AF not involved; copy-number robust",
            "af_resolved": "several candidates, read-AF picked one; on CN-neutral ground this selection is not copy-number confounded",
            "still_tied": "no single tree; candidate set / shape only",
            "determined": "structure_unique OR af_resolved",
            "deep_edge": "an edge whose parent is not ROOT, i.e. it expresses an ordering between two acquired mutations rather than a first acquisition",
        },
        "cn_neutral": {
            "all": profile(neu, "CN-neutral, all ranked"),
            "k_ge_2": profile([u for u in neu if u["k"] >= 2], "CN-neutral, k>=2"),
            "k_ge_3": profile([u for u in neu if u["k"] >= 3], "CN-neutral, k>=3"),
            "by_k": {
                b: profile([u for u in neu if kband(u["k"]) == b], f"CN-neutral {b}")
                for b in ["k=1", "k=2", "k=3", "k>=4"]
            },
        },
        "cn_altered_reference": {
            "k_ge_2": profile([u for u in alt if u["k"] >= 2], "CN-altered, k>=2"),
            "k_ge_3": profile([u for u in alt if u["k"] >= 3], "CN-altered, k>=3"),
        },
        "clean_k_ge_3_detail": {
            "units": [
                {
                    "region_id": u["rid"],
                    "chrom": u["chrom"],
                    "k": u["k"],
                    "state": u["state"],
                    "candidate_trees": u["ttc"],
                    "edges": u["n_edges"],
                    "branching": u["branching"],
                    "shapes": u["shapes"],
                }
                for u in sorted(
                    [x for x in neu if x["k"] >= 3], key=lambda x: (-x["k"], x["chrom"])
                )
            ],
            "chrom_distribution": dict(
                Counter(u["chrom"] for u in neu if u["k"] >= 3).most_common()
            ),
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    for key in ["all", "k_ge_2", "k_ge_3"]:
        p = out["cn_neutral"][key]
        print(f"\n=== CN-neutral {key} (n={p.get('n_units')}) ===")
        for f in ["structure_unique", "af_resolved", "still_tied", "determined_total"]:
            v = p.get(f)
            if v:
                print(f"  {f:20s} {v['count']:4d}/{v['n']:4d} = {v['percent']:6.2f}% ({v['lo']}-{v['hi']})")
        print(f"  determined edges     : {p.get('determined_edges')}")
        print(f"  determined deep edges: {p.get('determined_deep_edges')}")
        print(f"  chromosomes covered  : {p.get('chromosomes_covered')}")
    print("\n=== by k (CN-neutral) ===")
    for b, p in out["cn_neutral"]["by_k"].items():
        if p.get("n_units"):
            print(
                f"  {b:6s} n={p['n_units']:4d} determined={p['determined_total']['percent']:6.2f}% "
                f"edges={p['determined_edges']:5d} deep={p['determined_deep_edges']:4d}"
            )
    print("\n=== CN-altered reference ===")
    for key in ["k_ge_2", "k_ge_3"]:
        p = out["cn_altered_reference"][key]
        print(
            f"  {key:8s} n={p['n_units']:5d} determined={p['determined_total']['percent']:6.2f}% "
            f"edges={p['determined_edges']:6d} deep={p['determined_deep_edges']:5d}"
        )
    print("\nclean k>=3 chrom:", out["clean_k_ge_3_detail"]["chrom_distribution"])
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
