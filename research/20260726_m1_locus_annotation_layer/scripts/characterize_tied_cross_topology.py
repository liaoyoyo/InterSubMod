#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Characterise WHY 11.7421% of ranked units land in TIED_CROSS_TOPOLOGY.

The census stores, per unit, both halves of the ambiguity:

  best_vertex_set_count > 1  -> the data does not determine WHICH minimal set of states
                                (including hidden/Steiner states) to realise
  best_vertex_set_count == 1 -> the state set is determined, but a child can attach to
                                more than one equally-scoring parent

Those are different phenomena with different remedies, so they must never be reported as
one lump. This script splits the 8,449 cross-topology units along that line and against
active k, tie count, and dataset.

Read-only over the frozen census; no authority file is modified.
"""
import argparse
import collections
import glob
import json
import os
import sys

CENSUS = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1"
)
CLASSES = ("UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY")


def bucket_k(k):
    if k <= 1:
        return "k<=1"
    if k <= 4:
        return f"k={k}"
    if k <= 6:
        return "k=5-6"
    if k <= 9:
        return "k=7-9"
    return "k>=10"


def bucket_ties(n):
    if n <= 1:
        return "1"
    if n <= 2:
        return "2"
    if n <= 4:
        return "3-4"
    if n <= 8:
        return "5-8"
    if n <= 16:
        return "9-16"
    return ">16"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=CENSUS)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    total = collections.Counter()
    by_sample = collections.defaultdict(collections.Counter)
    # cross-topology breakdown
    cross_mechanism = collections.Counter()
    cross_by_k = collections.Counter()
    cross_by_ties = collections.Counter()
    cross_coarse_span = collections.Counter()
    cross_sig_count = collections.Counter()
    # denominators per k so we can report a RATE, not just a count
    ranked_by_k = collections.Counter()
    cross_by_k_and_mech = collections.Counter()
    # same comparison for the two resolved classes, to show what differs
    mech_by_class = collections.defaultdict(collections.Counter)
    k_by_class = collections.defaultdict(collections.Counter)
    ties_by_class = collections.defaultdict(collections.Counter)
    vertexset_by_class = collections.defaultdict(collections.Counter)
    examples = []

    files = sorted(glob.glob(os.path.join(a.census, "*.census.jsonl")))
    if not files:
        sys.exit(f"ERROR: no census files under {a.census}")

    for path in files:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                cls = rec["resolution_class"]
                if cls not in CLASSES:
                    continue
                sample = rec["sample"]
                k = int(rec["active_bit_count"])
                ties = int(rec["best_tree_tie_count"])
                bvs = int(rec["best_vertex_set_count"])
                mvs = int(rec["minimum_vertex_set_count"])
                sigs = int(rec["topology_signature_count"])
                coarse_n = int(rec["coarse_class_count"])

                mech = ("VERTEX_SET_AMBIGUOUS" if bvs > 1
                        else "PARENT_CHOICE_ONLY")

                total[cls] += 1
                by_sample[sample][cls] += 1
                ranked_by_k[bucket_k(k)] += 1
                mech_by_class[cls][mech] += 1
                k_by_class[cls][bucket_k(k)] += 1
                ties_by_class[cls][bucket_ties(ties)] += 1
                vertexset_by_class[cls]["bvs=1" if bvs == 1 else
                                        "bvs=2" if bvs == 2 else "bvs>=3"] += 1

                if cls == "TIED_CROSS_TOPOLOGY":
                    cross_mechanism[mech] += 1
                    cross_by_k[bucket_k(k)] += 1
                    cross_by_ties[bucket_ties(ties)] += 1
                    cross_coarse_span[f"coarse_classes={coarse_n}"] += 1
                    cross_sig_count[f"signatures={min(sigs,6)}" if sigs < 6
                                    else "signatures>=6"] += 1
                    cross_by_k_and_mech[(bucket_k(k), mech)] += 1
                    if len(examples) < 6 and sigs >= 2:
                        examples.append({
                            "sample": sample, "region_id": rec["region_id"],
                            "active_k": k, "minimum_vertex_set_count": mvs,
                            "best_vertex_set_count": bvs,
                            "best_tree_tie_count": ties,
                            "topology_signature_count": sigs,
                            "coarse_class_count": coarse_n,
                            "coarse_class_tree_counts": rec.get("coarse_class_tree_counts"),
                            "mechanism": mech,
                        })

    ranked = sum(total.values())
    cross = total["TIED_CROSS_TOPOLOGY"]

    def rate_by_k():
        out = {}
        for kb, denom in sorted(ranked_by_k.items()):
            num = cross_by_k.get(kb, 0)
            out[kb] = {"cross": num, "ranked": denom,
                       "rate": round(num / denom, 6) if denom else None}
        return out

    out = {
        "schema_name": "intersubmod.tied_cross_topology.characterisation",
        "schema_version": "1.0.0",
        "census_root": a.census,
        "ranked_units": ranked,
        "resolution_class_counts": dict(total),
        "cross_topology_units": cross,
        "cross_topology_fraction": round(cross / ranked, 6) if ranked else None,
        "mechanism_split_of_cross": dict(cross_mechanism),
        "mechanism_split_by_class": {c: dict(v) for c, v in mech_by_class.items()},
        "best_vertex_set_count_by_class": {c: dict(v) for c, v in vertexset_by_class.items()},
        "cross_by_active_k": dict(cross_by_k),
        "cross_rate_by_active_k": rate_by_k(),
        "active_k_by_class": {c: dict(v) for c, v in k_by_class.items()},
        "cross_by_tie_count": dict(cross_by_ties),
        "tie_count_by_class": {c: dict(v) for c, v in ties_by_class.items()},
        "cross_signature_multiplicity": dict(cross_sig_count),
        "cross_coarse_class_span": dict(cross_coarse_span),
        "cross_by_k_and_mechanism": {f"{k}|{m}": v
                                     for (k, m), v in sorted(cross_by_k_and_mech.items())},
        "per_sample": {s: dict(v) for s, v in by_sample.items()},
        "per_sample_cross_rate": {
            s: round(v["TIED_CROSS_TOPOLOGY"] / sum(v.values()), 6)
            for s, v in by_sample.items() if sum(v.values())
        },
        "examples": examples,
        "interpretation_guard": (
            "VERTEX_SET_AMBIGUOUS means the minimal state set itself is not determined "
            "(competing hidden/Steiner completions). PARENT_CHOICE_ONLY means the state "
            "set is fixed but at least one child has several equally-scoring parents. "
            "Neither is an error; both are honest non-identifiability under the coverage "
            "axiom and the enumerate-not-optimise rule."),
    }
    with open(a.out_json, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
