#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test whether the AF-tied units are tied because the read-AF axis carries no information.

The read-AF score is a whole-tree sum of (newly acquired mutation AF - ancestor mutation
AF) terms. If every active site in a unit shares the same AF, every term is zero, the
score is zero for *every* candidate tree, and the ranking step cannot discriminate at all.
The extreme and most common version of that is AF = 1.0 everywhere: no REF reads survive
at any active site.

Per unit we record:
  af_values_distinct  how many distinct active-site AFs exist (1 => no discrimination)
  all_af_one          every active site has ref_reads == 0
  best_score_zero     the winning whole-tree score is exactly 0
and cross-tabulate against the resolution class.

Read-only over the frozen canonical topology + signature census.
"""
import argparse
import collections
import json
import os
import sys
from fractions import Fraction

TOPOLOGY = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples"
)
CENSUS = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1"
)
SAMPLES = ("HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009",
           "HCC1937", "HCC1954")
CLASSES = ("UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY")


def parse_fraction(text):
    try:
        num, den = str(text).split("/")
        return Fraction(int(num), int(den))
    except (ValueError, ZeroDivisionError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    stats = {c: collections.Counter() for c in CLASSES}
    distinct_hist = {c: collections.Counter() for c in CLASSES}
    per_sample = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    cross_no_info_by_k = collections.Counter()
    cross_total_by_k = collections.Counter()
    examples_informative_cross = []

    for sample in SAMPLES:
        census_path = os.path.join(CENSUS, f"{sample}.census.jsonl")
        topo_path = os.path.join(TOPOLOGY, sample, f"{sample}.topology.jsonl")
        if not (os.path.exists(census_path) and os.path.exists(topo_path)):
            sys.exit(f"ERROR: missing inputs for {sample}")
        cls_by_region = {}
        with open(census_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                cls_by_region[rec["region_id"]] = rec["resolution_class"]

        with open(topo_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                cls = cls_by_region.get(rec.get("region_id"))
                if cls not in CLASSES:
                    continue
                active = set(int(p) for p in (rec.get("active_positions") or []))
                cov = [c for c in (rec.get("af_coverage") or [])
                       if int(c["position"]) in active]
                if not cov:
                    stats[cls]["no_active_coverage"] += 1
                    continue
                afs = []
                all_one = True
                for c in cov:
                    ref = int(c["ref_reads"])
                    alt = int(c["alt_reads"])
                    if ref + alt == 0:
                        afs.append(None)
                        all_one = False
                        continue
                    afs.append(Fraction(alt, ref + alt))
                    if ref != 0:
                        all_one = False
                distinct = len({x for x in afs if x is not None})
                score = parse_fraction(rec.get("best_score_fraction"))
                score_zero = (score == 0) if score is not None else False

                stats[cls]["units"] += 1
                distinct_hist[cls][min(distinct, 5)] += 1
                if distinct <= 1:
                    stats[cls]["af_uninformative_distinct_le_1"] += 1
                if all_one:
                    stats[cls]["all_active_sites_af_1"] += 1
                if score_zero:
                    stats[cls]["best_score_zero"] += 1
                if distinct <= 1 and all_one:
                    stats[cls]["af1_and_uninformative"] += 1
                per_sample[sample][cls]["units"] += 1
                if all_one:
                    per_sample[sample][cls]["all_af_1"] += 1

                if cls == "TIED_CROSS_TOPOLOGY":
                    k = int(rec["active_bit_count"])
                    kb = f"k={k}" if k <= 4 else "k>=5"
                    cross_total_by_k[kb] += 1
                    if distinct <= 1:
                        cross_no_info_by_k[kb] += 1
                    elif len(examples_informative_cross) < 5:
                        examples_informative_cross.append({
                            "sample": sample, "region_id": rec["region_id"],
                            "active_k": k,
                            "best_score_fraction": rec.get("best_score_fraction"),
                            "distinct_af": distinct,
                            "af_by_site": [
                                {"pos": c["position"], "ref": c["ref_reads"],
                                 "alt": c["alt_reads"]} for c in cov],
                            "ties": rec.get("best_tree_tie_count"),
                        })

    def rate(cls, key):
        n = stats[cls]["units"]
        return round(stats[cls][key] / n, 6) if n else None

    out = {
        "schema_name": "intersubmod.tie_af_degeneracy.v1",
        "schema_version": "1.0.0",
        "per_class": {
            c: {
                "units": stats[c]["units"],
                "af_uninformative_distinct_le_1": stats[c]["af_uninformative_distinct_le_1"],
                "rate_af_uninformative": rate(c, "af_uninformative_distinct_le_1"),
                "all_active_sites_af_1": stats[c]["all_active_sites_af_1"],
                "rate_all_af_1": rate(c, "all_active_sites_af_1"),
                "best_score_zero": stats[c]["best_score_zero"],
                "rate_best_score_zero": rate(c, "best_score_zero"),
                "af1_and_uninformative": stats[c]["af1_and_uninformative"],
                "no_active_coverage": stats[c]["no_active_coverage"],
            } for c in CLASSES
        },
        "distinct_af_histogram": {
            c: {str(k): v for k, v in sorted(distinct_hist[c].items())} for c in CLASSES},
        "cross_af_uninformative_by_k": {
            kb: {"uninformative": cross_no_info_by_k.get(kb, 0),
                 "total": cross_total_by_k[kb],
                 "rate": round(cross_no_info_by_k.get(kb, 0) / cross_total_by_k[kb], 6)}
            for kb in sorted(cross_total_by_k)},
        "per_sample_all_af_1": {
            s: {c: {"units": v[c]["units"], "all_af_1": v[c]["all_af_1"],
                    "rate": round(v[c]["all_af_1"] / v[c]["units"], 6) if v[c]["units"] else None}
                for c in CLASSES} for s, v in per_sample.items()},
        "examples_cross_with_informative_af": examples_informative_cross,
        "interpretation_guard": (
            "distinct AF <= 1 means the read-AF axis is constant across the unit's active "
            "sites, so every candidate tree scores identically and the ranking step cannot "
            "discriminate. AF = 1 at every active site is the dominant special case: no REF "
            "read survives, which is expected under LOH, high purity, or copy-number loss of "
            "the reference allele. This is an observation about the AF axis, not a CN call."),
    }
    with open(a.out_json, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
