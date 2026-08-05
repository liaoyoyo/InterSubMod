#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Measure whether the deterministic representative tree is shape-neutral.

For every TIED_CROSS_TOPOLOGY unit the canonical pipeline still emits ONE
`representative_best_morphology`. That representative is only defensible for display and
joining if it does not systematically favour one coarse geometry over the others -- any
bias would propagate straight into the coarse-topology composition that carries the only
cross-sample signal with exact-permutation support.

Comparison per unit:
  observed  = the coarse class of the representative
  expected  = that unit's own tied-tree class distribution (census coarse_class_tree_counts)

A neutral rule reproduces the tied distribution in aggregate. Any gap is the bias.

Read-only over frozen canonical + census outputs.
"""
import argparse
import collections
import json
import os
import sys

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
CLASSES = ("Single-only", "Sister-only", "Direct-only", "Sister+direct")


def morphology_class(morph):
    if not isinstance(morph, dict):
        return None
    direct = bool(morph.get("has_direct"))
    sister = bool(morph.get("has_sister"))
    if direct and sister:
        return "Sister+direct"
    if direct:
        return "Direct-only"
    if sister:
        return "Sister-only"
    return "Single-only"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    observed = collections.Counter()          # representative's class, cross units only
    expected = collections.Counter()          # share-weighted tied distribution
    observed_same = collections.Counter()     # sanity: TIED_SAME units
    per_sample_gap = {}
    units = 0
    missing_census = 0
    depth_first_agreement = collections.Counter()

    for sample in SAMPLES:
        census_path = os.path.join(CENSUS, f"{sample}.census.jsonl")
        topo_path = os.path.join(TOPOLOGY, sample, f"{sample}.topology.jsonl")
        if not (os.path.exists(census_path) and os.path.exists(topo_path)):
            sys.exit(f"ERROR: missing inputs for {sample}")

        census = {}
        with open(census_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                census[rec["region_id"]] = rec

        s_obs = collections.Counter()
        s_exp = collections.Counter()
        with open(topo_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                rid = rec.get("region_id")
                cen = census.get(rid)
                if cen is None:
                    continue
                cls = cen["resolution_class"]
                rep = morphology_class(rec.get("representative_best_morphology"))
                if rep is None:
                    missing_census += 1
                    continue
                if cls == "TIED_SAME_TOPOLOGY":
                    observed_same[rep] += 1
                if cls != "TIED_CROSS_TOPOLOGY":
                    continue
                units += 1
                observed[rep] += 1
                s_obs[rep] += 1
                counts = cen.get("coarse_class_tree_counts") or {}
                total = sum(int(v) for v in counts.values())
                if total:
                    for key, value in counts.items():
                        share = int(value) / total
                        expected[key] += share
                        s_exp[key] += share
                # would a "deepest chain" rule agree with the representative?
                # Direct-only / Sister+direct are the chain-bearing classes.
                if counts:
                    deepest = "Sister+direct" if "Sister+direct" in counts else (
                        "Direct-only" if "Direct-only" in counts else sorted(counts)[0])
                    depth_first_agreement[
                        "agrees_with_current_rep" if deepest == rep
                        else "differs_from_current_rep"] += 1

        denom_o = sum(s_obs.values())
        denom_e = sum(s_exp.values())
        if denom_o and denom_e:
            per_sample_gap[sample] = {
                cls: round(s_obs.get(cls, 0) / denom_o - s_exp.get(cls, 0) / denom_e, 6)
                for cls in CLASSES
            }

    tot_o = sum(observed.values())
    tot_e = sum(expected.values())
    out = {
        "schema_name": "intersubmod.representative_shape_bias.v1",
        "schema_version": "1.0.0",
        "cross_units_evaluated": units,
        "representative_class_counts": dict(observed),
        "representative_class_share": {
            c: round(observed.get(c, 0) / tot_o, 6) for c in CLASSES} if tot_o else {},
        "tied_tree_class_share_expected": {
            c: round(expected.get(c, 0) / tot_e, 6) for c in CLASSES} if tot_e else {},
        "bias_representative_minus_expected": {
            c: round(observed.get(c, 0) / tot_o - expected.get(c, 0) / tot_e, 6)
            for c in CLASSES} if tot_o and tot_e else {},
        "per_sample_bias": per_sample_gap,
        "tied_same_representative_classes": dict(observed_same),
        "deepest_chain_rule_vs_current_representative": dict(depth_first_agreement),
        "records_without_morphology": missing_census,
        "interpretation_guard": (
            "expected = each unit's own tied-tree coarse distribution, share-weighted so "
            "every unit counts once. A shape-neutral tie-break reproduces it in aggregate. "
            "A positive gap means the representative over-selects that geometry."),
    }
    with open(a.out_json, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
