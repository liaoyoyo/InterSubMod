#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aggregate the authority tumor-REF control table over all M1 sites.

Answers, with explicit denominators:
  1. how many M1 sites have a testable tumor-REF counterpart at all;
  2. how many have a stable joint (ALT+REF) clustering;
  3. how many have the methyl partition aligned with the ALT/REF allele axis;
  4. how those interact with the HP / technical confound flags.

Also cross-validates the derived methyl_K (from the read-assignment stream) against the
authority `coarse_ng` on the same (dataset, chrom, pos) key -- a mismatch is reported,
never silently reconciled.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
import sys

CONTROL = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
    "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel/"
    "all_ssnv_tumor_ref_control_site_results.tsv.gz"
)


def truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes")


def as_int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--control", default=CONTROL)
    ap.add_argument("--m1-table", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-aligned-tsv", required=True)
    a = ap.parse_args()

    # derived K from the read-assignment stream, keyed by site
    derived = {}
    with open(a.m1_table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            derived[(row["dataset"], row["chrom"], row["pos"])] = (
                as_int(row["methyl_K"]), row["hp_presence"])

    n = 0
    per_dataset = collections.Counter()
    ref_status = collections.Counter()
    ref_not_testable = collections.Counter()
    joint_status = collections.Counter()
    joint_not_testable = collections.Counter()
    allele_testable = 0
    allele_aligned = 0
    evidence_tier = collections.Counter()
    hp_conf = collections.Counter()
    tech_conf = collections.Counter()
    residual = collections.Counter()
    aligned_by_dataset = collections.Counter()
    testable_by_dataset = collections.Counter()
    ref_reads_hist = collections.Counter()
    k_match = k_mismatch = k_absent = 0
    hp_presence_of_aligned = collections.Counter()
    aligned_rows = []

    keep_cols = ["dataset", "chrom", "pos", "ref", "alt", "truth_label",
                 "n_tumor_alt", "n_tumor_ref", "coarse_ng", "cluster_sizes",
                 "joint_allele_v", "joint_allele_p_perm", "joint_allele_n",
                 "ref_status", "joint_status", "evidence_tier",
                 "hp_axis_confound", "technical_axis_confound",
                 "residual_unexplained_multigroup",
                 "phase_anchored_robust_epigenetic_candidate",
                 "alt_hp_family_counts"]

    with gzip.open(a.control, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            n += 1
            ds = row["dataset"]
            per_dataset[ds] += 1
            ref_status[row.get("ref_status", "")] += 1
            if not truthy(row.get("ref_evaluable")):
                ref_not_testable[row.get("ref_not_testable_reason", "") or "UNSPECIFIED"] += 1
            joint_status[row.get("joint_status", "")] += 1
            if not truthy(row.get("joint_evaluable")):
                joint_not_testable[row.get("joint_not_testable_reason", "") or "UNSPECIFIED"] += 1
            evidence_tier[row.get("evidence_tier", "")] += 1
            hp_conf[str(row.get("hp_axis_confound", ""))] += 1
            tech_conf[str(row.get("technical_axis_confound", ""))] += 1
            residual[str(row.get("residual_unexplained_multigroup", ""))] += 1

            nref = as_int(row.get("n_tumor_ref"))
            if nref is None:
                ref_reads_hist["NA"] += 1
            elif nref == 0:
                ref_reads_hist["0"] += 1
            elif nref < 3:
                ref_reads_hist["1-2"] += 1
            elif nref < 6:
                ref_reads_hist["3-5"] += 1
            elif nref < 11:
                ref_reads_hist["6-10"] += 1
            elif nref < 21:
                ref_reads_hist["11-20"] += 1
            else:
                ref_reads_hist[">=21"] += 1

            if truthy(row.get("joint_allele_testable")):
                allele_testable += 1
                testable_by_dataset[ds] += 1
            if truthy(row.get("joint_allele_axis_aligned")):
                allele_aligned += 1
                aligned_by_dataset[ds] += 1
                key = (ds, row["chrom"], row["pos"])
                hp_presence_of_aligned[derived.get(key, (None, "NOT_IN_M1_TABLE"))[1]] += 1
                aligned_rows.append({c: row.get(c, "") for c in keep_cols})

            # cross-validate derived K vs authority coarse_ng
            key = (ds, row["chrom"], row["pos"])
            dk = derived.get(key)
            ang = as_int(row.get("coarse_ng"))
            if dk is None:
                k_absent += 1
            elif dk[0] == ang:
                k_match += 1
            else:
                k_mismatch += 1

    os.makedirs(os.path.dirname(a.out_aligned_tsv), exist_ok=True)
    with open(a.out_aligned_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keep_cols, delimiter="\t")
        w.writeheader()
        for r in sorted(aligned_rows, key=lambda r: (r["dataset"], r["chrom"], as_int(r["pos"]) or 0)):
            w.writerow(r)

    sha = hashlib.sha256()
    with open(a.out_aligned_tsv, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha.update(chunk)

    out = {
        "schema_name": "intersubmod.m1_locus_annotation.alt_ref_axis_aggregate",
        "schema_version": "1.0.0",
        "control_source": a.control,
        "control_bytes": os.path.getsize(a.control),
        "m1_table": a.m1_table,
        "rows_read": n,
        "per_dataset": dict(per_dataset),
        "tumor_ref_read_count_bins": dict(ref_reads_hist),
        "ref_status": dict(ref_status),
        "ref_not_testable_reason_top": dict(ref_not_testable.most_common(10)),
        "joint_status": dict(joint_status),
        "joint_not_testable_reason_top": dict(joint_not_testable.most_common(10)),
        "joint_allele_testable": allele_testable,
        "joint_allele_axis_aligned": allele_aligned,
        "joint_allele_axis_aligned_by_dataset": dict(aligned_by_dataset),
        "joint_allele_testable_by_dataset": dict(testable_by_dataset),
        "aligned_hp_presence_from_m1_table": dict(hp_presence_of_aligned),
        "evidence_tier": dict(evidence_tier),
        "hp_axis_confound": dict(hp_conf),
        "technical_axis_confound": dict(tech_conf),
        "residual_unexplained_multigroup": dict(residual),
        "derived_K_vs_authority_coarse_ng": {
            "match": k_match, "mismatch": k_mismatch, "site_absent_from_m1_table": k_absent},
        "out_aligned_tsv": a.out_aligned_tsv,
        "out_aligned_tsv_sha256": sha.hexdigest(),
        "guardrails_from_authority_summary": {
            "ref_nonreplication": "Tumor-REF nonreplication does not prove a subclone.",
            "ref_replication": "Tumor-REF replication weakens ALT-specificity.",
            "joint_allele_association": ("posthoc Cramer's V with 499 label permutations "
                                         "after stable joint clustering"),
        },
    }
    with open(a.out_json, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
