#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Decompose the methylation screen into the two views asked for.

View 1  every LongPhase-S PASS sSNV -> not evaluable / evaluable-no-difference /
        evaluable-with-stable-methyl-difference (M1).

View 2  the M1 sites only -> which measured axis the methyl split runs along:
          HP axis        the split follows the germline haplotype (HP1 vs HP2)
          allele axis    the split follows ALT vs REF at the focal site
          technical axis strand / start / end / length / MAPQ / CpG-called
          residual       no measured axis explains it -> within-ALT substructure

The axes overlap, so both the raw per-axis counts and a precedence-ordered mutually
exclusive partition are emitted; the page must never show only one of the two.
Precedence puts HP first because a haplotype-aligned split is the strongest confound:
it is indistinguishable from germline allele-specific methylation.
"""
import argparse
import collections
import csv
import gzip
import json
import os
import sys

BASE = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
        "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation")
ALL_SITES = f"{BASE}/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz"
ALL_SUMMARY = f"{BASE}/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_summary.json"
CONTROL = f"{BASE}/all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel/all_ssnv_tumor_ref_control_site_results.tsv.gz"


def truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()
    for path in (ALL_SITES, ALL_SUMMARY, CONTROL):
        if not os.path.exists(path):
            sys.exit(f"FAIL CLOSED: missing authority artifact {path}")

    summary = json.load(open(ALL_SUMMARY))
    pooled = summary["pooled_site_weighted"]

    # ---------- view 1: whole screen, per dataset
    v1 = collections.defaultdict(collections.Counter)
    total = 0
    with gzip.open(ALL_SITES, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            total += 1
            ds = row["dataset"]
            if not truthy(row.get("m1_evaluable")):
                v1[ds]["not_evaluable"] += 1
            elif truthy(row.get("stable_null_multigroup")):
                v1[ds]["m1_difference"] += 1
            else:
                v1[ds]["no_difference"] += 1
    if total != int(pooled["n_sites"]):
        sys.exit(f"FAIL CLOSED: streamed {total} sites, summary says {pooled['n_sites']}")
    cohort_v1 = collections.Counter()
    for counts in v1.values():
        cohort_v1.update(counts)
    if cohort_v1["m1_difference"] != int(pooled["n_stable_null_multigroup"]):
        sys.exit("FAIL CLOSED: M1 count disagrees with authority summary")
    if cohort_v1["not_evaluable"] != total - int(pooled["n_evaluable"]):
        sys.exit("FAIL CLOSED: evaluable count disagrees with authority summary")

    # ---------- view 2: the M1 sites, by which axis explains the split
    v2 = collections.defaultdict(collections.Counter)
    overlap = collections.Counter()
    m1_rows = 0
    with gzip.open(CONTROL, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            m1_rows += 1
            ds = row["dataset"]
            hp = truthy(row.get("hp_family_aligned")) or truthy(row.get("hp_exact_aligned")) \
                or truthy(row.get("hp_axis_confound"))
            allele = truthy(row.get("joint_allele_axis_aligned"))
            tech = truthy(row.get("technical_axis_confound"))
            residual = truthy(row.get("residual_unexplained_multigroup"))
            overlap[(hp, allele, tech)] += 1
            # precedence: HP > allele > technical > residual
            if hp:
                bucket = "hp_axis"
            elif allele:
                bucket = "allele_axis"
            elif tech:
                bucket = "technical_axis"
            elif residual:
                bucket = "residual_within_alt"
            else:
                bucket = "unclassified"
            v2[ds][bucket] += 1
            v2[ds]["_raw_hp"] += int(hp)
            v2[ds]["_raw_allele"] += int(allele)
            v2[ds]["_raw_tech"] += int(tech)
            v2[ds]["_raw_residual"] += int(residual)
    if m1_rows != cohort_v1["m1_difference"]:
        sys.exit(f"FAIL CLOSED: control table has {m1_rows} rows, screen says "
                 f"{cohort_v1['m1_difference']}")
    cohort_v2 = collections.Counter()
    for counts in v2.values():
        cohort_v2.update(counts)
    partition = sum(cohort_v2[k] for k in
                    ("hp_axis", "allele_axis", "technical_axis",
                     "residual_within_alt", "unclassified"))
    if partition != m1_rows:
        sys.exit("FAIL CLOSED: axis partition does not cover every M1 site")

    out = {
        "schema_name": "intersubmod.methyl_axis_decomposition.v1",
        "schema_version": "1.0.0",
        "sources": {"all_sites": ALL_SITES, "all_summary": ALL_SUMMARY, "control": CONTROL},
        "view1_definition": {
            "not_evaluable": "m1_evaluable=false（ALT read 不足或距離矩陣不完整）",
            "no_difference": "可評估但未通過穩定多群門檻",
            "m1_difference": "stable_null_multigroup=true（M1）",
            "contract": summary["clustering_contract"]["m1_stability_gate_contract"],
        },
        "view1_cohort": dict(cohort_v1),
        "view1_per_dataset": {k: dict(v) for k, v in v1.items()},
        "view1_authority_check": {
            "n_sites": int(pooled["n_sites"]),
            "n_evaluable": int(pooled["n_evaluable"]),
            "n_stable_null_multigroup": int(pooled["n_stable_null_multigroup"]),
        },
        "view2_definition": {
            "hp_axis": "分群沿 germline haplotype（HP1 vs HP2）— 與 germline ASM 不可分辨",
            "allele_axis": "分群沿焦點位點的 ALT vs REF",
            "technical_axis": "分群沿 strand/start/end/length/MAPQ/CpG-called 之一",
            "residual_within_alt": "已測軸皆無法解釋 — 單 HP 的 ALT read 內部殘餘結構",
            "precedence": "hp_axis > allele_axis > technical_axis > residual_within_alt",
        },
        "view2_cohort_exclusive": {k: cohort_v2[k] for k in
                                   ("hp_axis", "allele_axis", "technical_axis",
                                    "residual_within_alt", "unclassified")},
        "view2_cohort_raw_overlapping": {
            "hp_axis": cohort_v2["_raw_hp"], "allele_axis": cohort_v2["_raw_allele"],
            "technical_axis": cohort_v2["_raw_tech"],
            "residual_within_alt": cohort_v2["_raw_residual"]},
        "view2_per_dataset": {k: dict(v) for k, v in v2.items()},
        "view2_overlap_hp_allele_tech": {f"hp={h}|allele={al}|tech={t}": n
                                         for (h, al, t), n in sorted(overlap.items())},
        "claim_ceiling": (
            "M1 是 operational screen yield，不是 subclone prevalence。allele_axis 與 "
            "germline allele-specific methylation 在 ALT 僅落單一 HP 時不可分辨，且 "
            "germline-het null 尚未執行。residual 只代表已測軸無法解釋，不代表已排除未測 confound。"),
    }
    with open(a.out_json, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("view1_per_dataset", "view2_per_dataset", "sources")},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
