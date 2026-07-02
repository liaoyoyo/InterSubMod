#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regenerate multisample_summary.json by aggregating each sample's CURRENT per-sample
topology_per_region.json + candidate_scoring.json + region_gene_annotation.json.

Fills the gap found 2026-07-02 audit: multisample_summary.json had no committed builder,
so the 06-29 aggregate never picked up the 07-01/07-02 C2/C3/D4 topology fix while the
per-sample JSONs (07-02) already had it. This script is pure aggregation (reads existing
per-sample JSONs, no BAM, no C++) so it is safe/fast to rerun any time the per-sample
outputs change.

Sources:
  HCC1395  -> canonical 4axis data dir (../data)
  6 others -> MSROOT/<sample>/
Determinacy fields come from topology (post-fix); gene-annotation counts are
fix-independent (depend on genomic coords, not tree topology).

Usage: python3 gen_multisample_summary.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA4AXIS = os.path.normpath(os.path.join(HERE, "..", "data"))
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"

# preserve original sample order
SAMPLES = ["HCC1395", "COLO829", "H1437", "H2009", "HCC1395_DORADO", "HCC1937", "HCC1954"]


def load(path):
    return json.load(open(path, encoding="utf-8")) if os.path.exists(path) else None


def prefix_get(d, prefix, default=0):
    """determinacy/topology_type keys carry Chinese suffixes, e.g. 'A_determined(單分子向量)'."""
    for k, v in d.items():
        if k.startswith(prefix):
            return v
    return default


def sample_dir(s):
    return DATA4AXIS if s == "HCC1395" else os.path.join(MSROOT, s)


def summarize(s):
    dr = sample_dir(s)
    tp = load(os.path.join(dr, "topology_per_region.json"))
    cs = load(os.path.join(dr, "candidate_scoring.json"))
    ga = load(os.path.join(dr, "region_gene_annotation.json"))
    if tp is None:
        return None
    stats = tp["stats"]
    det = stats["determinacy"]
    tt = stats["topology_type"]
    with_vector = (
        prefix_get(det, "A_determined")
        + prefix_get(det, "A_ambiguous")
        + det.get("B_pairwise_structure", 0)
        + det.get("C_underdetermined", 0)
        + det.get("incompatible", 0)
        + det.get("other", 0)
    )
    gsum = (ga or {}).get("summary", {})
    return {
        "with_vector": with_vector,
        "total_regions": sum(tt.values()),
        "A_determined": prefix_get(det, "A_determined"),
        "A_ambiguous": prefix_get(det, "A_ambiguous"),
        "B_pairwise": det.get("B_pairwise_structure", 0),
        "C_underdetermined": det.get("C_underdetermined", 0),
        "branched": prefix_get(tt, "branched"),
        "linear": prefix_get(tt, "linear"),
        "single": prefix_get(tt, "single"),
        "has_cycle": det.get("incompatible", 0),
        "cancer_gene_regions": gsum.get("regions_with_cancer_gene", 0),
        "druggable_regions": gsum.get("regions_with_druggable", 0),
        "promoter_regions": gsum.get("regions_with_promoter", 0),
        "queue_n": (cs or {}).get("n_need_confirm", 0),
    }


def main():
    out = {}
    for s in SAMPLES:
        r = summarize(s)
        if r is None:
            print(f"  WARN {s}: topology_per_region.json missing, skipped")
            continue
        out[s] = r
    out["_provenance"] = {
        "regenerated": "2026-07-02 post C2/C3/D4 fix",
        "note": "aggregated from per-sample topology_per_region.json (07-01/07-02 post-fix) "
                "+ candidate_scoring.json (queue_n) + region_gene_annotation.json (fix-independent). "
                "HCC1395 from 4axis canonical; 6 others from MSROOT.",
        "determinacy_denominator": "with_vector",
    }
    path = os.path.join(DATA4AXIS, "multisample_summary.json")
    json.dump(out, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    # human-readable table
    cols = ["total_regions", "with_vector", "A_determined", "branched", "linear", "has_cycle",
            "cancer_gene_regions", "druggable_regions", "queue_n"]
    print(f"{'sample':16}" + "".join(f"{c[:9]:>10}" for c in cols))
    for s in SAMPLES:
        if s not in out:
            continue
        r = out[s]
        print(f"{s:16}" + "".join(f"{r[c]:>10}" for c in cols))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
