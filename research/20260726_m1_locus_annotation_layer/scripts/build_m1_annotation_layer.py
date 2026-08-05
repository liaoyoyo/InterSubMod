#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the M1 locus-level annotation layer from authority artifacts only.

Supersedes the v1 derived table: the coarse group count is taken from the authority
`coarse_ng`, never re-derived from read labels (the v1 derivation disagreed with authority
on 50,531 / 102,842 sites because it counted peeled reads as their own group).

Three orthogonal axes, each with its own claim ceiling:

  A  methyl multigroup   coarse_ng / cluster_sizes / evidence_tier / confound flags
                         -> "this locus's focal-ALT reads split into K stable methyl groups"
  B  ALT/REF allele axis joint_allele_axis_aligned / V / permutation p / REF-ALT depths
                         -> "the methyl split runs along the ALT vs REF axis"
                            NOT "the mutation caused it" - germline-het null not run
  C  HP anchoring        alt_hp_family_counts -> HP1_SIDE / HP2_SIDE / CROSS_HP / AMBIGUOUS
                         -> somatic-ALT phasing quality. NOT an LOH call: LOH needs the
                            full (REF+ALT) read-set haplotype composition, which this
                            artifact does not carry.

Also joins every site to the exact-PS final groups via MLHP `positions`, so the workstation
can colour regions and plot loci that belong to no region (k=1 abstain sites).

Fail-closed: unknown axis classes, unparsable coarse_ng, or a row-count drift abort.
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
MLHP_ROOT = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples"
)
EXPECTED_ROWS = 102842

TIER_CODES = ["E2_STABLE_NULL_MULTIGROUP",
              "E3_UNEXPLAINED_AFTER_MEASURED_AXES",
              "E4_PHASE_ANCHORED_ROBUST_EPIGENETIC_CANDIDATE"]
HP_CODES = ["HP1_SIDE", "HP2_SIDE", "CROSS_HP", "AMBIGUOUS_ONLY", "NO_HP"]
ALLELE_CODES = ["ALIGNED", "TESTED_NOT_ALIGNED", "NOT_TESTABLE"]

TSV_COLUMNS = [
    "dataset", "chrom", "pos", "ref", "alt", "truth_label",
    # axis A
    "coarse_ng", "cluster_sizes", "core_read_n", "min_group_n",
    "evidence_tier", "residual_unexplained_multigroup",
    "hp_axis_confound", "technical_axis_confound",
    # axis B
    "allele_class", "joint_allele_v", "joint_allele_p_perm", "joint_allele_n",
    "ref_status", "n_tumor_alt", "n_tumor_ref", "caller_af", "normal_af", "normal_dp",
    # axis C
    "hp_anchor_class", "hp_ambiguous_fraction", "alt_hp_family_counts",
    # region join
    "in_final_group", "region_ids",
]


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


def classify_hp(raw):
    """Axis C from the authority `alt_hp_family_counts` JSON object."""
    if not raw or not raw.strip():
        return "NO_HP", 0.0, {}
    try:
        counts = json.loads(raw)
    except json.JSONDecodeError:
        return "NO_HP", 0.0, {}
    if not counts:
        return "NO_HP", 0.0, {}
    total = sum(int(v) for v in counts.values())
    hp1 = int(counts.get("HP1-side", 0))
    hp2 = int(counts.get("HP2-side", 0))
    amb = total - hp1 - hp2
    frac = (amb / total) if total else 0.0
    if hp1 and hp2:
        cls = "CROSS_HP"
    elif hp1:
        cls = "HP1_SIDE"
    elif hp2:
        cls = "HP2_SIDE"
    elif total:
        cls = "AMBIGUOUS_ONLY"
    else:
        cls = "NO_HP"
    return cls, frac, counts


def classify_allele(row):
    """Axis B."""
    if truthy(row.get("joint_allele_axis_aligned")):
        return "ALIGNED"
    if truthy(row.get("joint_allele_testable")):
        return "TESTED_NOT_ALIGNED"
    return "NOT_TESTABLE"


def cluster_stats(raw):
    """core read count and smallest group size from the authority `cluster_sizes` map."""
    if not raw or not raw.strip():
        return None, None
    try:
        sizes = json.loads(raw)
    except json.JSONDecodeError:
        return None, None
    vals = [int(v) for v in sizes.values()]
    if not vals:
        return None, None
    return sum(vals), min(vals)


def load_region_index(datasets):
    """(dataset, chrom, pos) -> [region_id, ...] from every sample's MLHP groups."""
    idx = collections.defaultdict(list)
    per_ds_groups = {}
    for ds in datasets:
        path = os.path.join(MLHP_ROOT, ds, f"{ds}.exact_ps_mlhp.json")
        if not os.path.exists(path):
            sys.exit(f"ERROR: MLHP missing for {ds}: {path}")
        with open(path) as fh:
            doc = json.load(fh)
        groups = doc.get("groups") or []
        per_ds_groups[ds] = len(groups)
        for g in groups:
            rid = g.get("region_id")
            chrom = g.get("chrom")
            for pos in (g.get("positions") or []):
                idx[(ds, chrom, int(pos))].append(rid)
    return idx, per_ds_groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--control", default=CONTROL)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    data_dir = os.path.join(a.out_dir, "data")
    side_dir = os.path.join(data_dir, "sidecar")
    res_dir = os.path.join(a.out_dir, "results")
    for d in (data_dir, side_dir, res_dir):
        os.makedirs(d, exist_ok=True)

    # pass 1: which datasets appear
    datasets = []
    seen = set()
    with gzip.open(a.control, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ds = row["dataset"]
            if ds not in seen:
                seen.add(ds)
                datasets.append(ds)
    region_index, per_ds_groups = load_region_index(datasets)

    tsv_path = os.path.join(data_dir, "m1_annotation_site_table_v2.tsv")
    per_dataset = collections.Counter()
    tier_c = collections.Counter()
    hp_c = collections.Counter()
    allele_c = collections.Counter()
    ng_c = collections.Counter()
    in_group_c = collections.Counter()
    cross = collections.Counter()          # (allele_class, hp_anchor_class)
    ref_absent = collections.Counter()
    per_ds_payload = collections.defaultdict(lambda: collections.defaultdict(list))
    per_ds_regions = collections.defaultdict(list)
    per_ds_region_pos = collections.defaultdict(dict)
    per_ds_chroms = collections.defaultdict(list)
    per_ds_chrom_pos = collections.defaultdict(dict)
    n = 0

    with gzip.open(a.control, "rt") as fh, open(tsv_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=TSV_COLUMNS, delimiter="\t")
        w.writeheader()
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1
            ds, chrom = row["dataset"], row["chrom"]
            pos = as_int(row["pos"])
            if pos is None:
                sys.exit(f"ERROR: row {n} unparsable pos -> fail closed")

            ng = as_int(row.get("coarse_ng"))
            if ng is None or ng < 2:
                sys.exit(f"ERROR: row {n} coarse_ng={row.get('coarse_ng')!r} violates M1 "
                         f"contract (>=2) -> fail closed")

            tier = row.get("evidence_tier", "")
            if tier not in TIER_CODES:
                sys.exit(f"ERROR: row {n} unknown evidence_tier {tier!r} -> fail closed")

            hp_cls, amb_frac, _ = classify_hp(row.get("alt_hp_family_counts", ""))
            al_cls = classify_allele(row)
            core_n, min_n = cluster_stats(row.get("cluster_sizes", ""))
            rids = region_index.get((ds, chrom, pos), [])

            w.writerow({
                "dataset": ds, "chrom": chrom, "pos": pos,
                "ref": row.get("ref", ""), "alt": row.get("alt", ""),
                "truth_label": row.get("truth_label", ""),
                "coarse_ng": ng, "cluster_sizes": row.get("cluster_sizes", ""),
                "core_read_n": core_n if core_n is not None else "",
                "min_group_n": min_n if min_n is not None else "",
                "evidence_tier": tier,
                "residual_unexplained_multigroup": row.get("residual_unexplained_multigroup", ""),
                "hp_axis_confound": row.get("hp_axis_confound", ""),
                "technical_axis_confound": row.get("technical_axis_confound", ""),
                "allele_class": al_cls,
                "joint_allele_v": row.get("joint_allele_v", ""),
                "joint_allele_p_perm": row.get("joint_allele_p_perm", ""),
                "joint_allele_n": row.get("joint_allele_n", ""),
                "ref_status": row.get("ref_status", ""),
                "n_tumor_alt": row.get("n_tumor_alt", ""),
                "n_tumor_ref": row.get("n_tumor_ref", ""),
                "caller_af": row.get("caller_af", ""),
                "normal_af": row.get("normal_af", ""),
                "normal_dp": row.get("normal_dp", ""),
                "hp_anchor_class": hp_cls,
                "hp_ambiguous_fraction": f"{amb_frac:.6f}",
                "alt_hp_family_counts": row.get("alt_hp_family_counts", ""),
                "in_final_group": bool(rids),
                "region_ids": ";".join(rids),
            })

            per_dataset[ds] += 1
            tier_c[tier] += 1
            hp_c[hp_cls] += 1
            allele_c[al_cls] += 1
            ng_c[min(ng, 9)] += 1
            in_group_c["in_group" if rids else "locus_only"] += 1
            cross[(al_cls, hp_cls)] += 1
            nref = as_int(row.get("n_tumor_ref"))
            ref_absent["ref_0" if nref == 0 else "ref_gt0"] += 1

            # compact per-dataset sidecar payload
            p = per_ds_payload[ds]
            cmap = per_ds_chrom_pos[ds]
            if chrom not in cmap:
                cmap[chrom] = len(per_ds_chroms[ds])
                per_ds_chroms[ds].append(chrom)
            # a physical locus can belong to several exact-PS x HP containers, so every
            # region it lands in must be recorded -- keeping only the first undercounts
            # the regions that carry an M1 site (28,983 vs the true 45,008).
            rmap = per_ds_region_pos[ds]
            r_idx = []
            for rid in rids:
                if rid not in rmap:
                    rmap[rid] = len(per_ds_regions[ds])
                    per_ds_regions[ds].append(rid)
                r_idx.append(rmap[rid])
            p["c"].append(cmap[chrom])
            p["p"].append(pos)
            p["ng"].append(ng)
            p["t"].append(TIER_CODES.index(tier))
            p["h"].append(HP_CODES.index(hp_cls))
            p["a"].append(ALLELE_CODES.index(al_cls))
            v = as_float(row.get("joint_allele_v"))
            p["v"].append(round(v, 4) if v is not None else None)
            p["na"].append(as_int(row.get("n_tumor_alt")))
            p["nr"].append(as_int(row.get("n_tumor_ref")))
            p["r"].append(r_idx)

    if n != EXPECTED_ROWS:
        sys.exit(f"ERROR: row count {n} != expected {EXPECTED_ROWS} -> fail closed")

    sidecar_files = {}
    for ds in datasets:
        p = per_ds_payload[ds]
        doc = {
            "schema_name": "intersubmod.m1_locus_annotation.sidecar",
            "schema_version": "2.0.0",
            "dataset": ds,
            "n_sites": len(p["p"]),
            "codes": {"tier": TIER_CODES, "hp": HP_CODES, "allele": ALLELE_CODES},
            "chroms": per_ds_chroms[ds],
            "regions": per_ds_regions[ds],
            "fields": {"c": "chrom index", "p": "position", "ng": "authority coarse_ng",
                       "t": "evidence_tier index", "h": "hp_anchor_class index",
                       "a": "allele_class index", "v": "joint_allele Cramer V",
                       "na": "n_tumor_alt", "nr": "n_tumor_ref",
                       "r": ("list of region indices this locus belongs to; empty list = "
                             "not in any exact-PS final group. A locus can sit in several "
                             "exact PS x HP containers, so this is a list, not a scalar.")},
            "sites": {k: p[k] for k in ("c", "p", "ng", "t", "h", "a", "v", "na", "nr", "r")},
            "claim_ceiling": (
                "Axis A is an operational methyl-multigroup screen, not a cellular subclone. "
                "Axis B alignment with the ALT/REF axis does not establish that the mutation "
                "caused the methylation difference; the germline-het null has not been run. "
                "Axis C is somatic-ALT phasing quality, not an LOH call."),
        }
        fp = os.path.join(side_dir, f"{ds}.m1_annotation.json")
        with open(fp, "w") as fh:
            json.dump(doc, fh, ensure_ascii=False, separators=(",", ":"))
        h = hashlib.sha256()
        with open(fp, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        sidecar_files[ds] = {"path": fp, "bytes": os.path.getsize(fp),
                             "sha256": h.hexdigest(), "n_sites": doc["n_sites"],
                             "n_regions_referenced": len(per_ds_regions[ds])}

    h = hashlib.sha256()
    with open(tsv_path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)

    receipt = {
        "schema_name": "intersubmod.m1_locus_annotation.build_receipt",
        "schema_version": "2.0.0",
        "supersedes": "v1 derived-K table (50,531/102,842 disagreed with authority coarse_ng)",
        "authority_control": a.control,
        "authority_control_bytes": os.path.getsize(a.control),
        "mlhp_root": MLHP_ROOT,
        "mlhp_groups_per_dataset": per_ds_groups,
        "rows": n,
        "per_dataset": dict(per_dataset),
        "axis_A_evidence_tier": dict(tier_c),
        "axis_A_coarse_ng_capped9": {str(k): v for k, v in sorted(ng_c.items())},
        "axis_B_allele_class": dict(allele_c),
        "axis_B_tumor_ref_absent": dict(ref_absent),
        "axis_C_hp_anchor_class": dict(hp_c),
        "axis_B_x_axis_C": {f"{a_}|{h_}": v for (a_, h_), v in sorted(cross.items())},
        "region_join": dict(in_group_c),
        "out_tsv": tsv_path,
        "out_tsv_bytes": os.path.getsize(tsv_path),
        "out_tsv_sha256": h.hexdigest(),
        "sidecars": sidecar_files,
        "checks": {
            "row_count_equals_expected": n == EXPECTED_ROWS,
            "all_coarse_ng_ge_2": True,
            "all_evidence_tier_known": True,
            "hp_classes_closed_set": set(hp_c) <= set(HP_CODES),
            "allele_classes_closed_set": set(allele_c) <= set(ALLELE_CODES),
            "region_join_partition": (in_group_c["in_group"] + in_group_c["locus_only"]) == n,
        },
        "claim_ceiling": {
            "axis_A": "operational methyl-multigroup screen; not a cellular subclone or clone count.",
            "axis_B": ("ALT/REF axis alignment. Confounded with germline allele-specific "
                       "methylation whenever the somatic ALT is confined to one HP side; "
                       "the germline-het null has NOT been run."),
            "axis_C": ("somatic-ALT phasing quality (HP1_SIDE / HP2_SIDE / CROSS_HP / "
                       "AMBIGUOUS_ONLY). This is NOT an LOH determination - LOH requires the "
                       "full REF+ALT read-set haplotype composition."),
            "region_join": ("region_ids come from exact-PS MLHP `positions`. locus_only sites "
                            "are not k=1 topology results; they are simply outside every "
                            "final group."),
        },
    }
    receipt["checks"]["hp_classes_closed_set"] = bool(receipt["checks"]["hp_classes_closed_set"])
    receipt["checks"]["allele_classes_closed_set"] = bool(receipt["checks"]["allele_classes_closed_set"])
    receipt["all_checks_pass"] = all(receipt["checks"].values())

    rp = os.path.join(res_dir, "m1_annotation_receipt_v2.json")
    with open(rp, "w") as fh:
        json.dump(receipt, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in receipt.items() if k != "sidecars"},
                     ensure_ascii=False, indent=1))
    print("\nsidecars:")
    for ds, meta in sidecar_files.items():
        print(f"  {ds:<16} {meta['n_sites']:>7,} sites  {meta['bytes']/1e6:>6.2f} MB  "
              f"{meta['n_regions_referenced']:>6,} regions  {meta['sha256'][:12]}")


if __name__ == "__main__":
    main()
