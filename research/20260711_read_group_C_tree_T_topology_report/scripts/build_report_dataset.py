#!/usr/bin/env python3
"""Build the corrected C/T/topology/VAF report dataset for seven datasets."""

from __future__ import annotations

import argparse
import collections
import csv
import json
import math
from datetime import datetime
from pathlib import Path


def count_bins(values, cap=6):
    counts = collections.Counter()
    for value in values:
        counts[str(value) if value <= cap else f">{cap}"] += 1
    return {key: counts.get(key, 0) for key in [*(str(i) for i in range(cap + 1)), f">{cap}"]}


def t_bins(values):
    counts = collections.Counter()
    for value in values:
        key = str(value) if value <= 6 else ">6"
        counts[key] += 1
    return {key: counts.get(key, 0) for key in ["1", "2", "3", "4", "5", "6", ">6"]}


def region_key(group):
    return f"{group['chrom']}:{group['start']}-{group['end']}"


def load_groups(sample_dir):
    groups = {}
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        document = json.loads(path.read_text(encoding="utf-8"))
        for group in document.get("groups", []):
            key = region_key(group)
            if key in groups:
                raise RuntimeError(f"duplicate group {sample_dir.name} {key}")
            groups[key] = group
    return groups


def c_alt_for_family(group, family):
    full = ((group.get("populations_by_hp") or {}).get(family) or {})
    return sum("A" in genotype for genotype in full)


def read_af_available(group, family):
    positions = group.get("positions", [])
    coverage = ((group.get("col_coverage_by_hp") or {}).get(family) or {})
    for position in positions:
        counts = coverage.get(str(position)) or coverage.get(position)
        if not counts or sum(counts[:2]) == 0:
            return False
    return True


def cn_class(value):
    return value if value in {"neutral", "gain", "loss", "loh"} else "unavailable"


def load_layered(sample_dir, sample):
    path = sample_dir / f"layered_reconstruction_{sample}.json"
    return path, json.loads(path.read_text(encoding="utf-8"))


def analyze_sample(meta, run_root, topo_sample, regional_sample, read_af_sample):
    sample = meta["sample"]
    sample_dir = run_root / "samples" / sample
    layered_path, layered = load_layered(sample_dir, sample)
    groups = load_groups(sample_dir)
    detail = layered["detail"]
    tree_regions = {unit["region"] for unit in detail}
    primary_units = [unit for unit in detail if unit.get("is_primary_lineage")]
    primary_by_region = collections.defaultdict(list)
    for unit in primary_units:
        primary_by_region[unit["region"]].append(unit)
    h3_regions = {
        unit["region"]
        for unit in detail
        if unit.get("is_h3_auxiliary") and unit.get("mutation_bearing")
    }

    pooled_values = []
    pooled_contract_ok = True
    support_contract_ok = True
    for group in groups.values():
        populations = group.get("populations") or {}
        expected_alt = sum("A" in genotype for genotype in populations)
        pooled_contract_ok &= (
            group.get("n_populations") == len(populations)
            and group.get("n_populations_with_ALT") == expected_alt
        )
        support_contract_ok &= all(count >= 3 for count in populations.values())
        pooled_values.append(expected_alt)

    c_unit = {}
    c_unit_values = []
    for unit in primary_units:
        group = groups.get(unit["region"])
        if group is None:
            raise RuntimeError(f"missing group for primary unit {sample} {unit['region']}")
        value = c_alt_for_family(group, unit["family"])
        c_unit[(unit["region"], unit["family"])] = value
        c_unit_values.append(value)

    hp_counts = collections.Counter()
    hp_h3 = collections.Counter()
    c_region_sum_values = []
    single_c = []
    double_c_sum = []
    double_pairs = collections.Counter()
    t_joint_values = []
    topology_alternative_values = []
    topology_classes = collections.Counter()
    nested_region_cross = collections.Counter()
    topology_by_hp = collections.defaultdict(collections.Counter)
    hidden_cross = collections.defaultdict(collections.Counter)
    hidden_zero = 0
    hidden_positive = 0
    incomplete_regions = 0

    for region in sorted(tree_regions):
        units = sorted(primary_by_region.get(region, []), key=lambda unit: unit["family"])
        if len(units) == 1:
            hp_class = "single_HP1_xor_HP2"
        elif len(units) == 2:
            hp_class = "double_HP1_and_HP2"
        elif len(units) == 0:
            hp_class = "no_primary"
        else:
            raise RuntimeError(f"unexpected primary multiplicity {sample} {region}: {len(units)}")
        hp_counts[hp_class] += 1
        hp_h3[f"{hp_class}|{'with_H3' if region in h3_regions else 'without_H3'}"] += 1
        if not units:
            continue

        c_values = [c_unit[(region, unit["family"])] for unit in units]
        c_sum = sum(c_values)
        c_region_sum_values.append(c_sum)
        if len(units) == 1:
            single_c.append(c_sum)
            c_state = f"HP{units[0]['family']}={c_values[0]}"
        else:
            by_family = {unit["family"]: c_unit[(region, unit["family"])] for unit in units}
            c1, c2 = by_family.get("1", 0), by_family.get("2", 0)
            double_pairs[f"HP1={c1};HP2={c2}"] += 1
            double_c_sum.append(c1 + c2)
            c_state = f"HP1={c1};HP2={c2}"

        if not all(unit.get("analysis_candidate_set_complete") is True for unit in units):
            incomplete_regions += 1
            topology_classes["incomplete"] += 1
            topology_by_hp[hp_class]["incomplete"] += 1
            nested_region_cross[(hp_class, c_state, "incomplete", "not_evaluated")] += 1
            continue
        t_joint = math.prod(int(unit["n_trees"]) for unit in units)
        topo_joint = math.prod(int(unit["n_distinct_shapes_exact"]) for unit in units)
        t_joint_values.append(t_joint)
        topology_alternative_values.append(topo_joint)
        hidden = sum(int(unit["n_hidden"]) for unit in units)
        if t_joint == 1 and topo_joint == 1:
            category = "T=1|Topo=1"
        elif t_joint > 1 and topo_joint == 1:
            category = "T>1|Topo=1"
        elif t_joint > 1 and topo_joint > 1:
            category = "T>1|Topo>1"
        else:
            category = "T=1|Topo>1_impossible"
        topology_classes[category] += 1
        topology_by_hp[hp_class][category] += 1
        hidden_cross[category]["hidden=0" if hidden == 0 else "hidden>0"] += 1
        nested_region_cross[
            (hp_class, c_state, category, "hidden=0" if hidden == 0 else "hidden>0")
        ] += 1
        if hidden == 0:
            hidden_zero += 1
        else:
            hidden_positive += 1

    ambiguous_units = [
        unit
        for unit in primary_units
        if not unit.get("capped") and unit.get("L1_base_class") == "ambiguous"
    ]
    prepared_by_cn = collections.Counter()
    for unit in ambiguous_units:
        group = groups[unit["region"]]
        if read_af_available(group, unit["family"]):
            prepared_by_cn[cn_class(unit.get("cn"))] += 1
    default_rank = read_af_sample.get("default") or {}
    reached_by_cn = collections.Counter(default_rank.get("by_cn") or {})

    funnel = layered.get("input_funnel") or {}
    w_k1 = int(funnel.get("n_positional_singleton", 0))
    w_k_gt1 = int(funnel.get("n_multilocus_pre_cap_groups", 0))
    w_ret = len(groups)
    w_tree = len(tree_regions)
    w_primary = len(primary_by_region)
    topology_complete = sum(
        topology_classes[key]
        for key in ("T=1|Topo=1", "T>1|Topo=1", "T>1|Topo>1")
    )

    checks = {
        "scope_conservation": bool(funnel.get("check_scope_conservation")),
        "W_kgt1_conservation": w_k_gt1
        == int(funnel.get("n_groups_retained", -1))
        + int(funnel.get("n_groups_read_unsupported", 0)),
        "W_ret_matches_groups": w_ret == int(funnel.get("n_groups_retained", -1)),
        "W_tree_HP_conservation": w_tree == sum(hp_counts.values()),
        "primary_region_conservation": w_primary
        == hp_counts["single_HP1_xor_HP2"] + hp_counts["double_HP1_and_HP2"],
        "primary_unit_conservation": len(primary_units)
        == hp_counts["single_HP1_xor_HP2"] + 2 * hp_counts["double_HP1_and_HP2"],
        "pooled_C_conservation": len(pooled_values) == w_ret,
        "primary_unit_C_conservation": len(c_unit_values) == len(primary_units),
        "primary_region_C_conservation": len(c_region_sum_values) == w_primary,
        "topology_region_conservation": topology_complete + incomplete_regions == w_primary,
        "impossible_topology_zero": topology_classes["T=1|Topo>1_impossible"] == 0,
        "hidden_conservation": hidden_zero + hidden_positive == topology_complete,
        "nested_region_cross_conservation": sum(nested_region_cross.values()) == w_primary,
        "pooled_field_contract": pooled_contract_ok,
        "pooled_MINREAD_contract": support_contract_ok,
        "topology_catalog_checks": all(topo_sample["checks"].values()),
        "topology_primary_units_match": topo_sample["primary_units"] == len(primary_units),
        "regional_catalog_checks": all(regional_sample["checks"].values()),
        "regional_primary_regions_match": regional_sample["primary_regions"] == w_primary,
        "regional_complete_regions_match": (
            regional_sample["fully_complete_regions"] == topology_complete
        ),
        "regional_incomplete_regions_match": (
            regional_sample["incomplete_regions"] == incomplete_regions
        ),
        "regional_topology_alternatives_match": (
            regional_sample["sum_region_topology_alternatives"]
            == sum(topology_alternative_values)
        ),
        "regional_exact_joint_T_match": (
            regional_sample["sum_exact_joint_tree_candidates"] == sum(t_joint_values)
        ),
        "read_AF_prepared_matches": sum(prepared_by_cn.values())
        == read_af_sample["n_units_analyzed_all_candidates"],
        "read_AF_missing_matches": len(ambiguous_units) - sum(prepared_by_cn.values())
        == read_af_sample["n_missing_read_af"],
        "read_AF_candidate_mismatch_zero": read_af_sample["n_candidate_mismatch_or_incomplete"] == 0,
    }

    return {
        "sample": sample,
        "biological_id": meta["biological_id"],
        "source_layered": str(layered_path),
        "S": {
            "all_contigs": meta["vcf"]["snv_total"],
            "autosomal_chr1_22": meta["vcf"]["autosomal_chr1_22"],
            "out_of_scope": meta["vcf"]["out_of_scope"],
        },
        "W": {
            "all_pre": w_k1 + w_k_gt1,
            "k1": w_k1,
            "k_gt1": w_k_gt1,
            "retained": w_ret,
            "tree_view": w_tree,
            "primary": w_primary,
            "primary_units": len(primary_units),
        },
        "HP": {
            "classes": dict(hp_counts),
            "by_H3": dict(hp_h3),
            "H3_regions": len(h3_regions),
        },
        "C": {
            "definition": "ALT-containing MINREAD-supported full R/A groups; ROOT/all-reference/hidden/partial excluded",
            "pooled_region": count_bins(pooled_values),
            "primary_HP_unit": count_bins(c_unit_values),
            "primary_region_sum": count_bins(c_region_sum_values),
            "single_HP_region": count_bins(single_c),
            "double_HP_region_sum": count_bins(double_c_sum),
            "double_HP_pairs": dict(sorted(double_pairs.items())),
            "pooled_vs_primary_HP_sum_discordant_primary_regions": sum(
                int(groups[region].get("n_populations_with_ALT", 0))
                != sum(c_unit[(region, unit["family"])] for unit in units)
                for region, units in primary_by_region.items()
            ),
            "pooled_vs_HP1_HP2_sum_discordant_all_retained_regions": sum(
                int(group.get("n_populations_with_ALT", 0))
                != c_alt_for_family(group, "1") + c_alt_for_family(group, "2")
                for group in groups.values()
            ),
            "pooled_vs_all_family_sum_discordant_all_retained_regions": sum(
                int(group.get("n_populations_with_ALT", 0))
                != sum(
                    sum("A" in genotype for genotype in populations)
                    for populations in (group.get("populations_by_hp") or {}).values()
                )
                for group in groups.values()
            ),
        },
        "T_and_Topology": {
            "T_joint_bins": {
                **t_bins(t_joint_values),
                "incomplete": incomplete_regions,
            },
            "classes": dict(topology_classes),
            "classes_by_HP": {
                key: dict(value) for key, value in topology_by_hp.items()
            },
            "nested_region_cross": [
                {
                    "HP_scope": hp_class,
                    "C_state": c_state,
                    "T_Topo_state": category,
                    "H_state": hidden_state,
                    "regions": count,
                }
                for (hp_class, c_state, category, hidden_state), count in sorted(
                    nested_region_cross.items()
                )
            ],
            "hidden": {
                "hidden=0": hidden_zero,
                "hidden>0": hidden_positive,
                "incomplete": incomplete_regions,
                "by_topology": {
                    key: dict(value) for key, value in hidden_cross.items()
                },
            },
            "shape_catalog": {
                "complete_primary_units": topo_sample["complete_units"],
                "incomplete_primary_units": topo_sample["incomplete_units"],
                "exact_candidate_trees": topo_sample["candidate_trees"],
                "unit_shape_incidence": topo_sample["unit_shape_incidence"],
                "distinct_shapes": topo_sample["distinct_shapes"],
                "rerun_units": topo_sample["rerun_units"],
            },
            "regional_ordered_forest": regional_sample,
        },
        "VAF_ranking": {
            "status": "exploratory family-specific read-AF ordering; not independent confirmation",
            "ambiguous_primary_units": read_af_sample["n_ambiguous_primary_units"],
            "prepared_all_candidates": read_af_sample["n_units_analyzed_all_candidates"],
            "missing_read_AF": read_af_sample["n_missing_read_af"],
            "candidate_mismatch": read_af_sample["n_candidate_mismatch_or_incomplete"],
            "prepared_by_CN": dict(prepared_by_cn),
            "top_weight_ge_0_60": default_rank.get("reached", 0),
            "top_weight_lt_0_60": read_af_sample["n_units_analyzed_all_candidates"]
            - default_rank.get("reached", 0),
            "winner_order_consistent": default_rank.get("winner_order_consistent", 0),
            "reached_by_CN": dict(reached_by_cn),
            "neutral_reach_rate": (
                reached_by_cn["neutral"] / prepared_by_cn["neutral"]
                if prepared_by_cn["neutral"]
                else None
            ),
        },
        "checks": checks,
    }


def build_global_shapes(topology, sample_order):
    global_rows = {}
    shape_only_by_sample = collections.defaultdict(collections.Counter)
    for sample in topology["samples"]:
        sample_name = sample["sample"]
        for unit in sample["unit_rows"]:
            shape_ids = list(unit["shape_candidate_counts"])
            if len(shape_ids) == 1:
                shape_only_by_sample[sample_name][shape_ids[0]] += 1
        for row in sample["catalog"]:
            stable_id = row["shape_id"]
            target = global_rows.setdefault(
                stable_id,
                {
                    key: row[key]
                    for key in (
                        "shape_id",
                        "signature",
                        "n_nodes_including_root",
                        "n_mutation_state_nodes",
                        "n_edges",
                        "n_leaves",
                        "max_depth",
                        "root_degree",
                        "n_internal_branch_nodes",
                        "max_outdegree",
                        "coarse_shape",
                        "example_region",
                        "example_hp",
                        "example_edges",
                    )
                },
            )
            target.setdefault("candidate_trees", 0)
            target.setdefault("unit_incidence", 0)
            target.setdefault("region_incidence", 0)
            target.setdefault("shape_only_units", 0)
            target.setdefault("by_sample_units", {})
            target.setdefault("by_sample_candidates", {})
            target["candidate_trees"] += row["candidate_trees"]
            target["unit_incidence"] += row["unit_incidence"]
            target["region_incidence"] += row["region_incidence"]
            target["shape_only_units"] += shape_only_by_sample[sample_name][stable_id]
            target["by_sample_units"][sample_name] = row["unit_incidence"]
            target["by_sample_candidates"][sample_name] = row["candidate_trees"]
    rows = sorted(
        global_rows.values(),
        key=lambda row: (-row["unit_incidence"], row["n_nodes_including_root"], row["signature"]),
    )
    for rank, row in enumerate(rows, 1):
        row["display_id"] = f"Topo_{rank:02d}"
        for sample in sample_order:
            row["by_sample_units"].setdefault(sample, 0)
            row["by_sample_candidates"].setdefault(sample, 0)
    return rows


def aggregate_samples(samples, shapes):
    def sum_path(*path):
        return sum(
            next_value(sample, path)
            for sample in samples
        )

    aggregate = {
        "dataset_rows": len(samples),
        "biological_samples": len({sample["biological_id"] for sample in samples}),
        "S_all_contigs": sum_path("S", "all_contigs"),
        "S_autosomal": sum_path("S", "autosomal_chr1_22"),
        "W_all_pre": sum_path("W", "all_pre"),
        "W_k1": sum_path("W", "k1"),
        "W_k_gt1": sum_path("W", "k_gt1"),
        "W_retained": sum_path("W", "retained"),
        "W_tree": sum_path("W", "tree_view"),
        "W_primary": sum_path("W", "primary"),
        "primary_units": sum_path("W", "primary_units"),
        "exact_candidate_trees": sum_path("T_and_Topology", "shape_catalog", "exact_candidate_trees"),
        "unit_shape_incidence": sum_path("T_and_Topology", "shape_catalog", "unit_shape_incidence"),
        "distinct_shapes_global": len(shapes),
        "region_topology_alternatives": sum_path(
            "T_and_Topology", "regional_ordered_forest", "sum_region_topology_alternatives"
        ),
        "joint_exact_T_candidates": sum_path(
            "T_and_Topology", "regional_ordered_forest", "sum_exact_joint_tree_candidates"
        ),
        "VAF_ambiguous_units": sum_path("VAF_ranking", "ambiguous_primary_units"),
        "VAF_prepared": sum_path("VAF_ranking", "prepared_all_candidates"),
        "VAF_missing": sum_path("VAF_ranking", "missing_read_AF"),
        "VAF_top_weight_ge_0_60": sum_path("VAF_ranking", "top_weight_ge_0_60"),
        "VAF_winner_order_consistent": sum_path("VAF_ranking", "winner_order_consistent"),
    }
    for section, field in (
        ("C", "pooled_region"),
        ("C", "primary_HP_unit"),
        ("C", "primary_region_sum"),
        ("T_and_Topology", "T_joint_bins"),
        ("T_and_Topology", "classes"),
    ):
        total = collections.Counter()
        for sample in samples:
            total.update(sample[section][field])
        aggregate[f"{section}.{field}"] = dict(total)
    return aggregate


def next_value(document, path):
    value = document
    for key in path:
        value = value[key]
    return value


def write_tsv(path, rows, fieldnames):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--regional-topology", required=True, type=Path)
    parser.add_argument("--read-af", required=True, type=Path)
    parser.add_argument("--vaf-example", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    topology = json.loads(args.topology.read_text(encoding="utf-8"))
    regional_topology = json.loads(args.regional_topology.read_text(encoding="utf-8"))
    read_af = json.loads(args.read_af.read_text(encoding="utf-8"))
    vaf_example = json.loads(args.vaf_example.read_text(encoding="utf-8"))
    topo_by_sample = {sample["sample"]: sample for sample in topology["samples"]}
    regional_by_sample = {
        sample["sample"]: sample for sample in regional_topology["samples"]
    }
    read_af_by_sample = {sample["sample"]: sample for sample in read_af["samples"]}
    sample_order = [sample["sample"] for sample in manifest["samples"]]

    samples = [
        analyze_sample(
            meta,
            args.run_root,
            topo_by_sample[meta["sample"]],
            regional_by_sample[meta["sample"]],
            read_af_by_sample[meta["sample"]],
        )
        for meta in manifest["samples"]
    ]
    shapes = build_global_shapes(topology, sample_order)
    aggregate = aggregate_samples(samples, shapes)

    hcc = next(sample for sample in samples if sample["sample"] == "HCC1395")
    hcc_dir = args.run_root / "samples" / "HCC1395"
    hcc_groups = load_groups(hcc_dir)
    _, hcc_layered = load_layered(hcc_dir, "HCC1395")
    example_region = "chr1:3294434-3310766"
    example_group = hcc_groups[example_region]
    example_unit = next(
        unit
        for unit in hcc_layered["detail"]
        if unit["region"] == example_region and unit.get("is_primary_lineage")
    )
    example_full = (
        (example_group.get("populations_by_hp") or {}).get(example_unit["family"]) or {}
    )
    example = {
        "region": example_region,
        "k": example_group["n_sSNV"],
        "HP": example_unit["family"],
        "read_groups": example_full,
        "C": sum("A" in genotype for genotype in example_full),
        "T": example_unit["n_trees"],
        "Topo": example_unit["n_distinct_shapes_exact"],
        "hidden": example_unit["n_hidden"],
        "edges": example_unit["trees"][0]["edges"],
        "site_read_AF": {
            str(position): (
                ((example_group["col_coverage_by_hp"][example_unit["family"]][str(position)][1])
                 / sum(example_group["col_coverage_by_hp"][example_unit["family"]][str(position)][:2]))
            )
            for position in example_group["positions"]
        },
    }

    validation_failures = [
        f"{sample['sample']}:{check}"
        for sample in samples
        for check, passed in sample["checks"].items()
        if not passed
    ]
    output = {
        "schema_version": "1.0",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "partial_historical_engineering_snapshot",
        "claim_ceiling": "regional mutation-state candidate trees; not confirmed clone/subclone phylogeny",
        "definitions": {
            "C_unit": "ALT-containing, MINREAD-supported full R/A read groups for one primary HP; excludes ROOT/reference/hidden/partial",
            "C_region": "sum of C_unit across primary HP1/HP2; double-HP pair retained",
            "T_unit": "number of exact candidate Steiner arborescences for one primary HP unit",
            "T_joint": "product of T_unit across primary HP1/HP2 units in one region",
            "Topo_shape": "rooted directed unlabeled unordered AHU-style canonical tree shape",
            "VAF_winner": "family-specific read-AF softmax top weight >=0.60 at temperature 0.05; exploratory",
        },
        "sources": {
            "run_root": str(args.run_root),
            "manifest": str(args.manifest),
            "topology_catalog": str(args.topology),
            "regional_topology_catalog": str(args.regional_topology),
            "read_AF_ordering": str(args.read_af),
            "VAF_teaching_example": str(args.vaf_example),
        },
        "HCC1395_example": example,
        "HCC1395_VAF_multi_T_example": vaf_example,
        "samples": samples,
        "aggregate": aggregate,
        "shape_catalog": shapes,
        "validation": {
            "status": "PASS" if not validation_failures else "FAIL",
            "failures": validation_failures,
        },
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "c_t_topology_report_data.json").write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    sample_rows = []
    for sample in samples:
        sample_rows.append(
            {
                "sample": sample["sample"],
                "biological_id": sample["biological_id"],
                "S": sample["S"]["all_contigs"],
                "S_autosomal": sample["S"]["autosomal_chr1_22"],
                "W_all": sample["W"]["all_pre"],
                "W_k1": sample["W"]["k1"],
                "W_k_gt1": sample["W"]["k_gt1"],
                "W_tree": sample["W"]["tree_view"],
                "W_primary": sample["W"]["primary"],
                "primary_units": sample["W"]["primary_units"],
                "C_pooled_0": sample["C"]["pooled_region"]["0"],
                "C_pooled_1": sample["C"]["pooled_region"]["1"],
                "C_pooled_2": sample["C"]["pooled_region"]["2"],
                "C_pooled_3plus": sum(
                    sample["C"]["pooled_region"].get(key, 0)
                    for key in ("3", "4", "5", "6", ">6")
                ),
                "T1_Topo1": sample["T_and_Topology"]["classes"].get("T=1|Topo=1", 0),
                "Tn_Topo1": sample["T_and_Topology"]["classes"].get("T>1|Topo=1", 0),
                "Tn_Topon": sample["T_and_Topology"]["classes"].get("T>1|Topo>1", 0),
                "T_incomplete": sample["T_and_Topology"]["classes"].get("incomplete", 0),
                "distinct_shapes": sample["T_and_Topology"]["shape_catalog"]["distinct_shapes"],
                "exact_candidate_trees": sample["T_and_Topology"]["shape_catalog"]["exact_candidate_trees"],
                "ordered_regional_shapes": sample["T_and_Topology"]["regional_ordered_forest"]["distinct_ordered_regional_shape_signatures"],
                "region_topology_alternatives": sample["T_and_Topology"]["regional_ordered_forest"]["sum_region_topology_alternatives"],
                "joint_exact_T_candidates": sample["T_and_Topology"]["regional_ordered_forest"]["sum_exact_joint_tree_candidates"],
                "VAF_ambiguous": sample["VAF_ranking"]["ambiguous_primary_units"],
                "VAF_prepared": sample["VAF_ranking"]["prepared_all_candidates"],
                "VAF_top_ge_0_60": sample["VAF_ranking"]["top_weight_ge_0_60"],
                "VAF_missing": sample["VAF_ranking"]["missing_read_AF"],
            }
        )
    write_tsv(
        args.output_dir / "all_sample_summary.tsv",
        sample_rows,
        list(sample_rows[0]),
    )

    shape_rows = []
    for row in shapes:
        shape_rows.append(
            {
                "display_id": row["display_id"],
                "stable_id": row["shape_id"],
                "signature": row["signature"],
                "coarse_shape": row["coarse_shape"],
                "mutation_state_nodes": row["n_mutation_state_nodes"],
                "max_depth": row["max_depth"],
                "root_degree": row["root_degree"],
                "leaves": row["n_leaves"],
                "internal_branch_nodes": row["n_internal_branch_nodes"],
                "unit_incidence": row["unit_incidence"],
                "shape_only_units": row["shape_only_units"],
                "candidate_trees": row["candidate_trees"],
                **{f"{sample}_units": row["by_sample_units"][sample] for sample in sample_order},
            }
        )
    write_tsv(
        args.output_dir / "topology_shape_catalog.tsv",
        shape_rows,
        list(shape_rows[0]),
    )

    pair_rows = [
        {"C_pair": key, "regions": value}
        for key, value in hcc["C"]["double_HP_pairs"].items()
    ]
    write_tsv(
        args.output_dir / "HCC1395_double_HP_C_pairs.tsv",
        pair_rows,
        ["C_pair", "regions"],
    )
    validation_rows = [
        {"sample": sample["sample"], "check": check, "pass": passed}
        for sample in samples
        for check, passed in sample["checks"].items()
    ]
    write_tsv(
        args.output_dir / "validation_checks.tsv",
        validation_rows,
        ["sample", "check", "pass"],
    )
    nested_rows = [
        {"sample": sample["sample"], **row}
        for sample in samples
        for row in sample["T_and_Topology"]["nested_region_cross"]
    ]
    write_tsv(
        args.output_dir / "primary_region_HP_C_T_Topo_H_cross.tsv",
        nested_rows,
        ["sample", "HP_scope", "C_state", "T_Topo_state", "H_state", "regions"],
    )
    print(
        f"REPORT DATA -> {args.output_dir}; validation={output['validation']['status']}; "
        f"samples={len(samples)} shapes={len(shapes)}"
    )
    if validation_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
