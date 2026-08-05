#!/usr/bin/env python3
"""Build the Markdown source and canonical portable-HTML artifact input."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
import subprocess
from datetime import datetime
from pathlib import Path


SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]

COARSE_ZH = {
    "single": "單節點",
    "linear": "鏈狀",
    "root_star": "根部分岔",
    "internal_branch": "內部分岔",
    "mixed_root_and_internal_branch": "根與內部皆分岔",
}


def pct(numerator, denominator):
    return numerator / denominator if denominator else None


def fmt(value):
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.2%}"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def md_escape(value):
    return fmt(value).replace("|", "\\|").replace("\n", " ")


def md_table(rows, columns):
    labels = [label for _, label in columns]
    output = [
        "| " + " | ".join(labels) + " |",
        "|" + "|".join("---" for _ in labels) + "|",
    ]
    for row in rows:
        output.append(
            "| " + " | ".join(md_escape(row.get(field)) for field, _ in columns) + " |"
        )
    return "\n".join(output)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_context():
    def run(*args):
        return subprocess.check_output(args, text=True).strip()

    try:
        return {
            "branch": run("git", "branch", "--show-current"),
            "commit": run("git", "rev-parse", "HEAD"),
            "dirty": bool(run("git", "status", "--porcelain")),
        }
    except (OSError, subprocess.CalledProcessError):
        return {"branch": "unknown", "commit": "unknown", "dirty": None}


def production_state(path, expected_samples):
    latest = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            latest[row["sample"]] = {
                "status": row["status"],
                "timestamp": row["timestamp"],
                "detail": row["detail"],
            }
    passes = sum(value["status"] == "PASS" for value in latest.values())
    active = [sample for sample, value in latest.items() if value["status"] == "START"]
    return {
        "latest": latest,
        "passes": passes,
        "expected": len(expected_samples),
        "active": active,
        "not_started": [sample for sample in expected_samples if sample not in latest],
        "mtime": datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds"),
    }


def ordered_edges(edges):
    children = {}
    for parent, child in edges:
        children.setdefault(parent, []).append(child)
    ordered = []

    def visit(parent):
        for child in sorted(children.get(parent, [])):
            ordered.append(f"{parent}→{child}")
            visit(child)

    visit("ROOT")
    return "; ".join(ordered)


def columns(*items):
    return [{"field": field, "label": label, **options} for field, label, options in items]


def sql_literal(value):
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (int, float)):
        return repr(value)
    return "'" + str(value).replace("'", "''") + "'"


def freeze_dataset_with_sql(dataset_name, rows):
    """Execute a literal SQLite query and return both the query and its rows."""
    if not rows:
        raise RuntimeError(f"cannot freeze empty dataset: {dataset_name}")
    fields = list(rows[0])
    if any(list(row) != fields for row in rows):
        raise RuntimeError(f"dataset columns are not stable: {dataset_name}")
    identifiers = ", ".join(f'"{field}"' for field in fields)
    values = ",\n    ".join(
        "(" + ", ".join(sql_literal(row[field]) for field in fields) + ")"
        for row in rows
    )
    query = (
        f'WITH "{dataset_name}" ({identifiers}) AS (\n'
        f"  VALUES\n    {values}\n"
        f')\nSELECT {identifiers} FROM "{dataset_name}"'
    )
    connection = sqlite3.connect(":memory:")
    try:
        cursor = connection.execute(query)
        frozen = [dict(zip(fields, result)) for result in cursor.fetchall()]
    finally:
        connection.close()
    if len(frozen) != len(rows):
        raise RuntimeError(f"SQL row-count mismatch: {dataset_name}")
    return query, frozen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    parser.add_argument("--nested-cross", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--regional-topology", required=True, type=Path)
    parser.add_argument("--read-af", required=True, type=Path)
    parser.add_argument("--vaf-example", required=True, type=Path)
    parser.add_argument("--current-status", required=True, type=Path)
    parser.add_argument("--report-md", required=True, type=Path)
    parser.add_argument("--artifact", required=True, type=Path)
    args = parser.parse_args()

    data = json.loads(args.data.read_text(encoding="utf-8"))
    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    by_sample = {row["sample"]: row for row in data["samples"]}
    samples = [by_sample[sample] for sample in SAMPLE_ORDER]
    hcc = by_sample["HCC1395"]
    aggregate = data["aggregate"]
    simple_example = data["HCC1395_example"]
    vaf_example = data["HCC1395_VAF_multi_T_example"]
    current = production_state(args.current_status, SAMPLE_ORDER)
    git = git_context()

    hcc_core = [{
        "S_autosomal": hcc["S"]["autosomal_chr1_22"],
        "W_retained": hcc["W"]["retained"],
        "W_primary": hcc["W"]["primary"],
        "primary_units": hcc["W"]["primary_units"],
        "exact_T": hcc["T_and_Topology"]["shape_catalog"]["exact_candidate_trees"],
        "shape_types": hcc["T_and_Topology"]["shape_catalog"]["distinct_shapes"],
    }]
    global_core = [{
        "dataset_rows": aggregate["dataset_rows"],
        "biological_samples": aggregate["biological_samples"],
        "primary_units": aggregate["primary_units"],
        "exact_T": aggregate["exact_candidate_trees"],
        "shape_types": aggregate["distinct_shapes_global"],
        "validation_failures": len(data["validation"]["failures"]),
    }]

    hcc_scope = [
        {"stage": "S：輸入 sSNV（全 contigs）", "count": hcc["S"]["all_contigs"], "meaning": "VCF 中全部 sSNV"},
        {"stage": "S：chr1–22", "count": hcc["S"]["autosomal_chr1_22"], "meaning": "本次 W/C/T/Topo 分析母體"},
        {"stage": "W：k=1", "count": hcc["W"]["k1"], "meaning": "單點區域，不做多位點樹"},
        {"stage": "W：k>1 pre-cap", "count": hcc["W"]["k_gt1"], "meaning": "相鄰 gap≤50 kb 的多位點區域"},
        {"stage": "W：retained", "count": hcc["W"]["retained"], "meaning": "通過讀段支援並保留"},
        {"stage": "W：tree view", "count": hcc["W"]["tree_view"], "meaning": "至少產生一個 lineage unit"},
        {"stage": "W：primary HP", "count": hcc["W"]["primary"], "meaning": "HP1 xor HP2 或 HP1 and HP2"},
        {"stage": "primary HP units", "count": hcc["W"]["primary_units"], "meaning": "region×HP，雙 HP 區域算兩個 unit"},
    ]

    hp_labels = {
        "single_HP1_xor_HP2": "HP1 xor HP2",
        "double_HP1_and_HP2": "HP1 and HP2",
        "no_primary": "無 primary HP1/HP2",
    }
    hp_order = ["single_HP1_xor_HP2", "double_HP1_and_HP2", "no_primary"]
    hcc_hp = []
    hcc_hp_long = []
    for key in hp_order:
        without = hcc["HP"]["by_H3"].get(f"{key}|without_H3", 0)
        with_h3 = hcc["HP"]["by_H3"].get(f"{key}|with_H3", 0)
        hcc_hp.append({"class": hp_labels[key], "without_H3": without, "with_H3": with_h3, "total": without + with_h3})
        hcc_hp_long.extend([
            {"class": hp_labels[key], "H3": "not HP3", "count": without},
            {"class": hp_labels[key], "H3": "with HP3", "count": with_h3},
        ])

    c_keys = [str(value) for value in range(7)] + [">6"]
    hcc_c_unit = [{"C": key, "count": hcc["C"]["primary_HP_unit"].get(key, 0)} for key in c_keys]
    hcc_c_pooled = [{"C": key, "count": hcc["C"]["pooled_region"].get(key, 0)} for key in c_keys]
    hcc_c_pairs = [
        {"C_pair": key.replace(";", ", "), "regions": value}
        for key, value in hcc["C"]["double_HP_pairs"].items()
    ]

    topology_labels = {
        "T=1|Topo=1": "T=1；Topo=1",
        "T>1|Topo=1": "T>1；Topo=1",
        "T>1|Topo>1": "T>1；Topo>1",
        "incomplete": "候選集不完整",
    }
    topology_order = ["T=1|Topo=1", "T>1|Topo=1", "T>1|Topo>1", "incomplete"]
    hcc_topology = [
        {"class": topology_labels[key], "regions": hcc["T_and_Topology"]["classes"].get(key, 0)}
        for key in topology_order
    ]
    hcc_hidden = []
    for key in topology_order[:3]:
        values = hcc["T_and_Topology"]["hidden"]["by_topology"].get(key, {})
        hcc_hidden.append({
            "class": topology_labels[key],
            "hidden_0": values.get("hidden=0", 0),
            "hidden_positive": values.get("hidden>0", 0),
            "total": values.get("hidden=0", 0) + values.get("hidden>0", 0),
        })

    simple_groups = [
        {"genotype": genotype, "reads": reads, "counts_as_C": "是（含 A）" if "A" in genotype else "否"}
        for genotype, reads in simple_example["read_groups"].items()
    ]
    simple_tree = [{
        "region": simple_example["region"],
        "HP": f"HP{simple_example['HP']}",
        "k": simple_example["k"],
        "C": simple_example["C"],
        "T": simple_example["T"],
        "Topo": simple_example["Topo"],
        "H": simple_example["hidden"],
        "edges": ordered_edges(simple_example["edges"]),
    }]

    vaf_sites = [
        {
            "site": f"位點 {row['index'] + 1}",
            "position": row["position"],
            "REF_reads": row["REF_reads"],
            "ALT_reads": row["ALT_reads"],
            "read_AF": row["read_AF"],
        }
        for row in vaf_example["sites"]
    ]
    vaf_candidates = [
        {
            "candidate": row["candidate"],
            "tree": ordered_edges(row["edges"]),
            "shape": row["shape_signature"],
            "score": row["score"],
            "score_text": f"{row['score']:.6f}",
            "weight": row["softmax_weight"],
            "selected": "VAF-supported top" if row["selected"] else "否",
        }
        for row in vaf_example["candidates"]
    ]

    shapes = []
    for row in data["shape_catalog"]:
        shapes.append({
            "Topo_ID": row["display_id"],
            "stable_id": row["shape_id"],
            "signature": row["signature"],
            "shape_family": COARSE_ZH.get(row["coarse_shape"], row["coarse_shape"]),
            "nodes_with_ROOT": row["n_nodes_including_root"],
            "mutation_state_nodes": row["n_mutation_state_nodes"],
            "max_depth": row["max_depth"],
            "root_degree": row["root_degree"],
            "leaves": row["n_leaves"],
            "unit_support": row["unit_incidence"],
            "shape_only_units": row["shape_only_units"],
            "candidate_T_occurrences": row["candidate_trees"],
            "HCC1395": row["by_sample_units"]["HCC1395"],
            "HCC1395_DORADO": row["by_sample_units"]["HCC1395_DORADO"],
            "COLO829": row["by_sample_units"]["COLO829"],
            "H1437": row["by_sample_units"]["H1437"],
            "H2009": row["by_sample_units"]["H2009"],
            "HCC1937": row["by_sample_units"]["HCC1937"],
            "HCC1954": row["by_sample_units"]["HCC1954"],
        })
    hcc_top_shapes = sorted(
        [row for row in shapes if row["HCC1395"]],
        key=lambda row: (-row["HCC1395"], row["Topo_ID"]),
    )[:10]
    hcc_top_shapes_chart = [
        {
            "shape": f"{row['Topo_ID']} · {row['signature']}",
            "unit_support": row["HCC1395"],
            "global_unit_support": row["unit_support"],
        }
        for row in hcc_top_shapes
    ]

    all_scope = []
    all_hp = []
    all_c_pooled = []
    all_c_unit = []
    all_topology = []
    all_regional = []
    all_vaf = []
    for sample in samples:
        name = sample["sample"]
        all_scope.append({
            "sample": name,
            "S_all": sample["S"]["all_contigs"],
            "S_chr1_22": sample["S"]["autosomal_chr1_22"],
            "W_k1": sample["W"]["k1"],
            "W_k_gt1": sample["W"]["k_gt1"],
            "W_retained": sample["W"]["retained"],
            "W_tree": sample["W"]["tree_view"],
            "W_primary": sample["W"]["primary"],
            "primary_units": sample["W"]["primary_units"],
        })
        hp_row = {"sample": name}
        for key, prefix in (("single_HP1_xor_HP2", "single"), ("double_HP1_and_HP2", "double"), ("no_primary", "none")):
            hp_row[f"{prefix}_no_H3"] = sample["HP"]["by_H3"].get(f"{key}|without_H3", 0)
            hp_row[f"{prefix}_with_H3"] = sample["HP"]["by_H3"].get(f"{key}|with_H3", 0)
        all_hp.append(hp_row)
        pooled = {"sample": name}
        unit = {"sample": name}
        for key in c_keys:
            safe = "gt6" if key == ">6" else key
            pooled[f"C_{safe}"] = sample["C"]["pooled_region"].get(key, 0)
            unit[f"C_{safe}"] = sample["C"]["primary_HP_unit"].get(key, 0)
        all_c_pooled.append(pooled)
        all_c_unit.append(unit)
        classes = sample["T_and_Topology"]["classes"]
        shape_catalog = sample["T_and_Topology"]["shape_catalog"]
        all_topology.append({
            "sample": name,
            "T1_Topo1": classes.get("T=1|Topo=1", 0),
            "Tn_Topo1": classes.get("T>1|Topo=1", 0),
            "Tn_Topon": classes.get("T>1|Topo>1", 0),
            "incomplete_regions": classes.get("incomplete", 0),
            "hidden_0": sample["T_and_Topology"]["hidden"]["hidden=0"],
            "hidden_positive": sample["T_and_Topology"]["hidden"]["hidden>0"],
            "complete_units": shape_catalog["complete_primary_units"],
            "incomplete_units": shape_catalog["incomplete_primary_units"],
            "exact_T": shape_catalog["exact_candidate_trees"],
            "unit_shape_incidences": shape_catalog["unit_shape_incidence"],
            "shape_types": shape_catalog["distinct_shapes"],
        })
        regional = sample["T_and_Topology"]["regional_ordered_forest"]
        all_regional.append({
            "sample": name,
            "primary_regions": regional["primary_regions"],
            "complete_regions": regional["fully_complete_regions"],
            "incomplete_regions": regional["incomplete_regions"],
            "topology_alternatives": regional["sum_region_topology_alternatives"],
            "joint_exact_T": regional["sum_exact_joint_tree_candidates"],
            "ordered_shape_signatures": regional["distinct_ordered_regional_shape_signatures"],
        })
        vaf = sample["VAF_ranking"]
        neutral_prepared = vaf["prepared_by_CN"].get("neutral", 0)
        neutral_top = vaf["reached_by_CN"].get("neutral", 0)
        all_vaf.append({
            "sample": name,
            "ambiguous_units": vaf["ambiguous_primary_units"],
            "prepared": vaf["prepared_all_candidates"],
            "missing_AF": vaf["missing_read_AF"],
            "top_ge_0_60": vaf["top_weight_ge_0_60"],
            "reach_rate_prepared": pct(vaf["top_weight_ge_0_60"], vaf["prepared_all_candidates"]),
            "winner_consistent": vaf["winner_order_consistent"],
            "neutral_prepared": neutral_prepared,
            "neutral_top": neutral_top,
            "neutral_reach_rate": pct(neutral_top, neutral_prepared),
        })

    hp_scope_labels = {
        "single_HP1_xor_HP2": "單 HP（HP1 xor HP2）",
        "double_HP1_and_HP2": "雙 HP（HP1 and HP2）",
    }
    nested_state_labels = {
        "T=1|Topo=1": "T=1；Topo=1",
        "T>1|Topo=1": "T>1；Topo=1",
        "T>1|Topo>1": "T>1；Topo>1",
        "incomplete": "候選集不完整",
    }
    all_nested_cross = []
    for sample in samples:
        for row in sample["T_and_Topology"]["nested_region_cross"]:
            all_nested_cross.append(
                {
                    "sample": sample["sample"],
                    "HP_scope": hp_scope_labels[row["HP_scope"]],
                    "C_state": row["C_state"].replace(";", ", "),
                    "T_Topo_state": nested_state_labels[row["T_Topo_state"]],
                    "H_state": row["H_state"].replace("hidden", "H"),
                    "regions": row["regions"],
                }
            )
    hcc_nested_cross = [row for row in all_nested_cross if row["sample"] == "HCC1395"]

    validation_rows = []
    for sample in samples:
        checks = sample["checks"]
        validation_rows.append({
            "sample": sample["sample"],
            "checks": len(checks),
            "passed": sum(checks.values()),
            "failed": sum(not value for value in checks.values()),
            "status": "PASS" if all(checks.values()) else "FAIL",
        })

    provenance_paths = [
        ("corrected report data", args.data),
        ("validation checks", args.validation),
        ("nested HP/C/T/Topo/H cross", args.nested_cross),
        ("exact topology catalog", args.topology),
        ("regional topology composition", args.regional_topology),
        ("read-AF ordering", args.read_af),
        ("HCC1395 VAF example", args.vaf_example),
    ]
    provenance = [
        {"artifact": label, "path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}
        for label, path in provenance_paths
    ]

    datasets = {
        "hcc_core": hcc_core,
        "global_core": global_core,
        "hcc_scope": hcc_scope,
        "hcc_hp": hcc_hp,
        "hcc_hp_long": hcc_hp_long,
        "hcc_c_unit": hcc_c_unit,
        "hcc_c_pooled": hcc_c_pooled,
        "hcc_c_pairs": hcc_c_pairs,
        "hcc_topology": hcc_topology,
        "hcc_hidden": hcc_hidden,
        "simple_groups": simple_groups,
        "simple_tree": simple_tree,
        "vaf_sites": vaf_sites,
        "vaf_candidates": vaf_candidates,
        "hcc_top_shapes": hcc_top_shapes_chart,
        "all_scope": all_scope,
        "all_hp": all_hp,
        "all_c_pooled": all_c_pooled,
        "all_c_unit": all_c_unit,
        "all_topology": all_topology,
        "all_regional": all_regional,
        "all_vaf": all_vaf,
        "hcc_nested_cross": hcc_nested_cross,
        "all_nested_cross": all_nested_cross,
        "shape_catalog": shapes,
        "validation": validation_rows,
        "provenance": provenance,
    }

    dataset_source_paths = {
        "hcc_top_shapes": args.topology,
        "shape_catalog": args.topology,
        "all_regional": args.regional_topology,
        "vaf_sites": args.vaf_example,
        "vaf_candidates": args.vaf_example,
        "all_vaf": args.read_af,
        "validation": args.validation,
        "hcc_nested_cross": args.nested_cross,
        "all_nested_cross": args.nested_cross,
    }
    dataset_sources = []
    for dataset_name, rows in list(datasets.items()):
        query, frozen_rows = freeze_dataset_with_sql(dataset_name, rows)
        datasets[dataset_name] = frozen_rows
        upstream = dataset_source_paths.get(dataset_name, args.data)
        dataset_sources.append(
            {
                "id": f"dataset_{dataset_name}",
                "label": f"Frozen SQL snapshot: {dataset_name}",
                "path": str(upstream),
                "query": {
                    "engine": "sqlite",
                    "language": "sql",
                    "description": (
                        "Literal frozen snapshot executed by build_report_artifact.py; "
                        f"upstream calculation remains traceable to {upstream}."
                    ),
                    "executed_at": generated_at,
                    "sql": query,
                    "tables_used": [str(upstream)],
                },
            }
        )

    source_root = "research/20260711_read_group_C_tree_T_topology_report"
    sources = [
        {
            "id": "report_data",
            "label": "C/T/Topo corrected report dataset",
            "path": f"{source_root}/data/c_t_topology_report_data.json",
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Recomputes read-group C, HP classes, region-level T/Topo classes, hidden-state crosses, and joins exact topology and read-AF summaries.",
                "executed_at": data["generated_at"],
                "sql": "python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_report_dataset.py --run-root output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 --manifest research/20260710_layered_reconstruction_v2/input_manifest_lock.json --topology research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json --regional-topology research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json --read-af research/20260711_read_group_C_tree_T_topology_report/data/read_af_tree_ordering_historical.json --vaf-example research/20260711_read_group_C_tree_T_topology_report/data/HCC1395_VAF_multi_T_example.json --output-dir research/20260711_read_group_C_tree_T_topology_report/data",
                "tables_used": [
                    "output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/*/mlhp_part_*.json",
                    "output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/*/layered_reconstruction_*.json",
                ],
                "metric_definitions": [
                    "C_unit=count of ALT-containing MINREAD-supported full R/A groups within one primary HP family; ROOT, partial, and H_ nodes excluded",
                    "T_joint=product of exact candidate-tree counts across HP1/HP2 primary units in a region",
                    "Topo count=product of exact rooted-shape counts across HP1/HP2 primary units in a region",
                ],
            },
        },
        {
            "id": "topology_catalog",
            "label": "Exact rooted topology catalog",
            "path": f"{source_root}/data/exact_topology_catalog.json",
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Re-enumerates every complete primary HP candidate set with tree_cap=0 and canonicalizes rooted unordered shapes.",
                "executed_at": data["generated_at"],
                "sql": "python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_exact_topology_catalog.py --run-root output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 --output research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json",
                "tables_used": ["output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/*/layered_reconstruction_*.json"],
            },
        },
        {
            "id": "regional_topology",
            "label": "Ordered HP1/HP2 regional forest catalog",
            "path": f"{source_root}/data/regional_topology_composition.json",
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Forms ordered HP1/HP2 topology pairs and Cartesian products without joining the two haplotypes into one tree.",
                "executed_at": data["generated_at"],
                "sql": "python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_regional_topology_composition.py --topology research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json --output research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json",
                "tables_used": [f"{source_root}/data/exact_topology_catalog.json"],
            },
        },
        {
            "id": "read_af",
            "label": "Historical family-specific read-AF ranking rerun",
            "path": f"{source_root}/data/read_af_tree_ordering_historical.json",
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Scores all exact candidates using HP-specific per-site read-AF and reports sensitivity-grid reach.",
                "executed_at": data["generated_at"],
                "sql": "python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/read_af_tree_ordering_multisample.py --run-root output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 --input-manifest output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/input_manifest.json --output research/20260711_read_group_C_tree_T_topology_report/data/read_af_tree_ordering_historical.json",
                "tables_used": ["output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/*/mlhp_part_*.json"],
                "metric_definitions": ["read_AF=ALT_reads/(REF_reads+ALT_reads) within one HP family", "top candidate requires softmax weight >=0.60 at temperature=0.05"],
            },
        },
        {
            "id": "vaf_example",
            "label": "HCC1395 multi-T VAF teaching example",
            "path": f"{source_root}/data/HCC1395_VAF_multi_T_example.json",
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Re-enumerates one CN-neutral HCC1395 T=2 unit and preserves site counts, candidate edges, scores, and weights.",
                "executed_at": data["generated_at"],
                "sql": "python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_teaching_example.py --sample-dir output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/HCC1395 --sample HCC1395 --region chr4:40979546-40998095 --hp 2 --method-script-dir docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts --output research/20260711_read_group_C_tree_T_topology_report/data/HCC1395_VAF_multi_T_example.json",
                "tables_used": ["output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/HCC1395/mlhp_part_*.json"],
            },
        },
        {
            "id": "validation",
            "label": "Per-sample conservation checks",
            "path": f"{source_root}/data/validation_checks.tsv",
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Emits every per-sample boolean conservation check from the corrected report-data build.",
                "executed_at": data["generated_at"],
                "sql": "python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_report_dataset.py --run-root output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 --manifest research/20260710_layered_reconstruction_v2/input_manifest_lock.json --topology research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json --regional-topology research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json --read-af research/20260711_read_group_C_tree_T_topology_report/data/read_af_tree_ordering_historical.json --vaf-example research/20260711_read_group_C_tree_T_topology_report/data/HCC1395_VAF_multi_T_example.json --output-dir research/20260711_read_group_C_tree_T_topology_report/data",
                "tables_used": [f"{source_root}/data/c_t_topology_report_data.json"],
            },
        },
        {
            "id": "nested_cross",
            "label": "Primary-region HP × C × T/Topo × H long table",
            "path": f"{source_root}/data/primary_region_HP_C_T_Topo_H_cross.tsv",
        },
        {"id": "current_status", "label": "Normalized raw-all producer status", "path": "output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/run_status.tsv"},
        {"id": "tree_solver", "label": "Exact tree enumeration solver", "path": "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py"},
        {"id": "vaf_method", "label": "Family-specific read-AF ordering method", "path": "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/read_af_tree_ordering_multisample.py"},
    ]
    # Widget validation requires actual SQL. Preserve Python commands as prose/file
    # provenance, while the dataset_* sources below are the SQL actually executed
    # to materialize the portable snapshot rows.
    for source in sources:
        query = source.get("query")
        if query and query.get("language") != "sql":
            query.pop("sql", None)
    sources.extend(dataset_sources)

    cards = [
        {
            "id": "hcc_core_card",
            "description": "HCC1395 chr1–22 historical engineering snapshot.",
            "dataset": "hcc_core",
            "sourceId": "report_data",
            "metrics": [
                {"label": "chr1–22 sSNV (S)", "field": "S_autosomal", "format": "number"},
                {"label": "retained W", "field": "W_retained", "format": "number"},
                {"label": "primary regions", "field": "W_primary", "format": "number"},
                {"label": "primary HP units", "field": "primary_units", "format": "number"},
                {"label": "exact candidate T", "field": "exact_T", "format": "number"},
                {"label": "HCC topology shapes", "field": "shape_types", "format": "number"},
            ],
        },
        {
            "id": "global_core_card",
            "description": "Seven dataset rows; HCC1395 and DORADO share one biological sample.",
            "dataset": "global_core",
            "sourceId": "report_data",
            "metrics": [
                {"label": "dataset rows", "field": "dataset_rows", "format": "number"},
                {"label": "biological samples", "field": "biological_samples", "format": "number"},
                {"label": "primary HP units", "field": "primary_units", "format": "number"},
                {"label": "exact candidate T", "field": "exact_T", "format": "number"},
                {"label": "global shapes", "field": "shape_types", "format": "number"},
                {"label": "validation failures", "field": "validation_failures", "format": "number"},
            ],
        },
    ]

    charts = [
        {
            "id": "hcc_hp_chart",
            "title": "HCC1395 的 HP 分類與 HP3 交叉",
            "subtitle": "分母為 7,927 個 tree-view regions；HP3 是額外交叉標記。",
            "type": "stackedBar",
            "dataset": "hcc_hp_long",
            "sourceId": "report_data",
            "intent": "composition",
            "comparisonContext": {"denominator": "7,927 tree-view regions", "grain": "region", "unit": "regions"},
            "encodings": {
                "x": {"field": "class", "type": "nominal", "label": "Primary HP class"},
                "y": {"field": "count", "type": "quantitative", "label": "Regions"},
                "color": {"field": "H3", "type": "nominal", "label": "HP3"},
            },
            "valueFormat": "number",
        },
        {
            "id": "hcc_c_chart",
            "title": "HCC1395 每個 primary HP tree unit 的 C",
            "subtitle": "C=read 直接支持的 ALT-containing full genotype groups；ROOT、partial 與 H 節點不算。",
            "type": "bar",
            "dataset": "hcc_c_unit",
            "sourceId": "report_data",
            "intent": "distribution",
            "comparisonContext": {"denominator": "10,355 primary HP units", "grain": "region×HP", "unit": "units"},
            "encodings": {
                "x": {"field": "C", "type": "ordinal", "label": "C"},
                "y": {"field": "count", "type": "quantitative", "label": "Primary HP units"},
            },
            "valueFormat": "number",
        },
        {
            "id": "hcc_topology_chart",
            "title": "HCC1395 region-level T / Topo identifiability",
            "subtitle": "雙 HP region 以 ordered forest 的 Cartesian product 計算 joint T 與 joint Topo。",
            "type": "bar",
            "dataset": "hcc_topology",
            "sourceId": "report_data",
            "intent": "composition",
            "comparisonContext": {"denominator": "7,590 primary regions", "grain": "region", "unit": "regions"},
            "encodings": {
                "x": {"field": "class", "type": "nominal", "label": "Identifiability class"},
                "y": {"field": "regions", "type": "quantitative", "label": "Regions"},
            },
            "valueFormat": "number",
        },
        {
            "id": "vaf_candidates_chart",
            "title": "HCC1395 多樹例子的 read-AF relative weight",
            "subtitle": "Softmax weight 是固定 heuristic 下的相對權重，不是校準 posterior。",
            "type": "bar",
            "dataset": "vaf_candidates",
            "sourceId": "vaf_example",
            "intent": "comparison",
            "comparisonContext": {"grain": "candidate exact T", "normalization": "softmax at temperature 0.05", "unit": "relative weight"},
            "encodings": {
                "x": {"field": "candidate", "type": "nominal", "label": "Candidate"},
                "y": {"field": "weight", "type": "quantitative", "label": "Relative weight", "format": "percent"},
            },
            "valueFormat": "percent",
        },
        {
            "id": "hcc_shapes_chart",
            "title": "HCC1395 最常見的 10 種 unit topology shapes",
            "subtitle": "以 unit support 排序；同一 unit 可支持多種 shape，因此不是互斥組成。",
            "type": "horizontalBar",
            "dataset": "hcc_top_shapes",
            "sourceId": "topology_catalog",
            "intent": "comparison",
            "comparisonContext": {"denominator": "9,702 complete HCC primary HP units", "grain": "region×HP", "unit": "unit supports"},
            "encodings": {
                "x": {"field": "shape", "type": "nominal", "label": "Topology shape"},
                "y": {"field": "unit_support", "type": "quantitative", "label": "HCC1395 unit support"},
            },
            "valueFormat": "number",
        },
        {
            "id": "all_vaf_chart",
            "title": "各 dataset 的 VAF-supported top reach",
            "subtitle": "分母為有完整逐位點 read-AF 的 ambiguous primary units；不是全部 regions。",
            "type": "horizontalBar",
            "dataset": "all_vaf",
            "sourceId": "read_af",
            "intent": "comparison",
            "comparisonContext": {"denominator": "prepared ambiguous primary units", "grain": "region×HP", "unit": "rate"},
            "encodings": {
                "x": {"field": "sample", "type": "nominal", "label": "Dataset"},
                "y": {"field": "reach_rate_prepared", "type": "quantitative", "label": "Top weight ≥0.60", "format": "percent"},
            },
            "valueFormat": "percent",
            "referenceLines": [{"axis": "y", "value": 0.60, "label": "60% descriptive reference", "lineStyle": "dashed", "color": "neutral"}],
        },
    ]

    tables = [
        {"id": "hcc_scope_table", "title": "HCC1395：S → W → primary units", "subtitle": "不同列的 grain 不同，不能直接當同一漏斗比例。", "dataset": "hcc_scope", "sourceId": "report_data", "density": "dense", "columns": columns(("stage", "層級", {"type": "text"}), ("count", "數量", {"format": "number"}), ("meaning", "意義", {"type": "text"}))},
        {"id": "simple_groups_table", "title": "教學例 1：read groups 如何得到 C", "subtitle": simple_example["region"], "dataset": "simple_groups", "sourceId": "report_data", "density": "dense", "columns": columns(("genotype", "R/A group", {"type": "text"}), ("reads", "支持 reads", {"format": "number"}), ("counts_as_C", "是否計入 C", {"type": "text"}))},
        {"id": "simple_tree_table", "title": "教學例 1：唯一樹", "subtitle": "ROOT → AR → AA；C=2、T=1、Topo=1、H=0。", "dataset": "simple_tree", "sourceId": "report_data", "density": "dense", "columns": columns(("region", "Region", {"type": "text"}), ("HP", "HP", {"type": "text"}), ("k", "k", {"format": "number"}), ("C", "C", {"format": "number"}), ("T", "T", {"format": "number"}), ("Topo", "Topo", {"format": "number"}), ("H", "H", {"format": "number"}), ("edges", "Edges", {"type": "text"}))},
        {"id": "vaf_sites_table", "title": "教學例 2：每個位點的 HP-specific read-AF", "subtitle": f"{vaf_example['region']} · HP{vaf_example['HP']} · CN={vaf_example['CN']}", "dataset": "vaf_sites", "sourceId": "vaf_example", "density": "dense", "columns": columns(("site", "位點", {"type": "text"}), ("position", "座標", {"format": "number"}), ("REF_reads", "REF reads", {"format": "number"}), ("ALT_reads", "ALT reads", {"format": "number"}), ("read_AF", "read-AF", {"format": "percent"}))},
        {"id": "vaf_candidates_table", "title": "教學例 2：兩棵 T 的分數與相對權重", "subtitle": "兩棵 exact T 都是同一 chain topology；VAF 決定哪個 mutation 先。", "dataset": "vaf_candidates", "sourceId": "vaf_example", "density": "dense", "columns": columns(("candidate", "Candidate", {"type": "text"}), ("tree", "Exact T", {"type": "text"}), ("shape", "Shape", {"type": "text"}), ("score", "Score", {"format": "number"}), ("weight", "Relative weight", {"format": "percent"}), ("selected", "判定", {"type": "text"}))},
        {"id": "hcc_hp_table", "title": "HCC1395 HP × HP3 完整交叉表", "subtitle": "single + double + no primary = 7,927 tree-view regions。", "dataset": "hcc_hp", "sourceId": "report_data", "density": "dense", "columns": columns(("class", "Primary HP class", {"type": "text"}), ("without_H3", "not HP3", {"format": "number"}), ("with_H3", "with HP3", {"format": "number"}), ("total", "合計", {"format": "number"}))},
        {"id": "hcc_c_pairs_table", "title": "HCC1395 雙 HP 的 ordered C pair", "subtitle": "保留 HP1/HP2 身分；不可只看 C_total。", "dataset": "hcc_c_pairs", "sourceId": "report_data", "density": "dense", "columns": columns(("C_pair", "(C_HP1, C_HP2)", {"type": "text"}), ("regions", "Regions", {"format": "number"}))},
        {"id": "hcc_nested_cross_table", "title": "HCC1395：單/雙 HP × C × H × T/Topo 完整巢狀交叉", "subtitle": "Region grain；單 HP 保留 HP1/HP2，雙 HP 保留 ordered (C_HP1,C_HP2)。雙 HP 的 H=0 表示兩棵皆 H=0；H>0 表示至少一棵 H>0。71 rows 合計 7,590 primary regions。", "dataset": "hcc_nested_cross", "sourceId": "nested_cross", "density": "dense", "columns": columns(("HP_scope", "HP scope", {"type": "text"}), ("C_state", "C / ordered C pair", {"type": "text"}), ("H_state", "H", {"type": "text"}), ("T_Topo_state", "T / Topo", {"type": "text"}), ("regions", "Regions", {"format": "number"}))},
        {"id": "hcc_hidden_table", "title": "HCC1395 Topo class × extra-state H", "subtitle": "H_ 同時混合未觀測 intermediate 與 partial-supported completion，不能直接叫 hidden clone。", "dataset": "hcc_hidden", "sourceId": "report_data", "density": "dense", "columns": columns(("class", "T / Topo class", {"type": "text"}), ("hidden_0", "H=0", {"format": "number"}), ("hidden_positive", "H>0", {"format": "number"}), ("total", "合計", {"format": "number"}))},
        {"id": "all_scope_table", "title": "全 7 dataset：S 與 W", "subtitle": "HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 dataset rows。", "dataset": "all_scope", "sourceId": "report_data", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), ("S_all", "S all", {"format": "number"}), ("S_chr1_22", "S chr1–22", {"format": "number"}), ("W_k1", "W k=1", {"format": "number"}), ("W_k_gt1", "W k>1", {"format": "number"}), ("W_retained", "W retained", {"format": "number"}), ("W_tree", "W tree", {"format": "number"}), ("W_primary", "W primary", {"format": "number"}), ("primary_units", "Primary units", {"format": "number"}))},
        {"id": "all_hp_table", "title": "全 7 dataset：HP 分類 × HP3", "subtitle": "Region grain。", "dataset": "all_hp", "sourceId": "report_data", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), ("single_no_H3", "xor / not H3", {"format": "number"}), ("single_with_H3", "xor / H3", {"format": "number"}), ("double_no_H3", "and / not H3", {"format": "number"}), ("double_with_H3", "and / H3", {"format": "number"}), ("none_no_H3", "no primary / not H3", {"format": "number"}), ("none_with_H3", "no primary / H3", {"format": "number"}))},
        {"id": "all_c_pooled_table", "title": "全 7 dataset：pooled region C（raw census）", "subtitle": "跨 HP pooled；不可直接對應某一棵 per-HP tree。", "dataset": "all_c_pooled", "sourceId": "report_data", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), *((f"C_{key}", f"C={key if key != 'gt6' else '>6'}", {"format": "number"}) for key in ["0", "1", "2", "3", "4", "5", "6", "gt6"]))},
        {"id": "all_c_unit_table", "title": "全 7 dataset：per-primary-HP unit C（tree-aligned）", "subtitle": "這是應與 T 對齊的 C 主表。", "dataset": "all_c_unit", "sourceId": "report_data", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), *((f"C_{key}", f"C={key if key != 'gt6' else '>6'}", {"format": "number"}) for key in ["0", "1", "2", "3", "4", "5", "6", "gt6"]))},
        {"id": "all_topology_table", "title": "全 7 dataset：T / Topo 與 unit shape catalog", "subtitle": "Region identifiability 與 unit shape census 同表但不同 grain，欄位明示。", "dataset": "all_topology", "sourceId": "report_data", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), ("T1_Topo1", "T1/Topo1 regions", {"format": "number"}), ("Tn_Topo1", "Tn/Topo1 regions", {"format": "number"}), ("Tn_Topon", "Tn/Topon regions", {"format": "number"}), ("incomplete_regions", "Incomplete regions", {"format": "number"}), ("hidden_0", "H=0 regions", {"format": "number"}), ("hidden_positive", "H>0 regions", {"format": "number"}), ("complete_units", "Complete units", {"format": "number"}), ("incomplete_units", "Incomplete units", {"format": "number"}), ("exact_T", "Exact T (unit candidates)", {"format": "number"}), ("unit_shape_incidences", "Unit-shape incidences", {"format": "number"}), ("shape_types", "Shape types", {"format": "number"}))},
        {"id": "all_nested_cross_table", "title": "全 7 dataset：單/雙 HP × C × H × T/Topo 長表", "subtitle": "457 rows；每個 dataset 的 regions 合計各自回到 W_primary。", "dataset": "all_nested_cross", "sourceId": "nested_cross", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), ("HP_scope", "HP scope", {"type": "text"}), ("C_state", "C / ordered C pair", {"type": "text"}), ("H_state", "H", {"type": "text"}), ("T_Topo_state", "T / Topo", {"type": "text"}), ("regions", "Regions", {"format": "number"}))},
        {"id": "all_regional_table", "title": "全 7 dataset：ordered HP1/HP2 regional forest", "subtitle": "雙 HP 不是一棵連通樹；F_region=(T_HP1,T_HP2)。", "dataset": "all_regional", "sourceId": "regional_topology", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), ("primary_regions", "Primary regions", {"format": "number"}), ("complete_regions", "Complete", {"format": "number"}), ("incomplete_regions", "Incomplete", {"format": "number"}), ("topology_alternatives", "Ordered Topo alternatives", {"format": "number"}), ("joint_exact_T", "Joint exact-T candidates", {"format": "number"}), ("ordered_shape_signatures", "Distinct ordered signatures", {"format": "number"}))},
        {"id": "all_vaf_table", "title": "全 7 dataset：VAF-supported top candidate", "subtitle": "Primary HP unit grain；僅是固定 read-AF heuristic。", "dataset": "all_vaf", "sourceId": "read_af", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), ("ambiguous_units", "Ambiguous", {"format": "number"}), ("prepared", "有完整 AF", {"format": "number"}), ("missing_AF", "缺 AF", {"format": "number"}), ("top_ge_0_60", "Top≥0.60", {"format": "number"}), ("reach_rate_prepared", "Reach / prepared", {"format": "percent"}), ("winner_consistent", "方向自洽", {"format": "number"}), ("neutral_prepared", "CN-neutral prepared", {"format": "number"}), ("neutral_top", "CN-neutral top", {"format": "number"}), ("neutral_reach_rate", "Neutral reach", {"format": "percent"}))},
        {"id": "shape_catalog_table", "title": "完整 46 種 rooted topology shapes", "subtitle": "Topo_ID 依全域 unit support 排名；candidate occurrences 會被高度歧義 unit 放大。", "dataset": "shape_catalog", "sourceId": "topology_catalog", "density": "dense", "defaultSort": {"field": "unit_support", "direction": "desc"}, "columns": columns(("Topo_ID", "TopoShape-ID", {"type": "text"}), ("stable_id", "Stable hash", {"type": "text"}), ("signature", "Canonical signature", {"type": "text"}), ("shape_family", "形狀族", {"type": "text"}), ("nodes_with_ROOT", "Nodes incl. ROOT", {"format": "number"}), ("max_depth", "Depth", {"format": "number"}), ("root_degree", "Root degree", {"format": "number"}), ("leaves", "Leaves", {"format": "number"}), ("unit_support", "Global unit support", {"format": "number"}), ("shape_only_units", "Shape-only units", {"format": "number"}), ("candidate_T_occurrences", "Candidate T occurrences", {"format": "number"}), ("HCC1395", "HCC1395", {"format": "number"}), ("HCC1395_DORADO", "DORADO", {"format": "number"}), ("COLO829", "COLO829", {"format": "number"}), ("H1437", "H1437", {"format": "number"}), ("H2009", "H2009", {"format": "number"}), ("HCC1937", "HCC1937", {"format": "number"}), ("HCC1954", "HCC1954", {"format": "number"}))},
        {"id": "validation_table", "title": "7/7 dataset 守恆檢查", "subtitle": "C field contract、MINREAD、HP、T/Topo、regional forest、read-AF candidate count。", "dataset": "validation", "sourceId": "validation", "density": "dense", "columns": columns(("sample", "Dataset", {"type": "text"}), ("checks", "Checks", {"format": "number"}), ("passed", "Passed", {"format": "number"}), ("failed", "Failed", {"format": "number"}), ("status", "Status", {"type": "text"}))},
        {"id": "provenance_table", "title": "資料產物與 SHA-256", "subtitle": "本報告引用的 frozen engineering artifacts。", "dataset": "provenance", "sourceId": "report_data", "density": "dense", "columns": columns(("artifact", "Artifact", {"type": "text"}), ("path", "Path", {"type": "text"}), ("sha256", "SHA-256", {"type": "text"}), ("bytes", "Bytes", {"format": "number"}))},
    ]

    for spec in [*cards, *charts, *tables]:
        spec["sourceId"] = f"dataset_{spec['dataset']}"

    scope_body = f"""## 🔴 PARTIAL / SCIENTIFIC NO-GO：完整 7-row 歷史工程 census，不是 clean-v3 終版

截至 `{current['mtime']}`，normalized raw-all producer 為 **{current['passes']}/{current['expected']} PASS**；目前 status ledger 的非終態樣本為 **{', '.join(current['active']) or '無'}**，尚未開始為 **{', '.join(current['not_started']) or '無'}**。clean layered-v3 的 7/7 aggregate 驗證仍未產生。

因此本 HTML 可用於 **C/T/Topo 定義、HCC1395 教學、歷史 7-row 工程數量與拓撲形狀 census**；不可升級成 clean-v3 biological validation 或真實 clone phylogeny 結論。"""

    summary_body = f"""## TL;DR

- `C` 已修正為 **reads 支持的 mutation-bearing full genotype 群數**，不再表示候選樹組合數。
- `T` 表示 exact rooted directed Steiner tree；雙 HP 是 ordered forest，joint exact-T 候選數為 Cartesian product。
- HCC1395 有 **{hcc['W']['primary']:,} primary regions / {hcc['W']['primary_units']:,} primary HP units**；complete units 中共枚舉 **{hcc['T_and_Topology']['shape_catalog']['exact_candidate_trees']:,} exact T**，看到 **{hcc['T_and_Topology']['shape_catalog']['distinct_shapes']} 種 unit tree shapes**。
- 全 7 dataset rows 有 **{aggregate['exact_candidate_trees']:,} exact unit-tree candidates、{aggregate['distinct_shapes_global']} 種 rooted shapes**；HCC1395 的 ordered regional forest 另有 **{hcc['T_and_Topology']['regional_ordered_forest']['distinct_ordered_regional_shape_signatures']} 種 signatures**。
- HCC1395 的 {hcc['VAF_ranking']['ambiguous_primary_units']:,} ambiguous units 中，{hcc['VAF_ranking']['prepared_all_candidates']:,} 有完整逐位點 read-AF；**{hcc['VAF_ranking']['top_weight_ge_0_60']:,}/{hcc['VAF_ranking']['prepared_all_candidates']:,} = {pct(hcc['VAF_ranking']['top_weight_ge_0_60'], hcc['VAF_ranking']['prepared_all_candidates']):.2%}** 達 top weight≥0.60。這是「VAF-supported most-likely candidate」，不是獨立確認的真樹。"""

    notation_body = """## 先把 C、T、Topo 三件事分清楚

`C_(r,h)`：region `r`、primary HP family `h` 中，`populations_by_hp[h]` 內 **含 A、完整跨越 k 個位點、且 read count≥MINREAD=3** 的 R/A genotype keys 數。`ROOT=RR…R`、partial R/A/X constraints、以及額外 `H_` 節點都不算 C。

對一棵具體 HP tree：

```text
|V(T)| = 1 ROOT + C + H_T
```

`T`：一棵 exact rooted directed Steiner arborescence candidate。`N_T,unit` 是一個 region×HP unit 的 exact T 數；雙 HP region 是保留 HP 身分的 ordered forest：

```text
F_region = (T_HP1, T_HP2)
N_T,region = N_T,HP1 × N_T,HP2
```

`τ(T)`：移除 mutation label 與 sibling 順序、但保留 ROOT 與方向後的形狀。`Topo_1 / Topo_n` 是「候選集合內有 1 / >1 種形狀」；`TopoShape-ID` 則是全資料實際看到哪一種 canonical shape，兩者不可混稱。邏輯上 `T=1 ⇒ Topo=1`，所以 `T=1 / Topo>1` 不可能。"""

    vaf_body = f"""## 為何逐位點 VAF 可以在多棵 T 中找 top candidate？

在無 copy-number／purity／multiplicity confounding、沒有回突變的理想情況，descendant mutation 所在細胞集合應是 ancestor mutation 的子集合，因此通常預期 `VAF_ancestor ≥ VAF_descendant`。本流程用同一 HP family 內 BAM reads 重算：

```text
r_j = ALT_reads_j / (REF_reads_j + ALT_reads_j)
Score(T) = Σ(parent→child) Σ(ancestor i in parent) (r_i − r_new-child)
weight(T) ∝ exp(Score(T) / 0.05)
```

教學例 2 是 HCC1395 `{vaf_example['region']}`、HP{vaf_example['HP']}、CN={vaf_example['CN']}：兩棵 exact T 的 topology 都是 chain，但一棵假設位點 1 先發生，另一棵假設位點 2 先發生。位點 read-AF 為 **{vaf_sites[0]['read_AF']:.4f}** 與 **{vaf_sites[1]['read_AF']:.4f}**，因此 `ROOT→H_AR→AA` 得到正的 VAF 差，relative weight={vaf_example['winner_weight']:.9f}。

這能在**固定 heuristic 與候選集合內**指出最有支持的相對 ordering；但 softmax weight 未校準成 Bayesian posterior，而且 raw read-AF 不是 purity/CN-corrected CCF，所以安全用語是 **VAF-supported most-likely candidate**。"""

    shape_body = """## 拓撲形狀怎麼看？

Canonical rule 是 `τ(v) = "(" + sort(τ(children)) + ")"`：

```text
(())        ROOT──●             單一 mutation-state node
((()))      ROOT──●──●          linear chain
(()())      ROOT──┬─●           root star（兩個 leaves）
                  └─●
((()()))    ROOT──●──┬─●        internal branch
                     └─●
```

完整 catalog 以 **unit support** 排名，因為 candidate occurrences 會讓高歧義 unit 被重複加權；同一 unit 可支持多種 shape，所以 unit support 跨 shape 不互斥。`Shape-only units` 表示該 unit 的完整候選集合只有這一種 topology shape。相同 graph shape 可以跨不同 k 出現，只代表 branching skeleton 相同，不代表相同突變或相同生物歷史。Incomplete/capped units 未列入 exact catalog，因此 46 種不是理論上限。"""

    limitations_body = """## 科學界線與限制

1. 這批數字來自 historical layered-v2 engineering snapshot；latest clean-v3 仍缺 7/7 aggregate validation。
2. `C=0` 只表示沒有達 MINREAD 的 full ALT genotype group；仍可能有 ALT partial-read subcube constraint 並產生 tree。
3. `H_` 目前同時包含真正未觀測 intermediate 與 partial-supported completion state，不能直接解讀為 hidden clone。
4. VAF ranking 是 unit-level、同資料 heuristic；winner direction consistency 不是獨立驗證。
5. Copy number、LOH、purity、mutation multiplicity 與 read-depth uncertainty 都可能改變 VAF ordering。HCC1395 的 top units 多數為 CN-altered，因此報告另列 CN-neutral 子集。
6. HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 dataset rows；7 rows 只代表 6 個 biological samples。"""

    reproduction_body = f"""## 可重現命令與版本

Repo：branch `{git['branch']}`；commit `{git['commit']}`；worktree dirty=`{git['dirty']}`。本任務只新增自己的 report workspace 與 in-progress report，不覆寫既有使用者變更。

```bash
python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_exact_topology_catalog.py \\
  --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 \\
  --output research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json

python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_regional_topology_composition.py \\
  --topology research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json \\
  --output research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json

python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_report_dataset.py [...]
```

Historical engineering snapshot 內部守恆：`{data['validation']['status']}`，failures=`{len(data['validation']['failures'])}`；clean-v3 scientific validation 仍 pending。"""

    blocks = [
        {"id": "scope", "type": "markdown", "body": scope_body, "sourceId": "current_status"},
        {"id": "summary", "type": "markdown", "body": summary_body, "sourceId": "report_data"},
        {"id": "hcc_metrics", "type": "metric-strip", "cardIds": ["hcc_core_card"]},
        {"id": "notation", "type": "markdown", "body": notation_body},
        {"id": "simple_groups", "type": "table", "tableId": "simple_groups_table"},
        {"id": "simple_tree", "type": "table", "tableId": "simple_tree_table"},
        {"id": "vaf_reason", "type": "markdown", "body": vaf_body, "sourceId": "vaf_example"},
        {"id": "vaf_sites", "type": "table", "tableId": "vaf_sites_table"},
        {"id": "vaf_candidates", "type": "table", "tableId": "vaf_candidates_table"},
        {"id": "vaf_chart", "type": "chart", "chartId": "vaf_candidates_chart"},
        {"id": "hcc_scope_heading", "type": "markdown", "body": "## HCC1395 完整數據"},
        {"id": "hcc_scope", "type": "table", "tableId": "hcc_scope_table"},
        {"id": "hcc_hp_chart_block", "type": "chart", "chartId": "hcc_hp_chart"},
        {"id": "hcc_hp_table_block", "type": "table", "tableId": "hcc_hp_table"},
        {"id": "hcc_c_chart_block", "type": "chart", "chartId": "hcc_c_chart"},
        {"id": "hcc_c_pairs_block", "type": "table", "tableId": "hcc_c_pairs_table"},
        {"id": "hcc_nested_cross_block", "type": "table", "tableId": "hcc_nested_cross_table"},
        {"id": "hcc_topology_chart_block", "type": "chart", "chartId": "hcc_topology_chart"},
        {"id": "hcc_hidden_block", "type": "table", "tableId": "hcc_hidden_table"},
        {"id": "shape_explain", "type": "markdown", "body": shape_body},
        {"id": "hcc_shapes_chart_block", "type": "chart", "chartId": "hcc_shapes_chart"},
        {"id": "all_heading", "type": "markdown", "body": "## 全部 7 dataset rows 數據"},
        {"id": "global_metrics", "type": "metric-strip", "cardIds": ["global_core_card"]},
        {"id": "all_scope_block", "type": "table", "tableId": "all_scope_table"},
        {"id": "all_hp_block", "type": "table", "tableId": "all_hp_table"},
        {"id": "all_c_pooled_block", "type": "table", "tableId": "all_c_pooled_table"},
        {"id": "all_c_unit_block", "type": "table", "tableId": "all_c_unit_table"},
        {"id": "all_topology_block", "type": "table", "tableId": "all_topology_table"},
        {"id": "all_nested_cross_block", "type": "table", "tableId": "all_nested_cross_table"},
        {"id": "all_regional_block", "type": "table", "tableId": "all_regional_table"},
        {"id": "all_vaf_chart_block", "type": "chart", "chartId": "all_vaf_chart"},
        {"id": "all_vaf_block", "type": "table", "tableId": "all_vaf_table"},
        {"id": "shape_catalog_block", "type": "table", "tableId": "shape_catalog_table"},
        {"id": "validation_heading", "type": "markdown", "body": "## 方法、驗證與限制"},
        {"id": "validation_block", "type": "table", "tableId": "validation_table"},
        {"id": "limitations", "type": "markdown", "body": limitations_body},
        {"id": "provenance_block", "type": "table", "tableId": "provenance_table"},
        {"id": "reproduction", "type": "markdown", "body": reproduction_body},
    ]

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "HCC1395 read 群 C、Steiner 樹 T 與拓撲形狀：全樣本數據報告",
            "description": "HCC1395 先行教學、C/T/Topo 修正版、逐位點 VAF 排序說明與 7 dataset rows 完整工程 census。",
            "generatedAt": generated_at,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": sources,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "partial",
            "datasets": datasets,
            "accessIssues": [
                {
                    "id": "clean_v3_pending",
                    "scope": "scientific validation",
                    "sourceId": "current_status",
                    "message": f"Normalized raw-all producer is {current['passes']}/{current['expected']} PASS at cutoff; clean layered-v3 7/7 aggregate output is unavailable. Historical v2 counts are complete for 7 dataset rows but remain engineering evidence.",
                },
                {
                    "id": "read_af_incomplete",
                    "dataset": "all_vaf",
                    "sourceId": "read_af",
                    "message": f"Read-AF is missing for {aggregate['VAF_missing']:,}/{aggregate['VAF_ambiguous_units']:,} ambiguous primary units; ranking completeness is false.",
                },
            ],
        },
        "sources": sources,
    }

    # Markdown companion: same answer-first path with complete tables.
    report = []
    report.append(
        "<!--\n"
        f"建立時間: {generated_at}\n"
        "目標: 校正 C/T/Topo 定義，以 HCC1395 教學後提供全 7 dataset rows 完整工程 census\n"
        "處理範圍: chr1–22；historical layered-v2；PARTIAL（clean layered-v3 7/7 pending）\n"
        "服務目標: G3 / G4 / G5\n"
        f"Git: branch={git['branch']}; commit={git['commit']}; dirty={git['dirty']}\n"
        "關聯檔案: InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json\n"
        "-->"
    )
    report.append("# HCC1395 read 群 C、Steiner 樹 T 與拓撲形狀：全樣本數據報告")
    report.append("用 SCQA + Feynman：先直答定義與數量，再用兩個 HCC1395 例子拆解，最後給全樣本表、方法與限制。")
    report.extend([scope_body, summary_body, notation_body])
    report.append("## HCC1395 教學例 1：C=2、T=1、Topo=1")
    report.append(md_table(simple_groups, [("genotype", "R/A group"), ("reads", "支持 reads"), ("counts_as_C", "計入 C")]))
    report.append(md_table(simple_tree, [("region", "Region"), ("HP", "HP"), ("k", "k"), ("C", "C"), ("T", "T"), ("Topo", "Topo"), ("H", "H"), ("edges", "Edges")]))
    report.append(vaf_body)
    report.append(md_table(vaf_sites, [("site", "位點"), ("position", "座標"), ("REF_reads", "REF reads"), ("ALT_reads", "ALT reads"), ("read_AF", "read-AF")]))
    report.append(md_table(vaf_candidates, [("candidate", "Candidate"), ("tree", "Exact T"), ("shape", "Shape"), ("score_text", "Score"), ("weight", "Weight"), ("selected", "判定")]))
    report.append("## HCC1395 完整數據")
    report.append(md_table(hcc_scope, [("stage", "層級"), ("count", "數量"), ("meaning", "意義")]))
    report.append("### HP × HP3")
    report.append(md_table(hcc_hp, [("class", "Primary HP class"), ("without_H3", "not H3"), ("with_H3", "with H3"), ("total", "合計")]))
    report.append("### C：tree-aligned per-HP unit")
    report.append(md_table(hcc_c_unit, [("C", "C"), ("count", "Primary HP units")]))
    report.append("### C：pooled region raw census")
    report.append(md_table(hcc_c_pooled, [("C", "C"), ("count", "Regions")]))
    report.append(f"Primary regions 內 pooled C 與 primary HP1/2 C_sum 不同：**{hcc['C']['pooled_vs_primary_HP_sum_discordant_primary_regions']:,}/{hcc['W']['primary']:,}**；全部 retained regions 對 HP1+HP2 為 **{hcc['C']['pooled_vs_HP1_HP2_sum_discordant_all_retained_regions']:,}/{hcc['W']['retained']:,}**；若把 HP3/none 等全部 family 都納入，則為 **{hcc['C']['pooled_vs_all_family_sum_discordant_all_retained_regions']:,}/{hcc['W']['retained']:,}**。差異來自跨 family pooling 的 collapse／threshold crossing，因此 pooled C 不可直接對 per-HP T。")
    report.append("### 雙 HP ordered C pair")
    report.append(md_table(hcc_c_pairs, [("C_pair", "(C_HP1,C_HP2)"), ("regions", "Regions")]))
    report.append("### T / Topo identifiability")
    report.append(md_table(hcc_topology, [("class", "Class"), ("regions", "Regions")]))
    report.append(md_table(hcc_hidden, [("class", "Class"), ("hidden_0", "H=0"), ("hidden_positive", "H>0"), ("total", "Total")]))
    report.append("### 單/雙 HP × C × H × T/Topo 完整巢狀交叉")
    report.append("此表是最初問題的直接答案：單 HP 顯示 `HP1=n` 或 `HP2=n`；雙 HP 顯示 ordered `HP1=n, HP2=m`，不以 C_sum 取代。雙 HP 的 `H=0` 表示兩棵都無 extra state；`H>0` 表示至少一棵有 extra state。")
    report.append(md_table(hcc_nested_cross, [("HP_scope", "HP scope"), ("C_state", "C / ordered pair"), ("H_state", "H"), ("T_Topo_state", "T / Topo"), ("regions", "Regions")]))
    regional = hcc["T_and_Topology"]["regional_ordered_forest"]
    report.append(f"HCC1395 complete primary regions={regional['fully_complete_regions']:,}；ordered topology alternatives 加總={regional['sum_region_topology_alternatives']:,}；joint exact-T candidate combinations 加總={regional['sum_exact_joint_tree_candidates']:,}；candidate catalog 枚舉到 ordered regional signatures={regional['distinct_ordered_regional_shape_signatures']:,}。")
    report.append(shape_body)
    report.append("### HCC1395 unit support 前 10 種 shape")
    report.append(md_table(hcc_top_shapes, [("Topo_ID", "TopoShape-ID"), ("signature", "Signature"), ("shape_family", "形狀族"), ("HCC1395", "HCC units"), ("unit_support", "Global units"), ("candidate_T_occurrences", "Global candidate T")]))
    report.append("## 全部 7 dataset rows 數據")
    report.append("### S / W")
    report.append(md_table(all_scope, [("sample", "Dataset"), ("S_all", "S all"), ("S_chr1_22", "S chr1–22"), ("W_k1", "W k=1"), ("W_k_gt1", "W k>1"), ("W_retained", "W retained"), ("W_tree", "W tree"), ("W_primary", "W primary"), ("primary_units", "Primary units")]))
    report.append("### HP × HP3")
    report.append(md_table(all_hp, [("sample", "Dataset"), ("single_no_H3", "xor/not H3"), ("single_with_H3", "xor/H3"), ("double_no_H3", "and/not H3"), ("double_with_H3", "and/H3"), ("none_no_H3", "none/not H3"), ("none_with_H3", "none/H3")]))
    report.append("### Pooled region C")
    report.append(md_table(all_c_pooled, [("sample", "Dataset"), *[(f"C_{key}", f"C={'>6' if key == 'gt6' else key}") for key in ["0", "1", "2", "3", "4", "5", "6", "gt6"]]]))
    report.append("### Per-primary-HP unit C")
    report.append(md_table(all_c_unit, [("sample", "Dataset"), *[(f"C_{key}", f"C={'>6' if key == 'gt6' else key}") for key in ["0", "1", "2", "3", "4", "5", "6", "gt6"]]]))
    report.append("### T / Topo / unit shapes")
    report.append(md_table(all_topology, [("sample", "Dataset"), ("T1_Topo1", "T1/Topo1 (regions)"), ("Tn_Topo1", "Tn/Topo1 (regions)"), ("Tn_Topon", "Tn/Topon (regions)"), ("incomplete_regions", "Incomplete (regions)"), ("hidden_0", "H=0 (regions)"), ("hidden_positive", "H>0 (regions)"), ("complete_units", "Complete units"), ("incomplete_units", "Incomplete units"), ("exact_T", "Exact T (unit candidates)"), ("unit_shape_incidences", "Unit-shape incidences"), ("shape_types", "Shape types")]))
    report.append("### 全 7 dataset 單/雙 HP × C × H × T/Topo 長表")
    report.append(md_table(all_nested_cross, [("sample", "Dataset"), ("HP_scope", "HP scope"), ("C_state", "C / ordered pair"), ("H_state", "H"), ("T_Topo_state", "T / Topo"), ("regions", "Regions")]))
    report.append("### Ordered regional forest")
    report.append(md_table(all_regional, [("sample", "Dataset"), ("primary_regions", "Primary R"), ("complete_regions", "Complete"), ("incomplete_regions", "Incomplete"), ("topology_alternatives", "Topo alternatives"), ("joint_exact_T", "Joint exact T"), ("ordered_shape_signatures", "Ordered signatures")]))
    report.append("### VAF-supported top candidate")
    report.append(md_table(all_vaf, [("sample", "Dataset"), ("ambiguous_units", "Ambiguous"), ("prepared", "Prepared"), ("missing_AF", "Missing AF"), ("top_ge_0_60", "Top≥.60"), ("reach_rate_prepared", "Reach/prepared"), ("winner_consistent", "Direction OK"), ("neutral_prepared", "Neutral prep"), ("neutral_top", "Neutral top"), ("neutral_reach_rate", "Neutral reach")]))
    report.append("## 完整 46 種 TopoShape catalog")
    report.append(md_table(shapes, [("Topo_ID", "TopoShape-ID"), ("stable_id", "Stable ID"), ("signature", "Signature"), ("shape_family", "形狀族"), ("nodes_with_ROOT", "Nodes incl. ROOT"), ("max_depth", "Depth"), ("root_degree", "Root degree"), ("leaves", "Leaves"), ("unit_support", "Global unit support"), ("shape_only_units", "Shape-only units"), ("candidate_T_occurrences", "Candidate T"), ("HCC1395", "HCC1395"), ("HCC1395_DORADO", "DORADO"), ("COLO829", "COLO829"), ("H1437", "H1437"), ("H2009", "H2009"), ("HCC1937", "HCC1937"), ("HCC1954", "HCC1954")]))
    report.append("## 驗證")
    report.append(md_table(validation_rows, [("sample", "Dataset"), ("checks", "Checks"), ("passed", "Passed"), ("failed", "Failed"), ("status", "Status")]))
    report.append(limitations_body)
    report.append("## Provenance")
    report.append(md_table(provenance, [("artifact", "Artifact"), ("path", "Path"), ("sha256", "SHA-256"), ("bytes", "Bytes")]))
    report.append(reproduction_body)
    report.append("---\n\nPARTIAL — historical 7-row engineering snapshot；clean layered-v3 7/7 pending。本文件不可作 final biological validation evidence。")

    args.report_md.parent.mkdir(parents=True, exist_ok=True)
    args.report_md.write_text("\n\n".join(report) + "\n", encoding="utf-8")
    args.artifact.parent.mkdir(parents=True, exist_ok=True)
    args.artifact.write_text(
        json.dumps(artifact, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"MARKDOWN -> {args.report_md}")
    print(f"ARTIFACT -> {args.artifact}")
    print(
        f"CURRENT PRODUCER -> pass={current['passes']}/{current['expected']} "
        f"active={','.join(current['active']) or 'none'} "
        f"not_started={','.join(current['not_started']) or 'none'}"
    )


if __name__ == "__main__":
    main()
