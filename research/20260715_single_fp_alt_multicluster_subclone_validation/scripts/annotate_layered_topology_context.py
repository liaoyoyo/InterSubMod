#!/usr/bin/env python3
"""Attach latest layered regional-tree context without upgrading single-site evidence."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def truth(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def tree_shape(tree: dict[str, Any]) -> str:
    edges = tree.get("edges") or []
    if tree.get("recurrence"):
        return "recurrence"
    if not edges:
        return "no_mutation_edge"
    if len(edges) == 1:
        return "trivial_one_edge"
    outdegree = Counter(parent for parent, _ in edges)
    if max(outdegree.values(), default=0) > 1:
        return "branching"
    return "linear_chain"


def load_detail_index(path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    index: dict[tuple[str, str], dict[str, Any]] = {}
    for detail in data["detail"]:
        key = (detail["region"], str(detail["family"]))
        if key in index:
            raise RuntimeError(f"Duplicate layered detail unit: {key}")
        index[key] = detail
    return index


def annotate_row(row: dict[str, str], index: dict[tuple[str, str], dict[str, Any]]) -> dict[str, Any]:
    output: dict[str, Any] = dict(row)
    region = row.get("selected_group_id") or ""
    try:
        alt_families = json.loads(row.get("layered_alt_families") or "[]")
    except json.JSONDecodeError:
        alt_families = []
    units = [index[(region, str(family))] for family in alt_families if (region, str(family)) in index]
    unit_records: list[dict[str, Any]] = []
    for unit in units:
        shapes = [tree_shape(tree) for tree in unit.get("trees") or []]
        unit_records.append(
            {
                "family": str(unit["family"]),
                "fam_label": unit.get("fam_label"),
                "unit_role": unit.get("unit_role"),
                "is_primary_lineage": unit.get("is_primary_lineage"),
                "is_h3_auxiliary": unit.get("is_h3_auxiliary"),
                "is_h4_auxiliary": unit.get("is_h4_auxiliary"),
                "L1_class": unit.get("L1_class"),
                "candidate_set_complete": unit.get("analysis_candidate_set_complete"),
                "capped": unit.get("capped"),
                "n_sSNV": unit.get("n_sSNV"),
                "n_trees": unit.get("n_trees"),
                "n_trees_stored": unit.get("n_trees_stored"),
                "stored_tree_shapes": shapes,
            }
        )
    primary = [unit for unit in unit_records if unit["is_primary_lineage"]]
    if row.get("ssnv_branch") != "retained" or not region:
        context = "NO_REGIONAL_TREE_FOR_SITE"
    elif not units:
        context = "NO_ALT_FAMILY_UNIT_MATCH"
    elif not primary:
        context = "AUXILIARY_OR_UNPHASED_ONLY"
    elif any(unit["capped"] or not unit["candidate_set_complete"] for unit in primary):
        context = "PRIMARY_CANDIDATE_SET_INCOMPLETE"
    else:
        shapes = {shape for unit in primary for shape in unit["stored_tree_shapes"]}
        if not shapes or shapes <= {"no_mutation_edge", "trivial_one_edge"}:
            context = "PRIMARY_TRIVIAL_ONE_EDGE_ONLY"
        elif "branching" in shapes or "recurrence" in shapes:
            context = "PRIMARY_INCLUDES_BRANCHING_OR_RECURRENCE"
        elif shapes <= {"linear_chain", "trivial_one_edge"}:
            context = "PRIMARY_ALL_STORED_TREES_LINEAR"
        else:
            context = "PRIMARY_TOPOLOGY_MIXED_OR_UNRESOLVED"
    output.update(
        {
            "layered_alt_family_units": json.dumps(unit_records, ensure_ascii=False, separators=(",", ":")),
            "layered_n_alt_family_units": len(unit_records),
            "layered_n_primary_alt_family_units": len(primary),
            "layered_topology_context": context,
            "focal_single_site_linear_topology_identifiable": False,
            "focal_single_site_subclone_confirmed_by_layered_tree": False,
        }
    )
    return output


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    subsets = {
        "all_sites": rows,
        "stable_null_multigroup": [row for row in rows if truth(row.get("stable_null_multigroup"))],
        "phase_anchored_robust_epigenetic_candidate": [
            row for row in rows if truth(row.get("phase_anchored_robust_epigenetic_candidate"))
        ],
    }
    return {
        name: {
            "n": len(values),
            "topology_context_counts": dict(Counter(row["layered_topology_context"] for row in values)),
            "focal_single_site_linear_topology_identifiable": 0,
            "focal_single_site_subclone_confirmed": 0,
        }
        for name, values in subsets.items()
    }


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-results", type=Path, required=True)
    parser.add_argument("--layered-root", type=Path, default=Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
        "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5"
    ))
    parser.add_argument("--output-tsv", type=Path, default=root / "results" / "focal_alt_multicluster" / "latest_site_results_with_topology.tsv")
    parser.add_argument("--summary", type=Path, default=root / "results" / "focal_alt_multicluster" / "latest_topology_context_summary.json")
    args = parser.parse_args()

    with args.site_results.open(encoding="utf-8") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in source_rows:
        by_sample[row["sample"]].append(row)
    output_rows: list[dict[str, Any]] = []
    for sample in DATASETS:
        if sample not in by_sample:
            continue
        path = args.layered_root / "samples" / sample / f"layered_reconstruction_{sample}.json"
        index = load_detail_index(path)
        annotated = [annotate_row(row, index) for row in by_sample[sample]]
        output_rows.extend(annotated)
        print(f"[{sample}] sites={len(annotated)} detail_units={len(index)}", flush=True)
    fields = sorted({key for row in output_rows for key in row})
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    summary = {
        "schema_name": "intersubmod.focal_alt_layered_topology_context",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "site_results": str(args.site_results),
        "layered_root": str(args.layered_root),
        "counts": summarize(output_rows),
        "guardrail": (
            "Regional tree context is same-pipeline corroboration. One focal sSNV supplies at most a "
            "trivial ROOT-to-ALT edge and cannot identify linear versus branching cellular ancestry."
        ),
        "pass": len(output_rows) == len(source_rows),
    }
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output_tsv), "summary": str(args.summary), "pass": summary["pass"]}, indent=2))


if __name__ == "__main__":
    main()
