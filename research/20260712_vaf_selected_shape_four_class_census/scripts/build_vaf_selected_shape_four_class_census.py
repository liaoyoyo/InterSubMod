#!/usr/bin/env python3
"""Materialize and classify final VAF-selected regional topology shapes.

The classifier operates on rooted-unlabeled mutation-state graph shapes.  It
does not infer clone counts.  Primary HP components remain ordered and are
collapsed to a region class by OR-combining sister/direct graph features.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
CATEGORIES = ["single_only", "sister_only", "direct_only", "sister_and_direct"]
CATEGORY_LABELS = {
    "single_only": "Single/no-within-HP-relation",
    "sister_only": "Sister-pattern",
    "direct_only": "Direct-chain-pattern",
    "sister_and_direct": "Sister+direct-pattern",
}
STRUCTURAL_TOPO1 = {"T=1|Topo=1", "T>1|Topo=1"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def truth(value: str | bool | None) -> bool:
    return value is True or value == "True"


def pct(numerator: int, denominator: int) -> float:
    return 100.0 * numerator / denominator if denominator else 0.0


def shape_category(shape: dict) -> str:
    sister = int(shape["max_outdegree"]) >= 2
    direct = int(shape["max_depth"]) >= 2
    if sister and direct:
        return "sister_and_direct"
    if sister:
        return "sister_only"
    if direct:
        return "direct_only"
    return "single_only"


def collapse_component_categories(categories: list[str]) -> tuple[str, bool, bool, bool]:
    sister = any(value in {"sister_only", "sister_and_direct"} for value in categories)
    direct = any(value in {"direct_only", "sister_and_direct"} for value in categories)
    split_across_hp = sister and direct and "sister_and_direct" not in categories
    if sister and direct:
        return "sister_and_direct", sister, direct, split_across_hp
    if sister:
        return "sister_only", sister, direct, split_across_hp
    if direct:
        return "direct_only", sister, direct, split_across_hp
    return "single_only", sister, direct, split_across_hp


def unique_rows(rows: list[dict[str, str]], keys: tuple[str, ...], name: str) -> dict[tuple[str, ...], dict]:
    result: dict[tuple[str, ...], dict] = {}
    for row in rows:
        key = tuple(row[value] for value in keys)
        if key in result:
            raise RuntimeError(f"duplicate {name} key: {key}")
        result[key] = row
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-catalog", type=Path, required=True)
    parser.add_argument("--ct-report", type=Path, required=True)
    parser.add_argument("--coarse-regions", type=Path, required=True)
    parser.add_argument("--vaf-units", type=Path, required=True)
    parser.add_argument("--vaf-regions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    exact_catalog = json.loads(args.exact_catalog.read_text(encoding="utf-8"))
    ct_report = json.loads(args.ct_report.read_text(encoding="utf-8"))
    coarse_rows = read_tsv(args.coarse_regions)
    vaf_unit_rows = read_tsv(args.vaf_units)
    vaf_region_rows = read_tsv(args.vaf_regions)

    shape_meta = {row["shape_id"]: row for row in ct_report["shape_catalog"]}
    shape_class = {shape_id: shape_category(row) for shape_id, row in shape_meta.items()}
    if len(shape_meta) != int(ct_report["aggregate"]["distinct_shapes_global"]):
        raise RuntimeError("shape catalog cardinality mismatch")

    coarse_by_key = unique_rows(coarse_rows, ("sample", "region"), "coarse region")
    vaf_region_by_key = unique_rows(vaf_region_rows, ("sample", "region"), "VAF region")
    vaf_unit_by_key = unique_rows(vaf_unit_rows, ("sample", "region", "family"), "VAF unit")
    if set(coarse_by_key) != set(vaf_region_by_key):
        raise RuntimeError("coarse/VAF region key sets differ")

    catalog_units: dict[tuple[str, str, str], dict] = {}
    catalog_units_by_region: dict[tuple[str, str], list[dict]] = collections.defaultdict(list)
    for sample in exact_catalog["samples"]:
        sample_name = sample["sample"]
        for unit in sample["unit_rows"]:
            row = {
                "sample": sample_name,
                "region": unit["region"],
                "family": str(unit["hp"]),
                "n_trees": int(unit["n_trees"]),
                "shape_ids": tuple(sorted(unit["shape_candidate_counts"])),
            }
            key = (row["sample"], row["region"], row["family"])
            if key in catalog_units:
                raise RuntimeError(f"duplicate exact-catalog unit: {key}")
            if not row["shape_ids"] or not set(row["shape_ids"]) <= set(shape_meta):
                raise RuntimeError(f"missing/unknown exact-catalog shape: {key}")
            catalog_units[key] = row
            catalog_units_by_region[(row["sample"], row["region"])].append(row)
    for rows in catalog_units_by_region.values():
        rows.sort(key=lambda row: int(row["family"]))

    selected_rows: list[dict] = []
    status_counts: dict[str, collections.Counter] = {
        sample: collections.Counter() for sample in SAMPLE_ORDER
    }
    structural_category_mismatches = 0
    resolved_exact_class_counts: dict[str, collections.Counter] = {
        sample: collections.Counter() for sample in SAMPLE_ORDER
    }
    direction_counts: dict[str, collections.Counter] = {
        sample: collections.Counter() for sample in SAMPLE_ORDER
    }

    for coarse in coarse_rows:
        sample, region = coarse["sample"], coarse["region"]
        if sample not in status_counts:
            raise RuntimeError(f"unexpected sample: {sample}")
        key = (sample, region)
        status_counts[sample]["primary_regions"] += 1
        if not truth(coarse["complete"]):
            status_counts[sample]["incomplete"] += 1
            continue

        status_counts[sample]["complete"] += 1
        vaf_region = vaf_region_by_key[key]
        structural = coarse["structural_class"]
        if structural != vaf_region["structural_class"]:
            raise RuntimeError(f"structural class mismatch: {key}")
        units = catalog_units_by_region.get(key, [])
        if len(units) != int(coarse["primary_hp_units"]):
            raise RuntimeError(f"primary HP unit count mismatch: {key}")

        if structural in STRUCTURAL_TOPO1:
            source = "structural_topo1"
            status_counts[sample][source] += 1
        elif (
            structural == "T>1|Topo>1"
            and truth(vaf_region["evaluable"])
            and int(vaf_region["n_top_shapes_joint_exact"]) == 1
        ):
            source = "vaf_resolved_topogt1"
            status_counts[sample][source] += 1
            resolved_exact_class_counts[sample][vaf_region["exact_top_class"]] += 1
            direction_counts[sample][vaf_region["top_direction_state"] or "blank"] += 1
        else:
            if vaf_region["final_class"] == "tied_first_different_topology":
                status_counts[sample]["unresolved_vaf_tie"] += 1
            else:
                status_counts[sample]["vaf_not_evaluable"] += 1
            continue

        component_rows = []
        component_categories = []
        for unit in units:
            family = unit["family"]
            if source == "structural_topo1":
                selected_shapes = unit["shape_ids"]
            elif unit["n_trees"] == 1:
                selected_shapes = unit["shape_ids"]
            else:
                vaf_unit = vaf_unit_by_key.get((sample, region, family))
                if not vaf_unit or vaf_unit["status"] != "ranked":
                    raise RuntimeError(f"VAF-resolved region lacks ranked unit: {(sample, region, family)}")
                selected_shapes = tuple(sorted(json.loads(vaf_unit["exact_top_shape_ids"])))
                if not set(selected_shapes) <= set(unit["shape_ids"]):
                    raise RuntimeError(f"selected shape outside structural candidate set: {(sample, region, family)}")
            if len(selected_shapes) != 1:
                raise RuntimeError(f"selected region component does not have one shape: {(sample, region, family)}")
            shape_id = selected_shapes[0]
            category = shape_class[shape_id]
            component_categories.append(category)
            component_rows.append((family, shape_id, category))

        category, sister, direct, split_across_hp = collapse_component_categories(component_categories)
        ordered_shape_tuple = "|".join(f"HP{family}={shape_id}" for family, shape_id, _ in component_rows)
        ordered_category_tuple = "|".join(
            f"HP{family}={component_category}" for family, _, component_category in component_rows
        )
        region_shape_id = "RSH-" + hashlib.sha1(ordered_shape_tuple.encode()).hexdigest()[:12]
        original_category = coarse["coarse_category"]
        if original_category == "no_within_hp_relation":
            original_category = "single_only"
        if source == "structural_topo1" and original_category != category:
            structural_category_mismatches += 1

        selected_rows.append(
            {
                "sample": sample,
                "biological_id": vaf_region["biological_id"],
                "region": region,
                "primary_hp_units": int(coarse["primary_hp_units"]),
                "structural_class": structural,
                "selection_source": source,
                "vaf_exact_top_class": vaf_region["exact_top_class"],
                "vaf_top_direction_state": vaf_region["top_direction_state"],
                "selected_region_shape_id": region_shape_id,
                "selected_ordered_shape_tuple": ordered_shape_tuple,
                "selected_ordered_component_categories": ordered_category_tuple,
                "final_shape_category": category,
                "final_shape_category_label": CATEGORY_LABELS[category],
                "has_sister_pattern": sister,
                "has_direct_chain_pattern": direct,
                "sister_direct_split_across_hp": split_across_hp,
            }
        )

    selected_by_sample: dict[str, list[dict]] = collections.defaultdict(list)
    for row in selected_rows:
        selected_by_sample[row["sample"]].append(row)

    summary_rows: list[dict] = []
    source_rows: list[dict] = []
    direction_rows: list[dict] = []
    checks: list[dict] = []

    def add_check(scope: str, check: str, observed: int | str, expected: int | str) -> None:
        checks.append(
            {
                "scope": scope,
                "check": check,
                "observed": observed,
                "expected": expected,
                "pass": str(observed == expected),
            }
        )

    aggregate_status = collections.Counter()
    for sample in SAMPLE_ORDER:
        aggregate_status.update(status_counts[sample])

    for sample in SAMPLE_ORDER + ["aggregate"]:
        rows = selected_rows if sample == "aggregate" else selected_by_sample[sample]
        status = aggregate_status if sample == "aggregate" else status_counts[sample]
        category_counts = collections.Counter(row["final_shape_category"] for row in rows)
        source_counts = collections.Counter(row["selection_source"] for row in rows)
        final_total = len(rows)
        complete = int(status["complete"])
        unresolved = int(status["unresolved_vaf_tie"] + status["vaf_not_evaluable"])
        row = {
            "sample": sample,
            "primary_regions": int(status["primary_regions"]),
            "complete_regions": complete,
            "incomplete_regions": int(status["incomplete"]),
            "structural_topo1_regions": int(source_counts["structural_topo1"]),
            "vaf_resolved_topogt1_regions": int(source_counts["vaf_resolved_topogt1"]),
            "final_single_shape_regions": final_total,
            "final_single_shape_pct_complete": f"{pct(final_total, complete):.2f}",
            "unresolved_regions": unresolved,
            "unresolved_pct_complete": f"{pct(unresolved, complete):.2f}",
            "unresolved_vaf_tie": int(status["unresolved_vaf_tie"]),
            "vaf_not_evaluable": int(status["vaf_not_evaluable"]),
        }
        for category in CATEGORIES:
            count = int(category_counts[category])
            row[category] = count
            row[f"{category}_pct_final_shape"] = f"{pct(count, final_total):.2f}"
            row[f"{category}_pct_complete"] = f"{pct(count, complete):.2f}"
        row["single_only_primary_hp1"] = sum(
            item["final_shape_category"] == "single_only" and int(item["primary_hp_units"]) == 1
            for item in rows
        )
        row["single_only_primary_hp2"] = sum(
            item["final_shape_category"] == "single_only" and int(item["primary_hp_units"]) == 2
            for item in rows
        )
        row["sister_and_direct_split_across_hp"] = sum(
            truth(item["sister_direct_split_across_hp"]) for item in rows
        )
        summary_rows.append(row)

        for source in ("structural_topo1", "vaf_resolved_topogt1"):
            source_subset = [item for item in rows if item["selection_source"] == source]
            source_total = len(source_subset)
            source_category_counts = collections.Counter(
                item["final_shape_category"] for item in source_subset
            )
            for category in CATEGORIES:
                count = int(source_category_counts[category])
                source_rows.append(
                    {
                        "sample": sample,
                        "selection_source": source,
                        "final_shape_category": category,
                        "regions": count,
                        "source_total": source_total,
                        "pct_within_source": f"{pct(count, source_total):.2f}",
                        "pct_of_final_single_shape": f"{pct(count, final_total):.2f}",
                    }
                )

        direction_counter = collections.Counter(
            item["vaf_top_direction_state"] or "blank"
            for item in rows
            if item["selection_source"] == "vaf_resolved_topogt1"
        )
        direction_total = sum(direction_counter.values())
        for direction_state, count in sorted(direction_counter.items()):
            direction_rows.append(
                {
                    "sample": sample,
                    "top_direction_state": direction_state,
                    "regions": count,
                    "vaf_resolved_total": direction_total,
                    "pct_vaf_resolved": f"{pct(count, direction_total):.2f}",
                }
            )

        add_check(sample, "selected_plus_unresolved_conserve_complete", final_total + unresolved, complete)
        add_check(sample, "four_categories_conserve_selected", sum(category_counts.values()), final_total)
        add_check(
            sample,
            "sources_conserve_selected",
            source_counts["structural_topo1"] + source_counts["vaf_resolved_topogt1"],
            final_total,
        )

    corrected_classes = ct_report["aggregate"]["T_and_Topology.classes"]
    aggregate_summary = summary_rows[-1]
    add_check("aggregate", "coarse_and_vaf_region_key_count", len(coarse_by_key), len(vaf_region_by_key))
    add_check("aggregate", "all_structural_categories_reproduced", structural_category_mismatches, 0)
    add_check("aggregate", "primary_matches_W_primary", len(coarse_rows), int(ct_report["aggregate"]["W_primary"]))
    add_check(
        "aggregate",
        "complete_matches_corrected_report",
        int(aggregate_summary["complete_regions"]),
        sum(int(corrected_classes[key]) for key in STRUCTURAL_TOPO1 | {"T>1|Topo>1"}),
    )
    add_check(
        "aggregate",
        "structural_topo1_matches_corrected_report",
        int(aggregate_summary["structural_topo1_regions"]),
        int(corrected_classes["T=1|Topo=1"]) + int(corrected_classes["T>1|Topo=1"]),
    )
    add_check("aggregate", "vaf_resolved_matches_frozen_endpoint", int(aggregate_summary["vaf_resolved_topogt1_regions"]), 15063)
    add_check("aggregate", "final_single_shape_matches_frozen_endpoint", int(aggregate_summary["final_single_shape_regions"]), 37039)
    add_check("aggregate", "unresolved_vaf_tie_matches_frozen_endpoint", int(aggregate_summary["unresolved_vaf_tie"]), 2205)
    add_check("aggregate", "vaf_not_evaluable_matches_frozen_endpoint", int(aggregate_summary["vaf_not_evaluable"]), 641)
    add_check(
        "aggregate",
        "selected_component_shape_tuples_complete",
        sum(len(row["selected_ordered_shape_tuple"].split("|")) for row in selected_rows),
        sum(int(row["primary_hp_units"]) for row in selected_rows),
    )
    add_check(
        "aggregate",
        "vaf_resolved_exact_top_classes_conserve",
        sum(sum(counter.values()) for counter in resolved_exact_class_counts.values()),
        int(aggregate_summary["vaf_resolved_topogt1_regions"]),
    )
    add_check(
        "aggregate",
        "vaf_resolved_direction_states_conserve",
        sum(sum(counter.values()) for counter in direction_counts.values()),
        int(aggregate_summary["vaf_resolved_topogt1_regions"]),
    )
    add_check("aggregate", "all_checks_pass", sum(row["pass"] == "True" for row in checks), len(checks))
    if not all(row["pass"] == "True" for row in checks):
        failed = [row for row in checks if row["pass"] != "True"]
        raise SystemExit(f"validation failed: {failed[:5]}")

    output_dir = args.output_dir
    region_output = output_dir / "20260712_vaf_final_single_shape_regions.tsv"
    summary_output = output_dir / "20260712_vaf_final_single_shape_four_class_summary.tsv"
    source_output = output_dir / "20260712_vaf_final_single_shape_four_class_by_source.tsv"
    direction_output = output_dir / "20260712_vaf_final_single_shape_direction_states.tsv"
    checks_output = output_dir / "20260712_vaf_final_single_shape_checks.tsv"
    json_output = output_dir / "20260712_vaf_final_single_shape_four_class_census.json"

    write_tsv(region_output, selected_rows, list(selected_rows[0]))
    write_tsv(summary_output, summary_rows, list(summary_rows[0]))
    write_tsv(source_output, source_rows, list(source_rows[0]))
    write_tsv(direction_output, direction_rows, list(direction_rows[0]))
    write_tsv(checks_output, checks, ["scope", "check", "observed", "expected", "pass"])

    sources = {
        name: {"path": str(path.resolve()), "sha256": sha256(path)}
        for name, path in {
            "exact_catalog": args.exact_catalog,
            "ct_report": args.ct_report,
            "coarse_regions": args.coarse_regions,
            "vaf_units": args.vaf_units,
            "vaf_regions": args.vaf_regions,
        }.items()
    }
    canonical = {
        "schema_version": "1.0",
        "generated_at": datetime.now(ZoneInfo("Asia/Taipei")).isoformat(timespec="seconds"),
        "status": "PASS",
        "task_type": "B_comprehensive_validation",
        "scope": "chr1-22; 7 dataset rows; 6 biological samples; historical layered-v2 engineering snapshot",
        "claim_ceiling": "VAF-selected rooted-unlabeled mutation-state graph patterns; not confirmed clones/subclones",
        "definitions": {
            "single_only": "No primary HP component has max_outdegree>=2 or max_depth>=2",
            "sister_only": "At least one primary HP component has max_outdegree>=2; none has max_depth>=2",
            "direct_only": "At least one primary HP component has max_depth>=2; none has max_outdegree>=2",
            "sister_and_direct": "Region-level OR has both features; they may occur in different HP components",
            "selection": "All structural Topo=1 regions plus original Topo>1 regions whose exact VAF first-rank set maps to one ordered HP shape tuple",
        },
        "sources": sources,
        "samples": summary_rows[:-1],
        "aggregate": summary_rows[-1],
        "vaf_resolved_exact_top_classes": {
            sample: dict(resolved_exact_class_counts[sample]) for sample in SAMPLE_ORDER
        },
        "vaf_resolved_direction_states": {
            sample: dict(direction_counts[sample]) for sample in SAMPLE_ORDER
        },
        "validation": {
            "checks": len(checks),
            "passed": sum(row["pass"] == "True" for row in checks),
        },
        "outputs": {
            "regions": str(region_output.resolve()),
            "summary": str(summary_output.resolve()),
            "by_source": str(source_output.resolve()),
            "direction_states": str(direction_output.resolve()),
            "checks": str(checks_output.resolve()),
        },
    }
    write_json(json_output, canonical)

    print(f"INPUT exact catalog -> {args.exact_catalog.resolve()}")
    print(f"INPUT VAF region/unit -> {args.vaf_regions.resolve()} | {args.vaf_units.resolve()}")
    print(f"OUTPUT summary -> {summary_output.resolve()}")
    print(f"OUTPUT regions -> {region_output.resolve()}")
    print(f"OUTPUT checks -> {checks_output.resolve()}")
    print(
        "RESULT pooled -> "
        + ", ".join(
            f"{category}={aggregate_summary[category]}"
            for category in CATEGORIES
        )
        + f", final={aggregate_summary['final_single_shape_regions']}"
    )
    print(f"STATUS -> PASS ({len(checks)}/{len(checks)} checks)")


if __name__ == "__main__":
    main()
