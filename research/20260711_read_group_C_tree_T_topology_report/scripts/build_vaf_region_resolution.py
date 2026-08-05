#!/usr/bin/env python3
"""Derive region-level VAF-prioritized candidate convergence for seven datasets.

This analysis does not alter the structural candidate set and does not claim a
confirmed biological topology.  A complete structurally ambiguous region is
called VAF-prioritized only when every ambiguous primary HP unit in that region
has a top candidate with relative weight >= the configured threshold.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
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
STRUCTURAL_CLASSES = ("T=1|Topo=1", "T>1|Topo=1", "T>1|Topo>1", "incomplete")
AMBIGUOUS_CLASSES = ("T>1|Topo=1", "T>1|Topo>1")
VAF_REGION_STATUSES = (
    "vaf_top_consistent",
    "vaf_top_direction_inconsistent",
    "below_threshold",
    "missing_or_mismatch",
    "not_evaluable_recurrence",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dump_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def structural_class(units: list[dict]) -> str:
    if not all(unit.get("analysis_candidate_set_complete") is True for unit in units):
        return "incomplete"
    t_joint = math.prod(int(unit["n_trees"]) for unit in units)
    topo_joint = math.prod(int(unit["n_distinct_shapes_exact"]) for unit in units)
    if t_joint == 1 and topo_joint == 1:
        return "T=1|Topo=1"
    if t_joint > 1 and topo_joint == 1:
        return "T>1|Topo=1"
    if t_joint > 1 and topo_joint > 1:
        return "T>1|Topo>1"
    raise RuntimeError(f"impossible structural state T={t_joint} Topo={topo_joint}")


def analyze_sample(meta: dict, sample_dir: Path, read_af_ref: dict, corrected_ref: dict, ra, solver, params: dict):
    sample = meta["sample"]
    layered_path = sample_dir / f"layered_reconstruction_{sample}.json"
    layered = json.loads(layered_path.read_text(encoding="utf-8"))
    groups = ra.load_groups(sample_dir)
    primary_units = [unit for unit in layered["detail"] if unit.get("is_primary_lineage")]
    primary_by_region: dict[str, list[dict]] = defaultdict(list)
    for unit in primary_units:
        primary_by_region[unit["region"]].append(unit)
    for units in primary_by_region.values():
        units.sort(key=lambda unit: unit["family"])

    ambiguous_units = [
        unit
        for unit in primary_units
        if not unit.get("capped") and unit.get("L1_base_class") == "ambiguous"
    ]
    unit_results = {}
    for unit in ambiguous_units:
        key = (unit["region"], unit["family"])
        if key in unit_results:
            raise RuntimeError(f"duplicate primary unit key: {sample} {key}")
        group = groups.get(unit["region"])
        if not group:
            unit_results[key] = {"status": "missing_group", "cn": ra.cn_class(unit.get("cn"))}
            continue
        family = unit["family"]
        full = (group.get("populations_by_hp", {}) or {}).get(family, {}) or {}
        partial = list(
            ((group.get("subread_groups_by_hp", {}) or {}).get(family, {}) or {}).keys()
        )
        k = len(group.get("positions", []))
        result = solver.enumerate_min_trees(full, partial, k, tree_cap=0)
        if (
            result.get("capped")
            or not result.get("trees_complete")
            or int(result["n_trees"]) != int(unit["n_trees"])
        ):
            unit_results[key] = {
                "status": "candidate_mismatch",
                "cn": ra.cn_class(unit.get("cn")),
            }
            continue
        read_af = ra.read_af_from_colcov(
            ((group.get("col_coverage_by_hp", {}) or {}).get(family, {})),
            group.get("positions", []),
            k,
        )
        if read_af is None:
            unit_results[key] = {"status": "missing_vaf", "cn": ra.cn_class(unit.get("cn"))}
            continue
        scored = [ra.ordering_score(tree["edges"], read_af) for tree in result["trees"]]
        weights = ra.posterior([item[0] for item in scored], params["temperature"])
        top_weight = max(weights)
        winner_index = weights.index(top_weight)
        direction_consistent = all(
            delta >= -params["violation_margin"]
            for delta in scored[winner_index][1]
        )
        if top_weight < params["posterior_threshold"]:
            status = "below_threshold"
        elif direction_consistent:
            status = "vaf_top_consistent"
        else:
            status = "vaf_top_direction_inconsistent"
        unit_results[key] = {
            "status": status,
            "cn": ra.cn_class(unit.get("cn")),
            "top_weight": top_weight,
            "direction_consistent": direction_consistent,
            "n_trees": int(result["n_trees"]),
        }

    unit_counts = Counter(result["status"] for result in unit_results.values())
    prepared = sum(
        unit_counts[key]
        for key in (
            "vaf_top_consistent",
            "vaf_top_direction_inconsistent",
            "below_threshold",
        )
    )
    reached = unit_counts["vaf_top_consistent"] + unit_counts["vaf_top_direction_inconsistent"]

    structural_counts = Counter()
    ambiguous_region_status = Counter()
    cross = defaultdict(Counter)
    complete_recurrence_tn_units = 0
    all_tn_units_have_supported_class = True
    for region, units in sorted(primary_by_region.items()):
        category = structural_class(units)
        structural_counts[category] += 1
        if category not in AMBIGUOUS_CLASSES:
            continue
        t_multiple_units = [unit for unit in units if int(unit["n_trees"]) > 1]
        recurrence_units = [
            unit
            for unit in t_multiple_units
            if unit.get("L1_base_class") == "recurrence_required"
        ]
        complete_recurrence_tn_units += len(recurrence_units)
        unsupported = [
            unit
            for unit in t_multiple_units
            if unit.get("L1_base_class") not in {"ambiguous", "recurrence_required"}
        ]
        if unsupported:
            all_tn_units_have_supported_class = False
            raise RuntimeError(
                f"T>1 unit has unsupported L1 class: {sample} {region} "
                f"{[unit.get('L1_base_class') for unit in unsupported]}"
            )
        if recurrence_units:
            region_status = "not_evaluable_recurrence"
            ambiguous_region_status[region_status] += 1
            cross[category][region_status] += 1
            continue
        ambiguous_keys = [(region, unit["family"]) for unit in t_multiple_units]
        if not ambiguous_keys:
            raise RuntimeError(f"T>1 region lacks T>1 units: {sample} {region}")
        statuses = [unit_results[key]["status"] for key in ambiguous_keys]
        if any(status in {"missing_group", "missing_vaf", "candidate_mismatch"} for status in statuses):
            region_status = "missing_or_mismatch"
        elif any(status == "below_threshold" for status in statuses):
            region_status = "below_threshold"
        elif any(status == "vaf_top_direction_inconsistent" for status in statuses):
            region_status = "vaf_top_direction_inconsistent"
        elif all(status == "vaf_top_consistent" for status in statuses):
            region_status = "vaf_top_consistent"
        else:
            raise RuntimeError(f"unclassified VAF statuses: {sample} {region} {statuses}")
        ambiguous_region_status[region_status] += 1
        cross[category][region_status] += 1

    expected_classes = corrected_ref["T_and_Topology"]["classes"]
    expected_vaf = read_af_ref
    checks = {
        "structural_classes_match_corrected_report": all(
            structural_counts[key] == int(expected_classes.get(key, 0))
            for key in STRUCTURAL_CLASSES
        ),
        "ambiguous_unit_count_matches_read_af_report": len(ambiguous_units)
        == int(expected_vaf["n_ambiguous_primary_units"]),
        "prepared_unit_count_matches_read_af_report": prepared
        == int(expected_vaf["n_units_analyzed_all_candidates"]),
        "missing_unit_count_matches_read_af_report": (
            unit_counts["missing_group"] + unit_counts["missing_vaf"]
        )
        == int(expected_vaf["n_missing_read_af"]),
        "candidate_mismatch_matches_read_af_report": unit_counts["candidate_mismatch"]
        == int(expected_vaf["n_candidate_mismatch_or_incomplete"]),
        "top_threshold_count_matches_read_af_report": reached
        == int(expected_vaf["default"]["reached"]),
        "winner_consistency_matches_read_af_report": unit_counts["vaf_top_consistent"]
        == int(expected_vaf["default"]["winner_order_consistent"]),
        "ambiguous_region_status_conservation": sum(ambiguous_region_status.values())
        == structural_counts["T>1|Topo=1"] + structural_counts["T>1|Topo>1"],
        "cross_conservation_Tn_Topo1": sum(cross["T>1|Topo=1"].values())
        == structural_counts["T>1|Topo=1"],
        "cross_conservation_Tn_Topon": sum(cross["T>1|Topo>1"].values())
        == structural_counts["T>1|Topo>1"],
        "all_Tn_units_have_supported_L1_class": all_tn_units_have_supported_class,
    }

    return {
        "sample": sample,
        "biological_id": meta["biological_id"],
        "source_layered": str(layered_path),
        "structural_region_classes": {
            key: structural_counts[key] for key in STRUCTURAL_CLASSES
        },
        "ambiguous_complete_regions": sum(ambiguous_region_status.values()),
        "region_vaf_status": {
            key: ambiguous_region_status[key] for key in VAF_REGION_STATUSES
        },
        "by_original_structural_class": {
            category: {key: cross[category][key] for key in VAF_REGION_STATUSES}
            for category in AMBIGUOUS_CLASSES
        },
        "unit_vaf_status": {
            "ambiguous_primary_HP_units": len(ambiguous_units),
            "prepared": prepared,
            "missing_vaf": unit_counts["missing_group"] + unit_counts["missing_vaf"],
            "candidate_mismatch": unit_counts["candidate_mismatch"],
            "top_weight_ge_threshold": reached,
            "top_and_direction_consistent": unit_counts["vaf_top_consistent"],
            "top_direction_inconsistent": unit_counts["vaf_top_direction_inconsistent"],
            "below_threshold": unit_counts["below_threshold"],
            "complete_region_recurrence_required_Tn_units": complete_recurrence_tn_units,
        },
        "checks": checks,
    }


def sum_nested(samples: list[dict], *path: str) -> int:
    total = 0
    for sample in samples:
        value = sample
        for key in path:
            value = value[key]
        total += int(value)
    return total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--read-af-summary", required=True, type=Path)
    parser.add_argument("--corrected-report", required=True, type=Path)
    parser.add_argument("--method-script-dir", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-summary-tsv", required=True, type=Path)
    parser.add_argument("--output-checks-tsv", required=True, type=Path)
    parser.add_argument("--temperature", type=float, default=0.05)
    parser.add_argument("--posterior-threshold", type=float, default=0.60)
    parser.add_argument("--violation-margin", type=float, default=0.05)
    args = parser.parse_args()

    sys.path.insert(0, str(args.method_script_dir.resolve()))
    import read_af_tree_ordering_multisample as ra  # noqa: PLC0415
    import tree_enumeration_solver as solver  # noqa: PLC0415

    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    read_af = json.loads(args.read_af_summary.read_text(encoding="utf-8"))
    corrected = json.loads(args.corrected_report.read_text(encoding="utf-8"))
    meta_by_name = {item["sample"]: item for item in manifest["samples"]}
    read_af_by_name = {item["sample"]: item for item in read_af["samples"]}
    corrected_by_name = {item["sample"]: item for item in corrected["samples"]}
    if set(meta_by_name) != set(SAMPLE_ORDER) or set(read_af_by_name) != set(SAMPLE_ORDER) or set(corrected_by_name) != set(SAMPLE_ORDER):
        raise SystemExit("seven-row sample scope mismatch")

    params = {
        "temperature": args.temperature,
        "posterior_threshold": args.posterior_threshold,
        "violation_margin": args.violation_margin,
    }
    samples = [
        analyze_sample(
            meta_by_name[name],
            args.run_root / "samples" / name,
            read_af_by_name[name],
            corrected_by_name[name],
            ra,
            solver,
            params,
        )
        for name in SAMPLE_ORDER
    ]

    aggregate_structural = {
        key: sum_nested(samples, "structural_region_classes", key)
        for key in STRUCTURAL_CLASSES
    }
    aggregate_region_vaf = {
        key: sum_nested(samples, "region_vaf_status", key)
        for key in VAF_REGION_STATUSES
    }
    aggregate_cross = {
        category: {
            key: sum_nested(samples, "by_original_structural_class", category, key)
            for key in VAF_REGION_STATUSES
        }
        for category in AMBIGUOUS_CLASSES
    }
    aggregate_unit_vaf = {
        key: sum_nested(samples, "unit_vaf_status", key)
        for key in samples[0]["unit_vaf_status"]
    }
    all_checks = [value for sample in samples for value in sample["checks"].values()]
    aggregate_checks = {
        "all_sample_checks_pass": all(all_checks),
        "structural_complete_conservation": sum(
            aggregate_structural[key]
            for key in ("T=1|Topo=1", "T>1|Topo=1", "T>1|Topo>1")
        )
        == int(corrected["aggregate"]["W_primary"])
        - int(corrected["aggregate"]["T_and_Topology.classes"]["incomplete"]),
        "ambiguous_region_conservation": sum(aggregate_region_vaf.values())
        == aggregate_structural["T>1|Topo=1"] + aggregate_structural["T>1|Topo>1"],
        "aggregate_unit_vaf_matches_report": (
            aggregate_unit_vaf["ambiguous_primary_HP_units"]
            == int(corrected["aggregate"]["VAF_ambiguous_units"])
            and aggregate_unit_vaf["prepared"] == int(corrected["aggregate"]["VAF_prepared"])
            and aggregate_unit_vaf["missing_vaf"] == int(corrected["aggregate"]["VAF_missing"])
            and aggregate_unit_vaf["top_weight_ge_threshold"]
            == int(corrected["aggregate"]["VAF_top_weight_ge_0_60"])
            and aggregate_unit_vaf["top_and_direction_consistent"]
            == int(corrected["aggregate"]["VAF_winner_order_consistent"])
        ),
    }
    validation_pass = all(aggregate_checks.values())
    generated_at = datetime.now(ZoneInfo("Asia/Taipei")).isoformat(timespec="seconds")
    output = {
        "schema_version": "1.0",
        "generated_at": generated_at,
        "status": "PASS" if validation_pass else "FAIL",
        "scope": "chr1-22; 7 dataset rows / 6 biological samples; historical layered-v2 engineering snapshot",
        "claim_ceiling": "VAF-prioritized candidate selection within the existing exact-T set; not corrected or confirmed biological topology",
        "unit_of_analysis": {
            "structural": "complete primary region",
            "ranking_input": "ambiguous primary HP unit",
            "region_resolution": "structurally ambiguous complete primary region",
        },
        "rule": (
            "A structurally ambiguous complete region is VAF-prioritized only when every ambiguous "
            "primary HP unit has a top exact-T relative weight >= posterior_threshold. Direction "
            "consistency is reported separately and uses the same read-VAF evidence."
        ),
        "parameters": params,
        "aggregate": {
            "structural_region_classes": aggregate_structural,
            "ambiguous_complete_regions": sum(aggregate_region_vaf.values()),
            "region_vaf_status": aggregate_region_vaf,
            "by_original_structural_class": aggregate_cross,
            "unit_vaf_status": aggregate_unit_vaf,
        },
        "samples": samples,
        "validation": {
            "status": "PASS" if validation_pass else "FAIL",
            "aggregate_checks": aggregate_checks,
            "failures": [key for key, value in aggregate_checks.items() if not value]
            + [
                f"{sample['sample']}:{key}"
                for sample in samples
                for key, value in sample["checks"].items()
                if not value
            ],
        },
        "sources": {
            "run_root": str(args.run_root.resolve()),
            "input_manifest": {"path": str(args.input_manifest.resolve()), "sha256": sha256(args.input_manifest)},
            "read_af_summary": {"path": str(args.read_af_summary.resolve()), "sha256": sha256(args.read_af_summary)},
            "corrected_report": {"path": str(args.corrected_report.resolve()), "sha256": sha256(args.corrected_report)},
            "method_script": {
                "path": str((args.method_script_dir / "read_af_tree_ordering_multisample.py").resolve()),
                "sha256": sha256(args.method_script_dir / "read_af_tree_ordering_multisample.py"),
            },
        },
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    dump_json(args.output_json, output)
    with args.output_summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "ambiguous_complete_regions",
                "vaf_top_consistent",
                "vaf_top_direction_inconsistent",
                "below_threshold",
                "missing_or_mismatch",
                "not_evaluable_recurrence",
            ]
        )
        for sample in samples:
            status = sample["region_vaf_status"]
            writer.writerow(
                [
                    sample["sample"],
                    sample["ambiguous_complete_regions"],
                    status["vaf_top_consistent"],
                    status["vaf_top_direction_inconsistent"],
                    status["below_threshold"],
                    status["missing_or_mismatch"],
                    status["not_evaluable_recurrence"],
                ]
            )
    with args.output_checks_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["scope", "check", "pass"])
        for sample in samples:
            for key, value in sample["checks"].items():
                writer.writerow([sample["sample"], key, str(bool(value))])
        for key, value in aggregate_checks.items():
            writer.writerow(["aggregate", key, str(bool(value))])

    print(f"INPUT RUN ROOT -> {args.run_root.resolve()}")
    print(f"INPUT READ-AF SUMMARY -> {args.read_af_summary.resolve()}")
    print(f"OUTPUT JSON -> {args.output_json.resolve()}")
    print(f"OUTPUT SUMMARY TSV -> {args.output_summary_tsv.resolve()}")
    print(f"OUTPUT CHECKS TSV -> {args.output_checks_tsv.resolve()}")
    print(f"STATUS -> {'PASS' if validation_pass else 'FAIL'}")
    if not validation_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
