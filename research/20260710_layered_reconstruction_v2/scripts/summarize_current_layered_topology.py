#!/usr/bin/env python3
"""Summarize current canonical and sensitivity regional topology outputs."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
from pathlib import Path


EXPECTED_SAMPLES = {
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
}

TOPOLOGY_FIELDS = (
    "exact_and_topology_unique",
    "topology_unique_exact_multiple",
    "topology_multiple_exact_multiple",
    "incomplete",
    "impossible_exact_unique_topology_multiple",
)


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str):
    if not condition:
        raise ValueError(message)


def load_topology_function(module_path: Path):
    spec = importlib.util.spec_from_file_location("layered_observation_report", module_path)
    require(spec is not None and spec.loader is not None, f"cannot import topology module: {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.region_candidate_topology


def sum_nested(rows, key, fields):
    return {field: sum(int(row[key][field]) for row in rows) for field in fields}


def summarize_sample(run_root: Path, verification_row, topology_function):
    sample = verification_row["sample"]
    sample_root = run_root / "samples" / sample
    layered_path = sample_root / f"layered_reconstruction_{sample}.json"
    region_path = sample_root / f"layered_region_view_{sample}.json"
    require(layered_path.is_file(), f"missing layered reconstruction: {layered_path}")
    require(region_path.is_file(), f"missing region view: {region_path}")

    layered = load_json(layered_path)
    region_view = load_json(region_path)
    require(layered.get("sample") == sample, f"layered sample mismatch: {sample}")
    require(region_view.get("sample") == sample, f"region sample mismatch: {sample}")

    detail = layered.get("detail", [])
    regions = region_view.get("regions", [])
    detail_by_key = {
        (unit.get("region"), str(unit.get("family"))): unit
        for unit in detail
    }
    topology = topology_function(regions, detail_by_key)
    census = region_view["census"]
    l1 = census["L1"]
    funnel = census["funnel"]
    tags = census["read_tag_census"]

    metrics = verification_row.get("metrics", {})
    require(verification_row.get("pass") is True, f"sample verification failed: {sample}")
    require(len(regions) == int(metrics.get("n_regions", -1)), f"region count mismatch: {sample}")
    require(len(regions) == int(metrics.get("n_groups", -1)), f"group/region mismatch: {sample}")
    require(len(detail) == int(metrics.get("n_detail_units", -1)), f"detail count mismatch: {sample}")
    require(l1.get("all_V1V7_pass") is True, f"V1-V7 aggregate failed: {sample}")
    require(l1.get("all_eligible_V1V7_pass") is True, f"eligible V1-V7 failed: {sample}")
    require(int(l1.get("n_verification_fail", -1)) == 0, f"verification failures present: {sample}")
    require(int(l1.get("n_eligible_skipped_V4V5", -1)) == 0, f"eligible V4/V5 skipped: {sample}")
    require(tags.get("check_exact_sidecar_join") is True, f"read-tag exact join failed: {sample}")
    require(int(tags.get("alignment_group_exposures", -1)) == int(tags.get("sidecar_exact_matches", -2)),
            f"read-tag exposure/match mismatch: {sample}")
    for field in ("sidecar_missing", "sidecar_conflicts", "sidecar_extra", "sidecar_malformed",
                  "alignment_identity_allele_conflicts"):
        require(int(tags.get(field, -1)) == 0, f"{field} is non-zero: {sample}")

    w_tree = len(regions)
    w_primary = sum(1 for region in regions if int(region.get("n_primary_lineages", 0)) > 0)
    no_primary = w_tree - w_primary
    topology_total = sum(int(topology["topology_classes"].get(field, 0)) for field in TOPOLOGY_FIELDS)
    require(topology_total == w_primary, f"topology classes do not conserve primary regions: {sample}")
    require(topology["complete_regions"] + topology["incomplete_regions"] == w_primary,
            f"complete/incomplete does not conserve primary regions: {sample}")
    require(int(topology["topology_classes"]["impossible_exact_unique_topology_multiple"]) == 0,
            f"impossible topology state observed: {sample}")
    require(not topology["invalid_states"], f"invalid topology states observed: {sample}")
    require(int(funnel["n_groups_retained"]) == w_tree, f"funnel/region mismatch: {sample}")

    primary_units = [unit for unit in detail if unit.get("is_primary_lineage")]
    full_and_partial = sum(
        1 for unit in primary_units
        if int(unit.get("n_full_pops", 0)) > 0 and int(unit.get("n_partial", 0)) > 0
    )
    partial_only = sum(
        1 for unit in primary_units
        if int(unit.get("n_full_pops", 0)) == 0 and int(unit.get("n_partial", 0)) > 0
    )
    full_only = sum(
        1 for unit in primary_units
        if int(unit.get("n_full_pops", 0)) > 0 and int(unit.get("n_partial", 0)) == 0
    )

    return {
        "sample": sample,
        "pass": True,
        "tree_backbone_source": census.get("U1_backbone_source"),
        "tree_input_records": int(census["U1_sSNV_somatic_total"]),
        "autosomal_biallelic_sSNV": int(census["U1_sSNV_scope_chr1_22"]),
        "retained_sSNV": int(funnel["n_sSNV_retained"]),
        "W_tree": w_tree,
        "W_primary": w_primary,
        "no_primary_lineage": no_primary,
        "primary_units": len(primary_units),
        "primary_units_full_and_partial": full_and_partial,
        "primary_units_partial_only": partial_only,
        "primary_units_full_only": full_only,
        "complete_regions": int(topology["complete_regions"]),
        "incomplete_regions": int(topology["incomplete_regions"]),
        "topology_classes": topology["topology_classes"],
        "hidden_classes": topology["hidden_classes"],
        "hp_h3": topology["hp_h3"],
        "C_bins": topology["C_bins"]["all"],
        "joint_exact_combinations_total": int(topology["joint_exact_combinations_total"]),
        "joint_topology_combinations_total": int(topology["joint_topology_combinations_total"]),
        "max_C_region": int(topology["max_C_region"]),
        "max_Topo_region": int(topology["max_Topo_region"]),
        "legacy_has_capped_label_count": int(topology["legacy_has_capped_label_count"]),
        "capped_label_discordance_count": len(topology["capped_label_discordance"]),
        "read_tag_exposures": int(tags["alignment_group_exposures"]),
        "read_tag_exact_matches": int(tags["sidecar_exact_matches"]),
        "mixed_PS_regions": int(tags.get("phase_set_region_counts", {}).get("multiple", 0)),
        "all_invariants_pass": True,
        "paths": {
            "layered_reconstruction": str(layered_path),
            "layered_region_view": str(region_path),
        },
        "sha256": {
            "layered_reconstruction": sha256(layered_path),
            "layered_region_view": sha256(region_path),
        },
    }


def summarize_root(label: str, run_root: Path, topology_function):
    success_path = run_root / "_SUCCESS"
    verification_path = run_root / "verification_summary.json"
    require(success_path.is_file(), f"missing _SUCCESS: {run_root}")
    require(verification_path.is_file(), f"missing verification summary: {run_root}")
    verification = load_json(verification_path)
    require(verification.get("all_pass") is True, f"run verification failed: {run_root}")
    require(int(verification.get("n_pass", -1)) == 7, f"run is not 7/7 PASS: {run_root}")
    require(int(verification.get("n_fail", -1)) == 0, f"run has failed datasets: {run_root}")
    sample_names = [row.get("sample") for row in verification.get("samples", [])]
    require(len(sample_names) == 7 and set(sample_names) == EXPECTED_SAMPLES,
            f"unexpected sample set: {sample_names}")

    samples = [summarize_sample(run_root, row, topology_function) for row in verification["samples"]]
    topology_classes = sum_nested(samples, "topology_classes", TOPOLOGY_FIELDS)
    hidden_fields = ("hidden_zero", "hidden_positive", "incomplete")
    hp_h3_fields = tuple(samples[0]["hp_h3"].keys())
    c_fields = tuple(samples[0]["C_bins"].keys())
    aggregate = {
        "dataset_count": 7,
        "tree_input_records": sum(row["tree_input_records"] for row in samples),
        "autosomal_biallelic_sSNV": sum(row["autosomal_biallelic_sSNV"] for row in samples),
        "retained_sSNV": sum(row["retained_sSNV"] for row in samples),
        "W_tree": sum(row["W_tree"] for row in samples),
        "W_primary": sum(row["W_primary"] for row in samples),
        "no_primary_lineage": sum(row["no_primary_lineage"] for row in samples),
        "primary_units": sum(row["primary_units"] for row in samples),
        "primary_units_full_and_partial": sum(row["primary_units_full_and_partial"] for row in samples),
        "primary_units_partial_only": sum(row["primary_units_partial_only"] for row in samples),
        "primary_units_full_only": sum(row["primary_units_full_only"] for row in samples),
        "complete_regions": sum(row["complete_regions"] for row in samples),
        "incomplete_regions": sum(row["incomplete_regions"] for row in samples),
        "topology_classes": topology_classes,
        "hidden_classes": sum_nested(samples, "hidden_classes", hidden_fields),
        "hp_h3": sum_nested(samples, "hp_h3", hp_h3_fields),
        "C_bins": sum_nested(samples, "C_bins", c_fields),
        "joint_exact_combinations_total": sum(row["joint_exact_combinations_total"] for row in samples),
        "joint_topology_combinations_total": sum(row["joint_topology_combinations_total"] for row in samples),
        "max_C_region": max(row["max_C_region"] for row in samples),
        "max_Topo_region": max(row["max_Topo_region"] for row in samples),
        "legacy_has_capped_label_count": sum(row["legacy_has_capped_label_count"] for row in samples),
        "capped_label_discordance_count": sum(row["capped_label_discordance_count"] for row in samples),
        "read_tag_exposures": sum(row["read_tag_exposures"] for row in samples),
        "read_tag_exact_matches": sum(row["read_tag_exact_matches"] for row in samples),
        "mixed_PS_regions": sum(row["mixed_PS_regions"] for row in samples),
        "all_invariants_pass": True,
    }
    require(aggregate["W_primary"] == aggregate["complete_regions"] + aggregate["incomplete_regions"],
            f"aggregate complete/incomplete mismatch: {label}")
    require(aggregate["W_primary"] == sum(topology_classes.values()),
            f"aggregate topology-class mismatch: {label}")
    require(aggregate["read_tag_exposures"] == aggregate["read_tag_exact_matches"],
            f"aggregate read-tag mismatch: {label}")
    return {
        "label": label,
        "run_root": str(run_root),
        "success_sha256": sha256(success_path),
        "verification_path": str(verification_path),
        "verification_sha256": sha256(verification_path),
        "aggregate": aggregate,
        "samples": samples,
    }


def write_tsv(path: Path, summaries):
    fields = [
        "backbone", "sample", "tree_input_records", "autosomal_biallelic_sSNV", "retained_sSNV",
        "W_tree", "W_primary", "no_primary_lineage", "primary_units", "complete_regions",
        "incomplete_regions", "C1_T1", "Cgt1_T1", "Cgt1_Tgt1", "impossible", "hidden_zero",
        "hidden_positive", "read_tag_exposures", "read_tag_exact_matches", "mixed_PS_regions",
        "all_invariants_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for summary in summaries:
            for row in summary["samples"]:
                topology = row["topology_classes"]
                hidden = row["hidden_classes"]
                writer.writerow({
                    "backbone": summary["label"],
                    "sample": row["sample"],
                    "tree_input_records": row["tree_input_records"],
                    "autosomal_biallelic_sSNV": row["autosomal_biallelic_sSNV"],
                    "retained_sSNV": row["retained_sSNV"],
                    "W_tree": row["W_tree"],
                    "W_primary": row["W_primary"],
                    "no_primary_lineage": row["no_primary_lineage"],
                    "primary_units": row["primary_units"],
                    "complete_regions": row["complete_regions"],
                    "incomplete_regions": row["incomplete_regions"],
                    "C1_T1": topology["exact_and_topology_unique"],
                    "Cgt1_T1": topology["topology_unique_exact_multiple"],
                    "Cgt1_Tgt1": topology["topology_multiple_exact_multiple"],
                    "impossible": topology["impossible_exact_unique_topology_multiple"],
                    "hidden_zero": hidden["hidden_zero"],
                    "hidden_positive": hidden["hidden_positive"],
                    "read_tag_exposures": row["read_tag_exposures"],
                    "read_tag_exact_matches": row["read_tag_exact_matches"],
                    "mixed_PS_regions": row["mixed_PS_regions"],
                    "all_invariants_pass": str(row["all_invariants_pass"]).lower(),
                })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-root", required=True, type=Path)
    parser.add_argument("--sensitivity-root", required=True, type=Path)
    parser.add_argument("--backbone-comparison", required=True, type=Path)
    parser.add_argument("--topology-module", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    args = parser.parse_args()

    require(args.backbone_comparison.is_file(), "missing backbone comparison")
    require(args.topology_module.is_file(), "missing topology module")
    topology_function = load_topology_function(args.topology_module)
    canonical = summarize_root("longphase_s_recalibrated_FILTER_PASS", args.canonical_root, topology_function)
    sensitivity = summarize_root("clairs_FILTER_PASS_sensitivity", args.sensitivity_root, topology_function)
    comparison = load_json(args.backbone_comparison)
    require(comparison.get("aggregate", {}).get("verdict") == "backbone_sensitive",
            "unexpected backbone comparison verdict")
    require(len(comparison.get("comparisons", [])) == 7, "backbone comparison is not 7 datasets")

    result = {
        "schema_name": "intersubmod.current_layered_topology_summary",
        "schema_version": "1.0.0",
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "scope": "7 datasets / 6 biological samples / chr1-22",
        "claim_scope": "regional mutation-state candidate trees; not confirmed cell clones",
        "definitions": {
            "C_region": "product of n_trees across analysis-complete primary HP1/HP2 units",
            "Topo_region": "product of n_distinct_shapes_exact across those units",
            "complete": "all primary units are non-capped, full-pass, and candidate-complete",
            "PS": "phase-block QC context only; not a topology edge or lineage label",
        },
        "all_pass": True,
        "canonical": canonical,
        "sensitivity": sensitivity,
        "backbone_comparison": {
            "path": str(args.backbone_comparison),
            "sha256": sha256(args.backbone_comparison),
            "aggregate": comparison["aggregate"],
        },
        "code_provenance": {
            "summarizer": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
            "topology_module": {"path": str(args.topology_module.resolve()), "sha256": sha256(args.topology_module)},
        },
    }
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_tsv(args.output_tsv, (canonical, sensitivity))
    print(f"SUMMARY_JSON={args.output_json}")
    print(f"SUMMARY_TSV={args.output_tsv}")
    print("ALL_PASS=true")
    print(f"CANONICAL_W_TREE={canonical['aggregate']['W_tree']}")
    print(f"CANONICAL_W_PRIMARY={canonical['aggregate']['W_primary']}")
    print(f"CANONICAL_COMPLETE={canonical['aggregate']['complete_regions']}")
    print(f"CANONICAL_INCOMPLETE={canonical['aggregate']['incomplete_regions']}")
    print(f"CANONICAL_TOPOLOGY={json.dumps(canonical['aggregate']['topology_classes'], sort_keys=True)}")


if __name__ == "__main__":
    main()
