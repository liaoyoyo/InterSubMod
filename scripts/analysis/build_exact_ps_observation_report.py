#!/usr/bin/env python3
"""Build one authority-verified, offline Exact-PS cohort observation report."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable


SCHEMA_NAME = "intersubmod.exact_ps_observation_report.data"
SCHEMA_VERSION = "1.0.0"
BUILD_RECEIPT_SCHEMA = "intersubmod.exact_ps_observation_report.build_receipt"
BUILD_RECEIPT_VERSION = "1.0.0"
AUTHORITY_SCHEMA = "intersubmod.exact_ps_tree_research.ai_handoff_authority"
STRICT_READY_SCHEMA = "intersubmod.strict_region_all7_bundle_ready"
REPORT_STATUS = "validated-derived-observation"
EXPECTED_ARTIFACT_IDS = {
    "strict_linkage_ready",
    "all7_cohort_receipt",
    "all7_summary",
    "all7_cpp_reader_report",
    "topology_receipt",
    "topology_summary",
    "topology_reader_report",
    "candidate_factorization_receipt",
    "read_af_decision_report",
    "methyl_manifest",
    "methyl_sidecar_receipt",
    "methyl_report_data",
    "methyl_reader_report",
}
REQUIRED_JSON_ARTIFACT_IDS = {
    "strict_linkage_ready",
    "all7_cohort_receipt",
    "all7_summary",
    "topology_receipt",
    "topology_summary",
    "candidate_factorization_receipt",
    "methyl_sidecar_receipt",
    "methyl_report_data",
}
EXPECTED_DENOMINATORS = {
    "ssnv_dataset_records": (469849, 469849),
    "strict_components": (255752, 255752),
    "k1_linkage_abstain": (170131, 255752),
    "strict_read_linked_w": (85621, 255752),
    "final_groups": (98955, 98955),
    "no_active_alt": (13014, 98955),
    "mutation_bearing_units": (85941, 98955),
    "family_complete": (75224, 85941),
    "resource_abstain": (10717, 85941),
    "ranked_complete": (71955, 75224),
    "zero_denominator": (3224, 75224),
    "af_recurrence_screen_abstain": (45, 75224),
    "unique_best_tree": (39648, 71955),
    "tied_same_rooted_unlabeled_topology": (23858, 71955),
    "tied_cross_topology": (8449, 71955),
    "one_rooted_unlabeled_topology": (63506, 71955),
    "methyl_formal_units": (1045, 1045),
    "methyl_evaluable_units": (811, 1045),
    "methyl_robust_associations": (3, 811),
}
RESOLUTION_KEYS = (
    "UNIQUE_TREE",
    "TIED_SAME_TOPOLOGY",
    "TIED_CROSS_TOPOLOGY",
)
COARSE_KEYS = (
    "Single-only",
    "Sister-only",
    "Direct-only",
    "Sister+direct",
)


class BuildError(RuntimeError):
    """A fail-closed report build error."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BuildError(message)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing input: {path}")
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_path(resolved),
    }


def read_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing JSON: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise BuildError(f"cannot read JSON {path}: {exc}") from exc
    require(isinstance(payload, dict), f"JSON root is not an object: {path}")
    return payload


def write_text_exclusive(path: Path, text: str, encoding: str = "utf-8") -> None:
    try:
        with path.open("x", encoding=encoding) as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise BuildError(f"refusing to overwrite output: {path}") from exc


def write_json_exclusive(path: Path, payload: dict[str, Any]) -> None:
    text = json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        indent=2,
        sort_keys=True,
    ) + "\n"
    write_text_exclusive(path, text)


def write_sidecar_exclusive(path: Path) -> Path:
    sidecar = path.with_name(f"{path.name}.sha256")
    write_text_exclusive(
        sidecar,
        f"{sha256_path(path)}  {path.name}\n",
        encoding="ascii",
    )
    return sidecar


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--authority-manifest", type=Path, required=True)
    parser.add_argument("--denominator-registry", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def as_int(value: Any, label: str) -> int:
    require(not isinstance(value, bool), f"{label} must be an integer")
    try:
        converted = int(value)
    except (TypeError, ValueError) as exc:
        raise BuildError(f"{label} must be an integer: {value!r}") from exc
    require(str(converted) == str(value).strip(), f"{label} is not canonical integer")
    return converted


def index_unique(
    rows: Iterable[dict[str, Any]],
    key: str,
    label: str,
) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for row in rows:
        require(isinstance(row, dict), f"{label} row is not an object")
        value = row.get(key)
        require(isinstance(value, str) and value, f"{label} row lacks {key}")
        require(value not in result, f"duplicate {label} {key}: {value}")
        result[value] = row
    return result


def load_authority_manifest(
    path: Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]], list[dict[str, Any]]]:
    manifest = read_json(path)
    require(
        manifest.get("schema_name") == AUTHORITY_SCHEMA,
        "authority manifest schema_name drift",
    )
    require(manifest.get("schema_version") == "1.0.0", "authority schema drift")
    artifacts_raw = manifest.get("artifacts")
    require(isinstance(artifacts_raw, list), "authority artifacts is not an array")
    artifacts = index_unique(artifacts_raw, "artifact_id", "authority artifact")
    require(
        set(artifacts) == EXPECTED_ARTIFACT_IDS,
        "authority artifact set drift: "
        f"missing={sorted(EXPECTED_ARTIFACT_IDS - set(artifacts))}, "
        f"extra={sorted(set(artifacts) - EXPECTED_ARTIFACT_IDS)}",
    )
    verified: list[dict[str, Any]] = []
    for artifact_id in sorted(artifacts):
        record = artifacts[artifact_id]
        artifact_path = Path(str(record.get("path", ""))).expanduser()
        observed = identity(artifact_path)
        expected_sha = record.get("sha256")
        require(
            isinstance(expected_sha, str) and len(expected_sha) == 64,
            f"authority SHA malformed: {artifact_id}",
        )
        require(
            observed["sha256"] == expected_sha,
            f"authority SHA mismatch: {artifact_id}",
        )
        verified.append(
            {
                "artifact_id": artifact_id,
                "role": record.get("role"),
                "status": record.get("status"),
                **observed,
            }
        )
    for artifact_id in REQUIRED_JSON_ARTIFACT_IDS:
        read_json(Path(artifacts[artifact_id]["path"]))
    return manifest, artifacts, verified


def load_denominator_registry(path: Path) -> list[dict[str, Any]]:
    require(path.is_file(), f"missing denominator registry: {path}")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise BuildError(f"cannot read denominator registry: {exc}") from exc
    require(len(rows) == len(EXPECTED_DENOMINATORS), "denominator row count drift")
    indexed = index_unique(rows, "metric", "denominator")
    require(set(indexed) == set(EXPECTED_DENOMINATORS), "denominator metric set drift")
    normalized: list[dict[str, Any]] = []
    for metric, expected in EXPECTED_DENOMINATORS.items():
        row = indexed[metric]
        numerator = as_int(row.get("numerator"), f"{metric}.numerator")
        denominator = as_int(row.get("denominator"), f"{metric}.denominator")
        require((numerator, denominator) == expected, f"{metric} frozen count drift")
        require(row.get("status") == "current", f"{metric} is not current")
        observed_pct = float(row["percent"])
        expected_pct = 100.0 * numerator / denominator
        require(
            math.isclose(observed_pct, expected_pct, abs_tol=0.0001),
            f"{metric} percent is inconsistent",
        )
        normalized.append(
            {
                "metric": metric,
                "numerator": numerator,
                "denominator": denominator,
                "percent": observed_pct,
                "scope": row["scope"],
                "status": row["status"],
                "authority_artifact_id": row["authority_artifact_id"],
                "parent_metric": None
                if row["parent_metric"] == "NA"
                else row["parent_metric"],
                "excluded_count": as_int(
                    row["excluded_count"], f"{metric}.excluded_count"
                ),
                "unknown_count": as_int(
                    row["unknown_count"], f"{metric}.unknown_count"
                ),
                "exclusion_reason": None
                if row["exclusion_reason"] == "NA"
                else row["exclusion_reason"],
                "interpretation": row["interpretation"],
            }
        )
    return normalized


def nested_strict_data(
    ready_path: Path,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    ready = read_json(ready_path)
    require(ready.get("schema_name") == STRICT_READY_SCHEMA, "strict READY schema drift")
    require(ready.get("schema_version") == "1.0.0", "strict READY version drift")
    require(ready.get("all_pass") is True, "strict READY is not all-pass")
    ready_checks = ready.get("checks")
    require(
        isinstance(ready_checks, dict)
        and ready_checks
        and all(value is True for value in ready_checks.values()),
        "strict READY checks are not all-pass",
    )
    ready_scope = ready.get("scope")
    require(isinstance(ready_scope, dict), "strict READY scope missing")
    require(
        ready_scope.get("dataset_chromosome_rows") == 154,
        "strict READY dataset-chromosome scope drift",
    )
    require(
        ready_scope.get("task_type") == "B_comprehensive_validation",
        "strict READY task type drift",
    )
    bundle_files = ready.get("bundle_files")
    require(isinstance(bundle_files, dict), "strict READY bundle_files missing")
    expected_bundle_names = {
        "artifact",
        "artifact_receipt",
        "artifact_sidecar",
        "data",
        "data_receipt",
        "data_sidecar",
        "delivery_receipt",
        "html",
        "visual_qa_receipt",
    }
    require(
        set(bundle_files) == expected_bundle_names,
        "strict READY bundle identity set drift",
    )
    verified_bundle: list[dict[str, Any]] = []
    for bundle_name in sorted(bundle_files):
        record = bundle_files[bundle_name]
        require(isinstance(record, dict), f"strict bundle record invalid: {bundle_name}")
        bundle_path = Path(str(record.get("path", "")))
        observed_bundle = identity(bundle_path)
        require(
            observed_bundle["sha256"] == record.get("sha256"),
            f"strict nested bundle SHA mismatch: {bundle_name}",
        )
        require(
            observed_bundle["size_bytes"] == record.get("size_bytes"),
            f"strict nested bundle size mismatch: {bundle_name}",
        )
        verified_bundle.append({"bundle_name": bundle_name, **observed_bundle})
    data_record = bundle_files.get("data")
    require(isinstance(data_record, dict), "strict READY data identity missing")
    data_path = Path(str(data_record.get("path", "")))
    observed = identity(data_path)
    data = read_json(data_path)
    require(
        data.get("schema_name") == "intersubmod.strict_region_all7_report_data",
        "strict report-data schema drift",
    )
    return data, observed, verified_bundle


def count_methyl_by_sample(
    cases: list[dict[str, Any]], samples: list[str]
) -> dict[str, dict[str, int]]:
    result = {
        sample: {
            "formal_units": 0,
            "evaluable_units": 0,
            "not_evaluable_units": 0,
            "robust_associations": 0,
            "evaluable_no_robust": 0,
            "confounded": 0,
        }
        for sample in samples
    }
    allowed_assessments = {
        "NOT_EVALUABLE",
        "ROBUST_ASSOCIATION",
        "EVALUABLE_NO_ROBUST_ASSOCIATION",
        "CONFOUNDED",
        "LOCAL_CIS_COMPATIBLE",
        "TAG_DEPENDENT",
    }
    for case in cases:
        require(isinstance(case, dict), "methyl case is not an object")
        sample = case.get("dataset")
        assessment = case.get("assessment")
        require(sample in result, f"methyl case sample outside scope: {sample}")
        require(assessment in allowed_assessments, f"unknown methyl assessment: {assessment}")
        target = result[sample]
        target["formal_units"] += 1
        if assessment == "NOT_EVALUABLE":
            target["not_evaluable_units"] += 1
        else:
            target["evaluable_units"] += 1
        if assessment == "ROBUST_ASSOCIATION":
            target["robust_associations"] += 1
        elif assessment == "EVALUABLE_NO_ROBUST_ASSOCIATION":
            target["evaluable_no_robust"] += 1
        elif assessment == "CONFOUNDED":
            target["confounded"] += 1
    return result


def n(record: dict[str, Any], key: str) -> int:
    value = record.get(key)
    require(isinstance(value, int) and not isinstance(value, bool), f"{key} not integer")
    return value


def build_report_data(
    manifest_path: Path,
    registry_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    manifest, artifacts, verified_artifacts = load_authority_manifest(manifest_path)
    registry = load_denominator_registry(registry_path)
    require(
        all(row["authority_artifact_id"] in artifacts for row in registry),
        "denominator registry references an unknown authority artifact",
    )
    registry_by_metric = {row["metric"]: row for row in registry}
    for metric, denominator in manifest.get("denominators", {}).items():
        if metric not in registry_by_metric:
            continue
        row = registry_by_metric[metric]
        require(row["numerator"] == denominator.get("count"), f"{metric} manifest drift")
        require(
            row["denominator"] == denominator.get("denominator", denominator.get("count")),
            f"{metric} manifest denominator drift",
        )

    strict_data, strict_data_identity, strict_bundle_identities = nested_strict_data(
        Path(artifacts["strict_linkage_ready"]["path"])
    )
    cohort_receipt = read_json(Path(artifacts["all7_cohort_receipt"]["path"]))
    all7 = read_json(Path(artifacts["all7_summary"]["path"]))
    topology = read_json(Path(artifacts["topology_summary"]["path"]))
    factorization = read_json(
        Path(artifacts["candidate_factorization_receipt"]["path"])
    )
    methyl = read_json(Path(artifacts["methyl_report_data"]["path"]))

    scope = manifest.get("scope")
    require(isinstance(scope, dict), "authority scope missing")
    samples = scope.get("technical_datasets")
    require(
        isinstance(samples, list)
        and len(samples) == 7
        and all(isinstance(value, str) for value in samples),
        "technical dataset scope drift",
    )
    require(len(set(samples)) == 7, "duplicate technical datasets")

    strict_by_sample = index_unique(strict_data.get("datasets", []), "dataset", "strict")
    all7_by_sample = index_unique(all7.get("samples", []), "sample", "all7")
    topology_by_sample = index_unique(topology.get("samples", []), "sample", "topology")
    require(set(strict_by_sample) == set(samples), "strict sample scope mismatch")
    require(set(all7_by_sample) == set(samples), "all7 sample scope mismatch")
    require(set(topology_by_sample) == set(samples), "topology sample scope mismatch")
    require(
        cohort_receipt.get("technical_all_pass") is True,
        "cohort receipt technical status failed",
    )
    require(
        cohort_receipt.get("validation_evidence_eligible") is False,
        "cohort validation-evidence status drift",
    )
    require(all7.get("technical_all_pass") is True, "all7 summary technical status failed")
    require(
        topology.get("cohort", {}).get("checks", {}).get("all_sample_checks_pass")
        is True,
        "topology cohort checks failed",
    )
    require(factorization.get("all_pass") is True, "factorization receipt failed")
    require(
        all(value is True for value in factorization.get("checks", {}).values()),
        "factorization checks failed",
    )

    methyl_cases = methyl.get("cases")
    require(isinstance(methyl_cases, list), "methyl cases missing")
    methyl_by_sample = count_methyl_by_sample(methyl_cases, samples)
    aggregate_methyl = methyl.get("aggregate")
    require(isinstance(aggregate_methyl, dict), "methyl aggregate missing")
    require(
        sum(row["formal_units"] for row in methyl_by_sample.values())
        == n(aggregate_methyl, "analysis_units")
        == 1045,
        "methyl formal-unit conservation failed",
    )
    require(
        sum(row["evaluable_units"] for row in methyl_by_sample.values())
        == n(aggregate_methyl, "analysis_evaluable")
        == 811,
        "methyl evaluable conservation failed",
    )
    require(
        sum(row["robust_associations"] for row in methyl_by_sample.values()) == 3,
        "methyl robust conservation failed",
    )

    per_sample: list[dict[str, Any]] = []
    for sample in samples:
        strict = strict_by_sample[sample]
        structural = all7_by_sample[sample]
        topo = topology_by_sample[sample]
        require(
            n(strict, "tree_eligible_W") == n(strict, "HP1_W") + n(strict, "HP2_W"),
            f"{sample} HP1/HP2 W conservation failed",
        )
        require(
            n(structural, "groups_total")
            == n(structural, "no_active_alt_units")
            + n(structural, "mutation_bearing_units"),
            f"{sample} group partition failed",
        )
        require(
            n(structural, "mutation_bearing_units")
            == n(structural, "mutation_family_complete_units")
            + n(structural, "mutation_family_abstain_units"),
            f"{sample} family partition failed",
        )
        require(
            n(structural, "mutation_family_complete_units")
            == n(structural, "ranked_units")
            + n(structural, "zero_denominator_units")
            + n(structural, "recurrence_required_units"),
            f"{sample} AF outcome partition failed",
        )
        resolution = topo.get("resolution")
        require(isinstance(resolution, dict), f"{sample} topology resolution missing")
        resolution_counts = {
            key: n(resolution[key], "n") for key in RESOLUTION_KEYS
        }
        ranked = n(structural, "ranked_units")
        require(sum(resolution_counts.values()) == ranked, f"{sample} resolution drift")
        require(
            resolution_counts["UNIQUE_TREE"]
            == n(structural, "best_tree_unique_units"),
            f"{sample} unique-tree drift",
        )
        require(
            resolution_counts["TIED_SAME_TOPOLOGY"]
            + resolution_counts["TIED_CROSS_TOPOLOGY"]
            == n(structural, "best_tree_tied_units"),
            f"{sample} tied-tree drift",
        )
        coarse_raw = topo.get("resolved_coarse_class")
        require(isinstance(coarse_raw, dict), f"{sample} coarse class missing")
        coarse = {key: n(coarse_raw[key], "n") for key in COARSE_KEYS}
        one_topology = n(topo.get("one_exact_topology", {}), "n")
        require(
            one_topology
            == resolution_counts["UNIQUE_TREE"]
            + resolution_counts["TIED_SAME_TOPOLOGY"],
            f"{sample} one-topology identity failed",
        )
        require(sum(coarse.values()) >= one_topology, f"{sample} coarse class underflow")
        active_all = {
            str(key): int(value)
            for key, value in structural.get("active_k_distribution", {}).items()
        }
        active_ranked = {
            str(key): n(value, "n") for key, value in topo.get("active_k", {}).items()
        }
        require(sum(active_all.values()) == n(structural, "groups_total"), f"{sample} k drift")
        require(sum(active_ranked.values()) == ranked, f"{sample} ranked k drift")
        stage = structural.get("stage_wall_seconds")
        require(isinstance(stage, dict), f"{sample} stage timing missing")
        per_sample.append(
            {
                "sample": sample,
                "biological_id": "HCC1395"
                if sample == "HCC1395_DORADO"
                else sample,
                "technical_replicate": sample == "HCC1395_DORADO",
                "input": {
                    "candidate_loci_s": n(strict, "candidate_loci_S"),
                    "strict_components": n(strict, "all_components"),
                    "k1_components": n(strict, "singleton_k1_components"),
                    "strict_read_linked_w": n(strict, "tree_eligible_W"),
                    "hp1_w": n(strict, "HP1_W"),
                    "hp2_w": n(strict, "HP2_W"),
                },
                "solver": {
                    "final_groups": n(structural, "groups_total"),
                    "no_active_alt": n(structural, "no_active_alt_units"),
                    "mutation_bearing": n(structural, "mutation_bearing_units"),
                    "family_complete": n(
                        structural, "mutation_family_complete_units"
                    ),
                    "resource_abstain": n(
                        structural, "mutation_family_abstain_units"
                    ),
                    "ranked": ranked,
                    "zero_denominator": n(structural, "zero_denominator_units"),
                    "recurrence_af_abstain": n(
                        structural, "recurrence_required_units"
                    ),
                },
                "topology": {
                    "resolution": resolution_counts,
                    "one_exact_topology": one_topology,
                    "one_exact_topology_pct_ranked": 100.0
                    * one_topology
                    / ranked,
                    "resolved_coarse_class": coarse,
                    "best_tree_total": n(topo, "best_tree_total"),
                },
                "active_k": {
                    "all_groups": active_all,
                    "ranked": active_ranked,
                },
                "runtime_seconds": {
                    "adapter": float(stage["exact_ps_partition_to_mlhp"]),
                    "cpp_topology": float(stage["cpp_topology"]),
                    "stage_sum": float(stage["total"]),
                },
                "methyl": methyl_by_sample[sample],
            }
        )

    totals = all7.get("totals")
    cohort_topology = topology.get("cohort")
    require(isinstance(totals, dict), "all7 totals missing")
    require(isinstance(cohort_topology, dict), "topology cohort missing")
    expected_all7_totals = {
        "groups_total": 98955,
        "no_active_alt_units": 13014,
        "mutation_bearing_units": 85941,
        "mutation_family_complete_units": 75224,
        "mutation_family_abstain_units": 10717,
        "ranked_units": 71955,
        "zero_denominator_units": 3224,
        "recurrence_required_units": 45,
        "best_tree_unique_units": 39648,
        "best_tree_tied_units": 32307,
    }
    for field, expected in expected_all7_totals.items():
        require(n(totals, field) == expected, f"all7 total drift: {field}")
    sample_total_observed = {
        "groups_total": sum(row["solver"]["final_groups"] for row in per_sample),
        "no_active_alt_units": sum(row["solver"]["no_active_alt"] for row in per_sample),
        "mutation_bearing_units": sum(
            row["solver"]["mutation_bearing"] for row in per_sample
        ),
        "mutation_family_complete_units": sum(
            row["solver"]["family_complete"] for row in per_sample
        ),
        "mutation_family_abstain_units": sum(
            row["solver"]["resource_abstain"] for row in per_sample
        ),
        "ranked_units": sum(row["solver"]["ranked"] for row in per_sample),
        "zero_denominator_units": sum(
            row["solver"]["zero_denominator"] for row in per_sample
        ),
        "recurrence_required_units": sum(
            row["solver"]["recurrence_af_abstain"] for row in per_sample
        ),
        "best_tree_unique_units": sum(
            row["topology"]["resolution"]["UNIQUE_TREE"] for row in per_sample
        ),
        "best_tree_tied_units": sum(
            row["topology"]["resolution"]["TIED_SAME_TOPOLOGY"]
            + row["topology"]["resolution"]["TIED_CROSS_TOPOLOGY"]
            for row in per_sample
        ),
    }
    require(
        sample_total_observed == expected_all7_totals,
        "per-sample totals do not reproduce all7 totals",
    )
    cohort_resolution = {
        key: n(cohort_topology["resolution"][key], "n") for key in RESOLUTION_KEYS
    }
    cohort_coarse = {
        key: n(cohort_topology["resolved_coarse_class"][key], "n")
        for key in COARSE_KEYS
    }
    factor_totals = factorization.get("totals")
    require(isinstance(factor_totals, dict), "factorization totals missing")

    checks: list[dict[str, Any]] = []

    def add_check(check_id: str, observed: int | bool, expected: int | bool) -> None:
        passed = observed == expected
        checks.append(
            {
                "check_id": check_id,
                "observed": observed,
                "expected": expected,
                "status": "PASS" if passed else "FAIL",
            }
        )
        require(passed, f"conservation check failed: {check_id}")

    add_check(
        "candidate_loci_sum_samples",
        sum(row["input"]["candidate_loci_s"] for row in per_sample),
        469849,
    )
    add_check(
        "strict_components_partition",
        sum(row["input"]["strict_components"] for row in per_sample),
        255752,
    )
    add_check(
        "strict_k1_sum_samples",
        sum(row["input"]["k1_components"] for row in per_sample),
        170131,
    )
    add_check(
        "strict_w_sum_samples",
        sum(row["input"]["strict_read_linked_w"] for row in per_sample),
        85621,
    )
    add_check(
        "strict_components_k1_plus_w",
        170131 + 85621,
        255752,
    )
    add_check(
        "final_groups_sum_samples",
        sum(row["solver"]["final_groups"] for row in per_sample),
        98955,
    )
    add_check("final_groups_no_alt_plus_mutation", 13014 + 85941, 98955)
    add_check("mutation_family_complete_plus_abstain", 75224 + 10717, 85941)
    add_check(
        "complete_family_af_outcomes",
        71955 + 3224 + 45,
        75224,
    )
    add_check(
        "ranked_resolution",
        sum(cohort_resolution.values()),
        71955,
    )
    add_check(
        "one_topology_resolution",
        cohort_resolution["UNIQUE_TREE"]
        + cohort_resolution["TIED_SAME_TOPOLOGY"],
        63506,
    )
    add_check(
        "coarse_resolution",
        sum(cohort_coarse.values()),
        63511,
    )
    add_check(
        "cross_exact_same_coarse_gap",
        sum(cohort_coarse.values()) - n(cohort_topology["one_exact_topology"], "n"),
        5,
    )
    add_check(
        "active_k_all_groups",
        sum(int(value) for value in totals["active_k_distribution"].values()),
        98955,
    )
    add_check(
        "active_k_ranked",
        sum(n(value, "n") for value in cohort_topology["active_k"].values()),
        71955,
    )
    add_check(
        "cohort_best_tree_total",
        n(cohort_topology, "best_tree_total"),
        n(factor_totals, "global_best_trees"),
    )
    add_check(
        "minimum_tree_count",
        n(factor_totals, "minimum_trees"),
        972592,
    )
    add_check(
        "methyl_assessment_partition",
        sum(aggregate_methyl["assessment_counts"].values()),
        1045,
    )
    add_check(
        "methyl_formal_partition",
        n(aggregate_methyl, "analysis_evaluable")
        + n(aggregate_methyl["assessment_counts"], "NOT_EVALUABLE"),
        1045,
    )
    add_check(
        "methyl_evaluable_partition",
        n(aggregate_methyl["assessment_counts"], "ROBUST_ASSOCIATION")
        + n(aggregate_methyl["assessment_counts"], "EVALUABLE_NO_ROBUST_ASSOCIATION")
        + n(aggregate_methyl["assessment_counts"], "CONFOUNDED"),
        811,
    )
    add_check(
        "authority_hashes",
        len(verified_artifacts),
        13,
    )
    add_check(
        "strict_nested_bundle_hashes",
        len(strict_bundle_identities),
        9,
    )
    add_check(
        "cohort_technical_status",
        cohort_receipt.get("technical_all_pass") is True,
        True,
    )

    registry_records = {
        row["metric"]: {
            "count": row["numerator"],
            "denominator": row["denominator"],
            "percent": row["percent"],
            "interpretation": row["interpretation"],
        }
        for row in registry
    }
    relationship_types = {
        "ssnv_dataset_records": "SOURCE_GRAIN",
        "strict_components": "GRAIN_CHANGE",
        "k1_linkage_abstain": "SAME_GRAIN_PARTITION",
        "strict_read_linked_w": "SAME_GRAIN_PARTITION",
        "final_groups": "ONE_TO_MANY_FAN_OUT",
        "no_active_alt": "SAME_GRAIN_PARTITION",
        "mutation_bearing_units": "SAME_GRAIN_PARTITION",
        "family_complete": "SAME_GRAIN_PARTITION",
        "resource_abstain": "SAME_GRAIN_PARTITION",
        "ranked_complete": "SAME_GRAIN_PARTITION",
        "zero_denominator": "SAME_GRAIN_PARTITION",
        "af_recurrence_screen_abstain": "SAME_GRAIN_PARTITION",
        "unique_best_tree": "SAME_GRAIN_PARTITION",
        "tied_same_rooted_unlabeled_topology": "SAME_GRAIN_PARTITION",
        "tied_cross_topology": "SAME_GRAIN_PARTITION",
        "one_rooted_unlabeled_topology": "DERIVED_UNION",
        "methyl_formal_units": "INDEPENDENT_SIDECAR_GRAIN",
        "methyl_evaluable_units": "SAME_GRAIN_PARTITION",
        "methyl_robust_associations": "SAME_GRAIN_SUBSET",
    }
    unit_names = {
        "ssnv_dataset_records": "sSNV dataset-record",
        "strict_components": "strict component",
        "k1_linkage_abstain": "strict component",
        "strict_read_linked_w": "strict component W",
        "final_groups": "bounded final group",
        "no_active_alt": "bounded final group",
        "mutation_bearing_units": "mutation-bearing group",
        "family_complete": "mutation-bearing group",
        "resource_abstain": "mutation-bearing group",
        "ranked_complete": "family-complete group",
        "zero_denominator": "family-complete group",
        "af_recurrence_screen_abstain": "family-complete group",
        "unique_best_tree": "ranked group",
        "tied_same_rooted_unlabeled_topology": "ranked group",
        "tied_cross_topology": "ranked group",
        "one_rooted_unlabeled_topology": "ranked group",
        "methyl_formal_units": "methyl formal unit",
        "methyl_evaluable_units": "methyl formal unit",
        "methyl_robust_associations": "methyl evaluable unit",
    }
    funnel_stages = [
        {
            "metric": row["metric"],
            "count": row["numerator"],
            "denominator": row["denominator"],
            "unit": unit_names[row["metric"]],
            "authority_artifact_id": row["authority_artifact_id"],
            "relationship_type": relationship_types[row["metric"]],
            "parent_metric": row["parent_metric"],
        }
        for row in registry
    ]
    data = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "as_of_date": manifest.get("as_of_date"),
        "report_status": REPORT_STATUS,
        "title": "Exact-PS × HP 區域突變狀態樹：全 7 技術資料集觀察",
        "subtitle": "從 read 串聯、最小候選族、read-AF 排序到拓撲可辨識度",
        "scope": {
            "task_type": "B_comprehensive_validation",
            "technical_datasets": samples,
            "technical_dataset_count": 7,
            "biological_id_count": scope.get("biological_ids"),
            "chromosomes": scope.get("chromosomes"),
            "chromosome_scope": "chr1–22",
            "primary_grain": scope.get("primary_grain"),
            "linkage_min_distinct_molecules": scope.get(
                "linkage_min_distinct_molecules"
            ),
            "maximum_original_block_k": scope.get("maximum_original_block_k"),
            "maximum_active_k": scope.get("maximum_active_k"),
            "technical_pair_note": scope.get("technical_pair_note"),
        },
        "claim_boundary": manifest.get("claim_boundary"),
        "model_contract": manifest.get("model_contract"),
        "headline": {
            "final_groups": 98955,
            "mutation_bearing": 85941,
            "resource_abstain": 10717,
            "ranked": 71955,
            "unique_best_tree": 39648,
            "one_rooted_unlabeled_topology": 63506,
            "tied_cross_topology": 8449,
            "methyl_robust_associations": 3,
        },
        "cohort_funnel": {
            "records": registry_records,
            "stages": funnel_stages,
            "grain_rules": [
                {
                    "from": "ssnv_dataset_records",
                    "to": "strict_components",
                    "relation": "grain_change",
                    "warning": "sSNV dataset-record 與 graph component 不是同一計數單位，不能作守恆流量。",
                },
                {
                    "from": "strict_read_linked_w",
                    "to": "final_groups",
                    "relation": "one_to_many_partition",
                    "warning": "一個 W 可切成多個 bounded blocks；98,955/85,621 是 fan-out，不是比例。",
                },
            ],
        },
        "tree_decision": {
            "policy": manifest.get("tree_decision_policy"),
            "materialization_status": "POLICY_ONLY_NOT_MATERIALIZED",
            "safe_empirical_buckets": {
                "no_mutation_tree_inference": 13014,
                "abstain_incomplete": 10717,
                "af_abstain_zero_denominator": 3224,
                "af_abstain_recurrence_screen": 45,
                "raw_af_unique_representative": 39648,
                "topology_representative_only": 23858,
                "candidate_set_only_cross_topology": 8449,
            },
            "warning": "目前 frozen row output 尚未把完整 tree_decision enum 逐列物化；此處只呈現可由權威統計直接守恆的桶。",
        },
        "topology": {
            "ranked_units": 71955,
            "resolution": cohort_resolution,
            "one_exact_topology": n(cohort_topology["one_exact_topology"], "n"),
            "one_exact_topology_pct_ranked": float(
                cohort_topology["one_exact_topology"]["pct_ranked"]
            ),
            "resolved_coarse_class": cohort_coarse,
            "cross_exact_but_same_coarse_class": n(
                cohort_topology, "cross_exact_topology_but_same_coarse_class"
            ),
            "maximum_best_tree_tie_count": n(
                cohort_topology, "maximum_best_tree_tie_count"
            ),
            "candidate_counts": {
                "minimum_vertex_sets": n(factor_totals, "minimum_vertex_sets"),
                "minimum_trees": n(factor_totals, "minimum_trees"),
                "global_af_best_trees": n(factor_totals, "global_best_trees"),
            },
            "classification_definitions": {
                "Single-only": "沒有分枝，且從根到任何節點都沒有兩層以上的非根路徑。",
                "Sister-only": "存在分枝，但沒有兩層以上的祖先—後代路徑。",
                "Direct-only": "存在兩層以上的路徑，但沒有分枝。",
                "Sister+direct": "同時存在分枝與兩層以上的路徑。",
            },
            "interpretation_ceiling": "這四類是 rooted graph geometry；不能直接翻成已確認的 clone 或 subclone。",
        },
        "active_k": {
            "all_groups": {
                str(key): int(value)
                for key, value in totals.get("active_k_distribution", {}).items()
            },
            "ranked": {
                str(key): n(value, "n")
                for key, value in cohort_topology.get("active_k", {}).items()
            },
            "guard": {
                "maximum_active_k": 12,
                "maximum_search_nodes": cohort_receipt["parameters"][
                    "max_search_nodes"
                ],
                "maximum_family_size": cohort_receipt["parameters"][
                    "max_family_size"
                ],
                "resource_abstain_is_nonrandom": True,
                "note": "10,717 個 canonical abstain 全部為 search-node guard；集中在較大的 active k。",
            },
        },
        "runtime": {
            "cohort_runner_wall_seconds": float(cohort_receipt["wall_seconds"]),
            "sample_stage_sum_seconds": float(totals["stage_wall_seconds"]["total"]),
            "solver_internal_seconds": float(
                totals["solver_elapsed_microseconds"]["sum"]
            )
            / 1_000_000.0,
            "clock_warning": "三個數值來自不同計時範圍，不能相加或互相當作百分比。",
            "parameters": cohort_receipt.get("parameters"),
        },
        "cn_loh_readiness": {
            "status": "NOT_INTEGRATED",
            "availability": None,
            "consequence": "read-AF 是 molecule-level、copy-unadjusted 排序；不得翻成 cell fraction、CCF 或已確認生物祖先關係。",
            "required_before_biological_single_tree": [
                "allele-specific copy number",
                "LOH / copy-neutral LOH screen",
                "purity and mutation multiplicity model",
                "AF uncertainty or likelihood calibration",
            ],
        },
        "methyl_auxiliary": {
            "status": "association-only",
            "formal_units": 1045,
            "evaluable_units": 811,
            "not_evaluable_units": 234,
            "robust_associations": 3,
            "evaluable_no_robust": 627,
            "confounded": 181,
            "assessment_counts": aggregate_methyl.get("assessment_counts"),
            "claim_ceiling": methyl.get("claim_ceiling"),
            "topology_rescoring": False,
            "node_clustering_materialized": False,
        },
        "per_sample": per_sample,
        "conservation_checks": checks,
        "provenance": {
            "authority_manifest": identity(manifest_path),
            "denominator_registry": identity(registry_path),
            "strict_nested_data": strict_data_identity,
            "strict_nested_bundle_identities": strict_bundle_identities,
            "strict_nested_bundle_count": len(strict_bundle_identities),
            "verified_authority_artifacts": verified_artifacts,
            "authority_artifact_count": len(verified_artifacts),
        },
        "forbidden_claims": manifest.get("claim_boundary", {}).get("forbidden"),
    }
    build_context = {
        "manifest": manifest,
        "verified_artifacts": verified_artifacts,
        "registry": registry,
    }
    return data, build_context


def fmt_int(value: int) -> str:
    return f"{value:,}"


def fmt_pct(numerator: int, denominator: int, digits: int = 2) -> str:
    return f"{100.0 * numerator / denominator:.{digits}f}%"


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def topology_svg(data: dict[str, Any]) -> str:
    rows = [
        ("全部", data["topology"]["resolution"], data["topology"]["ranked_units"]),
        *[
            (
                row["sample"],
                row["topology"]["resolution"],
                row["solver"]["ranked"],
            )
            for row in data["per_sample"]
        ],
    ]
    colors = {
        "UNIQUE_TREE": "#176b65",
        "TIED_SAME_TOPOLOGY": "#d49a34",
        "TIED_CROSS_TOPOLOGY": "#b44d3a",
    }
    labels = {
        "UNIQUE_TREE": "唯一最佳樹",
        "TIED_SAME_TOPOLOGY": "同拓撲並列",
        "TIED_CROSS_TOPOLOGY": "跨拓撲未解",
    }
    left = 170
    width = 850
    y0 = 92
    row_height = 46
    chunks = [
        '<svg class="chart-svg" viewBox="0 0 1160 500" role="img" '
        'aria-labelledby="topology-chart-title topology-chart-desc">',
        '<title id="topology-chart-title">各樣本 read-AF 最佳候選解析狀態</title>',
        '<desc id="topology-chart-desc">每一列以可排名區域為分母，顯示唯一最佳樹、同拓撲並列與跨拓撲未解比例。</desc>',
        '<g class="chart-legend">',
    ]
    x = left
    for key in RESOLUTION_KEYS:
        chunks.append(
            f'<rect x="{x}" y="28" width="18" height="18" rx="3" fill="{colors[key]}"/>'
            f'<text x="{x + 26}" y="43">{esc(labels[key])}</text>'
        )
        x += 230
    chunks.append("</g>")
    for index, (sample, counts, total) in enumerate(rows):
        y = y0 + index * row_height
        chunks.append(
            f'<text x="{left - 18}" y="{y + 17}" text-anchor="end" '
            f'class="chart-row-label">{esc(sample)}</text>'
        )
        cursor = left
        for key in RESOLUTION_KEYS:
            segment = width * counts[key] / total
            chunks.append(
                f'<rect x="{cursor:.2f}" y="{y}" width="{segment:.2f}" '
                f'height="24" fill="{colors[key]}"><title>'
                f'{esc(sample)} · {esc(labels[key])}: {fmt_int(counts[key])} '
                f'({fmt_pct(counts[key], total)})</title></rect>'
            )
            cursor += segment
        chunks.append(
            f'<text x="{left + width + 16}" y="{y + 17}" class="chart-total">'
            f'n={fmt_int(total)}</text>'
        )
    chunks.extend(
        [
            f'<line x1="{left}" y1="470" x2="{left + width}" y2="470" '
            'stroke="#918a80" stroke-width="1"/>',
            f'<text x="{left}" y="490">0%</text>',
            f'<text x="{left + width / 2}" y="490" text-anchor="middle">50%</text>',
            f'<text x="{left + width}" y="490" text-anchor="end">100%</text>',
            "</svg>",
        ]
    )
    return "".join(chunks)


def active_k_svg(data: dict[str, Any]) -> str:
    all_groups = data["active_k"]["all_groups"]
    ranked = data["active_k"]["ranked"]
    values = [all_groups.get(str(k), 0) for k in range(13)]
    max_value = max(values)
    chart_left = 75
    chart_bottom = 315
    chart_height = 235
    bar_width = 55
    gap = 19
    chunks = [
        '<svg class="chart-svg" viewBox="0 0 1080 370" role="img" '
        'aria-labelledby="k-chart-title k-chart-desc">',
        '<title id="k-chart-title">final groups 的 active k 分布與可排名子集</title>',
        '<desc id="k-chart-desc">每個 active k 顯示全部 final groups，深色內條顯示其中可用 read-AF 排名者；k 越大資源 abstain 越集中。</desc>',
        f'<line x1="{chart_left}" y1="{chart_bottom}" x2="1045" '
        f'y2="{chart_bottom}" stroke="#635f58" stroke-width="1"/>',
    ]
    for k in range(13):
        x = chart_left + k * (bar_width + gap)
        total = all_groups.get(str(k), 0)
        ranked_n = ranked.get(str(k), 0)
        height = chart_height * total / max_value
        ranked_height = chart_height * ranked_n / max_value
        y = chart_bottom - height
        ranked_y = chart_bottom - ranked_height
        chunks.append(
            f'<rect x="{x}" y="{y:.2f}" width="{bar_width}" height="{height:.2f}" '
            'rx="5" fill="#d9d2c5"><title>'
            f'active k={k}: 全部 {fmt_int(total)}</title></rect>'
        )
        chunks.append(
            f'<rect x="{x + 9}" y="{ranked_y:.2f}" width="{bar_width - 18}" '
            f'height="{ranked_height:.2f}" rx="4" fill="#176b65"><title>'
            f'active k={k}: 可排名 {fmt_int(ranked_n)}</title></rect>'
        )
        chunks.append(
            f'<text x="{x + bar_width / 2}" y="{chart_bottom + 24}" '
            f'text-anchor="middle">{k}</text>'
        )
        if total >= 1000:
            chunks.append(
                f'<text x="{x + bar_width / 2}" y="{max(18, y - 8):.2f}" '
                f'text-anchor="middle" class="chart-small">{fmt_int(total)}</text>'
            )
    chunks.extend(
        [
            '<rect x="760" y="30" width="18" height="18" rx="3" fill="#d9d2c5"/>',
            '<text x="786" y="44">全部 final groups</text>',
            '<rect x="905" y="30" width="18" height="18" rx="3" fill="#176b65"/>',
            '<text x="931" y="44">可排名</text>',
            '<text x="560" y="359" text-anchor="middle" class="axis-label">active k</text>',
            "</svg>",
        ]
    )
    return "".join(chunks)


def funnel_svg(data: dict[str, Any]) -> str:
    return f"""
<svg class="chart-svg funnel-svg" viewBox="0 0 1180 430" role="img"
     aria-labelledby="funnel-title funnel-desc">
  <title id="funnel-title">Exact-PS 分母與決策流程</title>
  <desc id="funnel-desc">圖中明確區分 sSNV、component、W 與 final group 的不同計數單位；W 到 final group 是一對多切割，不是守恆漏斗。</desc>
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5"
            markerWidth="7" markerHeight="7" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#635f58"/>
    </marker>
    <marker id="arrow-rust" viewBox="0 0 10 10" refX="9" refY="5"
            markerWidth="7" markerHeight="7" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#b44d3a"/>
    </marker>
  </defs>
  <g class="flow-box">
    <rect x="18" y="32" width="190" height="80" rx="10"/>
    <text x="113" y="60" text-anchor="middle" class="flow-number">469,849</text>
    <text x="113" y="86" text-anchor="middle">sSNV dataset-records</text>
  </g>
  <line x1="208" y1="72" x2="275" y2="72" stroke="#635f58" stroke-dasharray="6 5"
        marker-end="url(#arrow)"/>
  <text x="242" y="52" text-anchor="middle" class="flow-note">粒度改變</text>
  <g class="flow-box">
    <rect x="280" y="32" width="190" height="80" rx="10"/>
    <text x="375" y="60" text-anchor="middle" class="flow-number">255,752</text>
    <text x="375" y="86" text-anchor="middle">strict components</text>
  </g>
  <path d="M470 72 H515 V172 H570" fill="none" stroke="#635f58"
        marker-end="url(#arrow)"/>
  <path d="M470 72 H515 V72 H570" fill="none" stroke="#635f58"
        marker-end="url(#arrow)"/>
  <g class="flow-box muted">
    <rect x="575" y="32" width="190" height="80" rx="10"/>
    <text x="670" y="60" text-anchor="middle" class="flow-number">170,131</text>
    <text x="670" y="86" text-anchor="middle">k=1；不建多點樹</text>
  </g>
  <g class="flow-box">
    <rect x="575" y="132" width="190" height="80" rx="10"/>
    <text x="670" y="160" text-anchor="middle" class="flow-number">85,621</text>
    <text x="670" y="186" text-anchor="middle">read-linked W（k≥2）</text>
  </g>
  <line x1="765" y1="172" x2="835" y2="172" stroke="#b44d3a"
        stroke-width="2" stroke-dasharray="7 5" marker-end="url(#arrow-rust)"/>
  <text x="800" y="144" text-anchor="middle" class="flow-note rust">一對多切割</text>
  <g class="flow-box accent">
    <rect x="840" y="132" width="210" height="80" rx="10"/>
    <text x="945" y="160" text-anchor="middle" class="flow-number">98,955</text>
    <text x="945" y="186" text-anchor="middle">final topology groups</text>
  </g>
  <path d="M945 212 V258 H710" fill="none" stroke="#635f58"
        marker-end="url(#arrow)"/>
  <path d="M945 212 V258 H970" fill="none" stroke="#635f58"
        marker-end="url(#arrow)"/>
  <g class="flow-box muted">
    <rect x="495" y="268" width="220" height="78" rx="10"/>
    <text x="605" y="296" text-anchor="middle" class="flow-number">13,014</text>
    <text x="605" y="321" text-anchor="middle">no active ALT</text>
  </g>
  <g class="flow-box">
    <rect x="865" y="268" width="220" height="78" rx="10"/>
    <text x="975" y="296" text-anchor="middle" class="flow-number">85,941</text>
    <text x="975" y="321" text-anchor="middle">mutation-bearing</text>
  </g>
  <path d="M975 346 V376 H765" fill="none" stroke="#635f58"
        marker-end="url(#arrow)"/>
  <path d="M975 346 V376 H1085" fill="none" stroke="#635f58"
        marker-end="url(#arrow)"/>
  <text x="680" y="408" text-anchor="middle" class="flow-final">
    75,224 family complete → 71,955 ranked
  </text>
  <text x="1025" y="408" text-anchor="middle" class="flow-final rust">
    10,717 resource abstain
  </text>
</svg>
"""


def sample_table_rows(data: dict[str, Any]) -> str:
    rows = []
    for row in data["per_sample"]:
        resolution = row["topology"]["resolution"]
        solver = row["solver"]
        methyl = row["methyl"]
        rows.append(
            f"""
<tr data-sample="{esc(row['sample'])}">
  <th scope="row">
    <span class="sample-name">{esc(row['sample'])}</span>
    {"<span class='tech-tag'>technical pair</span>" if row["technical_replicate"] else ""}
  </th>
  <td data-field="strict_read_linked_w" data-value="{row["input"]["strict_read_linked_w"]}">{fmt_int(row["input"]["strict_read_linked_w"])}</td>
  <td data-field="final_groups" data-value="{solver["final_groups"]}">{fmt_int(solver["final_groups"])}</td>
  <td data-field="resource_abstain" data-value="{solver["resource_abstain"]}">{fmt_int(solver["resource_abstain"])}
      <span class="subvalue">{fmt_pct(solver["resource_abstain"], solver["mutation_bearing"])}</span></td>
  <td data-field="ranked" data-value="{solver["ranked"]}">{fmt_int(solver["ranked"])}</td>
  <td data-field="unique_tree" data-value="{resolution["UNIQUE_TREE"]}">{fmt_int(resolution["UNIQUE_TREE"])}
      <span class="subvalue">{fmt_pct(resolution["UNIQUE_TREE"], solver["ranked"])}</span></td>
  <td data-field="tied_same_topology" data-value="{resolution["TIED_SAME_TOPOLOGY"]}">{fmt_int(resolution["TIED_SAME_TOPOLOGY"])}
      <span class="subvalue">{fmt_pct(resolution["TIED_SAME_TOPOLOGY"], solver["ranked"])}</span></td>
  <td data-field="tied_cross_topology" data-value="{resolution["TIED_CROSS_TOPOLOGY"]}">{fmt_int(resolution["TIED_CROSS_TOPOLOGY"])}
      <span class="subvalue">{fmt_pct(resolution["TIED_CROSS_TOPOLOGY"], solver["ranked"])}</span></td>
  <td data-field="one_exact_topology" data-value="{row["topology"]["one_exact_topology"]}">{fmt_int(row["topology"]["one_exact_topology"])}
      <span class="subvalue">{row["topology"]["one_exact_topology_pct_ranked"]:.2f}%</span></td>
  <td data-field="methyl_formal" data-value="{methyl["formal_units"]}">{fmt_int(methyl["formal_units"])} /
      <span data-field="methyl_robust" data-value="{methyl["robust_associations"]}">{fmt_int(methyl["robust_associations"])}</span></td>
</tr>
"""
        )
    return "".join(rows)


def provenance_rows(data: dict[str, Any]) -> str:
    chunks = []
    for artifact in data["provenance"]["verified_authority_artifacts"]:
        chunks.append(
            f"""
<tr>
  <th scope="row"><code>{esc(artifact["artifact_id"])}</code></th>
  <td>{esc(artifact["status"] or "documented")}</td>
  <td><code>{esc(artifact["sha256"][:16])}…</code></td>
  <td class="path-cell">{esc(artifact["path"])}</td>
</tr>
"""
        )
    return "".join(chunks)


def coarse_cards(data: dict[str, Any]) -> str:
    definitions = data["topology"]["classification_definitions"]
    counts = data["topology"]["resolved_coarse_class"]
    icons = {
        "Single-only": "○",
        "Sister-only": "⑂",
        "Direct-only": "↳",
        "Sister+direct": "⑂↳",
    }
    return "".join(
        f"""
<article class="definition-card">
  <div class="definition-icon" aria-hidden="true">{esc(icons[key])}</div>
  <h3>{esc(key)}</h3>
  <p class="definition-count">{fmt_int(counts[key])}</p>
  <p>{esc(definitions[key])}</p>
</article>
"""
        for key in COARSE_KEYS
    )


def render_html(data: dict[str, Any]) -> str:
    embedded = json.dumps(
        data,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).replace("</", "<\\/")
    resolution = data["topology"]["resolution"]
    coarse_resolved = sum(data["topology"]["resolved_coarse_class"].values())
    runtime = data["runtime"]
    return f"""<!doctype html>
<html lang="zh-Hant">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta http-equiv="Content-Security-Policy"
        content="default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'; img-src data:; base-uri 'none'; form-action 'none'; object-src 'none'">
  <meta name="color-scheme" content="light">
  <title>Exact-PS × HP｜全 7 資料集觀察</title>
  <style>
    :root {{
      --paper: #f5f0e6;
      --paper-deep: #ece3d3;
      --ink: #1d2826;
      --muted: #655f57;
      --line: #cbc0ae;
      --teal: #176b65;
      --teal-dark: #0f4e4a;
      --teal-soft: #dcebe7;
      --rust: #b44d3a;
      --rust-soft: #f1ddd6;
      --gold: #d49a34;
      --gold-soft: #f2e5c7;
      --white: #fffdfa;
      --shadow: 0 14px 34px rgba(49, 43, 35, .09);
      --serif: Georgia, "Times New Roman", "Noto Serif TC", serif;
      --sans: -apple-system, BlinkMacSystemFont, "Segoe UI", "Noto Sans TC", sans-serif;
    }}
    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; background: var(--paper); }}
    body {{
      margin: 0;
      min-width: 0;
      color: var(--ink);
      background:
        linear-gradient(90deg, rgba(23,107,101,.035) 1px, transparent 1px) 0 0 / 34px 34px,
        linear-gradient(rgba(23,107,101,.025) 1px, transparent 1px) 0 0 / 34px 34px,
        var(--paper);
      font-family: var(--sans);
      line-height: 1.65;
      overflow-x: hidden;
    }}
    a {{ color: var(--teal-dark); text-underline-offset: 3px; }}
    code {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: .88em;
      overflow-wrap: anywhere;
    }}
    .skip-link {{
      position: fixed;
      top: .5rem;
      left: .5rem;
      transform: translateY(-160%);
      background: var(--ink);
      color: white;
      padding: .6rem .8rem;
      z-index: 30;
    }}
    .skip-link:focus {{ transform: none; }}
    .topbar {{
      position: sticky;
      top: 0;
      z-index: 20;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 1rem;
      min-height: 58px;
      padding: .55rem clamp(1rem, 4vw, 3.5rem);
      border-bottom: 1px solid rgba(29,40,38,.14);
      background: rgba(245,240,230,.94);
      backdrop-filter: blur(12px);
    }}
    .brand {{
      font-family: var(--serif);
      font-weight: 700;
      letter-spacing: .02em;
      white-space: nowrap;
    }}
    nav {{
      display: flex;
      gap: .25rem;
      overflow-x: auto;
      scrollbar-width: thin;
    }}
    nav a {{
      display: inline-block;
      flex: 0 0 auto;
      padding: .35rem .58rem;
      color: var(--muted);
      font-size: .83rem;
      text-decoration: none;
      border-radius: 999px;
    }}
    nav a:hover, nav a:focus-visible {{ background: var(--teal-soft); color: var(--teal-dark); }}
    main {{ width: min(1180px, calc(100% - 2rem)); margin: 0 auto; }}
    .hero {{
      position: relative;
      display: grid;
      grid-template-columns: minmax(0, 1.45fr) minmax(270px, .55fr);
      gap: clamp(2rem, 6vw, 6rem);
      align-items: end;
      padding: clamp(4rem, 10vw, 8.5rem) 0 3.5rem;
    }}
    .hero::before {{
      content: "";
      position: absolute;
      width: 220px;
      height: 220px;
      right: 5%;
      top: 14%;
      border: 1px solid rgba(23,107,101,.22);
      border-radius: 50%;
      box-shadow: 45px 25px 0 -1px var(--paper), 45px 25px 0 rgba(180,77,58,.18);
      pointer-events: none;
    }}
    .eyebrow {{
      margin: 0 0 .8rem;
      color: var(--teal-dark);
      font-size: .78rem;
      font-weight: 800;
      letter-spacing: .14em;
      text-transform: uppercase;
    }}
    h1, h2, h3 {{ font-family: var(--serif); line-height: 1.12; }}
    h1 {{
      max-width: 900px;
      margin: 0;
      font-size: clamp(2.35rem, 6vw, 5.8rem);
      letter-spacing: -.045em;
      text-wrap: balance;
    }}
    .lede {{
      max-width: 760px;
      margin: 1.35rem 0 0;
      color: var(--muted);
      font-size: clamp(1.03rem, 1.5vw, 1.3rem);
    }}
    .hero-status {{
      position: relative;
      z-index: 1;
      padding: 1.2rem;
      border-left: 4px solid var(--teal);
      background: rgba(255,253,250,.8);
      box-shadow: var(--shadow);
    }}
    .hero-status strong {{ display: block; font-family: var(--serif); font-size: 1.2rem; }}
    .hero-status p {{ margin: .35rem 0 0; color: var(--muted); font-size: .92rem; }}
    .claim-strip {{
      margin: 0 0 2rem;
      padding: 1rem 1.2rem;
      border: 1px solid var(--rust);
      background: var(--rust-soft);
      color: #743124;
    }}
    .claim-strip strong {{ margin-right: .45rem; }}
    .metric-grid {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 1px;
      margin: 0 0 4rem;
      background: var(--line);
      border: 1px solid var(--line);
    }}
    .metric {{
      min-width: 0;
      padding: 1.35rem;
      background: var(--white);
    }}
    .metric-value {{
      display: block;
      font-family: var(--serif);
      font-size: clamp(1.75rem, 3vw, 2.7rem);
      line-height: 1;
      font-variant-numeric: tabular-nums;
    }}
    .metric-label {{ display: block; margin-top: .55rem; color: var(--muted); font-size: .86rem; }}
    section {{
      scroll-margin-top: 75px;
      padding: clamp(3.4rem, 7vw, 6rem) 0;
      border-top: 1px solid var(--line);
    }}
    .section-head {{
      display: grid;
      grid-template-columns: minmax(220px, .72fr) minmax(0, 1.28fr);
      gap: clamp(1.4rem, 5vw, 5rem);
      margin-bottom: 2.3rem;
      align-items: start;
    }}
    h2 {{ margin: 0; font-size: clamp(2rem, 4vw, 3.8rem); letter-spacing: -.035em; }}
    .section-copy {{ margin: .2rem 0 0; color: var(--muted); font-size: 1.02rem; }}
    .assertion {{
      border-left: 5px solid var(--teal);
      padding: .9rem 1.1rem;
      background: var(--teal-soft);
      font-weight: 700;
    }}
    .warning {{
      border-left-color: var(--rust);
      background: var(--rust-soft);
      color: #6f3025;
    }}
    .chart-frame {{
      padding: clamp(.7rem, 2vw, 1.4rem);
      border: 1px solid var(--line);
      background: rgba(255,253,250,.86);
      box-shadow: var(--shadow);
      overflow: hidden;
    }}
    .chart-svg {{ display: block; width: 100%; height: auto; overflow: visible; }}
    .chart-svg text {{ fill: var(--ink); font: 14px var(--sans); }}
    .chart-svg .chart-row-label {{ font-weight: 700; }}
    .chart-svg .chart-total, .chart-svg .chart-small {{ fill: var(--muted); font-size: 12px; }}
    .chart-svg .flow-box rect {{ fill: var(--white); stroke: var(--line); stroke-width: 1.5; }}
    .chart-svg .flow-box.muted rect {{ fill: var(--paper-deep); }}
    .chart-svg .flow-box.accent rect {{ fill: var(--teal-soft); stroke: var(--teal); }}
    .chart-svg .flow-number {{ font: 700 25px var(--serif); }}
    .chart-svg .flow-note {{ fill: var(--muted); font-size: 12px; }}
    .chart-svg .rust {{ fill: var(--rust); }}
    .chart-svg .flow-final {{ font-weight: 750; }}
    .chart-caption {{
      display: flex;
      justify-content: space-between;
      gap: 1rem;
      margin: .9rem .2rem 0;
      color: var(--muted);
      font-size: .83rem;
    }}
    .definition-grid {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 1rem;
      margin: 2rem 0;
    }}
    .definition-card {{
      min-width: 0;
      padding: 1.2rem;
      border: 1px solid var(--line);
      background: var(--white);
    }}
    .definition-card h3 {{ margin: .8rem 0 .25rem; font-size: 1.2rem; }}
    .definition-card p {{ margin: .4rem 0 0; color: var(--muted); font-size: .9rem; }}
    .definition-icon {{
      display: grid;
      place-items: center;
      width: 44px;
      height: 44px;
      border-radius: 50%;
      background: var(--teal-soft);
      color: var(--teal-dark);
      font: 700 1.2rem var(--serif);
    }}
    .definition-card .definition-count {{
      color: var(--ink);
      font: 700 1.7rem var(--serif);
    }}
    .sample-filters {{
      display: flex;
      flex-wrap: wrap;
      gap: .5rem;
      margin: 0 0 1rem;
    }}
    .sample-filters button {{
      appearance: none;
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: .45rem .75rem;
      color: var(--ink);
      background: var(--white);
      font: 650 .82rem var(--sans);
      cursor: pointer;
    }}
    .sample-filters button:hover,
    .sample-filters button:focus-visible,
    .sample-filters button[aria-pressed="true"] {{
      border-color: var(--teal);
      background: var(--teal);
      color: white;
    }}
    .table-wrap {{
      width: 100%;
      overflow-x: auto;
      border: 1px solid var(--line);
      background: var(--white);
      box-shadow: var(--shadow);
    }}
    table {{
      width: 100%;
      min-width: 980px;
      border-collapse: collapse;
      font-variant-numeric: tabular-nums;
    }}
    caption {{
      padding: 1rem 1.1rem;
      text-align: left;
      color: var(--muted);
      font-size: .84rem;
    }}
    th, td {{
      padding: .8rem .72rem;
      border-bottom: 1px solid var(--line);
      text-align: right;
      vertical-align: top;
    }}
    thead th {{
      position: sticky;
      top: 0;
      z-index: 1;
      color: var(--muted);
      background: var(--paper-deep);
      font-size: .76rem;
      line-height: 1.3;
    }}
    th:first-child, td:first-child {{ text-align: left; }}
    tbody tr:last-child th, tbody tr:last-child td {{ border-bottom: 0; }}
    tbody tr {{ transition: background .16s ease, box-shadow .16s ease; }}
    tbody tr.is-focus, tbody tr[data-highlighted="true"] {{
      background: var(--gold-soft);
      box-shadow: inset 5px 0 var(--gold);
    }}
    .sample-name {{ display: block; white-space: nowrap; }}
    .subvalue {{ display: block; color: var(--muted); font-size: .76rem; }}
    .tech-tag {{
      display: inline-block;
      margin-top: .25rem;
      padding: .1rem .35rem;
      border-radius: 4px;
      color: #6f3025;
      background: var(--rust-soft);
      font-size: .68rem;
      font-weight: 700;
    }}
    .decision-grid, .readiness-grid, .runtime-grid {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 1rem;
    }}
    .decision-card, .readiness-card, .runtime-card {{
      min-width: 0;
      padding: 1.25rem;
      border-top: 4px solid var(--teal);
      background: var(--white);
      box-shadow: var(--shadow);
    }}
    .decision-card.warn, .readiness-card.warn {{ border-top-color: var(--rust); }}
    .decision-card h3, .readiness-card h3, .runtime-card h3 {{
      margin: 0;
      font-size: 1.12rem;
      overflow-wrap: anywhere;
    }}
    .decision-card .big, .readiness-card .big, .runtime-card .big {{
      display: block;
      margin: .65rem 0;
      font: 700 clamp(1.6rem, 3vw, 2.5rem) var(--serif);
    }}
    .decision-card p, .readiness-card p, .runtime-card p {{
      margin: .5rem 0 0;
      color: var(--muted);
      font-size: .9rem;
    }}
    .status-pill {{
      display: inline-flex;
      align-items: center;
      gap: .4rem;
      padding: .3rem .6rem;
      border-radius: 999px;
      background: var(--rust-soft);
      color: #743124;
      font: 800 .76rem ui-monospace, monospace;
    }}
    .status-pill.good {{ background: var(--teal-soft); color: var(--teal-dark); }}
    ul.compact {{ margin: .8rem 0 0; padding-left: 1.15rem; color: var(--muted); }}
    ul.compact li + li {{ margin-top: .35rem; }}
    details {{
      margin-top: 1rem;
      border: 1px solid var(--line);
      background: rgba(255,253,250,.75);
    }}
    summary {{
      cursor: pointer;
      padding: .9rem 1rem;
      font-weight: 750;
    }}
    details > div {{ padding: 0 1rem 1rem; }}
    .provenance-table table {{ min-width: 1050px; font-size: .8rem; }}
    .path-cell {{ max-width: 540px; overflow-wrap: anywhere; text-align: left; }}
    .checks {{
      display: flex;
      flex-wrap: wrap;
      gap: .45rem;
      margin-top: 1rem;
    }}
    .check {{
      padding: .26rem .52rem;
      border: 1px solid #a8c8c3;
      border-radius: 999px;
      color: var(--teal-dark);
      background: var(--teal-soft);
      font-size: .74rem;
    }}
    footer {{
      margin-top: 4rem;
      padding: 2.5rem max(1rem, calc((100vw - 1180px)/2));
      color: #dbe7e4;
      background: var(--ink);
    }}
    footer p {{ margin: .3rem 0; max-width: 900px; }}
    .no-js-note {{ display: block; }}
    @media (max-width: 900px) {{
      .hero, .section-head {{ grid-template-columns: 1fr; }}
      .hero-status {{ max-width: 600px; }}
      .metric-grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
      .definition-grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
      .decision-grid, .readiness-grid, .runtime-grid {{ grid-template-columns: 1fr; }}
    }}
    @media (max-width: 560px) {{
      main {{ width: min(100% - 1rem, 1180px); }}
      .topbar {{ align-items: flex-start; padding-inline: .75rem; }}
      .brand {{ font-size: .85rem; padding-top: .25rem; }}
      .hero {{ padding-top: 3.2rem; }}
      .hero::before {{
        width: 130px;
        height: 130px;
        right: 12%;
        box-shadow: 20px 14px 0 -1px var(--paper), 20px 14px 0 rgba(180,77,58,.18);
      }}
      h1 {{ font-size: clamp(2.15rem, 13vw, 3.4rem); }}
      .metric-grid {{ grid-template-columns: 1fr 1fr; }}
      .metric {{ padding: 1rem .85rem; }}
      .metric-value {{ font-size: 1.65rem; }}
      .definition-grid {{ grid-template-columns: 1fr; }}
      .chart-frame {{ padding: .25rem; }}
      .chart-caption {{ display: block; }}
      .funnel-svg {{ min-width: 760px; }}
      .chart-frame:has(.funnel-svg) {{ overflow-x: auto; }}
    }}
    @media print {{
      @page {{ size: A4; margin: 13mm; }}
      html, body {{ background: white !important; }}
      body {{ color: black; font-size: 9.5pt; }}
      .topbar, .sample-filters, .skip-link {{ display: none !important; }}
      main {{ width: 100%; }}
      .hero {{ padding: 0 0 12mm; grid-template-columns: 1fr; }}
      .hero::before {{ display: none; }}
      h1 {{ font-size: 30pt; }}
      section {{ padding: 12mm 0; break-before: page; }}
      .metric-grid, .definition-grid {{ break-inside: avoid; }}
      .metric, .definition-card, .decision-card, .readiness-card, .runtime-card,
      .chart-frame, .table-wrap {{ box-shadow: none; break-inside: avoid; }}
      details {{ break-inside: auto; }}
      details > div {{ display: block !important; }}
      .table-wrap {{ overflow: visible; }}
      table {{ min-width: 0; font-size: 7pt; }}
      th, td {{ padding: 4pt 3pt; }}
      .provenance-table {{ display: none; }}
      footer {{ margin-top: 0; padding: 10mm 0; background: white; color: black; }}
    }}
  </style>
</head>
<body data-report-status="{REPORT_STATUS}">
<a class="skip-link" href="#main">跳到主要內容</a>
<header class="topbar">
  <div class="brand">InterSubMod · Exact-PS Observatory</div>
  <nav aria-label="報告段落">
    <a href="#overview">結論</a>
    <a href="#funnel">流程</a>
    <a href="#topology">拓撲</a>
    <a href="#samples">樣本</a>
    <a href="#cn-loh">CN/LOH</a>
    <a href="#methyl">甲基</a>
    <a href="#provenance">來源</a>
  </nav>
</header>
<main id="main">
  <section class="hero" id="overview">
    <div>
      <p class="eyebrow">All 7 technical datasets · chr1–22 · as of {esc(data["as_of_date"])}</p>
      <h1>從 read 共現，走到能說多遠的區域樹。</h1>
      <p class="lede">這份觀察把 strict PS×HP read-linked blocks、最小 mutation-state 候選族、
      family-specific read-AF 排序與 exact topology census 接成一條可稽核的證據鏈；所有無法辨識之處保留 abstain。</p>
    </div>
    <aside class="hero-status" aria-label="報告狀態">
      <strong>VALIDATED DERIVED OBSERVATION</strong>
      <p>13/13 權威輸入 SHA 驗證；19 個分母列守恆。此 HTML 是衍生觀察，不取代 canonical JSON／receipts。</p>
    </aside>
  </section>

  <div class="claim-strip"><strong>Claim ceiling</strong>
  {esc(data["claim_boundary"]["canonical_sentence"])}</div>

  <div class="metric-grid" aria-label="核心數字">
    <article class="metric" data-metric="final_groups" data-value="{data["headline"]["final_groups"]}"><span class="metric-value">{fmt_int(data["headline"]["final_groups"])}</span><span class="metric-label">final topology groups</span></article>
    <article class="metric" data-metric="ranked" data-value="{data["headline"]["ranked"]}"><span class="metric-value">{fmt_int(data["headline"]["ranked"])}</span><span class="metric-label">可用 read-AF 排名</span></article>
    <article class="metric" data-metric="one_rooted_unlabeled_topology" data-value="{data["headline"]["one_rooted_unlabeled_topology"]}"><span class="metric-value">{fmt_int(data["headline"]["one_rooted_unlabeled_topology"])}</span><span class="metric-label">第一順位僅一種 rooted-unlabeled topology</span></article>
    <article class="metric" data-metric="resource_abstain" data-value="{data["headline"]["resource_abstain"]}"><span class="metric-value">{fmt_int(data["headline"]["resource_abstain"])}</span><span class="metric-label">resource abstain；不可隱藏</span></article>
  </div>

  <section id="funnel">
    <div class="section-head">
      <div>
        <p class="eyebrow">01 · Denominator discipline</p>
        <h2>這不是一條單純漏斗。</h2>
      </div>
      <div>
        <p class="section-copy">sSNV、component、W 與 final group 是不同資料粒度。
        85,621 個 W 經 bounded partition 可產生 98,955 個 groups；因此數字變大是合法的一對多切割，
        不能解讀成比例增加。</p>
        <p class="assertion">結論：所有比例都必須跟著自己的分母；拓撲結果只以 71,955 個 ranked units 為分母。</p>
      </div>
    </div>
    <figure class="chart-frame">
      {funnel_svg(data)}
      <figcaption class="chart-caption">
        <span>虛線＝粒度改變或一對多切割；實線＝同粒度內的狀態分流。</span>
        <span>來源：strict-linkage READY + all7 summary</span>
      </figcaption>
    </figure>
  </section>

  <section id="topology">
    <div class="section-head">
      <div>
        <p class="eyebrow">02 · Identifiability</p>
        <h2>88.26% 只剩一種 shape，不等於 88.26% 真實生物樹。</h2>
      </div>
      <div>
        <p class="section-copy">在 71,955 個 read-AF 可排名區域中，39,648 個只有一棵最佳樹，
        23,858 個雖有多棵並列樹但共享同一 rooted-unlabeled topology；合計 63,506。
        另外 8,449 個仍跨 topology 並列。</p>
        <p class="assertion">read-AF 是 frozen candidate family 內的條件式排序，
        不是 likelihood、posterior，也尚未做 CN/LOH 校正。</p>
      </div>
    </div>
    <figure class="chart-frame">
      {topology_svg(data)}
      <figcaption class="chart-caption">
        <span>每列分母＝該資料集的 ranked units；HCC1395 與 DORADO 是同一生物樣本的技術資料。</span>
        <span>唯一 {fmt_int(resolution["UNIQUE_TREE"])} · 同拓撲並列 {fmt_int(resolution["TIED_SAME_TOPOLOGY"])} · 跨拓撲 {fmt_int(resolution["TIED_CROSS_TOPOLOGY"])}</span>
      </figcaption>
    </figure>

    <div class="definition-grid" aria-label="粗拓撲定義">
      {coarse_cards(data)}
    </div>
    <p class="assertion warning">粗分類的分母是 coarse-resolved {fmt_int(coarse_resolved)} 個 ranked units。
    hidden nodes 會參與 depth 與 branching 判斷；分類描述 graph geometry，不能直接等同「單 clone／多 clone／subclone」。</p>

    <figure class="chart-frame" style="margin-top:2rem">
      {active_k_svg(data)}
      <figcaption class="chart-caption">
        <span>淺色＝98,955 全部 groups；深色＝71,955 ranked。active k=0 不建 mutation tree。</span>
        <span>k=9–12 仍有 groups，但 canonical solver 在大 active k 大量 fail-closed abstain。</span>
      </figcaption>
    </figure>
  </section>

  <section id="samples">
    <div class="section-head">
      <div>
        <p class="eyebrow">03 · Cohort table</p>
        <h2>七個技術資料集，同一個分母規格。</h2>
      </div>
      <p class="section-copy">點選樣本只做視覺聚焦，不改變數據或分母。表中「拓撲單一」＝唯一樹＋同拓撲並列；
      methyl 欄為 formal units／robust associations，並非拿來重排樹。</p>
    </div>
    <div class="sample-filters" role="group" aria-label="聚焦樣本">
      <button type="button" data-sample-filter="ALL" aria-pressed="true">全部</button>
      {"".join(f'<button type="button" data-sample-filter="{esc(sample)}" aria-pressed="false">{esc(sample)}</button>' for sample in data["scope"]["technical_datasets"])}
    </div>
    <div class="table-wrap">
      <table id="sample-table">
        <caption>Exact-PS × HP 全 7 技術資料集；數字為區域／group 計數，非細胞數。</caption>
        <thead>
          <tr>
            <th scope="col">資料集</th>
            <th scope="col">W</th>
            <th scope="col">final groups</th>
            <th scope="col">resource<br>abstain</th>
            <th scope="col">ranked</th>
            <th scope="col">唯一樹</th>
            <th scope="col">同拓撲<br>並列</th>
            <th scope="col">跨拓撲<br>未解</th>
            <th scope="col">拓撲單一</th>
            <th scope="col">methyl<br>formal / robust</th>
          </tr>
        </thead>
        <tbody>{sample_table_rows(data)}</tbody>
      </table>
    </div>
    <noscript><p class="no-js-note">JavaScript 已關閉；完整七列仍保留，只有樣本聚焦按鈕不作用。</p></noscript>
  </section>

  <section id="tree-decision">
    <div class="section-head">
      <div>
        <p class="eyebrow">04 · Output policy</p>
        <h2>多棵候選，不應被壓成一棵「真樹」。</h2>
      </div>
      <div>
        <p class="section-copy">目前可安全發布的是候選集合、共同 topology 或 deterministic display representative。
        完整 tree_decision enum 已定義，但 frozen row output 尚未逐列物化，所以報告只列可由權威統計直接重建的桶。</p>
      </div>
    </div>
    <div class="decision-grid">
      <article class="decision-card">
        <h3>RAW_AF_UNIQUE_REPRESENTATIVE</h3>
        <span class="big">39,648</span>
        <p>read-AF 分數唯一第一名；仍是 CN/LOH 未校正的局部數學候選，不是 confirmed biological tree。</p>
      </article>
      <article class="decision-card">
        <h3>TOPOLOGY_REPRESENTATIVE_ONLY</h3>
        <span class="big">23,858</span>
        <p>多棵最佳 labeled trees，但只剩一種 rooted-unlabeled shape；可畫代表樹，必須標明非唯一。</p>
      </article>
      <article class="decision-card warn">
        <h3>CANDIDATE_SET_ONLY</h3>
        <span class="big">8,449</span>
        <p>AF 最高分仍跨 topology；保留完整候選或 consensus backbone，不任選一棵。</p>
      </article>
      <article class="decision-card warn">
        <h3>ABSTAIN_INCOMPLETE</h3>
        <span class="big">10,717</span>
        <p>search-node guard 前未完整列完 minimum family；不得發布 winner 或 topology。</p>
      </article>
      <article class="decision-card warn">
        <h3>AF_ABSTAIN</h3>
        <span class="big">3,224 + 45</span>
        <p>分別為 AF denominator 缺失與 recurrence screen；結構證據可保留，AF 方向不可下結論。</p>
      </article>
      <article class="decision-card">
        <h3>NO_MUTATION_TREE</h3>
        <span class="big">13,014</span>
        <p>active k=0 的 root-only 技術狀態；不能稱為「一棵單突變樹」或「單 clone」。</p>
      </article>
    </div>
  </section>

  <section id="cn-loh">
    <div class="section-head">
      <div>
        <p class="eyebrow">05 · Biological applicability</p>
        <h2>CN/LOH 是目前最重要的未整合層。</h2>
      </div>
      <div>
        <span class="status-pill">NOT_INTEGRATED</span>
        <p class="section-copy">同一 HP 仍可能有 allele-specific amplification、subclonal LOH 或 mutation multiplicity。
        因此 molecule fraction 不一定等於 cell fraction；read-state 也不必然代表不同細胞群。</p>
      </div>
    </div>
    <div class="readiness-grid">
      <article class="readiness-card warn">
        <h3>現在可以說</h3>
        <span class="big">local candidate</span>
        <p>在 exact PS×HP read evidence 與 frozen mutation model 下相容的局部 mutation-state 候選樹。</p>
      </article>
      <article class="readiness-card warn">
        <h3>現在不能說</h3>
        <span class="big">cellular truth</span>
        <p>不能把 AF winner 直接稱為真實 clone/subclone、CCF、全樣本演化樹或確定祖先—後代。</p>
      </article>
      <article class="readiness-card">
        <h3>升格前需要</h3>
        <span class="big">CN + LOH + purity</span>
        <p>allele-specific copy number、LOH、mutation multiplicity、purity 與不確定性／likelihood 校準。</p>
      </article>
    </div>
  </section>

  <section id="methyl">
    <div class="section-head">
      <div>
        <p class="eyebrow">06 · Epigenetic sidecar</p>
        <h2>甲基提供關聯訊號，尚未改寫 topology。</h2>
      </div>
      <div>
        <span class="status-pill good">association-only</span>
        <p class="section-copy">1,045 個 formal units 中 811 個可評估；3 個通過 robust association gates。
        這表示 methylation 與已知 genetic-state label 有穩健關聯，不等於一個 node 內已確認兩個 clone，
        也不能判斷平行或垂直演化。</p>
      </div>
    </div>
    <div class="metric-grid">
      <article class="metric"><span class="metric-value">1,045</span><span class="metric-label">formal units</span></article>
      <article class="metric"><span class="metric-value">811</span><span class="metric-label">evaluable（77.61%）</span></article>
      <article class="metric"><span class="metric-value">3</span><span class="metric-label">robust associations（0.37% evaluable）</span></article>
      <article class="metric"><span class="metric-value">181</span><span class="metric-label">confounded；不升格</span></article>
    </div>
    <p class="assertion warning">目前沒有 same-node unsupervised methyl cluster estimator，
    也沒有 methyl edge scorer 進入 canonical AF 排序；任何用甲基拆 node 或改 tree 的結果都屬新方法，需另立模擬與驗證契約。</p>
  </section>

  <section id="runtime">
    <div class="section-head">
      <div>
        <p class="eyebrow">07 · Runtime</p>
        <h2>三個時間，三種計時範圍。</h2>
      </div>
      <p class="section-copy">{esc(runtime["clock_warning"])}
      這些數字用來定位瓶頸，不能相加或直接互算「solver 占整體幾成」。</p>
    </div>
    <div class="runtime-grid">
      <article class="runtime-card"><h3>cohort runner wall</h3>
        <span class="big">{runtime["cohort_runner_wall_seconds"] / 60:.2f} min</span>
        <p>從 runner 外層量到的實際 wall clock：{runtime["cohort_runner_wall_seconds"]:.3f} s。</p></article>
      <article class="runtime-card"><h3>sample stage sum</h3>
        <span class="big">{runtime["sample_stage_sum_seconds"] / 60:.2f} min</span>
        <p>七個樣本 adapter＋C++ stage 個別 wall time 的加總：{runtime["sample_stage_sum_seconds"]:.3f} s。</p></article>
      <article class="runtime-card"><h3>solver internal</h3>
        <span class="big">{runtime["solver_internal_seconds"]:.2f} s</span>
        <p>solver 內部 microseconds 加總；不含 adapter、IO、process orchestration。</p></article>
    </div>
  </section>

  <section id="provenance">
    <div class="section-head">
      <div>
        <p class="eyebrow">08 · Provenance</p>
        <h2>報告可重生，但不建立第二個真相來源。</h2>
      </div>
      <div>
        <p class="section-copy">Python builder 先核對 authority manifest 中 13 個 path＋SHA，
        再核對 denominator registry、nested strict-linkage data 與 cohort 守恆，最後才輸出 HTML。
        canonical JSON／receipts 仍是數值權威。</p>
        <div class="checks">
          {"".join(f'<span class="check">{esc(check["check_id"])} · PASS</span>' for check in data["conservation_checks"])}
        </div>
      </div>
    </div>
    <details>
      <summary>顯示 13 個權威 artifact 身分</summary>
      <div class="table-wrap provenance-table">
        <table>
          <caption>每一筆都在 build 當下重新計算 SHA-256；路徑只作追溯，不由 HTML 修改。</caption>
          <thead><tr><th scope="col">artifact ID</th><th scope="col">status</th><th scope="col">SHA-256</th><th scope="col">absolute path</th></tr></thead>
          <tbody>{provenance_rows(data)}</tbody>
        </table>
      </div>
    </details>
    <details>
      <summary>方法與宣稱邊界</summary>
      <div>
        <p><strong>Model：</strong>{esc(data["model_contract"]["model_name"])}</p>
        <p><strong>Mutation model：</strong><code>{esc(data["model_contract"]["mutation_model"])}</code></p>
        <p><strong>Forbidden claims：</strong></p>
        <ul class="compact">{"".join(f"<li>{esc(item)}</li>" for item in data["forbidden_claims"])}</ul>
      </div>
    </details>
  </section>
</main>
<footer>
  <p><strong>InterSubMod Exact-PS Observatory</strong> · generated {esc(data["created_at_utc"])}</p>
  <p>Derived observation only. Canonical data were read and hash-verified, not modified.</p>
</footer>
<script id="report-data" type="application/json">{embedded}</script>
<script>
(() => {{
  "use strict";
  const buttons = Array.from(document.querySelectorAll("[data-sample-filter]"));
  const rows = Array.from(document.querySelectorAll("#sample-table tbody tr[data-sample]"));
  buttons.forEach((button) => {{
    button.addEventListener("click", () => {{
      const target = button.dataset.sampleFilter;
      buttons.forEach((candidate) => {{
        candidate.setAttribute("aria-pressed", candidate === button ? "true" : "false");
      }});
      rows.forEach((row) => {{
        const active = target !== "ALL" && row.dataset.sample === target;
        row.classList.toggle("is-focus", active);
        row.dataset.highlighted = active ? "true" : "false";
      }});
    }});
  }});
}})();
</script>
</body>
</html>
"""


def validate_rendered_html(document: str, data: dict[str, Any]) -> list[dict[str, Any]]:
    checks = [
        ("one_h1", document.count("<h1>") == 1),
        ("offline_no_http", "http://" not in document and "https://" not in document),
        ("status_marker", f'data-report-status="{REPORT_STATUS}"' in document),
        ("sample_rows", document.count("<tr data-sample=") == 7),
        ("svg_count", document.count('class="chart-svg') >= 3),
        ("embedded_schema", SCHEMA_NAME in document),
        ("cn_status", "NOT_INTEGRATED" in document),
        ("methyl_status", "association-only" in document),
        ("ranked_sentinel", "71,955" in document),
        ("cross_topology_sentinel", "8,449" in document),
        (
            "all_data_conservation_pass",
            all(check["status"] == "PASS" for check in data["conservation_checks"]),
        ),
    ]
    normalized = [
        {"check_id": check_id, "status": "PASS" if passed else "FAIL"}
        for check_id, passed in checks
    ]
    failed = [check["check_id"] for check in normalized if check["status"] != "PASS"]
    require(not failed, f"render contract failed: {failed}")
    return normalized


def prepare_output_dir(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.exists():
        require(resolved.is_dir(), f"output path is not a directory: {resolved}")
        require(not any(resolved.iterdir()), f"output directory is not empty: {resolved}")
    else:
        resolved.mkdir(parents=True, exist_ok=False)
    return resolved


def main() -> int:
    args = parse_args()
    manifest_path = args.authority_manifest.expanduser().resolve()
    registry_path = args.denominator_registry.expanduser().resolve()
    require(manifest_path.is_file(), f"missing authority manifest: {manifest_path}")
    require(registry_path.is_file(), f"missing denominator registry: {registry_path}")
    data, context = build_report_data(manifest_path, registry_path)
    document = render_html(data)
    render_checks = validate_rendered_html(document, data)
    output_dir = prepare_output_dir(args.output_dir)
    report_data_path = output_dir / "report_data.json"
    report_html_path = output_dir / "report.standalone.html"
    write_json_exclusive(report_data_path, data)
    write_text_exclusive(report_html_path, document)
    data_sidecar = write_sidecar_exclusive(report_data_path)
    html_sidecar = write_sidecar_exclusive(report_html_path)
    receipt = {
        "schema_name": BUILD_RECEIPT_SCHEMA,
        "schema_version": BUILD_RECEIPT_VERSION,
        "created_at_utc": utc_now(),
        "inputs": {
            "authority_manifest": identity(manifest_path),
            "denominator_registry": identity(registry_path),
            "verified_authority_artifacts": context["verified_artifacts"],
            "strict_nested_bundle_identities": data["provenance"][
                "strict_nested_bundle_identities"
            ],
        },
        "outputs": {
            "report_data": identity(report_data_path),
            "report_html": identity(report_html_path),
            "report_data_sidecar": identity(data_sidecar),
            "report_html_sidecar": identity(html_sidecar),
        },
        "checks": {
            "authority_hash_count": len(context["verified_artifacts"]),
            "strict_nested_bundle_hash_count": data["provenance"][
                "strict_nested_bundle_count"
            ],
            "denominator_row_count": len(context["registry"]),
            "conservation": data["conservation_checks"],
            "render": render_checks,
            "canonical_data_unmodified": True,
            "cn_loh_status": "NOT_INTEGRATED",
            "methyl_status": "association-only",
        },
    }
    receipt["all_pass"] = (
        receipt["checks"]["authority_hash_count"] == 13
        and receipt["checks"]["strict_nested_bundle_hash_count"] == 9
        and receipt["checks"]["denominator_row_count"] == 19
        and all(
            check["status"] == "PASS"
            for check in receipt["checks"]["conservation"]
        )
        and all(
            check["status"] == "PASS" for check in receipt["checks"]["render"]
        )
        and receipt["checks"]["canonical_data_unmodified"] is True
        and receipt["checks"]["cn_loh_status"] == "NOT_INTEGRATED"
        and receipt["checks"]["methyl_status"] == "association-only"
    )
    require(receipt["all_pass"], "build receipt checks are not all-pass")
    receipt_path = output_dir / "report_build_receipt.json"
    write_json_exclusive(receipt_path, receipt)
    receipt_sidecar = write_sidecar_exclusive(receipt_path)
    print(f"[INPUT] authority_manifest={manifest_path}")
    print(f"[INPUT] denominator_registry={registry_path}")
    print(f"[OUTPUT] report_data={report_data_path}")
    print(f"[OUTPUT] html={report_html_path}")
    print(f"[OUTPUT] receipt={receipt_path}")
    print(f"[OUTPUT] receipt_sidecar={receipt_sidecar}")
    print(
        json.dumps(
            {
                "all_pass": True,
                "authority_hashes": len(context["verified_artifacts"]),
                "denominator_rows": len(context["registry"]),
                "final_groups": data["headline"]["final_groups"],
                "ranked": data["headline"]["ranked"],
                "one_topology": data["headline"]["one_rooted_unlabeled_topology"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BuildError as exc:
        print(f"[FAIL] {exc}", file=os.sys.stderr)
        raise SystemExit(2)
