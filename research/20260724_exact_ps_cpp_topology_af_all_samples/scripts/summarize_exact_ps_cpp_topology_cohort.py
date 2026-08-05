#!/usr/bin/env python3
"""Build a fail-closed compact summary of completed exact-PS topology samples.

The default mode requires a completed ``cohort_receipt.json`` and summarizes
exactly the samples declared there.  A deliberately explicit partial mode is
available for pilot inspection:

    --samples HCC1395,HCC1395_DORADO --allow-partial

A sample is never opened until its completed ``pipeline_receipt.json`` exists,
passes its SHA256 sidecar, and declares a technical PASS.  Family/resource
ABSTAINs are retained as results and are kept separate from technical validity.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import io
import json
import os
from pathlib import Path
import re
import sys
from typing import Any, Iterable, Mapping, Sequence


SUMMARY_SCHEMA = "intersubmod.exact_ps_cpp_topology_af.cohort_summary"
SUMMARY_VERSION = "1.0.0"
PIPELINE_SCHEMA = "intersubmod.exact_ps_cpp_topology_cohort.sample_pipeline_receipt"
COHORT_SCHEMA = "intersubmod.exact_ps_cpp_topology_cohort.cohort_receipt"
MLHP_RECEIPT_SCHEMA = "intersubmod.exact_ps_layered_tree_input.receipt"
TOPOLOGY_RECEIPT_SCHEMA = "intersubmod.exact_ps_cpp_topology_af.receipt"
TOPOLOGY_UNIT_SCHEMA = "intersubmod.exact_ps_cpp_topology_af.unit"
SOURCE_SCHEMA_VERSION = "1.0.0"
SAMPLE_RE = re.compile(r"^[A-Za-z0-9_.-]+$")

TSV_FIELDS = (
    "sample",
    "groups_total",
    "mutation_bearing_units",
    "mutation_family_complete_units",
    "mutation_family_abstain_units",
    "objective_certified_units",
    "objective_abstain_units",
    "ranked_units",
    "zero_denominator_units",
    "recurrence_required_units",
    "no_active_alt_units",
    "best_tree_unique_units",
    "best_tree_tied_units",
    "morphology_single_only_unique",
    "morphology_sister_only_unique",
    "morphology_direct_only_unique",
    "morphology_sister_plus_direct_unique",
    "morphology_unresolved_tied",
    "active_k_distribution_json",
    "search_nodes_sum",
    "search_nodes_max",
    "search_nodes_at_or_above_guard_units",
    "solver_elapsed_microseconds_sum",
    "solver_elapsed_microseconds_max",
    "adapter_wall_seconds",
    "cpp_topology_wall_seconds",
    "pipeline_stage_wall_seconds",
    "topology_runtime_seconds",
    "all_units_objective_certified",
    "all_mutation_bearing_families_complete",
)


class SummaryError(RuntimeError):
    """Raised when an input identity or conservation check fails."""


@dataclass(frozen=True)
class FileIdentity:
    path: str
    size_bytes: int
    sha256: str

    def as_json(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        }


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SummaryError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> FileIdentity:
    require(path.is_file(), f"required file is missing: {path}")
    resolved = path.resolve(strict=True)
    return FileIdentity(
        path=str(resolved),
        size_bytes=resolved.stat().st_size,
        sha256=sha256_path(resolved),
    )


def parse_identity(value: Any, label: str) -> FileIdentity:
    require(isinstance(value, Mapping), f"{label} identity is not an object")
    path = value.get("path")
    size = value.get("size_bytes")
    digest = value.get("sha256")
    require(isinstance(path, str) and path, f"{label} identity path is invalid")
    require(isinstance(size, int) and size >= 0, f"{label} identity size is invalid")
    require(
        isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
        f"{label} identity SHA256 is invalid",
    )
    return FileIdentity(path=path, size_bytes=size, sha256=digest)


def require_identity(
    actual: FileIdentity,
    expected_value: Any,
    label: str,
    *,
    expected_path: Path | None = None,
) -> None:
    expected = parse_identity(expected_value, label)
    if expected_path is not None:
        require(
            Path(expected.path).resolve(strict=True) == expected_path.resolve(strict=True),
            f"{label} identity points to an unexpected path",
        )
    require(actual.size_bytes == expected.size_bytes, f"{label} size mismatch")
    require(actual.sha256 == expected.sha256, f"{label} SHA256 mismatch")


def read_json(path: Path, label: str) -> Mapping[str, Any]:
    require(path.is_file(), f"{label} is missing: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SummaryError(f"cannot read {label}: {path}: {exc}") from exc
    require(isinstance(value, Mapping), f"{label} root is not an object: {path}")
    return value


def verify_sha256_sidecar(path: Path, label: str) -> None:
    sidecar = path.with_name(f"{path.name}.sha256")
    require(sidecar.is_file(), f"{label} SHA256 sidecar is missing: {sidecar}")
    fields = sidecar.read_text(encoding="ascii").strip().split()
    require(len(fields) == 2, f"{label} SHA256 sidecar format is invalid")
    require(fields[1] == path.name, f"{label} SHA256 sidecar filename mismatch")
    require(fields[0] == sha256_path(path), f"{label} SHA256 sidecar mismatch")


def schema_is(document: Mapping[str, Any], name: str, label: str) -> None:
    require(document.get("schema_name") == name, f"{label} schema_name mismatch")
    require(
        document.get("schema_version") == SOURCE_SCHEMA_VERSION,
        f"{label} schema_version mismatch",
    )


def expected_sample_paths(sample_root: Path, sample: str) -> dict[str, Path]:
    return {
        "mlhp": sample_root / f"{sample}.exact_ps_mlhp.json",
        "mlhp_receipt": sample_root / f"{sample}.exact_ps_mlhp.json.receipt.json",
        "topology": sample_root / f"{sample}.topology.jsonl",
        "topology_receipt": sample_root / f"{sample}.topology.receipt.json",
    }


def pipeline_bound_identities(
    pipeline: Mapping[str, Any],
    sample_root: Path,
    sample: str,
) -> dict[Path, FileIdentity]:
    values = pipeline.get("bound_artifacts")
    require(isinstance(values, list) and values, "pipeline bound_artifacts is invalid")
    result: dict[Path, FileIdentity] = {}
    for index, value in enumerate(values):
        spec = parse_identity(value, f"pipeline bound_artifacts[{index}]")
        path = Path(spec.path)
        require(path.is_absolute(), "pipeline bound artifact path must be absolute")
        resolved = path.resolve(strict=True)
        require(
            resolved == sample_root or sample_root in resolved.parents,
            f"pipeline bound artifact escapes sample root: {resolved}",
        )
        require(resolved not in result, f"duplicate pipeline bound artifact: {resolved}")
        result[resolved] = spec

    for label, expected in expected_sample_paths(sample_root, sample).items():
        resolved = expected.resolve(strict=True)
        require(resolved in result, f"pipeline does not bind expected {label}: {resolved}")
    return result


def verify_eager_pipeline_artifacts(
    bound: Mapping[Path, FileIdentity],
    topology_path: Path,
) -> list[dict[str, Any]]:
    """Verify every bound file except topology JSONL, which is hashed while parsed."""
    identities: list[dict[str, Any]] = []
    topology_resolved = topology_path.resolve(strict=True)
    for path, expected in sorted(bound.items(), key=lambda item: str(item[0])):
        if path == topology_resolved:
            continue
        actual = identity(path)
        require_identity(actual, expected.as_json(), f"pipeline bound artifact {path}")
        identities.append(actual.as_json())
    return identities


def counter_json(counter: Counter[Any]) -> dict[str, int]:
    return {
        str(key): int(counter[key])
        for key in sorted(counter, key=lambda item: (str(type(item)), str(item)))
    }


def require_census(
    observed: Counter[str],
    expected_value: Any,
    label: str,
    total: int,
) -> None:
    require(isinstance(expected_value, Mapping), f"{label} census is missing")
    expected: dict[str, int] = {}
    for key, value in expected_value.items():
        require(isinstance(key, str), f"{label} census key is invalid")
        require(isinstance(value, int) and value >= 0, f"{label} census value is invalid")
        expected[key] = value
    require(sum(expected.values()) == total, f"{label} census does not conserve total")
    require(dict(observed) == expected, f"{label} JSONL/receipt census mismatch")


def coarse_morphology(value: Any) -> str:
    require(isinstance(value, Mapping), "unique best tree lacks morphology")
    sister = value.get("has_sister")
    direct = value.get("has_direct")
    require(isinstance(sister, bool), "morphology has_sister is invalid")
    require(isinstance(direct, bool), "morphology has_direct is invalid")
    if sister and direct:
        return "sister_plus_direct"
    if sister:
        return "sister_only"
    if direct:
        return "direct_only"
    return "single_only"


def parse_nonnegative_int(value: Any, label: str) -> int:
    require(isinstance(value, int) and not isinstance(value, bool) and value >= 0, f"{label} is invalid")
    return value


def parse_positive_decimal(value: Any, label: str) -> int:
    require(
        isinstance(value, str) and re.fullmatch(r"[1-9][0-9]*", value) is not None,
        f"{label} is not a positive decimal string",
    )
    return int(value)


def summarize_topology_jsonl(
    path: Path,
    *,
    sample: str,
    expected_identity: FileIdentity,
    topology_receipt: Mapping[str, Any],
    max_search_nodes: int,
) -> tuple[dict[str, Any], FileIdentity]:
    digest = hashlib.sha256()
    size = 0
    rows = 0
    active_k: Counter[int] = Counter()
    unit_status: Counter[str] = Counter()
    family_status: Counter[str] = Counter()
    objective_status: Counter[str] = Counter()
    read_af_status: Counter[str] = Counter()
    morphology: Counter[str] = Counter()
    mutation_bearing = 0
    mutation_family_complete = 0
    ranked = 0
    zero_denominator = 0
    recurrence_required = 0
    no_active_alt = 0
    best_unique = 0
    best_tied = 0
    search_nodes_sum = 0
    search_nodes_max = 0
    search_nodes_at_guard = 0
    solver_us_sum = 0
    solver_us_max = 0

    with path.open("rb") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            digest.update(raw_line)
            size += len(raw_line)
            require(raw_line.endswith(b"\n"), f"topology line {line_number} lacks newline")
            try:
                row = json.loads(raw_line)
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                raise SummaryError(f"invalid topology JSONL line {line_number}: {exc}") from exc
            require(isinstance(row, Mapping), f"topology line {line_number} is not an object")
            schema_is(row, TOPOLOGY_UNIT_SCHEMA, f"topology line {line_number}")
            require(row.get("sample") == sample, f"topology line {line_number} sample mismatch")
            k = parse_nonnegative_int(
                row.get("active_bit_count"), f"topology line {line_number} active_bit_count"
            )
            require(k <= 12, f"topology line {line_number} exceeds active k=12")
            active_k[k] += 1
            rows += 1
            if k > 0:
                mutation_bearing += 1

            statuses: list[tuple[str, Counter[str]]] = [
                ("unit_status", unit_status),
                ("family_status", family_status),
                ("objective_status", objective_status),
                ("read_af_status", read_af_status),
            ]
            for field, counter in statuses:
                value = row.get(field)
                require(isinstance(value, str) and value, f"topology line {line_number} {field} is invalid")
                counter[value] += 1

            family = row["family_status"]
            require(
                family == "FAMILY_COMPLETE" or family.startswith("ABSTAIN_"),
                f"topology line {line_number} has unknown family status",
            )
            objective = row["objective_status"]
            require(
                objective == "OBJECTIVE_CERTIFIED" or objective.startswith("ABSTAIN_"),
                f"topology line {line_number} has unknown objective status",
            )
            if k > 0 and family == "FAMILY_COMPLETE":
                mutation_family_complete += 1

            nodes = parse_nonnegative_int(
                row.get("search_nodes"), f"topology line {line_number} search_nodes"
            )
            solver_us = parse_nonnegative_int(
                row.get("solver_elapsed_microseconds"),
                f"topology line {line_number} solver_elapsed_microseconds",
            )
            search_nodes_sum += nodes
            search_nodes_max = max(search_nodes_max, nodes)
            solver_us_sum += solver_us
            solver_us_max = max(solver_us_max, solver_us)
            if max_search_nodes > 0 and nodes >= max_search_nodes:
                search_nodes_at_guard += 1

            status = row["unit_status"]
            if status == "ranked":
                ranked += 1
                require(row["read_af_status"] == "ranked_complete", "ranked/read-AF status mismatch")
                unique = row.get("best_tree_unique")
                require(isinstance(unique, bool), "ranked row best_tree_unique is invalid")
                tie_count = parse_positive_decimal(
                    row.get("best_tree_tie_count"), "ranked row best_tree_tie_count"
                )
                require(unique == (tie_count == 1), "best tree unique/tie count mismatch")
                if unique:
                    best_unique += 1
                    morphology[coarse_morphology(row.get("representative_best_morphology"))] += 1
                else:
                    best_tied += 1
                    morphology["unresolved_tied"] += 1
            else:
                require(
                    row.get("best_tree_unique") is None,
                    f"non-ranked line {line_number} carries best_tree_unique",
                )

            zero_denominator += int(status == "zero_denominator")
            recurrence_required += int(status == "recurrence_required")
            no_active_alt += int(status == "no_active_alt")

    actual_identity = FileIdentity(
        path=str(path.resolve(strict=True)),
        size_bytes=size,
        sha256=digest.hexdigest(),
    )
    require_identity(actual_identity, expected_identity.as_json(), "topology JSONL")
    counts = topology_receipt.get("counts")
    census = topology_receipt.get("status_census")
    require(isinstance(counts, Mapping), "topology receipt counts is invalid")
    require(isinstance(census, Mapping), "topology receipt status_census is invalid")
    expected_total = parse_nonnegative_int(counts.get("total_units"), "receipt total_units")
    require(rows == expected_total, "topology JSONL row count differs from receipt")
    require(
        mutation_bearing
        == parse_nonnegative_int(counts.get("mutation_bearing_units"), "receipt mutation units"),
        "mutation-bearing count mismatch",
    )
    require(
        mutation_family_complete
        == parse_nonnegative_int(
            counts.get("family_complete_mutation_bearing_units"),
            "receipt complete mutation families",
        ),
        "complete mutation-family count mismatch",
    )
    require(
        ranked
        == parse_nonnegative_int(counts.get("ranking_complete_count"), "receipt ranked count"),
        "ranked count mismatch",
    )
    require(
        objective_status["OBJECTIVE_CERTIFIED"]
        == parse_nonnegative_int(
            counts.get("objective_certified_units"), "receipt certified objective count"
        ),
        "objective-certified count mismatch",
    )
    require_census(unit_status, census.get("unit_status"), "unit_status", rows)
    require_census(family_status, census.get("family_status"), "family_status", rows)
    require_census(objective_status, census.get("objective_status"), "objective_status", rows)
    require_census(read_af_status, census.get("read_af_status"), "read_af_status", rows)
    require(best_unique + best_tied == ranked, "unique/tied ranked conservation failed")
    require(
        sum(morphology[key] for key in (
            "single_only",
            "sister_only",
            "direct_only",
            "sister_plus_direct",
        ))
        == best_unique,
        "unique morphology conservation failed",
    )
    require(morphology["unresolved_tied"] == best_tied, "tied morphology conservation failed")
    require(sum(active_k.values()) == rows, "active-k conservation failed")

    family_abstain = mutation_bearing - mutation_family_complete
    objective_certified = objective_status["OBJECTIVE_CERTIFIED"]
    result = {
        "groups_total": rows,
        "mutation_bearing_units": mutation_bearing,
        "mutation_family_complete_units": mutation_family_complete,
        "mutation_family_abstain_units": family_abstain,
        "objective_certified_units": objective_certified,
        "objective_abstain_units": rows - objective_certified,
        "ranked_units": ranked,
        "zero_denominator_units": zero_denominator,
        "recurrence_required_units": recurrence_required,
        "no_active_alt_units": no_active_alt,
        "best_tree_unique_units": best_unique,
        "best_tree_tied_units": best_tied,
        "unique_best_morphology": {
            "single_only": morphology["single_only"],
            "sister_only": morphology["sister_only"],
            "direct_only": morphology["direct_only"],
            "sister_plus_direct": morphology["sister_plus_direct"],
        },
        "morphology_unresolved_tied": morphology["unresolved_tied"],
        "morphology_scope": (
            "Coarse graph geometry is counted only when the globally best tree is unique; "
            "all AF-tied best trees are unresolved and their deterministic representative "
            "is not treated as a topology conclusion."
        ),
        "active_k_distribution": {
            str(k): active_k.get(k, 0) for k in range(0, 13)
        },
        "search_nodes": {
            "sum": search_nodes_sum,
            "max": search_nodes_max,
            "at_or_above_configured_guard_units": search_nodes_at_guard,
        },
        "solver_elapsed_microseconds": {
            "sum": solver_us_sum,
            "max": solver_us_max,
        },
        "status_census": {
            "unit_status": counter_json(unit_status),
            "family_status": counter_json(family_status),
            "objective_status": counter_json(objective_status),
            "read_af_status": counter_json(read_af_status),
        },
    }
    return result, actual_identity


def stage_wall_seconds(pipeline: Mapping[str, Any]) -> dict[str, float]:
    stages = pipeline.get("stages")
    require(isinstance(stages, list), "pipeline stages is invalid")
    values: dict[str, float] = {}
    total = 0.0
    for index, stage in enumerate(stages):
        require(isinstance(stage, Mapping), f"pipeline stage {index} is invalid")
        name = stage.get("stage")
        wall = stage.get("wall_seconds")
        require(isinstance(name, str) and name, f"pipeline stage {index} name is invalid")
        require(
            isinstance(wall, (int, float)) and not isinstance(wall, bool) and wall >= 0,
            f"pipeline stage {index} wall_seconds is invalid",
        )
        require(name not in values, f"duplicate pipeline stage: {name}")
        values[name] = float(wall)
        total += float(wall)
    require("exact_ps_partition_to_mlhp" in values, "adapter stage is missing")
    require("cpp_topology" in values, "C++ topology stage is missing")
    return {
        "exact_ps_partition_to_mlhp": values["exact_ps_partition_to_mlhp"],
        "cpp_topology": values["cpp_topology"],
        "total": total,
    }


def summarize_sample(cohort_root: Path, sample: str) -> dict[str, Any]:
    require(SAMPLE_RE.fullmatch(sample) is not None, f"invalid sample name: {sample}")
    sample_root = cohort_root / "samples" / sample
    require(sample_root.is_dir(), f"sample directory is missing: {sample_root}")

    # Completion gate: do not inspect any downstream artifact before this exists.
    pipeline_path = sample_root / "pipeline_receipt.json"
    require(
        pipeline_path.is_file(),
        f"sample is unfinished (pipeline receipt absent); refusing to read it: {sample}",
    )
    verify_sha256_sidecar(pipeline_path, f"{sample} pipeline receipt")
    pipeline = read_json(pipeline_path, f"{sample} pipeline receipt")
    schema_is(pipeline, PIPELINE_SCHEMA, f"{sample} pipeline receipt")
    require(pipeline.get("all_pass") is True, f"{sample} pipeline technical all_pass is false")
    require(
        pipeline.get("technical_all_pass") is True,
        f"{sample} pipeline technical_all_pass is false",
    )
    require(pipeline.get("sample") == sample, f"{sample} pipeline sample mismatch")

    paths = expected_sample_paths(sample_root, sample)
    bound = pipeline_bound_identities(pipeline, sample_root.resolve(strict=True), sample)
    source_identities = verify_eager_pipeline_artifacts(bound, paths["topology"])

    mlhp_receipt = read_json(paths["mlhp_receipt"], f"{sample} MLHP receipt")
    schema_is(mlhp_receipt, MLHP_RECEIPT_SCHEMA, f"{sample} MLHP receipt")
    require(mlhp_receipt.get("all_pass") is True, f"{sample} MLHP receipt is not PASS")
    mlhp_identity = identity(paths["mlhp"])
    require_identity(
        mlhp_identity,
        mlhp_receipt.get("output"),
        f"{sample} MLHP output",
        expected_path=paths["mlhp"],
    )

    topology_receipt = read_json(paths["topology_receipt"], f"{sample} topology receipt")
    schema_is(topology_receipt, TOPOLOGY_RECEIPT_SCHEMA, f"{sample} topology receipt")
    require(topology_receipt.get("all_pass") is True, f"{sample} topology technical PASS is false")
    require_identity(
        mlhp_identity,
        topology_receipt.get("input"),
        f"{sample} topology input",
        expected_path=paths["mlhp"],
    )
    parameters = topology_receipt.get("parameters")
    require(isinstance(parameters, Mapping), f"{sample} topology parameters is invalid")
    max_search_nodes = parse_nonnegative_int(
        parameters.get("max_search_nodes"), f"{sample} max_search_nodes"
    )
    expected_topology_identity = parse_identity(
        topology_receipt.get("output"), f"{sample} topology output"
    )
    require(
        Path(expected_topology_identity.path).resolve(strict=True)
        == paths["topology"].resolve(strict=True),
        f"{sample} topology output path mismatch",
    )
    bound_topology = bound[paths["topology"].resolve(strict=True)]
    require(
        expected_topology_identity == bound_topology,
        f"{sample} topology receipt/pipeline output identity mismatch",
    )
    summary, topology_identity = summarize_topology_jsonl(
        paths["topology"],
        sample=sample,
        expected_identity=expected_topology_identity,
        topology_receipt=topology_receipt,
        max_search_nodes=max_search_nodes,
    )
    source_identities.append(topology_identity.as_json())

    mlhp_counts = mlhp_receipt.get("counts")
    pipeline_counts = pipeline.get("counts")
    require(isinstance(mlhp_counts, Mapping), f"{sample} MLHP counts is invalid")
    require(isinstance(pipeline_counts, Mapping), f"{sample} pipeline counts is invalid")
    mlhp_groups = parse_nonnegative_int(
        mlhp_counts.get("groups_tree_input"), f"{sample} MLHP groups_tree_input"
    )
    require(mlhp_groups == summary["groups_total"], f"{sample} MLHP/topology groups mismatch")
    require(
        parse_nonnegative_int(pipeline_counts.get("mlhp_groups"), f"{sample} pipeline MLHP groups")
        == summary["groups_total"],
        f"{sample} pipeline MLHP group count mismatch",
    )
    require(
        parse_nonnegative_int(
            pipeline_counts.get("topology_jsonl_rows"), f"{sample} pipeline topology rows"
        )
        == summary["groups_total"],
        f"{sample} pipeline topology row count mismatch",
    )

    all_family_complete = summary["mutation_family_abstain_units"] == 0
    all_objective_certified = summary["objective_abstain_units"] == 0
    require(
        topology_receipt.get("all_mutation_bearing_families_complete")
        is all_family_complete,
        f"{sample} topology family-completeness flag mismatch",
    )
    require(
        topology_receipt.get("all_units_objective_certified")
        is all_objective_certified,
        f"{sample} topology objective-completeness flag mismatch",
    )
    require(
        pipeline.get("all_mutation_bearing_families_complete") is all_family_complete,
        f"{sample} pipeline family-completeness flag mismatch",
    )
    require(
        pipeline.get("all_units_objective_certified") is all_objective_certified,
        f"{sample} pipeline objective-completeness flag mismatch",
    )

    runtime = topology_receipt.get("runtime")
    require(isinstance(runtime, Mapping), f"{sample} topology runtime is invalid")
    elapsed = runtime.get("elapsed_seconds")
    require(
        isinstance(elapsed, (int, float)) and not isinstance(elapsed, bool) and elapsed >= 0,
        f"{sample} topology elapsed_seconds is invalid",
    )
    summary.update(
        {
            "sample": sample,
            "technical_all_pass": True,
            "all_units_objective_certified": all_objective_certified,
            "all_mutation_bearing_families_complete": all_family_complete,
            "configured_guards": {
                "max_family_size": parse_nonnegative_int(
                    parameters.get("max_family_size"), f"{sample} max_family_size"
                ),
                "max_search_nodes": max_search_nodes,
                "max_active_bits": 12,
            },
            "stage_wall_seconds": stage_wall_seconds(pipeline),
            "topology_runtime_seconds": float(elapsed),
            "source_identities": sorted(source_identities, key=lambda item: item["path"]),
            "pipeline_receipt_identity": identity(pipeline_path).as_json(),
        }
    )
    return summary


def validate_cohort_receipt(cohort_root: Path) -> tuple[Mapping[str, Any], list[str], FileIdentity]:
    receipt_path = cohort_root / "cohort_receipt.json"
    require(receipt_path.is_file(), f"completed cohort receipt is missing: {receipt_path}")
    verify_sha256_sidecar(receipt_path, "cohort receipt")
    receipt = read_json(receipt_path, "cohort receipt")
    schema_is(receipt, COHORT_SCHEMA, "cohort receipt")
    require(receipt.get("all_pass") is True, "cohort receipt technical all_pass is false")
    require(receipt.get("technical_all_pass") is True, "cohort receipt technical_all_pass is false")
    scope = receipt.get("scope")
    require(isinstance(scope, Mapping), "cohort scope is invalid")
    samples = scope.get("samples")
    require(isinstance(samples, list) and samples, "cohort sample scope is invalid")
    require(all(isinstance(sample, str) for sample in samples), "cohort sample name is invalid")
    require(len(samples) == len(set(samples)), "cohort sample scope contains duplicates")
    sample_receipts = receipt.get("sample_receipts")
    require(isinstance(sample_receipts, Mapping), "cohort sample_receipts is invalid")
    require(set(sample_receipts) == set(samples), "cohort sample receipt grid mismatch")
    for sample in samples:
        block = sample_receipts[sample]
        require(isinstance(block, Mapping), f"cohort receipt block is invalid: {sample}")
        pipeline_path = cohort_root / "samples" / sample / "pipeline_receipt.json"
        actual = identity(pipeline_path)
        require_identity(
            actual,
            block.get("pipeline"),
            f"cohort pipeline identity {sample}",
            expected_path=pipeline_path,
        )
    return receipt, list(samples), identity(receipt_path)


def sum_sample_field(samples: Iterable[Mapping[str, Any]], field: str) -> int:
    return sum(int(sample[field]) for sample in samples)


def build_cohort_summary(
    cohort_root: Path,
    samples: Sequence[str],
    *,
    partial: bool,
    cohort_receipt: Mapping[str, Any] | None,
    cohort_receipt_identity: FileIdentity | None,
) -> dict[str, Any]:
    require(samples, "no samples were requested")
    require(len(samples) == len(set(samples)), "sample list contains duplicates")
    sample_summaries = [summarize_sample(cohort_root, sample) for sample in samples]

    totals: dict[str, Any] = {}
    additive_fields = (
        "groups_total",
        "mutation_bearing_units",
        "mutation_family_complete_units",
        "mutation_family_abstain_units",
        "objective_certified_units",
        "objective_abstain_units",
        "ranked_units",
        "zero_denominator_units",
        "recurrence_required_units",
        "no_active_alt_units",
        "best_tree_unique_units",
        "best_tree_tied_units",
        "morphology_unresolved_tied",
    )
    for field in additive_fields:
        totals[field] = sum_sample_field(sample_summaries, field)
    totals["unique_best_morphology"] = {
        label: sum(
            sample["unique_best_morphology"][label] for sample in sample_summaries
        )
        for label in ("single_only", "sister_only", "direct_only", "sister_plus_direct")
    }
    totals["active_k_distribution"] = {
        str(k): sum(sample["active_k_distribution"][str(k)] for sample in sample_summaries)
        for k in range(0, 13)
    }
    totals["search_nodes"] = {
        "sum": sum(sample["search_nodes"]["sum"] for sample in sample_summaries),
        "max": max(sample["search_nodes"]["max"] for sample in sample_summaries),
        "at_or_above_configured_guard_units": sum(
            sample["search_nodes"]["at_or_above_configured_guard_units"]
            for sample in sample_summaries
        ),
    }
    totals["solver_elapsed_microseconds"] = {
        "sum": sum(
            sample["solver_elapsed_microseconds"]["sum"] for sample in sample_summaries
        ),
        "max": max(
            sample["solver_elapsed_microseconds"]["max"] for sample in sample_summaries
        ),
    }
    totals["stage_wall_seconds"] = {
        "exact_ps_partition_to_mlhp": sum(
            sample["stage_wall_seconds"]["exact_ps_partition_to_mlhp"]
            for sample in sample_summaries
        ),
        "cpp_topology": sum(
            sample["stage_wall_seconds"]["cpp_topology"] for sample in sample_summaries
        ),
        "total": sum(sample["stage_wall_seconds"]["total"] for sample in sample_summaries),
    }
    totals["topology_runtime_seconds"] = sum(
        sample["topology_runtime_seconds"] for sample in sample_summaries
    )

    checks = {
        "groups_equal_active_k_sum": (
            totals["groups_total"] == sum(totals["active_k_distribution"].values())
        ),
        "mutation_complete_plus_abstain_conserved": (
            totals["mutation_bearing_units"]
            == totals["mutation_family_complete_units"]
            + totals["mutation_family_abstain_units"]
        ),
        "objective_certified_plus_abstain_conserved": (
            totals["groups_total"]
            == totals["objective_certified_units"] + totals["objective_abstain_units"]
        ),
        "ranked_unique_plus_tied_conserved": (
            totals["ranked_units"]
            == totals["best_tree_unique_units"] + totals["best_tree_tied_units"]
        ),
        "unique_morphology_conserved": (
            totals["best_tree_unique_units"]
            == sum(totals["unique_best_morphology"].values())
        ),
        "tied_morphology_is_unresolved": (
            totals["best_tree_tied_units"] == totals["morphology_unresolved_tied"]
        ),
        "all_source_identity_gates_passed": True,
        "only_completed_pipeline_samples_read": True,
    }
    require(all(checks.values()), "cohort conservation check failed")

    all_family_complete = totals["mutation_family_abstain_units"] == 0
    all_objective_certified = totals["objective_abstain_units"] == 0
    if cohort_receipt is not None:
        require(
            cohort_receipt.get("all_mutation_bearing_families_complete")
            is all_family_complete,
            "cohort family-completeness flag mismatch",
        )
        require(
            cohort_receipt.get("all_units_objective_certified")
            is all_objective_certified,
            "cohort objective-completeness flag mismatch",
        )

    scope: dict[str, Any] = {
        "cohort_root": str(cohort_root),
        "samples": list(samples),
        "partial": partial,
    }
    if cohort_receipt is not None:
        scope["runner_scope"] = cohort_receipt.get("scope")
    result = {
        "schema_name": SUMMARY_SCHEMA,
        "schema_version": SUMMARY_VERSION,
        "all_pass": True,
        "technical_all_pass": True,
        "all_mutation_bearing_families_complete": all_family_complete,
        "all_units_objective_certified": all_objective_certified,
        "validation_evidence_eligible": (
            not partial
            and cohort_receipt is not None
            and bool((cohort_receipt.get("scope") or {}).get("canonical_all7"))
            and bool((cohort_receipt.get("scope") or {}).get("autosomes_complete"))
            and all_family_complete
            and all_objective_certified
        ),
        "scope": scope,
        "samples": sample_summaries,
        "totals": totals,
        "checks": checks,
        "cohort_receipt_identity": (
            cohort_receipt_identity.as_json()
            if cohort_receipt_identity is not None
            else None
        ),
        "claim_ceiling": (
            "Technical PASS proves receipt/schema/SHA256 identities and arithmetic "
            "conservation for completed samples. Resource-limit families remain explicit "
            "ABSTAINs. read-AF ranking orders recurrence-allowed minimum graph models; "
            "it is not calibrated proof of a unique biological clone tree. Sister/direct "
            "labels are graph geometry, and are reported only for a unique global AF winner."
        ),
    }
    return result


def sample_to_tsv(sample: Mapping[str, Any]) -> dict[str, Any]:
    morphology = sample["unique_best_morphology"]
    return {
        "sample": sample["sample"],
        "groups_total": sample["groups_total"],
        "mutation_bearing_units": sample["mutation_bearing_units"],
        "mutation_family_complete_units": sample["mutation_family_complete_units"],
        "mutation_family_abstain_units": sample["mutation_family_abstain_units"],
        "objective_certified_units": sample["objective_certified_units"],
        "objective_abstain_units": sample["objective_abstain_units"],
        "ranked_units": sample["ranked_units"],
        "zero_denominator_units": sample["zero_denominator_units"],
        "recurrence_required_units": sample["recurrence_required_units"],
        "no_active_alt_units": sample["no_active_alt_units"],
        "best_tree_unique_units": sample["best_tree_unique_units"],
        "best_tree_tied_units": sample["best_tree_tied_units"],
        "morphology_single_only_unique": morphology["single_only"],
        "morphology_sister_only_unique": morphology["sister_only"],
        "morphology_direct_only_unique": morphology["direct_only"],
        "morphology_sister_plus_direct_unique": morphology["sister_plus_direct"],
        "morphology_unresolved_tied": sample["morphology_unresolved_tied"],
        "active_k_distribution_json": json.dumps(
            sample["active_k_distribution"], sort_keys=True, separators=(",", ":")
        ),
        "search_nodes_sum": sample["search_nodes"]["sum"],
        "search_nodes_max": sample["search_nodes"]["max"],
        "search_nodes_at_or_above_guard_units": sample["search_nodes"][
            "at_or_above_configured_guard_units"
        ],
        "solver_elapsed_microseconds_sum": sample["solver_elapsed_microseconds"]["sum"],
        "solver_elapsed_microseconds_max": sample["solver_elapsed_microseconds"]["max"],
        "adapter_wall_seconds": sample["stage_wall_seconds"][
            "exact_ps_partition_to_mlhp"
        ],
        "cpp_topology_wall_seconds": sample["stage_wall_seconds"]["cpp_topology"],
        "pipeline_stage_wall_seconds": sample["stage_wall_seconds"]["total"],
        "topology_runtime_seconds": sample["topology_runtime_seconds"],
        "all_units_objective_certified": str(
            sample["all_units_objective_certified"]
        ).lower(),
        "all_mutation_bearing_families_complete": str(
            sample["all_mutation_bearing_families_complete"]
        ).lower(),
    }


def tsv_payload(samples: Sequence[Mapping[str, Any]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=TSV_FIELDS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for sample in samples:
        writer.writerow(sample_to_tsv(sample))
    return buffer.getvalue().encode("utf-8")


def write_outputs_exclusive(
    summary: dict[str, Any],
    output_json: Path,
    output_tsv: Path,
) -> None:
    require(output_json != output_tsv, "JSON and TSV output paths must differ")
    require(not output_json.exists(), f"JSON output already exists: {output_json}")
    require(not output_tsv.exists(), f"TSV output already exists: {output_tsv}")
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    tsv_bytes = tsv_payload(summary["samples"])
    summary["tsv_output"] = {
        "path": str(output_tsv.absolute()),
        "size_bytes": len(tsv_bytes),
        "sha256": hashlib.sha256(tsv_bytes).hexdigest(),
    }
    json_bytes = (
        json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")
    with output_tsv.open("xb") as handle:
        handle.write(tsv_bytes)
    with output_json.open("xb") as handle:
        handle.write(json_bytes)


def parse_samples(value: str) -> tuple[str, ...]:
    result = tuple(part.strip() for part in value.split(",") if part.strip())
    require(result, "--samples is empty")
    for sample in result:
        require(SAMPLE_RE.fullmatch(sample) is not None, f"invalid sample name: {sample}")
    require(len(result) == len(set(result)), "--samples contains duplicates")
    return result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-root", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument(
        "--samples",
        help=(
            "comma-separated completed sample subset; requires --allow-partial. "
            "Omit for the completed cohort_receipt.json scope"
        ),
    )
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="explicitly permit a sample subset without a completed cohort receipt",
    )
    return parser.parse_args(argv)


def execute(args: argparse.Namespace) -> dict[str, Any]:
    cohort_root = args.cohort_root.expanduser().resolve(strict=True)
    require(cohort_root.is_dir(), f"cohort root is not a directory: {cohort_root}")
    if args.samples is None:
        require(
            not args.allow_partial,
            "--allow-partial requires an explicit --samples subset",
        )
        cohort_receipt, samples, receipt_identity = validate_cohort_receipt(cohort_root)
        partial = False
    else:
        require(args.allow_partial, "explicit --samples requires --allow-partial")
        samples = list(parse_samples(args.samples))
        cohort_receipt = None
        receipt_identity = None
        partial = True
    summary = build_cohort_summary(
        cohort_root,
        samples,
        partial=partial,
        cohort_receipt=cohort_receipt,
        cohort_receipt_identity=receipt_identity,
    )
    write_outputs_exclusive(summary, args.output_json, args.output_tsv)
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        summary = execute(args)
    except (SummaryError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        "PASS "
        f"samples={len(summary['samples'])} "
        f"groups={summary['totals']['groups_total']} "
        f"ranked={summary['totals']['ranked_units']} "
        f"family_abstain={summary['totals']['mutation_family_abstain_units']} "
        f"output_json={args.output_json} output_tsv={args.output_tsv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
