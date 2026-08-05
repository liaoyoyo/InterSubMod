#!/usr/bin/env python3
"""Deterministically merge split pattern-methyl analysis runs.

The confirmatory and secondary analyses are executed with different
permutation budgets but belong to one frozen report universe.  This utility
validates both run bundles, optionally overlays higher-permutation
confirmatory rows on the same exact unit key, then recomputes multiplicity and
the reader-facing frozen assessment gates with the authoritative analyzer.

No input assessment, q-value, or adjusted p-value is trusted during the
merge.  A failure leaves no authoritative output directory and moves any
partial staging directory into ``_failed_staging_archive`` for audit.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence


MERGE_SCHEMA_NAME = "intersubmod.pattern_methyl_run_merge"
MERGE_SCHEMA_VERSION = "1.0.0"
KEY_FIELDS = ("dataset", "chrom", "region_id", "hp_raw")
CONFIRMATORY_FAMILY = "CONFIRMATORY_FULL4_OR_LONG"
SECONDARY_FAMILY = "SECONDARY_PAIR_CONTRAST"
CONFIRMATORY_PERMUTATIONS = 999
SECONDARY_PERMUTATIONS = 199
REFINED_PERMUTATIONS = 49999
OUTPUT_EVIDENCE = "pattern_methyl_evidence.combined.v1.tsv.gz"
OUTPUT_DETAILS = "pattern_methyl_details.combined.v1.jsonl.gz"
OUTPUT_SUMMARY = "analysis_summary.combined.json"
OUTPUT_RECEIPT = "merge_receipt.json"
OUTPUT_RECEIPT_SHA256 = "merge_receipt.json.sha256"
AUTHORITATIVE_RESULT_STATUS = "AUTHORITATIVE_WITHIN_DECLARED_FAMILY"
PROVISIONAL_RESULT_STATUS = (
    "PROVISIONAL_SUBSET_REFINEMENT_REQUIRES_FULL_FAMILY_MERGE"
)
REQUIRED_SOURCE_FIELDS = frozenset(
    {
        "pattern_counts",
        "artifact_catalog",
        "candidate_shards",
        "candidate_manifest",
        "unit_key_file",
        "unit_key_receipt",
        "screen_evidence",
        "topology_root",
        "topology_jsonl",
    }
)
REFINEMENT_SOURCE_FIELDS = frozenset(
    {"unit_key_file", "unit_key_receipt", "screen_evidence"}
)

# Only these evidence fields may differ in a higher-permutation refinement.
# Multiplicity and assessment fields are discarded and recomputed.
REFINEMENT_MUTABLE_FIELDS = frozenset(
    {
        "permanova_p",
        "permanova_permutations_requested",
        "permanova_permutations_realized",
        "permdisp_p",
        "multiplicity_family",
        "q_by",
        "p_holm",
        "assessment",
    }
)


class MergeContractError(RuntimeError):
    """Raised when a merge input violates a frozen identity or schema."""


def _load_analysis_module() -> Any:
    script = Path(__file__).resolve().with_name("analyze_pattern_methylation.py")
    if not script.is_file():
        raise MergeContractError(f"authoritative analyzer missing: {script}")
    module_name = "_intersubmod_pattern_methyl_analysis_for_merge"
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        return loaded
    spec = importlib.util.spec_from_file_location(module_name, script)
    if spec is None or spec.loader is None:
        raise MergeContractError(f"cannot import authoritative analyzer: {script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


ANALYSIS = _load_analysis_module()


@dataclass(frozen=True)
class RunBundle:
    label: str
    family: str
    evidence_path: Path
    details_path: Path
    summary_path: Path
    rows: tuple[dict[str, str], ...]
    row_by_key: Mapping[tuple[str, str, str, str], dict[str, str]]
    details: tuple[dict[str, Any], ...]
    detail_by_key: Mapping[tuple[str, str, str, str], dict[str, Any]]
    summary: Mapping[str, Any]
    identities: Mapping[str, Mapping[str, Any]]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _require_regular_nonempty(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise MergeContractError(f"{label} missing or empty: {resolved}")
    return resolved


def file_identity(path: Path, *, output_name_only: bool = False) -> dict[str, Any]:
    resolved = _require_regular_nonempty(path, "artifact")
    return {
        "path": resolved.name if output_name_only else str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


@contextlib.contextmanager
def open_text(path: Path) -> Iterator[io.TextIOBase]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8", newline="") as handle:
            yield handle


@contextlib.contextmanager
def deterministic_gzip_text_writer(
    path: Path, *, newline: str | None = None
) -> Iterator[io.TextIOWrapper]:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, compresslevel=6, mtime=0
        ) as gzip_handle:
            with io.TextIOWrapper(
                gzip_handle, encoding="utf-8", newline=newline, write_through=True
            ) as text_handle:
                yield text_handle


def exact_key(row: Mapping[str, Any]) -> tuple[str, str, str, str]:
    key = tuple(str(row.get(field, "")).strip() for field in KEY_FIELDS)
    if any(not value for value in key):
        raise MergeContractError(f"empty exact-key field: {dict(zip(KEY_FIELDS, key))}")
    return key  # type: ignore[return-value]


def row_sort_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        str(row["dataset"]),
        ANALYSIS.natural_chrom_key(str(row["chrom"])),
        str(row["region_id"]),
        str(row["hp_raw"]),
    )


def _parse_int(
    value: object, *, label: str, minimum: int | None = None
) -> int:
    try:
        number = int(str(value))
    except (TypeError, ValueError) as exc:
        raise MergeContractError(f"{label} is not an integer: {value!r}") from exc
    if minimum is not None and number < minimum:
        raise MergeContractError(f"{label} must be >= {minimum}: {number}")
    return number


def _parse_probability(value: object, *, label: str) -> float:
    try:
        number = float(str(value))
    except (TypeError, ValueError) as exc:
        raise MergeContractError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(number) or not 0.0 <= number <= 1.0:
        raise MergeContractError(f"{label} is outside [0,1]: {number}")
    return number


def _validate_evaluable_row(
    row: Mapping[str, str], key: tuple[str, str, str, str], label: str
) -> None:
    requested = _parse_int(
        row["permanova_permutations_requested"],
        label=f"{label} {key} requested permutations",
        minimum=1,
    )
    realized = _parse_int(
        row["permanova_permutations_realized"],
        label=f"{label} {key} realized permutations",
        minimum=1,
    )
    if realized > requested:
        raise MergeContractError(
            f"{label} {key} realized permutations exceed requested: "
            f"{realized}>{requested}"
        )
    _parse_probability(row["permanova_p"], label=f"{label} {key} permanova_p")
    _parse_probability(row["permdisp_p"], label=f"{label} {key} permdisp_p")
    numeric_gates: dict[str, float] = {}
    for field in (
        "permanova_r2",
        "best_pair_distance_contrast",
        "best_pair_standardized_effect",
        "max_geometry_smd",
    ):
        try:
            number = float(row[field])
        except (TypeError, ValueError) as exc:
            raise MergeContractError(
                f"{label} {key} missing numeric gate field {field}"
            ) from exc
        if not math.isfinite(number):
            raise MergeContractError(
                f"{label} {key} non-finite gate field {field}: {number}"
            )
        numeric_gates[field] = number
    for field in ("equal_n_retention", "rarefaction_retention"):
        raw = str(row[field]).strip()
        if not raw and numeric_gates["permanova_r2"] == 0.0:
            continue
        try:
            number = float(raw)
        except (TypeError, ValueError) as exc:
            raise MergeContractError(
                f"{label} {key} missing numeric gate field {field}"
            ) from exc
        if not math.isfinite(number):
            raise MergeContractError(
                f"{label} {key} non-finite gate field {field}: {number}"
            )


def load_evidence(
    path: Path, *, label: str, expected_family: str
) -> tuple[
    tuple[dict[str, str], ...],
    dict[tuple[str, str, str, str], dict[str, str]],
]:
    resolved = _require_regular_nonempty(path, f"{label} evidence")
    rows: list[dict[str, str]] = []
    by_key: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with open_text(resolved) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_header = list(ANALYSIS.RESULT_FIELDS)
        if reader.fieldnames != expected_header:
            raise MergeContractError(
                f"{label} evidence header mismatch: "
                f"observed={reader.fieldnames!r} expected={expected_header!r}"
            )
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise MergeContractError(
                    f"{label} evidence line {line_number} has extra cells"
                )
            row = {field: str(raw.get(field, "")) for field in expected_header}
            if row["schema_version"] != ANALYSIS.SCHEMA_VERSION:
                raise MergeContractError(
                    f"{label} evidence schema version mismatch at line {line_number}"
                )
            key = exact_key(row)
            if key in by_key:
                raise MergeContractError(f"duplicate exact key in {label} evidence: {key}")
            computed_family = ANALYSIS.assign_multiplicity_family(row)
            if computed_family != expected_family:
                raise MergeContractError(
                    f"{label} evidence key {key} belongs to {computed_family}, "
                    f"not {expected_family}"
                )
            if row["multiplicity_family"] != computed_family:
                raise MergeContractError(
                    f"{label} evidence key {key} declares inconsistent "
                    f"multiplicity_family={row['multiplicity_family']!r}"
                )
            status = row["evaluation_status"]
            if status not in {"EVALUABLE", "NOT_EVALUABLE"}:
                raise MergeContractError(
                    f"{label} evidence key {key} has invalid evaluation_status={status!r}"
                )
            if status == "EVALUABLE":
                _validate_evaluable_row(row, key, label)
            elif not row["invalid_reason"]:
                raise MergeContractError(
                    f"{label} non-evaluable key {key} lacks invalid_reason"
                )
            rows.append(row)
            by_key[key] = row
    if not rows:
        raise MergeContractError(f"{label} evidence has no rows")
    return tuple(rows), by_key


def load_details(
    path: Path,
    *,
    label: str,
    evidence_by_key: Mapping[tuple[str, str, str, str], Mapping[str, str]],
) -> tuple[
    tuple[dict[str, Any], ...],
    dict[tuple[str, str, str, str], dict[str, Any]],
]:
    resolved = _require_regular_nonempty(path, f"{label} details")
    details: list[dict[str, Any]] = []
    by_key: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    with open_text(resolved) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                raise MergeContractError(
                    f"{label} details contains blank line {line_number}"
                )
            try:
                detail = json.loads(line)
            except json.JSONDecodeError as exc:
                raise MergeContractError(
                    f"{label} details invalid JSON at line {line_number}"
                ) from exc
            if not isinstance(detail, dict):
                raise MergeContractError(
                    f"{label} details line {line_number} is not an object"
                )
            if detail.get("schema_name") != f"{ANALYSIS.SCHEMA_NAME}.detail":
                raise MergeContractError(
                    f"{label} details schema mismatch at line {line_number}"
                )
            if detail.get("schema_version") != ANALYSIS.SCHEMA_VERSION:
                raise MergeContractError(
                    f"{label} details version mismatch at line {line_number}"
                )
            key = exact_key(detail)
            if key in by_key:
                raise MergeContractError(f"duplicate exact key in {label} details: {key}")
            evidence = evidence_by_key.get(key)
            if evidence is None:
                raise MergeContractError(
                    f"{label} detail key absent from evidence: {key}"
                )
            if detail.get("assessment") != evidence["assessment"]:
                raise MergeContractError(
                    f"{label} detail/evidence assessment mismatch: {key}"
                )
            details.append(detail)
            by_key[key] = detail
    return tuple(details), by_key


def _counter_dict(counter: Counter[str]) -> dict[str, int]:
    return dict(sorted((str(key), int(value)) for key, value in counter.items()))


def expected_counts(
    rows: Sequence[Mapping[str, str]], details_n: int
) -> dict[str, Any]:
    return {
        "units_total": len(rows),
        "units_evaluable": sum(
            row["evaluation_status"] == "EVALUABLE" for row in rows
        ),
        "detail_records": details_n,
        "assessment": _counter_dict(Counter(row["assessment"] for row in rows)),
        "multiplicity_family": _counter_dict(
            Counter(row["multiplicity_family"] for row in rows)
        ),
        "invalid_reason": _counter_dict(
            Counter(
                row["invalid_reason"]
                for row in rows
                if row["evaluation_status"] != "EVALUABLE"
            )
        ),
    }


def _binding_path(summary_path: Path, binding: Mapping[str, Any]) -> Path:
    raw = binding.get("path")
    if not isinstance(raw, str) or not raw:
        raise MergeContractError(f"summary output binding lacks path: {binding!r}")
    path = Path(raw)
    if not path.is_absolute():
        path = summary_path.parent / path
    return path.resolve()


def _validate_summary_output_binding(
    summary_path: Path,
    binding: object,
    actual_path: Path,
    *,
    label: str,
) -> None:
    if not isinstance(binding, dict):
        raise MergeContractError(f"{label} summary output binding is not an object")
    if _binding_path(summary_path, binding) != actual_path.resolve():
        raise MergeContractError(f"{label} summary output path does not bind input")
    identity = file_identity(actual_path)
    if binding.get("sha256") != identity["sha256"]:
        raise MergeContractError(f"{label} summary output sha256 mismatch")
    if int(binding.get("size_bytes", -1)) != identity["size_bytes"]:
        raise MergeContractError(f"{label} summary output size mismatch")


def load_and_validate_summary(
    path: Path,
    *,
    label: str,
    rows: Sequence[Mapping[str, str]],
    details_n: int,
    evidence_path: Path,
    details_path: Path,
) -> dict[str, Any]:
    resolved = _require_regular_nonempty(path, f"{label} summary")
    try:
        summary = json.loads(resolved.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise MergeContractError(f"{label} summary is invalid JSON") from exc
    if not isinstance(summary, dict):
        raise MergeContractError(f"{label} summary root is not an object")
    if summary.get("schema_name") != f"{ANALYSIS.SCHEMA_NAME}.summary":
        raise MergeContractError(f"{label} summary schema_name mismatch")
    if summary.get("schema_version") != ANALYSIS.SCHEMA_VERSION:
        raise MergeContractError(f"{label} summary schema_version mismatch")
    config = summary.get("config")
    if not isinstance(config, dict):
        raise MergeContractError(f"{label} summary config missing")
    configured_permutations = _parse_int(
        config.get("permutations"),
        label=f"{label} summary configured permutations",
        minimum=1,
    )
    requested_permutations = {
        _parse_int(
            row["permanova_permutations_requested"],
            label=f"{label} evidence requested permutations",
            minimum=1,
        )
        for row in rows
        if row["evaluation_status"] == "EVALUABLE"
    }
    if requested_permutations and requested_permutations != {
        configured_permutations
    }:
        raise MergeContractError(
            f"{label} summary/evidence permutation budget mismatch: "
            f"summary={configured_permutations} evidence={sorted(requested_permutations)}"
        )
    expected = expected_counts(rows, details_n)
    observed = summary.get("counts")
    if observed != expected:
        raise MergeContractError(
            f"{label} summary counts mismatch: observed={observed!r} expected={expected!r}"
        )
    outputs = summary.get("outputs")
    if not isinstance(outputs, dict):
        raise MergeContractError(f"{label} summary outputs missing")
    _validate_summary_output_binding(
        resolved,
        outputs.get("evidence"),
        evidence_path,
        label=f"{label} evidence",
    )
    _validate_summary_output_binding(
        resolved,
        outputs.get("details"),
        details_path,
        label=f"{label} details",
    )
    return summary


def load_bundle(
    *,
    label: str,
    family: str,
    evidence_path: Path,
    details_path: Path,
    summary_path: Path,
) -> RunBundle:
    evidence_path = _require_regular_nonempty(evidence_path, f"{label} evidence")
    details_path = _require_regular_nonempty(details_path, f"{label} details")
    summary_path = _require_regular_nonempty(summary_path, f"{label} summary")
    rows, row_by_key = load_evidence(
        evidence_path, label=label, expected_family=family
    )
    details, detail_by_key = load_details(
        details_path, label=label, evidence_by_key=row_by_key
    )
    summary = load_and_validate_summary(
        summary_path,
        label=label,
        rows=rows,
        details_n=len(details),
        evidence_path=evidence_path,
        details_path=details_path,
    )
    return RunBundle(
        label=label,
        family=family,
        evidence_path=evidence_path,
        details_path=details_path,
        summary_path=summary_path,
        rows=rows,
        row_by_key=row_by_key,
        details=details,
        detail_by_key=detail_by_key,
        summary=summary,
        identities={
            "evidence": file_identity(evidence_path),
            "details": file_identity(details_path),
            "summary": file_identity(summary_path),
        },
    )


def validate_bundle_budget(bundle: RunBundle, expected_permutations: int) -> None:
    config = bundle.summary.get("config")
    if not isinstance(config, dict):
        raise MergeContractError(f"{bundle.label} summary config missing")
    configured = _parse_int(
        config.get("permutations"),
        label=f"{bundle.label} configured permutations",
        minimum=1,
    )
    if configured != expected_permutations:
        raise MergeContractError(
            f"{bundle.label} must use frozen permutation budget "
            f"{expected_permutations}, observed {configured}"
        )
    for key, row in bundle.row_by_key.items():
        if row["evaluation_status"] != "EVALUABLE":
            continue
        requested = _parse_int(
            row["permanova_permutations_requested"],
            label=f"{bundle.label} requested permutations {key}",
            minimum=1,
        )
        realized = _parse_int(
            row["permanova_permutations_realized"],
            label=f"{bundle.label} realized permutations {key}",
            minimum=1,
        )
        if requested != expected_permutations or realized != expected_permutations:
            raise MergeContractError(
                f"{bundle.label} evidence does not realize frozen permutation "
                f"budget {expected_permutations} at {key}: "
                f"{requested}/{realized}"
            )


def _detail_identity_payload(detail: Mapping[str, Any]) -> dict[str, Any]:
    """Return the detail payload unaffected by multiplicity classification."""
    return {
        field: value
        for field, value in detail.items()
        if field != "assessment"
    }


def _validate_source_file_binding(
    bundle: RunBundle,
    binding: object,
    *,
    label: str,
    expected_path: Path | None = None,
    binding_document_path: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    if not isinstance(binding, dict):
        raise MergeContractError(f"{bundle.label} source {label} is not bound")
    path = _binding_path(binding_document_path or bundle.summary_path, binding)
    if expected_path is not None and path != expected_path.resolve():
        raise MergeContractError(
            f"{bundle.label} source {label} path does not bind expected input"
        )
    identity = file_identity(path)
    if binding.get("sha256") != identity["sha256"]:
        raise MergeContractError(
            f"{bundle.label} source {label} sha256 mismatch"
        )
    try:
        bound_size = int(binding.get("size_bytes", -1))
    except (TypeError, ValueError) as exc:
        raise MergeContractError(
            f"{bundle.label} source {label} size is invalid"
        ) from exc
    if bound_size != identity["size_bytes"]:
        raise MergeContractError(f"{bundle.label} source {label} size mismatch")
    return path, identity


def _adaptive_eligible_keys(
    confirmatory: RunBundle,
) -> set[tuple[str, str, str, str]]:
    floor = 1.0 / (CONFIRMATORY_PERMUTATIONS + 1)
    selected: set[tuple[str, str, str, str]] = set()
    for key, row in confirmatory.row_by_key.items():
        if row["evaluation_status"] != "EVALUABLE":
            continue
        try:
            eligible = (
                math.isclose(
                    float(row["permanova_p"]),
                    floor,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                )
                and float(row["permanova_r2"]) >= 0.10
                and float(row["best_pair_distance_contrast"]) >= 0.10
                and float(row["best_pair_standardized_effect"]) >= 0.50
                and float(row["permdisp_p"]) >= 0.05
                and float(row["max_geometry_smd"]) < 0.50
                and ANALYSIS.bool_value(row["all_states_n8"])
                and float(row["equal_n_retention"]) >= 0.50
                and float(row["rarefaction_retention"]) >= 0.50
            )
        except (TypeError, ValueError) as exc:
            raise MergeContractError(
                f"confirmatory adaptive gate is invalid at {key}"
            ) from exc
        if eligible:
            selected.add(key)
    return selected


def validate_refinement_receipt(
    refined_bundle: RunBundle,
    confirmatory: RunBundle,
    refined_keys: set[tuple[str, str, str, str]],
) -> None:
    sources = refined_bundle.summary.get("sources")
    if not isinstance(sources, dict):
        raise MergeContractError("refined confirmatory summary sources missing")
    unit_key_path, unit_key_identity = _validate_source_file_binding(
        refined_bundle,
        sources.get("unit_key_file"),
        label="unit_key_file",
    )
    receipt_path, _receipt_identity = _validate_source_file_binding(
        refined_bundle,
        sources.get("unit_key_receipt"),
        label="unit_key_receipt",
    )
    screen_path, screen_identity = _validate_source_file_binding(
        refined_bundle,
        sources.get("screen_evidence"),
        label="screen_evidence",
    )
    confirmatory_identity = file_identity(confirmatory.evidence_path)
    if screen_identity != confirmatory_identity:
        raise MergeContractError(
            "refined receipt screen_evidence path/size/SHA does not equal "
            "confirmatory evidence"
        )
    if screen_path != confirmatory.evidence_path.resolve():
        raise MergeContractError(
            "refined receipt screen_evidence path does not equal "
            "confirmatory evidence"
        )

    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise MergeContractError("refined unit-key receipt is invalid JSON") from exc
    if not isinstance(receipt, dict):
        raise MergeContractError("refined unit-key receipt root is not an object")
    if receipt.get("schema_name") != ANALYSIS.ADAPTIVE_SELECTION_SCHEMA_NAME:
        raise MergeContractError("refined unit-key receipt schema_name mismatch")
    if receipt.get("schema_version") != ANALYSIS.ADAPTIVE_SELECTION_SCHEMA_VERSION:
        raise MergeContractError("refined unit-key receipt schema_version mismatch")
    if receipt.get("all_pass") is not True:
        raise MergeContractError("refined unit-key receipt is not all_pass")
    contract = receipt.get("contract")
    if not isinstance(contract, dict):
        raise MergeContractError("refined unit-key receipt contract missing")
    if contract.get("family") != CONFIRMATORY_FAMILY:
        raise MergeContractError("refined unit-key receipt family mismatch")
    if (
        _parse_int(
            contract.get("screen_permutations"),
            label="refined receipt screen permutations",
            minimum=1,
        )
        != CONFIRMATORY_PERMUTATIONS
        or _parse_int(
            contract.get("refinement_permutations"),
            label="refined receipt refinement permutations",
            minimum=1,
        )
        != REFINED_PERMUTATIONS
    ):
        raise MergeContractError("refined unit-key receipt budget mismatch")
    try:
        screen_floor = float(contract.get("screen_floor", float("nan")))
    except (TypeError, ValueError) as exc:
        raise MergeContractError(
            "refined unit-key receipt screen floor is invalid"
        ) from exc
    if not math.isclose(
        screen_floor,
        1.0 / (CONFIRMATORY_PERMUTATIONS + 1),
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise MergeContractError("refined unit-key receipt screen floor mismatch")
    counts = receipt.get("counts")
    if not isinstance(counts, dict):
        raise MergeContractError("refined unit-key receipt counts missing")
    selected_count = _parse_int(
        counts.get("selected_for_refinement"),
        label="refined receipt selected count",
        minimum=0,
    )
    if selected_count != len(refined_keys):
        raise MergeContractError(
            "refined receipt selected count does not equal refined keys: "
            f"{selected_count}!={len(refined_keys)}"
        )

    outputs = receipt.get("outputs")
    inputs = receipt.get("inputs")
    if not isinstance(outputs, dict) or not isinstance(inputs, dict):
        raise MergeContractError("refined unit-key receipt bindings missing")
    receipt_unit_path, receipt_unit_identity = _validate_source_file_binding(
        refined_bundle,
        outputs.get("unit_keys"),
        label="receipt.outputs.unit_keys",
        expected_path=unit_key_path,
        binding_document_path=receipt_path,
    )
    receipt_screen_path, receipt_screen_identity = _validate_source_file_binding(
        refined_bundle,
        inputs.get("screen_evidence"),
        label="receipt.inputs.screen_evidence",
        expected_path=confirmatory.evidence_path,
        binding_document_path=receipt_path,
    )
    if (
        receipt_unit_path != unit_key_path
        or receipt_unit_identity != unit_key_identity
        or receipt_screen_path != screen_path
        or receipt_screen_identity != screen_identity
    ):
        raise MergeContractError(
            "refined analyzer sources drift from unit-key receipt bindings"
        )
    try:
        unit_keys = ANALYSIS.load_unit_keys(unit_key_path)
    except Exception as exc:
        raise MergeContractError(
            f"refined unit-key TSV is invalid: {type(exc).__name__}: {exc}"
        ) from exc
    if unit_keys != refined_keys:
        raise MergeContractError(
            "refined unit-key TSV does not exactly equal refined evidence keys"
        )


def validate_refinement(
    refined_bundle: RunBundle,
    confirmatory: RunBundle,
) -> None:
    if (
        refined_bundle.summary.get("result_status")
        != PROVISIONAL_RESULT_STATUS
    ):
        raise MergeContractError(
            "refined confirmatory run is not marked provisional subset refinement"
        )
    refined_keys = set(refined_bundle.row_by_key)
    expected_refined_keys = _adaptive_eligible_keys(confirmatory)
    if refined_keys != expected_refined_keys:
        raise MergeContractError(
            "refined exact keys do not equal frozen adaptive-gate selection: "
            f"missing={sorted(expected_refined_keys - refined_keys)[:3]} "
            f"extra={sorted(refined_keys - expected_refined_keys)[:3]}"
        )
    refined_detail_keys = set(refined_bundle.detail_by_key)
    if refined_detail_keys != refined_keys:
        missing = sorted(refined_keys - refined_detail_keys)
        extra = sorted(refined_detail_keys - refined_keys)
        raise MergeContractError(
            "refined confirmatory details must exactly cover refined evidence: "
            f"missing={missing[:3]} extra={extra[:3]}"
        )

    base_config = dict(confirmatory.summary.get("config", {}))
    refined_config = dict(refined_bundle.summary.get("config", {}))
    base_config.pop("permutations", None)
    refined_config.pop("permutations", None)
    if refined_config != base_config:
        raise MergeContractError(
            "refined confirmatory config drift outside permutations: "
            f"base={base_config!r} refined={refined_config!r}"
        )
    validate_source_compatibility(
        confirmatory,
        refined_bundle,
        allow_refinement_source_difference=True,
    )
    refined_sources = refined_bundle.summary.get("sources", {})
    if not isinstance(refined_sources.get("unit_key_file"), dict):
        raise MergeContractError(
            "refined confirmatory source contract lacks unit_key_file binding"
        )
    validate_refinement_receipt(refined_bundle, confirmatory, refined_keys)

    for key, refined in refined_bundle.row_by_key.items():
        base = confirmatory.row_by_key.get(key)
        if base is None:
            raise MergeContractError(
                f"refined confirmatory exact key absent from base run: {key}"
            )
        if base["evaluation_status"] != "EVALUABLE" or refined[
            "evaluation_status"
        ] != "EVALUABLE":
            raise MergeContractError(
                f"refinement requires evaluable base and refined rows: {key}"
            )
        for field in ANALYSIS.RESULT_FIELDS:
            if field in REFINEMENT_MUTABLE_FIELDS:
                continue
            if refined[field] != base[field]:
                raise MergeContractError(
                    f"refinement identity/statistic drift at {key} field={field}: "
                    f"base={base[field]!r} refined={refined[field]!r}"
                )
        base_requested = _parse_int(
            base["permanova_permutations_requested"],
            label=f"base refinement requested {key}",
            minimum=1,
        )
        refined_requested = _parse_int(
            refined["permanova_permutations_requested"],
            label=f"refined requested {key}",
            minimum=1,
        )
        base_realized = _parse_int(
            base["permanova_permutations_realized"],
            label=f"base refinement realized {key}",
            minimum=1,
        )
        refined_realized = _parse_int(
            refined["permanova_permutations_realized"],
            label=f"refined realized {key}",
            minimum=1,
        )
        if refined_requested <= base_requested or refined_realized <= base_realized:
            raise MergeContractError(
                f"refinement permutations did not increase at {key}: "
                f"requested {base_requested}->{refined_requested}, "
                f"realized {base_realized}->{refined_realized}"
            )
        if key not in confirmatory.detail_by_key:
            raise MergeContractError(
                f"refined confirmatory key lacks base detail record: {key}"
            )
        refined_detail = refined_bundle.detail_by_key[key]
        if _detail_identity_payload(refined_detail) != _detail_identity_payload(
            confirmatory.detail_by_key[key]
        ):
            raise MergeContractError(
                f"refinement detail identity/statistic drift at {key}"
            )


def source_contract(
    bundle: RunBundle, *, ignore_refinement_sources: bool
) -> dict[str, Any]:
    sources = bundle.summary.get("sources")
    if not isinstance(sources, dict):
        raise MergeContractError(f"{bundle.label} summary sources missing")
    missing = REQUIRED_SOURCE_FIELDS - set(sources)
    if missing:
        raise MergeContractError(
            f"{bundle.label} summary sources missing fields: {sorted(missing)}"
        )
    if not isinstance(sources.get("candidate_manifest"), dict):
        raise MergeContractError(
            f"{bundle.label} source contract requires candidate_manifest binding"
        )
    normalized = dict(sources)
    if ignore_refinement_sources:
        for field in REFINEMENT_SOURCE_FIELDS:
            normalized.pop(field, None)
    return normalized


def validate_source_compatibility(
    first: RunBundle,
    second: RunBundle,
    *,
    allow_refinement_source_difference: bool = False,
) -> None:
    first_sources = source_contract(
        first, ignore_refinement_sources=allow_refinement_source_difference
    )
    second_sources = source_contract(
        second, ignore_refinement_sources=allow_refinement_source_difference
    )
    if first_sources != second_sources:
        raise MergeContractError(
            f"source contract drift between {first.label} and {second.label}"
        )


def validate_base_run_compatibility(
    confirmatory: RunBundle, secondary: RunBundle
) -> None:
    for bundle in (confirmatory, secondary):
        if bundle.summary.get("result_status") != AUTHORITATIVE_RESULT_STATUS:
            raise MergeContractError(
                f"{bundle.label} run is not authoritative within its family"
            )
        sources = bundle.summary.get("sources")
        if not isinstance(sources, dict):
            raise MergeContractError(f"{bundle.label} summary sources missing")
        non_null = [
            field for field in REFINEMENT_SOURCE_FIELDS
            if sources.get(field) is not None
        ]
        if non_null:
            raise MergeContractError(
                f"{bundle.label} base sources must set refinement bindings to "
                f"None: {sorted(non_null)}"
            )
    confirmatory_config = dict(confirmatory.summary.get("config", {}))
    secondary_config = dict(secondary.summary.get("config", {}))
    confirmatory_config.pop("permutations", None)
    secondary_config.pop("permutations", None)
    if confirmatory_config != secondary_config:
        raise MergeContractError(
            "confirmatory/secondary config drift outside permutations: "
            f"confirmatory={confirmatory_config!r} "
            f"secondary={secondary_config!r}"
        )
    validate_source_compatibility(confirmatory, secondary)


def merge_rows_and_details(
    confirmatory: RunBundle,
    secondary: RunBundle,
    refined: RunBundle | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    overlap = set(confirmatory.row_by_key).intersection(secondary.row_by_key)
    if overlap:
        raise MergeContractError(
            f"confirmatory/secondary exact keys overlap: {sorted(overlap)[:3]}"
        )
    combined: dict[tuple[str, str, str, str], dict[str, Any]] = {
        key: dict(row)
        for key, row in {
            **confirmatory.row_by_key,
            **secondary.row_by_key,
        }.items()
    }
    if refined is not None:
        for key, row in refined.row_by_key.items():
            combined[key] = dict(row)
    rows = sorted(combined.values(), key=row_sort_key)
    for row in rows:
        row["multiplicity_family"] = ""
        row["q_by"] = ""
        row["p_holm"] = ""
        row["assessment"] = "PENDING_MULTIPLICITY"
    ANALYSIS.adjust_p_values(rows)
    ANALYSIS.classify_rows(rows)

    detail_overlap = set(confirmatory.detail_by_key).intersection(
        secondary.detail_by_key
    )
    if detail_overlap:
        raise MergeContractError(
            f"confirmatory/secondary detail keys overlap: {sorted(detail_overlap)[:3]}"
        )
    evidence_assessment = {exact_key(row): row["assessment"] for row in rows}
    combined_details: dict[
        tuple[str, str, str, str], Mapping[str, Any]
    ] = {
        exact_key(detail): detail
        for detail in list(confirmatory.details) + list(secondary.details)
    }
    if refined is not None:
        for key, detail in refined.detail_by_key.items():
            combined_details[key] = detail

    details: list[dict[str, Any]] = []
    for detail in combined_details.values():
        copied = dict(detail)
        key = exact_key(copied)
        copied["assessment"] = evidence_assessment[key]
        details.append(copied)
    details.sort(key=row_sort_key)
    return rows, details


def write_evidence(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    with deterministic_gzip_text_writer(path, newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ANALYSIS.RESULT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: ANALYSIS.output_value(row.get(field, ""))
                    for field in ANALYSIS.RESULT_FIELDS
                }
            )


def write_details(path: Path, details: Sequence[Mapping[str, Any]]) -> None:
    with deterministic_gzip_text_writer(path) as handle:
        for detail in details:
            handle.write(
                json.dumps(
                    detail,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                )
            )
            handle.write("\n")


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )


def _output_identity(path: Path) -> dict[str, Any]:
    return file_identity(path, output_name_only=True)


def _combined_counts(
    rows: Sequence[Mapping[str, Any]],
    details_n: int,
    *,
    confirmatory_n: int,
    secondary_n: int,
    refined_n: int,
) -> dict[str, Any]:
    base = expected_counts(rows, details_n)  # type: ignore[arg-type]
    base.update(
        {
            "confirmatory_input_units": confirmatory_n,
            "secondary_input_units": secondary_n,
            "refined_confirmatory_units": refined_n,
        }
    )
    return base


def merge_to_staging(
    staging: Path,
    *,
    confirmatory: RunBundle,
    secondary: RunBundle,
    refined: RunBundle | None,
) -> dict[str, Any]:
    rows, details = merge_rows_and_details(confirmatory, secondary, refined)
    evidence_path = staging / OUTPUT_EVIDENCE
    details_path = staging / OUTPUT_DETAILS
    summary_path = staging / OUTPUT_SUMMARY
    receipt_path = staging / OUTPUT_RECEIPT
    receipt_sidecar_path = staging / OUTPUT_RECEIPT_SHA256

    write_evidence(evidence_path, rows)
    write_details(details_path, details)
    counts = _combined_counts(
        rows,
        len(details),
        confirmatory_n=len(confirmatory.rows),
        secondary_n=len(secondary.rows),
        refined_n=len(refined.rows) if refined is not None else 0,
    )
    input_identities: dict[str, Any] = {
        "confirmatory": confirmatory.identities,
        "secondary": secondary.identities,
        "refined_confirmatory": (
            refined.identities if refined is not None else None
        ),
    }
    summary = {
        "schema_name": f"{ANALYSIS.SCHEMA_NAME}.combined_summary",
        "schema_version": MERGE_SCHEMA_VERSION,
        "claim_ceiling": (
            "pattern-conditioned regional methylation association only; "
            "not ancestry, clone, causality, or topology rescoring"
        ),
        "merge_contract": {
            "schema_name": MERGE_SCHEMA_NAME,
            "schema_version": MERGE_SCHEMA_VERSION,
            "exact_key_fields": list(KEY_FIELDS),
            "refinement_mutable_fields": sorted(REFINEMENT_MUTABLE_FIELDS),
            "multiplicity": (
                "BY and Holm recomputed independently within the frozen "
                "confirmatory and secondary families"
            ),
            "assessment": (
                "recomputed with authoritative "
                "analyze_pattern_methylation.classify_rows"
            ),
            "details": (
                "refined detail records overlay matching base keys; "
                "assessment refreshed"
            ),
        },
        "config_by_family": {
            CONFIRMATORY_FAMILY: confirmatory.summary.get("config"),
            SECONDARY_FAMILY: secondary.summary.get("config"),
            "refined_confirmatory": (
                refined.summary.get("config") if refined is not None else None
            ),
        },
        "counts": counts,
        "inputs": input_identities,
        "outputs": {
            "evidence": _output_identity(evidence_path),
            "details": _output_identity(details_path),
        },
    }
    write_json(summary_path, summary)
    outputs = {
        "evidence": _output_identity(evidence_path),
        "details": _output_identity(details_path),
        "summary": _output_identity(summary_path),
    }
    receipt = {
        "schema_name": f"{MERGE_SCHEMA_NAME}.receipt",
        "schema_version": MERGE_SCHEMA_VERSION,
        "pass": True,
        "contracts": summary["merge_contract"],
        "counts": counts,
        "inputs": input_identities,
        "outputs": outputs,
    }
    write_json(receipt_path, receipt)
    receipt_sha256 = sha256_file(receipt_path)
    receipt_sidecar_path.write_text(
        f"{receipt_sha256}  {OUTPUT_RECEIPT}\n", encoding="utf-8"
    )
    return {
        "counts": counts,
        "outputs": outputs,
        "receipt": {
            **_output_identity(receipt_path),
            "sidecar": _output_identity(receipt_sidecar_path),
        },
    }


def archive_failed_staging(staging: Path, output_parent: Path) -> Path:
    """Move a failed staging directory into a retained audit location."""
    archive_root = output_parent / "_failed_staging_archive"
    archive_root.mkdir(parents=True, exist_ok=True)
    destination = archive_root / staging.name
    suffix = 1
    while os.path.lexists(destination):
        destination = archive_root / f"{staging.name}.{suffix}"
        suffix += 1
    os.rename(staging, destination)
    return destination


def execute(args: argparse.Namespace) -> dict[str, Any]:
    output_dir = args.output_dir.expanduser().absolute()
    if os.path.lexists(output_dir):
        raise MergeContractError(f"refusing to overwrite output: {output_dir}")
    refined_paths = (
        args.refined_confirmatory_evidence,
        args.refined_confirmatory_details,
        args.refined_confirmatory_summary,
    )
    refined_presence = tuple(path is not None for path in refined_paths)
    if any(refined_presence) and not all(refined_presence):
        raise MergeContractError(
            "refined confirmatory evidence/details/summary must be provided together"
        )

    confirmatory = load_bundle(
        label="confirmatory",
        family=CONFIRMATORY_FAMILY,
        evidence_path=args.confirmatory_evidence,
        details_path=args.confirmatory_details,
        summary_path=args.confirmatory_summary,
    )
    secondary = load_bundle(
        label="secondary",
        family=SECONDARY_FAMILY,
        evidence_path=args.secondary_evidence,
        details_path=args.secondary_details,
        summary_path=args.secondary_summary,
    )
    validate_bundle_budget(confirmatory, CONFIRMATORY_PERMUTATIONS)
    validate_bundle_budget(secondary, SECONDARY_PERMUTATIONS)
    validate_base_run_compatibility(confirmatory, secondary)
    expected_refined_keys = _adaptive_eligible_keys(confirmatory)
    refined: RunBundle | None = None
    if all(refined_presence):
        refined = load_bundle(
            label="refined_confirmatory",
            family=CONFIRMATORY_FAMILY,
            evidence_path=args.refined_confirmatory_evidence,
            details_path=args.refined_confirmatory_details,
            summary_path=args.refined_confirmatory_summary,
        )
        validate_bundle_budget(refined, REFINED_PERMUTATIONS)
        validate_refinement(refined, confirmatory)
    elif expected_refined_keys:
        raise MergeContractError(
            "confirmatory screen contains frozen adaptive-gate candidates but "
            "no refined confirmatory bundle was provided: "
            f"{sorted(expected_refined_keys)[:3]}"
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.staging.",
            dir=output_dir.parent,
        )
    )
    try:
        result = merge_to_staging(
            staging,
            confirmatory=confirmatory,
            secondary=secondary,
            refined=refined,
        )
        if os.path.lexists(output_dir):
            raise MergeContractError(
                f"output appeared during staging; refusing publish: {output_dir}"
            )
        os.rename(staging, output_dir)
        return result
    except Exception as exc:
        if os.path.lexists(staging):
            try:
                archived = archive_failed_staging(staging, output_dir.parent)
            except Exception as archive_exc:
                raise MergeContractError(
                    "merge failed and staging archive failed; "
                    f"staging retained at {staging}; "
                    f"archive_error={type(archive_exc).__name__}: {archive_exc}"
                ) from exc
            raise MergeContractError(
                f"merge failed; staging archived at {archived}: "
                f"{type(exc).__name__}: {exc}"
            ) from exc
        raise


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--confirmatory-evidence", type=Path, required=True)
    parser.add_argument("--confirmatory-details", type=Path, required=True)
    parser.add_argument("--confirmatory-summary", type=Path, required=True)
    parser.add_argument("--secondary-evidence", type=Path, required=True)
    parser.add_argument("--secondary-details", type=Path, required=True)
    parser.add_argument("--secondary-summary", type=Path, required=True)
    parser.add_argument("--refined-confirmatory-evidence", type=Path)
    parser.add_argument("--refined-confirmatory-details", type=Path)
    parser.add_argument("--refined-confirmatory-summary", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    result = execute(parse_args(argv))
    print(
        json.dumps(
            {
                "pass": True,
                "counts": result["counts"],
                "receipt": result["receipt"],
            },
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MergeContractError as exc:
        print(f"FAIL_CLOSED: {exc}", file=sys.stderr)
        raise SystemExit(2)
