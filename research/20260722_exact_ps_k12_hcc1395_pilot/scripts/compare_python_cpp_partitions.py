#!/usr/bin/env python3
"""Fail-closed parity comparison for Python and C++ exact-PS partitions."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import sys
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCHEMA_NAME = "intersubmod.exact_ps_k12_python_cpp_comparison"
SCHEMA_VERSION = "1.0.0"
PYTHON_SCHEMA_NAME = "intersubmod.exact_ps_k12_partition"
PYTHON_SCHEMA_VERSION = "0.1.0"
CPP_SCHEMA_NAME = "intersubmod.exact_ps_k12_cpp_partition_summary"
CPP_SCHEMA_VERSION = "1.0.0"

RETAINED = "retained"
CUT = "cut"
UNAVOIDABLE = "unavoidable_span_gt_max_block_size"
DISPOSITIONS = (RETAINED, CUT, UNAVOIDABLE)
MAX_MISMATCH_SAMPLES = 20

UNITS_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "k", "positions1",
)
CONSTRAINT_FIELDS = (
    "dataset", "chrom", "unit_id", "constraint_id", "hp_family",
    "phase_set", "positions1", "call_codes", "molecule_weight",
    "pattern_count",
)
PYTHON_BLOCK_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "block_id", "block_index",
    "start1", "end1", "k", "positions1", "retained_molecule_weight",
    "retained_pattern_count",
)
PYTHON_MEMBERSHIP_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "site_index", "pos1",
    "unit_local_index", "block_id", "block_index",
)
PYTHON_DISPOSITION_FIELDS = (
    "dataset", "chrom", "unit_id", "constraint_id", "hp_family",
    "phase_set", "positions1", "call_codes", "molecule_weight",
    "pattern_count", "disposition", "span_sites", "crossed_cut_count",
    "retained_block_index",
)
CPP_BLOCK_FIELDS = PYTHON_BLOCK_FIELDS + (
    "start_index_zero_based", "end_index_exclusive_zero_based",
    "retained_constraint_count", "unit_cut_indices_zero_based",
    "unit_cut_gaps", "unit_cut_gap_sum", "unit_retained_molecule_weight",
    "unit_retained_pattern_count", "unit_total_molecule_weight",
    "unit_total_pattern_count",
)
CPP_MEMBERSHIP_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "unit_local_index", "pos1",
    "block_id", "block_index", "site_index_in_unit_zero_based",
    "site_index_in_block_zero_based",
)
CPP_DISPOSITION_FIELDS = PYTHON_DISPOSITION_FIELDS + ("retained_block_id",)

MEMBERSHIP_COMMON_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "unit_local_index", "pos1",
    "block_id", "block_index",
)

CHECK_NAMES = (
    "normalized_input_payload_bytes_equal",
    "python_receipt_valid",
    "cpp_summary_valid",
    "headers_exact",
    "keys_unique",
    "blocks_exact",
    "membership_exact",
    "dispositions_exact",
    "aggregate_exact",
    "declared_aggregates_exact",
)

PYTHON_REQUIRED_PARAMETERS = {
    "threshold": 3,
    "max_block_size": 12,
    "accepted_inference_role": "PRIMARY_PS_AWARE",
    "accepted_linkage_basis": ["PS_HP1", "PS_HP2"],
    "phase_set_required": True,
    "fixed_linkage_calls": ["A", "R"],
    "non_linking_calls": ["D", "L", "O", "S", "X"],
}
PYTHON_REQUIRED_CHECKS = (
    "unit_sites_assigned_exactly_once",
    "block_k_at_most_12",
    "k_at_most_12_has_one_block",
    "cross_ps_zero",
    "cross_hp_zero",
    "constraint_ids_disposed_exactly_once",
    "constraint_mass_conserved",
)
CPP_REQUIRED_PARAMETERS = {
    "max_block_size": 12,
    "partition_type": "contiguous_nonoverlapping",
    "molecule_weight_type": "arbitrary_precision_nonnegative_integer",
    "objective_order": [
        "max_retained_molecule_weight",
        "max_retained_pattern_count",
        "min_blocks",
        "max_cut_gap_sum",
        "lexicographically_smaller_cut_tuple",
    ],
}
CPP_REQUIRED_CHECKS = (
    "exact_nonmissing_ps_hp_units",
    "strictly_increasing_unique_positions",
    "constraint_positions_within_unit",
    "unique_unit_and_constraint_ids",
    "max_block_size_respected",
    "site_membership_conserved",
    "constraint_dispositions_mutually_exclusive_and_conserved",
    "objective_matches_dispositions",
)

CANONICAL_NONNEGATIVE_INTEGER = re.compile(r"(?:0|[1-9][0-9]*)\Z")


class StopComparison(Exception):
    """Stop after a hard gate while preserving a diagnostic receipt."""


class ComparisonState:
    """Mutable comparison evidence collected into one deterministic receipt."""

    def __init__(self) -> None:
        self.checks: Dict[str, bool] = {name: False for name in CHECK_NAMES}
        self.artifacts: Dict[str, Any] = {
            "python": {},
            "cpp_inputs": {},
            "cpp_outputs": {},
            "normalized_payloads": {},
        }
        self.counts: Dict[str, Any] = {
            "python_rows": {},
            "cpp_rows": {},
            "aggregate": {},
        }
        self.mismatch_count = 0
        self.mismatch_samples: List[Dict[str, Any]] = []
        self.failure: Optional[str] = None

    def mismatch(
        self,
        category: str,
        key: Any,
        python_value: Any,
        cpp_value: Any,
        message: str,
    ) -> None:
        self.mismatch_count += 1
        if len(self.mismatch_samples) < MAX_MISMATCH_SAMPLES:
            self.mismatch_samples.append(
                {
                    "category": category,
                    "key": key,
                    "python": python_value,
                    "cpp": cpp_value,
                    "message": message,
                }
            )

    def receipt(self) -> Dict[str, Any]:
        all_pass = self.mismatch_count == 0 and all(self.checks.values())
        return {
            "schema_name": SCHEMA_NAME,
            "schema_version": SCHEMA_VERSION,
            "all_pass": all_pass,
            "checks": self.checks,
            "artifacts": self.artifacts,
            "counts": self.counts,
            "mismatch_count": self.mismatch_count,
            "mismatch_samples": self.mismatch_samples,
            "failure": self.failure,
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Python exact-PS normalized partition outputs with the "
            "independent C++ partition outputs."
        )
    )
    parser.add_argument("--python-dir", required=True, type=Path)
    parser.add_argument("--cpp-input-units", required=True, type=Path)
    parser.add_argument("--cpp-input-constraints", required=True, type=Path)
    parser.add_argument("--cpp-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def read_regular_file(path: Path, label: str) -> bytes:
    if not path.is_file():
        raise ValueError(f"{label} is not a regular file: {path}")
    return path.read_bytes()


def identity(path: Path, payload: Optional[bytes] = None) -> Dict[str, Any]:
    observed = read_regular_file(path, str(path)) if payload is None else payload
    return {
        "path": str(path.resolve()),
        "bytes": len(observed),
        "sha256": sha256(observed),
    }


def load_json(path: Path, label: str) -> Tuple[Dict[str, Any], bytes]:
    payload = read_regular_file(path, label)

    def reject_duplicate_keys(pairs: Iterable[Tuple[str, Any]]) -> Dict[str, Any]:
        result: Dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{label} contains duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        decoded = payload.decode("utf-8")
        document = json.loads(decoded, object_pairs_hook=reject_duplicate_keys)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot parse {label}: {error}") from error
    if not isinstance(document, dict):
        raise ValueError(f"{label} root must be a JSON object")
    return document, payload


def parse_tsv_payload(
    payload: bytes,
    expected_header: Sequence[str],
    label: str,
) -> List[Dict[str, str]]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} is not valid UTF-8: {error}") from error
    reader = csv.reader(io.StringIO(text, newline=""), delimiter="\t")
    try:
        header = next(reader)
    except StopIteration as error:
        raise ValueError(f"{label} is empty") from error
    if tuple(header) != tuple(expected_header):
        raise ValueError(
            f"{label} header mismatch; expected {list(expected_header)!r}, "
            f"observed {header!r}"
        )
    rows: List[Dict[str, str]] = []
    for line_number, fields in enumerate(reader, start=2):
        if not fields or (len(fields) == 1 and fields[0] == ""):
            raise ValueError(f"{label} contains a blank row at line {line_number}")
        if len(fields) != len(expected_header):
            raise ValueError(
                f"{label} line {line_number} has {len(fields)} fields; "
                f"expected {len(expected_header)}"
            )
        rows.append(dict(zip(expected_header, fields)))
    return rows


def read_tsv(
    path: Path,
    expected_header: Sequence[str],
    label: str,
    compressed: bool,
) -> Tuple[List[Dict[str, str]], bytes]:
    file_payload = read_regular_file(path, label)
    if compressed:
        try:
            semantic_payload = gzip.decompress(file_payload)
        except (OSError, EOFError) as error:
            raise ValueError(f"cannot decompress {label}: {error}") from error
    else:
        semantic_payload = file_payload
    return parse_tsv_payload(semantic_payload, expected_header, label), file_payload


def canonical_nonnegative_integer(value: Any, label: str) -> int:
    if not isinstance(value, str) or CANONICAL_NONNEGATIVE_INTEGER.fullmatch(value) is None:
        raise ValueError(f"{label} must be a canonical non-negative decimal string")
    return int(value)


def positive_integer(value: Any, label: str) -> int:
    parsed = canonical_nonnegative_integer(value, label)
    if parsed == 0:
        raise ValueError(f"{label} must be positive")
    return parsed


def bind_normalized_inputs(
    state: ComparisonState,
    python_dir: Path,
    cpp_units_path: Path,
    cpp_constraints_path: Path,
) -> Tuple[bytes, bytes]:
    """Bind C++ inputs byte-for-byte before any semantic result comparison."""
    specifications = (
        ("units", python_dir / "units.tsv.gz", cpp_units_path),
        ("constraints", python_dir / "constraints.tsv.gz", cpp_constraints_path),
    )
    decompressed_payloads: List[bytes] = []
    all_equal = True
    for label, gzip_path, plain_path in specifications:
        compressed = read_regular_file(gzip_path, f"Python {label}.tsv.gz")
        plain = read_regular_file(plain_path, f"C++ input {label}.tsv")
        try:
            decompressed = gzip.decompress(compressed)
        except (OSError, EOFError) as error:
            raise ValueError(f"cannot decompress Python {label}.tsv.gz: {error}") from error
        state.artifacts["python"][f"{label}.tsv.gz"] = identity(gzip_path, compressed)
        state.artifacts["cpp_inputs"][f"{label}.tsv"] = identity(plain_path, plain)
        state.artifacts["normalized_payloads"][label] = {
            "python_decompressed_bytes": len(decompressed),
            "python_decompressed_sha256": sha256(decompressed),
            "cpp_plain_bytes": len(plain),
            "cpp_plain_sha256": sha256(plain),
            "bytes_equal": decompressed == plain,
        }
        decompressed_payloads.append(decompressed)
        if decompressed != plain:
            all_equal = False
            state.mismatch(
                "input_payload",
                label,
                {
                    "bytes": len(decompressed),
                    "sha256": sha256(decompressed),
                },
                {"bytes": len(plain), "sha256": sha256(plain)},
                "C++ plain input is not byte-identical to the gzip-decompressed Python payload",
            )
    state.checks["normalized_input_payload_bytes_equal"] = all_equal
    if not all_equal:
        state.failure = "normalized input payload byte binding failed; semantic comparison skipped"
        raise StopComparison(state.failure)
    return decompressed_payloads[0], decompressed_payloads[1]


def validate_python_receipt(
    state: ComparisonState,
    python_dir: Path,
    output_identities: Mapping[str, Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    path = python_dir / "receipt.json"
    try:
        receipt, payload = load_json(path, "Python receipt")
        state.artifacts["python"]["receipt.json"] = identity(path, payload)
    except Exception as error:
        state.mismatch("python_receipt", "receipt.json", "valid JSON object", None, str(error))
        return None

    valid = True
    expected_scalars = {
        "schema_name": PYTHON_SCHEMA_NAME,
        "schema_version": PYTHON_SCHEMA_VERSION,
        "all_pass": True,
    }
    for key, expected in expected_scalars.items():
        observed = receipt.get(key)
        if observed != expected or type(observed) is not type(expected):
            valid = False
            state.mismatch(
                "python_receipt", key, expected, observed,
                "Python receipt contract value differs",
            )
    valid = validate_parameters_and_checks(
        state,
        "python_receipt",
        receipt,
        PYTHON_REQUIRED_PARAMETERS,
        PYTHON_REQUIRED_CHECKS,
    ) and valid
    outputs = receipt.get("outputs")
    if not isinstance(outputs, dict):
        valid = False
        state.mismatch(
            "python_receipt", "outputs", "object", outputs,
            "Python receipt outputs must be an object",
        )
    else:
        for name, expected_identity in sorted(output_identities.items()):
            observed = outputs.get(name)
            if observed != expected_identity:
                valid = False
                state.mismatch(
                    "python_receipt", f"outputs.{name}", expected_identity, observed,
                    "Python receipt output identity differs from current artifact",
                )
    state.checks["python_receipt_valid"] = valid
    return receipt


def validate_cpp_summary(state: ComparisonState, cpp_dir: Path) -> Optional[Dict[str, Any]]:
    path = cpp_dir / "summary.json"
    try:
        summary, payload = load_json(path, "C++ summary")
        state.artifacts["cpp_outputs"]["summary.json"] = identity(path, payload)
    except Exception as error:
        state.mismatch("cpp_summary", "summary.json", "valid JSON object", None, str(error))
        return None
    valid = True
    expected_scalars = {
        "schema_name": CPP_SCHEMA_NAME,
        "schema_version": CPP_SCHEMA_VERSION,
        "all_pass": True,
    }
    for key, expected in expected_scalars.items():
        observed = summary.get(key)
        if observed != expected or type(observed) is not type(expected):
            valid = False
            state.mismatch(
                "cpp_summary", key, expected, observed,
                "C++ summary contract value differs",
            )
    valid = validate_parameters_and_checks(
        state,
        "cpp_summary",
        summary,
        CPP_REQUIRED_PARAMETERS,
        CPP_REQUIRED_CHECKS,
    ) and valid
    state.checks["cpp_summary_valid"] = valid
    return summary


def validate_parameters_and_checks(
    state: ComparisonState,
    category: str,
    document: Mapping[str, Any],
    required_parameters: Mapping[str, Any],
    required_checks: Sequence[str],
) -> bool:
    valid = True
    parameters = document.get("parameters")
    if not isinstance(parameters, dict):
        state.mismatch(
            category, "parameters", required_parameters, parameters,
            "required algorithm parameters are absent or not an object",
        )
        valid = False
    else:
        for key, expected in required_parameters.items():
            observed = parameters.get(key)
            if observed != expected or type(observed) is not type(expected):
                state.mismatch(
                    category, f"parameters.{key}", expected, observed,
                    "algorithm parameter differs from the fixed comparison contract",
                )
                valid = False

    checks = document.get("checks")
    if not isinstance(checks, dict):
        state.mismatch(
            category, "checks", list(required_checks), checks,
            "required internal checks are absent or not an object",
        )
        valid = False
    else:
        missing = sorted(set(required_checks) - set(checks))
        for key in missing:
            state.mismatch(
                category, f"checks.{key}", True, None,
                "required internal check is missing",
            )
            valid = False
        for key, observed in sorted(checks.items()):
            if observed is not True:
                state.mismatch(
                    category, f"checks.{key}", True, observed,
                    "every declared internal check must be boolean true",
                )
                valid = False
    return valid


def key_text(key: Tuple[str, ...]) -> str:
    return "|".join(key)


def index_rows(
    state: ComparisonState,
    category: str,
    side: str,
    rows: Sequence[Mapping[str, str]],
    key_fields: Sequence[str],
) -> Tuple[Dict[Tuple[str, ...], Mapping[str, str]], bool]:
    indexed: Dict[Tuple[str, ...], Mapping[str, str]] = {}
    unique = True
    for row_number, row in enumerate(rows, start=2):
        key = tuple(row[field] for field in key_fields)
        if any(value == "" for value in key):
            unique = False
            state.mismatch(
                category,
                f"{side}:line:{row_number}",
                "non-empty key",
                list(key),
                f"{side} row contains an empty key field",
            )
            continue
        if key in indexed:
            unique = False
            state.mismatch(
                category,
                key_text(key),
                "unique key",
                side,
                f"duplicate {side} row key",
            )
            continue
        indexed[key] = row
    return indexed, unique


def compare_rows(
    state: ComparisonState,
    category: str,
    python_rows: Sequence[Mapping[str, str]],
    cpp_rows: Sequence[Mapping[str, str]],
    fields: Sequence[str],
    key_fields: Sequence[str],
) -> Tuple[bool, bool]:
    before = state.mismatch_count
    python_index, python_unique = index_rows(
        state, category, "python", python_rows, key_fields
    )
    cpp_index, cpp_unique = index_rows(state, category, "cpp", cpp_rows, key_fields)
    for key in sorted(set(python_index) | set(cpp_index)):
        python_row = python_index.get(key)
        cpp_row = cpp_index.get(key)
        if python_row is None or cpp_row is None:
            state.mismatch(
                category,
                key_text(key),
                None if python_row is None else {field: python_row[field] for field in fields},
                None if cpp_row is None else {field: cpp_row[field] for field in fields},
                "sorted key exists on only one side",
            )
            continue
        python_projection = {field: python_row[field] for field in fields}
        cpp_projection = {field: cpp_row[field] for field in fields}
        if python_projection != cpp_projection:
            state.mismatch(
                category,
                key_text(key),
                python_projection,
                cpp_projection,
                "sorted common-field row differs",
            )
    return state.mismatch_count == before, python_unique and cpp_unique


def compute_aggregate(
    state: ComparisonState,
    side: str,
    units: Sequence[Mapping[str, str]],
    constraints: Sequence[Mapping[str, str]],
    blocks: Sequence[Mapping[str, str]],
    membership: Sequence[Mapping[str, str]],
    dispositions: Sequence[Mapping[str, str]],
) -> Dict[str, Any]:
    sites = 0
    constraint_weight: Dict[Tuple[str, str, str, str], int] = {}
    constraint_pattern_count: Dict[Tuple[str, str, str, str], int] = {}
    for row in units:
        sites += positive_integer(row["k"], f"{side} unit {row['unit_id']} k")
    for row in constraints:
        key = (row["dataset"], row["chrom"], row["unit_id"], row["constraint_id"])
        constraint_weight[key] = canonical_nonnegative_integer(
            row["molecule_weight"], f"{side} constraint {row['constraint_id']} weight"
        )
        constraint_pattern_count[key] = canonical_nonnegative_integer(
            row["pattern_count"],
            f"{side} constraint {row['constraint_id']} pattern_count",
        )

    status_counts = {status: 0 for status in DISPOSITIONS}
    status_weights = {status: 0 for status in DISPOSITIONS}
    status_pattern_counts = {status: 0 for status in DISPOSITIONS}
    disposition_keys = set()
    for row in dispositions:
        key = (row["dataset"], row["chrom"], row["unit_id"], row["constraint_id"])
        disposition_keys.add(key)
        status = row["disposition"]
        if status not in status_counts:
            raise ValueError(
                f"{side} disposition {row['constraint_id']} has unknown status {status!r}"
            )
        weight = canonical_nonnegative_integer(
            row["molecule_weight"], f"{side} disposition {row['constraint_id']} weight"
        )
        pattern_count = canonical_nonnegative_integer(
            row["pattern_count"],
            f"{side} disposition {row['constraint_id']} pattern_count",
        )
        if key not in constraint_weight:
            state.mismatch(
                "aggregate", key_text(key), "known constraint", side,
                f"{side} disposition references an unknown constraint",
            )
        elif constraint_weight[key] != weight:
            state.mismatch(
                "aggregate", key_text(key), str(constraint_weight[key]), str(weight),
                f"{side} disposition weight differs from normalized constraint",
            )
        if key in constraint_pattern_count and constraint_pattern_count[key] != pattern_count:
            state.mismatch(
                "aggregate", key_text(key), str(constraint_pattern_count[key]),
                str(pattern_count),
                f"{side} disposition pattern_count differs from normalized constraint",
            )
        status_counts[status] += 1
        status_weights[status] += weight
        status_pattern_counts[status] += pattern_count

    expected_constraint_keys = set(constraint_weight)
    if disposition_keys != expected_constraint_keys:
        state.mismatch(
            "aggregate",
            f"{side}.constraint_disposition_keys",
            sorted(key_text(key) for key in expected_constraint_keys),
            sorted(key_text(key) for key in disposition_keys),
            f"{side} dispositions do not cover normalized constraints exactly once",
        )
    if len(membership) != sites:
        state.mismatch(
            "aggregate", f"{side}.sites", sites, len(membership),
            f"{side} membership row count differs from declared unit-site count",
        )

    total_weight = sum(constraint_weight.values())
    retained_weight = status_weights[RETAINED]
    cut_weight = status_weights[CUT]
    unavoidable_weight = status_weights[UNAVOIDABLE]
    lost_weight = cut_weight + unavoidable_weight
    if total_weight != retained_weight + lost_weight:
        state.mismatch(
            "aggregate", f"{side}.weight_conservation", str(total_weight),
            str(retained_weight + lost_weight),
            f"{side} total weight is not retained + cut + unavoidable",
        )
    total_pattern_count = sum(constraint_pattern_count.values())
    retained_pattern_count = status_pattern_counts[RETAINED]
    cut_pattern_count = status_pattern_counts[CUT]
    unavoidable_pattern_count = status_pattern_counts[UNAVOIDABLE]
    lost_pattern_count = cut_pattern_count + unavoidable_pattern_count
    if total_pattern_count != retained_pattern_count + lost_pattern_count:
        state.mismatch(
            "aggregate", f"{side}.pattern_count_conservation",
            str(total_pattern_count), str(retained_pattern_count + lost_pattern_count),
            f"{side} total pattern_count is not retained + cut + unavoidable",
        )
    partition_cuts = len(blocks) - len(units)
    if partition_cuts < 0:
        raise ValueError(f"{side} has fewer blocks than units")
    return {
        "units": len(units),
        "sites": sites,
        "constraints": len(constraints),
        "blocks": len(blocks),
        "partition_cuts": partition_cuts,
        "retained": status_counts[RETAINED],
        "cut": status_counts[CUT],
        "unavoidable": status_counts[UNAVOIDABLE],
        "weights": {
            "total": str(total_weight),
            "retained": str(retained_weight),
            "cut": str(cut_weight),
            "unavoidable": str(unavoidable_weight),
            "lost": str(lost_weight),
        },
        "pattern_counts": {
            "total": str(total_pattern_count),
            "retained": str(retained_pattern_count),
            "cut": str(cut_pattern_count),
            "unavoidable": str(unavoidable_pattern_count),
            "lost": str(lost_pattern_count),
        },
    }


def json_count(document: Mapping[str, Any], path: str) -> Any:
    current: Any = document
    for component in path.split("."):
        if not isinstance(current, dict) or component not in current:
            return None
        current = current[component]
    return current


def compare_declared_value(
    state: ComparisonState,
    category: str,
    path: str,
    expected: Any,
    observed: Any,
    expect_string: bool,
) -> bool:
    valid_type = (
        isinstance(observed, str)
        and CANONICAL_NONNEGATIVE_INTEGER.fullmatch(observed) is not None
        if expect_string
        else isinstance(observed, int) and not isinstance(observed, bool) and observed >= 0
    )
    if not valid_type or observed != expected:
        state.mismatch(
            category, path, expected, observed,
            "declared aggregate differs from aggregates recomputed from rows",
        )
        return False
    return True


def validate_declared_aggregates(
    state: ComparisonState,
    python_receipt: Optional[Mapping[str, Any]],
    cpp_summary: Optional[Mapping[str, Any]],
    aggregate: Mapping[str, Any],
) -> bool:
    if python_receipt is None or cpp_summary is None:
        state.mismatch(
            "declared_aggregates", "documents", "both valid JSON documents", None,
            "cannot validate declared aggregates without both documents",
        )
        return False
    checks: List[bool] = []
    python_mapping = {
        "counts.eligible_units": aggregate["units"],
        "counts.eligible_unit_sites": aggregate["sites"],
        "counts.constraints": aggregate["constraints"],
        "counts.blocks": aggregate["blocks"],
        "counts.retained_constraints": aggregate["retained"],
        "counts.cut_constraints": aggregate["cut"],
        "counts.unavoidable_constraints": aggregate["unavoidable"],
    }
    for path, expected in python_mapping.items():
        checks.append(
            compare_declared_value(
                state, "declared_aggregates", f"python.{path}", expected,
                json_count(python_receipt, path), False,
            )
        )
    for name in ("total", "retained", "cut", "unavoidable"):
        path = f"constraint_mass.{name}"
        checks.append(
            compare_declared_value(
                state, "declared_aggregates", f"python.{path}",
                aggregate["weights"][name], json_count(python_receipt, path), True,
            )
        )

    cpp_mapping = {
        "counts.units": aggregate["units"],
        "counts.sites": aggregate["sites"],
        "counts.constraints_total": aggregate["constraints"],
        "counts.blocks": aggregate["blocks"],
        "counts.cuts": aggregate["partition_cuts"],
        "counts.constraints_retained": aggregate["retained"],
        "counts.constraints_cut": aggregate["cut"],
        "counts.constraints_unavoidable": aggregate["unavoidable"],
    }
    for path, expected in cpp_mapping.items():
        checks.append(
            compare_declared_value(
                state, "declared_aggregates", f"cpp.{path}", expected,
                json_count(cpp_summary, path), False,
            )
        )
    for name in ("total", "retained", "lost"):
        path = f"weights.{name}"
        checks.append(
            compare_declared_value(
                state, "declared_aggregates", f"cpp.{path}",
                aggregate["weights"][name], json_count(cpp_summary, path), True,
            )
        )
    for name in ("total", "retained", "lost"):
        path = f"pattern_counts.{name}"
        checks.append(
            compare_declared_value(
                state, "declared_aggregates", f"cpp.{path}",
                aggregate["pattern_counts"][name], json_count(cpp_summary, path), True,
            )
        )
    return all(checks)


def run_comparison(args: argparse.Namespace, state: ComparisonState) -> None:
    python_dir = args.python_dir
    cpp_dir = args.cpp_dir

    units_payload, constraints_payload = bind_normalized_inputs(
        state, python_dir, args.cpp_input_units, args.cpp_input_constraints
    )
    python_units = parse_tsv_payload(units_payload, UNITS_FIELDS, "Python units.tsv.gz payload")
    python_constraints = parse_tsv_payload(
        constraints_payload, CONSTRAINT_FIELDS, "Python constraints.tsv.gz payload"
    )

    python_tables: Dict[str, List[Dict[str, str]]] = {
        "units": python_units,
        "constraints": python_constraints,
    }
    python_specs = (
        ("blocks", PYTHON_BLOCK_FIELDS),
        ("membership", PYTHON_MEMBERSHIP_FIELDS),
        ("dispositions", PYTHON_DISPOSITION_FIELDS),
    )
    python_output_identities: Dict[str, Mapping[str, Any]] = {
        "units.tsv.gz": state.artifacts["python"]["units.tsv.gz"],
        "constraints.tsv.gz": state.artifacts["python"]["constraints.tsv.gz"],
    }
    for name, fields in python_specs:
        path = python_dir / f"{name}.tsv.gz"
        rows, file_payload = read_tsv(path, fields, f"Python {name}.tsv.gz", True)
        python_tables[name] = rows
        file_identity = identity(path, file_payload)
        state.artifacts["python"][f"{name}.tsv.gz"] = file_identity
        python_output_identities[f"{name}.tsv.gz"] = file_identity

    cpp_tables: Dict[str, List[Dict[str, str]]] = {}
    cpp_specs = (
        ("blocks", CPP_BLOCK_FIELDS),
        ("membership", CPP_MEMBERSHIP_FIELDS),
        ("dispositions", CPP_DISPOSITION_FIELDS),
    )
    for name, fields in cpp_specs:
        path = cpp_dir / f"{name}.tsv"
        rows, file_payload = read_tsv(path, fields, f"C++ {name}.tsv", False)
        cpp_tables[name] = rows
        state.artifacts["cpp_outputs"][f"{name}.tsv"] = identity(path, file_payload)
    state.checks["headers_exact"] = True

    python_receipt = validate_python_receipt(
        state, python_dir, python_output_identities
    )
    cpp_summary = validate_cpp_summary(state, cpp_dir)

    state.counts["python_rows"] = {
        name: len(rows) for name, rows in sorted(python_tables.items())
    }
    state.counts["cpp_rows"] = {
        "blocks": len(cpp_tables["blocks"]),
        "membership": len(cpp_tables["membership"]),
        "dispositions": len(cpp_tables["dispositions"]),
        "input_units": len(python_units),
        "input_constraints": len(python_constraints),
    }

    keys_unique = True
    blocks_exact, blocks_unique = compare_rows(
        state,
        "blocks",
        python_tables["blocks"],
        cpp_tables["blocks"],
        PYTHON_BLOCK_FIELDS,
        ("dataset", "chrom", "unit_id", "block_index"),
    )
    keys_unique = keys_unique and blocks_unique
    membership_exact, membership_unique = compare_rows(
        state,
        "membership",
        python_tables["membership"],
        cpp_tables["membership"],
        MEMBERSHIP_COMMON_FIELDS,
        ("dataset", "chrom", "unit_id", "unit_local_index"),
    )
    keys_unique = keys_unique and membership_unique
    dispositions_exact, disposition_unique = compare_rows(
        state,
        "dispositions",
        python_tables["dispositions"],
        cpp_tables["dispositions"],
        PYTHON_DISPOSITION_FIELDS,
        ("dataset", "chrom", "unit_id", "constraint_id"),
    )
    keys_unique = keys_unique and disposition_unique

    # The normalized inputs also require unique stable keys.
    _, units_unique = index_rows(
        state, "units", "normalized", python_units,
        ("dataset", "chrom", "unit_id"),
    )
    _, constraints_unique = index_rows(
        state, "constraints", "normalized", python_constraints,
        ("dataset", "chrom", "unit_id", "constraint_id"),
    )
    keys_unique = keys_unique and units_unique and constraints_unique
    state.checks["keys_unique"] = keys_unique
    state.checks["blocks_exact"] = blocks_exact
    state.checks["membership_exact"] = membership_exact
    state.checks["dispositions_exact"] = dispositions_exact

    aggregate_before = state.mismatch_count
    python_aggregate = compute_aggregate(
        state,
        "python",
        python_units,
        python_constraints,
        python_tables["blocks"],
        python_tables["membership"],
        python_tables["dispositions"],
    )
    cpp_aggregate = compute_aggregate(
        state,
        "cpp",
        python_units,
        python_constraints,
        cpp_tables["blocks"],
        cpp_tables["membership"],
        cpp_tables["dispositions"],
    )
    state.counts["python_aggregate"] = python_aggregate
    state.counts["cpp_aggregate"] = cpp_aggregate
    state.counts["aggregate"] = python_aggregate
    if python_aggregate != cpp_aggregate:
        state.mismatch(
            "aggregate", "python_vs_cpp", python_aggregate, cpp_aggregate,
            "aggregates recomputed from Python and C++ rows differ",
        )
    state.checks["aggregate_exact"] = (
        state.mismatch_count == aggregate_before and python_aggregate == cpp_aggregate
    )
    state.checks["declared_aggregates_exact"] = validate_declared_aggregates(
        state, python_receipt, cpp_summary, python_aggregate
    )


def write_receipt_exclusive(path: Path, receipt: Mapping[str, Any]) -> None:
    payload = (
        json.dumps(receipt, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    ).encode("utf-8")
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        raise


def main() -> int:
    args = parse_args()
    if os.path.lexists(args.output):
        print(f"ERROR: output already exists; refusing to overwrite: {args.output}", file=sys.stderr)
        return 1

    state = ComparisonState()
    try:
        run_comparison(args, state)
    except StopComparison:
        pass
    except Exception as error:
        state.failure = str(error)
        if not state.checks["headers_exact"] and "header" in str(error).lower():
            category = "header"
        else:
            category = "contract"
        state.mismatch(category, "comparison", "valid comparison contract", None, str(error))

    receipt = state.receipt()
    try:
        write_receipt_exclusive(args.output, receipt)
    except FileExistsError:
        print(f"ERROR: output already exists; refusing to overwrite: {args.output}", file=sys.stderr)
        return 1
    except Exception as error:
        print(f"ERROR: cannot write comparison receipt: {error}", file=sys.stderr)
        return 1

    summary = {
        "output": str(args.output.resolve()),
        "all_pass": receipt["all_pass"],
        "mismatch_count": receipt["mismatch_count"],
    }
    stream = sys.stdout if receipt["all_pass"] else sys.stderr
    print(json.dumps(summary, sort_keys=True), file=stream)
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
