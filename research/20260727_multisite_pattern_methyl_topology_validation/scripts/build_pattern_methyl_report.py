#!/usr/bin/env python3
"""Build a standalone exact-raw-HP pattern × methylation technical report.

The report is an association-only reader over frozen Task-B outputs.  It
validates receipt bindings, injects all reader-facing values from the supplied
TSV/JSON artifacts, and emits one offline HTML file plus the exact embedded
report-data JSON.  It never maps partial X patterns to topology vertices, and
Hamming-distance greater than one is rendered only as a pair band.
"""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence, TextIO


SCHEMA_NAME = "intersubmod.pattern_methyl_standalone_report"
SCHEMA_VERSION = "1.1.0"
TITLE = "多 sSNV pattern × 區域甲基關聯驗證"
CLAIM_CEILING_ZH = (
    "本報告只支援 pattern-conditioned regional methylation association；"
    "不建立細胞克隆、演化次序、因果方向，也不以甲基重算 topology。"
)
RESOLUTION_LIMITED_FAMILY = "CONFIRMATORY_FULL4_OR_LONG"
RESOLUTION_LIMITED_PERMUTATIONS = 49_999
RESOLUTION_LIMITED_P = 1.0 / 50_000.0
RESOLUTION_LIMITED_Q_THRESHOLD = 0.05
ANALYSIS_SUMMARY_SCHEMAS = frozenset(
    {
        "intersubmod.pattern_methyl_evidence.summary",
        "intersubmod.pattern_methyl_evidence.combined_summary",
    }
)
SENTINEL_TARGETS = (
    {
        "id": "h2009-chr5-18096980",
        "dataset": "H2009",
        "chrom": "chr5",
        "position": 18096980,
        "hp_raw": "2-1",
        "label": "H2009 chr5:18,096,980",
    },
    {
        "id": "hcc1395-chr22-46257699",
        "dataset": "HCC1395",
        "chrom": "chr22",
        "position": 46257699,
        "hp_raw": "1-1",
        "label": "HCC1395 chr22:46,257,699",
    },
)
ASSESSMENTS = (
    "ROBUST_ASSOCIATION",
    "LOCAL_CIS_COMPATIBLE",
    "TAG_DEPENDENT",
    "CONFOUNDED",
    "EVALUABLE_NO_ROBUST_ASSOCIATION",
    "NOT_EVALUABLE",
)
ASSESSMENT_PRIORITY = {
    "ROBUST_ASSOCIATION": 0,
    "LOCAL_CIS_COMPATIBLE": 1,
    "TAG_DEPENDENT": 2,
    "CONFOUNDED": 3,
    "EVALUABLE_NO_ROBUST_ASSOCIATION": 4,
    "NOT_EVALUABLE": 5,
}
PATTERN_COUNT_REQUIRED = {
    "dataset",
    "chrom",
    "region_id",
    "unit_id",
    "phase_set",
    "hp_family",
    "hp_raw",
    "active_positions",
    "n_active_bits",
    "n_complete",
    "n_partial",
    "state_count_json",
    "partial_state_count_json",
    "formal_n5",
}
CATALOG_REQUIRED = {
    "dataset",
    "chrom",
    "position1",
    "status",
}
EVIDENCE_REQUIRED = {
    "schema_version",
    "dataset",
    "chrom",
    "region_id",
    "unit_id",
    "phase_set",
    "hp_family",
    "hp_raw",
    "active_positions",
    "n_active_bits",
    "pair_full4",
    "k_ge_3",
    "input_n_complete",
    "input_state_counts_json",
    "analysis_n",
    "analysis_state_counts_json",
    "n_common_cpg",
    "n_distal_cpg",
    "permanova_r2",
    "permanova_p",
    "permanova_permutations_requested",
    "permdisp_p",
    "best_pair",
    "best_pair_hamming",
    "best_pair_distance_contrast",
    "best_pair_standardized_effect",
    "best_pair_topology_relation",
    "max_geometry_smd",
    "all_states_n8",
    "all_states_n10",
    "equal_n_retention",
    "rarefaction_retention",
    "distal_retention",
    "multiplicity_family",
    "q_by",
    "p_holm",
    "assessment",
    "evaluation_status",
    "invalid_reason",
}
MAX_PROFILE_CPG = 72
MAX_LAYER_CASES = 8
MAX_SENTINEL_MATCHES = 3
READ_CLUSTER_MATRIX_MAX_N = 160
READ_CLUSTER_DISTANCE_LEVELS = 254
READ_CLUSTER_NA_CODE = 255
READ_CLUSTER_SYMMETRY_TOLERANCE = 1e-6
READ_CLUSTER_DIAGONAL_TOLERANCE = 1e-6
READ_DIGEST_PATTERN = re.compile(r"^[0-9a-f]{64}$")
HAMMING1_RELATIONS = frozenset(
    {
        "HAMMING1_GLOBAL_BEST_UNANIMOUS",
        "HAMMING1_NOT_UNANIMOUS",
        "HAMMING1_NOT_IN_GLOBAL_BEST",
        "HAMMING1_TOPOLOGY_UNAVAILABLE",
    }
)
HAMMING_GT1_RELATION = "PAIR_BAND_ONLY_HAMMING_GT1"


class ReportContractError(RuntimeError):
    """Raised when a report input violates its frozen contract."""


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise ReportContractError(f"{label} is missing or not a file: {path}")
    if path.stat().st_size <= 0:
        raise ReportContractError(f"{label} is empty: {path}")
    return path


def load_json(path: Path, label: str) -> dict[str, Any]:
    require_file(path, label)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ReportContractError(f"{label} is not valid JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise ReportContractError(f"{label} must contain a JSON object")
    return payload


def load_catalog_summary(path: Path) -> dict[str, Any]:
    """Read a compact summary or stop before a large pretty-printed records array."""
    require_file(path, "artifact catalog summary")
    if path.stat().st_size <= 64 * 1024 * 1024:
        return load_json(path, "artifact catalog summary")
    prefix: list[str] = []
    found_records = False
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith('  "records": ['):
                found_records = True
                break
            prefix.append(line)
            if sum(len(value) for value in prefix) > 64 * 1024 * 1024:
                raise ReportContractError(
                    "artifact catalog summary prefix exceeds 64 MiB"
                )
    if not found_records:
        raise ReportContractError(
            "large artifact catalog JSON has no top-level records boundary"
        )
    text = "".join(prefix).rstrip()
    if text.endswith(","):
        text = text[:-1]
    try:
        payload = json.loads(text + "\n}\n")
    except json.JSONDecodeError as exc:
        raise ReportContractError(
            "cannot parse artifact catalog summary before records"
        ) from exc
    if not isinstance(payload, dict):
        raise ReportContractError("artifact catalog summary must be an object")
    return payload


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def int_value(value: object) -> int | None:
    text = str(value or "").strip()
    if text == "":
        return None
    try:
        return int(text)
    except ValueError as exc:
        raise ReportContractError(f"expected integer, received {value!r}") from exc


def float_value(value: object) -> float | None:
    text = str(value or "").strip()
    if text == "":
        return None
    try:
        number = float(text)
    except ValueError as exc:
        raise ReportContractError(f"expected number, received {value!r}") from exc
    return number if math.isfinite(number) else None


def parse_positions(value: object) -> list[int]:
    text = str(value or "").strip()
    if not text:
        return []
    try:
        positions = [int(token) for token in text.split(",") if token]
    except ValueError as exc:
        raise ReportContractError(f"invalid active_positions: {value!r}") from exc
    if (
        positions != sorted(positions)
        or len(positions) != len(set(positions))
        or any(position <= 0 for position in positions)
    ):
        raise ReportContractError(
            f"active_positions must be positive, sorted, and unique: {value!r}"
        )
    return positions


def parse_count_json(value: object, label: str) -> dict[str, int]:
    text = str(value or "{}").strip() or "{}"
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ReportContractError(f"{label} is not valid JSON") from exc
    if not isinstance(payload, dict):
        raise ReportContractError(f"{label} must be a JSON object")
    output: dict[str, int] = {}
    for key, raw_count in payload.items():
        if not isinstance(key, str) or not key:
            raise ReportContractError(f"{label} has an invalid state key")
        if isinstance(raw_count, bool) or not isinstance(raw_count, int) or raw_count < 0:
            raise ReportContractError(f"{label}[{key!r}] must be a non-negative integer")
        output[key] = raw_count
    return dict(sorted(output.items()))


def validate_pattern_states(
    states: Mapping[str, int],
    *,
    n_active_bits: int,
    complete: bool,
    label: str,
) -> None:
    for state in states:
        if len(state) != n_active_bits or any(code not in {"R", "A", "X"} for code in state):
            raise ReportContractError(
                f"{label} contains a state outside the R/A/X k={n_active_bits} contract"
            )
        if complete and "X" in state:
            raise ReportContractError(f"{label} contains X in a complete state")
        if not complete and "X" not in state:
            raise ReportContractError(
                f"{label} contains a non-X state in partial evidence"
            )


def require_fields(
    path: Path, fieldnames: Sequence[str] | None, required: Iterable[str]
) -> None:
    available = set(fieldnames or ())
    missing = sorted(set(required) - available)
    if missing:
        raise ReportContractError(f"{path}: missing columns {missing}")


def verify_bound_output(
    receipt: Mapping[str, Any], output_label: str, path: Path, receipt_label: str
) -> None:
    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping):
        raise ReportContractError(f"{receipt_label} has no outputs mapping")
    binding = outputs.get(output_label)
    if not isinstance(binding, Mapping):
        raise ReportContractError(
            f"{receipt_label} has no {output_label!r} output binding"
        )
    expected_sha = binding.get("sha256")
    if not isinstance(expected_sha, str) or len(expected_sha) != 64:
        raise ReportContractError(
            f"{receipt_label} {output_label!r} has no valid SHA-256"
        )
    actual_sha = sha256_path(path)
    if actual_sha != expected_sha:
        raise ReportContractError(
            f"{receipt_label} SHA mismatch for {output_label}: "
            f"{actual_sha} != {expected_sha}"
        )


def identity_key(row: Mapping[str, object]) -> tuple[str, str, str, str]:
    return (
        str(row["dataset"]),
        str(row["chrom"]),
        str(row["region_id"]),
        str(row["hp_raw"]),
    )


def case_id_for(key: tuple[str, str, str, str]) -> str:
    return hashlib.sha256("\x1f".join(key).encode("utf-8")).hexdigest()[:16]


def load_pattern_counts(
    path: Path, census_receipt: Mapping[str, Any]
) -> tuple[
    dict[tuple[str, str, str, str], dict[str, Any]],
    dict[str, Any],
    list[dict[str, Any]],
]:
    require_file(path, "pattern counts")
    verify_bound_output(census_receipt, "pattern_counts", path, "census receipt")
    formal: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    hp_raw_values: set[str] = set()
    rows_total = 0
    formal_n5 = 0
    complete_total = 0
    partial_total = 0
    sentinel_matches: list[dict[str, Any]] = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, PATTERN_COUNT_REQUIRED)
        for line_number, row in enumerate(reader, start=2):
            rows_total += 1
            key = identity_key(row)
            positions = parse_positions(row["active_positions"])
            n_active_bits = int_value(row["n_active_bits"])
            if n_active_bits is None or n_active_bits != len(positions):
                raise ReportContractError(
                    f"{path}:{line_number}: n_active_bits/active_positions mismatch"
                )
            n_complete = int_value(row["n_complete"])
            n_partial = int_value(row["n_partial"])
            if n_complete is None or n_complete < 0 or n_partial is None or n_partial < 0:
                raise ReportContractError(
                    f"{path}:{line_number}: invalid complete/partial counts"
                )
            complete_states = parse_count_json(
                row["state_count_json"], f"{path}:{line_number}:state_count_json"
            )
            partial_states = parse_count_json(
                row["partial_state_count_json"],
                f"{path}:{line_number}:partial_state_count_json",
            )
            validate_pattern_states(
                complete_states,
                n_active_bits=n_active_bits,
                complete=True,
                label=f"{path}:{line_number}:state_count_json",
            )
            validate_pattern_states(
                partial_states,
                n_active_bits=n_active_bits,
                complete=False,
                label=f"{path}:{line_number}:partial_state_count_json",
            )
            if sum(complete_states.values()) != n_complete:
                raise ReportContractError(
                    f"{path}:{line_number}: complete state mass mismatch"
                )
            if sum(partial_states.values()) != n_partial:
                raise ReportContractError(
                    f"{path}:{line_number}: partial state mass mismatch"
                )
            normalized = {
                "dataset": row["dataset"],
                "chrom": row["chrom"],
                "region_id": row["region_id"],
                "unit_id": row["unit_id"],
                "phase_set": row["phase_set"],
                "hp_family": row["hp_family"],
                "hp_raw": row["hp_raw"],
                "active_positions": positions,
                "n_active_bits": n_active_bits,
                "n_complete": n_complete,
                "n_partial": n_partial,
                "state_counts": complete_states,
                "partial_state_counts": partial_states,
                "formal_n5": bool_value(row["formal_n5"]),
            }
            hp_raw_values.add(row["hp_raw"])
            complete_total += n_complete
            partial_total += n_partial
            if normalized["formal_n5"]:
                if key in formal:
                    raise ReportContractError(f"duplicate formal count key: {key}")
                formal[key] = normalized
                formal_n5 += 1
            for target in SENTINEL_TARGETS:
                if (
                    row["dataset"] == target["dataset"]
                    and row["chrom"] == target["chrom"]
                    and int(target["position"]) in positions
                ):
                    sentinel_matches.append(normalized)
    receipt_counts = census_receipt.get("counts")
    if not isinstance(receipt_counts, Mapping):
        raise ReportContractError("census receipt has no counts mapping")
    for field, actual in (
        ("pattern_count_rows", rows_total),
        ("formal_n5", formal_n5),
    ):
        expected = receipt_counts.get(field)
        if expected is not None and int(expected) != actual:
            raise ReportContractError(
                f"census receipt {field} mismatch: {expected} != {actual}"
            )
    aggregate = {
        "pattern_count_rows": rows_total,
        "formal_n5": formal_n5,
        "complete_read_projections": complete_total,
        "partial_read_projections": partial_total,
        "hp_raw_values": sorted(hp_raw_values),
    }
    return formal, aggregate, sentinel_matches


def load_catalog(
    summary_path: Path,
    tsv_path: Path,
    expected_marker_keys: set[tuple[str, str, int]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    summary = load_catalog_summary(summary_path)
    if summary.get("schema_name") != "intersubmod.multisite_pattern_methyl_artifact_catalog":
        raise ReportContractError("unexpected artifact catalog schema_name")
    if summary.get("schema_version") != "1.0.0":
        raise ReportContractError("unexpected artifact catalog schema_version")
    if summary.get("pass") is not True or summary.get("status") != "PASS":
        raise ReportContractError("artifact catalog is not PASS")
    require_file(tsv_path, "artifact catalog TSV")
    rows = 0
    passed = 0
    failed = 0
    by_dataset: Counter[str] = Counter()
    by_chrom: Counter[str] = Counter()
    pass_keys: set[tuple[str, str, int]] = set()
    with open_text(tsv_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(tsv_path, reader.fieldnames, CATALOG_REQUIRED)
        for line_number, row in enumerate(reader, start=2):
            rows += 1
            try:
                position = int(row["position1"])
            except ValueError as exc:
                raise ReportContractError(
                    f"{tsv_path}:{line_number}: invalid position1"
                ) from exc
            if position <= 0:
                raise ReportContractError(
                    f"{tsv_path}:{line_number}: non-positive position1"
                )
            status = row["status"]
            if status == "PASS":
                marker_key = (row["dataset"], row["chrom"], position)
                if marker_key in pass_keys:
                    raise ReportContractError(
                        f"{tsv_path}:{line_number}: duplicate PASS marker key "
                        f"{marker_key}"
                    )
                pass_keys.add(marker_key)
                passed += 1
            else:
                failed += 1
            by_dataset[row["dataset"]] += 1
            by_chrom[row["chrom"]] += 1
    declared = summary.get("summary")
    if not isinstance(declared, Mapping):
        raise ReportContractError("artifact catalog has no summary")
    comparisons = {
        "markers_total": rows,
        "markers_pass": passed,
        "markers_fail": failed,
    }
    for field, actual in comparisons.items():
        if int(declared.get(field, -1)) != actual:
            raise ReportContractError(
                f"artifact catalog {field} mismatch: {declared.get(field)} != {actual}"
            )
    if failed:
        raise ReportContractError(f"artifact catalog TSV contains {failed} failed rows")
    if pass_keys != expected_marker_keys:
        missing = sorted(
            expected_marker_keys - pass_keys,
            key=lambda key: (key[0], natural_chrom_key(key[1]), key[2]),
        )
        unexpected = sorted(
            pass_keys - expected_marker_keys,
            key=lambda key: (key[0], natural_chrom_key(key[1]), key[2]),
        )
        raise ReportContractError(
            "artifact catalog PASS marker set does not exactly match the formal "
            f"marker universe: missing={missing[:5]!r} "
            f"unexpected={unexpected[:5]!r}"
        )
    return summary, {
        **comparisons,
        "dataset_marker_counts": dict(sorted(by_dataset.items())),
        "chromosome_marker_counts": dict(
            sorted(by_chrom.items(), key=lambda item: natural_chrom_key(item[0]))
        ),
    }


def load_bernoulli_parity(
    path: Path,
    catalog_tsv_path: Path,
    catalog_aggregate: Mapping[str, Any],
) -> dict[str, Any]:
    summary = load_json(path, "Bernoulli parity summary")
    if summary.get("schema_name") != "intersubmod.bernoulli_artifact_parity":
        raise ReportContractError("unexpected Bernoulli parity schema_name")
    if summary.get("schema_version") != "1.0.0":
        raise ReportContractError("unexpected Bernoulli parity schema_version")
    if summary.get("all_pass") is not True:
        raise ReportContractError("Bernoulli parity summary is not all_pass")

    counts = summary.get("counts")
    if not isinstance(counts, Mapping):
        raise ReportContractError("Bernoulli parity summary has no counts mapping")

    def count_value(field: str) -> int:
        value = counts.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ReportContractError(
                f"Bernoulli parity {field} must be a non-negative integer"
            )
        return value

    markers_total = count_value("markers_total")
    markers_pass = count_value("markers_pass")
    markers_fail = count_value("markers_fail")
    pair_cells_checked = count_value("pair_cells_checked")
    invalid_mask_mismatches = count_value("invalid_mask_mismatches")
    catalog_total = int(catalog_aggregate["markers_total"])
    catalog_pass = int(catalog_aggregate["markers_pass"])
    if markers_total != catalog_total or markers_total != catalog_pass:
        raise ReportContractError(
            "Bernoulli parity markers_total does not match catalog "
            f"markers_total/markers_pass: {markers_total} != "
            f"{catalog_total}/{catalog_pass}"
        )
    if markers_pass != markers_total or markers_fail != 0:
        raise ReportContractError(
            "Bernoulli parity PASS/FAIL counts contradict all_pass"
        )
    if invalid_mask_mismatches != 0:
        raise ReportContractError(
            "Bernoulli parity invalid-mask mismatches contradict all_pass"
        )

    contract = summary.get("contract")
    if not isinstance(contract, Mapping):
        raise ReportContractError("Bernoulli parity summary has no contract mapping")
    max_reads = contract.get("max_reads_per_marker")
    min_common = contract.get("min_common_cpg")
    tolerance = contract.get("absolute_tolerance")
    if isinstance(max_reads, bool) or not isinstance(max_reads, int) or max_reads < 2:
        raise ReportContractError(
            "Bernoulli parity max_reads_per_marker must be at least two"
        )
    if isinstance(min_common, bool) or not isinstance(min_common, int) or min_common < 1:
        raise ReportContractError(
            "Bernoulli parity min_common_cpg must be a positive integer"
        )
    if contract.get("all_cpgs_retained") is not True:
        raise ReportContractError("Bernoulli parity did not retain all CpGs")
    if isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)):
        raise ReportContractError(
            "Bernoulli parity absolute_tolerance must be numeric"
        )
    tolerance_float = float(tolerance)
    if not math.isfinite(tolerance_float) or tolerance_float <= 0.0:
        raise ReportContractError(
            "Bernoulli parity absolute_tolerance must be finite and positive"
        )
    max_error = summary.get("max_absolute_error")
    if isinstance(max_error, bool) or not isinstance(max_error, (int, float)):
        raise ReportContractError(
            "Bernoulli parity max_absolute_error must be numeric"
        )
    max_error_float = float(max_error)
    if (
        not math.isfinite(max_error_float)
        or max_error_float < 0.0
        or max_error_float > tolerance_float
    ):
        raise ReportContractError(
            "Bernoulli parity max_absolute_error exceeds its declared tolerance"
        )

    inputs = summary.get("inputs")
    catalog_binding = (
        inputs.get("artifact_catalog") if isinstance(inputs, Mapping) else None
    )
    if not isinstance(catalog_binding, Mapping):
        raise ReportContractError(
            "Bernoulli parity summary has no artifact catalog binding"
        )
    declared_sha = str(catalog_binding.get("sha256", ""))
    observed_sha = sha256_path(catalog_tsv_path)
    if declared_sha != observed_sha:
        raise ReportContractError(
            "Bernoulli parity artifact catalog SHA mismatch"
        )

    return {
        "all_pass": True,
        "markers_total": markers_total,
        "markers_pass": markers_pass,
        "markers_fail": markers_fail,
        "pair_cells_checked": pair_cells_checked,
        "invalid_mask_mismatches": invalid_mask_mismatches,
        "max_absolute_error": max_error_float,
        "absolute_tolerance": tolerance_float,
        "max_reads_per_marker": max_reads,
        "all_cpgs_retained": True,
        "min_common_cpg": min_common,
    }


def natural_chrom_key(chrom: str) -> tuple[int, str]:
    suffix = chrom.removeprefix("chr")
    return (int(suffix), chrom) if suffix.isdigit() else (10**9, chrom)


def normalize_case(
    row: Mapping[str, str], count_row: Mapping[str, Any]
) -> dict[str, Any]:
    n_active_bits = int_value(row["n_active_bits"])
    if n_active_bits is None or n_active_bits < 2:
        raise ReportContractError("formal evidence must have n_active_bits >= 2")
    analysis_states = parse_count_json(
        row["analysis_state_counts_json"], "analysis_state_counts_json"
    )
    input_states = parse_count_json(
        row["input_state_counts_json"], "input_state_counts_json"
    )
    validate_pattern_states(
        input_states,
        n_active_bits=n_active_bits,
        complete=True,
        label="input_state_counts_json",
    )
    validate_pattern_states(
        analysis_states,
        n_active_bits=n_active_bits,
        complete=True,
        label="analysis_state_counts_json",
    )
    states = analysis_states or input_states
    analysis_n = int_value(row["analysis_n"])
    if row["evaluation_status"] == "EVALUABLE":
        if analysis_n is None or analysis_n != sum(analysis_states.values()):
            raise ReportContractError(
                "analysis_n must equal analysis_state_counts_json mass for "
                f"evaluable evidence {identity_key(row)}"
            )
        if not analysis_states:
            raise ReportContractError(
                f"evaluable evidence has no analysis states: {identity_key(row)}"
            )
    elif analysis_states or analysis_n not in (None, 0):
        raise ReportContractError(
            "non-evaluable evidence must not contain analyzed state mass: "
            f"{identity_key(row)}"
        )
    key = identity_key(row)
    permanova_p = float_value(row["permanova_p"])
    permutations_requested = int_value(
        row["permanova_permutations_requested"]
    )
    q_by = float_value(row["q_by"])
    resolution_limited = is_resolution_limited(
        family=row["multiplicity_family"],
        permutations_requested=permutations_requested,
        p_value=permanova_p,
        q_by=q_by,
    )
    return {
        "id": case_id_for(key),
        "dataset": row["dataset"],
        "chrom": row["chrom"],
        "region_id": row["region_id"],
        "unit_id": row["unit_id"],
        "phase_set": row["phase_set"],
        "hp_family": row["hp_family"],
        "hp_raw": row["hp_raw"],
        "active_positions": parse_positions(row["active_positions"]),
        "n_active_bits": n_active_bits,
        "pair_full4": bool_value(row["pair_full4"]),
        "k_ge_3": bool_value(row["k_ge_3"]),
        "input_n_complete": int_value(row["input_n_complete"]),
        "input_state_counts": input_states,
        "partial_state_counts": dict(count_row.get("partial_state_counts", {})),
        "analysis_n": analysis_n,
        "analysis_state_counts": analysis_states,
        "patterns": sorted(states),
        "n_common_cpg": int_value(row["n_common_cpg"]),
        "n_distal_cpg": int_value(row["n_distal_cpg"]),
        "permanova_r2": float_value(row["permanova_r2"]),
        "permanova_p": permanova_p,
        "permanova_permutations_requested": permutations_requested,
        "permdisp_p": float_value(row["permdisp_p"]),
        "best_pair": row["best_pair"],
        "best_pair_hamming": int_value(row["best_pair_hamming"]),
        "best_pair_distance_contrast": float_value(
            row["best_pair_distance_contrast"]
        ),
        "best_pair_standardized_effect": float_value(
            row["best_pair_standardized_effect"]
        ),
        "best_pair_topology_relation": row["best_pair_topology_relation"],
        "max_geometry_smd": float_value(row["max_geometry_smd"]),
        "all_states_n8": bool_value(row["all_states_n8"]),
        "all_states_n10": bool_value(row["all_states_n10"]),
        "equal_n_retention": float_value(row["equal_n_retention"]),
        "rarefaction_retention": float_value(row["rarefaction_retention"]),
        "distal_retention": float_value(row["distal_retention"]),
        "multiplicity_family": row["multiplicity_family"],
        "q_by": q_by,
        "p_holm": float_value(row["p_holm"]),
        "resolution_limited": resolution_limited,
        "assessment": row["assessment"],
        "evaluation_status": row["evaluation_status"],
        "invalid_reason": row["invalid_reason"],
        "detail": None,
    }


def is_resolution_limited(
    *,
    family: str,
    permutations_requested: int | None,
    p_value: float | None,
    q_by: float | None,
) -> bool:
    return (
        family == RESOLUTION_LIMITED_FAMILY
        and permutations_requested == RESOLUTION_LIMITED_PERMUTATIONS
        and p_value == RESOLUTION_LIMITED_P
        and q_by is not None
        and q_by > RESOLUTION_LIMITED_Q_THRESHOLD
    )


def load_evidence(
    path: Path,
    analysis_summary: Mapping[str, Any],
    formal_counts: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
) -> list[dict[str, Any]]:
    require_file(path, "merged evidence")
    verify_bound_output(analysis_summary, "evidence", path, "analysis summary")
    cases: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str]] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, EVIDENCE_REQUIRED)
        for line_number, row in enumerate(reader, start=2):
            if row["schema_version"] != "1.0.0":
                raise ReportContractError(
                    f"{path}:{line_number}: unexpected evidence schema_version"
                )
            key = identity_key(row)
            if key in seen:
                raise ReportContractError(f"duplicate evidence key: {key}")
            seen.add(key)
            if key not in formal_counts:
                raise ReportContractError(
                    f"{path}:{line_number}: evidence key is absent from formal counts"
                )
            count_row = formal_counts[key]
            for field in ("unit_id", "phase_set", "hp_family"):
                if row[field] != str(count_row[field]):
                    raise ReportContractError(
                        f"{path}:{line_number}: evidence/count {field} mismatch "
                        f"for {key}: {row[field]!r} != {count_row[field]!r}"
                    )
            evidence_positions = parse_positions(row["active_positions"])
            if evidence_positions != list(count_row["active_positions"]):
                raise ReportContractError(
                    f"{path}:{line_number}: evidence/count active_positions "
                    f"mismatch for {key}"
                )
            evidence_k = int_value(row["n_active_bits"])
            if evidence_k != int(count_row["n_active_bits"]):
                raise ReportContractError(
                    f"{path}:{line_number}: evidence/count n_active_bits "
                    f"mismatch for {key}"
                )
            evidence_input_n = int_value(row["input_n_complete"])
            if evidence_input_n != int(count_row["n_complete"]):
                raise ReportContractError(
                    f"{path}:{line_number}: evidence/count input_n_complete "
                    f"mismatch for {key}"
                )
            evidence_input_states = parse_count_json(
                row["input_state_counts_json"],
                f"{path}:{line_number}:input_state_counts_json",
            )
            if evidence_input_states != dict(count_row["state_counts"]):
                raise ReportContractError(
                    f"{path}:{line_number}: evidence/count input_state_counts "
                    f"mismatch for {key}"
                )
            if row["assessment"] not in ASSESSMENTS:
                raise ReportContractError(
                    f"{path}:{line_number}: unsupported assessment "
                    f"{row['assessment']!r}"
                )
            cases.append(normalize_case(row, formal_counts[key]))
    missing_evidence = sorted(set(formal_counts) - seen)
    if missing_evidence:
        raise ReportContractError(
            f"formal counts have {len(missing_evidence)} keys absent from evidence"
        )
    return sorted(
        cases,
        key=lambda row: (
            row["dataset"],
            natural_chrom_key(row["chrom"]),
            row["region_id"],
            row["hp_raw"],
        ),
    )


def evenly_spaced_indices(length: int, limit: int) -> list[int]:
    if length <= limit:
        return list(range(length))
    indices = {
        round(index * (length - 1) / (limit - 1))
        for index in range(limit)
    }
    return sorted(indices)


def pattern_hamming(first: str, second: str) -> int:
    if len(first) != len(second):
        raise ReportContractError(
            f"cannot compare patterns of unequal width: {first!r}, {second!r}"
        )
    return sum(left != right for left, right in zip(first, second))


def encode_u8(values: Sequence[int]) -> str:
    if any(
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < 0
        or value > 255
        for value in values
    ):
        raise ReportContractError("uint8 payload contains a value outside 0..255")
    return base64.b64encode(bytes(values)).decode("ascii")


def encode_categorical(values: Sequence[str]) -> tuple[list[str], str]:
    levels = sorted(set(values))
    if len(levels) > 256:
        raise ReportContractError(
            f"categorical read annotation has {len(levels)} levels; maximum is 256"
        )
    level_index = {value: index for index, value in enumerate(levels)}
    return levels, encode_u8([level_index[value] for value in values])


def pseudonymize_read_group_levels(levels: Sequence[str]) -> list[str]:
    """Keep RG category membership while omitting run-level source strings."""
    labels: list[str] = []
    next_ordinal = 1
    for level in levels:
        if level == ".":
            labels.append(".")
        else:
            labels.append(str(next_ordinal))
            next_ordinal += 1
    return labels


def validate_read_order(
    raw: object,
    *,
    case: Mapping[str, Any],
    state_counts: Mapping[str, int],
) -> list[dict[str, str]]:
    if not isinstance(raw, list):
        raise ReportContractError("detail read_order must be a list")
    analysis_n = int(case["analysis_n"])
    if len(raw) != analysis_n:
        raise ReportContractError(
            f"detail read_order length mismatch for {identity_key(case)}: "
            f"{len(raw)} != {analysis_n}"
        )
    reads: list[dict[str, str]] = []
    qnames: set[str] = set()
    observed_patterns: Counter[str] = Counter()
    for index, item in enumerate(raw):
        if not isinstance(item, Mapping):
            raise ReportContractError(
                f"detail read_order[{index}] must be an object"
            )
        qname = item.get("qname_sha256")
        pattern = item.get("pattern")
        strand = item.get("strand")
        read_group = item.get("read_group")
        if not isinstance(qname, str) or READ_DIGEST_PATTERN.fullmatch(qname) is None:
            raise ReportContractError(
                f"detail read_order[{index}] has an invalid qname SHA-256"
            )
        if qname in qnames:
            raise ReportContractError(
                f"detail read_order contains duplicate qname digest for "
                f"{identity_key(case)}"
            )
        qnames.add(qname)
        if (
            not isinstance(pattern, str)
            or len(pattern) != int(case["n_active_bits"])
            or any(code not in {"R", "A"} for code in pattern)
            or pattern not in state_counts
        ):
            raise ReportContractError(
                f"detail read_order[{index}] has an invalid complete pattern"
            )
        if not isinstance(strand, str) or not strand:
            raise ReportContractError(
                f"detail read_order[{index}] has an invalid strand"
            )
        if not isinstance(read_group, str) or not read_group:
            raise ReportContractError(
                f"detail read_order[{index}] has an invalid read_group"
            )
        observed_patterns[pattern] += 1
        reads.append(
            {
                "pattern": pattern,
                "strand": strand,
                "read_group": read_group,
            }
        )
    if dict(sorted(observed_patterns.items())) != dict(sorted(state_counts.items())):
        raise ReportContractError(
            f"detail read_order/state count mismatch for {identity_key(case)}"
        )
    return reads


def validate_distance_matrix(
    raw: object,
    *,
    n_reads: int,
    case: Mapping[str, Any],
) -> list[list[float]] | None:
    if raw is None:
        if n_reads <= READ_CLUSTER_MATRIX_MAX_N:
            raise ReportContractError(
                f"detail distance_matrix is null at N={n_reads} for "
                f"{identity_key(case)}"
            )
        return None
    if n_reads > READ_CLUSTER_MATRIX_MAX_N:
        raise ReportContractError(
            f"detail embeds distance_matrix above N={READ_CLUSTER_MATRIX_MAX_N} "
            f"for {identity_key(case)}"
        )
    if not isinstance(raw, list) or len(raw) != n_reads:
        raise ReportContractError(
            f"detail distance_matrix height mismatch for {identity_key(case)}"
        )
    matrix: list[list[float]] = []
    for row_index, raw_row in enumerate(raw):
        if not isinstance(raw_row, list) or len(raw_row) != n_reads:
            raise ReportContractError(
                f"detail distance_matrix width mismatch at row {row_index} for "
                f"{identity_key(case)}"
            )
        row: list[float] = []
        for column_index, raw_value in enumerate(raw_row):
            if (
                isinstance(raw_value, bool)
                or not isinstance(raw_value, (int, float))
            ):
                raise ReportContractError(
                    f"detail distance_matrix[{row_index}][{column_index}] "
                    "must be numeric"
                )
            value = float(raw_value)
            if not math.isfinite(value) or value < 0.0 or value > 1.0:
                raise ReportContractError(
                    f"detail distance_matrix[{row_index}][{column_index}] "
                    "must be finite and within 0..1"
                )
            row.append(value)
        matrix.append(row)
    for index in range(n_reads):
        if abs(matrix[index][index]) > READ_CLUSTER_DIAGONAL_TOLERANCE:
            raise ReportContractError(
                f"detail distance_matrix has non-zero diagonal for "
                f"{identity_key(case)}"
            )
        for other in range(index + 1, n_reads):
            if (
                abs(matrix[index][other] - matrix[other][index])
                > READ_CLUSTER_SYMMETRY_TOLERANCE
            ):
                raise ReportContractError(
                    f"detail distance_matrix is asymmetric for "
                    f"{identity_key(case)}"
                )
    return matrix


def average_linkage_order(
    matrix: Sequence[Sequence[float]],
) -> tuple[list[int], list[list[int | float]]]:
    """Return deterministic label-independent UPGMA leaf order and linkage."""
    n_reads = len(matrix)
    if n_reads < 2:
        raise ReportContractError("UPGMA read clustering requires at least two reads")
    active = set(range(n_reads))
    sizes = {index: 1 for index in active}
    minimum_leaf = {index: index for index in active}
    leaf_orders = {index: [index] for index in active}
    distances: dict[tuple[int, int], float] = {
        (first, second): float(matrix[first][second])
        for first in range(n_reads)
        for second in range(first + 1, n_reads)
    }
    merges: list[list[int | float]] = []
    next_cluster = n_reads
    while len(active) > 1:
        if not distances:
            raise ReportContractError("UPGMA distance map became empty")
        (first, second), height = min(
            distances.items(),
            key=lambda item: (
                item[1],
                min(minimum_leaf[item[0][0]], minimum_leaf[item[0][1]]),
                max(minimum_leaf[item[0][0]], minimum_leaf[item[0][1]]),
                item[0],
            ),
        )
        remaining = sorted(active - {first, second})
        first_size = sizes[first]
        second_size = sizes[second]
        new_distances: dict[int, float] = {}
        for other in remaining:
            first_key = tuple(sorted((first, other)))
            second_key = tuple(sorted((second, other)))
            new_distances[other] = (
                first_size * distances[first_key]
                + second_size * distances[second_key]
            ) / (first_size + second_size)

        first_order = leaf_orders[first]
        second_order = leaf_orders[second]
        first_variants = (first_order, list(reversed(first_order)))
        second_variants = (second_order, list(reversed(second_order)))
        candidates: list[tuple[float, tuple[int, ...], list[int]]] = []
        for left in first_variants:
            for right in second_variants:
                candidates.append(
                    (
                        float(matrix[left[-1]][right[0]]),
                        tuple(left + right),
                        left + right,
                    )
                )
                candidates.append(
                    (
                        float(matrix[right[-1]][left[0]]),
                        tuple(right + left),
                        right + left,
                    )
                )
        chosen_order = min(candidates, key=lambda item: (item[0], item[1]))[2]

        stale_keys = [
            key for key in distances if first in key or second in key
        ]
        for key in stale_keys:
            del distances[key]
        active.remove(first)
        active.remove(second)
        active.add(next_cluster)
        sizes[next_cluster] = first_size + second_size
        minimum_leaf[next_cluster] = min(
            minimum_leaf[first], minimum_leaf[second]
        )
        leaf_orders[next_cluster] = chosen_order
        for other, value in new_distances.items():
            distances[tuple(sorted((next_cluster, other)))] = value
        merges.append(
            [
                first,
                second,
                round(float(height), 6),
                first_size + second_size,
            ]
        )
        next_cluster += 1
    root = next(iter(active))
    return leaf_orders[root], merges


def normalize_read_cluster(
    payload: Mapping[str, Any],
    *,
    case: Mapping[str, Any],
    state_counts: Mapping[str, int],
) -> dict[str, Any]:
    reads = validate_read_order(
        payload.get("read_order"),
        case=case,
        state_counts=state_counts,
    )
    matrix = validate_distance_matrix(
        payload.get("distance_matrix"),
        n_reads=len(reads),
        case=case,
    )
    if matrix is None:
        return {
            "available": False,
            "source_status": "SOURCE_NULL_N_GT_160",
            "n_reads": len(reads),
            "source_limit_n": READ_CLUSTER_MATRIX_MAX_N,
        }

    order, raw_merges = average_linkage_order(matrix)
    order_position = {source_index: index for index, source_index in enumerate(order)}
    n_reads = len(order)
    reordered_matrix = [
        matrix[source_row][source_column]
        for source_row in order
        for source_column in order
    ]
    matrix_codes = [
        min(
            READ_CLUSTER_DISTANCE_LEVELS,
            max(0, round(value * READ_CLUSTER_DISTANCE_LEVELS)),
        )
        for value in reordered_matrix
    ]
    matrix_bytes = bytes(matrix_codes)
    ordered_patterns = [reads[index]["pattern"] for index in order]
    ordered_strands = [reads[index]["strand"] for index in order]
    ordered_read_groups = [reads[index]["read_group"] for index in order]
    pattern_levels, pattern_codes = encode_categorical(ordered_patterns)
    strand_levels, strand_codes = encode_categorical(ordered_strands)
    raw_read_group_levels, read_group_codes = encode_categorical(
        ordered_read_groups
    )
    read_group_levels = pseudonymize_read_group_levels(raw_read_group_levels)

    remapped_merges: list[list[int | float]] = []
    for left, right, height, size in raw_merges:
        left_id = int(left)
        right_id = int(right)
        remapped_merges.append(
            [
                order_position[left_id] if left_id < n_reads else left_id,
                order_position[right_id] if right_id < n_reads else right_id,
                height,
                size,
            ]
        )

    within: list[float] = []
    between: list[float] = []
    for first in range(n_reads):
        for second in range(first + 1, n_reads):
            target = (
                within
                if reads[first]["pattern"] == reads[second]["pattern"]
                else between
            )
            target.append(matrix[first][second])
    within_mean = sum(within) / len(within) if within else None
    between_mean = sum(between) / len(between) if between else None
    block_count = 1 + sum(
        first != second
        for first, second in zip(ordered_patterns, ordered_patterns[1:])
    )
    adjacent_same_fraction = (
        sum(
            first == second
            for first, second in zip(ordered_patterns, ordered_patterns[1:])
        )
        / (n_reads - 1)
    )
    return {
        "available": True,
        "source_status": "EMBEDDED_QUANTIZED",
        "n_reads": n_reads,
        "order_method": "UPGMA_AVERAGE_LINKAGE_LABEL_INDEPENDENT_V1",
        "order_uses_pattern_labels": False,
        "order_tie_breaker": "SOURCE_LEAF_ORDINAL_FOR_EXACT_DISTANCE_TIES",
        "matrix_encoding": "UINT8_ROW_MAJOR_BASE64_V1",
        "matrix_u8_b64": base64.b64encode(matrix_bytes).decode("ascii"),
        "matrix_u8_sha256": hashlib.sha256(matrix_bytes).hexdigest(),
        "distance_levels": READ_CLUSTER_DISTANCE_LEVELS,
        "na_code": READ_CLUSTER_NA_CODE,
        "quantization_max_abs_error": 0.5 / READ_CLUSTER_DISTANCE_LEVELS,
        "pattern_levels": pattern_levels,
        "pattern_codes_u8_b64": pattern_codes,
        "strand_levels": strand_levels,
        "strand_codes_u8_b64": strand_codes,
        "read_group_levels": read_group_levels,
        "read_group_codes_u8_b64": read_group_codes,
        "read_group_levels_pseudonymized": True,
        "dendrogram_merges": remapped_merges,
        "within_pattern_mean": within_mean,
        "between_pattern_mean": between_mean,
        "between_minus_within": (
            between_mean - within_mean
            if within_mean is not None and between_mean is not None
            else None
        ),
        "pattern_block_count": block_count,
        "adjacent_same_pattern_fraction": adjacent_same_fraction,
    }


def normalize_detail(
    payload: Mapping[str, Any], case: Mapping[str, Any]
) -> dict[str, Any]:
    positions_raw = payload.get("cpg_positions")
    profiles_raw = payload.get("state_mean_profiles")
    effects_raw = payload.get("pairwise_effects")
    active_positions_raw = payload.get("active_positions")
    if active_positions_raw != case["active_positions"]:
        raise ReportContractError(
            f"detail/evidence active_positions mismatch for {identity_key(case)}"
        )
    if not isinstance(positions_raw, list) or not all(
        isinstance(position, int) for position in positions_raw
    ):
        raise ReportContractError("detail cpg_positions must be an integer list")
    if not isinstance(profiles_raw, Mapping) or not profiles_raw:
        raise ReportContractError("detail state_mean_profiles must be non-empty")
    if not isinstance(effects_raw, list):
        raise ReportContractError("detail pairwise_effects must be a list")
    indices = evenly_spaced_indices(len(positions_raw), MAX_PROFILE_CPG)
    state_counts_raw = payload.get("state_counts")
    if not isinstance(state_counts_raw, Mapping):
        raise ReportContractError("detail state_counts must be an object")
    state_counts: dict[str, int] = {}
    for state, raw_count in state_counts_raw.items():
        if (
            not isinstance(state, str)
            or not state
            or isinstance(raw_count, bool)
            or not isinstance(raw_count, int)
            or raw_count < 0
        ):
            raise ReportContractError("detail state_counts is malformed")
        state_counts[state] = raw_count
    state_counts = dict(sorted(state_counts.items()))
    if state_counts != case["analysis_state_counts"]:
        raise ReportContractError(
            f"detail/evidence state_counts mismatch for {identity_key(case)}"
        )
    if set(profiles_raw) != set(state_counts):
        raise ReportContractError(
            f"detail profile/state key mismatch for {identity_key(case)}"
        )
    profiles: dict[str, list[float | None]] = {}
    for state, values in profiles_raw.items():
        if not isinstance(state, str) or not isinstance(values, list):
            raise ReportContractError("detail state profile is malformed")
        if len(values) != len(positions_raw):
            raise ReportContractError(
                f"detail profile width mismatch for state {state!r}"
            )
        if (
            len(state) != int(case["n_active_bits"])
            or any(code not in {"R", "A"} for code in state)
        ):
            raise ReportContractError(
                f"detail profile has invalid complete state {state!r}"
            )
        profiles[state] = [
            None if values[index] is None else float_value(values[index])
            for index in indices
        ]
    effects: list[dict[str, Any]] = []
    seen_pairs: set[frozenset[str]] = set()
    for raw in effects_raw:
        if not isinstance(raw, Mapping):
            raise ReportContractError("detail pairwise effect must be an object")
        first = str(raw.get("first", ""))
        second = str(raw.get("second", ""))
        if not first or not second or "X" in first or "X" in second:
            raise ReportContractError(
                "detail pairwise effects must contain complete R/A states only"
            )
        if first == second or first not in state_counts or second not in state_counts:
            raise ReportContractError(
                "detail pairwise effect references an unknown or identical state"
            )
        pair_key = frozenset((first, second))
        if pair_key in seen_pairs:
            raise ReportContractError("detail contains a duplicate state pair")
        seen_pairs.add(pair_key)
        observed_hamming = int(raw["hamming"])
        expected_hamming = pattern_hamming(first, second)
        if observed_hamming != expected_hamming:
            raise ReportContractError(
                f"detail pair hamming mismatch for {first}|{second}: "
                f"{observed_hamming} != {expected_hamming}"
            )
        relation = str(raw.get("topology_relation", ""))
        if expected_hamming > 1 and relation != HAMMING_GT1_RELATION:
            raise ReportContractError(
                f"detail Hamming>1 pair has directional/invalid relation: {relation!r}"
            )
        if expected_hamming == 1 and relation not in HAMMING1_RELATIONS:
            raise ReportContractError(
                f"detail Hamming=1 pair has invalid relation: {relation!r}"
            )
        effects.append(
            {
                "first": first,
                "second": second,
                "hamming": observed_hamming,
                "between_mean": float_value(raw.get("between_mean")),
                "pooled_within_mean": float_value(raw.get("pooled_within_mean")),
                "distance_contrast": float_value(raw.get("distance_contrast")),
                "standardized_effect": float_value(raw.get("standardized_effect")),
                "topology_relation": relation,
            }
        )
    expected_pairs = {
        frozenset((first, second))
        for index, first in enumerate(state_counts)
        for second in list(state_counts)[index + 1 :]
    }
    if seen_pairs != expected_pairs:
        raise ReportContractError(
            f"detail pairwise state-pair coverage mismatch for {identity_key(case)}"
        )
    read_cluster = normalize_read_cluster(
        payload,
        case=case,
        state_counts=state_counts,
    )
    best_pair = str(case["best_pair"])
    best_parts = best_pair.split("|")
    if len(best_parts) != 2 or not all(best_parts):
        raise ReportContractError(
            f"detail-bearing evidence has malformed best_pair: {best_pair!r}"
        )
    best_key = frozenset(best_parts)
    best_effect = next(
        (effect for effect in effects if frozenset((effect["first"], effect["second"])) == best_key),
        None,
    )
    if best_effect is None:
        raise ReportContractError(
            f"evidence best_pair is absent from detail: {best_pair!r}"
        )
    if (
        best_effect["hamming"] != case["best_pair_hamming"]
        or best_effect["topology_relation"] != case["best_pair_topology_relation"]
    ):
        raise ReportContractError(
            f"detail/evidence best-pair hamming/relation mismatch for "
            f"{identity_key(case)}"
        )
    topology_raw = payload.get("topology")
    topology: dict[str, Any] | None = None
    if isinstance(topology_raw, Mapping):
        undirected_edges: list[dict[str, Any]] = []
        raw_edges = topology_raw.get("representative_best_edges", [])
        if not isinstance(raw_edges, list):
            raise ReportContractError(
                "detail topology representative_best_edges must be a list"
            )
        for edge in raw_edges:
            if not isinstance(edge, Mapping):
                raise ReportContractError("detail topology edge must be an object")
            endpoints = [
                {
                    "vertex": edge.get("parent_vertex"),
                    "label": edge.get("parent_label"),
                },
                {
                    "vertex": edge.get("child_vertex"),
                    "label": edge.get("child_label"),
                },
            ]
            endpoints.sort(
                key=lambda endpoint: (
                    str(endpoint.get("label", "")),
                    str(endpoint.get("vertex", "")),
                )
            )
            undirected_edges.append({"endpoints": endpoints})
        topology = {
            "best_tree_unique": bool(topology_raw.get("best_tree_unique")),
            "best_tree_tie_count": topology_raw.get("best_tree_tie_count"),
            "representative_edge_memberships": undirected_edges,
        }
    return {
        "cpg_positions_total": len(positions_raw),
        "cpg_positions_displayed": [positions_raw[index] for index in indices],
        "state_mean_profiles": dict(sorted(profiles.items())),
        "state_counts": state_counts,
        "pairwise_effects": effects,
        "topology": topology,
        "balanced_n_per_state": payload.get("balanced_n_per_state"),
        "read_cluster": read_cluster,
    }


def load_details(
    path: Path,
    analysis_summary: Mapping[str, Any],
    cases: Sequence[dict[str, Any]],
) -> None:
    require_file(path, "analysis details")
    verify_bound_output(analysis_summary, "details", path, "analysis summary")
    by_key = {identity_key(case): case for case in cases}
    seen: set[tuple[str, str, str, str]] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                raise ReportContractError(f"{path}:{line_number}: blank JSONL row")
            try:
                payload = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ReportContractError(
                    f"{path}:{line_number}: malformed detail JSON"
                ) from exc
            if not isinstance(payload, dict):
                raise ReportContractError(
                    f"{path}:{line_number}: detail row must be an object"
                )
            if payload.get("schema_name") != (
                "intersubmod.pattern_methyl_evidence.detail"
            ):
                raise ReportContractError(
                    f"{path}:{line_number}: unexpected detail schema_name"
                )
            if payload.get("schema_version") != "1.0.0":
                raise ReportContractError(
                    f"{path}:{line_number}: unexpected detail schema_version"
                )
            key = identity_key(payload)
            if key in seen:
                raise ReportContractError(f"duplicate detail key: {key}")
            seen.add(key)
            case = by_key.get(key)
            if case is None:
                raise ReportContractError(f"detail key absent from evidence: {key}")
            if payload.get("assessment") != case["assessment"]:
                raise ReportContractError(f"detail assessment mismatch for {key}")
            case["detail"] = normalize_detail(payload, case)
    declared_counts = analysis_summary.get("counts")
    if not isinstance(declared_counts, Mapping):
        raise ReportContractError("analysis summary has no counts mapping")
    expected = int(declared_counts.get("detail_records", -1))
    if len(seen) != expected:
        raise ReportContractError(
            f"analysis detail count mismatch: {len(seen)} != {expected}"
        )


def validate_analysis_summary(
    summary: Mapping[str, Any], cases: Sequence[Mapping[str, Any]]
) -> dict[str, int]:
    if summary.get("schema_name") not in ANALYSIS_SUMMARY_SCHEMAS:
        raise ReportContractError("unexpected analysis summary schema_name")
    if summary.get("schema_version") != "1.0.0":
        raise ReportContractError("unexpected analysis summary schema_version")
    counts = summary.get("counts")
    if not isinstance(counts, Mapping):
        raise ReportContractError("analysis summary has no counts mapping")
    if int(counts.get("units_total", -1)) != len(cases):
        raise ReportContractError(
            f"analysis unit count mismatch: {counts.get('units_total')} != {len(cases)}"
        )
    actual = Counter(str(case["assessment"]) for case in cases)
    declared = counts.get("assessment")
    if not isinstance(declared, Mapping):
        raise ReportContractError("analysis summary assessment counts are missing")
    normalized_declared = {
        str(key): int(value) for key, value in declared.items()
    }
    if dict(sorted(actual.items())) != dict(sorted(normalized_declared.items())):
        raise ReportContractError(
            "analysis summary assessment counts do not match evidence"
        )
    evaluable = sum(case["evaluation_status"] == "EVALUABLE" for case in cases)
    if int(counts.get("units_evaluable", -1)) != evaluable:
        raise ReportContractError("analysis evaluable count does not match evidence")
    return {assessment: int(actual[assessment]) for assessment in ASSESSMENTS}


def extreme_sort_key(case: Mapping[str, Any]) -> tuple[object, ...]:
    q_value = case.get("q_by")
    r_squared = case.get("permanova_r2")
    contrast = case.get("best_pair_distance_contrast")
    return (
        ASSESSMENT_PRIORITY[str(case["assessment"])],
        case.get("detail") is None,
        math.inf if q_value is None else float(q_value),
        -(-1.0 if r_squared is None else float(r_squared)),
        -(-1.0 if contrast is None else float(contrast)),
        case["dataset"],
        natural_chrom_key(str(case["chrom"])),
        case["region_id"],
        case["hp_raw"],
    )


def explained_sort_key(case: Mapping[str, Any]) -> tuple[object, ...]:
    assessment_order = {
        "CONFOUNDED": 0,
        "NOT_EVALUABLE": 1,
        "EVALUABLE_NO_ROBUST_ASSOCIATION": 2,
        "TAG_DEPENDENT": 3,
        "LOCAL_CIS_COMPATIBLE": 4,
        "ROBUST_ASSOCIATION": 5,
    }
    geometry = case.get("max_geometry_smd")
    return (
        assessment_order[str(case["assessment"])],
        -(0.0 if geometry is None else float(geometry)),
        case["dataset"],
        natural_chrom_key(str(case["chrom"])),
        case["region_id"],
        case["hp_raw"],
    )


def build_sentinels(
    matches: Sequence[Mapping[str, Any]],
    cases: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    case_by_key = {identity_key(case): case for case in cases}
    output: list[dict[str, Any]] = []
    for target in SENTINEL_TARGETS:
        target_matches = [
            dict(row)
            for row in matches
            if row["dataset"] == target["dataset"]
            and row["chrom"] == target["chrom"]
            and int(target["position"]) in row["active_positions"]
        ]
        target_matches.sort(
            key=lambda row: (
                row["hp_raw"] != target["hp_raw"],
                not bool(row["formal_n5"]),
                -int(row["n_complete"]),
                row["region_id"],
                row["hp_raw"],
            )
        )
        rendered_matches: list[dict[str, Any]] = []
        for row in target_matches[:MAX_SENTINEL_MATCHES]:
            key = identity_key(row)
            evidence = case_by_key.get(key)
            rendered_matches.append(
                {
                    **row,
                    "case_id": evidence["id"] if evidence else None,
                    "assessment": (
                        evidence["assessment"] if evidence else "NOT_IN_FORMAL_EVIDENCE"
                    ),
                    "analysis_n": evidence["analysis_n"] if evidence else None,
                    "permanova_r2": (
                        evidence["permanova_r2"] if evidence else None
                    ),
                    "q_by": evidence["q_by"] if evidence else None,
                }
            )
        output.append(
            {
                **target,
                "found_in_pattern_counts": bool(target_matches),
                "match_count": len(target_matches),
                "preferred_hp_match_count": sum(
                    row["hp_raw"] == target["hp_raw"] for row in target_matches
                ),
                "display_limit": MAX_SENTINEL_MATCHES,
                "displayed_match_count": len(rendered_matches),
                "matches_truncated": len(target_matches) > len(rendered_matches),
                "truncated_match_count": len(target_matches) - len(rendered_matches),
                "matches": rendered_matches,
            }
        )
    return output


def source_binding(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def technical_headline(assessment_counts: Mapping[str, int], full_scope: bool) -> str:
    scope_text = "完整 7-dataset、chr1–22 frozen family" if full_scope else "目前輸入範圍"
    robust = int(assessment_counts.get("ROBUST_ASSOCIATION", 0))
    local = int(assessment_counts.get("LOCAL_CIS_COMPATIBLE", 0))
    confounded = int(assessment_counts.get("CONFOUNDED", 0))
    if robust:
        return (
            f"{scope_text} 找到 {robust} 個通過正式與 robustness gates 的"
            f"區域甲基關聯；另有 {local} 個僅與局部 cis 相容、"
            f"{confounded} 個未通過既定混雜 gate（此分類不證明混雜造成結果）。"
        )
    return (
        f"{scope_text} 尚無通過全部正式與 robustness gates 的區域甲基關聯；"
        f"{local} 個僅與局部 cis 相容，{confounded} 個未通過既定混雜 gate"
        "（此分類不證明混雜造成結果）。"
    )


def strongest_secondary_case(
    cases: Sequence[Mapping[str, Any]],
) -> Mapping[str, Any] | None:
    secondary = [
        case
        for case in cases
        if case["multiplicity_family"] == "SECONDARY_PAIR_CONTRAST"
        and case["evaluation_status"] == "EVALUABLE"
        and case["permanova_r2"] is not None
    ]
    return max(
        secondary,
        key=lambda case: (
            float(case["permanova_r2"]),
            float(case["best_pair_distance_contrast"] or -math.inf),
        ),
        default=None,
    )


def direct_answer(cases: Sequence[Mapping[str, Any]]) -> str:
    full4 = sum(bool(case["pair_full4"]) for case in cases)
    robust = [
        case for case in cases if case["assessment"] == "ROBUST_ASSOCIATION"
    ]
    robust_long = sum(bool(case["k_ge_3"]) for case in robust)
    strongest = strongest_secondary_case(cases)
    first = (
        f"同一 exact raw HP 下，formal full-four RR/RA/AR/AA 單元為 {full4}；"
        f"{len(robust)} 個 robust 關聯中有 {robust_long} 個來自 k≥3 的較長 signature。"
    )
    if strongest is None:
        return first + "目前沒有可評估的 secondary 二位點案例。"
    q_value = strongest.get("q_by")
    q_text = "NA" if q_value is None else f"{float(q_value):.3f}"
    return (
        first
        + f"最強 secondary 二位點案例是 {strongest['dataset']} "
        f"{strongest['chrom']} {strongest['best_pair']}（R² "
        f"{float(strongest['permanova_r2']):.3f}、BY q {q_text}），"
        "仍屬 evaluable · no robust association。"
    )


def build_report_data(
    *,
    census_receipt_path: Path,
    pattern_counts_path: Path,
    catalog_summary_path: Path,
    catalog_tsv_path: Path,
    bernoulli_parity_summary_path: Path,
    evidence_path: Path,
    details_path: Path,
    analysis_summary_path: Path,
) -> dict[str, Any]:
    census_receipt = load_json(census_receipt_path, "census receipt")
    if census_receipt.get("schema_name") != (
        "intersubmod.exact_raw_hp_pattern_census"
    ):
        raise ReportContractError("unexpected census receipt schema_name")
    if census_receipt.get("schema_version") != "1.0.0":
        raise ReportContractError("unexpected census receipt schema_version")
    if census_receipt.get("all_pass") is not True:
        raise ReportContractError("census receipt is not all_pass")

    analysis_summary = load_json(analysis_summary_path, "analysis summary")
    formal_counts, census_aggregate, sentinel_matches = load_pattern_counts(
        pattern_counts_path, census_receipt
    )
    formal_marker_keys = {
        (str(row["dataset"]), str(row["chrom"]), int(position))
        for row in formal_counts.values()
        for position in row["active_positions"]
    }
    catalog_summary, catalog_aggregate = load_catalog(
        catalog_summary_path, catalog_tsv_path, formal_marker_keys
    )
    bernoulli_parity = load_bernoulli_parity(
        bernoulli_parity_summary_path,
        catalog_tsv_path,
        catalog_aggregate,
    )
    cases = load_evidence(evidence_path, analysis_summary, formal_counts)
    assessment_counts = validate_analysis_summary(analysis_summary, cases)
    load_details(details_path, analysis_summary, cases)

    scope_raw = census_receipt.get("scope")
    if not isinstance(scope_raw, Mapping):
        raise ReportContractError("census receipt has no scope")
    datasets = [str(value) for value in scope_raw.get("datasets", [])]
    chromosomes = [str(value) for value in scope_raw.get("chromosomes", [])]
    if not datasets or not chromosomes:
        raise ReportContractError("census scope has no datasets or chromosomes")
    full_scope = (
        len(datasets) == 7
        and set(chromosomes) == {f"chr{index}" for index in range(1, 23)}
    )

    extreme = sorted(cases, key=extreme_sort_key)[:MAX_LAYER_CASES]
    explained = sorted(cases, key=explained_sort_key)[:MAX_LAYER_CASES]
    sentinels = build_sentinels(sentinel_matches, cases)
    detail_records = sum(case["detail"] is not None for case in cases)
    read_cluster_records = sum(
        bool(case["detail"]["read_cluster"]["available"])
        for case in cases
        if case["detail"] is not None
    )
    read_cluster_source_null_records = sum(
        case["detail"]["read_cluster"]["source_status"]
        == "SOURCE_NULL_N_GT_160"
        for case in cases
        if case["detail"] is not None
    )
    read_cluster_evaluable_no_detail_records = sum(
        case["evaluation_status"] == "EVALUABLE" and case["detail"] is None
        for case in cases
    )
    read_cluster_non_evaluable_records = sum(
        case["evaluation_status"] != "EVALUABLE"
        and (
            case["detail"] is None
            or not case["detail"]["read_cluster"]["available"]
        )
        for case in cases
    )
    if (
        read_cluster_records
        + read_cluster_source_null_records
        + read_cluster_evaluable_no_detail_records
        + read_cluster_non_evaluable_records
        != len(cases)
    ):
        raise ReportContractError(
            "read-cluster availability partition does not equal formal cases"
        )
    read_cluster_reads_total = sum(
        int(case["detail"]["read_cluster"]["n_reads"])
        for case in cases
        if case["detail"] is not None
        and case["detail"]["read_cluster"]["available"]
    )
    read_cluster_cells_total = sum(
        int(case["detail"]["read_cluster"]["n_reads"]) ** 2
        for case in cases
        if case["detail"] is not None
        and case["detail"]["read_cluster"]["available"]
    )
    strongest_secondary = strongest_secondary_case(cases)
    data: dict[str, Any] = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "title": TITLE,
        "claim_ceiling": CLAIM_CEILING_ZH,
        "technical_headline": technical_headline(assessment_counts, full_scope),
        "direct_answer": direct_answer(cases),
        "scope": {
            "task_type": str(scope_raw.get("task_type", "B_comprehensive_validation")),
            "datasets": datasets,
            "chromosomes": chromosomes,
            "technical_dataset_count": len(datasets),
            "biological_sample_count": 6 if full_scope else None,
            "full_task_b_scope": full_scope,
        },
        "aggregate": {
            **census_aggregate,
            **catalog_aggregate,
            "analysis_units": len(cases),
            "analysis_evaluable": sum(
                case["evaluation_status"] == "EVALUABLE" for case in cases
            ),
            "detail_records": detail_records,
            "read_cluster_records": read_cluster_records,
            "read_cluster_source_null_records": read_cluster_source_null_records,
            "read_cluster_evaluable_no_detail_records": (
                read_cluster_evaluable_no_detail_records
            ),
            "read_cluster_non_evaluable_records": (
                read_cluster_non_evaluable_records
            ),
            "read_cluster_unavailable_records": len(cases) - read_cluster_records,
            "read_cluster_reads_total": read_cluster_reads_total,
            "read_cluster_cells_total": read_cluster_cells_total,
            "resolution_limited_units": sum(
                bool(case["resolution_limited"]) for case in cases
            ),
            "full4_units": sum(bool(case["pair_full4"]) for case in cases),
            "k_ge_3_units": sum(bool(case["k_ge_3"]) for case in cases),
            "robust_k_ge_3_units": sum(
                case["assessment"] == "ROBUST_ASSOCIATION"
                and bool(case["k_ge_3"])
                for case in cases
            ),
            "dataset_formal_unit_counts": {
                dataset: sum(case["dataset"] == dataset for case in cases)
                for dataset in datasets
            },
            "assessment_counts": assessment_counts,
        },
        "bernoulli_parity": bernoulli_parity,
        "sentinels": sentinels,
        "extreme_case_ids": [case["id"] for case in extreme],
        "explained_case_ids": [case["id"] for case in explained],
        "strongest_secondary_case_id": (
            strongest_secondary["id"] if strongest_secondary is not None else None
        ),
        "cases": cases,
        "method_contract": {
            "primary_grain": (
                "dataset × chromosome × exact PS × exact raw HP × topology region"
            ),
            "formal_gate": (
                "active k≥2、complete R/A N≥40、至少兩個 complete states 各 n≥5"
            ),
            "partial_pattern_rule": (
                "含 X 的 pattern 只保留為 subcube evidence，不映射 topology vertex"
            ),
            "pair_relation_rule": (
                "Hamming>1 只畫 pair band；不得投影為 topology edge"
            ),
            "resolution_limited_rule": (
                "confirmatory family、49,999 requested permutations、"
                "p=1/50,000 且 BY q>0.05；只加旗標，不改 assessment"
            ),
            "methyl_role": "association overlay only; topology 保持 frozen",
            "read_cluster_rule": (
                "報告側以 Bernoulli read×read distance 做 label-independent "
                "UPGMA 排序；pattern/strand/read-group 只在排序後疊加，"
                "exact-distance ties 以 source leaf ordinal 決勝；"
                "RG 僅輸出匿名類別；不切 cluster、不建立 clone 或演化關係"
            ),
            "read_cluster_coverage_rule": (
                "source detail gate 為 primary p≤0.05 或 R²≥0.05；"
                "distance matrix 僅在 analysis N≤160 時存在"
            ),
            "read_cluster_encoding": (
                "0..1 distance 量化為 uint8 0..254；255 保留 NA；"
                f"最大顯示量化誤差 {0.5 / READ_CLUSTER_DISTANCE_LEVELS:.7f}"
            ),
        },
        "chart_map": [
            {
                "section": "aggregate",
                "question": "各 assessment 類別的 formal unit 數量",
                "family": "comparison",
                "type": "horizontal bars",
                "fields": ["assessment", "count"],
            },
            {
                "section": "case detail",
                "question": "各 complete state 的 CpG 平均 profile",
                "family": "matrix",
                "type": "profile heatmap",
                "fields": ["state", "CpG position", "mean methylation"],
            },
            {
                "section": "case detail",
                "question": (
                    "同一 exact PS × raw HP 內 reads 依完整 R/A pattern "
                    "呈現局部聚集或混合"
                ),
                "family": "matrix and cohort",
                "type": "clustered read×read Bernoulli distance heatmap",
                "fields": [
                    "anonymous read ordinal",
                    "complete pattern",
                    "strand",
                    "read group",
                    "Bernoulli distance",
                ],
            },
            {
                "section": "case detail",
                "question": "state pair 的 Bernoulli distance contrast",
                "family": "matrix",
                "type": "state-pair matrix",
                "fields": ["state A", "state B", "distance contrast", "Hamming"],
            },
            {
                "section": "case detail",
                "question": "最佳 state pair 與 frozen topology 的關係",
                "family": "relationship",
                "type": "topology relation band",
                "fields": ["best pair", "Hamming", "topology relation"],
            },
        ],
        "sources": {
            "census_receipt": source_binding(census_receipt_path),
            "pattern_counts": source_binding(pattern_counts_path),
            "artifact_catalog_summary": source_binding(catalog_summary_path),
            "artifact_catalog_tsv": source_binding(catalog_tsv_path),
            "bernoulli_parity_summary": source_binding(
                bernoulli_parity_summary_path
            ),
            "evidence": source_binding(evidence_path),
            "details": source_binding(details_path),
            "analysis_summary": source_binding(analysis_summary_path),
        },
        "catalog_contract": {
            "status": catalog_summary.get("status"),
            "reads_hp_authoritative": False,
            "exact_raw_hp_authority": (
                "frozen LongPhase-S sparse read-call sidecar"
            ),
        },
    }
    # Reject non-finite values before the same object is written and embedded.
    json.dumps(data, ensure_ascii=False, allow_nan=False)
    return data


HTML_TEMPLATE = r"""<!doctype html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<title>多 sSNV pattern × 區域甲基關聯驗證</title>
<style>
:root{
  --paper:#f4f0e7;--paper-2:#fffdf8;--ink:#142327;--muted:#58676a;
  --line:#c9c4b8;--navy:#153b52;--blue:#1e6585;--cyan:#d7edf2;
  --copper:#9f4d2f;--gold:#c48a2b;--olive:#63734e;--pink:#a94f68;
  --shadow:0 16px 40px rgba(20,35,39,.10);--radius:18px;
  --display:"Noto Serif TC","Source Han Serif TC","Songti TC","PMingLiU",serif;
  --body:"Noto Sans TC","PingFang TC","Microsoft JhengHei",sans-serif;
  --mono:"IBM Plex Mono","SFMono-Regular",Consolas,monospace;
}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{
  margin:0;color:var(--ink);font-family:var(--body);line-height:1.65;
  background-color:var(--paper);
  background-image:radial-gradient(rgba(21,59,82,.10) .65px,transparent .65px);
  background-size:13px 13px;
}
button,select{font:inherit}
button,select,a{min-height:44px}
:focus-visible{outline:3px solid var(--gold);outline-offset:3px}
.skip-link{position:absolute;left:12px;top:-80px;background:var(--ink);color:#fff;padding:10px 14px;z-index:20}
.skip-link:focus{top:12px}
.page{width:min(1440px,100%);margin:0 auto;background:rgba(244,240,231,.94)}
.hero{position:relative;padding:clamp(28px,5vw,78px) clamp(18px,5vw,76px) 34px;border-bottom:1px solid var(--line);overflow:hidden}
.hero::after{content:"";position:absolute;right:-70px;top:-90px;width:260px;height:260px;border:42px solid rgba(30,101,133,.10);border-radius:50%;pointer-events:none}
.eyebrow{font:700 .76rem/1.2 var(--mono);letter-spacing:.16em;text-transform:uppercase;color:var(--copper)}
h1,h2,h3{font-family:var(--display);letter-spacing:.01em}
h1{max-width:940px;margin:.55rem 0 1rem;font-size:clamp(2.2rem,5.4vw,5.5rem);line-height:1.04}
.headline{max-width:930px;margin:0;font-size:clamp(1.04rem,2vw,1.35rem);font-weight:650}
.ceiling{display:grid;grid-template-columns:auto 1fr;gap:12px;align-items:start;max-width:1040px;margin:24px 0 0;padding:16px 18px;background:var(--navy);color:#fff;border-left:7px solid var(--gold);border-radius:0 12px 12px 0}
.ceiling b{font:700 .72rem/1.4 var(--mono);letter-spacing:.08em;color:#f5ca78}
.answer-band{display:grid;grid-template-columns:auto 1fr;gap:12px;align-items:start;max-width:1040px;margin-top:12px;padding:14px 18px;border:1px solid rgba(159,77,47,.45);border-left:7px solid var(--copper);border-radius:0 12px 12px 0;background:#fffaf0}
.answer-band b{font:700 .72rem/1.4 var(--mono);letter-spacing:.08em;color:var(--copper)}
.scope-chips{display:flex;flex-wrap:wrap;gap:8px;margin-top:18px}
.chip,.badge{display:inline-flex;align-items:center;gap:6px;border:1px solid var(--line);border-radius:999px;padding:5px 10px;background:var(--paper-2);font:650 .78rem/1.3 var(--mono)}
.resolution-flag{display:inline-flex;align-items:center;width:max-content;margin-top:7px;padding:4px 8px;border:1px dashed var(--copper);border-radius:6px;background:#f7ecd2;color:var(--copper);font:700 .68rem/1.3 var(--mono)}
.parity-strip{display:grid;grid-template-columns:auto repeat(3,minmax(0,1fr));gap:12px;align-items:center;max-width:1040px;margin-top:18px;padding:13px 15px;border:1px solid rgba(21,59,82,.35);border-radius:12px;background:rgba(255,253,248,.88)}
.parity-strip>strong{font:700 .72rem/1.3 var(--mono);letter-spacing:.08em;color:var(--blue)}.parity-stat{display:grid;gap:2px}.parity-stat span{font-size:.7rem;color:var(--muted)}.parity-stat b{font:700 .88rem/1.2 var(--mono)}
.badge{white-space:nowrap}
.badge[data-assessment="ROBUST_ASSOCIATION"]{border-color:var(--blue);background:var(--cyan)}
.badge[data-assessment="LOCAL_CIS_COMPATIBLE"]{border-color:var(--olive);background:#e8eddc}
.badge[data-assessment="TAG_DEPENDENT"]{border-color:var(--gold);background:#f8e8c3}
.badge[data-assessment="CONFOUNDED"]{border-color:var(--copper);background:#f4ddd2}
.badge[data-assessment="NOT_EVALUABLE"]{border-style:dashed;background:#ece9e1}
.filter-shell{position:sticky;top:0;z-index:10;padding:12px clamp(14px,4vw,56px);background:rgba(255,253,248,.96);border-bottom:1px solid var(--line);box-shadow:0 7px 18px rgba(20,35,39,.07);backdrop-filter:blur(10px)}
.filters{display:grid;grid-template-columns:repeat(5,minmax(120px,1fr)) auto;gap:10px;align-items:end}
.field{display:grid;gap:4px}.field label{font:700 .69rem/1.2 var(--mono);letter-spacing:.07em;color:var(--muted)}
select{width:100%;border:1px solid #9da8a7;border-radius:9px;background:#fff;color:var(--ink);padding:8px 34px 8px 10px}
.reset{border:1px solid var(--navy);border-radius:9px;background:var(--navy);color:#fff;padding:8px 16px;cursor:pointer}
.filter-status{grid-column:1/-1;margin:0;font-size:.82rem;color:var(--muted)}
main{padding:0 clamp(14px,4vw,56px) 70px}
.layer{padding:clamp(44px,7vw,92px) 0;border-bottom:1px solid var(--line)}
.section-kicker{margin:0;color:var(--copper);font:700 .72rem/1.2 var(--mono);letter-spacing:.14em}
h2{margin:.45rem 0 .8rem;font-size:clamp(1.65rem,3.5vw,3.15rem);line-height:1.12}
.section-lede{max-width:900px;margin:0 0 24px;color:var(--muted)}
.metric-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin:22px 0 30px}
.metric{padding:18px;background:var(--paper-2);border:1px solid var(--line);border-radius:var(--radius);box-shadow:var(--shadow)}
.metric span{display:block;color:var(--muted);font-size:.78rem}.metric strong{display:block;margin-top:5px;font:700 clamp(1.65rem,3vw,2.6rem)/1 var(--mono)}
.chart-card,.case-detail,.sentinel,.case-card{background:var(--paper-2);border:1px solid var(--line);border-radius:var(--radius);box-shadow:var(--shadow)}
.chart-card{padding:22px;margin-top:20px}.chart-card h3{margin:0;font-size:1.22rem}.chart-note{margin:4px 0 18px;color:var(--muted);font-size:.88rem}
.assessment-bars{display:grid;gap:10px}.bar-row{display:grid;grid-template-columns:minmax(190px,1fr) 4fr 58px;gap:10px;align-items:center}
.bar-label{font:650 .77rem/1.3 var(--mono)}.bar-track{height:16px;border:1px solid #b8b5ad;background:#ece9e1}.bar-fill{display:block;height:100%;background:var(--blue)}.bar-value{text-align:right;font:700 .8rem/1 var(--mono)}
.sentinel-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}
.sentinel{padding:22px}.sentinel h3{margin:0}.sentinel .coord{font:700 .8rem/1.4 var(--mono);color:var(--blue)}
.empty-note{padding:16px;border:1px dashed #9da8a7;background:#f2efe8;color:var(--muted)}
.match{margin-top:14px;padding-top:14px;border-top:1px solid var(--line)}.match dl{display:grid;grid-template-columns:auto 1fr;gap:3px 12px;margin:9px 0 0}.match dt{color:var(--muted)}.match dd{margin:0;font-family:var(--mono)}
.cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}.case-card{display:flex;flex-direction:column;padding:18px;min-width:0}
.case-card h3{margin:.5rem 0 .3rem;font:700 .92rem/1.35 var(--mono);overflow-wrap:anywhere}.case-card p{margin:0 0 12px;color:var(--muted);font-size:.84rem}
.case-card .mini{display:grid;grid-template-columns:1fr 1fr;gap:8px;margin:auto 0 13px}.mini span{padding:7px;background:#f0ede5;font:650 .74rem/1.2 var(--mono)}
.secondary-highlight{margin-top:24px;padding:18px;border:1px solid rgba(159,77,47,.45);border-radius:var(--radius);background:#fffaf0}.secondary-highlight>h3{margin:0}.secondary-highlight>p{max-width:900px;color:var(--muted)}.secondary-highlight .cards{grid-template-columns:minmax(0,420px)}
.open-case{border:1px solid var(--navy);border-radius:8px;background:transparent;color:var(--navy);font-weight:700;cursor:pointer}
.explorer{padding-top:58px}.explorer h2{font-size:clamp(1.55rem,3vw,2.5rem)}
.case-table-wrap{overflow:auto;border:1px solid var(--line);border-radius:14px;background:var(--paper-2)}
table{border-collapse:collapse;width:100%}.case-table th,.case-table td{padding:10px 12px;border-bottom:1px solid #dedad1;text-align:left;vertical-align:top}.case-table th{position:sticky;top:0;background:#e8e3d8;font:700 .72rem/1.3 var(--mono);z-index:1}.case-table td{font-size:.82rem}.case-table code{font-family:var(--mono);overflow-wrap:anywhere}
.table-action{min-height:36px;border:1px solid var(--blue);border-radius:7px;background:#fff;color:var(--blue);cursor:pointer}
.case-detail{max-width:100%;min-width:0;margin-top:22px;padding:clamp(16px,3vw,32px);overflow:hidden;scroll-margin-top:132px}.detail-head{display:flex;justify-content:space-between;gap:20px;align-items:flex-start}.detail-head h3{margin:.25rem 0;font:700 clamp(1rem,2vw,1.35rem)/1.35 var(--mono);overflow-wrap:anywhere}
.detail-metrics{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:8px;margin:18px 0}.detail-metrics div{padding:9px;background:#f0ede5}.detail-metrics span{display:block;color:var(--muted);font-size:.7rem}.detail-metrics b{font:700 .88rem/1.2 var(--mono)}
.visual-grid{display:grid;min-width:0;gap:18px}.visual{min-width:0;padding-top:22px;border-top:1px solid var(--line)}.visual h3{margin:0;font-size:1.2rem}.visual p{margin:5px 0 13px;color:var(--muted);font-size:.88rem}
.matrix-scroll{width:100%;max-width:100%;min-width:0;overflow:auto;padding-bottom:5px}.profile-svg,.relation-svg{display:block;min-width:720px;width:100%;height:auto;background:#fff}
.read-cluster-shell{display:grid;gap:12px}.cluster-summary-grid{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:8px}.cluster-summary-grid div{padding:9px;border:1px solid #dedad1;background:#f8f5ee}.cluster-summary-grid span{display:block;color:var(--muted);font-size:.68rem}.cluster-summary-grid b{font:700 .82rem/1.25 var(--mono)}
.read-cluster-scroll{position:relative;width:100%;max-width:100%;overflow:auto;border:1px solid #dedad1;background:#fff}.read-cluster-stage{position:relative;width:max-content;min-width:760px;padding:12px;background:#fff}.read-cluster-canvas{display:block;background:#fff}.read-cluster-canvas:focus-visible{outline:4px solid var(--gold);outline-offset:-4px}
.cluster-tooltip{position:absolute;z-index:4;display:none;max-width:330px;padding:8px 10px;border:1px solid var(--ink);border-radius:7px;background:rgba(20,35,39,.94);color:#fff;font:650 .72rem/1.45 var(--mono);pointer-events:none;box-shadow:0 8px 20px rgba(20,35,39,.24)}
.cluster-legend{display:flex;flex-wrap:wrap;align-items:center;gap:8px 14px;padding:10px 12px;border:1px solid #dedad1;background:#f8f5ee}.cluster-legend strong{font:700 .7rem/1.2 var(--mono);letter-spacing:.06em}.distance-key{display:inline-flex;align-items:center;gap:7px;font-size:.75rem}.distance-gradient{width:132px;height:12px;border:1px solid #8d9695;background:linear-gradient(90deg,#fffdf8,#d7edf2,#1e6585,#153b52)}
.annotation-key{display:inline-flex;align-items:center;gap:6px;font:650 .72rem/1.2 var(--mono)}.annotation-swatch{width:12px;height:12px;border:1px solid rgba(20,35,39,.35)}.cluster-live{min-height:1.4em;margin:0;color:var(--muted);font:650 .75rem/1.4 var(--mono)}
.state-matrix{width:auto;min-width:520px}.state-matrix th,.state-matrix td{min-width:76px;height:62px;padding:6px;border:2px solid #fff;text-align:center;font:650 .72rem/1.25 var(--mono)}.state-matrix th{background:#e8e3d8}.state-matrix td[data-hamming="gt1"]{outline:2px dashed var(--copper);outline-offset:-4px}.state-matrix small{display:block;color:var(--muted)}
.rule-band{display:grid;grid-template-columns:auto 1fr;gap:10px;margin:18px 0 0;padding:13px;border-left:5px solid var(--gold);background:#f7ecd2}.rule-band b{font-family:var(--mono)}
.sr-only{position:absolute!important;width:1px!important;max-width:1px!important;height:1px!important;padding:0!important;margin:-1px!important;overflow:hidden!important;clip:rect(0,0,0,0)!important;white-space:nowrap!important;border:0!important;contain:strict!important;table-layout:fixed!important}
.footer{padding:30px clamp(18px,5vw,76px);background:var(--ink);color:#e9f0ef}.footer p{max-width:960px;margin:0}.footer code{font-family:var(--mono);color:#f5ca78}
@media (max-width:1050px){
  .filters{grid-template-columns:repeat(3,minmax(120px,1fr))}.metric-grid,.cards{grid-template-columns:repeat(2,minmax(0,1fr))}.detail-metrics{grid-template-columns:repeat(3,minmax(0,1fr))}.cluster-summary-grid{grid-template-columns:repeat(3,minmax(0,1fr))}
}
@media (max-width:700px){
  body{background-image:none}.hero{padding:28px 18px}.filter-shell{position:relative;padding:12px 14px}.filters{grid-template-columns:1fr 1fr}.reset{grid-column:1/-1}
  .parity-strip{grid-template-columns:1fr 1fr}.parity-strip>strong{grid-column:1/-1}
  main{padding-inline:14px}.metric-grid,.cards,.sentinel-grid{grid-template-columns:1fr}.bar-row{grid-template-columns:1fr 2fr 44px}.detail-head{display:block}.detail-metrics{grid-template-columns:1fr 1fr}
  .case-detail{scroll-margin-top:16px}.cluster-summary-grid{grid-template-columns:1fr 1fr}
  .case-table thead{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0,0,0,0)}
  .case-table,.case-table tbody,.case-table tr,.case-table td{display:block;width:100%}.case-table tr{padding:10px;border-bottom:2px solid var(--line)}.case-table td{display:grid;grid-template-columns:105px 1fr;gap:8px;border:0;padding:5px}.case-table td::before{content:attr(data-label);font:700 .68rem/1.3 var(--mono);color:var(--muted)}
}
@media (prefers-reduced-motion:reduce){html{scroll-behavior:auto}*{animation:none!important;transition:none!important}}
@media print{
  @page{size:A4;margin:14mm}body{background:#fff;font-size:9.5pt}.page{width:100%;background:#fff}.filter-shell,.open-case,.table-action,.skip-link,.cluster-tooltip{display:none!important}.hero{padding:0 0 10mm}.layer{padding:10mm 0}.metric,.chart-card,.case-detail,.sentinel,.case-card{box-shadow:none;break-inside:avoid}.cards{grid-template-columns:1fr 1fr}.profile-svg,.relation-svg{min-width:0}.read-cluster-scroll{overflow:visible}.read-cluster-stage{min-width:0;width:100%!important;padding:0}.read-cluster-canvas{max-width:100%;height:auto!important}.case-table th{position:static}.footer{background:#fff;color:#000;border-top:1px solid #000;padding:8mm 0}.footer code{color:#000}
}
</style>
</head>
<body>
<a class="skip-link" href="#main">跳到主要內容</a>
<div class="page" id="report-root" data-testid="report-root">
  <header class="hero">
    <p class="eyebrow">Task B · comprehensive validation · exact raw HP</p>
    <h1>多 sSNV pattern × 區域甲基關聯驗證</h1>
    <p class="headline" id="technicalHeadline" data-testid="technical-summary"></p>
    <div class="ceiling" data-testid="claim-ceiling"><b>CLAIM CEILING</b><span id="claimCeiling"></span></div>
    <div class="answer-band" data-testid="direct-answer"><b>DIRECT ANSWER</b><span id="directAnswer"></span></div>
    <div class="scope-chips" id="scopeChips" data-testid="task-scope"></div>
    <div class="parity-strip" id="parityStrip" data-testid="bernoulli-parity" aria-label="Bernoulli distance parity QA"></div>
  </header>

  <aside class="filter-shell" data-testid="filter-bar" aria-label="案例篩選">
    <div class="filters">
      <div class="field"><label for="filterDataset">DATASET</label><select id="filterDataset" data-testid="filter-dataset"></select></div>
      <div class="field"><label for="filterChrom">CHROMOSOME</label><select id="filterChrom" data-testid="filter-chrom"></select></div>
      <div class="field"><label for="filterHp">EXACT RAW HP</label><select id="filterHp" data-testid="filter-hp"></select></div>
      <div class="field"><label for="filterPattern">COMPLETE PATTERN</label><select id="filterPattern" data-testid="filter-pattern"></select></div>
      <div class="field"><label for="filterAssessment">ASSESSMENT</label><select id="filterAssessment" data-testid="filter-assessment"></select></div>
      <button class="reset" id="resetFilters" type="button">重設篩選</button>
      <p class="filter-status" id="filterStatus" aria-live="polite"></p>
    </div>
  </aside>

  <main id="main">
    <section class="layer" id="aggregateLayer" data-testid="aggregate-layer">
      <p class="section-kicker">01 / AGGREGATE</p>
      <h2>先看完整母體，再解讀任何單一位點</h2>
      <p class="section-lede">正式單元只來自 complete R/A gate；含 X 的 reads 保留在 partial 分母，但不進正式 state 檢定。標示「目前篩選」的數量會隨篩選更新；標示「全域」者固定為完整 1,045-unit evidence universe。</p>
      <div class="metric-grid" id="metricGrid"></div>
      <div class="chart-card">
        <h3>Assessment 數量</h3>
        <p class="chart-note">同一 frozen family 的互斥 reader-facing assessment；長條從零開始，數值為 formal units。</p>
        <div class="assessment-bars" id="assessmentBars" role="img" aria-label="Assessment unit counts"></div>
      </div>
    </section>

    <section class="layer" id="canonicalLayer" data-testid="canonical-layer">
      <p class="section-kicker">02 / CANONICAL SENTINELS</p>
      <h2>兩個固定哨點不因結果排名而消失</h2>
      <p class="section-lede">H2009 chr5:18,096,980 與 HCC1395 chr22:46,257,699 是預先指定的 canonical sentinels。若未進 formal evidence，畫面會保留「未進正式分析」而不是捏造陰性。</p>
      <div class="sentinel-grid" id="sentinelGrid"></div>
    </section>

    <section class="layer" id="extremeLayer" data-testid="extreme-layer">
      <p class="section-kicker">03 / EXTREME</p>
      <h2>最強訊號仍須逐項通過混雜與 robustness gates</h2>
      <p class="section-lede">排序先看 assessment，再看 BY q、R² 與 state-pair distance contrast；它是審閱順序，不是重新定義顯著門檻。</p>
      <div class="cards" id="extremeCards"></div>
      <aside class="secondary-highlight" data-testid="strongest-secondary">
        <h3>最強 secondary exploratory</h3>
        <p>固定顯示 secondary family 中 R² 最高的案例；即使 effect 很大，只要 BY q 未過門檻，就不能提升為 robust association。</p>
        <div class="cards" id="secondaryCard"></div>
      </aside>
    </section>

    <section class="layer" id="explainedLayer" data-testid="explained-layer">
      <p class="section-kicker">04 / WELL-EXPLAINED</p>
      <h2>負結果與混雜案例同樣是主要證據</h2>
      <p class="section-lede">這一層優先顯示 CONFOUNDED、NOT_EVALUABLE 與無 robust association 的單元，讓 strand/read-group、geometry、dispersion 或資料不足不被誤寫成生物差異。</p>
      <div class="cards" id="explainedCards"></div>
    </section>

    <section class="explorer" aria-labelledby="explorerTitle">
      <p class="section-kicker">AUDITABLE CASE EXPLORER</p>
      <h2 id="explorerTitle">逐單元檢視 read 聚集、profile 與 frozen topology 關係</h2>
      <p class="section-lede">所有表格數字、heatmap cell 與 pair relation 均由輸入 evidence/details 注入。點選案例後，四個視覺使用同一 exact PS × exact raw HP 單元；read dendrogram 只用 Bernoulli distance 排序，pattern 標籤不參與建樹。</p>
      <div class="case-table-wrap" data-testid="case-list"><table class="case-table">
        <thead><tr><th>Dataset / chr</th><th>Region</th><th>HP</th><th>Patterns</th><th>N</th><th>R² / q</th><th>Assessment</th><th>操作</th></tr></thead>
        <tbody id="caseRows"></tbody>
      </table></div>
      <article class="case-detail" id="caseDetail" data-testid="case-detail" tabindex="-1"></article>
    </section>

    <div class="rule-band"><b>固定視覺規則</b><span><strong>X pattern 僅為 subcube evidence</strong>，不映射 topology vertex；<strong>Hamming&gt;1 僅畫 pair band</strong>，絕不畫成 edge。</span></div>
  </main>

  <footer class="footer"><p><code>association-only</code> · exact raw HP 是分析 strata；artifact 內舊 HP 欄位不具權威性。HCC1395 與 HCC1395_DORADO 是 technical datasets，跨樣本分母不得重複計為兩個 biological samples。</p></footer>
</div>
<script id="report-data" type="application/json">__REPORT_DATA_JSON__</script>
<script>
"use strict";
const DATA=JSON.parse(document.getElementById("report-data").textContent);
const MAX_PROFILE_CPG=72;
const $=id=>document.getElementById(id);
const esc=value=>String(value??"").replace(/[&<>"']/g,ch=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[ch]));
const fmt=(value,digits=3)=>value===null||value===undefined||value===""?"NA":Number(value).toLocaleString("zh-TW",{maximumFractionDigits:digits});
const pct=value=>value===null||value===undefined?"NA":`${(Number(value)*100).toFixed(1)}%`;
const assessmentZh={
  ROBUST_ASSOCIATION:"Robust association",LOCAL_CIS_COMPATIBLE:"Local cis compatible",
  TAG_DEPENDENT:"Raw-HP tag dependent",CONFOUNDED:"Confound gate not passed",
  EVALUABLE_NO_ROBUST_ASSOCIATION:"Evaluable · no robust association",NOT_EVALUABLE:"Not evaluable",
  NOT_IN_FORMAL_EVIDENCE:"未進 formal evidence"
};
const relationZh={
  PAIR_BAND_ONLY_HAMMING_GT1:"Hamming>1：僅 pair band，非 edge",
  HAMMING1_GLOBAL_BEST_UNANIMOUS:"Hamming=1：global-best edge membership 一致（無方向）",
  HAMMING1_NOT_UNANIMOUS:"Hamming=1：best tree 不唯一，edge membership 不一致（無方向）",
  HAMMING1_NOT_IN_GLOBAL_BEST:"Hamming=1：不在 global-best edge set（無方向）",
  HAMMING1_TOPOLOGY_UNAVAILABLE:"Hamming=1：topology edge membership 不可用（無方向）"
};
const allCases=DATA.cases;
const caseById=new Map(allCases.map(row=>[row.id,row]));
const filterIds=["filterDataset","filterChrom","filterHp","filterPattern","filterAssessment"];
function chromSort(a,b){const x=Number(a.replace("chr","")),y=Number(b.replace("chr",""));return (x-y)||a.localeCompare(b)}
function setOptions(id,values,label){
  const select=$(id);select.innerHTML=`<option value="ALL">${esc(label)}</option>`+values.map(value=>`<option value="${esc(value)}">${esc(id==="filterAssessment"?(assessmentZh[value]||value):value)}</option>`).join("");
}
function init(){
  $("technicalHeadline").textContent=DATA.technical_headline;
  $("claimCeiling").textContent=DATA.claim_ceiling;
  $("directAnswer").textContent=DATA.direct_answer;
  const s=DATA.scope;
  const datasetsWithFormal=Object.values(DATA.aggregate.dataset_formal_unit_counts).filter(value=>value>0).length;
  $("scopeChips").innerHTML=[
    `${s.task_type}`,`${s.technical_dataset_count} technical datasets`,`${s.chromosomes.length} chromosomes`,
    s.biological_sample_count?`${s.biological_sample_count} biological samples`:"subset / synthetic scope",
    `${datasetsWithFormal}/${s.technical_dataset_count} datasets with formal units`,
    "exact PS × exact raw HP",
    `${fmt(DATA.aggregate.read_cluster_records,0)} read cluster matrices`
  ].map(value=>`<span class="chip">${esc(value)}</span>`).join("");
  const parity=DATA.bernoulli_parity;
  $("parityStrip").innerHTML=[
    `<strong>BERNOULLI PARITY · ${parity.all_pass?"PASS":"FAIL"}</strong>`,
    `<span class="parity-stat"><span>Pair cells checked</span><b>${fmt(parity.pair_cells_checked,0)}</b></span>`,
    `<span class="parity-stat"><span>Invalid-mask mismatch</span><b>${fmt(parity.invalid_mask_mismatches,0)}</b></span>`,
    `<span class="parity-stat"><span>Max absolute error</span><b>${Number(parity.max_absolute_error).toExponential(2)}</b></span>`
  ].join("");
  setOptions("filterDataset",[...DATA.scope.datasets].sort(),"全部 dataset");
  setOptions("filterChrom",[...new Set(allCases.map(row=>row.chrom))].sort(chromSort),"全部 chromosome");
  setOptions("filterHp",[...new Set(allCases.map(row=>row.hp_raw))].sort(),"全部 exact HP");
  setOptions("filterPattern",[...new Set(allCases.flatMap(row=>row.patterns))].sort(),"全部 complete pattern");
  setOptions("filterAssessment",[...new Set(allCases.map(row=>row.assessment))].sort(),"全部 assessment");
  filterIds.forEach(id=>$(id).addEventListener("change",render));
  $("resetFilters").addEventListener("click",()=>{filterIds.forEach(id=>$(id).value="ALL");render()});
  document.addEventListener("click",event=>{
    const button=event.target.closest("[data-case-id]");
    if(button){selectCase(button.dataset.caseId,true)}
  });
  renderSentinels();renderSecondaryHighlight();render();
  const first=DATA.extreme_case_ids.find(id=>caseById.get(id)?.detail)||allCases.find(row=>row.detail)?.id||allCases[0]?.id;
  if(first)selectCase(first,false);
}
function filteredCases(){
  const f={dataset:$("filterDataset").value,chrom:$("filterChrom").value,hp:$("filterHp").value,pattern:$("filterPattern").value,assessment:$("filterAssessment").value};
  return allCases.filter(row=>(f.dataset==="ALL"||row.dataset===f.dataset)&&(f.chrom==="ALL"||row.chrom===f.chrom)&&(f.hp==="ALL"||row.hp_raw===f.hp)&&(f.pattern==="ALL"||row.patterns.includes(f.pattern))&&(f.assessment==="ALL"||row.assessment===f.assessment));
}
function metric(label,value,note=""){return `<div class="metric"><span>${esc(label)}</span><strong>${esc(value)}</strong>${note?`<span>${esc(note)}</span>`:""}</div>`}
function resolutionFlag(row){return row.resolution_limited?`<span class="resolution-flag" data-testid="resolution-limited-flag">PERMUTATION-RESOLUTION LIMITED</span>`:""}
function renderAggregate(rows){
  const counts={};rows.forEach(row=>counts[row.assessment]=(counts[row.assessment]||0)+1);
  $("metricGrid").innerHTML=[
    metric("目前篩選 formal units",fmt(rows.length,0),`完整 evidence：${fmt(DATA.aggregate.analysis_units,0)}`),
    metric("可評估 units",fmt(rows.filter(row=>row.evaluation_status==="EVALUABLE").length,0)),
    metric("全域 Formal full-four",fmt(DATA.aggregate.full4_units,0),"RR / RA / AR / AA"),
    metric("全域 k≥3 formal units",fmt(DATA.aggregate.k_ge_3_units,0),`robust：${fmt(DATA.aggregate.robust_k_ge_3_units,0)}`),
    metric("全域 Artifact markers",fmt(DATA.aggregate.markers_pass,0),"catalog PASS"),
    metric("全域 Partial read projections",fmt(DATA.aggregate.partial_read_projections,0),"X 不進正式 state test"),
    metric("全域 Read 聚集圖",fmt(DATA.aggregate.read_cluster_records,0),`${fmt(DATA.aggregate.read_cluster_non_evaluable_records,0)} 不可評估 · ${fmt(DATA.aggregate.read_cluster_evaluable_no_detail_records,0)} 未達 detail gate · ${fmt(DATA.aggregate.read_cluster_source_null_records,0)} N>160`),
    metric("Resolution-limited",fmt(rows.filter(row=>row.resolution_limited).length,0),"旗標；不改 assessment")
  ].join("");
  const max=Math.max(1,...Object.values(counts));
  $("assessmentBars").innerHTML=Object.entries(assessmentZh).filter(([key])=>key!=="NOT_IN_FORMAL_EVIDENCE"&&(counts[key]||DATA.aggregate.assessment_counts[key])).map(([key,label])=>{
    const count=counts[key]||0;return `<div class="bar-row"><span class="bar-label">${esc(label)}</span><span class="bar-track"><span class="bar-fill" style="width:${100*count/max}%"></span></span><span class="bar-value">${count}</span></div>`;
  }).join("")||`<p class="empty-note">目前篩選沒有案例。</p>`;
}
function renderSentinels(){
  $("sentinelGrid").innerHTML=DATA.sentinels.map(target=>{
    const matches=target.matches.length?target.matches.map(row=>`<div class="match"><span class="badge" data-assessment="${esc(row.assessment)}">${esc(assessmentZh[row.assessment]||row.assessment)}</span><dl><dt>Exact HP</dt><dd>${esc(row.hp_raw)}</dd><dt>Pattern states</dt><dd>${esc(Object.entries(row.state_counts).map(([k,v])=>`${k}:${v}`).join(" · "))}</dd><dt>Complete / partial</dt><dd>${fmt(row.n_complete,0)} / ${fmt(row.n_partial,0)}</dd><dt>R² / BY q</dt><dd>${fmt(row.permanova_r2)} / ${fmt(row.q_by)}</dd></dl>${row.case_id?`<button class="open-case" type="button" data-case-id="${esc(row.case_id)}">開啟正式案例</button>`:""}</div>`).join(""):`<p class="empty-note">此座標未出現在輸入 pattern-count active marker universe；不可據此宣稱陰性。</p>`;
    const truncation=target.matches_truncated?`<p class="empty-note" data-testid="sentinel-truncation">固定座標的 canonical exact HP ${esc(target.hp_raw)} 優先；僅顯示 ${fmt(target.displayed_match_count,0)} / ${fmt(target.match_count,0)} 個 strata，另有 ${fmt(target.truncated_match_count,0)} 個未展開。</p>`:"";
    return `<article class="sentinel"><span class="coord">${esc(target.label)} · canonical exact HP ${esc(target.hp_raw)}</span><h3>${target.found_in_pattern_counts?`${fmt(target.match_count,0)} 個 exact-HP strata 命中`:"未命中 active marker universe"}</h3>${truncation}${matches}</article>`;
  }).join("");
}
function caseCard(row){
  return `<article class="case-card"><span class="badge" data-assessment="${esc(row.assessment)}">${esc(assessmentZh[row.assessment]||row.assessment)}</span>${resolutionFlag(row)}<h3>${esc(row.dataset)} · ${esc(row.chrom)} · ${esc(row.region_id)}</h3><p>HP ${esc(row.hp_raw)} · ${esc(row.patterns.join(" / ")||"無 analysis states")}</p><div class="mini"><span>N ${fmt(row.analysis_n,0)}</span><span>R² ${fmt(row.permanova_r2)}</span><span>q ${fmt(row.q_by)}</span><span>Δdist ${fmt(row.best_pair_distance_contrast)}</span></div><button class="open-case" type="button" data-case-id="${esc(row.id)}">查看完整證據</button></article>`;
}
function renderSecondaryHighlight(){
  const row=caseById.get(DATA.strongest_secondary_case_id);
  $("secondaryCard").innerHTML=row?caseCard(row):`<p class="empty-note">目前沒有可評估的 secondary 案例。</p>`;
}
function renderLayer(containerId,ids,rows){
  const allowed=new Set(rows.map(row=>row.id));const selected=ids.map(id=>caseById.get(id)).filter(row=>row&&allowed.has(row.id));
  $(containerId).innerHTML=selected.map(caseCard).join("")||`<p class="empty-note" data-testid="empty-state">目前篩選沒有此層案例。</p>`;
}
function renderTable(rows){
  $("caseRows").innerHTML=rows.map(row=>`<tr><td data-label="Dataset / chr">${esc(row.dataset)} · ${esc(row.chrom)}</td><td data-label="Region"><code>${esc(row.region_id)}</code></td><td data-label="HP">${esc(row.hp_raw)}</td><td data-label="Patterns">${esc(row.patterns.join(" / ")||"NA")}</td><td data-label="N">${fmt(row.analysis_n,0)}</td><td data-label="R² / q">${fmt(row.permanova_r2)} / ${fmt(row.q_by)}</td><td data-label="Assessment"><span class="badge" data-assessment="${esc(row.assessment)}">${esc(assessmentZh[row.assessment]||row.assessment)}</span>${resolutionFlag(row)}</td><td data-label="操作"><button class="table-action" type="button" data-case-id="${esc(row.id)}">開啟</button></td></tr>`).join("")||`<tr><td colspan="8">目前篩選沒有案例。</td></tr>`;
}
function render(){
  const rows=filteredCases();$("filterStatus").textContent=`顯示 ${rows.length.toLocaleString("zh-TW")} / ${allCases.length.toLocaleString("zh-TW")} 個 formal units`;
  renderAggregate(rows);renderLayer("extremeCards",DATA.extreme_case_ids,rows);renderLayer("explainedCards",DATA.explained_case_ids,rows);renderTable(rows);
}
function colorFor(value){
  if(value===null||value===undefined)return "#dedbd3";const t=Math.max(0,Math.min(1,Number(value)));
  const a=[239,233,216],b=[30,101,133];return `rgb(${a.map((v,i)=>Math.round(v+(b[i]-v)*t)).join(",")})`;
}
const patternPalette=["#1e6585","#9f4d2f","#63734e","#c48a2b","#a94f68","#326f73","#6f5b91","#7a6c58","#357a4f","#9b6280","#4f6e93","#9a7038"];
const strandPalette=["#153b52","#c48a2b","#8d9695","#63734e"];
const readGroupPalette=["#58676a","#9f4d2f","#1e6585","#63734e","#a94f68","#c48a2b"];
function decodeU8(encoded){
  const binary=atob(encoded||""),values=new Uint8Array(binary.length);
  for(let index=0;index<binary.length;index+=1)values[index]=binary.charCodeAt(index);
  return values;
}
function categoryColor(kind,index){
  const palette=kind==="pattern"?patternPalette:(kind==="strand"?strandPalette:readGroupPalette);
  return palette[index%palette.length];
}
function distanceColor(code,levels,naCode){
  if(code===naCode)return "#c8c8c3";
  const t=Math.max(0,Math.min(1,code/levels));
  const stops=[[255,253,248],[215,237,242],[30,101,133],[21,59,82]];
  const scaled=t*(stops.length-1),low=Math.floor(scaled),high=Math.min(stops.length-1,low+1),fraction=scaled-low;
  return `rgb(${stops[low].map((value,index)=>Math.round(value+(stops[high][index]-value)*fraction)).join(",")})`;
}
function legendItems(label,levels,kind){
  return levels.map((value,index)=>`<span class="annotation-key"><span class="annotation-swatch" style="background:${categoryColor(kind,index)}"></span>${esc(label)} ${esc(value)}</span>`).join("");
}
function readClusterPanel(row){
  const detail=row.detail;
  if(!detail){
    const explanation=row.evaluation_status==="EVALUABLE"
      ?"此可評估單元未達 source detail trigger（primary p>0.05 且 R²<0.05），因此沒有保存 read matrix；這是圖形不可用，不等於證明 reads 沒有聚集。"
      :"此單元未通過正式可評估條件，source 沒有產生 read matrix。";
    return `<p class="empty-note" data-testid="read-cluster-unavailable">${esc(explanation)}</p>`;
  }
  const cluster=detail.read_cluster;
  if(!cluster||!cluster.available){
    const explanation=cluster?.source_status==="SOURCE_NULL_N_GT_160"
      ?`Source detail 有 ${fmt(cluster.n_reads,0)} reads，超過固定 N≤${fmt(cluster.source_limit_n,0)} matrix 保存上限；不以 state mean 補造 read×read matrix。`
      :"Source detail 沒有可用 read×read matrix；不做推測性重建。";
    return `<p class="empty-note" data-testid="read-cluster-unavailable">${esc(explanation)}</p>`;
  }
  const patternLegend=legendItems("Pattern",cluster.pattern_levels,"pattern");
  const strandLegend=legendItems("Strand",cluster.strand_levels,"strand");
  const groupLegend=legendItems("RG",cluster.read_group_levels,"group");
  return `<div class="read-cluster-shell" data-testid="read-cluster-heatmap">
    <div class="cluster-summary-grid">
      <div><span>Reads</span><b>${fmt(cluster.n_reads,0)}</b></div>
      <div><span>Within-pattern mean</span><b>${fmt(cluster.within_pattern_mean)}</b></div>
      <div><span>Between-pattern mean</span><b>${fmt(cluster.between_pattern_mean)}</b></div>
      <div><span>Between − within</span><b>${fmt(cluster.between_minus_within)}</b></div>
      <div><span>Pattern blocks after order</span><b>${fmt(cluster.pattern_block_count,0)}</b></div>
    </div>
    <div class="cluster-legend"><strong>DISTANCE</strong><span class="distance-key">0 相似 <span class="distance-gradient" aria-hidden="true"></span> 1 不同</span><strong>ANNOTATIONS</strong>${patternLegend}${strandLegend}${groupLegend}</div>
    <div class="read-cluster-scroll">
      <div class="read-cluster-stage" id="readClusterStage">
        <canvas class="read-cluster-canvas" id="readClusterCanvas" data-testid="read-cluster-canvas" role="img" tabindex="0" aria-describedby="readClusterDescription readClusterLive"></canvas>
        <div class="cluster-tooltip" id="readClusterTooltip" role="tooltip"></div>
      </div>
    </div>
    <p id="readClusterDescription">排序法：label-independent UPGMA average linkage；完整 pattern、strand 與匿名 RG 類別只在排序完成後疊加。exact-distance ties 以 source leaf ordinal 作 deterministic tie-break，因此不主張任意 row permutation 完全 invariant。矩陣是描述性視覺，不設定 cluster cut，也不代表 clone、祖源或演化方向。</p>
    <p class="cluster-live" id="readClusterLive" aria-live="polite"></p>
  </div>`;
}
function drawReadCluster(row){
  const cluster=row.detail?.read_cluster,canvas=$("readClusterCanvas"),stage=$("readClusterStage"),tooltip=$("readClusterTooltip"),live=$("readClusterLive");
  if(!cluster?.available||!canvas||!stage||!tooltip||!live)return;
  const n=cluster.n_reads,matrix=decodeU8(cluster.matrix_u8_b64),patterns=decodeU8(cluster.pattern_codes_u8_b64),strands=decodeU8(cluster.strand_codes_u8_b64),groups=decodeU8(cluster.read_group_codes_u8_b64);
  if(matrix.length!==n*n||patterns.length!==n||strands.length!==n||groups.length!==n)throw new Error("read cluster payload length mismatch");
  const matrixSize=Math.max(420,Math.min(720,n*6)),left=164,top=122,right=34,bottom=58,width=Math.ceil(left+matrixSize+right),height=Math.ceil(top+matrixSize+bottom),cell=matrixSize/n,dpr=Math.max(1,window.devicePixelRatio||1);
  canvas.width=Math.ceil(width*dpr);canvas.height=Math.ceil(height*dpr);canvas.style.width=`${width}px`;canvas.style.height=`${height}px`;stage.style.width=`${width+24}px`;
  canvas.setAttribute("aria-label",`${row.dataset} ${row.chrom} HP ${row.hp_raw}；${n} reads 的 Bernoulli distance UPGMA 聚集圖；patterns ${cluster.pattern_levels.join("、")}`);
  let keyboardRow=0,keyboardColumn=0;
  const describe=(r,c)=>{
    const code=matrix[r*n+c],distance=code===cluster.na_code?"NA":(code/cluster.distance_levels).toFixed(3),first=cluster.pattern_levels[patterns[r]],second=cluster.pattern_levels[patterns[c]];
    return `r${String(r+1).padStart(3,"0")} (${first}) × r${String(c+1).padStart(3,"0")} (${second})；Bernoulli distance ${distance}；${first===second?"同 pattern":"不同 pattern"}`;
  };
  const paint=(cursor=null)=>{
    const context=canvas.getContext("2d");context.setTransform(dpr,0,0,dpr,0,0);context.clearRect(0,0,width,height);context.fillStyle="#fff";context.fillRect(0,0,width,height);
    const dendrogramTop=18,dendrogramBottom=82,maxHeight=Math.max(...cluster.dendrogram_merges.map(merge=>Number(merge[2])),1e-9),nodes=new Map();
    for(let index=0;index<n;index+=1)nodes.set(index,{x:left+(index+.5)*cell,y:dendrogramBottom});
    context.strokeStyle="#58676a";context.lineWidth=.8;
    cluster.dendrogram_merges.forEach((merge,index)=>{
      const first=nodes.get(Number(merge[0])),second=nodes.get(Number(merge[1])),mergeY=dendrogramBottom-(Number(merge[2])/maxHeight)*(dendrogramBottom-dendrogramTop);
      if(!first||!second)throw new Error("invalid dendrogram child id");
      context.beginPath();context.moveTo(first.x,first.y);context.lineTo(first.x,mergeY);context.lineTo(second.x,mergeY);context.lineTo(second.x,second.y);context.stroke();
      nodes.set(n+index,{x:(first.x+second.x)/2,y:mergeY});
    });
    for(let column=0;column<n;column+=1){context.fillStyle=categoryColor("pattern",patterns[column]);context.fillRect(left+column*cell,top-24,Math.ceil(cell+.2),10)}
    const stripX=[left-120,left-80,left-40],stripValues=[patterns,strands,groups],stripKinds=["pattern","strand","group"],stripWidth=24;
    for(let rowIndex=0;rowIndex<n;rowIndex+=1){
      stripValues.forEach((values,index)=>{context.fillStyle=categoryColor(stripKinds[index],values[rowIndex]);context.fillRect(stripX[index],top+rowIndex*cell,stripWidth,Math.ceil(cell+.2))});
      for(let column=0;column<n;column+=1){const code=matrix[rowIndex*n+column];context.fillStyle=distanceColor(code,cluster.distance_levels,cluster.na_code);context.fillRect(left+column*cell,top+rowIndex*cell,Math.ceil(cell+.25),Math.ceil(cell+.25))}
    }
    context.strokeStyle="rgba(159,77,47,.42)";context.lineWidth=1;
    for(let index=1;index<n;index+=1){
      if(patterns[index]!==patterns[index-1]){const position=left+index*cell;context.beginPath();context.moveTo(position,top);context.lineTo(position,top+matrixSize);context.stroke();context.beginPath();context.moveTo(left,top+index*cell);context.lineTo(left+matrixSize,top+index*cell);context.stroke()}
    }
    context.strokeStyle="#142327";context.lineWidth=1;context.strokeRect(left,top,matrixSize,matrixSize);
    context.fillStyle="#142327";context.font="600 12px sans-serif";context.textAlign="center";context.fillText("Reads（UPGMA ordered）",left+matrixSize/2,height-17);
    context.save();context.translate(width-10,top+matrixSize/2);context.rotate(Math.PI/2);context.fillText("Reads（same order）",0,0);context.restore();
    context.font="600 10px monospace";context.textAlign="center";["R/A","Strand","RG"].forEach((label,index)=>context.fillText(label,stripX[index]+stripWidth/2,top-7));
    context.textAlign="left";context.fillStyle="#58676a";context.fillText(`UPGMA · N=${n}`,14,16);
    if(cursor){context.strokeStyle="#c48a2b";context.lineWidth=2;context.strokeRect(left+cursor.column*cell,top+cursor.row*cell,Math.max(2,cell),Math.max(2,cell))}
  };
  const locate=event=>{
    const rect=canvas.getBoundingClientRect(),x=(event.clientX-rect.left)*(width/rect.width),y=(event.clientY-rect.top)*(height/rect.height),column=Math.floor((x-left)/cell),rowIndex=Math.floor((y-top)/cell);
    return rowIndex>=0&&rowIndex<n&&column>=0&&column<n?{row:rowIndex,column,x,y}:null;
  };
  canvas.addEventListener("mousemove",event=>{const cellInfo=locate(event);if(!cellInfo){tooltip.style.display="none";return}tooltip.textContent=describe(cellInfo.row,cellInfo.column);tooltip.style.display="block";tooltip.style.left=`${Math.min(width-310,Math.max(4,cellInfo.x+18))}px`;tooltip.style.top=`${Math.min(height-64,Math.max(4,cellInfo.y+18))}px`});
  canvas.addEventListener("mouseleave",()=>{tooltip.style.display="none"});
  canvas.addEventListener("keydown",event=>{
    const moves={ArrowUp:[-1,0],ArrowDown:[1,0],ArrowLeft:[0,-1],ArrowRight:[0,1]};if(!moves[event.key])return;event.preventDefault();keyboardRow=Math.max(0,Math.min(n-1,keyboardRow+moves[event.key][0]));keyboardColumn=Math.max(0,Math.min(n-1,keyboardColumn+moves[event.key][1]));paint({row:keyboardRow,column:keyboardColumn});live.textContent=describe(keyboardRow,keyboardColumn);
  });
  canvas.addEventListener("focus",()=>{paint({row:keyboardRow,column:keyboardColumn});live.textContent=`方向鍵瀏覽矩陣。${describe(keyboardRow,keyboardColumn)}`});
  canvas.addEventListener("blur",()=>{paint();live.textContent=""});
  paint();
}
function ordinalMarkerX(position,positions,left,cell){
  if(!positions.length)return left;
  const center=index=>left+index*cell+(cell-1)/2;
  if(position<positions[0])return center(0)-12;
  if(position>positions[positions.length-1])return center(positions.length-1)+12;
  if(position===positions[0])return center(0);
  for(let index=1;index<positions.length;index+=1){
    if(position<=positions[index]){
      const low=positions[index-1],high=positions[index];
      const fraction=high===low?0:(position-low)/(high-low);
      return center(index-1)+fraction*(center(index)-center(index-1));
    }
  }
  return center(positions.length-1);
}
function profileHeatmap(row){
  const detail=row.detail;if(!detail)return `<p class="empty-note">此單元未產生 details record；不顯示推測 profile。</p>`;
  const positions=detail.cpg_positions_displayed,states=Object.keys(detail.state_mean_profiles),cell=18,left=112,top=58,width=Math.max(720,left+positions.length*cell+24),height=top+states.length*34+52;
  let svg=`<svg class="profile-svg" viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="profileTitle profileDesc"><title id="profileTitle">State mean methylation profile heatmap</title><desc id="profileDesc">${states.length} states across ${positions.length} displayed CpGs; values range from zero to one.</desc>`;
  row.active_positions.forEach(position=>{const x=ordinalMarkerX(position,positions,left,cell);const outside=position<(positions[0]??position)||position>(positions[positions.length-1]??position);svg+=`<path d="M${x} 27 l-5 9 h10 z" fill="#9f4d2f"><title>active sSNV ${position}${outside?" · outside displayed CpG span":""}</title></path>`});
  states.forEach((state,r)=>{const y=top+r*34;svg+=`<text x="${left-10}" y="${y+14}" text-anchor="end" font-family="monospace" font-size="12" fill="#142327">${esc(state)} · n=${fmt(detail.state_counts[state],0)}</text>`;detail.state_mean_profiles[state].forEach((value,c)=>{const x=left+c*cell;svg+=`<rect x="${x}" y="${y}" width="${cell-1}" height="24" fill="${colorFor(value)}" stroke="#fff"><title>${esc(state)} · CpG ${positions[c]} · mean ${fmt(value)}</title></rect>`})});
  [0,Math.floor((positions.length-1)/2),positions.length-1].filter((v,i,a)=>positions.length&&a.indexOf(v)===i).forEach(index=>{svg+=`<text x="${left+index*cell+(cell-1)/2}" y="${height-14}" text-anchor="middle" font-family="monospace" font-size="10" fill="#58676a">${positions[index].toLocaleString("zh-TW")}</text>`});
  svg+=`<rect x="${left}" y="${height-38}" width="70" height="8" fill="#efe9d8"/><rect x="${left+70}" y="${height-38}" width="70" height="8" fill="#1e6585"/><text x="${left+150}" y="${height-30}" font-size="10" fill="#58676a">mean methylation 0 → 1</text></svg>`;
  const table=`<table class="sr-only"><caption>State mean methylation values</caption><thead><tr><th>State</th>${positions.map(p=>`<th>${p}</th>`).join("")}</tr></thead><tbody>${states.map(state=>`<tr><th>${esc(state)}</th>${detail.state_mean_profiles[state].map(value=>`<td>${fmt(value)}</td>`).join("")}</tr>`).join("")}</tbody></table>`;
  return `<div class="matrix-scroll" data-testid="profile-heatmap">${svg}${table}</div>`;
}
function pairMatrix(row){
  const detail=row.detail;if(!detail)return `<p class="empty-note">此單元沒有可注入的 pairwise detail。</p>`;
  const states=Object.keys(detail.state_mean_profiles),lookup=new Map();
  detail.pairwise_effects.forEach(effect=>{lookup.set(`${effect.first}|${effect.second}`,effect);lookup.set(`${effect.second}|${effect.first}`,effect)});
  const body=states.map(a=>`<tr><th scope="row">${esc(a)}</th>${states.map(b=>{if(a===b)return `<td>—<small>同 state</small></td>`;const effect=lookup.get(`${a}|${b}`);if(!effect)return `<td>NA</td>`;const relation=relationZh[effect.topology_relation]||effect.topology_relation;return `<td data-hamming="${effect.hamming>1?"gt1":"one"}" style="background:${effect.hamming>1?"#f4ddd2":"#d7edf2"}"><strong>${fmt(effect.distance_contrast)}</strong><small>H=${effect.hamming} · ${esc(relation)}</small></td>`}).join("")}</tr>`).join("");
  return `<div class="matrix-scroll" data-testid="state-pair-matrix"><table class="state-matrix"><caption class="sr-only">State pair Bernoulli distance contrast matrix</caption><thead><tr><th>State</th>${states.map(state=>`<th>${esc(state)}</th>`).join("")}</tr></thead><tbody>${body}</tbody></table></div>`;
}
function topologyRelation(row){
  const pair=(row.best_pair||"|").split("|"),relation=row.best_pair_topology_relation||"HAMMING1_TOPOLOGY_UNAVAILABLE",gt=Number(row.best_pair_hamming)>1;
  if(!pair[0]||!pair[1])return `<p class="empty-note">沒有可顯示的 best state pair。</p>`;
  const dashed=gt||relation!=="HAMMING1_GLOBAL_BEST_UNANIMOUS";
  const label=relationZh[relation]||relation;
  return `<div class="matrix-scroll" data-testid="topology-relation"><svg class="relation-svg" viewBox="0 0 760 190" role="img" aria-labelledby="relationTitle relationDesc"><title id="relationTitle">Best state pair and frozen topology undirected membership</title><desc id="relationDesc">${esc(label)}</desc><rect x="60" y="58" width="210" height="64" rx="9" fill="#d7edf2" stroke="#153b52" stroke-width="2"/><rect x="490" y="58" width="210" height="64" rx="9" fill="#f7ecd2" stroke="#9f4d2f" stroke-width="2"/><text x="165" y="97" text-anchor="middle" font-family="monospace" font-size="22" fill="#142327">${esc(pair[0])}</text><text x="595" y="97" text-anchor="middle" font-family="monospace" font-size="22" fill="#142327">${esc(pair[1])}</text><line x1="270" y1="90" x2="490" y2="90" stroke="${gt?"#9f4d2f":"#1e6585"}" stroke-width="4" ${dashed?'stroke-dasharray="10 8"':""}/><text x="380" y="70" text-anchor="middle" font-family="monospace" font-size="12" fill="#58676a">Hamming ${fmt(row.best_pair_hamming,0)} · undirected</text><text x="380" y="154" text-anchor="middle" font-size="13" fill="#142327">${esc(label)}</text></svg></div>`;
}
function selectCase(id,scroll){
  const row=caseById.get(id);if(!row)return;
  const reason=row.invalid_reason?`<p class="empty-note"><strong>未評估原因：</strong>${esc(row.invalid_reason)}</p>`:"";
  $("caseDetail").innerHTML=`<div class="detail-head"><div><span class="badge" data-assessment="${esc(row.assessment)}">${esc(assessmentZh[row.assessment]||row.assessment)}</span>${resolutionFlag(row)}<h3>${esc(row.dataset)} · ${esc(row.chrom)} · ${esc(row.region_id)} · HP ${esc(row.hp_raw)}</h3><p>${esc(row.patterns.join(" / ")||"No analyzed complete states")}</p></div><span class="chip">active k=${fmt(row.n_active_bits,0)}</span></div>${reason}<div class="detail-metrics"><div><span>Analysis N</span><b>${fmt(row.analysis_n,0)}</b></div><div><span>Common CpGs</span><b>${fmt(row.n_common_cpg,0)}</b></div><div><span>PERMANOVA R²</span><b>${fmt(row.permanova_r2)}</b></div><div><span>BY q</span><b>${fmt(row.q_by)}</b></div><div><span>Permutations</span><b>${fmt(row.permanova_permutations_requested,0)}</b></div><div><span>PERMDISP p</span><b>${fmt(row.permdisp_p)}</b></div><div><span>Geometry max SMD</span><b>${fmt(row.max_geometry_smd)}</b></div></div><div class="visual-grid"><section class="visual"><h3>Read × read Bernoulli 聚集關係</h3><p>延續舊觀察圖的 UPGMA dendrogram＋對稱 distance heatmap＋annotation strips。排序只使用 read 距離；完整 R/A pattern、strand、匿名 RG 類別在排序後疊加。淺色=相似、深色=不同；不切群、不推斷 clone 或演化方向。</p>${readClusterPanel(row)}</section><section class="visual"><h3>State mean methylation profile</h3><p>每列是一個 complete R/A state；最多顯示 ${MAX_PROFILE_CPG} 個等距取樣 CpG。三角以 genomic position 在同一 CpG ordinal 軸分段內插；超出顯示 span 時置於邊界外。色階是 0–1 機率，不代表效應顯著。</p>${profileHeatmap(row)}</section><section class="visual"><h3>State-pair Bernoulli distance contrast</h3><p>格內為 between − pooled within；虛線框表示 Hamming&gt;1，只能解讀為 pair band。</p>${pairMatrix(row)}</section><section class="visual"><h3>Best pair × frozen topology relation</h3><p>所有 state-pair topology relation 一律以無方向 edge membership 或 pair band 顯示，不推斷先後或演化方向。</p>${topologyRelation(row)}</section></div>`;
  drawReadCluster(row);
  if(scroll){$("caseDetail").focus({preventScroll:true});$("caseDetail").scrollIntoView({behavior:"smooth",block:"start"})}
}
init();
</script>
</body>
</html>
"""


def render_html(report_data: Mapping[str, Any]) -> str:
    embedded = json.dumps(
        report_data,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).replace("</", "<\\/")
    return HTML_TEMPLATE.replace("__REPORT_DATA_JSON__", embedded)


def write_outputs(
    report_data: Mapping[str, Any], output_html: Path, output_data_json: Path
) -> None:
    if output_html.resolve() == output_data_json.resolve():
        raise ReportContractError("HTML and report-data JSON paths must differ")
    for path in (output_html, output_data_json):
        if path.exists():
            raise ReportContractError(f"refusing to overwrite output: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
    html_partial = output_html.with_suffix(output_html.suffix + ".partial")
    data_partial = output_data_json.with_suffix(output_data_json.suffix + ".partial")
    for path in (html_partial, data_partial):
        if path.exists():
            raise ReportContractError(f"staging output already exists: {path}")
    data_text = json.dumps(
        report_data,
        ensure_ascii=False,
        sort_keys=True,
        indent=2,
        allow_nan=False,
    ) + "\n"
    html_text = render_html(report_data)
    with data_partial.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(data_text)
    with html_partial.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(html_text)
    data_partial.rename(output_data_json)
    html_partial.rename(output_html)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census-receipt", type=Path, required=True)
    parser.add_argument("--pattern-counts", type=Path, required=True)
    parser.add_argument("--artifact-catalog-summary", type=Path, required=True)
    parser.add_argument("--artifact-catalog-tsv", type=Path, required=True)
    parser.add_argument("--bernoulli-parity-summary", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--details", type=Path, required=True)
    parser.add_argument("--analysis-summary", type=Path, required=True)
    parser.add_argument("--output-html", type=Path, required=True)
    parser.add_argument("--output-data-json", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    data = build_report_data(
        census_receipt_path=args.census_receipt.resolve(),
        pattern_counts_path=args.pattern_counts.resolve(),
        catalog_summary_path=args.artifact_catalog_summary.resolve(),
        catalog_tsv_path=args.artifact_catalog_tsv.resolve(),
        bernoulli_parity_summary_path=args.bernoulli_parity_summary.resolve(),
        evidence_path=args.evidence.resolve(),
        details_path=args.details.resolve(),
        analysis_summary_path=args.analysis_summary.resolve(),
    )
    write_outputs(data, args.output_html.resolve(), args.output_data_json.resolve())
    print(
        json.dumps(
            {
                "status": "PASS",
                "output_html": str(args.output_html.resolve()),
                "output_data_json": str(args.output_data_json.resolve()),
                "analysis_units": data["aggregate"]["analysis_units"],
                "detail_records": data["aggregate"]["detail_records"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
