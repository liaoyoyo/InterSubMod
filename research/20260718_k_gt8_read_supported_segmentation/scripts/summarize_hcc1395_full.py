#!/usr/bin/env python3
"""Validate and summarize the canonical HCC1395 chr1-22 segmentation run.

This is an independent, read-only validator for the full-run evidence tree.
It authenticates the runner contract/receipt, every chromosome stage receipt,
all 21 partition child receipts, and every partition TSV identity before it
creates a new summary directory.  Metrics are re-derived from the TSV rows;
receipt counters are used as cross-checks, not as the sole source of truth.

The output directory must not exist.  The script never overwrites source or
summary artifacts.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from math import floor, isfinite
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


DATASET = "HCC1395"
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
RUN_SCHEMA = "intersubmod.hcc1395_full_k_gt8_segmentation"
RUN_SCHEMA_VERSION = "1.0.0"
PARTITION_SCHEMA = "intersubmod.k_gt8_read_supported_segmentation"
PARTITION_SCHEMA_VERSION = "0.1.0"
SUMMARY_SCHEMA = "intersubmod.hcc1395_full_k_gt8_segmentation.summary"
SUMMARY_SCHEMA_VERSION = "1.1.0"
OUTER_STAGE_WALL_TOLERANCE_SECONDS = 0.25
UNAVOIDABLE = "unavoidable_span_gt_max_block_size"
RETAINED = "retained"
CUT = "cut"

CANONICAL_EXPECTED: dict[str, tuple[int, int, int, int]] = {
    "chr1": (3045, 13, 146, 18),
    "chr2": (3074, 17, 253, 49),
    "chr3": (2381, 9, 103, 23),
    "chr4": (2598, 21, 245, 23),
    "chr5": (2305, 13, 146, 17),
    "chr6": (27099, 83, 25657, 3574),
    "chr7": (2845, 46, 523, 23),
    "chr8": (3178, 59, 1205, 153),
    "chr9": (1523, 6, 80, 22),
    "chr10": (1504, 9, 106, 15),
    "chr11": (1328, 2, 28, 18),
    "chr12": (1573, 13, 160, 43),
    "chr13": (1019, 5, 52, 12),
    "chr14": (1343, 10, 156, 56),
    "chr15": (1154, 9, 123, 28),
    "chr16": (18819, 63, 18162, 2674),
    "chr17": (1093, 10, 188, 65),
    "chr18": (991, 5, 53, 13),
    "chr19": (761, 3, 27, 9),
    "chr20": (1075, 6, 59, 11),
    "chr21": (436, 0, 0, 0),
    "chr22": (543, 6, 98, 40),
}

COMPONENT_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "legacy_group_ordinal",
    "start1",
    "end1",
    "span_bp",
    "pre_cap_k",
    "old_densest8_selected",
    "old_cap_excluded",
    "new_block_count",
    "new_site_retained",
    "primary_active_site_count",
    "primary_active_site_fraction",
    "exact_pattern_count",
    "raw_total_molecule_weight",
    "raw_retained_molecule_weight",
    "raw_lost_molecule_weight",
    "raw_retention_ratio",
    "old_densest8_retained_molecule_weight",
    "old_densest8_retention_ratio",
    "retained_exact_pattern_count",
    "lost_exact_pattern_count",
    "unavoidable_pattern_count",
    "old_densest8_retained_pattern_count",
    "weight_stable",
    "raw_cut_indices",
    "equal_pattern_cut_indices",
    "log1p_cut_indices",
    "cut_gap_sum_bp",
    "status",
    "positions_sha256",
)

BLOCK_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "block_index",
    "block_id",
    "start1",
    "end1",
    "k",
    "positions",
    "component_status",
    "raw_retained_molecule_weight",
    "retained_exact_pattern_count",
    "primary_active_site_count",
)

BLOCK_ALL_FIELDS = (*BLOCK_FIELDS, "evidence_gate")

MEMBERSHIP_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "site_index",
    "component_local_index",
    "pos1",
    "ref",
    "alt",
    "block_id",
    "primary_linkage_observed",
    "old_densest8_disposition",
)

CONSTRAINT_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "constraint_id",
    "unit_id",
    "hp_family",
    "phase_set",
    "positions",
    "call_codes",
    "n_fixed_ra",
    "span_sites",
    "molecule_weight",
    "disposition",
    "crossed_cut_count",
    "retained_block_index",
)

RUNNER_SUMMARY_FIELDS = (
    "chrom",
    "ssnv_sites",
    "k_gt8_components",
    "k_gt8_sites",
    "k_gt8_max_k",
    "extraction_status",
    "extraction_wall_seconds",
    "extraction_stage_receipt",
    "extraction_stage_receipt_sha256",
    "partition_status",
    "partition_wall_seconds",
    "partition_stage_receipt",
    "partition_stage_receipt_sha256",
)

SUMMARY_FIELDS = (
    "chrom",
    "partition_stage_status",
    "ssnv_sites",
    "components",
    "sites",
    "old_selected_sites",
    "old_excluded_sites",
    "new_blocks",
    "new_retained_sites",
    "primary_active_sites",
    "exact_patterns",
    "raw_total_molecule_weight",
    "new_retained_molecule_weight",
    "new_lost_molecule_weight",
    "new_retention_ratio",
    "old_densest8_retained_molecule_weight",
    "old_densest8_retention_ratio",
    "retained_exact_patterns",
    "lost_exact_patterns",
    "unavoidable_patterns",
    "unavoidable_molecule_weight",
    "unavoidable_n_fixed_ra_gt8_patterns",
    "unavoidable_n_fixed_ra_gt8_molecule_weight",
    "unavoidable_n_fixed_ra_lte8_span_gt8_patterns",
    "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight",
    "old_densest8_retained_patterns",
    "weight_stable_components",
    "weight_sensitive_components",
    "zero_evidence_blocks",
    "zero_evidence_block_sites",
    "tree_ready_blocks",
    "tree_ready_block_sites",
    "abstain_blocks",
    "abstain_block_sites",
    "status_counts_json",
    "block_k_distribution_json",
    "block_span_bp_distribution_json",
    "component_retention_delta_distribution_json",
    "raw_overlapping_alignments",
    "eligible_alignments_pre_identity_collapse",
    "duplicate_alignment_identities_collapsed",
    "canonical_unique_molecules",
    "molecule_sparse_rows",
    "sparse_site_calls",
    "known_ps_hp12_molecule_rows",
    "known_ps_hp1_molecule_rows",
    "known_ps_hp2_molecule_rows",
    "extraction_failed_checks",
    "duplicate_identity_conflicts",
    "extraction_wall_seconds",
    "partition_wall_seconds",
    "partition_pattern_load_aggregate_seconds",
    "partition_ordered_hypergraph_dp_seconds",
    "extraction_max_rss_kb",
    "partition_max_rss_kb",
)

EXTRACTION_METRIC_FIELDS = (
    "raw_overlapping_alignments",
    "eligible_alignments_pre_identity_collapse",
    "duplicate_alignment_identities_collapsed",
    "canonical_unique_molecules",
    "molecule_sparse_rows",
    "sparse_site_calls",
    "known_ps_hp12_molecule_rows",
    "known_ps_hp1_molecule_rows",
    "known_ps_hp2_molecule_rows",
    "extraction_failed_checks",
    "duplicate_identity_conflicts",
)


class ContractError(RuntimeError):
    """Raised when source evidence or a conservation invariant fails."""


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def semantic_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def require_regular_file(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_file():
        raise ContractError(f"{label} is missing, nonregular, or a symlink: {path}")


def require_real_directory(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_dir():
        raise ContractError(f"{label} is missing, not a directory, or a symlink: {path}")


def file_identity(path: Path) -> dict[str, Any]:
    require_regular_file(path, "file")
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def strict_json_load(path: Path) -> Any:
    require_regular_file(path, "JSON")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ContractError(f"duplicate JSON key in {path}: {key}")
            result[key] = value
        return result

    try:
        return json.loads(
            path.read_text(encoding="utf-8", errors="strict"),
            object_pairs_hook=reject_duplicates,
            parse_constant=lambda value: (_ for _ in ()).throw(
                ContractError(f"non-finite JSON value in {path}: {value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"cannot parse strict JSON {path}: {exc}") from exc


def verify_sha256_sidecar(path: Path) -> str:
    require_regular_file(path, "SHA-protected file")
    sidecar = path.with_name(f"{path.name}.sha256")
    require_regular_file(sidecar, "SHA-256 sidecar")
    try:
        text = sidecar.read_text(encoding="ascii", errors="strict")
    except (OSError, UnicodeError) as exc:
        raise ContractError(f"cannot read SHA-256 sidecar {sidecar}: {exc}") from exc
    match = re.fullmatch(r"([0-9a-f]{64})  ([^\n]+)\n?", text)
    if not match or match.group(2) != path.name:
        raise ContractError(f"malformed SHA-256 sidecar: {sidecar}")
    observed = sha256_path(path)
    if match.group(1) != observed:
        raise ContractError(
            f"SHA-256 mismatch for {path}: sidecar={match.group(1)} observed={observed}"
        )
    return observed


def verify_identity(
    raw: Any,
    expected_path: Path,
    label: str,
    *,
    require_sidecar: bool = False,
) -> dict[str, Any]:
    if not isinstance(raw, Mapping) or set(raw) != {"path", "size_bytes", "sha256"}:
        raise ContractError(f"{label} identity must contain path/size_bytes/sha256")
    observed = file_identity(expected_path)
    if observed != dict(raw):
        raise ContractError(
            f"{label} identity drift: expected={dict(raw)} observed={observed}"
        )
    if require_sidecar:
        sidecar_sha = verify_sha256_sidecar(expected_path)
        if sidecar_sha != observed["sha256"]:
            raise ContractError(f"{label} sidecar disagrees with identity")
    return observed


def _int(value: str, label: str, *, minimum: int | None = None) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{label} is not an integer: {value!r}") from exc
    if minimum is not None and result < minimum:
        raise ContractError(f"{label} must be >= {minimum}: {result}")
    return result


def _decimal(value: str, label: str, *, minimum: Decimal | None = None) -> Decimal:
    try:
        result = Decimal(value)
    except Exception as exc:
        raise ContractError(f"{label} is not a finite decimal: {value!r}") from exc
    if not result.is_finite():
        raise ContractError(f"{label} is not finite: {value!r}")
    if minimum is not None and result < minimum:
        raise ContractError(f"{label} must be >= {minimum}: {result}")
    return result


def _csv_ints(value: str, label: str) -> tuple[int, ...]:
    if not value:
        return ()
    result = tuple(_int(token, label) for token in value.split(","))
    if len(result) != len(set(result)):
        raise ContractError(f"{label} contains duplicate positions")
    return result


def _ratio(numerator: Decimal, denominator: Decimal) -> str:
    if denominator == 0:
        return ""
    return format(numerator / denominator, "f")


def _json_compact(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _distribution(
    values: Sequence[int | Decimal],
    *,
    total_count: int | None = None,
) -> dict[str, Any]:
    """Return deterministic R-7 linear-interpolated p10/p50/p90 summaries."""

    ordered = sorted(Decimal(value) for value in values)
    total = len(ordered) if total_count is None else total_count
    if total < len(ordered):
        raise ContractError("distribution total_count is smaller than defined values")
    if not ordered:
        return {
            "n_total": total,
            "n_defined": 0,
            "n_undefined": total,
            "min": "",
            "p10": "",
            "median": "",
            "p90": "",
            "max": "",
            "mean": "",
        }

    def percentile(probability: Decimal) -> Decimal:
        if len(ordered) == 1:
            return ordered[0]
        coordinate = Decimal(len(ordered) - 1) * probability
        lower = floor(coordinate)
        fraction = coordinate - lower
        if fraction == 0:
            return ordered[lower]
        return ordered[lower] + (ordered[lower + 1] - ordered[lower]) * fraction

    return {
        "n_total": total,
        "n_defined": len(ordered),
        "n_undefined": total - len(ordered),
        "min": format(ordered[0], "f"),
        "p10": format(percentile(Decimal("0.10")), "f"),
        "median": format(percentile(Decimal("0.50")), "f"),
        "p90": format(percentile(Decimal("0.90")), "f"),
        "max": format(ordered[-1], "f"),
        "mean": format(sum(ordered, Decimal(0)) / len(ordered), "f"),
    }


def _read_tsv_rows(path: Path, expected_fields: Sequence[str]) -> Iterable[dict[str, str]]:
    require_regular_file(path, "gzip TSV")
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != tuple(expected_fields):
                raise ContractError(
                    f"TSV header drift in {path}: expected={tuple(expected_fields)} "
                    f"observed={tuple(reader.fieldnames or ())}"
                )
            for row_number, row in enumerate(reader, start=2):
                if None in row or set(row) != set(expected_fields):
                    raise ContractError(f"malformed TSV row {row_number} in {path}")
                if any(value is None for value in row.values()):
                    raise ContractError(f"missing TSV value at row {row_number} in {path}")
                yield dict(row)
    except (OSError, UnicodeError, csv.Error, gzip.BadGzipFile) as exc:
        raise ContractError(f"cannot read gzip TSV {path}: {exc}") from exc


def _validate_all_true(raw: Any, label: str) -> None:
    if not isinstance(raw, Mapping) or not raw:
        raise ContractError(f"{label} must be a nonempty check map")
    failures = [name for name, value in raw.items() if value is not True]
    if failures:
        raise ContractError(f"{label} contains failed checks: {failures}")


def parse_gnu_time_v(path: Path) -> dict[str, Any]:
    """Parse command identity and stable fields from GNU ``/usr/bin/time -v``."""

    require_regular_file(path, "GNU time resource file")
    fields: dict[str, str] = {}
    try:
        for raw_line in path.read_text(
            encoding="utf-8", errors="strict"
        ).splitlines():
            line = raw_line.strip()
            if ": " not in line:
                continue
            key, value = line.split(": ", 1)
            if key in fields:
                raise ContractError(f"duplicate GNU time field in {path}: {key}")
            fields[key] = value
    except (OSError, UnicodeError) as exc:
        raise ContractError(f"cannot read GNU time output {path}: {exc}") from exc

    required = {
        "Command being timed",
        "User time (seconds)",
        "System time (seconds)",
        "Percent of CPU this job got",
        "Elapsed (wall clock) time (h:mm:ss or m:ss)",
        "Maximum resident set size (kbytes)",
        "File system inputs",
        "File system outputs",
        "Exit status",
    }
    missing = sorted(required - set(fields))
    if missing:
        raise ContractError(f"GNU time output lacks fields in {path}: {missing}")

    recorded_command = fields["Command being timed"].strip()
    if (
        len(recorded_command) >= 2
        and recorded_command[0] == '"'
        and recorded_command[-1] == '"'
    ):
        command_text = recorded_command[1:-1]
    else:
        command_text = recorded_command
    try:
        command_argv = shlex.split(command_text)
    except ValueError as exc:
        raise ContractError(
            f"invalid GNU time command field in {path}: {recorded_command}"
        ) from exc
    if not command_argv:
        raise ContractError(f"GNU time command field is empty in {path}")

    elapsed_text = fields["Elapsed (wall clock) time (h:mm:ss or m:ss)"]
    parts = elapsed_text.split(":")
    try:
        if len(parts) == 2:
            elapsed_seconds = int(parts[0]) * 60 + float(parts[1])
        elif len(parts) == 3:
            elapsed_seconds = (
                int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
            )
        else:
            raise ValueError("wrong field count")
    except ValueError as exc:
        raise ContractError(f"invalid GNU elapsed time in {path}: {elapsed_text}") from exc

    cpu_text = fields["Percent of CPU this job got"]
    if not cpu_text.endswith("%"):
        raise ContractError(f"invalid CPU percentage in {path}: {cpu_text}")
    result = {
        "path": str(path.resolve()),
        "sha256": sha256_path(path),
        "command_text": command_text,
        "command_argv": command_argv,
        "user_seconds": float(fields["User time (seconds)"]),
        "system_seconds": float(fields["System time (seconds)"]),
        "cpu_percent": float(cpu_text[:-1]),
        "elapsed_text": elapsed_text,
        "elapsed_seconds": elapsed_seconds,
        "max_rss_kb": _int(
            fields["Maximum resident set size (kbytes)"],
            f"{path} maximum RSS",
            minimum=0,
        ),
        "filesystem_inputs": _int(
            fields["File system inputs"], f"{path} filesystem inputs", minimum=0
        ),
        "filesystem_outputs": _int(
            fields["File system outputs"], f"{path} filesystem outputs", minimum=0
        ),
        "exit_status": _int(fields["Exit status"], f"{path} exit status"),
    }
    if result["exit_status"] != 0:
        raise ContractError(f"GNU time reports nonzero exit status in {path}")
    if any(
        result[key] < 0
        for key in ("user_seconds", "system_seconds", "cpu_percent", "elapsed_seconds")
    ):
        raise ContractError(f"GNU time contains negative resource value in {path}")
    return result


def _recorded_absolute_path(token: str, label: str) -> Path:
    """Resolve a command token only when its provenance is location-complete."""

    candidate = Path(token)
    if not candidate.is_absolute():
        raise ContractError(f"{label} must be an absolute path in GNU time command")
    return candidate.resolve(strict=False)


def _single_runner_option(
    runner_args: Sequence[str],
    option: str,
) -> str:
    """Read one path-valued runner option while rejecting ambiguity."""

    values: list[str] = []
    for index, token in enumerate(runner_args):
        if token == option:
            if index + 1 >= len(runner_args):
                raise ContractError(f"outer GNU time command has no value for {option}")
            values.append(runner_args[index + 1])
        elif token.startswith(f"{option}="):
            values.append(token.split("=", 1)[1])
    if len(values) != 1 or not values[0]:
        raise ContractError(
            f"outer GNU time command must contain exactly one {option}; "
            f"observed={len(values)}"
        )
    return values[0]


def validate_outer_time_contract(
    outer_resource: Mapping[str, Any],
    *,
    full_root: Path,
    runner_path: Path,
    manifest_path: Path,
    sequential_stage_wall_seconds: float,
    tolerance_seconds: float = OUTER_STAGE_WALL_TOLERANCE_SECONDS,
) -> dict[str, Any]:
    """Bind outer GNU time to this fresh run and enforce its wall-time lower bound."""

    argv = outer_resource.get("command_argv")
    if (
        not isinstance(argv, list)
        or not argv
        or any(not isinstance(token, str) or not token for token in argv)
    ):
        raise ContractError("outer GNU time command argv is missing or invalid")

    expected_runner = runner_path.resolve(strict=False)
    runner_indices = [
        index
        for index, token in enumerate(argv)
        if Path(token).is_absolute()
        and Path(token).resolve(strict=False) == expected_runner
    ]
    if len(runner_indices) != 1:
        raise ContractError(
            "outer GNU time command must bind exactly once to the contracted runner; "
            f"observed={len(runner_indices)}"
        )
    runner_args = argv[runner_indices[0] + 1 :]
    if "--resume" in runner_args or any(
        token.startswith("--resume=") for token in runner_args
    ):
        raise ContractError(
            "formal outer GNU time must cover a fresh runner command without --resume"
        )

    recorded_manifest = _recorded_absolute_path(
        _single_runner_option(runner_args, "--manifest"),
        "outer GNU time --manifest",
    )
    recorded_output_root = _recorded_absolute_path(
        _single_runner_option(runner_args, "--output-root"),
        "outer GNU time --output-root",
    )
    if recorded_manifest != manifest_path.resolve(strict=False):
        raise ContractError(
            "outer GNU time --manifest does not match the frozen run contract"
        )
    if recorded_output_root != full_root.resolve(strict=False):
        raise ContractError(
            "outer GNU time --output-root does not match the validated full-run root"
        )

    elapsed_raw = outer_resource.get("elapsed_seconds")
    if (
        isinstance(elapsed_raw, bool)
        or not isinstance(elapsed_raw, (int, float))
        or not isfinite(float(elapsed_raw))
        or float(elapsed_raw) < 0
    ):
        raise ContractError("outer GNU time elapsed_seconds is missing or invalid")
    if (
        isinstance(sequential_stage_wall_seconds, bool)
        or not isinstance(sequential_stage_wall_seconds, (int, float))
        or not isfinite(float(sequential_stage_wall_seconds))
        or float(sequential_stage_wall_seconds) < 0
    ):
        raise ContractError("sequential stage wall sum is invalid")
    if (
        isinstance(tolerance_seconds, bool)
        or not isinstance(tolerance_seconds, (int, float))
        or not isfinite(float(tolerance_seconds))
        or float(tolerance_seconds) < 0
    ):
        raise ContractError("outer-stage wall tolerance is invalid")

    elapsed_seconds = float(elapsed_raw)
    stage_wall_seconds = float(sequential_stage_wall_seconds)
    tolerance = float(tolerance_seconds)
    outer_minus_stage = elapsed_seconds - stage_wall_seconds
    if outer_minus_stage < -tolerance:
        raise ContractError(
            "outer GNU time elapsed is shorter than the sequential extraction + "
            "partition stage wall sum beyond tolerance: "
            f"outer={elapsed_seconds:.6f}s stages={stage_wall_seconds:.6f}s "
            f"tolerance={tolerance:.6f}s"
        )
    return {
        "command_binding_verified": True,
        "fresh_non_resume_command_verified": True,
        "sequential_stage_wall_seconds": stage_wall_seconds,
        "lower_bound_tolerance_seconds": tolerance,
        "outer_minus_sequential_stage_wall_seconds": outer_minus_stage,
        "runner_overhead_seconds": max(0.0, outer_minus_stage),
    }


def _verify_output_map(
    raw_outputs: Any,
    expected_paths: Mapping[str, Path],
    label: str,
) -> dict[str, dict[str, Any]]:
    if not isinstance(raw_outputs, Mapping) or set(raw_outputs) != set(expected_paths):
        raise ContractError(
            f"{label} output keys drift: expected={sorted(expected_paths)} "
            f"observed={sorted(raw_outputs) if isinstance(raw_outputs, Mapping) else raw_outputs}"
        )
    return {
        name: verify_identity(raw_outputs[name], path, f"{label} output {name}")
        for name, path in expected_paths.items()
    }


def _expected_status(
    *,
    pattern_count: int,
    lost_pattern_count: int,
    active_site_count: int,
    component_k: int,
    weight_stable: bool,
) -> str:
    if pattern_count == 0:
        return "ABSTAIN_NO_PRIMARY_LINKAGE"
    if lost_pattern_count:
        return (
            "PARTIAL_LOCAL_ONLY_WEIGHT_STABLE"
            if weight_stable
            else "PARTIAL_LOCAL_ONLY_WEIGHT_SENSITIVE"
        )
    coverage = (
        "FULL_SITE_COVERAGE"
        if active_site_count == component_k
        else "PARTIAL_SITE_COVERAGE"
    )
    stability = "WEIGHT_STABLE" if weight_stable else "WEIGHT_SENSITIVE"
    return f"OBSERVED_CONSTRAINT_NO_LOSS_{coverage}_{stability}"


def block_evidence_gate(
    retained_exact_pattern_count: int,
    raw_retained_molecule_weight: int,
    primary_active_site_count: int,
) -> str:
    """Return the conservative gate into local candidate-tree inference."""

    if (
        retained_exact_pattern_count > 0
        and raw_retained_molecule_weight > 0
        and primary_active_site_count >= 2
    ):
        return "TREE_READY_LOCAL"
    return "ABSTAIN_ZERO_LOCAL_EVIDENCE"


def validate_partition_internal_timings(
    receipt: Mapping[str, Any],
    chrom: str,
) -> dict[str, float]:
    """Fail closed on pattern loading and the full partition component-loop timer."""

    timings = receipt.get("timings_seconds")
    if not isinstance(timings, Mapping):
        raise ContractError(f"partition child timings_seconds missing for {chrom}")
    source_to_summary = {
        "load_and_aggregate_primary_patterns": (
            "partition_pattern_load_aggregate_seconds"
        ),
        "ordered_hypergraph_dp": "partition_ordered_hypergraph_dp_seconds",
    }
    result: dict[str, float] = {}
    for source_field, summary_field in source_to_summary.items():
        raw = timings.get(source_field)
        if (
            isinstance(raw, bool)
            or not isinstance(raw, (int, float))
            or not isfinite(float(raw))
            or float(raw) < 0
        ):
            raise ContractError(
                f"partition child timing {source_field} is missing/invalid for "
                f"{chrom}: {raw!r}"
            )
        result[summary_field] = float(raw)
    return result


def validate_partition_output(
    partition_dir: Path,
    chrom: str,
    expected_inventory: Mapping[str, Any],
    *,
    expected_input_paths: Mapping[str, Path] | None = None,
) -> tuple[dict[str, Any], list[dict[str, str]], dict[str, Any]]:
    """Authenticate one partition child and independently rederive its metrics."""

    require_real_directory(partition_dir, f"{chrom} partition directory")
    if chrom not in AUTOSOMES:
        raise ContractError(f"invalid chromosome: {chrom}")
    prefix = f"{DATASET}.{chrom}"
    expected_paths = {
        "legacy_components": partition_dir / f"{prefix}.legacy_components.tsv.gz",
        "blocks": partition_dir / f"{prefix}.blocks.tsv.gz",
        "site_membership": partition_dir / f"{prefix}.site_membership.tsv.gz",
        "cut_constraints": partition_dir / f"{prefix}.cut_constraints.tsv.gz",
    }
    receipt_path = partition_dir / "receipt.json"
    receipt_sha = verify_sha256_sidecar(receipt_path)
    receipt = strict_json_load(receipt_path)
    if not isinstance(receipt, Mapping):
        raise ContractError(f"partition receipt root is not an object: {receipt_path}")
    if receipt.get("schema_name") != PARTITION_SCHEMA:
        raise ContractError(f"partition schema drift for {chrom}")
    if receipt.get("schema_version") != PARTITION_SCHEMA_VERSION:
        raise ContractError(f"partition schema version drift for {chrom}")
    if receipt.get("all_pass") is not True:
        raise ContractError(f"partition receipt is not all_pass for {chrom}")
    _validate_all_true(receipt.get("checks"), f"{chrom} partition checks")
    scope = receipt.get("scope")
    if not isinstance(scope, Mapping):
        raise ContractError(f"partition scope missing for {chrom}")
    if (
        scope.get("dataset") != DATASET
        or scope.get("chrom") != chrom
        or scope.get("site_catalog_sites") != expected_inventory["ssnv_sites"]
    ):
        raise ContractError(f"partition scope mismatch for {chrom}: {scope}")
    parameters = receipt.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ContractError(f"partition parameters missing for {chrom}")
    expected_parameter_values = {
        "legacy_gap_bp": 50_000,
        "max_block_size": 8,
        "primary_hp_families": ["1", "2"],
        "require_exact_known_phase_set": True,
        "fixed_call_codes": ["A", "R"],
        "primary_weighting": "raw_molecule",
        "sensitivity_weightings": ["equal_pattern", "log1p_molecule"],
    }
    for key, expected in expected_parameter_values.items():
        if parameters.get(key) != expected:
            raise ContractError(
                f"partition parameter mismatch for {chrom} {key}: "
                f"expected={expected!r} observed={parameters.get(key)!r}"
            )
    outputs = _verify_output_map(
        receipt.get("outputs"), expected_paths, f"{chrom} partition"
    )
    internal_timings = validate_partition_internal_timings(receipt, chrom)

    inputs = receipt.get("inputs")
    if not isinstance(inputs, Mapping) or set(inputs) != {"site_catalog", "molecule_calls"}:
        raise ContractError(f"partition input identities missing for {chrom}")
    if expected_input_paths is not None and set(expected_input_paths) != set(inputs):
        raise ContractError(f"expected partition input paths incomplete for {chrom}")
    verified_inputs: dict[str, dict[str, Any]] = {}
    for name, raw in inputs.items():
        path = expected_input_paths[name] if expected_input_paths else Path(str(raw["path"]))
        verified_inputs[name] = verify_identity(raw, path, f"{chrom} partition input {name}")

    component_rows = list(_read_tsv_rows(expected_paths["legacy_components"], COMPONENT_FIELDS))
    components: dict[str, dict[str, str]] = {}
    for row in component_rows:
        component_id = row["legacy_component_id"]
        if row["dataset"] != DATASET or row["chrom"] != chrom:
            raise ContractError(f"component scope mismatch for {component_id}")
        if component_id in components:
            raise ContractError(f"duplicate component ID for {chrom}: {component_id}")
        components[component_id] = row

    membership_by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen_site_indices: set[int] = set()
    for row in _read_tsv_rows(expected_paths["site_membership"], MEMBERSHIP_FIELDS):
        component_id = row["legacy_component_id"]
        if (
            row["dataset"] != DATASET
            or row["chrom"] != chrom
            or component_id not in components
        ):
            raise ContractError(f"membership scope/component mismatch for {chrom}")
        site_index = _int(row["site_index"], f"{chrom} site_index", minimum=0)
        if site_index in seen_site_indices:
            raise ContractError(f"duplicate target site index for {chrom}: {site_index}")
        seen_site_indices.add(site_index)
        if row["primary_linkage_observed"] not in {"true", "false"}:
            raise ContractError(f"invalid primary linkage flag for {chrom}")
        if row["old_densest8_disposition"] not in {"selected", "cap_excluded"}:
            raise ContractError(f"invalid old-cap disposition for {chrom}")
        membership_by_component[component_id].append(row)

    component_positions: dict[str, tuple[int, ...]] = {}
    selected_positions: dict[str, set[int]] = {}
    active_positions: dict[str, set[int]] = {}
    membership_block_by_position: dict[str, dict[int, str]] = {}
    for component_id, component in components.items():
        rows = membership_by_component.get(component_id, [])
        rows.sort(key=lambda item: _int(item["component_local_index"], "local index"))
        indices = tuple(_int(row["component_local_index"], "local index") for row in rows)
        if indices != tuple(range(len(rows))):
            raise ContractError(f"noncontiguous component-local indices: {component_id}")
        positions = tuple(_int(row["pos1"], "membership position", minimum=1) for row in rows)
        if any(right <= left for left, right in zip(positions, positions[1:])):
            raise ContractError(f"component positions are not strictly ordered: {component_id}")
        component_positions[component_id] = positions
        selected_positions[component_id] = {
            _int(row["pos1"], "membership position")
            for row in rows
            if row["old_densest8_disposition"] == "selected"
        }
        active_positions[component_id] = {
            _int(row["pos1"], "membership position")
            for row in rows
            if row["primary_linkage_observed"] == "true"
        }
        membership_block_by_position[component_id] = {
            _int(row["pos1"], "membership position"): row["block_id"] for row in rows
        }
        if len(membership_block_by_position[component_id]) != len(rows):
            raise ContractError(f"duplicate component position in membership: {component_id}")

    block_rows_by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
    block_agg_by_component: dict[str, Counter[str]] = defaultdict(Counter)
    block_k_distribution: Counter[int] = Counter()
    block_span_bp_values: list[int] = []
    block_all_rows: list[dict[str, str]] = []
    zero_evidence_blocks = 0
    zero_evidence_block_sites = 0
    tree_ready_blocks = 0
    tree_ready_block_sites = 0
    abstain_blocks = 0
    abstain_block_sites = 0
    total_block_sites = 0
    for row in _read_tsv_rows(expected_paths["blocks"], BLOCK_FIELDS):
        component_id = row["legacy_component_id"]
        if (
            row["dataset"] != DATASET
            or row["chrom"] != chrom
            or component_id not in components
        ):
            raise ContractError(f"block scope/component mismatch for {chrom}")
        positions = _csv_ints(row["positions"], f"{component_id} block positions")
        k = _int(row["k"], f"{component_id} block k", minimum=1)
        if len(positions) != k or k > 8:
            raise ContractError(f"invalid block k/positions for {row['block_id']}")
        if tuple(sorted(positions)) != positions:
            raise ContractError(f"unordered block positions for {row['block_id']}")
        if row["start1"] != str(positions[0]) or row["end1"] != str(positions[-1]):
            raise ContractError(f"block endpoints disagree for {row['block_id']}")
        for position in positions:
            if membership_block_by_position[component_id].get(position) != row["block_id"]:
                raise ContractError(f"block/membership disagreement for {row['block_id']}")
        if row["component_status"] != components[component_id]["status"]:
            raise ContractError(f"block status disagreement for {row['block_id']}")
        retained_weight = _int(
            row["raw_retained_molecule_weight"],
            f"{row['block_id']} retained molecule weight",
            minimum=0,
        )
        retained_patterns = _int(
            row["retained_exact_pattern_count"],
            f"{row['block_id']} retained patterns",
            minimum=0,
        )
        active_count = _int(
            row["primary_active_site_count"],
            f"{row['block_id']} active sites",
            minimum=0,
        )
        if active_count > k:
            raise ContractError(f"block active-site count exceeds k: {row['block_id']}")
        if (retained_weight == 0) != (retained_patterns == 0):
            raise ContractError(
                f"zero-evidence block weight/pattern mismatch: {row['block_id']}"
            )
        span_bp = positions[-1] - positions[0]
        block_span_bp_values.append(span_bp)
        aggregate = block_agg_by_component[component_id]
        aggregate["retained_weight"] += retained_weight
        aggregate["retained_patterns"] += retained_patterns
        aggregate["active_sites"] += active_count
        if retained_weight == 0:
            zero_evidence_blocks += 1
            zero_evidence_block_sites += k
        evidence_gate = block_evidence_gate(
            retained_patterns,
            retained_weight,
            active_count,
        )
        tree_ready = evidence_gate == "TREE_READY_LOCAL"
        if tree_ready:
            tree_ready_blocks += 1
            tree_ready_block_sites += k
        else:
            abstain_blocks += 1
            abstain_block_sites += k
        block_all_rows.append({**row, "evidence_gate": evidence_gate})
        block_rows_by_component[component_id].append(row)
        block_k_distribution[k] += 1
        total_block_sites += k

    block_index_by_position: dict[str, dict[int, int]] = {}
    for component_id, rows in block_rows_by_component.items():
        rows.sort(key=lambda item: _int(item["block_index"], "block index"))
        indices = tuple(_int(row["block_index"], "block index") for row in rows)
        if indices != tuple(range(1, len(rows) + 1)):
            raise ContractError(f"noncontiguous block indices: {component_id}")
        position_map: dict[int, int] = {}
        flattened: list[int] = []
        for row in rows:
            block_index = _int(row["block_index"], "block index")
            positions = _csv_ints(row["positions"], "block positions")
            for position in positions:
                if position in position_map:
                    raise ContractError(f"position assigned to multiple blocks: {component_id}")
                position_map[position] = block_index
                flattened.append(position)
        if tuple(flattened) != component_positions[component_id]:
            raise ContractError(f"blocks do not partition component positions: {component_id}")
        block_index_by_position[component_id] = position_map
    if set(block_index_by_position) != set(components):
        missing = sorted(set(components) - set(block_index_by_position))
        raise ContractError(f"components without block partitions for {chrom}: {missing}")

    constraint_agg: dict[str, Counter[str]] = {
        component_id: Counter() for component_id in components
    }
    seen_constraints: set[str] = set()
    for row in _read_tsv_rows(expected_paths["cut_constraints"], CONSTRAINT_FIELDS):
        component_id = row["legacy_component_id"]
        if (
            row["dataset"] != DATASET
            or row["chrom"] != chrom
            or component_id not in components
        ):
            raise ContractError(f"constraint scope/component mismatch for {chrom}")
        constraint_id = row["constraint_id"]
        if constraint_id in seen_constraints:
            raise ContractError(f"duplicate constraint ID for {chrom}: {constraint_id}")
        seen_constraints.add(constraint_id)
        if row["hp_family"] not in {"1", "2"} or not row["phase_set"]:
            raise ContractError(f"invalid HP/PS unit for constraint {constraint_id}")
        if row["unit_id"] != f"HP{row['hp_family']}|PS{row['phase_set']}":
            raise ContractError(f"unit ID mismatch for constraint {constraint_id}")
        positions = _csv_ints(row["positions"], f"{constraint_id} positions")
        call_codes = row["call_codes"]
        n_fixed = _int(row["n_fixed_ra"], f"{constraint_id} n_fixed_ra", minimum=1)
        if len(positions) != n_fixed or len(call_codes) != n_fixed:
            raise ContractError(f"fixed-call cardinality mismatch for {constraint_id}")
        if any(code not in {"R", "A"} for code in call_codes):
            raise ContractError(f"non-fixed call code for {constraint_id}")
        component_position_tuple = component_positions[component_id]
        local_index = {position: index for index, position in enumerate(component_position_tuple)}
        if any(position not in local_index for position in positions):
            raise ContractError(f"constraint position outside component: {constraint_id}")
        expected_span = local_index[positions[-1]] - local_index[positions[0]] + 1
        if _int(row["span_sites"], f"{constraint_id} span_sites") != expected_span:
            raise ContractError(f"constraint span mismatch for {constraint_id}")
        weight = _int(row["molecule_weight"], f"{constraint_id} weight", minimum=1)
        disposition = row["disposition"]
        if disposition not in {RETAINED, CUT, UNAVOIDABLE}:
            raise ContractError(f"unknown constraint disposition: {disposition}")
        involved_blocks = {block_index_by_position[component_id][pos] for pos in positions}
        crossed_cut_count = _int(
            row["crossed_cut_count"], f"{constraint_id} crossed_cut_count", minimum=0
        )
        if crossed_cut_count != max(involved_blocks) - min(involved_blocks):
            raise ContractError(f"crossed-cut count mismatch for {constraint_id}")
        if disposition == RETAINED:
            if len(involved_blocks) != 1:
                raise ContractError(f"retained constraint crosses a cut: {constraint_id}")
            if row["retained_block_index"] != str(next(iter(involved_blocks))):
                raise ContractError(f"retained block mismatch for {constraint_id}")
        else:
            if len(involved_blocks) <= 1 or row["retained_block_index"] != "":
                raise ContractError(f"lost constraint disposition mismatch: {constraint_id}")
        if disposition == UNAVOIDABLE and expected_span <= 8:
            raise ContractError(f"avoidable constraint labeled unavoidable: {constraint_id}")
        if disposition == CUT and expected_span > 8:
            raise ContractError(f"long-span constraint not labeled unavoidable: {constraint_id}")

        aggregate = constraint_agg[component_id]
        aggregate["patterns"] += 1
        aggregate["total_weight"] += weight
        if disposition == RETAINED:
            aggregate["retained_patterns"] += 1
            aggregate["retained_weight"] += weight
        else:
            aggregate["lost_patterns"] += 1
            aggregate["lost_weight"] += weight
        if disposition == UNAVOIDABLE:
            aggregate["unavoidable_patterns"] += 1
            aggregate["unavoidable_weight"] += weight
            if n_fixed > 8:
                aggregate["unavoidable_n_fixed_ra_gt8_patterns"] += 1
                aggregate["unavoidable_n_fixed_ra_gt8_weight"] += weight
            else:
                aggregate[
                    "unavoidable_n_fixed_ra_lte8_span_gt8_patterns"
                ] += 1
                aggregate[
                    "unavoidable_n_fixed_ra_lte8_span_gt8_weight"
                ] += weight
        if set(positions).issubset(selected_positions[component_id]):
            aggregate["old_retained_patterns"] += 1
            aggregate["old_retained_weight"] += weight

    status_counts: Counter[str] = Counter()
    totals: Counter[str] = Counter()
    weight_stable_components = 0
    component_retention_delta_values: list[Decimal] = []
    for component_id, row in components.items():
        positions = component_positions.get(component_id, ())
        k = _int(row["pre_cap_k"], f"{component_id} pre_cap_k", minimum=9)
        if len(positions) != k:
            raise ContractError(f"component k disagrees with membership: {component_id}")
        old_selected = _int(row["old_densest8_selected"], "old selected")
        old_excluded = _int(row["old_cap_excluded"], "old excluded")
        if old_selected != 8 or len(selected_positions[component_id]) != 8:
            raise ContractError(f"old densest-8 selection is not exactly 8: {component_id}")
        expected_dense_start = min(
            range(k - 8 + 1),
            key=lambda start: positions[start + 7] - positions[start],
        )
        expected_dense_positions = set(
            positions[expected_dense_start : expected_dense_start + 8]
        )
        if selected_positions[component_id] != expected_dense_positions:
            raise ContractError(f"old densest-8 window drift: {component_id}")
        if old_selected + old_excluded != k:
            raise ContractError(f"old selected/excluded conservation fails: {component_id}")
        blocks = block_rows_by_component.get(component_id, [])
        if _int(row["new_block_count"], "new block count") != len(blocks):
            raise ContractError(f"new block count mismatch: {component_id}")
        if _int(row["new_site_retained"], "new retained sites") != k:
            raise ContractError(f"new site conservation mismatch: {component_id}")
        active_count = len(active_positions[component_id])
        if _int(row["primary_active_site_count"], "active site count") != active_count:
            raise ContractError(f"active site count mismatch: {component_id}")
        if row["primary_active_site_fraction"] != f"{active_count / k:.12f}":
            raise ContractError(f"active site fraction mismatch: {component_id}")
        if row["start1"] != str(positions[0]) or row["end1"] != str(positions[-1]):
            raise ContractError(f"component endpoints mismatch: {component_id}")
        if _int(row["span_bp"], "component span") != positions[-1] - positions[0]:
            raise ContractError(f"component span mismatch: {component_id}")
        if row["positions_sha256"] != semantic_sha256(positions):
            raise ContractError(f"component position digest mismatch: {component_id}")

        aggregate = constraint_agg[component_id]
        component_expectations = {
            "exact_pattern_count": aggregate["patterns"],
            "raw_total_molecule_weight": aggregate["total_weight"],
            "raw_retained_molecule_weight": aggregate["retained_weight"],
            "raw_lost_molecule_weight": aggregate["lost_weight"],
            "old_densest8_retained_molecule_weight": aggregate["old_retained_weight"],
            "retained_exact_pattern_count": aggregate["retained_patterns"],
            "lost_exact_pattern_count": aggregate["lost_patterns"],
            "unavoidable_pattern_count": aggregate["unavoidable_patterns"],
            "old_densest8_retained_pattern_count": aggregate["old_retained_patterns"],
        }
        for field, expected in component_expectations.items():
            if _int(row[field], f"{component_id} {field}", minimum=0) != expected:
                raise ContractError(
                    f"component constraint aggregate mismatch: {component_id} {field}"
                )
        block_aggregate = block_agg_by_component[component_id]
        if (
            block_aggregate["retained_weight"] != aggregate["retained_weight"]
            or block_aggregate["retained_patterns"] != aggregate["retained_patterns"]
            or block_aggregate["active_sites"] != active_count
        ):
            raise ContractError(f"block/component evidence aggregate mismatch: {component_id}")
        if row["raw_retention_ratio"] != _ratio(
            Decimal(aggregate["retained_weight"]), Decimal(aggregate["total_weight"])
        ):
            raise ContractError(f"new retention ratio mismatch: {component_id}")
        if row["old_densest8_retention_ratio"] != _ratio(
            Decimal(aggregate["old_retained_weight"]), Decimal(aggregate["total_weight"])
        ):
            raise ContractError(f"old retention ratio mismatch: {component_id}")
        if aggregate["total_weight"] > 0:
            component_retention_delta_values.append(
                Decimal(aggregate["retained_weight"])
                / Decimal(aggregate["total_weight"])
                - Decimal(aggregate["old_retained_weight"])
                / Decimal(aggregate["total_weight"])
            )

        raw_cuts = _csv_ints(row["raw_cut_indices"], f"{component_id} raw cuts")
        expected_cuts: list[int] = []
        cumulative = 0
        for block in blocks[:-1]:
            cumulative += _int(block["k"], "block k")
            expected_cuts.append(cumulative)
        if raw_cuts != tuple(expected_cuts):
            raise ContractError(f"raw cut indices disagree with blocks: {component_id}")
        expected_gap_sum = sum(
            positions[index] - positions[index - 1] for index in raw_cuts
        )
        if _int(row["cut_gap_sum_bp"], "cut gap sum") != expected_gap_sum:
            raise ContractError(f"cut-gap sum mismatch: {component_id}")
        weight_stable = row["weight_stable"] == "true"
        if row["weight_stable"] not in {"true", "false"}:
            raise ContractError(f"invalid weight_stable flag: {component_id}")
        sensitivity_equal = (
            row["raw_cut_indices"]
            == row["equal_pattern_cut_indices"]
            == row["log1p_cut_indices"]
        )
        if weight_stable != sensitivity_equal:
            raise ContractError(f"weight-stability flag mismatch: {component_id}")
        expected_status = _expected_status(
            pattern_count=aggregate["patterns"],
            lost_pattern_count=aggregate["lost_patterns"],
            active_site_count=active_count,
            component_k=k,
            weight_stable=weight_stable,
        )
        if row["status"] != expected_status:
            raise ContractError(f"component status mismatch: {component_id}")
        status_counts[row["status"]] += 1
        weight_stable_components += int(weight_stable)

        totals["components"] += 1
        totals["sites"] += k
        totals["old_selected_sites"] += old_selected
        totals["old_excluded_sites"] += old_excluded
        totals["new_blocks"] += len(blocks)
        totals["new_retained_sites"] += k
        totals["primary_active_sites"] += active_count
        totals["exact_patterns"] += aggregate["patterns"]
        totals["raw_total_molecule_weight"] += aggregate["total_weight"]
        totals["new_retained_molecule_weight"] += aggregate["retained_weight"]
        totals["new_lost_molecule_weight"] += aggregate["lost_weight"]
        totals["old_densest8_retained_molecule_weight"] += aggregate[
            "old_retained_weight"
        ]
        totals["retained_exact_patterns"] += aggregate["retained_patterns"]
        totals["lost_exact_patterns"] += aggregate["lost_patterns"]
        totals["unavoidable_patterns"] += aggregate["unavoidable_patterns"]
        totals["unavoidable_molecule_weight"] += aggregate["unavoidable_weight"]
        totals["unavoidable_n_fixed_ra_gt8_patterns"] += aggregate[
            "unavoidable_n_fixed_ra_gt8_patterns"
        ]
        totals["unavoidable_n_fixed_ra_gt8_molecule_weight"] += aggregate[
            "unavoidable_n_fixed_ra_gt8_weight"
        ]
        totals[
            "unavoidable_n_fixed_ra_lte8_span_gt8_patterns"
        ] += aggregate["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"]
        totals[
            "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight"
        ] += aggregate["unavoidable_n_fixed_ra_lte8_span_gt8_weight"]
        totals["old_densest8_retained_patterns"] += aggregate[
            "old_retained_patterns"
        ]

    totals["zero_evidence_blocks"] = zero_evidence_blocks
    totals["zero_evidence_block_sites"] = zero_evidence_block_sites
    totals["tree_ready_blocks"] = tree_ready_blocks
    totals["tree_ready_block_sites"] = tree_ready_block_sites
    totals["abstain_blocks"] = abstain_blocks
    totals["abstain_block_sites"] = abstain_block_sites

    if total_block_sites != totals["sites"] or len(seen_site_indices) != totals["sites"]:
        raise ContractError(f"partition site conservation fails for {chrom}")
    if (
        tree_ready_blocks + abstain_blocks != totals["new_blocks"]
        or tree_ready_block_sites + abstain_block_sites != totals["sites"]
    ):
        raise ContractError(f"tree-ready/abstain block conservation fails for {chrom}")
    if (
        totals["unavoidable_n_fixed_ra_gt8_patterns"]
        + totals["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"]
        != totals["unavoidable_patterns"]
        or totals["unavoidable_n_fixed_ra_gt8_molecule_weight"]
        + totals["unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight"]
        != totals["unavoidable_molecule_weight"]
    ):
        raise ContractError(f"unavoidable mechanism conservation fails for {chrom}")
    counts = receipt.get("counts")
    if not isinstance(counts, Mapping):
        raise ContractError(f"partition receipt counts missing for {chrom}")
    receipt_count_expectations = {
        "target_components": totals["components"],
        "target_sites": totals["sites"],
        "old_selected_sites": totals["old_selected_sites"],
        "old_cap_excluded_sites": totals["old_excluded_sites"],
        "new_blocks": totals["new_blocks"],
        "primary_active_sites_component_sum": totals["primary_active_sites"],
        "exact_patterns": totals["exact_patterns"],
        "raw_total_molecule_weight": totals["raw_total_molecule_weight"],
        "raw_retained_molecule_weight": totals["new_retained_molecule_weight"],
        "raw_lost_molecule_weight": totals["new_lost_molecule_weight"],
        "unavoidable_patterns": totals["unavoidable_patterns"],
        "weight_stable_components": weight_stable_components,
    }
    for key, expected in receipt_count_expectations.items():
        if counts.get(key) != expected:
            raise ContractError(
                f"partition receipt count mismatch for {chrom} {key}: "
                f"expected={expected} observed={counts.get(key)}"
            )
    expected_components = expected_inventory["k_gt8_components"]
    expected_sites = expected_inventory["k_gt8_sites"]
    if totals["components"] != expected_components or totals["sites"] != expected_sites:
        raise ContractError(
            f"canonical target mismatch for {chrom}: "
            f"components={totals['components']}/{expected_components} "
            f"sites={totals['sites']}/{expected_sites}"
        )
    if receipt.get("status_counts") != dict(sorted(status_counts.items())):
        raise ContractError(f"partition status counts mismatch for {chrom}")
    if (
        totals["raw_total_molecule_weight"]
        != totals["new_retained_molecule_weight"] + totals["new_lost_molecule_weight"]
    ):
        raise ContractError(f"molecule-weight conservation fails for {chrom}")

    metrics = {
        "chrom": chrom,
        "ssnv_sites": expected_inventory["ssnv_sites"],
        **dict(totals),
        **internal_timings,
        "new_retention_ratio": _ratio(
            Decimal(totals["new_retained_molecule_weight"]),
            Decimal(totals["raw_total_molecule_weight"]),
        ),
        "old_densest8_retention_ratio": _ratio(
            Decimal(totals["old_densest8_retained_molecule_weight"]),
            Decimal(totals["raw_total_molecule_weight"]),
        ),
        "weight_stable_components": weight_stable_components,
        "weight_sensitive_components": totals["components"] - weight_stable_components,
        "status_counts": dict(sorted(status_counts.items())),
        "block_k_distribution": {
            str(k): count for k, count in sorted(block_k_distribution.items())
        },
        "block_span_bp_distribution": _distribution(block_span_bp_values),
        "component_retention_delta_distribution": _distribution(
            component_retention_delta_values,
            total_count=totals["components"],
        ),
        "_block_span_bp_values": block_span_bp_values,
        "_component_retention_delta_values": component_retention_delta_values,
    }
    evidence = {
        "partition_receipt": {
            "path": str(receipt_path.resolve()),
            "size_bytes": receipt_path.stat().st_size,
            "sha256": receipt_sha,
        },
        "inputs": verified_inputs,
        "outputs": outputs,
        "semantic_result_sha256": receipt.get("semantic_result_sha256"),
        "pattern_diagnostics": receipt.get("pattern_diagnostics"),
        "primary_rows_by_hp_family": receipt.get("primary_rows_by_hp_family"),
        "timings_seconds": dict(receipt["timings_seconds"]),
        "_block_all_rows": block_all_rows,
    }
    return metrics, component_rows, evidence


def _verify_stage_receipt(
    stage_path: Path,
    *,
    chrom: str,
    stage: str,
    contract_sha: str,
    expected_identity: Any,
    expected_tool: Any,
) -> tuple[dict[str, Any], dict[str, Any]]:
    identity = verify_identity(
        expected_identity,
        stage_path,
        f"{chrom} {stage} stage receipt",
        require_sidecar=True,
    )
    receipt = strict_json_load(stage_path)
    expected_scalars = {
        "schema_name": f"{RUN_SCHEMA}.stage_receipt",
        "schema_version": RUN_SCHEMA_VERSION,
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "stage": stage,
        "status": "COMPLETED",
        "exit_code": 0,
        "contract_sha256": contract_sha,
    }
    for key, expected in expected_scalars.items():
        if receipt.get(key) != expected:
            raise ContractError(f"{chrom} {stage} stage receipt mismatch for {key}")
    if receipt.get("tool") != expected_tool:
        raise ContractError(f"{chrom} {stage} tool identity drift")
    command = receipt.get("command")
    if not isinstance(command, list) or not command:
        raise ContractError(f"{chrom} {stage} stage command missing")
    expected_resource = stage_path.parent / f"{stage}.time_v.txt"
    expected_logs = {
        "stdout": stage_path.parent / f"{stage}.stdout.log",
        "stderr": stage_path.parent / f"{stage}.stderr.log",
        "resource": expected_resource,
    }
    _verify_output_map(receipt.get("logs"), expected_logs, f"{chrom} {stage} logs")
    expected_timed = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(expected_resource),
        "--",
        *command,
    ]
    if receipt.get("timed_command") != expected_timed:
        raise ContractError(f"{chrom} {stage} timed command drift")
    wall_seconds = receipt.get("wall_seconds")
    if not isinstance(wall_seconds, (int, float)) or wall_seconds < 0:
        raise ContractError(f"{chrom} {stage} invalid wall_seconds")
    return receipt, parse_gnu_time_v(expected_resource)


def _verify_extract_child(
    extract_dir: Path,
    *,
    chrom: str,
    inventory: Mapping[str, Any],
    stage_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    receipt_path = extract_dir / "receipt.json"
    child_identity = stage_receipt.get("child_receipt")
    identity = verify_identity(
        child_identity,
        receipt_path,
        f"{chrom} extraction child receipt",
        require_sidecar=True,
    )
    receipt = strict_json_load(receipt_path)
    if (
        receipt.get("schema_name")
        != "intersubmod.lossless_read_linkage_chromosome_receipt"
        or receipt.get("schema_version") != "1.3.0"
        or receipt.get("all_pass") is not True
    ):
        raise ContractError(f"invalid extraction child receipt for {chrom}")
    expected_scope = {
        "dataset": DATASET,
        "chrom": chrom,
        "n_sSNV": inventory["ssnv_sites"],
    }
    if receipt.get("scope") != expected_scope:
        raise ContractError(f"extraction child scope mismatch for {chrom}")
    _validate_all_true(receipt.get("checks"), f"{chrom} extraction checks")
    prefix = f"{DATASET}.{chrom}"
    names = (
        f"{prefix}.targets.bed",
        f"{prefix}.molecule_sparse_calls.tsv.gz",
        f"{prefix}.site_catalog.tsv.gz",
        f"{prefix}.cut_support.tsv.gz",
        f"{prefix}.components.tsv.gz",
        f"{prefix}.site_component_membership.tsv.gz",
        "samtools_view.stderr.log",
    )
    outputs = _verify_output_map(
        receipt.get("outputs"),
        {name: extract_dir / name for name in names},
        f"{chrom} extraction",
    )
    counts = receipt.get("counts")
    if not isinstance(counts, Mapping):
        raise ContractError(f"extraction counts missing for {chrom}")
    required_count_names = {
        "raw_overlapping_alignments",
        "canonical_eligible_alignments_pre_identity_collapse",
        "duplicate_alignment_identities_collapsed",
        "duplicate_alignment_identities_seen",
        "canonical_eligible_alignments",
        "unique_alignment_ids",
        "unique_molecule_ids",
        "molecule_sparse_rows_written",
        "site_call_rows_sparse",
        "known_ps_hp12_molecule_rows",
        "known_ps_hp1_molecule_rows",
        "known_ps_hp2_molecule_rows",
        "sidecar_missing",
    }
    missing = sorted(required_count_names - set(counts))
    if missing:
        raise ContractError(f"extraction counts missing for {chrom}: {missing}")
    numeric_counts = {
        key: _int(str(counts[key]), f"{chrom} extraction {key}", minimum=0)
        for key in required_count_names
    }
    conflicts = _int(
        str(counts.get("duplicate_alignment_identity_conflicts", 0)),
        f"{chrom} duplicate identity conflicts",
        minimum=0,
    )
    pre_collapse = numeric_counts[
        "canonical_eligible_alignments_pre_identity_collapse"
    ]
    collapsed = numeric_counts["duplicate_alignment_identities_collapsed"]
    canonical = numeric_counts["canonical_eligible_alignments"]
    if pre_collapse != canonical + collapsed:
        raise ContractError(f"duplicate-collapse conservation fails for {chrom}")
    if numeric_counts["duplicate_alignment_identities_seen"] != collapsed:
        raise ContractError(f"duplicate seen/collapsed count mismatch for {chrom}")
    for key in (
        "unique_alignment_ids",
        "unique_molecule_ids",
        "molecule_sparse_rows_written",
    ):
        if numeric_counts[key] != canonical:
            raise ContractError(f"canonical molecule conservation fails for {chrom}: {key}")
    known_hp1 = numeric_counts["known_ps_hp1_molecule_rows"]
    known_hp2 = numeric_counts["known_ps_hp2_molecule_rows"]
    if numeric_counts["known_ps_hp12_molecule_rows"] != known_hp1 + known_hp2:
        raise ContractError(f"known-PS HP1/HP2 row conservation fails for {chrom}")
    if conflicts != 0 or numeric_counts["sidecar_missing"] != 0:
        raise ContractError(f"failed extraction identity joins for {chrom}")
    checks = receipt.get("checks")
    failed_checks = (
        sum(value is not True for value in checks.values())
        if isinstance(checks, Mapping)
        else 1
    )
    if failed_checks != 0:
        raise ContractError(f"extraction failed-check count is nonzero for {chrom}")
    diagnostics = {
        "raw_overlapping_alignments": numeric_counts["raw_overlapping_alignments"],
        "eligible_alignments_pre_identity_collapse": pre_collapse,
        "duplicate_alignment_identities_collapsed": collapsed,
        "canonical_unique_molecules": canonical,
        "molecule_sparse_rows": numeric_counts["molecule_sparse_rows_written"],
        "sparse_site_calls": numeric_counts["site_call_rows_sparse"],
        "known_ps_hp12_molecule_rows": numeric_counts[
            "known_ps_hp12_molecule_rows"
        ],
        "known_ps_hp1_molecule_rows": known_hp1,
        "known_ps_hp2_molecule_rows": known_hp2,
        "extraction_failed_checks": failed_checks,
        "duplicate_identity_conflicts": conflicts,
    }
    return {"identity": identity, "outputs": outputs, "diagnostics": diagnostics}


def _parse_runner_summary(path: Path) -> list[dict[str, str]]:
    require_regular_file(path, "runner chromosome summary")
    try:
        with path.open("r", encoding="utf-8", errors="strict", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != RUNNER_SUMMARY_FIELDS:
                raise ContractError("runner chromosome summary header drift")
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ContractError(f"cannot read runner chromosome summary: {exc}") from exc
    if tuple(row["chrom"] for row in rows) != AUTOSOMES:
        raise ContractError("runner chromosome summary is not exact chr1-22 order")
    return rows


def validate_full_run(
    full_root: Path,
    outer_time_path: Path | None = None,
) -> tuple[
    dict[str, Any],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    """Validate the production runner tree and return summary payload + components."""

    require_real_directory(full_root, "full-run root")
    contract_path = full_root / "run_contract.json"
    contract_sha = verify_sha256_sidecar(contract_path)
    contract = strict_json_load(contract_path)
    if not isinstance(contract, Mapping):
        raise ContractError("run contract root is not an object")
    if (
        contract.get("schema_name") != f"{RUN_SCHEMA}.run_contract"
        or contract.get("schema_version") != RUN_SCHEMA_VERSION
        or contract.get("sample") != DATASET
    ):
        raise ContractError("run contract schema/sample drift")
    expected_contract_scope = {
        "chromosomes": list(AUTOSOMES),
        "comprehensive_chr1_22": True,
        "test_mode": False,
    }
    if contract.get("scope") != expected_contract_scope:
        raise ContractError("run contract is not the comprehensive production scope")
    if contract.get("canonical_expected_inventory_enforced") is not True:
        raise ContractError("canonical inventory was not enforced by the runner")
    parameters = contract.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ContractError("run parameters missing")
    if (
        parameters.get("legacy_gap_bp") != 50_000
        or parameters.get("max_block_size") != 8
        or parameters.get("baseq_min") != 20
    ):
        raise ContractError("runner gap/block/BASEQ parameter drift")
    tools = contract.get("tools")
    if not isinstance(tools, Mapping) or not {
        "runner",
        "extractor",
        "partitioner",
    }.issubset(tools):
        raise ContractError("runner tool identities missing")
    runner_identity = tools["runner"]
    if not isinstance(runner_identity, Mapping) or not isinstance(
        runner_identity.get("path"), str
    ):
        raise ContractError("runner tool identity path is missing")
    runner_path = Path(runner_identity["path"])
    verify_identity(runner_identity, runner_path, "runner tool identity")
    manifest_identity = contract.get("frozen_manifest")
    if not isinstance(manifest_identity, Mapping) or not isinstance(
        manifest_identity.get("path"), str
    ):
        raise ContractError("frozen manifest identity path is missing")
    manifest_path = Path(manifest_identity["path"])
    verify_identity(manifest_identity, manifest_path, "frozen manifest identity")
    inventory_rows = contract.get("vcf_derived_inventory")
    if not isinstance(inventory_rows, list):
        raise ContractError("VCF-derived inventory missing")
    if tuple(row.get("chrom") for row in inventory_rows) != AUTOSOMES:
        raise ContractError("VCF-derived inventory is not exact chr1-22 order")
    inventory_by_chrom: dict[str, Mapping[str, Any]] = {}
    for row in inventory_rows:
        chrom = row["chrom"]
        observed = tuple(
            row.get(key)
            for key in ("ssnv_sites", "k_gt8_components", "k_gt8_sites", "k_gt8_max_k")
        )
        if observed != CANONICAL_EXPECTED[chrom]:
            raise ContractError(
                f"canonical VCF inventory mismatch for {chrom}: "
                f"expected={CANONICAL_EXPECTED[chrom]} observed={observed}"
            )
        inventory_by_chrom[chrom] = row

    receipt_path = full_root / "receipt.json"
    receipt_sha = verify_sha256_sidecar(receipt_path)
    receipt = strict_json_load(receipt_path)
    if not isinstance(receipt, Mapping):
        raise ContractError("runner receipt root is not an object")
    expected_receipt_scalars = {
        "schema_name": f"{RUN_SCHEMA}.run_receipt",
        "schema_version": RUN_SCHEMA_VERSION,
        "all_pass": True,
        "comprehensive_all_pass": True,
        "sample": DATASET,
        "scope": {"chromosomes": list(AUTOSOMES), "test_mode": False},
    }
    for key, expected in expected_receipt_scalars.items():
        if receipt.get(key) != expected:
            raise ContractError(f"runner receipt mismatch for {key}")
    verify_identity(
        receipt.get("contract"),
        contract_path,
        "runner contract identity",
        require_sidecar=True,
    )
    expected_totals = {
        "chromosomes": 22,
        "ssnv_sites": 79_687,
        "k_gt8_components": 408,
        "k_gt8_sites": 47_570,
        "k_gt8_max_k": 3_574,
        "partitioned_chromosomes": 21,
        "zero_target_skipped_chromosomes": 1,
    }
    if receipt.get("totals") != expected_totals:
        raise ContractError(f"runner canonical totals mismatch: {receipt.get('totals')}")

    marker_path = full_root / "_SUCCESS"
    marker = strict_json_load(marker_path)
    expected_marker_values = {
        "schema_name": f"{RUN_SCHEMA}.success_marker",
        "schema_version": RUN_SCHEMA_VERSION,
        "all_pass": True,
        "comprehensive": True,
        "sample": DATASET,
        "scope": list(AUTOSOMES),
        "run_receipt": {
            "path": str(receipt_path.resolve()),
            "sha256": receipt_sha,
        },
    }
    for key, expected in expected_marker_values.items():
        if marker.get(key) != expected:
            raise ContractError(f"_SUCCESS marker mismatch for {key}")

    runner_summary_path = full_root / "chromosome_summary.tsv"
    runner_summary_sha = verify_sha256_sidecar(runner_summary_path)
    runner_summary_identity = verify_identity(
        receipt.get("outputs", {}).get("chromosome_summary"),
        runner_summary_path,
        "runner chromosome summary",
        require_sidecar=True,
    )
    runner_summary_rows = _parse_runner_summary(runner_summary_path)
    runner_summary_by_chrom = {row["chrom"]: row for row in runner_summary_rows}
    chromosome_receipts = receipt.get("chromosomes")
    if (
        not isinstance(chromosome_receipts, list)
        or tuple(row.get("chrom") for row in chromosome_receipts) != AUTOSOMES
    ):
        raise ContractError("runner chromosome receipt list is not exact chr1-22")

    per_chrom: list[dict[str, Any]] = []
    component_all: list[dict[str, str]] = []
    block_all: list[dict[str, str]] = []
    evidence_by_chrom: dict[str, Any] = {}
    extraction_resources: list[dict[str, Any]] = []
    partition_resources: list[dict[str, Any]] = []
    for chrom, chromosome_receipt in zip(AUTOSOMES, chromosome_receipts):
        inventory = inventory_by_chrom[chrom]
        if chromosome_receipt.get("inventory") != inventory:
            raise ContractError(f"runner inventory receipt drift for {chrom}")
        summary_row = runner_summary_by_chrom[chrom]
        expected_inventory_values = CANONICAL_EXPECTED[chrom]
        for field, expected in zip(
            ("ssnv_sites", "k_gt8_components", "k_gt8_sites", "k_gt8_max_k"),
            expected_inventory_values,
        ):
            if _int(summary_row[field], f"{chrom} runner summary {field}") != expected:
                raise ContractError(f"runner summary inventory mismatch for {chrom} {field}")

        chrom_root = full_root / "chromosomes" / chrom
        require_real_directory(chrom_root, f"{chrom} stage root")
        extract_stage_path = chrom_root / "extract.stage_receipt.json"
        extract_stage, extract_resource = _verify_stage_receipt(
            extract_stage_path,
            chrom=chrom,
            stage="extract",
            contract_sha=contract_sha,
            expected_identity=chromosome_receipt.get("extraction"),
            expected_tool=tools["extractor"],
        )
        extract_child = _verify_extract_child(
            chrom_root / "extract",
            chrom=chrom,
            inventory=inventory,
            stage_receipt=extract_stage,
        )
        extraction_resources.append(extract_resource)
        if (
            summary_row["extraction_status"] != "COMPLETED"
            or summary_row["extraction_stage_receipt"] != str(extract_stage_path.resolve())
            or summary_row["extraction_stage_receipt_sha256"]
            != chromosome_receipt["extraction"]["sha256"]
            or float(summary_row["extraction_wall_seconds"])
            != float(extract_stage["wall_seconds"])
        ):
            raise ContractError(f"runner extraction summary drift for {chrom}")

        if chrom == "chr21":
            skip_path = chrom_root / "partition.skip_receipt.json"
            skip_identity = verify_identity(
                chromosome_receipt.get("partition"),
                skip_path,
                "chr21 partition skip receipt",
                require_sidecar=True,
            )
            skip = strict_json_load(skip_path)
            expected_skip = {
                "schema_name": f"{RUN_SCHEMA}.partition_skip_receipt",
                "schema_version": RUN_SCHEMA_VERSION,
                "all_pass": True,
                "sample": DATASET,
                "chrom": "chr21",
                "stage": "partition",
                "status": "SKIP_NO_K_GT8_TARGET",
                "contract_sha256": contract_sha,
                "expected_target_components": 0,
                "expected_target_sites": 0,
                "extraction_receipt": chromosome_receipt["extraction"],
                "inventory": inventory,
            }
            for key, expected in expected_skip.items():
                if skip.get(key) != expected:
                    raise ContractError(f"chr21 skip receipt mismatch for {key}")
            if (
                chromosome_receipt.get("partition_status")
                != "SKIP_NO_K_GT8_TARGET"
                or summary_row["partition_status"] != "SKIP_NO_K_GT8_TARGET"
                or summary_row["partition_wall_seconds"] != ""
                or summary_row["partition_stage_receipt"] != str(skip_path.resolve())
                or summary_row["partition_stage_receipt_sha256"] != skip_identity["sha256"]
            ):
                raise ContractError("chr21 runner partition summary drift")
            per_chrom.append(
                {
                    "chrom": chrom,
                    "partition_stage_status": "SKIP_NO_K_GT8_TARGET",
                    "ssnv_sites": inventory["ssnv_sites"],
                    **extract_child["diagnostics"],
                    "components": 0,
                    "sites": 0,
                    "old_selected_sites": 0,
                    "old_excluded_sites": 0,
                    "new_blocks": 0,
                    "new_retained_sites": 0,
                    "primary_active_sites": 0,
                    "exact_patterns": 0,
                    "raw_total_molecule_weight": 0,
                    "new_retained_molecule_weight": 0,
                    "new_lost_molecule_weight": 0,
                    "new_retention_ratio": "",
                    "old_densest8_retained_molecule_weight": 0,
                    "old_densest8_retention_ratio": "",
                    "retained_exact_patterns": 0,
                    "lost_exact_patterns": 0,
                    "unavoidable_patterns": 0,
                    "unavoidable_molecule_weight": 0,
                    "unavoidable_n_fixed_ra_gt8_patterns": 0,
                    "unavoidable_n_fixed_ra_gt8_molecule_weight": 0,
                    "unavoidable_n_fixed_ra_lte8_span_gt8_patterns": 0,
                    "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight": 0,
                    "old_densest8_retained_patterns": 0,
                    "weight_stable_components": 0,
                    "weight_sensitive_components": 0,
                    "zero_evidence_blocks": 0,
                    "zero_evidence_block_sites": 0,
                    "tree_ready_blocks": 0,
                    "tree_ready_block_sites": 0,
                    "abstain_blocks": 0,
                    "abstain_block_sites": 0,
                    "status_counts": {},
                    "block_k_distribution": {},
                    "block_span_bp_distribution": _distribution(()),
                    "component_retention_delta_distribution": _distribution(
                        (), total_count=0
                    ),
                    "_block_span_bp_values": [],
                    "_component_retention_delta_values": [],
                    "extraction_wall_seconds": extract_stage["wall_seconds"],
                    "partition_wall_seconds": "",
                    "partition_pattern_load_aggregate_seconds": 0.0,
                    "partition_ordered_hypergraph_dp_seconds": 0.0,
                    "extraction_max_rss_kb": extract_resource["max_rss_kb"],
                    "partition_max_rss_kb": "",
                }
            )
            evidence_by_chrom[chrom] = {
                "extraction": extract_child,
                "partition_skip_receipt": skip_identity,
            }
            continue

        partition_stage_path = chrom_root / "partition.stage_receipt.json"
        partition_stage, partition_resource = _verify_stage_receipt(
            partition_stage_path,
            chrom=chrom,
            stage="partition",
            contract_sha=contract_sha,
            expected_identity=chromosome_receipt.get("partition"),
            expected_tool=tools["partitioner"],
        )
        partition_dir = chrom_root / "partition"
        partition_receipt_path = partition_dir / "receipt.json"
        verify_identity(
            partition_stage.get("child_receipt"),
            partition_receipt_path,
            f"{chrom} partition child receipt",
            require_sidecar=True,
        )
        metrics, component_rows, partition_evidence = validate_partition_output(
            partition_dir,
            chrom,
            inventory,
            expected_input_paths={
                "site_catalog": chrom_root
                / "extract"
                / f"{DATASET}.{chrom}.site_catalog.tsv.gz",
                "molecule_calls": chrom_root
                / "extract"
                / f"{DATASET}.{chrom}.molecule_sparse_calls.tsv.gz",
            },
        )
        block_all.extend(partition_evidence.pop("_block_all_rows"))
        pattern_diagnostics = partition_evidence["pattern_diagnostics"]
        primary_by_family = partition_evidence["primary_rows_by_hp_family"]
        if not isinstance(pattern_diagnostics, Mapping) or not isinstance(
            primary_by_family, Mapping
        ):
            raise ContractError(f"partition pattern diagnostics missing for {chrom}")
        extraction_diagnostics = extract_child["diagnostics"]
        expected_partition_extraction_counts = {
            "molecule_rows_seen": extraction_diagnostics["molecule_sparse_rows"],
            "primary_known_ps_rows": extraction_diagnostics[
                "known_ps_hp12_molecule_rows"
            ],
        }
        for key, expected in expected_partition_extraction_counts.items():
            if pattern_diagnostics.get(key) != expected:
                raise ContractError(
                    f"partition/extraction diagnostic mismatch for {chrom} {key}"
                )
        if primary_by_family != {
            "1": extraction_diagnostics["known_ps_hp1_molecule_rows"],
            "2": extraction_diagnostics["known_ps_hp2_molecule_rows"],
        }:
            raise ContractError(
                f"partition/extraction known-PS family counts mismatch for {chrom}"
            )
        internal_partition_seconds = (
            metrics["partition_pattern_load_aggregate_seconds"]
            + metrics["partition_ordered_hypergraph_dp_seconds"]
        )
        if internal_partition_seconds > float(partition_stage["wall_seconds"]) + 0.1:
            raise ContractError(
                f"partition child internal timing exceeds stage wall for {chrom}: "
                f"internal={internal_partition_seconds} "
                f"stage={partition_stage['wall_seconds']}"
            )
        partition_resources.append(partition_resource)
        if (
            chromosome_receipt.get("partition_status") != "COMPLETED"
            or summary_row["partition_status"] != "COMPLETED"
            or summary_row["partition_stage_receipt"]
            != str(partition_stage_path.resolve())
            or summary_row["partition_stage_receipt_sha256"]
            != chromosome_receipt["partition"]["sha256"]
            or float(summary_row["partition_wall_seconds"])
            != float(partition_stage["wall_seconds"])
        ):
            raise ContractError(f"runner partition summary drift for {chrom}")
        metrics.update(
            {
                "partition_stage_status": "COMPLETED",
                **extract_child["diagnostics"],
                "extraction_wall_seconds": extract_stage["wall_seconds"],
                "partition_wall_seconds": partition_stage["wall_seconds"],
                "extraction_max_rss_kb": extract_resource["max_rss_kb"],
                "partition_max_rss_kb": partition_resource["max_rss_kb"],
            }
        )
        per_chrom.append(metrics)
        component_all.extend(component_rows)
        evidence_by_chrom[chrom] = {
            "extraction": extract_child,
            "partition": partition_evidence,
        }

    if len(per_chrom) != 22 or len(component_all) != 408:
        raise ContractError(
            "validated cardinality is not 22 chromosomes/408 components"
        )

    aggregate: Counter[str] = Counter()
    aggregate_status: Counter[str] = Counter()
    aggregate_block_k: Counter[int] = Counter()
    aggregate_block_spans: list[int] = []
    aggregate_component_deltas: list[Decimal] = []
    for row in per_chrom:
        for field in (
            "components",
            "sites",
            "old_selected_sites",
            "old_excluded_sites",
            "new_blocks",
            "new_retained_sites",
            "primary_active_sites",
            "exact_patterns",
            "raw_total_molecule_weight",
            "new_retained_molecule_weight",
            "new_lost_molecule_weight",
            "old_densest8_retained_molecule_weight",
            "retained_exact_patterns",
            "lost_exact_patterns",
            "unavoidable_patterns",
            "unavoidable_molecule_weight",
            "unavoidable_n_fixed_ra_gt8_patterns",
            "unavoidable_n_fixed_ra_gt8_molecule_weight",
            "unavoidable_n_fixed_ra_lte8_span_gt8_patterns",
            "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight",
            "old_densest8_retained_patterns",
            "weight_stable_components",
            "weight_sensitive_components",
            "zero_evidence_blocks",
            "zero_evidence_block_sites",
            "tree_ready_blocks",
            "tree_ready_block_sites",
            "abstain_blocks",
            "abstain_block_sites",
            "partition_pattern_load_aggregate_seconds",
            "partition_ordered_hypergraph_dp_seconds",
            *EXTRACTION_METRIC_FIELDS,
        ):
            aggregate[field] += row[field]
        aggregate_status.update(row["status_counts"])
        aggregate_block_k.update(
            {int(k): value for k, value in row["block_k_distribution"].items()}
        )
        aggregate_block_spans.extend(row.pop("_block_span_bp_values"))
        aggregate_component_deltas.extend(
            row.pop("_component_retention_delta_values")
        )
    if aggregate["components"] != 408 or aggregate["sites"] != 47_570:
        raise ContractError("aggregate target totals differ from 408 components/47,570 sites")
    if aggregate["old_selected_sites"] != 408 * 8:
        raise ContractError("aggregate old densest-8 selection is not 3,264 sites")
    if (
        aggregate["old_selected_sites"] + aggregate["old_excluded_sites"]
        != aggregate["sites"]
        or aggregate["new_retained_sites"] != aggregate["sites"]
        or sum(k * count for k, count in aggregate_block_k.items())
        != aggregate["sites"]
    ):
        raise ContractError("aggregate site/block conservation failed")
    if (
        aggregate["raw_total_molecule_weight"]
        != aggregate["new_retained_molecule_weight"]
        + aggregate["new_lost_molecule_weight"]
    ):
        raise ContractError("aggregate molecule-weight conservation failed")
    if (
        aggregate["unavoidable_n_fixed_ra_gt8_patterns"]
        + aggregate["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"]
        != aggregate["unavoidable_patterns"]
        or aggregate["unavoidable_n_fixed_ra_gt8_molecule_weight"]
        + aggregate["unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight"]
        != aggregate["unavoidable_molecule_weight"]
    ):
        raise ContractError("aggregate unavoidable mechanism conservation failed")
    if (
        sum(aggregate_block_k.values()) != aggregate["new_blocks"]
        or aggregate["zero_evidence_blocks"] > aggregate["new_blocks"]
        or aggregate["zero_evidence_block_sites"] > aggregate["sites"]
    ):
        raise ContractError("aggregate block/zero-evidence conservation failed")
    if (
        aggregate["tree_ready_blocks"] + aggregate["abstain_blocks"]
        != aggregate["new_blocks"]
        or aggregate["tree_ready_block_sites"] + aggregate["abstain_block_sites"]
        != aggregate["sites"]
        or len(block_all) != aggregate["new_blocks"]
        or sum(_int(row["k"], "block_all k", minimum=1) for row in block_all)
        != aggregate["sites"]
    ):
        raise ContractError("aggregate tree-ready/abstain conservation failed")
    if (
        aggregate["eligible_alignments_pre_identity_collapse"]
        != aggregate["canonical_unique_molecules"]
        + aggregate["duplicate_alignment_identities_collapsed"]
        or aggregate["molecule_sparse_rows"]
        != aggregate["canonical_unique_molecules"]
        or aggregate["known_ps_hp12_molecule_rows"]
        != aggregate["known_ps_hp1_molecule_rows"]
        + aggregate["known_ps_hp2_molecule_rows"]
        or aggregate["extraction_failed_checks"] != 0
        or aggregate["duplicate_identity_conflicts"] != 0
    ):
        raise ContractError("aggregate extraction/collapse conservation failed")

    totals = {
        "chrom": "ALL",
        "partition_stage_status": "21_COMPLETED_1_ZERO_TARGET_SKIPPED",
        "ssnv_sites": 79_687,
        **dict(aggregate),
        "new_retention_ratio": _ratio(
            Decimal(aggregate["new_retained_molecule_weight"]),
            Decimal(aggregate["raw_total_molecule_weight"]),
        ),
        "old_densest8_retention_ratio": _ratio(
            Decimal(aggregate["old_densest8_retained_molecule_weight"]),
            Decimal(aggregate["raw_total_molecule_weight"]),
        ),
        "status_counts": dict(sorted(aggregate_status.items())),
        "block_k_distribution": {
            str(k): count for k, count in sorted(aggregate_block_k.items())
        },
        "block_span_bp_distribution": _distribution(aggregate_block_spans),
        "component_retention_delta_distribution": _distribution(
            aggregate_component_deltas,
            total_count=aggregate["components"],
        ),
        "extraction_wall_seconds": sum(
            float(row["extraction_wall_seconds"]) for row in per_chrom
        ),
        "partition_wall_seconds": sum(
            float(row["partition_wall_seconds"])
            for row in per_chrom
            if row["partition_wall_seconds"] != ""
        ),
        "extraction_max_rss_kb": max(
            resource["max_rss_kb"] for resource in extraction_resources
        ),
        "partition_max_rss_kb": max(
            resource["max_rss_kb"] for resource in partition_resources
        ),
    }

    if outer_time_path is None:
        outer_time_path = full_root.with_name(f"{full_root.name}.outer_time.txt")
    outer_resource = parse_gnu_time_v(outer_time_path)
    sequential_stage_wall_seconds = (
        totals["extraction_wall_seconds"] + totals["partition_wall_seconds"]
    )
    outer_resource.update(
        validate_outer_time_contract(
            outer_resource,
            full_root=full_root,
            runner_path=runner_path,
            manifest_path=manifest_path,
            sequential_stage_wall_seconds=sequential_stage_wall_seconds,
        )
    )
    resource_summary = {
        "extraction": {
            "stage_count": len(extraction_resources),
            "wall_seconds_sum": totals["extraction_wall_seconds"],
            "gnu_elapsed_seconds_sum": sum(
                item["elapsed_seconds"] for item in extraction_resources
            ),
            "user_seconds_sum": sum(item["user_seconds"] for item in extraction_resources),
            "system_seconds_sum": sum(
                item["system_seconds"] for item in extraction_resources
            ),
            "max_rss_kb_max": totals["extraction_max_rss_kb"],
            "filesystem_inputs_sum": sum(
                item["filesystem_inputs"] for item in extraction_resources
            ),
            "filesystem_outputs_sum": sum(
                item["filesystem_outputs"] for item in extraction_resources
            ),
        },
        "partition": {
            "stage_count": len(partition_resources),
            "wall_seconds_sum": totals["partition_wall_seconds"],
            "pattern_load_aggregate_seconds_sum": totals[
                "partition_pattern_load_aggregate_seconds"
            ],
            "ordered_hypergraph_dp_seconds_sum": totals[
                "partition_ordered_hypergraph_dp_seconds"
            ],
            "gnu_elapsed_seconds_sum": sum(
                item["elapsed_seconds"] for item in partition_resources
            ),
            "user_seconds_sum": sum(item["user_seconds"] for item in partition_resources),
            "system_seconds_sum": sum(
                item["system_seconds"] for item in partition_resources
            ),
            "max_rss_kb_max": totals["partition_max_rss_kb"],
            "filesystem_inputs_sum": sum(
                item["filesystem_inputs"] for item in partition_resources
            ),
            "filesystem_outputs_sum": sum(
                item["filesystem_outputs"] for item in partition_resources
            ),
        },
        "outer": outer_resource,
    }
    source = {
        "full_root": str(full_root.resolve()),
        "run_contract": {
            "path": str(contract_path.resolve()),
            "size_bytes": contract_path.stat().st_size,
            "sha256": contract_sha,
        },
        "run_receipt": {
            "path": str(receipt_path.resolve()),
            "size_bytes": receipt_path.stat().st_size,
            "sha256": receipt_sha,
        },
        "runner_chromosome_summary": {
            **runner_summary_identity,
            "sha256_sidecar_verified": runner_summary_sha,
        },
        "success_marker": file_identity(marker_path),
        "partition_receipts_verified": 21,
        "partition_tsv_identities_verified": 84,
        "chromosome_evidence": evidence_by_chrom,
    }
    payload = {
        "schema_name": SUMMARY_SCHEMA,
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "all_pass": True,
        "comprehensive_all_pass": True,
        "sample": DATASET,
        "scope": {
            "chromosomes": list(AUTOSOMES),
            "partitioned_chromosomes": [
                chrom for chrom in AUTOSOMES if chrom != "chr21"
            ],
            "zero_target_skipped_chromosome": "chr21",
        },
        "definitions": {
            "component": (
                "legacy adjacent-gap positional component with pre-cap k>8"
            ),
            "read_retention_denominator": (
                "component-level HP1/HP2 × known-PS exact-pattern molecule-count "
                "incidence; the same physical molecule may contribute once in each "
                "component it overlaps, so the aggregate is not a genome-wide unique-"
                "read count"
            ),
            "new_retention": (
                "within each component, molecule-count incidence whose complete set "
                "of fixed R/A target sites for that component constraint is contained "
                "in one selected k<=8 block, divided by the component-level incidence "
                "denominator"
            ),
            "old_densest8_retention": (
                "densest-contiguous-8 counterfactual recomputed from the same new "
                "BASEQ20 extraction and component-level incidence denominator; retained "
                "only when every fixed R/A target site in that component constraint is "
                "inside the selected eight sites; it is not an empirical retention "
                "measurement from the historical v5 run"
            ),
            "unavoidable": (
                "exact pattern whose ordered component span exceeds max block "
                "size 8; no valid disjoint k<=8 partition can retain it"
            ),
            "unavoidable_n_fixed_ra_gt8": (
                "unavoidable pattern with more than 8 fixed R/A calls in the "
                "read pattern itself"
            ),
            "unavoidable_n_fixed_ra_lte8_span_gt8": (
                "unavoidable pattern with at most 8 fixed R/A calls but endpoints "
                "spanning more than 8 ordered target sites because intervening "
                "targets must remain in contiguous blocks"
            ),
            "weight_stable": (
                "raw-molecule, equal-pattern, and log1p-molecule objectives "
                "select identical cut indices"
            ),
            "zero_evidence_block": (
                "selected k<=8 block whose raw_retained_molecule_weight is zero"
            ),
            "evidence_gate": (
                "TREE_READY_LOCAL only when retained_exact_pattern_count>0, "
                "raw_retained_molecule_weight>0, and "
                "primary_active_site_count>=2; otherwise "
                "ABSTAIN_ZERO_LOCAL_EVIDENCE. TREE_READY_LOCAL permits entry "
                "to candidate-tree inference only and is not evidence of a "
                "true or unique evolutionary tree"
            ),
            "tree_ready_block_sites": (
                "sum of k across TREE_READY_LOCAL blocks; this is block membership "
                "coverage and does not assert that every counted site has direct "
                "read support"
            ),
            "block_span_bp_distribution": (
                "block end1-start1; R-7 linear-interpolated p10/median/p90"
            ),
            "component_retention_delta_distribution": (
                "per-component new retention ratio minus the densest-8 counterfactual "
                "recomputed on the same BASEQ20 extraction and incidence denominator, "
                "excluding components with zero total molecule-count incidence; "
                "R-7 linear-interpolated p10/median/p90"
            ),
            "partition_pattern_load_aggregate_seconds": (
                "child receipt timings_seconds.load_and_aggregate_primary_patterns; "
                "loads sparse molecule calls and aggregates exact patterns, excluding "
                "BAM extraction and other partition-stage overhead"
            ),
            "partition_ordered_hypergraph_dp_seconds": (
                "child receipt timings_seconds.ordered_hypergraph_dp; three-weight "
                "partition component loop including raw/equal/log DP passes, the "
                "densest-8 baseline, diagnostics, row materialization, and aggregation; "
                "excludes pattern loading, BAM extraction, final output writes, and "
                "process overhead"
            ),
            "outer_runner_overhead_seconds": (
                "outer GNU time elapsed minus the sequential sum of extraction and "
                "partition stage walls; clipped at zero only for the reported overhead "
                "field, with a 0.25-second lower-bound tolerance for GNU time display "
                "rounding"
            ),
        },
        "totals": totals,
        "per_chromosome": per_chrom,
        "resources": resource_summary,
        "source_evidence": source,
        "checks": {
            "runner_receipt_and_success_marker_authenticated": True,
            "exact_canonical_408_components_47570_sites": True,
            "all_21_partition_receipts_authenticated": True,
            "all_84_partition_tsv_sha256_identities_verified": True,
            "outer_command_bound_to_runner_manifest_output_root": True,
            "outer_fresh_command_excludes_resume": True,
            "outer_elapsed_covers_sequential_stage_wall_sum": True,
            "component_metrics_rederived_from_tsv": True,
            "old_selected_plus_excluded_equals_target": True,
            "new_retained_sites_equal_target": True,
            "block_k_lte_8_and_sites_conserved": True,
            "raw_weight_equals_retained_plus_lost": True,
            "chr21_zero_target_skip_authenticated": True,
            "resource_exit_status_zero": True,
            "extraction_duplicate_collapse_conserved": True,
            "extraction_failed_checks_zero": True,
            "known_ps_hp1_hp2_rows_conserved": True,
            "zero_evidence_blocks_rederived": True,
            "block_span_distribution_rederived": True,
            "component_retention_delta_distribution_rederived": True,
            "tree_ready_plus_abstain_equals_all_blocks": True,
            "tree_ready_plus_abstain_sites_equals_all_block_sites": True,
            "partition_internal_timings_authenticated_and_aggregated": True,
            "partition_internal_timings_do_not_exceed_stage_wall": True,
            "unavoidable_mechanism_patterns_conserved": True,
            "unavoidable_mechanism_molecule_weight_conserved": True,
        },
        "claim_ceiling": (
            "validates deterministic engineering conservation and read-supported "
            "local k<=8 segmentation for the frozen HCC1395 run. "
            "TREE_READY_LOCAL only admits a block to candidate-tree inference; "
            "it does not identify or validate a true or unique evolutionary tree"
        ),
    }
    return payload, component_all, block_all


def _write_text_exclusive(path: Path, text: str, *, encoding: str = "utf-8") -> None:
    with path.open("x", encoding=encoding, newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def _write_sha_sidecar(path: Path) -> None:
    sidecar = path.with_name(f"{path.name}.sha256")
    _write_text_exclusive(sidecar, f"{sha256_path(path)}  {path.name}\n", encoding="ascii")


def _deterministic_gzip_writer(path: Path):
    raw = path.open("xb")
    compressed = gzip.GzipFile(
        filename="",
        mode="wb",
        compresslevel=6,
        fileobj=raw,
        mtime=0,
    )
    return raw, io.TextIOWrapper(compressed, encoding="utf-8", newline="")


def write_summary_outputs(
    output_dir: Path,
    payload: Mapping[str, Any],
    component_rows: Sequence[Mapping[str, str]],
    block_rows: Sequence[Mapping[str, str]],
) -> dict[str, dict[str, Any]]:
    if output_dir.exists() or output_dir.is_symlink():
        raise ContractError(f"refusing to overwrite output directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(exist_ok=False)

    summary_json = output_dir / "summary.json"
    _write_text_exclusive(
        summary_json,
        json.dumps(
            payload,
            ensure_ascii=False,
            allow_nan=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )

    summary_tsv = output_dir / "summary.tsv"
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=SUMMARY_FIELDS,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    summary_rows = [*payload["per_chromosome"], payload["totals"]]
    for raw in summary_rows:
        row = {field: raw.get(field, "") for field in SUMMARY_FIELDS}
        row["status_counts_json"] = _json_compact(raw["status_counts"])
        row["block_k_distribution_json"] = _json_compact(
            raw["block_k_distribution"]
        )
        row["block_span_bp_distribution_json"] = _json_compact(
            raw["block_span_bp_distribution"]
        )
        row["component_retention_delta_distribution_json"] = _json_compact(
            raw["component_retention_delta_distribution"]
        )
        writer.writerow(row)
    _write_text_exclusive(summary_tsv, buffer.getvalue())

    component_path = output_dir / "component_all.tsv.gz"
    raw, text = _deterministic_gzip_writer(component_path)
    try:
        component_writer = csv.DictWriter(
            text,
            fieldnames=COMPONENT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        component_writer.writeheader()
        component_writer.writerows(component_rows)
    finally:
        text.close()
        if not raw.closed:
            raw.close()

    block_path = output_dir / "block_all.tsv.gz"
    raw, text = _deterministic_gzip_writer(block_path)
    try:
        block_writer = csv.DictWriter(
            text,
            fieldnames=BLOCK_ALL_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        block_writer.writeheader()
        block_writer.writerows(block_rows)
    finally:
        text.close()
        if not raw.closed:
            raw.close()

    for path in (summary_json, summary_tsv, component_path, block_path):
        _write_sha_sidecar(path)
    return {
        path.name: file_identity(path)
        for path in (summary_json, summary_tsv, component_path, block_path)
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--outer-time",
        type=Path,
        help="GNU time -v file; defaults to <full-root>.outer_time.txt",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        full_root = args.full_root.resolve(strict=True)
        output_candidate = args.output_dir.resolve(strict=False)
        if output_candidate == full_root or full_root in output_candidate.parents:
            raise ContractError("summary output must be outside the immutable full-run root")
        payload, component_rows, block_rows = validate_full_run(
            full_root, args.outer_time
        )
        outputs = write_summary_outputs(
            args.output_dir, payload, component_rows, block_rows
        )
    except (ContractError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "all_pass": True,
                "comprehensive_all_pass": True,
                "components": payload["totals"]["components"],
                "sites": payload["totals"]["sites"],
                "outputs": outputs,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
