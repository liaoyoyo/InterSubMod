#!/usr/bin/env python3
"""Build a canonical Data Analytics report artifact from HCC1395 summary outputs.

The builder is intentionally read-only with respect to the full-run evidence tree.
It accepts immutable summarizer products, the separately audited historical
runtime baseline, exact-product log remediation, and observed-constraint HP/PS
unit retention evidence; it cross-checks their identities and conservation
contracts before writing one new ``artifact.json`` file.

Incomplete/subset inputs are rejected unless ``--fixture`` is explicit.  Fixture
artifacts are visibly labelled and must never be used as comprehensive evidence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Mapping, Sequence


SUMMARY_SCHEMA = "intersubmod.hcc1395_full_k_gt8_segmentation.summary"
SUMMARY_SCHEMA_VERSION = "1.1.0"
DATASET = "HCC1395"
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
EXPECTED_COMPONENTS = 408
EXPECTED_TARGET_SITES = 47_570
EXPECTED_SSNV_SITES = 79_687
SPAN_GRID_SCHEMA = "intersubmod.hcc1395_cached_span_sensitivity.run_receipt"
SPAN_GRID_SCHEMA_VERSION = "1.0.0"
EXACT_LOG_SCHEMA = "intersubmod.k_gt8_exact_log_sensitivity_remediation"
EXACT_LOG_SCHEMA_VERSION = "0.1.0"
HP_PS_UNIT_SCHEMA = "intersubmod.hp_ps_observed_constraint_retention_audit"
HP_PS_UNIT_SCHEMA_VERSION = "1.0.0"
CANONICAL_SPAN_CAPS = (50_000, 100_000, 200_000)
MAX_ROWS_PER_DATASET = 2_000
MAX_SNAPSHOT_BYTES = 3 * 1024 * 1024
MAX_INLINE_SOURCE_CHARS = 200_000
MAX_HP_PS_WORST_UNITS = 25

EXACT_LOG_REQUIRED_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "pre_cap_k",
    "exact_pattern_count",
    "legacy_log_matches_exact",
    "legacy_weight_stable_reported",
    "legacy_weight_stable_reconstructed",
    "corrected_weight_stable",
    "correction_changed_stability",
    "remediation_class",
)

HP_PS_UNIT_REQUIRED_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "hp_family",
    "phase_set",
    "unit_id",
    "component_k",
    "component_primary_active_site_count",
    "unit_active_site_count",
    "unit_active_site_fraction",
    "total_pattern_rows",
    "retained_pattern_rows",
    "cut_lost_pattern_rows",
    "unavoidable_pattern_rows",
    "nonretained_pattern_rows",
    "total_molecule_component_incidence_weight",
    "retained_molecule_component_incidence_weight",
    "cut_lost_molecule_component_incidence_weight",
    "unavoidable_molecule_component_incidence_weight",
    "nonretained_molecule_component_incidence_weight",
    "retention_ratio",
    "ratio_status",
    "support_stratum",
    "eligible_headline",
)

HP_PS_PAIR_REQUIRED_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "hp1_total_pattern_rows",
    "hp1_total_molecule_component_incidence_weight",
    "hp1_retained_molecule_component_incidence_weight",
    "hp1_retention_ratio",
    "hp2_total_pattern_rows",
    "hp2_total_molecule_component_incidence_weight",
    "hp2_retained_molecule_component_incidence_weight",
    "hp2_retention_ratio",
    "hp1_minus_hp2_retention_delta",
    "absolute_retention_delta",
    "both_hp_headline_eligible",
)

HP_PS_SOURCE_CHECKS = frozenset(
    {
        "component_active_site_sum_matches_receipt",
        "component_count_matches_receipt",
        "constraint_active_sites_match_membership_flags",
        "constraint_rows_match_receipt",
        "every_component_k_gt_max_block_size",
        "lost_weight_matches_receipt",
        "per_component_legacy_metrics_match_constraints",
        "retained_weight_matches_receipt",
        "site_count_matches_receipt",
        "total_weight_matches_receipt",
        "unavoidable_patterns_match_receipt",
    }
)

SUMMARY_REQUIRED_FIELDS = (
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
    "raw_total_molecule_weight",
    "new_retained_molecule_weight",
    "new_lost_molecule_weight",
    "unavoidable_patterns",
    "unavoidable_molecule_weight",
    "unavoidable_n_fixed_ra_gt8_patterns",
    "unavoidable_n_fixed_ra_gt8_molecule_weight",
    "unavoidable_n_fixed_ra_lte8_span_gt8_patterns",
    "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight",
    "new_retention_ratio",
    "old_densest8_retained_molecule_weight",
    "old_densest8_retention_ratio",
    "weight_stable_components",
    "weight_sensitive_components",
    "zero_evidence_blocks",
    "zero_evidence_block_sites",
    "tree_ready_blocks",
    "tree_ready_block_sites",
    "abstain_blocks",
    "abstain_block_sites",
    "extraction_wall_seconds",
    "partition_wall_seconds",
    "partition_pattern_load_aggregate_seconds",
    "partition_ordered_hypergraph_dp_seconds",
)

SUMMARY_TSV_CROSSCHECK_FIELDS = (
    "ssnv_sites",
    "components",
    "sites",
    "old_selected_sites",
    "old_excluded_sites",
    "new_blocks",
    "new_retained_sites",
    "primary_active_sites",
    "raw_total_molecule_weight",
    "new_retained_molecule_weight",
    "new_lost_molecule_weight",
    "unavoidable_patterns",
    "unavoidable_molecule_weight",
    "unavoidable_n_fixed_ra_gt8_patterns",
    "unavoidable_n_fixed_ra_gt8_molecule_weight",
    "unavoidable_n_fixed_ra_lte8_span_gt8_patterns",
    "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight",
    "new_retention_ratio",
    "old_densest8_retained_molecule_weight",
    "old_densest8_retention_ratio",
    "weight_stable_components",
    "weight_sensitive_components",
    "zero_evidence_blocks",
    "zero_evidence_block_sites",
    "tree_ready_blocks",
    "tree_ready_block_sites",
    "abstain_blocks",
    "abstain_block_sites",
    "extraction_wall_seconds",
    "partition_wall_seconds",
    "partition_pattern_load_aggregate_seconds",
    "partition_ordered_hypergraph_dp_seconds",
)

COMPONENT_REQUIRED_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
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
    "raw_total_molecule_weight",
    "raw_retained_molecule_weight",
    "raw_lost_molecule_weight",
    "raw_retention_ratio",
    "old_densest8_retained_molecule_weight",
    "old_densest8_retention_ratio",
    "weight_stable",
    "status",
    "positions_sha256",
)

SPAN_GRID_REQUIRED_FIELDS = (
    "span_cap_bp",
    "chrom",
    "status",
    "ssnv_sites",
    "k_gt8_components",
    "k_gt8_sites",
    "k_gt8_max_k",
    "wall_seconds",
    "new_blocks",
    "exact_patterns",
    "raw_total_molecule_weight",
    "raw_retained_molecule_weight",
    "raw_lost_molecule_weight",
    "unavoidable_patterns",
    "unavoidable_size_patterns",
    "unavoidable_span_cap_patterns",
    "unavoidable_both_limits_patterns",
)


class ReportInputError(RuntimeError):
    """Raised when source evidence cannot support the requested report state."""


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def require_regular_file(path: Path, label: str) -> Path:
    if path.is_symlink() or not path.is_file():
        raise ReportInputError(f"{label} is missing, nonregular, or a symlink: {path}")
    return path


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportInputError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ReportInputError(f"JSON root must be an object: {path}")
    return value


def read_tsv(path: Path, *, compressed: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if compressed else open
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or ())
            rows = [dict(row) for row in reader]
    except (OSError, csv.Error) as exc:
        raise ReportInputError(f"cannot read TSV {path}: {exc}") from exc
    if not fields:
        raise ReportInputError(f"TSV has no header: {path}")
    return fields, rows


def as_decimal(value: Any, label: str, *, blank_ok: bool = False) -> Decimal | None:
    if value is None or value == "":
        if blank_ok:
            return None
        raise ReportInputError(f"{label} is blank")
    if isinstance(value, bool):
        raise ReportInputError(f"{label} must be numeric, not boolean")
    try:
        parsed = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise ReportInputError(f"{label} is not numeric: {value!r}") from exc
    if not parsed.is_finite():
        raise ReportInputError(f"{label} must be finite")
    return parsed


def as_int(value: Any, label: str) -> int:
    parsed = as_decimal(value, label)
    assert parsed is not None
    if parsed != parsed.to_integral_value():
        raise ReportInputError(f"{label} must be an integer: {value!r}")
    return int(parsed)


def as_float(value: Any, label: str, *, blank_ok: bool = False) -> float | None:
    parsed = as_decimal(value, label, blank_ok=blank_ok)
    return None if parsed is None else float(parsed)


def as_bool(value: Any, label: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    raise ReportInputError(f"{label} is not boolean-like: {value!r}")


def ratio(numerator: int | float, denominator: int | float) -> float | None:
    if denominator == 0:
        return None
    return float(numerator) / float(denominator)


def validate_component_retention(
    row: Mapping[str, Any], label: str
) -> tuple[float | None, float | None, float | None]:
    total = as_int(row["raw_total_molecule_weight"], f"{label}.raw_total")
    new_retained = as_int(
        row["raw_retained_molecule_weight"], f"{label}.new_retained"
    )
    old_retained = as_int(
        row["old_densest8_retained_molecule_weight"], f"{label}.old_retained"
    )
    new_ratio = as_float(
        row["raw_retention_ratio"], f"{label}.raw_retention_ratio", blank_ok=True
    )
    old_ratio = as_float(
        row["old_densest8_retention_ratio"],
        f"{label}.old_densest8_retention_ratio",
        blank_ok=True,
    )
    expected_new = ratio(new_retained, total)
    expected_old = ratio(old_retained, total)
    for name, observed, expected in (
        ("raw_retention_ratio", new_ratio, expected_new),
        ("old_densest8_retention_ratio", old_ratio, expected_old),
    ):
        if observed is None or expected is None:
            if observed is not expected:
                raise ReportInputError(
                    f"{label}.{name} must be blank exactly when the denominator is zero"
                )
        elif not math.isclose(observed, expected, rel_tol=0, abs_tol=1e-12):
            raise ReportInputError(f"{label}.{name} does not match retained/total")
    gain = None if new_ratio is None else new_ratio - old_ratio
    return new_ratio, old_ratio, gain


def numeric_equal(left: Any, right: Any, label: str) -> bool:
    left_value = as_decimal(left, f"{label} JSON", blank_ok=True)
    right_value = as_decimal(right, f"{label} TSV", blank_ok=True)
    if left_value is None or right_value is None:
        return left_value is right_value
    return abs(left_value - right_value) <= Decimal("1e-12")


def chrom_number(chrom: str) -> int:
    match = re.fullmatch(r"chr([1-9]|1[0-9]|2[0-2])", chrom)
    if not match:
        raise ReportInputError(f"unsupported chromosome label: {chrom!r}")
    return int(match.group(1))


def safe_source_path(path: Path) -> str:
    """Return a package-safe local filename, never a machine-local absolute path."""
    name = path.name
    if not name or name in {".", ".."} or "/" in name or "\\" in name:
        raise ReportInputError(f"unsafe source filename: {path}")
    return name


def parse_generated_at(value: str | None) -> str:
    if value is None:
        return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
            "+00:00", "Z"
        )
    normalized = value[:-1] + "+00:00" if value.endswith("Z") else value
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError as exc:
        raise ReportInputError("--generated-at must be ISO-8601") from exc
    if parsed.tzinfo is None:
        raise ReportInputError("--generated-at must include a timezone")
    return parsed.astimezone(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def parse_baseline_runtime(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    patterns = (
        r"\|\s*\*\*HCC1395 downstream total\*\*\s*\|\s*\*\*`?([0-9][0-9,.]+)`?\*\*",
        r"高精度時間差為\s*([0-9][0-9,.]+)\s*秒",
    )
    seconds: float | None = None
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            seconds = float(match.group(1).replace(",", ""))
            break
    if seconds is None or not math.isfinite(seconds) or seconds <= 0:
        raise ReportInputError(
            "baseline-runtime-audit.md does not expose the audited HCC1395 seconds"
        )
    is_proxy = bool(
        re.search(r"filesystem birth-timestamp proxy|historical wall proxy", text, re.I)
    )
    if not is_proxy:
        raise ReportInputError(
            "baseline runtime lacks its required historical wall-proxy qualification"
        )
    return {
        "seconds": seconds,
        "minutes": seconds / 60.0,
        "is_historical_wall_proxy": True,
    }


def resolve_span_grid_paths(candidate: Path) -> tuple[Path, Path]:
    """Resolve a span-grid root, summary.tsv, or receipt.json to both inputs."""
    if candidate.is_symlink():
        raise ReportInputError(f"span-grid input must not be a symlink: {candidate}")
    if candidate.is_dir():
        summary_path = candidate / "summary.tsv"
        receipt_path = candidate / "receipt.json"
    elif candidate.name == "summary.tsv":
        summary_path = candidate
        receipt_path = candidate.parent / "receipt.json"
    elif candidate.name == "receipt.json":
        receipt_path = candidate
        summary_path = candidate.parent / "summary.tsv"
    else:
        raise ReportInputError(
            "--span-grid-summary must point to the run root, summary.tsv, or receipt.json"
        )
    return (
        require_regular_file(summary_path, "span-grid summary.tsv"),
        require_regular_file(receipt_path, "span-grid receipt.json"),
    )


def verify_file_identity(
    identity: Any,
    *,
    expected_path: Path,
    label: str,
) -> None:
    if not isinstance(identity, Mapping):
        raise ReportInputError(f"{label} identity must be an object")
    if set(identity) != {"path", "size_bytes", "sha256"}:
        raise ReportInputError(f"{label} identity keys differ")
    try:
        identity_path = Path(str(identity["path"])).resolve()
    except OSError as exc:
        raise ReportInputError(f"{label} path cannot be resolved") from exc
    if identity_path != expected_path.resolve():
        raise ReportInputError(f"{label} identity path does not match supplied input")
    if as_int(identity["size_bytes"], f"{label}.size_bytes") != expected_path.stat().st_size:
        raise ReportInputError(f"{label} size identity mismatch")
    observed_sha = sha256_path(expected_path)
    if identity["sha256"] != observed_sha:
        raise ReportInputError(f"{label} SHA-256 identity mismatch")


def verify_recorded_file_identity(identity: Any, *, label: str) -> Path:
    """Verify a receipt-recorded immutable file without trusting its pathname."""
    if not isinstance(identity, Mapping):
        raise ReportInputError(f"{label} identity must be an object")
    if set(identity) != {"path", "size_bytes", "sha256"}:
        raise ReportInputError(f"{label} identity keys differ")
    path = require_regular_file(Path(str(identity["path"])), label)
    if as_int(identity["size_bytes"], f"{label}.size_bytes") != path.stat().st_size:
        raise ReportInputError(f"{label} recorded size identity mismatch")
    if not re.fullmatch(r"[0-9a-f]{64}", str(identity["sha256"])):
        raise ReportInputError(f"{label} recorded SHA-256 is malformed")
    if sha256_path(path) != identity["sha256"]:
        raise ReportInputError(f"{label} recorded SHA-256 identity mismatch")
    return path


def validate_producer_identity_shape(identity: Any, *, label: str) -> None:
    """Keep producer provenance without requiring the live script to stay frozen."""
    if not isinstance(identity, Mapping) or set(identity) != {
        "path",
        "size_bytes",
        "sha256",
    }:
        raise ReportInputError(f"{label} identity keys differ")
    if as_int(identity["size_bytes"], f"{label}.size_bytes") <= 0:
        raise ReportInputError(f"{label}.size_bytes must be positive")
    if not re.fullmatch(r"[0-9a-f]{64}", str(identity["sha256"])):
        raise ReportInputError(f"{label}.sha256 is malformed")


def require_all_true(value: Any, *, label: str) -> Mapping[str, Any]:
    if (
        not isinstance(value, Mapping)
        or not value
        or any(item is not True for item in value.values())
    ):
        raise ReportInputError(f"{label} must be a nonempty all-true object")
    return value


def quantile_type7(values: Sequence[float], probability: float) -> float:
    """R/NumPy-compatible type-7 linear quantile for audited unit distributions."""
    if not values:
        raise ReportInputError("cannot compute a quantile from zero values")
    if not 0 <= probability <= 1:
        raise ReportInputError("quantile probability must be between zero and one")
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def validate_hp_ps_summary_slice(
    recorded: Any,
    units: Sequence[Mapping[str, Any]],
    *,
    label: str,
) -> None:
    """Cross-check one chromosome/HP/support-stratum summary against unit rows."""
    if not isinstance(recorded, Mapping):
        raise ReportInputError(f"{label} summary slice must be an object")
    expected_eligible = sum(bool(unit["eligible_headline"]) for unit in units)
    if (
        as_int(recorded.get("observed_constraint_units"), f"{label}.units")
        != len(units)
        or as_int(recorded.get("eligible_headline_units"), f"{label}.eligible")
        != expected_eligible
    ):
        raise ReportInputError(f"{label} unit counts differ")
    incidence = recorded.get("molecule_component_incidences")
    patterns = recorded.get("pattern_rows")
    if not isinstance(incidence, Mapping) or not isinstance(patterns, Mapping):
        raise ReportInputError(f"{label} mass summaries are required")
    incidence_expected = {
        "total": sum(int(unit["total_weight"]) for unit in units),
        "retained": sum(int(unit["retained_weight"]) for unit in units),
        "cut_lost": sum(int(unit["cut_lost_weight"]) for unit in units),
        "unavoidable": sum(int(unit["unavoidable_weight"]) for unit in units),
    }
    pattern_expected = {
        "total": sum(int(unit["total_pattern_rows"]) for unit in units),
        "retained": sum(int(unit["retained_pattern_rows"]) for unit in units),
        "cut_lost": sum(int(unit["cut_lost_pattern_rows"]) for unit in units),
        "unavoidable": sum(int(unit["unavoidable_pattern_rows"]) for unit in units),
    }
    for name, observed in incidence_expected.items():
        if as_int(incidence.get(name), f"{label}.incidence.{name}") != observed:
            raise ReportInputError(f"{label} incidence mass differs: {name}")
    for name, observed in pattern_expected.items():
        if as_int(patterns.get(name), f"{label}.patterns.{name}") != observed:
            raise ReportInputError(f"{label} pattern mass differs: {name}")
    weighted = as_float(
        incidence.get("weighted_retention_ratio"), f"{label}.weighted_retention"
    )
    expected_weighted = (
        incidence_expected["retained"] / incidence_expected["total"]
    )
    if weighted is None or not math.isclose(
        weighted, expected_weighted, rel_tol=0, abs_tol=1e-12
    ):
        raise ReportInputError(f"{label} weighted retention differs")


def verify_sha256_sidecar(sidecar: Path, target: Path, *, label: str) -> None:
    text = sidecar.read_text(encoding="utf-8").strip()
    match = re.fullmatch(r"([0-9a-f]{64})(?:\s+\*?([^\s]+))?", text)
    if not match or match.group(1) != sha256_path(target):
        raise ReportInputError(f"{label} sidecar digest mismatch")
    if match.group(2) is not None and Path(match.group(2)).name != target.name:
        raise ReportInputError(f"{label} sidecar filename mismatch")


def resolve_exact_log_paths(candidate: Path) -> dict[str, Path]:
    """Resolve an exact-log audit root, receipt.json, or summary.json."""
    if candidate.is_symlink():
        raise ReportInputError(f"exact-log audit input must not be a symlink: {candidate}")
    if candidate.is_dir():
        root = candidate
    elif candidate.name in {"receipt.json", "summary.json"}:
        root = candidate.parent
    else:
        raise ReportInputError(
            "--exact-log-audit must point to an audit root, receipt.json, or summary.json"
        )
    return {
        "root": root,
        "receipt": require_regular_file(root / "receipt.json", "exact-log receipt"),
        "receipt_sha": require_regular_file(
            root / "receipt.json.sha256", "exact-log receipt SHA sidecar"
        ),
        "summary": require_regular_file(root / "summary.json", "exact-log summary"),
        "detail": require_regular_file(
            root / "exact_log_sensitivity.tsv.gz", "exact-log detail"
        ),
        "success": require_regular_file(root / "_SUCCESS", "exact-log success marker"),
    }


def resolve_hp_ps_unit_paths(candidate: Path) -> dict[str, Path]:
    """Resolve an HP/PS unit-audit root, receipt.json, or summary.json."""
    if candidate.is_symlink():
        raise ReportInputError(f"HP/PS unit audit input must not be a symlink: {candidate}")
    if candidate.is_dir():
        root = candidate
    elif candidate.name in {"receipt.json", "summary.json"}:
        root = candidate.parent
    else:
        raise ReportInputError(
            "--hp-ps-unit-audit must point to an audit root, receipt.json, or summary.json"
        )
    return {
        "root": root,
        "receipt": require_regular_file(root / "receipt.json", "HP/PS unit receipt"),
        "receipt_sha": require_regular_file(
            root / "receipt.json.sha256", "HP/PS receipt SHA sidecar"
        ),
        "summary": require_regular_file(root / "summary.json", "HP/PS unit summary"),
        "summary_tsv": require_regular_file(root / "summary.tsv", "HP/PS unit summary TSV"),
        "units": require_regular_file(
            root / f"{DATASET}.hp_ps_observed_units.tsv.gz",
            "HP/PS observed-unit table",
        ),
        "pairs": require_regular_file(
            root / f"{DATASET}.hp1_hp2_paired_components.tsv.gz",
            "HP1/HP2 paired-component table",
        ),
    }


def _validate_exact_source_identities(receipt: Mapping[str, Any]) -> None:
    sources = receipt.get("sources")
    if not isinstance(sources, Mapping):
        raise ReportInputError("exact-log receipt.sources must be an object")
    full_root = sources.get("full_root")
    if full_root is not None:
        if not isinstance(full_root, Mapping) or set(full_root) != {
            "marker",
            "receipt",
            "totals",
        }:
            raise ReportInputError("exact-log full-root identity contract differs")
        verify_recorded_file_identity(
            full_root["marker"], label="exact-log full-root success marker"
        )
        verify_recorded_file_identity(
            full_root["receipt"], label="exact-log full-root receipt"
        )
        expected_totals = {
            "chromosomes": 22,
            "ssnv_sites": EXPECTED_SSNV_SITES,
            "k_gt8_components": EXPECTED_COMPONENTS,
            "k_gt8_sites": EXPECTED_TARGET_SITES,
            "k_gt8_max_k": 3_574,
            "partitioned_chromosomes": 21,
            "zero_target_skipped_chromosomes": 1,
        }
        if full_root["totals"] != expected_totals:
            raise ReportInputError("exact-log full-root canonical totals differ")
    partitions = sources.get("partitions")
    if not isinstance(partitions, list) or not partitions:
        raise ReportInputError("exact-log receipt must bind source partitions")
    required = {
        "partition_receipt",
        "legacy_components",
        "cut_constraints",
        "site_membership_coordinate_witness",
    }
    for index, source in enumerate(partitions, start=1):
        if not isinstance(source, Mapping):
            raise ReportInputError(f"exact-log source[{index}] must be an object")
        if source.get("dataset") != DATASET:
            raise ReportInputError(f"exact-log source[{index}] dataset differs")
        chrom_number(str(source.get("chrom", "")))
        missing = required - set(source)
        if missing:
            raise ReportInputError(
                f"exact-log source[{index}] is missing identities: {sorted(missing)}"
            )
        for name in sorted(required):
            verify_recorded_file_identity(
                source[name], label=f"exact-log source[{index}].{name}"
            )


def validate_exact_log_audits(
    candidates: Sequence[Path],
    *,
    fixture: bool,
    component_rows: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    """Authenticate exact-product remediation and bind it to every report component."""
    if not candidates:
        raise ReportInputError("at least one exact-log audit is required")
    component_keys = {
        (str(row["chrom"]), str(row["legacy_component_id"])) for row in component_rows
    }
    if len(component_keys) != len(component_rows):
        raise ReportInputError("report component IDs are not unique")

    merged_rows: dict[tuple[str, str], dict[str, Any]] = {}
    aggregate_counts = {
        "components": 0,
        "sites": 0,
        "patterns": 0,
        "legacy_log_matches_exact": 0,
        "legacy_log_differs_exact": 0,
        "legacy_weight_stable": 0,
        "corrected_weight_stable": 0,
        "correction_changed_stability": 0,
    }
    chromosome_counts: dict[str, dict[str, int]] = {}
    paths: list[dict[str, Path]] = []
    modes: list[str] = []

    for audit_index, candidate in enumerate(candidates, start=1):
        resolved = resolve_exact_log_paths(candidate)
        receipt = load_json(resolved["receipt"])
        summary = load_json(resolved["summary"])
        verify_sha256_sidecar(
            resolved["receipt_sha"],
            resolved["receipt"],
            label="exact-log receipt",
        )
        success = load_json(resolved["success"])
        success_identity = success.get("receipt")
        if (
            success.get("schema_name") != f"{EXACT_LOG_SCHEMA}.success_marker"
            or success.get("schema_version") != EXACT_LOG_SCHEMA_VERSION
            or success.get("all_pass") is not True
            or not isinstance(success_identity, Mapping)
            or Path(str(success_identity.get("path", ""))).resolve()
            != resolved["receipt"].resolve()
            or success_identity.get("sha256") != sha256_path(resolved["receipt"])
        ):
            raise ReportInputError("exact-log _SUCCESS does not bind receipt.json")
        if receipt.get("schema_name") != EXACT_LOG_SCHEMA:
            raise ReportInputError("unexpected exact-log receipt schema")
        if receipt.get("schema_version") != EXACT_LOG_SCHEMA_VERSION:
            raise ReportInputError("unexpected exact-log receipt schema version")
        if receipt.get("all_pass") is not True:
            raise ReportInputError("exact-log receipt must be all_pass")
        require_all_true(receipt.get("checks"), label="exact-log receipt.checks")
        validate_producer_identity_shape(receipt.get("script"), label="exact-log script")
        _validate_exact_source_identities(receipt)

        outputs = receipt.get("outputs")
        if not isinstance(outputs, Mapping):
            raise ReportInputError("exact-log receipt.outputs must be an object")
        for output_name, path_key in (
            ("exact_log_sensitivity", "detail"),
            ("summary", "summary"),
        ):
            verify_file_identity(
                outputs.get(output_name),
                expected_path=resolved[path_key],
                label=f"exact-log {output_name}",
            )
        if "summary_sha256" not in outputs:
            raise ReportInputError("exact-log receipt does not bind summary SHA sidecar")
        verify_recorded_file_identity(
            outputs["summary_sha256"], label="exact-log summary SHA sidecar"
        )

        mode = str(receipt.get("mode", ""))
        if mode not in {"single_chromosome", "comprehensive_chr1_22"}:
            raise ReportInputError(f"unsupported exact-log mode: {mode!r}")
        if (
            mode == "comprehensive_chr1_22"
            and not isinstance(receipt.get("sources", {}).get("full_root"), Mapping)
        ):
            raise ReportInputError(
                "comprehensive exact-log audit must bind the authenticated full root"
            )
        modes.append(mode)
        scope = receipt.get("scope")
        counts = receipt.get("counts")
        if not isinstance(scope, Mapping) or not isinstance(counts, Mapping):
            raise ReportInputError("exact-log scope/counts must be objects")
        scope_chroms = tuple(str(chrom) for chrom in scope.get("chromosomes", ()))
        if not scope_chroms or len(scope_chroms) != len(set(scope_chroms)):
            raise ReportInputError("exact-log chromosome scope must be nonempty and unique")
        for chrom in scope_chroms:
            chrom_number(chrom)
        if tuple(scope.get("datasets", ())) != (DATASET,):
            raise ReportInputError("exact-log scope dataset must be HCC1395")

        if summary.get("schema_name") != f"{EXACT_LOG_SCHEMA}.summary":
            raise ReportInputError("unexpected exact-log summary schema")
        if (
            summary.get("schema_version") != EXACT_LOG_SCHEMA_VERSION
            or summary.get("all_pass") is not True
            or summary.get("mode") != mode
            or summary.get("semantic_result_sha256")
            != receipt.get("semantic_result_sha256")
        ):
            raise ReportInputError("exact-log summary/receipt binding differs")
        if summary.get("counts") != counts:
            raise ReportInputError("exact-log summary/receipt counts differ")
        if summary.get("chromosome_counts") != receipt.get("chromosome_counts"):
            raise ReportInputError("exact-log chromosome counts differ")
        interpretation = summary.get("interpretation")
        parameters = receipt.get("parameters")
        if (
            not isinstance(interpretation, Mapping)
            or interpretation.get("primary_raw_partition_changed") is not False
            or not isinstance(parameters, Mapping)
            or parameters.get("primary_raw_partition_changed") is not False
        ):
            raise ReportInputError(
                "exact-log remediation must state primary raw cuts are unchanged"
            )

        fields, rows = read_tsv(resolved["detail"], compressed=True)
        missing = [field for field in EXACT_LOG_REQUIRED_FIELDS if field not in fields]
        if missing:
            raise ReportInputError(
                "exact-log detail is missing fields: " + ", ".join(missing)
            )
        local_counts = {name: 0 for name in aggregate_counts}
        local_by_chrom: dict[str, dict[str, int]] = {}
        for row_index, row in enumerate(rows, start=1):
            label = f"exact-log audit[{audit_index}] row[{row_index}]"
            if row["dataset"] != DATASET:
                raise ReportInputError(f"{label} dataset differs")
            chrom = row["chrom"]
            chrom_number(chrom)
            if chrom not in scope_chroms:
                raise ReportInputError(f"{label} chromosome escaped receipt scope")
            key = (chrom, row["legacy_component_id"])
            if key in merged_rows:
                raise ReportInputError(f"duplicate exact-log component: {key}")
            pre_cap_k = as_int(row["pre_cap_k"], f"{label}.pre_cap_k")
            patterns = as_int(row["exact_pattern_count"], f"{label}.exact_pattern_count")
            if pre_cap_k <= 8 or patterns < 0:
                raise ReportInputError(f"{label} has invalid k/pattern count")
            legacy_matches = as_bool(
                row["legacy_log_matches_exact"], f"{label}.legacy_log_matches_exact"
            )
            legacy_reported = as_bool(
                row["legacy_weight_stable_reported"],
                f"{label}.legacy_weight_stable_reported",
            )
            legacy_reconstructed = as_bool(
                row["legacy_weight_stable_reconstructed"],
                f"{label}.legacy_weight_stable_reconstructed",
            )
            corrected = as_bool(
                row["corrected_weight_stable"], f"{label}.corrected_weight_stable"
            )
            changed = as_bool(
                row["correction_changed_stability"],
                f"{label}.correction_changed_stability",
            )
            if legacy_reported != legacy_reconstructed:
                raise ReportInputError(f"{label} legacy stability reconstruction differs")
            if changed != (corrected != legacy_reported):
                raise ReportInputError(f"{label} correction-change flag differs")
            normalized = {
                "legacy_log_matches_exact": legacy_matches,
                "legacy_log_differs_exact": not legacy_matches,
                "legacy_weight_stable": legacy_reported,
                "corrected_weight_stable": corrected,
                "correction_changed_stability": changed,
                "remediation_class": row["remediation_class"],
            }
            merged_rows[key] = normalized
            row_counts = {
                "components": 1,
                "sites": pre_cap_k,
                "patterns": patterns,
                "legacy_log_matches_exact": int(legacy_matches),
                "legacy_log_differs_exact": int(not legacy_matches),
                "legacy_weight_stable": int(legacy_reported),
                "corrected_weight_stable": int(corrected),
                "correction_changed_stability": int(changed),
            }
            chrom_counts = local_by_chrom.setdefault(
                chrom, {name: 0 for name in aggregate_counts}
            )
            for name, value in row_counts.items():
                local_counts[name] += value
                chrom_counts[name] += value

        for name, observed in local_counts.items():
            if as_int(counts.get(name), f"exact-log counts.{name}") != observed:
                raise ReportInputError(f"exact-log count mismatch: {name}")
        if (
            as_int(scope.get("components"), "exact-log scope.components")
            != local_counts["components"]
            or as_int(scope.get("sites"), "exact-log scope.sites")
            != local_counts["sites"]
        ):
            raise ReportInputError("exact-log scope totals differ from detail")
        receipt_by_chrom = receipt.get("chromosome_counts")
        if not isinstance(receipt_by_chrom, Mapping):
            raise ReportInputError("exact-log chromosome_counts must be an object")
        for chrom, observed_counts in local_by_chrom.items():
            expected_counts = receipt_by_chrom.get(chrom)
            if not isinstance(expected_counts, Mapping):
                raise ReportInputError(f"exact-log chromosome count missing: {chrom}")
            for name, observed in observed_counts.items():
                if as_int(
                    expected_counts.get(name), f"exact-log {chrom}.{name}"
                ) != observed:
                    raise ReportInputError(
                        f"exact-log chromosome count mismatch: {chrom}.{name}"
                    )
            if chrom in chromosome_counts:
                raise ReportInputError(f"duplicate exact-log chromosome: {chrom}")
            chromosome_counts[chrom] = observed_counts
        for name, observed in local_counts.items():
            aggregate_counts[name] += observed
        paths.append(resolved)

    if set(merged_rows) != component_keys:
        missing = len(component_keys - set(merged_rows))
        extra = len(set(merged_rows) - component_keys)
        raise ReportInputError(
            f"exact-log component coverage differs from report: missing={missing}, extra={extra}"
        )
    if not fixture:
        if len(candidates) != 1 or modes != ["comprehensive_chr1_22"]:
            raise ReportInputError(
                "formal exact-log audit must be one comprehensive_chr1_22 receipt"
            )
        if (
            aggregate_counts["components"] != EXPECTED_COMPONENTS
            or aggregate_counts["sites"] != EXPECTED_TARGET_SITES
            or tuple(chromosome_counts)
            != tuple(chrom for chrom in AUTOSOMES if chrom != "chr21")
        ):
            raise ReportInputError(
                "formal exact-log audit must conserve 408 components / 47,570 sites "
                "across the 21 nonempty autosomes, with chr21 authenticated by full-root"
            )
    return {
        "component_rows": merged_rows,
        "counts": aggregate_counts,
        "chromosome_counts": chromosome_counts,
        "paths": paths,
        "modes": modes,
        "primary_raw_partition_changed": False,
    }


def _validate_hp_ps_source_identities(receipt: Mapping[str, Any]) -> None:
    inputs = receipt.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ReportInputError("HP/PS receipt.inputs must be a nonempty array")
    required = {
        "partition_receipt",
        "legacy_components",
        "site_membership",
        "cut_constraints",
    }
    for index, source in enumerate(inputs, start=1):
        if not isinstance(source, Mapping) or source.get("dataset") != DATASET:
            raise ReportInputError(f"HP/PS input[{index}] scope differs")
        chrom_number(str(source.get("chrom", "")))
        missing = required - set(source)
        if missing:
            raise ReportInputError(
                f"HP/PS input[{index}] is missing identities: {sorted(missing)}"
            )
        for name in sorted(required):
            verify_recorded_file_identity(
                source[name], label=f"HP/PS input[{index}].{name}"
            )


def validate_hp_ps_unit_audit(
    candidate: Path,
    *,
    fixture: bool,
    component_rows: Sequence[Mapping[str, str]],
    report_chromosomes: Sequence[str],
) -> dict[str, Any]:
    """Authenticate the v5 observed-constraint HP/PS unit retention contract."""
    resolved = resolve_hp_ps_unit_paths(candidate)
    receipt = load_json(resolved["receipt"])
    summary = load_json(resolved["summary"])
    verify_sha256_sidecar(
        resolved["receipt_sha"], resolved["receipt"], label="HP/PS receipt"
    )
    if receipt.get("schema_name") != HP_PS_UNIT_SCHEMA:
        raise ReportInputError("unexpected HP/PS unit receipt schema")
    if receipt.get("schema_version") != HP_PS_UNIT_SCHEMA_VERSION:
        raise ReportInputError("unexpected HP/PS unit receipt schema version")
    if receipt.get("all_pass") is not True:
        raise ReportInputError("HP/PS unit receipt must be all_pass")
    checks = require_all_true(receipt.get("checks"), label="HP/PS receipt.checks")
    required_final_checks = {
        "full_scope_is_exact_autosomes",
        "no_unobserved_unit_rows_synthesized",
        "source_checks_all_pass",
        "unit_molecule_incidence_mass_conserved",
        "unit_pattern_mass_conserved",
        "unit_rows_have_observed_denominator",
    }
    if set(checks) != required_final_checks:
        raise ReportInputError("HP/PS final-check contract differs")
    validate_producer_identity_shape(receipt.get("tool"), label="HP/PS producer")
    _validate_hp_ps_source_identities(receipt)

    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping):
        raise ReportInputError("HP/PS receipt.outputs must be an object")
    for output_name, path_key in (
        ("summary_json", "summary"),
        ("summary_tsv", "summary_tsv"),
        ("unit_table", "units"),
        ("paired_component_table", "pairs"),
    ):
        verify_file_identity(
            outputs.get(output_name),
            expected_path=resolved[path_key],
            label=f"HP/PS {output_name}",
        )
    if receipt.get("summary") != summary:
        raise ReportInputError("HP/PS receipt.summary differs from authenticated summary.json")
    summary_scope = summary.get("scope_contract")
    if (
        not isinstance(summary_scope, Mapping)
        or summary_scope.get("aggregation_weight")
        != "molecule_x_component_incidence"
        or summary_scope.get("scope_ceiling")
        != "observed_constraint_units_only"
        or "synthetic" not in str(
            summary_scope.get("unobserved_opportunity_policy", "")
        ).lower()
    ):
        raise ReportInputError("HP/PS summary scope/no-synthetic contract differs")

    scope = receipt.get("scope")
    definitions = receipt.get("definitions")
    parameters = receipt.get("parameters")
    if not all(isinstance(value, Mapping) for value in (scope, definitions, parameters)):
        raise ReportInputError("HP/PS scope/definitions/parameters must be objects")
    assert isinstance(scope, Mapping)
    assert isinstance(definitions, Mapping)
    assert isinstance(parameters, Mapping)
    if (
        scope.get("dataset") != DATASET
        or scope.get("scope_ceiling") != "observed_constraint_units_only"
        or definitions.get("aggregation_weight") != "molecule_x_component_incidence"
        or as_int(
            parameters.get("eligible_min_total_weight"),
            "HP/PS eligible_min_total_weight",
        )
        != 20
        or as_int(
            parameters.get("eligible_min_pattern_rows"),
            "HP/PS eligible_min_pattern_rows",
        )
        != 5
    ):
        raise ReportInputError("HP/PS scope or eligibility contract differs")
    target_chroms = tuple(
        str(chrom) for chrom in scope.get("chromosomes_with_partition_targets", ())
    )
    skipped_chroms = tuple(
        str(chrom) for chrom in scope.get("skipped_zero_target_chromosomes", ())
    )
    if len(set((*target_chroms, *skipped_chroms))) != len(
        (*target_chroms, *skipped_chroms)
    ):
        raise ReportInputError("HP/PS chromosome scope contains duplicates")
    for chrom in (*target_chroms, *skipped_chroms):
        chrom_number(chrom)
    if not fixture:
        if (
            scope.get("mode") != "full"
            or tuple(sorted((*target_chroms, *skipped_chroms), key=chrom_number))
            != AUTOSOMES
        ):
            raise ReportInputError(
                "formal HP/PS audit requires mode=full and exact chr1-chr22 scope"
            )
    elif not set(target_chroms).issubset(set(report_chromosomes)):
        raise ReportInputError("fixture HP/PS target chromosomes escaped report scope")

    source_checks = receipt.get("source_checks")
    if not isinstance(source_checks, Mapping) or set(source_checks) != set(target_chroms):
        raise ReportInputError("HP/PS source-check chromosome scope differs")
    for chrom, chrom_checks in source_checks.items():
        if (
            not isinstance(chrom_checks, Mapping)
            or set(chrom_checks) != HP_PS_SOURCE_CHECKS
            or any(value is not True for value in chrom_checks.values())
        ):
            raise ReportInputError(
                f"HP/PS {chrom} must pass the 11-source-check v5 contract"
            )

    summary_tsv_fields, summary_tsv_rows = read_tsv(resolved["summary_tsv"])
    required_summary_tsv_fields = {
        "scope_type",
        "scope_id",
        "metric",
        "value",
        "unit",
        "denominator_note",
    }
    if not required_summary_tsv_fields.issubset(summary_tsv_fields):
        raise ReportInputError("HP/PS summary.tsv v5 fields differ")
    unit_fields, raw_unit_rows = read_tsv(resolved["units"], compressed=True)
    pair_fields, raw_pair_rows = read_tsv(resolved["pairs"], compressed=True)
    missing_units = [
        field for field in HP_PS_UNIT_REQUIRED_FIELDS if field not in unit_fields
    ]
    missing_pairs = [
        field for field in HP_PS_PAIR_REQUIRED_FIELDS if field not in pair_fields
    ]
    if missing_units or missing_pairs:
        raise ReportInputError(
            "HP/PS detail fields differ: "
            f"units={missing_units}, pairs={missing_pairs}"
        )

    main_components = {
        (str(row["chrom"]), str(row["legacy_component_id"])): row
        for row in component_rows
    }
    main_keys = set(main_components)
    main_by_chrom: dict[str, set[str]] = {}
    for chrom, component_id in main_keys:
        main_by_chrom.setdefault(chrom, set()).add(component_id)
    normalized_units: list[dict[str, Any]] = []
    unit_keys: set[tuple[str, str, int, str]] = set()
    observed_components: set[tuple[str, str]] = set()
    hp_components = {1: set(), 2: set()}
    totals = {
        "total": 0,
        "retained": 0,
        "cut_lost": 0,
        "unavoidable": 0,
        "pattern_total": 0,
        "pattern_retained": 0,
        "pattern_cut_lost": 0,
        "pattern_unavoidable": 0,
    }
    for index, row in enumerate(raw_unit_rows, start=1):
        label = f"HP/PS unit[{index}]"
        if row["dataset"] != DATASET:
            raise ReportInputError(f"{label} dataset differs")
        chrom = row["chrom"]
        component_key = (chrom, row["legacy_component_id"])
        if component_key not in main_keys:
            raise ReportInputError(f"{label} component escaped report scope")
        hp = as_int(row["hp_family"], f"{label}.hp_family")
        if hp not in {1, 2}:
            raise ReportInputError(f"{label}.hp_family must be 1 or 2")
        unit_key = (chrom, row["legacy_component_id"], hp, row["phase_set"])
        if unit_key in unit_keys:
            raise ReportInputError(f"duplicate HP/PS unit: {unit_key}")
        unit_keys.add(unit_key)
        observed_components.add(component_key)
        hp_components[hp].add(component_key)

        component_k = as_int(row["component_k"], f"{label}.component_k")
        component_active_sites = as_int(
            row["component_primary_active_site_count"],
            f"{label}.component_primary_active_site_count",
        )
        unit_active_sites = as_int(
            row["unit_active_site_count"], f"{label}.unit_active_site_count"
        )
        unit_active_fraction = as_float(
            row["unit_active_site_fraction"], f"{label}.unit_active_site_fraction"
        )
        assert unit_active_fraction is not None
        if (
            component_k <= 8
            or not 0 <= unit_active_sites <= component_k
            or not 0 <= component_active_sites <= component_k
            or not math.isclose(
                unit_active_fraction,
                unit_active_sites / component_k,
                rel_tol=0,
                abs_tol=1e-12,
            )
        ):
            raise ReportInputError(f"{label} component/unit active-site fields differ")
        main_component = main_components[component_key]
        if "pre_cap_k" in main_component and as_int(
            main_component["pre_cap_k"], f"{label}.report_pre_cap_k"
        ) != component_k:
            raise ReportInputError(f"{label} component k differs from report")
        pattern_total = as_int(row["total_pattern_rows"], f"{label}.pattern_total")
        pattern_retained = as_int(
            row["retained_pattern_rows"], f"{label}.pattern_retained"
        )
        pattern_cut = as_int(row["cut_lost_pattern_rows"], f"{label}.pattern_cut")
        pattern_unavoidable = as_int(
            row["unavoidable_pattern_rows"], f"{label}.pattern_unavoidable"
        )
        pattern_nonretained = as_int(
            row["nonretained_pattern_rows"], f"{label}.pattern_nonretained"
        )
        weight_total = as_int(
            row["total_molecule_component_incidence_weight"], f"{label}.weight_total"
        )
        weight_retained = as_int(
            row["retained_molecule_component_incidence_weight"],
            f"{label}.weight_retained",
        )
        weight_cut = as_int(
            row["cut_lost_molecule_component_incidence_weight"],
            f"{label}.weight_cut",
        )
        weight_unavoidable = as_int(
            row["unavoidable_molecule_component_incidence_weight"],
            f"{label}.weight_unavoidable",
        )
        weight_nonretained = as_int(
            row["nonretained_molecule_component_incidence_weight"],
            f"{label}.weight_nonretained",
        )
        if (
            min(
                pattern_total,
                pattern_retained,
                pattern_cut,
                pattern_unavoidable,
                pattern_nonretained,
                weight_total,
                weight_retained,
                weight_cut,
                weight_unavoidable,
                weight_nonretained,
            )
            < 0
            or
            pattern_total
            != pattern_retained + pattern_cut + pattern_unavoidable
            or pattern_nonretained != pattern_cut + pattern_unavoidable
            or weight_total != weight_retained + weight_cut + weight_unavoidable
            or weight_nonretained != weight_cut + weight_unavoidable
            or weight_total <= 0
        ):
            raise ReportInputError(f"{label} incidence mass conservation failed")
        retention = as_float(row["retention_ratio"], f"{label}.retention_ratio")
        assert retention is not None
        if (
            not 0 <= retention <= 1
            or not math.isclose(
                retention, weight_retained / weight_total, rel_tol=0, abs_tol=1e-12
            )
            or row["ratio_status"] != "OBSERVED_CONSTRAINT_DENOMINATOR"
        ):
            raise ReportInputError(f"{label} observed retention denominator differs")
        eligible = as_bool(row["eligible_headline"], f"{label}.eligible_headline")
        if eligible != (weight_total >= 20 and pattern_total >= 5):
            raise ReportInputError(f"{label} eligibility differs")
        expected_stratum = (
            "1-4"
            if weight_total < 5
            else "5-19"
            if weight_total < 20
            else "20-49"
            if weight_total < 50
            else ">=50"
        )
        if row["support_stratum"] != expected_stratum:
            raise ReportInputError(f"{label} support stratum differs")
        normalized_units.append(
            {
                "dataset": DATASET,
                "chrom": chrom,
                "component_id": row["legacy_component_id"],
                "hp_family": f"HP{hp}",
                "phase_set": row["phase_set"],
                "unit_id": row["unit_id"],
                "component_k": component_k,
                "component_primary_active_site_count": component_active_sites,
                "unit_active_site_count": unit_active_sites,
                "unit_active_site_fraction": unit_active_fraction,
                "total_pattern_rows": pattern_total,
                "retained_pattern_rows": pattern_retained,
                "cut_lost_pattern_rows": pattern_cut,
                "unavoidable_pattern_rows": pattern_unavoidable,
                "nonretained_pattern_rows": pattern_nonretained,
                "total_weight": weight_total,
                "retained_weight": weight_retained,
                "cut_lost_weight": weight_cut,
                "unavoidable_weight": weight_unavoidable,
                "nonretained_weight": weight_nonretained,
                "retention_ratio": retention,
                "support_stratum": expected_stratum,
                "eligible_headline": eligible,
            }
        )
        for name, value in (
            ("total", weight_total),
            ("retained", weight_retained),
            ("cut_lost", weight_cut),
            ("unavoidable", weight_unavoidable),
            ("pattern_total", pattern_total),
            ("pattern_retained", pattern_retained),
            ("pattern_cut_lost", pattern_cut),
            ("pattern_unavoidable", pattern_unavoidable),
        ):
            totals[name] += value

    slice_specs = (
        (
            "by_chromosome",
            lambda unit: str(unit["chrom"]),
        ),
        (
            "by_hp_family",
            lambda unit: str(unit["hp_family"]),
        ),
        (
            "by_support_stratum",
            lambda unit: str(unit["support_stratum"]),
        ),
    )
    for summary_name, key_function in slice_specs:
        recorded_slices = summary.get(summary_name)
        if not isinstance(recorded_slices, Mapping):
            raise ReportInputError(f"HP/PS {summary_name} summaries are required")
        expected_slice_ids = {key_function(unit) for unit in normalized_units}
        if set(recorded_slices) != expected_slice_ids:
            raise ReportInputError(f"HP/PS {summary_name} scope differs")
        for slice_id in sorted(expected_slice_ids):
            validate_hp_ps_summary_slice(
                recorded_slices[slice_id],
                [
                    unit
                    for unit in normalized_units
                    if key_function(unit) == slice_id
                ],
                label=f"HP/PS {summary_name}.{slice_id}",
            )

    counts = summary.get("counts")
    coverage = summary.get("component_hp_unit_coverage")
    if not isinstance(counts, Mapping) or not isinstance(coverage, Mapping):
        raise ReportInputError("HP/PS v5 summary counts/component coverage are required")
    component_scope = len(main_keys)
    if (
        as_int(counts.get("components_in_partition_scope"), "HP/PS component scope")
        != component_scope
        or as_int(
            counts.get("components_with_observed_constraint_units"),
            "HP/PS observed components",
        )
        != len(observed_components)
        or as_int(
            counts.get("components_without_observed_constraint_units"),
            "HP/PS no-constraint components",
        )
        != component_scope - len(observed_components)
        or as_int(counts.get("observed_constraint_units"), "HP/PS observed units")
        != len(normalized_units)
    ):
        raise ReportInputError("HP/PS component/unit coverage counts differ")
    if not fixture and component_scope != EXPECTED_COMPONENTS:
        raise ReportInputError("formal HP/PS audit must cover 408 report components")

    coverage_expected = {
        "components_in_partition_scope": component_scope,
        "components_with_any_observed_unit": len(observed_components),
        "components_without_observed_unit": component_scope - len(observed_components),
        "components_hp1_only": len(hp_components[1] - hp_components[2]),
        "components_hp2_only": len(hp_components[2] - hp_components[1]),
        "components_hp1_and_hp2": len(hp_components[1] & hp_components[2]),
    }
    for name, observed in coverage_expected.items():
        if as_int(coverage.get(name), f"HP/PS coverage.{name}") != observed:
            raise ReportInputError(f"HP/PS v5 component coverage differs: {name}")
    overall_coverage_rows = {
        row["metric"]: row
        for row in summary_tsv_rows
        if row.get("scope_type") == "component_coverage"
        and row.get("scope_id") == "overall"
    }
    if set(overall_coverage_rows) != set(coverage_expected):
        raise ReportInputError("HP/PS summary.tsv overall component coverage differs")
    for name, observed in coverage_expected.items():
        if as_int(
            overall_coverage_rows[name]["value"],
            f"HP/PS summary.tsv coverage.{name}",
        ) != observed:
            raise ReportInputError(f"HP/PS summary.tsv coverage differs: {name}")
    coverage_by_chrom = coverage.get("by_chromosome")
    if not isinstance(coverage_by_chrom, Mapping):
        raise ReportInputError("HP/PS v5 per-chromosome coverage is required")
    for chrom, component_ids in main_by_chrom.items():
        chrom_row = coverage_by_chrom.get(chrom)
        if not isinstance(chrom_row, Mapping):
            raise ReportInputError(f"HP/PS coverage missing chromosome: {chrom}")
        chrom_hp1 = {key for key in hp_components[1] if key[0] == chrom}
        chrom_hp2 = {key for key in hp_components[2] if key[0] == chrom}
        chrom_observed = chrom_hp1 | chrom_hp2
        expected = {
            "components_in_partition_scope": len(component_ids),
            "components_with_any_observed_unit": len(chrom_observed),
            "components_without_observed_unit": len(component_ids)
            - len(chrom_observed),
            "components_hp1_only": len(chrom_hp1 - chrom_hp2),
            "components_hp2_only": len(chrom_hp2 - chrom_hp1),
            "components_hp1_and_hp2": len(chrom_hp1 & chrom_hp2),
        }
        for name, observed in expected.items():
            if as_int(chrom_row.get(name), f"HP/PS {chrom}.{name}") != observed:
                raise ReportInputError(f"HP/PS chromosome coverage differs: {chrom}.{name}")

    incidence = summary.get("molecule_component_incidence_totals")
    pattern_totals = summary.get("pattern_row_totals")
    if not isinstance(incidence, Mapping) or not isinstance(pattern_totals, Mapping):
        raise ReportInputError("HP/PS incidence totals are required")
    for summary_name, observed in (
        ("total", totals["total"]),
        ("retained", totals["retained"]),
        ("cut_lost", totals["cut_lost"]),
        ("unavoidable", totals["unavoidable"]),
    ):
        if as_int(incidence.get(summary_name), f"HP/PS incidence.{summary_name}") != observed:
            raise ReportInputError(f"HP/PS molecule incidence differs: {summary_name}")
    for summary_name, observed in (
        ("total", totals["pattern_total"]),
        ("retained", totals["pattern_retained"]),
        ("cut_lost", totals["pattern_cut_lost"]),
        ("unavoidable", totals["pattern_unavoidable"]),
    ):
        if (
            as_int(pattern_totals.get(summary_name), f"HP/PS patterns.{summary_name}")
            != observed
        ):
            raise ReportInputError(f"HP/PS pattern mass differs: {summary_name}")
    weighted_retention = totals["retained"] / totals["total"]
    recorded_weighted = as_float(
        incidence.get("weighted_retention_ratio"), "HP/PS weighted retention"
    )
    if recorded_weighted is None or not math.isclose(
        recorded_weighted, weighted_retention, rel_tol=0, abs_tol=1e-12
    ):
        raise ReportInputError("HP/PS weighted retention ratio differs")

    eligible_units = [row for row in normalized_units if row["eligible_headline"]]
    eligible_values = [float(row["retention_ratio"]) for row in eligible_units]
    if (
        as_int(counts.get("eligible_headline_units"), "HP/PS eligible units")
        != len(eligible_units)
        or not eligible_values
    ):
        raise ReportInputError("HP/PS eligible-unit count differs or is zero")
    distribution = summary.get("retention_distribution_eligible_headline_units")
    if not isinstance(distribution, Mapping):
        raise ReportInputError("HP/PS eligible retention distribution is required")
    distribution_counts = distribution.get("cumulative_threshold_counts")
    distribution_quantiles = distribution.get("quantiles")
    if not isinstance(distribution_counts, Mapping) or not isinstance(
        distribution_quantiles, Mapping
    ):
        raise ReportInputError("HP/PS eligible distribution fields are required")
    eligible_lt_05 = sum(value < 0.5 for value in eligible_values)
    eligible_lt_08 = sum(value < 0.8 for value in eligible_values)
    if (
        as_int(distribution.get("n_units"), "HP/PS eligible distribution n")
        != len(eligible_values)
        or as_int(distribution_counts.get("lt_0_5"), "HP/PS eligible lt0.5")
        != eligible_lt_05
        or as_int(distribution_counts.get("lt_0_8"), "HP/PS eligible lt0.8")
        != eligible_lt_08
    ):
        raise ReportInputError("HP/PS eligible tail counts differ")
    eligible_quantiles = {
        "p10": quantile_type7(eligible_values, 0.10),
        "p25": quantile_type7(eligible_values, 0.25),
        "median": quantile_type7(eligible_values, 0.50),
    }
    for name in ("p25", "median"):
        recorded = as_float(
            distribution_quantiles.get(name), f"HP/PS eligible {name}"
        )
        if recorded is None or not math.isclose(
            recorded, eligible_quantiles[name], rel_tol=0, abs_tol=1e-12
        ):
            raise ReportInputError(f"HP/PS eligible quantile differs: {name}")

    normalized_pairs: list[dict[str, Any]] = []
    pair_components: set[tuple[str, str]] = set()
    for index, row in enumerate(raw_pair_rows, start=1):
        label = f"HP/PS pair[{index}]"
        component_key = (row["chrom"], row["legacy_component_id"])
        if row["dataset"] != DATASET or component_key not in main_keys:
            raise ReportInputError(f"{label} scope differs")
        if component_key in pair_components:
            raise ReportInputError(f"duplicate HP/PS paired component: {component_key}")
        pair_components.add(component_key)
        hp1_total = as_int(
            row["hp1_total_molecule_component_incidence_weight"],
            f"{label}.hp1_total",
        )
        hp1_retained = as_int(
            row["hp1_retained_molecule_component_incidence_weight"],
            f"{label}.hp1_retained",
        )
        hp2_total = as_int(
            row["hp2_total_molecule_component_incidence_weight"],
            f"{label}.hp2_total",
        )
        hp2_retained = as_int(
            row["hp2_retained_molecule_component_incidence_weight"],
            f"{label}.hp2_retained",
        )
        hp1_ratio = as_float(row["hp1_retention_ratio"], f"{label}.hp1_ratio")
        hp2_ratio = as_float(row["hp2_retention_ratio"], f"{label}.hp2_ratio")
        assert hp1_ratio is not None and hp2_ratio is not None
        if (
            hp1_total <= 0
            or hp2_total <= 0
            or not math.isclose(hp1_ratio, hp1_retained / hp1_total, abs_tol=1e-12)
            or not math.isclose(hp2_ratio, hp2_retained / hp2_total, abs_tol=1e-12)
        ):
            raise ReportInputError(f"{label} retention ratio differs")
        delta = hp1_ratio - hp2_ratio
        absolute_delta = abs(delta)
        recorded_delta = as_float(
            row["hp1_minus_hp2_retention_delta"], f"{label}.delta"
        )
        recorded_absolute = as_float(
            row["absolute_retention_delta"], f"{label}.absolute_delta"
        )
        if (
            recorded_delta is None
            or recorded_absolute is None
            or not math.isclose(recorded_delta, delta, abs_tol=1e-12)
            or not math.isclose(recorded_absolute, absolute_delta, abs_tol=1e-12)
        ):
            raise ReportInputError(f"{label} delta differs")
        both_eligible = as_bool(
            row["both_hp_headline_eligible"], f"{label}.both_eligible"
        )
        expected_eligible = (
            hp1_total >= 20
            and as_int(row["hp1_total_pattern_rows"], f"{label}.hp1_patterns") >= 5
            and hp2_total >= 20
            and as_int(row["hp2_total_pattern_rows"], f"{label}.hp2_patterns") >= 5
        )
        if both_eligible != expected_eligible:
            raise ReportInputError(f"{label} paired eligibility differs")
        normalized_pairs.append(
            {
                "chrom": row["chrom"],
                "component_id": row["legacy_component_id"],
                "hp1_retention_ratio": hp1_ratio,
                "hp2_retention_ratio": hp2_ratio,
                "hp1_minus_hp2_retention_delta": delta,
                "absolute_retention_delta": absolute_delta,
                "both_hp_headline_eligible": both_eligible,
            }
        )
    if pair_components != hp_components[1] & hp_components[2]:
        raise ReportInputError("HP/PS paired-component coverage differs")
    pair_summary = summary.get("hp1_hp2_component_paired")
    if not isinstance(pair_summary, Mapping):
        raise ReportInputError("HP/PS paired-component summary is required")
    eligible_pairs = [
        row for row in normalized_pairs if row["both_hp_headline_eligible"]
    ]
    pair_deltas = [
        float(row["absolute_retention_delta"]) for row in normalized_pairs
    ]
    if (
        as_int(pair_summary.get("n_pairs"), "HP/PS paired components")
        != len(normalized_pairs)
        or as_int(
            pair_summary.get("n_both_headline_eligible"),
            "HP/PS eligible paired components",
        )
        != len(eligible_pairs)
        or not pair_deltas
    ):
        raise ReportInputError("HP/PS paired-component counts differ")
    pair_delta_quantiles = pair_summary.get("absolute_delta_quantiles")
    if not isinstance(pair_delta_quantiles, Mapping):
        raise ReportInputError("HP/PS pair-delta quantiles are required")
    pair_delta_median = quantile_type7(pair_deltas, 0.5)
    pair_delta_ge_025 = sum(value >= 0.25 for value in pair_deltas)
    recorded_pair_median = as_float(
        pair_delta_quantiles.get("median"), "HP/PS pair delta median"
    )
    if (
        recorded_pair_median is None
        or not math.isclose(
            recorded_pair_median,
            pair_delta_median,
            rel_tol=0,
            abs_tol=1e-12,
        )
        or as_int(
            pair_summary.get("absolute_delta_ge_0_25"),
            "HP/PS pair delta >=0.25",
        )
        != pair_delta_ge_025
    ):
        raise ReportInputError("HP/PS pair-delta summary differs")

    hp_weighted: dict[str, float] = {}
    for hp_label in ("HP1", "HP2"):
        hp_rows = [row for row in normalized_units if row["hp_family"] == hp_label]
        hp_total = sum(row["total_weight"] for row in hp_rows)
        hp_retained = sum(row["retained_weight"] for row in hp_rows)
        hp_weighted[hp_label] = hp_retained / hp_total
        recorded_hp = as_float(
            summary.get("by_hp_family", {})
            .get(hp_label, {})
            .get("molecule_component_incidences", {})
            .get("weighted_retention_ratio"),
            f"HP/PS {hp_label} weighted retention",
        )
        if recorded_hp is None or not math.isclose(
            recorded_hp, hp_weighted[hp_label], abs_tol=1e-12
        ):
            raise ReportInputError(f"HP/PS {hp_label} weighted retention differs")

    return {
        "paths": resolved,
        "summary": summary,
        "units": normalized_units,
        "pairs": normalized_pairs,
        "metrics": {
            "components_in_partition_scope": component_scope,
            "components_with_observed_constraint_units": len(observed_components),
            "components_without_observed_constraint_units": component_scope
            - len(observed_components),
            "components_hp1_only": coverage_expected["components_hp1_only"],
            "components_hp2_only": coverage_expected["components_hp2_only"],
            "components_hp1_and_hp2": coverage_expected["components_hp1_and_hp2"],
            "observed_constraint_units": len(normalized_units),
            "eligible_headline_units": len(eligible_units),
            "weighted_retention_ratio": weighted_retention,
            "eligible_retention_p10": eligible_quantiles["p10"],
            "eligible_retention_p25": eligible_quantiles["p25"],
            "eligible_retention_median": eligible_quantiles["median"],
            "eligible_retention_lt_0_5": eligible_lt_05,
            "eligible_retention_lt_0_8": eligible_lt_08,
            "hp1_weighted_retention_ratio": hp_weighted["HP1"],
            "hp2_weighted_retention_ratio": hp_weighted["HP2"],
            "hp1_minus_hp2_weighted_retention_delta": hp_weighted["HP1"]
            - hp_weighted["HP2"],
            "paired_components": len(normalized_pairs),
            "paired_components_both_eligible": len(eligible_pairs),
            "paired_component_absolute_delta_median": pair_delta_median,
            "paired_component_absolute_delta_ge_0_25": pair_delta_ge_025,
        },
    }


def validate_span_grid(
    candidate: Path,
    *,
    fixture: bool,
    main_totals: Mapping[str, Any],
) -> dict[str, Any]:
    """Authenticate and aggregate the cached 50/100/200 kb sensitivity grid."""
    summary_path, receipt_path = resolve_span_grid_paths(candidate)
    receipt = load_json(receipt_path)
    if receipt.get("schema_name") != SPAN_GRID_SCHEMA:
        raise ReportInputError(
            f"unexpected span-grid receipt schema: {receipt.get('schema_name')!r}"
        )
    if receipt.get("schema_version") != SPAN_GRID_SCHEMA_VERSION:
        raise ReportInputError(
            f"unexpected span-grid schema version: {receipt.get('schema_version')!r}"
        )
    if receipt.get("sample") != DATASET or receipt.get("all_pass") is not True:
        raise ReportInputError("span-grid receipt must be all-pass HCC1395")

    scope = receipt.get("scope")
    if not isinstance(scope, Mapping):
        raise ReportInputError("span-grid receipt.scope must be an object")
    scope_chroms = tuple(scope.get("chromosomes", ()))
    scope_caps = tuple(scope.get("span_caps_bp", ()))
    if scope_caps != CANONICAL_SPAN_CAPS:
        raise ReportInputError(
            "span-grid caps must be exactly 50000,100000,200000 in canonical order"
        )
    if not scope_chroms or len(scope_chroms) != len(set(scope_chroms)):
        raise ReportInputError("span-grid chromosome scope must be nonempty and unique")
    for chrom in scope_chroms:
        chrom_number(str(chrom))
    if not fixture:
        if (
            receipt.get("comprehensive_all_pass") is not True
            or scope_chroms != AUTOSOMES
            or scope.get("test_mode") is not False
        ):
            raise ReportInputError(
                "formal span-grid requires comprehensive chr1-chr22 production receipt"
            )

    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping) or "summary" not in outputs:
        raise ReportInputError("span-grid receipt does not bind summary.tsv")
    verify_file_identity(
        outputs["summary"],
        expected_path=summary_path,
        label="span-grid summary.tsv",
    )

    fields, rows = read_tsv(summary_path)
    missing = [field for field in SPAN_GRID_REQUIRED_FIELDS if field not in fields]
    if missing:
        raise ReportInputError(
            "span-grid summary.tsv is missing fields: " + ", ".join(missing)
        )
    expected_pairs = [
        (cap, str(chrom)) for cap in CANONICAL_SPAN_CAPS for chrom in scope_chroms
    ]
    observed_pairs = [
        (
            as_int(row["span_cap_bp"], f"span-grid row[{index}].span_cap_bp"),
            row["chrom"],
        )
        for index, row in enumerate(rows, start=1)
    ]
    if observed_pairs != expected_pairs:
        raise ReportInputError(
            "span-grid summary row order/scope must be cap-major canonical grid"
        )

    additive_fields = (
        "ssnv_sites",
        "k_gt8_components",
        "k_gt8_sites",
        "new_blocks",
        "exact_patterns",
        "raw_total_molecule_weight",
        "raw_retained_molecule_weight",
        "raw_lost_molecule_weight",
        "unavoidable_patterns",
        "unavoidable_size_patterns",
        "unavoidable_span_cap_patterns",
        "unavoidable_both_limits_patterns",
    )
    aggregate_by_cap: dict[int, dict[str, Any]] = {}
    for cap in CANONICAL_SPAN_CAPS:
        cap_rows = [row for row in rows if as_int(row["span_cap_bp"], "span_cap_bp") == cap]
        aggregate: dict[str, Any] = {
            "span_cap_bp": cap,
            "chromosomes": len(cap_rows),
            "completed_partition_chromosomes": 0,
            "zero_target_skipped_chromosomes": 0,
            "cached_partition_wall_seconds": 0.0,
        }
        for field in additive_fields:
            aggregate[field] = 0
        for index, row in enumerate(cap_rows, start=1):
            label = f"span-grid {cap} row[{index}]"
            chrom_number(row["chrom"])
            status = row["status"]
            if status not in {"COMPLETED", "SKIP_NO_K_GT8_TARGET"}:
                raise ReportInputError(f"{label} has unsupported status: {status!r}")
            wall_seconds = as_float(row["wall_seconds"], f"{label}.wall_seconds")
            assert wall_seconds is not None
            if wall_seconds < 0:
                raise ReportInputError(f"{label}.wall_seconds must be nonnegative")
            if status == "COMPLETED":
                aggregate["completed_partition_chromosomes"] += 1
                aggregate["cached_partition_wall_seconds"] += wall_seconds
            else:
                aggregate["zero_target_skipped_chromosomes"] += 1
            for field in additive_fields:
                value = as_int(row[field], f"{label}.{field}")
                if value < 0:
                    raise ReportInputError(f"{label}.{field} must be nonnegative")
                aggregate[field] += value
            retained = as_int(
                row["raw_retained_molecule_weight"], f"{label}.retained"
            )
            lost = as_int(row["raw_lost_molecule_weight"], f"{label}.lost")
            total = as_int(row["raw_total_molecule_weight"], f"{label}.total")
            if retained + lost != total:
                raise ReportInputError(f"{label} molecule-weight conservation failed")
            if status == "SKIP_NO_K_GT8_TARGET" and any(
                as_int(row[field], f"{label}.{field}") != 0
                for field in (
                    "k_gt8_components",
                    "k_gt8_sites",
                    "new_blocks",
                    "raw_total_molecule_weight",
                )
            ):
                raise ReportInputError(f"{label} zero-target skip contains target metrics")
        aggregate_by_cap[cap] = aggregate

    cap_receipts = receipt.get("caps")
    if not isinstance(cap_receipts, list) or [
        item.get("span_cap_bp") if isinstance(item, Mapping) else None
        for item in cap_receipts
    ] != list(CANONICAL_SPAN_CAPS):
        raise ReportInputError("span-grid receipt.caps must follow canonical cap order")
    receipt_by_cap = {
        as_int(item["span_cap_bp"], "receipt cap"): item for item in cap_receipts
    }
    for cap, aggregate in aggregate_by_cap.items():
        item = receipt_by_cap[cap]
        totals = item.get("totals")
        if not isinstance(totals, Mapping):
            raise ReportInputError(f"span-grid cap {cap} totals must be an object")
        for field, observed in aggregate.items():
            if field in {"span_cap_bp", "ssnv_sites"}:
                continue
            if field == "cached_partition_wall_seconds":
                expected = as_float(totals.get(field), f"cap {cap}.{field}")
                if expected is None or not math.isclose(
                    observed, expected, rel_tol=0, abs_tol=1e-9
                ):
                    raise ReportInputError(f"span-grid cap {cap} runtime total mismatch")
            elif as_int(totals.get(field), f"cap {cap}.{field}") != observed:
                raise ReportInputError(
                    f"span-grid cap {cap} receipt mismatch for {field}"
                )

    main_components = as_int(main_totals["components"], "ALL.components")
    main_sites = as_int(main_totals["sites"], "ALL.sites")
    main_raw_weight = as_int(
        main_totals["raw_total_molecule_weight"], "ALL.raw_total_molecule_weight"
    )
    for cap, aggregate in aggregate_by_cap.items():
        if (
            aggregate["k_gt8_components"] != main_components
            or aggregate["k_gt8_sites"] != main_sites
            or aggregate["raw_total_molecule_weight"] != main_raw_weight
        ):
            raise ReportInputError(
                f"span-grid cap {cap} does not share the k-only target/denominator"
            )

    run_totals = receipt.get("totals")
    if not isinstance(run_totals, Mapping):
        raise ReportInputError("span-grid run totals must be an object")
    expected_run_totals = {
        "span_caps": len(CANONICAL_SPAN_CAPS),
        "chromosome_cap_tasks": len(rows),
        "completed_partition_tasks": sum(
            aggregate["completed_partition_chromosomes"]
            for aggregate in aggregate_by_cap.values()
        ),
        "zero_target_skipped_tasks": sum(
            aggregate["zero_target_skipped_chromosomes"]
            for aggregate in aggregate_by_cap.values()
        ),
    }
    for field, observed in expected_run_totals.items():
        if as_int(run_totals.get(field), f"span-grid totals.{field}") != observed:
            raise ReportInputError(f"span-grid run total mismatch for {field}")
    observed_wall = sum(
        aggregate["cached_partition_wall_seconds"]
        for aggregate in aggregate_by_cap.values()
    )
    expected_wall = as_float(
        run_totals.get("cached_partition_wall_seconds"),
        "span-grid totals.cached_partition_wall_seconds",
    )
    if expected_wall is None or not math.isclose(
        observed_wall, expected_wall, rel_tol=0, abs_tol=1e-9
    ):
        raise ReportInputError("span-grid aggregate runtime total mismatch")
    expected_timing_scope = (
        "cached_partition_wall_seconds excludes source identity "
        "validation and excludes upstream BAM extraction"
    )
    if receipt.get("timing_scope") != expected_timing_scope:
        raise ReportInputError("span-grid timing_scope differs from runner contract")

    return {
        "summary_path": summary_path,
        "receipt_path": receipt_path,
        "receipt": receipt,
        "aggregates": [aggregate_by_cap[cap] for cap in CANONICAL_SPAN_CAPS],
    }


def validate_summary_row(row: Mapping[str, Any], label: str) -> None:
    missing = [field for field in SUMMARY_REQUIRED_FIELDS if field not in row]
    if missing:
        raise ReportInputError(f"{label} is missing fields: {', '.join(missing)}")

    nonnegative_integer_fields = (
        "ssnv_sites",
        "components",
        "sites",
        "old_selected_sites",
        "old_excluded_sites",
        "new_blocks",
        "new_retained_sites",
        "primary_active_sites",
        "raw_total_molecule_weight",
        "new_retained_molecule_weight",
        "new_lost_molecule_weight",
        "unavoidable_patterns",
        "unavoidable_molecule_weight",
        "unavoidable_n_fixed_ra_gt8_patterns",
        "unavoidable_n_fixed_ra_gt8_molecule_weight",
        "unavoidable_n_fixed_ra_lte8_span_gt8_patterns",
        "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight",
        "unavoidable_patterns",
        "unavoidable_molecule_weight",
        "unavoidable_n_fixed_ra_gt8_patterns",
        "unavoidable_n_fixed_ra_gt8_molecule_weight",
        "unavoidable_n_fixed_ra_lte8_span_gt8_patterns",
        "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight",
        "old_densest8_retained_molecule_weight",
        "weight_stable_components",
        "weight_sensitive_components",
        "zero_evidence_blocks",
        "zero_evidence_block_sites",
        "tree_ready_blocks",
        "tree_ready_block_sites",
        "abstain_blocks",
        "abstain_block_sites",
    )
    values: dict[str, int] = {}
    for field in nonnegative_integer_fields:
        values[field] = as_int(row[field], f"{label}.{field}")
        if values[field] < 0:
            raise ReportInputError(f"{label}.{field} must be nonnegative")

    if values["old_selected_sites"] + values["old_excluded_sites"] != values["sites"]:
        raise ReportInputError(f"{label} old selected/excluded site conservation failed")
    if values["new_retained_sites"] != values["sites"]:
        raise ReportInputError(f"{label} new retained-site conservation failed")
    if (
        values["new_retained_molecule_weight"] + values["new_lost_molecule_weight"]
        != values["raw_total_molecule_weight"]
    ):
        raise ReportInputError(f"{label} molecule-weight conservation failed")
    if (
        values["unavoidable_n_fixed_ra_gt8_patterns"]
        + values["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"]
        != values["unavoidable_patterns"]
    ):
        raise ReportInputError(
            f"{label} unavoidable pattern-mechanism conservation failed"
        )
    if (
        values["unavoidable_n_fixed_ra_gt8_molecule_weight"]
        + values["unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight"]
        != values["unavoidable_molecule_weight"]
    ):
        raise ReportInputError(
            f"{label} unavoidable weight-mechanism conservation failed"
        )
    if values["unavoidable_molecule_weight"] > values["new_lost_molecule_weight"]:
        raise ReportInputError(f"{label} unavoidable weight exceeds all lost weight")
    if (
        values["weight_stable_components"] + values["weight_sensitive_components"]
        != values["components"]
    ):
        raise ReportInputError(f"{label} weight-stability component conservation failed")
    if values["primary_active_sites"] > values["sites"]:
        raise ReportInputError(f"{label} primary-active sites exceed target sites")
    if values["zero_evidence_blocks"] > values["new_blocks"]:
        raise ReportInputError(f"{label} zero-evidence blocks exceed new blocks")
    if values["zero_evidence_block_sites"] > values["sites"]:
        raise ReportInputError(f"{label} zero-evidence block sites exceed target sites")
    if values["zero_evidence_blocks"] > values["abstain_blocks"]:
        raise ReportInputError(f"{label} zero-evidence blocks exceed ABSTAIN blocks")
    if values["zero_evidence_block_sites"] > values["abstain_block_sites"]:
        raise ReportInputError(f"{label} zero-evidence sites exceed ABSTAIN sites")
    if values["tree_ready_blocks"] + values["abstain_blocks"] != values["new_blocks"]:
        raise ReportInputError(f"{label} tree-ready/ABSTAIN block conservation failed")
    if (
        values["tree_ready_block_sites"] + values["abstain_block_sites"]
        != values["sites"]
    ):
        raise ReportInputError(f"{label} tree-ready/ABSTAIN site conservation failed")

    new_ratio = as_float(
        row["new_retention_ratio"],
        f"{label}.new_retention_ratio",
        blank_ok=True,
    )
    old_ratio = as_float(
        row["old_densest8_retention_ratio"],
        f"{label}.old_densest8_retention_ratio",
        blank_ok=True,
    )
    expected_new = ratio(
        values["new_retained_molecule_weight"], values["raw_total_molecule_weight"]
    )
    expected_old = ratio(
        values["old_densest8_retained_molecule_weight"],
        values["raw_total_molecule_weight"],
    )
    if expected_new is None:
        if new_ratio not in {None, 0}:
            raise ReportInputError(
                f"{label} zero-denominator new ratio must be blank or zero"
            )
    elif new_ratio is None:
        raise ReportInputError(f"{label} nonzero-denominator new ratio is blank")
    elif not math.isclose(new_ratio, expected_new, rel_tol=0, abs_tol=1e-12):
        raise ReportInputError(f"{label} new retention ratio does not match weights")
    if expected_old is None:
        if old_ratio not in {None, 0}:
            raise ReportInputError(
                f"{label} zero-denominator old ratio must be blank or zero"
            )
    elif old_ratio is None:
        raise ReportInputError(f"{label} nonzero-denominator old ratio is blank")
    elif not math.isclose(old_ratio, expected_old, rel_tol=0, abs_tol=1e-12):
        raise ReportInputError(f"{label} old retention ratio does not match weights")

    timing_fields = (
        "extraction_wall_seconds",
        "partition_wall_seconds",
        "partition_pattern_load_aggregate_seconds",
        "partition_ordered_hypergraph_dp_seconds",
    )
    timing_values: dict[str, float | None] = {}
    for field in timing_fields:
        value = as_float(row[field], f"{label}.{field}", blank_ok=True)
        if value is not None and value < 0:
            raise ReportInputError(f"{label}.{field} must be nonnegative")
        timing_values[field] = value
    pattern_seconds = timing_values["partition_pattern_load_aggregate_seconds"]
    component_loop_seconds = timing_values[
        "partition_ordered_hypergraph_dp_seconds"
    ]
    partition_seconds = timing_values["partition_wall_seconds"]
    if pattern_seconds is None or component_loop_seconds is None:
        raise ReportInputError(f"{label} partition internal timings must be explicit")
    if partition_seconds is None:
        if pattern_seconds != 0 or component_loop_seconds != 0:
            raise ReportInputError(
                f"{label} skipped partition must have zero internal timings"
            )
    elif pattern_seconds + component_loop_seconds > partition_seconds + 0.1:
        raise ReportInputError(
            f"{label} partition internal timings exceed partition stage wall"
        )


def validate_inputs(
    summary: Mapping[str, Any],
    summary_tsv_rows: Sequence[Mapping[str, str]],
    summary_tsv_fields: Sequence[str],
    component_rows: Sequence[Mapping[str, str]],
    component_fields: Sequence[str],
    *,
    fixture: bool,
) -> tuple[list[Mapping[str, Any]], Mapping[str, Any]]:
    if summary.get("schema_name") != SUMMARY_SCHEMA:
        raise ReportInputError(f"unexpected summary schema: {summary.get('schema_name')!r}")
    if summary.get("schema_version") != SUMMARY_SCHEMA_VERSION:
        raise ReportInputError(
            f"unexpected summary schema version: {summary.get('schema_version')!r}"
        )
    if summary.get("sample") != DATASET:
        raise ReportInputError(f"summary sample must be {DATASET}")

    per_chrom = summary.get("per_chromosome")
    totals = summary.get("totals")
    if not isinstance(per_chrom, list) or not per_chrom:
        raise ReportInputError("summary.per_chromosome must be a nonempty array")
    if not isinstance(totals, dict):
        raise ReportInputError("summary.totals must be an object")
    if any(not isinstance(row, dict) for row in per_chrom):
        raise ReportInputError("summary.per_chromosome rows must be objects")

    chroms = [str(row.get("chrom", "")) for row in per_chrom]
    if len(set(chroms)) != len(chroms):
        raise ReportInputError("summary.per_chromosome contains duplicate chromosomes")
    for row in per_chrom:
        chrom_number(str(row["chrom"]))
        validate_summary_row(row, str(row["chrom"]))
    if totals.get("chrom") != "ALL":
        raise ReportInputError("summary.totals.chrom must be ALL")
    validate_summary_row(totals, "ALL")

    for field in SUMMARY_REQUIRED_FIELDS:
        if field not in summary_tsv_fields:
            raise ReportInputError(f"summary.tsv is missing required field: {field}")
    tsv_by_chrom = {row.get("chrom", ""): row for row in summary_tsv_rows}
    expected_tsv_keys = {*chroms, "ALL"}
    if set(tsv_by_chrom) != expected_tsv_keys or len(summary_tsv_rows) != len(
        expected_tsv_keys
    ):
        raise ReportInputError("summary.tsv chromosome rows do not match summary.json")
    for json_row in [*per_chrom, totals]:
        chrom = str(json_row["chrom"])
        tsv_row = tsv_by_chrom[chrom]
        for field in SUMMARY_TSV_CROSSCHECK_FIELDS:
            if not numeric_equal(json_row[field], tsv_row[field], f"{chrom}.{field}"):
                raise ReportInputError(
                    f"summary.json/summary.tsv mismatch: {chrom}.{field}"
                )

    additive_fields = (
        "ssnv_sites",
        "components",
        "sites",
        "old_selected_sites",
        "old_excluded_sites",
        "new_blocks",
        "new_retained_sites",
        "primary_active_sites",
        "raw_total_molecule_weight",
        "new_retained_molecule_weight",
        "new_lost_molecule_weight",
        "old_densest8_retained_molecule_weight",
        "weight_stable_components",
        "weight_sensitive_components",
        "zero_evidence_blocks",
        "zero_evidence_block_sites",
        "tree_ready_blocks",
        "tree_ready_block_sites",
        "abstain_blocks",
        "abstain_block_sites",
    )
    for field in additive_fields:
        observed = sum(as_int(row[field], f"{row['chrom']}.{field}") for row in per_chrom)
        expected = as_int(totals[field], f"ALL.{field}")
        if observed != expected:
            raise ReportInputError(f"per-chromosome sum mismatch for {field}")
    for field in (
        "partition_pattern_load_aggregate_seconds",
        "partition_ordered_hypergraph_dp_seconds",
    ):
        observed = sum(
            as_float(row[field], f"{row['chrom']}.{field}") or 0.0
            for row in per_chrom
        )
        expected = as_float(totals[field], f"ALL.{field}")
        if expected is None or not math.isclose(
            observed, expected, rel_tol=0, abs_tol=1e-9
        ):
            raise ReportInputError(f"per-chromosome sum mismatch for {field}")

    missing_component_fields = [
        field for field in COMPONENT_REQUIRED_FIELDS if field not in component_fields
    ]
    if missing_component_fields:
        raise ReportInputError(
            "component_all.tsv.gz is missing fields: "
            + ", ".join(missing_component_fields)
        )
    if len(component_rows) != as_int(totals["components"], "ALL.components"):
        raise ReportInputError("component row count does not match summary total")
    component_chroms = {row["chrom"] for row in component_rows}
    if not component_chroms.issubset(set(chroms)):
        raise ReportInputError("component rows contain chromosomes absent from summary")
    if any(row["dataset"] != DATASET for row in component_rows):
        raise ReportInputError(f"all component rows must have dataset={DATASET}")

    component_sums = {
        "new_site_retained": 0,
        "old_densest8_selected": 0,
        "old_cap_excluded": 0,
        "primary_active_site_count": 0,
        "raw_total_molecule_weight": 0,
        "raw_retained_molecule_weight": 0,
        "raw_lost_molecule_weight": 0,
        "old_densest8_retained_molecule_weight": 0,
    }
    for index, row in enumerate(component_rows, start=1):
        label = f"component[{index}]"
        chrom_number(row["chrom"])
        if as_int(row["pre_cap_k"], f"{label}.pre_cap_k") <= 8:
            raise ReportInputError(f"{label}.pre_cap_k must be greater than 8")
        for field in component_sums:
            value = as_int(row[field], f"{label}.{field}")
            if value < 0:
                raise ReportInputError(f"{label}.{field} must be nonnegative")
            component_sums[field] += value
        if (
            as_int(row["raw_retained_molecule_weight"], f"{label}.retained")
            + as_int(row["raw_lost_molecule_weight"], f"{label}.lost")
            != as_int(row["raw_total_molecule_weight"], f"{label}.total")
        ):
            raise ReportInputError(f"{label} molecule-weight conservation failed")
        as_bool(row["weight_stable"], f"{label}.weight_stable")
        primary_active = as_int(
            row["primary_active_site_count"], f"{label}.primary_active_site_count"
        )
        component_sites = as_int(
            row["new_site_retained"], f"{label}.new_site_retained"
        )
        if primary_active > component_sites:
            raise ReportInputError(f"{label} primary-active sites exceed component sites")
        primary_fraction = as_float(
            row["primary_active_site_fraction"],
            f"{label}.primary_active_site_fraction",
        )
        expected_fraction = ratio(primary_active, component_sites)
        if (
            primary_fraction is None
            or expected_fraction is None
            or not math.isclose(
                primary_fraction, expected_fraction, rel_tol=0, abs_tol=1e-12
            )
        ):
            raise ReportInputError(
                f"{label} primary-active fraction does not match site counts"
            )
        if not row["status"]:
            raise ReportInputError(f"{label}.status must be nonempty")
        if not re.fullmatch(r"[0-9a-fA-F]{64}", row["positions_sha256"]):
            raise ReportInputError(f"{label}.positions_sha256 is not SHA-256")

    expected_component_sums = {
        "new_site_retained": "new_retained_sites",
        "old_densest8_selected": "old_selected_sites",
        "old_cap_excluded": "old_excluded_sites",
        "primary_active_site_count": "primary_active_sites",
        "raw_total_molecule_weight": "raw_total_molecule_weight",
        "raw_retained_molecule_weight": "new_retained_molecule_weight",
        "raw_lost_molecule_weight": "new_lost_molecule_weight",
        "old_densest8_retained_molecule_weight": (
            "old_densest8_retained_molecule_weight"
        ),
    }
    for component_field, summary_field in expected_component_sums.items():
        if component_sums[component_field] != as_int(
            totals[summary_field], f"ALL.{summary_field}"
        ):
            raise ReportInputError(
                f"component sum mismatch: {component_field} vs {summary_field}"
            )

    if not fixture:
        if summary.get("all_pass") is not True or summary.get(
            "comprehensive_all_pass"
        ) is not True:
            raise ReportInputError(
                "formal report requires all_pass=true and comprehensive_all_pass=true"
            )
        if tuple(chroms) != AUTOSOMES:
            raise ReportInputError("formal report requires ordered chr1-chr22 scope")
        scope = summary.get("scope")
        if not isinstance(scope, dict) or tuple(scope.get("chromosomes", ())) != AUTOSOMES:
            raise ReportInputError("formal report scope.chromosomes must be chr1-chr22")
        expected = {
            "components": EXPECTED_COMPONENTS,
            "sites": EXPECTED_TARGET_SITES,
            "ssnv_sites": EXPECTED_SSNV_SITES,
        }
        for field, value in expected.items():
            if as_int(totals[field], f"ALL.{field}") != value:
                raise ReportInputError(
                    f"formal report canonical total mismatch: {field} != {value}"
                )
        checks = summary.get("checks")
        if not isinstance(checks, dict) or not checks or any(
            value is not True for value in checks.values()
        ):
            raise ReportInputError("formal report requires every summarizer check to pass")
    return list(per_chrom), totals


def fmt_integer(value: int | float) -> str:
    return f"{int(value):,}"


def fmt_percent(value: float | None) -> str:
    return "無分母" if value is None else f"{value * 100:.2f}%"


def fmt_minutes(value: float | None) -> str:
    return "未提供" if value is None else f"{value:,.2f} 分鐘"


def source_spec(
    *,
    source_id: str,
    label: str,
    path: str | None,
    sql: str,
    description: str,
    tables: Sequence[str],
    filters: Sequence[str],
    metric_definitions: Sequence[str],
    generated_at: str,
) -> dict[str, Any]:
    source: dict[str, Any] = {
        "id": source_id,
        "label": label,
        "query": {
            "engine": "SQLite",
            "language": "sql",
            "sql": sql,
            "description": description,
            "executed_at": generated_at,
            "tables_used": list(tables),
            "filters": list(filters),
            "metric_definitions": list(metric_definitions),
        },
    }
    if path is not None:
        source["path"] = path
    return source


def build_artifact(
    *,
    summary: Mapping[str, Any],
    per_chrom: Sequence[Mapping[str, Any]],
    totals: Mapping[str, Any],
    component_rows: Sequence[Mapping[str, str]],
    baseline: Mapping[str, Any],
    span_grid: Mapping[str, Any] | None,
    exact_log_audit: Mapping[str, Any] | None,
    hp_ps_unit_audit: Mapping[str, Any] | None,
    input_paths: Mapping[str, Path],
    generated_at: str,
    fixture: bool,
    top_components: int,
) -> dict[str, Any]:
    per_chrom_dataset: list[dict[str, Any]] = []
    retention_dataset: list[dict[str, Any]] = []
    runtime_dataset: list[dict[str, Any]] = []
    exact_by_component = (
        exact_log_audit.get("component_rows", {})
        if exact_log_audit is not None
        else {}
    )
    exact_by_chrom = (
        exact_log_audit.get("chromosome_counts", {})
        if exact_log_audit is not None
        else {}
    )
    abstain_components_by_chrom = {chrom: 0 for chrom in AUTOSOMES}
    evidence_components_by_chrom = {chrom: 0 for chrom in AUTOSOMES}
    for component in component_rows:
        chrom = component["chrom"]
        status = component["status"]
        retained_weight = as_int(
            component["raw_retained_molecule_weight"],
            f"{chrom}.{component['legacy_component_id']}.raw_retained_molecule_weight",
        )
        if status.startswith("ABSTAIN"):
            abstain_components_by_chrom[chrom] = (
                abstain_components_by_chrom.get(chrom, 0) + 1
            )
        if retained_weight > 0 and not status.startswith("ABSTAIN"):
            evidence_components_by_chrom[chrom] = (
                evidence_components_by_chrom.get(chrom, 0) + 1
            )

    for raw in sorted(per_chrom, key=lambda row: chrom_number(str(row["chrom"]))):
        chrom = str(raw["chrom"])
        order = chrom_number(chrom)
        sites = as_int(raw["sites"], f"{chrom}.sites")
        raw_weight = as_int(
            raw["raw_total_molecule_weight"], f"{chrom}.raw_total_molecule_weight"
        )
        new_weight = as_int(
            raw["new_retained_molecule_weight"],
            f"{chrom}.new_retained_molecule_weight",
        )
        old_weight = as_int(
            raw["old_densest8_retained_molecule_weight"],
            f"{chrom}.old_densest8_retained_molecule_weight",
        )
        extraction_seconds = as_float(
            raw["extraction_wall_seconds"], f"{chrom}.extraction_wall_seconds"
        )
        partition_seconds = as_float(
            raw["partition_wall_seconds"],
            f"{chrom}.partition_wall_seconds",
            blank_ok=True,
        )
        pattern_load_seconds = as_float(
            raw["partition_pattern_load_aggregate_seconds"],
            f"{chrom}.partition_pattern_load_aggregate_seconds",
        )
        component_loop_seconds = as_float(
            raw["partition_ordered_hypergraph_dp_seconds"],
            f"{chrom}.partition_ordered_hypergraph_dp_seconds",
        )
        assert extraction_seconds is not None
        assert (
            pattern_load_seconds is not None
            and component_loop_seconds is not None
        )
        new_read_ratio = ratio(new_weight, raw_weight)
        old_read_ratio = ratio(old_weight, raw_weight)
        new_blocks = as_int(raw["new_blocks"], f"{chrom}.new_blocks")
        zero_evidence_blocks = as_int(
            raw["zero_evidence_blocks"], f"{chrom}.zero_evidence_blocks"
        )
        zero_evidence_block_sites = as_int(
            raw["zero_evidence_block_sites"],
            f"{chrom}.zero_evidence_block_sites",
        )
        primary_active_sites = as_int(
            raw["primary_active_sites"], f"{chrom}.primary_active_sites"
        )
        tree_ready_blocks = as_int(
            raw["tree_ready_blocks"], f"{chrom}.tree_ready_blocks"
        )
        tree_ready_block_sites = as_int(
            raw["tree_ready_block_sites"], f"{chrom}.tree_ready_block_sites"
        )
        abstain_blocks = as_int(
            raw["abstain_blocks"], f"{chrom}.abstain_blocks"
        )
        abstain_block_sites = as_int(
            raw["abstain_block_sites"], f"{chrom}.abstain_block_sites"
        )
        legacy_weight_stable_components = as_int(
            raw["weight_stable_components"],
            f"{chrom}.weight_stable_components",
        )
        corrected_weight_stable_components = (
            as_int(
                exact_by_chrom[chrom]["corrected_weight_stable"],
                f"exact-log {chrom}.corrected_weight_stable",
            )
            if chrom in exact_by_chrom
            else legacy_weight_stable_components
        )
        corrected_weight_sensitive_components = (
            as_int(raw["components"], f"{chrom}.components")
            - corrected_weight_stable_components
        )
        row = {
            "chrom": chrom,
            "chromosome_number": order,
            "partition_stage_status": str(raw["partition_stage_status"]),
            "ssnv_sites": as_int(raw["ssnv_sites"], f"{chrom}.ssnv_sites"),
            "target_components": as_int(raw["components"], f"{chrom}.components"),
            "target_sites": sites,
            "old_selected_sites": as_int(
                raw["old_selected_sites"], f"{chrom}.old_selected_sites"
            ),
            "old_excluded_sites": as_int(
                raw["old_excluded_sites"], f"{chrom}.old_excluded_sites"
            ),
            "new_blocks": new_blocks,
            "new_retained_sites": as_int(
                raw["new_retained_sites"], f"{chrom}.new_retained_sites"
            ),
            "primary_active_sites": primary_active_sites,
            "primary_active_site_ratio": ratio(primary_active_sites, sites),
            "zero_evidence_blocks": zero_evidence_blocks,
            "zero_evidence_block_sites": zero_evidence_block_sites,
            "tree_ready_local_blocks": tree_ready_blocks,
            "tree_ready_local_sites": tree_ready_block_sites,
            "abstain_blocks": abstain_blocks,
            "abstain_block_sites": abstain_block_sites,
            "abstain_components": abstain_components_by_chrom.get(chrom, 0),
            "evidence_positive_components": evidence_components_by_chrom.get(chrom, 0),
            "is_chr6_chr16_extreme": chrom in {"chr6", "chr16"},
            "old_site_retention_ratio": ratio(
                as_int(raw["old_selected_sites"], f"{chrom}.old_selected_sites"),
                sites,
            ),
            "new_site_retention_ratio": ratio(
                as_int(raw["new_retained_sites"], f"{chrom}.new_retained_sites"),
                sites,
            ),
            "raw_total_molecule_weight": raw_weight,
            "old_retained_molecule_weight": old_weight,
            "new_retained_molecule_weight": new_weight,
            "new_lost_molecule_weight": as_int(
                raw["new_lost_molecule_weight"],
                f"{chrom}.new_lost_molecule_weight",
            ),
            "unavoidable_patterns": as_int(
                raw["unavoidable_patterns"], f"{chrom}.unavoidable_patterns"
            ),
            "unavoidable_molecule_weight": as_int(
                raw["unavoidable_molecule_weight"],
                f"{chrom}.unavoidable_molecule_weight",
            ),
            "unavoidable_n_fixed_ra_gt8_patterns": as_int(
                raw["unavoidable_n_fixed_ra_gt8_patterns"],
                f"{chrom}.unavoidable_n_fixed_ra_gt8_patterns",
            ),
            "unavoidable_n_fixed_ra_gt8_molecule_weight": as_int(
                raw["unavoidable_n_fixed_ra_gt8_molecule_weight"],
                f"{chrom}.unavoidable_n_fixed_ra_gt8_molecule_weight",
            ),
            "unavoidable_n_fixed_ra_lte8_span_gt8_patterns": as_int(
                raw["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"],
                f"{chrom}.unavoidable_n_fixed_ra_lte8_span_gt8_patterns",
            ),
            "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight": as_int(
                raw["unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight"],
                f"{chrom}.unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight",
            ),
            "old_read_retention_ratio": old_read_ratio,
            "new_read_retention_ratio": new_read_ratio,
            "read_retention_gain": (
                None
                if new_read_ratio is None or old_read_ratio is None
                else new_read_ratio - old_read_ratio
            ),
            "legacy_weight_stable_components": legacy_weight_stable_components,
            "corrected_weight_stable_components": (
                corrected_weight_stable_components
                if exact_log_audit is not None
                else None
            ),
            "weight_stable_components": corrected_weight_stable_components,
            "legacy_weight_sensitive_components": as_int(
                raw["weight_sensitive_components"],
                f"{chrom}.weight_sensitive_components",
            ),
            "weight_sensitive_components": corrected_weight_sensitive_components,
            "legacy_log_differs_exact": (
                as_int(
                    exact_by_chrom[chrom]["legacy_log_differs_exact"],
                    f"exact-log {chrom}.legacy_log_differs_exact",
                )
                if chrom in exact_by_chrom
                else None
            ),
            "correction_changed_stability": (
                as_int(
                    exact_by_chrom[chrom]["correction_changed_stability"],
                    f"exact-log {chrom}.correction_changed_stability",
                )
                if chrom in exact_by_chrom
                else None
            ),
            "stability_basis": (
                "exact_product corrected"
                if exact_log_audit is not None
                else "legacy log1p (exact audit not provided)"
            ),
            "extraction_minutes": extraction_seconds / 60.0,
            "partition_minutes": (
                None if partition_seconds is None else partition_seconds / 60.0
            ),
            "pattern_load_aggregate_minutes": pattern_load_seconds / 60.0,
            "three_weight_partition_component_loop_minutes": (
                component_loop_seconds / 60.0
            ),
            "stage_wall_sum_minutes": (
                extraction_seconds + (partition_seconds or 0.0)
            )
            / 60.0,
        }
        per_chrom_dataset.append(row)

        if raw_weight > 0:
            for method_order, method, retained, retained_ratio in (
                (1, "densest-8 counterfactual", old_weight, old_read_ratio),
                (2, "新 read-supported", new_weight, new_read_ratio),
            ):
                retention_dataset.append(
                    {
                        "chrom": chrom,
                        "chromosome_number": order,
                        "method": method,
                        "method_order": method_order,
                        "retention_ratio": retained_ratio,
                        "retained_molecule_weight": retained,
                        "raw_total_molecule_weight": raw_weight,
                        "target_components": row["target_components"],
                        "target_sites": sites,
                        "old_selected_sites": row["old_selected_sites"],
                        "new_retained_sites": row["new_retained_sites"],
                    }
                )

        runtime_dataset.append(
            {
                "chrom": chrom,
                "chromosome_number": order,
                "stage": "Extract",
                "stage_detail": "read-linkage extraction",
                "stage_order": 1,
                "wall_minutes": extraction_seconds / 60.0,
                "timing_scope": (
                    "fixture synthetic stage wall"
                    if fixture
                    else "runner stage wall"
                ),
                "is_partition_child_internal": False,
                "fixture_value": fixture,
                "target_components": row["target_components"],
                "target_sites": sites,
                "max_rss_kb": as_int(
                    raw.get("extraction_max_rss_kb", 0),
                    f"{chrom}.extraction_max_rss_kb",
                ),
            }
        )
        runtime_dataset.extend(
            [
                {
                    "chrom": chrom,
                    "chromosome_number": order,
                    "stage": "Partition",
                    "stage_detail": "partition stage total",
                    "stage_order": 2,
                    "wall_minutes": (partition_seconds or 0.0) / 60.0,
                    "timing_scope": (
                        "fixture synthetic partition stage wall"
                        if fixture
                        else "runner stage wall; zero for skipped chr21"
                    ),
                    "is_partition_child_internal": False,
                    "fixture_value": fixture,
                    "target_components": row["target_components"],
                    "target_sites": sites,
                    "max_rss_kb": as_int(
                        raw.get("partition_max_rss_kb", 0) or 0,
                        f"{chrom}.partition_max_rss_kb",
                    ),
                },
                {
                    "chrom": chrom,
                    "chromosome_number": order,
                    "stage": "Load",
                    "stage_detail": "pattern load + aggregate",
                    "stage_order": 3,
                    "wall_minutes": pattern_load_seconds / 60.0,
                    "timing_scope": (
                        "fixture synthetic partition child internal timer"
                        if fixture
                        else "authenticated partition child internal timer"
                    ),
                    "is_partition_child_internal": True,
                    "fixture_value": fixture,
                    "target_components": row["target_components"],
                    "target_sites": sites,
                    "max_rss_kb": None,
                },
                {
                    "chrom": chrom,
                    "chromosome_number": order,
                    "stage": "Loop",
                    "stage_detail": "three-weight partition component loop",
                    "stage_order": 4,
                    "wall_minutes": component_loop_seconds / 60.0,
                    "timing_scope": (
                        "fixture synthetic component loop (DP + diagnostics)"
                        if fixture
                        else (
                            "authenticated partition child component loop "
                            "(DP + diagnostics)"
                        )
                    ),
                    "is_partition_child_internal": True,
                    "fixture_value": fixture,
                    "target_components": row["target_components"],
                    "target_sites": sites,
                    "max_rss_kb": None,
                },
            ]
        )

    component_diagnostic_rows: list[dict[str, Any]] = []
    for raw in component_rows:
        component_label = f"{raw['chrom']}.{raw['legacy_component_id']}"
        new_ratio, old_ratio, retention_gain = validate_component_retention(
            raw,
            component_label,
        )
        start1 = as_int(raw["start1"], "component.start1")
        end1 = as_int(raw["end1"], "component.end1")
        span_bp = as_int(raw["span_bp"], "component.span_bp")
        legacy_weight_stable = as_bool(
            raw["weight_stable"], "component.weight_stable"
        )
        exact_component = exact_by_component.get(
            (raw["chrom"], raw["legacy_component_id"])
        )
        corrected_weight_stable = (
            bool(exact_component["corrected_weight_stable"])
            if exact_component is not None
            else legacy_weight_stable
        )
        primary_active_sites = as_int(
            raw["primary_active_site_count"], "component.primary_active_site_count"
        )
        component_sites = as_int(
            raw["new_site_retained"], "component.new_site_retained"
        )
        retained_weight = as_int(
            raw["raw_retained_molecule_weight"], "component.new_retained"
        )
        is_abstain = raw["status"].startswith("ABSTAIN")
        if retained_weight == 0 or is_abstain:
            evidence_status = "ABSTAIN / zero-retained-weight"
            evidence_status_label = "ABSTAIN"
        elif corrected_weight_stable:
            evidence_status = "Non-ABSTAIN / weight-stable"
            evidence_status_label = "Stable"
        else:
            evidence_status = "Non-ABSTAIN / weight-sensitive"
            evidence_status_label = "Sensitive"
        component_diagnostic_rows.append(
            {
                "chrom": raw["chrom"],
                "chromosome_number": chrom_number(raw["chrom"]),
                "component_id": raw["legacy_component_id"],
                "start1": start1,
                "end1": end1,
                "span_bp": span_bp,
                "start_mb": start1 / 1_000_000.0,
                "end_mb": end1 / 1_000_000.0,
                "span_mb": span_bp / 1_000_000.0,
                "pre_cap_k": as_int(raw["pre_cap_k"], "component.pre_cap_k"),
                "old_selected_sites": as_int(
                    raw["old_densest8_selected"], "component.old_selected"
                ),
                "old_excluded_sites": as_int(
                    raw["old_cap_excluded"], "component.old_excluded"
                ),
                "new_block_count": as_int(
                    raw["new_block_count"], "component.new_block_count"
                ),
                "new_retained_sites": as_int(
                    raw["new_site_retained"], "component.new_site_retained"
                ),
                "primary_active_sites": primary_active_sites,
                "primary_active_site_ratio": ratio(
                    primary_active_sites, component_sites
                ),
                "raw_total_molecule_weight": as_int(
                    raw["raw_total_molecule_weight"], "component.raw_total"
                ),
                "new_retained_molecule_weight": retained_weight,
                "new_lost_molecule_weight": as_int(
                    raw["raw_lost_molecule_weight"], "component.new_lost"
                ),
                "old_retained_molecule_weight": as_int(
                    raw["old_densest8_retained_molecule_weight"],
                    "component.old_retained",
                ),
                "old_read_retention_ratio": old_ratio,
                "new_read_retention_ratio": new_ratio,
                "read_retention_gain": retention_gain,
                "legacy_weight_stable": legacy_weight_stable,
                "corrected_weight_stable": (
                    corrected_weight_stable
                    if exact_component is not None
                    else None
                ),
                "weight_stable": corrected_weight_stable,
                "legacy_log_differs_exact": (
                    exact_component["legacy_log_differs_exact"]
                    if exact_component is not None
                    else None
                ),
                "correction_changed_stability": (
                    exact_component["correction_changed_stability"]
                    if exact_component is not None
                    else None
                ),
                "stability_basis": (
                    "exact_product corrected"
                    if exact_component is not None
                    else "legacy log1p (exact audit not provided)"
                ),
                "weight_stability_status": evidence_status_label,
                "evidence_status_detail": evidence_status,
                "is_abstain": is_abstain,
                "has_retained_component_weight": retained_weight > 0,
                "is_chr6_chr16_extreme": raw["chrom"] in {"chr6", "chr16"},
                "status": raw["status"],
                "positions_sha256_prefix": raw["positions_sha256"][:12],
            }
        )
    genomic_component_rows = sorted(
        (dict(row) for row in component_diagnostic_rows),
        key=lambda row: (
            row["chromosome_number"],
            row["start1"],
            row["end1"],
            row["component_id"],
        ),
    )
    top_component_rows = sorted(
        (
            dict(row)
            for row in component_diagnostic_rows
            if row["read_retention_gain"] is not None
        ),
        key=lambda row: (
            -row["read_retention_gain"],
            -row["raw_total_molecule_weight"],
            row["chromosome_number"],
            row["start1"],
        ),
    )
    top_component_rows = top_component_rows[:top_components]
    for rank, row in enumerate(top_component_rows, start=1):
        row["rank"] = rank

    total_sites = as_int(totals["sites"], "ALL.sites")
    total_raw_weight = as_int(
        totals["raw_total_molecule_weight"], "ALL.raw_total_molecule_weight"
    )
    total_new_weight = as_int(
        totals["new_retained_molecule_weight"],
        "ALL.new_retained_molecule_weight",
    )
    total_old_weight = as_int(
        totals["old_densest8_retained_molecule_weight"],
        "ALL.old_densest8_retained_molecule_weight",
    )
    new_site_ratio = ratio(
        as_int(totals["new_retained_sites"], "ALL.new_retained_sites"), total_sites
    )
    old_site_ratio = ratio(
        as_int(totals["old_selected_sites"], "ALL.old_selected_sites"), total_sites
    )
    new_read_ratio = ratio(total_new_weight, total_raw_weight)
    old_read_ratio = ratio(total_old_weight, total_raw_weight)
    if any(
        value is None
        for value in (new_site_ratio, old_site_ratio, new_read_ratio, old_read_ratio)
    ):
        raise ReportInputError("aggregate report ratios require nonzero denominators")
    assert (
        new_site_ratio is not None
        and old_site_ratio is not None
        and new_read_ratio is not None
        and old_read_ratio is not None
    )

    stage_wall_seconds = (
        as_float(totals["extraction_wall_seconds"], "ALL.extraction_wall_seconds")
        or 0.0
    ) + (
        as_float(
            totals["partition_wall_seconds"],
            "ALL.partition_wall_seconds",
            blank_ok=True,
        )
        or 0.0
    )
    resources = summary.get("resources")
    outer_seconds: float | None = None
    outer_minus_stage_seconds: float | None = None
    runner_overhead_seconds: float | None = None
    outer_resource: Mapping[str, Any] | None = None
    if isinstance(resources, dict) and isinstance(resources.get("outer"), dict):
        outer_resource = resources["outer"]
        outer_seconds = as_float(
            outer_resource.get("elapsed_seconds"),
            "resources.outer.elapsed_seconds",
            blank_ok=True,
        )
    if not fixture and outer_seconds is None:
        raise ReportInputError(
            "formal report requires the validated runner outer GNU time"
        )
    if outer_seconds is not None:
        derived_outer_minus_stage = outer_seconds - stage_wall_seconds
        if fixture:
            outer_minus_stage_seconds = derived_outer_minus_stage
            runner_overhead_seconds = max(0.0, derived_outer_minus_stage)
        else:
            assert outer_resource is not None
            for field in (
                "command_binding_verified",
                "fresh_non_resume_command_verified",
            ):
                if outer_resource.get(field) is not True:
                    raise ReportInputError(
                        f"formal report requires resources.outer.{field}=true"
                    )
            recorded_stage_sum = as_float(
                outer_resource.get("sequential_stage_wall_seconds"),
                "resources.outer.sequential_stage_wall_seconds",
            )
            outer_minus_stage_seconds = as_float(
                outer_resource.get(
                    "outer_minus_sequential_stage_wall_seconds"
                ),
                "resources.outer.outer_minus_sequential_stage_wall_seconds",
            )
            runner_overhead_seconds = as_float(
                outer_resource.get("runner_overhead_seconds"),
                "resources.outer.runner_overhead_seconds",
            )
            if (
                recorded_stage_sum is None
                or outer_minus_stage_seconds is None
                or runner_overhead_seconds is None
                or not math.isclose(
                    recorded_stage_sum,
                    stage_wall_seconds,
                    rel_tol=0,
                    abs_tol=1e-6,
                )
                or not math.isclose(
                    outer_minus_stage_seconds,
                    derived_outer_minus_stage,
                    rel_tol=0,
                    abs_tol=1e-6,
                )
                or not math.isclose(
                    runner_overhead_seconds,
                    max(0.0, derived_outer_minus_stage),
                    rel_tol=0,
                    abs_tol=1e-6,
                )
            ):
                raise ReportInputError(
                    "formal outer GNU time/stage-wall/runner-overhead cross-check failed"
                )
    report_wall_seconds = outer_seconds if outer_seconds is not None else stage_wall_seconds
    runtime_basis = (
        "runner outer /usr/bin/time wall"
        if outer_seconds is not None
        else "per-chromosome stage wall sum"
    )
    total_primary_active_sites = as_int(
        totals["primary_active_sites"], "ALL.primary_active_sites"
    )
    total_new_blocks = as_int(totals["new_blocks"], "ALL.new_blocks")
    total_zero_evidence_blocks = as_int(
        totals["zero_evidence_blocks"], "ALL.zero_evidence_blocks"
    )
    total_zero_evidence_block_sites = as_int(
        totals["zero_evidence_block_sites"], "ALL.zero_evidence_block_sites"
    )
    total_tree_ready_blocks = as_int(
        totals["tree_ready_blocks"], "ALL.tree_ready_blocks"
    )
    total_tree_ready_block_sites = as_int(
        totals["tree_ready_block_sites"], "ALL.tree_ready_block_sites"
    )
    total_abstain_blocks = as_int(
        totals["abstain_blocks"], "ALL.abstain_blocks"
    )
    total_abstain_block_sites = as_int(
        totals["abstain_block_sites"], "ALL.abstain_block_sites"
    )
    total_pattern_load_seconds = as_float(
        totals["partition_pattern_load_aggregate_seconds"],
        "ALL.partition_pattern_load_aggregate_seconds",
    )
    total_component_loop_seconds = as_float(
        totals["partition_ordered_hypergraph_dp_seconds"],
        "ALL.partition_ordered_hypergraph_dp_seconds",
    )
    assert (
        total_pattern_load_seconds is not None
        and total_component_loop_seconds is not None
    )
    total_unavoidable_patterns = as_int(
        totals["unavoidable_patterns"], "ALL.unavoidable_patterns"
    )
    total_unavoidable_weight = as_int(
        totals["unavoidable_molecule_weight"], "ALL.unavoidable_molecule_weight"
    )
    total_unavoidable_dense_patterns = as_int(
        totals["unavoidable_n_fixed_ra_gt8_patterns"],
        "ALL.unavoidable_n_fixed_ra_gt8_patterns",
    )
    total_unavoidable_dense_weight = as_int(
        totals["unavoidable_n_fixed_ra_gt8_molecule_weight"],
        "ALL.unavoidable_n_fixed_ra_gt8_molecule_weight",
    )
    total_unavoidable_sparse_patterns = as_int(
        totals["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"],
        "ALL.unavoidable_n_fixed_ra_lte8_span_gt8_patterns",
    )
    total_unavoidable_sparse_weight = as_int(
        totals["unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight"],
        "ALL.unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight",
    )
    total_abstain_components = sum(abstain_components_by_chrom.values())
    total_evidence_components = sum(evidence_components_by_chrom.values())
    total_legacy_weight_stable = as_int(
        totals["weight_stable_components"], "ALL.weight_stable_components"
    )
    total_corrected_weight_stable = (
        as_int(
            exact_log_audit["counts"]["corrected_weight_stable"],
            "exact-log corrected_weight_stable",
        )
        if exact_log_audit is not None
        else total_legacy_weight_stable
    )
    headline_metrics = [
        {
            "target_components": as_int(totals["components"], "ALL.components"),
            "target_sites": total_sites,
            "old_selected_sites": as_int(
                totals["old_selected_sites"], "ALL.old_selected_sites"
            ),
            "new_retained_sites": as_int(
                totals["new_retained_sites"], "ALL.new_retained_sites"
            ),
            "primary_active_sites": total_primary_active_sites,
            "primary_active_site_ratio": ratio(
                total_primary_active_sites, total_sites
            ),
            "zero_evidence_blocks": total_zero_evidence_blocks,
            "zero_evidence_block_sites": total_zero_evidence_block_sites,
            "tree_ready_local_blocks": total_tree_ready_blocks,
            "tree_ready_local_sites": total_tree_ready_block_sites,
            "abstain_blocks": total_abstain_blocks,
            "abstain_block_sites": total_abstain_block_sites,
            "tree_ready_local_block_ratio": ratio(
                total_tree_ready_blocks, total_new_blocks
            ),
            "abstain_components": total_abstain_components,
            "abstain_component_ratio": ratio(
                total_abstain_components,
                as_int(totals["components"], "ALL.components"),
            ),
            "evidence_positive_components": total_evidence_components,
            "legacy_weight_stable_components": total_legacy_weight_stable,
            "corrected_weight_stable_components": (
                total_corrected_weight_stable
                if exact_log_audit is not None
                else None
            ),
            "weight_stable_components": total_corrected_weight_stable,
            "weight_sensitive_components": (
                as_int(totals["components"], "ALL.components")
                - total_corrected_weight_stable
            ),
            "legacy_log_differs_exact": (
                as_int(
                    exact_log_audit["counts"]["legacy_log_differs_exact"],
                    "exact-log legacy_log_differs_exact",
                )
                if exact_log_audit is not None
                else None
            ),
            "correction_changed_stability": (
                as_int(
                    exact_log_audit["counts"]["correction_changed_stability"],
                    "exact-log correction_changed_stability",
                )
                if exact_log_audit is not None
                else None
            ),
            "primary_raw_partition_changed": (
                False if exact_log_audit is not None else None
            ),
            "stability_basis": (
                "exact_product corrected"
                if exact_log_audit is not None
                else "legacy log1p (exact audit not provided)"
            ),
            "old_site_retention_ratio": old_site_ratio,
            "new_site_retention_ratio": new_site_ratio,
            "site_retention_gain": new_site_ratio - old_site_ratio,
            "old_read_retention_ratio": old_read_ratio,
            "new_read_retention_ratio": new_read_ratio,
            "read_retention_gain": new_read_ratio - old_read_ratio,
            "new_blocks": total_new_blocks,
            "stage_wall_sum_minutes": stage_wall_seconds / 60.0,
            "report_wall_minutes": report_wall_seconds / 60.0,
            "outer_minus_sequential_stage_wall_minutes": (
                None
                if outer_minus_stage_seconds is None
                else outer_minus_stage_seconds / 60.0
            ),
            "runner_overhead_minutes": (
                None
                if runner_overhead_seconds is None
                else runner_overhead_seconds / 60.0
            ),
            "baseline_wall_proxy_minutes": float(baseline["minutes"]),
            "runtime_ratio_vs_baseline_proxy": (
                report_wall_seconds / float(baseline["seconds"])
            ),
            "runtime_basis": runtime_basis,
            "partition_pattern_load_aggregate_minutes": (
                total_pattern_load_seconds / 60.0
            ),
            "three_weight_partition_component_loop_minutes": (
                total_component_loop_seconds / 60.0
            ),
            "unavoidable_patterns": total_unavoidable_patterns,
            "unavoidable_molecule_weight": total_unavoidable_weight,
            "unavoidable_n_fixed_ra_gt8_patterns": (
                total_unavoidable_dense_patterns
            ),
            "unavoidable_n_fixed_ra_gt8_molecule_weight": (
                total_unavoidable_dense_weight
            ),
            "unavoidable_n_fixed_ra_lte8_span_gt8_patterns": (
                total_unavoidable_sparse_patterns
            ),
            "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight": (
                total_unavoidable_sparse_weight
            ),
        }
    ]

    hp_ps_summary_dataset: list[dict[str, Any]] = []
    hp_ps_unit_distribution: list[dict[str, Any]] = []
    hp_ps_worst_units: list[dict[str, Any]] = []
    if hp_ps_unit_audit is not None:
        hp_metrics = dict(hp_ps_unit_audit["metrics"])
        hp_ps_summary_dataset.append(hp_metrics)
        bucket_specs = (
            ("0–<10%", 0.0, 0.1),
            ("10–<20%", 0.1, 0.2),
            ("20–<30%", 0.2, 0.3),
            ("30–<40%", 0.3, 0.4),
            ("40–<50%", 0.4, 0.5),
            ("50–<60%", 0.5, 0.6),
            ("60–<70%", 0.6, 0.7),
            ("70–<80%", 0.7, 0.8),
            ("80–<90%", 0.8, 0.9),
            ("90–<100%", 0.9, 1.0),
            ("100%", 1.0, 1.0),
        )
        for bucket_order, (label, lower, upper) in enumerate(bucket_specs, start=1):
            for eligibility_order, (eligibility, eligible_value) in enumerate(
                (
                    ("Eligible（weight≥20 且 patterns≥5）", True),
                    ("未達 headline eligibility", False),
                ),
                start=1,
            ):
                count = 0
                for unit in hp_ps_unit_audit["units"]:
                    value = float(unit["retention_ratio"])
                    in_bucket = (
                        value == 1.0
                        if lower == upper == 1.0
                        else lower <= value < upper
                    )
                    if (
                        in_bucket
                        and bool(unit["eligible_headline"]) is eligible_value
                    ):
                        count += 1
                hp_ps_unit_distribution.append(
                    {
                        "bucket_order": bucket_order,
                        "retention_bucket": label,
                        "eligibility_order": eligibility_order,
                        "eligibility_status": eligibility,
                        "unit_count": count,
                    }
                )
        hp_ps_worst_units = sorted(
            (
                dict(unit)
                for unit in hp_ps_unit_audit["units"]
                if unit["eligible_headline"]
            ),
            key=lambda unit: (
                unit["retention_ratio"],
                -unit["total_weight"],
                -unit["total_pattern_rows"],
                chrom_number(unit["chrom"]),
                unit["component_id"],
                unit["hp_family"],
                unit["phase_set"],
            ),
        )[:MAX_HP_PS_WORST_UNITS]
        for rank, unit in enumerate(hp_ps_worst_units, start=1):
            unit["rank"] = rank

    span_sensitivity_dataset: list[dict[str, Any]] = []
    if span_grid is not None:
        no_cap_partition_seconds = (
            as_float(
                totals["partition_wall_seconds"],
                "ALL.partition_wall_seconds",
                blank_ok=True,
            )
            or 0.0
        )
        no_cap_blocks = as_int(totals["new_blocks"], "ALL.new_blocks")
        no_cap_retained = as_int(
            totals["new_retained_molecule_weight"],
            "ALL.new_retained_molecule_weight",
        )
        no_cap_lost = as_int(
            totals["new_lost_molecule_weight"], "ALL.new_lost_molecule_weight"
        )
        no_cap_ratio = ratio(no_cap_retained, total_raw_weight)
        assert no_cap_ratio is not None
        span_sensitivity_dataset.append(
            {
                "cap_order": 0,
                "cap_label": "No cap（k≤8 only）",
                "span_cap_bp": None,
                "hard_cap_enabled": False,
                "target_components": as_int(totals["components"], "ALL.components"),
                "target_sites": total_sites,
                "new_blocks": no_cap_blocks,
                "raw_total_molecule_weight": total_raw_weight,
                "raw_retained_molecule_weight": no_cap_retained,
                "raw_lost_molecule_weight": no_cap_lost,
                "read_retention_ratio": no_cap_ratio,
                "retention_delta_vs_no_cap": 0.0,
                "blocks_delta_vs_no_cap": 0,
                "cached_partition_wall_seconds": no_cap_partition_seconds,
                "cached_partition_wall_minutes": no_cap_partition_seconds / 60.0,
                "runtime_ratio_vs_no_cap": 1.0,
                "unavoidable_patterns": as_int(
                    totals.get("unavoidable_patterns", 0),
                    "ALL.unavoidable_patterns",
                ),
                "unavoidable_size_patterns": as_int(
                    totals.get("unavoidable_patterns", 0),
                    "ALL.unavoidable_size_patterns",
                ),
                "unavoidable_span_cap_patterns": 0,
                "unavoidable_both_limits_patterns": 0,
                "timing_scope": "k-only full-run partition stage wall sum",
            }
        )
        for cap_order, aggregate in enumerate(
            span_grid["aggregates"], start=1
        ):
            cap = as_int(aggregate["span_cap_bp"], "span-grid span_cap_bp")
            retained = as_int(
                aggregate["raw_retained_molecule_weight"],
                f"span-grid {cap}.raw_retained_molecule_weight",
            )
            lost = as_int(
                aggregate["raw_lost_molecule_weight"],
                f"span-grid {cap}.raw_lost_molecule_weight",
            )
            retained_ratio = ratio(retained, total_raw_weight)
            assert retained_ratio is not None
            wall_seconds = as_float(
                aggregate["cached_partition_wall_seconds"],
                f"span-grid {cap}.cached_partition_wall_seconds",
            )
            assert wall_seconds is not None
            blocks = as_int(aggregate["new_blocks"], f"span-grid {cap}.new_blocks")
            span_sensitivity_dataset.append(
                {
                    "cap_order": cap_order,
                    "cap_label": f"{cap // 1000} kb",
                    "span_cap_bp": cap,
                    "hard_cap_enabled": True,
                    "target_components": as_int(
                        aggregate["k_gt8_components"],
                        f"span-grid {cap}.k_gt8_components",
                    ),
                    "target_sites": as_int(
                        aggregate["k_gt8_sites"],
                        f"span-grid {cap}.k_gt8_sites",
                    ),
                    "new_blocks": blocks,
                    "raw_total_molecule_weight": total_raw_weight,
                    "raw_retained_molecule_weight": retained,
                    "raw_lost_molecule_weight": lost,
                    "read_retention_ratio": retained_ratio,
                    "retention_delta_vs_no_cap": retained_ratio - no_cap_ratio,
                    "blocks_delta_vs_no_cap": blocks - no_cap_blocks,
                    "cached_partition_wall_seconds": wall_seconds,
                    "cached_partition_wall_minutes": wall_seconds / 60.0,
                    "runtime_ratio_vs_no_cap": (
                        None
                        if no_cap_partition_seconds == 0
                        else wall_seconds / no_cap_partition_seconds
                    ),
                    "unavoidable_patterns": as_int(
                        aggregate["unavoidable_patterns"],
                        f"span-grid {cap}.unavoidable_patterns",
                    ),
                    "unavoidable_size_patterns": as_int(
                        aggregate["unavoidable_size_patterns"],
                        f"span-grid {cap}.unavoidable_size_patterns",
                    ),
                    "unavoidable_span_cap_patterns": as_int(
                        aggregate["unavoidable_span_cap_patterns"],
                        f"span-grid {cap}.unavoidable_span_cap_patterns",
                    ),
                    "unavoidable_both_limits_patterns": as_int(
                        aggregate["unavoidable_both_limits_patterns"],
                        f"span-grid {cap}.unavoidable_both_limits_patterns",
                    ),
                    "timing_scope": (
                        "cached partition only; excludes source validation and BAM extraction"
                    ),
                }
            )

    hashes = {name: sha256_path(path) for name, path in input_paths.items()}
    hash_filters = {
        name: f"{path.name} sha256={hashes[name]}"
        for name, path in input_paths.items()
    }
    sources = [
        source_spec(
            source_id="src_headline_summary",
            label="HCC1395 summarizer headline totals",
            path=safe_source_path(input_paths["summary_json"]),
            sql=(
                'SELECT target_components, target_sites, old_selected_sites, '
                "new_retained_sites, old_site_retention_ratio, "
                "new_site_retention_ratio, old_read_retention_ratio, "
                "new_read_retention_ratio, primary_active_sites, "
                "tree_ready_local_blocks, tree_ready_local_sites, "
                "zero_evidence_blocks, zero_evidence_block_sites, "
                "abstain_blocks, abstain_block_sites, abstain_components, "
                "unavoidable_patterns, unavoidable_molecule_weight, "
                "unavoidable_n_fixed_ra_gt8_patterns, "
                "unavoidable_n_fixed_ra_lte8_span_gt8_patterns "
                'FROM "headline_metrics";'
            ),
            description=(
                "Selects report headline counts and retention rates after the builder "
                "cross-checks summary.json against summary.tsv and component_all.tsv.gz."
            ),
            tables=("artifact_snapshot.headline_metrics",),
            filters=(
                f"sample={DATASET}",
                "pre-cap k>8 components",
                hash_filters["summary_json"],
                hash_filters["summary_tsv"],
                hash_filters["component_all"],
            ),
            metric_definitions=(
                "Target components = count of legacy adjacent-gap components with pre-cap k>8.",
                "Target sites = sum of sites in those pre-cap k>8 components.",
                "Site retention = retained target sites / all target sites.",
                "Read-retention denominator = component-level HP1/HP2 × known-PS exact-pattern molecule-count incidence; one physical molecule may be counted in multiple components, so it is not a genome-wide unique-read count.",
                "TREE_READY_LOCAL requires retained_exact_pattern_count>0, raw_retained_molecule_weight>0, and primary_active_site_count>=2.",
                "ABSTAIN_ZERO_LOCAL_EVIDENCE is every block that fails the TREE_READY_LOCAL gate.",
                "Zero-evidence blocks are the raw-retained-weight=0 subset and do not cover all ABSTAIN causes.",
                "Unavoidable A = n_fixed_ra>8 in the read pattern itself.",
                "Unavoidable B = n_fixed_ra<=8 but span_sites>8 because ordered intervening target sites keep the contiguous block too wide.",
            ),
            generated_at=generated_at,
        ),
        source_spec(
            source_id="src_per_chrom",
            label="HCC1395 per-chromosome reviewed summary",
            path=safe_source_path(input_paths["summary_tsv"]),
            sql=(
                'SELECT * FROM "per_chromosome_metrics" '
                "ORDER BY chromosome_number ASC;"
            ),
            description=(
                "Selects the reviewed chromosome-grain metrics materialized from the "
                "JSON/TSV cross-check."
            ),
            tables=("artifact_snapshot.per_chromosome_metrics",),
            filters=(
                f"sample={DATASET}",
                "chromosome scope supplied by authenticated summarizer",
                hash_filters["summary_json"],
                hash_filters["summary_tsv"],
            ),
            metric_definitions=(
                "Stage wall sum = extraction wall + complete partition-stage wall for one chromosome.",
                (
                    "Weight stable = raw, equal-pattern, and exact-product objectives "
                    "choose identical cuts; corrected exact-product status is used."
                    if exact_log_audit is not None
                    else "Weight stable = legacy raw/equal/log1p status; exact audit not provided."
                ),
            ),
            generated_at=generated_at,
        ),
        source_spec(
            source_id="src_retention",
            label="New method vs same-extraction densest-8 counterfactual",
            path=safe_source_path(input_paths["summary_json"]),
            sql=(
                'SELECT * FROM "retention_by_chrom_method" '
                "ORDER BY chromosome_number ASC, method_order ASC;"
            ),
            description=(
                "Reshapes the same reviewed component-level molecule-count incidence "
                "denominator into a tidy densest-8 counterfactual versus new "
                "read-supported chromosome comparison."
            ),
            tables=("artifact_snapshot.retention_by_chrom_method",),
            filters=(
                "raw_total_molecule_weight > 0",
                f"sample={DATASET}",
                hash_filters["summary_json"],
            ),
            metric_definitions=(
                "Densest-8 retention is a counterfactual recomputed on the new BASEQ20 extraction, not empirical retention from the historical v5 run.",
                "Old counterfactual retention = component-level molecule-count incidence whose complete fixed R/A target-site set for that component constraint falls in the densest contiguous eight sites / the same component-level incidence denominator.",
                "New retention = component-level molecule-count incidence whose complete fixed R/A target-site set for that component constraint is contained in one selected disjoint k<=8 block / the same denominator.",
            ),
            generated_at=generated_at,
        ),
        source_spec(
            source_id="src_runtime",
            label="New run runtime and audited historical proxy",
            path=None,
            sql=(
                'SELECT h.report_wall_minutes, h.baseline_wall_proxy_minutes, '
                'h.runtime_ratio_vs_baseline_proxy, r.* FROM "headline_metrics" h '
                'JOIN "runtime_by_chrom_stage" r ON 1=1 '
                "ORDER BY r.chromosome_number ASC, r.stage_order ASC;"
            ),
            description=(
                "Combines the new summarizer runtime receipts with the separately "
                "audited old HCC1395 historical wall proxy. The old value is not a "
                "monotonic benchmark. Partition-child internal timers are nested "
                "within, not additive to, the partition stage wall."
            ),
            tables=(
                "artifact_snapshot.headline_metrics",
                "artifact_snapshot.runtime_by_chrom_stage",
                "baseline-runtime-audit.md",
            ),
            filters=(
                f"sample={DATASET}",
                hash_filters["summary_json"],
                hash_filters["baseline_audit"],
            ),
            metric_definitions=(
                "New report wall = runner outer GNU time when available; otherwise explicitly labelled sum of per-chromosome stage walls.",
                "Historical baseline = 5,086.484135464-second filesystem birth-timestamp wall proxy audited from the terminal-success v5 run.",
                "Runtime ratio = new report wall seconds / audited historical wall-proxy seconds.",
                "Runner overhead = max(0, outer GNU time elapsed - sequential extraction/partition stage-wall sum); the raw signed difference is retained separately.",
                "Pattern load/aggregate = authenticated child receipt timings_seconds.load_and_aggregate_primary_patterns.",
                "Three-weight partition component loop (DP + diagnostics) = authenticated child receipt timings_seconds.ordered_hypergraph_dp; includes raw/equal/log DP passes, the densest-8 baseline, diagnostics, row materialization, and aggregation; excludes BAM extraction, pattern loading, final writes, and process overhead.",
            ),
            generated_at=generated_at,
        ),
        source_spec(
            source_id="src_genomic_components",
            label="All k>8 component genomic regions",
            path=safe_source_path(input_paths["component_all"]),
            sql=(
                'SELECT * FROM "genomic_components" '
                "ORDER BY chromosome_number ASC, start1 ASC, end1 ASC;"
            ),
            description=(
                "Selects every reviewed component at component grain for the "
                "genome-region scatter; no top-N filtering is applied."
            ),
            tables=("artifact_snapshot.genomic_components",),
            filters=(
                f"sample={DATASET}",
                "pre_cap_k > 8",
                "all authenticated component rows",
                hash_filters["component_all"],
            ),
            metric_definitions=(
                "Start Mb = 1-based component start / 1,000,000.",
                "Bubble size = pre-cap component k.",
                (
                    "Color class combines component ABSTAIN/zero-retained-weight "
                    "status with corrected exact-product weight stability."
                    if exact_log_audit is not None
                    else "Color class uses legacy weight stability because exact audit is absent."
                ),
                "Block-level TREE_READY_LOCAL is reported separately and additionally requires a retained exact pattern and at least two primary-active sites.",
            ),
            generated_at=generated_at,
        ),
        source_spec(
            source_id="src_components",
            label="Ranked k>8 component diagnostics",
            path=safe_source_path(input_paths["component_all"]),
            sql=(
                'SELECT * FROM "top_components" '
                "ORDER BY read_retention_gain DESC, raw_total_molecule_weight DESC "
                f"LIMIT {top_components};"
            ),
            description=(
                "Ranks reviewed component rows by new minus same-extraction densest-8 "
                "counterfactual retention gain; the snapshot retains only the bounded "
                "top diagnostic rows."
            ),
            tables=("artifact_snapshot.top_components",),
            filters=(
                f"sample={DATASET}",
                "pre_cap_k > 8",
                f"top {top_components} by read-retention gain",
                hash_filters["component_all"],
            ),
            metric_definitions=(
                "Component retention gain = new retention ratio - the densest-8 counterfactual recomputed on the same BASEQ20 extraction and component-level incidence denominator.",
                "Span bp = component end1 - start1 as emitted by the summarizer.",
            ),
            generated_at=generated_at,
        ),
        {
            "id": "src_baseline_audit",
            "label": "HCC1395 old densest-8 runtime baseline audit",
            "path": safe_source_path(input_paths["baseline_audit"]),
            "query": {
                "description": (
                    "Read-only methodology audit establishing the 5,086.484135464-second "
                    "historical downstream wall proxy and its comparison limits."
                ),
                "executed_at": generated_at,
                "filters": [hash_filters["baseline_audit"]],
            },
        },
    ]
    if exact_log_audit is not None:
        exact_filter_keys = sorted(
            name for name in hash_filters if name.startswith("exact_")
        )
        sources.append(
            source_spec(
                source_id="src_exact_log_audit",
                label="Exact-product log-sensitivity remediation audit",
                path=safe_source_path(input_paths["exact_receipt_1"]),
                sql=(
                    'SELECT legacy_weight_stable_components, '
                    "corrected_weight_stable_components, "
                    "legacy_log_differs_exact, correction_changed_stability "
                    'FROM "headline_metrics";'
                ),
                description=(
                    "Authenticates the exact integer-product replacement for the "
                    "legacy floating log1p comparison and materializes corrected "
                    "stability without changing primary raw-molecule cuts."
                ),
                tables=(
                    "artifact_snapshot.headline_metrics",
                    "artifact_snapshot.per_chromosome_metrics",
                    "artifact_snapshot.genomic_components",
                ),
                filters=tuple(
                    [
                        f"sample={DATASET}",
                        "every report component has one authenticated exact-log row",
                    ]
                    + [hash_filters[name] for name in exact_filter_keys]
                ),
                metric_definitions=(
                    "Corrected weight stable = raw, equal-pattern, and exact-product objectives choose identical cuts.",
                    "Legacy-log differs exact = the legacy floating log1p cut "
                    "differs from the exact integer-product cut.",
                    "Correction changed stability = corrected exact-product "
                    "stability differs from the legacy reported stability.",
                    "Primary raw partition changed = false; remediation changes only the log-sensitivity comparison.",
                ),
                generated_at=generated_at,
            )
        )
    if hp_ps_unit_audit is not None:
        hp_filter_keys = sorted(
            name for name in hash_filters if name.startswith("hp_ps_")
        )
        sources.append(
            source_spec(
                source_id="src_hp_ps_unit_audit",
                label="HP/PS observed-constraint unit retention audit v5 contract",
                path=safe_source_path(input_paths["hp_ps_receipt"]),
                sql='SELECT * FROM "hp_ps_summary_metrics";',
                description=(
                    "Reports observed HP/PS exact-phase-set unit retention at "
                    "molecule-by-component-incidence grain, including the full "
                    "eligible-unit distribution and its worst supported tail."
                ),
                tables=(
                    "artifact_snapshot.hp_ps_summary_metrics",
                    "artifact_snapshot.hp_ps_unit_distribution",
                    "artifact_snapshot.hp_ps_worst_units",
                ),
                filters=tuple(
                    [
                        f"sample={DATASET}",
                        "scope ceiling = observed_constraint_units_only",
                        "headline eligible = total weight >=20 and pattern rows >=5",
                        "no synthetic rows for components without observed constraints",
                    ]
                    + [hash_filters[name] for name in hp_filter_keys]
                ),
                metric_definitions=(
                    "Unit = dataset × chromosome × legacy component × HP1/HP2 × exact known phase set.",
                    "Retention = retained molecule×component incidence / total observed incidence for that unit.",
                    "Weighted retention = sum retained incidence / sum total "
                    "incidence; it is not unique-read coverage.",
                    "Eligible quantiles are unweighted type-7 quantiles across "
                    "units; p10 is recomputed from authenticated unit rows.",
                    "Components without observed constraint units have no "
                    "reconstructed denominator and are never assigned synthetic "
                    "0% or 100% retention.",
                    "HP1/HP2 component delta compares component-aggregated "
                    "retention only where both HP families are observed.",
                ),
                generated_at=generated_at,
            )
        )
    if span_grid is not None:
        sources.append(
            source_spec(
                source_id="src_span_sensitivity",
                label="HCC1395 cached physical-span sensitivity grid",
                path=safe_source_path(input_paths["span_summary"]),
                sql=(
                    'SELECT * FROM "span_sensitivity" '
                    "ORDER BY cap_order ASC;"
                ),
                description=(
                    "Compares the authenticated k-only no-cap baseline with cached "
                    "50/100/200 kb hard-span partitions using the same target and "
                    "component-level molecule-count-incidence denominator."
                ),
                tables=(
                    "artifact_snapshot.span_sensitivity",
                    "span-grid summary.tsv",
                    "span-grid receipt.json",
                ),
                filters=(
                    f"sample={DATASET}",
                    "span caps = 50000,100000,200000 bp",
                    hash_filters["summary_json"],
                    hash_filters["span_summary"],
                    hash_filters["span_receipt"],
                ),
                metric_definitions=(
                    "Hard span cap is inclusive: block end1-start1 must be <= cap while k remains <=8.",
                    "Read retention uses the common component-level HP1/HP2 × known-PS exact-pattern molecule-count-incidence denominator; it is not a genome-wide unique-read count.",
                    "New blocks = sum of selected local blocks across the supplied chromosome scope.",
                    "Cached partition wall excludes source validation and upstream BAM extraction.",
                    "No-cap runtime is the k-only full-run partition-stage wall sum; capped runtimes are cached re-partition wall sums.",
                ),
                generated_at=generated_at,
            )
        )

    fixture_note = (
        "**FIXTURE 狀態：本頁只用 chr1／chr22 架構資料驗證欄位、圖表與來源連結；"
        "所有 fixture 數字都不可當作 HCC1395 chr1–chr22 最終結果。**"
    )
    ready_note = (
        "**正式狀態：summarizer 的 chr1–chr22、JSON/TSV/component 守恆與所有 "
        "checks 均通過，且 exact-log 與 full HP/PS unit audits 已綁定；"
        "以下數字可作本輪 frozen run 的完整工程證據。**"
    )
    title = (
        "[FIXTURE] HCC1395 k>8 Read-supported 分割報告架構預覽"
        if fixture
        else "HCC1395 k>8 Read-supported 分割完整技術報告"
    )
    status_note = fixture_note if fixture else ready_note
    runtime_comparison = (
        f"新流程採 `{runtime_basis}`，為 {fmt_minutes(report_wall_seconds / 60.0)}；"
        f"相對舊流程 {float(baseline['minutes']):.2f} 分鐘歷史 wall proxy，"
        f"比值為 {report_wall_seconds / float(baseline['seconds']):.2f}×。"
        + (
            "外層 wall 減 sequential extraction＋partition stage wall 的 runner "
            f"overhead 為 {fmt_minutes(runner_overhead_seconds / 60.0)}。"
            if runner_overhead_seconds is not None
            else ""
        )
    )
    runtime_internal_summary = (
        "其中 authenticated partition child internal timers 顯示：pattern "
        f"load／aggregate 合計 {fmt_minutes(total_pattern_load_seconds / 60.0)}，"
        "three-weight partition component loop（含 DP + diagnostics）合計 "
        f"{fmt_minutes(total_component_loop_seconds / 60.0)}。兩者是 partition "
        "stage 內部子區段，不能與 stage wall 相加，也不包含 BAM extraction。"
    )
    unavoidable_summary = (
        f"Unavoidable 共 {fmt_integer(total_unavoidable_patterns)} patterns／"
        f"{fmt_integer(total_unavoidable_weight)} molecule weight；其中 read "
        f"本身 n_fixed_ra>8 為 {fmt_integer(total_unavoidable_dense_patterns)}／"
        f"{fmt_integer(total_unavoidable_dense_weight)}，n_fixed_ra≤8 但跨越 >8 "
        "ordered target sites 為 "
        f"{fmt_integer(total_unavoidable_sparse_patterns)}／"
        f"{fmt_integer(total_unavoidable_sparse_weight)}。"
    )
    evidence_summary = (
        f"工程上共指派 {fmt_integer(total_sites)} 個位點，但 primary-active 僅 "
        f"{fmt_integer(total_primary_active_sites)} 個；"
        f"{fmt_integer(total_zero_evidence_blocks)} 個區塊（"
        f"{fmt_integer(total_zero_evidence_block_sites)} 個位點）沒有 retained raw "
        f"weight，因此 tree-ready local blocks 為 "
        f"{fmt_integer(total_tree_ready_blocks)}／{fmt_integer(total_new_blocks)}；"
        f"其餘 {fmt_integer(total_abstain_blocks)} 個 blocks（"
        f"{fmt_integer(total_abstain_block_sites)} sites）為 ABSTAIN。另有 "
        f"{fmt_integer(total_abstain_components)} 個 component 標為 ABSTAIN。"
    )
    per_chrom_lookup = {row["chrom"]: row for row in per_chrom_dataset}
    if {"chr6", "chr16"}.issubset(per_chrom_lookup):
        chr6 = per_chrom_lookup["chr6"]
        chr16 = per_chrom_lookup["chr16"]
        extreme_components = chr6["target_components"] + chr16["target_components"]
        extreme_sites = chr6["target_sites"] + chr16["target_sites"]
        extreme_summary = (
            f"chr6 與 chr16 合計 {fmt_integer(extreme_components)} 個 component"
            f"（{fmt_percent(extreme_components / headline_metrics[0]['target_components'])}）、"
            f"{fmt_integer(extreme_sites)} 個目標位點"
            f"（{fmt_percent(extreme_sites / total_sites)}），是區域圖必須保留且不能用"
            "全基因組平均掩蓋的極端集中區。"
        )
    else:
        extreme_summary = (
            "目前 fixture 不含 chr6／chr16；正式 22-chrom 報告會在區域圖完整顯示兩者。"
        )
    if span_sensitivity_dataset:
        span_status_sentence = (
            "50／100／200 kb hard-cap grid 已通過 receipt、summary 與共同分母核對。"
        )
    else:
        span_status_sentence = (
            "未提供 `--span-grid-summary`；hard-cap sensitivity 尚未跑入本報告，"
            "但不阻擋 k-only 結果。"
        )
    if exact_log_audit is not None:
        exact_counts = exact_log_audit["counts"]
        exact_status_sentence = (
            "Exact-product remediation 已驗證全部 component：corrected "
            f"weight-stable {fmt_integer(exact_counts['corrected_weight_stable'])}／"
            f"{fmt_integer(exact_counts['components'])}；legacy log 與 exact cut "
            f"不同 {fmt_integer(exact_counts['legacy_log_differs_exact'])} 個，"
            "修正後改變 stability 分類 "
            f"{fmt_integer(exact_counts['correction_changed_stability'])} 個。"
            "**Primary raw-molecule cuts 完全未改動。**"
        )
        exact_audit_body = (
            "## Exact-product 修正後的 weight stability\n\n"
            f"{exact_status_sentence}\n\n"
            "Headline、染色體表與基因組 component 圖的 `weight_stable` 已改用 "
            "`corrected_weight_stable`；legacy 值仍保留供稽核。這項修正只把 "
            "`sum(log(m+1))` 的浮點比較換成數學等價的整數 `product(m+1)`，"
            "不重切 primary raw partition，也不驗證 clone、topology 或真實樹。"
        )
    else:
        exact_status_sentence = (
            "未提供 `--exact-log-audit`；fixture 僅顯示 legacy log1p stability，"
            "不可當 corrected stability 證據。"
        )
        exact_audit_body = (
            "## Exact-product stability audit 尚未提供\n\n"
            f"{exact_status_sentence}"
        )
    if hp_ps_unit_audit is not None:
        hp_metrics = hp_ps_summary_dataset[0]
        hp_status_sentence = (
            f"HP/PS observed-constraint audit 覆蓋 "
            f"{fmt_integer(hp_metrics['components_in_partition_scope'])} 個 "
            f"partition components；其中 "
            f"{fmt_integer(hp_metrics['components_with_observed_constraint_units'])} "
            "個有 observed constraint units、"
            f"{fmt_integer(hp_metrics['components_without_observed_constraint_units'])} "
            "個沒有可識別分母。Observed units 共 "
            f"{fmt_integer(hp_metrics['observed_constraint_units'])}，headline eligible "
            f"{fmt_integer(hp_metrics['eligible_headline_units'])}。"
        )
        hp_ps_body = (
            "## HP/PS unit retention：分布與高支持尾端\n\n"
            f"{hp_status_sentence} molecule×component-incidence weighted retention "
            f"為 {fmt_percent(hp_metrics['weighted_retention_ratio'])}；eligible unit "
            f"retention 的 median／p10／p25 分別為 "
            f"{fmt_percent(hp_metrics['eligible_retention_median'])}／"
            f"{fmt_percent(hp_metrics['eligible_retention_p10'])}／"
            f"{fmt_percent(hp_metrics['eligible_retention_p25'])}，低於 0.5／0.8 "
            f"分別有 {fmt_integer(hp_metrics['eligible_retention_lt_0_5'])}／"
            f"{fmt_integer(hp_metrics['eligible_retention_lt_0_8'])} units。\n\n"
            f"Component coverage 分為 HP1 only "
            f"{fmt_integer(hp_metrics['components_hp1_only'])}、HP2 only "
            f"{fmt_integer(hp_metrics['components_hp2_only'])}、HP1+HP2 "
            f"{fmt_integer(hp_metrics['components_hp1_and_hp2'])}；兩側都達 eligibility "
            f"的 paired components 為 "
            f"{fmt_integer(hp_metrics['paired_components_both_eligible'])}；所有 "
            f"paired components 的 |HP1−HP2 retention| median 為 "
            f"{fmt_percent(hp_metrics['paired_component_absolute_delta_median'])}，"
            f"≥25 percentage points 有 "
            f"{fmt_integer(hp_metrics['paired_component_absolute_delta_ge_0_25'])} 個。"
            "下圖顯示所有 observed units 的分布，表格則直接列出最低 retention "
            f"的最多 {MAX_HP_PS_WORST_UNITS} 個 eligible units，避免只看 aggregate "
            "weighted retention 掩蓋低尾。\n\n"
            "**沒有 observed constraint 的 component 不建立 unit row，也不填造 "
            "0% 或 100% retention。** 這裡的 weight 是 molecule×component "
            "incidence，不是 unique molecules。此工程 audit 不用 VAF，也不判斷"
            "姐妹／直系／混合拓撲；若下游有多棵候選樹，要選「最可能」樹，必須另以"
            "每個位點 VAF 的 likelihood／一致性做排序。"
        )
        hp_ps_worst_body = (
            "### Lowest-retention eligible units\n\n"
            "此表只納入已達 weight≥20 且 pattern rows≥5 的 units，依 retention "
            "由低到高，再依 total incidence 與 pattern rows 由高到低排序。"
            "因此它直接揭露高支持低尾，但不是生物學重要性、clone 大小或拓撲排名。"
        )
    else:
        hp_status_sentence = (
            "未提供 `--hp-ps-unit-audit`；fixture 不顯示 unit retention，也不對"
            "沒有 observed constraint 的 component 填造 0% 或 100%。"
        )
        hp_ps_body = (
            "## HP/PS unit retention audit 尚未提供\n\n"
            f"{hp_status_sentence}"
        )
        hp_ps_worst_body = ""
    executive_body = (
        "## Executive Summary\n\n"
        f"{status_note}\n\n"
        f"本報告範圍含 {fmt_integer(len(per_chrom_dataset))} 條染色體、"
        f"{fmt_integer(headline_metrics[0]['target_components'])} 個 k>8 區域與 "
        f"{fmt_integer(total_sites)} 個目標位點。新法工程位點指派率為 "
        f"{fmt_percent(new_site_ratio)}，densest-8 counterfactual 為 "
        f"{fmt_percent(old_site_ratio)}；以同一 component-level molecule-count "
        f"incidence 分母計算，新法 constraint read retention 為 "
        f"{fmt_percent(new_read_ratio)}，同一 BASEQ20 extraction 上重算的 "
        f"densest-8 counterfactual 為 {fmt_percent(old_read_ratio)}。"
        f"{evidence_summary} "
        f"{unavoidable_summary} {extreme_summary} {runtime_comparison} "
        f"{runtime_internal_summary} {exact_status_sentence} "
        f"{hp_status_sentence} {span_status_sentence}\n\n"
        "**結論界線：工程 assignment 只證明位點被分配到 k≤8 block；"
        "它不等於該 block 有 read-supported topology evidence，更不等於唯一演化樹。**"
    )
    concept_body = (
        "## 先讀懂 assignment、evidence 與 tree-ready\n\n"
        "### 工程位點指派率\n"
        "問的是「原本 k>8 區域中的 sSNV 位點，有多少被分配進 k≤8 區塊」。"
        "它只衡量工程覆蓋；100% assignment 不代表可分析、救回或已有 topology。\n\n"
        "### Primary-active sites 與 zero-evidence blocks\n"
        "Primary-active sites 至少參與 primary linkage；zero-evidence block 的 "
        "retained raw-molecule weight 為 0。後者即使已接收位點，也不能當作樹證據；"
        "但其中仍可能有曾出現在被切斷或 unavoidable constraint 的 primary-active "
        "sites，因此不等於完全沒有觀測。\n\n"
        "### Tree-ready local blocks\n"
        "嚴格 gate 同時要求 retained exact pattern > 0、retained raw-molecule "
        "weight > 0、primary-active sites ≥2；任一不滿足即 "
        "ABSTAIN_ZERO_LOCAL_EVIDENCE。通過只表示可進入後續局部候選樹推論，"
        "尚未證明 topology 唯一或樹結構唯一。\n\n"
        "### Unavoidable 不只有一種原因\n"
        f"{unavoidable_summary} 第一類是 read pattern 自己已有超過 8 個固定 "
        "R/A calls；第二類 calls 數沒有超過 8，但 endpoints 中間夾著其他 "
        "target sites，使 contiguous ordered block 的 span_sites 超過 8。"
        "兩類不可混稱成「read 都太長」。\n\n"
        "### Read retention\n"
        "分母是 component-level HP1/HP2 × known-PS exact-pattern "
        "molecule-count incidence；同一 physical molecule 橫跨多個 component "
        "時可在各 component 重複計數，因此不是全基因組 unique reads。該 "
        "component constraint 內全部 fixed R/A target sites 必須完整落在同一"
        "區塊才算新法保留。densest-8 是在相同 BASEQ20 extraction 上重算的"
        "counterfactual，不是舊 v5 的實測 retention。\n\n"
        "### Weight stable\n"
        "正式 headline 使用 raw-molecule、equal-pattern 與 exact-product "
        "三種目標選到相同 cut 的 `corrected_weight_stable`；legacy log1p "
        "結果只保留作稽核。sensitive 代表切點會受權重定義影響。ABSTAIN 則需與 "
        "zero-evidence 分開看，不能用 100% assignment 掩蓋。"
    )
    retention_body = (
        "## 新法 vs densest-8 counterfactual：每條染色體的 constraint 證據保留\n\n"
        "圖中每一對柱使用相同 component-level molecule-count-incidence 分母；"
        "差距表示新分割相對同一 BASEQ20 extraction 上重算的 densest-8 "
        "counterfactual 多保留多少 constraint 證據，而不是與舊 v5 實測 reads "
        "相比。零分母染色體不畫入比例圖。"
    )
    target_body = (
        "## k>8 工程目標位點分布\n\n"
        "柱高是各染色體進入本輪分割的目標位點數，用來辨識工作量集中位置；"
        "它不是突變負荷，也不是 evidence-positive 或 tree-ready 位點數。"
    )
    genomic_body = (
        "## 基因組區域圖：位置、k 與 evidence status\n\n"
        f"每一點是一個完整 component；x 為起點 Mb、y 為染色體序號、泡泡大小為 "
        f"pre-cap k，顏色區分 component ABSTAIN／zero-retained-weight 與 "
        + (
            "corrected exact-product weight stability。"
            if exact_log_audit is not None
            else "legacy weight stability（exact audit 未提供）。"
        )
        + " Tooltip 保留 start/end/span、primary-active、read retention 與 status。"
        + f" {extreme_summary} 這張圖顯示工程與證據分布，不是基因註解或演化樹圖。"
    )
    runtime_body = (
        "## 每條染色體的 runtime\n\n"
        "圖分成 BAM/read-linkage extraction、完整 partition stage，以及 "
        "partition child 內的 pattern load／aggregate 與 three-weight partition "
        "component loop（含 DP + diagnostics）。"
        "後兩者是完整 partition stage 的內含子區段，不可再次相加。這可避免"
        "把主要 I/O、read 掃描、Python 啟動或 output write 等等待誤認成 "
        f"component-loop 成本。舊 {float(baseline['minutes']):.2f} "
        "分鐘數字是 filesystem "
        "birth-timestamp 歷史 wall proxy，不是同機、同負載、單一 `/usr/bin/time` "
        "benchmark；舊 run 為 BASEQ_MIN=0 且與其他 samples 並行，新 run 為 "
        "BASEQ_MIN=20、chr1–22 sequential 且繼承 nice/ionice。"
        f"{runtime_comparison} {runtime_internal_summary}"
    )
    per_chrom_body = (
        "## 染色體層級稽核表\n\n"
        "表格保留位點守恆、read retention、"
        + (
            "corrected exact-product weight sensitivity"
            if exact_log_audit is not None
            else "legacy weight sensitivity"
        )
        + " 與 stage wall；"
        + "預設依目標位點數由大到小，先看對總工作量影響最大的染色體。"
    )
    components_body = (
        "## Read-retention 改善最大的元件\n\n"
        f"下表只顯示 gain 最大的 {len(top_component_rows)} 個 component，"
        "用途是診斷與人工 review，不代表演化重要性排名。"
    )
    if span_sensitivity_dataset:
        no_cap_span = span_sensitivity_dataset[0]
        capped_by_label = {
            row["cap_label"]: row for row in span_sensitivity_dataset[1:]
        }
        span_retention_body = (
            "## Hard-cap sensitivity：read retention\n\n"
            "Hard cap 同時要求 `k≤8` 且區塊內 `end1-start1≤cap`。"
            f"No-cap retention 為 {fmt_percent(no_cap_span['read_retention_ratio'])}；"
            + "、".join(
                f"{label} 為 {fmt_percent(capped_by_label[label]['read_retention_ratio'])}"
                for label in ("50 kb", "100 kb", "200 kb")
            )
            + "。更低 retention 表示跨越實體距離上限的 read pattern 無法完整留在同一 block，"
            "不是 variant 被刪除。"
        )
        span_blocks_body = (
            "### Hard cap 對 block 數的影響\n\n"
            "較嚴格的實體距離上限通常會把同一 component 切成更多 local blocks。"
            "Block 數增加是幾何切割結果，不能直接解讀為更多 subclones 或更多樹。"
        )
        span_runtime_body = (
            "### Cached partition runtime 的比較界線\n\n"
            "Capped 數字只含 cached partition wall，排除 source validation 與 BAM "
            "extraction；no-cap 使用 full-run 的 k-only partition-stage wall sum。"
            "兩者可用於 partition operational cost 比較，但不是 raw-to-final runtime。"
        )
    else:
        span_retention_body = (
            "## Hard-cap sensitivity 尚未提供\n\n"
            "本次沒有 `--span-grid-summary`，因此不顯示 no-cap 對 50／100／200 kb "
            "的 retention、block 數與 cached runtime 圖。這是可選證據缺口，不會把 "
            "k-only snapshot 標成 partial；提供通過驗證的 final summary／receipt 後即可補入。"
        )
        span_blocks_body = ""
        span_runtime_body = ""
    span_caveat = (
        "- Hard-cap grid 已驗證，但 capped runtime 只涵蓋 cached partition。\n"
        if span_sensitivity_dataset
        else "- Hard-cap grid 未提供；本報告不對 50／100／200 kb sensitivity 下結論。\n"
    )
    audit_methodology = (
        "6. Exact-log audit 逐 component 驗證 source/output identity、守恆與 "
        "`corrected_weight_stable`，並確認 primary raw cuts 未變。\n"
        if exact_log_audit is not None
        else "6. Exact-log audit 未提供；fixture 只保留 legacy stability。\n"
    ) + (
        "7. HP/PS v5 contract 驗證 receipt-summary binding、每條染色體 11 項 "
        "source checks、observed-unit mass、component coverage 與 HP1/HP2 pair delta。\n"
        if hp_ps_unit_audit is not None
        else "7. HP/PS unit audit 未提供；不產生 retention 分母或估計。\n"
    )
    methodology_body = (
        "## 方法與可重現性\n\n"
        "1. Builder 先逐欄交叉核對 `summary.json` 與 `summary.tsv`。\n"
        "2. 再以 `component_all.tsv.gz` 重算 component 數、位點與 molecule-weight "
        "守恆。\n"
        "3. Tree-ready 採 summarizer block gate：retained pattern >0、raw weight >0、"
        "primary-active sites ≥2；不以 `new_blocks-zero_evidence` 近似。\n"
        "4. 圖表只使用 snapshot 中已審核的 plain row arrays；每個 dataset "
        f"最多 {MAX_ROWS_PER_DATASET:,} 列。\n"
        "5. 所有 source 都保存輸入 SHA-256；報告不含 machine-local 絕對路徑。\n"
        + audit_methodology
        + (
            "8. Optional span grid 另核對 final receipt、summary identity、cap totals "
            "與共同 component-level molecule-count-incidence denominator。\n\n"
            if span_sensitivity_dataset
            else "8. Optional span grid 未提供，故只產生 k-only 報告。\n\n"
        )
        + "Read-supported segmentation 的結論上限是「在 frozen inputs 上，產生可重現、"
        "守恆的局部 k≤8 分割，並標出可進入 candidate-tree inference 的 blocks」。"
        "它不等於辨識唯一真實的腫瘤演化樹。"
    )
    caveat_body = (
        "## 資料狀態與限制\n\n"
        f"{status_note}\n\n"
        "- 舊 runtime 是歷史 wall proxy；只能做 operational cost 比較，不能宣稱精確加速比。\n"
        "- Runtime scope 並非同一實驗：舊流程 BASEQ_MIN=0、與其他 samples 並行且含後續"
        "樹／region／ledger；新流程 BASEQ_MIN=20、chr1–22 sequential，停在 sparse "
        "extraction＋k>8 segmentation。Read-retention 的 densest-8 是在新 BASEQ20 "
        "extraction 與同一 component-level incidence 分母上重算的 counterfactual，"
        "不是舊 v5 實測 retention，因此不受 runtime scope 差異混淆。\n"
        "- 新 full run 啟動時 load1=25.2（48 logical CPUs），但 07:31 spot check "
        "load1=67.42 且 `/big7_disk` 裝置在 1 秒 iostat 視窗為 100% util；本次正式 "
        "wall 是 shared-host、nice/ionice 低優先序的實際 operational time，不是隔離式"
        "微基準。\n"
        "- Read retention 是 component-level HP1/HP2 × known-PS exact-pattern "
        "molecule-count incidence；同一 physical molecule 可跨 component 重複，"
        "不是全基因組 unique-read count，也不是 VAF、coverage 或 variant-caller F1。\n"
        "- 100% 位點 assignment 只代表工程分配；primary-active、TREE_READY_LOCAL "
        "與 ABSTAIN 必須分開報告。\n"
        "- Zero-evidence 是 ABSTAIN 的子集；active sites <2 或無 retained pattern "
        "也會使 block ABSTAIN。Zero retained local weight 仍可能同時存在被切斷或 "
        "unavoidable constraint 的 active calls，不能改寫成完全沒有 read observation。\n"
        "- Component gain 高只表示相對同 extraction densest-8 counterfactual 保存"
        "更多同分母 constraint 證據。\n"
        "- HP/PS unit audit 的 scope ceiling 是 observed constraints；無 observed "
        "constraint 不等於 0% 或 100%，也不等於沒有 biological relation。\n"
        "- 此 contract 僅涵蓋 primary HP1／HP2，不能推導 HP3／not-HP3 分類；"
        "也未計算無 HP 關係、姐妹 only、直系 only、姐妹＋直系或 Topo>1。\n"
        "- VAF 未用於此工程 retention audit；多棵候選樹的最可能排序須在下游"
        "逐位點納入 VAF likelihood，不能由本圖自行宣稱。\n"
        "- 報告不推論唯一 topology、唯一 tree structure 或臨床效益。\n"
        + span_caveat
        + (
            "- 目前為 fixture；必須等正式 22-chrom summarizer 輸出通過後重建，"
            "才能移除 fixture 標籤。"
            if fixture
            else "- 本報告僅對 summary 所指向的 frozen HCC1395 run 有效。"
        )
    )
    next_body = (
        "## 下一步判讀\n\n"
        "先檢查 read-retention gain 最大且 weight-sensitive 的 component，再對照"
        "零證據區塊、span 與 cut constraints。若要談 topology 或樹結構，必須另接"
        "候選樹、VAF likelihood 與跨樣本一致性證據；本報告本身不越過該證據界線。"
    )

    cards = [
        {
            "id": "card_target_components",
            "description": "歷史相鄰距離 component 在 cap 前 k>8 的區域數。",
            "dataset": "headline_metrics",
            "sourceId": "src_headline_summary",
            "metrics": [
                {
                    "label": "k>8 目標區域",
                    "field": "target_components",
                    "format": "number",
                }
            ],
        },
        {
            "id": "card_target_sites",
            "description": "所有 k>8 目標區域包含的 sSNV 位點總數。",
            "dataset": "headline_metrics",
            "sourceId": "src_headline_summary",
            "metrics": [
                {
                    "label": "目標 sSNV 位點",
                    "field": "target_sites",
                    "format": "number",
                }
            ],
        },
        {
            "id": "card_site_retention",
            "description": (
                "工程 assignment 覆蓋；不代表位點已有 read-supported topology evidence。"
            ),
            "dataset": "headline_metrics",
            "sourceId": "src_headline_summary",
            "metrics": [
                {
                    "label": "工程位點指派率",
                    "field": "new_site_retention_ratio",
                    "format": "percent",
                },
                {
                    "label": "densest-8",
                    "field": "old_site_retention_ratio",
                    "format": "percent",
                },
                {
                    "label": "差值",
                    "field": "site_retention_gain",
                    "format": "percent",
                    "signed": True,
                },
            ],
        },
        {
            "id": "card_read_retention",
            "description": (
                "以同一 BASEQ20 extraction 與 component-level molecule-count-"
                "incidence 分母比較新法與 densest-8 counterfactual。"
            ),
            "dataset": "headline_metrics",
            "sourceId": "src_headline_summary",
            "metrics": [
                {
                    "label": "新法 read retention",
                    "field": "new_read_retention_ratio",
                    "format": "percent",
                },
                {
                    "label": "densest-8 counterfactual",
                    "field": "old_read_retention_ratio",
                    "format": "percent",
                },
                {
                    "label": "差值",
                    "field": "read_retention_gain",
                    "format": "percent",
                    "signed": True,
                },
            ],
        },
        {
            "id": "card_primary_active_sites",
            "description": "至少參與 primary linkage 的目標位點；比 assignment 更接近可用證據。",
            "dataset": "headline_metrics",
            "sourceId": "src_headline_summary",
            "metrics": [
                {
                    "label": "Primary-active sites",
                    "field": "primary_active_sites",
                    "format": "number",
                },
                {
                    "label": "占目標位點",
                    "field": "primary_active_site_ratio",
                    "format": "percent",
                },
            ],
        },
        {
            "id": "card_tree_ready_blocks",
            "description": (
                "同時要求 retained pattern、raw weight 與至少兩個 primary-active sites。"
            ),
            "dataset": "headline_metrics",
            "sourceId": "src_headline_summary",
            "metrics": [
                {
                    "label": "Tree-ready local blocks",
                    "field": "tree_ready_local_blocks",
                    "format": "number",
                },
                {
                    "label": "ABSTAIN blocks",
                    "field": "abstain_blocks",
                    "format": "number",
                },
                {
                    "label": "ABSTAIN sites",
                    "field": "abstain_block_sites",
                    "format": "number",
                },
            ],
        },
        {
            "id": "card_zero_evidence",
            "description": "Retained raw-molecule weight=0 的 ABSTAIN 子集。",
            "dataset": "headline_metrics",
            "sourceId": "src_headline_summary",
            "metrics": [
                {
                    "label": "Zero-evidence blocks",
                    "field": "zero_evidence_blocks",
                    "format": "number",
                },
                {
                    "label": "涉及位點",
                    "field": "zero_evidence_block_sites",
                    "format": "number",
                },
            ],
        },
        {
            "id": "card_abstain_components",
            "description": "Component status 以 ABSTAIN 開頭者；不可由 100% assignment 覆蓋掉。",
            "dataset": "headline_metrics",
            "sourceId": "src_headline_summary",
            "metrics": [
                {
                    "label": "ABSTAIN components",
                    "field": "abstain_components",
                    "format": "number",
                },
                {
                    "label": "占 k>8 區域",
                    "field": "abstain_component_ratio",
                    "format": "percent",
                },
                {
                    "label": "Evidence-positive",
                    "field": "evidence_positive_components",
                    "format": "number",
                },
            ],
        },
        {
            "id": "card_runtime",
            "description": (
                "新流程 wall 與舊 historical wall proxy 的 operational 比較；"
                "兩者不是同質 benchmark。"
            ),
            "dataset": "headline_metrics",
            "sourceId": "src_runtime",
            "metrics": [
                {
                    "label": "新流程 wall（min）",
                    "field": "report_wall_minutes",
                    "format": "number",
                },
                {
                    "label": "舊 proxy（min）",
                    "field": "baseline_wall_proxy_minutes",
                    "format": "number",
                },
                {
                    "label": "新／舊",
                    "field": "runtime_ratio_vs_baseline_proxy",
                    "format": "number",
                },
                {
                    "label": "Runner overhead（min）",
                    "field": "runner_overhead_minutes",
                    "format": "number",
                },
            ],
        },
    ]
    if exact_log_audit is not None:
        cards.append(
            {
                "id": "card_exact_stability",
                "description": (
                    "Exact integer-product correction; primary raw-molecule cuts remain unchanged."
                ),
                "dataset": "headline_metrics",
                "sourceId": "src_exact_log_audit",
                "metrics": [
                    {
                        "label": "Corrected stable",
                        "field": "corrected_weight_stable_components",
                        "format": "number",
                    },
                    {
                        "label": "Legacy≠exact cuts",
                        "field": "legacy_log_differs_exact",
                        "format": "number",
                    },
                    {
                        "label": "Stability changed",
                        "field": "correction_changed_stability",
                        "format": "number",
                    },
                ],
            }
        )
    if hp_ps_unit_audit is not None:
        cards.extend(
            [
                {
                    "id": "card_hp_ps_coverage",
                    "description": (
                        "Observed constraint coverage only; no synthetic denominator for unobserved components."
                    ),
                    "dataset": "hp_ps_summary_metrics",
                    "sourceId": "src_hp_ps_unit_audit",
                    "metrics": [
                        {
                            "label": "Observed components",
                            "field": "components_with_observed_constraint_units",
                            "format": "number",
                        },
                        {
                            "label": "No observed constraint",
                            "field": "components_without_observed_constraint_units",
                            "format": "number",
                        },
                        {
                            "label": "Eligible units",
                            "field": "eligible_headline_units",
                            "format": "number",
                        },
                    ],
                },
                {
                    "id": "card_hp_ps_retention",
                    "description": (
                        "Molecule×component-incidence weighted aggregate plus unweighted eligible-unit tail."
                    ),
                    "dataset": "hp_ps_summary_metrics",
                    "sourceId": "src_hp_ps_unit_audit",
                    "metrics": [
                        {
                            "label": "Weighted retention",
                            "field": "weighted_retention_ratio",
                            "format": "percent",
                        },
                        {
                            "label": "Eligible median",
                            "field": "eligible_retention_median",
                            "format": "percent",
                        },
                        {
                            "label": "Eligible p10",
                            "field": "eligible_retention_p10",
                            "format": "percent",
                        },
                        {
                            "label": "Eligible <0.5",
                            "field": "eligible_retention_lt_0_5",
                            "format": "number",
                        },
                    ],
                },
            ]
        )

    charts = [
        {
            "id": "chart_read_retention",
            "title": "每條染色體的新法與 densest-8 counterfactual retention",
            "subtitle": (
                "同一 BASEQ20 extraction 與 component-level incidence denominator；"
                "零分母染色體排除。"
            ),
            "intent": "comparison",
            "question": "新 read-supported 分割相對 densest-8 保存多少 read 證據？",
            "rationale": "Grouped bars keep the common denominator comparison explicit at chromosome grain.",
            "comparisonContext": {
                "baseline": "densest-8 counterfactual",
                "denominator": (
                    "component-level HP1/HP2 × known-PS exact-pattern "
                    "molecule-count incidence"
                ),
                "grain": "chromosome × method",
                "unit": "fraction",
            },
            "type": "bar",
            "dataset": "retention_by_chrom_method",
            "sourceId": "src_retention",
            "encodings": {
                "x": {
                    "field": "chrom",
                    "type": "ordinal",
                    "label": "染色體",
                },
                "y": {
                    "field": "retention_ratio",
                    "type": "quantitative",
                    "format": "percent",
                    "label": "Read retention",
                },
                "color": {
                    "field": "method",
                    "type": "nominal",
                    "label": "方法",
                },
                "tooltip": [
                    {"field": "retained_molecule_weight", "type": "quantitative"},
                    {"field": "raw_total_molecule_weight", "type": "quantitative"},
                    {"field": "target_components", "type": "quantitative"},
                    {"field": "target_sites", "type": "quantitative"},
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "legend": {"position": "bottom", "sort": "spec", "title": "方法"},
            "palette": {"kind": "categorical"},
            "settings": {
                "orientation": "vertical",
                "groupMode": "grouped",
                "categoryLabelPolicy": "rotate",
                "sort": "custom",
            },
            "surface": {"surface": "card", "interactiveLegend": True},
        },
        {
            "id": "chart_target_sites",
            "title": "每條染色體的 k>8 目標位點",
            "subtitle": "顯示分割工作量分布，不等同全染色體突變負荷。",
            "intent": "comparison",
            "question": "k>8 目標位點集中在哪些染色體？",
            "rationale": "A single-series bar chart makes chromosome workload concentration easy to compare.",
            "comparisonContext": {
                "grain": "chromosome",
                "unit": "sites",
                "denominator": "none; exact target-site count",
            },
            "type": "bar",
            "dataset": "per_chromosome_metrics",
            "sourceId": "src_per_chrom",
            "encodings": {
                "x": {
                    "field": "chrom",
                    "type": "ordinal",
                    "label": "染色體",
                },
                "y": {
                    "field": "target_sites",
                    "type": "quantitative",
                    "format": "compact",
                    "label": "目標位點",
                },
                "tooltip": [
                    {"field": "stage_detail", "type": "text"},
                    {"field": "target_components", "type": "quantitative"},
                    {"field": "ssnv_sites", "type": "quantitative"},
                    {"field": "new_blocks", "type": "quantitative"},
                ],
            },
            "valueFormat": "compact",
            "layout": "full",
            "labels": {"values": "auto"},
            "palette": {"kind": "sequential", "name": "blue"},
            "settings": {
                "orientation": "vertical",
                "groupMode": "single",
                "categoryLabelPolicy": "rotate",
                "sort": "custom",
            },
            "surface": {"surface": "card"},
        },
        {
            "id": "chart_genomic_components",
            "title": "k>8 components 的全基因組區域位置",
            "subtitle": (
                "x=起點 Mb；y=染色體序號；大小=pre-cap k；顏色=evidence／"
                + (
                    "corrected exact-product stability。"
                    if exact_log_audit is not None
                    else "legacy stability（exact audit absent）。"
                )
            ),
            "intent": "relationship",
            "question": "大型 k>8 components 位於哪些基因組區域，證據狀態如何？",
            "rationale": (
                "A component-grain genomic scatter preserves all regions and makes "
                "chr6/chr16 concentration, high-k bubbles, and zero-evidence status visible."
            ),
            "comparisonContext": {
                "grain": "one k>8 component",
                "unit": "start Mb × chromosome number",
                "denominator": "all authenticated component_all rows",
            },
            "type": "scatter",
            "dataset": "genomic_components",
            "sourceId": "src_genomic_components",
            "encodings": {
                "x": {
                    "field": "start_mb",
                    "type": "quantitative",
                    "format": "number",
                    "label": "Component start（Mb）",
                },
                "y": {
                    "field": "chromosome_number",
                    "type": "quantitative",
                    "format": "number",
                    "label": "Chromosome number",
                },
                "size": {
                    "field": "pre_cap_k",
                    "type": "quantitative",
                    "label": "Pre-cap k",
                },
                "color": {
                    "field": "weight_stability_status",
                    "type": "nominal",
                    "label": "Evidence／weight status",
                },
                "tooltip": [
                    {"field": "evidence_status_detail", "type": "text"},
                    {"field": "component_id", "type": "text"},
                    {"field": "chrom", "type": "text"},
                    {"field": "start1", "type": "quantitative"},
                    {"field": "end1", "type": "quantitative"},
                    {"field": "span_bp", "type": "quantitative"},
                    {"field": "pre_cap_k", "type": "quantitative"},
                    {"field": "primary_active_sites", "type": "quantitative"},
                    {"field": "raw_total_molecule_weight", "type": "quantitative"},
                    {
                        "field": "new_read_retention_ratio",
                        "type": "quantitative",
                        "format": "percent",
                    },
                    {
                        "field": "old_read_retention_ratio",
                        "type": "quantitative",
                        "format": "percent",
                    },
                    {"field": "status", "type": "text"},
                ],
            },
            "layout": "full",
            "legend": {
                "position": "bottom",
                "sort": "spec",
                "title": "Evidence／weight status",
            },
            "palette": {"kind": "categorical"},
            "maxRows": 500,
            "compatibleTypes": ["scatter"],
            "surface": {"surface": "card", "interactiveLegend": True},
        },
        {
            "id": "chart_runtime",
            "title": "每條染色體的 extraction、partition 與內部演算法時間",
            "subtitle": (
                "Pattern load 與 three-weight partition component loop"
                "（含 DP + diagnostics）是 stage 內含子區段，不可重複相加。"
            ),
            "intent": "comparison",
            "question": (
                "每條染色體的時間主要花在 BAM extraction、partition overhead、"
                "pattern aggregation 還是 three-weight partition component loop？"
            ),
            "rationale": (
                "Grouped bars separate runner stage wall from authenticated "
                "partition-child internal timers."
            ),
            "comparisonContext": {
                "grain": "chromosome × stage",
                "unit": "minutes",
                "normalization": "none",
            },
            "type": "bar",
            "dataset": "runtime_by_chrom_stage",
            "sourceId": "src_runtime",
            "encodings": {
                "x": {
                    "field": "chrom",
                    "type": "ordinal",
                    "label": "染色體",
                },
                "y": {
                    "field": "wall_minutes",
                    "type": "quantitative",
                    "format": "number",
                    "label": "Wall time（min）",
                },
                "color": {
                    "field": "stage",
                    "type": "nominal",
                    "label": "階段",
                },
                "tooltip": [
                    {"field": "target_components", "type": "quantitative"},
                    {"field": "target_sites", "type": "quantitative"},
                    {"field": "max_rss_kb", "type": "quantitative"},
                    {"field": "timing_scope", "type": "text"},
                ],
            },
            "valueFormat": "number",
            "unit": "min",
            "layout": "full",
            "legend": {"position": "bottom", "sort": "spec", "title": "階段"},
            "palette": {"kind": "categorical"},
            "settings": {
                "orientation": "vertical",
                "groupMode": "grouped",
                "categoryLabelPolicy": "rotate",
                "sort": "custom",
            },
            "surface": {"surface": "card", "interactiveLegend": True},
        },
    ]
    if hp_ps_unit_audit is not None:
        charts.append(
            {
                "id": "chart_hp_ps_unit_distribution",
                "title": "Observed HP/PS unit retention 分布",
                "subtitle": (
                    "所有 observed units；顏色區分 headline eligibility，低尾另以 unit 表逐列揭露。"
                ),
                "intent": "distribution",
                "question": "Observed HP/PS units 的 retention 分布與低尾有多寬？",
                "rationale": (
                    "Fixed-width bins preserve the complete observed-unit count while "
                    "the adjacent worst-unit table prevents aggregate masking."
                ),
                "comparisonContext": {
                    "grain": "retention bin × eligibility",
                    "unit": "observed HP/PS units",
                    "denominator": "observed_constraint_units_only",
                },
                "type": "bar",
                "dataset": "hp_ps_unit_distribution",
                "sourceId": "src_hp_ps_unit_audit",
                "encodings": {
                    "x": {
                        "field": "retention_bucket",
                        "type": "ordinal",
                        "label": "Retention bucket",
                    },
                    "y": {
                        "field": "unit_count",
                        "type": "quantitative",
                        "format": "number",
                        "label": "Unit count",
                    },
                    "color": {
                        "field": "eligibility_status",
                        "type": "nominal",
                        "label": "Eligibility",
                    },
                    "tooltip": [
                        {"field": "unit_count", "type": "quantitative"},
                        {"field": "eligibility_status", "type": "text"},
                    ],
                },
                "valueFormat": "number",
                "layout": "full",
                "legend": {
                    "position": "bottom",
                    "sort": "spec",
                    "title": "Eligibility",
                },
                "palette": {"kind": "categorical"},
                "settings": {
                    "orientation": "vertical",
                    "groupMode": "grouped",
                    "categoryLabelPolicy": "rotate",
                    "sort": "custom",
                },
                "surface": {"surface": "card", "interactiveLegend": True},
            }
        )
    if span_sensitivity_dataset:
        charts.extend(
            [
                {
                    "id": "chart_span_retention",
                    "title": "No-cap 與 hard-cap read retention",
                    "subtitle": (
                        "共同 component-level molecule-count-incidence denominator；"
                        "50／100／200 kb 為 inclusive span caps。"
                    ),
                    "intent": "comparison",
                    "question": "實體 span cap 對 constraint read retention 有何影響？",
                    "rationale": "Four ordered bars show the direct no-cap versus hard-cap retention sensitivity.",
                    "comparisonContext": {
                        "baseline": "No cap（k≤8 only）",
                        "denominator": "common raw_total_molecule_weight",
                        "grain": "partition policy",
                        "unit": "fraction",
                    },
                    "type": "bar",
                    "dataset": "span_sensitivity",
                    "sourceId": "src_span_sensitivity",
                    "encodings": {
                        "x": {
                            "field": "cap_label",
                            "type": "ordinal",
                            "label": "Span policy",
                        },
                        "y": {
                            "field": "read_retention_ratio",
                            "type": "quantitative",
                            "format": "percent",
                            "label": "Read retention",
                        },
                        "tooltip": [
                            {
                                "field": "raw_retained_molecule_weight",
                                "type": "quantitative",
                            },
                            {
                                "field": "raw_lost_molecule_weight",
                                "type": "quantitative",
                            },
                            {
                                "field": "raw_total_molecule_weight",
                                "type": "quantitative",
                            },
                            {
                                "field": "retention_delta_vs_no_cap",
                                "type": "quantitative",
                                "format": "percent",
                            },
                            {
                                "field": "unavoidable_span_cap_patterns",
                                "type": "quantitative",
                            },
                        ],
                    },
                    "valueFormat": "percent",
                    "layout": "full",
                    "labels": {"values": "all"},
                    "palette": {"kind": "sequential", "name": "blue"},
                    "settings": {
                        "orientation": "vertical",
                        "groupMode": "single",
                        "sort": "custom",
                    },
                    "surface": {"surface": "card"},
                },
                {
                    "id": "chart_span_blocks",
                    "title": "No-cap 與 hard-cap local block 數",
                    "subtitle": "同一組 k>8 components；block 增加是幾何切割，不是 subclone 數。",
                    "intent": "comparison",
                    "question": "Hard span cap 會增加多少 local blocks？",
                    "rationale": "Ordered bars compare partition fragmentation under the four policies.",
                    "comparisonContext": {
                        "baseline": "No cap（k≤8 only）",
                        "grain": "partition policy",
                        "unit": "blocks",
                    },
                    "type": "bar",
                    "dataset": "span_sensitivity",
                    "sourceId": "src_span_sensitivity",
                    "encodings": {
                        "x": {
                            "field": "cap_label",
                            "type": "ordinal",
                            "label": "Span policy",
                        },
                        "y": {
                            "field": "new_blocks",
                            "type": "quantitative",
                            "format": "compact",
                            "label": "Local blocks",
                        },
                        "tooltip": [
                            {
                                "field": "blocks_delta_vs_no_cap",
                                "type": "quantitative",
                            },
                            {"field": "target_components", "type": "quantitative"},
                            {"field": "target_sites", "type": "quantitative"},
                            {
                                "field": "unavoidable_patterns",
                                "type": "quantitative",
                            },
                        ],
                    },
                    "valueFormat": "compact",
                    "layout": "full",
                    "labels": {"values": "all"},
                    "palette": {"kind": "sequential", "name": "orange"},
                    "settings": {
                        "orientation": "vertical",
                        "groupMode": "single",
                        "sort": "custom",
                    },
                    "surface": {"surface": "card"},
                },
                {
                    "id": "chart_span_runtime",
                    "title": "No-cap 與 hard-cap cached partition wall",
                    "subtitle": "Partition-only 分鐘；排除 source validation、BAM extraction 與下游報告。",
                    "intent": "comparison",
                    "question": "各 span policy 的 cached partition operational cost 為何？",
                    "rationale": "Ordered bars compare partition-only wall time without mixing extraction cost.",
                    "comparisonContext": {
                        "baseline": "No-cap k-only partition-stage wall sum",
                        "grain": "partition policy",
                        "unit": "minutes",
                    },
                    "type": "bar",
                    "dataset": "span_sensitivity",
                    "sourceId": "src_span_sensitivity",
                    "encodings": {
                        "x": {
                            "field": "cap_label",
                            "type": "ordinal",
                            "label": "Span policy",
                        },
                        "y": {
                            "field": "cached_partition_wall_minutes",
                            "type": "quantitative",
                            "format": "number",
                            "label": "Partition wall（min）",
                        },
                        "tooltip": [
                            {
                                "field": "cached_partition_wall_seconds",
                                "type": "quantitative",
                            },
                            {
                                "field": "runtime_ratio_vs_no_cap",
                                "type": "quantitative",
                            },
                            {"field": "new_blocks", "type": "quantitative"},
                            {
                                "field": "read_retention_ratio",
                                "type": "quantitative",
                                "format": "percent",
                            },
                        ],
                    },
                    "valueFormat": "number",
                    "unit": "min",
                    "layout": "full",
                    "labels": {"values": "all"},
                    "palette": {"kind": "sequential", "name": "olive"},
                    "settings": {
                        "orientation": "vertical",
                        "groupMode": "single",
                        "sort": "custom",
                    },
                    "surface": {"surface": "card"},
                },
            ]
        )

    tables = [
        {
            "id": "table_per_chrom",
            "title": "染色體層級 metrics",
            "subtitle": (
                "預設按目標位點數由大到小；完整 timing、active-site 與 stability "
                "欄位保留在同一 snapshot dataset。"
            ),
            "dataset": "per_chromosome_metrics",
            "sourceId": "src_per_chrom",
            "defaultSort": {"field": "target_sites", "direction": "desc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "chrom", "label": "染色體", "type": "text"},
                {
                    "field": "target_components",
                    "label": "k>8 區域",
                    "format": "number",
                },
                {"field": "target_sites", "label": "目標位點", "format": "number"},
                {"field": "new_blocks", "label": "新區塊", "format": "number"},
                {
                    "field": "tree_ready_local_blocks",
                    "label": "Tree-ready blocks",
                    "format": "number",
                },
                {
                    "field": "tree_ready_local_sites",
                    "label": "Tree-ready sites",
                    "format": "number",
                },
                {
                    "field": "old_read_retention_ratio",
                    "label": "Densest-8 retention",
                    "format": "percent",
                },
                {
                    "field": "new_read_retention_ratio",
                    "label": "新 read retention",
                    "format": "percent",
                },
            ],
        },
        {
            "id": "table_top_components",
            "title": "Read-retention gain 最大的 components",
            "subtitle": (
                f"Top {len(top_component_rows)}；診斷排序，不是生物學重要性排名。"
                "完整 retained/lost、active-site 與 exact-log 欄位保留在 snapshot。"
            ),
            "dataset": "top_components",
            "sourceId": "src_components",
            "defaultSort": {"field": "read_retention_gain", "direction": "desc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "rank", "label": "Rank", "format": "number"},
                {"field": "chrom", "label": "染色體", "type": "text"},
                {"field": "component_id", "label": "Component", "type": "text"},
                {"field": "pre_cap_k", "label": "原始 k", "format": "number"},
                {"field": "span_bp", "label": "Span（bp）", "format": "compact"},
                {
                    "field": "raw_total_molecule_weight",
                    "label": "Raw weight",
                    "format": "compact",
                },
                {
                    "field": "read_retention_gain",
                    "label": "Gain",
                    "format": "percent",
                    "movement": True,
                },
                {
                    "field": "weight_stable",
                    "label": (
                        "Corrected stable"
                        if exact_log_audit is not None
                        else "Legacy stable"
                    ),
                    "type": "text",
                },
            ],
        },
    ]
    if hp_ps_unit_audit is not None:
        tables.append(
            {
                "id": "table_hp_ps_worst_units",
                "title": "Lowest-retention eligible HP/PS units",
                "subtitle": (
                    f"最多 {MAX_HP_PS_WORST_UNITS} 列；eligible=weight≥20 且 pattern "
                    "rows≥5。完整 retained/cut/unavoidable 與 active-site 欄位保留在 "
                    "snapshot。"
                ),
                "dataset": "hp_ps_worst_units",
                "sourceId": "src_hp_ps_unit_audit",
                "defaultSort": {"field": "retention_ratio", "direction": "asc"},
                "density": "dense",
                "layout": "full",
                "columns": [
                    {"field": "rank", "label": "Rank", "format": "number"},
                    {"field": "chrom", "label": "染色體", "type": "text"},
                    {
                        "field": "component_id",
                        "label": "Legacy component",
                        "type": "text",
                    },
                    {"field": "hp_family", "label": "HP", "type": "text"},
                    {"field": "phase_set", "label": "Exact PS", "type": "text"},
                    {
                        "field": "component_k",
                        "label": "Component k",
                        "format": "number",
                    },
                    {
                        "field": "retention_ratio",
                        "label": "Retention",
                        "format": "percent",
                    },
                    {
                        "field": "total_weight",
                        "label": "Total incidence",
                        "format": "number",
                    },
                ],
            }
        )

    if span_sensitivity_dataset:
        span_blocks = [
            {
                "id": "span_retention_intro",
                "type": "markdown",
                "body": span_retention_body,
                "layout": "full",
                "sourceId": "src_span_sensitivity",
            },
            {
                "id": "span_retention_chart",
                "type": "chart",
                "chartId": "chart_span_retention",
                "layout": "full",
            },
            {
                "id": "span_blocks_intro",
                "type": "markdown",
                "body": span_blocks_body,
                "layout": "full",
                "sourceId": "src_span_sensitivity",
            },
            {
                "id": "span_blocks_chart",
                "type": "chart",
                "chartId": "chart_span_blocks",
                "layout": "full",
            },
            {
                "id": "span_runtime_intro",
                "type": "markdown",
                "body": span_runtime_body,
                "layout": "full",
                "sourceId": "src_span_sensitivity",
            },
            {
                "id": "span_runtime_chart",
                "type": "chart",
                "chartId": "chart_span_runtime",
                "layout": "full",
            },
        ]
    else:
        span_blocks = [
            {
                "id": "span_not_provided",
                "type": "markdown",
                "body": span_retention_body,
                "layout": "full",
            }
        ]
    exact_audit_blocks = [
        {
            "id": "exact_log_audit",
            "type": "markdown",
            "body": exact_audit_body,
            "layout": "full",
            **(
                {"sourceId": "src_exact_log_audit"}
                if exact_log_audit is not None
                else {}
            ),
        }
    ]
    hp_ps_blocks: list[dict[str, Any]] = [
        {
            "id": "hp_ps_unit_intro",
            "type": "markdown",
            "body": hp_ps_body,
            "layout": "full",
            **(
                {"sourceId": "src_hp_ps_unit_audit"}
                if hp_ps_unit_audit is not None
                else {}
            ),
        }
    ]
    if hp_ps_unit_audit is not None:
        hp_ps_blocks.extend(
            [
                {
                    "id": "hp_ps_unit_distribution_chart",
                    "type": "chart",
                    "chartId": "chart_hp_ps_unit_distribution",
                    "layout": "full",
                },
                {
                    "id": "hp_ps_worst_units_intro",
                    "type": "markdown",
                    "body": hp_ps_worst_body,
                    "layout": "full",
                    "sourceId": "src_hp_ps_unit_audit",
                },
                {
                    "id": "hp_ps_worst_units_table",
                    "type": "table",
                    "tableId": "table_hp_ps_worst_units",
                    "layout": "full",
                },
            ]
        )

    blocks = [
        {"id": "title", "type": "markdown", "body": f"# {title}", "layout": "full"},
        {
            "id": "executive_summary",
            "type": "markdown",
            "body": executive_body,
            "layout": "full",
        },
        {
            "id": "headline_metrics",
            "type": "metric-strip",
            "cardIds": [card["id"] for card in cards],
            "layout": "full",
        },
        {
            "id": "concepts",
            "type": "markdown",
            "body": concept_body,
            "layout": "full",
        },
        *exact_audit_blocks,
        *hp_ps_blocks,
        {
            "id": "retention_intro",
            "type": "markdown",
            "body": retention_body,
            "layout": "full",
            "sourceId": "src_retention",
        },
        {
            "id": "retention_chart",
            "type": "chart",
            "chartId": "chart_read_retention",
            "layout": "full",
        },
        {
            "id": "target_intro",
            "type": "markdown",
            "body": target_body,
            "layout": "full",
            "sourceId": "src_per_chrom",
        },
        {
            "id": "target_chart",
            "type": "chart",
            "chartId": "chart_target_sites",
            "layout": "full",
        },
        {
            "id": "genomic_intro",
            "type": "markdown",
            "body": genomic_body,
            "layout": "full",
            "sourceId": "src_genomic_components",
        },
        {
            "id": "genomic_chart",
            "type": "chart",
            "chartId": "chart_genomic_components",
            "layout": "full",
        },
        {
            "id": "runtime_intro",
            "type": "markdown",
            "body": runtime_body,
            "layout": "full",
            "sourceId": "src_runtime",
        },
        {
            "id": "runtime_chart",
            "type": "chart",
            "chartId": "chart_runtime",
            "layout": "full",
        },
        *span_blocks,
        {
            "id": "per_chrom_intro",
            "type": "markdown",
            "body": per_chrom_body,
            "layout": "full",
            "sourceId": "src_per_chrom",
        },
        {
            "id": "per_chrom_table",
            "type": "table",
            "tableId": "table_per_chrom",
            "layout": "full",
        },
        {
            "id": "components_intro",
            "type": "markdown",
            "body": components_body,
            "layout": "full",
            "sourceId": "src_components",
        },
        {
            "id": "components_table",
            "type": "table",
            "tableId": "table_top_components",
            "layout": "full",
        },
        {
            "id": "methodology",
            "type": "markdown",
            "body": methodology_body,
            "layout": "full",
        },
        {
            "id": "caveats",
            "type": "markdown",
            "body": caveat_body,
            "layout": "full",
        },
        {
            "id": "next_steps",
            "type": "markdown",
            "body": next_body,
            "layout": "full",
        },
    ]

    snapshot = {
        "version": 1,
        "generatedAt": generated_at,
        "status": "fixture" if fixture else "ready",
        "datasets": {
            "headline_metrics": headline_metrics,
            "per_chromosome_metrics": per_chrom_dataset,
            "retention_by_chrom_method": retention_dataset,
            "runtime_by_chrom_stage": runtime_dataset,
            "genomic_components": genomic_component_rows,
            "top_components": top_component_rows,
            **(
                {
                    "hp_ps_summary_metrics": hp_ps_summary_dataset,
                    "hp_ps_unit_distribution": hp_ps_unit_distribution,
                    "hp_ps_worst_units": hp_ps_worst_units,
                }
                if hp_ps_unit_audit is not None
                else {}
            ),
            **(
                {"span_sensitivity": span_sensitivity_dataset}
                if span_sensitivity_dataset
                else {}
            ),
        },
    }
    manifest = {
        "version": 1,
        "surface": "report",
        "title": title,
        "description": (
            "HCC1395 k>8 read-supported segmentation 的位點、read retention、"
            "染色體 runtime 與 component 診斷報告。"
        ),
        "generatedAt": generated_at,
        "cards": cards,
        "charts": charts,
        "tables": tables,
        "sources": sources,
        "blocks": blocks,
    }
    artifact = {
        "surface": "report",
        "manifest": manifest,
        "snapshot": snapshot,
        "sources": sources,
    }
    validate_artifact_shape(artifact)
    return artifact


def validate_artifact_shape(artifact: Mapping[str, Any]) -> None:
    """Run local structural checks without invoking the MCP artifact validator."""
    if artifact.get("surface") != "report":
        raise ReportInputError("artifact surface must be report")
    manifest = artifact.get("manifest")
    snapshot = artifact.get("snapshot")
    if not isinstance(manifest, dict) or not isinstance(snapshot, dict):
        raise ReportInputError("artifact manifest/snapshot must be objects")
    title = manifest.get("title")
    blocks = manifest.get("blocks")
    if not isinstance(title, str) or not isinstance(blocks, list) or not blocks:
        raise ReportInputError("manifest requires title and nonempty blocks")
    first = blocks[0]
    if (
        not isinstance(first, dict)
        or first.get("type") != "markdown"
        or first.get("body") != f"# {title}"
    ):
        raise ReportInputError("first report block must be the matching # title")

    datasets = snapshot.get("datasets")
    if not isinstance(datasets, dict):
        raise ReportInputError("snapshot.datasets must be an object")
    for dataset_id, rows in datasets.items():
        if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
            raise ReportInputError(f"snapshot dataset must be a plain row array: {dataset_id}")
        if len(rows) > MAX_ROWS_PER_DATASET:
            raise ReportInputError(
                f"snapshot dataset exceeds {MAX_ROWS_PER_DATASET} rows: {dataset_id}"
            )
    snapshot_json = json.dumps(
        snapshot, ensure_ascii=False, allow_nan=False, separators=(",", ":")
    )
    if len(snapshot_json.encode("utf-8")) > MAX_SNAPSHOT_BYTES:
        raise ReportInputError("snapshot exceeds 3 MB")

    source_ids = {
        source.get("id")
        for source in manifest.get("sources", [])
        if isinstance(source, dict)
    }
    inline_source_chars = 0
    for source in manifest.get("sources", []):
        if not isinstance(source, dict) or not isinstance(source.get("id"), str):
            raise ReportInputError("every manifest source must have an id")
        path = source.get("path")
        if isinstance(path, str) and (Path(path).is_absolute() or ".." in Path(path).parts):
            raise ReportInputError(f"unsafe source path in artifact: {path}")
        inline_source_chars += len(json.dumps(source, ensure_ascii=False))
    if inline_source_chars > MAX_INLINE_SOURCE_CHARS:
        raise ReportInputError("inline source metadata exceeds 200k characters")

    for collection_name in ("cards", "charts", "tables"):
        collection = manifest.get(collection_name, [])
        if not isinstance(collection, list):
            raise ReportInputError(f"manifest.{collection_name} must be an array")
        for item in collection:
            if item.get("dataset") not in datasets:
                raise ReportInputError(
                    f"{collection_name} item references missing dataset: {item.get('id')}"
                )
            if item.get("sourceId") not in source_ids and "source" not in item:
                raise ReportInputError(
                    f"{collection_name} item lacks valid provenance: {item.get('id')}"
                )

    for chart in manifest.get("charts", []):
        encodings = chart.get("encodings")
        if (
            not isinstance(encodings, dict)
            or not isinstance(encodings.get("x"), dict)
            or not isinstance(encodings.get("y"), dict)
            or not encodings["x"].get("field")
            or not (
                encodings["y"].get("field")
                or isinstance(encodings["y"].get("fields"), list)
            )
        ):
            raise ReportInputError(f"chart lacks canonical x/y encodings: {chart['id']}")
    for table in manifest.get("tables", []):
        columns = table.get("columns")
        sort = table.get("defaultSort")
        fields = {
            column.get("field")
            for column in columns
            if isinstance(column, dict) and isinstance(column.get("field"), str)
        }
        if (
            not isinstance(sort, dict)
            or sort.get("field") not in fields
            or sort.get("direction") not in {"asc", "desc"}
        ):
            raise ReportInputError(f"table defaultSort is invalid: {table['id']}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--summary-tsv", required=True, type=Path)
    parser.add_argument("--component-all", required=True, type=Path)
    parser.add_argument("--baseline-runtime-audit", required=True, type=Path)
    parser.add_argument(
        "--span-grid-summary",
        type=Path,
        help=(
            "Optional cached span-sensitivity run root, summary.tsv, or receipt.json; "
            "adds no-cap versus 50/100/200 kb evidence."
        ),
    )
    parser.add_argument(
        "--exact-log-audit",
        type=Path,
        action="append",
        default=[],
        help=(
            "Exact-product remediation root, receipt.json, or summary.json. "
            "Repeat for multiple single-chromosome fixture/probe audits; formal "
            "reports require one comprehensive chr1-chr22 audit."
        ),
    )
    parser.add_argument(
        "--hp-ps-unit-audit",
        type=Path,
        help=(
            "Observed-constraint HP/PS unit audit root, receipt.json, or summary.json. "
            "Formal reports require an authenticated mode=full audit."
        ),
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--fixture",
        action="store_true",
        help="Allow subset input and force a visible non-evidence fixture report.",
    )
    parser.add_argument(
        "--top-components",
        type=int,
        default=25,
        help="Bounded number of component diagnostics to include (default: 25).",
    )
    parser.add_argument(
        "--generated-at",
        help="ISO-8601 timestamp; defaults to current UTC time.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if not 1 <= args.top_components <= MAX_ROWS_PER_DATASET:
            raise ReportInputError(
                f"--top-components must be between 1 and {MAX_ROWS_PER_DATASET}"
            )
        input_paths = {
            "summary_json": require_regular_file(args.summary_json, "summary.json"),
            "summary_tsv": require_regular_file(args.summary_tsv, "summary.tsv"),
            "component_all": require_regular_file(
                args.component_all, "component_all.tsv.gz"
            ),
            "baseline_audit": require_regular_file(
                args.baseline_runtime_audit, "baseline-runtime-audit.md"
            ),
        }
        summary = load_json(input_paths["summary_json"])
        summary_tsv_fields, summary_tsv_rows = read_tsv(input_paths["summary_tsv"])
        component_fields, component_rows = read_tsv(
            input_paths["component_all"], compressed=True
        )
        baseline = parse_baseline_runtime(input_paths["baseline_audit"])
        generated_at = parse_generated_at(args.generated_at)
        per_chrom, totals = validate_inputs(
            summary,
            summary_tsv_rows,
            summary_tsv_fields,
            component_rows,
            component_fields,
            fixture=args.fixture,
        )
        if not args.fixture and not args.exact_log_audit:
            raise ReportInputError(
                "formal report requires --exact-log-audit comprehensive evidence"
            )
        if not args.fixture and args.hp_ps_unit_audit is None:
            raise ReportInputError(
                "formal report requires --hp-ps-unit-audit full evidence"
            )
        exact_log_audit = None
        if args.exact_log_audit:
            exact_log_audit = validate_exact_log_audits(
                args.exact_log_audit,
                fixture=args.fixture,
                component_rows=component_rows,
            )
            for index, paths in enumerate(exact_log_audit["paths"], start=1):
                input_paths[f"exact_receipt_{index}"] = paths["receipt"]
                input_paths[f"exact_receipt_sha_{index}"] = paths["receipt_sha"]
                input_paths[f"exact_summary_{index}"] = paths["summary"]
                input_paths[f"exact_detail_{index}"] = paths["detail"]
                input_paths[f"exact_success_{index}"] = paths["success"]
        hp_ps_unit_audit = None
        if args.hp_ps_unit_audit is not None:
            hp_ps_unit_audit = validate_hp_ps_unit_audit(
                args.hp_ps_unit_audit,
                fixture=args.fixture,
                component_rows=component_rows,
                report_chromosomes=[str(row["chrom"]) for row in per_chrom],
            )
            hp_paths = hp_ps_unit_audit["paths"]
            input_paths["hp_ps_receipt"] = hp_paths["receipt"]
            input_paths["hp_ps_receipt_sha"] = hp_paths["receipt_sha"]
            input_paths["hp_ps_summary"] = hp_paths["summary"]
            input_paths["hp_ps_summary_tsv"] = hp_paths["summary_tsv"]
            input_paths["hp_ps_units"] = hp_paths["units"]
            input_paths["hp_ps_pairs"] = hp_paths["pairs"]
        span_grid = None
        if args.span_grid_summary is not None:
            span_grid = validate_span_grid(
                args.span_grid_summary,
                fixture=args.fixture,
                main_totals=totals,
            )
            input_paths["span_summary"] = span_grid["summary_path"]
            input_paths["span_receipt"] = span_grid["receipt_path"]
        artifact = build_artifact(
            summary=summary,
            per_chrom=per_chrom,
            totals=totals,
            component_rows=component_rows,
            baseline=baseline,
            span_grid=span_grid,
            exact_log_audit=exact_log_audit,
            hp_ps_unit_audit=hp_ps_unit_audit,
            input_paths=input_paths,
            generated_at=generated_at,
            fixture=args.fixture,
            top_components=args.top_components,
        )
        output = args.output
        if output.exists() or output.is_symlink():
            raise ReportInputError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        serialized = (
            json.dumps(
                artifact,
                ensure_ascii=False,
                allow_nan=False,
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        with output.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(serialized)
        receipt = {
            "ok": True,
            "surface": "report",
            "snapshot_status": artifact["snapshot"]["status"],
            "output": str(output),
            "output_sha256": sha256_path(output),
            "dataset_rows": {
                name: len(rows)
                for name, rows in artifact["snapshot"]["datasets"].items()
            },
            "artifact_profile": {
                "datasets": len(artifact["snapshot"]["datasets"]),
                "rows": sum(
                    len(rows) for rows in artifact["snapshot"]["datasets"].values()
                ),
                "blocks": len(artifact["manifest"]["blocks"]),
                "cards": len(artifact["manifest"]["cards"]),
                "charts": len(artifact["manifest"]["charts"]),
                "tables": len(artifact["manifest"]["tables"]),
                "sources": len(artifact["manifest"]["sources"]),
            },
            "portable_contract": {
                "snapshot_datasets_shape": "object_of_plain_row_arrays",
                "max_rows_per_dataset": MAX_ROWS_PER_DATASET,
                "max_snapshot_bytes": MAX_SNAPSHOT_BYTES,
                "first_markdown_matches_title": True,
                "all_tables_have_default_sort": True,
                "exact_log_audit": (
                    "authenticated" if exact_log_audit is not None else "not_provided"
                ),
                "hp_ps_unit_audit": (
                    "authenticated" if hp_ps_unit_audit is not None else "not_provided"
                ),
                "formal_profile_without_optional_span": {
                    "datasets": 9,
                    "row_formula": (
                        "584 + min(top_components,408) + "
                        "min(eligible_headline_units,25)"
                    ),
                    "default_max_rows": 634,
                    "blocks": 25,
                    "charts": 5,
                    "tables": 3,
                },
                "formal_profile_with_optional_span": {
                    "datasets": 10,
                    "row_formula": (
                        "588 + min(top_components,408) + "
                        "min(eligible_headline_units,25)"
                    ),
                    "default_max_rows": 638,
                    "blocks": 30,
                    "charts": 8,
                    "tables": 3,
                },
            },
        }
        print(json.dumps(receipt, ensure_ascii=False, sort_keys=True))
        return 0
    except (ReportInputError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
