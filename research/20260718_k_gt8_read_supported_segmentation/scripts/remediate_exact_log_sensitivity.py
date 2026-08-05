#!/usr/bin/env python3
"""Fail-closed exact-product remediation for legacy log1p sensitivity cuts.

For positive integer molecule support ``m``, maximizing

    sum(log(m + 1))

over retained patterns is mathematically equivalent to maximizing the exact
integer product

    product(m + 1).

The latter comparison cannot be reversed by Decimal approximation or
aggregation order.  This standalone audit never modifies the source
partition.  It reads its immutable TSV outputs, independently recomputes the
raw/equal objectives, and replaces only the *log sensitivity comparison* with
the exact-product objective.  The primary raw-molecule cut is never changed.

The legacy component and constraint tables do not contain coordinates for
inactive sites.  Therefore the sibling site-membership table is a mandatory
coordinate witness: without it, local indices and the genomic-gap tie-break
cannot be reconstructed.  All three files and the source receipt are verified
before any output directory is created.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_NAME = "intersubmod.k_gt8_exact_log_sensitivity_remediation"
SCHEMA_VERSION = "0.1.0"
SOURCE_SCHEMA = "intersubmod.k_gt8_read_supported_segmentation"
FULL_SOURCE_SCHEMA = "intersubmod.hcc1395_full_k_gt8_segmentation"
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))

COMPONENT_REQUIRED = {
    "dataset",
    "chrom",
    "legacy_component_id",
    "start1",
    "end1",
    "pre_cap_k",
    "exact_pattern_count",
    "raw_total_molecule_weight",
    "raw_retained_molecule_weight",
    "raw_lost_molecule_weight",
    "retained_exact_pattern_count",
    "lost_exact_pattern_count",
    "weight_stable",
    "raw_cut_indices",
    "equal_pattern_cut_indices",
    "log1p_cut_indices",
    "positions_sha256",
}
MEMBERSHIP_REQUIRED = {
    "dataset",
    "chrom",
    "legacy_component_id",
    "site_index",
    "component_local_index",
    "pos1",
}
CONSTRAINT_REQUIRED = {
    "dataset",
    "chrom",
    "legacy_component_id",
    "constraint_id",
    "positions",
    "n_fixed_ra",
    "span_sites",
    "molecule_weight",
}

OUTPUT_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "start1",
    "end1",
    "pre_cap_k",
    "exact_pattern_count",
    "raw_cut_indices",
    "equal_pattern_cut_indices",
    "legacy_log1p_cut_indices",
    "exact_product_cut_indices",
    "raw_recomputed_matches_source",
    "equal_recomputed_matches_source",
    "legacy_log_matches_exact",
    "legacy_weight_stable_reported",
    "legacy_weight_stable_reconstructed",
    "corrected_weight_stable",
    "correction_changed_stability",
    "legacy_log_exact_product_relation",
    "legacy_log_retained_pattern_count",
    "exact_product_retained_pattern_count",
    "legacy_log_block_count",
    "exact_product_block_count",
    "legacy_log_cut_gap_sum_bp",
    "exact_product_cut_gap_sum_bp",
    "legacy_log_product_bit_length",
    "exact_product_bit_length",
    "legacy_log_product_sha256",
    "exact_product_sha256",
    "remediation_class",
)


class ContractError(RuntimeError):
    """Raised when an immutable input or conservation contract is violated."""


@dataclass(frozen=True)
class Constraint:
    """One aggregated read-pattern constraint."""

    constraint_id: str
    lo: int
    hi: int
    molecule_weight: int


@dataclass(frozen=True)
class State:
    """One DP prefix state for one objective."""

    objective: int
    retained_patterns: int
    block_count: int
    gap_sum: int
    cuts: tuple[int, ...]


@dataclass(frozen=True)
class ObjectiveResult:
    """Optimal cuts and the associated exact objective diagnostics."""

    mode: str
    cuts: tuple[int, ...]
    objective: int
    retained_patterns: int
    block_count: int
    gap_sum: int


@dataclass(frozen=True)
class SourceBundle:
    """Verified path set for one chromosome partition."""

    components: Path
    constraints: Path
    membership: Path
    receipt: Path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def semantic_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise ContractError(f"required regular file is missing: {path}")
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def strict_json_load(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"invalid JSON: {path}: {exc}") from exc


def verify_sha256_sidecar(path: Path) -> str:
    sidecar = Path(f"{path}.sha256")
    if not sidecar.is_file() or sidecar.is_symlink():
        raise ContractError(f"missing SHA-256 sidecar: {sidecar}")
    tokens = sidecar.read_text(encoding="utf-8").strip().split()
    if len(tokens) != 2 or tokens[1] != path.name:
        raise ContractError(f"malformed SHA-256 sidecar: {sidecar}")
    observed = sha256_path(path)
    if tokens[0] != observed:
        raise ContractError(
            f"SHA-256 sidecar mismatch for {path}: "
            f"expected={tokens[0]} observed={observed}"
        )
    return observed


def parse_nonnegative_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{label} is not an integer: {value!r}") from exc
    if parsed < 0:
        raise ContractError(f"{label} must be non-negative: {parsed}")
    return parsed


def parse_positive_int(value: str, label: str) -> int:
    parsed = parse_nonnegative_int(value, label)
    if parsed < 1:
        raise ContractError(f"{label} must be positive: {parsed}")
    return parsed


def parse_bool(value: str, label: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ContractError(f"{label} must be lowercase true/false: {value!r}")


def parse_csv_ints(value: str, label: str) -> tuple[int, ...]:
    if value == "":
        return ()
    try:
        parsed = tuple(int(token) for token in value.split(","))
    except ValueError as exc:
        raise ContractError(f"{label} is not a comma-separated integer list") from exc
    return parsed


def format_cuts(cuts: Sequence[int]) -> str:
    return ",".join(map(str, cuts))


def read_tsv_gz(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ContractError(f"required gzip TSV is missing: {path}")
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                missing = sorted(required - set(reader.fieldnames or ()))
                raise ContractError(f"TSV missing columns {missing}: {path}")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ContractError(f"cannot read gzip TSV {path}: {exc}") from exc
    return rows


def _state_is_better(candidate: State, incumbent: State | None) -> bool:
    if incumbent is None:
        return True
    candidate_key = (
        candidate.objective,
        candidate.retained_patterns,
        -candidate.block_count,
        candidate.gap_sum,
    )
    incumbent_key = (
        incumbent.objective,
        incumbent.retained_patterns,
        -incumbent.block_count,
        incumbent.gap_sum,
    )
    if candidate_key != incumbent_key:
        return candidate_key > incumbent_key
    return candidate.cuts < incumbent.cuts


def validate_cuts(cuts: Sequence[int], n_sites: int, max_block_size: int) -> None:
    normalized = tuple(cuts)
    if any(isinstance(cut, bool) or not isinstance(cut, int) for cut in normalized):
        raise ContractError("cut indices must be integers")
    if tuple(sorted(set(normalized))) != normalized:
        raise ContractError(f"cut indices must be strictly increasing: {normalized}")
    if any(cut <= 0 or cut >= n_sites for cut in normalized):
        raise ContractError(f"cut index outside component: n={n_sites} cuts={normalized}")
    bounds = (0,) + normalized + (n_sites,)
    if any(
        right - left > max_block_size
        for left, right in zip(bounds, bounds[1:])
    ):
        raise ContractError(
            f"cuts violate max_block_size={max_block_size}: "
            f"n={n_sites} cuts={normalized}"
        )


def evaluate_cuts(
    positions: Sequence[int],
    constraints: Sequence[Constraint],
    cuts: Sequence[int],
    *,
    max_block_size: int,
) -> dict[str, int]:
    """Evaluate all exact objectives for one valid cut tuple."""

    validate_cuts(cuts, len(positions), max_block_size)
    raw = 0
    patterns = 0
    product = 1
    for constraint in constraints:
        retained = not any(constraint.lo < cut <= constraint.hi for cut in cuts)
        if retained:
            raw += constraint.molecule_weight
            patterns += 1
            product *= constraint.molecule_weight + 1
    return {
        "raw": raw,
        "patterns": patterns,
        "product": product,
        "block_count": len(cuts) + 1,
        "gap_sum": sum(
            positions[cut] - positions[cut - 1] for cut in cuts
        ),
    }


def solve_ordered_partition(
    positions: Sequence[int],
    constraints: Sequence[Constraint],
    *,
    max_block_size: int,
    mode: str,
) -> ObjectiveResult:
    """Run exact ordered contiguous bounded DP for one integer objective."""

    if mode not in {"raw", "equal", "exact_product"}:
        raise ValueError(f"unknown objective mode: {mode}")
    ordered_positions = tuple(positions)
    if not ordered_positions:
        raise ContractError("positions must not be empty")
    if any(
        right <= left
        for left, right in zip(ordered_positions, ordered_positions[1:])
    ):
        raise ContractError("positions must be strictly increasing")
    if max_block_size < 1:
        raise ContractError("max_block_size must be positive")

    n_sites = len(ordered_positions)
    interval_raw: dict[tuple[int, int], int] = defaultdict(int)
    interval_patterns: dict[tuple[int, int], int] = defaultdict(int)
    interval_product: dict[tuple[int, int], int] = defaultdict(lambda: 1)
    seen_ids: set[str] = set()
    for constraint in constraints:
        if constraint.constraint_id in seen_ids:
            raise ContractError(f"duplicate constraint ID: {constraint.constraint_id}")
        seen_ids.add(constraint.constraint_id)
        if not (0 <= constraint.lo <= constraint.hi < n_sites):
            raise ContractError(
                f"constraint interval outside positions: {constraint}"
            )
        if constraint.molecule_weight < 1:
            raise ContractError(
                f"constraint support must be positive: {constraint}"
            )
        if constraint.hi - constraint.lo + 1 > max_block_size:
            continue
        key = (constraint.lo, constraint.hi)
        interval_raw[key] += constraint.molecule_weight
        interval_patterns[key] += 1
        interval_product[key] *= constraint.molecule_weight + 1

    block_gain: dict[tuple[int, int], tuple[int, int, int]] = {}
    for start in range(n_sites):
        raw = 0
        patterns = 0
        product = 1
        for end_inclusive in range(
            start, min(n_sites, start + max_block_size)
        ):
            for lo in range(start, end_inclusive + 1):
                key = (lo, end_inclusive)
                raw += interval_raw[key]
                patterns += interval_patterns[key]
                product *= interval_product[key]
            block_gain[(start, end_inclusive + 1)] = (
                raw,
                patterns,
                product,
            )

    initial_objective = 1 if mode == "exact_product" else 0
    states: list[State | None] = [None] * (n_sites + 1)
    states[0] = State(initial_objective, 0, 0, 0, ())
    for end_exclusive in range(1, n_sites + 1):
        best: State | None = None
        for start in range(
            max(0, end_exclusive - max_block_size), end_exclusive
        ):
            previous = states[start]
            if previous is None:
                raise AssertionError("missing DP prefix state")
            raw_gain, pattern_gain, product_gain = block_gain[
                (start, end_exclusive)
            ]
            if mode == "raw":
                objective = previous.objective + raw_gain
            elif mode == "equal":
                objective = previous.objective + pattern_gain
            else:
                objective = previous.objective * product_gain
            if start:
                cuts = previous.cuts + (start,)
                gap_sum = (
                    previous.gap_sum
                    + ordered_positions[start]
                    - ordered_positions[start - 1]
                )
            else:
                cuts = previous.cuts
                gap_sum = previous.gap_sum
            candidate = State(
                objective=objective,
                retained_patterns=previous.retained_patterns + pattern_gain,
                block_count=previous.block_count + 1,
                gap_sum=gap_sum,
                cuts=cuts,
            )
            if _state_is_better(candidate, best):
                best = candidate
        states[end_exclusive] = best

    final = states[n_sites]
    if final is None:
        raise AssertionError("DP produced no terminal state")
    evaluation = evaluate_cuts(
        ordered_positions,
        constraints,
        final.cuts,
        max_block_size=max_block_size,
    )
    expected_objective = {
        "raw": evaluation["raw"],
        "equal": evaluation["patterns"],
        "exact_product": evaluation["product"],
    }[mode]
    if final.objective != expected_objective:
        raise AssertionError(
            f"DP objective does not match disposition audit: {mode}"
        )
    if final.retained_patterns != evaluation["patterns"]:
        raise AssertionError("DP pattern count does not match disposition audit")
    return ObjectiveResult(
        mode=mode,
        cuts=final.cuts,
        objective=final.objective,
        retained_patterns=final.retained_patterns,
        block_count=final.block_count,
        gap_sum=final.gap_sum,
    )


def product_sha256(value: int) -> str:
    return hashlib.sha256(str(value).encode("ascii")).hexdigest()


def relation(left: int, right: int) -> str:
    if left < right:
        return "LT"
    if left > right:
        return "GT"
    return "EQ"


def remediation_class(
    legacy: Mapping[str, int],
    exact: ObjectiveResult,
    legacy_cuts: Sequence[int],
) -> str:
    if tuple(legacy_cuts) == exact.cuts:
        return "LEGACY_LOG_MATCHES_EXACT"
    if legacy["product"] < exact.objective:
        return "LEGACY_LOG_SUBOPTIMAL_EXACT_PRODUCT"
    if legacy["patterns"] < exact.retained_patterns:
        return "LEGACY_LOG_SUBOPTIMAL_PATTERN_TIEBREAK"
    if legacy["block_count"] > exact.block_count:
        return "LEGACY_LOG_SUBOPTIMAL_BLOCK_COUNT_TIEBREAK"
    if legacy["gap_sum"] < exact.gap_sum:
        return "LEGACY_LOG_SUBOPTIMAL_GAP_TIEBREAK"
    if tuple(legacy_cuts) > exact.cuts:
        return "LEGACY_LOG_SUBOPTIMAL_LEX_TIEBREAK"
    raise ContractError(
        "legacy log cut differs from exact optimum but is not worse under "
        "the declared objective/tie-break"
    )


def verify_source_bundle(bundle: SourceBundle) -> tuple[dict[str, Any], int]:
    receipt_file_identity = file_identity(bundle.receipt)
    receipt_sha = verify_sha256_sidecar(bundle.receipt)
    receipt = strict_json_load(bundle.receipt)
    if not isinstance(receipt, Mapping):
        raise ContractError(f"partition receipt is not an object: {bundle.receipt}")
    if receipt.get("schema_name") != SOURCE_SCHEMA:
        raise ContractError(f"unexpected source schema: {bundle.receipt}")
    if receipt.get("all_pass") is not True:
        raise ContractError(f"source partition is not all_pass: {bundle.receipt}")
    parameters = receipt.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ContractError("source receipt parameters are missing")
    max_block_size = parameters.get("max_block_size")
    if max_block_size != 8:
        raise ContractError(
            f"exact-log remediation expects source max_block_size=8: "
            f"observed={max_block_size}"
        )
    expected_tie = [
        "retained_weight_desc",
        "retained_pattern_count_desc",
        "block_count_asc",
        "cut_gap_sum_desc",
        "cut_indices_lexicographic_asc",
    ]
    if parameters.get("tie_break") != expected_tie:
        raise ContractError("source tie-break contract differs")
    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping):
        raise ContractError("source receipt outputs are missing")
    expected_paths = {
        "legacy_components": bundle.components,
        "cut_constraints": bundle.constraints,
        "site_membership": bundle.membership,
    }
    for name, path in expected_paths.items():
        expected = outputs.get(name)
        if not isinstance(expected, Mapping):
            raise ContractError(f"source output identity is missing: {name}")
        observed = file_identity(path)
        if dict(expected) != observed:
            raise ContractError(
                f"source output identity mismatch for {name}: "
                f"expected={expected} observed={observed}"
            )
    return {
        "path": receipt_file_identity["path"],
        "size_bytes": receipt_file_identity["size_bytes"],
        "sha256": receipt_sha,
    }, int(max_block_size)


def infer_single_bundle(args: argparse.Namespace) -> SourceBundle:
    if args.components is None or args.cut_constraints is None:
        raise ContractError(
            "single-chromosome mode requires --components and --cut-constraints"
        )
    components = args.components.resolve()
    constraints = args.cut_constraints.resolve()
    suffix = ".legacy_components.tsv.gz"
    if not components.name.endswith(suffix):
        raise ContractError(
            f"component filename must end with {suffix}: {components}"
        )
    prefix = components.name[: -len(suffix)]
    expected_constraints = components.parent / f"{prefix}.cut_constraints.tsv.gz"
    if constraints != expected_constraints.resolve():
        raise ContractError(
            "components and constraints are not the matching sibling pair"
        )
    membership = (
        args.site_membership.resolve()
        if args.site_membership is not None
        else (components.parent / f"{prefix}.site_membership.tsv.gz").resolve()
    )
    receipt = (
        args.partition_receipt.resolve()
        if args.partition_receipt is not None
        else (components.parent / "receipt.json").resolve()
    )
    return SourceBundle(components, constraints, membership, receipt)


def verify_full_root(full_root: Path) -> dict[str, Any]:
    marker_path = full_root / "_SUCCESS"
    marker_identity = file_identity(marker_path)
    marker = strict_json_load(marker_path)
    if not isinstance(marker, Mapping):
        raise ContractError("full-root _SUCCESS is not an object")
    if marker.get("all_pass") is not True or marker.get("comprehensive") is not True:
        raise ContractError("full-root _SUCCESS is not comprehensive all_pass")
    if tuple(marker.get("scope", ())) != AUTOSOMES:
        raise ContractError("full-root _SUCCESS scope is not chr1-22")
    receipt_path = full_root / "receipt.json"
    receipt_file_identity = file_identity(receipt_path)
    receipt_sha = verify_sha256_sidecar(receipt_path)
    receipt = strict_json_load(receipt_path)
    if not isinstance(receipt, Mapping):
        raise ContractError("full-root receipt is not an object")
    if receipt.get("schema_name") != f"{FULL_SOURCE_SCHEMA}.run_receipt":
        raise ContractError("unexpected full-root receipt schema")
    if (
        receipt.get("all_pass") is not True
        or receipt.get("comprehensive_all_pass") is not True
    ):
        raise ContractError("full-root receipt is not comprehensive all_pass")
    run_identity = marker.get("run_receipt")
    expected_identity = {
        "path": str(receipt_path.resolve()),
        "sha256": receipt_sha,
    }
    if run_identity != expected_identity:
        raise ContractError("_SUCCESS does not identify the verified run receipt")
    totals = receipt.get("totals")
    if not isinstance(totals, Mapping):
        raise ContractError("full-root totals are missing")
    expected_totals = {
        "chromosomes": 22,
        "ssnv_sites": 79_687,
        "k_gt8_components": 408,
        "k_gt8_sites": 47_570,
        "k_gt8_max_k": 3_574,
        "partitioned_chromosomes": 21,
        "zero_target_skipped_chromosomes": 1,
    }
    if dict(totals) != expected_totals:
        raise ContractError(
            f"full-root totals differ: expected={expected_totals} observed={totals}"
        )

    chromosome_records = receipt.get("chromosomes")
    if (
        not isinstance(chromosome_records, list)
        or tuple(
            record.get("chrom") if isinstance(record, Mapping) else None
            for record in chromosome_records
        )
        != AUTOSOMES
    ):
        raise ContractError("full-root chromosome receipt chain is incomplete")
    for record in chromosome_records:
        chrom = str(record["chrom"])
        partition_status = record.get("partition_status")
        partition_identity = record.get("partition")
        if not isinstance(partition_identity, Mapping):
            raise ContractError(f"partition stage identity is missing: {chrom}")
        stage_path = Path(str(partition_identity.get("path", "")))
        observed_stage_identity = file_identity(stage_path)
        verify_sha256_sidecar(stage_path)
        if dict(partition_identity) != observed_stage_identity:
            raise ContractError(f"partition stage identity differs: {chrom}")
        expected_parent = full_root / "chromosomes" / chrom
        if stage_path.resolve().parent != expected_parent.resolve():
            raise ContractError(f"partition stage receipt escaped scope: {chrom}")
        stage_receipt = strict_json_load(stage_path)
        if not isinstance(stage_receipt, Mapping):
            raise ContractError(f"partition stage receipt is not an object: {chrom}")
        if stage_receipt.get("all_pass") is not True:
            raise ContractError(f"partition stage is not all_pass: {chrom}")
        if chrom == "chr21":
            if (
                partition_status != "SKIP_NO_K_GT8_TARGET"
                or stage_receipt.get("status") != "SKIP_NO_K_GT8_TARGET"
            ):
                raise ContractError("chr21 is not the authenticated zero-target skip")
            continue
        if (
            partition_status != "COMPLETED"
            or stage_receipt.get("stage") != "partition"
            or stage_receipt.get("status") != "COMPLETED"
            or stage_receipt.get("chrom") != chrom
        ):
            raise ContractError(f"partition stage contract differs: {chrom}")
        child_path = expected_parent / "partition" / "receipt.json"
        child_file_identity = file_identity(child_path)
        child_sha = verify_sha256_sidecar(child_path)
        expected_child_identity = {
            "path": child_file_identity["path"],
            "size_bytes": child_file_identity["size_bytes"],
            "sha256": child_sha,
        }
        if stage_receipt.get("child_receipt") != expected_child_identity:
            raise ContractError(
                f"partition stage does not authenticate current child: {chrom}"
            )
    return {
        "marker": marker_identity,
        "receipt": {
            "path": receipt_file_identity["path"],
            "size_bytes": receipt_file_identity["size_bytes"],
            "sha256": receipt_sha,
        },
        "totals": expected_totals,
    }


def full_root_bundles(full_root: Path) -> tuple[list[SourceBundle], dict[str, Any]]:
    full_identity = verify_full_root(full_root)
    bundles = []
    for chrom in AUTOSOMES:
        partition = full_root / "chromosomes" / chrom / "partition"
        if not partition.is_dir():
            if chrom == "chr21":
                continue
            raise ContractError(f"partition directory is missing for {chrom}")
        prefix = f"HCC1395.{chrom}"
        bundles.append(
            SourceBundle(
                partition / f"{prefix}.legacy_components.tsv.gz",
                partition / f"{prefix}.cut_constraints.tsv.gz",
                partition / f"{prefix}.site_membership.tsv.gz",
                partition / "receipt.json",
            )
        )
    return bundles, full_identity


def process_bundle(
    bundle: SourceBundle,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    receipt_identity, max_block_size = verify_source_bundle(bundle)
    component_rows = read_tsv_gz(bundle.components, COMPONENT_REQUIRED)
    membership_rows = read_tsv_gz(bundle.membership, MEMBERSHIP_REQUIRED)
    constraint_rows = read_tsv_gz(bundle.constraints, CONSTRAINT_REQUIRED)
    if not component_rows:
        raise ContractError(f"source has no target components: {bundle.components}")

    components: dict[str, dict[str, str]] = {}
    source_dataset_chrom: set[tuple[str, str]] = set()
    for row in component_rows:
        component_id = row["legacy_component_id"]
        if not component_id or component_id in components:
            raise ContractError(f"duplicate/empty component ID: {component_id!r}")
        components[component_id] = row
        source_dataset_chrom.add((row["dataset"], row["chrom"]))
    if len(source_dataset_chrom) != 1:
        raise ContractError("component table has multiple dataset/chrom scopes")
    dataset, chrom = next(iter(source_dataset_chrom))

    membership_by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen_site_indices: set[int] = set()
    for row in membership_rows:
        if (row["dataset"], row["chrom"]) != (dataset, chrom):
            raise ContractError("membership dataset/chrom differs from components")
        component_id = row["legacy_component_id"]
        if component_id not in components:
            raise ContractError(f"membership references unknown component: {component_id}")
        site_index = parse_nonnegative_int(row["site_index"], "site_index")
        if site_index in seen_site_indices:
            raise ContractError(f"duplicate site_index in target membership: {site_index}")
        seen_site_indices.add(site_index)
        membership_by_component[component_id].append(row)

    positions_by_component: dict[str, tuple[int, ...]] = {}
    for component_id, component in components.items():
        rows = sorted(
            membership_by_component.get(component_id, ()),
            key=lambda row: parse_nonnegative_int(
                row["component_local_index"], "component_local_index"
            ),
        )
        k = parse_positive_int(component["pre_cap_k"], "pre_cap_k")
        local_indices = tuple(
            parse_nonnegative_int(row["component_local_index"], "component_local_index")
            for row in rows
        )
        if local_indices != tuple(range(k)):
            raise ContractError(
                f"membership does not provide every local index once: {component_id}"
            )
        positions = tuple(
            parse_positive_int(row["pos1"], "membership pos1") for row in rows
        )
        if any(right <= left for left, right in zip(positions, positions[1:])):
            raise ContractError(f"component positions are not increasing: {component_id}")
        if positions[0] != int(component["start1"]) or positions[-1] != int(
            component["end1"]
        ):
            raise ContractError(f"component endpoints differ from membership: {component_id}")
        if semantic_sha256(positions) != component["positions_sha256"]:
            raise ContractError(f"positions SHA differs: {component_id}")
        positions_by_component[component_id] = positions

    constraints_by_component: dict[str, list[Constraint]] = defaultdict(list)
    seen_constraint_ids: set[str] = set()
    for row in constraint_rows:
        if (row["dataset"], row["chrom"]) != (dataset, chrom):
            raise ContractError("constraint dataset/chrom differs from components")
        component_id = row["legacy_component_id"]
        if component_id not in components:
            raise ContractError(f"constraint references unknown component: {component_id}")
        constraint_id = row["constraint_id"]
        if not constraint_id or constraint_id in seen_constraint_ids:
            raise ContractError(f"duplicate/empty constraint ID: {constraint_id!r}")
        seen_constraint_ids.add(constraint_id)
        constraint_positions = parse_csv_ints(
            row["positions"], f"constraint positions {constraint_id}"
        )
        if (
            len(constraint_positions) < 2
            or tuple(sorted(set(constraint_positions))) != constraint_positions
        ):
            raise ContractError(
                f"constraint positions must be >=2, sorted, unique: {constraint_id}"
            )
        positions = positions_by_component[component_id]
        position_index = {position: index for index, position in enumerate(positions)}
        unknown = tuple(
            position
            for position in constraint_positions
            if position not in position_index
        )
        if unknown:
            raise ContractError(
                f"constraint references sites outside component: "
                f"{constraint_id}: {unknown}"
            )
        indices = tuple(position_index[position] for position in constraint_positions)
        lo, hi = min(indices), max(indices)
        span_sites = parse_positive_int(row["span_sites"], "span_sites")
        if span_sites != hi - lo + 1:
            raise ContractError(f"constraint span_sites differs: {constraint_id}")
        if parse_positive_int(row["n_fixed_ra"], "n_fixed_ra") != len(
            constraint_positions
        ):
            raise ContractError(f"constraint n_fixed_ra differs: {constraint_id}")
        weight = parse_positive_int(row["molecule_weight"], "molecule_weight")
        constraints_by_component[component_id].append(
            Constraint(constraint_id, lo, hi, weight)
        )

    output_rows = []
    counters: Counter[str] = Counter()
    for component_id, component in components.items():
        positions = positions_by_component[component_id]
        constraints = tuple(
            sorted(
                constraints_by_component.get(component_id, ()),
                key=lambda item: item.constraint_id,
            )
        )
        k = len(positions)
        source_pattern_count = parse_nonnegative_int(
            component["exact_pattern_count"], "exact_pattern_count"
        )
        source_total_weight = parse_nonnegative_int(
            component["raw_total_molecule_weight"], "raw_total_molecule_weight"
        )
        if len(constraints) != source_pattern_count:
            raise ContractError(f"constraint count differs: {component_id}")
        if sum(item.molecule_weight for item in constraints) != source_total_weight:
            raise ContractError(f"constraint molecule mass differs: {component_id}")

        raw_source = parse_csv_ints(component["raw_cut_indices"], "raw cuts")
        equal_source = parse_csv_ints(
            component["equal_pattern_cut_indices"], "equal cuts"
        )
        legacy_log = parse_csv_ints(component["log1p_cut_indices"], "log1p cuts")
        for cuts in (raw_source, equal_source, legacy_log):
            validate_cuts(cuts, k, max_block_size)

        raw_result = solve_ordered_partition(
            positions, constraints, max_block_size=max_block_size, mode="raw"
        )
        equal_result = solve_ordered_partition(
            positions, constraints, max_block_size=max_block_size, mode="equal"
        )
        exact_result = solve_ordered_partition(
            positions,
            constraints,
            max_block_size=max_block_size,
            mode="exact_product",
        )
        raw_matches = raw_result.cuts == raw_source
        equal_matches = equal_result.cuts == equal_source
        if not raw_matches or not equal_matches:
            raise ContractError(
                f"independent primary objective differs from source: {component_id}; "
                f"raw={raw_result.cuts}/{raw_source} "
                f"equal={equal_result.cuts}/{equal_source}"
            )

        raw_evaluation = evaluate_cuts(
            positions,
            constraints,
            raw_source,
            max_block_size=max_block_size,
        )
        if raw_evaluation["raw"] != parse_nonnegative_int(
            component["raw_retained_molecule_weight"],
            "raw_retained_molecule_weight",
        ):
            raise ContractError(f"raw retained mass differs: {component_id}")
        if source_total_weight - raw_evaluation["raw"] != parse_nonnegative_int(
            component["raw_lost_molecule_weight"], "raw_lost_molecule_weight"
        ):
            raise ContractError(f"raw lost mass differs: {component_id}")
        if raw_evaluation["patterns"] != parse_nonnegative_int(
            component["retained_exact_pattern_count"],
            "retained_exact_pattern_count",
        ):
            raise ContractError(f"raw retained pattern count differs: {component_id}")
        if source_pattern_count - raw_evaluation["patterns"] != parse_nonnegative_int(
            component["lost_exact_pattern_count"], "lost_exact_pattern_count"
        ):
            raise ContractError(f"raw lost pattern count differs: {component_id}")

        legacy_evaluation = evaluate_cuts(
            positions,
            constraints,
            legacy_log,
            max_block_size=max_block_size,
        )
        exact_evaluation = evaluate_cuts(
            positions,
            constraints,
            exact_result.cuts,
            max_block_size=max_block_size,
        )
        legacy_matches_exact = legacy_log == exact_result.cuts
        legacy_stable_reconstructed = (
            raw_source == equal_source == legacy_log
        )
        legacy_stable_reported = parse_bool(
            component["weight_stable"], "weight_stable"
        )
        if legacy_stable_reconstructed != legacy_stable_reported:
            raise ContractError(f"source weight_stable differs: {component_id}")
        corrected_stable = (
            raw_result.cuts == equal_result.cuts == exact_result.cuts
        )
        class_name = remediation_class(
            legacy_evaluation, exact_result, legacy_log
        )

        output_rows.append(
            {
                "dataset": dataset,
                "chrom": chrom,
                "legacy_component_id": component_id,
                "start1": component["start1"],
                "end1": component["end1"],
                "pre_cap_k": k,
                "exact_pattern_count": source_pattern_count,
                "raw_cut_indices": format_cuts(raw_source),
                "equal_pattern_cut_indices": format_cuts(equal_source),
                "legacy_log1p_cut_indices": format_cuts(legacy_log),
                "exact_product_cut_indices": format_cuts(exact_result.cuts),
                "raw_recomputed_matches_source": "true",
                "equal_recomputed_matches_source": "true",
                "legacy_log_matches_exact": str(legacy_matches_exact).lower(),
                "legacy_weight_stable_reported": str(
                    legacy_stable_reported
                ).lower(),
                "legacy_weight_stable_reconstructed": str(
                    legacy_stable_reconstructed
                ).lower(),
                "corrected_weight_stable": str(corrected_stable).lower(),
                "correction_changed_stability": str(
                    corrected_stable != legacy_stable_reported
                ).lower(),
                "legacy_log_exact_product_relation": relation(
                    legacy_evaluation["product"], exact_result.objective
                ),
                "legacy_log_retained_pattern_count": legacy_evaluation[
                    "patterns"
                ],
                "exact_product_retained_pattern_count": exact_result.retained_patterns,
                "legacy_log_block_count": legacy_evaluation["block_count"],
                "exact_product_block_count": exact_result.block_count,
                "legacy_log_cut_gap_sum_bp": legacy_evaluation["gap_sum"],
                "exact_product_cut_gap_sum_bp": exact_result.gap_sum,
                "legacy_log_product_bit_length": legacy_evaluation[
                    "product"
                ].bit_length(),
                "exact_product_bit_length": exact_result.objective.bit_length(),
                "legacy_log_product_sha256": product_sha256(
                    legacy_evaluation["product"]
                ),
                "exact_product_sha256": product_sha256(exact_result.objective),
                "remediation_class": class_name,
            }
        )
        counters["components"] += 1
        counters["sites"] += k
        counters["patterns"] += source_pattern_count
        counters["legacy_log_matches_exact"] += int(legacy_matches_exact)
        counters["legacy_log_differs_exact"] += int(not legacy_matches_exact)
        counters["legacy_weight_stable"] += int(legacy_stable_reported)
        counters["corrected_weight_stable"] += int(corrected_stable)
        counters["correction_changed_stability"] += int(
            corrected_stable != legacy_stable_reported
        )
        counters[
            f"remediation_class::{class_name}"
        ] += 1

    output_rows.sort(
        key=lambda row: (
            int(row["chrom"][3:]),
            int(row["start1"]),
            row["legacy_component_id"],
        )
    )
    source = {
        "dataset": dataset,
        "chrom": chrom,
        "partition_receipt": receipt_identity,
        "legacy_components": file_identity(bundle.components),
        "cut_constraints": file_identity(bundle.constraints),
        "site_membership_coordinate_witness": file_identity(bundle.membership),
    }
    return output_rows, dict(sorted(counters.items())), source


def deterministic_gzip_writer(path: Path):
    raw = path.open("xb")
    compressed = gzip.GzipFile(
        filename="",
        mode="wb",
        compresslevel=6,
        fileobj=raw,
        mtime=0,
    )
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")

    class Context:
        def __enter__(self):
            return text

        def __exit__(self, exc_type, exc, traceback):
            try:
                text.close()
            finally:
                if not raw.closed:
                    raw.close()
            return False

    return Context()


def write_tsv_gz(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    with deterministic_gzip_writer(path) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json_exclusive(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def write_sha256_sidecar(path: Path) -> Path:
    sidecar = Path(f"{path}.sha256")
    sidecar.write_text(
        f"{sha256_path(path)}  {path.name}\n",
        encoding="utf-8",
    )
    return sidecar


def run(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() or output_dir.is_symlink():
        raise ContractError(f"output directory must not exist: {output_dir}")

    full_identity = None
    if args.full_root is not None:
        if any(
            value is not None
            for value in (
                args.components,
                args.cut_constraints,
                args.site_membership,
                args.partition_receipt,
            )
        ):
            raise ContractError(
                "--full-root cannot be combined with single-chromosome inputs"
            )
        bundles, full_identity = full_root_bundles(args.full_root.resolve())
        mode = "comprehensive_chr1_22"
    else:
        bundles = [infer_single_bundle(args)]
        mode = "single_chromosome"

    rows: list[dict[str, Any]] = []
    aggregate: Counter[str] = Counter()
    source_bundles = []
    chromosome_counts: dict[str, dict[str, int]] = {}
    for bundle in bundles:
        bundle_rows, counts, source = process_bundle(bundle)
        rows.extend(bundle_rows)
        aggregate.update(counts)
        source_bundles.append(source)
        chromosome_counts[source["chrom"]] = counts
    rows.sort(
        key=lambda row: (
            int(row["chrom"][3:]),
            int(row["start1"]),
            row["legacy_component_id"],
        )
    )

    if args.expected_components is not None and len(rows) != args.expected_components:
        raise ContractError(
            f"component total differs: expected={args.expected_components} "
            f"observed={len(rows)}"
        )
    if mode == "comprehensive_chr1_22":
        if len(rows) != 408 or aggregate["sites"] != 47_570:
            raise ContractError(
                "comprehensive source does not conserve canonical "
                f"408 components / 47,570 sites: "
                f"components={len(rows)} sites={aggregate['sites']}"
            )

    semantic_rows = [
        {
            "component_id": row["legacy_component_id"],
            "raw": row["raw_cut_indices"],
            "equal": row["equal_pattern_cut_indices"],
            "legacy_log": row["legacy_log1p_cut_indices"],
            "exact_product": row["exact_product_cut_indices"],
            "corrected_weight_stable": row["corrected_weight_stable"],
            "class": row["remediation_class"],
        }
        for row in rows
    ]
    semantic_digest = semantic_sha256(semantic_rows)
    elapsed_before_write = time.perf_counter() - started
    summary = {
        "schema_name": f"{SCHEMA_NAME}.summary",
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "mode": mode,
        "counts": dict(sorted(aggregate.items())),
        "chromosome_counts": dict(sorted(chromosome_counts.items())),
        "semantic_result_sha256": semantic_digest,
        "timings_seconds": {
            "validate_and_recompute_before_write": elapsed_before_write,
        },
        "interpretation": {
            "primary_raw_partition_changed": False,
            "exact_product_scope": (
                "log1p sensitivity comparison only; exact product(m+1) is "
                "mathematically equivalent to sum(log(m+1))"
            ),
            "claim_ceiling": (
                "exact correction of log-weight cut comparison under the "
                "existing ordered contiguous k<=8 objective; not a tree, "
                "clone, or biological truth validation"
            ),
        },
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir()
    table_path = output_dir / "exact_log_sensitivity.tsv.gz"
    summary_path = output_dir / "summary.json"
    write_tsv_gz(table_path, OUTPUT_FIELDS, rows)
    write_json_exclusive(summary_path, summary)
    summary_sha_path = write_sha256_sidecar(summary_path)

    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "mode": mode,
        "scope": {
            "datasets": sorted({row["dataset"] for row in rows}),
            "chromosomes": sorted(
                {row["chrom"] for row in rows},
                key=lambda chrom: int(chrom[3:]),
            ),
            "components": len(rows),
            "sites": aggregate["sites"],
        },
        "parameters": {
            "max_block_size": 8,
            "remediated_sensitivity_objective": "exact_product(molecule_weight+1)",
            "mathematical_equivalence": (
                "argmax sum(log(m+1)) == argmax product(m+1)"
            ),
            "tie_break": [
                "retained_exact_product_desc",
                "retained_pattern_count_desc",
                "block_count_asc",
                "cut_gap_sum_desc",
                "cut_indices_lexicographic_asc",
            ],
            "primary_raw_partition_changed": False,
        },
        "checks": {
            "source_receipts_and_sha_verified": True,
            "coordinate_witnesses_verified": True,
            "positions_semantic_sha_verified": True,
            "raw_objective_independently_reproduced": True,
            "equal_objective_independently_reproduced": True,
            "constraint_count_and_mass_conserved": True,
            "all_exact_product_cuts_k_lte_8": True,
            "output_directory_was_absent": True,
        },
        "counts": dict(sorted(aggregate.items())),
        "chromosome_counts": dict(sorted(chromosome_counts.items())),
        "sources": {
            "full_root": full_identity,
            "partitions": source_bundles,
        },
        "outputs": {
            "exact_log_sensitivity": file_identity(table_path),
            "summary": file_identity(summary_path),
            "summary_sha256": file_identity(summary_sha_path),
        },
        "script": file_identity(Path(__file__).resolve()),
        "semantic_result_sha256": semantic_digest,
        "timings_seconds": {
            "validate_and_recompute_before_write": elapsed_before_write,
            "total_before_receipt": time.perf_counter() - started,
        },
        "command": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "claim_ceiling": (
            "exact correction of log-weight sensitivity cuts; the source "
            "primary raw-molecule partition is unchanged, and no evolutionary "
            "tree truth claim is made"
        ),
    }
    receipt_path = output_dir / "receipt.json"
    write_json_exclusive(receipt_path, receipt)
    receipt_sha_path = write_sha256_sidecar(receipt_path)
    success = {
        "schema_name": f"{SCHEMA_NAME}.success_marker",
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "receipt": {
            "path": str(receipt_path.resolve()),
            "sha256": sha256_path(receipt_path),
        },
    }
    write_json_exclusive(output_dir / "_SUCCESS", success)
    return {
        "all_pass": True,
        "mode": mode,
        "scope": receipt["scope"],
        "counts": receipt["counts"],
        "timings_seconds": receipt["timings_seconds"],
        "semantic_result_sha256": semantic_digest,
        "receipt": str(receipt_path),
        "receipt_sha256": sha256_path(receipt_path),
        "receipt_sidecar": str(receipt_sha_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Independently remediate legacy Decimal log1p sensitivity cuts "
            "with an exact integer-product DP"
        )
    )
    parser.add_argument("--components", type=Path)
    parser.add_argument("--cut-constraints", type=Path)
    parser.add_argument("--site-membership", type=Path)
    parser.add_argument("--partition-receipt", type=Path)
    parser.add_argument("--full-root", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-components", type=int)
    return parser.parse_args()


def main() -> None:
    try:
        result = run(parse_args())
    except ContractError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
