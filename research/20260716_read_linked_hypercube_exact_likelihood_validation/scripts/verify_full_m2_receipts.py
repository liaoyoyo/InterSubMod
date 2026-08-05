#!/usr/bin/env python3
"""Independently verify the complete M2 extraction/ranking evidence bundle.

This verifier deliberately does not import either full-run orchestrator or the
chromosome ranker.  It reads the persisted child receipts/tables, recomputes
the task index and additive aggregates, and streams an independent
reconstruction of the canonical candidate table.  It also independently
recomputes the BQ-aware mixture objective and its simplex/KKT certificate from
persisted pattern, state, and mixture-weight data.  This separation prevents a
bug in production ranking or aggregation code from automatically certifying
itself.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import stat
import sys
from collections import Counter, defaultdict
from itertools import zip_longest
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
VERIFIER = Path(__file__).resolve()
RELEASE_SCHEMA_NAME = "intersubmod.m2_release_run_manifest"
RELEASE_SCHEMA_VERSION = "1.0.0"
RELEASE_FREEZER_RELATIVE = (
    "research/20260716_read_linked_hypercube_exact_likelihood_validation/"
    "scripts/freeze_m2_release_contract.py"
)
RELEASE_VERIFIER_RELATIVE = (
    "research/20260716_read_linked_hypercube_exact_likelihood_validation/"
    "scripts/verify_full_m2_receipts.py"
)
CANONICAL_MANIFEST_SHA256 = "16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45"
CANONICAL_SCHEMA_SHA256 = "47f54d86159beae66b68719e034c7755eba59b2d1719d4546153fa0ae74ff19f"
CANONICAL_SCHEMA_RELATIVE = (
    "docs/methodology/_assets/20260627_subclone_4axis_teaching/"
    "schemas/layered_input_manifest_v3.schema.json"
)
RELEASE_RESEARCH_RELATIVE = (
    "research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts"
)
RELEASE_SOURCE_PATHS = {
    "extractor": f"{RELEASE_RESEARCH_RELATIVE}/extract_lossless_read_linkage.py",
    "lossless_read_contract": f"{RELEASE_RESEARCH_RELATIVE}/lossless_read_contract.py",
    "full_extraction_runner": f"{RELEASE_RESEARCH_RELATIVE}/run_full_m2_extraction.py",
    "ranker": f"{RELEASE_RESEARCH_RELATIVE}/build_m2_patterns_and_rank.py",
    "hypercube_solver": f"{RELEASE_RESEARCH_RELATIVE}/hypercube_exact.py",
    "full_ranking_runner": f"{RELEASE_RESEARCH_RELATIVE}/run_full_m2_ranking.py",
    "pilot_verifier": f"{RELEASE_RESEARCH_RELATIVE}/verify_m2_single_task_pilot.py",
    "full_verifier": f"{RELEASE_RESEARCH_RELATIVE}/verify_full_m2_receipts.py",
    "input_identity_verifier": f"{RELEASE_RESEARCH_RELATIVE}/verify_m2_input_identity.py",
    "release_contract_freezer": f"{RELEASE_RESEARCH_RELATIVE}/freeze_m2_release_contract.py",
    "canonical_manifest_schema": CANONICAL_SCHEMA_RELATIVE,
}
RELEASE_CHECK_NAMES = {
    "canonical_or_explicit_synthetic_test_authority_selected",
    "supplied_pre_receipt_authenticated_and_semantically_exact",
    "fresh_input_verifier_rerun_passed_and_eligible",
    "fresh_snapshot_exactly_equals_supplied_pre_snapshot",
    "manifest_sha_matches_selected_authority", "exact_eleven_file_allowlist_frozen",
    "all_sources_regular_non_symlink_single_link_and_non_aliased",
    "all_copies_regular_single_link_non_aliased_and_sha_equal_source",
    "canonical_manifest_and_pre_receipt_copied_immutably",
    "publish_boundary_recheck_completed", "bootstrap20_and_seed20260716_frozen",
}
RELEASE_TOP_KEYS = {
    "schema_name", "schema_version", "task_type", "authority_mode",
    "validation_evidence_eligible", "created_at_utc", "scope", "entrypoints",
    "canonical_manifest", "pre_input_identity_receipt",
    "fresh_input_identity_verification", "source_snapshot", "canonical_schema",
    "runtime", "runtime_semantic_sha256", "parameters", "producer", "assurance",
    "checks", "all_pass", "receipt_integrity",
}
RELEASE_FILE_IDENTITY_KEYS = {
    "path", "st_dev", "st_ino", "st_nlink", "size_bytes", "mtime_ns", "ctime_ns",
    "mode_octal", "sha256",
}
FROZEN_RELEASE_PARAMETERS = {
    "extraction": {
        "mapq_min": 20, "baseq_min": 20, "bridge_thresholds": [1, 2, 3, 5],
        "workers": 2, "samtools_threads": 1,
    },
    "ranking": {
        "thresholds": [1, 2, 3, 5], "component_bases": ["PS_HP1", "PS_HP2"],
        "hp_families": ["1", "2"],
        "structural_exact_pattern_minread_grid": [1, 2, 3, 5],
        "primary_structural_exact_pattern_minread": 3, "exact_k_max": 12,
        "max_vertex_sets": 256, "solver_time_limit_seconds_per_milp": 30.0,
        "fixed_error_grid": [0.005, 0.01, 0.02, 0.05],
        "minimum_bq_error_rate": 0.000001, "maximum_bq_error_rate": 0.25,
        "conditional_candidate_ranking_bootstrap_replicates": 20,
        "conditional_candidate_ranking_bootstrap_seed": 20260716,
        "tie_tolerance_log_likelihood": 0.000001, "workers": 2,
    },
    "scheduler": {
        "initial_batch_tasks": 8, "subsequent_batch_tasks": 16,
        "task_timeout_seconds": 28800, "timeout_grace_seconds": 300,
    },
}
RELEASE_ASSURANCE = {
    "code_snapshot": "exact source bytes copied to regular non-aliased files; files 0444; directories 0555",
    "input_identity": "exact manifest-derived PRE plus mandatory fresh verifier snapshot equality",
    "immutable_inputs": "canonical manifest and supplied PRE receipt/sidecar copied read-only into contract",
    "publication": "same-filesystem renameat2(RENAME_NOREPLACE) after complete staging",
}
_DEEP_RELEASE_VERIFY_CACHE: dict[tuple[str, str], dict[str, Any]] = {}
DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
AUTOSOMES = tuple(f"chr{value}" for value in range(1, 23))
LINKAGE_BASES = (
    "pooled", "HP1", "HP2", "HP3", "HP4", "unphased",
    "PS_HP1", "PS_HP2", "MISSING_PS_HP1", "MISSING_PS_HP2",
)
PRIMARY_BASES = ("PS_HP1", "PS_HP2")
REQUIRED_EXTRACTION_SUFFIXES = (
    ".molecule_sparse_calls.tsv.gz",
    ".site_catalog.tsv.gz",
    ".components.tsv.gz",
    ".site_component_membership.tsv.gz",
)
CANDIDATE_SOURCE_NAME = "m2_compressed_vertex_set_candidates.tsv.gz"
PATTERN_SOURCE_NAME = "m2_symbolic_pattern_counts.tsv.gz"
RUNTIME_SOURCE_NAME = "m2_unit_runtime_diagnostics.tsv.gz"
RUNTIME_METRICS = (
    "candidate_generation_elapsed_seconds",
    "likelihood_fit_elapsed_seconds",
    "unit_total_elapsed_seconds",
)
RUNTIME_DIAGNOSTIC_FIELDS = (
    "dataset", "chrom", "threshold", "component_basis", "phase_set", "component_id", "family",
    "structural_exact_pattern_minread", "structural_minread_role",
    "candidate_generation_invoked", "likelihood_fit_invoked", *RUNTIME_METRICS,
)
CANDIDATE_TABLE_COLUMNS = (
    "unit_key", "dataset", "chrom", "component_id", "threshold", "hp_family", "ps",
    "candidate_id", "vertex_set_id", "vertex_states", "vertex_roles", "parent_choice_count",
    "profile_log_likelihood", "relative_log_likelihood", "mixture_weights_pi", "winner_status",
    "tie_group", "coarse_topology_class", "candidate_set_complete",
)
PROFILE_PATTERN_REQUIRED_COLUMNS = (
    "dataset", "chrom", "threshold", "component_basis", "phase_set", "component_id", "family",
    "structural_exact_pattern_minread", "structural_minread_role", "k", "pattern_rax",
    "fixed_base_qualities", "n_molecules", "scoring_eligible",
)
PROFILE_CANDIDATE_REQUIRED_COLUMNS = (
    "dataset", "chrom", "component_basis", "phase_set", "threshold", "component_id", "family",
    "structural_exact_pattern_minread", "vertex_set_id", "states_json", "primary_log_likelihood",
    "relative_likelihood_weight", "mixture_weights_json", "fit_converged", "fit_monotone",
    "optimizer_status", "global_log_likelihood_gap_bound", "simplex_residual", "is_winner",
    "is_tied_winner",
)
EXPECTED_OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE = 1e-8
EXPECTED_SIMPLEX_RESIDUAL_TOLERANCE = 1e-12
ORCHESTRATION_DIRECTORY = "_orchestration"
ORCHESTRATION_SCHEMA_VERSION = "1.1.0"
RESOURCE_GATE_REQUIRED_RESERVE_BYTES = 300 * 1024**3
RESOURCE_GATE_SCHEMA_NAME = "intersubmod.m2_resource_gate_receipt"
RESOURCE_GATE_SCHEMA_VERSION = "1.0.0"
RESOURCE_GATE_RECEIPT_KEYS = {
    "schema_name", "schema_version", "stage", "gate_scope", "gate_id", "gate_nonce",
    "target", "process_snapshot", "filesystem_snapshot", "producer_source",
    "observed_at_utc", "observed_monotonic_ns", "host_boot_id", "checks",
    "receipt_integrity",
}
RESOURCE_GATE_IDENTITY_KEYS = {
    "path", "sha256", "semantic_sha256", "gate_id", "sidecar_path", "sidecar_sha256",
}
ORCHESTRATION_SCHEMA_NAMES = {
    "session": "intersubmod.m2_orchestration_session",
    "batch": "intersubmod.m2_orchestration_batch_start",
    "grant": "intersubmod.m2_orchestration_child_grant",
    "completion": "intersubmod.m2_orchestration_child_completion",
}
ORCHESTRATION_SESSION_KEYS = {
    "schema_name", "schema_version", "stage", "session_id", "session_nonce",
    "created_at_utc", "created_monotonic_ns", "host_boot_id", "release_manifest",
    "release_binding_semantic_sha256", "run_contract_semantic_sha256", "scope",
    "output_root", "producer_sources", "scheduler_policy", "parent_extraction",
    "resource_gate", "receipt_integrity",
}
ORCHESTRATION_BATCH_KEYS = {
    "schema_name", "schema_version", "stage", "session_id", "session_sha256",
    "batch_index", "batch_id", "batch_nonce", "previous_chain_head", "before_count",
    "max_new_tasks", "effective_count", "selected_task_ids",
    "run_contract_semantic_sha256", "runner_source", "resource_gate", "created_at_utc",
    "created_monotonic_ns", "host_boot_id", "receipt_integrity",
}
ORCHESTRATION_GRANT_KEYS = {
    "schema_name", "schema_version", "stage", "session_id", "session_sha256",
    "batch_id", "batch_start_sha256", "task_id", "dataset", "chrom", "task_ordinal",
    "child_nonce", "command", "command_semantic_sha256", "producer_sources",
    "input_identity", "parameters_semantic_sha256", "expected_output_dir",
    "issued_at_utc", "issued_monotonic_ns", "host_boot_id", "receipt_integrity",
}
ORCHESTRATION_COMPLETION_KEYS = {
    "schema_name", "schema_version", "stage", "session_id", "session_sha256",
    "batch_id", "grant_identity", "task_id", "dataset", "chrom", "task_ordinal",
    "child_receipt_identity", "child_outputs_semantic_sha256", "command_semantic_sha256",
    "status", "returncode", "timed_out", "process_group_id", "started_monotonic_ns",
    "completed_at_utc", "completed_monotonic_ns", "host_boot_id", "runner_source",
    "receipt_integrity",
}
ORCHESTRATION_RECEIPT_KEYS = {
    "session_identity", "batch_start_identity", "previous_chain_head",
    "batch_completion_attestations", "cumulative_attested_tasks",
}
ORCHESTRATION_COMPLETED_COUNTS = (8, 24, 40, 56, 72, 88, 104, 120, 136, 152, 154)
POST_INPUT_IDENTITY_SCHEMA_NAME = "intersubmod.m2_input_identity_verification"
POST_INPUT_IDENTITY_SCHEMA_VERSION = "1.0.0"
POST_INPUT_IDENTITY_CHECK_NAMES = {
    "manifest_matches_selected_authority_sha256",
    "schema_matches_selected_authority_sha256",
    "manifest_passes_complete_canonical_v3_json_schema",
    "manifest_scope_is_exactly_seven_datasets",
    "direct_input_artifact_count_is_42",
    "all_42_direct_input_roles_resolve_to_distinct_files",
    "seven_bams_use_storage_identity_v1",
    "thirty_five_non_bam_artifact_roles_use_full_sha256",
    "twenty_one_manifest_defined_bam_chunks_recomputed",
    "all_manifest_identities_match_observed_files",
    "post_snapshot_exactly_matches_authenticated_pre",
}

# Independent literal contract.  Do not import this from the production
# ranker or full-run orchestrator: exact comparison must be capable of finding
# producer-side semantic drift.
EXPECTED_METHOD_CONTRACT = {
    "schema_name": "intersubmod.m2_partial_read_likelihood_method_contract",
    "schema_version": "1.0.0",
    "structural_group_scope": "DISTINCT_EXACT_RAX_COUNT_GE_THRESHOLD",
    "active_compatible_vertex_indices_materialized_for_sparse_rows": True,
    "completion_wise_tree_worlds_materialized": False,
    "cross_read_cartesian_products_materialized": False,
    "likelihood_primary_term": (
        "BQ_SYMMETRIC_SUBSTITUTION_CONDITIONAL_RA_READ_PATTERN_MIXTURE_PROFILE_LIKELIHOOD"
    ),
    "same_molecule_vaf_added_as_separate_term": False,
    "parent_edges_or_transitions_scored": False,
    "missing_calls_marginalized": True,
}
SUM_FIELDS = (
    "n_component_hp_units",
    "n_components",
    "molecule_component_projections",
    "informative_scoring_molecules",
    "all_x_excluded_molecules",
    "structural_retained_molecules",
    "below_minread_scoring_molecules",
    "bq_scoring_molecules",
    "non_ra_cells_marginalized",
    "raw_tree_candidates_T_complete_units",
    "distinct_vertex_sets_V_complete_units",
    "solver_complete_units",
    "solver_incomplete_or_not_run_units",
    "quality_primary_unique_vertex_units",
    "quality_primary_tied_vertex_units",
    "rank_abstain_units",
    "fixed_error_grid_stable_units",
    "fixed_error_grid_evaluated_units",
    "conditional_candidate_ranking_bootstrap_complete_units",
    "conditional_candidate_ranking_bootstrap_not_run_units",
    "conditional_candidate_ranking_bootstrap_nonconverged_units",
    "coarse_topology_class_unique_units",
    "coarse_topology_multiple_class_units",
    "parent_edge_assignment_unique_units",
    "exact_topology_proven_unique_units",
    "topology_evaluated_units",
    "topology_class_inclusion_counts_denominator",
    "k_component_sites_total",
    "k_observed_alt_active_total",
    "k_scoring_alt_observed_total",
    "not_structural_alt_active_sites_total",
    "structural_partial_pattern_groups",
    "partial_group_coverage_denominator",
    "partial_groups_covered",
    "partial_groups_unsatisfied",
)
COUNTER_FIELDS = (
    "selection_status_counts",
    "candidate_generation_status_counts",
    "k_route_counts",
    "projected_call_class_counts",
    "conditional_candidate_ranking_bootstrap_status_counts",
    "topology_class_inclusion_counts",
    "coarse_topology_unique_class_counts",
    "coarse_topology_ambiguous_class_set_counts",
    "topology_derivation_status_counts",
    "exact_topology_uniqueness_status_counts",
)


class VerificationError(RuntimeError):
    """Raised when persisted evidence violates the independent contract."""


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def semantic_json_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def _exact_keys(value: Any, expected: set[str], label: str) -> Mapping[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    require(set(value) == expected, f"{label} exact-key schema mismatch")
    return value


def _nonnegative_int(value: Any, label: str) -> int:
    require(
        isinstance(value, int) and not isinstance(value, bool) and value >= 0,
        f"{label} is not a nonnegative integer",
    )
    return int(value)


def _positive_int(value: Any, label: str) -> int:
    result = _nonnegative_int(value, label)
    require(result > 0, f"{label} is not positive")
    return result


def _hex_sha256(value: Any, label: str) -> str:
    require(
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value),
        f"{label} is not a lowercase SHA-256 hex digest",
    )
    return value


def _nonempty_string(value: Any, label: str) -> str:
    require(isinstance(value, str) and bool(value), f"{label} is not a nonempty string")
    return value


def _strict_json_bytes(payload: bytes, label: Path | str) -> dict[str, Any]:
    def reject_duplicate(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise VerificationError(f"duplicate JSON key {key!r}: {label}")
            result[key] = value
        return result

    try:
        value = json.loads(payload.decode("utf-8"), object_pairs_hook=reject_duplicate)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f"invalid JSON {label}: {exc}") from exc
    require(isinstance(value, dict), f"JSON root is not an object: {label}")
    return value


def _same_stat_state(left: os.stat_result, right: os.stat_result) -> bool:
    return (
        left.st_dev == right.st_dev
        and left.st_ino == right.st_ino
        and left.st_nlink == right.st_nlink
        and left.st_mode == right.st_mode
        and left.st_size == right.st_size
        and left.st_mtime_ns == right.st_mtime_ns
        and left.st_ctime_ns == right.st_ctime_ns
    )


def _stable_regular_bytes(
    path: Path, label: str, *, immutable_mode: int | None = None
) -> tuple[bytes, os.stat_result]:
    raw = path.absolute()
    try:
        before_path = os.lstat(raw)
    except OSError as exc:
        raise VerificationError(f"{label} is unavailable: {exc}") from exc
    require(
        stat.S_ISREG(before_path.st_mode)
        and not stat.S_ISLNK(before_path.st_mode)
        and before_path.st_nlink == 1,
        f"{label} is not a single-link regular non-symlink file",
    )
    if immutable_mode is not None:
        require(
            stat.S_IMODE(before_path.st_mode) == immutable_mode,
            f"{label} mode is not {immutable_mode:04o}",
        )
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(raw, flags)
    except OSError as exc:
        raise VerificationError(f"cannot open {label}: {exc}") from exc
    try:
        before_fd = os.fstat(descriptor)
        require(_same_stat_state(before_path, before_fd), f"{label} changed before read")
        chunks: list[bytes] = []
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            chunks.append(block)
        after_fd = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        after_path = os.lstat(raw)
    except OSError as exc:
        raise VerificationError(f"{label} disappeared after read: {exc}") from exc
    require(
        _same_stat_state(before_fd, after_fd)
        and _same_stat_state(after_fd, after_path),
        f"{label} changed while being authenticated",
    )
    return b"".join(chunks), after_fd


def _authenticate_json_with_sidecar(
    path: Path, label: str, *, immutable: bool = False
) -> tuple[dict[str, Any], dict[str, Any]]:
    mode = 0o444 if immutable else None
    payload, observed = _stable_regular_bytes(path, label, immutable_mode=mode)
    document = _strict_json_bytes(payload, path)
    integrity = _exact_keys(
        document.get("receipt_integrity"),
        {"scheme", "sidecar_name", "covers"},
        f"{label}.receipt_integrity",
    )
    require(
        integrity["scheme"] == "external_sha256_sidecar_v1"
        and integrity["sidecar_name"] == f"{path.name}.sha256"
        and integrity["covers"] == path.name,
        f"{label} receipt-integrity contract mismatch",
    )
    sidecar_path = path.with_name(f"{path.name}.sha256")
    sidecar_payload, sidecar_stat = _stable_regular_bytes(
        sidecar_path, f"{label} sidecar", immutable_mode=mode
    )
    try:
        fields = sidecar_payload.decode("ascii", errors="strict").strip().split()
    except UnicodeError as exc:
        raise VerificationError(f"{label} sidecar is not ASCII") from exc
    raw_sha256 = hashlib.sha256(payload).hexdigest()
    require(
        fields == [raw_sha256, path.name],
        f"{label} sidecar mismatch: does not authenticate the exact JSON bytes",
    )
    return document, {
        "path": str(path.resolve(strict=True)),
        "size_bytes": int(observed.st_size),
        "sha256": raw_sha256,
        "semantic_sha256": semantic_json_sha256(document),
        "st_dev": int(observed.st_dev),
        "st_ino": int(observed.st_ino),
        "sidecar_sha256": hashlib.sha256(sidecar_payload).hexdigest(),
        "sidecar_st_dev": int(sidecar_stat.st_dev),
        "sidecar_st_ino": int(sidecar_stat.st_ino),
    }


def _resource_gate_path(root: Path, gate_scope: str, batch_index: int | None = None) -> Path:
    if gate_scope == "session":
        name = "session.json"
    elif gate_scope == "batch" and isinstance(batch_index, int) and batch_index > 0:
        name = f"batch_{batch_index:03d}.json"
    else:
        raise VerificationError("resource gate scope/index is invalid")
    return root / ORCHESTRATION_DIRECTORY / "resource_gates" / name


def _verify_resource_gate(
    value: Any,
    *,
    root: Path,
    stage: str,
    gate_scope: str,
    expected_target: Mapping[str, Any],
    expected_producer_source: Mapping[str, Any],
    label: str,
    batch_index: int | None = None,
) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    identity = _exact_keys(value, RESOURCE_GATE_IDENTITY_KEYS, f"{label} identity")
    expected_path = _resource_gate_path(root, gate_scope, batch_index).resolve(strict=True)
    require(
        Path(str(identity["path"])).resolve(strict=True) == expected_path,
        f"{label} fixed path mismatch/path swap",
    )
    gate, authenticated = _authenticate_json_with_sidecar(
        expected_path, label, immutable=True
    )
    sidecar = expected_path.with_name(f"{expected_path.name}.sha256")
    expected_identity = {
        "path": str(expected_path),
        "sha256": authenticated["sha256"],
        "semantic_sha256": semantic_json_sha256(gate),
        "gate_id": gate.get("gate_id"),
        "sidecar_path": str(sidecar.resolve(strict=True)),
        "sidecar_sha256": sha256_path(sidecar),
    }
    require(identity == expected_identity, f"{label} identity/sidecar mismatch")
    _exact_keys(gate, RESOURCE_GATE_RECEIPT_KEYS, label)
    require(
        gate["schema_name"] == RESOURCE_GATE_SCHEMA_NAME
        and gate["schema_version"] == RESOURCE_GATE_SCHEMA_VERSION
        and gate["stage"] == stage
        and gate["gate_scope"] == gate_scope,
        f"{label} schema/stage/scope mismatch",
    )
    _hex_sha256(gate["gate_nonce"], f"{label}.gate_nonce")
    _hex_sha256(gate["gate_id"], f"{label}.gate_id")
    expected_gate_id = semantic_json_sha256({
        key: value for key, value in gate.items()
        if key not in {"gate_id", "receipt_integrity"}
    })
    require(gate["gate_id"] == expected_gate_id, f"{label}.gate_id mismatch")
    require(gate["target"] == dict(expected_target), f"{label} target/reuse mismatch")
    require(
        gate["producer_source"] == dict(expected_producer_source),
        f"{label} producer source mismatch",
    )
    process = _exact_keys(
        gate["process_snapshot"],
        {"process_count", "root_count", "representatives", "semantic_sha256"},
        f"{label}.process_snapshot",
    )
    process_core = {
        "process_count": process["process_count"],
        "root_count": process["root_count"],
        "representatives": process["representatives"],
    }
    require(
        process_core == {"process_count": 0, "root_count": 0, "representatives": []},
        f"{label} is not a zero-conflict observation",
    )
    require(
        process["semantic_sha256"] == semantic_json_sha256(process_core),
        f"{label} process snapshot SHA mismatch",
    )
    filesystem = _exact_keys(
        gate["filesystem_snapshot"],
        {
            "probe_path", "target_output_root", "st_dev", "f_bavail", "f_frsize",
            "available_bytes", "required_reserve_bytes", "disk_pass", "semantic_sha256",
        },
        f"{label}.filesystem_snapshot",
    )
    filesystem_core = {
        key: filesystem[key]
        for key in (
            "probe_path", "target_output_root", "st_dev", "f_bavail", "f_frsize",
            "available_bytes", "required_reserve_bytes", "disk_pass",
        )
    }
    root_stat = os.lstat(root.resolve(strict=True))
    _nonnegative_int(filesystem["f_bavail"], f"{label}.filesystem_snapshot.f_bavail")
    _nonnegative_int(filesystem["f_frsize"], f"{label}.filesystem_snapshot.f_frsize")
    _nonnegative_int(filesystem["available_bytes"], f"{label}.filesystem_snapshot.available_bytes")
    require(
        filesystem["probe_path"] == str(root.resolve(strict=True))
        and filesystem["target_output_root"] == str(root.resolve(strict=True))
        and filesystem["st_dev"] == int(root_stat.st_dev),
        f"{label} filesystem path/device identity mismatch",
    )
    require(
        filesystem["available_bytes"] == filesystem["f_bavail"] * filesystem["f_frsize"],
        f"{label} available-byte calculation mismatch",
    )
    require(
        filesystem["required_reserve_bytes"] == RESOURCE_GATE_REQUIRED_RESERVE_BYTES
        and filesystem["disk_pass"] is True
        and filesystem["available_bytes"] >= RESOURCE_GATE_REQUIRED_RESERVE_BYTES,
        f"{label} did not attest the exact 300 GiB reserve",
    )
    require(
        filesystem["semantic_sha256"] == semantic_json_sha256(filesystem_core),
        f"{label} filesystem snapshot SHA mismatch",
    )
    require(
        gate["checks"] == {"zero_conflict": True, "disk_reserve": True, "all_pass": True},
        f"{label} checks are not all PASS",
    )
    _nonempty_string(gate["observed_at_utc"], f"{label}.observed_at_utc")
    _nonnegative_int(gate["observed_monotonic_ns"], f"{label}.observed_monotonic_ns")
    _nonempty_string(gate["host_boot_id"], f"{label}.host_boot_id")
    return gate, expected_identity


def verify_post_input_identity_receipt(
    post_receipt_path: Path, release_binding: Mapping[str, Any]
) -> dict[str, Any]:
    """Authenticate POST and prove exact equality to the immutable frozen PRE.

    This is intentionally independent of the input-identity producer.  It does
    not re-read BAMs itself; it authenticates the producer receipt, independently
    verifies its exact authority/scope/role/count contract, and compares the
    complete persisted snapshot object to the PRE copy frozen in the release.
    """
    release_path = Path(str((release_binding.get("release_manifest") or {}).get("path", "")))
    release_document, release_identity = _authenticate_json_with_sidecar(
        release_path, "release manifest for POST gate", immutable=True
    )
    require(
        release_identity["sha256"]
        == (release_binding.get("release_manifest") or {}).get("sha256"),
        "POST gate release manifest differs from the authenticated release binding",
    )
    pre_contract = _exact_keys(
        release_document.get("pre_input_identity_receipt"),
        {
            "origin", "immutable_copy", "input_identity_snapshot_sha256",
            "receipt_semantic_sha256", "authority_mode", "validation_evidence_eligible",
            "artifact_roles",
        },
        "release.pre_input_identity_receipt",
    )
    pre_copy = pre_contract.get("immutable_copy")
    require(isinstance(pre_copy, dict), "release immutable PRE copy is malformed")
    pre_relative = Path(str(pre_copy.get("path", "")))
    require(
        bool(pre_relative.parts) and not pre_relative.is_absolute() and ".." not in pre_relative.parts,
        "release immutable PRE path is unsafe",
    )
    pre_path = release_path.parent / pre_relative
    pre_document, pre_identity = _authenticate_json_with_sidecar(
        pre_path, "immutable PRE receipt", immutable=True
    )
    require(pre_identity["sha256"] == pre_copy.get("sha256"), "immutable PRE raw SHA mismatch")
    require(
        pre_identity["semantic_sha256"] == pre_contract["receipt_semantic_sha256"],
        "immutable PRE semantic SHA mismatch",
    )
    require(
        pre_document.get("schema_name") == POST_INPUT_IDENTITY_SCHEMA_NAME
        and pre_document.get("schema_version") == POST_INPUT_IDENTITY_SCHEMA_VERSION
        and pre_document.get("mode") == "PRE"
        and pre_document.get("authority_mode") == "CANONICAL_V5_FROZEN"
        and pre_document.get("validation_evidence_eligible") is True
        and pre_document.get("all_pass") is True,
        "immutable PRE authority/schema/pass contract mismatch",
    )
    pre_snapshot = pre_document.get("input_identity_snapshot")
    require(isinstance(pre_snapshot, dict), "immutable PRE has no input snapshot")
    pre_snapshot_sha = semantic_json_sha256(pre_snapshot)
    require(
        pre_document.get("input_identity_snapshot_sha256") == pre_snapshot_sha
        == pre_contract["input_identity_snapshot_sha256"],
        "immutable PRE snapshot digest mismatch",
    )

    post_document, post_identity = _authenticate_json_with_sidecar(
        post_receipt_path, "POST input-identity receipt"
    )
    expected_post_keys = {
        "schema_name", "schema_version", "task_type", "mode", "authority_mode",
        "validation_evidence_eligible", "authority", "scope", "manifest",
        "canonical_schema", "verifier", "assurance", "input_identity_snapshot",
        "input_identity_snapshot_sha256", "artifact_audit", "compare_to", "checks",
        "all_pass", "receipt_integrity",
    }
    _exact_keys(post_document, expected_post_keys, "POST input-identity receipt")
    require(
        post_document["schema_name"] == POST_INPUT_IDENTITY_SCHEMA_NAME
        and post_document["schema_version"] == POST_INPUT_IDENTITY_SCHEMA_VERSION
        and post_document["task_type"] == "B_COMPREHENSIVE_VALIDATION"
        and post_document["mode"] == "POST_COMPARE"
        and post_document["authority_mode"] == "CANONICAL_V5_FROZEN"
        and post_document["validation_evidence_eligible"] is True
        and post_document["all_pass"] is True,
        "POST authority/schema/pass contract mismatch",
    )
    expected_scope = {
        "technical_datasets": 7,
        "biological_samples": 6,
        "chromosomes": "chr1-chr22",
        "tasks": 154,
        "datasets": list(DATASETS),
        "direct_input_artifacts": 42,
    }
    require(post_document["scope"] == expected_scope, "POST exact comprehensive scope mismatch")
    authority = _exact_keys(
        post_document["authority"],
        {
            "canonical_manifest_sha256", "canonical_schema_sha256",
            "canonical_schema_path", "selected_authority_is_canonical",
            "test_only_override",
        },
        "POST authority",
    )
    expected_schema = (release_binding.get("snapshot_sources") or {}).get(
        "canonical_manifest_schema"
    ) or {}
    require(
        authority["canonical_manifest_sha256"]
        == release_binding["canonical_input_manifest"]["sha256"]
        and authority["canonical_schema_sha256"] == expected_schema.get("sha256")
        and Path(str(authority["canonical_schema_path"])).resolve(strict=True)
        == Path(str(expected_schema.get("path", ""))).resolve(strict=True)
        and authority["selected_authority_is_canonical"] is True
        and authority["test_only_override"] is False,
        "POST canonical authority identity/selection mismatch",
    )
    checks = _exact_keys(
        post_document["checks"], POST_INPUT_IDENTITY_CHECK_NAMES, "POST checks"
    )
    require(all(checks[name] is True for name in POST_INPUT_IDENTITY_CHECK_NAMES), "POST named check failed")
    expected_verifier = (release_binding.get("snapshot_sources") or {}).get(
        "input_identity_verifier"
    )
    require(
        isinstance(expected_verifier, dict)
        and post_document["verifier"] == expected_verifier,
        "POST verifier path/SHA differs from frozen input-identity verifier",
    )
    canonical_input = release_binding.get("canonical_input_manifest") or {}
    require(
        (post_document.get("manifest") or {}).get("path") == canonical_input.get("path")
        and (post_document.get("manifest") or {}).get("sha256") == canonical_input.get("sha256")
        and (post_document.get("manifest") or {}).get("expected_sha256")
        == canonical_input.get("sha256"),
        "POST canonical manifest identity differs from release",
    )
    canonical_schema = post_document.get("canonical_schema") or {}
    require(
        canonical_schema.get("path") == expected_schema.get("path")
        and canonical_schema.get("sha256") == expected_schema.get("sha256")
        and canonical_schema.get("expected_sha256") == expected_schema.get("sha256"),
        "POST canonical schema identity differs from release",
    )
    assurance = _exact_keys(
        post_document["assurance"],
        {"bam", "other_direct_inputs", "pre_post", "temporal_immutability_proven"},
        "POST assurance",
    )
    require(
        assurance["temporal_immutability_proven"] is False,
        "POST overclaims temporal immutability",
    )
    compare_to = _exact_keys(
        post_document["compare_to"],
        {"path", "sha256", "pre_snapshot_sha256", "exact_snapshot_equal"},
        "POST compare_to",
    )
    require(
        Path(str(compare_to["path"])).resolve(strict=True) == pre_path.resolve(strict=True)
        and compare_to["sha256"] == pre_identity["sha256"]
        and compare_to["pre_snapshot_sha256"] == pre_snapshot_sha
        and compare_to["exact_snapshot_equal"] is True,
        "POST compare_to does not bind the immutable PRE",
    )
    post_snapshot = post_document["input_identity_snapshot"]
    require(
        post_snapshot == pre_snapshot
        and post_document["input_identity_snapshot_sha256"] == pre_snapshot_sha,
        "POST snapshot is not exactly equal to immutable PRE",
    )
    artifacts = post_snapshot.get("artifacts") if isinstance(post_snapshot, dict) else None
    require(isinstance(artifacts, list) and len(artifacts) == 42, "POST snapshot does not contain 42 roles")
    require(
        len({(row.get("dataset"), row.get("role")) for row in artifacts if isinstance(row, dict)})
        == 42,
        "POST snapshot artifact roles are not unique",
    )
    audit = post_document.get("artifact_audit") or {}
    require(
        audit.get("n_artifacts") == 42
        and audit.get("n_unique_resolved_files") == 42
        and audit.get("n_storage_identity_v1") == 7
        and audit.get("n_full_sha256") == 35
        and audit.get("n_sampled_bam_chunks") == 21
        and audit.get("n_mismatches") == 0,
        "POST artifact-audit count contract mismatch",
    )
    return {
        "post_receipt": {
            "path": post_identity["path"],
            "sha256": post_identity["sha256"],
            "semantic_sha256": post_identity["semantic_sha256"],
        },
        "immutable_pre_receipt": {
            "path": pre_identity["path"],
            "sha256": pre_identity["sha256"],
            "semantic_sha256": pre_identity["semantic_sha256"],
        },
        "input_identity_snapshot_sha256": pre_snapshot_sha,
        "n_artifact_roles": 42,
        "exact_snapshot_equal": True,
        "claim_boundary": (
            "Persistent PRE-to-POST input identity equality is verified; transient mutation "
            "followed by restoration is not proven absent."
        ),
    }


def _expected_orchestration_sources(
    stage: str, release_binding: Mapping[str, Any]
) -> dict[str, Any]:
    sources = release_binding["snapshot_sources"]
    if stage == "extraction":
        return {
            "runner": dict(sources["full_extraction_runner"]),
            "child_producer": dict(sources["extractor"]),
            "dependencies": {
                "lossless_read_contract": dict(sources["lossless_read_contract"]),
            },
        }
    require(stage == "ranking", f"unknown orchestration stage: {stage}")
    return {
        "runner": dict(sources["full_ranking_runner"]),
        "child_producer": dict(sources["ranker"]),
        "dependencies": {
            "hypercube_solver": dict(sources["hypercube_solver"]),
        },
    }


def _expected_scheduler_policy() -> dict[str, Any]:
    return {
        "initial_batch_tasks": 8,
        "subsequent_batch_tasks": 16,
        "selection_policy": "canonical DATASETS order then chr1-22; first N pending",
        "unattested_or_orphan_child_policy": "FAIL_CLOSED",
    }


def _canonical_task_id(dataset: str, chrom: str) -> str:
    return f"{dataset}:{chrom}"


def _canonical_task_ordinal(
    dataset: str, chrom: str, datasets: Sequence[str], chromosomes: Sequence[str]
) -> int:
    return datasets.index(dataset) * len(chromosomes) + chromosomes.index(chrom) + 1


def _identity_path_equals(value: Any, expected: Path, label: str) -> None:
    if isinstance(value, os.PathLike):
        value = os.fspath(value)
    _nonempty_string(value, label)
    try:
        actual = Path(str(value)).resolve(strict=True)
        wanted = expected.resolve(strict=True)
    except OSError as exc:
        raise VerificationError(f"{label} cannot be resolved: {exc}") from exc
    require(actual == wanted, f"{label} path mismatch")


def _require_real_directory_chain(root: Path, directory: Path, label: str) -> None:
    """Reject symlinked/escaping ancestors beneath a formal output root."""
    root_absolute = root.absolute()
    directory_absolute = directory.absolute()
    try:
        relative = directory_absolute.relative_to(root_absolute)
    except ValueError as exc:
        raise VerificationError(f"{label} escapes the output root") from exc
    current = root_absolute
    for part in relative.parts:
        current = current / part
        try:
            observed = os.lstat(current)
        except OSError as exc:
            raise VerificationError(f"{label} ancestor is unavailable: {current}: {exc}") from exc
        require(
            stat.S_ISDIR(observed.st_mode) and not stat.S_ISLNK(observed.st_mode),
            f"{label} ancestor is not a real directory: {current}",
        )
    require(
        directory_absolute.resolve(strict=True).is_relative_to(root_absolute.resolve(strict=True)),
        f"{label} resolves outside the output root",
    )


def _expected_child_parameters(
    stage: str, run_contract: Mapping[str, Any]
) -> Mapping[str, Any]:
    if stage == "extraction":
        return {
            "mapq_min": run_contract.get("mapq_min"),
            "baseq_min": run_contract.get("baseq_min"),
            "samtools_threads": run_contract.get("samtools_threads"),
            "bridge_thresholds": run_contract.get("bridge_thresholds"),
            "component_linkage_bases": run_contract.get("component_linkage_bases"),
        }
    parameters = run_contract.get("parameters")
    require(isinstance(parameters, dict), "ranking run contract has no child parameters")
    return parameters


def _verify_granted_command(
    command: Any,
    *,
    stage: str,
    dataset: str,
    chrom: str,
    output_dir: Path,
    run_contract: Mapping[str, Any],
    producer_sources: Mapping[str, Any],
    input_identity: Mapping[str, Any],
    extraction_root: Path | None,
    frozen_python_executable: Path,
) -> list[str]:
    require(
        isinstance(command, list)
        and command
        and all(isinstance(token, str) and token for token in command),
        f"{stage} grant command is not a nonempty argv list: {dataset}/{chrom}",
    )
    observed = list(command)
    require(len(observed) >= 2, f"{stage} grant command is truncated: {dataset}/{chrom}")
    _identity_path_equals(
        observed[0], frozen_python_executable, f"{stage} grant Python interpreter"
    )
    _identity_path_equals(
        observed[1], Path(str(producer_sources["child_producer"]["path"])),
        f"{stage} grant child producer",
    )
    if stage == "extraction":
        require(len(observed) == 18, f"extraction grant argv length mismatch: {dataset}/{chrom}")
        manifest_path = Path(str(input_identity["manifest_path"]))
        _identity_path_equals(observed[3], manifest_path, "extraction command manifest")
        expected = [
            observed[0], str(Path(str(producer_sources["child_producer"]["path"]))),
            "--manifest", str(manifest_path), "--sample", dataset, "--chrom", chrom,
            "--output-dir", str(output_dir.resolve()), "--mapq-min",
            str(run_contract["mapq_min"]), "--baseq-min", str(run_contract["baseq_min"]),
            "--bridge-thresholds",
            ",".join(str(value) for value in run_contract["bridge_thresholds"]),
            "--samtools-threads", str(run_contract["samtools_threads"]),
        ]
    else:
        require(extraction_root is not None, "ranking command verification needs extraction root")
        parameters = run_contract.get("parameters") or {}
        extraction_dir = extraction_root / "samples" / dataset / chrom
        expected = [
            observed[0], str(Path(str(producer_sources["child_producer"]["path"]))),
            "--input-dir", str(extraction_dir.resolve()),
            "--output-dir", str(output_dir.resolve()),
            "--structural-exact-pattern-minread-grid",
            ",".join(str(value) for value in parameters["structural_exact_pattern_minread_grid"]),
            "--primary-structural-exact-pattern-minread",
            str(parameters["primary_structural_exact_pattern_minread"]),
            "--exact-k-max", str(parameters["exact_k_max"]),
            "--max-vertex-sets", str(parameters["max_vertex_sets"]),
            "--solver-time-limit-seconds",
            str(parameters["solver_time_limit_seconds_per_milp"]),
            "--fixed-error-grid",
            ",".join(
                format(value, ".12g")
                for value in parameters[
                    "fixed_error_grid_conditional_binary_flip_probability"
                ]
            ),
            "--minimum-bq-error-rate", str(parameters["minimum_bq_error_rate"]),
            "--maximum-bq-error-rate", str(parameters["maximum_bq_error_rate"]),
            "--conditional-candidate-ranking-bootstrap-replicates",
            str(parameters["conditional_candidate_ranking_bootstrap_replicates"]),
            "--conditional-candidate-ranking-bootstrap-seed",
            str(parameters["conditional_candidate_ranking_bootstrap_base_seed"]),
            "--tie-tolerance", str(parameters["tie_tolerance_log_likelihood"]),
        ]
        for threshold in run_contract["thresholds"]:
            expected.extend(("--threshold", str(threshold)))
        for basis in run_contract["component_bases"]:
            expected.extend(("--component-basis", str(basis)))
        for family in run_contract["hp_families"]:
            expected.extend(("--hp-family", str(family)))
    require(observed == expected, f"{stage} grant command contract mismatch: {dataset}/{chrom}")
    return observed


def _authenticate_orchestration_session(
    root: Path,
    *,
    stage: str,
    full_receipt: Mapping[str, Any],
    release_binding: Mapping[str, Any],
    datasets: Sequence[str],
    chromosomes: Sequence[str],
    parent_extraction: Mapping[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    try:
        root_lstat = os.lstat(root)
    except OSError as exc:
        raise VerificationError(f"{stage} output root is unavailable: {exc}") from exc
    require(
        stat.S_ISDIR(root_lstat.st_mode) and not stat.S_ISLNK(root_lstat.st_mode),
        f"{stage} output root is not a real directory",
    )
    session_path = root / ORCHESTRATION_DIRECTORY / "session.json"
    _require_real_directory_chain(
        root, session_path.parent, f"{stage} orchestration session directory"
    )
    session, identity = _authenticate_json_with_sidecar(
        session_path, f"{stage} orchestration session", immutable=True
    )
    _exact_keys(session, ORCHESTRATION_SESSION_KEYS, f"{stage} orchestration session")
    require(
        session["schema_name"] == ORCHESTRATION_SCHEMA_NAMES["session"]
        and session["schema_version"] == ORCHESTRATION_SCHEMA_VERSION
        and session["stage"] == stage,
        f"{stage} orchestration session schema/stage mismatch",
    )
    _hex_sha256(session["session_nonce"], f"{stage} session_nonce")
    _nonempty_string(session["created_at_utc"], f"{stage} session created_at_utc")
    _nonnegative_int(session["created_monotonic_ns"], f"{stage} session monotonic time")
    _nonempty_string(session["host_boot_id"], f"{stage} session host_boot_id")
    require(
        session["release_manifest"] == release_binding["release_manifest"],
        f"{stage} session release manifest mismatch",
    )
    require(
        session["release_binding_semantic_sha256"]
        == semantic_json_sha256(release_binding),
        f"{stage} session release-binding semantic SHA mismatch",
    )
    run_contract = full_receipt.get("run_contract")
    require(isinstance(run_contract, dict), f"{stage} full receipt has no run contract")
    require(
        session["run_contract_semantic_sha256"] == semantic_json_sha256(run_contract),
        f"{stage} session run-contract semantic SHA mismatch",
    )
    require(
        session["scope"] == {
            "datasets": list(datasets), "chromosomes": list(chromosomes),
            "expected_tasks": len(datasets) * len(chromosomes),
        },
        f"{stage} orchestration session scope mismatch",
    )
    require(
        session["output_root"] == {
            "path": str(root.resolve(strict=True)),
            "st_dev": int(root_lstat.st_dev), "st_ino": int(root_lstat.st_ino),
        },
        f"{stage} orchestration output-root path/inode mismatch",
    )
    expected_sources = _expected_orchestration_sources(stage, release_binding)
    require(
        session["producer_sources"] == expected_sources,
        f"{stage} orchestration producer source mismatch",
    )
    require(
        session["scheduler_policy"] == _expected_scheduler_policy(),
        f"{stage} scheduler policy mismatch",
    )
    if stage == "extraction":
        require(session["parent_extraction"] is None, "extraction session has a parent extraction")
        session_target = {
            "output_root": dict(session["output_root"]),
            "release_manifest": dict(release_binding["release_manifest"]),
            "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
        }
        session_gate, _ = _verify_resource_gate(
            session["resource_gate"], root=root, stage="extraction",
            gate_scope="session", expected_target=session_target,
            expected_producer_source=session["producer_sources"]["runner"],
            label="extraction session resource_gate",
        )
        require(
            session["host_boot_id"] == session_gate["host_boot_id"]
            and session["created_monotonic_ns"] >= session_gate["observed_monotonic_ns"],
            "extraction session is not temporally bound to its resource gate",
        )
    else:
        require(parent_extraction is not None, "ranking session has no verified extraction parent")
        expected_parent = {
            "session_id": parent_extraction["session_id"],
            "terminal_receipt_path": parent_extraction["terminal_receipt"]["path"],
            "terminal_receipt_sha256": parent_extraction["terminal_receipt"]["sha256"],
        }
        require(
            session["parent_extraction"] == expected_parent,
            "ranking session is not bound to the verified extraction session/terminal",
        )
        session_target = {
            "output_root": dict(session["output_root"]),
            "release_manifest": dict(release_binding["release_manifest"]),
            "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
            "parent_extraction": expected_parent,
        }
        session_gate, _ = _verify_resource_gate(
            session["resource_gate"], root=root, stage="ranking",
            gate_scope="session", expected_target=session_target,
            expected_producer_source=session["producer_sources"]["runner"],
            label="ranking session resource_gate",
        )
        require(
            session["host_boot_id"] == session_gate["host_boot_id"]
            and session["created_monotonic_ns"] >= session_gate["observed_monotonic_ns"],
            "ranking session is not temporally bound to its resource gate",
        )
    expected_session_id = semantic_json_sha256({
        key: value for key, value in session.items()
        if key not in {"session_id", "receipt_integrity"}
    })
    require(session["session_id"] == expected_session_id, f"{stage} session_id mismatch")
    identity["session_id"] = expected_session_id
    return session, identity


def _verify_orchestration_grant_completion(
    root: Path,
    *,
    stage: str,
    dataset: str,
    chrom: str,
    task_ordinal: int,
    session: Mapping[str, Any],
    session_identity: Mapping[str, Any],
    batch: Mapping[str, Any],
    batch_identity: Mapping[str, Any],
    run_contract: Mapping[str, Any],
    release_binding: Mapping[str, Any],
    extraction_root: Path | None,
    extraction_children: Mapping[tuple[str, str], Mapping[str, Any]] | None,
    frozen_python_executable: Path,
) -> tuple[dict[str, Any], dict[str, str]]:
    task_id = _canonical_task_id(dataset, chrom)
    task_attestation_dir = root / ORCHESTRATION_DIRECTORY / "tasks" / dataset / chrom
    _require_real_directory_chain(root, task_attestation_dir, f"{stage} task attestation directory")
    grant_path = task_attestation_dir / "grant.json"
    completion_path = task_attestation_dir / "completion.json"
    grant, grant_identity = _authenticate_json_with_sidecar(
        grant_path, f"{stage} child grant {task_id}", immutable=True
    )
    _exact_keys(grant, ORCHESTRATION_GRANT_KEYS, f"{stage} child grant {task_id}")
    require(
        grant["schema_name"] == ORCHESTRATION_SCHEMA_NAMES["grant"]
        and grant["schema_version"] == ORCHESTRATION_SCHEMA_VERSION
        and grant["stage"] == stage
        and grant["session_id"] == session["session_id"]
        and grant["session_sha256"] == session_identity["sha256"]
        and grant["batch_id"] == batch["batch_id"]
        and grant["batch_start_sha256"] == batch_identity["sha256"]
        and grant["task_id"] == task_id
        and grant["dataset"] == dataset
        and grant["chrom"] == chrom
        and grant["task_ordinal"] == task_ordinal,
        f"{stage} grant session/batch/task binding mismatch: {task_id}",
    )
    _hex_sha256(grant["child_nonce"], f"{stage} grant child_nonce {task_id}")
    require(
        grant["producer_sources"] == session["producer_sources"],
        f"{stage} grant producer source mismatch: {task_id}",
    )
    expected_output_dir = root / "samples" / dataset / chrom
    _require_real_directory_chain(root, expected_output_dir, f"{stage} child output directory")
    require(
        grant["expected_output_dir"] == str(expected_output_dir.resolve(strict=True)),
        f"{stage} grant expected output directory mismatch: {task_id}",
    )
    expected_parameters = _expected_child_parameters(stage, run_contract)
    require(
        grant["parameters_semantic_sha256"] == semantic_json_sha256(expected_parameters),
        f"{stage} grant parameter semantic SHA mismatch: {task_id}",
    )
    input_identity = grant["input_identity"]
    if stage == "extraction":
        _exact_keys(
            input_identity, {"manifest_path", "manifest_sha256"},
            f"extraction grant input identity {task_id}",
        )
        canonical = release_binding["canonical_input_manifest"]
        require(
            input_identity == {
                "manifest_path": str(Path(canonical["path"]).resolve(strict=True)),
                "manifest_sha256": canonical["sha256"],
            },
            f"extraction grant canonical manifest identity mismatch: {task_id}",
        )
    else:
        require(
            extraction_root is not None and extraction_children is not None,
            "ranking grant verification lacks extraction evidence",
        )
        _exact_keys(
            input_identity,
            {
                "extraction_child_receipt_path", "extraction_child_receipt_sha256",
                "extraction_child_receipt_semantic_sha256",
                "extraction_outputs_semantic_sha256",
            },
            f"ranking grant input identity {task_id}",
        )
        extraction_path = extraction_root / "samples" / dataset / chrom / "receipt.json"
        extraction_document = extraction_children[(dataset, chrom)]
        require(
            input_identity == {
                "extraction_child_receipt_path": str(extraction_path.resolve(strict=True)),
                "extraction_child_receipt_sha256": sha256_path(extraction_path),
                "extraction_child_receipt_semantic_sha256": semantic_json_sha256(extraction_document),
                "extraction_outputs_semantic_sha256": semantic_json_sha256(
                    extraction_document.get("outputs") or {}
                ),
            },
            f"ranking grant extraction child identity mismatch: {task_id}",
        )
    command = _verify_granted_command(
        grant["command"], stage=stage, dataset=dataset, chrom=chrom,
        output_dir=expected_output_dir, run_contract=run_contract,
        producer_sources=session["producer_sources"], input_identity=input_identity,
        extraction_root=extraction_root,
        frozen_python_executable=frozen_python_executable,
    )
    require(
        grant["command_semantic_sha256"] == semantic_json_sha256(command),
        f"{stage} grant command semantic SHA mismatch: {task_id}",
    )
    _nonempty_string(grant["issued_at_utc"], f"{stage} grant issued_at_utc {task_id}")
    issued_ns = _nonnegative_int(
        grant["issued_monotonic_ns"], f"{stage} grant issued_monotonic_ns {task_id}"
    )
    _nonempty_string(grant["host_boot_id"], f"{stage} grant host_boot_id {task_id}")
    require(
        grant["host_boot_id"] == batch["host_boot_id"]
        and issued_ns >= batch["created_monotonic_ns"],
        f"{stage} grant is not temporally bound to its batch: {task_id}",
    )

    completion, completion_identity = _authenticate_json_with_sidecar(
        completion_path, f"{stage} child completion {task_id}", immutable=True
    )
    _exact_keys(
        completion, ORCHESTRATION_COMPLETION_KEYS,
        f"{stage} child completion {task_id}",
    )
    expected_grant_identity = {
        "path": grant_identity["path"], "sha256": grant_identity["sha256"],
        "semantic_sha256": grant_identity["semantic_sha256"],
    }
    require(
        completion["schema_name"] == ORCHESTRATION_SCHEMA_NAMES["completion"]
        and completion["schema_version"] == ORCHESTRATION_SCHEMA_VERSION
        and completion["stage"] == stage
        and completion["session_id"] == session["session_id"]
        and completion["session_sha256"] == session_identity["sha256"]
        and completion["batch_id"] == batch["batch_id"]
        and completion["grant_identity"] == expected_grant_identity
        and completion["task_id"] == task_id
        and completion["dataset"] == dataset
        and completion["chrom"] == chrom
        and completion["task_ordinal"] == task_ordinal
        and completion["command_semantic_sha256"] == grant["command_semantic_sha256"]
        and completion["status"] == "PASS"
        and completion["returncode"] == 0
        and completion["timed_out"] is False
        and completion["runner_source"] == session["producer_sources"]["runner"],
        f"{stage} completion session/grant/task/status binding mismatch: {task_id}",
    )
    _positive_int(completion["process_group_id"], f"{stage} completion process_group_id {task_id}")
    started_ns = _nonnegative_int(
        completion["started_monotonic_ns"],
        f"{stage} completion started_monotonic_ns {task_id}",
    )
    completed_ns = _nonnegative_int(
        completion["completed_monotonic_ns"],
        f"{stage} completion completed_monotonic_ns {task_id}",
    )
    _nonempty_string(completion["completed_at_utc"], f"{stage} completion completed_at_utc {task_id}")
    require(
        completion["host_boot_id"] == grant["host_boot_id"]
        and issued_ns <= started_ns <= completed_ns,
        f"{stage} completion timing/boot order mismatch: {task_id}",
    )
    child_path = root / "samples" / dataset / chrom / "receipt.json"
    child_document, child_identity = _authenticate_json_with_sidecar(
        child_path, f"{stage} attested child receipt {task_id}", immutable=True
    )
    declared_child = _exact_keys(
        completion["child_receipt_identity"],
        {"path", "size_bytes", "sha256", "semantic_sha256"},
        f"{stage} completion child receipt identity {task_id}",
    )
    require(
        declared_child == {
            "path": child_identity["path"], "size_bytes": child_identity["size_bytes"],
            "sha256": child_identity["sha256"],
            "semantic_sha256": child_identity["semantic_sha256"],
        },
        f"{stage} completion child receipt physical/semantic identity mismatch: {task_id}",
    )
    outputs = child_document.get("outputs") or {}
    require(isinstance(outputs, dict), f"{stage} child outputs are malformed: {task_id}")
    require(
        completion["child_outputs_semantic_sha256"] == semantic_json_sha256(outputs),
        f"{stage} completion child-output semantic SHA mismatch: {task_id}",
    )
    for output_identity in outputs.values():
        verify_identity(output_identity, child_path.parent)
    require(
        not (child_path.parent / "runner_task_status.json").exists(),
        f"{stage} timeout/task-status marker exists beside an attested child: {task_id}",
    )
    return child_document, {
        "path": completion_identity["path"], "sha256": completion_identity["sha256"],
    }


def verify_orchestration_stage(
    root: Path,
    full_receipt: Mapping[str, Any],
    *,
    stage: str,
    datasets: Sequence[str],
    chromosomes: Sequence[str],
    release_binding: Mapping[str, Any],
    extraction_root: Path | None = None,
    extraction_children: Mapping[tuple[str, str], Mapping[str, Any]] | None = None,
    parent_extraction: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    expected_pairs = expected_keys(datasets, chromosomes)
    require(len(expected_pairs) == 154, f"{stage} orchestration scope is not exactly 154 tasks")
    session, session_identity = _authenticate_orchestration_session(
        root, stage=stage, full_receipt=full_receipt, release_binding=release_binding,
        datasets=datasets, chromosomes=chromosomes, parent_extraction=parent_extraction,
    )
    run_contract = full_receipt["run_contract"]
    release_document, _ = _authenticate_json_with_sidecar(
        Path(str(release_binding["release_manifest"]["path"])),
        f"{stage} frozen runtime release manifest", immutable=True,
    )
    frozen_python_executable = Path(
        str(((release_document.get("runtime") or {}).get("python") or {}).get("executable", ""))
    )
    _nonempty_string(str(frozen_python_executable), f"{stage} frozen Python executable")
    expected_task_ids = [_canonical_task_id(*pair) for pair in expected_pairs]
    previous_chain_head: dict[str, str] | None = None
    completion_identities: dict[str, dict[str, str]] = {}
    authenticated_children: dict[tuple[str, str], dict[str, Any]] = {}
    expected_checkpoint_paths: set[Path] = set()
    expected_batch_paths: set[Path] = set()
    expected_grant_paths: set[Path] = set()
    expected_completion_paths: set[Path] = set()
    previous_count = 0
    terminal_identity: dict[str, str] | None = None
    for batch_index, cumulative_count in enumerate(ORCHESTRATION_COMPLETED_COUNTS, start=1):
        effective_count = cumulative_count - previous_count
        batch_path = (
            root / ORCHESTRATION_DIRECTORY / "batches" /
            f"batch_{batch_index:03d}_start.json"
        )
        _require_real_directory_chain(
            root, batch_path.parent, f"{stage} batch attestation directory"
        )
        expected_batch_paths.add(batch_path.resolve(strict=True))
        batch, batch_identity = _authenticate_json_with_sidecar(
            batch_path, f"{stage} batch {batch_index} start", immutable=True
        )
        _exact_keys(batch, ORCHESTRATION_BATCH_KEYS, f"{stage} batch {batch_index} start")
        selected_ids = expected_task_ids[previous_count:cumulative_count]
        require(
            batch["schema_name"] == ORCHESTRATION_SCHEMA_NAMES["batch"]
            and batch["schema_version"] == ORCHESTRATION_SCHEMA_VERSION
            and batch["stage"] == stage
            and batch["session_id"] == session["session_id"]
            and batch["session_sha256"] == session_identity["sha256"]
            and batch["batch_index"] == batch_index
            and batch["previous_chain_head"] == previous_chain_head
            and batch["before_count"] == previous_count
            and batch["max_new_tasks"] == (8 if batch_index == 1 else 16)
            and batch["effective_count"] == effective_count
            and batch["selected_task_ids"] == selected_ids
            and batch["run_contract_semantic_sha256"] == semantic_json_sha256(run_contract)
            and batch["runner_source"] == session["producer_sources"]["runner"],
            f"{stage} batch {batch_index} chain/count/order/source mismatch",
        )
        _hex_sha256(batch["batch_nonce"], f"{stage} batch {batch_index} nonce")
        expected_batch_id = semantic_json_sha256({
            key: value for key, value in batch.items()
            if key not in {"batch_id", "receipt_integrity"}
        })
        require(batch["batch_id"] == expected_batch_id, f"{stage} batch {batch_index} id mismatch")
        _nonempty_string(batch["created_at_utc"], f"{stage} batch {batch_index} created_at_utc")
        _nonnegative_int(batch["created_monotonic_ns"], f"{stage} batch {batch_index} monotonic time")
        _nonempty_string(batch["host_boot_id"], f"{stage} batch {batch_index} host_boot_id")
        gate_target = {
            "output_root": dict(session["output_root"]),
            "session_id": session["session_id"],
            "session_sha256": session_identity["sha256"],
            "batch_index": batch_index,
            "before_count": previous_count,
            "max_new_tasks": 8 if batch_index == 1 else 16,
            "effective_count": effective_count,
            "selected_task_ids": selected_ids,
            "previous_chain_head": previous_chain_head,
        }
        batch_gate, _ = _verify_resource_gate(
            batch["resource_gate"], root=root, stage=stage, gate_scope="batch",
            batch_index=batch_index, expected_target=gate_target,
            expected_producer_source=session["producer_sources"]["runner"],
            label=f"{stage} batch {batch_index} resource_gate",
        )
        require(
            batch["host_boot_id"] == batch_gate["host_boot_id"]
            and batch["created_monotonic_ns"] >= batch_gate["observed_monotonic_ns"],
            f"{stage} batch {batch_index} is not temporally bound to its resource gate",
        )

        current_completion_rows: list[dict[str, str]] = []
        for offset, pair in enumerate(expected_pairs[previous_count:cumulative_count], start=previous_count + 1):
            dataset, chrom = pair
            child, completion_identity = _verify_orchestration_grant_completion(
                root, stage=stage, dataset=dataset, chrom=chrom, task_ordinal=offset,
                session=session, session_identity=session_identity, batch=batch,
                batch_identity=batch_identity, run_contract=run_contract,
                release_binding=release_binding, extraction_root=extraction_root,
                extraction_children=extraction_children,
                frozen_python_executable=frozen_python_executable,
            )
            task_id = _canonical_task_id(dataset, chrom)
            require(task_id not in completion_identities, f"duplicate completion: {stage}/{task_id}")
            completion_identities[task_id] = completion_identity
            authenticated_children[(dataset, chrom)] = child
            current_completion_rows.append({"task_id": task_id, **completion_identity})
            task_dir = root / ORCHESTRATION_DIRECTORY / "tasks" / dataset / chrom
            expected_grant_paths.add((task_dir / "grant.json").resolve(strict=True))
            expected_completion_paths.add((task_dir / "completion.json").resolve(strict=True))

        if cumulative_count == 154:
            chain_receipt_path = root / (
                "full_extraction_receipt.json" if stage == "extraction"
                else "full_ranking_receipt.json"
            )
        else:
            chain_receipt_path = (
                root / "checkpoints" / f"checkpoint_{cumulative_count:03d}_of_154.json"
            )
            _require_real_directory_chain(
                root, chain_receipt_path.parent, f"{stage} checkpoint directory"
            )
            expected_checkpoint_paths.add(chain_receipt_path.resolve(strict=True))
        chain_receipt, chain_identity = _authenticate_json_with_sidecar(
            chain_receipt_path, f"{stage} orchestration chain receipt {cumulative_count}",
            immutable=True,
        )
        require(
            chain_receipt.get("schema_name")
            == ("intersubmod.m2_full_extraction_receipt" if stage == "extraction"
                else "intersubmod.m2_full_ranking_receipt")
            and chain_receipt.get("schema_version") == ("1.2.0" if stage == "extraction" else "2.0.0"),
            f"{stage} chain receipt schema mismatch at {cumulative_count}",
        )
        orchestration = _exact_keys(
            chain_receipt.get("orchestration"), ORCHESTRATION_RECEIPT_KEYS,
            f"{stage} chain receipt orchestration {cumulative_count}",
        )
        expected_session_identity = {
            "path": session_identity["path"], "sha256": session_identity["sha256"],
            "session_id": session["session_id"],
        }
        expected_batch_identity = {
            "path": batch_identity["path"], "sha256": batch_identity["sha256"],
            "batch_id": batch["batch_id"], "batch_index": batch_index,
        }
        require(
            orchestration["session_identity"] == expected_session_identity
            and orchestration["batch_start_identity"] == expected_batch_identity
            and orchestration["previous_chain_head"] == previous_chain_head
            and orchestration["batch_completion_attestations"] == current_completion_rows
            and orchestration["cumulative_attested_tasks"] == cumulative_count,
            f"{stage} chain receipt orchestration mismatch at {cumulative_count}",
        )
        results = chain_receipt.get("results")
        require(
            isinstance(results, list) and len(results) == cumulative_count,
            f"{stage} chain receipt result count mismatch at {cumulative_count}",
        )
        result_index: dict[str, Mapping[str, Any]] = {}
        for result in results:
            require(isinstance(result, dict), f"{stage} chain result is malformed")
            task_id = _canonical_task_id(str(result.get("dataset")), str(result.get("chrom")))
            require(task_id not in result_index, f"{stage} duplicate chain result {task_id}")
            result_index[task_id] = result
        require(
            set(result_index) == set(expected_task_ids[:cumulative_count]),
            f"{stage} chain result task set mismatch at {cumulative_count}",
        )
        for task_id in expected_task_ids[:cumulative_count]:
            result = result_index[task_id]
            require(
                result.get("status") in {"PASS", "REUSED_PASS"}
                and result.get("orchestration_completion") == completion_identities[task_id],
                f"{stage} chain result lacks exact completion identity: {task_id}",
            )
            dataset, chrom = task_id.split(":", 1)
            if stage == "extraction":
                require(
                    result.get("receipt") == authenticated_children[(dataset, chrom)],
                    f"extraction chain embeds a different child receipt: {task_id}",
                )
            else:
                rank_identity = result.get("rank_receipt") or {}
                rank_path = root / "samples" / dataset / chrom / "receipt.json"
                require(
                    Path(str(rank_identity.get("path", ""))).resolve(strict=True)
                    == rank_path.resolve(strict=True)
                    and rank_identity.get("sha256") == sha256_path(rank_path),
                    f"ranking chain child receipt identity mismatch: {task_id}",
                )
        previous_chain_head = {
            "path": chain_identity["path"], "sha256": chain_identity["sha256"],
        }
        previous_count = cumulative_count
        if cumulative_count == 154:
            terminal_identity = dict(previous_chain_head)

    require(len(completion_identities) == 154, f"{stage} does not have exactly 154 completions")
    actual_batches = {
        path.resolve(strict=True)
        for path in (root / ORCHESTRATION_DIRECTORY / "batches").glob("batch_*_start.json")
    }
    require(actual_batches == expected_batch_paths, f"{stage} has an open/orphan/extra batch")
    expected_gate_paths = {_resource_gate_path(root, "session").resolve(strict=True)}
    expected_gate_paths.update(
        _resource_gate_path(root, "batch", index).resolve(strict=True)
        for index in range(1, len(ORCHESTRATION_COMPLETED_COUNTS) + 1)
    )
    actual_gate_paths = {
        path.resolve(strict=True)
        for path in (root / ORCHESTRATION_DIRECTORY / "resource_gates").glob("*.json")
    }
    require(
        actual_gate_paths == expected_gate_paths,
        f"{stage} has missing/orphan/extra/reused resource gate receipts",
    )
    actual_grants = {
        path.resolve(strict=True)
        for path in (root / ORCHESTRATION_DIRECTORY / "tasks").glob("*/*/grant.json")
    }
    actual_completions = {
        path.resolve(strict=True)
        for path in (root / ORCHESTRATION_DIRECTORY / "tasks").glob("*/*/completion.json")
    }
    require(actual_grants == expected_grant_paths, f"{stage} has orphan/missing/extra child grants")
    require(
        actual_completions == expected_completion_paths,
        f"{stage} has orphan/missing/extra child completions",
    )
    checkpoint_dir = root / "checkpoints"
    actual_checkpoints = {
        path.resolve(strict=True) for path in checkpoint_dir.glob("checkpoint_*_of_154.json")
    }
    require(
        actual_checkpoints == expected_checkpoint_paths,
        f"{stage} checkpoint count chain has a gap/reorder/extra receipt",
    )
    ending_root_stat = os.lstat(root)
    require(
        stat.S_ISDIR(ending_root_stat.st_mode)
        and not stat.S_ISLNK(ending_root_stat.st_mode)
        and ending_root_stat.st_dev == session["output_root"]["st_dev"]
        and ending_root_stat.st_ino == session["output_root"]["st_ino"],
        f"{stage} output root changed during orchestration verification",
    )
    require(terminal_identity is not None, f"{stage} terminal identity was not reconstructed")
    return {
        "stage": stage,
        "session_id": session["session_id"],
        "session_receipt": {
            "path": session_identity["path"], "sha256": session_identity["sha256"],
        },
        "terminal_receipt": terminal_identity,
        "legal_cumulative_counts": list(ORCHESTRATION_COMPLETED_COUNTS),
        "n_batches": len(ORCHESTRATION_COMPLETED_COUNTS),
        "n_authenticated_resource_gates": len(expected_gate_paths),
        "n_attested_children": len(completion_identities),
        "all_child_receipts_and_outputs_rehashed": True,
        "all_grants_and_completions_session_bound": True,
        "all_checkpoints_and_terminal_chain_bound": True,
        "no_open_orphan_or_preseeded_child_accepted": True,
        "claim_boundary": (
            "Authenticated parent grants/completions and immutable chain receipts prove the "
            "persisted orchestration protocol under the non-hostile same-UID assumption; they "
            "do not constitute an operating-system process-ancestry proof against a malicious "
            "same-UID actor."
        ),
    }


def _reauthenticate_orchestration_publication_boundary(
    audit: Mapping[str, Any], label: str
) -> None:
    for key in ("session_receipt", "terminal_receipt"):
        expected = audit.get(key) or {}
        path = Path(str(expected.get("path", "")))
        _, observed = _authenticate_json_with_sidecar(
            path, f"{label} publication-boundary {key}", immutable=True
        )
        require(
            observed["path"] == expected.get("path")
            and observed["sha256"] == expected.get("sha256"),
            f"{label} {key} drifted before final publication",
        )


def write_immutable_receipt_exclusive(path: Path, receipt: Mapping[str, Any]) -> None:
    """Publish JSON+sidecar without overwriting any prior PASS or FAIL evidence."""
    sidecar = path.with_name(f"{path.name}.sha256")
    require(
        not os.path.lexists(path) and not os.path.lexists(sidecar),
        "final output or sidecar already exists",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        receipt, ensure_ascii=False, allow_nan=False, indent=2
    ).encode("utf-8") + b"\n"
    with path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    digest = hashlib.sha256(payload).hexdigest()
    sidecar_payload = f"{digest}  {path.name}\n".encode("ascii")
    with sidecar.open("xb") as handle:
        handle.write(sidecar_payload)
        handle.flush()
        os.fsync(handle.fileno())
    path.chmod(0o444)
    sidecar.chmod(0o444)


def _deep_verify_frozen_release(
    release_path: Path,
    release_sha256: str,
    identities: Mapping[str, Mapping[str, str]],
    *,
    force_fresh: bool = False,
) -> dict[str, Any]:
    """Execute only the freezer adjacent to this authenticated verifier."""
    anchor_relative = Path(RELEASE_VERIFIER_RELATIVE)
    anchor = VERIFIER.resolve(strict=True)
    if tuple(anchor.parts[-len(anchor_relative.parts):]) != anchor_relative.parts:
        raise VerificationError("release verifier is not at its hardcoded snapshot path")
    snapshot_root = anchor
    for _ in anchor_relative.parts:
        snapshot_root = snapshot_root.parent
    contract_root = release_path.parent.resolve(strict=True)
    require(
        snapshot_root.parent == contract_root,
        "release manifest is not adjacent to this verifier's code_snapshot",
    )
    freezer = snapshot_root / RELEASE_FREEZER_RELATIVE
    expected = identities.get("release_contract_freezer") or {}
    try:
        observed = os.lstat(freezer)
        resolved_freezer = freezer.resolve(strict=True)
    except OSError as exc:
        raise VerificationError(f"anchored frozen freezer unavailable: {exc}") from exc
    require(
        str(resolved_freezer) == expected.get("path")
        and stat.S_ISREG(observed.st_mode)
        and not stat.S_ISLNK(observed.st_mode)
        and stat.S_IMODE(observed.st_mode) == 0o444
        and observed.st_nlink == 1
        and sha256_path(freezer) == expected.get("sha256"),
        "anchored frozen freezer identity mismatch",
    )
    freezer_sha256 = sha256_path(freezer)
    cache_key = (release_sha256, freezer_sha256)
    if not force_fresh and cache_key in _DEEP_RELEASE_VERIFY_CACHE:
        return _DEEP_RELEASE_VERIFY_CACHE[cache_key]
    module_name = f"_m2_frozen_freezer_{freezer_sha256[:16]}"
    spec = importlib.util.spec_from_file_location(module_name, freezer)
    require(spec is not None and spec.loader is not None, "cannot load frozen release verifier")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
        verification = module.verify_release_contract(contract_root)
    except Exception as exc:
        raise VerificationError(f"frozen freezer deep verification failed: {exc}") from exc
    finally:
        sys.modules.pop(module_name, None)
    identity = verification.get("run_manifest") or {}
    require(
        verification.get("all_pass") is True
        and verification.get("validation_evidence_eligible") is True
        and identity.get("sha256") == release_sha256,
        "frozen freezer did not certify this release manifest",
    )
    summary = {
        "mode": "FROZEN_FREEZER_VERIFY_RELEASE_CONTRACT",
        "freezer_path": str(freezer.resolve()), "freezer_sha256": freezer_sha256,
        "release_manifest_sha256": release_sha256,
        "verified_snapshot_files": (verification.get("snapshot") or {}).get("n_files"),
        "all_pass": True,
    }
    _DEEP_RELEASE_VERIFY_CACHE[cache_key] = summary
    return summary


def _release_file_record(
    path_value: str, observed: os.stat_result, sha256: str
) -> dict[str, Any]:
    return {
        "path": path_value,
        "st_dev": int(observed.st_dev), "st_ino": int(observed.st_ino),
        "st_nlink": int(observed.st_nlink), "size_bytes": int(observed.st_size),
        "mtime_ns": int(observed.st_mtime_ns), "ctime_ns": int(observed.st_ctime_ns),
        "mode_octal": format(stat.S_IMODE(observed.st_mode), "04o"), "sha256": sha256,
    }


def _release_validate_copy(
    record: Any, physical: Path, declared_path: str, label: str
) -> dict[str, Any]:
    row = _exact_keys(record, RELEASE_FILE_IDENTITY_KEYS, label)
    payload, observed = _stable_regular_bytes(physical, label, immutable_mode=0o444)
    expected = _release_file_record(
        declared_path, observed, hashlib.sha256(payload).hexdigest()
    )
    require(row == expected, f"{label} identity/stat/SHA mismatch")
    return dict(row)


def _release_validate_source(
    record: Any, repo_root: Path, repo_relative: str, label: str
) -> tuple[dict[str, Any], tuple[int, int]]:
    row = _exact_keys(record, RELEASE_FILE_IDENTITY_KEYS, label)
    physical = repo_root / repo_relative
    _require_real_directory_chain(repo_root, physical.parent, f"{label} parent")
    payload, observed = _stable_regular_bytes(physical, label)
    expected = _release_file_record(
        str(physical.resolve(strict=True)), observed, hashlib.sha256(payload).hexdigest()
    )
    require(row == expected, f"{label} origin identity/stat/SHA mismatch")
    return dict(row), (int(observed.st_dev), int(observed.st_ino))


def _release_validate_runtime(document: Mapping[str, Any]) -> None:
    runtime = _exact_keys(
        document.get("runtime"), {"python", "packages", "samtools", "platform"},
        "release runtime",
    )
    python = _exact_keys(
        runtime["python"], {"executable", "version", "implementation"},
        "release runtime.python",
    )
    packages = _exact_keys(
        runtime["packages"], {"numpy", "scipy", "pysam"},
        "release runtime.packages",
    )
    samtools = _exact_keys(
        runtime["samtools"], {"path", "version_line", "htslib_version_line"},
        "release runtime.samtools",
    )
    require(
        all(
            isinstance(value, str) and bool(value)
            for value in [
                *python.values(), *packages.values(), *samtools.values(), runtime["platform"],
            ]
        ),
        "release runtime contains an empty/non-string field",
    )
    require(
        document.get("runtime_semantic_sha256") == semantic_json_sha256(runtime),
        "release runtime semantic digest mismatch",
    )


def load_release_contract_binding(
    manifest_path: Path,
    *,
    required_sources: Mapping[str, Path],
    _skip_deep_verification_for_test: bool = False,
    _force_deep_reverification: bool = False,
) -> dict[str, Any]:
    """Statically authenticate the exact release before executing its freezer.

    In particular, every origin source is independently re-read and matched to
    its frozen copy before the release-provided freezer is imported.  This
    ordering prevents a self-signed manifest from turning the deep verifier
    hook into a code-execution primitive.
    """
    release_path = manifest_path.absolute()
    document, release_identity = _authenticate_json_with_sidecar(
        release_path, "release manifest", immutable=True
    )
    _exact_keys(document, RELEASE_TOP_KEYS, "release manifest")
    expected_scope = {
        "technical_datasets": 7, "biological_samples": 6, "chromosomes": 22,
        "chromosome_names": list(AUTOSOMES), "tasks": len(DATASETS) * len(AUTOSOMES),
        "datasets": list(DATASETS),
    }
    require(
        document["schema_name"] == RELEASE_SCHEMA_NAME
        and document["schema_version"] == RELEASE_SCHEMA_VERSION
        and document["task_type"] == "B_COMPREHENSIVE_VALIDATION"
        and document["authority_mode"] == "CANONICAL_V5_FROZEN"
        and document["validation_evidence_eligible"] is True
        and document["all_pass"] is True
        and isinstance(document["created_at_utc"], str)
        and document["created_at_utc"].endswith("Z")
        and document["scope"] == expected_scope,
        "release canonical authority/task/scope metadata mismatch",
    )
    checks = _exact_keys(document["checks"], RELEASE_CHECK_NAMES, "release checks")
    require(all(value is True for value in checks.values()), "release named check failed")
    require(
        document["parameters"] == FROZEN_RELEASE_PARAMETERS,
        "release frozen parameters differ from canonical M2 parameters",
    )
    require(document["assurance"] == RELEASE_ASSURANCE, "release assurance mismatch")
    _release_validate_runtime(document)

    snapshot = _exact_keys(
        document["source_snapshot"],
        {
            "repo_root", "snapshot_root", "n_files", "entries", "entrypoints",
            "exact_allowlist_semantic_sha256",
        },
        "release source_snapshot",
    )
    repo_root = Path(str(snapshot["repo_root"]))
    require(repo_root.is_absolute(), "release source repo_root is not absolute")
    try:
        repo_stat = os.lstat(repo_root)
    except OSError as exc:
        raise VerificationError(f"release source repo_root is unavailable: {exc}") from exc
    require(
        stat.S_ISDIR(repo_stat.st_mode) and not stat.S_ISLNK(repo_stat.st_mode),
        "release source repo_root is not a real directory",
    )
    require(snapshot["snapshot_root"] == "code_snapshot", "release snapshot root mismatch")
    snapshot_root = release_path.parent / "code_snapshot"
    _require_real_directory_chain(
        release_path.parent, snapshot_root, "release code_snapshot root"
    )
    entries = snapshot["entries"]
    require(
        isinstance(entries, list) and len(entries) == 11 and snapshot["n_files"] == 11,
        "release snapshot cardinality is not exactly eleven",
    )
    expected_code_entrypoints = {
        role: (Path("code_snapshot") / relative).as_posix()
        for role, relative in RELEASE_SOURCE_PATHS.items()
        if role != "canonical_manifest_schema"
    }
    require(
        snapshot["entrypoints"] == expected_code_entrypoints,
        "release snapshot exact entrypoint map mismatch",
    )
    expected_allowlist = [
        {"role": role, "repo_relative_path": relative}
        for role, relative in RELEASE_SOURCE_PATHS.items()
    ]
    require(
        snapshot["exact_allowlist_semantic_sha256"]
        == semantic_json_sha256(expected_allowlist),
        "release snapshot allowlist semantic digest mismatch",
    )
    identities: dict[str, dict[str, str]] = {}
    producer_records: dict[str, dict[str, Any]] = {}
    source_inodes: set[tuple[int, int]] = set()
    copy_inodes: set[tuple[int, int]] = set()
    for entry in entries:
        _exact_keys(
            entry, {"role", "repo_relative_path", "source", "snapshot"},
            "release snapshot entry",
        )
        role = entry["role"]
        require(
            role in RELEASE_SOURCE_PATHS and role not in identities,
            f"release duplicate/unknown snapshot role: {role!r}",
        )
        repo_relative = RELEASE_SOURCE_PATHS[role]
        require(
            entry["repo_relative_path"] == repo_relative,
            f"release snapshot repo path mismatch: {role}",
        )
        source, source_inode = _release_validate_source(
            entry["source"], repo_root, repo_relative, f"release source/{role}"
        )
        relative = Path("code_snapshot") / repo_relative
        physical = release_path.parent / relative
        _require_real_directory_chain(
            release_path.parent, physical.parent, f"release snapshot/{role} parent"
        )
        copied = _release_validate_copy(
            entry["snapshot"], physical, relative.as_posix(), f"release snapshot/{role}"
        )
        require(
            source["sha256"] == copied["sha256"]
            and source["st_nlink"] == copied["st_nlink"] == 1,
            f"release source/copy bytes or single-link contract mismatch: {role}",
        )
        copy_inode = (int(copied["st_dev"]), int(copied["st_ino"]))
        require(
            source_inode not in source_inodes and copy_inode not in copy_inodes
            and source_inode != copy_inode,
            f"release source/copy inode alias detected: {role}",
        )
        source_inodes.add(source_inode)
        copy_inodes.add(copy_inode)
        identities[role] = {
            "path": str(physical.resolve(strict=True)), "sha256": copied["sha256"],
        }
        producer_records[role] = dict(entry)
    require(
        set(identities) == set(RELEASE_SOURCE_PATHS)
        and len(source_inodes) == len(copy_inodes) == 11,
        "release snapshot role/inode set is not the exact canonical eleven",
    )
    for role, actual_path in required_sources.items():
        require(
            role in identities and role in expected_code_entrypoints,
            f"required frozen source role is missing: {role}",
        )
        _identity_path_equals(actual_path, Path(identities[role]["path"]), f"current source/{role}")
        require(
            sha256_path(Path(actual_path)) == identities[role]["sha256"],
            f"current source SHA differs from frozen role: {role}",
        )

    input_root = release_path.parent / "input_contract"
    _require_real_directory_chain(release_path.parent, input_root, "release input_contract")
    canonical = _exact_keys(
        document["canonical_manifest"], {"expected_sha256", "origin", "immutable_copy"},
        "release canonical_manifest",
    )
    require(
        canonical["expected_sha256"] == CANONICAL_MANIFEST_SHA256,
        "release canonical manifest authority SHA mismatch",
    )
    _exact_keys(canonical["origin"], RELEASE_FILE_IDENTITY_KEYS, "canonical manifest origin")
    canonical_relative = Path("input_contract/canonical_manifest.json")
    canonical_path = release_path.parent / canonical_relative
    canonical_copy = _release_validate_copy(
        canonical["immutable_copy"], canonical_path, canonical_relative.as_posix(),
        "release canonical manifest copy",
    )
    require(
        canonical_copy["sha256"] == CANONICAL_MANIFEST_SHA256,
        "release canonical manifest copy SHA mismatch",
    )
    copy_identity = {
        "path": str(canonical_path.resolve(strict=True)), "sha256": canonical_copy["sha256"],
    }
    require(
        document["canonical_schema"] == {
            "role": "canonical_manifest_schema",
            "repo_relative_path": CANONICAL_SCHEMA_RELATIVE,
            "sha256": CANONICAL_SCHEMA_SHA256,
        },
        "release canonical schema authority mismatch",
    )

    pre = _exact_keys(
        document["pre_input_identity_receipt"],
        {
            "origin", "immutable_copy", "input_identity_snapshot_sha256",
            "receipt_semantic_sha256", "authority_mode",
            "validation_evidence_eligible", "artifact_roles",
        },
        "release PRE contract",
    )
    require(
        pre["authority_mode"] == "CANONICAL_V5_FROZEN"
        and pre["validation_evidence_eligible"] is True
        and pre["artifact_roles"] == 42,
        "release PRE authority/scope mismatch",
    )
    _exact_keys(
        pre["origin"], RELEASE_FILE_IDENTITY_KEYS | {"sidecar"}, "release PRE origin"
    )
    pre_copy = _exact_keys(
        pre["immutable_copy"], RELEASE_FILE_IDENTITY_KEYS | {"sidecar"},
        "release PRE immutable copy",
    )
    pre_relative = Path(str(pre_copy["path"]))
    require(
        not pre_relative.is_absolute() and pre_relative.parts
        and pre_relative.parts[0] == "input_contract" and ".." not in pre_relative.parts,
        "release PRE copy path escapes input_contract",
    )
    pre_path = release_path.parent / pre_relative
    require(
        pre_path.resolve(strict=True).is_relative_to(input_root.resolve(strict=True)),
        "release PRE copy resolves outside input_contract",
    )
    _release_validate_copy(
        {key: pre_copy[key] for key in RELEASE_FILE_IDENTITY_KEYS},
        pre_path, pre_relative.as_posix(), "release PRE copy",
    )
    pre_side = _exact_keys(
        pre_copy["sidecar"], RELEASE_FILE_IDENTITY_KEYS, "release PRE copy sidecar"
    )
    pre_side_relative = Path(str(pre_side["path"]))
    require(
        not pre_side_relative.is_absolute() and pre_side_relative.parts
        and pre_side_relative.parts[0] == "input_contract"
        and ".." not in pre_side_relative.parts,
        "release PRE sidecar path escapes input_contract",
    )
    pre_side_path = release_path.parent / pre_side_relative
    require(
        pre_side_path.resolve(strict=True).is_relative_to(input_root.resolve(strict=True)),
        "release PRE sidecar resolves outside input_contract",
    )
    _release_validate_copy(
        pre_side, pre_side_path, pre_side_relative.as_posix(), "release PRE copy sidecar"
    )
    pre_document, pre_identity = _authenticate_json_with_sidecar(
        pre_path, "release immutable PRE", immutable=True
    )
    pre_snapshot = pre_document.get("input_identity_snapshot")
    require(
        pre_document.get("schema_name") == POST_INPUT_IDENTITY_SCHEMA_NAME
        and pre_document.get("schema_version") == POST_INPUT_IDENTITY_SCHEMA_VERSION
        and pre_document.get("mode") == "PRE"
        and pre_document.get("authority_mode") == "CANONICAL_V5_FROZEN"
        and pre_document.get("validation_evidence_eligible") is True
        and pre_document.get("all_pass") is True
        and isinstance(pre_snapshot, dict)
        and pre_document.get("input_identity_snapshot_sha256")
        == semantic_json_sha256(pre_snapshot)
        == pre["input_identity_snapshot_sha256"]
        and pre_identity["semantic_sha256"] == pre["receipt_semantic_sha256"],
        "release immutable PRE semantic/authority/snapshot mismatch",
    )

    fresh = _exact_keys(
        document["fresh_input_identity_verification"],
        {
            "execution_mode", "verifier_path", "verifier_sha256", "receipt_sha256",
            "receipt_semantic_sha256", "checks_semantic_sha256",
            "artifact_audit_semantic_sha256", "input_identity_snapshot_sha256",
            "all_pass", "validation_evidence_eligible",
            "exactly_equals_supplied_pre_snapshot",
        },
        "release fresh input identity verification",
    )
    input_verifier_source = producer_records["input_identity_verifier"]["source"]
    require(
        fresh["execution_mode"] == "production_subprocess_required"
        and fresh["all_pass"] is True
        and fresh["validation_evidence_eligible"] is True
        and fresh["exactly_equals_supplied_pre_snapshot"] is True
        and fresh["input_identity_snapshot_sha256"] == pre["input_identity_snapshot_sha256"]
        and fresh["verifier_path"] == input_verifier_source["path"]
        and fresh["verifier_sha256"] == input_verifier_source["sha256"],
        "release fresh input identity verification summary mismatch",
    )
    for digest_name in (
        "receipt_sha256", "receipt_semantic_sha256", "checks_semantic_sha256",
        "artifact_audit_semantic_sha256", "input_identity_snapshot_sha256",
    ):
        _hex_sha256(fresh[digest_name], f"release fresh.{digest_name}")

    freezer_entry = producer_records["release_contract_freezer"]
    require(
        document["producer"] == {
            "role": "release_contract_freezer",
            "repo_relative_path": RELEASE_SOURCE_PATHS["release_contract_freezer"],
            "source_sha256": freezer_entry["source"]["sha256"],
            "immutable_copy_path": freezer_entry["snapshot"]["path"],
            "immutable_copy_sha256": freezer_entry["snapshot"]["sha256"],
        },
        "release producer is not bound to the authenticated freezer source/copy",
    )
    expected_top_entrypoints = {
        **expected_code_entrypoints,
        "canonical_manifest_copy": canonical_relative.as_posix(),
        "pre_input_identity_receipt_copy": pre_relative.as_posix(),
        "canonical_schema_copy": (
            Path("code_snapshot") / CANONICAL_SCHEMA_RELATIVE
        ).as_posix(),
    }
    require(
        document["entrypoints"] == expected_top_entrypoints,
        "release top-level exact entrypoint map mismatch",
    )

    if _skip_deep_verification_for_test:
        deep_verification = {"mode": "TEST_ONLY_SKIPPED", "all_pass": False}
    else:
        # This is deliberately last: no release-provided Python is imported
        # until every static byte/path/source/PRE invariant above has passed.
        deep_verification = _deep_verify_frozen_release(
            release_path.resolve(strict=True), release_identity["sha256"], identities,
            force_fresh=_force_deep_reverification,
        )
    return {
        "schema_name": "intersubmod.m2_release_binding",
        "schema_version": "1.0.0",
        "release_manifest": {
            "path": str(release_path.resolve(strict=True)),
            "sha256": release_identity["sha256"],
            "semantic_sha256": release_identity["semantic_sha256"],
            "sidecar": {
                "path": str(release_path.with_name(f"{release_path.name}.sha256").resolve()),
                "sha256": release_identity["sidecar_sha256"],
            },
        },
        "authority_mode": document["authority_mode"],
        "validation_evidence_eligible": deep_verification["all_pass"] is True,
        "canonical_input_manifest": copy_identity,
        "snapshot_sources": dict(sorted(identities.items())),
        "frozen_parameters": document["parameters"],
        "frozen_parameters_semantic_sha256": semantic_json_sha256(document["parameters"]),
        "deep_release_verification": deep_verification,
    }


def _require_source_identity(
    declared: Mapping[str, Any], expected: Mapping[str, str], label: str
) -> None:
    try:
        actual_path = Path(str(declared.get("path", ""))).resolve(strict=True)
    except OSError as exc:
        raise VerificationError(f"{label} source unavailable: {exc}") from exc
    require(str(actual_path) == expected["path"], f"{label} source path differs from release")
    require(declared.get("sha256") == expected["sha256"], f"{label} declared SHA differs from release")
    require(sha256_path(actual_path) == expected["sha256"], f"{label} physical SHA differs from release")


def verify_extraction_release_binding(
    full: Mapping[str, Any], binding: Mapping[str, Any]
) -> None:
    run_contract = full.get("run_contract") or {}
    require(run_contract.get("release_binding") == binding, "full extraction release binding mismatch")
    sources = binding["snapshot_sources"]
    _require_source_identity(run_contract.get("runner") or {}, sources["full_extraction_runner"], "extraction runner")
    _require_source_identity(
        run_contract.get("extractor") or {}, sources["extractor"], "extractor"
    )
    _require_source_identity(
        run_contract.get("lossless_read_contract") or {},
        sources["lossless_read_contract"],
        "lossless read contract",
    )
    require(
        run_contract.get("extractor_sha256") == sources["extractor"]["sha256"],
        "extraction producer SHA differs from release",
    )
    require(
        run_contract.get("manifest_sha256") == binding["canonical_input_manifest"]["sha256"],
        "extraction input manifest SHA differs from release",
    )
    frozen = binding["frozen_parameters"]
    extraction = frozen.get("extraction") or {}
    scheduler = frozen.get("scheduler") or {}
    require(
        {
            "mapq_min": run_contract.get("mapq_min"),
            "baseq_min": run_contract.get("baseq_min"),
            "bridge_thresholds": run_contract.get("bridge_thresholds"),
            "workers": run_contract.get("workers"),
            "samtools_threads": run_contract.get("samtools_threads"),
        } == extraction,
        "full extraction frozen parameter mismatch",
    )
    require(
        run_contract.get("task_timeout_seconds") == scheduler.get("task_timeout_seconds")
        and run_contract.get("timeout_grace_seconds") == scheduler.get("timeout_grace_seconds"),
        "full extraction frozen scheduler mismatch",
    )
    invocation = full.get("invocation") or {}
    require(
        invocation.get("max_new_tasks") in {
            scheduler.get("initial_batch_tasks"), scheduler.get("subsequent_batch_tasks")
        },
        "full extraction operational batch is not release-approved",
    )


def verify_ranking_release_binding(
    full: Mapping[str, Any], extraction_full: Mapping[str, Any], binding: Mapping[str, Any]
) -> None:
    run_contract = full.get("run_contract") or {}
    require(run_contract.get("release_binding") == binding, "full ranking release binding mismatch")
    require(
        ((extraction_full.get("run_contract") or {}).get("release_binding")) == binding,
        "extraction/ranking release contracts differ",
    )
    sources = binding["snapshot_sources"]
    _require_source_identity(run_contract.get("runner") or {}, sources["full_ranking_runner"], "ranking runner")
    _require_source_identity(run_contract.get("ranker") or {}, sources["ranker"], "ranker")
    _require_source_identity(
        run_contract.get("hypercube_solver") or {},
        sources["hypercube_solver"],
        "hypercube solver",
    )
    frozen = binding["frozen_parameters"]
    ranking = frozen.get("ranking") or {}
    scheduler = frozen.get("scheduler") or {}
    parameters = run_contract.get("parameters") or {}
    expected_parameters = {
        "structural_exact_pattern_minread_grid": ranking.get("structural_exact_pattern_minread_grid"),
        "primary_structural_exact_pattern_minread": ranking.get("primary_structural_exact_pattern_minread"),
        "scoring_minread": 1,
        "exact_k_max": ranking.get("exact_k_max"),
        "max_vertex_sets": ranking.get("max_vertex_sets"),
        "solver_time_limit_seconds_per_milp": ranking.get("solver_time_limit_seconds_per_milp"),
        "minimum_bq_error_rate": ranking.get("minimum_bq_error_rate"),
        "maximum_bq_error_rate": ranking.get("maximum_bq_error_rate"),
        "fixed_error_grid_conditional_binary_flip_probability": ranking.get("fixed_error_grid"),
        "conditional_candidate_ranking_bootstrap_replicates": ranking.get(
            "conditional_candidate_ranking_bootstrap_replicates"
        ),
        "conditional_candidate_ranking_bootstrap_base_seed": ranking.get(
            "conditional_candidate_ranking_bootstrap_seed"
        ),
        "tie_tolerance_log_likelihood": ranking.get("tie_tolerance_log_likelihood"),
    }
    require(parameters == expected_parameters, "full ranking frozen scientific parameters mismatch")
    require(
        run_contract.get("thresholds") == ranking.get("thresholds")
        and run_contract.get("component_bases") == ranking.get("component_bases")
        and run_contract.get("hp_families") == ranking.get("hp_families")
        and run_contract.get("workers") == ranking.get("workers"),
        "full ranking frozen scope/worker parameters mismatch",
    )
    require(
        run_contract.get("task_timeout_seconds") == scheduler.get("task_timeout_seconds")
        and run_contract.get("timeout_grace_seconds") == scheduler.get("timeout_grace_seconds"),
        "full ranking frozen scheduler mismatch",
    )
    invocation = full.get("invocation") or {}
    require(
        invocation.get("max_new_tasks") in {
            scheduler.get("initial_batch_tasks"), scheduler.get("subsequent_batch_tasks")
        }
        and invocation.get("aggregate_only") is False,
        "full ranking operational invocation is not release-approved",
    )


def verify_child_method_contract_and_ranker_source(
    child: Mapping[str, Any],
    expected_ranker_path: Path,
    expected_ranker_sha256: str,
    label: str,
) -> None:
    """Exact-compare one child contract and bind it to physical ranker bytes."""
    contract = (child.get("parameters") or {}).get("method_contract")
    require(contract == EXPECTED_METHOD_CONTRACT, f"method contract mismatch: {label}")
    ranker = (child.get("provenance") or {}).get("ranker") or {}
    try:
        declared_path = Path(str(ranker.get("path", ""))).resolve(strict=True)
        expected_path = expected_ranker_path.resolve(strict=True)
    except OSError as exc:
        raise VerificationError(f"ranker source path unavailable: {label}: {exc}") from exc
    require(declared_path == expected_path, f"child ranker path mismatch: {label}")
    require(ranker.get("sha256") == expected_ranker_sha256, f"child ranker SHA mismatch: {label}")
    require(
        sha256_path(expected_path) == expected_ranker_sha256,
        f"actual ranker source bytes mismatch declared SHA: {label}",
    )


def independent_runtime_summary(values: Iterable[float]) -> dict[str, Any]:
    """Independently compute exact nearest-rank timing statistics by full sort."""
    data = np.fromiter((float(value) for value in values), dtype=np.float64)
    n = int(data.size)
    if n == 0:
        return {"n": 0, "sum": 0.0, "max": None, "p50": None, "p95": None, "p99": None}
    require(bool(np.isfinite(data).all()) and not bool((data < 0.0).any()), "invalid runtime value")
    total = math.fsum(float(value) for value in data)
    data.sort()

    def nearest(probability: float) -> float:
        return float(data[min(n - 1, math.ceil(probability * n) - 1)])

    return {
        "n": n,
        "sum": total,
        "max": float(data[-1]),
        "p50": nearest(0.50),
        "p95": nearest(0.95),
        "p99": nearest(0.99),
    }


def verify_child_runtime_diagnostics(
    child: Mapping[str, Any], expected_dataset: str, expected_chrom: str
) -> dict[str, dict[str, list[float]]]:
    runtime = child.get("runtime_diagnostics") or {}
    require(runtime.get("schema_name") == "intersubmod.m2_unit_runtime_diagnostics", "runtime schema mismatch")
    require(runtime.get("schema_version") == "1.0.0", "runtime schema version mismatch")
    require(runtime.get("clock") == "time.monotonic_ns", "runtime clock mismatch")
    require(runtime.get("unit") == "seconds", "runtime unit mismatch")
    require(runtime.get("per_unit_output") == RUNTIME_SOURCE_NAME, "runtime output name mismatch")
    identity = (child.get("outputs") or {}).get(RUNTIME_SOURCE_NAME) or {}
    path = Path(str(identity.get("path")))
    values = {metric: [] for metric in RUNTIME_METRICS}
    all_values = {metric: [] for metric in RUNTIME_METRICS}
    invoked_values = {
        metric: []
        for metric in (
            "candidate_generation_elapsed_seconds", "likelihood_fit_elapsed_seconds",
        )
    }
    n_all = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(tuple(reader.fieldnames or ()) == RUNTIME_DIAGNOSTIC_FIELDS, "runtime columns mismatch")
        for row in reader:
            n_all += 1
            require(row["dataset"] == expected_dataset and row["chrom"] == expected_chrom, "runtime scope row mismatch")
            parsed = {metric: float(row[metric]) for metric in RUNTIME_METRICS}
            require(
                all(math.isfinite(value) and value >= 0.0 for value in parsed.values())
                and parsed["candidate_generation_elapsed_seconds"]
                + parsed["likelihood_fit_elapsed_seconds"]
                <= parsed["unit_total_elapsed_seconds"] + 1e-9,
                "runtime segment conservation failure",
            )
            for metric, value in parsed.items():
                all_values[metric].append(value)
            if row["structural_minread_role"] == "PRIMARY":
                for metric, value in parsed.items():
                    values[metric].append(value)
                invocation = {
                    "candidate_generation_elapsed_seconds": row["candidate_generation_invoked"],
                    "likelihood_fit_elapsed_seconds": row["likelihood_fit_invoked"],
                }
                require(all(value in {"true", "false"} for value in invocation.values()), "runtime invocation flag invalid")
                for metric, value in invocation.items():
                    if value == "true":
                        invoked_values[metric].append(parsed[metric])
    scopes = runtime.get("scopes") or {}
    primary = scopes.get("primary_unit_evaluations") or {}
    all_scope = scopes.get("all_structural_minread_unit_evaluations") or {}
    require(int(primary.get("n_unit_evaluations", -1)) == len(values[RUNTIME_METRICS[0]]), "primary runtime row count mismatch")
    require(int(all_scope.get("n_unit_evaluations", -1)) == n_all, "all-minread runtime row count mismatch")
    for metric in RUNTIME_METRICS:
        compare_exact(primary.get(metric) or {}, independent_runtime_summary(values[metric]), f"runtime/primary/{metric}")
        compare_exact(
            all_scope.get(metric) or {},
            independent_runtime_summary(all_values[metric]),
            f"runtime/all_minread/{metric}",
        )
    declared_invoked = runtime.get("primary_invoked_segment_scopes") or {}
    require(set(declared_invoked) == set(invoked_values), "runtime invoked-segment scope mismatch")
    for metric, metric_values in invoked_values.items():
        compare_exact(
            declared_invoked.get(metric) or {},
            independent_runtime_summary(metric_values),
            f"runtime/primary_invoked/{metric}",
        )
    return {"all_primary": values, "when_invoked": invoked_values}


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing JSON: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise VerificationError(f"invalid JSON {path}: {exc}") from exc
    require(isinstance(payload, dict), f"JSON root is not an object: {path}")
    return payload


def verify_sidecar(path: Path, receipt: Mapping[str, Any]) -> str:
    integrity = receipt.get("receipt_integrity") or {}
    require(integrity.get("scheme") == "external_sha256_sidecar_v1", f"bad sidecar scheme: {path}")
    sidecar = path.parent / str(integrity.get("sidecar_name") or f"{path.name}.sha256")
    require(sidecar.is_file(), f"missing sidecar: {sidecar}")
    try:
        fields = sidecar.read_text(encoding="ascii", errors="strict").strip().split()
    except (OSError, UnicodeError) as exc:
        raise VerificationError(f"unreadable sidecar {sidecar}: {exc}") from exc
    actual = sha256_path(path)
    require(len(fields) == 2 and fields[0] == actual and fields[1] == path.name, f"sidecar mismatch: {path}")
    return actual


def verify_identity(identity: Mapping[str, Any], expected_parent: Path | None = None) -> dict[str, Any]:
    try:
        path = Path(str(identity["path"])).absolute()
        declared_size = int(identity["size_bytes"])
        declared_sha = str(identity["sha256"])
    except (KeyError, TypeError, ValueError) as exc:
        raise VerificationError(f"malformed file identity: {identity}") from exc
    if expected_parent is not None:
        require(
            path.parent.resolve(strict=True) == expected_parent.resolve(strict=True),
            f"output escapes task directory: {path}",
        )
    try:
        before_path = os.lstat(path)
    except OSError as exc:
        raise VerificationError(f"missing declared output: {path}: {exc}") from exc
    require(
        stat.S_ISREG(before_path.st_mode)
        and not stat.S_ISLNK(before_path.st_mode)
        and before_path.st_nlink == 1,
        f"declared output is not a single-link regular non-symlink file: {path}",
    )
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise VerificationError(f"cannot open declared output: {path}: {exc}") from exc
    digest = hashlib.sha256()
    actual_size = 0
    try:
        before_fd = os.fstat(descriptor)
        require(_same_stat_state(before_path, before_fd), f"declared output changed before read: {path}")
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            actual_size += len(block)
            digest.update(block)
        after_fd = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = os.lstat(path)
    require(
        _same_stat_state(before_fd, after_fd) and _same_stat_state(after_fd, after_path),
        f"declared output changed while being hashed: {path}",
    )
    actual_sha = digest.hexdigest()
    require(actual_size == declared_size, f"size mismatch: {path}")
    require(actual_sha == declared_sha, f"SHA-256 mismatch: {path}")
    return {
        "path": str(path.resolve(strict=True)), "size_bytes": actual_size,
        "sha256": actual_sha,
    }


def authenticate_receipt(
    path: Path, schema_name: str, schema_version: str, *, require_all_checks: bool = True
) -> tuple[dict[str, Any], str]:
    receipt = load_json(path)
    require(receipt.get("schema_name") == schema_name, f"schema_name mismatch: {path}")
    require(receipt.get("schema_version") == schema_version, f"schema_version mismatch: {path}")
    require(receipt.get("all_pass") is True, f"receipt is not PASS: {path}")
    if require_all_checks:
        checks = receipt.get("checks") or {}
        require(bool(checks), f"receipt has no checks: {path}")
        failed = sorted(key for key, value in checks.items() if value is not True)
        require(not failed, f"receipt contains failed/non-boolean checks {failed}: {path}")
    return receipt, verify_sidecar(path, receipt)


def expected_keys(datasets: Sequence[str], chromosomes: Sequence[str]) -> list[tuple[str, str]]:
    return [(dataset, chrom) for dataset in datasets for chrom in chromosomes]


def validate_result_index(
    results: Sequence[Mapping[str, Any]], datasets: Sequence[str], chromosomes: Sequence[str], label: str
) -> dict[tuple[str, str], Mapping[str, Any]]:
    expected = expected_keys(datasets, chromosomes)
    require(len(results) == len(expected), f"{label}: expected {len(expected)} results, got {len(results)}")
    index: dict[tuple[str, str], Mapping[str, Any]] = {}
    for result in results:
        key = (str(result.get("dataset")), str(result.get("chrom")))
        require(key not in index, f"{label}: duplicate task {key}")
        index[key] = result
    require(set(index) == set(expected), f"{label}: task set differs from declared Cartesian scope")
    return index


def _integer_mapping(payload: Mapping[str, Any]) -> dict[str, int]:
    return {
        str(key): int(value)
        for key, value in payload.items()
        if isinstance(value, int) and not isinstance(value, bool)
    }


def extraction_count_checks(counts: Mapping[str, Any], label: str) -> dict[str, bool]:
    def value(key: str) -> int:
        return int(counts.get(key, 0))

    checks = {
        "raw_alignment_class_conserved": value("raw_overlapping_alignments")
        == value("alignment_class_primary") + value("alignment_class_secondary")
        + value("alignment_class_supplementary") + value("alignment_class_unmapped"),
        "raw_filter_funnel_conserved": value("raw_overlapping_alignments")
        == value("excluded_by_flag") + value("mapq_rejected_after_flag")
        + value("canonical_eligible_alignments"),
        "eligible_rows_conserved": value("canonical_eligible_alignments")
        == value("molecule_sparse_rows_written") == value("unique_molecule_ids"),
        "fixed_calls_not_greater_than_sparse_calls": value("fixed_ra_calls") <= value("site_call_rows_sparse"),
        "alt_calls_not_greater_than_fixed_calls": value("alt_calls") <= value("fixed_ra_calls"),
        "sidecar_matches_equal_eligible": value("sidecar_exact_matches") == value("canonical_eligible_alignments"),
        "sidecar_missing_zero": value("sidecar_missing") == 0,
    }
    failed = sorted(key for key, passed in checks.items() if not passed)
    require(not failed, f"{label}: independent count conservation failed: {failed}")
    return checks


def verify_component_cell(cell: Mapping[str, Any], label: str) -> dict[str, Any]:
    distribution = {int(key): int(value) for key, value in (cell.get("k_component_sites_distribution") or {}).items()}
    n_components = int(cell.get("n_components", 0))
    checks = {
        "distribution_sums_to_components": sum(distribution.values()) == n_components,
        "singletons_match_distribution": int(cell.get("n_singletons_k1", 0)) == distribution.get(1, 0),
        "multisite_matches_distribution": int(cell.get("n_multisite_k_gt1", 0))
        == sum(value for key, value in distribution.items() if key > 1),
        "max_k_matches_distribution": int(cell.get("max_k_component_sites", 0))
        == (max(distribution) if distribution else 0),
    }
    failed = sorted(key for key, passed in checks.items() if not passed)
    require(not failed, f"{label}: component distribution failed: {failed}")
    return checks


def aggregate_extraction_children(
    child_receipts: Mapping[tuple[str, str], Mapping[str, Any]],
    datasets: Sequence[str],
    thresholds: Sequence[int],
    task_statuses: Mapping[tuple[str, str], str],
) -> dict[str, Any]:
    totals: Counter[str] = Counter()
    by_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    task_counts: dict[str, Counter[str]] = defaultdict(Counter)
    component_sums: dict[str, dict[str, Counter[str]]] = {
        basis: {str(threshold): Counter() for threshold in thresholds} for basis in LINKAGE_BASES
    }
    component_max: dict[str, dict[str, dict[str, int]]] = {
        basis: {
            str(threshold): {"max_k_component_sites": 0, "max_k": 0} for threshold in thresholds
        }
        for basis in LINKAGE_BASES
    }
    k_distributions: dict[str, dict[str, Counter[int]]] = {
        basis: {str(threshold): Counter() for threshold in thresholds} for basis in LINKAGE_BASES
    }
    phase_totals = {
        "known_phase_set_chromosome_units_by_hp_family": Counter(),
        "known_ps_active_site_memberships_by_hp_family": Counter(),
        "missing_ps_active_sites_by_hp_family": Counter(),
    }
    legacy: dict[str, dict[str, Counter[str]]] = defaultdict(lambda: defaultdict(Counter))

    for (dataset, chrom), receipt in child_receipts.items():
        counts = _integer_mapping(receipt.get("counts") or {})
        totals.update(counts)
        by_dataset[dataset].update(counts)
        task_counts[dataset][task_statuses[(dataset, chrom)]] += 1
        phase = receipt.get("phase_set_contract_counts") or {}
        for target, source in (
            ("known_phase_set_chromosome_units_by_hp_family", "known_phase_sets_by_hp_family"),
            ("known_ps_active_site_memberships_by_hp_family", "known_ps_active_site_memberships_by_hp_family"),
            ("missing_ps_active_sites_by_hp_family", "missing_ps_active_sites_by_hp_family"),
        ):
            phase_totals[target].update(_integer_mapping(phase.get(source) or {}))
        for family, threshold_values in (receipt.get("legacy_cross_phase_set_aggregation_audit") or {}).items():
            for threshold, metrics in threshold_values.items():
                legacy[str(family)][str(threshold)].update(_integer_mapping(metrics))
        summaries = receipt.get("component_summary_by_linkage_basis") or {}
        require(set(summaries) == set(LINKAGE_BASES), f"{dataset}: linkage-basis set mismatch")
        for basis in LINKAGE_BASES:
            require(set(summaries[basis]) == {str(value) for value in thresholds}, f"{dataset}/{basis}: threshold set mismatch")
            for threshold in thresholds:
                key = str(threshold)
                cell = summaries[basis][key]
                verify_component_cell(cell, f"{dataset}/{basis}/{key}")
                for name, value in cell.items():
                    if (
                        isinstance(value, int) and not isinstance(value, bool)
                        and name not in {"max_k_component_sites", "max_k"}
                    ):
                        component_sums[basis][key][name] += value
                for name in ("max_k_component_sites", "max_k"):
                    component_max[basis][key][name] = max(
                        component_max[basis][key][name], int(cell.get(name, 0))
                    )
                k_distributions[basis][key].update({
                    int(k): int(value)
                    for k, value in (cell.get("k_component_sites_distribution") or {}).items()
                })

    components = {
        basis: {
            key: {
                **dict(component_sums[basis][key]),
                **component_max[basis][key],
                "k_component_sites_distribution": {
                    str(k): value for k, value in sorted(k_distributions[basis][key].items())
                },
                "k_distribution": {
                    str(k): value for k, value in sorted(k_distributions[basis][key].items())
                },
            }
            for key in component_sums[basis]
        }
        for basis in LINKAGE_BASES
    }
    return {
        "counts": dict(totals),
        "component_summary_by_linkage_basis": components,
        "component_summary_by_threshold": components["pooled"],
        "phase_set_contract_totals": {
            key: dict(sorted(counter.items())) for key, counter in phase_totals.items()
        },
        "legacy_cross_phase_set_aggregation_audit": {
            family: {threshold: dict(metrics) for threshold, metrics in sorted(values.items())}
            for family, values in sorted(legacy.items())
        },
        "by_dataset": {
            dataset: {
                "task_status_counts": dict(task_counts[dataset]),
                "counts": dict(by_dataset[dataset]),
            }
            for dataset in datasets
        },
    }


def compare_exact(actual: Any, expected: Any, label: str) -> None:
    require(actual == expected, f"aggregate mismatch at {label}")


def verify_extraction_bundle(
    extraction_root: Path,
    datasets: Sequence[str],
    chromosomes: Sequence[str],
    release_binding: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    full_path = extraction_root / "full_extraction_receipt.json"
    full, full_sha = authenticate_receipt(
        full_path, "intersubmod.m2_full_extraction_receipt", "1.2.0"
    )
    if release_binding is not None:
        verify_extraction_release_binding(full, release_binding)
    scope = full.get("scope") or {}
    require(scope.get("datasets") == list(datasets), "full extraction dataset scope mismatch")
    require(scope.get("chromosomes") == list(chromosomes), "full extraction chromosome scope mismatch")
    require(int(scope.get("expected_tasks", -1)) == len(datasets) * len(chromosomes), "full extraction expected_tasks mismatch")
    results = full.get("results") or []
    result_index = validate_result_index(results, datasets, chromosomes, "full extraction")
    orchestration = None
    if release_binding is not None:
        orchestration = verify_orchestration_stage(
            extraction_root, full, stage="extraction", datasets=datasets,
            chromosomes=chromosomes, release_binding=release_binding,
        )
    thresholds = tuple(int(value) for value in (full.get("run_contract") or {}).get("bridge_thresholds") or ())
    require(thresholds, "full extraction has no bridge thresholds")
    children: dict[tuple[str, str], dict[str, Any]] = {}
    task_index = []
    for dataset, chrom in expected_keys(datasets, chromosomes):
        result = result_index[(dataset, chrom)]
        require(result.get("status") in {"PASS", "REUSED_PASS"}, f"extraction task not passing: {dataset}/{chrom}")
        child_path = extraction_root / "samples" / dataset / chrom / "receipt.json"
        child, child_sha = authenticate_receipt(
            child_path, "intersubmod.lossless_read_linkage_chromosome_receipt", "1.2.0"
        )
        child_scope = child.get("scope") or {}
        require(child_scope.get("dataset") == dataset and child_scope.get("chrom") == chrom, f"child extraction scope mismatch: {dataset}/{chrom}")
        require(result.get("receipt") == child, f"embedded extraction child differs from disk: {dataset}/{chrom}")
        if release_binding is not None:
            producer = (child.get("provenance") or {}).get("extractor") or {}
            _require_source_identity(
                producer,
                release_binding["snapshot_sources"]["extractor"],
                f"child extractor {dataset}/{chrom}",
            )
            child_manifest = (child.get("provenance") or {}).get("manifest") or {}
            _require_source_identity(
                child_manifest,
                release_binding["canonical_input_manifest"],
                f"child manifest {dataset}/{chrom}",
            )
        extraction_count_checks(child.get("counts") or {}, f"{dataset}/{chrom}")
        outputs = child.get("outputs") or {}
        required_matches = {
            suffix: [name for name in outputs if name.endswith(suffix)] for suffix in REQUIRED_EXTRACTION_SUFFIXES
        }
        require(all(len(matches) == 1 for matches in required_matches.values()), f"required extraction outputs missing/duplicated: {dataset}/{chrom}")
        verified_outputs = {
            name: verify_identity(identity, child_path.parent)
            for name, identity in outputs.items()
        }
        children[(dataset, chrom)] = child
        task_index.append({
            "dataset": dataset,
            "chrom": chrom,
            "child_receipt_sha256": child_sha,
            "n_outputs_verified": len(verified_outputs),
            "n_sSNV": int(child_scope.get("n_sSNV", 0)),
            "canonical_eligible_alignments": int((child.get("counts") or {}).get("canonical_eligible_alignments", 0)),
            "all_pass": True,
        })
    recomputed = aggregate_extraction_children(
        children,
        datasets,
        thresholds,
        {key: str(result_index[key]["status"]) for key in expected_keys(datasets, chromosomes)},
    )
    compare_exact(full.get("aggregate"), recomputed, "full_extraction.aggregate")
    return {
        "receipt_path": str(full_path.resolve()),
        "receipt_sha256": full_sha,
        "n_tasks": len(task_index),
        "task_index": task_index,
        "recomputed_aggregate": recomputed,
        "children": children,
        "full_receipt": full,
        "thresholds": list(thresholds),
        "orchestration": orchestration,
    }


def new_rank_cell() -> dict[str, Any]:
    return {
        "sums": Counter(),
        "counters": {field: Counter() for field in COUNTER_FIELDS},
        "partial_pattern_funnel": {},
    }


def merge_numeric_tree(target: dict[str, Any], source: Mapping[str, Any], path: str = "") -> None:
    for key, value in source.items():
        current_path = f"{path}/{key}" if path else str(key)
        if key == "definitions":
            if key in target:
                require(target[key] == value, f"partial-funnel definition drift at {current_path}")
            else:
                target[key] = dict(value)
        elif isinstance(value, Mapping):
            destination = target.setdefault(key, {})
            require(isinstance(destination, dict), f"partial-funnel type drift at {current_path}")
            merge_numeric_tree(destination, value, current_path)
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            target[key] = target.get(key, 0) + value
        else:
            if key in target:
                require(target[key] == value, f"partial-funnel scalar drift at {current_path}")
            else:
                target[key] = value


def add_rank_summary(cell: dict[str, Any], summary: Mapping[str, Any], partial: Mapping[str, Any] | None) -> None:
    for field in SUM_FIELDS:
        cell["sums"][field] += int(summary.get(field, 0))
    for field in COUNTER_FIELDS:
        cell["counters"][field].update(_integer_mapping(summary.get(field) or {}))
    if partial is not None:
        merge_numeric_tree(cell["partial_pattern_funnel"], partial)


def freeze_rank_cell(cell: Mapping[str, Any]) -> dict[str, Any]:
    return {
        **dict(cell["sums"]),
        **{field: dict(sorted(cell["counters"][field].items())) for field in COUNTER_FIELDS},
        "partial_pattern_funnel": cell["partial_pattern_funnel"],
    }


def rank_conservation(cell: Mapping[str, Any]) -> dict[str, bool]:
    n_units = int(cell.get("n_component_hp_units", 0))
    return {
        "selection_status_sum_equals_units": sum((cell.get("selection_status_counts") or {}).values()) == n_units,
        "unique_tied_abstain_partition_units": int(cell.get("quality_primary_unique_vertex_units", 0))
        + int(cell.get("quality_primary_tied_vertex_units", 0))
        + int(cell.get("rank_abstain_units", 0)) == n_units,
        "solver_partition_units": int(cell.get("solver_complete_units", 0))
        + int(cell.get("solver_incomplete_or_not_run_units", 0)) == n_units,
        "projection_funnel": int(cell.get("molecule_component_projections", 0))
        == int(cell.get("informative_scoring_molecules", 0)) + int(cell.get("all_x_excluded_molecules", 0)),
        "structural_scoring_funnel": int(cell.get("informative_scoring_molecules", 0))
        == int(cell.get("structural_retained_molecules", 0)) + int(cell.get("below_minread_scoring_molecules", 0)),
        "bq_equals_informative": int(cell.get("bq_scoring_molecules", 0))
        == int(cell.get("informative_scoring_molecules", 0)),
        "raw_T_not_less_than_V": int(cell.get("raw_tree_candidates_T_complete_units", 0))
        >= int(cell.get("distinct_vertex_sets_V_complete_units", 0)),
        "topology_partition": int(cell.get("coarse_topology_class_unique_units", 0))
        + int(cell.get("coarse_topology_multiple_class_units", 0))
        == int(cell.get("topology_evaluated_units", 0)),
        "partial_coverage": int(cell.get("partial_groups_covered", 0))
        + int(cell.get("partial_groups_unsatisfied", 0))
        == int(cell.get("partial_group_coverage_denominator", 0)),
        "partial_unsatisfied_zero": int(cell.get("partial_groups_unsatisfied", 0)) == 0,
        "k_route_partition": sum((cell.get("k_route_counts") or {}).values()) == n_units,
        "effective_k_mass": int(cell.get("k_component_sites_total", 0))
        == int(cell.get("k_observed_alt_active_total", 0))
        + int(cell.get("not_structural_alt_active_sites_total", 0)),
    }


def compare_rank_cell(declared: Mapping[str, Any], recomputed: Mapping[str, Any], label: str) -> dict[str, bool]:
    for field in SUM_FIELDS:
        compare_exact(int(declared.get(field, 0)), int(recomputed.get(field, 0)), f"{label}/{field}")
    for field in COUNTER_FIELDS:
        compare_exact(declared.get(field) or {}, recomputed.get(field) or {}, f"{label}/{field}")
    compare_exact(
        declared.get("partial_pattern_funnel") or {},
        recomputed.get("partial_pattern_funnel") or {},
        f"{label}/partial_pattern_funnel",
    )
    checks = rank_conservation(recomputed)
    failed = sorted(key for key, value in checks.items() if not value)
    require(not failed, f"{label}: independent rank conservation failed: {failed}")
    return checks


def add_child_rank_aggregates(
    child: Mapping[str, Any], dataset: str,
    primary_global: dict[tuple[str, str], dict[str, Any]],
    primary_dataset: dict[str, dict[tuple[str, str], dict[str, Any]]],
    sensitivity_global: dict[str, dict[tuple[str, str], dict[str, Any]]],
    sensitivity_dataset: dict[str, dict[str, dict[tuple[str, str], dict[str, Any]]]],
) -> None:
    partial_primary = child.get("partial_pattern_funnel_by_linkage_basis_threshold") or {}
    for basis, thresholds in (child.get("aggregate_by_linkage_basis_threshold") or {}).items():
        for threshold, summary in thresholds.items():
            key = (str(basis), str(threshold))
            partial = (partial_primary.get(str(basis)) or {}).get(str(threshold))
            add_rank_summary(primary_global.setdefault(key, new_rank_cell()), summary, partial)
            add_rank_summary(primary_dataset.setdefault(dataset, {}).setdefault(key, new_rank_cell()), summary, partial)
    for minread, payload in (child.get("sensitivity_by_structural_exact_pattern_minread") or {}).items():
        partial_nested = payload.get("partial_pattern_funnel_by_linkage_basis_threshold") or {}
        for basis, thresholds in (payload.get("by_linkage_basis_threshold") or {}).items():
            for threshold, summary in thresholds.items():
                key = (str(basis), str(threshold))
                partial = (partial_nested.get(str(basis)) or {}).get(str(threshold))
                add_rank_summary(
                    sensitivity_global.setdefault(str(minread), {}).setdefault(key, new_rank_cell()),
                    summary,
                    partial,
                )
                add_rank_summary(
                    sensitivity_dataset.setdefault(str(minread), {}).setdefault(dataset, {}).setdefault(key, new_rank_cell()),
                    summary,
                    partial,
                )


def nested_rank_cells(cells: Mapping[tuple[str, str], Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = defaultdict(dict)
    for (basis, threshold), cell in sorted(cells.items()):
        output[basis][threshold] = freeze_rank_cell(cell)
    return dict(output)


def compare_rank_cell_set(
    declared: Mapping[str, Any], recomputed: Mapping[str, Any], label: str
) -> tuple[int, list[dict[str, Any]]]:
    require(set(declared) == set(recomputed), f"{label}: linkage-basis set mismatch")
    checks = []
    n = 0
    for basis in sorted(recomputed):
        require(set(declared[basis]) == set(recomputed[basis]), f"{label}/{basis}: threshold set mismatch")
        for threshold in sorted(recomputed[basis], key=int):
            cell_checks = compare_rank_cell(
                declared[basis][threshold], recomputed[basis][threshold], f"{label}/{basis}/{threshold}"
            )
            checks.append({"path": f"{label}/{basis}/{threshold}", "checks": cell_checks, "all_pass": True})
            n += 1
    return n, checks


def _canonical_unit_key(row: Mapping[str, str], primary_minread: int) -> str:
    chrom_number = int(str(row["chrom"]).removeprefix("chr"))
    return (
        f"{row['dataset']}|chr{chrom_number:02d}|{row['component_basis']}|PS={row['phase_set']}|"
        f"B{int(row['threshold']):03d}|{row['component_id']}|HP{row['family']}|M{primary_minread}"
    )


def _strict_bool(value: str, label: str) -> bool:
    require(value in {"true", "false"}, f"invalid boolean {label}: {value!r}")
    return value == "true"


def _serialized_12g_tolerance(*values: float) -> float:
    """Bound normal decimal serialization plus floating summation noise.

    Numeric scalar columns are written with ``.12g``.  JSON mixture weights
    retain Python's full round-trip precision, so only the persisted scalar is
    rounded.  Four half-units at the twelfth significant digit is deliberately
    conservative while remaining far below any scientifically meaningful
    likelihood adjustment.
    """
    scale = max((abs(float(value)) for value in values), default=1.0)
    return 2e-11 * max(1.0, scale)


def _profile_source_groups(
    path: Path,
    required_columns: Sequence[str],
    dataset: str,
    chrom: str,
    primary_minread: int,
    *,
    source_kind: str,
) -> Iterator[tuple[tuple[str, ...], str, list[dict[str, str]]]]:
    """Stream contiguous primary-unit rows without retaining a whole child."""
    current_order: tuple[str, ...] | None = None
    current_unit: str | None = None
    current_rows: list[dict[str, str]] = []
    previous_order: tuple[str, ...] | None = None
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        missing = sorted(set(required_columns) - set(fields))
        require(not missing, f"{source_kind} columns missing {missing}: {path}")
        for row in reader:
            try:
                row_minread = int(row["structural_exact_pattern_minread"])
            except (TypeError, ValueError) as exc:
                raise VerificationError(f"invalid minread in {source_kind}: {path}") from exc
            if row_minread != primary_minread or row["component_basis"] not in PRIMARY_BASES:
                continue
            require(
                row["dataset"] == dataset and row["chrom"] == chrom,
                f"{source_kind} scope mismatch: {path}",
            )
            require(
                row["phase_set"] not in {"", ".", "None", "NA", "UNKNOWN"},
                f"{source_kind} primary row lacks known PS: {path}",
            )
            if source_kind == "pattern":
                require(
                    row["structural_minread_role"] == "PRIMARY",
                    f"primary pattern row has non-primary role: {path}",
                )
            # This is the producer's physical order: the values are stringified
            # before sorting.  Keeping the same literal order permits a lockstep
            # merge while the canonical key remains the human-readable identity.
            order = (
                row["component_basis"], row["phase_set"], row["threshold"],
                row["component_id"], row["family"],
            )
            unit = _canonical_unit_key(row, primary_minread)
            if current_order is None:
                current_order, current_unit = order, unit
            elif order != current_order:
                require(
                    previous_order is None or current_order > previous_order,
                    f"{source_kind} unit order regression: {path}",
                )
                yield current_order, str(current_unit), current_rows
                previous_order = current_order
                current_order, current_unit, current_rows = order, unit, []
            require(unit == current_unit, f"{source_kind} unit identity drift: {path}")
            current_rows.append(row)
        if current_order is not None:
            require(
                previous_order is None or current_order > previous_order,
                f"{source_kind} final unit order regression: {path}",
            )
            yield current_order, str(current_unit), current_rows


def _parse_profile_patterns(
    rows: Sequence[Mapping[str, str]], unit_key: str
) -> tuple[int, list[tuple[str, tuple[int, ...], int]], int]:
    require(bool(rows), f"profile unit has no pattern rows: {unit_key}")
    parsed: list[tuple[str, tuple[int, ...], int]] = []
    seen: set[tuple[str, tuple[int, ...]]] = set()
    k_values: set[int] = set()
    scoring_molecules = 0
    for row in rows:
        try:
            k = int(row["k"])
            count = int(row["n_molecules"])
        except (TypeError, ValueError) as exc:
            raise VerificationError(f"invalid pattern k/count: {unit_key}") from exc
        pattern = row["pattern_rax"]
        require(k > 0 and len(pattern) == k, f"pattern k mismatch: {unit_key}")
        require(set(pattern) <= {"R", "A", "X"}, f"invalid R/A/X pattern: {unit_key}")
        require(count > 0, f"non-positive pattern count: {unit_key}")
        quality_cells = row["fixed_base_qualities"].split(",")
        require(len(quality_cells) == k, f"BQ vector length mismatch: {unit_key}")
        qualities: list[int] = []
        for symbol, cell in zip(pattern, quality_cells):
            if symbol == "X":
                require(cell == "", f"X call unexpectedly has BQ: {unit_key}")
                qualities.append(-1)
            else:
                try:
                    quality = int(cell)
                except (TypeError, ValueError) as exc:
                    raise VerificationError(f"fixed R/A call lacks integer BQ: {unit_key}") from exc
                require(quality >= 0, f"negative fixed-call BQ: {unit_key}")
                qualities.append(quality)
        quality_tuple = tuple(qualities)
        duplicate_key = (pattern, quality_tuple)
        require(duplicate_key not in seen, f"duplicate pattern/BQ group: {unit_key}")
        seen.add(duplicate_key)
        informative = set(pattern) != {"X"}
        require(
            _strict_bool(row["scoring_eligible"], f"{unit_key}/scoring_eligible") == informative,
            f"scoring-eligible flag disagrees with R/A/X pattern: {unit_key}",
        )
        if informative:
            parsed.append((pattern, quality_tuple, count))
            scoring_molecules += count
        k_values.add(k)
    require(len(k_values) == 1, f"mixed k within profile unit: {unit_key}")
    require(bool(parsed), f"candidate unit has no informative scoring patterns: {unit_key}")
    return next(iter(k_values)), parsed, scoring_molecules


def _candidate_profile_at_persisted_pi(
    row: Mapping[str, str],
    patterns: Sequence[tuple[str, tuple[int, ...], int]],
    k: int,
    minimum_error_rate: float,
    maximum_error_rate: float,
    unit_key: str,
) -> dict[str, Any]:
    """Independently recompute conditional BQ LL and its simplex KKT gap."""
    try:
        states_payload = json.loads(row["states_json"])
        weights_payload = json.loads(row["mixture_weights_json"])
    except json.JSONDecodeError as exc:
        raise VerificationError(f"invalid states/pi JSON: {unit_key}") from exc
    require(isinstance(states_payload, list) and states_payload, f"empty candidate states: {unit_key}")
    require(isinstance(weights_payload, dict), f"candidate pi is not a mapping: {unit_key}")
    states: list[int] = []
    for state in states_payload:
        require(isinstance(state, Mapping), f"invalid state record: {unit_key}")
        bitmask = state.get("bitmask")
        require(
            isinstance(bitmask, int) and not isinstance(bitmask, bool) and 0 <= bitmask < (1 << k),
            f"state bitmask outside k-cube: {unit_key}",
        )
        expected_state = "".join("A" if bitmask & (1 << bit) else "R" for bit in range(k))
        require(state.get("state_rax") == expected_state, f"state_rax disagrees with bitmask: {unit_key}")
        states.append(bitmask)
    require(len(states) == len(set(states)), f"duplicate candidate state: {unit_key}")
    require(0 in states, f"candidate state set omits root: {unit_key}")
    expected_vertex_set_id = hashlib.sha256(json.dumps(
        {"k": int(k), "vertices": sorted(states)},
        sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")).hexdigest()
    require(
        row["vertex_set_id"] == expected_vertex_set_id,
        f"vertex_set_id disagrees with candidate states: {unit_key}",
    )
    try:
        parsed_weights = {int(key): float(value) for key, value in weights_payload.items()}
    except (TypeError, ValueError) as exc:
        raise VerificationError(f"invalid candidate pi: {unit_key}") from exc
    require(set(parsed_weights) == set(states), f"pi/state key mismatch: {unit_key}")
    weights = [parsed_weights[state] for state in states]
    require(
        all(math.isfinite(value) and value >= 0.0 for value in weights),
        f"candidate pi is non-finite or negative: {unit_key}",
    )
    emission_rows: list[list[float]] = []
    count_values: list[float] = []
    for pattern, qualities, count in patterns:
        emissions: list[float] = []
        for state in states:
            probability = 1.0
            for bit, (symbol, quality) in enumerate(zip(pattern, qualities)):
                if symbol == "X":
                    continue
                error = min(max(10.0 ** (-quality / 10.0), minimum_error_rate), maximum_error_rate)
                denominator = 1.0 - (2.0 * error / 3.0)
                match_probability = (1.0 - error) / denominator
                flip_probability = (error / 3.0) / denominator
                expected_alt = bool(state & (1 << bit))
                probability *= (
                    match_probability if expected_alt == (symbol == "A") else flip_probability
                )
            emissions.append(probability)
        emission_rows.append(emissions)
        count_values.append(float(count))
    # This matrix is bounded to one candidate in one unit and immediately
    # released.  NumPy's operation order mirrors the persisted objective while
    # the implementation remains independent of the production ranker/solver.
    emission = np.asarray(emission_rows, dtype=np.float64)
    counts = np.asarray(count_values, dtype=np.float64)
    pi = np.asarray(weights, dtype=np.float64)
    mixture_probabilities = emission @ pi
    require(
        bool(np.isfinite(mixture_probabilities).all())
        and bool((mixture_probabilities > 0.0).all()),
        f"non-positive profile likelihood denominator: {unit_key}",
    )
    log_likelihood = float(np.dot(counts, np.log(mixture_probabilities)))
    gradients = emission.T @ (counts / mixture_probabilities)
    weighted_gradient = float(np.dot(pi, gradients))
    gap = max(0.0, float(np.max(gradients) - weighted_gradient))
    simplex_residual = max(
        abs(math.fsum(weights) - 1.0),
        max(0.0, -min(weights)),
    )
    return {
        "log_likelihood": log_likelihood,
        "global_log_likelihood_gap_bound": gap,
        "simplex_residual": simplex_residual,
        "n_states": len(states),
        "n_scoring_molecules": int(counts.sum()),
        "n_emission_cells": int(emission.size),
    }


def _require_recomputed_float(
    stored_raw: str, recomputed: float, label: str
) -> tuple[float, float]:
    try:
        stored = float(stored_raw)
    except (TypeError, ValueError) as exc:
        raise VerificationError(f"invalid persisted float {label}: {stored_raw!r}") from exc
    require(math.isfinite(stored) and math.isfinite(recomputed), f"non-finite value {label}")
    delta = abs(stored - recomputed)
    require(
        delta <= _serialized_12g_tolerance(stored, recomputed),
        f"{label} differs from independent recomputation: delta={delta:.17g}",
    )
    return stored, delta


def _verify_profile_unit(
    unit_key: str,
    pattern_rows: Sequence[Mapping[str, str]],
    candidate_rows: Sequence[Mapping[str, str]],
    *,
    minimum_error_rate: float,
    maximum_error_rate: float,
    tie_tolerance: float,
    max_vertex_sets: int,
) -> dict[str, Any]:
    require(bool(candidate_rows), f"empty candidate group: {unit_key}")
    require(
        len(candidate_rows) <= max_vertex_sets,
        f"candidate group exceeds declared max_vertex_sets: {unit_key}",
    )
    vertex_ids = [row["vertex_set_id"] for row in candidate_rows]
    require(len(vertex_ids) == len(set(vertex_ids)), f"duplicate vertex_set_id: {unit_key}")
    k, patterns, scoring_molecules = _parse_profile_patterns(pattern_rows, unit_key)
    audits: list[dict[str, Any]] = []
    max_ll_delta = 0.0
    max_gap_delta = 0.0
    max_simplex_delta = 0.0
    peak_states = 0
    peak_emission_cells = 0
    for row in candidate_rows:
        recomputed = _candidate_profile_at_persisted_pi(
            row, patterns, k, minimum_error_rate, maximum_error_rate, unit_key
        )
        stored_ll, ll_delta = _require_recomputed_float(
            row["primary_log_likelihood"], recomputed["log_likelihood"],
            f"{unit_key}/{row['vertex_set_id']}/primary_log_likelihood",
        )
        stored_gap, gap_delta = _require_recomputed_float(
            row["global_log_likelihood_gap_bound"],
            recomputed["global_log_likelihood_gap_bound"],
            f"{unit_key}/{row['vertex_set_id']}/global_log_likelihood_gap_bound",
        )
        stored_simplex, simplex_delta = _require_recomputed_float(
            row["simplex_residual"], recomputed["simplex_residual"],
            f"{unit_key}/{row['vertex_set_id']}/simplex_residual",
        )
        monotone = _strict_bool(row["fit_monotone"], f"{unit_key}/fit_monotone")
        converged = _strict_bool(row["fit_converged"], f"{unit_key}/fit_converged")
        independently_certified = (
            recomputed["global_log_likelihood_gap_bound"]
            <= EXPECTED_OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE
            and recomputed["simplex_residual"] <= EXPECTED_SIMPLEX_RESIDUAL_TOLERANCE
            and monotone
        )
        require(
            converged == independently_certified,
            f"fit_converged disagrees with independently recomputed certificate: {unit_key}",
        )
        require(
            row["optimizer_status"].startswith("CERTIFIED_") == converged,
            f"optimizer_status disagrees with certificate: {unit_key}",
        )
        audits.append({
            "row": row,
            "log_likelihood": recomputed["log_likelihood"],
            "stored_log_likelihood": stored_ll,
            "converged": converged,
            "monotone": monotone,
        })
        max_ll_delta = max(max_ll_delta, ll_delta)
        max_gap_delta = max(max_gap_delta, gap_delta)
        max_simplex_delta = max(max_simplex_delta, simplex_delta)
        peak_states = max(peak_states, int(recomputed["n_states"]))
        peak_emission_cells = max(peak_emission_cells, int(recomputed["n_emission_cells"]))

    best_ll = max(audit["log_likelihood"] for audit in audits)
    likelihood_terms = [
        math.exp(max(-745.0, audit["log_likelihood"] - best_ll)) for audit in audits
    ]
    normalizer = math.fsum(likelihood_terms)
    require(normalizer > 0.0 and math.isfinite(normalizer), f"invalid candidate normalizer: {unit_key}")
    top_indices = {
        index for index, audit in enumerate(audits)
        if best_ll - audit["log_likelihood"] <= tie_tolerance
    }
    optimizer_pass = all(audit["converged"] and audit["monotone"] for audit in audits)
    max_relative_delta = 0.0
    for index, (audit, term) in enumerate(zip(audits, likelihood_terms)):
        row = audit["row"]
        _, relative_delta = _require_recomputed_float(
            row["relative_likelihood_weight"], term / normalizer,
            f"{unit_key}/{row['vertex_set_id']}/relative_likelihood_weight",
        )
        max_relative_delta = max(max_relative_delta, relative_delta)
        expected_winner = optimizer_pass and index in top_indices
        expected_tied = expected_winner and len(top_indices) > 1
        require(
            _strict_bool(row["is_winner"], f"{unit_key}/is_winner") == expected_winner,
            f"winner flag differs from recomputed LL ranking: {unit_key}",
        )
        require(
            _strict_bool(row["is_tied_winner"], f"{unit_key}/is_tied_winner") == expected_tied,
            f"tie flag differs from recomputed LL ranking: {unit_key}",
        )
    return {
        "n_pattern_rows": len(pattern_rows),
        "n_scoring_molecules": scoring_molecules,
        "n_candidates": len(candidate_rows),
        "max_abs_ll_delta": max_ll_delta,
        "max_abs_relative_weight_delta": max_relative_delta,
        "max_abs_gap_delta": max_gap_delta,
        "max_abs_simplex_residual_delta": max_simplex_delta,
        "peak_states_per_candidate": peak_states,
        "peak_emission_cells_per_candidate": peak_emission_cells,
    }


def verify_child_profile_likelihood(
    child: Mapping[str, Any],
    dataset: str,
    chrom: str,
    expected_parameters: Mapping[str, Any],
) -> dict[str, Any]:
    """Stream/join one child's primary pattern and candidate tables."""
    parameters = child.get("parameters") or {}
    parameter_keys = (
        "minimum_bq_error_rate", "maximum_bq_error_rate",
        "tie_tolerance_log_likelihood", "primary_structural_exact_pattern_minread",
        "max_vertex_sets",
    )
    for key in parameter_keys:
        require(key in parameters and key in expected_parameters, f"missing profile parameter: {key}")
        compare_exact(parameters[key], expected_parameters[key], f"{dataset}/{chrom}/parameters/{key}")
    optimizer_contract = parameters.get("optimizer_contract") or {}
    require(
        float(optimizer_contract.get("global_log_likelihood_gap_tolerance", math.nan))
        == EXPECTED_OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE,
        f"optimizer KKT tolerance drift: {dataset}/{chrom}",
    )
    require(
        float(optimizer_contract.get("simplex_residual_tolerance", math.nan))
        == EXPECTED_SIMPLEX_RESIDUAL_TOLERANCE,
        f"optimizer simplex tolerance drift: {dataset}/{chrom}",
    )
    minimum_error_rate = float(parameters["minimum_bq_error_rate"])
    maximum_error_rate = float(parameters["maximum_bq_error_rate"])
    tie_tolerance = float(parameters["tie_tolerance_log_likelihood"])
    primary_minread = int(parameters["primary_structural_exact_pattern_minread"])
    max_vertex_sets = int(parameters["max_vertex_sets"])
    require(
        0.0 < minimum_error_rate <= maximum_error_rate < 0.5
        and tie_tolerance >= 0.0 and primary_minread > 0 and max_vertex_sets > 0,
        f"invalid profile parameter range: {dataset}/{chrom}",
    )
    outputs = child.get("outputs") or {}
    pattern_identity = outputs.get(PATTERN_SOURCE_NAME) or {}
    candidate_identity = outputs.get(CANDIDATE_SOURCE_NAME) or {}
    pattern_path = Path(str(pattern_identity.get("path", "")))
    candidate_path = Path(str(candidate_identity.get("path", "")))
    require(pattern_path.is_file(), f"missing child pattern table: {dataset}/{chrom}")
    require(candidate_path.is_file(), f"missing child candidate table: {dataset}/{chrom}")
    pattern_groups = iter(_profile_source_groups(
        pattern_path, PROFILE_PATTERN_REQUIRED_COLUMNS, dataset, chrom, primary_minread,
        source_kind="pattern",
    ))
    candidate_groups = _profile_source_groups(
        candidate_path, PROFILE_CANDIDATE_REQUIRED_COLUMNS, dataset, chrom, primary_minread,
        source_kind="candidate",
    )
    current_pattern = next(pattern_groups, None)
    totals = Counter()
    maxima = {
        "max_abs_ll_delta": 0.0,
        "max_abs_relative_weight_delta": 0.0,
        "max_abs_gap_delta": 0.0,
        "max_abs_simplex_residual_delta": 0.0,
        "peak_pattern_rows_per_unit": 0,
        "peak_candidates_per_unit": 0,
        "peak_states_per_candidate": 0,
        "peak_emission_cells_per_candidate": 0,
    }
    n_units = 0
    for candidate_order, candidate_unit, candidate_rows in candidate_groups:
        while current_pattern is not None and current_pattern[0] < candidate_order:
            current_pattern = next(pattern_groups, None)
        require(
            current_pattern is not None and current_pattern[0] == candidate_order
            and current_pattern[1] == candidate_unit,
            f"candidate unit has no matching primary pattern group: {candidate_unit}",
        )
        audit = _verify_profile_unit(
            candidate_unit, current_pattern[2], candidate_rows,
            minimum_error_rate=minimum_error_rate,
            maximum_error_rate=maximum_error_rate,
            tie_tolerance=tie_tolerance,
            max_vertex_sets=max_vertex_sets,
        )
        n_units += 1
        totals["n_pattern_rows"] += int(audit["n_pattern_rows"])
        totals["n_scoring_molecules"] += int(audit["n_scoring_molecules"])
        totals["n_candidates"] += int(audit["n_candidates"])
        for key in (
            "max_abs_ll_delta", "max_abs_relative_weight_delta", "max_abs_gap_delta",
            "max_abs_simplex_residual_delta",
        ):
            maxima[key] = max(maxima[key], float(audit[key]))
        maxima["peak_states_per_candidate"] = max(
            maxima["peak_states_per_candidate"], int(audit["peak_states_per_candidate"])
        )
        maxima["peak_emission_cells_per_candidate"] = max(
            maxima["peak_emission_cells_per_candidate"],
            int(audit["peak_emission_cells_per_candidate"]),
        )
        maxima["peak_pattern_rows_per_unit"] = max(
            maxima["peak_pattern_rows_per_unit"], len(current_pattern[2])
        )
        maxima["peak_candidates_per_unit"] = max(
            maxima["peak_candidates_per_unit"], len(candidate_rows)
        )
        current_pattern = next(pattern_groups, None)
    return {
        "dataset": dataset,
        "chrom": chrom,
        "n_units": n_units,
        **dict(totals),
        **maxima,
        "all_profile_likelihoods_and_certificates_match": True,
        "all_relative_weights_match": True,
        "all_winner_tie_partitions_match": True,
    }


def aggregate_profile_likelihood_audits(
    children: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    sum_fields = ("n_units", "n_candidates", "n_pattern_rows", "n_scoring_molecules")
    max_float_fields = (
        "max_abs_ll_delta", "max_abs_relative_weight_delta", "max_abs_gap_delta",
        "max_abs_simplex_residual_delta",
    )
    max_integer_fields = (
        "peak_pattern_rows_per_unit", "peak_candidates_per_unit", "peak_states_per_candidate",
        "peak_emission_cells_per_candidate",
    )
    return {
        "n_children": len(children),
        **{key: sum(int(child.get(key, 0)) for child in children) for key in sum_fields},
        **{
            key: max((float(child.get(key, 0.0)) for child in children), default=0.0)
            for key in max_float_fields
        },
        **{
            key: max((int(child.get(key, 0)) for child in children), default=0)
            for key in max_integer_fields
        },
        "all_profile_likelihoods_and_certificates_match": all(
            child.get("all_profile_likelihoods_and_certificates_match") is True for child in children
        ),
        "all_relative_weights_match": all(
            child.get("all_relative_weights_match") is True for child in children
        ),
        "all_winner_tie_partitions_match": all(
            child.get("all_winner_tie_partitions_match") is True for child in children
        ),
        "count_unit_definitions": {
            "n_children": "authenticated dataset-chromosome child receipts",
            "n_units": "primary HP-family x known-PS x read-linked component x threshold units with candidates",
            "n_candidates": "distinct complete candidate vertex sets evaluated within those units",
            "n_pattern_rows": "distinct R/A/X plus fixed-BQ groups joined to candidate-bearing units",
            "n_scoring_molecules": (
                "sum of n_molecules over joined pattern groups; a molecule projected to multiple units or "
                "thresholds contributes once per projection and is not a globally unique molecule count"
            ),
        },
        "numeric_contract": {
            "persisted_scalar_format": ".12g",
            "independent_comparison_absolute_tolerance": "2e-11 * max(1, abs(stored), abs(recomputed))",
            "global_log_likelihood_gap_tolerance": EXPECTED_OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE,
            "simplex_residual_tolerance": EXPECTED_SIMPLEX_RESIDUAL_TOLERANCE,
        },
        "streaming_memory_contract": (
            "one child opened at a time; pattern and candidate TSVs lockstep-streamed; only one "
            "primary unit's pattern rows and candidate rows plus one candidate's pattern-by-state "
            "emission matrix retained; no whole-child or whole-run table materialization"
        ),
        "child_summaries": list(children),
    }


def transformed_candidate_rows(
    child_receipts: Mapping[tuple[str, str], Mapping[str, Any]],
    datasets: Sequence[str], chromosomes: Sequence[str], primary_minread: int,
) -> Iterator[dict[str, str]]:
    previous_global_unit: str | None = None
    for dataset, chrom in expected_keys(datasets, chromosomes):
        child = child_receipts[(dataset, chrom)]
        identity = (child.get("outputs") or {}).get(CANDIDATE_SOURCE_NAME)
        require(identity is not None, f"missing child candidate table: {dataset}/{chrom}")
        source_path = Path(str(identity["path"]))
        current_unit: str | None = None
        group: list[dict[str, str]] = []

        def emit(unit_key: str, rows: list[dict[str, str]]) -> Iterator[dict[str, str]]:
            nonlocal previous_global_unit
            require(previous_global_unit is None or unit_key >= previous_global_unit, "child candidate global unit order regression")
            previous_global_unit = unit_key
            ordered = sorted(rows, key=lambda item: item["vertex_set_id"])
            require(len({row["vertex_set_id"] for row in ordered}) == len(ordered), f"duplicate vertex_set_id: {unit_key}")
            optimizer_pass = all(
                row["fit_converged"].lower() == "true" and row["fit_monotone"].lower() == "true"
                for row in ordered
            )
            best_ll = max(float(row["primary_log_likelihood"]) for row in ordered)
            for index, row in enumerate(ordered, start=1):
                states = json.loads(row["states_json"])
                state_map = {str(item["bitmask"]): item["state_rax"] for item in states}
                role_map = {str(item["bitmask"]): item["roles"] for item in states}
                winner = row["is_winner"].lower() == "true"
                tied = row["is_tied_winner"].lower() == "true"
                winner_status = (
                    "ABSTAIN_UNIT_OPTIMIZER" if not optimizer_pass
                    else "TIED_WINNER" if winner and tied
                    else "UNIQUE_WINNER" if winner
                    else "NON_WINNER"
                )
                yield {
                    "unit_key": unit_key,
                    "dataset": row["dataset"],
                    "chrom": row["chrom"],
                    "component_id": row["component_id"],
                    "threshold": row["threshold"],
                    "hp_family": row["family"],
                    "ps": row["phase_set"],
                    "candidate_id": f"C{index:06d}",
                    "vertex_set_id": row["vertex_set_id"],
                    "vertex_states": json.dumps(state_map, sort_keys=True, separators=(",", ":")),
                    "vertex_roles": json.dumps(role_map, sort_keys=True, separators=(",", ":")),
                    "parent_choice_count": row["parent_choice_count"],
                    "profile_log_likelihood": row["primary_log_likelihood"],
                    "relative_log_likelihood": format(float(row["primary_log_likelihood"]) - best_ll, ".12g"),
                    "mixture_weights_pi": row["mixture_weights_json"],
                    "winner_status": winner_status,
                    "tie_group": "TOP_TIE" if tied else "",
                    "coarse_topology_class": row["coarse_topology_classes_json"],
                    "candidate_set_complete": "true",
                }

        with gzip.open(source_path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if int(row["structural_exact_pattern_minread"]) != primary_minread:
                    continue
                require(row["dataset"] == dataset and row["chrom"] == chrom, f"candidate source scope mismatch: {source_path}")
                require(row["component_basis"] in PRIMARY_BASES, f"non-primary basis in candidate source: {source_path}")
                require(row["phase_set"] not in {"", ".", "None", "NA", "UNKNOWN"}, f"missing PS in candidate source: {source_path}")
                unit_key = _canonical_unit_key(row, primary_minread)
                if current_unit is None:
                    current_unit = unit_key
                elif unit_key != current_unit:
                    require(unit_key > current_unit, f"child candidate unit order regression: {source_path}")
                    yield from emit(current_unit, group)
                    current_unit, group = unit_key, []
                group.append(row)
            if current_unit is not None:
                yield from emit(current_unit, group)


def verify_candidate_table(
    ranking_root: Path,
    metadata: Mapping[str, Any],
    child_receipts: Mapping[tuple[str, str], Mapping[str, Any]],
    datasets: Sequence[str], chromosomes: Sequence[str], primary_minread: int,
) -> dict[str, Any]:
    require(tuple(metadata.get("columns") or ()) == CANDIDATE_TABLE_COLUMNS, "candidate metadata columns mismatch")
    identity = verify_identity(metadata, ranking_root)
    path = Path(identity["path"])
    expected = transformed_candidate_rows(child_receipts, datasets, chromosomes, primary_minread)
    semantic = hashlib.sha256()
    n_rows = 0
    n_units = 0
    previous_unit: str | None = None
    candidate_index = 0
    unit_statuses: list[str] = []

    def close_unit() -> None:
        nonlocal n_units
        if not unit_statuses:
            return
        statuses = Counter(unit_statuses)
        valid = (
            set(statuses) == {"ABSTAIN_UNIT_OPTIMIZER"}
            or (statuses.get("UNIQUE_WINNER", 0) == 1 and statuses.get("TIED_WINNER", 0) == 0 and statuses.get("ABSTAIN_UNIT_OPTIMIZER", 0) == 0)
            or (statuses.get("TIED_WINNER", 0) >= 2 and statuses.get("UNIQUE_WINNER", 0) == 0 and statuses.get("ABSTAIN_UNIT_OPTIMIZER", 0) == 0)
        )
        require(valid, f"candidate winner partition invalid for {previous_unit}: {dict(statuses)}")
        n_units += 1

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(tuple(reader.fieldnames or ()) == CANDIDATE_TABLE_COLUMNS, "candidate table header mismatch")
        sentinel = object()
        for actual, expected_row in zip_longest(reader, expected, fillvalue=sentinel):
            require(actual is not sentinel and expected_row is not sentinel, "candidate table row count differs from child reconstruction")
            require(actual == expected_row, f"candidate row differs from independently reconstructed source at row {n_rows + 1}")
            unit = actual["unit_key"]
            if previous_unit is None or unit != previous_unit:
                if previous_unit is not None:
                    close_unit()
                require(previous_unit is None or unit > previous_unit, "candidate table unit sort regression")
                previous_unit = unit
                candidate_index = 0
                unit_statuses = []
            candidate_index += 1
            require(actual["candidate_id"] == f"C{candidate_index:06d}", f"candidate_id sequence mismatch: {unit}")
            require(actual["candidate_set_complete"] == "true", f"candidate_set_complete is not true: {unit}")
            unit_statuses.append(actual["winner_status"])
            semantic.update(
                json.dumps(actual, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8") + b"\n"
            )
            n_rows += 1
        close_unit()
    require(n_rows == int(metadata.get("n_rows", -1)), "candidate metadata n_rows mismatch")
    require(n_units == int(metadata.get("n_units", -1)), "candidate metadata n_units mismatch")
    require(semantic.hexdigest() == metadata.get("semantic_sha256"), "candidate semantic SHA-256 mismatch")
    return {
        **identity,
        "semantic_sha256": semantic.hexdigest(),
        "n_rows": n_rows,
        "n_units": n_units,
        "all_rows_match_independent_child_reconstruction": True,
        "winner_partitions_conserved": True,
    }


def verify_ranking_bundle(
    ranking_root: Path,
    extraction_root: Path,
    extraction_children: Mapping[tuple[str, str], Mapping[str, Any]],
    datasets: Sequence[str], chromosomes: Sequence[str],
    extraction_full: Mapping[str, Any] | None = None,
    release_binding: Mapping[str, Any] | None = None,
    extraction_orchestration: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    full_path = ranking_root / "full_ranking_receipt.json"
    full, full_sha = authenticate_receipt(full_path, "intersubmod.m2_full_ranking_receipt", "2.0.0")
    extraction_full_path = extraction_root / "full_extraction_receipt.json"
    extraction_full_sha = sha256_path(extraction_full_path)
    declared_full_extraction = (full.get("run_contract") or {}).get(
        "full_extraction_receipt"
    ) or {}
    require(
        Path(str(declared_full_extraction.get("path", ""))).resolve()
        == extraction_full_path.resolve()
        and declared_full_extraction.get("sha256") == extraction_full_sha,
        "ranking run contract full-extraction receipt identity mismatch",
    )
    declared_upstream = full.get("upstream_extraction_receipt") or {}
    require(
        Path(str(declared_upstream.get("path", ""))).resolve()
        == extraction_full_path.resolve()
        and declared_upstream.get("sha256") == extraction_full_sha,
        "ranking top-level full-extraction receipt identity mismatch",
    )
    if release_binding is not None:
        require(extraction_full is not None, "release verification requires the full extraction receipt")
        verify_ranking_release_binding(full, extraction_full, release_binding)
        orchestration = verify_orchestration_stage(
            ranking_root, full, stage="ranking", datasets=datasets,
            chromosomes=chromosomes, release_binding=release_binding,
            extraction_root=extraction_root,
            extraction_children=extraction_children,
            parent_extraction=extraction_orchestration,
        )
    else:
        orchestration = None
    scope = full.get("scope") or {}
    require(scope.get("datasets") == list(datasets), "full ranking dataset scope mismatch")
    require(scope.get("chromosomes") == list(chromosomes), "full ranking chromosome scope mismatch")
    require(int(scope.get("expected_tasks", -1)) == len(datasets) * len(chromosomes), "full ranking expected_tasks mismatch")
    run_contract = full.get("run_contract") or {}
    require(
        run_contract.get("method_contract") == EXPECTED_METHOD_CONTRACT,
        "full ranking method contract mismatch",
    )
    declared_ranker = run_contract.get("ranker") or {}
    try:
        actual_ranker_path = Path(str(declared_ranker.get("path", ""))).resolve(strict=True)
    except OSError as exc:
        raise VerificationError(f"full ranking ranker source unavailable: {exc}") from exc
    actual_ranker_sha256 = sha256_path(actual_ranker_path)
    require(
        declared_ranker.get("sha256") == actual_ranker_sha256,
        "full ranking ranker SHA does not match actual source bytes",
    )
    result_index = validate_result_index(full.get("results") or [], datasets, chromosomes, "full ranking")
    task_index_declared = validate_result_index(full.get("task_index") or [], datasets, chromosomes, "full ranking task_index")
    child_receipts: dict[tuple[str, str], dict[str, Any]] = {}
    task_index = []
    primary_global: dict[tuple[str, str], dict[str, Any]] = {}
    primary_dataset: dict[str, dict[tuple[str, str], dict[str, Any]]] = {}
    sensitivity_global: dict[str, dict[tuple[str, str], dict[str, Any]]] = {}
    sensitivity_dataset: dict[str, dict[str, dict[tuple[str, str], dict[str, Any]]]] = {}
    input_totals: Counter[str] = Counter()
    call_totals: Counter[str] = Counter()
    input_by_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    calls_by_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    hp_by_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    runtime_values = {metric: [] for metric in RUNTIME_METRICS}
    runtime_invoked_values = {
        metric: []
        for metric in (
            "candidate_generation_elapsed_seconds", "likelihood_fit_elapsed_seconds",
        )
    }
    profile_child_audits: list[dict[str, Any]] = []
    expected_profile_parameters = run_contract.get("parameters") or {}

    for dataset, chrom in expected_keys(datasets, chromosomes):
        result = result_index[(dataset, chrom)]
        require(result.get("status") in {"PASS", "REUSED_PASS"}, f"ranking task not passing: {dataset}/{chrom}")
        child_path = ranking_root / "samples" / dataset / chrom / "receipt.json"
        child, child_sha = authenticate_receipt(
            child_path, "intersubmod.m2_symbolic_patterns_vertex_rank_receipt", "2.0.0"
        )
        child_scope = child.get("scope") or {}
        require(child_scope.get("dataset") == [dataset] and child_scope.get("chrom") == [chrom], f"child ranking scope mismatch: {dataset}/{chrom}")
        require(sorted(child_scope.get("component_bases") or []) == sorted(PRIMARY_BASES), f"child ranking basis mismatch: {dataset}/{chrom}")
        require(child_scope.get("component_basis_mode") == "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY", f"child ranking mode mismatch: {dataset}/{chrom}")
        verify_child_method_contract_and_ranker_source(
            child, actual_ranker_path, actual_ranker_sha256, f"{dataset}/{chrom}"
        )
        declared_rank = result.get("rank_receipt") or {}
        require(Path(str(declared_rank.get("path"))).resolve() == child_path.resolve(), f"full ranking child path mismatch: {dataset}/{chrom}")
        require(declared_rank.get("sha256") == child_sha, f"full ranking child SHA mismatch: {dataset}/{chrom}")
        extraction_path = extraction_root / "samples" / dataset / chrom / "receipt.json"
        upstream = (child.get("provenance") or {}).get("upstream_extraction_receipt") or {}
        require(Path(str(upstream.get("path"))).resolve() == extraction_path.resolve(), f"ranking upstream extraction path mismatch: {dataset}/{chrom}")
        require(upstream.get("sha256") == sha256_path(extraction_path), f"ranking upstream extraction SHA mismatch: {dataset}/{chrom}")
        extraction_outputs = extraction_children[(dataset, chrom)].get("outputs") or {}
        input_files = child.get("input_files") or {}
        require(set(input_files) == {
            name for name in extraction_outputs if name.endswith(REQUIRED_EXTRACTION_SUFFIXES)
        }, f"ranking input file set differs from extraction receipt: {dataset}/{chrom}")
        for name, identity in input_files.items():
            compare_exact(identity, extraction_outputs[name], f"{dataset}/{chrom}/input_files/{name}")
            verify_identity(identity, extraction_path.parent)
        for identity in (child.get("outputs") or {}).values():
            verify_identity(identity, child_path.parent)
        profile_child_audits.append(
            verify_child_profile_likelihood(
                child, dataset, chrom, expected_profile_parameters
            )
        )
        child_runtime = verify_child_runtime_diagnostics(child, dataset, chrom)
        for metric in RUNTIME_METRICS:
            runtime_values[metric].extend(child_runtime["all_primary"][metric])
        for metric in runtime_invoked_values:
            runtime_invoked_values[metric].extend(child_runtime["when_invoked"][metric])
        declared_task = task_index_declared[(dataset, chrom)]
        for flag in (
            "all_pass", "parameters_match_extraction", "input_hashes_match_extraction",
            "upstream_outputs_verified", "no_cross_ps_pattern_pooling", "known_ps_never_mixed",
            "missing_ps_separate_diagnostic", "runtime_diagnostics_contract_valid",
            "method_contract_matches", "ranker_source_bound",
        ):
            require(declared_task.get(flag) is True, f"declared task_index flag false: {dataset}/{chrom}/{flag}")
        require(declared_task.get("schema_version") == "2.0.0", f"task_index schema mismatch: {dataset}/{chrom}")
        child_receipts[(dataset, chrom)] = child
        add_child_rank_aggregates(
            child, dataset, primary_global, primary_dataset, sensitivity_global, sensitivity_dataset
        )
        input_counts = child.get("input_counts") or {}
        integers = _integer_mapping(input_counts)
        input_totals.update(integers)
        input_by_dataset[dataset].update(integers)
        calls = _integer_mapping(input_counts.get("selected_sparse_call_code_counts") or {})
        call_totals.update(calls)
        calls_by_dataset[dataset].update(calls)
        hp_by_dataset[dataset].update(_integer_mapping(input_counts.get("hp_family_rows") or {}))
        task_index.append({
            "dataset": dataset,
            "chrom": chrom,
            "rank_receipt_sha256": child_sha,
            "n_outputs_verified": len(child.get("outputs") or {}),
            "n_component_hp_units": int((child.get("aggregate") or {}).get("n_component_hp_units", 0)),
            "all_pass": True,
        })

    recomputed_primary = nested_rank_cells(primary_global)
    aggregate = full.get("aggregate") or {}
    declared_primary = aggregate.get("by_linkage_basis_threshold") or {}
    n_cells, conservation_cells = compare_rank_cell_set(
        declared_primary, recomputed_primary, "aggregate/by_linkage_basis_threshold"
    )
    compare_exact(aggregate.get("input_call_funnel") or {}, dict(input_totals), "aggregate/input_call_funnel")
    compare_exact(aggregate.get("input_sparse_call_code_counts") or {}, dict(sorted(call_totals.items())), "aggregate/input_sparse_call_code_counts")
    for dataset in datasets:
        declared_dataset = (aggregate.get("by_dataset") or {}).get(dataset) or {}
        expected_dataset_primary = nested_rank_cells(primary_dataset.get(dataset, {}))
        count, cells = compare_rank_cell_set(
            declared_dataset.get("by_linkage_basis_threshold") or {},
            expected_dataset_primary,
            f"aggregate/by_dataset/{dataset}/by_linkage_basis_threshold",
        )
        n_cells += count
        conservation_cells.extend(cells)
        compare_exact(declared_dataset.get("input_call_funnel") or {}, dict(input_by_dataset[dataset]), f"aggregate/by_dataset/{dataset}/input_call_funnel")
        compare_exact(declared_dataset.get("input_sparse_call_code_counts") or {}, dict(sorted(calls_by_dataset[dataset].items())), f"aggregate/by_dataset/{dataset}/input_sparse_call_code_counts")
        compare_exact(declared_dataset.get("input_hp_family_rows") or {}, dict(sorted(hp_by_dataset[dataset].items())), f"aggregate/by_dataset/{dataset}/input_hp_family_rows")

    declared_sensitivity = aggregate.get("by_structural_exact_pattern_minread") or {}
    require(set(declared_sensitivity) == set(sensitivity_global), "sensitivity minread grid mismatch")
    recomputed_sensitivity: dict[str, Any] = {}
    for minread in sorted(sensitivity_global, key=int):
        expected_global = nested_rank_cells(sensitivity_global[minread])
        count, cells = compare_rank_cell_set(
            (declared_sensitivity[minread] or {}).get("by_linkage_basis_threshold") or {},
            expected_global,
            f"aggregate/by_structural_exact_pattern_minread/{minread}",
        )
        n_cells += count
        conservation_cells.extend(cells)
        recomputed_sensitivity[minread] = {"by_linkage_basis_threshold": expected_global, "by_dataset": {}}
        declared_minread_datasets = (declared_sensitivity[minread] or {}).get("by_dataset") or {}
        require(set(declared_minread_datasets) == set(datasets), f"sensitivity dataset set mismatch: minread={minread}")
        for dataset in datasets:
            expected_cells = nested_rank_cells(sensitivity_dataset[minread].get(dataset, {}))
            count, cells = compare_rank_cell_set(
                (declared_minread_datasets[dataset] or {}).get("by_linkage_basis_threshold") or {},
                expected_cells,
                f"aggregate/by_structural_exact_pattern_minread/{minread}/by_dataset/{dataset}",
            )
            n_cells += count
            conservation_cells.extend(cells)
            recomputed_sensitivity[minread]["by_dataset"][dataset] = {
                "by_linkage_basis_threshold": expected_cells
            }

    primary_minread = int(((full.get("run_contract") or {}).get("parameters") or {}).get("primary_structural_exact_pattern_minread", -1))
    require(primary_minread > 0, "missing primary structural exact-pattern minread")
    candidate_audit = verify_candidate_table(
        ranking_root,
        full.get("candidate_table") or {},
        child_receipts,
        datasets,
        chromosomes,
        primary_minread,
    )
    rankable_units = sum(
        int(cell.get("solver_complete_units", 0))
        for thresholds in recomputed_primary.values() for cell in thresholds.values()
    )
    require(candidate_audit["n_units"] == rankable_units, "candidate table does not cover all solver-complete units")
    declared_runtime = full.get("runtime_diagnostics") or {}
    require(
        declared_runtime.get("schema_name") == "intersubmod.m2_full_primary_runtime_diagnostics"
        and declared_runtime.get("schema_version") == "1.0.0"
        and declared_runtime.get("clock") == "time.monotonic_ns"
        and declared_runtime.get("unit") == "seconds",
        "full runtime metadata mismatch",
    )
    require(
        int(declared_runtime.get("n_child_runtime_files", -1)) == len(datasets) * len(chromosomes),
        "full runtime child count mismatch",
    )
    n_runtime_units = len(runtime_values[RUNTIME_METRICS[0]])
    require(
        int(declared_runtime.get("n_unit_evaluations", -1)) == n_runtime_units,
        "full runtime unit count mismatch",
    )
    recomputed_runtime = {
        metric: independent_runtime_summary(values)
        for metric, values in runtime_values.items()
    }
    for metric in RUNTIME_METRICS:
        compare_exact(
            (declared_runtime.get("metrics") or {}).get(metric) or {},
            recomputed_runtime[metric],
            f"full_runtime/{metric}",
        )
    recomputed_runtime_invoked = {
        metric: independent_runtime_summary(values)
        for metric, values in runtime_invoked_values.items()
    }
    for metric, summary in recomputed_runtime_invoked.items():
        compare_exact(
            (declared_runtime.get("metrics_when_invoked") or {}).get(metric) or {},
            summary,
            f"full_runtime_when_invoked/{metric}",
        )
    profile_likelihood_audit = aggregate_profile_likelihood_audits(profile_child_audits)
    require(
        profile_likelihood_audit["n_children"] == len(datasets) * len(chromosomes),
        "profile likelihood child count mismatch",
    )
    require(
        profile_likelihood_audit["n_units"] == candidate_audit["n_units"]
        and profile_likelihood_audit["n_candidates"] == candidate_audit["n_rows"],
        "profile likelihood audit does not cover the full primary candidate table",
    )
    return {
        "receipt_path": str(full_path.resolve()),
        "receipt_sha256": full_sha,
        "n_tasks": len(task_index),
        "task_index": task_index,
        "n_aggregate_cells_recomputed": n_cells,
        "all_aggregate_cells_conserved": all(row["all_pass"] for row in conservation_cells),
        "aggregate_cell_checks": conservation_cells,
        "recomputed": {
            "input_call_funnel": dict(input_totals),
            "input_sparse_call_code_counts": dict(sorted(call_totals.items())),
            "by_linkage_basis_threshold": recomputed_primary,
            "by_structural_exact_pattern_minread": recomputed_sensitivity,
        },
        "candidate_table": candidate_audit,
        "profile_likelihood_independent_recomputation": profile_likelihood_audit,
        "runtime_diagnostics": {
            "n_child_runtime_files": len(datasets) * len(chromosomes),
            "n_unit_evaluations": n_runtime_units,
            "metrics": recomputed_runtime,
            "metrics_when_invoked": recomputed_runtime_invoked,
            "all_child_and_full_runtime_summaries_independently_recomputed": True,
        },
        "method_contract_verification": {
            "contract": EXPECTED_METHOD_CONTRACT,
            "ranker_source_path": str(actual_ranker_path),
            "ranker_source_sha256": actual_ranker_sha256,
            "n_children_exactly_matched_and_source_bound": len(task_index),
            "all_children_exactly_matched_and_source_bound": (
                len(task_index) == len(datasets) * len(chromosomes)
            ),
            "verification_mode": "STATIC_EXACT_CONTRACT_AND_ACTUAL_SOURCE_SHA_BINDING",
        },
        "orchestration": orchestration,
    }


def verify_bundle(
    extraction_root: Path,
    ranking_root: Path,
    datasets: Sequence[str] = DATASETS,
    chromosomes: Sequence[str] = AUTOSOMES,
    release_contract_manifest: Path | None = None,
    post_input_identity_receipt: Path | None = None,
) -> dict[str, Any]:
    release_binding = None
    if release_contract_manifest is not None:
        release_binding = load_release_contract_binding(
            release_contract_manifest,
            required_sources={"full_verifier": VERIFIER},
        )
        require(
            post_input_identity_receipt is not None,
            "formal release verification requires --post-input-identity-receipt",
        )
        post_input_identity = verify_post_input_identity_receipt(
            post_input_identity_receipt, release_binding
        )
    else:
        post_input_identity = None
    extraction = verify_extraction_bundle(
        extraction_root, datasets, chromosomes, release_binding
    )
    ranking = verify_ranking_bundle(
        ranking_root,
        extraction_root,
        extraction["children"],
        datasets,
        chromosomes,
        extraction.get("full_receipt"),
        release_binding,
        extraction.get("orchestration"),
    )
    # Child payloads are needed only during verification and would make the
    # independent receipt unnecessarily large.
    extraction_public = {
        key: value for key, value in extraction.items()
        if key not in {"children", "full_receipt"}
    }
    receipt = {
        "schema_name": "intersubmod.m2_full_independent_verification",
        "schema_version": "1.0.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "scope": {
            "datasets": list(datasets),
            "chromosomes": list(chromosomes),
            "expected_tasks": len(datasets) * len(chromosomes),
            "n_technical_datasets": len(datasets),
        },
        "verification_independence": {
            "imports_production_aggregator": False,
            "imports_production_ranker": False,
            "reads_bam": False,
            "recomputes_from_child_receipts_and_tables": True,
            "candidate_table_compared_rowwise_to_independent_child_reconstruction": True,
            "method_contract_verification_mode": (
                "STATIC_EXACT_CONTRACT_AND_ACTUAL_SOURCE_SHA_BINDING"
            ),
            "profile_likelihood_recomputed_from_patterns_states_pi": True,
            "profile_likelihood_claim_boundary": (
                "Independently streams persisted pattern R/A/X+BQ+count groups and candidate states/pi; "
                "recomputes the conditional symmetric-substitution BQ mixture LL at persisted pi, the "
                "simplex and Frank-Wolfe/KKT certificate, relative weights, and winner/tie partition. "
                "The certificate exploits objective concavity; the verifier does not rerun SLSQP."
            ),
        },
        "release_binding": release_binding,
        "post_input_identity": post_input_identity,
        "extraction": extraction_public,
        "ranking": ranking,
        "checks": {
            "all_expected_extraction_tasks_verified": extraction_public["n_tasks"] == len(datasets) * len(chromosomes),
            "all_expected_ranking_tasks_verified": ranking["n_tasks"] == len(datasets) * len(chromosomes),
            "extraction_full_aggregate_exactly_recomputed": True,
            "ranking_core_aggregates_exactly_recomputed": True,
            "ranking_all_recomputed_cells_conserve": ranking["all_aggregate_cells_conserved"],
            "candidate_table_physical_and_semantic_hashes_verified": True,
            "candidate_table_rows_match_child_sources": ranking["candidate_table"]["all_rows_match_independent_child_reconstruction"],
            "candidate_table_winner_partitions_conserve": ranking["candidate_table"]["winner_partitions_conserved"],
            "runtime_diagnostics_independently_recomputed": ranking["runtime_diagnostics"][
                "all_child_and_full_runtime_summaries_independently_recomputed"
            ],
            "all_child_method_contracts_exactly_compared_and_source_bound": (
                ranking["method_contract_verification"]
                ["all_children_exactly_matched_and_source_bound"]
            ),
            "profile_likelihood_recomputed_independently": (
                ranking["profile_likelihood_independent_recomputation"]["n_units"]
                == ranking["candidate_table"]["n_units"]
                and ranking["profile_likelihood_independent_recomputation"]["n_candidates"]
                == ranking["candidate_table"]["n_rows"]
            ),
            "profile_likelihood_certificates_match": ranking[
                "profile_likelihood_independent_recomputation"
            ]["all_profile_likelihoods_and_certificates_match"],
            "profile_relative_weights_match": ranking[
                "profile_likelihood_independent_recomputation"
            ]["all_relative_weights_match"],
            "profile_winner_tie_partitions_match": ranking[
                "profile_likelihood_independent_recomputation"
            ]["all_winner_tie_partitions_match"],
        },
    }
    if release_binding is not None:
        receipt["checks"].update({
            "release_contract_authenticated_and_eligible": True,
            "release_contract_all_snapshot_sources_rehashed": True,
            "extraction_and_ranking_bound_to_same_release": True,
            "full_runner_dependency_paths_and_shas_match_release": True,
            "frozen_scientific_and_scheduler_parameters_match": True,
            "post_input_identity_authenticated_and_exactly_equals_frozen_pre": (
                post_input_identity is not None
                and post_input_identity["exact_snapshot_equal"] is True
            ),
            "extraction_orchestration_session_batch_grant_completion_chain_verified": (
                extraction_public["orchestration"] is not None
                and extraction_public["orchestration"]["n_attested_children"] == 154
            ),
            "extraction_session_plus_11_resource_gates_authenticated": (
                extraction_public["orchestration"] is not None
                and extraction_public["orchestration"].get("n_authenticated_resource_gates") == 12
            ),
            "ranking_orchestration_session_batch_grant_completion_chain_verified": (
                ranking["orchestration"] is not None
                and ranking["orchestration"]["n_attested_children"] == 154
            ),
            "ranking_session_plus_11_resource_gates_authenticated": (
                ranking["orchestration"] is not None
                and ranking["orchestration"].get("n_authenticated_resource_gates") == 12
            ),
            "ranking_orchestration_bound_to_verified_extraction_session": (
                ranking["orchestration"] is not None
                and extraction_public["orchestration"] is not None
            ),
        })
        # Publication-boundary verification is intentionally fresh and last.
        # A long 154-task/table scan must not certify release/source/POST or
        # terminal bytes that drifted after the initial authentication.
        current_binding = load_release_contract_binding(
            release_contract_manifest,
            required_sources={"full_verifier": VERIFIER},
            _force_deep_reverification=True,
        )
        require(
            current_binding == release_binding,
            "release contract or a frozen source drifted before final publication",
        )
        current_post = verify_post_input_identity_receipt(
            post_input_identity_receipt, current_binding
        )
        require(
            current_post == post_input_identity,
            "POST input-identity evidence drifted before final publication",
        )
        _reauthenticate_orchestration_publication_boundary(
            extraction_public["orchestration"], "extraction"
        )
        _reauthenticate_orchestration_publication_boundary(
            ranking["orchestration"], "ranking"
        )
        receipt["checks"].update({
            "publication_boundary_release_and_all_frozen_sources_force_reverified": True,
            "publication_boundary_post_receipt_reauthenticated": True,
            "publication_boundary_extraction_session_and_terminal_reauthenticated": True,
            "publication_boundary_ranking_session_and_terminal_reauthenticated": True,
        })
    receipt["all_pass"] = bool(receipt["checks"]) and all(
        value is True for value in receipt["checks"].values()
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extraction-root", required=True, type=Path)
    parser.add_argument("--ranking-root", required=True, type=Path)
    parser.add_argument("--release-contract-manifest", required=True, type=Path)
    parser.add_argument("--post-input-identity-receipt", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    failure: str | None = None
    try:
        receipt = verify_bundle(
            args.extraction_root.resolve(),
            args.ranking_root.resolve(),
            release_contract_manifest=args.release_contract_manifest,
            post_input_identity_receipt=args.post_input_identity_receipt,
        )
    except (VerificationError, OSError, ValueError, KeyError, TypeError, csv.Error, gzip.BadGzipFile) as exc:
        failure = f"{type(exc).__name__}: {exc}"
        receipt = {
            "schema_name": "intersubmod.m2_full_independent_verification",
            "schema_version": "1.0.0",
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "scope": {"datasets": list(DATASETS), "chromosomes": list(AUTOSOMES), "expected_tasks": 154},
            "inputs": {
                "extraction_root": str(args.extraction_root.resolve()),
                "ranking_root": str(args.ranking_root.resolve()),
                "release_contract_manifest": str(args.release_contract_manifest.resolve()),
                "post_input_identity_receipt": str(args.post_input_identity_receipt.resolve()),
            },
            "failure": failure,
            "checks": {"verification_completed_without_contract_violation": False},
            "all_pass": False,
        }
    receipt["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": f"{args.output.name}.sha256",
        "covers": args.output.name,
    }
    write_immutable_receipt_exclusive(args.output, receipt)
    print(json.dumps({
        "output": str(args.output.resolve()),
        "sha256": sha256_path(args.output),
        "all_pass": receipt["all_pass"],
        "failure": failure,
    }, ensure_ascii=False))
    raise SystemExit(0 if receipt["all_pass"] else 1)


if __name__ == "__main__":
    main()
