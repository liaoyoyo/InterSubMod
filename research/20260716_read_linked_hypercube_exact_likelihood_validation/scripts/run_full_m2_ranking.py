#!/usr/bin/env python3
"""Resumable 7-dataset × chr1-22 M2 ranking orchestrator and aggregator.

The runner fails closed unless the upstream full extraction receipt covers all
154 tasks.  Every reusable chromosome result is revalidated against the exact
extraction receipt, ranker hash, parameters, output hashes, and external receipt
checksum.  It does not weaken a k>12 or optimizer abstention into a result.
"""

from __future__ import annotations

import argparse
from array import array
import concurrent.futures
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import secrets
import shlex
import signal
import stat
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
RUNNER = Path(__file__).resolve()
RANKER = SCRIPT_DIR / "build_m2_patterns_and_rank.py"
HYPERCUBE_SOLVER = SCRIPT_DIR / "hypercube_exact.py"
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
EXPECTED_TASKS = len(DATASETS) * len(AUTOSOMES)
ORCHESTRATION_POLICY = {
    "initial_batch_tasks": 8,
    "subsequent_batch_tasks": 16,
    "selection_policy": "canonical DATASETS order then chr1-22; first N pending",
    "unattested_or_orphan_child_policy": "FAIL_CLOSED",
}
ORCHESTRATION_SCHEMA_VERSION = "1.1.0"
LEGAL_CUMULATIVE_COUNTS = (8, 24, 40, 56, 72, 88, 104, 120, 136, 152, 154)
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
RELEASE_SCHEMA_NAME = "intersubmod.m2_release_run_manifest"
RELEASE_SCHEMA_VERSION = "1.0.0"
CANONICAL_MANIFEST_SHA256 = "16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45"
CANONICAL_SCHEMA_SHA256 = "47f54d86159beae66b68719e034c7755eba59b2d1719d4546153fa0ae74ff19f"
CANONICAL_SCHEMA_RELATIVE = (
    "docs/methodology/_assets/20260627_subclone_4axis_teaching/"
    "schemas/layered_input_manifest_v3.schema.json"
)
RELEASE_RESEARCH_RELATIVE = "research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts"
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
    "canonical_manifest_and_pre_receipt_copied_immutably", "publish_boundary_recheck_completed",
    "bootstrap20_and_seed20260716_frozen",
}
RELEASE_TOP_KEYS = {
    "schema_name", "schema_version", "task_type", "authority_mode", "validation_evidence_eligible",
    "created_at_utc", "scope", "entrypoints", "canonical_manifest",
    "pre_input_identity_receipt", "fresh_input_identity_verification", "source_snapshot",
    "canonical_schema", "runtime", "runtime_semantic_sha256", "parameters", "producer",
    "assurance", "checks", "all_pass", "receipt_integrity",
}
RELEASE_FILE_IDENTITY_KEYS = {
    "path", "st_dev", "st_ino", "st_nlink", "size_bytes", "mtime_ns", "ctime_ns",
    "mode_octal", "sha256",
}
FROZEN_RELEASE_PARAMETERS = {
    "extraction": {"mapq_min": 20, "baseq_min": 20, "bridge_thresholds": [1, 2, 3, 5], "workers": 2, "samtools_threads": 1},
    "ranking": {
        "thresholds": [1, 2, 3, 5], "component_bases": ["PS_HP1", "PS_HP2"],
        "hp_families": ["1", "2"], "structural_exact_pattern_minread_grid": [1, 2, 3, 5],
        "primary_structural_exact_pattern_minread": 3, "exact_k_max": 12,
        "max_vertex_sets": 256, "solver_time_limit_seconds_per_milp": 30.0,
        "fixed_error_grid": [0.005, 0.01, 0.02, 0.05], "minimum_bq_error_rate": 0.000001,
        "maximum_bq_error_rate": 0.25, "conditional_candidate_ranking_bootstrap_replicates": 20,
        "conditional_candidate_ranking_bootstrap_seed": 20260716,
        "tie_tolerance_log_likelihood": 0.000001, "workers": 2,
    },
    "scheduler": {"initial_batch_tasks": 8, "subsequent_batch_tasks": 16, "task_timeout_seconds": 28800, "timeout_grace_seconds": 300},
}
RELEASE_ASSURANCE = {
    "code_snapshot": "exact source bytes copied to regular non-aliased files; files 0444; directories 0555",
    "input_identity": "exact manifest-derived PRE plus mandatory fresh verifier snapshot equality",
    "immutable_inputs": "canonical manifest and supplied PRE receipt/sidecar copied read-only into contract",
    "publication": "same-filesystem renameat2(RENAME_NOREPLACE) after complete staging",
}
_DEEP_RELEASE_VERIFY_CACHE: dict[tuple[str, str], dict[str, Any]] = {}
TASK_STATUS_NAME = "runner_task_status.json"
BIOLOGICAL_IDS = {
    "COLO829": "COLO829", "H1437": "H1437", "H2009": "H2009",
    "HCC1395": "HCC1395", "HCC1395_DORADO": "HCC1395",
    "HCC1937": "HCC1937", "HCC1954": "HCC1954",
}
CANDIDATE_TABLE_COLUMNS = (
    "unit_key", "dataset", "chrom", "component_id", "threshold", "hp_family", "ps",
    "candidate_id", "vertex_set_id", "vertex_states", "vertex_roles", "parent_choice_count",
    "profile_log_likelihood", "relative_log_likelihood", "mixture_weights_pi", "winner_status",
    "tie_group", "coarse_topology_class", "candidate_set_complete",
)

# Deliberately duplicated rather than imported from the chromosome ranker: the
# orchestrator must reject a child whose declared method drifts even if the
# producer changed its own constant at the same time.
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

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_m2_patterns_and_rank import (  # noqa: E402
    RUNTIME_DIAGNOSTIC_FIELDS,
    RUNTIME_METRICS,
    sha256_path,
    summarize_runtime_values,
)


def verify_sha256_sidecar(path: Path, sidecar_name: str | None = None) -> bool:
    sidecar = path.parent / (sidecar_name or f"{path.name}.sha256")
    if not sidecar.is_file():
        return False
    fields = sidecar.read_text(encoding="ascii", errors="strict").strip().split()
    return len(fields) == 2 and fields[0] == sha256_path(path) and fields[1] == path.name


def write_json_and_sha256_exclusive(path: Path, payload: dict[str, Any]) -> Path:
    """Create an immutable checkpoint/receipt and checksum without overwriting."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    checksum_path = path.with_name(f"{path.name}.sha256")
    with checksum_path.open("x", encoding="ascii") as handle:
        handle.write(f"{sha256_path(path)}  {path.name}\n")
    return checksum_path


def semantic_json_sha256(payload: Any) -> str:
    import hashlib

    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    ).hexdigest()


SESSION_KEYS = {
    "schema_name", "schema_version", "stage", "session_id", "session_nonce",
    "created_at_utc", "created_monotonic_ns", "host_boot_id", "release_manifest",
    "release_binding_semantic_sha256", "run_contract_semantic_sha256", "scope",
    "output_root", "producer_sources", "scheduler_policy", "parent_extraction",
    "resource_gate", "receipt_integrity",
}
BATCH_START_KEYS = {
    "schema_name", "schema_version", "stage", "session_id", "session_sha256",
    "batch_index", "batch_id", "batch_nonce", "previous_chain_head", "before_count",
    "max_new_tasks", "effective_count", "selected_task_ids",
    "run_contract_semantic_sha256", "runner_source", "resource_gate", "created_at_utc",
    "created_monotonic_ns", "host_boot_id", "receipt_integrity",
}
GRANT_KEYS = {
    "schema_name", "schema_version", "stage", "session_id", "session_sha256",
    "batch_id", "batch_start_sha256", "task_id", "dataset", "chrom", "task_ordinal",
    "child_nonce", "command", "command_semantic_sha256", "producer_sources",
    "input_identity", "parameters_semantic_sha256", "expected_output_dir",
    "issued_at_utc", "issued_monotonic_ns", "host_boot_id", "receipt_integrity",
}
COMPLETION_KEYS = {
    "schema_name", "schema_version", "stage", "session_id", "session_sha256",
    "batch_id", "grant_identity", "task_id", "dataset", "chrom", "task_ordinal",
    "child_receipt_identity", "child_outputs_semantic_sha256", "command_semantic_sha256",
    "status", "returncode", "timed_out", "process_group_id", "started_monotonic_ns",
    "completed_at_utc", "completed_monotonic_ns", "host_boot_id", "runner_source",
    "receipt_integrity",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def host_boot_id() -> str:
    try:
        return Path("/proc/sys/kernel/random/boot_id").read_text(encoding="ascii").strip()
    except OSError:
        return "unavailable"


def assert_no_symlink_ancestors(path: Path) -> None:
    absolute = path.absolute()
    current = Path(absolute.anchor)
    for part in absolute.parts[1:]:
        current = current / part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        if stat.S_ISLNK(observed.st_mode):
            raise RuntimeError(f"release path contains a symlink component: {current}")
        if current != absolute and not stat.S_ISDIR(observed.st_mode):
            raise RuntimeError(f"release path ancestor is not a directory: {current}")


def _receipt_integrity(path: Path) -> dict[str, str]:
    return {"scheme": "external_sha256_sidecar_v1", "sidecar_name": f"{path.name}.sha256", "covers": path.name}


def write_immutable_json_exclusive(path: Path, payload: Mapping[str, Any]) -> dict[str, str]:
    assert_no_symlink_ancestors(path.parent)
    write_json_and_sha256_exclusive(path, dict(payload))
    sidecar = path.with_name(f"{path.name}.sha256")
    path.chmod(0o444)
    sidecar.chmod(0o444)
    return {"path": str(path.resolve()), "sha256": sha256_path(path)}


def load_immutable_json(path: Path, label: str) -> dict[str, Any]:
    raw = path.absolute()
    assert_no_symlink_ancestors(raw)
    try:
        observed, side = os.lstat(raw), os.lstat(raw.with_name(f"{raw.name}.sha256"))
    except OSError as exc:
        raise RuntimeError(f"{label} is missing: {exc}") from exc
    if (
        not stat.S_ISREG(observed.st_mode) or stat.S_ISLNK(observed.st_mode)
        or stat.S_IMODE(observed.st_mode) != 0o444 or observed.st_nlink != 1
        or not stat.S_ISREG(side.st_mode) or stat.S_ISLNK(side.st_mode)
        or stat.S_IMODE(side.st_mode) != 0o444 or side.st_nlink != 1
        or not verify_sha256_sidecar(raw)
    ):
        raise RuntimeError(f"{label} immutable file/sidecar contract failed")
    try:
        value = json.loads(raw.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"{label} unreadable: {exc}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} is not a JSON object")
    if "receipt_integrity" in value and value["receipt_integrity"] != _receipt_integrity(raw):
        raise RuntimeError(f"{label} receipt_integrity filename binding mismatch")
    return value


def _resource_gate_path(output_root: Path, gate_scope: str, batch_index: int | None = None) -> Path:
    if gate_scope == "session":
        name = "session.json"
    elif gate_scope == "batch" and isinstance(batch_index, int) and batch_index > 0:
        name = f"batch_{batch_index:03d}.json"
    else:
        raise ValueError("resource gate scope/index is invalid")
    return output_root / "_orchestration" / "resource_gates" / name


def _resource_gate_identity(path: Path, payload: Mapping[str, Any]) -> dict[str, Any]:
    sidecar = path.with_name(f"{path.name}.sha256")
    return {
        "path": str(path.resolve()),
        "sha256": sha256_path(path),
        "semantic_sha256": semantic_json_sha256(payload),
        "gate_id": payload["gate_id"],
        "sidecar_path": str(sidecar.resolve()),
        "sidecar_sha256": sha256_path(sidecar),
    }


def resource_gate_preview(output_root: Path, conflicts: Mapping[str, Any]) -> dict[str, Any]:
    """Return a live, non-persisted process/disk decision for CLI preflight output."""
    probe = output_root if output_root.exists() else output_root.parent
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    probe = probe.resolve(strict=True)
    observed = os.lstat(probe)
    filesystem = os.statvfs(probe)
    available_bytes = int(filesystem.f_bavail) * int(filesystem.f_frsize)
    return {
        "process_count": int(conflicts["process_count"]),
        "root_count": int(conflicts["root_count"]),
        "representatives": list(conflicts["representatives"]),
        "filesystem_probe_path": str(probe),
        "filesystem_st_dev": int(observed.st_dev),
        "filesystem_f_bavail": int(filesystem.f_bavail),
        "filesystem_f_frsize": int(filesystem.f_frsize),
        "available_bytes": available_bytes,
        "required_reserve_bytes": RESOURCE_GATE_REQUIRED_RESERVE_BYTES,
        "zero_conflict_pass": int(conflicts["process_count"]) == 0,
        "disk_pass": available_bytes >= RESOURCE_GATE_REQUIRED_RESERVE_BYTES,
    }


def create_resource_gate_receipt(
    output_root: Path,
    *,
    stage: str,
    gate_scope: str,
    target: Mapping[str, Any],
    producer_source: Mapping[str, Any],
    conflicts: Mapping[str, Any] | None = None,
    batch_index: int | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Persist a one-use, source-bound zero-conflict + 300 GiB gate receipt."""
    if stage not in {"extraction", "ranking"}:
        raise ValueError("resource gate stage is invalid")
    root = output_root.resolve(strict=True)
    root_stat = os.lstat(root)
    if not stat.S_ISDIR(root_stat.st_mode) or stat.S_ISLNK(root_stat.st_mode):
        raise RuntimeError("resource gate output root is not a real directory")
    conflicts = active_conflicts() if conflicts is None else dict(conflicts)
    process_core = {
        "process_count": int(conflicts["process_count"]),
        "root_count": int(conflicts["root_count"]),
        "representatives": list(conflicts["representatives"]),
    }
    process_snapshot = {**process_core, "semantic_sha256": semantic_json_sha256(process_core)}
    filesystem = os.statvfs(root)
    available_bytes = int(filesystem.f_bavail) * int(filesystem.f_frsize)
    filesystem_core = {
        "probe_path": str(root),
        "target_output_root": str(root),
        "st_dev": int(root_stat.st_dev),
        "f_bavail": int(filesystem.f_bavail),
        "f_frsize": int(filesystem.f_frsize),
        "available_bytes": available_bytes,
        "required_reserve_bytes": RESOURCE_GATE_REQUIRED_RESERVE_BYTES,
        "disk_pass": available_bytes >= RESOURCE_GATE_REQUIRED_RESERVE_BYTES,
    }
    filesystem_snapshot = {
        **filesystem_core, "semantic_sha256": semantic_json_sha256(filesystem_core),
    }
    checks = {
        "zero_conflict": process_core
        == {"process_count": 0, "root_count": 0, "representatives": []},
        "disk_reserve": filesystem_core["disk_pass"],
        "all_pass": False,
    }
    checks["all_pass"] = checks["zero_conflict"] and checks["disk_reserve"]
    if not checks["all_pass"]:
        raise RuntimeError(
            "resource gate failed: "
            f"process_count={process_core['process_count']} "
            f"available_bytes={available_bytes} "
            f"required_reserve_bytes={RESOURCE_GATE_REQUIRED_RESERVE_BYTES}"
        )
    gate_path = _resource_gate_path(output_root, gate_scope, batch_index)
    payload: dict[str, Any] = {
        "schema_name": RESOURCE_GATE_SCHEMA_NAME,
        "schema_version": RESOURCE_GATE_SCHEMA_VERSION,
        "stage": stage,
        "gate_scope": gate_scope,
        "gate_id": "",
        "gate_nonce": secrets.token_hex(32),
        "target": dict(target),
        "process_snapshot": process_snapshot,
        "filesystem_snapshot": filesystem_snapshot,
        "producer_source": dict(producer_source),
        "observed_at_utc": utc_now(),
        "observed_monotonic_ns": time.monotonic_ns(),
        "host_boot_id": host_boot_id(),
        "checks": checks,
        "receipt_integrity": _receipt_integrity(gate_path),
    }
    payload["gate_id"] = semantic_json_sha256({
        key: value for key, value in payload.items()
        if key not in {"gate_id", "receipt_integrity"}
    })
    write_immutable_json_exclusive(gate_path, payload)
    return payload, _resource_gate_identity(gate_path, payload)


def load_resource_gate_receipt(
    identity: Mapping[str, Any],
    *,
    expected_path: Path,
    expected_stage: str,
    expected_scope: str,
    expected_target: Mapping[str, Any],
    expected_producer_source: Mapping[str, Any],
) -> dict[str, Any]:
    """Runner-side fail-closed validation; the final verifier reimplements this independently."""
    if set(identity) != RESOURCE_GATE_IDENTITY_KEYS:
        raise RuntimeError("resource gate identity exact-key schema mismatch")
    path = Path(str(identity.get("path", "")))
    if path.resolve(strict=True) != expected_path.resolve(strict=True):
        raise RuntimeError("resource gate identity path swap detected")
    payload = load_immutable_json(path, "resource gate receipt")
    if set(payload) != RESOURCE_GATE_RECEIPT_KEYS:
        raise RuntimeError("resource gate receipt exact-key schema mismatch")
    if _resource_gate_identity(path, payload) != dict(identity):
        raise RuntimeError("resource gate identity/sidecar mismatch")
    expected_gate_id = semantic_json_sha256({
        key: value for key, value in payload.items()
        if key not in {"gate_id", "receipt_integrity"}
    })
    process = payload.get("process_snapshot") or {}
    process_core = {
        "process_count": process.get("process_count"),
        "root_count": process.get("root_count"),
        "representatives": process.get("representatives"),
    }
    filesystem = payload.get("filesystem_snapshot") or {}
    filesystem_core = {
        key: filesystem.get(key)
        for key in (
            "probe_path", "target_output_root", "st_dev", "f_bavail", "f_frsize",
            "available_bytes", "required_reserve_bytes", "disk_pass",
        )
    }
    root_stat = os.lstat(expected_path.parents[2])
    if (
        payload.get("schema_name") != RESOURCE_GATE_SCHEMA_NAME
        or payload.get("schema_version") != RESOURCE_GATE_SCHEMA_VERSION
        or payload.get("stage") != expected_stage
        or payload.get("gate_scope") != expected_scope
        or payload.get("target") != dict(expected_target)
        or payload.get("producer_source") != dict(expected_producer_source)
        or payload.get("gate_id") != expected_gate_id
        or process_core != {"process_count": 0, "root_count": 0, "representatives": []}
        or process.get("semantic_sha256") != semantic_json_sha256(process_core)
        or not isinstance(filesystem_core["f_bavail"], int)
        or not isinstance(filesystem_core["f_frsize"], int)
        or filesystem_core["available_bytes"]
        != filesystem_core["f_bavail"] * filesystem_core["f_frsize"]
        or filesystem_core["required_reserve_bytes"] != RESOURCE_GATE_REQUIRED_RESERVE_BYTES
        or filesystem_core["disk_pass"] is not True
        or filesystem_core["available_bytes"] < RESOURCE_GATE_REQUIRED_RESERVE_BYTES
        or filesystem.get("semantic_sha256") != semantic_json_sha256(filesystem_core)
        or filesystem_core["probe_path"] != str(expected_path.parents[2].resolve())
        or filesystem_core["target_output_root"] != str(expected_path.parents[2].resolve())
        or filesystem_core["st_dev"] != int(root_stat.st_dev)
        or payload.get("checks")
        != {"zero_conflict": True, "disk_reserve": True, "all_pass": True}
    ):
        raise RuntimeError("resource gate semantic/process/disk/source binding mismatch")
    return payload


def release_producer_sources(binding: Mapping[str, Any]) -> dict[str, Any]:
    sources = binding["snapshot_sources"]
    return {
        "runner": dict(sources["full_ranking_runner"]),
        "child_producer": dict(sources["ranker"]),
        "dependencies": {"hypercube_solver": dict(sources["hypercube_solver"])},
    }


def extraction_parent_identity(extraction_root: Path, terminal_path: Path) -> dict[str, Any]:
    session_path = extraction_root / "_orchestration" / "session.json"
    session = load_immutable_json(session_path, "parent extraction session")
    terminal = load_immutable_json(terminal_path, "parent extraction terminal receipt")
    session_identity = (terminal.get("orchestration") or {}).get("session_identity") or {}
    if (
        session.get("schema_name") != "intersubmod.m2_orchestration_session"
        or session.get("stage") != "extraction"
        or session_identity != {
            "path": str(session_path.resolve()), "sha256": sha256_path(session_path),
            "session_id": session.get("session_id"),
        }
    ):
        raise RuntimeError("ranking parent extraction session/terminal mismatch")
    return {
        "session_id": session["session_id"],
        "terminal_receipt_path": str(terminal_path.resolve()),
        "terminal_receipt_sha256": sha256_path(terminal_path),
    }


def _session_identity(output_root: Path, session: Mapping[str, Any]) -> dict[str, Any]:
    path = output_root / "_orchestration" / "session.json"
    return {"path": str(path.resolve()), "sha256": sha256_path(path), "session_id": session["session_id"]}


def ensure_release_orchestration_session(
    output_root: Path,
    binding: Mapping[str, Any],
    run_contract: Mapping[str, Any],
    parent_extraction: Mapping[str, Any],
    gate: Mapping[str, Any] | None,
) -> dict[str, Any]:
    output_root = output_root.absolute()
    assert_no_symlink_ancestors(output_root)
    session_path = output_root / "_orchestration" / "session.json"
    if session_path.exists():
        if output_root.is_symlink() or not output_root.is_dir():
            raise RuntimeError("release ranking output root is not a real directory")
        session = load_immutable_json(session_path, "ranking orchestration session")
    else:
        if not output_root.exists():
            output_root.parent.mkdir(parents=True, exist_ok=True)
            output_root.mkdir()
        if output_root.is_symlink() or not output_root.is_dir():
            raise RuntimeError("release ranking output root is not a real directory")
        if gate is None:
            raise RuntimeError("new ranking session requires an authenticated resource gate")
        root_stat = os.lstat(output_root)
        if not stat.S_ISDIR(root_stat.st_mode) or root_stat.st_nlink < 2:
            raise RuntimeError("new ranking release output root identity is invalid")
        session = {
            "schema_name": "intersubmod.m2_orchestration_session",
            "schema_version": ORCHESTRATION_SCHEMA_VERSION, "stage": "ranking", "session_id": "",
            "session_nonce": secrets.token_hex(32), "created_at_utc": utc_now(),
            "created_monotonic_ns": time.monotonic_ns(), "host_boot_id": host_boot_id(),
            "release_manifest": dict(binding["release_manifest"]),
            "release_binding_semantic_sha256": semantic_json_sha256(binding),
            "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
            "scope": {"datasets": list(DATASETS), "chromosomes": list(AUTOSOMES), "expected_tasks": EXPECTED_TASKS},
            "output_root": {"path": str(output_root.resolve()), "st_dev": int(root_stat.st_dev), "st_ino": int(root_stat.st_ino)},
            "producer_sources": release_producer_sources(binding),
            "scheduler_policy": dict(ORCHESTRATION_POLICY),
            "parent_extraction": dict(parent_extraction), "resource_gate": dict(gate),
            "receipt_integrity": _receipt_integrity(session_path),
        }
        session["session_id"] = semantic_json_sha256({
            key: value for key, value in session.items()
            if key not in {"session_id", "receipt_integrity"}
        })
        write_immutable_json_exclusive(session_path, session)
    if set(session) != SESSION_KEYS:
        raise RuntimeError("ranking orchestration session exact-key schema mismatch")
    root_stat = os.lstat(output_root)
    if not stat.S_ISDIR(root_stat.st_mode) or root_stat.st_nlink < 2:
        raise RuntimeError("ranking release output root identity is invalid")
    static = {
        "schema_name": "intersubmod.m2_orchestration_session",
        "schema_version": ORCHESTRATION_SCHEMA_VERSION,
        "stage": "ranking", "release_manifest": dict(binding["release_manifest"]),
        "release_binding_semantic_sha256": semantic_json_sha256(binding),
        "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
        "scope": {"datasets": list(DATASETS), "chromosomes": list(AUTOSOMES), "expected_tasks": EXPECTED_TASKS},
        "output_root": {"path": str(output_root.resolve()), "st_dev": int(root_stat.st_dev), "st_ino": int(root_stat.st_ino)},
        "producer_sources": release_producer_sources(binding),
        "scheduler_policy": dict(ORCHESTRATION_POLICY), "parent_extraction": dict(parent_extraction),
    }
    if any(session.get(key) != value for key, value in static.items()):
        raise RuntimeError("ranking orchestration session binding/root/parent mismatch")
    session_target = {
        "output_root": static["output_root"],
        "release_manifest": dict(binding["release_manifest"]),
        "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
        "parent_extraction": dict(parent_extraction),
    }
    gate_payload = load_resource_gate_receipt(
        session.get("resource_gate") or {},
        expected_path=_resource_gate_path(output_root, "session"),
        expected_stage="ranking",
        expected_scope="session",
        expected_target=session_target,
        expected_producer_source=session["producer_sources"]["runner"],
    )
    if gate is not None and dict(gate) != session.get("resource_gate"):
        raise RuntimeError("new ranking session gate differs from persisted session")
    if (
        gate_payload.get("host_boot_id") != session.get("host_boot_id")
        or gate_payload.get("observed_monotonic_ns", -1) > session.get("created_monotonic_ns", -1)
    ):
        raise RuntimeError("ranking session is not temporally bound to its resource gate")
    expected_id = semantic_json_sha256({key: value for key, value in session.items() if key not in {"session_id", "receipt_integrity"}})
    if session.get("session_id") != expected_id:
        raise RuntimeError("ranking orchestration session semantic identity mismatch")
    return session


def _task_id(dataset: str, chrom: str) -> str:
    return f"{dataset}:{chrom}"


def _task_ordinal(dataset: str, chrom: str) -> int:
    return DATASETS.index(dataset) * len(AUTOSOMES) + AUTOSOMES.index(chrom) + 1


def canonical_sort_results(results: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        results,
        key=lambda item: (DATASETS.index(item["dataset"]), AUTOSOMES.index(item["chrom"])),
    )


def create_batch_start_and_grants(
    output_root: Path,
    session: Mapping[str, Any],
    run_contract: Mapping[str, Any],
    selected_specs: Sequence[tuple[Any, ...]],
    *, before_count: int, previous_chain_head: Mapping[str, str] | None,
    batch_index: int, max_new_tasks: int, gate: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    effective = len(selected_specs)
    expected_max = 8 if before_count == 0 else 16
    if (
        max_new_tasks != expected_max or effective != min(expected_max, EXPECTED_TASKS - before_count)
        or before_count not in (0, *LEGAL_CUMULATIVE_COUNTS[:-1])
    ):
        raise RuntimeError("illegal ranking release batch size/count chain")
    selected_ids = [_task_id(spec[0], spec[1]) for spec in selected_specs]
    session_path = output_root / "_orchestration" / "session.json"
    path = output_root / "_orchestration" / "batches" / f"batch_{batch_index:03d}_start.json"
    gate_target = {
        "output_root": dict(session["output_root"]),
        "session_id": session["session_id"],
        "session_sha256": sha256_path(session_path),
        "batch_index": batch_index,
        "before_count": before_count,
        "max_new_tasks": max_new_tasks,
        "effective_count": effective,
        "selected_task_ids": selected_ids,
        "previous_chain_head": None if previous_chain_head is None else dict(previous_chain_head),
    }
    gate_payload = load_resource_gate_receipt(
        gate,
        expected_path=_resource_gate_path(output_root, "batch", batch_index),
        expected_stage="ranking",
        expected_scope="batch",
        expected_target=gate_target,
        expected_producer_source=session["producer_sources"]["runner"],
    )
    batch: dict[str, Any] = {
        "schema_name": "intersubmod.m2_orchestration_batch_start",
        "schema_version": ORCHESTRATION_SCHEMA_VERSION,
        "stage": "ranking", "session_id": session["session_id"],
        "session_sha256": sha256_path(session_path), "batch_index": batch_index,
        "batch_id": "", "batch_nonce": secrets.token_hex(32),
        "previous_chain_head": None if previous_chain_head is None else dict(previous_chain_head),
        "before_count": before_count, "max_new_tasks": max_new_tasks,
        "effective_count": effective, "selected_task_ids": selected_ids,
        "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
        "runner_source": dict(session["producer_sources"]["runner"]),
        "resource_gate": dict(gate), "created_at_utc": utc_now(),
        "created_monotonic_ns": time.monotonic_ns(), "host_boot_id": host_boot_id(),
        "receipt_integrity": _receipt_integrity(path),
    }
    batch["batch_id"] = semantic_json_sha256({key: value for key, value in batch.items() if key not in {"batch_id", "receipt_integrity"}})
    if (
        gate_payload.get("host_boot_id") != batch["host_boot_id"]
        or gate_payload.get("observed_monotonic_ns", -1) > batch["created_monotonic_ns"]
    ):
        raise RuntimeError("ranking batch start is not temporally bound to its resource gate")
    write_immutable_json_exclusive(path, batch)
    batch_sha = sha256_path(path)
    grants: dict[str, dict[str, Any]] = {}
    for spec in selected_specs:
        dataset, chrom, extraction_dir, output_dir, command = spec[:5]
        extraction_receipt_path = extraction_dir / "receipt.json"
        extraction_receipt = json.loads(extraction_receipt_path.read_text(encoding="utf-8"))
        task_id = _task_id(dataset, chrom)
        grant_path = output_root / "_orchestration" / "tasks" / dataset / chrom / "grant.json"
        grant = {
            "schema_name": "intersubmod.m2_orchestration_child_grant",
            "schema_version": ORCHESTRATION_SCHEMA_VERSION,
            "stage": "ranking", "session_id": session["session_id"],
            "session_sha256": sha256_path(session_path), "batch_id": batch["batch_id"],
            "batch_start_sha256": batch_sha, "task_id": task_id, "dataset": dataset,
            "chrom": chrom, "task_ordinal": _task_ordinal(dataset, chrom),
            "child_nonce": secrets.token_hex(32), "command": list(command),
            "command_semantic_sha256": semantic_json_sha256(list(command)),
            "producer_sources": dict(session["producer_sources"]),
            "input_identity": {
                "extraction_child_receipt_path": str(extraction_receipt_path.resolve()),
                "extraction_child_receipt_sha256": sha256_path(extraction_receipt_path),
                "extraction_child_receipt_semantic_sha256": semantic_json_sha256(extraction_receipt),
                "extraction_outputs_semantic_sha256": semantic_json_sha256(extraction_receipt.get("outputs") or {}),
            },
            "parameters_semantic_sha256": semantic_json_sha256(spec[7]),
            "expected_output_dir": str(output_dir.resolve()), "issued_at_utc": utc_now(),
            "issued_monotonic_ns": time.monotonic_ns(), "host_boot_id": host_boot_id(),
            "receipt_integrity": _receipt_integrity(grant_path),
        }
        write_immutable_json_exclusive(grant_path, grant)
        grants[task_id] = {"path": str(grant_path.resolve()), "sha256": sha256_path(grant_path), "document": grant}
    return batch | {"path": str(path.resolve()), "sha256": batch_sha}, grants


def write_child_completion(
    output_root: Path, session: Mapping[str, Any], batch: Mapping[str, Any],
    grant_record: Mapping[str, Any], result: Mapping[str, Any],
) -> dict[str, str]:
    dataset, chrom = str(result["dataset"]), str(result["chrom"])
    task_id = _task_id(dataset, chrom)
    if result.get("status") != "PASS" or int(result.get("returncode", -1)) != 0:
        raise RuntimeError(f"cannot attest nonpassing ranking child: {task_id}")
    receipt_path = Path(str(result["receipt_path"]))
    if not verify_sha256_sidecar(receipt_path):
        raise RuntimeError(f"rank child receipt sidecar failed: {task_id}")
    receipt_path.chmod(0o444)
    receipt_path.with_name(f"{receipt_path.name}.sha256").chmod(0o444)
    receipt = load_immutable_json(receipt_path, "ranking child receipt")
    grant_path, grant = Path(str(grant_record["path"])), grant_record["document"]
    path = output_root / "_orchestration" / "tasks" / dataset / chrom / "completion.json"
    completion = {
        "schema_name": "intersubmod.m2_orchestration_child_completion",
        "schema_version": ORCHESTRATION_SCHEMA_VERSION,
        "stage": "ranking", "session_id": session["session_id"],
        "session_sha256": sha256_path(output_root / "_orchestration" / "session.json"),
        "batch_id": batch["batch_id"],
        "grant_identity": {"path": str(grant_path.resolve()), "sha256": sha256_path(grant_path), "semantic_sha256": semantic_json_sha256(grant)},
        "task_id": task_id, "dataset": dataset, "chrom": chrom,
        "task_ordinal": _task_ordinal(dataset, chrom),
        "child_receipt_identity": {"path": str(receipt_path.resolve()), "size_bytes": receipt_path.stat().st_size, "sha256": sha256_path(receipt_path), "semantic_sha256": semantic_json_sha256(receipt)},
        "child_outputs_semantic_sha256": semantic_json_sha256(receipt.get("outputs") or {}),
        "command_semantic_sha256": grant["command_semantic_sha256"], "status": "PASS",
        "returncode": 0, "timed_out": False, "process_group_id": result.get("process_group_id"),
        "started_monotonic_ns": result.get("started_monotonic_ns"), "completed_at_utc": utc_now(),
        "completed_monotonic_ns": time.monotonic_ns(), "host_boot_id": host_boot_id(),
        "runner_source": dict(session["producer_sources"]["runner"]),
        "receipt_integrity": _receipt_integrity(path),
    }
    return write_immutable_json_exclusive(path, completion)


def validate_attested_completion(
    path: Path, session: Mapping[str, Any], *, expected_task_id: str | None = None,
    expected_batch: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, str]]:
    completion = load_immutable_json(path, "ranking child completion")
    if set(completion) != COMPLETION_KEYS:
        raise RuntimeError("ranking child completion exact-key schema mismatch")
    task_id = completion.get("task_id")
    dataset, chrom = str(completion.get("dataset")), str(completion.get("chrom"))
    session_root = Path(str((session.get("output_root") or {}).get("path", "")))
    session_path = session_root / "_orchestration" / "session.json"
    if (
        completion.get("schema_name") != "intersubmod.m2_orchestration_child_completion"
        or completion.get("schema_version") != ORCHESTRATION_SCHEMA_VERSION
        or completion.get("stage") != "ranking"
        or completion.get("session_id") != session.get("session_id")
        or (expected_task_id is not None and task_id != expected_task_id)
        or completion.get("status") != "PASS" or completion.get("returncode") != 0
        or completion.get("timed_out") is not False
        or completion.get("session_sha256") != sha256_path(session_path)
        or task_id != _task_id(dataset, chrom)
        or completion.get("task_ordinal") != _task_ordinal(dataset, chrom)
        or completion.get("runner_source") != session["producer_sources"]["runner"]
        or (expected_batch is not None and completion.get("host_boot_id") != expected_batch.get("host_boot_id"))
        or path.resolve() != (
            session_root / "_orchestration" / "tasks" / dataset / chrom / "completion.json"
        ).resolve()
        or (expected_batch is not None and completion.get("batch_id") != expected_batch.get("batch_id"))
    ):
        raise RuntimeError("ranking completion session/task/status mismatch")
    grant_identity = completion.get("grant_identity") or {}
    grant_path = Path(str(grant_identity.get("path", "")))
    grant = load_immutable_json(grant_path, "ranking child grant")
    if (
        set(grant) != GRANT_KEYS
        or grant.get("schema_name") != "intersubmod.m2_orchestration_child_grant"
        or grant.get("schema_version") != ORCHESTRATION_SCHEMA_VERSION
        or grant.get("stage") != "ranking"
        or grant.get("task_id") != task_id
        or grant.get("session_id") != session.get("session_id")
        or sha256_path(grant_path) != grant_identity.get("sha256")
        or semantic_json_sha256(grant) != grant_identity.get("semantic_sha256")
        or grant.get("command_semantic_sha256") != completion.get("command_semantic_sha256")
        or grant.get("session_sha256") != sha256_path(session_path)
        or grant.get("dataset") != dataset or grant.get("chrom") != chrom
        or grant.get("task_ordinal") != completion.get("task_ordinal")
        or grant.get("producer_sources") != session.get("producer_sources")
        or grant.get("host_boot_id") != completion.get("host_boot_id")
        or grant_path.resolve() != (
            session_root / "_orchestration" / "tasks" / dataset / chrom / "grant.json"
        ).resolve()
        or Path(str(grant.get("expected_output_dir", ""))).resolve() != (
            session_root / "samples" / dataset / chrom
        ).resolve()
        or semantic_json_sha256(grant.get("command")) != grant.get("command_semantic_sha256")
        or (expected_batch is not None and (
            grant.get("batch_id") != expected_batch.get("batch_id")
            or grant.get("batch_start_sha256") != expected_batch.get("sha256")
            or not (
                expected_batch.get("created_monotonic_ns", -1)
                <= grant.get("issued_monotonic_ns", -1)
                <= completion.get("started_monotonic_ns", -1)
                <= completion.get("completed_monotonic_ns", -1)
            )
        ))
    ):
        raise RuntimeError("ranking grant/completion binding mismatch")
    child = completion.get("child_receipt_identity") or {}
    receipt_path = Path(str(child.get("path", "")))
    receipt = load_immutable_json(receipt_path, "attested ranking child receipt")
    if (
        (receipt_path.parent / TASK_STATUS_NAME).exists()
        or sha256_path(receipt_path) != child.get("sha256")
        or receipt_path.stat().st_size != child.get("size_bytes")
        or semantic_json_sha256(receipt) != child.get("semantic_sha256")
        or semantic_json_sha256(receipt.get("outputs") or {}) != completion.get("child_outputs_semantic_sha256")
    ):
        raise RuntimeError("attested ranking child receipt/output mismatch")
    for output_identity in (receipt.get("outputs") or {}).values():
        output_path = Path(str(output_identity.get("path", "")))
        output_stat = os.lstat(output_path)
        if (
            output_path.resolve().parent != receipt_path.parent.resolve()
            or not stat.S_ISREG(output_stat.st_mode) or stat.S_ISLNK(output_stat.st_mode)
            or output_stat.st_nlink != 1
            or output_path.stat().st_size != output_identity.get("size_bytes")
            or sha256_path(output_path) != output_identity.get("sha256")
        ):
            raise RuntimeError("attested ranking declared output identity mismatch")
    upstream = ((receipt.get("provenance") or {}).get("upstream_extraction_receipt") or {})
    extraction_path = Path(str(upstream.get("path", "")))
    extraction_receipt = json.loads(extraction_path.read_text(encoding="utf-8"))
    expected_input = {
        "extraction_child_receipt_path": str(extraction_path.resolve()),
        "extraction_child_receipt_sha256": sha256_path(extraction_path),
        "extraction_child_receipt_semantic_sha256": semantic_json_sha256(extraction_receipt),
        "extraction_outputs_semantic_sha256": semantic_json_sha256(extraction_receipt.get("outputs") or {}),
    }
    if grant.get("input_identity") != expected_input:
        raise RuntimeError("ranking grant upstream extraction identity mismatch")
    for identity in (extraction_receipt.get("outputs") or {}).values():
        source = Path(str(identity.get("path", "")))
        source_stat = os.lstat(source)
        if (
            not stat.S_ISREG(source_stat.st_mode) or stat.S_ISLNK(source_stat.st_mode)
            or source_stat.st_nlink != 1 or source.stat().st_size != identity.get("size_bytes")
            or sha256_path(source) != identity.get("sha256")
        ):
            raise RuntimeError("ranking grant upstream extraction output drift")
    expected_parameters = {
        key: value for key, value in (receipt.get("parameters") or {}).items()
        if key in {
            "structural_exact_pattern_minread_grid", "primary_structural_exact_pattern_minread",
            "scoring_minread", "exact_k_max", "max_vertex_sets",
            "solver_time_limit_seconds_per_milp", "minimum_bq_error_rate",
            "maximum_bq_error_rate", "fixed_error_grid_conditional_binary_flip_probability",
            "conditional_candidate_ranking_bootstrap_replicates",
            "conditional_candidate_ranking_bootstrap_base_seed", "tie_tolerance_log_likelihood",
        }
    }
    if grant.get("parameters_semantic_sha256") != semantic_json_sha256(expected_parameters):
        raise RuntimeError("ranking grant parameter identity differs from child receipt")
    return completion, {"path": str(path.resolve()), "sha256": sha256_path(path)}


def load_release_orchestration_state(output_root: Path, session: Mapping[str, Any]) -> dict[str, Any]:
    checkpoint_dir = output_root / "checkpoints"
    paths = list(checkpoint_dir.glob("checkpoint_*_of_154.json")) if checkpoint_dir.exists() else []
    terminal = output_root / "full_ranking_receipt.json"
    if terminal.exists():
        paths.append(terminal)
    rows = []
    for path in paths:
        receipt = load_immutable_json(path, "ranking chain receipt")
        count = (receipt.get("orchestration") or {}).get("cumulative_attested_tasks")
        if not isinstance(count, int):
            raise RuntimeError("ranking chain receipt lacks attested count")
        expected_path = (
            terminal
            if count == EXPECTED_TASKS
            else checkpoint_dir / f"checkpoint_{count:03d}_of_154.json"
        )
        if (
            path.resolve() != expected_path.resolve()
            or receipt.get("receipt_integrity") != _receipt_integrity(path)
        ):
            raise RuntimeError("ranking chain receipt exact path/integrity mismatch")
        if (count == EXPECTED_TASKS) != (path == terminal):
            raise RuntimeError("ranking terminal orchestration count/path mismatch")
        rows.append((count, path, receipt))
    rows.sort(key=lambda row: row[0])
    if [row[0] for row in rows] != list(LEGAL_CUMULATIVE_COUNTS[:len(rows)]):
        raise RuntimeError("ranking checkpoint count chain gap/reorder")
    previous = None
    completions: dict[str, dict[str, str]] = {}
    child_receipts: dict[str, dict[str, Any]] = {}
    child_receipt_paths: dict[str, str] = {}
    for index, (count, path, receipt) in enumerate(rows, start=1):
        orch = receipt.get("orchestration") or {}
        if set(orch) != {"session_identity", "batch_start_identity", "previous_chain_head", "batch_completion_attestations", "cumulative_attested_tasks"}:
            raise RuntimeError("ranking receipt orchestration exact-key mismatch")
        if orch["session_identity"] != _session_identity(output_root, session) or orch["previous_chain_head"] != previous:
            raise RuntimeError("ranking receipt session/previous head mismatch")
        start_identity = orch["batch_start_identity"]
        start_path = Path(str(start_identity.get("path", "")))
        start = load_immutable_json(start_path, "ranking batch start")
        expected_before = 0 if index == 1 else LEGAL_CUMULATIVE_COUNTS[index - 2]
        expected_effective = count - expected_before
        expected_selected_ids = [
            _task_id(dataset, chrom)
            for dataset in DATASETS for chrom in AUTOSOMES
        ][expected_before:count]
        if (
            set(start) != BATCH_START_KEYS or start.get("batch_index") != index
            or start.get("session_id") != session.get("session_id")
            or sha256_path(start_path) != start_identity.get("sha256")
            or start.get("batch_id") != start_identity.get("batch_id")
            or start.get("previous_chain_head") != previous
            or start.get("before_count") != expected_before
            or start.get("max_new_tasks") != (8 if index == 1 else 16)
            or start.get("effective_count") != expected_effective
            or start.get("session_sha256") != sha256_path(output_root / "_orchestration" / "session.json")
            or start.get("run_contract_semantic_sha256") != session.get("run_contract_semantic_sha256")
            or start.get("runner_source") != session["producer_sources"]["runner"]
            or start.get("batch_id") != semantic_json_sha256({
                key: value for key, value in start.items()
                if key not in {"batch_id", "receipt_integrity"}
            })
            or start.get("selected_task_ids") != expected_selected_ids
        ):
            raise RuntimeError("ranking batch start/count chain mismatch")
        gate_target = {
            "output_root": dict(session["output_root"]),
            "session_id": session["session_id"],
            "session_sha256": sha256_path(output_root / "_orchestration" / "session.json"),
            "batch_index": index,
            "before_count": expected_before,
            "max_new_tasks": 8 if index == 1 else 16,
            "effective_count": expected_effective,
            "selected_task_ids": expected_selected_ids,
            "previous_chain_head": previous,
        }
        gate_payload = load_resource_gate_receipt(
            start.get("resource_gate") or {},
            expected_path=_resource_gate_path(output_root, "batch", index),
            expected_stage="ranking",
            expected_scope="batch",
            expected_target=gate_target,
            expected_producer_source=session["producer_sources"]["runner"],
        )
        if (
            gate_payload.get("host_boot_id") != start.get("host_boot_id")
            or gate_payload.get("observed_monotonic_ns", -1) > start.get("created_monotonic_ns", -1)
        ):
            raise RuntimeError("ranking batch is not temporally bound to resource gate")
        attestations = orch["batch_completion_attestations"]
        if [row.get("task_id") for row in attestations] != start.get("selected_task_ids"):
            raise RuntimeError("ranking completion order differs from batch start")
        for row in attestations:
            task_id = row["task_id"]
            if task_id in completions:
                raise RuntimeError("duplicate ranking completion")
            completion, identity = validate_attested_completion(
                Path(row["path"]), session, expected_task_id=task_id,
                expected_batch=start | {"sha256": sha256_path(start_path)},
            )
            if identity != {"path": row["path"], "sha256": row["sha256"]}:
                raise RuntimeError("ranking completion identity differs from chain receipt")
            completions[task_id] = identity
            child_path = Path(str(completion["child_receipt_identity"]["path"]))
            child_receipts[task_id] = load_immutable_json(
                child_path, f"attested ranking child receipt {task_id}"
            )
            child_receipt_paths[task_id] = str(child_path.resolve())
        result_rows = receipt.get("results")
        expected_cumulative_task_ids = [
            _task_id(dataset, chrom)
            for dataset in DATASETS for chrom in AUTOSOMES
        ][:count]
        if not isinstance(result_rows, list) or len(result_rows) != count:
            raise RuntimeError(
                "ranking receipt results do not cover the cumulative attested task count"
            )
        observed_result_task_ids = [
            _task_id(result.get("dataset"), result.get("chrom"))
            for result in result_rows
            if isinstance(result, dict)
        ]
        if observed_result_task_ids != expected_cumulative_task_ids:
            raise RuntimeError(
                "ranking receipt results are not an exact ordered bijection of cumulative child completions"
            )
        for result, task_id in zip(result_rows, observed_result_task_ids):
            if result.get("orchestration_completion") != completions.get(task_id):
                raise RuntimeError("ranking result lacks exact completion identity")
        if count < EXPECTED_TASKS:
            run_contract = receipt.get("run_contract")
            invocation = receipt.get("invocation")
            elapsed_seconds = receipt.get("elapsed_seconds")
            if (
                not isinstance(run_contract, dict)
                or semantic_json_sha256(run_contract)
                != session.get("run_contract_semantic_sha256")
                or not isinstance(invocation, dict)
                or not isinstance(elapsed_seconds, (int, float))
                or isinstance(elapsed_seconds, bool)
                or not math.isfinite(float(elapsed_seconds))
                or float(elapsed_seconds) < 0.0
            ):
                raise RuntimeError("ranking checkpoint rebuild inputs are invalid")
            current_batch_ids = set(start["selected_task_ids"])
            rebuilt_results = []
            for task_id in expected_cumulative_task_ids:
                dataset, chrom = task_id.split(":", 1)
                rebuilt_results.append({
                    "dataset": dataset,
                    "chrom": chrom,
                    "status": "PASS" if task_id in current_batch_ids else "REUSED_PASS",
                    "receipt": child_receipts[task_id],
                    "receipt_path": child_receipt_paths[task_id],
                    "orchestration_completion": completions[task_id],
                })
            rebuilt = build_ranking_checkpoint(
                rebuilt_results,
                run_contract=run_contract,
                invocation=invocation,
                elapsed_seconds=float(elapsed_seconds),
            )
            rebuilt["orchestration"] = orch
            rebuilt["receipt_integrity"] = _receipt_integrity(path)
            if receipt != rebuilt:
                raise RuntimeError(
                    "ranking checkpoint aggregate, checks, or canonical rows drifted"
                )
        previous = {"path": str(path.resolve()), "sha256": sha256_path(path)}
    starts = list((output_root / "_orchestration" / "batches").glob("batch_*_start.json")) if (output_root / "_orchestration" / "batches").exists() else []
    grants = list((output_root / "_orchestration" / "tasks").glob("*/*/grant.json")) if (output_root / "_orchestration" / "tasks").exists() else []
    completion_paths = list((output_root / "_orchestration" / "tasks").glob("*/*/completion.json")) if (output_root / "_orchestration" / "tasks").exists() else []
    gate_dir = output_root / "_orchestration" / "resource_gates"
    gate_paths = set(gate_dir.glob("*.json")) if gate_dir.exists() else set()
    expected_gate_paths = {_resource_gate_path(output_root, "session")}
    expected_gate_paths.update(
        _resource_gate_path(output_root, "batch", index)
        for index in range(1, len(rows) + 1)
    )
    if len(starts) != len(rows) or len(grants) != len(completions) or len(completion_paths) != len(completions):
        raise RuntimeError("open/orphan ranking batch, grant, or completion detected")
    if {path.resolve() for path in gate_paths} != {path.resolve() for path in expected_gate_paths}:
        raise RuntimeError("missing/orphan/reused ranking resource gate detected")
    return {
        "count": len(completions),
        "completions": completions,
        "child_receipts": child_receipts,
        "child_receipt_paths": child_receipt_paths,
        "previous_chain_head": previous,
        "next_batch_index": len(rows) + 1,
    }


def _release_require_exact_keys(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        raise RuntimeError(f"{label} exact-key schema mismatch")
    return value


def _release_validate_copy(record: Any, physical: Path, declared_path: str, label: str) -> dict[str, Any]:
    row = _release_require_exact_keys(record, RELEASE_FILE_IDENTITY_KEYS, label)
    observed = os.lstat(physical.absolute())
    if not stat.S_ISREG(observed.st_mode) or stat.S_ISLNK(observed.st_mode) or stat.S_IMODE(observed.st_mode) != 0o444 or observed.st_nlink != 1:
        raise RuntimeError(f"{label} is not an immutable single-link 0444 regular file")
    expected = {
        "path": declared_path, "st_dev": int(observed.st_dev), "st_ino": int(observed.st_ino),
        "st_nlink": int(observed.st_nlink), "size_bytes": int(observed.st_size),
        "mtime_ns": int(observed.st_mtime_ns), "ctime_ns": int(observed.st_ctime_ns),
        "mode_octal": "0444", "sha256": sha256_path(physical),
    }
    if row != expected:
        raise RuntimeError(f"{label} identity/stat/SHA mismatch")
    return row


def _release_validate_runtime(document: Mapping[str, Any]) -> None:
    runtime = _release_require_exact_keys(document.get("runtime"), {"python", "packages", "samtools", "platform"}, "runtime")
    python = _release_require_exact_keys(runtime["python"], {"executable", "version", "implementation"}, "runtime.python")
    packages = _release_require_exact_keys(runtime["packages"], {"numpy", "scipy", "pysam"}, "runtime.packages")
    samtools = _release_require_exact_keys(runtime["samtools"], {"path", "version_line", "htslib_version_line"}, "runtime.samtools")
    if not all(isinstance(value, str) and value for value in [*python.values(), *packages.values(), *samtools.values(), runtime["platform"]]):
        raise RuntimeError("release runtime contains an empty/non-string field")
    if document.get("runtime_semantic_sha256") != semantic_json_sha256(runtime):
        raise RuntimeError("release runtime semantic digest mismatch")


def _deep_verify_frozen_release(
    release_path: Path,
    release_sha256: str,
    identities: Mapping[str, Mapping[str, str]],
    *,
    force_fresh: bool = False,
) -> dict[str, Any]:
    """Execute only the freezer adjacent to this already-authenticated runner."""
    anchor_relative = Path(RELEASE_SOURCE_PATHS["full_ranking_runner"])
    anchor = RUNNER.resolve(strict=True)
    if tuple(anchor.parts[-len(anchor_relative.parts):]) != anchor_relative.parts:
        raise RuntimeError("release runner is not located at its hardcoded snapshot path")
    snapshot_root = anchor
    for _ in anchor_relative.parts:
        snapshot_root = snapshot_root.parent
    contract_root = release_path.parent.resolve(strict=True)
    if snapshot_root.parent != contract_root:
        raise RuntimeError("release manifest is not adjacent to this runner's code_snapshot")
    freezer = snapshot_root / RELEASE_SOURCE_PATHS["release_contract_freezer"]
    expected = identities.get("release_contract_freezer") or {}
    observed = os.lstat(freezer)
    if (
        str(freezer.resolve(strict=True)) != expected.get("path")
        or not stat.S_ISREG(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or stat.S_IMODE(observed.st_mode) != 0o444
        or observed.st_nlink != 1
        or sha256_path(freezer) != expected.get("sha256")
    ):
        raise RuntimeError("anchored frozen freezer identity mismatch")
    freezer_sha256 = sha256_path(freezer)
    cache_key = (release_sha256, freezer_sha256)
    if not force_fresh and cache_key in _DEEP_RELEASE_VERIFY_CACHE:
        return _DEEP_RELEASE_VERIFY_CACHE[cache_key]
    module_name = f"_m2_frozen_freezer_{freezer_sha256[:16]}"
    spec = importlib.util.spec_from_file_location(module_name, freezer)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen release verifier")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
        verification = module.verify_release_contract(contract_root)
    finally:
        sys.modules.pop(module_name, None)
    identity = verification.get("run_manifest") or {}
    if verification.get("all_pass") is not True or verification.get("validation_evidence_eligible") is not True or identity.get("sha256") != release_sha256:
        raise RuntimeError("frozen freezer deep verification did not certify the release")
    summary = {
        "mode": "FROZEN_FREEZER_VERIFY_RELEASE_CONTRACT",
        "freezer_path": str(freezer.resolve()), "freezer_sha256": freezer_sha256,
        "release_manifest_sha256": release_sha256,
        "verified_snapshot_files": (verification.get("snapshot") or {}).get("n_files"),
        "all_pass": True,
    }
    _DEEP_RELEASE_VERIFY_CACHE[cache_key] = summary
    return summary


def load_release_contract_binding(
    manifest_path: Path, *, required_sources: Mapping[str, Path],
    _skip_deep_verification_for_test: bool = False,
    _force_deep_reverification: bool = False,
) -> dict[str, Any]:
    """Authenticate the exact canonical release schema and physical frozen bytes."""
    raw_release_path = manifest_path.absolute()
    release_stat = os.lstat(raw_release_path)
    if not stat.S_ISREG(release_stat.st_mode) or stat.S_IMODE(release_stat.st_mode) != 0o444 or release_stat.st_nlink != 1:
        raise RuntimeError("release manifest is not an immutable single-link 0444 regular file")
    release_path = raw_release_path.resolve(strict=True)
    release_sha256 = sha256_path(release_path)
    deep_verification = (
        {"mode": "TEST_ONLY_SKIPPED", "all_pass": False}
        if _skip_deep_verification_for_test else None
    )
    sidecar_path = release_path.with_name(f"{release_path.name}.sha256")
    sidecar_stat = os.lstat(sidecar_path)
    if not stat.S_ISREG(sidecar_stat.st_mode) or stat.S_IMODE(sidecar_stat.st_mode) != 0o444 or sidecar_stat.st_nlink != 1:
        raise RuntimeError("release sidecar is not an immutable single-link 0444 regular file")
    fields = sidecar_path.read_text(encoding="ascii", errors="strict").strip().split()
    if len(fields) != 2 or fields != [release_sha256, release_path.name]:
        raise RuntimeError("release-contract sidecar mismatch")
    document = json.loads(release_path.read_text(encoding="utf-8"))
    _release_require_exact_keys(document, RELEASE_TOP_KEYS, "release manifest")
    expected_scope = {
        "technical_datasets": 7, "biological_samples": 6, "chromosomes": 22,
        "chromosome_names": list(AUTOSOMES), "tasks": EXPECTED_TASKS,
        "datasets": list(DATASETS),
    }
    if (
        document["schema_name"] != RELEASE_SCHEMA_NAME or document["schema_version"] != RELEASE_SCHEMA_VERSION
        or document["task_type"] != "B_COMPREHENSIVE_VALIDATION"
        or document["authority_mode"] != "CANONICAL_V5_FROZEN"
        or document["validation_evidence_eligible"] is not True or document["all_pass"] is not True
        or not isinstance(document["created_at_utc"], str) or not document["created_at_utc"].endswith("Z")
        or document["scope"] != expected_scope
    ):
        raise RuntimeError("release canonical authority/task/scope metadata mismatch")
    integrity = _release_require_exact_keys(document["receipt_integrity"], {"scheme", "sidecar_name", "covers"}, "release receipt_integrity")
    if integrity != {"scheme": "external_sha256_sidecar_v1", "sidecar_name": f"{release_path.name}.sha256", "covers": release_path.name}:
        raise RuntimeError("release receipt_integrity filename binding mismatch")
    checks = _release_require_exact_keys(document["checks"], RELEASE_CHECK_NAMES, "release checks")
    if any(value is not True for value in checks.values()):
        raise RuntimeError("release contract contains a failing named check")
    if document["parameters"] != FROZEN_RELEASE_PARAMETERS:
        raise RuntimeError("release frozen parameters differ from canonical M2 parameters")
    parameters = document["parameters"]
    if document["assurance"] != RELEASE_ASSURANCE:
        raise RuntimeError("release assurance declaration mismatch")
    _release_validate_runtime(document)
    snapshot = _release_require_exact_keys(document["source_snapshot"], {"repo_root", "snapshot_root", "n_files", "entries", "entrypoints", "exact_allowlist_semantic_sha256"}, "source_snapshot")
    if snapshot.get("snapshot_root") != "code_snapshot":
        raise RuntimeError("release snapshot root mismatch")
    if not Path(str(snapshot["repo_root"])).is_absolute():
        raise RuntimeError("release snapshot repo_root is not absolute")
    snapshot_root_path = release_path.parent / "code_snapshot"
    if snapshot_root_path.is_symlink():
        raise RuntimeError("release snapshot root must not be a symlink")
    resolved_snapshot_root = snapshot_root_path.resolve(strict=True)
    entries = snapshot.get("entries")
    if not isinstance(entries, list) or len(entries) != 11 or snapshot.get("n_files") != 11:
        raise RuntimeError("release snapshot cardinality is malformed")
    expected_code_entrypoints = {
        role: (Path("code_snapshot") / relative).as_posix()
        for role, relative in RELEASE_SOURCE_PATHS.items() if role != "canonical_manifest_schema"
    }
    entrypoints = snapshot["entrypoints"]
    if entrypoints != expected_code_entrypoints:
        raise RuntimeError("release snapshot exact entrypoint map mismatch")
    identities: dict[str, dict[str, str]] = {}
    producer_records: dict[str, dict[str, Any]] = {}
    for entry in entries:
        _release_require_exact_keys(entry, {"role", "repo_relative_path", "source", "snapshot"}, "snapshot entry")
        role = entry["role"]
        if role not in RELEASE_SOURCE_PATHS or role in identities:
            raise RuntimeError(f"duplicate/empty release snapshot role: {role!r}")
        repo_relative = RELEASE_SOURCE_PATHS[role]
        if entry["repo_relative_path"] != repo_relative:
            raise RuntimeError(f"release snapshot repo path mismatch: {role}")
        source = _release_require_exact_keys(entry["source"], RELEASE_FILE_IDENTITY_KEYS, f"snapshot source/{role}")
        if Path(str(source["path"])).absolute() != Path(str(snapshot["repo_root"])) / repo_relative:
            raise RuntimeError(f"release source origin path mismatch: {role}")
        if source["sha256"] != entry["snapshot"].get("sha256") or source["st_nlink"] != 1:
            raise RuntimeError(f"release source/copy identity mismatch: {role}")
        relative = Path("code_snapshot") / repo_relative
        physical_raw = release_path.parent / relative
        physical = physical_raw.resolve(strict=True)
        if not physical.is_relative_to(resolved_snapshot_root):
            raise RuntimeError(f"release snapshot path escapes its root: {role}")
        copied = _release_validate_copy(entry["snapshot"], physical_raw, relative.as_posix(), f"snapshot/{role}")
        identities[role] = {"path": str(physical), "sha256": copied["sha256"]}
        producer_records[role] = entry
    if set(identities) != set(RELEASE_SOURCE_PATHS):
        raise RuntimeError("release snapshot role set is not the exact canonical eleven")
    expected_allowlist = [{"role": role, "repo_relative_path": relative} for role, relative in RELEASE_SOURCE_PATHS.items()]
    if snapshot["exact_allowlist_semantic_sha256"] != semantic_json_sha256(expected_allowlist):
        raise RuntimeError("release allowlist semantic digest mismatch")
    for role, actual_path in required_sources.items():
        if role not in identities:
            raise RuntimeError(f"required release snapshot role is missing: {role}")
        if role not in entrypoints:
            raise RuntimeError(f"required release snapshot entrypoint is missing: {role}")
        resolved_actual = actual_path.resolve(strict=True)
        expected = identities[role]
        if str(resolved_actual) != expected["path"] or sha256_path(resolved_actual) != expected["sha256"]:
            raise RuntimeError(f"current source does not match frozen snapshot role: {role}")
    canonical = _release_require_exact_keys(document["canonical_manifest"], {"expected_sha256", "origin", "immutable_copy"}, "canonical_manifest")
    if canonical["expected_sha256"] != CANONICAL_MANIFEST_SHA256:
        raise RuntimeError("release canonical manifest authority SHA mismatch")
    _release_require_exact_keys(canonical["origin"], RELEASE_FILE_IDENTITY_KEYS, "canonical_manifest.origin")
    relative = Path("input_contract/canonical_manifest.json")
    canonical_path_raw = release_path.parent / relative
    canonical_path = canonical_path_raw.resolve(strict=True)
    input_root = release_path.parent / "input_contract"
    if input_root.is_symlink() or not canonical_path.is_relative_to(input_root.resolve(strict=True)):
        raise RuntimeError("immutable canonical manifest escapes input_contract")
    copied_manifest = _release_validate_copy(canonical["immutable_copy"], canonical_path_raw, relative.as_posix(), "canonical_manifest.copy")
    copy_identity = {"path": str(canonical_path), "sha256": copied_manifest["sha256"]}
    if copy_identity["sha256"] != CANONICAL_MANIFEST_SHA256:
        raise RuntimeError("release canonical manifest copy SHA mismatch")
    schema = _release_require_exact_keys(document["canonical_schema"], {"role", "repo_relative_path", "sha256"}, "canonical_schema")
    if schema != {"role": "canonical_manifest_schema", "repo_relative_path": CANONICAL_SCHEMA_RELATIVE, "sha256": CANONICAL_SCHEMA_SHA256}:
        raise RuntimeError("canonical schema authority mismatch")
    pre = _release_require_exact_keys(document["pre_input_identity_receipt"], {"origin", "immutable_copy", "input_identity_snapshot_sha256", "receipt_semantic_sha256", "authority_mode", "validation_evidence_eligible", "artifact_roles"}, "pre_input_identity_receipt")
    if pre["authority_mode"] != "CANONICAL_V5_FROZEN" or pre["validation_evidence_eligible"] is not True or pre["artifact_roles"] != 42:
        raise RuntimeError("PRE authority/scope mismatch")
    _release_require_exact_keys(pre["origin"], RELEASE_FILE_IDENTITY_KEYS | {"sidecar"}, "PRE.origin")
    pre_copy = _release_require_exact_keys(pre["immutable_copy"], RELEASE_FILE_IDENTITY_KEYS | {"sidecar"}, "PRE.copy")
    pre_relative = Path(str(pre_copy["path"]))
    if (
        pre_relative.is_absolute() or ".." in pre_relative.parts
        or not pre_relative.parts or pre_relative.parts[0] != "input_contract"
    ):
        raise RuntimeError("PRE copy path escapes input_contract")
    pre_physical_raw = release_path.parent / pre_relative
    pre_physical = pre_physical_raw.resolve(strict=True)
    if not pre_physical.is_relative_to(input_root.resolve(strict=True)):
        raise RuntimeError("PRE copy path escapes input_contract")
    _release_validate_copy(
        {key: pre_copy[key] for key in RELEASE_FILE_IDENTITY_KEYS},
        pre_physical_raw,
        pre_relative.as_posix(),
        "PRE.copy",
    )
    pre_side = pre_copy["sidecar"]
    pre_side_relative = Path(str(pre_side.get("path", "")))
    if (
        pre_side_relative.is_absolute() or ".." in pre_side_relative.parts
        or not pre_side_relative.parts or pre_side_relative.parts[0] != "input_contract"
    ):
        raise RuntimeError("PRE sidecar path escapes input_contract")
    pre_side_raw = release_path.parent / pre_side_relative
    pre_side_physical = pre_side_raw.resolve(strict=True)
    if not pre_side_physical.is_relative_to(input_root.resolve(strict=True)):
        raise RuntimeError("PRE sidecar path escapes input_contract")
    _release_validate_copy(
        pre_side, pre_side_raw, pre_side_relative.as_posix(), "PRE.copy.sidecar"
    )
    pre_side_fields = pre_side_raw.read_text(
        encoding="ascii", errors="strict"
    ).strip().split()
    if pre_side_fields != [sha256_path(pre_physical_raw), pre_physical_raw.name]:
        raise RuntimeError("PRE copied receipt sidecar mismatch")
    fresh_keys = {"execution_mode", "verifier_path", "verifier_sha256", "receipt_sha256", "receipt_semantic_sha256", "checks_semantic_sha256", "artifact_audit_semantic_sha256", "input_identity_snapshot_sha256", "all_pass", "validation_evidence_eligible", "exactly_equals_supplied_pre_snapshot"}
    fresh = _release_require_exact_keys(document["fresh_input_identity_verification"], fresh_keys, "fresh_input_identity_verification")
    if fresh["execution_mode"] != "production_subprocess_required" or fresh["all_pass"] is not True or fresh["validation_evidence_eligible"] is not True or fresh["exactly_equals_supplied_pre_snapshot"] is not True or fresh["input_identity_snapshot_sha256"] != pre["input_identity_snapshot_sha256"] or fresh["verifier_path"] != producer_records["input_identity_verifier"]["source"]["path"] or fresh["verifier_sha256"] != producer_records["input_identity_verifier"]["source"]["sha256"]:
        raise RuntimeError("fresh input identity verification summary mismatch")
    producer = _release_require_exact_keys(document["producer"], {"role", "repo_relative_path", "source_sha256", "immutable_copy_path", "immutable_copy_sha256"}, "producer")
    freezer_entry = producer_records["release_contract_freezer"]
    if producer != {"role": "release_contract_freezer", "repo_relative_path": RELEASE_SOURCE_PATHS["release_contract_freezer"], "source_sha256": freezer_entry["source"]["sha256"], "immutable_copy_path": freezer_entry["snapshot"]["path"], "immutable_copy_sha256": freezer_entry["snapshot"]["sha256"]}:
        raise RuntimeError("release producer is not bound to frozen freezer source")
    expected_top_entrypoints = {**expected_code_entrypoints, "canonical_manifest_copy": "input_contract/canonical_manifest.json", "pre_input_identity_receipt_copy": pre_relative.as_posix(), "canonical_schema_copy": (Path("code_snapshot") / CANONICAL_SCHEMA_RELATIVE).as_posix()}
    if document["entrypoints"] != expected_top_entrypoints:
        raise RuntimeError("release top-level exact entrypoint map mismatch")
    if deep_verification is None:
        deep_verification = _deep_verify_frozen_release(
            release_path, release_sha256, identities,
            force_fresh=_force_deep_reverification,
        )
    return {
        "schema_name": "intersubmod.m2_release_binding",
        "schema_version": "1.0.0",
        "release_manifest": {
            "path": str(release_path),
            "sha256": release_sha256,
            "semantic_sha256": semantic_json_sha256(document),
            "sidecar": {"path": str(sidecar_path.resolve()), "sha256": sha256_path(sidecar_path)},
        },
        "authority_mode": document.get("authority_mode"),
        "validation_evidence_eligible": deep_verification.get("all_pass") is True,
        "canonical_input_manifest": copy_identity,
        "snapshot_sources": dict(sorted(identities.items())),
        "frozen_parameters": parameters,
        "frozen_parameters_semantic_sha256": semantic_json_sha256(parameters),
        "deep_release_verification": deep_verification,
    }


def validate_release_ranking_parameters(
    binding: Mapping[str, Any], args: argparse.Namespace
) -> None:
    frozen = binding["frozen_parameters"]
    ranking = frozen.get("ranking") or {}
    scheduler = frozen.get("scheduler") or {}
    expected = {
        "thresholds": list(args.thresholds),
        "component_bases": list(args.component_bases),
        "hp_families": list(args.hp_families),
        "structural_exact_pattern_minread_grid": list(args.structural_exact_pattern_minreads),
        "primary_structural_exact_pattern_minread": args.primary_structural_exact_pattern_minread,
        "exact_k_max": args.exact_k_max,
        "max_vertex_sets": args.max_vertex_sets,
        "solver_time_limit_seconds_per_milp": args.solver_time_limit_seconds,
        "fixed_error_grid": list(args.fixed_error_grid_values),
        "minimum_bq_error_rate": args.minimum_bq_error_rate,
        "maximum_bq_error_rate": args.maximum_bq_error_rate,
        "conditional_candidate_ranking_bootstrap_replicates": (
            args.conditional_candidate_ranking_bootstrap_replicates
        ),
        "conditional_candidate_ranking_bootstrap_seed": (
            args.conditional_candidate_ranking_bootstrap_seed
        ),
        "tie_tolerance_log_likelihood": args.tie_tolerance,
        "workers": args.workers,
    }
    if ranking != expected:
        raise RuntimeError("ranking CLI parameters drift from the frozen release contract")
    if (
        args.task_timeout_seconds != scheduler.get("task_timeout_seconds")
        or args.timeout_grace_seconds != scheduler.get("timeout_grace_seconds")
        or args.max_new_tasks not in {
            scheduler.get("initial_batch_tasks"), scheduler.get("subsequent_batch_tasks")
        }
    ):
        raise RuntimeError("ranking scheduler/batch parameters drift from the frozen release contract")
    if args.aggregate_only:
        raise RuntimeError("--aggregate-only is forbidden for a release-bound ranking run")


def require_matching_upstream_release_binding(
    full_extraction: Mapping[str, Any], release_binding: Mapping[str, Any]
) -> None:
    upstream_binding = (full_extraction.get("run_contract") or {}).get("release_binding")
    if not isinstance(upstream_binding, dict):
        raise RuntimeError("release-bound ranking requires a release-bound full extraction")
    if (
        (upstream_binding.get("release_manifest") or {}).get("sha256")
        != (release_binding.get("release_manifest") or {}).get("sha256")
        or upstream_binding != release_binding
    ):
        raise RuntimeError("extraction and ranking are bound to different release contracts")


def _runtime_summary_contract_valid(summary: Mapping[str, Any]) -> bool:
    try:
        n = int(summary.get("n", -1))
        if n < 0 or set(summary) != {"n", "sum", "max", "p50", "p95", "p99"}:
            return False
        if n == 0:
            return summary == {
                "n": 0, "sum": 0.0, "max": None, "p50": None, "p95": None, "p99": None,
            }
        values = [float(summary[key]) for key in ("sum", "max", "p50", "p95", "p99")]
        return all(math.isfinite(value) and value >= 0.0 for value in values)
    except (TypeError, ValueError):
        return False


def child_runtime_contract_valid(receipt: Mapping[str, Any]) -> bool:
    runtime = receipt.get("runtime_diagnostics") or {}
    scopes = runtime.get("scopes") or {}
    invoked = runtime.get("primary_invoked_segment_scopes") or {}
    if (
        runtime.get("schema_name") != "intersubmod.m2_unit_runtime_diagnostics"
        or runtime.get("schema_version") != "1.0.0"
        or runtime.get("clock") != "time.monotonic_ns"
        or runtime.get("unit") != "seconds"
        or set(scopes) != {
            "primary_unit_evaluations", "all_structural_minread_unit_evaluations",
        }
        or set(invoked) != {
            "candidate_generation_elapsed_seconds", "likelihood_fit_elapsed_seconds",
        }
    ):
        return False
    output_name = runtime.get("per_unit_output")
    identity = (receipt.get("outputs") or {}).get(output_name)
    if not output_name or not identity:
        return False
    try:
        output_path = Path(identity["path"])
        with gzip.open(output_path, "rt", encoding="utf-8", newline="") as handle:
            header = tuple(next(csv.reader(handle, delimiter="\t")))
        if header != RUNTIME_DIAGNOSTIC_FIELDS:
            return False
        for scope in scopes.values():
            n_units = int(scope.get("n_unit_evaluations", -1))
            if n_units < 0 or set(scope) != {"n_unit_evaluations", *RUNTIME_METRICS}:
                return False
            if any(
                not _runtime_summary_contract_valid(scope.get(metric) or {})
                or int(scope[metric]["n"]) != n_units
                for metric in RUNTIME_METRICS
            ):
                return False
        if any(not _runtime_summary_contract_valid(summary) for summary in invoked.values()):
            return False
    except (OSError, StopIteration, KeyError, TypeError, ValueError):
        return False
    return True


def child_method_contract_and_ranker_source_bound(
    receipt: Mapping[str, Any], ranker_path: Path, ranker_sha256: str
) -> bool:
    """Validate one child against the independent full-run method contract."""
    try:
        resolved_ranker = ranker_path.resolve(strict=True)
        provenance = (receipt.get("provenance") or {}).get("ranker") or {}
        declared_path = Path(str(provenance.get("path", ""))).resolve(strict=True)
        return (
            (receipt.get("parameters") or {}).get("method_contract")
            == EXPECTED_METHOD_CONTRACT
            and declared_path == resolved_ranker
            and provenance.get("sha256") == ranker_sha256
            and sha256_path(resolved_ranker) == ranker_sha256
        )
    except (OSError, TypeError, ValueError):
        return False


def load_full_extraction_receipt(extraction_root: Path) -> tuple[dict[str, Any], Path]:
    path = extraction_root / "full_extraction_receipt.json"
    if not path.is_file():
        raise FileNotFoundError(f"missing full extraction receipt: {path}")
    receipt = json.loads(path.read_text(encoding="utf-8"))
    scope = receipt.get("scope") or {}
    if (
        receipt.get("schema_name") != "intersubmod.m2_full_extraction_receipt"
        or receipt.get("schema_version") != "1.2.0"
        or receipt.get("all_pass") is not True
        or scope.get("datasets") != list(DATASETS)
        or scope.get("chromosomes") != list(AUTOSOMES)
        or scope.get("expected_tasks") != 154
        or receipt.get("n_results") != 154
    ):
        raise RuntimeError("upstream full extraction receipt is incomplete or has the wrong scope")
    integrity = receipt.get("receipt_integrity") or {}
    if integrity.get("scheme") != "external_sha256_sidecar_v1" or not verify_sha256_sidecar(
        path, integrity.get("sidecar_name")
    ):
        raise RuntimeError("upstream full extraction receipt checksum failed")

    results = receipt.get("results")
    if not isinstance(results, list) or len(results) != 154:
        raise RuntimeError("upstream full extraction receipt does not embed exactly 154 child results")
    expected_pairs = {(dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES}
    observed_pairs: list[tuple[str, str]] = []
    for row in results:
        if not isinstance(row, dict):
            raise RuntimeError("upstream full extraction result is not an object")
        dataset = row.get("dataset")
        chrom = row.get("chrom")
        child = row.get("receipt")
        if (
            row.get("status") not in {"PASS", "REUSED_PASS"}
            or not isinstance(dataset, str)
            or not isinstance(chrom, str)
            or not isinstance(child, dict)
        ):
            raise RuntimeError("upstream full extraction result is not a passing embedded child")
        child_scope = child.get("scope") or {}
        if (
            child.get("schema_name")
            != "intersubmod.lossless_read_linkage_chromosome_receipt"
            or child.get("schema_version") != "1.2.0"
            or child.get("all_pass") is not True
            or child_scope.get("dataset") != dataset
            or child_scope.get("chrom") != chrom
        ):
            raise RuntimeError("upstream full extraction child receipt contract or scope mismatch")
        observed_pairs.append((dataset, chrom))
    if len(set(observed_pairs)) != 154 or set(observed_pairs) != expected_pairs:
        raise RuntimeError("upstream full extraction child pairs are not the exact canonical 154 tasks")
    return receipt, path


def _rank_parameters(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "structural_exact_pattern_minread_grid": list(args.structural_exact_pattern_minreads),
        "primary_structural_exact_pattern_minread": args.primary_structural_exact_pattern_minread,
        "scoring_minread": 1,
        "exact_k_max": args.exact_k_max,
        "max_vertex_sets": args.max_vertex_sets,
        "solver_time_limit_seconds_per_milp": args.solver_time_limit_seconds,
        "minimum_bq_error_rate": args.minimum_bq_error_rate,
        "maximum_bq_error_rate": args.maximum_bq_error_rate,
        "fixed_error_grid_conditional_binary_flip_probability": list(args.fixed_error_grid_values),
        "conditional_candidate_ranking_bootstrap_replicates": args.conditional_candidate_ranking_bootstrap_replicates,
        "conditional_candidate_ranking_bootstrap_base_seed": args.conditional_candidate_ranking_bootstrap_seed,
        "tie_tolerance_log_likelihood": args.tie_tolerance,
    }


def load_passing_rank_receipt(
    path: Path,
    *,
    dataset: str,
    chrom: str,
    extraction_receipt_path: Path,
    ranker_path: Path,
    ranker_sha256: str,
    expected_parameters: Mapping[str, Any],
    expected_bases: Sequence[str],
    expected_thresholds: Sequence[int],
    expected_hp_families: Sequence[str],
    expected_component_basis_mode: str,
    expected_extraction_receipt_semantic_sha256: str,
) -> dict[str, Any] | None:
    if (path.parent / TASK_STATUS_NAME).exists():
        return None
    if not path.is_file():
        return None
    try:
        extraction_receipt = json.loads(extraction_receipt_path.read_text(encoding="utf-8"))
        extraction_integrity = extraction_receipt.get("receipt_integrity") or {}
        if (
            semantic_json_sha256(extraction_receipt) != expected_extraction_receipt_semantic_sha256
            or extraction_integrity.get("scheme") != "external_sha256_sidecar_v1"
            or not verify_sha256_sidecar(
                extraction_receipt_path, extraction_integrity.get("sidecar_name")
            )
        ):
            return None
        if (
            extraction_receipt.get("schema_name")
            != "intersubmod.lossless_read_linkage_chromosome_receipt"
            or extraction_receipt.get("schema_version") != "1.2.0"
            or extraction_receipt.get("all_pass") is not True
        ):
            return None
        extraction_dir = extraction_receipt_path.parent.resolve()
        required_suffixes = (
            ".molecule_sparse_calls.tsv.gz",
            ".site_catalog.tsv.gz",
            ".components.tsv.gz",
            ".site_component_membership.tsv.gz",
        )
        extraction_outputs = extraction_receipt.get("outputs") or {}
        for suffix in required_suffixes:
            matches = [identity for name, identity in extraction_outputs.items() if name.endswith(suffix)]
            if len(matches) != 1:
                return None
            identity = matches[0]
            output = Path(identity["path"])
            if (
                output.resolve().parent != extraction_dir
                or not output.is_file()
                or output.stat().st_size != identity["size_bytes"]
                or sha256_path(output) != identity["sha256"]
            ):
                return None
        receipt = json.loads(path.read_text(encoding="utf-8"))
        scope = receipt.get("scope") or {}
        if (
            receipt.get("schema_name") != "intersubmod.m2_symbolic_patterns_vertex_rank_receipt"
            or receipt.get("schema_version") != "2.0.0"
            or receipt.get("all_pass") is not True
            or scope.get("dataset") != [dataset]
            or scope.get("chrom") != [chrom]
            or sorted(scope.get("component_bases") or []) != sorted(expected_bases)
            or sorted(scope.get("thresholds") or []) != sorted(expected_thresholds)
            or sorted(scope.get("hp_families") or []) != sorted(expected_hp_families)
            or scope.get("component_basis_mode") != expected_component_basis_mode
        ):
            return None
        parameters = receipt.get("parameters") or {}
        if any(parameters.get(key) != value for key, value in expected_parameters.items()):
            return None
        provenance = receipt.get("provenance") or {}
        if not child_method_contract_and_ranker_source_bound(
            receipt, ranker_path, ranker_sha256
        ):
            return None
        upstream = provenance.get("upstream_extraction_receipt") or {}
        if Path(upstream.get("path", "")).resolve() != extraction_receipt_path.resolve():
            return None
        if upstream.get("sha256") != sha256_path(extraction_receipt_path):
            return None
        integrity = receipt.get("receipt_integrity") or {}
        if integrity.get("scheme") != "external_sha256_sidecar_v1" or not verify_sha256_sidecar(
            path, integrity.get("sidecar_name")
        ):
            return None
        receipt_dir = path.parent.resolve()
        for identity in (receipt.get("outputs") or {}).values():
            output = Path(identity["path"])
            if output.resolve().parent != receipt_dir or not output.is_file():
                return None
            if output.stat().st_size != identity["size_bytes"] or sha256_path(output) != identity["sha256"]:
                return None
        if not child_runtime_contract_valid(receipt):
            return None
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return None
    return receipt


def task_command(
    extraction_dir: Path,
    output_dir: Path,
    args: argparse.Namespace,
) -> list[str]:
    command = [
        sys.executable,
        str(RANKER),
        "--input-dir",
        str(extraction_dir),
        "--output-dir",
        str(output_dir),
        "--structural-exact-pattern-minread-grid",
        ",".join(str(value) for value in args.structural_exact_pattern_minreads),
        "--primary-structural-exact-pattern-minread",
        str(args.primary_structural_exact_pattern_minread),
        "--exact-k-max",
        str(args.exact_k_max),
        "--max-vertex-sets",
        str(args.max_vertex_sets),
        "--solver-time-limit-seconds",
        str(args.solver_time_limit_seconds),
        "--fixed-error-grid",
        ",".join(format(value, ".12g") for value in args.fixed_error_grid_values),
        "--minimum-bq-error-rate",
        str(args.minimum_bq_error_rate),
        "--maximum-bq-error-rate",
        str(args.maximum_bq_error_rate),
        "--conditional-candidate-ranking-bootstrap-replicates",
        str(args.conditional_candidate_ranking_bootstrap_replicates),
        "--conditional-candidate-ranking-bootstrap-seed",
        str(args.conditional_candidate_ranking_bootstrap_seed),
        "--tie-tolerance",
        str(args.tie_tolerance),
    ]
    for threshold in args.thresholds:
        command.extend(("--threshold", str(threshold)))
    for basis in args.component_bases:
        command.extend(("--component-basis", basis))
    for family in args.hp_families:
        command.extend(("--hp-family", family))
    return command


def run_task(spec: tuple[Any, ...]) -> dict[str, Any]:
    (
        dataset,
        chrom,
        extraction_dir,
        output_dir,
        command,
        ranker_path,
        ranker_sha256,
        expected_parameters,
        expected_bases,
        expected_thresholds,
        expected_hp_families,
        expected_component_basis_mode,
        expected_extraction_receipt_semantic_sha256,
        aggregate_only,
        task_timeout_seconds,
        timeout_grace_seconds,
    ) = spec
    extraction_receipt_path = extraction_dir / "receipt.json"
    receipt_path = output_dir / "receipt.json"
    existing = load_passing_rank_receipt(
        receipt_path,
        dataset=dataset,
        chrom=chrom,
        extraction_receipt_path=extraction_receipt_path,
        ranker_path=ranker_path,
        ranker_sha256=ranker_sha256,
        expected_parameters=expected_parameters,
        expected_bases=expected_bases,
        expected_thresholds=expected_thresholds,
        expected_hp_families=expected_hp_families,
        expected_component_basis_mode=expected_component_basis_mode,
        expected_extraction_receipt_semantic_sha256=expected_extraction_receipt_semantic_sha256,
    )
    if existing is not None:
        return {
            "dataset": dataset,
            "chrom": chrom,
            "status": "REUSED_PASS",
            "receipt": existing,
            "receipt_path": str(receipt_path),
        }
    if aggregate_only:
        return {"dataset": dataset, "chrom": chrom, "status": "FAIL_MISSING_OR_STALE_RANK_RECEIPT"}
    if output_dir.exists():
        return {"dataset": dataset, "chrom": chrom, "status": "FAIL_EXISTING_NONPASS_DIRECTORY"}
    started = time.monotonic()
    started_monotonic_ns = time.monotonic_ns()
    process = run_command_with_process_group_timeout(
        command,
        task_timeout_seconds=task_timeout_seconds,
        timeout_grace_seconds=timeout_grace_seconds,
    )
    if process["timed_out"]:
        output_dir.mkdir(parents=True, exist_ok=True)
        marker_path = output_dir / TASK_STATUS_NAME
        elapsed_seconds = time.monotonic() - started
        marker = {
            "schema_name": "intersubmod.m2_runner_task_status",
            "schema_version": "1.0.0",
            "dataset": dataset,
            "chrom": chrom,
            "status": "TIMEOUT",
            "returncode": process["returncode"],
            "elapsed_seconds": elapsed_seconds,
            "stdout_tail": process["stdout"][-4000:],
            "stderr_tail": process["stderr"][-4000:],
            "task_timeout_seconds": task_timeout_seconds,
            "timeout_grace_seconds": timeout_grace_seconds,
            "process_group_id": process["process_group_id"],
            "termination_stage": process["termination_stage"],
            "command": command,
        }
        with marker_path.open("x", encoding="utf-8") as handle:
            json.dump(marker, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        return {
            "dataset": dataset,
            "chrom": chrom,
            "status": "TIMEOUT",
            "returncode": process["returncode"],
            "elapsed_seconds": elapsed_seconds,
            "stdout_tail": process["stdout"][-4000:],
            "stderr_tail": process["stderr"][-4000:],
            "task_dir": str(output_dir),
            "task_status_path": str(marker_path),
            "command": command,
            "process_group_id": process["process_group_id"],
            "termination_stage": process["termination_stage"],
        }
    receipt = load_passing_rank_receipt(
        receipt_path,
        dataset=dataset,
        chrom=chrom,
        extraction_receipt_path=extraction_receipt_path,
        ranker_path=ranker_path,
        ranker_sha256=ranker_sha256,
        expected_parameters=expected_parameters,
        expected_bases=expected_bases,
        expected_thresholds=expected_thresholds,
        expected_hp_families=expected_hp_families,
        expected_component_basis_mode=expected_component_basis_mode,
        expected_extraction_receipt_semantic_sha256=expected_extraction_receipt_semantic_sha256,
    )
    status = "PASS" if process["returncode"] == 0 and receipt is not None else "FAIL"
    marker_path: Path | None = None
    if status == "FAIL":
        output_dir.mkdir(parents=True, exist_ok=True)
        marker_path = output_dir / TASK_STATUS_NAME
        marker = {
            "schema_name": "intersubmod.m2_runner_task_status",
            "schema_version": "1.0.0",
            "dataset": dataset,
            "chrom": chrom,
            "status": "CHILD_FAILED_OR_INVALID_RECEIPT",
            "returncode": process["returncode"],
            "elapsed_seconds": time.monotonic() - started,
            "stdout_tail": process["stdout"][-4000:],
            "stderr_tail": process["stderr"][-4000:],
            "command": command,
        }
        with marker_path.open("x", encoding="utf-8") as handle:
            json.dump(marker, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
    return {
        "dataset": dataset,
        "chrom": chrom,
        "status": status,
        "returncode": process["returncode"],
        "elapsed_seconds": time.monotonic() - started,
        "stdout_tail": process["stdout"][-4000:],
        "stderr_tail": process["stderr"][-4000:],
        "command": command,
        "receipt": receipt,
        "receipt_path": str(receipt_path) if receipt is not None else None,
        "task_status_path": None if marker_path is None else str(marker_path),
        "process_group_id": process["process_group_id"],
        "started_monotonic_ns": started_monotonic_ns,
    }


def run_command_with_process_group_timeout(
    command: Sequence[str], *, task_timeout_seconds: float, timeout_grace_seconds: float
) -> dict[str, Any]:
    """Run one child in a new process group and terminate the whole tree on timeout."""
    process = subprocess.Popen(
        list(command),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    timed_out = False
    termination_stage = "NONE"
    try:
        stdout, stderr = process.communicate(timeout=task_timeout_seconds)
    except subprocess.TimeoutExpired:
        timed_out = True
        termination_stage = "TERM"
        grace_deadline = time.monotonic() + timeout_grace_seconds
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            stdout, stderr = process.communicate(timeout=timeout_grace_seconds)
        except subprocess.TimeoutExpired:
            stdout = stderr = None
        while _process_group_exists(process.pid) and time.monotonic() < grace_deadline:
            time.sleep(min(0.05, max(0.0, grace_deadline - time.monotonic())))
        if _process_group_exists(process.pid):
            termination_stage = "KILL"
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        if stdout is None or stderr is None:
            stdout, stderr = process.communicate()
    return {
        "returncode": process.returncode,
        "stdout": stdout or "",
        "stderr": stderr or "",
        "timed_out": timed_out,
        "process_group_id": process.pid,
        "termination_stage": termination_stage,
    }


def _process_group_exists(process_group_id: int) -> bool:
    try:
        os.killpg(process_group_id, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def scan_existing_rank_specs(
    specs: Sequence[tuple[Any, ...]],
    *,
    attested_completions: Mapping[str, Mapping[str, str]] | None = None,
) -> tuple[list[dict[str, Any]], list[tuple[Any, ...]], list[dict[str, str]]]:
    """Validate all existing rank children before scheduling any pending task."""
    reused: list[dict[str, Any]] = []
    pending: list[tuple[Any, ...]] = []
    invalid: list[dict[str, str]] = []
    for spec in specs:
        (
            dataset,
            chrom,
            extraction_dir,
            output_dir,
            _,
            ranker_path,
            ranker_sha256,
            expected_parameters,
            expected_bases,
            expected_thresholds,
            expected_hp_families,
            expected_component_basis_mode,
            expected_extraction_receipt_semantic_sha256,
            aggregate_only,
            _,
            _,
        ) = spec
        if not output_dir.exists():
            if aggregate_only:
                invalid.append({
                    "dataset": dataset,
                    "chrom": chrom,
                    "task_dir": str(output_dir),
                    "status": "FAIL_MISSING_OR_STALE_RANK_RECEIPT",
                })
            else:
                pending.append(spec)
            continue
        receipt = load_passing_rank_receipt(
            output_dir / "receipt.json",
            dataset=dataset,
            chrom=chrom,
            extraction_receipt_path=extraction_dir / "receipt.json",
            ranker_path=ranker_path,
            ranker_sha256=ranker_sha256,
            expected_parameters=expected_parameters,
            expected_bases=expected_bases,
            expected_thresholds=expected_thresholds,
            expected_hp_families=expected_hp_families,
            expected_component_basis_mode=expected_component_basis_mode,
            expected_extraction_receipt_semantic_sha256=expected_extraction_receipt_semantic_sha256,
        )
        if receipt is None:
            invalid.append({
                "dataset": dataset,
                "chrom": chrom,
                "task_dir": str(output_dir),
                "status": "FAIL_EXISTING_NONPASS_DIRECTORY",
            })
        else:
            task_id = _task_id(dataset, chrom)
            completion = None if attested_completions is None else attested_completions.get(task_id)
            if attested_completions is not None and completion is None:
                invalid.append({
                    "dataset": dataset, "chrom": chrom, "task_dir": str(output_dir),
                    "status": "FAIL_UNATTESTED_OR_ORPHAN_RANK_CHILD",
                })
                continue
            reused.append({
                "dataset": dataset,
                "chrom": chrom,
                "status": "REUSED_PASS",
                "receipt": receipt,
                "receipt_path": str(output_dir / "receipt.json"),
                **({"orchestration_completion": dict(completion)} if completion else {}),
            })
    return reused, pending, invalid


def run_specs_rolling(
    specs: Sequence[tuple[Any, ...]],
    *,
    workers: int,
    task_runner: Callable[[tuple[Any, ...]], dict[str, Any]] = run_task,
    progress: Callable[[dict[str, Any]], None] | None = None,
) -> tuple[list[dict[str, Any]], int]:
    """Maintain only ``workers`` submitted futures instead of queuing all 154."""
    results: list[dict[str, Any]] = []
    max_inflight = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        iterator = iter(specs)
        pending: dict[concurrent.futures.Future, tuple[Any, ...]] = {}

        def submit_one() -> bool:
            nonlocal max_inflight
            try:
                spec = next(iterator)
            except StopIteration:
                return False
            pending[executor.submit(task_runner, spec)] = spec
            max_inflight = max(max_inflight, len(pending))
            return True

        for _ in range(workers):
            if not submit_one():
                break
        failed = False
        while pending and not failed:
            done, _ = concurrent.futures.wait(
                pending, return_when=concurrent.futures.FIRST_COMPLETED
            )
            completed_batch: list[dict[str, Any]] = []
            for future in sorted(done, key=lambda item: specs.index(pending[item])):
                pending.pop(future)
                result = future.result()
                results.append(result)
                completed_batch.append(result)
                if progress is not None:
                    progress(result)
            failed = any(
                result["status"] not in {"PASS", "REUSED_PASS"}
                for result in completed_batch
            )
            if not failed:
                for _ in completed_batch:
                    submit_one()
        if failed:
            for future in pending:
                future.cancel()
    return results, max_inflight


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


def _empty_group() -> dict[str, Any]:
    return {
        "sums": Counter(),
        "selection_status_counts": Counter(),
        "candidate_generation_status_counts": Counter(),
        "k_route_counts": Counter(),
        "projected_call_class_counts": Counter(),
        "conditional_candidate_ranking_bootstrap_status_counts": Counter(),
        "topology_class_inclusion_counts": Counter(),
        "coarse_topology_unique_class_counts": Counter(),
        "coarse_topology_ambiguous_class_set_counts": Counter(),
        "topology_derivation_status_counts": Counter(),
        "exact_topology_uniqueness_status_counts": Counter(),
        "partial_pattern_funnel": {},
    }


def _merge_numeric_tree(target: dict[str, Any], source: Mapping[str, Any]) -> None:
    for key, value in source.items():
        if key == "definitions":
            if key in target and target[key] != value:
                raise RuntimeError("partial funnel definition drift across child receipts")
            target[key] = dict(value)
        elif isinstance(value, Mapping):
            child = target.setdefault(key, {})
            if not isinstance(child, dict):
                raise RuntimeError(f"partial funnel type drift at {key}")
            _merge_numeric_tree(child, value)
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            target[key] = target.get(key, 0) + value
        else:
            if key in target and target[key] != value:
                raise RuntimeError(f"partial funnel scalar drift at {key}")
            target[key] = value


def _add_summary(
    target: dict[str, Any], summary: Mapping[str, Any], partial_funnel: Mapping[str, Any] | None = None
) -> None:
    target["sums"].update({key: int(summary.get(key, 0)) for key in SUM_FIELDS})
    for field in (
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
    ):
        target[field].update({key: int(value) for key, value in (summary.get(field) or {}).items()})
    if partial_funnel is not None:
        _merge_numeric_tree(target["partial_pattern_funnel"], partial_funnel)


def _freeze_group(group: Mapping[str, Any]) -> dict[str, Any]:
    frozen = {
        **dict(group["sums"]),
        "selection_status_counts": dict(sorted(group["selection_status_counts"].items())),
        "candidate_generation_status_counts": dict(sorted(group["candidate_generation_status_counts"].items())),
        "k_route_counts": dict(sorted(group["k_route_counts"].items())),
        "projected_call_class_counts": dict(sorted(group["projected_call_class_counts"].items())),
        "conditional_candidate_ranking_bootstrap_status_counts": dict(
            sorted(group["conditional_candidate_ranking_bootstrap_status_counts"].items())
        ),
        "topology_class_inclusion_counts": dict(sorted(group["topology_class_inclusion_counts"].items())),
        "coarse_topology_unique_class_counts": dict(
            sorted(group["coarse_topology_unique_class_counts"].items())
        ),
        "coarse_topology_ambiguous_class_set_counts": dict(
            sorted(group["coarse_topology_ambiguous_class_set_counts"].items())
        ),
        "topology_derivation_status_counts": dict(sorted(group["topology_derivation_status_counts"].items())),
        "exact_topology_uniqueness_status_counts": dict(
            sorted(group["exact_topology_uniqueness_status_counts"].items())
        ),
        "partial_pattern_funnel": group["partial_pattern_funnel"],
    }
    n_units = frozen.get("n_component_hp_units", 0)
    selection_total = sum(frozen["selection_status_counts"].values())
    checks = {
        "selection_status_sum_equals_units": selection_total == n_units,
        "unique_tied_abstain_partition_units": (
            frozen.get("quality_primary_unique_vertex_units", 0)
            + frozen.get("quality_primary_tied_vertex_units", 0)
            + frozen.get("rank_abstain_units", 0) == n_units
        ),
        "solver_complete_plus_notrun_equals_units": (
            frozen.get("solver_complete_units", 0)
            + frozen.get("solver_incomplete_or_not_run_units", 0) == n_units
        ),
        "projection_funnel_conserved": frozen.get("molecule_component_projections", 0)
        == frozen.get("informative_scoring_molecules", 0)
        + frozen.get("all_x_excluded_molecules", 0),
        "structural_scoring_funnel_conserved": frozen.get("informative_scoring_molecules", 0)
        == frozen.get("structural_retained_molecules", 0)
        + frozen.get("below_minread_scoring_molecules", 0),
        "bq_molecules_equal_informative": frozen.get("bq_scoring_molecules", 0)
        == frozen.get("informative_scoring_molecules", 0),
        "raw_T_not_less_than_distinct_V": frozen.get("raw_tree_candidates_T_complete_units", 0)
        >= frozen.get("distinct_vertex_sets_V_complete_units", 0),
        "topology_unique_plus_multiple_equals_evaluated": (
            frozen.get("coarse_topology_class_unique_units", 0)
            + frozen.get("coarse_topology_multiple_class_units", 0)
            == frozen.get("topology_evaluated_units", 0)
        ),
        "partial_coverage_conserved_and_zero_unsatisfied": (
            frozen.get("partial_groups_covered", 0)
            + frozen.get("partial_groups_unsatisfied", 0)
            == frozen.get("partial_group_coverage_denominator", 0)
            and frozen.get("partial_groups_unsatisfied", 0) == 0
        ),
        "k_route_partition_equals_units": sum(frozen["k_route_counts"].values()) == n_units,
    }
    frozen["denominator_map"] = {
        "unit_denominator": n_units,
        "molecule_projection_denominator": frozen.get("molecule_component_projections", 0),
        "informative_molecule_denominator": frozen.get("informative_scoring_molecules", 0),
        "solver_complete_unit_denominator": frozen.get("solver_complete_units", 0),
        "topology_evaluated_unit_denominator": frozen.get("topology_evaluated_units", 0),
        "partial_group_coverage_denominator": frozen.get("partial_group_coverage_denominator", 0),
        "topology_inclusion_denominator": frozen.get("topology_class_inclusion_counts_denominator", 0),
    }
    frozen["conservation_checks"] = checks
    frozen["all_conserved"] = all(checks.values())
    funnel = frozen.get("partial_pattern_funnel") or {}
    denominator_sources = {
        "unique_patterns": funnel.get("unique_rax_pattern_groups") or {},
        "quality_groups": funnel.get("bq_quality_pattern_groups") or {},
        "molecule_projections": funnel.get("molecule_projections") or {},
    }
    frozen["solver_complete_incomplete_counts"] = {
        "COMPLETE": frozen.get("solver_complete_units", 0),
        "INCOMPLETE_OR_ABSTAIN": frozen.get("solver_incomplete_or_not_run_units", 0),
    }
    evaluated_stability = frozen.get("fixed_error_grid_evaluated_units", 0)
    stable = frozen.get("fixed_error_grid_stable_units", 0)
    frozen["likelihood_stability_counts"] = {
        "STABLE": stable,
        "UNSTABLE": max(0, evaluated_stability - stable),
        "NOT_ASSESSED": max(0, n_units - evaluated_stability),
    }
    frozen["coarse_topology_unique_units"] = frozen.get(
        "coarse_topology_class_unique_units", 0
    )
    frozen["coarse_topology_ambiguous_units"] = frozen.get(
        "coarse_topology_multiple_class_units", 0
    )
    frozen["ambiguous_topology_class_set_counts"] = frozen.get(
        "coarse_topology_ambiguous_class_set_counts", {}
    )
    frozen["partial_pattern_denominators"] = {
        key: int(value.get("denominator", 0)) for key, value in denominator_sources.items()
    }
    frozen["partial_u_distribution"] = {
        key: dict(value.get("u_number_of_X_distribution") or {})
        for key, value in denominator_sources.items()
    }
    return frozen


def aggregate_rank_receipts(results: Iterable[dict[str, Any]]) -> dict[str, Any]:
    global_groups: dict[tuple[str, str], dict[str, Any]] = defaultdict(_empty_group)
    dataset_groups: dict[str, dict[tuple[str, str], dict[str, Any]]] = defaultdict(
        lambda: defaultdict(_empty_group)
    )
    minread_global_groups: dict[str, dict[tuple[str, str], dict[str, Any]]] = defaultdict(
        lambda: defaultdict(_empty_group)
    )
    minread_dataset_groups: dict[str, dict[str, dict[tuple[str, str], dict[str, Any]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(_empty_group))
    )
    task_status = Counter()
    input_funnel = Counter()
    input_call_codes = Counter()
    dataset_input_funnel: dict[str, Counter[str]] = defaultdict(Counter)
    dataset_input_call_codes: dict[str, Counter[str]] = defaultdict(Counter)
    dataset_hp_rows: dict[str, Counter[str]] = defaultdict(Counter)
    for result in results:
        task_status[result["status"]] += 1
        receipt = result.get("receipt")
        if not receipt:
            continue
        dataset = result["dataset"]
        input_counts = receipt.get("input_counts") or {}
        input_funnel.update({key: int(value) for key, value in input_counts.items() if isinstance(value, int)})
        dataset_input_funnel[dataset].update(
            {key: int(value) for key, value in input_counts.items() if isinstance(value, int)}
        )
        input_call_codes.update({
            key: int(value) for key, value in (input_counts.get("selected_sparse_call_code_counts") or {}).items()
        })
        dataset_input_call_codes[dataset].update({
            key: int(value) for key, value in (input_counts.get("selected_sparse_call_code_counts") or {}).items()
        })
        dataset_hp_rows[dataset].update({
            key: int(value) for key, value in (input_counts.get("hp_family_rows") or {}).items()
        })
        primary_partial = receipt.get("partial_pattern_funnel_by_linkage_basis_threshold") or {}
        for basis, thresholds in (receipt.get("aggregate_by_linkage_basis_threshold") or {}).items():
            for threshold, summary in thresholds.items():
                key = (basis, str(threshold))
                funnel = (primary_partial.get(basis) or {}).get(str(threshold))
                _add_summary(global_groups[key], summary, funnel)
                _add_summary(dataset_groups[dataset][key], summary, funnel)
        for minread, payload in (receipt.get("sensitivity_by_structural_exact_pattern_minread") or {}).items():
            partial_by_basis = payload.get("partial_pattern_funnel_by_linkage_basis_threshold") or {}
            for basis, thresholds in (payload.get("by_linkage_basis_threshold") or {}).items():
                for threshold, summary in thresholds.items():
                    key = (basis, str(threshold))
                    funnel = (partial_by_basis.get(basis) or {}).get(str(threshold))
                    _add_summary(minread_global_groups[str(minread)][key], summary, funnel)
                    _add_summary(minread_dataset_groups[str(minread)][dataset][key], summary, funnel)

    def freeze_nested(groups: Mapping[tuple[str, str], Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
        nested: dict[str, dict[str, Any]] = defaultdict(dict)
        for (basis, threshold), group in sorted(groups.items()):
            nested[basis][threshold] = _freeze_group(group)
        return dict(nested)

    def freeze_dataset_input(dataset: str) -> dict[str, Any]:
        molecule_denominator = dataset_input_funnel[dataset].get("sparse_molecule_rows_total", 0)
        call_denominator = sum(dataset_input_call_codes[dataset].values())
        hp_included = dataset_input_funnel[dataset].get(
            "sparse_molecule_rows_included_in_at_least_one_selected_unit", 0
        )
        return {
            "input_call_funnel": dict(dataset_input_funnel[dataset]),
            "input_call_funnel_fraction_of_sparse_molecule_rows": {
                key: (value / molecule_denominator if molecule_denominator else None)
                for key, value in sorted(dataset_input_funnel[dataset].items())
                if "molecule_rows" in key
            },
            "input_sparse_call_code_counts": dict(sorted(dataset_input_call_codes[dataset].items())),
            "input_sparse_call_code_fraction": {
                key: (value / call_denominator if call_denominator else None)
                for key, value in sorted(dataset_input_call_codes[dataset].items())
            },
            "input_hp_family_rows": dict(sorted(dataset_hp_rows[dataset].items())),
            "denominators": {
                "sparse_molecule_rows": molecule_denominator,
                "selected_sparse_call_codes": call_denominator,
                "selected_fixed_ra_calls": (
                    dataset_input_funnel[dataset].get("selected_fixed_ra_calls_with_bq", 0)
                    + dataset_input_funnel[dataset].get("selected_fixed_ra_calls_without_bq", 0)
                ),
            },
            "sample_funnel": {
                "biological_id": BIOLOGICAL_IDS.get(dataset, dataset),
                "raw_sparse_molecules": molecule_denominator,
                "ps_known_molecules": dataset_input_funnel[dataset].get(
                    "sparse_molecule_rows_known_ps", 0
                ),
                "ps_missing_molecules": dataset_input_funnel[dataset].get(
                    "sparse_molecule_rows_missing_ps", 0
                ),
                "hp_included_molecules": hp_included,
                "hp_excluded_molecules": molecule_denominator - hp_included,
                "call_class_counts": {
                    code: dataset_input_call_codes[dataset].get(code, 0)
                    for code in ("R", "A", "O", "D", "S", "L", "X")
                },
            },
        }

    return {
        "task_status_counts": dict(sorted(task_status.items())),
        "input_call_funnel": dict(input_funnel),
        "input_sparse_call_code_counts": dict(sorted(input_call_codes.items())),
        "by_linkage_basis_threshold": freeze_nested(global_groups),
        "by_linkage_basis_and_threshold": freeze_nested(global_groups),
        "by_dataset": {
            dataset: {
                "by_linkage_basis_threshold": freeze_nested(groups),
                "by_linkage_basis_and_threshold": freeze_nested(groups),
                **freeze_dataset_input(dataset),
            }
            for dataset, groups in sorted(dataset_groups.items())
        },
        "by_structural_exact_pattern_minread": {
            minread: {
                "by_linkage_basis_threshold": freeze_nested(groups),
                "by_dataset": {
                    dataset: {"by_linkage_basis_threshold": freeze_nested(dataset_cells)}
                    for dataset, dataset_cells in sorted(minread_dataset_groups[minread].items())
                },
            }
            for minread, groups in sorted(minread_global_groups.items(), key=lambda item: int(item[0]))
        },
    }


def aggregate_primary_runtime_diagnostics(
    results: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Recompute exact primary-unit wall-time quantiles from streamed child TSVs.

    The input tables are streamed one child at a time.  Five packed
    ``array('d')`` vectors are retained (three all-primary metrics plus two
    invoked-segment tails) so exact nearest-rank p50/p95/p99 can be reported
    without materializing TSV rows.
    """
    aggregate_values = {metric: array("d") for metric in RUNTIME_METRICS}
    aggregate_invoked_values = {
        metric: array("d")
        for metric in (
            "candidate_generation_elapsed_seconds", "likelihood_fit_elapsed_seconds",
        )
    }
    n_child_files = 0
    for result in results:
        receipt = result.get("receipt") or {}
        if not child_runtime_contract_valid(receipt):
            raise RuntimeError(
                f"invalid child runtime contract: {result.get('dataset')} {result.get('chrom')}"
            )
        runtime = receipt["runtime_diagnostics"]
        identity = receipt["outputs"][runtime["per_unit_output"]]
        local_values = {metric: array("d") for metric in RUNTIME_METRICS}
        local_invoked_values = {metric: array("d") for metric in aggregate_invoked_values}
        with gzip.open(identity["path"], "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != RUNTIME_DIAGNOSTIC_FIELDS:
                raise RuntimeError("child runtime diagnostic columns drifted")
            for row in reader:
                if row["structural_minread_role"] != "PRIMARY":
                    continue
                values = {metric: float(row[metric]) for metric in RUNTIME_METRICS}
                if (
                    any(not math.isfinite(value) or value < 0.0 for value in values.values())
                    or values["candidate_generation_elapsed_seconds"]
                    + values["likelihood_fit_elapsed_seconds"]
                    > values["unit_total_elapsed_seconds"] + 1e-9
                ):
                    raise RuntimeError("invalid per-unit monotonic runtime diagnostic")
                for metric, value in values.items():
                    local_values[metric].append(value)
                invocation = {
                    "candidate_generation_elapsed_seconds": row["candidate_generation_invoked"],
                    "likelihood_fit_elapsed_seconds": row["likelihood_fit_invoked"],
                }
                if any(value not in {"true", "false"} for value in invocation.values()):
                    raise RuntimeError("invalid runtime invocation flag")
                for metric, value in invocation.items():
                    if value == "true":
                        local_invoked_values[metric].append(values[metric])
        declared = runtime["scopes"]["primary_unit_evaluations"]
        local_summary = {
            "n_unit_evaluations": len(local_values[RUNTIME_METRICS[0]]),
            **{
                metric: summarize_runtime_values(local_values[metric])
                for metric in RUNTIME_METRICS
            },
        }
        if local_summary != declared:
            raise RuntimeError("child runtime summary does not match streamed per-unit values")
        local_invoked_summary = {
            metric: summarize_runtime_values(values)
            for metric, values in local_invoked_values.items()
        }
        if local_invoked_summary != runtime["primary_invoked_segment_scopes"]:
            raise RuntimeError("child invoked-segment summary does not match per-unit flags")
        for metric in RUNTIME_METRICS:
            aggregate_values[metric].extend(local_values[metric])
        for metric in aggregate_invoked_values:
            aggregate_invoked_values[metric].extend(local_invoked_values[metric])
        n_child_files += 1
    n_units = len(aggregate_values[RUNTIME_METRICS[0]])
    if any(len(values) != n_units for values in aggregate_values.values()):
        raise AssertionError("runtime metric row counts do not align")
    return {
        "schema_name": "intersubmod.m2_full_primary_runtime_diagnostics",
        "schema_version": "1.0.0",
        "clock": "time.monotonic_ns",
        "unit": "seconds",
        "scope": "primary structural-minread unit evaluations across passing child receipts",
        "quantile_definition": (
            "exact empirical nearest-rank: rank=ceil(p*n), one-based, for p in {0.50,0.95,0.99}"
        ),
        "interpretation": (
            "process-local monotonic wall-clock performance diagnostic; environment/load dependent, "
            "not a scientific result or cross-machine reproducibility claim"
        ),
        "aggregation_memory_contract": (
            "child TSVs streamed; five packed float64 vectors retained (three all-primary metrics plus "
            "candidate/likelihood invoked-segment tails; at most 40 bytes per primary unit)"
        ),
        "n_child_runtime_files": n_child_files,
        "n_unit_evaluations": n_units,
        "metrics": {
            metric: summarize_runtime_values(values)
            for metric, values in aggregate_values.items()
        },
        "metrics_when_invoked": {
            metric: summarize_runtime_values(values)
            for metric, values in aggregate_invoked_values.items()
        },
        "all_child_summaries_recomputed_from_per_unit_rows": True,
    }


def aggregate_conservation_audit(payload: Any) -> dict[str, Any]:
    cells = []

    def visit(value: Any, path: tuple[str, ...]) -> None:
        if isinstance(value, Mapping):
            if "conservation_checks" in value:
                cells.append({
                    "path": "/".join(path),
                    "all_conserved": bool(value.get("all_conserved")),
                    "failed": sorted(
                        key for key, passed in (value.get("conservation_checks") or {}).items()
                        if not passed
                    ),
                })
            for key, child in value.items():
                visit(child, path + (str(key),))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                visit(child, path + (str(index),))

    visit(payload, ())
    return {
        "n_aggregate_cells_checked": len(cells),
        "n_failed_cells": sum(not cell["all_conserved"] for cell in cells),
        "all_aggregate_cells_conserved": bool(cells) and all(cell["all_conserved"] for cell in cells),
        "failed_cells": [cell for cell in cells if not cell["all_conserved"]],
    }


def _canonical_candidate_unit_key(
    row: Mapping[str, str], primary_minread: int
) -> str:
    chrom_number = int(str(row["chrom"]).removeprefix("chr"))
    return (
        f"{row['dataset']}|chr{chrom_number:02d}|{row['component_basis']}|PS={row['phase_set']}|"
        f"B{int(row['threshold']):03d}|{row['component_id']}|HP{row['family']}|M{primary_minread}"
    )


def _canonical_candidate_group_rows(
    unit_key: str, group: Sequence[Mapping[str, str]]
) -> list[dict[str, str]]:
    ordered_group = sorted(group, key=lambda row: row["vertex_set_id"])
    unit_optimizer_pass = all(
        row["fit_converged"].lower() == "true"
        and row["fit_monotone"].lower() == "true"
        for row in ordered_group
    )
    best_ll = max(float(row["primary_log_likelihood"]) for row in ordered_group)
    outputs: list[dict[str, str]] = []
    for candidate_index, row in enumerate(ordered_group, start=1):
        states = json.loads(row["states_json"])
        state_map = {str(item["bitmask"]): item["state_rax"] for item in states}
        role_map = {str(item["bitmask"]): item["roles"] for item in states}
        winner = row["is_winner"].lower() == "true"
        tied = row["is_tied_winner"].lower() == "true"
        winner_status = (
            "ABSTAIN_UNIT_OPTIMIZER" if not unit_optimizer_pass
            else "TIED_WINNER" if winner and tied
            else "UNIQUE_WINNER" if winner
            else "NON_WINNER"
        )
        outputs.append({
            "unit_key": unit_key,
            "dataset": row["dataset"],
            "chrom": row["chrom"],
            "component_id": row["component_id"],
            "threshold": row["threshold"],
            "hp_family": row["family"],
            "ps": row["phase_set"],
            "candidate_id": f"C{candidate_index:06d}",
            "vertex_set_id": row["vertex_set_id"],
            "vertex_states": json.dumps(state_map, sort_keys=True, separators=(",", ":")),
            "vertex_roles": json.dumps(role_map, sort_keys=True, separators=(",", ":")),
            "parent_choice_count": row["parent_choice_count"],
            "profile_log_likelihood": row["primary_log_likelihood"],
            "relative_log_likelihood": format(
                float(row["primary_log_likelihood"]) - best_ll, ".12g"
            ),
            "mixture_weights_pi": row["mixture_weights_json"],
            "winner_status": winner_status,
            "tie_group": "TOP_TIE" if tied else "",
            "coarse_topology_class": row["coarse_topology_classes_json"],
            "candidate_set_complete": "true",
        })
    return outputs


def _iter_canonical_candidate_rows(
    results: Sequence[Mapping[str, Any]],
    primary_minread: int,
    *,
    verify_child_identities: bool,
) -> Iterable[dict[str, str]]:
    """Stream the exact consolidated rows implied by child candidate artifacts."""
    previous_unit: str | None = None
    for result in results:
        receipt = result.get("receipt") or {}
        identity = (receipt.get("outputs") or {}).get(
            "m2_compressed_vertex_set_candidates.tsv.gz"
        )
        if not isinstance(identity, Mapping):
            raise RuntimeError(
                f"child candidate table missing: {result.get('dataset')} {result.get('chrom')}"
            )
        source_path = Path(str(identity.get("path", "")))
        if verify_child_identities:
            observed = os.lstat(source_path)
            if (
                not stat.S_ISREG(observed.st_mode)
                or stat.S_ISLNK(observed.st_mode)
                or observed.st_nlink != 1
                or observed.st_size != identity.get("size_bytes")
                or sha256_path(source_path) != identity.get("sha256")
            ):
                raise RuntimeError("child candidate artifact identity drifted")
        with gzip.open(source_path, "rt", encoding="utf-8", newline="") as source:
            current_unit: str | None = None
            group: list[dict[str, str]] = []
            for row in csv.DictReader(source, delimiter="\t"):
                if int(row["structural_exact_pattern_minread"]) != primary_minread:
                    continue
                if row["phase_set"] in {"", ".", "None", "NA", "UNKNOWN"}:
                    raise RuntimeError("primary child candidate table contains missing PS")
                unit_key = _canonical_candidate_unit_key(row, primary_minread)
                if current_unit is None:
                    current_unit = unit_key
                elif unit_key != current_unit:
                    if unit_key < current_unit:
                        raise RuntimeError("child candidate rows are not unit-sorted")
                    if previous_unit is not None and current_unit < previous_unit:
                        raise RuntimeError("candidate table unit sort order regression")
                    previous_unit = current_unit
                    yield from _canonical_candidate_group_rows(current_unit, group)
                    current_unit = unit_key
                    group = []
                group.append(row)
            if current_unit is not None:
                if previous_unit is not None and current_unit < previous_unit:
                    raise RuntimeError("candidate table unit sort order regression")
                previous_unit = current_unit
                yield from _canonical_candidate_group_rows(current_unit, group)


def build_full_candidate_table(
    results: Sequence[Mapping[str, Any]], output_root: Path, primary_minread: int
) -> dict[str, Any]:
    """Stream primary-minread child candidates into one canonical, sorted table."""
    output_path = output_root / "m2_ps_aware_candidate_table.tsv.gz"
    semantic = hashlib.sha256()
    n_rows = 0
    n_units = 0
    max_buffered_candidate_rows = 0
    current_unit: str | None = None
    current_group_size = 0
    with gzip.open(output_path, "wt", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(
            target, CANDIDATE_TABLE_COLUMNS, delimiter="\t", extrasaction="raise"
        )
        writer.writeheader()
        for output in _iter_canonical_candidate_rows(
            results, primary_minread, verify_child_identities=True
        ):
            writer.writerow(output)
            semantic.update(
                json.dumps(
                    output, sort_keys=True, separators=(",", ":"), ensure_ascii=False
                ).encode()
                + b"\n"
            )
            if output["unit_key"] != current_unit:
                if current_unit is not None:
                    max_buffered_candidate_rows = max(
                        max_buffered_candidate_rows, current_group_size
                    )
                current_unit = output["unit_key"]
                current_group_size = 0
                n_units += 1
            current_group_size += 1
            n_rows += 1
    if current_unit is not None:
        max_buffered_candidate_rows = max(
            max_buffered_candidate_rows, current_group_size
        )
    return {
        "schema_version": "2.0.0",
        "format": "tsv.gz",
        "columns": list(CANDIDATE_TABLE_COLUMNS),
        "sort_order": "unit_key,candidate_id",
        "path": str(output_path.resolve()),
        "size_bytes": output_path.stat().st_size,
        "sha256": sha256_path(output_path),
        "semantic_sha256": semantic.hexdigest(),
        "n_rows": n_rows,
        "n_units": n_units,
        "max_buffered_candidate_rows": max_buffered_candidate_rows,
    }


def _ranking_terminal_scope() -> dict[str, Any]:
    return {
        "datasets": list(DATASETS),
        "chromosomes": list(AUTOSOMES),
        "expected_tasks": EXPECTED_TASKS,
        "n_results": EXPECTED_TASKS,
        "biological_samples": [
            "COLO829", "H1437", "H2009", "HCC1395", "HCC1937", "HCC1954",
        ],
        "n_biological_samples": 6,
        "technical_datasets": ["HCC1395", "HCC1395_DORADO"],
        "n_technical_datasets": 7,
        "child_schema_version": "2.0.0",
        "primary_unit": "HP_family×known_PS×read_linked_component×threshold",
        "missing_ps_policy": "SEPARATE_DIAGNOSTIC_NOT_PRIMARY",
        "primary_likelihood": "BQ_SYMMETRIC_SUBSTITUTION_CONDITIONAL_RA",
    }


def _build_ranking_task_index(
    results: Sequence[Mapping[str, Any]], ranker_sha256: str
) -> list[dict[str, Any]]:
    return [
        {
            "dataset": result["dataset"],
            "chrom": result["chrom"],
            "all_pass": result["status"] in {"PASS", "REUSED_PASS"},
            "schema_version": (result.get("receipt") or {}).get("schema_version"),
            "parameters_match_extraction": result["status"] in {"PASS", "REUSED_PASS"},
            "input_hashes_match_extraction": result["status"] in {"PASS", "REUSED_PASS"},
            "upstream_outputs_verified": result["status"] in {"PASS", "REUSED_PASS"},
            "no_cross_ps_pattern_pooling": (result.get("receipt") or {}).get(
                "checks", {}
            ).get("no_cross_ps_pattern_pooling") is True,
            "known_ps_never_mixed": (result.get("receipt") or {}).get(
                "checks", {}
            ).get("known_ps_never_mixed") is True,
            "missing_ps_separate_diagnostic": (result.get("receipt") or {}).get(
                "checks", {}
            ).get("missing_ps_separate_diagnostic") is True,
            "runtime_diagnostics_contract_valid": child_runtime_contract_valid(
                result.get("receipt") or {}
            ),
            "method_contract_matches": (
                (result.get("receipt") or {}).get("parameters", {}).get("method_contract")
                == EXPECTED_METHOD_CONTRACT
            ),
            "ranker_source_bound": child_method_contract_and_ranker_source_bound(
                result.get("receipt") or {}, RANKER, ranker_sha256
            ),
        }
        for result in results
    ]


def _compact_ranking_results(
    results: Sequence[Mapping[str, Any]], output_root: Path
) -> list[dict[str, Any]]:
    return [
        {
            "dataset": result["dataset"],
            "chrom": result["chrom"],
            "status": result["status"],
            "rank_receipt": (
                None
                if not result.get("receipt")
                else {
                    "path": str(
                        output_root / "samples" / result["dataset"] / result["chrom"]
                        / "receipt.json"
                    ),
                    "sha256": sha256_path(
                        output_root / "samples" / result["dataset"] / result["chrom"]
                        / "receipt.json"
                    ),
                }
            ),
            **(
                {"orchestration_completion": result["orchestration_completion"]}
                if result.get("orchestration_completion") else {}
            ),
        }
        for result in results
    ]


def _verify_existing_candidate_table(
    metadata: Any,
    output_root: Path,
    *,
    results: Sequence[Mapping[str, Any]],
    primary_minread: int,
    require_full_scope: bool = True,
) -> bool:
    """Rebuild all canonical rows from children and compare the table lockstep."""
    expected_keys = {
        "schema_version", "format", "columns", "sort_order", "path",
        "size_bytes", "sha256", "semantic_sha256", "n_rows", "n_units",
        "max_buffered_candidate_rows",
    }
    if not isinstance(metadata, dict) or set(metadata) != expected_keys:
        return False
    expected_path = output_root / "m2_ps_aware_candidate_table.tsv.gz"
    try:
        path = Path(str(metadata["path"]))
        if path.resolve(strict=True) != expected_path.resolve(strict=True):
            return False
        observed = os.lstat(path)
        if (
            not stat.S_ISREG(observed.st_mode)
            or path.is_symlink()
            or observed.st_nlink != 1
            or metadata["schema_version"] != "2.0.0"
            or metadata["format"] != "tsv.gz"
            or metadata["columns"] != list(CANDIDATE_TABLE_COLUMNS)
            or metadata["sort_order"] != "unit_key,candidate_id"
            or metadata["size_bytes"] != observed.st_size
            or metadata["sha256"] != sha256_path(path)
        ):
            return False
        semantic = hashlib.sha256()
        n_rows = 0
        n_units = 0
        current_unit: str | None = None
        current_group_size = 0
        max_group_size = 0
        previous_sort_key: tuple[str, str] | None = None
        expected_task_ids = [
            _task_id(dataset, chrom)
            for dataset in DATASETS for chrom in AUTOSOMES
        ]
        observed_task_ids = [
            _task_id(result.get("dataset"), result.get("chrom"))
            for result in results
            if isinstance(result, Mapping)
        ]
        if require_full_scope and observed_task_ids != expected_task_ids:
            return False
        expected_rows = iter(_iter_canonical_candidate_rows(
            results, primary_minread, verify_child_identities=True
        ))
        missing = object()
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != CANDIDATE_TABLE_COLUMNS:
                return False
            for row in reader:
                if set(row) != set(CANDIDATE_TABLE_COLUMNS):
                    return False
                expected = next(expected_rows, missing)
                if expected is missing or row != expected:
                    return False
                sort_key = (row["unit_key"], row["candidate_id"])
                if previous_sort_key is not None and sort_key <= previous_sort_key:
                    return False
                previous_sort_key = sort_key
                if row["unit_key"] != current_unit:
                    if current_unit is not None:
                        max_group_size = max(max_group_size, current_group_size)
                    current_unit = row["unit_key"]
                    current_group_size = 0
                    n_units += 1
                current_group_size += 1
                n_rows += 1
                semantic.update(
                    json.dumps(
                        row, sort_keys=True, separators=(",", ":"), ensure_ascii=False
                    ).encode()
                    + b"\n"
                )
        if next(expected_rows, missing) is not missing:
            return False
        if current_unit is not None:
            max_group_size = max(max_group_size, current_group_size)
        return (
            metadata["semantic_sha256"] == semantic.hexdigest()
            and metadata["n_rows"] == n_rows
            and metadata["n_units"] == n_units
            and metadata["max_buffered_candidate_rows"] == max_group_size
        )
    except (
        OSError, KeyError, TypeError, ValueError, RuntimeError, gzip.BadGzipFile,
        json.JSONDecodeError,
    ):
        return False


def build_ranking_terminal_receipt(
    results: Sequence[dict[str, Any]],
    *,
    output_root: Path,
    full_extraction_path: Path,
    run_contract: Mapping[str, Any],
    invocation: Mapping[str, Any],
    elapsed_seconds: float,
    aggregate: Mapping[str, Any],
    conservation_audit: Mapping[str, Any],
    candidate_table: Mapping[str, Any],
    runtime_diagnostics: Mapping[str, Any],
    ranker_sha256: str,
    release_orchestration: Mapping[str, Any] | None,
) -> dict[str, Any]:
    all_children_pass = len(results) == EXPECTED_TASKS and all(
        result["status"] in {"PASS", "REUSED_PASS"} for result in results
    )
    all_child_method_contracts_identical_and_source_bound = (
        all_children_pass
        and all(
            child_method_contract_and_ranker_source_bound(
                result.get("receipt") or {}, RANKER, ranker_sha256
            )
            for result in results
        )
    )
    task_index = _build_ranking_task_index(results, ranker_sha256)
    rankable_units = sum(
        int(cell.get("solver_complete_units", 0))
        for thresholds in aggregate.get("by_linkage_basis_and_threshold", {}).values()
        for cell in thresholds.values()
    )
    primary_unit_evaluations = sum(
        int(cell.get("n_component_hp_units", 0))
        for thresholds in aggregate.get("by_linkage_basis_and_threshold", {}).values()
        for cell in thresholds.values()
    )
    checks = {
        "all_154_results_present": len(results) == EXPECTED_TASKS,
        "all_154_task_pairs_are_unique_and_exactly_canonical": (
            len(task_index) == EXPECTED_TASKS
            and {(row["dataset"], row["chrom"]) for row in task_index}
            == {(dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES}
        ),
        "all_154_rank_receipts_pass_and_match_contract": (
            len(results) == EXPECTED_TASKS
            and all(result["status"] in {"PASS", "REUSED_PASS"} for result in results)
        ),
        "all_seven_datasets_present": set(aggregate["by_dataset"]) == set(DATASETS),
        "scope_n_results_is_154": len(results) == EXPECTED_TASKS,
        "all_dataset_basis_threshold_and_global_cells_conserve": conservation_audit[
            "all_aggregate_cells_conserved"
        ],
        "all_tasks_pass": all_children_pass,
        "all_datasets_present": set(aggregate.get("by_dataset", {})) == set(DATASETS),
        "all_autosomes_present": len(task_index) == EXPECTED_TASKS
        and {row["chrom"] for row in task_index} == set(AUTOSOMES),
        "all_child_receipts_schema_2_0": len(task_index) == EXPECTED_TASKS
        and all(row["schema_version"] == "2.0.0" for row in task_index),
        "all_upstream_output_hashes_verified": len(task_index) == EXPECTED_TASKS
        and all(row["upstream_outputs_verified"] for row in task_index),
        "parameters_match_extraction": len(task_index) == EXPECTED_TASKS
        and all(row["parameters_match_extraction"] for row in task_index),
        "inputs_match_extraction_receipt": len(task_index) == EXPECTED_TASKS
        and all(row["input_hashes_match_extraction"] for row in task_index),
        "aggregate_conserved": conservation_audit["all_aggregate_cells_conserved"],
        "all_154_child_method_contracts_identical_and_source_bound": (
            all_child_method_contracts_identical_and_source_bound
        ),
        "k_gt12_never_claimed_global_optimal": all_children_pass,
        "same_read_vaf_not_double_counted": (
            all_child_method_contracts_identical_and_source_bound
            and EXPECTED_METHOD_CONTRACT["same_molecule_vaf_added_as_separate_term"] is False
        ),
        "no_cross_ps_pattern_pooling": len(task_index) == EXPECTED_TASKS
        and all(row["no_cross_ps_pattern_pooling"] for row in task_index),
        "known_ps_never_mixed": len(task_index) == EXPECTED_TASKS
        and all(row["known_ps_never_mixed"] for row in task_index),
        "missing_ps_separate_diagnostic": len(task_index) == EXPECTED_TASKS
        and all(row["missing_ps_separate_diagnostic"] for row in task_index),
        "all_child_runtime_diagnostics_contracts_valid": len(task_index) == EXPECTED_TASKS
        and all(row["runtime_diagnostics_contract_valid"] for row in task_index),
        "full_runtime_aggregate_covers_all_154_children": bool(runtime_diagnostics)
        and runtime_diagnostics.get("n_child_runtime_files") == EXPECTED_TASKS,
        "full_runtime_aggregate_covers_all_primary_unit_evaluations": bool(
            runtime_diagnostics
        ) and runtime_diagnostics.get("n_unit_evaluations") == primary_unit_evaluations,
        "full_runtime_aggregate_recomputed_from_streamed_unit_rows": bool(
            runtime_diagnostics
        ) and runtime_diagnostics.get(
            "all_child_summaries_recomputed_from_per_unit_rows"
        ) is True,
        "conditional_ra_model_not_claimed_full_generative": all_children_pass,
        "bootstrap_is_fixed_candidate_set_ranking_stability": all_children_pass,
        "topology_inclusion_counts_not_used_as_composition": all_children_pass,
        "candidate_table_hash_verified": bool(candidate_table)
        and candidate_table.get("sha256") == sha256_path(Path(candidate_table["path"])),
        "candidate_table_row_schema_complete": bool(candidate_table)
        and tuple(candidate_table.get("columns") or ()) == CANDIDATE_TABLE_COLUMNS,
        "candidate_table_covers_all_rankable_units": bool(candidate_table)
        and candidate_table.get("n_units") == rankable_units,
    }
    if run_contract.get("release_binding") is not None:
        checks.update({
            "terminal_receipt_is_bound_to_authenticated_release_contract": True,
            "upstream_extraction_uses_identical_release_contract": True,
            "runner_ranker_and_solver_match_frozen_snapshot": True,
            "frozen_scientific_and_scheduler_parameters_match": True,
        })
    receipt = {
        "schema_name": "intersubmod.m2_full_ranking_receipt",
        "schema_version": "2.0.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "scope": _ranking_terminal_scope(),
        "upstream_extraction_receipt": {
            "path": str(full_extraction_path.resolve()),
            "sha256": sha256_path(full_extraction_path),
        },
        "elapsed_seconds": elapsed_seconds,
        "run_contract": dict(run_contract),
        "invocation": dict(invocation),
        "operational_parameters_excluded_from_run_contract": ["max_new_tasks"],
        "aggregate": dict(aggregate),
        "by_dataset": dict(aggregate.get("by_dataset", {})),
        "aggregate_conservation_audit": dict(conservation_audit),
        "candidate_table": dict(candidate_table),
        "runtime_diagnostics": dict(runtime_diagnostics),
        "task_index": task_index,
        "results": _compact_ranking_results(results, output_root),
        "checks": checks,
    }
    if release_orchestration is not None:
        receipt["orchestration"] = dict(release_orchestration)
    receipt["all_pass"] = all(checks.values())
    return receipt


def _bind_terminal_results_to_reusable_rank_children(
    receipt_results: Any,
    reusable_results: Sequence[Mapping[str, Any]],
    output_root: Path,
) -> list[dict[str, Any]] | None:
    if not isinstance(receipt_results, list) or len(receipt_results) != EXPECTED_TASKS:
        return None
    expected_task_ids = [
        _task_id(dataset, chrom)
        for dataset in DATASETS for chrom in AUTOSOMES
    ]
    reusable_by_task: dict[str, Mapping[str, Any]] = {}
    for result in reusable_results:
        if not isinstance(result, Mapping):
            return None
        task_id = _task_id(result.get("dataset"), result.get("chrom"))
        if task_id in reusable_by_task:
            return None
        reusable_by_task[task_id] = result
    if list(reusable_by_task) != expected_task_ids:
        return None

    observed_task_ids: list[str] = []
    bound: list[dict[str, Any]] = []
    for row in receipt_results:
        if not isinstance(row, dict):
            return None
        task_id = _task_id(row.get("dataset"), row.get("chrom"))
        observed_task_ids.append(task_id)
        reusable = reusable_by_task.get(task_id)
        if reusable is None or row.get("status") not in {"PASS", "REUSED_PASS"}:
            return None
        expected_row = _compact_ranking_results(
            [{**dict(reusable), "status": row["status"]}], output_root
        )[0]
        if row != expected_row:
            return None
        bound.append({**dict(reusable), "status": row["status"]})
    if observed_task_ids != expected_task_ids:
        return None
    return bound


def validated_existing_final(
    path: Path,
    run_contract: Mapping[str, Any],
    *,
    reusable_results: Sequence[Mapping[str, Any]],
    output_root: Path,
    full_extraction_path: Path,
    ranker_sha256: str,
) -> dict[str, Any] | None:
    if not path.is_file() or not verify_sha256_sidecar(path):
        return None
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if (
        receipt.get("schema_name") != "intersubmod.m2_full_ranking_receipt"
        or receipt.get("all_pass") is not True
        or (receipt.get("scope") or {}).get("n_results") != EXPECTED_TASKS
        or receipt.get("run_contract") != run_contract
        or receipt.get("receipt_integrity") != _receipt_integrity(path)
    ):
        return None
    bound_results = _bind_terminal_results_to_reusable_rank_children(
        receipt.get("results"), reusable_results, output_root
    )
    invocation = receipt.get("invocation")
    elapsed_seconds = receipt.get("elapsed_seconds")
    primary_minread = (run_contract.get("parameters") or {}).get(
        "primary_structural_exact_pattern_minread"
    )
    if (
        bound_results is None
        or not isinstance(invocation, dict)
        or not isinstance(elapsed_seconds, (int, float))
        or isinstance(elapsed_seconds, bool)
        or not math.isfinite(float(elapsed_seconds))
        or float(elapsed_seconds) < 0.0
        or not isinstance(primary_minread, int)
        or isinstance(primary_minread, bool)
        or primary_minread < 1
        or not _verify_existing_candidate_table(
            receipt.get("candidate_table"),
            output_root,
            results=bound_results or (),
            primary_minread=primary_minread if isinstance(primary_minread, int) else -1,
        )
    ):
        return None
    try:
        aggregate = aggregate_rank_receipts(bound_results)
        conservation_audit = aggregate_conservation_audit(aggregate)
        runtime_diagnostics = aggregate_primary_runtime_diagnostics(bound_results)
        expected = build_ranking_terminal_receipt(
            bound_results,
            output_root=output_root,
            full_extraction_path=full_extraction_path,
            run_contract=run_contract,
            invocation=invocation,
            elapsed_seconds=elapsed_seconds,
            aggregate=aggregate,
            conservation_audit=conservation_audit,
            candidate_table=receipt["candidate_table"],
            runtime_diagnostics=runtime_diagnostics,
            ranker_sha256=ranker_sha256,
            release_orchestration=(
                receipt.get("orchestration")
                if all(result.get("orchestration_completion") for result in reusable_results)
                else None
            ),
        )
    except (OSError, KeyError, TypeError, ValueError, RuntimeError):
        return None
    expected["receipt_integrity"] = _receipt_integrity(path)
    if (
        not all(result.get("orchestration_completion") for result in reusable_results)
        and "orchestration" in receipt
    ):
        return None
    return receipt if receipt == expected else None


def conflict_kind(command: str) -> str | None:
    """Classify actual argv tokens for all known high-load extraction/ranking families."""
    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()
    basenames = [Path(token).name for token in tokens]
    if any(
        name.startswith("analyze_all_ssnv_focal_alt_multigroup") and name.endswith(".py")
        for name in basenames
    ):
        return "all_ssnv_focal_alt_multigroup"
    if any(
        name.startswith("audit_cooccurrence_task_contract_preflight")
        and name.endswith(".py")
        for name in basenames
    ):
        return "all_ssnv_cooccurrence_audit"
    if "extract_lossless_read_linkage.py" in basenames:
        return "m2_extractor"
    if "run_full_m2_extraction.py" in basenames:
        return "m2_full_extraction_runner"
    if "build_m2_patterns_and_rank.py" in basenames:
        return "m2_ranker"
    if "run_full_m2_ranking.py" in basenames:
        return "m2_full_ranking_runner"
    return None


def parse_process_table(output: str) -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    for line in output.splitlines():
        fields = line.strip().split(None, 4)
        if len(fields) != 5:
            continue
        pid, ppid, elapsed, cpu, command = fields
        try:
            rows[int(pid)] = {
                "pid": int(pid), "ppid": int(ppid), "elapsed_seconds": int(elapsed),
                "cpu_percent": float(cpu), "command": command,
            }
        except ValueError:
            continue
    return rows


def process_family(pids: dict[int, dict[str, Any]], current: int) -> set[int]:
    """Return current process, every ancestor, and every descendant."""
    ancestors = {current}
    cursor = current
    while cursor in pids:
        parent = pids[cursor]["ppid"]
        if parent <= 0 or parent in ancestors:
            break
        ancestors.add(parent)
        cursor = parent
    descendants = {current}
    changed = True
    while changed:
        changed = False
        for pid, row in pids.items():
            if row["ppid"] in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    return ancestors | descendants


def summarize_conflicts(pids: dict[int, dict[str, Any]], excluded: set[int]) -> dict[str, Any]:
    matched = {
        pid: row | {"kind": kind}
        for pid, row in pids.items()
        if pid not in excluded and (kind := conflict_kind(row["command"])) is not None
    }
    groups: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for pid, row in matched.items():
        root = pid
        cursor = row["ppid"]
        while cursor in matched:
            root = cursor
            cursor = matched[cursor]["ppid"]
        groups[root].append(row)
    representatives = []
    for root, members in sorted(groups.items()):
        representative = matched[root]
        representatives.append({
            **representative,
            "group_process_count": len(members),
            "group_cpu_percent_sum": sum(member["cpu_percent"] for member in members),
            "group_kinds": sorted({member["kind"] for member in members}),
            "member_pids": sorted(member["pid"] for member in members),
        })
    return {
        "process_count": len(matched),
        "root_count": len(representatives),
        "representatives": representatives,
    }


def active_conflicts(*, process_table: str | None = None, current_pid: int | None = None) -> dict[str, Any]:
    current = os.getpid() if current_pid is None else current_pid
    output = process_table
    if output is None:
        output = subprocess.check_output(
            ["ps", "-eo", "pid=,ppid=,etimes=,%cpu=,cmd="],
            text=True,
            errors="replace",
        )
    pids = parse_process_table(output)
    return summarize_conflicts(pids, process_family(pids, current))


def ensure_preflight_contract(
    path: Path, preflight: Mapping[str, Any], run_contract: Mapping[str, Any]
) -> None:
    """Persist initial preflight and reject source/science/timeout drift on resume."""
    if path.exists():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing.get("run_contract") != run_contract:
            raise RuntimeError("existing ranking output root has a different run contract")
        return
    path.write_text(json.dumps(preflight, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def build_ranking_checkpoint(
    results: Sequence[dict[str, Any]],
    *,
    run_contract: Mapping[str, Any],
    invocation: Mapping[str, Any],
    elapsed_seconds: float,
) -> dict[str, Any]:
    ordered = sorted(
        results,
        key=lambda item: (DATASETS.index(item["dataset"]), AUTOSOMES.index(item["chrom"])),
    )
    if not ordered or any(
        result["status"] not in {"PASS", "REUSED_PASS"} for result in ordered
    ):
        raise ValueError("checkpoint requires at least one passing child and no failed child")
    aggregate = aggregate_rank_receipts(ordered)
    conservation = aggregate_conservation_audit(aggregate)
    passing = len(ordered)
    task_keys = [(result["dataset"], result["chrom"]) for result in ordered]
    expected_task_keys = {(dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES}
    checks = {
        "all_recorded_tasks_pass_or_reused": True,
        "passing_count_matches_results": passing == len(ordered),
        "remaining_count_conserved": passing + (EXPECTED_TASKS - passing) == EXPECTED_TASKS,
        "recorded_task_pairs_are_unique_and_canonical": (
            len(task_keys) == len(set(task_keys))
            and set(task_keys).issubset(expected_task_keys)
        ),
        "terminal_full_scope_complete": False,
        "recorded_aggregate_cells_conserve": conservation["all_aggregate_cells_conserved"],
    }
    if run_contract.get("release_binding") is not None:
        checks["checkpoint_is_bound_to_authenticated_release_contract"] = True
    return {
        "schema_name": "intersubmod.m2_full_ranking_checkpoint",
        "schema_version": "2.0.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "scope": {
            "datasets": list(DATASETS),
            "chromosomes": list(AUTOSOMES),
            "expected_tasks": EXPECTED_TASKS,
            "n_results": passing,
            "child_schema_version": "2.0.0",
        },
        "elapsed_seconds": elapsed_seconds,
        "run_contract": dict(run_contract),
        "invocation": dict(invocation),
        "operational_parameters_excluded_from_run_contract": ["max_new_tasks"],
        "task_status_counts": dict(Counter(result["status"] for result in ordered)),
        "passing_tasks": passing,
        "remaining_tasks": EXPECTED_TASKS - passing,
        "aggregate": aggregate,
        "by_dataset": aggregate.get("by_dataset", {}),
        "aggregate_conservation_audit": conservation,
        "results": [
            {
                "dataset": result["dataset"],
                "chrom": result["chrom"],
                "status": result["status"],
                "rank_receipt": {
                    "path": result.get("receipt_path"),
                    "sha256": (
                        sha256_path(Path(result["receipt_path"]))
                        if result.get("receipt_path") else None
                    ),
                },
                **(
                    {"orchestration_completion": result["orchestration_completion"]}
                    if result.get("orchestration_completion") else {}
                ),
            }
            for result in ordered
        ],
        "checks": checks,
        "checkpoint_complete": all(
            value for key, value in checks.items() if key != "terminal_full_scope_complete"
        ),
        "all_pass": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extraction-root", required=True, type=Path)
    parser.add_argument("--release-contract-manifest", type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--thresholds", default="1,2,3,5")
    parser.add_argument("--component-bases", default="PS_HP1,PS_HP2")
    parser.add_argument("--hp-families", default="1,2")
    parser.add_argument("--structural-exact-pattern-minread-grid", default="1,2,3,5")
    parser.add_argument("--primary-structural-exact-pattern-minread", type=int, default=3)
    parser.add_argument("--exact-k-max", type=int, default=12)
    parser.add_argument("--max-vertex-sets", type=int, default=256)
    parser.add_argument("--solver-time-limit-seconds", type=float, default=30.0)
    parser.add_argument("--fixed-error-grid", default="0.005,0.01,0.02,0.05")
    parser.add_argument("--minimum-bq-error-rate", type=float, default=1e-6)
    parser.add_argument("--maximum-bq-error-rate", type=float, default=0.25)
    parser.add_argument("--conditional-candidate-ranking-bootstrap-replicates", type=int, default=20)
    parser.add_argument("--conditional-candidate-ranking-bootstrap-seed", type=int, default=20260716)
    parser.add_argument("--tie-tolerance", type=float, default=1e-6)
    parser.add_argument("--max-new-tasks", type=int)
    parser.add_argument("--task-timeout-seconds", type=float, default=28800.0)
    parser.add_argument(
        "--task-timeout-grace-seconds", "--timeout-grace-seconds",
        dest="timeout_grace_seconds", type=float, default=300.0,
    )
    parser.add_argument("--aggregate-only", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.thresholds = tuple(sorted({int(value) for value in args.thresholds.split(",") if value}))
    args.component_bases = tuple(value for value in args.component_bases.split(",") if value)
    args.hp_families = tuple(value for value in args.hp_families.split(",") if value)
    args.structural_exact_pattern_minreads = tuple(sorted({
        int(value) for value in args.structural_exact_pattern_minread_grid.split(",") if value
    }))
    args.fixed_error_grid_values = tuple(
        float(value) for value in args.fixed_error_grid.split(",") if value
    )
    if (
        not 1 <= args.workers <= 4
        or not args.thresholds
        or min(args.thresholds) < 1
        or not args.component_bases
        or not args.hp_families
        or args.conditional_candidate_ranking_bootstrap_replicates < 0
        or not args.structural_exact_pattern_minreads
        or min(args.structural_exact_pattern_minreads) < 1
        or args.primary_structural_exact_pattern_minread not in args.structural_exact_pattern_minreads
        or tuple(args.component_bases) != ("PS_HP1", "PS_HP2")
        or tuple(args.hp_families) != ("1", "2")
        or not args.fixed_error_grid_values
        or any(not 0 < value < 0.5 for value in args.fixed_error_grid_values)
        or not 0 < args.minimum_bq_error_rate <= args.maximum_bq_error_rate < 0.5
        or (args.max_new_tasks is not None and args.max_new_tasks < 1)
        or args.task_timeout_seconds <= 0
        or args.timeout_grace_seconds < 0
        or not math.isfinite(args.task_timeout_seconds)
        or not math.isfinite(args.timeout_grace_seconds)
    ):
        raise ValueError("invalid full-ranking worker, basis, threshold, bootstrap, or likelihood contract")
    full_extraction, full_extraction_path = load_full_extraction_receipt(args.extraction_root)
    release_binding = None
    if args.release_contract_manifest is not None:
        release_binding = load_release_contract_binding(
            args.release_contract_manifest,
            required_sources={
                "full_ranking_runner": RUNNER,
                "ranker": RANKER,
                "hypercube_solver": HYPERCUBE_SOLVER,
            },
        )
        validate_release_ranking_parameters(release_binding, args)
        require_matching_upstream_release_binding(full_extraction, release_binding)
    ranker_sha256 = sha256_path(RANKER)
    hypercube_solver_sha256 = sha256_path(HYPERCUBE_SOLVER)
    runner_sha256 = sha256_path(RUNNER)
    expected_parameters = _rank_parameters(args)
    run_contract = {
        "full_extraction_receipt": {
            "path": str(full_extraction_path.resolve()),
            "sha256": sha256_path(full_extraction_path),
        },
        "ranker": {"path": str(RANKER.resolve()), "sha256": ranker_sha256},
        "hypercube_solver": {
            "path": str(HYPERCUBE_SOLVER.resolve()), "sha256": hypercube_solver_sha256,
        },
        "runner": {"path": str(RUNNER), "sha256": runner_sha256},
        "method_contract": EXPECTED_METHOD_CONTRACT,
        "thresholds": list(args.thresholds),
        "component_bases": list(args.component_bases),
        "hp_families": list(args.hp_families),
        "parameters": expected_parameters,
        "workers": args.workers,
        "task_timeout_seconds": args.task_timeout_seconds,
        "timeout_grace_seconds": args.timeout_grace_seconds,
    }
    if release_binding is not None:
        run_contract["release_binding"] = release_binding
        run_contract["orchestration_policy"] = dict(ORCHESTRATION_POLICY)
    invocation = {
        "max_new_tasks": args.max_new_tasks,
        "pending_selection_policy": "canonical DATASETS order then chr1-22; first N pending",
        "task_timeout_seconds": args.task_timeout_seconds,
        "timeout_grace_seconds": args.timeout_grace_seconds,
        "workers": args.workers,
        "aggregate_only": args.aggregate_only,
    }
    conflicts = active_conflicts()
    resource_preview = resource_gate_preview(args.output_root, conflicts)
    preflight = {
        "schema_name": "intersubmod.m2_full_ranking_preflight",
        "schema_version": "2.0.0",
        "scope": {
            "datasets": list(DATASETS),
            "chromosomes": list(AUTOSOMES),
            "expected_tasks": 154,
        },
        "upstream_full_extraction_all_pass": full_extraction.get("all_pass") is True,
        "runner": {"path": str(RUNNER), "sha256": runner_sha256},
        "run_contract": run_contract,
        "release_binding": release_binding,
        "invocation": invocation,
        "operational_parameters_excluded_from_run_contract": ["max_new_tasks"],
        "active_conflict_process_count": conflicts["process_count"],
        "active_conflict_root_count": conflicts["root_count"],
        "active_conflicts": conflicts["representatives"],
        "resource_gate_preview": resource_preview,
        "resource_gate_pass": (
            conflicts["process_count"] == 0 and resource_preview["disk_pass"]
        ),
    }
    print(json.dumps(preflight, ensure_ascii=False, indent=2), flush=True)
    if args.preflight_only:
        existing_preflight = args.output_root / "preflight.json"
        if existing_preflight.exists():
            existing = json.loads(existing_preflight.read_text(encoding="utf-8"))
            if existing.get("run_contract") != run_contract:
                raise RuntimeError("existing ranking output root has a different run contract")
        raise SystemExit(0 if preflight["resource_gate_pass"] else 2)
    if not preflight["resource_gate_pass"]:
        raise SystemExit(2)
    release_session = None
    orchestration_state = None
    parent_extraction = None
    if release_binding is not None:
        parent_extraction = extraction_parent_identity(
            args.extraction_root, full_extraction_path
        )
        if not args.output_root.exists():
            args.output_root.parent.mkdir(parents=True, exist_ok=True)
            args.output_root.mkdir()
        root_stat = os.lstat(args.output_root)
        session_target = {
            "output_root": {
                "path": str(args.output_root.resolve()),
                "st_dev": int(root_stat.st_dev),
                "st_ino": int(root_stat.st_ino),
            },
            "release_manifest": dict(release_binding["release_manifest"]),
            "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
            "parent_extraction": dict(parent_extraction),
        }
        session_gate_identity = None
        session_path = args.output_root / "_orchestration" / "session.json"
        if not session_path.exists():
            _, session_gate_identity = create_resource_gate_receipt(
                args.output_root,
                stage="ranking",
                gate_scope="session",
                target=session_target,
                producer_source=release_producer_sources(release_binding)["runner"],
                conflicts=active_conflicts(),
            )
        release_session = ensure_release_orchestration_session(
            args.output_root, release_binding, run_contract, parent_extraction,
            session_gate_identity,
        )
        orchestration_state = load_release_orchestration_state(
            args.output_root, release_session
        )
    else:
        args.output_root.mkdir(parents=True, exist_ok=True)
    contract_path = args.output_root / "preflight.json"
    ensure_preflight_contract(contract_path, preflight, run_contract)

    expected_extraction_receipts = {
        (result.get("dataset"), result.get("chrom")): result.get("receipt")
        for result in (full_extraction.get("results") or [])
        if result.get("status") in {"PASS", "REUSED_PASS"} and result.get("receipt")
    }
    if set(expected_extraction_receipts) != {
        (dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES
    }:
        raise RuntimeError("full extraction receipt does not contain 154 passing child receipts")
    specs = []
    for dataset in DATASETS:
        for chrom in AUTOSOMES:
            extraction_dir = args.extraction_root / "samples" / dataset / chrom
            output_dir = args.output_root / "samples" / dataset / chrom
            specs.append(
                (
                    dataset,
                    chrom,
                    extraction_dir,
                    output_dir,
                    task_command(extraction_dir, output_dir, args),
                    RANKER,
                    ranker_sha256,
                    expected_parameters,
                    args.component_bases,
                    args.thresholds,
                    args.hp_families,
                    "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
                    semantic_json_sha256(expected_extraction_receipts[(dataset, chrom)]),
                    args.aggregate_only,
                    args.task_timeout_seconds,
                    args.timeout_grace_seconds,
                )
            )
    reused, pending_specs, invalid = scan_existing_rank_specs(
        specs,
        attested_completions=(
            None if orchestration_state is None else orchestration_state["completions"]
        ),
    )
    if invalid:
        print(json.dumps({"all_pass": False, "existing_nonpass": invalid}, ensure_ascii=False))
        raise SystemExit(1)
    if orchestration_state is not None and orchestration_state["count"] != len(reused):
        raise RuntimeError("attested ranking task count differs from reusable children")

    final_path = args.output_root / "full_ranking_receipt.json"
    existing_final = validated_existing_final(
        final_path,
        run_contract,
        reusable_results=reused,
        output_root=args.output_root,
        full_extraction_path=full_extraction_path,
        ranker_sha256=ranker_sha256,
    )
    if not pending_specs and existing_final is not None:
        print(json.dumps({"receipt": str(final_path), "all_pass": True, "status": "REUSED_FINAL"}))
        raise SystemExit(0)
    if final_path.exists():
        raise RuntimeError("existing terminal full ranking receipt is invalid or scope is no longer complete")

    selected = pending_specs[
        : args.max_new_tasks if args.max_new_tasks is not None else len(pending_specs)
    ]
    release_batch = None
    release_grants = None
    if release_session is not None:
        selected_ids = [_task_id(spec[0], spec[1]) for spec in selected]
        batch_target = {
            "output_root": dict(release_session["output_root"]),
            "session_id": release_session["session_id"],
            "session_sha256": sha256_path(
                args.output_root / "_orchestration" / "session.json"
            ),
            "batch_index": orchestration_state["next_batch_index"],
            "before_count": len(reused),
            "max_new_tasks": int(args.max_new_tasks),
            "effective_count": len(selected),
            "selected_task_ids": selected_ids,
            "previous_chain_head": orchestration_state["previous_chain_head"],
        }
        _, batch_gate_identity = create_resource_gate_receipt(
            args.output_root,
            stage="ranking",
            gate_scope="batch",
            batch_index=orchestration_state["next_batch_index"],
            target=batch_target,
            producer_source=release_session["producer_sources"]["runner"],
            conflicts=active_conflicts(),
        )
        release_batch, release_grants = create_batch_start_and_grants(
            args.output_root, release_session, run_contract, selected,
            before_count=len(reused),
            previous_chain_head=orchestration_state["previous_chain_head"],
            batch_index=orchestration_state["next_batch_index"],
            max_new_tasks=int(args.max_new_tasks),
            gate=batch_gate_identity,
        )
    started = time.monotonic()
    new_results: list[dict[str, Any]] = []

    def collect_progress(result: dict[str, Any]) -> None:
        new_results.append(result)
        passing = len(reused) + sum(
            item["status"] in {"PASS", "REUSED_PASS"} for item in new_results
        )
        print(
            f"M2_RANK_PROGRESS task={result['dataset']}:{result['chrom']} status={result['status']} "
            f"passing={passing}/{EXPECTED_TASKS} elapsed_seconds={time.monotonic() - started:.1f}",
            flush=True,
        )

    executed, max_inflight = run_specs_rolling(
        selected, workers=args.workers, progress=collect_progress
    )
    if executed != new_results:
        raise AssertionError("rolling execution result accounting drifted")
    if any(result["status"] not in {"PASS", "REUSED_PASS"} for result in executed):
        print(json.dumps({
            "all_pass": False,
            "status": "TASK_FAILURE",
            "results": executed,
            "max_inflight_futures": max_inflight,
        }, ensure_ascii=False))
        raise SystemExit(1)
    executed = canonical_sort_results(executed)
    if release_session is not None:
        for result in executed:
            task_id = _task_id(result["dataset"], result["chrom"])
            result["orchestration_completion"] = write_child_completion(
                args.output_root, release_session, release_batch,
                release_grants[task_id], result,
            )
    if sha256_path(RUNNER) != runner_sha256:
        raise RuntimeError("M2 ranking runner source changed during the invocation")
    results = reused + executed
    invocation = invocation | {
        "reused_tasks": len(reused),
        "pending_tasks_before_invocation": len(pending_specs),
        "selected_new_tasks": len(selected),
        "selected_task_ids": [f"{spec[0]}:{spec[1]}" for spec in selected],
        "max_inflight_futures": max_inflight,
    }
    results.sort(key=lambda item: (DATASETS.index(item["dataset"]), AUTOSOMES.index(item["chrom"])))
    release_orchestration = None
    if release_session is not None:
        release_orchestration = {
            "session_identity": _session_identity(args.output_root, release_session),
            "batch_start_identity": {
                "path": release_batch["path"], "sha256": release_batch["sha256"],
                "batch_id": release_batch["batch_id"], "batch_index": release_batch["batch_index"],
            },
            "previous_chain_head": orchestration_state["previous_chain_head"],
            "batch_completion_attestations": [
                {"task_id": _task_id(result["dataset"], result["chrom"]),
                 **result["orchestration_completion"]}
                for result in executed
            ],
            "cumulative_attested_tasks": len(results),
        }
    if release_binding is not None:
        current_binding = load_release_contract_binding(
            args.release_contract_manifest,
            required_sources={
                "full_ranking_runner": RUNNER,
                "ranker": RANKER,
                "hypercube_solver": HYPERCUBE_SOLVER,
            },
        )
        if current_binding != release_binding:
            raise RuntimeError("release contract or a frozen source changed during ranking")
        upstream_binding = (
            json.loads(full_extraction_path.read_text(encoding="utf-8")).get("run_contract") or {}
        ).get("release_binding")
        if upstream_binding != release_binding:
            raise RuntimeError("upstream extraction release binding changed during ranking")
    aggregate = aggregate_rank_receipts(results)
    conservation_audit = aggregate_conservation_audit(aggregate)
    if len(results) < EXPECTED_TASKS:
        receipt = build_ranking_checkpoint(
            results,
            run_contract=run_contract,
            invocation=invocation,
            elapsed_seconds=time.monotonic() - started,
        )
        passing = receipt["passing_tasks"]
        receipt_path = (
            args.output_root / "checkpoints" /
            f"checkpoint_{passing:03d}_of_{EXPECTED_TASKS}.json"
        )
        receipt["receipt_integrity"] = {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{receipt_path.name}.sha256",
            "covers": receipt_path.name,
        }
        if release_orchestration is not None:
            receipt["orchestration"] = release_orchestration
        if receipt["checkpoint_complete"] is not True:
            print(json.dumps({
                "all_pass": False,
                "status": "CHECKPOINT_CONSERVATION_FAILURE",
                "checks": receipt["checks"],
            }, ensure_ascii=False))
            raise SystemExit(1)
        if (
            sha256_path(RUNNER) != runner_sha256
            or sha256_path(RANKER) != ranker_sha256
            or sha256_path(HYPERCUBE_SOLVER) != hypercube_solver_sha256
            or sha256_path(full_extraction_path) != run_contract["full_extraction_receipt"]["sha256"]
        ):
            raise RuntimeError("M2 ranking code or upstream receipt changed before checkpoint publication")
        if release_binding is not None:
            current_binding = load_release_contract_binding(
                args.release_contract_manifest,
                required_sources={
                    "full_ranking_runner": RUNNER,
                    "ranker": RANKER,
                    "hypercube_solver": HYPERCUBE_SOLVER,
                },
                _force_deep_reverification=True,
            )
            if current_binding != release_binding:
                raise RuntimeError("release contract changed before ranking checkpoint publication")
            require_matching_upstream_release_binding(
                json.loads(full_extraction_path.read_text(encoding="utf-8")), release_binding
            )
        if release_session is not None:
            write_immutable_json_exclusive(receipt_path, receipt)
        else:
            write_json_and_sha256_exclusive(receipt_path, receipt)
        print(json.dumps({
            "checkpoint": str(receipt_path),
            "all_pass": False,
            "checkpoint_complete": True,
            "passing": passing,
            "remaining": receipt["remaining_tasks"],
        }, ensure_ascii=False))
        raise SystemExit(3)
    all_children_pass = len(results) == EXPECTED_TASKS and all(
        result["status"] in {"PASS", "REUSED_PASS"} for result in results
    )
    candidate_table = (
        build_full_candidate_table(
            results, args.output_root, args.primary_structural_exact_pattern_minread
        )
        if all_children_pass
        else {}
    )
    runtime_diagnostics = (
        aggregate_primary_runtime_diagnostics(results) if all_children_pass else {}
    )
    receipt = build_ranking_terminal_receipt(
        results,
        output_root=args.output_root,
        full_extraction_path=full_extraction_path,
        run_contract=run_contract,
        invocation=invocation,
        elapsed_seconds=time.monotonic() - started,
        aggregate=aggregate,
        conservation_audit=conservation_audit,
        candidate_table=candidate_table,
        runtime_diagnostics=runtime_diagnostics,
        ranker_sha256=ranker_sha256,
        release_orchestration=release_orchestration,
    )
    receipt_path = final_path
    receipt["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": f"{receipt_path.name}.sha256",
        "covers": receipt_path.name,
    }
    if receipt["all_pass"] is not True:
        print(json.dumps({
            "all_pass": False,
            "status": "TERMINAL_AGGREGATE_CHECK_FAILURE",
            "checks": receipt["checks"],
        }, ensure_ascii=False))
        raise SystemExit(1)
    if (
        sha256_path(RUNNER) != runner_sha256
        or sha256_path(RANKER) != ranker_sha256
        or sha256_path(HYPERCUBE_SOLVER) != hypercube_solver_sha256
        or sha256_path(full_extraction_path) != run_contract["full_extraction_receipt"]["sha256"]
    ):
        raise RuntimeError("M2 ranking code or upstream receipt changed before final publication")
    if release_binding is not None:
        current_binding = load_release_contract_binding(
            args.release_contract_manifest,
            required_sources={
                "full_ranking_runner": RUNNER,
                "ranker": RANKER,
                "hypercube_solver": HYPERCUBE_SOLVER,
            },
            _force_deep_reverification=True,
        )
        if current_binding != release_binding:
            raise RuntimeError("release contract changed before final ranking publication")
        require_matching_upstream_release_binding(
            json.loads(full_extraction_path.read_text(encoding="utf-8")), release_binding
        )
    if release_session is not None:
        write_immutable_json_exclusive(receipt_path, receipt)
    else:
        write_json_and_sha256_exclusive(receipt_path, receipt)
    print(json.dumps({"receipt": str(receipt_path), "all_pass": receipt["all_pass"]}, ensure_ascii=False))
    raise SystemExit(0)


if __name__ == "__main__":
    main()
