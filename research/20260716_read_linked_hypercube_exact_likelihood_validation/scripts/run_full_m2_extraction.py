#!/usr/bin/env python3
"""Resumable chr1-22 × 7-dataset M2 extraction orchestrator."""

from __future__ import annotations

import argparse
import concurrent.futures
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
from typing import Any, Callable, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
RUNNER = Path(__file__).resolve()
EXTRACTOR = SCRIPT_DIR / "extract_lossless_read_linkage.py"
LOSSLESS_READ_CONTRACT = SCRIPT_DIR / "lossless_read_contract.py"
AUTOSOMES = tuple(f"chr{value}" for value in range(1, 23))
DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
CONFLICT_PATTERNS = (
    # Match both the canonical script name and source-locked copies such as
    # ``analyze_all_ssnv_focal_alt_multigroup.pinned_<sha>.py``.
    "analyze_all_ssnv_focal_alt_multigroup",
    "audit_cooccurrence_task_contract_preflight",
    "extract_lossless_read_linkage.py",
    "run_full_m2_extraction.py",
    "build_m2_patterns_and_rank.py",
    "run_full_m2_ranking.py",
)
LINKAGE_BASES = (
    "pooled", "HP1", "HP2", "HP3", "HP4", "unphased",
    "PS_HP1", "PS_HP2", "MISSING_PS_HP1", "MISSING_PS_HP2",
)
TASK_STATUS_NAME = "runner_task_status.json"
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
    "manifest_sha_matches_selected_authority",
    "exact_eleven_file_allowlist_frozen",
    "all_sources_regular_non_symlink_single_link_and_non_aliased",
    "all_copies_regular_single_link_non_aliased_and_sha_equal_source",
    "canonical_manifest_and_pre_receipt_copied_immutably",
    "publish_boundary_recheck_completed",
    "bootstrap20_and_seed20260716_frozen",
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


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def semantic_json_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _release_require_exact_keys(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        raise RuntimeError(f"{label} exact-key schema mismatch")
    return value


def _release_validate_copy(
    record: Any, physical: Path, declared_path: str, label: str
) -> dict[str, Any]:
    row = _release_require_exact_keys(record, RELEASE_FILE_IDENTITY_KEYS, label)
    raw = physical.absolute()
    observed = os.lstat(raw)
    if (
        not stat.S_ISREG(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or stat.S_IMODE(observed.st_mode) != 0o444
        or observed.st_nlink != 1
    ):
        raise RuntimeError(f"{label} is not an immutable single-link 0444 regular file")
    expected = {
        "path": declared_path,
        "st_dev": int(observed.st_dev), "st_ino": int(observed.st_ino),
        "st_nlink": int(observed.st_nlink), "size_bytes": int(observed.st_size),
        "mtime_ns": int(observed.st_mtime_ns), "ctime_ns": int(observed.st_ctime_ns),
        "mode_octal": "0444", "sha256": sha256_path(raw),
    }
    if row != expected:
        raise RuntimeError(f"{label} identity/stat/SHA mismatch")
    return row


def _release_validate_runtime(document: Mapping[str, Any]) -> None:
    runtime = _release_require_exact_keys(
        document.get("runtime"), {"python", "packages", "samtools", "platform"}, "runtime"
    )
    python = _release_require_exact_keys(
        runtime["python"], {"executable", "version", "implementation"}, "runtime.python"
    )
    packages = _release_require_exact_keys(
        runtime["packages"], {"numpy", "scipy", "pysam"}, "runtime.packages"
    )
    samtools = _release_require_exact_keys(
        runtime["samtools"], {"path", "version_line", "htslib_version_line"}, "runtime.samtools"
    )
    if not all(
        isinstance(value, str) and bool(value)
        for value in [*python.values(), *packages.values(), *samtools.values(), runtime["platform"]]
    ):
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
    anchor_relative = Path(RELEASE_SOURCE_PATHS["full_extraction_runner"])
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
    manifest_identity = verification.get("run_manifest") or {}
    if (
        verification.get("all_pass") is not True
        or verification.get("validation_evidence_eligible") is not True
        or manifest_identity.get("sha256") != release_sha256
    ):
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
    manifest_path: Path,
    *,
    required_sources: dict[str, Path],
    input_manifest: Path,
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
    scope = document["scope"]
    expected_scope = {
        "technical_datasets": 7, "biological_samples": 6, "chromosomes": 22,
        "chromosome_names": list(AUTOSOMES), "tasks": EXPECTED_TASKS,
        "datasets": list(DATASETS),
    }
    if (
        document["schema_name"] != RELEASE_SCHEMA_NAME
        or document["schema_version"] != RELEASE_SCHEMA_VERSION
        or document["task_type"] != "B_COMPREHENSIVE_VALIDATION"
        or document["authority_mode"] != "CANONICAL_V5_FROZEN"
        or document["all_pass"] is not True
        or document["validation_evidence_eligible"] is not True
        or not isinstance(document["created_at_utc"], str)
        or not document["created_at_utc"].endswith("Z")
        or scope != expected_scope
    ):
        raise RuntimeError("release canonical authority/task/scope metadata mismatch")
    integrity = _release_require_exact_keys(
        document["receipt_integrity"], {"scheme", "sidecar_name", "covers"},
        "release receipt_integrity",
    )
    if integrity != {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": f"{release_path.name}.sha256", "covers": release_path.name,
    }:
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
    snapshot = _release_require_exact_keys(
        document["source_snapshot"],
        {"repo_root", "snapshot_root", "n_files", "entries", "entrypoints", "exact_allowlist_semantic_sha256"},
        "source_snapshot",
    )
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
        for role, relative in RELEASE_SOURCE_PATHS.items()
        if role != "canonical_manifest_schema"
    }
    entrypoints = snapshot["entrypoints"]
    if entrypoints != expected_code_entrypoints:
        raise RuntimeError("release snapshot exact entrypoint map mismatch")
    identities: dict[str, dict[str, str]] = {}
    allowlist_rows: list[dict[str, str]] = []
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
        allowlist_rows.append({"role": role, "repo_relative_path": repo_relative})
    if set(identities) != set(RELEASE_SOURCE_PATHS):
        raise RuntimeError("release snapshot role set is not the exact canonical eleven")
    expected_allowlist = [
        {"role": role, "repo_relative_path": relative}
        for role, relative in RELEASE_SOURCE_PATHS.items()
    ]
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

    canonical = _release_require_exact_keys(
        document["canonical_manifest"], {"expected_sha256", "origin", "immutable_copy"},
        "canonical_manifest",
    )
    if canonical["expected_sha256"] != CANONICAL_MANIFEST_SHA256:
        raise RuntimeError("release canonical manifest authority SHA mismatch")
    _release_require_exact_keys(canonical["origin"], RELEASE_FILE_IDENTITY_KEYS, "canonical_manifest.origin")
    relative = Path("input_contract/canonical_manifest.json")
    physical_copy_raw = release_path.parent / relative
    physical_copy = physical_copy_raw.resolve(strict=True)
    input_root = release_path.parent / "input_contract"
    if input_root.is_symlink() or not physical_copy.is_relative_to(input_root.resolve(strict=True)):
        raise RuntimeError("immutable canonical manifest escapes input_contract")
    copied_manifest = _release_validate_copy(
        canonical["immutable_copy"], physical_copy_raw, relative.as_posix(), "canonical_manifest.copy"
    )
    copy_identity = {"path": str(physical_copy), "sha256": copied_manifest["sha256"]}
    if copy_identity["sha256"] != CANONICAL_MANIFEST_SHA256:
        raise RuntimeError("immutable canonical manifest copy SHA mismatch")
    schema = _release_require_exact_keys(
        document["canonical_schema"], {"role", "repo_relative_path", "sha256"}, "canonical_schema"
    )
    if schema != {
        "role": "canonical_manifest_schema", "repo_relative_path": CANONICAL_SCHEMA_RELATIVE,
        "sha256": CANONICAL_SCHEMA_SHA256,
    }:
        raise RuntimeError("canonical schema authority mismatch")
    pre = _release_require_exact_keys(
        document["pre_input_identity_receipt"],
        {"origin", "immutable_copy", "input_identity_snapshot_sha256", "receipt_semantic_sha256", "authority_mode", "validation_evidence_eligible", "artifact_roles"},
        "pre_input_identity_receipt",
    )
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
    fresh_keys = {
        "execution_mode", "verifier_path", "verifier_sha256", "receipt_sha256",
        "receipt_semantic_sha256", "checks_semantic_sha256", "artifact_audit_semantic_sha256",
        "input_identity_snapshot_sha256", "all_pass", "validation_evidence_eligible",
        "exactly_equals_supplied_pre_snapshot",
    }
    fresh = _release_require_exact_keys(document["fresh_input_identity_verification"], fresh_keys, "fresh_input_identity_verification")
    if (
        fresh["execution_mode"] != "production_subprocess_required"
        or fresh["all_pass"] is not True
        or fresh["validation_evidence_eligible"] is not True
        or fresh["exactly_equals_supplied_pre_snapshot"] is not True
        or fresh["input_identity_snapshot_sha256"] != pre["input_identity_snapshot_sha256"]
        or fresh["verifier_path"] != producer_records["input_identity_verifier"]["source"]["path"]
        or fresh["verifier_sha256"] != producer_records["input_identity_verifier"]["source"]["sha256"]
    ):
        raise RuntimeError("fresh input identity verification summary mismatch")
    producer = _release_require_exact_keys(
        document["producer"], {"role", "repo_relative_path", "source_sha256", "immutable_copy_path", "immutable_copy_sha256"}, "producer"
    )
    freezer_entry = producer_records["release_contract_freezer"]
    if producer != {
        "role": "release_contract_freezer",
        "repo_relative_path": RELEASE_SOURCE_PATHS["release_contract_freezer"],
        "source_sha256": freezer_entry["source"]["sha256"],
        "immutable_copy_path": freezer_entry["snapshot"]["path"],
        "immutable_copy_sha256": freezer_entry["snapshot"]["sha256"],
    }:
        raise RuntimeError("release producer is not bound to frozen freezer source")
    expected_top_entrypoints = {
        **expected_code_entrypoints,
        "canonical_manifest_copy": "input_contract/canonical_manifest.json",
        "pre_input_identity_receipt_copy": pre_relative.as_posix(),
        "canonical_schema_copy": (Path("code_snapshot") / CANONICAL_SCHEMA_RELATIVE).as_posix(),
    }
    if document["entrypoints"] != expected_top_entrypoints:
        raise RuntimeError("release top-level exact entrypoint map mismatch")
    actual_manifest = input_manifest.resolve(strict=True)
    if (
        str(actual_manifest) != copy_identity["path"]
        or sha256_path(actual_manifest) != copy_identity["sha256"]
        or copy_identity["sha256"] != CANONICAL_MANIFEST_SHA256
    ):
        raise RuntimeError("runner --manifest is not the contract-bound canonical manifest copy")
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


def validate_release_extraction_parameters(
    binding: dict[str, Any], args: argparse.Namespace, thresholds: tuple[int, ...]
) -> None:
    frozen = binding["frozen_parameters"]
    extraction = frozen.get("extraction") or {}
    scheduler = frozen.get("scheduler") or {}
    expected = {
        "mapq_min": args.mapq_min,
        "baseq_min": args.baseq_min,
        "bridge_thresholds": list(thresholds),
        "workers": args.workers,
        "samtools_threads": args.samtools_threads,
    }
    if extraction != expected:
        raise RuntimeError("extraction CLI parameters drift from the frozen release contract")
    if (
        args.task_timeout_seconds != scheduler.get("task_timeout_seconds")
        or args.timeout_grace_seconds != scheduler.get("timeout_grace_seconds")
        or args.max_new_tasks not in {
            scheduler.get("initial_batch_tasks"), scheduler.get("subsequent_batch_tasks")
        }
    ):
        raise RuntimeError("extraction scheduler/batch parameters drift from the frozen release contract")
    if args.ignore_resource_gate:
        raise RuntimeError("--ignore-resource-gate is forbidden for a release-bound extraction")


def embedded_extractor_producers_match_release(
    results: Sequence[dict[str, Any]], binding: dict[str, Any]
) -> bool:
    expected = binding["snapshot_sources"]["extractor"]
    expected_manifest = binding["canonical_input_manifest"]
    for result in results:
        provenance = (result.get("receipt") or {}).get("provenance") or {}
        producer = provenance.get("extractor") or {}
        manifest = provenance.get("manifest") or {}
        try:
            path = Path(str(producer.get("path", ""))).resolve(strict=True)
            manifest_path = Path(str(manifest.get("path", ""))).resolve(strict=True)
        except OSError:
            return False
        if str(path) != expected["path"] or producer.get("sha256") != expected["sha256"]:
            return False
        if (
            str(manifest_path) != expected_manifest["path"]
            or manifest.get("sha256") != expected_manifest["sha256"]
        ):
            return False
    return True


def write_sha256_sidecar(path: Path) -> Path:
    checksum_path = path.with_name(f"{path.name}.sha256")
    checksum_path.write_text(f"{sha256_path(path)}  {path.name}\n", encoding="ascii")
    return checksum_path


def write_json_and_sha256_exclusive(path: Path, payload: dict[str, Any]) -> Path:
    """Create an immutable receipt/checkpoint and checksum without overwriting."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    checksum_path = path.with_name(f"{path.name}.sha256")
    try:
        with checksum_path.open("x", encoding="ascii") as handle:
            handle.write(f"{sha256_path(path)}  {path.name}\n")
    except BaseException:
        # The JSON remains as an auditable incomplete write; never hide it by
        # overwriting on resume.
        raise
    return checksum_path


def verify_sha256_sidecar(path: Path, sidecar_name: str | None = None) -> bool:
    checksum_path = path.parent / (sidecar_name or f"{path.name}.sha256")
    if not checksum_path.is_file():
        return False
    fields = checksum_path.read_text(encoding="ascii", errors="strict").strip().split()
    return len(fields) == 2 and fields[1] == path.name and fields[0] == sha256_path(path)


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
    path = Path("/proc/sys/kernel/random/boot_id")
    try:
        value = path.read_text(encoding="ascii").strip()
    except OSError:
        value = "unavailable"
    return value


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
    return {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": f"{path.name}.sha256",
        "covers": path.name,
    }


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
        observed = os.lstat(raw)
        sidecar = raw.with_name(f"{raw.name}.sha256")
        side_observed = os.lstat(sidecar)
    except OSError as exc:
        raise RuntimeError(f"{label} is missing: {exc}") from exc
    if (
        not stat.S_ISREG(observed.st_mode) or stat.S_ISLNK(observed.st_mode)
        or stat.S_IMODE(observed.st_mode) != 0o444 or observed.st_nlink != 1
        or not stat.S_ISREG(side_observed.st_mode) or stat.S_ISLNK(side_observed.st_mode)
        or stat.S_IMODE(side_observed.st_mode) != 0o444 or side_observed.st_nlink != 1
        or not verify_sha256_sidecar(raw)
    ):
        raise RuntimeError(f"{label} immutable file/sidecar contract failed")
    try:
        payload = json.loads(raw.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"{label} is unreadable: {exc}") from exc
    if not isinstance(payload, dict):
        raise RuntimeError(f"{label} is not a JSON object")
    if "receipt_integrity" in payload and payload["receipt_integrity"] != _receipt_integrity(raw):
        raise RuntimeError(f"{label} receipt_integrity filename binding mismatch")
    return payload


def release_producer_sources(binding: Mapping[str, Any]) -> dict[str, Any]:
    sources = binding["snapshot_sources"]
    return {
        "runner": dict(sources["full_extraction_runner"]),
        "child_producer": dict(sources["extractor"]),
        "dependencies": {"lossless_read_contract": dict(sources["lossless_read_contract"])},
    }


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
    receipt_path: Path | None = None,
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
    process_snapshot = {
        **process_core,
        "semantic_sha256": semantic_json_sha256(process_core),
    }
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
        **filesystem_core,
        "semantic_sha256": semantic_json_sha256(filesystem_core),
    }
    checks = {
        "zero_conflict": process_core == {
            "process_count": 0, "root_count": 0, "representatives": [],
        },
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
    if gate_scope == "pilot":
        if receipt_path is None:
            raise ValueError("pilot resource gate requires an explicit receipt path")
        gate_path = receipt_path.absolute()
    else:
        if receipt_path is not None:
            raise ValueError("session/batch gate path is fixed and cannot be overridden")
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


def _session_identity(output_root: Path, session: Mapping[str, Any]) -> dict[str, Any]:
    path = output_root / "_orchestration" / "session.json"
    return {"path": str(path.resolve()), "sha256": sha256_path(path), "session_id": session["session_id"]}


def ensure_release_orchestration_session(
    output_root: Path,
    binding: Mapping[str, Any],
    run_contract: Mapping[str, Any],
    gate: Mapping[str, Any] | None,
) -> dict[str, Any]:
    output_root = output_root.absolute()
    assert_no_symlink_ancestors(output_root)
    session_path = output_root / "_orchestration" / "session.json"
    if session_path.exists():
        if output_root.is_symlink() or not output_root.is_dir():
            raise RuntimeError("release output root is not a real directory")
        session = load_immutable_json(session_path, "orchestration session")
    else:
        if not output_root.exists():
            output_root.parent.mkdir(parents=True, exist_ok=True)
            output_root.mkdir()
        if output_root.is_symlink() or not output_root.is_dir():
            raise RuntimeError("release output root is not a real directory")
        if gate is None:
            raise RuntimeError("new extraction session requires an authenticated resource gate")
        root_stat = os.lstat(output_root)
        if not stat.S_ISDIR(root_stat.st_mode) or root_stat.st_nlink < 2:
            raise RuntimeError("new release output root directory identity is invalid")
        payload: dict[str, Any] = {
            "schema_name": "intersubmod.m2_orchestration_session",
            "schema_version": ORCHESTRATION_SCHEMA_VERSION,
            "stage": "extraction",
            "session_id": "",
            "session_nonce": secrets.token_hex(32),
            "created_at_utc": utc_now(),
            "created_monotonic_ns": time.monotonic_ns(),
            "host_boot_id": host_boot_id(),
            "release_manifest": dict(binding["release_manifest"]),
            "release_binding_semantic_sha256": semantic_json_sha256(binding),
            "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
            "scope": {
                "datasets": list(DATASETS), "chromosomes": list(AUTOSOMES),
                "expected_tasks": EXPECTED_TASKS,
            },
            "output_root": {
                "path": str(output_root.resolve()), "st_dev": int(root_stat.st_dev),
                "st_ino": int(root_stat.st_ino),
            },
            "producer_sources": release_producer_sources(binding),
            "scheduler_policy": dict(ORCHESTRATION_POLICY),
            "parent_extraction": None,
            "resource_gate": dict(gate),
            "receipt_integrity": _receipt_integrity(session_path),
        }
        payload["session_id"] = semantic_json_sha256({
            key: value for key, value in payload.items()
            if key not in {"session_id", "receipt_integrity"}
        })
        write_immutable_json_exclusive(session_path, payload)
        session = payload
    if set(session) != SESSION_KEYS:
        raise RuntimeError("orchestration session exact-key schema mismatch")
    root_stat = os.lstat(output_root)
    if not stat.S_ISDIR(root_stat.st_mode) or root_stat.st_nlink < 2:
        raise RuntimeError("release output root directory identity is invalid")
    expected_static = {
        "schema_name": "intersubmod.m2_orchestration_session",
        "schema_version": ORCHESTRATION_SCHEMA_VERSION, "stage": "extraction",
        "release_manifest": dict(binding["release_manifest"]),
        "release_binding_semantic_sha256": semantic_json_sha256(binding),
        "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
        "scope": {"datasets": list(DATASETS), "chromosomes": list(AUTOSOMES), "expected_tasks": EXPECTED_TASKS},
        "output_root": {"path": str(output_root.resolve()), "st_dev": int(root_stat.st_dev), "st_ino": int(root_stat.st_ino)},
        "producer_sources": release_producer_sources(binding),
        "scheduler_policy": dict(ORCHESTRATION_POLICY),
        "parent_extraction": None,
    }
    if any(session.get(key) != value for key, value in expected_static.items()):
        raise RuntimeError("orchestration session binding/root inode mismatch")
    session_target = {
        "output_root": expected_static["output_root"],
        "release_manifest": dict(binding["release_manifest"]),
        "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
    }
    gate_payload = load_resource_gate_receipt(
        session.get("resource_gate") or {},
        expected_path=_resource_gate_path(output_root, "session"),
        expected_stage="extraction",
        expected_scope="session",
        expected_target=session_target,
        expected_producer_source=session["producer_sources"]["runner"],
    )
    if gate is not None and dict(gate) != session.get("resource_gate"):
        raise RuntimeError("new session resource gate identity differs from persisted session")
    if (
        gate_payload.get("host_boot_id") != session.get("host_boot_id")
        or gate_payload.get("observed_monotonic_ns", -1) > session.get("created_monotonic_ns", -1)
    ):
        raise RuntimeError("orchestration session is not temporally bound to its resource gate")
    expected_session_id = semantic_json_sha256({
        key: value for key, value in session.items()
        if key not in {"session_id", "receipt_integrity"}
    })
    if session.get("session_id") != expected_session_id:
        raise RuntimeError("orchestration session semantic identity mismatch")
    return session


def _task_id(dataset: str, chrom: str) -> str:
    return f"{dataset}:{chrom}"


def _canonical_task_ordinal(dataset: str, chrom: str) -> int:
    return DATASETS.index(dataset) * len(AUTOSOMES) + AUTOSOMES.index(chrom) + 1


def canonical_sort_results(results: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        results,
        key=lambda item: (DATASETS.index(item["dataset"]), AUTOSOMES.index(item["chrom"])),
    )


def _compact_extraction_results(
    results: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Return the deterministic scientific child rows persisted in checkpoints.

    Process diagnostics such as stdout, PID, and elapsed time are invocation-local
    and cannot be reproduced on resume.  They therefore must not be part of the
    persisted cumulative scientific result rows.  Child receipt bytes and, for a
    release run, the completion identity are reproducible and are retained.
    """
    compact: list[dict[str, Any]] = []
    for result in canonical_sort_results([dict(item) for item in results]):
        row = {
            "dataset": result["dataset"],
            "chrom": result["chrom"],
            "status": result["status"],
            "receipt": result.get("receipt"),
        }
        if result.get("orchestration_completion") is not None:
            row["orchestration_completion"] = result["orchestration_completion"]
        compact.append(row)
    return compact


def create_batch_start_and_grants(
    output_root: Path,
    session: Mapping[str, Any],
    run_contract: Mapping[str, Any],
    selected_specs: Sequence[tuple[Any, ...]],
    *,
    before_count: int,
    previous_chain_head: Mapping[str, str] | None,
    batch_index: int,
    max_new_tasks: int,
    gate: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    effective = len(selected_specs)
    expected_max = 8 if before_count == 0 else 16
    expected_effective = min(expected_max, EXPECTED_TASKS - before_count)
    if (
        max_new_tasks != expected_max or effective != expected_effective
        or before_count not in (0, *LEGAL_CUMULATIVE_COUNTS[:-1])
    ):
        raise RuntimeError("illegal release batch size/count chain")
    selected_ids = [_task_id(spec[0], spec[1]) for spec in selected_specs]
    session_path = output_root / "_orchestration" / "session.json"
    batch_path = output_root / "_orchestration" / "batches" / f"batch_{batch_index:03d}_start.json"
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
        expected_stage="extraction",
        expected_scope="batch",
        expected_target=gate_target,
        expected_producer_source=session["producer_sources"]["runner"],
    )
    payload: dict[str, Any] = {
        "schema_name": "intersubmod.m2_orchestration_batch_start",
        "schema_version": ORCHESTRATION_SCHEMA_VERSION, "stage": "extraction",
        "session_id": session["session_id"], "session_sha256": sha256_path(session_path),
        "batch_index": batch_index, "batch_id": "", "batch_nonce": secrets.token_hex(32),
        "previous_chain_head": None if previous_chain_head is None else dict(previous_chain_head),
        "before_count": before_count, "max_new_tasks": max_new_tasks,
        "effective_count": effective, "selected_task_ids": selected_ids,
        "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
        "runner_source": dict(session["producer_sources"]["runner"]),
        "resource_gate": dict(gate), "created_at_utc": utc_now(),
        "created_monotonic_ns": time.monotonic_ns(), "host_boot_id": host_boot_id(),
        "receipt_integrity": _receipt_integrity(batch_path),
    }
    payload["batch_id"] = semantic_json_sha256({
        key: value for key, value in payload.items()
        if key not in {"batch_id", "receipt_integrity"}
    })
    if (
        gate_payload.get("host_boot_id") != payload["host_boot_id"]
        or gate_payload.get("observed_monotonic_ns", -1) > payload["created_monotonic_ns"]
    ):
        raise RuntimeError("batch start is not temporally bound to its resource gate")
    write_immutable_json_exclusive(batch_path, payload)
    batch_sha = sha256_path(batch_path)
    grants: dict[str, dict[str, Any]] = {}
    for spec in selected_specs:
        dataset, chrom, task_dir, command = spec[:4]
        task_id = _task_id(dataset, chrom)
        grant_path = output_root / "_orchestration" / "tasks" / dataset / chrom / "grant.json"
        grant = {
            "schema_name": "intersubmod.m2_orchestration_child_grant",
            "schema_version": ORCHESTRATION_SCHEMA_VERSION, "stage": "extraction",
            "session_id": session["session_id"], "session_sha256": sha256_path(session_path),
            "batch_id": payload["batch_id"], "batch_start_sha256": batch_sha,
            "task_id": task_id, "dataset": dataset, "chrom": chrom,
            "task_ordinal": _canonical_task_ordinal(dataset, chrom),
            "child_nonce": secrets.token_hex(32), "command": list(command),
            "command_semantic_sha256": semantic_json_sha256(list(command)),
            "producer_sources": dict(session["producer_sources"]),
            "input_identity": {
                "manifest_path": str(Path(run_contract["release_binding"]["canonical_input_manifest"]["path"]).resolve()),
                "manifest_sha256": run_contract["manifest_sha256"],
            },
            "parameters_semantic_sha256": semantic_json_sha256(spec[4]),
            "expected_output_dir": str(task_dir.resolve()), "issued_at_utc": utc_now(),
            "issued_monotonic_ns": time.monotonic_ns(), "host_boot_id": host_boot_id(),
            "receipt_integrity": _receipt_integrity(grant_path),
        }
        write_immutable_json_exclusive(grant_path, grant)
        grants[task_id] = {"path": str(grant_path.resolve()), "sha256": sha256_path(grant_path), "document": grant}
    return payload | {"path": str(batch_path.resolve()), "sha256": batch_sha}, grants


def write_child_completion(
    output_root: Path,
    session: Mapping[str, Any],
    batch: Mapping[str, Any],
    grant_record: Mapping[str, Any],
    result: Mapping[str, Any],
) -> dict[str, str]:
    dataset, chrom = str(result["dataset"]), str(result["chrom"])
    task_id = _task_id(dataset, chrom)
    if result.get("status") != "PASS" or int(result.get("returncode", -1)) != 0:
        raise RuntimeError(f"cannot attest nonpassing child: {task_id}")
    receipt_path = output_root / "samples" / dataset / chrom / "receipt.json"
    if not verify_sha256_sidecar(receipt_path):
        raise RuntimeError(f"child receipt sidecar failed before attestation: {task_id}")
    receipt_path.chmod(0o444)
    receipt_path.with_name(f"{receipt_path.name}.sha256").chmod(0o444)
    receipt = load_immutable_json(receipt_path, f"child receipt {task_id}")
    grant_path = Path(str(grant_record["path"]))
    grant = grant_record["document"]
    completion_path = output_root / "_orchestration" / "tasks" / dataset / chrom / "completion.json"
    payload = {
        "schema_name": "intersubmod.m2_orchestration_child_completion",
        "schema_version": ORCHESTRATION_SCHEMA_VERSION, "stage": "extraction",
        "session_id": session["session_id"],
        "session_sha256": sha256_path(output_root / "_orchestration" / "session.json"),
        "batch_id": batch["batch_id"],
        "grant_identity": {
            "path": str(grant_path.resolve()), "sha256": sha256_path(grant_path),
            "semantic_sha256": semantic_json_sha256(grant),
        },
        "task_id": task_id, "dataset": dataset, "chrom": chrom,
        "task_ordinal": _canonical_task_ordinal(dataset, chrom),
        "child_receipt_identity": {
            "path": str(receipt_path.resolve()), "size_bytes": receipt_path.stat().st_size,
            "sha256": sha256_path(receipt_path), "semantic_sha256": semantic_json_sha256(receipt),
        },
        "child_outputs_semantic_sha256": semantic_json_sha256(receipt.get("outputs") or {}),
        "command_semantic_sha256": grant["command_semantic_sha256"],
        "status": "PASS", "returncode": 0, "timed_out": False,
        "process_group_id": result.get("process_group_id"),
        "started_monotonic_ns": result.get("started_monotonic_ns"),
        "completed_at_utc": utc_now(), "completed_monotonic_ns": time.monotonic_ns(),
        "host_boot_id": host_boot_id(),
        "runner_source": dict(session["producer_sources"]["runner"]),
        "receipt_integrity": _receipt_integrity(completion_path),
    }
    identity = write_immutable_json_exclusive(completion_path, payload)
    return {"path": identity["path"], "sha256": identity["sha256"]}


def validate_attested_completion(
    completion_path: Path,
    session: Mapping[str, Any],
    *,
    expected_task_id: str | None = None,
    expected_batch: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, str]]:
    completion = load_immutable_json(completion_path, "child completion")
    if set(completion) != COMPLETION_KEYS:
        raise RuntimeError("child completion exact-key schema mismatch")
    task_id = completion.get("task_id")
    dataset, chrom = str(completion.get("dataset")), str(completion.get("chrom"))
    session_root = Path(str((session.get("output_root") or {}).get("path", "")))
    session_path = session_root / "_orchestration" / "session.json"
    expected_completion_path = session_root / "_orchestration" / "tasks" / dataset / chrom / "completion.json"
    if (
        completion.get("schema_name") != "intersubmod.m2_orchestration_child_completion"
        or completion.get("schema_version") != ORCHESTRATION_SCHEMA_VERSION
        or completion.get("stage") != "extraction"
        or completion.get("session_id") != session.get("session_id")
        or (expected_task_id is not None and task_id != expected_task_id)
        or completion.get("status") != "PASS" or completion.get("returncode") != 0
        or completion.get("timed_out") is not False
        or completion.get("session_sha256") != sha256_path(session_path)
        or task_id != _task_id(dataset, chrom)
        or completion.get("task_ordinal") != _canonical_task_ordinal(dataset, chrom)
        or completion.get("runner_source") != session["producer_sources"]["runner"]
        or (expected_batch is not None and completion.get("host_boot_id") != expected_batch.get("host_boot_id"))
        or (expected_batch is not None and completion.get("batch_id") != expected_batch.get("batch_id"))
        or completion_path.resolve() != expected_completion_path.resolve()
    ):
        raise RuntimeError("child completion session/task/status mismatch")
    grant_identity = completion.get("grant_identity") or {}
    grant_path = Path(str(grant_identity.get("path", "")))
    grant = load_immutable_json(grant_path, "child grant")
    expected_output_dir = Path(str((session.get("output_root") or {}).get("path", ""))) / "samples" / dataset / chrom
    if (
        set(grant) != GRANT_KEYS
        or grant.get("schema_name") != "intersubmod.m2_orchestration_child_grant"
        or grant.get("schema_version") != ORCHESTRATION_SCHEMA_VERSION
        or grant.get("stage") != "extraction"
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
        or Path(str(grant.get("expected_output_dir", ""))).resolve() != expected_output_dir.resolve()
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
        raise RuntimeError("child grant/completion binding mismatch")
    child_identity = completion.get("child_receipt_identity") or {}
    receipt_path = Path(str(child_identity.get("path", "")))
    receipt = load_immutable_json(receipt_path, "attested child receipt")
    if (
        (receipt_path.parent / TASK_STATUS_NAME).exists()
        or sha256_path(receipt_path) != child_identity.get("sha256")
        or receipt_path.stat().st_size != child_identity.get("size_bytes")
        or semantic_json_sha256(receipt) != child_identity.get("semantic_sha256")
        or semantic_json_sha256(receipt.get("outputs") or {})
        != completion.get("child_outputs_semantic_sha256")
    ):
        raise RuntimeError("attested child receipt/output identity mismatch")
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
            raise RuntimeError("attested child declared output identity mismatch")
    provenance = receipt.get("provenance") or {}
    manifest_identity = provenance.get("manifest") or {}
    if grant.get("input_identity") != {
        "manifest_path": str(Path(str(manifest_identity.get("path", ""))).resolve()),
        "manifest_sha256": manifest_identity.get("sha256"),
    }:
        raise RuntimeError("child grant input manifest differs from child receipt")
    manifest_path = Path(str(grant["input_identity"]["manifest_path"]))
    if sha256_path(manifest_path) != grant["input_identity"]["manifest_sha256"]:
        raise RuntimeError("child grant input manifest physical SHA mismatch")
    child_parameters = receipt.get("parameters") or {}
    expected_grant_parameters = {
        key: child_parameters.get(key)
        for key in (
            "mapq_min", "baseq_min", "samtools_threads", "bridge_thresholds",
            "component_linkage_bases",
        )
    }
    if grant.get("parameters_semantic_sha256") != semantic_json_sha256(expected_grant_parameters):
        raise RuntimeError("child grant parameter identity differs from child receipt")
    identity = {"path": str(completion_path.resolve()), "sha256": sha256_path(completion_path)}
    return completion, identity


def load_release_orchestration_state(
    output_root: Path,
    session: Mapping[str, Any],
) -> dict[str, Any]:
    checkpoint_dir = output_root / "checkpoints"
    receipt_paths = list(checkpoint_dir.glob("checkpoint_*_of_154.json")) if checkpoint_dir.exists() else []
    terminal = output_root / "full_extraction_receipt.json"
    if terminal.exists():
        receipt_paths.append(terminal)
    chain_rows: list[tuple[int, Path, dict[str, Any]]] = []
    for path in receipt_paths:
        receipt = load_immutable_json(path, "orchestration chain receipt")
        orchestration = receipt.get("orchestration") or {}
        count = orchestration.get("cumulative_attested_tasks")
        if not isinstance(count, int):
            raise RuntimeError("chain receipt lacks cumulative orchestration count")
        expected_path = (
            terminal
            if count == EXPECTED_TASKS
            else checkpoint_dir / f"checkpoint_{count:03d}_of_154.json"
        )
        if (
            path.resolve() != expected_path.resolve()
            or receipt.get("receipt_integrity") != _receipt_integrity(path)
        ):
            raise RuntimeError("orchestration chain receipt exact path/integrity mismatch")
        if (count == EXPECTED_TASKS) != (path == terminal):
            raise RuntimeError("terminal orchestration count/path mismatch")
        chain_rows.append((count, path, receipt))
    chain_rows.sort(key=lambda row: row[0])
    if [row[0] for row in chain_rows] != list(LEGAL_CUMULATIVE_COUNTS[:len(chain_rows)]):
        raise RuntimeError("orchestration checkpoint count chain has a gap/reorder")
    previous: dict[str, str] | None = None
    completions: dict[str, dict[str, str]] = {}
    child_receipts: dict[str, dict[str, Any]] = {}
    for index, (count, path, receipt) in enumerate(chain_rows, start=1):
        orchestration = receipt.get("orchestration") or {}
        if set(orchestration) != {
            "session_identity", "batch_start_identity", "previous_chain_head",
            "batch_completion_attestations", "cumulative_attested_tasks",
        }:
            raise RuntimeError("receipt orchestration exact-key schema mismatch")
        if orchestration["session_identity"] != _session_identity(output_root, session):
            raise RuntimeError("receipt is bound to a different orchestration session")
        if orchestration["previous_chain_head"] != previous:
            raise RuntimeError("orchestration previous-chain-head mismatch")
        start_identity = orchestration["batch_start_identity"]
        start_path = Path(str(start_identity.get("path", "")))
        start = load_immutable_json(start_path, "batch start")
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
            or start.get("session_sha256") != sha256_path(output_root / "_orchestration" / "session.json")
            or start.get("run_contract_semantic_sha256") != session.get("run_contract_semantic_sha256")
            or start.get("runner_source") != session["producer_sources"]["runner"]
            or start.get("batch_id") != semantic_json_sha256({
                key: value for key, value in start.items()
                if key not in {"batch_id", "receipt_integrity"}
            })
        ):
            raise RuntimeError("batch start/receipt chain binding mismatch")
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
            expected_stage="extraction",
            expected_scope="batch",
            expected_target=gate_target,
            expected_producer_source=session["producer_sources"]["runner"],
        )
        if (
            gate_payload.get("host_boot_id") != start.get("host_boot_id")
            or gate_payload.get("observed_monotonic_ns", -1) > start.get("created_monotonic_ns", -1)
        ):
            raise RuntimeError("batch start is not temporally bound to resource gate")
        if (
            start.get("before_count") != expected_before
            or start.get("max_new_tasks") != (8 if index == 1 else 16)
            or start.get("effective_count") != expected_effective
            or start.get("selected_task_ids") != expected_selected_ids
        ):
            raise RuntimeError("batch scheduler count chain mismatch")
        attested = orchestration["batch_completion_attestations"]
        if not isinstance(attested, list) or len(attested) != expected_effective:
            raise RuntimeError("batch completion cardinality mismatch")
        if [row.get("task_id") for row in attested] != start.get("selected_task_ids"):
            raise RuntimeError("batch completion order differs from batch grant order")
        for row in attested:
            task_id = row["task_id"]
            if task_id in completions:
                raise RuntimeError("duplicate child completion across batch chain")
            completion, identity = validate_attested_completion(
                Path(row["path"]), session, expected_task_id=task_id,
                expected_batch=start | {"sha256": sha256_path(start_path)},
            )
            if identity != {"path": row["path"], "sha256": row["sha256"]}:
                raise RuntimeError("child completion identity differs from receipt")
            completions[task_id] = identity
            child_receipts[task_id] = load_immutable_json(
                Path(str(completion["child_receipt_identity"]["path"])),
                f"attested extraction child receipt {task_id}",
            )
        result_rows = receipt.get("results")
        expected_cumulative_task_ids = [
            _task_id(dataset, chrom)
            for dataset in DATASETS for chrom in AUTOSOMES
        ][:count]
        if not isinstance(result_rows, list) or len(result_rows) != count:
            raise RuntimeError("receipt results do not cover the cumulative attested task count")
        observed_result_task_ids = [
            _task_id(result.get("dataset"), result.get("chrom"))
            for result in result_rows
            if isinstance(result, dict)
        ]
        if observed_result_task_ids != expected_cumulative_task_ids:
            raise RuntimeError(
                "receipt results are not an exact ordered bijection of cumulative child completions"
            )
        for result, task_id in zip(result_rows, observed_result_task_ids):
            if result.get("orchestration_completion") != completions.get(task_id):
                raise RuntimeError("receipt result lacks its exact child completion identity")
        run_contract = receipt.get("run_contract")
        invocation = receipt.get("invocation")
        elapsed_seconds = receipt.get("elapsed_seconds")
        thresholds_raw = (
            run_contract.get("bridge_thresholds")
            if isinstance(run_contract, Mapping) else None
        )
        if (
            not isinstance(run_contract, dict)
            or semantic_json_sha256(run_contract)
            != session.get("run_contract_semantic_sha256")
            or not isinstance(invocation, dict)
            or not isinstance(elapsed_seconds, (int, float))
            or isinstance(elapsed_seconds, bool)
            or not math.isfinite(float(elapsed_seconds))
            or float(elapsed_seconds) < 0.0
            or not isinstance(thresholds_raw, list)
            or not thresholds_raw
            or any(not isinstance(value, int) or isinstance(value, bool) for value in thresholds_raw)
        ):
            raise RuntimeError("extraction checkpoint rebuild inputs are invalid")
        current_batch_ids = set(start["selected_task_ids"])
        rebuilt_results = []
        for task_id in expected_cumulative_task_ids:
            dataset, chrom = task_id.split(":", 1)
            rebuilt_results.append({
                "dataset": dataset,
                "chrom": chrom,
                "status": "PASS" if task_id in current_batch_ids else "REUSED_PASS",
                "receipt": child_receipts[task_id],
                "orchestration_completion": completions[task_id],
            })
        rebuilt = build_extraction_receipt(
            rebuilt_results,
            tuple(thresholds_raw),
            run_contract,
            invocation,
            float(elapsed_seconds),
        )
        rebuilt["orchestration"] = orchestration
        rebuilt["receipt_integrity"] = _receipt_integrity(path)
        if receipt != rebuilt:
            raise RuntimeError(
                "extraction checkpoint/terminal aggregate, checks, or canonical rows drifted"
            )
        previous = {"path": str(path.resolve()), "sha256": sha256_path(path)}
    starts = list((output_root / "_orchestration" / "batches").glob("batch_*_start.json")) if (output_root / "_orchestration" / "batches").exists() else []
    if len(starts) != len(chain_rows):
        raise RuntimeError("open/orphan orchestration batch detected")
    gate_dir = output_root / "_orchestration" / "resource_gates"
    gate_paths = set(gate_dir.glob("*.json")) if gate_dir.exists() else set()
    expected_gate_paths = {_resource_gate_path(output_root, "session")}
    expected_gate_paths.update(
        _resource_gate_path(output_root, "batch", index)
        for index in range(1, len(chain_rows) + 1)
    )
    if {path.resolve() for path in gate_paths} != {path.resolve() for path in expected_gate_paths}:
        raise RuntimeError("missing/orphan/reused orchestration resource gate detected")
    grant_paths = list((output_root / "_orchestration" / "tasks").glob("*/*/grant.json")) if (output_root / "_orchestration" / "tasks").exists() else []
    completion_paths = list((output_root / "_orchestration" / "tasks").glob("*/*/completion.json")) if (output_root / "_orchestration" / "tasks").exists() else []
    if len(grant_paths) != len(completions) or len(completion_paths) != len(completions):
        raise RuntimeError("unreferenced/orphan child grant or completion detected")
    return {
        "count": len(completions), "completions": completions,
        "child_receipts": child_receipts,
        "previous_chain_head": previous, "next_batch_index": len(chain_rows) + 1,
    }


def conflict_kind(command: str) -> str | None:
    """Classify actual script argv tokens, not incidental grep/sed text."""
    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()
    basenames = [Path(token).name for token in tokens]
    if any(name.startswith("analyze_all_ssnv_focal_alt_multigroup") and name.endswith(".py") for name in basenames):
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
        return "m2_full_runner"
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
            row = {
                "pid": int(pid),
                "ppid": int(ppid),
                "elapsed_seconds": int(elapsed),
                "cpu_percent": float(cpu),
                "command": command,
            }
        except ValueError:
            continue
        rows[row["pid"]] = row
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
        representatives.append(
            {
                **representative,
                "group_process_count": len(members),
                "group_cpu_percent_sum": sum(member["cpu_percent"] for member in members),
                "group_kinds": sorted({member["kind"] for member in members}),
                "member_pids": sorted(member["pid"] for member in members),
            }
        )
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
            ["ps", "-eo", "pid=,ppid=,etimes=,%cpu=,cmd="], text=True, errors="replace"
        )
    pids = parse_process_table(output)
    return summarize_conflicts(pids, process_family(pids, current))


def load_passing_receipt(
    path: Path,
    *,
    dataset: str,
    chrom: str,
    expected_parameters: dict[str, Any],
    manifest_sha256: str,
    extractor_sha256: str,
) -> dict[str, Any] | None:
    if (path.parent / TASK_STATUS_NAME).exists():
        return None
    if not path.is_file():
        return None
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    if receipt.get("schema_name") != "intersubmod.lossless_read_linkage_chromosome_receipt":
        return None
    if receipt.get("schema_version") != "1.2.0" or receipt.get("all_pass") is not True:
        return None
    if receipt.get("scope", {}).get("dataset") != dataset or receipt.get("scope", {}).get("chrom") != chrom:
        return None
    parameters = receipt.get("parameters", {})
    if any(parameters.get(key) != value for key, value in expected_parameters.items()):
        return None
    provenance = receipt.get("provenance", {})
    if (provenance.get("manifest") or {}).get("sha256") != manifest_sha256:
        return None
    if (provenance.get("extractor") or {}).get("sha256") != extractor_sha256:
        return None
    integrity = receipt.get("receipt_integrity", {})
    if integrity.get("scheme") != "external_sha256_sidecar_v1":
        return None
    try:
        if not verify_sha256_sidecar(path, integrity.get("sidecar_name")):
            return None
        receipt_dir = path.parent.resolve()
        for identity in receipt.get("outputs", {}).values():
            output = Path(identity["path"])
            if output.resolve().parent != receipt_dir or not output.is_file():
                return None
            if output.stat().st_size != identity.get("size_bytes") or sha256_path(output) != identity.get("sha256"):
                return None
    except (KeyError, OSError, UnicodeError, TypeError):
        return None
    return receipt


def task_command(
    manifest: Path,
    dataset: str,
    chrom: str,
    task_dir: Path,
    args: argparse.Namespace,
) -> list[str]:
    return [
        sys.executable,
        str(EXTRACTOR),
        "--manifest",
        str(manifest),
        "--sample",
        dataset,
        "--chrom",
        chrom,
        "--output-dir",
        str(task_dir),
        "--mapq-min",
        str(args.mapq_min),
        "--baseq-min",
        str(args.baseq_min),
        "--bridge-thresholds",
        args.bridge_thresholds,
        "--samtools-threads",
        str(args.samtools_threads),
    ]


def run_task(spec: tuple[Any, ...]) -> dict[str, Any]:
    (
        dataset,
        chrom,
        task_dir,
        command,
        expected_parameters,
        manifest_sha256,
        extractor_sha256,
        task_timeout_seconds,
        timeout_grace_seconds,
    ) = spec
    receipt_path = task_dir / "receipt.json"
    existing = load_passing_receipt(
        receipt_path,
        dataset=dataset,
        chrom=chrom,
        expected_parameters=expected_parameters,
        manifest_sha256=manifest_sha256,
        extractor_sha256=extractor_sha256,
    )
    if existing is not None:
        return {"dataset": dataset, "chrom": chrom, "status": "REUSED_PASS", "receipt": existing}
    if task_dir.exists():
        return {
            "dataset": dataset,
            "chrom": chrom,
            "status": "FAIL_EXISTING_NONPASS_DIRECTORY",
            "task_dir": str(task_dir),
        }
    started = time.monotonic()
    started_monotonic_ns = time.monotonic_ns()
    process = run_command_with_process_group_timeout(
        command,
        task_timeout_seconds=task_timeout_seconds,
        timeout_grace_seconds=timeout_grace_seconds,
    )
    if process["timed_out"]:
        task_dir.mkdir(parents=True, exist_ok=True)
        marker_path = task_dir / TASK_STATUS_NAME
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
            "task_dir": str(task_dir),
            "task_status_path": str(marker_path),
            "command": command,
            "process_group_id": process["process_group_id"],
            "termination_stage": process["termination_stage"],
        }
    receipt = load_passing_receipt(
        receipt_path,
        dataset=dataset,
        chrom=chrom,
        expected_parameters=expected_parameters,
        manifest_sha256=manifest_sha256,
        extractor_sha256=extractor_sha256,
    )
    status = "PASS" if process["returncode"] == 0 and receipt is not None else "FAIL"
    marker_path: Path | None = None
    if status == "FAIL":
        task_dir.mkdir(parents=True, exist_ok=True)
        marker_path = task_dir / TASK_STATUS_NAME
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
        "receipt": receipt,
        "command": command,
        "task_status_path": None if marker_path is None else str(marker_path),
        "process_group_id": process["process_group_id"],
        "started_monotonic_ns": started_monotonic_ns,
    }


def run_command_with_process_group_timeout(
    command: Sequence[str], *, task_timeout_seconds: float, timeout_grace_seconds: float
) -> dict[str, Any]:
    """Run one task in a new process group and kill the whole group on timeout."""
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
        # ``communicate`` can return when the direct child exits even though a
        # grandchild with redirected stdio still survives in the process group.
        # Honor the same grace deadline for that case, then kill the group.
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


def scan_existing_specs(
    specs: Sequence[tuple[Any, ...]],
    *,
    attested_completions: Mapping[str, Mapping[str, str]] | None = None,
) -> tuple[list[dict[str, Any]], list[tuple[Any, ...]], list[dict[str, str]]]:
    """Validate every existing child before any pending task is launched."""
    reused: list[dict[str, Any]] = []
    pending: list[tuple[Any, ...]] = []
    invalid: list[dict[str, str]] = []
    for spec in specs:
        dataset, chrom, task_dir = spec[:3]
        if not task_dir.exists():
            pending.append(spec)
            continue
        receipt = load_passing_receipt(
            task_dir / "receipt.json",
            dataset=dataset,
            chrom=chrom,
            expected_parameters=spec[4],
            manifest_sha256=spec[5],
            extractor_sha256=spec[6],
        )
        if receipt is None:
            invalid.append({
                "dataset": dataset,
                "chrom": chrom,
                "task_dir": str(task_dir),
                "status": "FAIL_EXISTING_NONPASS_DIRECTORY",
            })
        else:
            task_id = _task_id(dataset, chrom)
            completion = None if attested_completions is None else attested_completions.get(task_id)
            if attested_completions is not None and completion is None:
                invalid.append({
                    "dataset": dataset, "chrom": chrom, "task_dir": str(task_dir),
                    "status": "FAIL_UNATTESTED_OR_ORPHAN_CHILD",
                })
                continue
            reused.append({
                "dataset": dataset,
                "chrom": chrom,
                "status": "REUSED_PASS",
                "receipt": receipt,
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
    """Execute with at most ``workers`` submitted futures, fail closed."""
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
            for future in sorted(
                done,
                key=lambda item: specs.index(pending[item]),
            ):
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


def aggregate_results(results: list[dict[str, Any]], thresholds: tuple[int, ...]) -> dict[str, Any]:
    counts = Counter()
    component_by_basis = {
        basis: {str(threshold): Counter() for threshold in thresholds} for basis in LINKAGE_BASES
    }
    component_max_by_basis = {
        basis: {str(threshold): {"max_k_component_sites": 0, "max_k": 0} for threshold in thresholds}
        for basis in LINKAGE_BASES
    }
    k_by_basis = {
        basis: {str(threshold): Counter() for threshold in thresholds} for basis in LINKAGE_BASES
    }
    by_dataset = defaultdict(lambda: {"tasks": Counter(), "counts": Counter()})
    phase_contract_totals = {
        "known_phase_set_chromosome_units_by_hp_family": Counter(),
        "known_ps_active_site_memberships_by_hp_family": Counter(),
        "missing_ps_active_sites_by_hp_family": Counter(),
    }
    legacy_cross_ps: dict[str, dict[str, Counter[str]]] = defaultdict(
        lambda: defaultdict(Counter)
    )
    for result in results:
        dataset = result["dataset"]
        by_dataset[dataset]["tasks"][result["status"]] += 1
        receipt = result.get("receipt")
        if not receipt:
            continue
        counts.update({key: value for key, value in receipt.get("counts", {}).items() if isinstance(value, int)})
        by_dataset[dataset]["counts"].update(
            {key: value for key, value in receipt.get("counts", {}).items() if isinstance(value, int)}
        )
        phase_counts = receipt.get("phase_set_contract_counts") or {}
        for target_key, source_key in (
            ("known_phase_set_chromosome_units_by_hp_family", "known_phase_sets_by_hp_family"),
            ("known_ps_active_site_memberships_by_hp_family", "known_ps_active_site_memberships_by_hp_family"),
            ("missing_ps_active_sites_by_hp_family", "missing_ps_active_sites_by_hp_family"),
        ):
            phase_contract_totals[target_key].update({
                key: int(value) for key, value in (phase_counts.get(source_key) or {}).items()
            })
        for family, threshold_values in (receipt.get("legacy_cross_phase_set_aggregation_audit") or {}).items():
            for threshold, metrics in threshold_values.items():
                legacy_cross_ps[str(family)][str(threshold)].update({
                    key: int(value) for key, value in metrics.items()
                })
        summaries = receipt["component_summary_by_linkage_basis"]
        for basis in LINKAGE_BASES:
            for threshold in thresholds:
                key = str(threshold)
                summary = summaries[basis][key]
                component_by_basis[basis][key].update(
                    {
                        name: value
                        for name, value in summary.items()
                        if name not in {
                            "k_distribution",
                            "k_component_sites_distribution",
                            "max_k_component_sites",
                            "max_k",
                        }
                        and isinstance(value, int)
                    }
                )
                for max_field in ("max_k_component_sites", "max_k"):
                    component_max_by_basis[basis][key][max_field] = max(
                        component_max_by_basis[basis][key][max_field],
                        int(summary.get(max_field, 0)),
                    )
                k_by_basis[basis][key].update(
                    {int(k): value for k, value in summary["k_component_sites_distribution"].items()}
                )
    aggregate_components = {
        basis: {
            key: {
                **dict(component_by_basis[basis][key]),
                **component_max_by_basis[basis][key],
                "k_component_sites_distribution": {
                    str(k): value for k, value in sorted(k_by_basis[basis][key].items())
                },
                "k_distribution": {
                    str(k): value for k, value in sorted(k_by_basis[basis][key].items())
                },
            }
            for key in component_by_basis[basis]
        }
        for basis in LINKAGE_BASES
    }
    return {
        "counts": dict(counts),
        "component_summary_by_linkage_basis": aggregate_components,
        "component_summary_by_threshold": aggregate_components["pooled"],
        "phase_set_contract_totals": {
            key: dict(sorted(value.items())) for key, value in phase_contract_totals.items()
        },
        "legacy_cross_phase_set_aggregation_audit": {
            family: {threshold: dict(metrics) for threshold, metrics in sorted(values.items())}
            for family, values in sorted(legacy_cross_ps.items())
        },
        "by_dataset": {
            dataset: {"task_status_counts": dict(value["tasks"]), "counts": dict(value["counts"])}
            for dataset, value in sorted(by_dataset.items())
        },
    }


def build_extraction_receipt(
    results: list[dict[str, Any]],
    thresholds: tuple[int, ...],
    run_contract: dict[str, Any],
    invocation: dict[str, Any],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """Build either a partial checkpoint payload or the terminal full receipt."""
    ordered = _compact_extraction_results(results)
    aggregate = aggregate_results(ordered, thresholds)
    status_counts = Counter(result["status"] for result in ordered)
    all_recorded_pass = all(
        result["status"] in {"PASS", "REUSED_PASS"} for result in ordered
    )
    passing = len(ordered) if all_recorded_pass else sum(
        result["status"] in {"PASS", "REUSED_PASS"} for result in ordered
    )
    task_keys = [(result["dataset"], result["chrom"]) for result in ordered]
    expected_task_keys = {(dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES}
    exact_unique_task_scope = len(task_keys) == EXPECTED_TASKS and set(task_keys) == expected_task_keys
    complete = (
        passing == EXPECTED_TASKS
        and len(ordered) == EXPECTED_TASKS
        and exact_unique_task_scope
    )
    common = {
        "schema_version": "1.2.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "scope": {
            "datasets": list(DATASETS),
            "chromosomes": list(AUTOSOMES),
            "expected_tasks": EXPECTED_TASKS,
            "n_technical_datasets": 7,
            "n_biological_samples": 6,
            "technical_replicate_pair": ["HCC1395", "HCC1395_DORADO"],
            "child_schema_version": "1.2.0",
        },
        "elapsed_seconds": elapsed_seconds,
        "run_contract": run_contract,
        "invocation": invocation,
        "operational_parameters_excluded_from_run_contract": ["max_new_tasks"],
        "task_status_counts": dict(status_counts),
        "n_results": len(ordered),
        "passing_tasks": passing,
        "remaining_tasks": EXPECTED_TASKS - passing,
        "aggregate": aggregate,
        "results": ordered,
    }
    if not complete:
        release_binding = run_contract.get("release_binding")
        release_child_check = (
            release_binding is None
            or embedded_extractor_producers_match_release(ordered, release_binding)
        )
        return common | {
            "schema_name": "intersubmod.m2_full_extraction_checkpoint",
            "checkpoint_complete": all_recorded_pass,
            "all_pass": False,
            "checks": {
                "all_recorded_tasks_pass_or_reused": all_recorded_pass,
                "passing_count_matches_results": passing == len(ordered),
                "remaining_count_conserved": passing + (EXPECTED_TASKS - passing) == EXPECTED_TASKS,
                "recorded_task_pairs_are_unique_and_canonical": (
                    len(task_keys) == len(set(task_keys))
                    and set(task_keys).issubset(expected_task_keys)
                ),
                "recorded_child_producers_match_release_contract": release_child_check,
                "terminal_full_scope_complete": False,
            },
        }
    checks = {
        "all_154_tasks_returned": len(ordered) == EXPECTED_TASKS,
        "all_154_task_pairs_are_unique_and_exactly_canonical": exact_unique_task_scope,
        "all_tasks_pass_or_reused": all_recorded_pass,
        "all_datasets_present": set(aggregate["by_dataset"]) == set(DATASETS),
        "all_tasks_pass": all_recorded_pass,
        "all_autosomes_present": {result["chrom"] for result in ordered} == set(AUTOSOMES),
        "all_child_receipts_schema_1_2": all(
            (result.get("receipt") or {}).get("schema_version") == "1.2.0"
            for result in ordered
        ),
    }
    release_binding = run_contract.get("release_binding")
    if release_binding is not None:
        checks.update({
            "all_154_child_producers_match_release_contract": (
                embedded_extractor_producers_match_release(ordered, release_binding)
            ),
            "terminal_receipt_is_bound_to_authenticated_release_contract": True,
            "runner_extractor_and_lossless_dependency_match_frozen_snapshot": True,
            "frozen_scientific_and_scheduler_parameters_match": True,
            "resource_gate_bypass_was_forbidden": True,
        })
    return common | {
        "schema_name": "intersubmod.m2_full_extraction_receipt",
        "checks": checks,
        "all_pass": all(checks.values()),
    }


def _bind_terminal_results_to_reusable_children(
    receipt_results: Any,
    reusable_results: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]] | None:
    """Require one persisted terminal row for every canonical validated child.

    Runtime-only PASS diagnostics cannot be reproduced during resume, but the
    scientific child receipt and orchestration completion can and must match
    exactly.  The persisted status may be PASS or REUSED_PASS; both mean the
    same already-validated child for aggregation purposes.
    """
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
        if (
            reusable is None
            or row.get("status") not in {"PASS", "REUSED_PASS"}
            or row.get("receipt") != reusable.get("receipt")
            or row.get("orchestration_completion")
            != reusable.get("orchestration_completion")
        ):
            return None
        bound.append(row)
    if observed_task_ids != expected_task_ids:
        return None
    return bound


def validated_existing_final(
    path: Path,
    run_contract: dict[str, Any],
    *,
    reusable_results: Sequence[Mapping[str, Any]],
    thresholds: tuple[int, ...],
) -> dict[str, Any] | None:
    if not path.is_file() or not verify_sha256_sidecar(path):
        return None
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if (
        receipt.get("schema_name") != "intersubmod.m2_full_extraction_receipt"
        or receipt.get("all_pass") is not True
        or receipt.get("n_results") != EXPECTED_TASKS
        or receipt.get("run_contract") != run_contract
        or receipt.get("receipt_integrity") != _receipt_integrity(path)
    ):
        return None
    bound_results = _bind_terminal_results_to_reusable_children(
        receipt.get("results"), reusable_results
    )
    invocation = receipt.get("invocation")
    elapsed_seconds = receipt.get("elapsed_seconds")
    if (
        bound_results is None
        or not isinstance(invocation, dict)
        or not isinstance(elapsed_seconds, (int, float))
        or isinstance(elapsed_seconds, bool)
        or not math.isfinite(float(elapsed_seconds))
        or float(elapsed_seconds) < 0.0
    ):
        return None
    rebuilt = build_extraction_receipt(
        bound_results,
        thresholds,
        run_contract,
        invocation,
        elapsed_seconds,
    )
    rebuilt["receipt_integrity"] = _receipt_integrity(path)
    expected_has_orchestration = all(
        result.get("orchestration_completion") is not None
        for result in reusable_results
    )
    if expected_has_orchestration:
        if not isinstance(receipt.get("orchestration"), dict):
            return None
        rebuilt["orchestration"] = receipt["orchestration"]
    elif "orchestration" in receipt:
        return None
    return receipt if receipt == rebuilt else None


def ensure_preflight_contract(
    path: Path, preflight: dict[str, Any], run_contract: dict[str, Any]
) -> None:
    """Persist the first preflight and reject any resume contract drift."""
    if path.exists():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing.get("run_contract") != run_contract:
            raise RuntimeError("existing output root was created under a different M2 run contract")
        return
    path.write_text(json.dumps(preflight, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--release-contract-manifest", type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--samtools-threads", type=int, default=1)
    parser.add_argument("--mapq-min", type=int, default=20)
    parser.add_argument("--baseq-min", type=int, default=20)
    parser.add_argument("--bridge-thresholds", default="1,2,3,5")
    parser.add_argument("--max-new-tasks", type=int)
    parser.add_argument("--task-timeout-seconds", type=float, default=28800.0)
    parser.add_argument(
        "--task-timeout-grace-seconds", "--timeout-grace-seconds",
        dest="timeout_grace_seconds", type=float, default=300.0,
    )
    parser.add_argument("--ignore-resource-gate", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--preflight-receipt-output", type=Path)
    parser.add_argument("--preflight-pilot-stage", choices=("extraction", "ranking"))
    parser.add_argument("--preflight-pilot-label")
    parser.add_argument("--preflight-pilot-dataset", choices=DATASETS)
    parser.add_argument("--preflight-pilot-chrom", choices=AUTOSOMES)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    thresholds = tuple(sorted({int(value) for value in args.bridge_thresholds.split(",")}))
    if not thresholds or min(thresholds) < 1:
        raise ValueError("bridge thresholds must be positive")
    args.bridge_thresholds = ",".join(str(value) for value in thresholds)
    if args.workers < 1 or args.workers > 2:
        raise ValueError("M2 I/O workers must be 1 or 2")
    if args.max_new_tasks is not None and args.max_new_tasks < 1:
        raise ValueError("max-new-tasks must be a positive integer")
    if (
        not math.isfinite(args.task_timeout_seconds)
        or not math.isfinite(args.timeout_grace_seconds)
        or args.task_timeout_seconds <= 0
        or args.timeout_grace_seconds < 0
    ):
        raise ValueError("task timeout must be positive and grace must be nonnegative")
    if args.samtools_threads < 0 or args.workers * max(1, args.samtools_threads) > 4:
        raise ValueError("M2 decompression concurrency exceeds pre-registered ceiling")
    pilot_gate_values = (
        args.preflight_receipt_output,
        args.preflight_pilot_stage,
        args.preflight_pilot_label,
        args.preflight_pilot_dataset,
        args.preflight_pilot_chrom,
    )
    if any(value is not None for value in pilot_gate_values) and (
        not all(value is not None for value in pilot_gate_values)
        or not args.preflight_only
        or args.ignore_resource_gate
    ):
        raise ValueError(
            "pilot preflight receipt requires all pilot fields, --preflight-only, and no bypass"
        )
    if args.preflight_receipt_output is not None and (
        (
            args.preflight_pilot_stage == "extraction"
            and args.preflight_pilot_label != "extraction"
        )
        or (
            args.preflight_pilot_stage == "ranking"
            and args.preflight_pilot_label not in {"ranking_bootstrap0", "ranking_bootstrap20"}
        )
    ):
        raise ValueError("pilot gate stage/label must be extraction, ranking_bootstrap0, or ranking_bootstrap20")
    release_binding = None
    if args.release_contract_manifest is not None:
        release_binding = load_release_contract_binding(
            args.release_contract_manifest,
            required_sources={
                "full_extraction_runner": RUNNER,
                "extractor": EXTRACTOR,
                "lossless_read_contract": LOSSLESS_READ_CONTRACT,
            },
            input_manifest=args.manifest,
        )
        validate_release_extraction_parameters(release_binding, args, thresholds)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    manifest_samples = manifest.get("samples", [])
    if (
        manifest.get("dataset_count") != 7
        or len(manifest_samples) != 7
        or {row.get("sample") for row in manifest_samples} != set(DATASETS)
    ):
        raise RuntimeError("manifest is not the declared seven-dataset scope")
    conflicts = active_conflicts()
    resource_preview = resource_gate_preview(args.output_root, conflicts)
    manifest_sha256 = sha256_path(args.manifest)
    extractor_sha256 = sha256_path(EXTRACTOR)
    lossless_sha256 = sha256_path(LOSSLESS_READ_CONTRACT)
    runner_sha256 = sha256_path(RUNNER)
    run_contract = {
        "manifest_sha256": manifest_sha256,
        "extractor_sha256": extractor_sha256,
        "extractor": {"path": str(EXTRACTOR.resolve()), "sha256": extractor_sha256},
        "lossless_read_contract": {
            "path": str(LOSSLESS_READ_CONTRACT.resolve()), "sha256": lossless_sha256,
        },
        "runner": {"path": str(RUNNER), "sha256": runner_sha256},
        "workers": args.workers,
        "samtools_threads": args.samtools_threads,
        "mapq_min": args.mapq_min,
        "baseq_min": args.baseq_min,
        "bridge_thresholds": list(thresholds),
        "component_linkage_bases": list(LINKAGE_BASES),
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
    }
    preflight = {
        "schema_name": "intersubmod.m2_full_extraction_preflight",
        "schema_version": "1.2.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "scope": {
            "datasets": list(DATASETS), "chromosomes": list(AUTOSOMES), "n_tasks": 154,
            "n_technical_datasets": 7, "n_biological_samples": 6,
            "technical_replicate_pair": ["HCC1395", "HCC1395_DORADO"],
        },
        "manifest": {"path": str(args.manifest.resolve()), "sha256": manifest_sha256},
        "extractor": {"path": str(EXTRACTOR.resolve()), "sha256": extractor_sha256},
        "runner": {"path": str(RUNNER), "sha256": runner_sha256},
        "parameters": vars(args) | {
            "manifest": str(args.manifest), "output_root": str(args.output_root),
            "release_contract_manifest": (
                None if args.release_contract_manifest is None
                else str(args.release_contract_manifest)
            ),
            "preflight_receipt_output": (
                None if args.preflight_receipt_output is None
                else str(args.preflight_receipt_output)
            ),
        },
        "invocation": invocation,
        "operational_parameters_excluded_from_run_contract": ["max_new_tasks"],
        "run_contract": run_contract,
        "release_binding": release_binding,
        "active_conflict_process_count": conflicts["process_count"],
        "active_conflict_root_count": conflicts["root_count"],
        "active_conflicts": conflicts["representatives"],
        "resource_gate_preview": resource_preview,
        "resource_gate_pass": (
            (conflicts["process_count"] == 0 or args.ignore_resource_gate)
            and resource_preview["disk_pass"]
        ),
    }
    print(json.dumps(preflight, ensure_ascii=False, indent=2), flush=True)
    if args.preflight_only:
        existing_preflight = args.output_root / "preflight.json"
        if existing_preflight.exists():
            existing = json.loads(existing_preflight.read_text(encoding="utf-8"))
            if existing.get("run_contract") != run_contract:
                raise RuntimeError("existing output root was created under a different M2 run contract")
        if not preflight["resource_gate_pass"]:
            raise SystemExit(2)
        if args.preflight_receipt_output is not None:
            args.output_root.mkdir(parents=True, exist_ok=True)
            output_stat = os.lstat(args.output_root)
            target = {
                "task_type": "B_COMPREHENSIVE_VALIDATION_RELEASE_PILOT",
                "dataset": args.preflight_pilot_dataset,
                "chrom": args.preflight_pilot_chrom,
                "gate_label": args.preflight_pilot_label,
                "output_root": {
                    "path": str(args.output_root.resolve()),
                    "st_dev": int(output_stat.st_dev),
                    "st_ino": int(output_stat.st_ino),
                },
                "manifest": {
                    "path": str(args.manifest.resolve()),
                    "sha256": manifest_sha256,
                },
                "release_manifest": (
                    None if release_binding is None
                    else dict(release_binding["release_manifest"])
                ),
            }
            producer = (
                {"path": str(RUNNER.resolve()), "sha256": sha256_path(RUNNER)}
                if release_binding is None
                else release_producer_sources(release_binding)["runner"]
            )
            _, gate_identity = create_resource_gate_receipt(
                args.output_root,
                stage=args.preflight_pilot_stage,
                gate_scope="pilot",
                target=target,
                producer_source=producer,
                conflicts=active_conflicts(),
                receipt_path=args.preflight_receipt_output,
            )
            print(json.dumps({
                "pilot_resource_gate_receipt": gate_identity,
                "all_pass": True,
            }, ensure_ascii=False), flush=True)
        raise SystemExit(0)
    if not preflight["resource_gate_pass"]:
        raise SystemExit(2)
    release_session = None
    orchestration_state = None
    if release_binding is not None:
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
        }
        session_gate_identity = None
        session_path = args.output_root / "_orchestration" / "session.json"
        if not session_path.exists():
            _, session_gate_identity = create_resource_gate_receipt(
                args.output_root,
                stage="extraction",
                gate_scope="session",
                target=session_target,
                producer_source=release_producer_sources(release_binding)["runner"],
                conflicts=active_conflicts(),
            )
        release_session = ensure_release_orchestration_session(
            args.output_root, release_binding, run_contract, session_gate_identity
        )
        orchestration_state = load_release_orchestration_state(
            args.output_root, release_session
        )
    else:
        args.output_root.mkdir(parents=True, exist_ok=True)
    preflight_path = args.output_root / "preflight.json"
    ensure_preflight_contract(preflight_path, preflight, run_contract)

    specs = []
    expected_child_parameters = {
        "mapq_min": args.mapq_min,
        "baseq_min": args.baseq_min,
        "samtools_threads": args.samtools_threads,
        "bridge_thresholds": list(thresholds),
        "component_linkage_bases": list(LINKAGE_BASES),
    }
    for dataset in DATASETS:
        for chrom in AUTOSOMES:
            task_dir = args.output_root / "samples" / dataset / chrom
            specs.append(
                (
                    dataset,
                    chrom,
                    task_dir,
                    task_command(args.manifest, dataset, chrom, task_dir, args),
                    expected_child_parameters,
                    manifest_sha256,
                    extractor_sha256,
                    args.task_timeout_seconds,
                    args.timeout_grace_seconds,
                )
            )
    reused, pending_specs, invalid = scan_existing_specs(
        specs,
        attested_completions=(
            None if orchestration_state is None else orchestration_state["completions"]
        ),
    )
    if invalid:
        print(json.dumps({"all_pass": False, "existing_nonpass": invalid}, ensure_ascii=False))
        raise SystemExit(1)
    if (
        release_binding is not None
        and not embedded_extractor_producers_match_release(reused, release_binding)
    ):
        raise RuntimeError("a reusable extraction child is not bound to the release code/manifest")
    if orchestration_state is not None and orchestration_state["count"] != len(reused):
        raise RuntimeError("attested task count differs from reusable child count")

    final_path = args.output_root / "full_extraction_receipt.json"
    existing_final = validated_existing_final(
        final_path,
        run_contract,
        reusable_results=reused,
        thresholds=thresholds,
    )
    if not pending_specs and existing_final is not None:
        print(json.dumps({"receipt": str(final_path), "all_pass": True, "status": "REUSED_FINAL"}))
        raise SystemExit(0)
    if final_path.exists():
        raise RuntimeError("existing terminal full extraction receipt is invalid or scope is no longer complete")

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
            stage="extraction",
            gate_scope="batch",
            batch_index=orchestration_state["next_batch_index"],
            target=batch_target,
            producer_source=release_session["producer_sources"]["runner"],
            conflicts=active_conflicts(),
        )
        release_batch, release_grants = create_batch_start_and_grants(
            args.output_root,
            release_session,
            run_contract,
            selected,
            before_count=len(reused),
            previous_chain_head=orchestration_state["previous_chain_head"],
            batch_index=orchestration_state["next_batch_index"],
            max_new_tasks=int(args.max_new_tasks),
            gate=batch_gate_identity,
        )
    started = time.monotonic()

    def progress(result: dict[str, Any]) -> None:
        current = len(reused) + sum(
            item["status"] in {"PASS", "REUSED_PASS"} for item in new_results
        )
        print(
            f"M2_PROGRESS task={result['dataset']}:{result['chrom']} status={result['status']} "
            f"passing={current}/{EXPECTED_TASKS} elapsed_seconds={time.monotonic() - started:.1f}",
            flush=True,
        )

    new_results: list[dict[str, Any]] = []

    def collect_progress(result: dict[str, Any]) -> None:
        new_results.append(result)
        progress(result)

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
                args.output_root,
                release_session,
                release_batch,
                release_grants[task_id],
                result,
            )
    if sha256_path(RUNNER) != runner_sha256:
        raise RuntimeError("M2 extraction runner source changed during the invocation")

    results = reused + executed
    invocation = invocation | {
        "reused_tasks": len(reused),
        "pending_tasks_before_invocation": len(pending_specs),
        "selected_new_tasks": len(selected),
        "selected_task_ids": [f"{spec[0]}:{spec[1]}" for spec in selected],
        "max_inflight_futures": max_inflight,
    }
    receipt = build_extraction_receipt(
        results, thresholds, run_contract, invocation, time.monotonic() - started
    )
    if release_session is not None:
        receipt["orchestration"] = {
            "session_identity": _session_identity(args.output_root, release_session),
            "batch_start_identity": {
                "path": release_batch["path"], "sha256": release_batch["sha256"],
                "batch_id": release_batch["batch_id"],
                "batch_index": release_batch["batch_index"],
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
                "full_extraction_runner": RUNNER,
                "extractor": EXTRACTOR,
                "lossless_read_contract": LOSSLESS_READ_CONTRACT,
            },
            input_manifest=args.manifest,
            _force_deep_reverification=True,
        )
        if current_binding != release_binding:
            raise RuntimeError("release contract or a frozen source changed during extraction")
    if receipt["remaining_tasks"]:
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
        if (
            sha256_path(RUNNER) != runner_sha256
            or sha256_path(EXTRACTOR) != extractor_sha256
            or sha256_path(LOSSLESS_READ_CONTRACT) != lossless_sha256
            or sha256_path(args.manifest) != manifest_sha256
        ):
            raise RuntimeError("M2 extraction code or manifest changed before checkpoint publication")
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

    receipt["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": f"{final_path.name}.sha256",
        "covers": final_path.name,
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
        or sha256_path(EXTRACTOR) != extractor_sha256
        or sha256_path(LOSSLESS_READ_CONTRACT) != lossless_sha256
        or sha256_path(args.manifest) != manifest_sha256
    ):
        raise RuntimeError("M2 extraction code or manifest changed before final publication")
    if release_session is not None:
        write_immutable_json_exclusive(final_path, receipt)
    else:
        write_json_and_sha256_exclusive(final_path, receipt)
    print(json.dumps({
        "receipt": str(final_path),
        "all_pass": receipt["all_pass"],
        **receipt["task_status_counts"],
    }, ensure_ascii=False))
    raise SystemExit(0)


if __name__ == "__main__":
    main()
