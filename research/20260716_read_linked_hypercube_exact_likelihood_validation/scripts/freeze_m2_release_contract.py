#!/usr/bin/env python3
"""Freeze and independently verify the source/input contract for an M2 release.

Production mode is deliberately narrow: it accepts only the canonical v5
manifest and a passing canonical PRE receipt from ``verify_m2_input_identity``.
It copies an exact source allowlist plus the canonical JSON Schema into a new,
read-only, repo-layout-preserving snapshot and publishes an authenticated run
manifest with no-replace semantics.

The CLI has no test-authority switch.  Synthetic authority can only be supplied
through the private Python API used by the colocated unit tests, and every such
artifact is marked ineligible as validation evidence.
"""

from __future__ import annotations

import argparse
import ctypes
import dataclasses
import datetime as dt
import errno
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence


CANONICAL_MANIFEST_SHA256 = "16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45"
CANONICAL_SCHEMA_SHA256 = "47f54d86159beae66b68719e034c7755eba59b2d1719d4546153fa0ae74ff19f"
CANONICAL_SCHEMA_RELATIVE = Path(
    "docs/methodology/_assets/20260627_subclone_4axis_teaching/"
    "schemas/layered_input_manifest_v3.schema.json"
)
RESEARCH_RELATIVE = Path("research/20260716_read_linked_hypercube_exact_likelihood_validation")
SCRIPTS_RELATIVE = RESEARCH_RELATIVE / "scripts"

SOURCE_ALLOWLIST: tuple[tuple[str, Path], ...] = (
    ("extractor", SCRIPTS_RELATIVE / "extract_lossless_read_linkage.py"),
    ("lossless_read_contract", SCRIPTS_RELATIVE / "lossless_read_contract.py"),
    ("full_extraction_runner", SCRIPTS_RELATIVE / "run_full_m2_extraction.py"),
    ("ranker", SCRIPTS_RELATIVE / "build_m2_patterns_and_rank.py"),
    ("hypercube_solver", SCRIPTS_RELATIVE / "hypercube_exact.py"),
    ("full_ranking_runner", SCRIPTS_RELATIVE / "run_full_m2_ranking.py"),
    ("pilot_verifier", SCRIPTS_RELATIVE / "verify_m2_single_task_pilot.py"),
    ("full_verifier", SCRIPTS_RELATIVE / "verify_full_m2_receipts.py"),
    ("input_identity_verifier", SCRIPTS_RELATIVE / "verify_m2_input_identity.py"),
    ("release_contract_freezer", SCRIPTS_RELATIVE / "freeze_m2_release_contract.py"),
)
CANONICAL_SCHEMA_ROLE = "canonical_manifest_schema"
ALL_FROZEN_ROLES = SOURCE_ALLOWLIST + ((CANONICAL_SCHEMA_ROLE, CANONICAL_SCHEMA_RELATIVE),)

PRE_SCHEMA_NAME = "intersubmod.m2_input_identity_verification"
PRE_SCHEMA_VERSION = "1.0.0"
RUN_SCHEMA_NAME = "intersubmod.m2_release_run_manifest"
RUN_SCHEMA_VERSION = "1.0.0"
VERIFY_SCHEMA_NAME = "intersubmod.m2_release_contract_verification"
VERIFY_SCHEMA_VERSION = "1.0.0"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
DATASETS = (
    "COLO829", "H1437", "H2009", "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
)
DIRECT_ROLES = (
    "alignment_bam", "alignment_bai", "tree_vcf", "tree_vcf_index",
    "read_tags_sidecar", "read_tags_index",
)
PRE_CHECK_NAMES = (
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
)
PRE_ASSURANCE = {
    "bam": "storage_identity_v1 metadata plus three manifest-defined sampled chunks; no full BAM hash",
    "other_direct_inputs": "full SHA-256",
    "pre_post": (
        "canonical boundary snapshot exact equality; this detects persistent drift but does not prove "
        "that no transient mutation was made and restored between PRE and POST"
    ),
    "temporal_immutability_proven": False,
}
RUNTIME_KEYS = {"python", "packages", "samtools", "platform"}
RUN_CHECK_NAMES = (
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
)

FROZEN_PARAMETERS: Mapping[str, Any] = {
    "extraction": {
        "mapq_min": 20,
        "baseq_min": 20,
        "bridge_thresholds": [1, 2, 3, 5],
        "workers": 2,
        "samtools_threads": 1,
    },
    "ranking": {
        "thresholds": [1, 2, 3, 5],
        "component_bases": ["PS_HP1", "PS_HP2"],
        "hp_families": ["1", "2"],
        "structural_exact_pattern_minread_grid": [1, 2, 3, 5],
        "primary_structural_exact_pattern_minread": 3,
        "exact_k_max": 12,
        "max_vertex_sets": 256,
        "solver_time_limit_seconds_per_milp": 30.0,
        "fixed_error_grid": [0.005, 0.01, 0.02, 0.05],
        "minimum_bq_error_rate": 0.000001,
        "maximum_bq_error_rate": 0.25,
        "conditional_candidate_ranking_bootstrap_replicates": 20,
        "conditional_candidate_ranking_bootstrap_seed": 20260716,
        "tie_tolerance_log_likelihood": 0.000001,
        "workers": 2,
    },
    "scheduler": {
        "initial_batch_tasks": 8,
        "subsequent_batch_tasks": 16,
        "task_timeout_seconds": 28800,
        "timeout_grace_seconds": 300,
    },
}


class ReleaseContractError(RuntimeError):
    """Fail-closed release-contract, identity, or filesystem violation."""


@dataclasses.dataclass(frozen=True)
class _Authority:
    mode: str
    validation_evidence_eligible: bool
    expected_manifest_sha256: str
    schema_relative_path: Path
    expected_schema_sha256: str
    expected_pre_authority_mode: str


PRODUCTION_AUTHORITY = _Authority(
    mode="CANONICAL_V5_FROZEN",
    validation_evidence_eligible=True,
    expected_manifest_sha256=CANONICAL_MANIFEST_SHA256,
    schema_relative_path=CANONICAL_SCHEMA_RELATIVE,
    expected_schema_sha256=CANONICAL_SCHEMA_SHA256,
    expected_pre_authority_mode="CANONICAL_V5_FROZEN",
)


def _synthetic_test_authority(
    manifest_sha256: str, schema_relative_path: Path, schema_sha256: str
) -> _Authority:
    """Create an explicitly non-evidentiary authority for isolated unit tests."""
    return _Authority(
        mode="SYNTHETIC_TEST_ONLY_NOT_VALIDATION_EVIDENCE",
        validation_evidence_eligible=False,
        expected_manifest_sha256=manifest_sha256,
        schema_relative_path=schema_relative_path,
        expected_schema_sha256=schema_sha256,
        expected_pre_authority_mode="TEST_ONLY_UNFROZEN",
    )


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReleaseContractError(message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return sha256_bytes(payload)


def _strict_int(value: Any, label: str, minimum: int = 0) -> int:
    require(
        isinstance(value, int) and not isinstance(value, bool) and value >= minimum,
        f"{label} must be an integer >= {minimum}",
    )
    return int(value)


def _exact_keys(value: Any, expected: set[str], label: str) -> Mapping[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    require(set(value) == expected, f"{label} schema mismatch: {sorted(set(value) ^ expected)}")
    return value


def _validate_full_identity(value: Any, label: str) -> dict[str, Any]:
    identity = _exact_keys(value, {"policy", "size_bytes", "sha256"}, label)
    require(identity["policy"] == "full_sha256", f"{label} policy mismatch")
    _strict_int(identity["size_bytes"], f"{label}.size_bytes")
    require(
        isinstance(identity["sha256"], str) and SHA256_RE.fullmatch(identity["sha256"]) is not None,
        f"{label}.sha256 is invalid",
    )
    return dict(identity)


def _path_components_are_real_dirs(repo_root: Path, source: Path, label: str) -> None:
    """Reject any symlink/non-directory component between repo root and source."""
    root = repo_root.absolute()
    candidate = source.absolute()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ReleaseContractError(f"{label} is outside repository root: {candidate}") from exc
    current = root
    root_state = os.lstat(current)
    require(stat.S_ISDIR(root_state.st_mode) and not stat.S_ISLNK(root_state.st_mode), f"repo root is not a real directory: {root}")
    relative = candidate.relative_to(root)
    for component in relative.parts[:-1]:
        current = current / component
        observed = os.lstat(current)
        require(
            stat.S_ISDIR(observed.st_mode) and not stat.S_ISLNK(observed.st_mode),
            f"{label} has symlink/non-directory parent component: {current}",
        )
    resolved_root = root.resolve(strict=True)
    resolved_source = candidate.resolve(strict=True)
    try:
        resolved_source.relative_to(resolved_root)
    except ValueError as exc:
        raise ReleaseContractError(f"{label} symlink escape resolves outside repository root: {candidate}") from exc


def _reject_duplicate_key(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ReleaseContractError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def strict_json_load_bytes(payload: bytes, label: Path | str) -> Any:
    try:
        return json.loads(
            payload.decode("utf-8", errors="strict"),
            object_pairs_hook=_reject_duplicate_key,
            parse_constant=lambda value: (_ for _ in ()).throw(
                ReleaseContractError(f"non-finite JSON number in {label}: {value}")
            ),
        )
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ReleaseContractError(f"invalid JSON {label}: {exc}") from exc


def _stat_identity(value: os.stat_result) -> dict[str, int]:
    return {
        "st_dev": int(value.st_dev),
        "st_ino": int(value.st_ino),
        "st_nlink": int(value.st_nlink),
        "size_bytes": int(value.st_size),
        "mtime_ns": int(value.st_mtime_ns),
        "ctime_ns": int(value.st_ctime_ns),
    }


def _same_file_state(left: os.stat_result, right: os.stat_result) -> bool:
    return _stat_identity(left) == _stat_identity(right)


def _stable_regular_file(path: Path, label: str) -> tuple[bytes, dict[str, Any]]:
    """Read a regular non-symlink file while detecting path/inode mutation."""
    absolute = path.absolute()
    try:
        before_path = os.lstat(absolute)
    except OSError as exc:
        raise ReleaseContractError(f"cannot lstat {label}: {absolute}: {exc}") from exc
    require(stat.S_ISREG(before_path.st_mode), f"{label} is not a regular file: {absolute}")
    require(not stat.S_ISLNK(before_path.st_mode), f"{label} is a symlink: {absolute}")
    require(before_path.st_nlink == 1, f"{label} hardlink count is not one: {absolute}")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(absolute, flags)
    except OSError as exc:
        raise ReleaseContractError(f"cannot safely open {label}: {absolute}: {exc}") from exc
    try:
        before_fd = os.fstat(descriptor)
        require(before_fd.st_nlink == 1, f"{label} open file hardlink count is not one: {absolute}")
        require(_same_file_state(before_path, before_fd), f"{label} path changed before read: {absolute}")
        pieces: list[bytes] = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            pieces.append(block)
        after_fd = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        after_path = os.lstat(absolute)
    except OSError as exc:
        raise ReleaseContractError(f"cannot lstat {label} after read: {absolute}: {exc}") from exc
    require(
        _same_file_state(before_fd, after_fd) and _same_file_state(before_fd, after_path),
        f"{label} changed while being read: {absolute}",
    )
    require(after_fd.st_nlink == 1 and after_path.st_nlink == 1, f"{label} gained a hardlink: {absolute}")
    payload = b"".join(pieces)
    require(len(payload) == before_fd.st_size, f"short read for {label}: {absolute}")
    return payload, {
        "path": str(absolute.resolve(strict=True)),
        **_stat_identity(before_fd),
        "mode_octal": format(stat.S_IMODE(before_fd.st_mode), "04o"),
        "sha256": sha256_bytes(payload),
    }


def _authenticate_json_receipt(path: Path, label: str) -> tuple[dict[str, Any], dict[str, Any]]:
    payload, receipt_identity = _stable_regular_file(path, label)
    sidecar_path = path.with_name(f"{path.name}.sha256")
    sidecar_payload, sidecar_identity = _stable_regular_file(sidecar_path, f"{label} sidecar")
    require(
        (receipt_identity["st_dev"], receipt_identity["st_ino"])
        != (sidecar_identity["st_dev"], sidecar_identity["st_ino"]),
        f"{label} aliases its sidecar",
    )
    try:
        fields = sidecar_payload.decode("ascii", errors="strict").strip().split()
    except UnicodeError as exc:
        raise ReleaseContractError(f"{label} sidecar is not ASCII") from exc
    require(
        len(fields) == 2
        and fields[0] == receipt_identity["sha256"]
        and fields[1] == path.name,
        f"{label} sidecar mismatch",
    )
    document = strict_json_load_bytes(payload, path)
    require(isinstance(document, dict), f"{label} root is not an object")
    return document, {
        **receipt_identity,
        "sidecar": sidecar_identity,
    }


def _validate_storage_identity(
    value: Any, bam_path: Path, bai_path: Path, label: str
) -> dict[str, Any]:
    required = {
        "policy", "assurance", "is_full_content_hash", "requested_path", "realpath",
        "logical_is_symlink", "logical_size_bytes", "logical_mtime_ns", "st_dev", "st_ino",
        "size_bytes", "mtime_ns", "ctime_ns", "chunk_size_bytes", "chunks", "index",
        "identity_sha256",
    }
    identity = _exact_keys(value, required, label)
    require(identity["policy"] == "storage_identity_v1", f"{label} policy mismatch")
    require(
        identity["assurance"] == "metadata_plus_sampled_chunks_not_full_content_hash"
        and identity["is_full_content_hash"] is False,
        f"{label} assurance mismatch",
    )
    require(identity["requested_path"] == str(bam_path), f"{label}.requested_path mismatch")
    logical = os.lstat(bam_path)
    resolved = bam_path.resolve(strict=True)
    target = resolved.stat()
    require(stat.S_ISREG(target.st_mode), f"{label} BAM target is not regular")
    require(identity["realpath"] == str(resolved), f"{label}.realpath mismatch")
    require(identity["logical_is_symlink"] is stat.S_ISLNK(logical.st_mode), f"{label}.logical_is_symlink mismatch")
    observed_numbers = {
        "logical_size_bytes": logical.st_size,
        "logical_mtime_ns": logical.st_mtime_ns,
        "st_dev": target.st_dev,
        "st_ino": target.st_ino,
        "size_bytes": target.st_size,
        "mtime_ns": target.st_mtime_ns,
        "ctime_ns": target.st_ctime_ns,
    }
    for key, observed in observed_numbers.items():
        require(_strict_int(identity[key], f"{label}.{key}") == int(observed), f"{label}.{key} mismatch")
    require(_strict_int(identity["size_bytes"], f"{label}.size_bytes", 1) > 0, f"{label} BAM is empty")
    require(identity["chunk_size_bytes"] == 1024 * 1024, f"{label}.chunk_size_bytes mismatch")
    chunks = identity["chunks"]
    require(isinstance(chunks, list) and len(chunks) == 3, f"{label} chunk schema mismatch")
    length = min(1024 * 1024, int(target.st_size))
    offsets = (0, max(0, (int(target.st_size) - length) // 2), max(0, int(target.st_size) - length))
    for index, (chunk, expected_name, expected_offset) in enumerate(zip(chunks, ("first", "middle", "last"), offsets)):
        row = _exact_keys(chunk, {"label", "offset", "length", "sha256"}, f"{label}.chunks[{index}]")
        require(row["label"] == expected_name and row["offset"] == expected_offset and row["length"] == length, f"{label}.chunks[{index}] geometry mismatch")
        require(isinstance(row["sha256"], str) and SHA256_RE.fullmatch(row["sha256"]) is not None, f"{label}.chunks[{index}] SHA invalid")
    index = _exact_keys(identity["index"], {"path", "identity"}, f"{label}.index")
    require(index["path"] == str(bai_path), f"{label}.index.path mismatch")
    _validate_full_identity(index["identity"], f"{label}.index.identity")
    require(isinstance(identity["identity_sha256"], str) and SHA256_RE.fullmatch(identity["identity_sha256"]) is not None, f"{label}.identity_sha256 invalid")
    without_digest = dict(identity)
    digest = without_digest.pop("identity_sha256")
    require(canonical_json_sha256(without_digest) == digest, f"{label} semantic digest mismatch")
    return dict(identity)


def _derive_manifest_role_contract(manifest: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    samples = manifest.get("samples")
    require(isinstance(samples, list) and len(samples) == 7, "manifest samples are not seven rows")
    by_dataset: dict[str, Mapping[str, Any]] = {}
    for sample in samples:
        require(isinstance(sample, dict), "manifest sample row is not an object")
        dataset = sample.get("sample")
        require(isinstance(dataset, str) and dataset not in by_dataset, "manifest sample is duplicate/invalid")
        by_dataset[dataset] = sample
    require(set(by_dataset) == set(DATASETS), "manifest dataset set mismatch")
    artifacts: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    for dataset in DATASETS:
        sample = by_dataset[dataset]
        expected_biological = "HCC1395" if dataset == "HCC1395_DORADO" else dataset
        expected_platform = "ONT_DORADO" if dataset == "HCC1395_DORADO" else "ONT"
        expected_replicate = "platform_replica" if dataset == "HCC1395_DORADO" else "canonical"
        require(sample.get("biological_id") == expected_biological, f"{dataset} biological_id mismatch")
        require(sample.get("platform") == expected_platform, f"{dataset} platform mismatch")
        require(sample.get("replicate_role") == expected_replicate, f"{dataset} replicate_role mismatch")
        alignment = sample.get("alignment_payload")
        somatic = sample.get("somatic")
        read_tags = sample.get("read_tags")
        require(isinstance(alignment, dict) and isinstance(somatic, dict) and isinstance(read_tags, dict), f"{dataset} input blocks missing")
        require(alignment.get("identity_policy") == "storage_identity_v1" and alignment.get("embedded_hp_ps_policy") == "ignore", f"{dataset} alignment policy mismatch")
        require(somatic.get("backbone_role") == "longphase_s_recalibrated_filter_pass", f"{dataset} somatic backbone mismatch")
        require(read_tags.get("authority") == "external_coordinate_sidecar" and read_tags.get("identity_schema") == "coordinate_join_v1" and read_tags.get("fallback_policy") == "forbidden", f"{dataset} read-tag policy mismatch")
        bam = Path(str(alignment.get("path", "")))
        bai = Path(str(alignment.get("index_path", "")))
        storage = _validate_storage_identity(alignment.get("storage_identity_v1"), bam, bai, f"{dataset}.alignment_payload")
        tree = somatic.get("tree_vcf")
        sidecar = read_tags.get("sidecar")
        sidecar_index = read_tags.get("index")
        require(isinstance(tree, dict) and isinstance(sidecar, dict) and isinstance(sidecar_index, dict), f"{dataset} direct roles missing")
        tree_index = tree.get("index")
        require(isinstance(tree_index, dict), f"{dataset} tree index missing")
        role_specs = (
            ("alignment_bam", bam, storage, "storage_identity_v1"),
            ("alignment_bai", bai, storage["index"]["identity"], "full_sha256"),
            ("tree_vcf", Path(str(tree.get("path", ""))), _validate_full_identity(tree.get("identity"), f"{dataset}.tree_vcf"), "full_sha256"),
            ("tree_vcf_index", Path(str(tree_index.get("path", ""))), _validate_full_identity(tree_index.get("identity"), f"{dataset}.tree_vcf_index"), "full_sha256"),
            ("read_tags_sidecar", Path(str(sidecar.get("path", ""))), _validate_full_identity(sidecar.get("identity"), f"{dataset}.read_tags_sidecar"), "full_sha256"),
            ("read_tags_index", Path(str(sidecar_index.get("path", ""))), _validate_full_identity(sidecar_index.get("identity"), f"{dataset}.read_tags_index"), "full_sha256"),
        )
        for role, path, expected, policy in role_specs:
            resolved = path.resolve(strict=True)
            observed = resolved.stat()
            require(stat.S_ISREG(observed.st_mode), f"{dataset}.{role} is not a regular file")
            artifacts.append({"dataset": dataset, "role": role, "policy": policy, "path": str(path), "expected": expected})
            inventory.append({
                "dataset": dataset, "role": role, "path": str(path), "resolved_path": str(resolved),
                "st_dev": int(observed.st_dev), "st_ino": int(observed.st_ino),
            })
    require(len(artifacts) == 42 and len(inventory) == 42, "manifest role inventory is not 42")
    require(len({(row["st_dev"], row["st_ino"]) for row in inventory}) == 42, "manifest direct roles alias files")
    return artifacts, inventory


def _validate_pre_document(
    receipt: Mapping[str, Any], manifest_document: Mapping[str, Any], expected_manifest_path: Path,
    authority: _Authority, expected_schema_path: Path, expected_verifier_path: Path,
    expected_verifier_sha256: str,
) -> dict[str, Any]:
    top_keys = {
        "schema_name", "schema_version", "task_type", "mode", "authority_mode",
        "validation_evidence_eligible", "authority", "scope", "manifest", "canonical_schema",
        "verifier", "assurance", "input_identity_snapshot", "input_identity_snapshot_sha256",
        "artifact_audit", "compare_to", "checks", "all_pass", "receipt_integrity",
    }
    _exact_keys(receipt, top_keys, "PRE receipt")
    require(receipt["schema_name"] == PRE_SCHEMA_NAME and receipt["schema_version"] == PRE_SCHEMA_VERSION, "PRE schema mismatch")
    require(receipt["task_type"] == "B_COMPREHENSIVE_VALIDATION" and receipt["mode"] == "PRE", "PRE task/mode mismatch")
    require(receipt["authority_mode"] == authority.expected_pre_authority_mode, "PRE authority_mode mismatch")
    require(receipt["validation_evidence_eligible"] is authority.validation_evidence_eligible, "PRE eligibility mismatch")
    authority_row = _exact_keys(receipt["authority"], {"canonical_manifest_sha256", "canonical_schema_sha256", "canonical_schema_path", "selected_authority_is_canonical", "test_only_override"}, "PRE authority")
    require(authority_row["canonical_manifest_sha256"] == CANONICAL_MANIFEST_SHA256 and authority_row["canonical_schema_sha256"] == CANONICAL_SCHEMA_SHA256, "PRE frozen canonical authority declaration mismatch")
    require(Path(str(authority_row["canonical_schema_path"])).absolute() == expected_schema_path.absolute(), "PRE authority schema path mismatch")
    require(authority_row["selected_authority_is_canonical"] is authority.validation_evidence_eligible and authority_row["test_only_override"] is (not authority.validation_evidence_eligible), "PRE authority booleans mismatch")
    require(receipt["scope"] == {"technical_datasets": 7, "biological_samples": 6, "chromosomes": "chr1-chr22", "tasks": 154, "datasets": list(DATASETS), "direct_input_artifacts": 42}, "PRE scope mismatch")
    manifest = _exact_keys(receipt["manifest"], {"path", "sha256", "expected_sha256"}, "PRE manifest")
    require(Path(str(manifest["path"])).absolute() == expected_manifest_path.absolute(), "PRE manifest path mismatch")
    require(manifest["sha256"] == manifest["expected_sha256"] == authority.expected_manifest_sha256, "PRE manifest SHA mismatch")
    schema = _exact_keys(receipt["canonical_schema"], {"path", "sha256", "expected_sha256"}, "PRE canonical_schema")
    require(Path(str(schema["path"])).absolute() == expected_schema_path.absolute(), "PRE canonical schema path mismatch")
    require(schema["sha256"] == schema["expected_sha256"] == authority.expected_schema_sha256, "PRE canonical schema SHA mismatch")
    verifier = _exact_keys(receipt["verifier"], {"path", "sha256"}, "PRE verifier")
    require(Path(str(verifier["path"])).absolute() == expected_verifier_path.absolute(), "PRE verifier source path mismatch")
    require(verifier["sha256"] == expected_verifier_sha256, "PRE verifier source SHA mismatch")
    require(receipt["assurance"] == PRE_ASSURANCE, "PRE assurance mismatch")
    require(receipt["compare_to"] is None, "PRE compare_to must be null")
    checks = _exact_keys(receipt["checks"], set(PRE_CHECK_NAMES), "PRE checks")
    require(all(checks[name] is True for name in PRE_CHECK_NAMES) and receipt["all_pass"] is True, "PRE checks did not all pass")
    integrity = _exact_keys(receipt["receipt_integrity"], {"scheme", "sidecar_name", "covers"}, "PRE receipt_integrity")
    require(integrity["scheme"] == "external_sha256_sidecar_v1" and integrity["sidecar_name"] == f"{integrity['covers']}.sha256" and isinstance(integrity["covers"], str) and integrity["covers"], "PRE receipt_integrity mismatch")

    derived, inventory = _derive_manifest_role_contract(manifest_document)
    expected_artifacts = [
        {**row, "observed": row["expected"], "match": True}
        for row in derived
    ]
    expected_audit = {
        "artifacts": expected_artifacts,
        "role_inventory": inventory,
        "n_artifacts": 42,
        "n_unique_resolved_files": 42,
        "n_storage_identity_v1": 7,
        "n_full_sha256": 35,
        "n_sampled_bam_chunks": 21,
        "n_mismatches": 0,
    }
    require(receipt["artifact_audit"] == expected_audit, "PRE artifact audit differs from manifest-derived 7x6 contract")
    expected_snapshot = {
        "manifest_sha256": authority.expected_manifest_sha256,
        "schema_sha256": authority.expected_schema_sha256,
        "datasets": list(DATASETS),
        "artifacts": [
            {"dataset": row["dataset"], "role": row["role"], "policy": row["policy"], "path": row["path"], "observed": row["expected"]}
            for row in derived
        ],
    }
    require(receipt["input_identity_snapshot"] == expected_snapshot, "PRE snapshot differs from manifest-derived snapshot")
    require(receipt["input_identity_snapshot_sha256"] == canonical_json_sha256(expected_snapshot), "PRE snapshot semantic digest mismatch")
    return expected_snapshot


def _validate_pre_receipt(
    pre_receipt_path: Path, manifest_path: Path, manifest_document: Mapping[str, Any], authority: _Authority,
    expected_schema_path: Path, expected_verifier_path: Path, expected_verifier_sha256: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    receipt, identity = _authenticate_json_receipt(pre_receipt_path, "PRE input identity receipt")
    _validate_pre_document(
        receipt, manifest_document, manifest_path, authority, expected_schema_path,
        expected_verifier_path, expected_verifier_sha256,
    )
    require(receipt["receipt_integrity"]["covers"] == pre_receipt_path.name, "PRE receipt filename binding mismatch")
    return receipt, identity


def _discover_repo_root() -> Path:
    script = Path(__file__).resolve()
    for parent in script.parents:
        if (parent / CANONICAL_SCHEMA_RELATIVE).is_file():
            return parent
    raise ReleaseContractError("cannot discover InterSubMod repository root")


def _runtime_identity() -> dict[str, Any]:
    packages = {}
    for name in ("numpy", "scipy", "pysam"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError as exc:
            raise ReleaseContractError(f"required Python package is unavailable: {name}") from exc
    samtools = shutil.which("samtools")
    require(samtools is not None, "samtools is unavailable on PATH")
    try:
        result = subprocess.run(
            [samtools, "--version"], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=10, check=True,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise ReleaseContractError(f"cannot obtain samtools runtime version: {exc}") from exc
    version_lines = result.stdout.splitlines()
    require(version_lines, "samtools --version returned no output")
    return {
        "python": {
            "executable": str(Path(sys.executable).resolve(strict=True)),
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
        },
        "packages": packages,
        "samtools": {
            "path": str(Path(samtools).resolve(strict=True)),
            "version_line": version_lines[0],
            "htslib_version_line": next((line for line in version_lines if line.startswith("Using htslib")), None),
        },
        "platform": platform.platform(),
    }


def _validate_runtime_identity(value: Any, label: str) -> dict[str, Any]:
    runtime = _exact_keys(value, RUNTIME_KEYS, label)
    python = _exact_keys(runtime["python"], {"executable", "version", "implementation"}, f"{label}.python")
    packages = _exact_keys(runtime["packages"], {"numpy", "scipy", "pysam"}, f"{label}.packages")
    samtools = _exact_keys(runtime["samtools"], {"path", "version_line", "htslib_version_line"}, f"{label}.samtools")
    for name, item in (
        ("python.executable", python["executable"]), ("python.version", python["version"]),
        ("python.implementation", python["implementation"]), ("platform", runtime["platform"]),
        *((f"packages.{key}", packages[key]) for key in ("numpy", "scipy", "pysam")),
        *((f"samtools.{key}", samtools[key]) for key in ("path", "version_line", "htslib_version_line")),
    ):
        require(isinstance(item, str) and item, f"{label}.{name} must be a non-empty string")
    return {
        "python": dict(python), "packages": dict(packages), "samtools": dict(samtools),
        "platform": runtime["platform"],
    }


def _write_new_file(path: Path, payload: bytes, mode: int) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0)
    descriptor = os.open(path, flags, mode)
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.fchmod(descriptor, mode)
        observed = os.fstat(descriptor)
        require(observed.st_nlink == 1, f"new file hardlink count is not one: {path}")
    finally:
        os.close(descriptor)
    return {
        "path": str(path),
        "size_bytes": len(payload),
        "sha256": sha256_bytes(payload),
        "mode_octal": format(stat.S_IMODE(observed.st_mode), "04o"),
        "st_dev": int(observed.st_dev),
        "st_ino": int(observed.st_ino),
        "st_nlink": int(observed.st_nlink),
        "mtime_ns": int(observed.st_mtime_ns),
        "ctime_ns": int(observed.st_ctime_ns),
    }


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _set_snapshot_directory_modes(root: Path) -> None:
    directories = [path for path in root.rglob("*") if path.is_dir()]
    for directory in sorted(directories, key=lambda value: len(value.parts), reverse=True):
        os.chmod(directory, 0o555)
        _fsync_directory(directory)
    os.chmod(root, 0o555)
    _fsync_directory(root)


def _rename_noreplace(source: Path, destination: Path) -> None:
    """Atomically publish a directory and never replace an existing path."""
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    require(renameat2 is not None, "renameat2(RENAME_NOREPLACE) is unavailable")
    renameat2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    renameat2.restype = ctypes.c_int
    result = renameat2(
        -100, os.fsencode(source), -100, os.fsencode(destination), 1  # AT_FDCWD, RENAME_NOREPLACE
    )
    if result != 0:
        observed_errno = ctypes.get_errno()
        if observed_errno in {errno.EEXIST, errno.ENOTEMPTY}:
            raise ReleaseContractError(f"refusing to overwrite release contract root: {destination}")
        raise ReleaseContractError(
            f"atomic release-contract publish failed: {os.strerror(observed_errno)}"
        )


def _failed_archive_path(parent: Path, kind: str) -> Path:
    """Return a unique, human-auditable failed-publication archive path."""
    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    return parent / f".failed-{kind}.{timestamp}.{os.getpid()}"


def _archive_failed_tree(source: Path, kind: str) -> Path | None:
    """Atomically preserve one failed tree; leave it in place if archival fails."""
    if not os.path.lexists(source):
        return None
    destination = _failed_archive_path(source.parent, kind)
    try:
        _rename_noreplace(source, destination)
        _fsync_directory(source.parent)
    except (ReleaseContractError, OSError):
        return None
    return destination


def _archive_failed_publication(paths: Sequence[Path], parent: Path) -> Path | None:
    """Preserve owned receipt publication remnants without deleting any file."""
    owned = [path for path in paths if os.path.lexists(path)]
    if not owned:
        return None
    archive = _failed_archive_path(parent, "publication")
    try:
        os.mkdir(archive, 0o700)
        _fsync_directory(parent)
    except OSError:
        return None
    for index, source in enumerate(owned):
        destination = archive / f"{index:02d}.{source.name}"
        try:
            _rename_noreplace(source, destination)
        except (ReleaseContractError, OSError):
            # Every object remains either in the audit archive or at its
            # original path.  Never destroy evidence while handling failure.
            return archive
    try:
        _fsync_directory(archive)
        _fsync_directory(parent)
    except OSError:
        return archive
    return archive


def _write_authenticated_exclusive(path: Path, document: dict[str, Any]) -> str:
    path = path.absolute()
    sidecar = path.with_name(f"{path.name}.sha256")
    require(not os.path.lexists(path) and not os.path.lexists(sidecar), f"refusing to overwrite evidence: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    document["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": sidecar.name,
        "covers": path.name,
    }
    payload = json.dumps(document, ensure_ascii=False, allow_nan=False, indent=2).encode("utf-8") + b"\n"
    digest = sha256_bytes(payload)
    receipt_fd, receipt_tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    side_fd, side_tmp_name = tempfile.mkstemp(prefix=f".{sidecar.name}.", dir=str(path.parent))
    receipt_tmp, side_tmp = Path(receipt_tmp_name), Path(side_tmp_name)
    published_receipt = published_sidecar = False
    try:
        with os.fdopen(receipt_fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        with os.fdopen(side_fd, "wb") as handle:
            handle.write(f"{digest}  {path.name}\n".encode("ascii"))
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(receipt_tmp, 0o444)
        os.chmod(side_tmp, 0o444)
        _rename_noreplace(receipt_tmp, path)
        published_receipt = True
        _rename_noreplace(side_tmp, sidecar)
        published_sidecar = True
        _fsync_directory(path.parent)
        _, receipt_identity = _stable_regular_file(path, "published verification receipt")
        side_payload, side_identity = _stable_regular_file(
            sidecar, "published verification receipt sidecar"
        )
        require(receipt_identity["mode_octal"] == "0444", "published receipt mode is not 0444")
        require(side_identity["mode_octal"] == "0444", "published sidecar mode is not 0444")
        require(receipt_identity["sha256"] == digest, "published receipt SHA mismatch")
        require(
            side_payload == f"{digest}  {path.name}\n".encode("ascii"),
            "published receipt sidecar content mismatch",
        )
        require(
            (receipt_identity["st_dev"], receipt_identity["st_ino"])
            != (side_identity["st_dev"], side_identity["st_ino"]),
            "published receipt aliases its sidecar",
        )
    except BaseException as exc:
        remnants: list[Path] = []
        if published_receipt:
            remnants.append(path)
        if published_sidecar:
            remnants.append(sidecar)
        remnants.extend((receipt_tmp, side_tmp))
        archive = _archive_failed_publication(remnants, path.parent)
        if archive is not None and hasattr(exc, "add_note"):
            exc.add_note(f"failed publication evidence archived at {archive}")
        raise
    return digest


def _run_fresh_input_verifier(
    *,
    manifest_path: Path,
    manifest_document: Mapping[str, Any],
    schema_path: Path,
    verifier_path: Path,
    verifier_sha256: str,
    authority: _Authority,
    expected_snapshot: Mapping[str, Any],
    fresh_factory: Any | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Mandatory fresh verifier execution; injection exists only for synthetic unit tests."""
    manifest_path = manifest_path.resolve(strict=True)
    schema_path = schema_path.resolve(strict=True)
    verifier_path = verifier_path.resolve(strict=True)
    if authority.validation_evidence_eligible:
        require(fresh_factory is None, "production fresh verifier injection is forbidden")
        with tempfile.TemporaryDirectory(prefix="m2-release-fresh-input-") as temporary:
            output = Path(temporary) / "fresh_input_identity_pre.json"
            command = [
                sys.executable, str(verifier_path), "--manifest", str(manifest_path),
                "--output", str(output),
            ]
            try:
                result = subprocess.run(
                    command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    timeout=28800, check=False,
                )
            except (OSError, subprocess.SubprocessError) as exc:
                raise ReleaseContractError(f"fresh input verifier could not run: {exc}") from exc
            require(result.returncode == 0, f"fresh input verifier failed ({result.returncode}): {result.stdout[-2000:]}")
            fresh, identity = _authenticate_json_receipt(output, "fresh PRE input identity receipt")
            _validate_pre_document(
                fresh, manifest_document, manifest_path, authority, schema_path,
                verifier_path, verifier_sha256,
            )
            require(fresh["receipt_integrity"]["covers"] == output.name, "fresh PRE filename binding mismatch")
            execution_mode = "production_subprocess_required"
            receipt_sha256 = identity["sha256"]
    else:
        require(fresh_factory is not None, "synthetic tests must inject a private fresh verifier receipt")
        produced = fresh_factory(
            manifest_path=manifest_path,
            manifest_document=manifest_document,
            schema_path=schema_path,
            verifier_path=verifier_path,
            verifier_sha256=verifier_sha256,
            authority=authority,
        )
        require(isinstance(produced, dict), "synthetic fresh verifier factory did not return an object")
        fresh = dict(produced)
        _validate_pre_document(
            fresh, manifest_document, manifest_path, authority, schema_path,
            verifier_path, verifier_sha256,
        )
        execution_mode = "private_synthetic_test_injection_not_validation_evidence"
        receipt_sha256 = canonical_json_sha256(fresh)
    require(fresh["all_pass"] is True, "fresh input verifier did not pass")
    require(
        fresh["validation_evidence_eligible"] is authority.validation_evidence_eligible,
        "fresh input verifier eligibility mismatch",
    )
    require(fresh["input_identity_snapshot"] == expected_snapshot, "fresh input snapshot differs from supplied PRE snapshot")
    summary = {
        "execution_mode": execution_mode,
        "verifier_path": str(verifier_path.resolve(strict=True)),
        "verifier_sha256": verifier_sha256,
        "receipt_sha256": receipt_sha256,
        "receipt_semantic_sha256": canonical_json_sha256(fresh),
        "checks_semantic_sha256": canonical_json_sha256(fresh["checks"]),
        "artifact_audit_semantic_sha256": canonical_json_sha256(fresh["artifact_audit"]),
        "input_identity_snapshot_sha256": fresh["input_identity_snapshot_sha256"],
        "all_pass": True,
        "validation_evidence_eligible": fresh["validation_evidence_eligible"],
        "exactly_equals_supplied_pre_snapshot": True,
    }
    return fresh, summary


def freeze_release_contract(
    manifest_path: Path,
    pre_input_identity_receipt: Path,
    output_contract_root: Path,
    *,
    _repo_root: Path | None = None,
    _authority: _Authority | None = None,
    _runtime: Mapping[str, Any] | None = None,
    _fresh_pre_factory: Any | None = None,
) -> dict[str, Any]:
    authority = PRODUCTION_AUTHORITY if _authority is None else _authority
    if authority.validation_evidence_eligible:
        require(_authority is None and _repo_root is None and _runtime is None, "production authority forbids private test overrides")
    else:
        require(_authority is not None and _repo_root is not None and _runtime is not None, "synthetic authority requires explicit private fixtures")
    repo_root = (_discover_repo_root() if _repo_root is None else _repo_root).absolute()
    output_root = output_contract_root.absolute()
    require(not os.path.lexists(output_root), f"output contract root already exists: {output_root}")
    output_root.parent.mkdir(parents=True, exist_ok=True)

    manifest_payload, manifest_identity = _stable_regular_file(manifest_path, "canonical manifest")
    require(
        manifest_identity["sha256"] == authority.expected_manifest_sha256,
        "manifest SHA does not match frozen authority",
    )
    manifest_document = strict_json_load_bytes(manifest_payload, manifest_path)
    require(isinstance(manifest_document, dict), "manifest root is not an object")
    expected_roles = ALL_FROZEN_ROLES[:-1] + ((CANONICAL_SCHEMA_ROLE, authority.schema_relative_path),)
    source_rows: list[dict[str, Any]] = []
    source_payloads: dict[Path, bytes] = {}
    source_inodes: set[tuple[int, int]] = set()
    for role, relative in expected_roles:
        require(not relative.is_absolute() and ".." not in relative.parts, f"unsafe allowlist path: {relative}")
        source = repo_root / relative
        _path_components_are_real_dirs(repo_root, source, f"release source/{role}")
        payload, identity = _stable_regular_file(source, f"release source/{role}")
        inode = (identity["st_dev"], identity["st_ino"])
        require(inode not in source_inodes, f"release sources alias one inode: {role}")
        source_inodes.add(inode)
        if role == CANONICAL_SCHEMA_ROLE:
            require(identity["sha256"] == authority.expected_schema_sha256, "canonical schema SHA mismatch")
        source_payloads[relative] = payload
        source_rows.append({
            "role": role,
            "repo_relative_path": relative.as_posix(),
            "source": identity,
        })
    require(len(source_rows) == 11 and len(source_inodes) == 11, "release allowlist is not eleven unique files")
    source_by_role = {row["role"]: row for row in source_rows}
    input_verifier_source = source_by_role["input_identity_verifier"]["source"]
    schema_source = source_by_role[CANONICAL_SCHEMA_ROLE]["source"]
    if authority.validation_evidence_eligible:
        _, executing_identity = _stable_regular_file(Path(__file__), "executing release freezer")
        require(
            Path(executing_identity["path"]).resolve(strict=True)
            == Path(source_by_role["release_contract_freezer"]["source"]["path"]).resolve(strict=True)
            and executing_identity["sha256"] == source_by_role["release_contract_freezer"]["source"]["sha256"],
            "production freezer is not the allowlisted source",
        )
    pre_document, pre_identity = _validate_pre_receipt(
        pre_input_identity_receipt,
        Path(manifest_identity["path"]),
        manifest_document,
        authority,
        Path(schema_source["path"]),
        Path(input_verifier_source["path"]),
        str(input_verifier_source["sha256"]),
    )
    pre_payload, pre_file_identity = _stable_regular_file(pre_input_identity_receipt, "supplied PRE receipt bytes")
    pre_sidecar_path = pre_input_identity_receipt.with_name(f"{pre_input_identity_receipt.name}.sha256")
    pre_sidecar_payload, pre_sidecar_identity = _stable_regular_file(pre_sidecar_path, "supplied PRE sidecar bytes")
    require({**pre_file_identity, "sidecar": pre_sidecar_identity} == pre_identity, "supplied PRE changed after authentication")
    require(pre_input_identity_receipt.name not in {"", ".", ".."}, "unsafe PRE receipt basename")

    staging = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.staging.", dir=str(output_root.parent)))
    publication_moved = False
    try:
        snapshot_root = staging / "code_snapshot"
        snapshot_rows: list[dict[str, Any]] = []
        copy_inodes: set[tuple[int, int]] = set()
        for role, relative in expected_roles:
            destination = snapshot_root / relative
            copy_identity = _write_new_file(destination, source_payloads[relative], 0o444)
            copy_inode = (copy_identity["st_dev"], copy_identity["st_ino"])
            require(copy_inode not in copy_inodes, f"snapshot copies alias one inode: {role}")
            copy_inodes.add(copy_inode)
            require(
                copy_identity["sha256"] == source_by_role[role]["source"]["sha256"],
                f"snapshot copy SHA mismatch: {role}",
            )
            snapshot_rows.append({
                **source_by_role[role],
                "snapshot": {
                    **copy_identity,
                    "path": (Path("code_snapshot") / relative).as_posix(),
                },
            })
        require(len(copy_inodes) == 11, "snapshot copies are not eleven unique files")

        input_contract_root = staging / "input_contract"
        manifest_copy_relative = Path("input_contract/canonical_manifest.json")
        manifest_copy = _write_new_file(staging / manifest_copy_relative, manifest_payload, 0o444)
        manifest_copy["path"] = manifest_copy_relative.as_posix()
        pre_copy_relative = Path("input_contract") / pre_input_identity_receipt.name
        pre_sidecar_copy_relative = Path("input_contract") / pre_sidecar_path.name
        pre_copy = _write_new_file(staging / pre_copy_relative, pre_payload, 0o444)
        pre_copy["path"] = pre_copy_relative.as_posix()
        pre_sidecar_copy = _write_new_file(staging / pre_sidecar_copy_relative, pre_sidecar_payload, 0o444)
        pre_sidecar_copy["path"] = pre_sidecar_copy_relative.as_posix()
        require(
            len({
                *(copy_inodes),
                (manifest_copy["st_dev"], manifest_copy["st_ino"]),
                (pre_copy["st_dev"], pre_copy["st_ino"]),
                (pre_sidecar_copy["st_dev"], pre_sidecar_copy["st_ino"]),
            }) == 14,
            "immutable contract copies alias another file",
        )

        _, fresh_summary = _run_fresh_input_verifier(
            manifest_path=Path(manifest_identity["path"]),
            manifest_document=manifest_document,
            schema_path=Path(schema_source["path"]),
            verifier_path=Path(input_verifier_source["path"]),
            verifier_sha256=str(input_verifier_source["sha256"]),
            authority=authority,
            expected_snapshot=pre_document["input_identity_snapshot"],
            fresh_factory=_fresh_pre_factory,
        )

        # First publication-boundary audit: every origin and immutable copy is
        # reread after the fresh verifier and before the run manifest is made.
        manifest_payload_now, manifest_identity_now = _stable_regular_file(manifest_path, "publish-boundary canonical manifest")
        require(manifest_payload_now == manifest_payload and manifest_identity_now == manifest_identity, "canonical manifest drifted before publish")
        pre_now, pre_identity_now = _validate_pre_receipt(
            pre_input_identity_receipt, Path(manifest_identity["path"]), manifest_document, authority,
            Path(schema_source["path"]), Path(input_verifier_source["path"]),
            str(input_verifier_source["sha256"]),
        )
        require(pre_now == pre_document and pre_identity_now == pre_identity, "supplied PRE drifted before publish")
        for row in source_rows:
            relative = Path(row["repo_relative_path"])
            payload_now, identity_now = _stable_regular_file(repo_root / relative, f"publish-boundary source/{row['role']}")
            require(payload_now == source_payloads[relative] and identity_now == row["source"], f"release source drifted before publish: {row['role']}")
        for row in snapshot_rows:
            _, observed_copy = _stable_regular_file(staging / row["snapshot"]["path"], f"publish-boundary snapshot/{row['role']}")
            require(observed_copy["sha256"] == row["snapshot"]["sha256"] and observed_copy["st_nlink"] == 1, f"snapshot copy drifted before publish: {row['role']}")
        _, observed_manifest_copy = _stable_regular_file(staging / manifest_copy_relative, "publish-boundary manifest copy")
        require(observed_manifest_copy["sha256"] == manifest_identity["sha256"] and observed_manifest_copy["st_nlink"] == 1, "manifest copy drifted before publish")
        copied_pre_document, copied_pre_identity = _authenticate_json_receipt(staging / pre_copy_relative, "publish-boundary PRE copy")
        require(copied_pre_document == pre_document and copied_pre_identity["sha256"] == pre_identity["sha256"] and copied_pre_identity["st_nlink"] == 1 and copied_pre_identity["sidecar"]["st_nlink"] == 1, "PRE copy drifted before publish")

        runtime = _validate_runtime_identity(
            _runtime_identity() if _runtime is None else dict(_runtime), "freeze runtime"
        )
        producer_entry = next(row for row in snapshot_rows if row["role"] == "release_contract_freezer")
        code_entrypoints = {
            role: (Path("code_snapshot") / relative).as_posix()
            for role, relative in expected_roles if role != CANONICAL_SCHEMA_ROLE
        }
        entrypoints = {
            **code_entrypoints,
            "canonical_manifest_copy": manifest_copy_relative.as_posix(),
            "pre_input_identity_receipt_copy": pre_copy_relative.as_posix(),
            "canonical_schema_copy": (
                Path("code_snapshot") / authority.schema_relative_path
            ).as_posix(),
        }
        created_at = dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")
        run_manifest: dict[str, Any] = {
            "schema_name": RUN_SCHEMA_NAME,
            "schema_version": RUN_SCHEMA_VERSION,
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "authority_mode": authority.mode,
            "validation_evidence_eligible": authority.validation_evidence_eligible,
            "created_at_utc": created_at,
            "scope": {
                "technical_datasets": 7,
                "biological_samples": 6,
                "chromosomes": 22,
                "chromosome_names": [f"chr{index}" for index in range(1, 23)],
                "tasks": 154,
                "datasets": list(DATASETS),
            },
            "entrypoints": entrypoints,
            "canonical_manifest": {
                "expected_sha256": authority.expected_manifest_sha256,
                "origin": manifest_identity,
                "immutable_copy": manifest_copy,
            },
            "pre_input_identity_receipt": {
                "origin": pre_identity,
                "immutable_copy": {**pre_copy, "sidecar": pre_sidecar_copy},
                "input_identity_snapshot_sha256": pre_document["input_identity_snapshot_sha256"],
                "receipt_semantic_sha256": canonical_json_sha256(pre_document),
                "authority_mode": pre_document["authority_mode"],
                "validation_evidence_eligible": pre_document["validation_evidence_eligible"],
                "artifact_roles": 42,
            },
            "fresh_input_identity_verification": fresh_summary,
            "source_snapshot": {
                "repo_root": str(repo_root),
                "snapshot_root": "code_snapshot",
                "n_files": len(snapshot_rows),
                "entries": snapshot_rows,
                "entrypoints": code_entrypoints,
                "exact_allowlist_semantic_sha256": canonical_json_sha256([
                    {"role": role, "repo_relative_path": relative.as_posix()}
                    for role, relative in expected_roles
                ]),
            },
            "canonical_schema": {
                "role": CANONICAL_SCHEMA_ROLE,
                "repo_relative_path": authority.schema_relative_path.as_posix(),
                "sha256": authority.expected_schema_sha256,
            },
            "runtime": runtime,
            "runtime_semantic_sha256": canonical_json_sha256(runtime),
            "parameters": FROZEN_PARAMETERS,
            "producer": {
                "role": "release_contract_freezer",
                "repo_relative_path": producer_entry["repo_relative_path"],
                "source_sha256": producer_entry["source"]["sha256"],
                "immutable_copy_path": producer_entry["snapshot"]["path"],
                "immutable_copy_sha256": producer_entry["snapshot"]["sha256"],
            },
            "assurance": {
                "code_snapshot": "exact source bytes copied to regular non-aliased files; files 0444; directories 0555",
                "input_identity": "exact manifest-derived PRE plus mandatory fresh verifier snapshot equality",
                "immutable_inputs": "canonical manifest and supplied PRE receipt/sidecar copied read-only into contract",
                "publication": "same-filesystem renameat2(RENAME_NOREPLACE) after complete staging",
            },
            "checks": {
                "canonical_or_explicit_synthetic_test_authority_selected": True,
                "supplied_pre_receipt_authenticated_and_semantically_exact": True,
                "fresh_input_verifier_rerun_passed_and_eligible": fresh_summary["all_pass"] is True and fresh_summary["validation_evidence_eligible"] is authority.validation_evidence_eligible,
                "fresh_snapshot_exactly_equals_supplied_pre_snapshot": fresh_summary["exactly_equals_supplied_pre_snapshot"] is True,
                "manifest_sha_matches_selected_authority": True,
                "exact_eleven_file_allowlist_frozen": len(snapshot_rows) == 11,
                "all_sources_regular_non_symlink_single_link_and_non_aliased": len(source_inodes) == 11 and all(row["source"]["st_nlink"] == 1 for row in source_rows),
                "all_copies_regular_single_link_non_aliased_and_sha_equal_source": len(copy_inodes) == 11 and all(row["snapshot"]["st_nlink"] == 1 and row["snapshot"]["sha256"] == row["source"]["sha256"] for row in snapshot_rows),
                "canonical_manifest_and_pre_receipt_copied_immutably": manifest_copy["st_nlink"] == pre_copy["st_nlink"] == pre_sidecar_copy["st_nlink"] == 1,
                "publish_boundary_recheck_completed": True,
                "bootstrap20_and_seed20260716_frozen": (
                    FROZEN_PARAMETERS["ranking"]["conditional_candidate_ranking_bootstrap_replicates"] == 20
                    and FROZEN_PARAMETERS["ranking"]["conditional_candidate_ranking_bootstrap_seed"] == 20260716
                ),
            },
        }
        require(set(run_manifest["checks"]) == set(RUN_CHECK_NAMES), "internal run check schema mismatch")
        run_manifest["all_pass"] = all(value is True for value in run_manifest["checks"].values())
        manifest_name = "m2_run_manifest.json"
        sidecar_name = f"{manifest_name}.sha256"
        run_manifest["receipt_integrity"] = {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": sidecar_name,
            "covers": manifest_name,
        }
        run_payload = json.dumps(
            run_manifest, ensure_ascii=False, allow_nan=False, indent=2
        ).encode("utf-8") + b"\n"
        run_digest = sha256_bytes(run_payload)
        _write_new_file(staging / manifest_name, run_payload, 0o444)
        _write_new_file(
            staging / sidecar_name, f"{run_digest}  {manifest_name}\n".encode("ascii"), 0o444
        )
        _set_snapshot_directory_modes(staging)

        # Final no-write boundary immediately before the atomic no-replace
        # rename.  This also authenticates the staged run manifest itself.
        _, final_manifest_origin = _stable_regular_file(manifest_path, "final-boundary canonical manifest origin")
        require(final_manifest_origin == manifest_identity, "canonical manifest changed at final publish boundary")
        final_pre, final_pre_identity = _validate_pre_receipt(
            pre_input_identity_receipt, Path(manifest_identity["path"]), manifest_document, authority,
            Path(schema_source["path"]), Path(input_verifier_source["path"]), str(input_verifier_source["sha256"]),
        )
        require(final_pre == pre_document and final_pre_identity == pre_identity, "supplied PRE changed at final publish boundary")
        for row in source_rows:
            _, observed_source = _stable_regular_file(repo_root / row["repo_relative_path"], f"final-boundary source/{row['role']}")
            require(observed_source == row["source"], f"release source changed at final publish boundary: {row['role']}")
        for row in snapshot_rows:
            _, observed_copy = _stable_regular_file(staging / row["snapshot"]["path"], f"final-boundary snapshot/{row['role']}")
            require(observed_copy["sha256"] == row["snapshot"]["sha256"] and observed_copy["st_nlink"] == 1 and observed_copy["mode_octal"] == "0444", f"snapshot copy failed final boundary: {row['role']}")
        staged_run, _ = _authenticate_json_receipt(staging / manifest_name, "staged M2 run manifest")
        require(staged_run == run_manifest, "staged run manifest semantic mismatch")
        staged_pre, _ = _authenticate_json_receipt(staging / pre_copy_relative, "staged immutable PRE")
        require(staged_pre == pre_document, "staged immutable PRE semantic mismatch")
        _rename_noreplace(staging, output_root)
        publication_moved = True
        _fsync_directory(output_root.parent)
        return {
            "output_contract_root": str(output_root),
            "run_manifest": str(output_root / manifest_name),
            "run_manifest_sha256": run_digest,
            "authority_mode": authority.mode,
            "validation_evidence_eligible": authority.validation_evidence_eligible,
            "n_snapshot_files": len(snapshot_rows),
        }
    except BaseException as exc:
        failed_tree = output_root if publication_moved else staging
        archive = _archive_failed_tree(failed_tree, "staging")
        if archive is not None and hasattr(exc, "add_note"):
            exc.add_note(f"failed release staging archived at {archive}")
        raise


def _expected_authority_for_manifest(
    run_manifest: Mapping[str, Any], test_authority: _Authority | None
) -> _Authority:
    mode = run_manifest.get("authority_mode")
    if mode == PRODUCTION_AUTHORITY.mode:
        return PRODUCTION_AUTHORITY
    require(test_authority is not None, "synthetic release contract is forbidden in production verify-only mode")
    require(mode == test_authority.mode, "synthetic authority mode mismatch")
    return test_authority


FILE_IDENTITY_KEYS = {
    "path", "st_dev", "st_ino", "st_nlink", "size_bytes", "mtime_ns", "ctime_ns",
    "mode_octal", "sha256",
}


def _validate_file_identity_record(value: Any, label: str) -> Mapping[str, Any]:
    row = _exact_keys(value, FILE_IDENTITY_KEYS, label)
    require(isinstance(row["path"], str) and row["path"], f"{label}.path invalid")
    for key in ("st_dev", "st_ino", "st_nlink", "size_bytes", "mtime_ns", "ctime_ns"):
        _strict_int(row[key], f"{label}.{key}")
    require(row["st_nlink"] == 1, f"{label}.st_nlink is not one")
    require(isinstance(row["mode_octal"], str) and re.fullmatch(r"[0-7]{4}", row["mode_octal"]) is not None, f"{label}.mode_octal invalid")
    require(isinstance(row["sha256"], str) and SHA256_RE.fullmatch(row["sha256"]) is not None, f"{label}.sha256 invalid")
    return row


def _validate_fresh_summary(value: Any, authority: _Authority, expected_snapshot_sha: str) -> Mapping[str, Any]:
    keys = {
        "execution_mode", "verifier_path", "verifier_sha256", "receipt_sha256",
        "receipt_semantic_sha256", "checks_semantic_sha256", "artifact_audit_semantic_sha256",
        "input_identity_snapshot_sha256", "all_pass", "validation_evidence_eligible",
        "exactly_equals_supplied_pre_snapshot",
    }
    row = _exact_keys(value, keys, "fresh_input_identity_verification")
    expected_mode = (
        "production_subprocess_required" if authority.validation_evidence_eligible
        else "private_synthetic_test_injection_not_validation_evidence"
    )
    require(row["execution_mode"] == expected_mode, "fresh verifier execution mode mismatch")
    require(isinstance(row["verifier_path"], str) and row["verifier_path"], "fresh verifier path invalid")
    for key in ("verifier_sha256", "receipt_sha256", "receipt_semantic_sha256", "checks_semantic_sha256", "artifact_audit_semantic_sha256", "input_identity_snapshot_sha256"):
        require(isinstance(row[key], str) and SHA256_RE.fullmatch(row[key]) is not None, f"fresh verifier {key} invalid")
    require(row["input_identity_snapshot_sha256"] == expected_snapshot_sha, "fresh verifier snapshot digest mismatch")
    require(row["all_pass"] is True and row["validation_evidence_eligible"] is authority.validation_evidence_eligible and row["exactly_equals_supplied_pre_snapshot"] is True, "fresh verifier result mismatch")
    return row


def verify_release_contract(
    contract_root: Path,
    *,
    _test_authority: _Authority | None = None,
    _runtime: Mapping[str, Any] | None = None,
    _fresh_pre_factory: Any | None = None,
) -> dict[str, Any]:
    root = contract_root.absolute()
    root_lstat = os.lstat(root)
    require(stat.S_ISDIR(root_lstat.st_mode) and not stat.S_ISLNK(root_lstat.st_mode), "contract root is not a regular directory")
    require(stat.S_IMODE(root_lstat.st_mode) == 0o555, "contract root mode is not 0555")
    require(
        {path.name for path in root.iterdir()} == {
            "code_snapshot", "input_contract", "m2_run_manifest.json", "m2_run_manifest.json.sha256",
        },
        "contract root contains an unexpected top-level entry",
    )
    manifest_path = root / "m2_run_manifest.json"
    run_manifest, run_identity = _authenticate_json_receipt(manifest_path, "M2 run manifest")
    require(run_identity.get("mode_octal") == "0444", "run manifest mode is not 0444")
    require(
        (run_identity.get("sidecar") or {}).get("mode_octal") == "0444",
        "run manifest sidecar mode is not 0444",
    )
    run_keys = {
        "schema_name", "schema_version", "task_type", "authority_mode",
        "validation_evidence_eligible", "created_at_utc", "scope", "entrypoints",
        "canonical_manifest", "pre_input_identity_receipt", "fresh_input_identity_verification",
        "source_snapshot", "canonical_schema", "runtime", "runtime_semantic_sha256",
        "parameters", "producer", "assurance", "checks", "all_pass", "receipt_integrity",
    }
    _exact_keys(run_manifest, run_keys, "run manifest")
    require(run_manifest["schema_name"] == RUN_SCHEMA_NAME and run_manifest["schema_version"] == RUN_SCHEMA_VERSION, "run manifest schema mismatch")
    require(run_manifest["task_type"] == "B_COMPREHENSIVE_VALIDATION", "run manifest task_type mismatch")
    require(run_manifest["all_pass"] is True, "run manifest did not pass")
    try:
        created = dt.datetime.fromisoformat(str(run_manifest["created_at_utc"]).replace("Z", "+00:00"))
    except ValueError as exc:
        raise ReleaseContractError("run manifest created_at_utc is invalid") from exc
    require(created.tzinfo is not None, "run manifest created_at_utc lacks timezone")
    authority = _expected_authority_for_manifest(run_manifest, _test_authority)
    if authority.validation_evidence_eligible:
        require(_test_authority is None and _runtime is None and _fresh_pre_factory is None, "production verify-only forbids private overrides")
    else:
        require(_test_authority is not None and _runtime is not None and _fresh_pre_factory is not None, "synthetic verify requires explicit private fixtures")
    require(
        run_manifest.get("validation_evidence_eligible") is authority.validation_evidence_eligible,
        "run manifest validation eligibility mismatch",
    )
    checks = _exact_keys(run_manifest["checks"], set(RUN_CHECK_NAMES), "run manifest checks")
    require(all(checks[name] is True for name in RUN_CHECK_NAMES), "run manifest named checks are failing")
    expected_scope = {
        "technical_datasets": 7, "biological_samples": 6, "chromosomes": 22,
        "chromosome_names": [f"chr{index}" for index in range(1, 23)],
        "tasks": 154, "datasets": list(DATASETS),
    }
    require(run_manifest["scope"] == expected_scope, "run manifest exact scope mismatch")
    require(run_manifest.get("parameters") == FROZEN_PARAMETERS, "run manifest frozen parameters mismatch")
    require(run_manifest["assurance"] == {
        "code_snapshot": "exact source bytes copied to regular non-aliased files; files 0444; directories 0555",
        "input_identity": "exact manifest-derived PRE plus mandatory fresh verifier snapshot equality",
        "immutable_inputs": "canonical manifest and supplied PRE receipt/sidecar copied read-only into contract",
        "publication": "same-filesystem renameat2(RENAME_NOREPLACE) after complete staging",
    }, "run manifest assurance mismatch")
    runtime = _validate_runtime_identity(run_manifest["runtime"], "run manifest runtime")
    require(run_manifest["runtime_semantic_sha256"] == canonical_json_sha256(runtime), "run manifest runtime semantic digest mismatch")
    current_runtime = _validate_runtime_identity(
        _runtime_identity() if _runtime is None else dict(_runtime), "current verify runtime"
    )
    require(current_runtime == runtime, "current runtime differs from frozen runtime")

    expected_roles = ALL_FROZEN_ROLES[:-1] + ((CANONICAL_SCHEMA_ROLE, authority.schema_relative_path),)
    expected_by_role = {role: relative for role, relative in expected_roles}
    snapshot = _exact_keys(run_manifest["source_snapshot"], {"repo_root", "snapshot_root", "n_files", "entries", "entrypoints", "exact_allowlist_semantic_sha256"}, "source_snapshot")
    require(Path(str(snapshot["repo_root"])).is_absolute(), "snapshot repo_root is not absolute")
    require(snapshot["snapshot_root"] == "code_snapshot", "snapshot_root mismatch")
    entries = snapshot.get("entries")
    require(isinstance(entries, list) and len(entries) == 11, "run manifest snapshot entry count mismatch")
    require(snapshot.get("n_files") == 11, "run manifest snapshot n_files mismatch")
    entry_by_role: dict[str, Mapping[str, Any]] = {}
    for entry in entries:
        _exact_keys(entry, {"role", "repo_relative_path", "source", "snapshot"}, "snapshot entry")
        role = str(entry.get("role"))
        require(role not in entry_by_role, f"duplicate snapshot role: {role}")
        entry_by_role[role] = entry
    require(set(entry_by_role) == set(expected_by_role), "snapshot role allowlist mismatch")
    expected_entrypoints = {
        role: (Path("code_snapshot") / relative).as_posix()
        for role, relative in expected_roles
        if role != CANONICAL_SCHEMA_ROLE
    }
    require(snapshot.get("entrypoints") == expected_entrypoints, "snapshot entrypoint map mismatch")
    manifest_contract = _exact_keys(run_manifest["canonical_manifest"], {"expected_sha256", "origin", "immutable_copy"}, "canonical_manifest")
    require(manifest_contract["expected_sha256"] == authority.expected_manifest_sha256, "canonical manifest authority SHA mismatch")
    manifest_origin = _validate_file_identity_record(manifest_contract["origin"], "canonical_manifest.origin")
    manifest_copy_record = _validate_file_identity_record(manifest_contract["immutable_copy"], "canonical_manifest.immutable_copy")
    pre_contract = _exact_keys(run_manifest["pre_input_identity_receipt"], {"origin", "immutable_copy", "input_identity_snapshot_sha256", "receipt_semantic_sha256", "authority_mode", "validation_evidence_eligible", "artifact_roles"}, "pre_input_identity_receipt")
    pre_origin = _exact_keys(pre_contract["origin"], FILE_IDENTITY_KEYS | {"sidecar"}, "pre_input_identity_receipt.origin")
    _validate_file_identity_record({key: pre_origin[key] for key in FILE_IDENTITY_KEYS}, "pre_input_identity_receipt.origin receipt")
    _validate_file_identity_record(pre_origin["sidecar"], "pre_input_identity_receipt.origin sidecar")
    pre_copy_record = _exact_keys(pre_contract["immutable_copy"], FILE_IDENTITY_KEYS | {"sidecar"}, "pre_input_identity_receipt.immutable_copy")
    _validate_file_identity_record({key: pre_copy_record[key] for key in FILE_IDENTITY_KEYS}, "pre_input_identity_receipt.copy receipt")
    _validate_file_identity_record(pre_copy_record["sidecar"], "pre_input_identity_receipt.copy sidecar")
    require(pre_contract["authority_mode"] == authority.expected_pre_authority_mode and pre_contract["validation_evidence_eligible"] is authority.validation_evidence_eligible and pre_contract["artifact_roles"] == 42, "PRE contract authority/scope mismatch")
    frozen_fresh = _validate_fresh_summary(run_manifest["fresh_input_identity_verification"], authority, str(pre_contract["input_identity_snapshot_sha256"]))

    snapshot_root = root / "code_snapshot"
    expected_files = {snapshot_root / relative for relative in expected_by_role.values()}
    actual_files: set[Path] = set()
    actual_directories: set[Path] = {snapshot_root}
    for path in snapshot_root.rglob("*"):
        observed = os.lstat(path)
        require(not stat.S_ISLNK(observed.st_mode), f"snapshot contains symlink: {path}")
        if stat.S_ISDIR(observed.st_mode):
            actual_directories.add(path)
        elif stat.S_ISREG(observed.st_mode):
            actual_files.add(path)
        else:
            raise ReleaseContractError(f"snapshot contains non-file/non-directory: {path}")
    require(actual_files == expected_files, "snapshot file set differs from exact allowlist")
    expected_directories = {snapshot_root}
    for path in expected_files:
        parent = path.parent
        while parent != root:
            expected_directories.add(parent)
            parent = parent.parent
    require(actual_directories == expected_directories, "snapshot directory set differs from exact layout")
    for directory in actual_directories:
        require(stat.S_IMODE(os.lstat(directory).st_mode) == 0o555, f"snapshot directory mode is not 0555: {directory}")

    copy_inodes: set[tuple[int, int]] = set()
    verified_entries = []
    for role, relative in expected_roles:
        entry = entry_by_role[role]
        require(entry.get("repo_relative_path") == relative.as_posix(), f"snapshot relative path mismatch: {role}")
        source = _validate_file_identity_record(entry["source"], f"snapshot source/{role}")
        require(Path(str(source["path"])).absolute() == Path(str(snapshot["repo_root"])) / relative, f"snapshot source path mismatch: {role}")
        copy = _validate_file_identity_record(entry["snapshot"], f"snapshot copy/{role}")
        expected_copy_relative = (Path("code_snapshot") / relative).as_posix()
        require(copy.get("path") == expected_copy_relative, f"snapshot copy path mismatch: {role}")
        payload, observed = _stable_regular_file(root / expected_copy_relative, f"snapshot/{role}")
        del payload
        require(stat.S_IMODE(os.lstat(root / expected_copy_relative).st_mode) == 0o444, f"snapshot file mode is not 0444: {role}")
        require(observed["sha256"] == copy.get("sha256"), f"snapshot SHA differs from manifest: {role}")
        require(observed["size_bytes"] == copy.get("size_bytes"), f"snapshot size differs from manifest: {role}")
        require(observed["st_dev"] == copy["st_dev"] and observed["st_ino"] == copy["st_ino"] and observed["st_nlink"] == copy["st_nlink"] == 1, f"snapshot inode/link identity differs: {role}")
        require(observed["mtime_ns"] == copy["mtime_ns"] and observed["ctime_ns"] == copy["ctime_ns"] and observed["mode_octal"] == copy["mode_octal"] == "0444", f"snapshot metadata differs: {role}")
        require(copy.get("sha256") == source.get("sha256"), f"source/copy SHA mismatch in manifest: {role}")
        inode = (observed["st_dev"], observed["st_ino"])
        require(inode not in copy_inodes, f"snapshot files alias one inode: {role}")
        copy_inodes.add(inode)
        if role == CANONICAL_SCHEMA_ROLE:
            require(observed["sha256"] == authority.expected_schema_sha256, "snapshot canonical schema SHA mismatch")
        verified_entries.append({"role": role, "sha256": observed["sha256"], "mode_octal": "0444"})

    expected_allowlist_digest = canonical_json_sha256([
        {"role": role, "repo_relative_path": relative.as_posix()}
        for role, relative in expected_roles
    ])
    require(
        snapshot.get("exact_allowlist_semantic_sha256") == expected_allowlist_digest,
        "snapshot allowlist semantic digest mismatch",
    )
    expected_top_entrypoints = {
        **expected_entrypoints,
        "canonical_manifest_copy": str(manifest_copy_record["path"]),
        "pre_input_identity_receipt_copy": str(pre_copy_record["path"]),
        "canonical_schema_copy": (Path("code_snapshot") / authority.schema_relative_path).as_posix(),
    }
    require(run_manifest["entrypoints"] == expected_top_entrypoints, "run manifest immutable entrypoints mismatch")
    require(run_manifest["entrypoints"]["canonical_manifest_copy"] == "input_contract/canonical_manifest.json", "canonical_manifest_copy entrypoint is not fixed")
    input_root = root / "input_contract"
    require(stat.S_IMODE(os.lstat(input_root).st_mode) == 0o555 and not stat.S_ISLNK(os.lstat(input_root).st_mode), "input_contract directory is not immutable")
    expected_input_files = {
        root / str(manifest_copy_record["path"]), root / str(pre_copy_record["path"]),
        root / str(pre_copy_record["sidecar"]["path"]),
    }
    require({path for path in input_root.iterdir()} == expected_input_files, "input_contract file set mismatch")
    canonical_path = root / str(manifest_copy_record["path"])
    canonical_payload, observed_manifest = _stable_regular_file(canonical_path, "immutable canonical manifest copy")
    for key in FILE_IDENTITY_KEYS - {"path"}:
        require(observed_manifest[key] == manifest_copy_record[key], f"immutable canonical manifest {key} mismatch")
    require(observed_manifest["sha256"] == authority.expected_manifest_sha256 == manifest_origin["sha256"], "immutable canonical manifest SHA mismatch")
    manifest_document = strict_json_load_bytes(canonical_payload, canonical_path)
    require(isinstance(manifest_document, dict), "immutable canonical manifest is not an object")

    input_verifier_entry = entry_by_role["input_identity_verifier"]
    input_verifier_copy = input_verifier_entry["snapshot"]
    schema_entry = entry_by_role[CANONICAL_SCHEMA_ROLE]
    schema_copy = schema_entry["snapshot"]
    pre_path = root / str(pre_copy_record["path"])
    pre_receipt, observed_pre = _validate_pre_receipt(
        pre_path,
        Path(str(manifest_origin["path"])),
        manifest_document,
        authority,
        Path(str(schema_entry["source"]["path"])),
        Path(str(input_verifier_entry["source"]["path"])),
        str(input_verifier_copy.get("sha256", "")),
    )
    for key in FILE_IDENTITY_KEYS - {"path"}:
        require(observed_pre[key] == pre_copy_record[key], f"immutable PRE receipt {key} mismatch")
    for key in FILE_IDENTITY_KEYS - {"path"}:
        require(observed_pre["sidecar"][key] == pre_copy_record["sidecar"][key], f"immutable PRE sidecar {key} mismatch")
    require(observed_pre["sha256"] == pre_origin["sha256"], "immutable PRE differs from frozen origin")
    require(canonical_json_sha256(pre_receipt) == pre_contract["receipt_semantic_sha256"], "immutable PRE semantic digest mismatch")
    require(
        pre_receipt["input_identity_snapshot_sha256"]
        == pre_contract.get("input_identity_snapshot_sha256"),
        "PRE snapshot digest differs from run manifest",
    )
    require(
        frozen_fresh["verifier_path"] == str(input_verifier_entry["source"]["path"])
        and frozen_fresh["verifier_sha256"] == input_verifier_entry["source"]["sha256"],
        "freeze-time fresh verifier is not bound to frozen source",
    )

    producer = _exact_keys(run_manifest["producer"], {"role", "repo_relative_path", "source_sha256", "immutable_copy_path", "immutable_copy_sha256"}, "run manifest producer")
    freezer_entry = entry_by_role["release_contract_freezer"]
    require(producer == {
        "role": "release_contract_freezer",
        "repo_relative_path": freezer_entry["repo_relative_path"],
        "source_sha256": freezer_entry["source"]["sha256"],
        "immutable_copy_path": freezer_entry["snapshot"]["path"],
        "immutable_copy_sha256": freezer_entry["snapshot"]["sha256"],
    }, "run manifest producer is not the immutable freezer copy")

    expected_current_verifier = root / str(freezer_entry["snapshot"]["path"])
    if authority.validation_evidence_eligible:
        require(Path(__file__).resolve(strict=True) == expected_current_verifier.resolve(strict=True), "production verify-only must execute the immutable freezer snapshot")
    _, fresh_verify_summary = _run_fresh_input_verifier(
        manifest_path=canonical_path,
        manifest_document=manifest_document,
        schema_path=root / str(schema_copy["path"]),
        verifier_path=root / str(input_verifier_copy["path"]),
        verifier_sha256=str(input_verifier_copy["sha256"]),
        authority=authority,
        expected_snapshot=pre_receipt["input_identity_snapshot"],
        fresh_factory=_fresh_pre_factory,
    )
    verifier_payload, verifier_identity = _stable_regular_file(
        Path(__file__), "release-contract verify-only source"
    )
    del verifier_payload
    verification = {
        "schema_name": VERIFY_SCHEMA_NAME,
        "schema_version": VERIFY_SCHEMA_VERSION,
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "authority_mode": authority.mode,
        "validation_evidence_eligible": authority.validation_evidence_eligible,
        "verified_at_utc": dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z"),
        "scope": expected_scope,
        "run_manifest": run_identity,
        "canonical_manifest": observed_manifest,
        "pre_input_identity_receipt": observed_pre,
        "fresh_input_identity_verification": fresh_verify_summary,
        "runtime": current_runtime,
        "runtime_semantic_sha256": canonical_json_sha256(current_runtime),
        "snapshot": {"n_files": len(verified_entries), "entries": verified_entries},
        "verifier": verifier_identity,
        "checks": {
            "run_manifest_authenticated": True,
            "run_manifest_scope_and_parameters_exact": True,
            "canonical_manifest_sha_reverified": True,
            "pre_receipt_and_sidecar_reauthenticated": True,
            "immutable_pre_receipt_authority_and_42_roles_revalidated": True,
            "fresh_snapshot_verifier_rerun_and_exact_pre_equality": fresh_verify_summary["exactly_equals_supplied_pre_snapshot"] is True,
            "runtime_schema_digest_and_current_identity_exact": current_runtime == runtime,
            "producer_bound_to_immutable_freezer_copy": True,
            "exact_snapshot_file_and_directory_sets_verified": True,
            "all_snapshot_files_sha_and_mode_verified": True,
            "all_snapshot_directories_mode_verified": True,
            "snapshot_files_non_symlink_single_link_and_non_aliased": len(copy_inodes) == 11,
        },
    }
    verification["all_pass"] = all(value is True for value in verification["checks"].values())
    return verification


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    freeze_parser = subparsers.add_parser("freeze", help="create a new canonical release contract")
    freeze_parser.add_argument("--manifest", required=True, type=Path)
    freeze_parser.add_argument("--pre-input-identity-receipt", required=True, type=Path)
    freeze_parser.add_argument("--output-contract-root", required=True, type=Path)
    verify_parser = subparsers.add_parser("verify-only", help="verify a frozen contract without running M2")
    verify_parser.add_argument("--contract-root", required=True, type=Path)
    verify_parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.command == "freeze":
        try:
            summary = freeze_release_contract(
                args.manifest, args.pre_input_identity_receipt, args.output_contract_root
            )
        except (ReleaseContractError, OSError, ValueError, TypeError, KeyError) as exc:
            print(json.dumps({"all_pass": False, "failure": f"{type(exc).__name__}: {exc}"}), file=sys.stderr)
            return 1
        print(json.dumps({"all_pass": True, **summary}, ensure_ascii=False))
        return 0 if summary["validation_evidence_eligible"] is True else 1

    failure: str | None = None
    try:
        verification = verify_release_contract(args.contract_root)
    except (ReleaseContractError, OSError, ValueError, TypeError, KeyError) as exc:
        failure = f"{type(exc).__name__}: {exc}"
        verification = {
            "schema_name": VERIFY_SCHEMA_NAME,
            "schema_version": VERIFY_SCHEMA_VERSION,
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "authority_mode": "UNKNOWN_OR_REJECTED",
            "validation_evidence_eligible": False,
            "contract_root": str(args.contract_root.absolute()),
            "failure": failure,
            "checks": {"verification_completed_without_contract_violation": False},
            "all_pass": False,
        }
    digest = _write_authenticated_exclusive(args.output, verification)
    print(json.dumps({
        "output": str(args.output.absolute()),
        "sha256": digest,
        "all_pass": verification["all_pass"],
        "validation_evidence_eligible": verification.get("validation_evidence_eligible", False),
        "failure": failure,
    }, ensure_ascii=False))
    return 0 if verification["all_pass"] and verification.get("validation_evidence_eligible") is True else 1


if __name__ == "__main__":
    raise SystemExit(main())
