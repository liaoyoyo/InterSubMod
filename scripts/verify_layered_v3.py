#!/usr/bin/env python3
"""Independent fail-closed verifier for the real layered-v3 contracts.

Authorities are limited to artifacts created by ``layered_v3_lifecycle.py``:
``frozen_input_lock.json``, ``launch_receipt.json``, the source bundle,
``environment_lock.json``, the hash-chained state events/current state, and the
named outputs below the published run root.  The frozen lock must conform to
``intersubmod.layered_input_lock`` 1.0.0 as emitted by
``validate_layered_v3_inputs.py``.

This program never creates ``_SUCCESS``.  On exit 0 it publishes only atomic
JSON/TSV/SHA-256 verification artifacts; ``RunLifecycle.succeed()`` remains the
sole authority that transitions VERIFYING -> SUCCEEDED and creates _SUCCESS.

Exit codes: 0 PASS, 2 CLI/bootstrap error, 7 controlled verification failure,
70 unexpected internal failure.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import stat
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


VERIFIER_VERSION = "2.0.0"
LOCK_SCHEMA_NAME = "intersubmod.layered_input_lock"
LOCK_SCHEMA_VERSION = "1.0.0"
LOCK_SCHEMA_REFERENCE = (
    "InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/"
    "schemas/layered_input_lock_v1.schema.json"
)
RECEIPT_SCHEMA_NAME = "intersubmod.layered_run_receipt"
RECEIPT_SCHEMA_VERSION = "1.0.0"
STATE_SCHEMA_NAME = "intersubmod.layered_run_state"
STATE_SCHEMA_VERSION = "1.0.0"
SOURCE_BUNDLE_SCHEMA_NAME = "intersubmod.layered_source_bundle"
SOURCE_BUNDLE_SCHEMA_VERSION = "1.0.0"
ENVIRONMENT_SCHEMA_NAME = "intersubmod.layered_environment_lock"
ENVIRONMENT_SCHEMA_VERSION = "1.0.0"
OUTPUT_MANIFEST_SCHEMA_NAME = "intersubmod.layered_sample_output_manifest"
OUTPUT_MANIFEST_SCHEMA_VERSION = "1.0.0"
SUMMARY_SCHEMA_NAME = "intersubmod.layered_verification_summary"
SUMMARY_SCHEMA_VERSION = "1.0.0"

RUN_LOCK_FILENAME = "frozen_input_lock.json"
RUN_RECEIPT_FILENAME = "launch_receipt.json"
RUN_STATE_FILENAME = "run_state.json"
SOURCE_BUNDLE_MANIFEST = "source_bundle/source_bundle_manifest.json"
ENVIRONMENT_LOCK_FILENAME = "environment_lock.json"

EXPECTED_BINDING = {
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}
EXPECTED_DATASETS = frozenset(EXPECTED_BINDING)
EXPECTED_BIOLOGICAL_IDS = frozenset(EXPECTED_BINDING.values())
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
CANONICAL_TREE_INPUT_CONTRACT = "longphase_s_recalibrated_FILTER_PASS"
SENSITIVITY_TREE_INPUT_CONTRACT = "clairs_FILTER_PASS_sensitivity"
TREE_PROFILES = {
    CANONICAL_TREE_INPUT_CONTRACT: {
        "task_type": "comprehensive_validation",
        "backbone_role": "longphase_s_recalibrated_filter_pass",
        "ledger_tree_contract": "longphase_recalibrated_PASS",
    },
    SENSITIVITY_TREE_INPUT_CONTRACT: {
        "task_type": "backbone_sensitivity",
        "backbone_role": "clairs_filter_pass_sensitivity",
        "ledger_tree_contract": "clairs_PASS_input",
    },
}
PART_CHROMOSOMES = {
    1: ("chr1", "chr6", "chr11", "chr16", "chr21"),
    2: ("chr2", "chr7", "chr12", "chr17", "chr22"),
    3: ("chr3", "chr8", "chr13", "chr18"),
    4: ("chr4", "chr9", "chr14", "chr19"),
    5: ("chr5", "chr10", "chr15", "chr20"),
}
PART_ROLES = tuple(f"mlhp_part_{number}" for number in range(1, 6))
SCIENTIFIC_ROLES = PART_ROLES + (
    "layered_reconstruction",
    "layered_region_view",
    "site_ledger",
    "site_ledger_summary",
)
FUNNEL_FIELDS = (
    "n_sSNV_scope_input",
    "n_positional_singleton",
    "n_multilocus_pre_cap_groups",
    "n_multilocus_pre_cap_sSNV",
    "n_groups_retained",
    "n_groups_read_unsupported",
    "n_sSNV_retained",
    "n_sSNV_read_unsupported",
    "n_groups_capped_by_MAX_SNV",
    "n_sSNV_cap_excluded",
    "n_sSNV_accounted",
)
PROVENANCE_FIELDS = frozenset(
    {
        "frozen_lock_sha256",
        "launch_receipt_sha256",
        "environment_lock_sha256",
        "source_bundle_manifest_sha256",
        "source_bundle_content_sha256",
        "input_set_sha256",
    }
)
COORDINATE_IDENTITY_COLUMNS = ("QNAME", "CHROM", "START0", "END0", "FLAG", "CIGAR_B2")
SIDECAR_HEADER = ("#CHROM", "START0", "END0", "QNAME", "FLAG", "MAPQ", "CIGAR_B2", "HP", "PS")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
ZERO_DIGEST = "0" * 64


class ContractError(Exception):
    def __init__(self, code: str, message: str, *, exit_code: int = 7):
        super().__init__(message)
        self.code = code
        self.exit_code = exit_code


@dataclass
class Check:
    name: str
    passed: bool
    code: str
    observed: Any
    expected: Any
    sample: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "pass": self.passed,
            "code": None if self.passed else self.code,
            "sample": self.sample,
            "observed": self.observed,
            "expected": self.expected,
        }


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number is forbidden: {value}")


def _no_duplicate_pairs(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_json_strict(path: Path, code: str = "E_SCHEMA_INVALID") -> Any:
    try:
        return json.loads(
            path.read_text(encoding="utf-8", errors="strict"),
            object_pairs_hook=_no_duplicate_pairs,
            parse_constant=_reject_constant,
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as exc:
        raise ContractError(code, f"cannot parse strict JSON {path}: {exc}") from exc


def canonical_bytes(value: Any) -> bytes:
    try:
        return json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ContractError("E_SCHEMA_INVALID", f"not canonical JSON: {exc}") from exc


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise ContractError("E_ARTIFACT_MISSING", f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def sha256_slice(path: Path, offset: int, length: int) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(offset)
            payload = handle.read(length)
    except OSError as exc:
        raise ContractError("E_POST_INPUT_IDENTITY", f"cannot read {path} at {offset}: {exc}") from exc
    if len(payload) != length:
        raise ContractError(
            "E_POST_INPUT_IDENTITY",
            f"short read for {path}: offset={offset} expected={length} observed={len(payload)}",
        )
    return hashlib.sha256(payload).hexdigest()


def require_mapping(value: Any, where: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ContractError("E_SCHEMA_INVALID", f"{where} must be an object")
    return value


def require_list(value: Any, where: str) -> list[Any]:
    if not isinstance(value, list):
        raise ContractError("E_SCHEMA_INVALID", f"{where} must be an array")
    return value


def require_exact_keys(value: Mapping[str, Any], expected: Iterable[str], where: str) -> None:
    expected_set = set(expected)
    observed = set(value)
    if observed != expected_set:
        raise ContractError(
            "E_SCHEMA_INVALID",
            f"{where} fields differ: missing={sorted(expected_set-observed)} extra={sorted(observed-expected_set)}",
        )


def require_sha256(value: Any, where: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise ContractError("E_SCHEMA_INVALID", f"{where} must be lowercase SHA-256")
    return value


def require_nonnegative_int(value: Any, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ContractError("E_SCHEMA_INVALID", f"{where} must be a non-negative integer")
    return value


def regular_file(path: Path, where: str, *, reject_symlink: bool = True) -> Path:
    try:
        if reject_symlink and path.is_symlink():
            raise ContractError("E_PATH_ESCAPE", f"{where} may not be a symlink: {path}")
        resolved = path.resolve(strict=True)
        if not stat.S_ISREG(resolved.stat().st_mode):
            raise ContractError("E_ARTIFACT_MISSING", f"{where} is not a regular file: {path}")
        return resolved
    except OSError as exc:
        raise ContractError("E_ARTIFACT_MISSING", f"{where} is unavailable: {path}: {exc}") from exc


def contained_file(path_text: Any, root: Path, where: str) -> Path:
    if not isinstance(path_text, str) or not path_text or "\x00" in path_text:
        raise ContractError("E_SCHEMA_INVALID", f"{where} must be a relative path")
    relative = Path(path_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise ContractError("E_PATH_ESCAPE", f"{where} escapes root: {path_text}")
    resolved = regular_file(root / relative, where)
    try:
        resolved.relative_to(root.resolve(strict=True))
    except ValueError as exc:
        raise ContractError("E_PATH_ESCAPE", f"{where} escapes root: {path_text}") from exc
    return resolved


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        raise ContractError("E_OUTPUT_EXISTS", f"refusing to replace verifier artifact: {path}")
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
        directory = os.open(path.parent, flags)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except BaseException:
        # A .partial is retained for audit and is never an authoritative output.
        raise


def input_set_digest(sample_meta: Mapping[str, Any]) -> str:
    return canonical_sha256(
        {
            "sample": sample_meta["sample"],
            "biological_id": sample_meta["biological_id"],
            "alignment_payload": sample_meta["alignment_payload"],
            "somatic": sample_meta["somatic"],
            "read_tags": sample_meta["read_tags"],
            "copy_number": sample_meta["copy_number"],
        }
    )


def output_name(sample: str, role: str) -> str:
    if role.startswith("mlhp_part_"):
        return f"{role}.json"
    return {
        "layered_reconstruction": f"layered_reconstruction_{sample}.json",
        "layered_region_view": f"layered_region_view_{sample}.json",
        "site_ledger": f"ssnv_site_ledger_{sample}.tsv.gz",
        "site_ledger_summary": f"ssnv_site_ledger_{sample}.summary.json",
    }[role]


class Verifier:
    def __init__(self, run_root: Path, frozen_lock: Path, launch_receipt: Path):
        self.run_root = run_root.resolve(strict=True)
        self.lock_path = frozen_lock.resolve(strict=True)
        self.receipt_path = launch_receipt.resolve(strict=True)
        self.checks: list[Check] = []
        self.samples: list[dict[str, Any]] = []
        self.hash_cache: dict[Path, str] = {}
        self.lock_file_sha256 = ""
        self.receipt_sha256 = ""
        self.environment_sha256 = ""
        self.bundle_manifest_sha256 = ""
        self.bundle_content_sha256 = ""
        self.run_id = ""
        self.tree_contract = ""
        self.tree_profile: Mapping[str, str] = {}

    def add(
        self,
        name: str,
        passed: bool,
        code: str,
        observed: Any,
        expected: Any,
        sample: str | None = None,
    ) -> bool:
        self.checks.append(Check(name, bool(passed), code, observed, expected, sample))
        return bool(passed)

    def hash(self, path: Path) -> str:
        resolved = path.resolve(strict=True)
        if resolved not in self.hash_cache:
            self.hash_cache[resolved] = sha256_file(resolved)
        return self.hash_cache[resolved]

    def verify(self) -> dict[str, Any]:
        if not self.run_root.is_dir() or self.run_root.is_symlink():
            raise ContractError("E_PATH_ESCAPE", f"run root is not a real directory: {self.run_root}")
        for path, expected_name, label in (
            (self.lock_path, RUN_LOCK_FILENAME, "frozen lock"),
            (self.receipt_path, RUN_RECEIPT_FILENAME, "launch receipt"),
        ):
            regular_file(path, label)
            try:
                path.relative_to(self.run_root)
            except ValueError as exc:
                raise ContractError("E_PATH_ESCAPE", f"{label} is outside run root: {path}") from exc
            if path.parent != self.run_root or path.name != expected_name:
                raise ContractError("E_PATH_ESCAPE", f"{label} must be {self.run_root / expected_name}")

        lock = require_mapping(load_json_strict(self.lock_path), "frozen lock")
        receipt = require_mapping(load_json_strict(self.receipt_path), "launch receipt")
        self.lock_file_sha256 = self.hash(self.lock_path)
        self.receipt_sha256 = self.hash(self.receipt_path)
        self._verify_lock(lock)
        self._verify_receipt(receipt)
        self.run_id = receipt["run_id"]
        self.add(
            "receipt_lock_file_digest",
            receipt["frozen_lock_sha256"] == self.lock_file_sha256,
            "E_FROZEN_LOCK_MISMATCH",
            receipt["frozen_lock_sha256"],
            self.lock_file_sha256,
        )
        self.add(
            "run_id_root_binding",
            self.run_id == self.run_root.name,
            "E_STATE_INVALID",
            self.run_id,
            self.run_root.name,
        )
        self._verify_frozen_provenance(lock, receipt)
        self._verify_state(receipt)
        sample_by_id = self._verify_dataset_contract(lock)
        self._verify_sample_directories(sample_by_id)
        self._verify_output_manifest_index()

        for sample in sorted(EXPECTED_DATASETS):
            if sample not in sample_by_id:
                continue
            start = len(self.checks)
            try:
                metrics = self._verify_sample(sample_by_id[sample], lock["analysis_contract"])
            except ContractError as exc:
                self.add("sample_contract", False, exc.code, str(exc), "valid sample", sample)
                metrics = {}
            sample_checks = self.checks[start:]
            self.samples.append(
                {
                    "sample": sample,
                    "biological_id": sample_by_id[sample]["biological_id"],
                    "pass": all(check.passed for check in sample_checks),
                    "error_codes": sorted({check.code for check in sample_checks if not check.passed}),
                    "metrics": metrics,
                }
            )

        marker = self.run_root / "_SUCCESS"
        self.add(
            "success_marker_not_precreated",
            not marker.exists() and not marker.is_symlink(),
            "E_STATE_INVALID",
            str(marker) if marker.exists() else None,
            "absent until RunLifecycle.succeed()",
        )
        return self._summary(lock)

    def _verify_lock(self, lock: Mapping[str, Any]) -> None:
        require_exact_keys(
            lock,
            {
                "$schema",
                "schema_name",
                "schema_version",
                "lock_id",
                "created_at_utc",
                "source_manifest",
                "validator",
                "dataset_count",
                "biological_sample_count",
                "analysis_contract",
                "production_summary",
                "all_pass",
                "samples",
                "lock_sha256",
            },
            "frozen lock",
        )
        schema_ok = (
            lock["$schema"] == LOCK_SCHEMA_REFERENCE
            and lock["schema_name"] == LOCK_SCHEMA_NAME
            and lock["schema_version"] == LOCK_SCHEMA_VERSION
        )
        self.add(
            "frozen_lock_schema",
            schema_ok,
            "E_SCHEMA_VERSION",
            [lock["$schema"], lock["schema_name"], lock["schema_version"]],
            [LOCK_SCHEMA_REFERENCE, LOCK_SCHEMA_NAME, LOCK_SCHEMA_VERSION],
        )
        require_sha256(lock["lock_sha256"], "lock.lock_sha256")
        digest_payload = dict(lock)
        declared = digest_payload.pop("lock_sha256")
        observed = canonical_sha256(digest_payload)
        self.add(
            "frozen_lock_embedded_digest",
            observed == declared,
            "E_FROZEN_LOCK_MISMATCH",
            observed,
            declared,
        )
        if lock["all_pass"] is not True:
            self.add("frozen_lock_all_pass", False, "E_STATE_INVALID", lock["all_pass"], True)
        require_list(lock["samples"], "lock.samples")
        analysis = require_mapping(lock["analysis_contract"], "lock.analysis_contract")
        self.tree_contract = str(analysis.get("tree_input_contract", ""))
        self.tree_profile = TREE_PROFILES.get(self.tree_contract, {})
        if not self.tree_profile:
            raise ContractError("E_TREE_VCF_ROLE", f"unsupported tree contract: {self.tree_contract}")
        require_mapping(lock["production_summary"], "lock.production_summary")

    def _verify_receipt(self, receipt: Mapping[str, Any]) -> None:
        require_exact_keys(
            receipt,
            {
                "schema_name",
                "schema_version",
                "run_id",
                "created_at_utc",
                "launcher_pid",
                "launcher_pid_start_time",
                "frozen_lock_sha256",
                "environment_lock_sha256",
                "source_bundle_manifest_sha256",
                "source_bundle_content_sha256",
                "extra",
            },
            "launch receipt",
        )
        if (
            receipt["schema_name"] != RECEIPT_SCHEMA_NAME
            or receipt["schema_version"] != RECEIPT_SCHEMA_VERSION
        ):
            raise ContractError("E_SCHEMA_VERSION", "unsupported lifecycle receipt schema")
        if not isinstance(receipt["run_id"], str) or SAFE_ID_RE.fullmatch(receipt["run_id"]) is None:
            raise ContractError("E_SCHEMA_INVALID", "unsafe receipt run_id")
        require_nonnegative_int(receipt["launcher_pid"], "receipt.launcher_pid")
        if not isinstance(receipt["launcher_pid_start_time"], str):
            raise ContractError("E_SCHEMA_INVALID", "receipt launcher PID start time must be string")
        for field in (
            "frozen_lock_sha256",
            "environment_lock_sha256",
            "source_bundle_manifest_sha256",
            "source_bundle_content_sha256",
        ):
            require_sha256(receipt[field], f"receipt.{field}")
        extra = require_mapping(receipt["extra"], "receipt.extra")
        require_exact_keys(
            extra,
            {
                "mode",
                "task_type",
                "analysis_params",
                "source_manifest_snapshot",
                "process_observation",
                "validator_execution",
                "worker_input_authority",
                "tree_backbone_role",
                "ledger_roles",
            },
            "receipt.extra",
        )
        params = require_mapping(extra["analysis_params"], "receipt.extra.analysis_params")
        profile_ok = (
            extra["mode"] == "full"
            and extra["task_type"] == self.tree_profile["task_type"]
            and extra["validator_execution"] == "bundled_core_via_runner_adapter"
            and extra["worker_input_authority"] == RUN_LOCK_FILENAME
            and extra["tree_backbone_role"] == self.tree_profile["backbone_role"]
            and extra["ledger_roles"]
            == [
                "caller_raw_vcf",
                "longphase_input_vcf",
                "caller_pass_baseline_vcf",
                "longphase_recalibrated_all_vcf",
                "longphase_recalibrated_pass_vcf",
                "tree_vcf",
            ]
            and params.get("scope") == "chr1-22"
            and params.get("contigs") == list(AUTOSOMES)
            and params.get("dataset_count") == 7
            and params.get("biological_sample_count") == 6
            and params.get("VERIFY_EVERY") == 1
            and params.get("ANALYSIS_TREE_CAP") == 0
            and params.get("parallel_parts_per_sample") == 1
            and isinstance(params.get("parallel_samples"), int)
            and 1 <= params["parallel_samples"] <= 4
        )
        self.add(
            "receipt_comprehensive_profile",
            profile_ok,
            "E_STATE_INVALID",
            extra,
            "full comprehensive chr1-22/7/6; exhaustive tree verification; explicit VCF roles",
        )

    def _verify_full_identity(self, spec: Mapping[str, Any], label: str, sample: str | None = None) -> Path:
        require_exact_keys(spec, {"path", "identity"}, label)
        path_text = spec["path"]
        if not isinstance(path_text, str) or not path_text:
            raise ContractError("E_SCHEMA_INVALID", f"{label}.path is invalid")
        path = regular_file(Path(path_text), label, reject_symlink=False)
        identity = require_mapping(spec["identity"], f"{label}.identity")
        require_exact_keys(identity, {"policy", "size_bytes", "sha256"}, f"{label}.identity")
        if identity["policy"] != "full_sha256":
            raise ContractError("E_SCHEMA_INVALID", f"{label} is not full_sha256")
        expected_size = require_nonnegative_int(identity["size_bytes"], f"{label}.size_bytes")
        expected_sha = require_sha256(identity["sha256"], f"{label}.sha256")
        actual_size = path.stat().st_size
        actual_sha = self.hash(path)
        self.add(
            f"identity:{label}",
            actual_size == expected_size and actual_sha == expected_sha,
            "E_POST_INPUT_IDENTITY",
            {"size_bytes": actual_size, "sha256": actual_sha},
            {"size_bytes": expected_size, "sha256": expected_sha},
            sample,
        )
        return path

    def _verify_indexed_artifact(self, spec: Mapping[str, Any], label: str, sample: str) -> None:
        require_exact_keys(spec, {"path", "identity", "index"}, label)
        self._verify_full_identity({"path": spec["path"], "identity": spec["identity"]}, label, sample)
        self._verify_full_identity(require_mapping(spec["index"], f"{label}.index"), f"{label}.index", sample)

    def _verify_frozen_provenance(self, lock: Mapping[str, Any], receipt: Mapping[str, Any]) -> None:
        source = require_mapping(lock["source_manifest"], "lock.source_manifest")
        require_exact_keys(source, {"path", "byte_sha256", "canonical_sha256"}, "lock.source_manifest")
        source_path = regular_file(Path(source["path"]), "source manifest", reject_symlink=False)
        source_doc = load_json_strict(source_path, "E_SOURCE_MANIFEST_MISMATCH")
        source_bytes_sha = self.hash(source_path)
        source_canonical_sha = canonical_sha256(source_doc)
        self.add(
            "source_manifest_identity",
            source_bytes_sha == source["byte_sha256"] and source_canonical_sha == source["canonical_sha256"],
            "E_SOURCE_MANIFEST_MISMATCH",
            {"byte_sha256": source_bytes_sha, "canonical_sha256": source_canonical_sha},
            source,
        )
        snapshot_spec = receipt["extra"]["source_manifest_snapshot"]
        if not isinstance(snapshot_spec, dict) or set(snapshot_spec) != {"path", "sha256"}:
            raise ContractError("E_SOURCE_MANIFEST_MISMATCH", "source manifest snapshot spec is malformed")
        snapshot_path = contained_file(
            snapshot_spec["path"], self.run_root, "source manifest snapshot"
        )
        snapshot_sha = self.hash(snapshot_path)
        self.add(
            "source_manifest_snapshot_binding",
            snapshot_sha == snapshot_spec["sha256"] == source["byte_sha256"],
            "E_SOURCE_MANIFEST_MISMATCH",
            snapshot_sha,
            source["byte_sha256"],
        )
        process_spec = receipt["extra"]["process_observation"]
        if not isinstance(process_spec, dict) or set(process_spec) != {"path", "sha256", "conflict_count"}:
            raise ContractError("E_STATE_INVALID", "process observation receipt spec is malformed")
        process_path = contained_file(process_spec["path"], self.run_root, "process observation")
        process_document = require_mapping(
            load_json_strict(process_path, "E_STATE_INVALID"), "process observation"
        )
        process_lock = process_document.get("global_scope_lock")
        process_ok = (
            process_spec["conflict_count"] == 0
            and self.hash(process_path) == process_spec["sha256"]
            and process_document.get("schema_name") == "intersubmod.layered_process_observation"
            and process_document.get("schema_version") == "1.0.0"
            and process_document.get("observer_pid") == receipt["launcher_pid"]
            and process_document.get("observer_pid_start_time") == receipt["launcher_pid_start_time"]
            and isinstance(process_document.get("processes_inspected"), int)
            and process_document.get("processes_inspected") >= 0
            and isinstance(process_document.get("match_policy"), list)
            and bool(process_document["match_policy"])
            and process_document.get("conflicts") == []
            and process_document.get("pass") is True
            and isinstance(process_lock, dict)
            and process_lock.get("held") is True
            and process_lock.get("run_id") == receipt["run_id"]
            and Path(str(process_lock.get("path", ""))).name == ".layered_chr1_22_7dataset_full.lock"
        )
        lock_owner_ok = False
        if isinstance(process_lock, dict):
            lock_path = Path(str(process_lock.get("path", "")))
            try:
                lock_owner_ok = lock_path.is_file() and f"run_id={receipt['run_id']}" in lock_path.read_text(
                    encoding="utf-8", errors="strict"
                )
            except (OSError, UnicodeError):
                lock_owner_ok = False
        self.add(
            "process_exclusivity_observation",
            process_ok and lock_owner_ok,
            "E_STATE_INVALID",
            {"document": process_document, "lock_owner_ok": lock_owner_ok},
            "hash-bound zero-conflict observation under held full-scope lock",
        )

        validator = require_mapping(lock["validator"], "lock.validator")
        require_exact_keys(validator, {"path", "sha256", "schema_path", "schema_sha256"}, "lock.validator")
        for path_field, sha_field in (("path", "sha256"), ("schema_path", "schema_sha256")):
            path = regular_file(Path(validator[path_field]), f"validator.{path_field}", reject_symlink=False)
            actual = self.hash(path)
            self.add(
                f"validator_{path_field}_identity",
                actual == validator[sha_field],
                "E_SOURCE_BUNDLE_MISMATCH",
                actual,
                validator[sha_field],
            )

        bundle_manifest_path = regular_file(
            self.run_root / SOURCE_BUNDLE_MANIFEST, "source bundle manifest"
        )
        bundle_manifest = require_mapping(
            load_json_strict(bundle_manifest_path, "E_SOURCE_BUNDLE_MISMATCH"), "source bundle manifest"
        )
        require_exact_keys(
            bundle_manifest,
            {"schema_name", "schema_version", "created_at_utc", "file_count", "bundle_content_sha256", "files"},
            "source bundle manifest",
        )
        if (
            bundle_manifest["schema_name"] != SOURCE_BUNDLE_SCHEMA_NAME
            or bundle_manifest["schema_version"] != SOURCE_BUNDLE_SCHEMA_VERSION
        ):
            raise ContractError("E_SOURCE_BUNDLE_MISMATCH", "unsupported source bundle schema")
        files = require_list(bundle_manifest["files"], "source bundle files")
        require_nonnegative_int(bundle_manifest["file_count"], "source bundle file_count")
        roles: set[str] = set()
        content_identity = []
        expected_paths = {bundle_manifest_path.resolve()}
        bundle_root = bundle_manifest_path.parent.resolve()
        for index, raw in enumerate(files):
            item = require_mapping(raw, f"source bundle files[{index}]")
            require_exact_keys(
                item,
                {
                    "role", "source_path", "bundled_path", "mode", "size", "sha256",
                    "source_device", "source_inode", "source_mtime_ns", "bundled_device", "bundled_inode",
                },
                f"source bundle files[{index}]",
            )
            role = item["role"]
            if not isinstance(role, str) or role in roles:
                raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"duplicate/invalid bundle role: {role!r}")
            roles.add(role)
            relative = Path(item["bundled_path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"unsafe bundled path: {relative}")
            path = regular_file(bundle_root / relative, f"source bundle:{role}")
            try:
                path.relative_to(bundle_root)
            except ValueError as exc:
                raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"bundle path escapes: {relative}") from exc
            actual = {"size": path.stat().st_size, "sha256": self.hash(path), "mode": stat.S_IMODE(path.stat().st_mode)}
            expected = {"size": item["size"], "sha256": item["sha256"], "mode": item["mode"]}
            self.add(
                f"source_bundle_file:{role}",
                actual == expected,
                "E_SOURCE_BUNDLE_MISMATCH",
                actual,
                expected,
            )
            expected_paths.add(path)
            content_identity.append(
                {
                    "role": role,
                    "bundled_path": item["bundled_path"],
                    "mode": item["mode"],
                    "size": item["size"],
                    "sha256": item["sha256"],
                }
            )
        observed_paths = {path.resolve() for path in bundle_root.rglob("*") if path.is_file()}
        self.add(
            "source_bundle_file_set",
            observed_paths == expected_paths,
            "E_SOURCE_BUNDLE_MISMATCH",
            sorted(str(path.relative_to(bundle_root)) for path in observed_paths),
            sorted(str(path.relative_to(bundle_root)) for path in expected_paths),
        )
        required_roles = {"runner", "validator", "verifier"}
        roles_ok = required_roles.issubset(roles) and any(role.startswith("imported:") for role in roles)
        self.add("source_bundle_roles", roles_ok, "E_SOURCE_BUNDLE_MISMATCH", sorted(roles), "core+imported")
        content_sha = canonical_sha256(content_identity)
        manifest_sha = self.hash(bundle_manifest_path)
        self.bundle_manifest_sha256 = manifest_sha
        self.bundle_content_sha256 = content_sha
        bundle_ok = (
            len(files) == bundle_manifest["file_count"]
            and content_sha == bundle_manifest["bundle_content_sha256"]
            and manifest_sha == receipt["source_bundle_manifest_sha256"]
            and content_sha == receipt["source_bundle_content_sha256"]
        )
        self.add(
            "source_bundle_receipt_binding",
            bundle_ok,
            "E_SOURCE_BUNDLE_MISMATCH",
            {"manifest_sha256": manifest_sha, "content_sha256": content_sha, "file_count": len(files)},
            {
                "manifest_sha256": receipt["source_bundle_manifest_sha256"],
                "content_sha256": receipt["source_bundle_content_sha256"],
                "file_count": bundle_manifest["file_count"],
            },
        )

        environment_path = regular_file(self.run_root / ENVIRONMENT_LOCK_FILENAME, "environment lock")
        environment = require_mapping(
            load_json_strict(environment_path, "E_ENVIRONMENT_MISMATCH"), "environment lock"
        )
        environment_sha = self.hash(environment_path)
        self.environment_sha256 = environment_sha
        environment_ok = (
            environment.get("schema_name") == ENVIRONMENT_SCHEMA_NAME
            and environment.get("schema_version") == ENVIRONMENT_SCHEMA_VERSION
            and environment_sha == receipt["environment_lock_sha256"]
        )
        self.add(
            "environment_receipt_binding",
            environment_ok,
            "E_ENVIRONMENT_MISMATCH",
            {
                "schema_name": environment.get("schema_name"),
                "schema_version": environment.get("schema_version"),
                "sha256": environment_sha,
            },
            {
                "schema_name": ENVIRONMENT_SCHEMA_NAME,
                "schema_version": ENVIRONMENT_SCHEMA_VERSION,
                "sha256": receipt["environment_lock_sha256"],
            },
        )

    def _state_digest_valid(self, state: Mapping[str, Any]) -> bool:
        if "state_sha256" not in state or not isinstance(state["state_sha256"], str):
            return False
        payload = dict(state)
        declared = payload.pop("state_sha256")
        return canonical_sha256(payload) == declared

    def _verify_state(self, receipt: Mapping[str, Any]) -> None:
        state_path = regular_file(self.run_root / RUN_STATE_FILENAME, "run state")
        state = require_mapping(load_json_strict(state_path, "E_STATE_INVALID"), "run state")
        state_fields = {
            "schema_name", "schema_version", "run_id", "sequence", "previous_state_sha256",
            "timestamp_utc", "state", "stage", "sample", "launcher_pid",
            "launcher_pid_start_time", "reason", "error_code", "state_sha256",
        }
        require_exact_keys(state, state_fields, "run state")
        expected_bindings = {
            "schema_name": STATE_SCHEMA_NAME,
            "schema_version": STATE_SCHEMA_VERSION,
            "run_id": receipt["run_id"],
            "state": "VERIFYING",
            "stage": "verifying",
            "launcher_pid": receipt["launcher_pid"],
            "launcher_pid_start_time": receipt["launcher_pid_start_time"],
        }
        self.add(
            "run_state_bindings",
            all(state.get(key) == value for key, value in expected_bindings.items()) and self._state_digest_valid(state),
            "E_STATE_INVALID",
            {key: state.get(key) for key in expected_bindings} | {"digest_valid": self._state_digest_valid(state)},
            expected_bindings | {"digest_valid": True},
        )
        sequence = require_nonnegative_int(state["sequence"], "run state.sequence")
        events_root = self.run_root / "state_events"
        if events_root.is_symlink() or not events_root.is_dir():
            raise ContractError("E_STATE_INVALID", f"state_events is missing/symlink: {events_root}")
        entries = sorted(events_root.iterdir())
        expected_previous = ZERO_DIGEST
        event_documents = []
        chain_ok = len(entries) == sequence
        for number, event_path in enumerate(entries, start=1):
            if event_path.is_symlink() or not event_path.is_file():
                chain_ok = False
                continue
            event = require_mapping(load_json_strict(event_path, "E_STATE_INVALID"), f"state event {number}")
            event_documents.append(event)
            expected_name = f"{number:06d}_{event.get('state')}.json"
            chain_ok &= event_path.name == expected_name
            chain_ok &= event.get("sequence") == number
            chain_ok &= event.get("previous_state_sha256") == expected_previous
            chain_ok &= self._state_digest_valid(event)
            expected_previous = event.get("state_sha256", "")
        chain_ok &= bool(event_documents) and event_documents[-1] == state
        self.add(
            "run_state_event_chain",
            chain_ok,
            "E_STATE_INVALID",
            {"events": len(entries), "sequence": sequence, "last_digest": expected_previous},
            {"events": sequence, "sequence": sequence, "last_digest": state["state_sha256"]},
        )

    def _verify_dataset_contract(self, lock: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
        samples = require_list(lock["samples"], "lock.samples")
        ids = []
        biological = []
        for index, raw in enumerate(samples):
            sample = require_mapping(raw, f"lock.samples[{index}]")
            require_exact_keys(
                sample,
                {"sample", "biological_id", "pass", "alignment_payload", "somatic", "read_tags", "copy_number"},
                f"lock.samples[{index}]",
            )
            if not isinstance(sample["sample"], str) or SAFE_ID_RE.fullmatch(sample["sample"]) is None:
                raise ContractError("E_SAMPLE_ID_UNSAFE", f"unsafe sample: {sample.get('sample')!r}")
            if not isinstance(sample["biological_id"], str) or SAFE_ID_RE.fullmatch(sample["biological_id"]) is None:
                raise ContractError("E_SAMPLE_ID_UNSAFE", f"unsafe biological ID: {sample.get('biological_id')!r}")
            ids.append(sample["sample"])
            biological.append(sample["biological_id"])
            self.add(
                "sample_lock_pass",
                sample["pass"] is True,
                "E_STATE_INVALID",
                sample["pass"],
                True,
                sample["sample"],
            )
        self.add(
            "sample_ids_unique",
            len(ids) == len(set(ids)),
            "E_SAMPLE_DUPLICATE",
            ids,
            "unique",
        )
        self.add(
            "exact_dataset_set",
            len(ids) == 7 and set(ids) == EXPECTED_DATASETS,
            "E_DATASET_SET_MISMATCH",
            ids,
            sorted(EXPECTED_DATASETS),
        )
        self.add(
            "exact_biological_set",
            len(set(biological)) == 6 and set(biological) == EXPECTED_BIOLOGICAL_IDS,
            "E_DATASET_SET_MISMATCH",
            sorted(set(biological)),
            sorted(EXPECTED_BIOLOGICAL_IDS),
        )
        bindings = dict(zip(ids, biological))
        self.add(
            "dataset_biological_binding",
            all(EXPECTED_BINDING.get(sample) == bio for sample, bio in bindings.items()),
            "E_SAMPLE_OUTPUT_BINDING",
            bindings,
            EXPECTED_BINDING,
        )
        self.add(
            "declared_counts",
            lock["dataset_count"] == 7 == len(ids)
            and lock["biological_sample_count"] == 6 == len(set(biological)),
            "E_DATASET_SET_MISMATCH",
            [lock["dataset_count"], lock["biological_sample_count"], len(ids), len(set(biological))],
            [7, 6, 7, 6],
        )
        analysis = lock["analysis_contract"]
        scope = analysis.get("scope") if isinstance(analysis, dict) else None
        contract_ok = (
            analysis.get("task_type") == self.tree_profile["task_type"]
            and isinstance(scope, dict)
            and scope.get("name") == "whole_autosomes_chr1_22"
            and scope.get("contigs") == list(AUTOSOMES)
            and analysis.get("read_tag_mode") == "external_sidecar"
            and analysis.get("embedded_tag_policy") == "ignore"
            and analysis.get("require_exact_join") is True
            and analysis.get("sidecar_identity_schema") == "coordinate_join_v1"
            and analysis.get("sidecar_assurance") == "bounded_coordinate_equivalence_not_global_payload_identity"
            and analysis.get("bam_identity_policy") == "storage_identity_v1"
            and analysis.get("longphase_input_contract") == "normalized_ClairS_raw_all"
            and analysis.get("tree_input_contract") == self.tree_contract
            and analysis.get("tagging_semantics")
            == "longphase_s_raw_all_bidirectional_recalibration_tags"
            and analysis.get("duplicate_identity_policy")
            == "collapse_redundant_rows_with_identical_HP_PS"
        )
        self.add(
            "analysis_contract",
            contract_ok,
            "E_SIDECAR_CONTRACT_UNSUPPORTED",
            analysis,
            "comprehensive chr1-22 lower-assurance coordinate_join_v1/storage_identity_v1",
        )
        method = analysis.get("join_method_validation")
        method_ok = (
            isinstance(method, dict)
            and method.get("assurance") == "bounded_coordinate_equivalence_not_global_payload_identity"
            and method.get("claim_limit") == "join_method_only_not_per_sample_global_payload_identity"
        )
        if isinstance(method, dict):
            for role in ("real_data_bounded_receipt", "synthetic_three_case_receipt"):
                spec = method.get(role)
                if isinstance(spec, dict):
                    self._verify_full_identity(spec, f"join_method.{role}")
                else:
                    method_ok = False
        self.add(
            "join_method_receipts",
            method_ok,
            "E_SIDECAR_CONTRACT_UNSUPPORTED",
            method,
            "bounded real-data receipt + three-case synthetic receipt",
        )
        production = lock["production_summary"]
        production_datasets = production.get("datasets") if isinstance(production, dict) else None
        production_ids = {
            item.get("sample")
            for item in production_datasets or []
            if isinstance(item, dict)
        }
        production_rows_ok = isinstance(production_datasets, list) and all(
            isinstance(item, dict)
            and item.get("sample") in EXPECTED_DATASETS
            and item.get("pass") is True
            and item.get("filter_policy") == "production_default_filter"
            and item.get("truth_flags_present") is False
            and item.get("identity_schema") == "coordinate_join_v1"
            and item.get("longphase_input_contract") == "normalized_ClairS_raw_all"
            and item.get("tree_backbone_role") == "longphase_s_recalibrated_filter_pass"
            and isinstance(item.get("validation_sha256"), str)
            and SHA256_RE.fullmatch(item["validation_sha256"]) is not None
            for item in production_datasets
        )
        production_ok = (
            production.get("expected_dataset_count") == 7
            and production.get("completed_dataset_count") == 7
            and production.get("passed_dataset_count") == 7
            and production.get("all_pass") is True
            and isinstance(production.get("datasets"), list)
            and len(production["datasets"]) == 7
            and production_ids == EXPECTED_DATASETS
            and production_rows_ok
        )
        self.add(
            "production_summary",
            production_ok,
            "E_DATASET_SET_MISMATCH",
            production,
            "7/7 completed and passed",
        )
        return {sample["sample"]: sample for sample in samples if isinstance(sample, dict)}

    def _verify_sample_directories(self, sample_by_id: Mapping[str, Any]) -> None:
        samples_root = self.run_root / "samples"
        if samples_root.is_symlink() or not samples_root.is_dir():
            raise ContractError("E_DATASET_SET_MISMATCH", f"samples root missing/symlink: {samples_root}")
        directories = []
        invalid = []
        for child in samples_root.iterdir():
            if child.is_symlink() or not child.is_dir():
                invalid.append(child.name)
            else:
                directories.append(child.name)
        self.add(
            "sample_directories_exact",
            not invalid and set(directories) == EXPECTED_DATASETS == set(sample_by_id),
            "E_DATASET_SET_MISMATCH",
            {"directories": sorted(directories), "invalid": sorted(invalid)},
            sorted(EXPECTED_DATASETS),
        )

    def _verify_output_manifest_index(self) -> None:
        path = regular_file(self.run_root / "output_manifests.json", "output manifest index")
        document = require_mapping(load_json_strict(path), "output manifest index")
        require_exact_keys(document, {"dataset_count", "manifests"}, "output manifest index")
        manifests = require_list(document["manifests"], "output manifest index.manifests")
        observed: dict[str, str] = {}
        all_hashes_ok = True
        for index, raw in enumerate(manifests):
            item = require_mapping(raw, f"output manifest index[{index}]")
            require_exact_keys(item, {"sample", "path", "sha256"}, f"output manifest index[{index}]")
            sample = item["sample"]
            if not isinstance(sample, str) or sample in observed:
                all_hashes_ok = False
                continue
            expected_relative = f"samples/{sample}/output_manifest.json"
            manifest_path = contained_file(item["path"], self.run_root, f"output manifest:{sample}")
            actual_sha = self.hash(manifest_path)
            all_hashes_ok &= item["path"] == expected_relative
            all_hashes_ok &= item["sha256"] == actual_sha
            observed[sample] = actual_sha
        self.add(
            "output_manifest_index",
            document["dataset_count"] == 7
            and len(manifests) == 7
            and set(observed) == EXPECTED_DATASETS
            and all_hashes_ok,
            "E_OUTPUT_HASH_MISMATCH",
            {"dataset_count": document["dataset_count"], "samples": sorted(observed), "hashes_ok": all_hashes_ok},
            {"dataset_count": 7, "samples": sorted(EXPECTED_DATASETS), "hashes_ok": True},
        )

    def _verify_storage_identity(self, alignment: Mapping[str, Any], sample: str) -> tuple[Path, str]:
        require_exact_keys(alignment, {"path", "storage_identity_v1"}, f"{sample}.alignment_payload")
        expected = require_mapping(alignment["storage_identity_v1"], f"{sample}.storage_identity_v1")
        required = {
            "policy", "assurance", "is_full_content_hash", "requested_path", "realpath",
            "logical_is_symlink", "logical_size_bytes", "logical_mtime_ns", "st_dev", "st_ino",
            "size_bytes", "mtime_ns", "ctime_ns", "chunk_size_bytes", "chunks", "index", "identity_sha256",
        }
        require_exact_keys(expected, required, f"{sample}.storage_identity_v1")
        path = Path(alignment["path"])
        try:
            logical = path.lstat()
            resolved = path.resolve(strict=True)
            target = resolved.stat()
        except OSError as exc:
            raise ContractError("E_POST_INPUT_IDENTITY", f"cannot observe BAM identity: {path}: {exc}") from exc
        chunks = []
        for raw in require_list(expected["chunks"], f"{sample}.storage_identity_v1.chunks"):
            chunk = require_mapping(raw, "storage chunk")
            require_exact_keys(chunk, {"label", "offset", "length", "sha256"}, "storage chunk")
            offset = require_nonnegative_int(chunk["offset"], "storage chunk.offset")
            length = require_nonnegative_int(chunk["length"], "storage chunk.length")
            chunks.append(
                {
                    "label": chunk["label"],
                    "offset": offset,
                    "length": length,
                    "sha256": sha256_slice(resolved, offset, length),
                }
            )
        index = require_mapping(expected["index"], f"{sample}.storage index")
        self._verify_full_identity(index, f"{sample}.alignment_payload.index", sample)
        observed = {
            "policy": "storage_identity_v1",
            "assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
            "is_full_content_hash": False,
            "requested_path": str(path),
            "realpath": str(resolved),
            "logical_is_symlink": stat.S_ISLNK(logical.st_mode),
            "logical_size_bytes": logical.st_size,
            "logical_mtime_ns": logical.st_mtime_ns,
            "st_dev": target.st_dev,
            "st_ino": target.st_ino,
            "size_bytes": target.st_size,
            "mtime_ns": target.st_mtime_ns,
            "ctime_ns": target.st_ctime_ns,
            "chunk_size_bytes": expected["chunk_size_bytes"],
            "chunks": chunks,
            "index": {"path": index["path"], "identity": index["identity"]},
        }
        observed["identity_sha256"] = canonical_sha256(observed)
        self.add(
            "alignment_payload_storage_identity",
            observed == expected,
            "E_POST_INPUT_IDENTITY",
            observed,
            expected,
            sample,
        )
        return resolved, expected["identity_sha256"]

    def _scan_coordinate_sidecar(self, path: Path, sample: str) -> dict[str, int]:
        rows = duplicates = conflicts = malformed = 0
        current_chrom: str | None = None
        current_start = -1
        finished_contigs: set[str] = set()
        bucket: dict[tuple[str, ...], tuple[str, str]] = {}
        sorted_ok = True
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
                header = tuple(handle.readline().rstrip("\n").split("\t"))
                if header != SIDECAR_HEADER:
                    raise ContractError("E_SIDECAR_PARSE", f"unexpected sidecar header: {header}")
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 9:
                        malformed += 1
                        continue
                    try:
                        start = int(fields[1])
                        int(fields[2]); int(fields[4]); int(fields[5])
                    except ValueError:
                        malformed += 1
                        continue
                    chrom = fields[0]
                    if chrom != current_chrom:
                        if current_chrom is not None:
                            finished_contigs.add(current_chrom)
                        sorted_ok &= chrom not in finished_contigs
                        current_chrom, current_start, bucket = chrom, start, {}
                    elif start != current_start:
                        sorted_ok &= start >= current_start
                        current_start, bucket = start, {}
                    key = (fields[3], chrom, fields[1], fields[2], fields[4], fields[6])
                    tag = (fields[7], fields[8])
                    if key in bucket:
                        duplicates += 1
                        conflicts += int(bucket[key] != tag)
                    else:
                        bucket[key] = tag
                    rows += 1
        except (OSError, UnicodeError) as exc:
            raise ContractError("E_SIDECAR_PARSE", f"cannot scan sidecar {path}: {exc}") from exc
        self.add(
            "sidecar_coordinate_order",
            sorted_ok,
            "E_SIDECAR_PARSE",
            sorted_ok,
            True,
            sample,
        )
        return {
            "mapped_alignment_count": rows,
            "identity_unique_count": rows - duplicates,
            "duplicate_count": duplicates,
            "conflict_count": conflicts,
            "malformed_count": malformed,
        }

    def _verify_read_tags(
        self,
        meta: Mapping[str, Any],
        analysis: Mapping[str, Any],
        alignment_identity_sha: str,
        somatic_hashes: Mapping[str, str],
        sample: str,
    ) -> dict[str, Any]:
        require_exact_keys(
            meta,
            {
                "sidecar",
                "index",
                "validation",
                "producer_capture_receipt",
                "duplicate_identity_policy",
                "subject_binding",
                "producer_policy",
                "producer_evidence",
            },
            f"{sample}.read_tags",
        )
        sidecar_path = self._verify_full_identity(meta["sidecar"], f"{sample}.read_tags.sidecar", sample)
        self._verify_full_identity(meta["index"], f"{sample}.read_tags.index", sample)
        validation_path = self._verify_full_identity(
            meta["validation"], f"{sample}.read_tags.validation", sample
        )
        producer_receipt_path = self._verify_full_identity(
            meta["producer_capture_receipt"],
            f"{sample}.read_tags.producer_capture_receipt",
            sample,
        )
        sidecar_sha = meta["sidecar"]["identity"]["sha256"]
        index_sha = meta["index"]["identity"]["sha256"]
        validation_sha = meta["validation"]["identity"]["sha256"]
        producer_receipt_sha = meta["producer_capture_receipt"]["identity"]["sha256"]
        subject = require_mapping(meta["subject_binding"], f"{sample}.read_tags.subject_binding")
        evidence = require_mapping(meta["producer_evidence"], f"{sample}.read_tags.producer_evidence")
        scan = self._scan_coordinate_sidecar(sidecar_path, sample)
        expected_subject = {
            "schema_name": "intersubmod.coordinate_join_subject",
            "schema_version": "1.0.0",
            "sample": sample,
            "duplicate_identity_policy": "collapse_redundant_rows_with_identical_HP_PS",
            "coordinate_identity_columns": list(COORDINATE_IDENTITY_COLUMNS),
            "sidecar_sha256": sidecar_sha,
            "sidecar_index_sha256": index_sha,
            "validation_sha256": validation_sha,
            "producer_capture_receipt_sha256": producer_receipt_sha,
            "alignment_payload_storage_identity_sha256": alignment_identity_sha,
            "producer_command_argv_sha256": evidence.get("command_argv_sha256"),
            "producer_input_binding_sha256": evidence.get("input_binding_sha256"),
            "producer_effective_options_sha256": evidence.get("effective_options_sha256"),
            "caller_raw_vcf_sha256": somatic_hashes["caller_raw_vcf"],
            "longphase_input_vcf_sha256": somatic_hashes["longphase_input_vcf"],
            "caller_pass_baseline_vcf_sha256": somatic_hashes["caller_pass_baseline_vcf"],
            "longphase_recalibrated_all_vcf_sha256": somatic_hashes["longphase_recalibrated_all_vcf"],
            "longphase_recalibrated_pass_vcf_sha256": somatic_hashes[
                "longphase_recalibrated_pass_vcf"
            ],
            "mapped_alignment_count": scan["mapped_alignment_count"],
            "identity_unique_count": scan["identity_unique_count"],
            "duplicate_count": scan["duplicate_count"],
            "conflict_count": scan["conflict_count"],
        }
        multiplicity_ok = (
            scan["mapped_alignment_count"]
            == scan["identity_unique_count"] + scan["duplicate_count"]
            and scan["conflict_count"] == 0
        )
        subject_ok = subject == expected_subject and scan["malformed_count"] == 0 and multiplicity_ok
        self.add(
            "sidecar_subject_binding",
            subject_ok,
            "E_SIDECAR_SUBJECT_MISMATCH",
            {"subject": subject, "scan": scan},
            expected_subject,
            sample,
        )

        validation = require_mapping(
            load_json_strict(validation_path, "E_SIDECAR_SUBJECT_MISMATCH"),
            f"{sample} sidecar validation",
        )
        receipt = require_mapping(
            load_json_strict(producer_receipt_path, "E_PRODUCER_RECEIPT_UNSUPPORTED"),
            f"{sample} producer receipt",
        )
        argv = receipt.get("command_argv")
        argv_ok = isinstance(argv, list) and bool(argv) and all(isinstance(value, str) for value in argv)
        forbidden = ("--truth-vcf", "--truth-bed", "--disableFilter")
        forbidden_present = argv_ok and any(
            token == option or token.startswith(f"{option}=")
            for token in argv
            for option in forbidden
        )
        policy = meta["producer_policy"]
        required_native_checks = (
            "truth_flags_absent",
            "parser_count_matches_input",
            "sidecar_row_count_matches_capture",
            "tagged_count_matches_execution",
            "sidecar_no_unknown_HP",
            "sidecar_no_exact_identity_conflicts",
            "recalibrated_preserves_all_input_keys",
        )
        native_ok = (
            validation.get("sample") in (None, sample)
            and validation.get("region") == "all"
            and validation.get("pass") is True
            and validation.get("duplicate_exact_alignment_rows") == scan["duplicate_count"]
            and validation.get("duplicate_exact_alignment_conflicts") == 0 == scan["conflict_count"]
            and all(validation.get("checks", {}).get(name) is True for name in required_native_checks)
        )
        receipt_inputs = receipt.get("producer_inputs")
        receipt_counts = receipt.get("global_coordinate_counts")
        receipt_outputs = receipt.get("capture_outputs")
        if isinstance(receipt_inputs, dict):
            require_exact_keys(
                receipt_inputs,
                {
                    "germline_phased_vcf", "normal_bam", "tumor_bam", "caller_raw_vcf",
                    "longphase_input_vcf", "caller_pass_baseline_vcf", "reference",
                },
                f"{sample}.producer_inputs",
            )
            self._verify_indexed_artifact(
                receipt_inputs["germline_phased_vcf"], f"{sample}.producer.germline_phased_vcf", sample
            )
            for vcf_role in ("caller_raw_vcf", "longphase_input_vcf", "caller_pass_baseline_vcf"):
                self._verify_indexed_artifact(
                    receipt_inputs[vcf_role], f"{sample}.producer.{vcf_role}", sample
                )
            for bam_role in ("normal_bam", "tumor_bam"):
                storage_spec = require_mapping(receipt_inputs[bam_role], f"{sample}.producer.{bam_role}")
                require_exact_keys(storage_spec, {"path", "index_path", "storage_identity_v1"}, f"{sample}.{bam_role}")
                self._verify_storage_identity(
                    {"path": storage_spec["path"], "storage_identity_v1": storage_spec["storage_identity_v1"]},
                    sample,
                )
                self.add(
                    f"producer_index_path:{bam_role}",
                    storage_spec["index_path"] == storage_spec["storage_identity_v1"]["index"]["path"],
                    "E_POST_INPUT_IDENTITY",
                    storage_spec["index_path"],
                    storage_spec["storage_identity_v1"]["index"]["path"],
                    sample,
                )
            reference = require_mapping(receipt_inputs["reference"], f"{sample}.producer.reference")
            require_exact_keys(reference, {"path", "fai_path", "storage_identity_v1"}, f"{sample}.reference")
            self._verify_storage_identity(
                {"path": reference["path"], "storage_identity_v1": reference["storage_identity_v1"]},
                sample,
            )
            self.add(
                "producer_reference_index_path",
                reference["fai_path"] == reference["storage_identity_v1"]["index"]["path"],
                "E_POST_INPUT_IDENTITY",
                reference["fai_path"],
                reference["storage_identity_v1"]["index"]["path"],
                sample,
            )
        if isinstance(receipt_outputs, dict):
            require_exact_keys(
                receipt_outputs,
                {
                    "sidecar", "sidecar_index", "native_validation", "stream_capture_summary",
                    "normalization_audit", "filter_transition_audit", "sample_verification",
                    "longphase_recalibrated_all_vcf", "longphase_recalibrated_pass_vcf",
                },
                f"{sample}.capture_outputs",
            )
            for output_role in (
                "sidecar", "sidecar_index", "native_validation", "stream_capture_summary",
                "normalization_audit", "filter_transition_audit", "sample_verification",
            ):
                self._verify_full_identity(
                    receipt_outputs[output_role], f"{sample}.capture.{output_role}", sample
                )
            for output_role in ("longphase_recalibrated_all_vcf", "longphase_recalibrated_pass_vcf"):
                self._verify_indexed_artifact(
                    receipt_outputs[output_role], f"{sample}.capture.{output_role}", sample
                )
        longphase_executable = receipt.get("longphase_executable")
        if isinstance(longphase_executable, dict):
            require_exact_keys(longphase_executable, {"path", "identity", "version"}, "longphase executable")
            self._verify_full_identity(
                {"path": longphase_executable["path"], "identity": longphase_executable["identity"]},
                f"{sample}.longphase_executable",
                sample,
            )
        patch_evidence = receipt.get("patch_evidence")
        if isinstance(patch_evidence, dict):
            require_exact_keys(
                patch_evidence, {"approval_scope", "patch_receipt", "source_patch"}, "patch evidence"
            )
            self._verify_full_identity(
                patch_evidence["patch_receipt"], f"{sample}.patch_receipt", sample
            )
            self._verify_full_identity(
                patch_evidence["source_patch"], f"{sample}.source_patch", sample
            )
        producer_code = receipt.get("producer_code")
        if isinstance(producer_code, dict):
            require_exact_keys(
                producer_code, {"runner", "capture", "validator", "filter_auditor"}, "producer code"
            )
            for code_role, code_spec in producer_code.items():
                self._verify_full_identity(code_spec, f"{sample}.producer_code.{code_role}", sample)
        producer_environment = receipt.get("environment_lock")
        if isinstance(producer_environment, dict):
            require_exact_keys(
                producer_environment,
                {
                    "production_manifest",
                    "run_params",
                    "code_hash_manifest",
                    "sample_input_inventory",
                    "sample_input_hash_manifest",
                    "sample_output_hash_manifest",
                    "python_executable",
                    "python_version",
                    "platform",
                    "normalizer_source",
                    "normalizer_argv",
                },
                f"{sample}.producer_environment",
            )
            for artifact_role in (
                "production_manifest",
                "run_params",
                "code_hash_manifest",
                "sample_input_inventory",
                "sample_input_hash_manifest",
                "sample_output_hash_manifest",
                "python_executable",
                "normalizer_source",
            ):
                self._verify_full_identity(
                    producer_environment[artifact_role],
                    f"{sample}.producer_environment.{artifact_role}",
                    sample,
                )
            producer_environment_ok = (
                isinstance(producer_environment["python_version"], str)
                and bool(producer_environment["python_version"])
                and isinstance(producer_environment["platform"], str)
                and bool(producer_environment["platform"])
                and isinstance(producer_environment["normalizer_argv"], list)
                and bool(producer_environment["normalizer_argv"])
                and all(isinstance(token, str) and token for token in producer_environment["normalizer_argv"])
            )
        else:
            producer_environment_ok = False
        transition_summary = receipt.get("filter_transition_summary")
        transition_ok = False
        audit_chain_ok = False
        if isinstance(transition_summary, dict) and isinstance(receipt_outputs, dict):
            transition_counts = transition_summary.get("transition_counts")
            transition_ok = (
                transition_summary.get("pass") is True
                and isinstance(transition_summary.get("input_record_count"), int)
                and transition_summary.get("input_record_count") > 0
                and transition_summary.get("output_record_count")
                == transition_summary.get("input_record_count")
                and isinstance(transition_counts, dict)
                and all(isinstance(value, int) and value >= 0 for value in transition_counts.values())
                and sum(transition_counts.values()) == transition_summary.get("input_record_count")
            )
            normalization_doc = require_mapping(
                load_json_strict(Path(receipt_outputs["normalization_audit"]["path"])),
                f"{sample}.normalization_audit",
            )
            transition_doc = require_mapping(
                load_json_strict(Path(receipt_outputs["filter_transition_audit"]["path"])),
                f"{sample}.filter_transition_audit",
            )
            sample_doc = require_mapping(
                load_json_strict(Path(receipt_outputs["sample_verification"]["path"])),
                f"{sample}.sample_verification",
            )
            expected_transition = {
                "input_record_count": transition_doc.get("input", {}).get("record_count"),
                "output_record_count": transition_doc.get("output", {}).get("record_count"),
                "rescued_nonpass_to_pass": transition_doc.get("rescued_nonpass_to_pass"),
                "removed_pass_to_nonpass": transition_doc.get("removed_pass_to_nonpass"),
                "transition_counts": transition_doc.get("filter_transition_counts"),
                "pass": transition_doc.get("pass"),
            }
            audit_chain_ok = (
                transition_summary == expected_transition
                and normalization_doc.get("pass") is True
                and normalization_doc.get("input", {}).get("record_count")
                == transition_summary.get("input_record_count")
                and normalization_doc.get("output", {}).get("record_count")
                == transition_summary.get("input_record_count")
                and normalization_doc.get("rescued_nonpass_to_pass") == 0
                and normalization_doc.get("removed_pass_to_nonpass") == 0
                and sample_doc.get("sample") == sample
                and sample_doc.get("pass") is True
                and sample_doc.get("normalization") == normalization_doc
                and sample_doc.get("filter_transitions") == transition_doc
            )
        bam_policy = receipt.get("bam_output_policy")
        bam_policy_ok = False
        if isinstance(bam_policy, dict):
            fifo = Path(str(bam_policy.get("consumed_fifo_path", "")))
            try:
                fifo_is_named_pipe = stat.S_ISFIFO(fifo.stat().st_mode)
                regular_bams = [path for path in fifo.parent.glob("*.bam") if path.is_file()]
            except OSError:
                fifo_is_named_pipe = False
                regular_bams = [fifo]
            bam_policy_ok = (
                bam_policy.get("transport") == "named_fifo"
                and bam_policy.get("persisted_bam") is False
                and bam_policy.get("regular_bam_count") == 0
                and bam_policy.get("is_fifo_at_closeout") is True
                and fifo_is_named_pipe
                and not regular_bams
            )
        patch_ok = False
        if isinstance(patch_evidence, dict) and isinstance(longphase_executable, dict):
            patch_document = require_mapping(
                load_json_strict(Path(patch_evidence["patch_receipt"]["path"])),
                f"{sample}.patch_receipt",
            )
            patch_ok = (
                patch_evidence.get("approval_scope") == "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY"
                and patch_document.get("status") == "APPROVED_FOR_FAIL_CLOSED_7_DATASET_VALIDATION"
                and patch_document.get("approval", {}).get("scope")
                == patch_evidence.get("approval_scope")
                and patch_document.get("patched_binary_sha256")
                == longphase_executable.get("identity", {}).get("sha256")
            )
        receipt_ok = (
            receipt.get("schema_name") == "intersubmod.longphase_raw_all_capture_receipt"
            and receipt.get("schema_version") == "2.0.0"
            and receipt.get("sample") == sample
            and receipt.get("evidence_origin")
            == "post_run_normalization_from_frozen_raw_all_execution_artifacts"
            and receipt.get("identity_schema") == "coordinate_join_v1"
            and receipt.get("assurance") == "bounded_coordinate_equivalence_not_global_payload_identity"
            and receipt.get("longphase_input_contract") == "normalized_ClairS_raw_all"
            and receipt.get("tree_input_contract") == "longphase_s_recalibrated_FILTER_PASS"
            and receipt.get("duplicate_identity_policy")
            == "collapse_redundant_rows_with_identical_HP_PS"
            and meta.get("duplicate_identity_policy")
            == receipt.get("duplicate_identity_policy")
            and receipt.get("production_policy") == policy == analysis.get("production_filter_policy")
            and argv_ok
            and not forbidden_present
            and receipt.get("command_argv_sha256") == canonical_sha256(argv)
            and isinstance(receipt_inputs, dict)
            and receipt.get("producer_input_binding_sha256") == canonical_sha256(receipt_inputs)
            and isinstance(receipt_counts, dict)
            and receipt_counts.get("mapped_alignment_count") == scan["mapped_alignment_count"]
            and receipt_counts.get("identity_unique_count") == scan["identity_unique_count"]
            and receipt_counts.get("duplicate_count") == scan["duplicate_count"]
            and receipt_counts.get("conflict_count") == 0 == scan["conflict_count"]
            and isinstance(receipt_outputs, dict)
            and receipt_outputs.get("sidecar") == meta["sidecar"]
            and receipt_outputs.get("sidecar_index") == meta["index"]
            and receipt_outputs.get("native_validation") == meta["validation"]
            and producer_environment_ok
            and transition_ok
            and audit_chain_ok
            and bam_policy_ok
            and patch_ok
        )
        # Exact artifact-object comparisons keep raw, normalized, baseline,
        # recalibrated-all, and recalibrated-PASS roles from exchanging places.
        role_bindings_ok = (
            isinstance(receipt_outputs, dict)
            and isinstance(receipt_inputs, dict)
            and receipt_outputs.get("longphase_recalibrated_all_vcf", {}).get("identity", {}).get("sha256")
            == somatic_hashes["longphase_recalibrated_all_vcf"]
            and receipt_outputs.get("longphase_recalibrated_pass_vcf", {}).get("identity", {}).get("sha256")
            == somatic_hashes["longphase_recalibrated_pass_vcf"]
            and receipt_inputs.get("caller_raw_vcf", {}).get("identity", {}).get("sha256")
            == somatic_hashes["caller_raw_vcf"]
            and receipt_inputs.get("longphase_input_vcf", {}).get("identity", {}).get("sha256")
            == somatic_hashes["longphase_input_vcf"]
            and receipt_inputs.get("caller_pass_baseline_vcf", {}).get("identity", {}).get("sha256")
            == somatic_hashes["caller_pass_baseline_vcf"]
            and receipt_inputs.get("tumor_bam", {}).get("storage_identity_v1", {}).get("identity_sha256")
            == alignment_identity_sha
        )
        evidence_ok = (
            evidence.get("command_argv_sha256") == receipt.get("command_argv_sha256")
            and evidence.get("input_binding_sha256") == receipt.get("producer_input_binding_sha256")
            and evidence.get("effective_options_sha256") == canonical_sha256(receipt.get("effective_options"))
            and evidence.get("producer_receipt_sha256") == producer_receipt_sha
            and evidence.get("mapped_alignment_count") == scan["mapped_alignment_count"]
            and evidence.get("identity_unique_count") == scan["identity_unique_count"]
            and evidence.get("duplicate_count") == scan["duplicate_count"]
            and evidence.get("conflict_count") == 0 == scan["conflict_count"]
            and evidence.get("producer_inputs") == receipt_inputs
            and evidence.get("longphase_executable") == receipt.get("longphase_executable")
            and evidence.get("producer_code") == receipt.get("producer_code")
            and evidence.get("environment_lock") == receipt.get("environment_lock")
            and evidence.get("patch_evidence") == receipt.get("patch_evidence")
            and evidence.get("filter_transition_summary") == transition_summary
            and evidence.get("bam_output_policy") == bam_policy
        )
        self.add(
            "sidecar_validation_contract",
            native_ok and receipt_ok and role_bindings_ok and evidence_ok,
            "E_SIDECAR_CONTRACT_UNSUPPORTED",
            {
                "native_ok": native_ok,
                "receipt_ok": receipt_ok,
                "role_bindings_ok": role_bindings_ok,
                "evidence_ok": evidence_ok,
                "forbidden_present": forbidden_present,
            },
            "native PASS + frozen receipt; redundant identical tags collapse; no conflicts/truth/disableFilter",
            sample,
        )
        return {
            "identity_schema": "coordinate_join_v1",
            "sidecar_sha256": sidecar_sha,
            "sidecar_index_sha256": index_sha,
            "alignment_payload_identity_sha256": alignment_identity_sha,
            "duplicate_count": scan["duplicate_count"],
            "duplicate_identity_policy": meta["duplicate_identity_policy"],
        }

    def _verify_sample(self, meta: Mapping[str, Any], analysis: Mapping[str, Any]) -> dict[str, Any]:
        sample = meta["sample"]
        _, alignment_identity = self._verify_storage_identity(meta["alignment_payload"], sample)
        somatic = require_mapping(meta["somatic"], f"{sample}.somatic")
        require_exact_keys(
            somatic,
            {
                "caller_raw_vcf", "longphase_input_vcf", "caller_pass_baseline_vcf",
                "longphase_recalibrated_all_vcf", "longphase_recalibrated_pass_vcf", "tree_vcf",
            },
            f"{sample}.somatic",
        )
        somatic_hashes: dict[str, str] = {}
        for role, spec in somatic.items():
            self._verify_indexed_artifact(spec, f"{sample}.somatic.{role}", sample)
            somatic_hashes[role] = spec["identity"]["sha256"]
        selected_role = (
            "longphase_recalibrated_pass_vcf"
            if self.tree_contract == CANONICAL_TREE_INPUT_CONTRACT
            else "caller_pass_baseline_vcf"
        )
        selected_ok = (
            somatic["tree_vcf"] == somatic[selected_role]
            and somatic["longphase_recalibrated_pass_vcf"] != somatic["longphase_recalibrated_all_vcf"]
        )
        self.add(
            "selected_tree_role",
            selected_ok,
            "E_TREE_VCF_ROLE",
            {
                "tree": somatic["tree_vcf"].get("path"),
                "selected_role": selected_role,
                "selected": somatic[selected_role].get("path"),
            },
            self.tree_contract,
            sample,
        )
        tag_binding = self._verify_read_tags(
            require_mapping(meta["read_tags"], f"{sample}.read_tags"),
            analysis,
            alignment_identity,
            somatic_hashes,
            sample,
        )
        copy_number = require_mapping(meta["copy_number"], f"{sample}.copy_number")
        cn_fields = {
            "availability",
            "source",
            "semantics",
            "coordinate_system",
            "unlisted_position_semantics",
            "allowed_states",
            "overlap_policy",
            "reason",
            "cn_bed",
            "cn_int_gain",
            "cn_int_loss",
            "integration_json",
        }
        require_exact_keys(copy_number, cn_fields, f"{sample}.copy_number")
        if copy_number.get("availability") == "measured":
            measured_ok = (
                isinstance(copy_number.get("source"), str)
                and bool(copy_number["source"])
                and isinstance(copy_number.get("semantics"), str)
                and bool(copy_number["semantics"])
                and copy_number.get("coordinate_system") == "0_based_half_open"
                and copy_number.get("unlisted_position_semantics") == "neutral"
                and copy_number.get("allowed_states") == ["gain", "loss", "loh", "neutral"]
                and copy_number.get("overlap_policy") == "forbid"
                and copy_number.get("reason") is None
                and copy_number.get("cn_bed") is not None
            )
            self.add(
                "copy_number_measured_contract",
                measured_ok,
                "E_CN_CONTRACT",
                copy_number,
                "reviewed measured source, explicit 0-based neutral-unlisted semantics, cn_bed",
                sample,
            )
            for role in ("cn_bed", "cn_int_gain", "cn_int_loss", "integration_json"):
                spec = copy_number.get(role)
                if spec is not None:
                    self._verify_full_identity(spec, f"{sample}.copy_number.{role}", sample)
        else:
            unavailable = {
                "availability": "unavailable",
                "source": "unavailable",
                "semantics": "missing; never interpreted neutral",
                "coordinate_system": None,
                "unlisted_position_semantics": "unavailable",
                "allowed_states": [],
                "overlap_policy": "not_applicable",
                "reason": "No reviewed measured CN source is available",
                "cn_bed": None,
                "cn_int_gain": None,
                "cn_int_loss": None,
                "integration_json": None,
            }
            self.add(
                "copy_number_unavailable_contract",
                copy_number == unavailable,
                "E_CN_CONTRACT",
                copy_number,
                unavailable,
                sample,
            )

        input_sha = input_set_digest(meta)
        provenance = {
            "frozen_lock_sha256": self.lock_file_sha256,
            "launch_receipt_sha256": self.receipt_sha256,
            "environment_lock_sha256": self.environment_sha256,
            "source_bundle_manifest_sha256": self.bundle_manifest_sha256,
            "source_bundle_content_sha256": self.bundle_content_sha256,
            "input_set_sha256": input_sha,
        }
        somatic_roles = {
            "longphase_input_role": "normalized_clairs_raw_all",
            "longphase_input_vcf_sha256": somatic_hashes["longphase_input_vcf"],
            "caller_pass_baseline_role": (
                "clairs_filter_pass_selected_sensitivity_tree"
                if self.tree_contract == SENSITIVITY_TREE_INPUT_CONTRACT
                else "clairs_filter_pass_sensitivity_only"
            ),
            "caller_pass_baseline_vcf_sha256": somatic_hashes["caller_pass_baseline_vcf"],
            "tree_input_contract": self.tree_contract,
            "tree_backbone_role": self.tree_profile["backbone_role"],
            "tree_vcf_sha256": somatic_hashes["tree_vcf"],
            "ledger_role": "clairs_raw",
            "caller_raw_vcf_sha256": somatic_hashes["caller_raw_vcf"],
            "longphase_recalibrated_all_vcf_sha256": somatic_hashes["longphase_recalibrated_all_vcf"],
            "longphase_recalibrated_pass_vcf_sha256": somatic_hashes[
                "longphase_recalibrated_pass_vcf"
            ],
        }
        output_paths = self._verify_output_manifest(
            meta, provenance, somatic_roles, copy_number
        )
        if set(output_paths) != set(SCIENTIFIC_ROLES):
            return {"input_set_sha256": input_sha}
        _, funnel, group_keys = self._verify_parts(
            sample, output_paths, provenance, tag_binding, somatic_roles
        )
        layered = self._verify_layered(sample, output_paths["layered_reconstruction"], provenance, funnel)
        region = self._verify_region(
            sample,
            output_paths["layered_region_view"],
            provenance,
            funnel,
            group_keys,
            layered["l1"],
            copy_number,
        )
        ledger = self._verify_site_ledger(
            sample,
            output_paths["site_ledger"],
            output_paths["site_ledger_summary"],
            provenance,
            self.tree_profile["ledger_tree_contract"],
        )
        return {
            "input_set_sha256": input_sha,
            "read_tag_identity_schema": tag_binding["identity_schema"],
            "n_sSNV_scope_input": funnel["n_sSNV_scope_input"],
            "n_groups": len(group_keys),
            "n_detail_units": layered["n_detail_units"],
            "n_regions": region["n_regions"],
            "site_ledger_rows": ledger["rows"],
        }

    def _verify_output_manifest(
        self,
        meta: Mapping[str, Any],
        provenance: Mapping[str, str],
        somatic_roles: Mapping[str, str],
        copy_number: Mapping[str, Any],
    ) -> dict[str, Path]:
        sample = meta["sample"]
        sample_root = self.run_root / "samples" / sample
        manifest_path = regular_file(sample_root / "output_manifest.json", f"{sample} output manifest")
        manifest = require_mapping(load_json_strict(manifest_path), f"{sample} output manifest")
        require_exact_keys(
            manifest,
            {
                "schema_name", "schema_version", "sample", "biological_id", "run_id",
                *PROVENANCE_FIELDS, "somatic_roles", "copy_number_contract", "outputs",
            },
            f"{sample} output manifest",
        )
        expected_bindings = {
            "schema_name": OUTPUT_MANIFEST_SCHEMA_NAME,
            "schema_version": OUTPUT_MANIFEST_SCHEMA_VERSION,
            "sample": sample,
            "biological_id": meta["biological_id"],
            "run_id": self.run_id,
            "somatic_roles": dict(somatic_roles),
            "copy_number_contract": dict(copy_number),
            **provenance,
        }
        self.add(
            "output_manifest_bindings",
            all(manifest.get(key) == value for key, value in expected_bindings.items()),
            "E_SAMPLE_OUTPUT_BINDING",
            {key: manifest.get(key) for key in expected_bindings},
            expected_bindings,
            sample,
        )
        outputs = require_list(manifest["outputs"], f"{sample}.outputs")
        by_role: dict[str, Mapping[str, Any]] = {}
        for index, raw in enumerate(outputs):
            item = require_mapping(raw, f"{sample}.outputs[{index}]")
            require_exact_keys(item, {"role", "path", "size_bytes", "sha256"}, f"{sample}.outputs[{index}]")
            if item["role"] in by_role:
                raise ContractError("E_REQUIRED_OUTPUT", f"duplicate output role: {item['role']}")
            by_role[item["role"]] = item
        self.add(
            "output_role_set",
            set(by_role) == set(SCIENTIFIC_ROLES),
            "E_REQUIRED_OUTPUT",
            sorted(by_role),
            sorted(SCIENTIFIC_ROLES),
            sample,
        )
        paths: dict[str, Path] = {}
        for role in SCIENTIFIC_ROLES:
            if role not in by_role:
                continue
            item = by_role[role]
            expected_relative = output_name(sample, role)
            path = contained_file(item["path"], sample_root, f"{sample}:{role}")
            actual_size, actual_sha = path.stat().st_size, self.hash(path)
            expected_size = require_nonnegative_int(item["size_bytes"], f"{sample}:{role}.size")
            expected_sha = require_sha256(item["sha256"], f"{sample}:{role}.sha")
            self.add(
                f"output_identity:{role}",
                item["path"] == expected_relative
                and actual_size == expected_size
                and actual_sha == expected_sha,
                "E_OUTPUT_HASH_MISMATCH",
                {"path": item["path"], "size_bytes": actual_size, "sha256": actual_sha},
                {"path": expected_relative, "size_bytes": expected_size, "sha256": expected_sha},
                sample,
            )
            paths[role] = path
        return paths

    def _verify_provenance(
        self, document: Mapping[str, Any], expected: Mapping[str, str], name: str, sample: str
    ) -> None:
        observed = document.get("provenance")
        self.add(
            f"provenance:{name}",
            isinstance(observed, dict) and set(observed) == PROVENANCE_FIELDS and observed == expected,
            "E_SAMPLE_OUTPUT_BINDING",
            observed,
            expected,
            sample,
        )

    def _parse_funnel(self, value: Any, where: str) -> dict[str, int]:
        funnel = require_mapping(value, where)
        parsed = {field: require_nonnegative_int(funnel.get(field), f"{where}.{field}") for field in FUNNEL_FIELDS}
        if funnel.get("check_scope_conservation") is not True:
            raise ContractError("E_PART_CONSERVATION", f"{where}.check_scope_conservation is not true")
        return parsed

    @staticmethod
    def _funnel_conserved(funnel: Mapping[str, int]) -> bool:
        return (
            funnel["n_sSNV_scope_input"]
            == funnel["n_positional_singleton"] + funnel["n_multilocus_pre_cap_sSNV"]
            == funnel["n_sSNV_accounted"]
            and funnel["n_multilocus_pre_cap_sSNV"]
            == funnel["n_sSNV_retained"]
            + funnel["n_sSNV_read_unsupported"]
            + funnel["n_sSNV_cap_excluded"]
            and funnel["n_sSNV_accounted"]
            == funnel["n_positional_singleton"]
            + funnel["n_sSNV_cap_excluded"]
            + funnel["n_sSNV_read_unsupported"]
            + funnel["n_sSNV_retained"]
            and funnel["n_multilocus_pre_cap_groups"]
            == funnel["n_groups_retained"] + funnel["n_groups_read_unsupported"]
            and funnel["n_groups_capped_by_MAX_SNV"] <= funnel["n_groups_retained"]
            + funnel["n_groups_read_unsupported"]
        )

    def _verify_parts(
        self,
        sample: str,
        paths: Mapping[str, Path],
        provenance: Mapping[str, str],
        tag: Mapping[str, str],
        somatic_roles: Mapping[str, str],
    ) -> tuple[list[dict[str, Any]], dict[str, int], set[tuple[str, int, int]]]:
        parts = []
        aggregate = {field: 0 for field in FUNNEL_FIELDS}
        group_keys: set[tuple[str, int, int]] = set()
        duplicates = []
        for part, role in enumerate(PART_ROLES, start=1):
            document = require_mapping(load_json_strict(paths[role]), f"{sample}:{role}")
            self._verify_provenance(document, provenance, role, sample)
            identity_ok = (
                document.get("schema_version") == "2.0"
                and document.get("sample") == sample
                and document.get("part") == part
                and document.get("chromosomes") == list(PART_CHROMOSOMES[part])
                and document.get("somatic_roles") == dict(somatic_roles)
            )
            self.add(
                f"part_identity:{part}",
                identity_ok,
                "E_PART_CONSERVATION",
                [
                    document.get("schema_version"), document.get("sample"), document.get("part"),
                    document.get("chromosomes"), document.get("somatic_roles"),
                ],
                ["2.0", sample, part, list(PART_CHROMOSOMES[part]), dict(somatic_roles)],
                sample,
            )
            funnel = self._parse_funnel(document.get("input_funnel"), f"{sample}:{role}.funnel")
            self.add(
                f"part_funnel:{part}",
                self._funnel_conserved(funnel),
                "E_PART_CONSERVATION",
                funnel,
                "lossless funnel",
                sample,
            )
            for field in FUNNEL_FIELDS:
                aggregate[field] += funnel[field]
            groups = require_list(document.get("groups"), f"{sample}:{role}.groups")
            self.add(
                f"part_group_count:{part}",
                len(groups) == funnel["n_groups_retained"],
                "E_PART_CONSERVATION",
                len(groups),
                funnel["n_groups_retained"],
                sample,
            )
            for raw in groups:
                group = require_mapping(raw, "group")
                chrom, start, end = group.get("chrom"), group.get("start"), group.get("end")
                if (
                    chrom not in PART_CHROMOSOMES[part]
                    or isinstance(start, bool) or not isinstance(start, int)
                    or isinstance(end, bool) or not isinstance(end, int)
                    or start < 1 or end < start
                ):
                    self.add("part_group_coordinates", False, "E_PART_CONSERVATION", group, "valid", sample)
                    continue
                key = (chrom, start, end)
                if key in group_keys:
                    duplicates.append(key)
                group_keys.add(key)
            census = require_mapping(document.get("read_tag_census"), f"{sample}:{role}.read_tag_census")
            zero_fields = (
                "sidecar_missing", "sidecar_conflicts", "sidecar_extra", "sidecar_malformed",
                "alignment_identity_allele_conflicts",
            )
            values = {
                field: require_nonnegative_int(census.get(field), f"{sample}:{role}.{field}")
                for field in zero_fields
                + ("sidecar_duplicates", "sidecar_exact_matches", "alignment_group_exposures")
            }
            tag_ok = (
                census.get("identity_schema") == tag["identity_schema"]
                and census.get("sidecar_sha256") == tag["sidecar_sha256"]
                and census.get("sidecar_index_sha256") == tag["sidecar_index_sha256"]
                and census.get("alignment_payload_identity_sha256") == tag["alignment_payload_identity_sha256"]
                and census.get("duplicate_identity_policy") == tag["duplicate_identity_policy"]
                and values["sidecar_duplicates"] == tag["duplicate_count"]
                and all(values[field] == 0 for field in zero_fields)
                and values["sidecar_exact_matches"] == values["alignment_group_exposures"]
            )
            self.add(
                f"part_tag_binding:{part}",
                tag_ok,
                "E_PART_CONSERVATION",
                {**values, "identity_schema": census.get("identity_schema")},
                "frozen tag digests; expected redundant rows; zero loss/conflict; match=exposure",
                sample,
            )
            parts.append(document)
        self.add("five_parts_exact", len(parts) == 5, "E_PART_CONSERVATION", len(parts), 5, sample)
        self.add("part_groups_unique", not duplicates, "E_PART_CONSERVATION", duplicates, [], sample)
        return parts, aggregate, group_keys

    def _verify_layered(
        self,
        sample: str,
        path: Path,
        provenance: Mapping[str, str],
        aggregate: Mapping[str, int],
    ) -> dict[str, Any]:
        document = require_mapping(load_json_strict(path), f"{sample}:layered")
        self._verify_provenance(document, provenance, "layered", sample)
        self.add(
            "layered_identity",
            document.get("schema_version") == "2.0" and document.get("sample") == sample,
            "E_LAYERED_INVARIANT",
            [document.get("schema_version"), document.get("sample")],
            ["2.0", sample],
            sample,
        )
        funnel = self._parse_funnel(document.get("input_funnel"), f"{sample}:layered.funnel")
        self.add(
            "layered_part_conservation",
            funnel == dict(aggregate),
            "E_LAYERED_INVARIANT",
            funnel,
            aggregate,
            sample,
        )
        detail = require_list(document.get("detail"), f"{sample}:layered.detail")
        l1 = require_mapping(document.get("L1_ssnv_algorithm"), f"{sample}:layered.L1")
        primary = reference = invalid_primary = skipped = failures = incomplete = tree_mismatch = 0
        missing_shape = negative_hidden = duplicate_units = 0
        keys: set[tuple[Any, ...]] = set()
        for raw in detail:
            unit = require_mapping(raw, "layered unit")
            is_primary = unit.get("is_primary_lineage") is True
            is_reference = unit.get("reference_only") is True
            primary += int(is_primary); reference += int(is_reference)
            invalid_primary += int(is_primary and (str(unit.get("family")) not in {"1", "2"} or is_reference))
            capped = unit.get("capped") is True
            skipped += int(not capped and unit.get("verification_skipped") is True)
            failures += int(unit.get("verification_status") == "fail")
            if not capped:
                incomplete += int(unit.get("analysis_candidate_set_complete") is not True)
                tree_mismatch += int(unit.get("analysis_trees_generated") != unit.get("n_trees"))
                missing_shape += int(unit.get("n_distinct_shapes_exact") is None)
            for tree in unit.get("trees", []):
                hidden = tree.get("n_hidden") if isinstance(tree, dict) else None
                negative_hidden += int(isinstance(hidden, bool) or not isinstance(hidden, int) or hidden < 0)
            key = (unit.get("chrom"), unit.get("start"), unit.get("end"), unit.get("family"), unit.get("unit_role"))
            duplicate_units += int(key in keys); keys.add(key)
        observed = {
            "detail_count": len(detail), "primary": primary, "reference": reference,
            "invalid_primary": invalid_primary, "skipped": skipped, "failures": failures,
            "incomplete": incomplete, "tree_mismatch": tree_mismatch, "missing_shape": missing_shape,
            "negative_hidden": negative_hidden, "duplicate_units": duplicate_units,
        }
        ok = (
            document.get("n_detail_units") == len(detail)
            and all(observed[key] == 0 for key in (
                "invalid_primary", "skipped", "failures", "incomplete", "tree_mismatch",
                "missing_shape", "negative_hidden", "duplicate_units"
            ))
            and l1.get("all_eligible_V1V7_pass") is True
            and l1.get("n_verification_fail") == 0
            and l1.get("n_eligible_skipped_V4V5") == 0
            and l1.get("n_primary_lineage_units") == primary
            and l1.get("n_reference_only_controls") == reference
            and l1.get("n_units_total_including_unphased") == len(detail)
        )
        self.add("layered_scientific_invariants", ok, "E_LAYERED_INVARIANT", {**observed, "L1": l1}, "recomputed", sample)
        return {"n_detail_units": len(detail), "l1": l1}

    def _verify_region(
        self,
        sample: str,
        path: Path,
        provenance: Mapping[str, str],
        aggregate: Mapping[str, int],
        group_keys: set[tuple[str, int, int]],
        layered_l1: Mapping[str, Any],
        copy_number: Mapping[str, Any],
    ) -> dict[str, Any]:
        document = require_mapping(load_json_strict(path), f"{sample}:region")
        self._verify_provenance(document, provenance, "region", sample)
        census = require_mapping(document.get("census"), f"{sample}:region.census")
        funnel = self._parse_funnel(census.get("funnel"), f"{sample}:region.funnel")
        regions = require_list(document.get("regions"), f"{sample}:regions")
        keys: set[tuple[str, int, int]] = set()
        duplicates = multiplicity = invalid_primary = 0
        for raw in regions:
            region = require_mapping(raw, "region")
            key = (region.get("chrom"), region.get("start"), region.get("end"))
            duplicates += int(key in keys); keys.add(key)
            families = set()
            for lineage_raw in require_list(region.get("lineages"), "region.lineages"):
                lineage = require_mapping(lineage_raw, "lineage")
                if lineage.get("is_primary_lineage") is True:
                    family = str(lineage.get("family")); families.add(family)
                    invalid_primary += int(family not in {"1", "2"} or lineage.get("reference_only") is True)
            multiplicity += int(region.get("hp_multiplicity") != len(families))
        shared = ("n_units_total_including_unphased", "n_primary_lineage_units", "n_reference_only_controls", "n_verification_fail")
        region_l1 = census.get("L1")
        l1_ok = isinstance(region_l1, dict) and all(region_l1.get(key) == layered_l1.get(key) for key in shared)
        ok = (
            document.get("schema_version") == "2.0" and document.get("sample") == sample
            and document.get("copy_number_contract") == dict(copy_number)
            and funnel == dict(aggregate)
            and census.get("n_regions") == len(regions) == len(group_keys)
            and keys == group_keys and duplicates == 0 and multiplicity == 0 and invalid_primary == 0
            and funnel["n_groups_retained"] == len(regions) and l1_ok
        )
        self.add(
            "region_scientific_invariants",
            ok,
            "E_REGION_INVARIANT",
            {
                "n_regions": len(regions), "group_match": keys == group_keys, "duplicates": duplicates,
                "multiplicity_mismatch": multiplicity, "invalid_primary": invalid_primary, "l1_match": l1_ok,
            },
            "part groups/funnel and layered L1 recompute",
            sample,
        )
        return {"n_regions": len(regions)}

    def _verify_site_ledger(
        self,
        sample: str,
        ledger_path: Path,
        summary_path: Path,
        provenance: Mapping[str, str],
        expected_tree_contract: str,
    ) -> dict[str, Any]:
        summary = require_mapping(load_json_strict(summary_path), f"{sample}:ledger summary")
        self._verify_provenance(summary, provenance, "site_ledger_summary", sample)
        checks = summary.get("checks")
        branches = summary.get("branch_counts")
        transitions = summary.get("filter_transition_counts")
        duplicate_counts = summary.get("duplicate_record_key_excess")
        raw_count = summary.get("raw_clairs_records")
        shape_ok = (
            summary.get("schema_version") == "2.0" and summary.get("sample") == sample
            and summary.get("longphase_input_contract") == "clairs_raw_all"
            and summary.get("tree_contract") == expected_tree_contract
            and summary.get("pass") is True and isinstance(checks, dict) and bool(checks)
            and all(value is True for value in checks.values())
            and isinstance(branches, dict) and isinstance(raw_count, int) and not isinstance(raw_count, bool)
            and raw_count >= 0
            and summary.get("longphase_input_records") == raw_count
            and summary.get("longphase_recalibrated_records") == raw_count
            and isinstance(transitions, dict)
            and all(isinstance(value, int) and value >= 0 for value in transitions.values())
            and sum(transitions.values()) == raw_count
            and isinstance(duplicate_counts, dict)
            and set(duplicate_counts) == {"raw_clairs", "longphase_input", "longphase_recalibrated", "tree_input"}
            and all(value == 0 for value in duplicate_counts.values())
        )
        parsed = {}
        if isinstance(branches, dict):
            parsed = {key: require_nonnegative_int(value, f"branch_counts.{key}") for key, value in branches.items()}
        rows = wrong_sample = 0
        duplicate_columns = False
        try:
            with gzip.open(ledger_path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields = reader.fieldnames or []
                duplicate_columns = len(fields) != len(set(fields))
                if "sample" not in fields:
                    raise ContractError("E_SITE_LEDGER_INVARIANT", "ledger has no sample column")
                for row in reader:
                    rows += 1
                    wrong_sample += int(row.get("sample") != sample)
        except (OSError, UnicodeError, csv.Error) as exc:
            raise ContractError("E_SITE_LEDGER_INVARIANT", f"cannot stream ledger: {exc}") from exc
        ok = shape_ok and sum(parsed.values()) == raw_count == rows and wrong_sample == 0 and not duplicate_columns
        self.add(
            "site_ledger_scientific_invariants",
            ok,
            "E_SITE_LEDGER_INVARIANT",
            {
                "shape_ok": shape_ok, "branch_sum": sum(parsed.values()), "raw_count": raw_count,
                "rows": rows, "wrong_sample": wrong_sample, "duplicate_columns": duplicate_columns,
            },
            "branch sum=raw count=streamed rows; exact sample",
            sample,
        )
        return {"rows": rows}

    def _summary(self, lock: Mapping[str, Any]) -> dict[str, Any]:
        failures = [check for check in self.checks if not check.passed]
        ordered = sorted(self.samples, key=lambda item: item["sample"])
        all_pass = not failures and len(ordered) == 7 and all(item["pass"] for item in ordered)
        return {
            "schema_name": SUMMARY_SCHEMA_NAME,
            "schema_version": SUMMARY_SCHEMA_VERSION,
            "verifier_version": VERIFIER_VERSION,
            "verified_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "run_id": self.run_id,
            "run_root": str(self.run_root),
            "frozen_lock_path": str(self.lock_path),
            "frozen_lock_sha256": self.lock_file_sha256,
            "launch_receipt_path": str(self.receipt_path),
            "launch_receipt_sha256": self.receipt_sha256,
            "environment_lock_sha256": self.environment_sha256,
            "source_bundle_manifest_sha256": self.bundle_manifest_sha256,
            "source_bundle_content_sha256": self.bundle_content_sha256,
            "dataset_count": len(lock.get("samples", [])),
            "biological_sample_count": len({item.get("biological_id") for item in lock.get("samples", []) if isinstance(item, dict)}),
            "all_pass": all_pass,
            "n_pass": sum(item["pass"] for item in ordered),
            "n_fail": sum(not item["pass"] for item in ordered),
            "error_codes": sorted({check.code for check in failures}),
            "checks": [check.as_dict() for check in self.checks],
            "samples": ordered,
            "success_marker_created": False,
        }


def failure_summary(run_root: Path, lock: Path, receipt: Path, code: str, message: str) -> dict[str, Any]:
    return {
        "schema_name": SUMMARY_SCHEMA_NAME,
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "verifier_version": VERIFIER_VERSION,
        "verified_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "run_id": None,
        "run_root": str(run_root),
        "frozen_lock_path": str(lock),
        "frozen_lock_sha256": None,
        "launch_receipt_path": str(receipt),
        "launch_receipt_sha256": None,
        "environment_lock_sha256": None,
        "source_bundle_manifest_sha256": None,
        "source_bundle_content_sha256": None,
        "dataset_count": None,
        "biological_sample_count": None,
        "all_pass": False,
        "n_pass": 0,
        "n_fail": None,
        "error_codes": [code],
        "checks": [Check("fatal_contract", False, code, message, "valid lifecycle/lock contract").as_dict()],
        "samples": [],
        "success_marker_created": False,
    }


def render_tsv(summary: Mapping[str, Any]) -> bytes:
    rows = ["scope\tsample\tcheck\tpass\terror_code\tobserved_json\texpected_json"]
    for check in summary.get("checks", []):
        sample = check.get("sample") or "ALL"
        values = [
            "sample" if sample != "ALL" else "global",
            sample,
            str(check.get("name", "")),
            "PASS" if check.get("pass") else "FAIL",
            "" if check.get("pass") else str(check.get("code") or "E_UNKNOWN"),
            json.dumps(check.get("observed"), ensure_ascii=False, sort_keys=True, separators=(",", ":")),
            json.dumps(check.get("expected"), ensure_ascii=False, sort_keys=True, separators=(",", ":")),
        ]
        rows.append("\t".join(value.replace("\t", "\\t").replace("\n", "\\n") for value in values))
    return ("\n".join(rows) + "\n").encode("utf-8")


def publish(output: Path, summary: Mapping[str, Any]) -> tuple[Path, Path, Path]:
    json_path = output
    tsv_path = output.with_suffix(".tsv")
    checksum_path = output.with_suffix(".sha256")
    for path in (json_path, tsv_path, checksum_path):
        if path.exists() or path.is_symlink():
            raise ContractError("E_OUTPUT_EXISTS", f"verification output exists: {path}")
    json_payload = json.dumps(summary, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False).encode() + b"\n"
    tsv_payload = render_tsv(summary)
    json_sha = hashlib.sha256(json_payload).hexdigest()
    tsv_sha = hashlib.sha256(tsv_payload).hexdigest()
    checksum = f"{json_sha}  {json_path.name}\n{tsv_sha}  {tsv_path.name}\n".encode()
    atomic_write(json_path, json_payload)
    atomic_write(tsv_path, tsv_payload)
    if sha256_file(json_path) != json_sha or sha256_file(tsv_path) != tsv_sha:
        raise ContractError("E_OUTPUT_HASH_MISMATCH", "verifier artifact readback mismatch")
    atomic_write(checksum_path, checksum)
    return json_path, tsv_path, checksum_path


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--run-root", required=True, type=Path)
    result.add_argument("--frozen-lock", type=Path)
    result.add_argument("--launch-receipt", type=Path)
    result.add_argument("--output", type=Path)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    run_root = args.run_root
    lock = args.frozen_lock or run_root / RUN_LOCK_FILENAME
    receipt = args.launch_receipt or run_root / RUN_RECEIPT_FILENAME
    output = args.output or run_root / "verification_summary.json"
    exit_code = 0
    try:
        resolved_root = run_root.resolve(strict=True)
        output.parent.resolve(strict=True).relative_to(resolved_root)
        summary = Verifier(resolved_root, lock, receipt).verify()
        if not summary["all_pass"]:
            exit_code = 7
    except ContractError as exc:
        exit_code = exc.exit_code
        summary = failure_summary(run_root, lock, receipt, exc.code, str(exc))
    except (OSError, ValueError) as exc:
        exit_code = 2
        summary = failure_summary(run_root, lock, receipt, "E_BOOTSTRAP_PATH", str(exc))
    except BaseException as exc:  # pragma: no cover
        exit_code = 70
        summary = failure_summary(run_root, lock, receipt, "E_INTERNAL", repr(exc))
    try:
        json_path, tsv_path, checksum_path = publish(output, summary)
    except ContractError as exc:
        print(json.dumps({"pass": False, "code": exc.code, "message": str(exc)}), file=sys.stderr)
        return exc.exit_code
    except BaseException as exc:  # pragma: no cover
        print(json.dumps({"pass": False, "code": "E_INTERNAL", "message": repr(exc)}), file=sys.stderr)
        return 70
    print(
        f"VERIFY layered-v3: {'PASS' if summary['all_pass'] else 'FAIL'} "
        f"{summary.get('n_pass')}/{summary.get('dataset_count')} datasets -> {json_path}"
    )
    print(f"VERIFY artifacts: {tsv_path} {checksum_path}")
    if summary["error_codes"]:
        print(f"VERIFY error_codes: {','.join(summary['error_codes'])}", file=sys.stderr)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
