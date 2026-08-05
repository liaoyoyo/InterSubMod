#!/usr/bin/env python3
"""Fail-closed verifier for the HCC1395 Big7 PARTIAL pilot binding.

The verifier intentionally does not compute a full BAM hash. It binds the BAM
with exact filesystem metadata, file size, and fixed first/middle/last chunks;
the BAI is fully hashed. The frozen M2 manifest remains the authority for every
non-alignment field, including the VCF and HP/PS sidecar records.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
from typing import Any, Mapping, Sequence


TOPIC_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = TOPIC_ROOT / "input_contract" / "big7_hcc1395_pilot_manifest.json"
RECEIPT_SCHEMA_NAME = "intersubmod.big7_hcc1395_input_binding_receipt"
RECEIPT_SCHEMA_VERSION = "1.0.0"
PILOT_SCHEMA_NAME = "intersubmod.big7_hcc1395_partial_pilot_manifest"
PILOT_SCHEMA_VERSION = "1.0.0"
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))

PILOT_TOP_KEYS = {
    "schema_name",
    "schema_version",
    "manifest_id",
    "created_at_utc",
    "task_type",
    "validation_evidence_eligible",
    "scope",
    "source_authority",
    "extractor_contract",
    "dataset_count",
    "biological_sample_count",
    "samples",
}
SCOPE_KEYS = {
    "claim_status",
    "dataset",
    "biological_id",
    "contigs",
    "purpose",
    "claim_limit",
}
SOURCE_AUTHORITY_KEYS = {
    "manifest_path",
    "manifest_sha256",
    "sample",
    "sample_semantic_sha256",
    "alignment_transfer_assurance",
    "is_full_bam_content_hash",
}
EXTRACTOR_CONTRACT = {
    "loader": "json.loads_plus_input_entry",
    "required_sample_cardinality": 1,
    "canonical_v3_schema_reused": False,
    "reason": (
        "canonical layered_input_manifest v3 fixes seven samples; this HCC1395-only "
        "PARTIAL pilot uses a separate binding receipt contract"
    ),
}
ALIGNMENT_KEYS = {
    "path",
    "index_path",
    "embedded_hp_ps_policy",
    "identity_policy",
    "storage_identity_v1",
}
STORAGE_KEYS = {
    "policy",
    "assurance",
    "is_full_content_hash",
    "requested_path",
    "realpath",
    "logical_is_symlink",
    "logical_size_bytes",
    "logical_mtime_ns",
    "st_dev",
    "st_ino",
    "size_bytes",
    "mtime_ns",
    "ctime_ns",
    "chunk_size_bytes",
    "chunks",
    "index",
    "identity_sha256",
}


class VerificationError(RuntimeError):
    """A closed-contract verification failure."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class ChunkSpec:
    label: str
    offset: int
    length: int
    sha256: str


@dataclass(frozen=True)
class BindingContract:
    sample: str
    biological_id: str
    source_manifest_path: Path
    source_manifest_sha256: str
    source_sample_semantic_sha256: str
    bam_path: Path
    bai_path: Path
    logical_is_symlink: bool
    logical_size_bytes: int
    logical_mtime_ns: int
    st_dev: int
    st_ino: int
    bam_size_bytes: int
    mtime_ns: int
    ctime_ns: int
    chunk_size_bytes: int
    chunks: tuple[ChunkSpec, ...]
    bai_size_bytes: int
    bai_sha256: str
    storage_identity_sha256: str


PRODUCTION_CONTRACT = BindingContract(
    sample="HCC1395",
    biological_id="HCC1395",
    source_manifest_path=Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
        "20260716_read_linked_hypercube_exact_likelihood_validation/"
        "m2_frozen_release_v2/release_contract/input_contract/canonical_manifest.json"
    ),
    source_manifest_sha256="16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45",
    source_sample_semantic_sha256="29660d73461606f6cd9a1a38c5ee821cc1e70a37a80d01cc8136e2190f2e1e05",
    bam_path=Path("/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam"),
    bai_path=Path("/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam.bai"),
    logical_is_symlink=False,
    logical_size_bytes=292055926761,
    logical_mtime_ns=1775198635489054310,
    st_dev=2081,
    st_ino=613688725,
    bam_size_bytes=292055926761,
    mtime_ns=1775198635489054310,
    ctime_ns=1775198635489054310,
    chunk_size_bytes=1048576,
    chunks=(
        ChunkSpec(
            "first",
            0,
            1048576,
            "df999f92f396ceb5790584d8efc1480ea0d843a095a807b3e427ae8c4dda58de",
        ),
        ChunkSpec(
            "middle",
            146027439092,
            1048576,
            "6866c24ac09f498655aceb36d84514e1a468cace1e88cc88bb3224bb5ba992b5",
        ),
        ChunkSpec(
            "last",
            292054878185,
            1048576,
            "d72a059eff2888039414c07e14b0672756c1b95f54ba19ee11ce1010c545b54f",
        ),
    ),
    bai_size_bytes=119741784,
    bai_sha256="f7a23610899ff535de17ccb324f99fc8bede402f65b992d5601d34ffb62c49d3",
    storage_identity_sha256="dbb54185e37a374a89d3420980c1bf11e7d6bfdedd98d5768fb8abcc607ca4a2",
)


def canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def strict_json_load(path: Path) -> Any:
    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise VerificationError("E_JSON_DUPLICATE_KEY", f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise VerificationError("E_JSON_NONFINITE", f"non-finite JSON number: {value}")

    try:
        return json.loads(
            path.read_text(encoding="utf-8", errors="strict"),
            object_pairs_hook=pairs_hook,
            parse_constant=reject_constant,
        )
    except VerificationError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise VerificationError("E_JSON_READ", f"cannot read strict JSON {path}: {exc}") from exc


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(block_size), b""):
                digest.update(block)
    except OSError as exc:
        raise VerificationError("E_FILE_HASH", f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def expected_storage_identity(contract: BindingContract) -> dict[str, Any]:
    value: dict[str, Any] = {
        "policy": "storage_identity_v1",
        "assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
        "is_full_content_hash": False,
        "requested_path": str(contract.bam_path),
        "realpath": str(contract.bam_path),
        "logical_is_symlink": contract.logical_is_symlink,
        "logical_size_bytes": contract.logical_size_bytes,
        "logical_mtime_ns": contract.logical_mtime_ns,
        "st_dev": contract.st_dev,
        "st_ino": contract.st_ino,
        "size_bytes": contract.bam_size_bytes,
        "mtime_ns": contract.mtime_ns,
        "ctime_ns": contract.ctime_ns,
        "chunk_size_bytes": contract.chunk_size_bytes,
        "chunks": [
            {
                "label": chunk.label,
                "offset": chunk.offset,
                "length": chunk.length,
                "sha256": chunk.sha256,
            }
            for chunk in contract.chunks
        ],
        "index": {
            "path": str(contract.bai_path),
            "identity": {
                "policy": "full_sha256",
                "size_bytes": contract.bai_size_bytes,
                "sha256": contract.bai_sha256,
            },
        },
    }
    value["identity_sha256"] = canonical_sha256(value)
    if value["identity_sha256"] != contract.storage_identity_sha256:
        raise VerificationError(
            "E_INTERNAL_CONTRACT",
            "compiled storage identity digest does not match the binding contract",
        )
    return value


class CheckLog:
    def __init__(self) -> None:
        self.rows: list[dict[str, Any]] = []

    def require(
        self,
        condition: bool,
        name: str,
        code: str,
        message: str,
        *,
        expected: Any = None,
        observed: Any = None,
    ) -> None:
        row: dict[str, Any] = {"name": name, "pass": bool(condition)}
        if expected is not None:
            row["expected"] = expected
        if observed is not None:
            row["observed"] = observed
        self.rows.append(row)
        if not condition:
            raise VerificationError(code, message)


def require_exact_keys(value: Any, keys: set[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        raise VerificationError("E_SCHEMA", f"{label} exact-key schema mismatch")
    return value


def path_fields(value: Any, prefix: str = "$") -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    if isinstance(value, dict):
        for key, item in value.items():
            location = f"{prefix}.{key}"
            if key in {"path", "index_path", "requested_path", "realpath", "manifest_path"}:
                if isinstance(item, str):
                    rows.append((location, item))
                else:
                    rows.append((location, ""))
            rows.extend(path_fields(item, location))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            rows.extend(path_fields(item, f"{prefix}[{index}]"))
    return rows


def unique_sample(document: Mapping[str, Any], sample: str, label: str) -> dict[str, Any]:
    samples = document.get("samples")
    if not isinstance(samples, list):
        raise VerificationError("E_SCHEMA", f"{label}.samples is not a list")
    matches = [row for row in samples if isinstance(row, dict) and row.get("sample") == sample]
    if len(matches) != 1:
        raise VerificationError(
            "E_SAMPLE_CARDINALITY",
            f"{label} must contain exactly one {sample} entry; observed {len(matches)}",
        )
    return matches[0]


def file_stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        int(value.st_dev),
        int(value.st_ino),
        int(value.st_size),
        int(value.st_mtime_ns),
        int(value.st_ctime_ns),
    )


def open_regular_nonsymlink(path: Path, label: str) -> tuple[int, os.stat_result, os.stat_result]:
    if not path.is_absolute():
        raise VerificationError("E_PATH_ABSOLUTE", f"{label} path is not absolute: {path}")
    try:
        logical = path.lstat()
    except OSError as exc:
        raise VerificationError("E_FILE_MISSING", f"cannot stat {label} {path}: {exc}") from exc
    if stat.S_ISLNK(logical.st_mode) or not stat.S_ISREG(logical.st_mode):
        raise VerificationError("E_FILE_TYPE", f"{label} is not a non-symlink regular file: {path}")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
        target = os.fstat(descriptor)
    except OSError as exc:
        raise VerificationError("E_FILE_OPEN", f"cannot safely open {label} {path}: {exc}") from exc
    if not stat.S_ISREG(target.st_mode) or (target.st_dev, target.st_ino) != (logical.st_dev, logical.st_ino):
        os.close(descriptor)
        raise VerificationError("E_FILE_RACE", f"{label} identity changed during open: {path}")
    return descriptor, logical, target


def observe_bam(contract: BindingContract) -> dict[str, Any]:
    descriptor, logical, before = open_regular_nonsymlink(contract.bam_path, "BAM")
    try:
        chunks: list[dict[str, Any]] = []
        for chunk in contract.chunks:
            payload = os.pread(descriptor, chunk.length, chunk.offset)
            if len(payload) != chunk.length:
                raise VerificationError(
                    "E_BAM_CHUNK_SHORT",
                    f"short BAM chunk {chunk.label}: {len(payload)} != {chunk.length}",
                )
            chunks.append(
                {
                    "label": chunk.label,
                    "offset": chunk.offset,
                    "length": chunk.length,
                    "sha256": hashlib.sha256(payload).hexdigest(),
                }
            )
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if file_stat_identity(before) != file_stat_identity(after):
        raise VerificationError("E_FILE_RACE", "BAM changed while sampled chunks were read")
    return {
        "logical_is_symlink": stat.S_ISLNK(logical.st_mode),
        "logical_size_bytes": int(logical.st_size),
        "logical_mtime_ns": int(logical.st_mtime_ns),
        "st_dev": int(before.st_dev),
        "st_ino": int(before.st_ino),
        "size_bytes": int(before.st_size),
        "mtime_ns": int(before.st_mtime_ns),
        "ctime_ns": int(before.st_ctime_ns),
        "chunks": chunks,
    }


def observe_bai(contract: BindingContract) -> dict[str, Any]:
    descriptor, _logical, before = open_regular_nonsymlink(contract.bai_path, "BAI")
    digest = hashlib.sha256()
    first = b""
    try:
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            if not first:
                first = block[:4]
            digest.update(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if file_stat_identity(before) != file_stat_identity(after):
        raise VerificationError("E_FILE_RACE", "BAI changed while its full SHA-256 was computed")
    return {
        "size_bytes": int(before.st_size),
        "sha256": digest.hexdigest(),
        "magic_hex": first.hex(),
    }


def observe_full_sha256(path: Path, label: str) -> dict[str, Any]:
    descriptor, _logical, before = open_regular_nonsymlink(path, label)
    digest = hashlib.sha256()
    try:
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            digest.update(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if file_stat_identity(before) != file_stat_identity(after):
        raise VerificationError("E_FILE_RACE", f"{label} changed while its full SHA-256 was computed")
    return {"size_bytes": int(before.st_size), "sha256": digest.hexdigest()}


def resolve_executable(command: str) -> Path:
    candidate = shutil.which(command)
    if candidate is None:
        raise VerificationError("E_SAMTOOLS_MISSING", f"samtools executable not found: {command}")
    path = Path(candidate).resolve(strict=True)
    if not path.is_file() or not os.access(path, os.X_OK):
        raise VerificationError("E_SAMTOOLS_EXECUTABLE", f"samtools is not executable: {path}")
    return path


def run_quickcheck(samtools: str, bam: Path) -> dict[str, Any]:
    executable = resolve_executable(samtools)
    try:
        version = subprocess.run(
            [str(executable), "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        completed = subprocess.run(
            [str(executable), "quickcheck", "-v", str(bam)],
            check=False,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise VerificationError("E_SAMTOOLS_RUN", f"cannot run samtools quickcheck: {exc}") from exc
    return {
        "executable": str(executable),
        "version_line": (version.stdout.splitlines() or version.stderr.splitlines() or [""])[0],
        "command": [str(executable), "quickcheck", "-v", str(bam)],
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def verify_authority_artifact_metadata(
    artifact: Mapping[str, Any], label: str, checks: CheckLog
) -> dict[str, Any]:
    path = Path(str(artifact.get("path", "")))
    identity = artifact.get("identity")
    checks.require(path.is_absolute(), f"{label}_path_absolute", "E_PATH_ABSOLUTE", f"{label} path is not absolute")
    checks.require(
        isinstance(identity, dict)
        and identity.get("policy") == "full_sha256"
        and isinstance(identity.get("size_bytes"), int)
        and isinstance(identity.get("sha256"), str),
        f"{label}_identity_record_well_formed",
        "E_AUTHORITY_IDENTITY",
        f"{label} authority identity record is malformed",
    )
    observed = observe_full_sha256(path, label)
    observed_size = observed["size_bytes"]
    checks.require(
        observed_size == identity["size_bytes"],
        f"{label}_size_matches_frozen_authority",
        "E_AUTHORITY_SIZE",
        f"{label} size differs from frozen authority",
        expected=identity["size_bytes"],
        observed=observed_size,
    )
    observed_sha256 = observed["sha256"]
    checks.require(
        observed_sha256 == identity["sha256"],
        f"{label}_fresh_full_sha256_matches_frozen_authority",
        "E_AUTHORITY_SHA256",
        f"{label} full SHA-256 differs from frozen authority",
        expected=identity["sha256"],
        observed=observed_sha256,
    )
    return {
        "path": str(path),
        "size_bytes": observed_size,
        "sha256": observed_sha256,
        "verification_level": "fresh_full_sha256",
    }


def verify_binding(
    manifest_path: Path,
    *,
    contract: BindingContract = PRODUCTION_CONTRACT,
    samtools: str = "samtools",
) -> dict[str, Any]:
    checks = CheckLog()
    failure: dict[str, str] | None = None
    manifest_identity: dict[str, Any] | None = None
    source_identity: dict[str, Any] | None = None
    bam_observed: dict[str, Any] | None = None
    bai_observed: dict[str, Any] | None = None
    quickcheck: dict[str, Any] | None = None
    authority_artifacts: dict[str, Any] = {}

    try:
        checks.require(
            manifest_path.is_absolute(),
            "pilot_manifest_path_absolute",
            "E_PATH_ABSOLUTE",
            "pilot manifest path must be absolute",
            observed=str(manifest_path),
        )
        manifest = strict_json_load(manifest_path)
        checks.require(isinstance(manifest, dict), "pilot_manifest_is_object", "E_SCHEMA", "pilot manifest is not an object")
        require_exact_keys(manifest, PILOT_TOP_KEYS, "pilot manifest")
        checks.require(
            manifest.get("schema_name") == PILOT_SCHEMA_NAME
            and manifest.get("schema_version") == PILOT_SCHEMA_VERSION,
            "pilot_manifest_schema_exact",
            "E_SCHEMA",
            "pilot manifest schema name/version mismatch",
        )
        checks.require(
            manifest.get("task_type") == "exploratory_pilot"
            and manifest.get("validation_evidence_eligible") is False,
            "partial_pilot_claim_boundary",
            "E_CLAIM_BOUNDARY",
            "manifest must remain an ineligible exploratory PARTIAL pilot",
        )
        checks.require(
            manifest.get("dataset_count") == 1 and manifest.get("biological_sample_count") == 1,
            "hcc1395_only_counts",
            "E_SAMPLE_CARDINALITY",
            "pilot count fields must both equal one",
        )

        scope = require_exact_keys(manifest.get("scope"), SCOPE_KEYS, "scope")
        expected_scope = {
            "claim_status": "PARTIAL",
            "dataset": contract.sample,
            "biological_id": contract.biological_id,
            "contigs": list(AUTOSOMES),
            "purpose": "exact_ps_x_hp_k_le_12_input_binding",
            "claim_limit": (
                "HCC1395 chr1-22 input binding only; not seven-sample, production, "
                "or paper-final validation"
            ),
        }
        checks.require(
            dict(scope) == expected_scope,
            "scope_is_exact_hcc1395_autosome_partial",
            "E_SCOPE",
            "pilot scope or claim limit changed",
        )
        checks.require(
            manifest.get("extractor_contract") == EXTRACTOR_CONTRACT,
            "extractor_loader_contract_exact",
            "E_EXTRACTOR_CONTRACT",
            "extractor compatibility contract changed",
        )

        source_authority = require_exact_keys(
            manifest.get("source_authority"), SOURCE_AUTHORITY_KEYS, "source_authority"
        )
        expected_source_authority = {
            "manifest_path": str(contract.source_manifest_path),
            "manifest_sha256": contract.source_manifest_sha256,
            "sample": contract.sample,
            "sample_semantic_sha256": contract.source_sample_semantic_sha256,
            "alignment_transfer_assurance": (
                "same_size_plus_fixed_first_middle_last_chunk_sha256_plus_full_bai_sha256"
            ),
            "is_full_bam_content_hash": False,
        }
        checks.require(
            dict(source_authority) == expected_source_authority,
            "source_authority_contract_exact",
            "E_SOURCE_AUTHORITY",
            "source authority record changed",
        )

        absolute_fields = path_fields(manifest)
        non_absolute = [location for location, value in absolute_fields if not value or not Path(value).is_absolute()]
        checks.require(
            not non_absolute,
            "all_declared_input_paths_absolute",
            "E_PATH_ABSOLUTE",
            f"non-absolute declared input paths: {non_absolute}",
            observed=non_absolute,
        )

        manifest_identity = {
            "path": str(manifest_path),
            "size_bytes": manifest_path.stat().st_size,
            "sha256": sha256_path(manifest_path),
            "semantic_sha256": canonical_sha256(manifest),
        }
        source_sha = sha256_path(contract.source_manifest_path)
        source_document = strict_json_load(contract.source_manifest_path)
        source_identity = {
            "path": str(contract.source_manifest_path),
            "size_bytes": contract.source_manifest_path.stat().st_size,
            "sha256": source_sha,
            "semantic_sha256": canonical_sha256(source_document),
        }
        checks.require(
            source_sha == contract.source_manifest_sha256,
            "frozen_source_manifest_full_sha256",
            "E_SOURCE_SHA256",
            "frozen source manifest SHA-256 mismatch",
            expected=contract.source_manifest_sha256,
            observed=source_sha,
        )

        pilot_sample = unique_sample(manifest, contract.sample, "pilot manifest")
        checks.require(
            isinstance(manifest.get("samples"), list) and len(manifest["samples"]) == 1,
            "pilot_has_exactly_one_sample",
            "E_SAMPLE_CARDINALITY",
            "pilot must contain only the HCC1395 sample",
        )
        source_sample = unique_sample(source_document, contract.sample, "frozen source manifest")
        source_sample_sha = canonical_sha256(source_sample)
        checks.require(
            source_sample_sha == contract.source_sample_semantic_sha256,
            "frozen_hcc1395_sample_semantic_sha256",
            "E_SOURCE_SAMPLE_SHA256",
            "frozen HCC1395 sample semantic SHA-256 mismatch",
            expected=contract.source_sample_semantic_sha256,
            observed=source_sample_sha,
        )

        pilot_non_alignment = copy.deepcopy(pilot_sample)
        source_non_alignment = copy.deepcopy(source_sample)
        pilot_alignment = pilot_non_alignment.pop("alignment_payload", None)
        source_alignment = source_non_alignment.pop("alignment_payload", None)
        checks.require(
            pilot_non_alignment == source_non_alignment,
            "all_non_alignment_fields_equal_frozen_hcc1395_authority",
            "E_FROZEN_FIELDS",
            "VCF, sidecar, or another non-alignment HCC1395 field differs from frozen authority",
        )
        checks.require(
            isinstance(pilot_alignment, dict) and isinstance(source_alignment, dict),
            "alignment_records_present",
            "E_ALIGNMENT_SCHEMA",
            "pilot/source alignment payload is missing",
        )
        require_exact_keys(pilot_alignment, ALIGNMENT_KEYS, "pilot alignment_payload")
        expected_alignment = copy.deepcopy(source_alignment)
        expected_alignment["path"] = str(contract.bam_path)
        expected_alignment["index_path"] = str(contract.bai_path)
        expected_alignment["storage_identity_v1"] = expected_storage_identity(contract)
        checks.require(
            pilot_alignment == expected_alignment,
            "pilot_alignment_record_exact",
            "E_ALIGNMENT_BINDING",
            "pilot alignment payload differs from the compiled Big7 binding contract",
        )
        storage = require_exact_keys(
            pilot_alignment.get("storage_identity_v1"), STORAGE_KEYS, "storage_identity_v1"
        )
        checks.require(
            storage.get("is_full_content_hash") is False,
            "bam_not_misrepresented_as_full_content_hash",
            "E_ASSURANCE",
            "BAM identity must explicitly remain a sampled-chunk identity",
        )
        checks.require(
            storage["chunks"][1]["offset"] == contract.chunks[1].offset,
            "middle_chunk_offset_exactly_recorded",
            "E_BAM_MIDDLE_OFFSET",
            "middle BAM chunk offset differs from the fixed contract",
            expected=contract.chunks[1].offset,
            observed=storage["chunks"][1]["offset"],
        )

        authority_artifacts["tree_vcf"] = verify_authority_artifact_metadata(
            pilot_sample["somatic"]["tree_vcf"], "tree_vcf", checks
        )
        authority_artifacts["tree_vcf_index"] = verify_authority_artifact_metadata(
            pilot_sample["somatic"]["tree_vcf"]["index"], "tree_vcf_index", checks
        )
        authority_artifacts["read_tag_sidecar"] = verify_authority_artifact_metadata(
            pilot_sample["read_tags"]["sidecar"], "read_tag_sidecar", checks
        )
        authority_artifacts["read_tag_sidecar_index"] = verify_authority_artifact_metadata(
            pilot_sample["read_tags"]["index"], "read_tag_sidecar_index", checks
        )

        bam_observed = observe_bam(contract)
        expected_bam_observed = {
            "logical_is_symlink": contract.logical_is_symlink,
            "logical_size_bytes": contract.logical_size_bytes,
            "logical_mtime_ns": contract.logical_mtime_ns,
            "st_dev": contract.st_dev,
            "st_ino": contract.st_ino,
            "size_bytes": contract.bam_size_bytes,
            "mtime_ns": contract.mtime_ns,
            "ctime_ns": contract.ctime_ns,
            "chunks": expected_storage_identity(contract)["chunks"],
        }
        checks.require(
            bam_observed == expected_bam_observed,
            "bam_metadata_size_and_three_fixed_chunks_exact",
            "E_BAM_IDENTITY",
            "Big7 BAM metadata, size, or sampled chunk SHA-256 differs",
            expected=expected_bam_observed,
            observed=bam_observed,
        )

        bai_observed = observe_bai(contract)
        expected_bai_observed = {
            "size_bytes": contract.bai_size_bytes,
            "sha256": contract.bai_sha256,
            "magic_hex": "42414901",
        }
        checks.require(
            bai_observed == expected_bai_observed,
            "bai_size_magic_and_full_sha256_exact",
            "E_BAI_IDENTITY",
            "Big7 BAI size, magic, or full SHA-256 differs",
            expected=expected_bai_observed,
            observed=bai_observed,
        )

        quickcheck = run_quickcheck(samtools, contract.bam_path)
        checks.require(
            quickcheck["returncode"] == 0,
            "samtools_quickcheck_passes",
            "E_SAMTOOLS_QUICKCHECK",
            "samtools quickcheck failed for the bound Big7 BAM",
            expected=0,
            observed=quickcheck["returncode"],
        )
        post_quickcheck = observe_bam(contract)
        checks.require(
            post_quickcheck == bam_observed,
            "bam_unchanged_across_quickcheck",
            "E_FILE_RACE",
            "Big7 BAM changed across samtools quickcheck",
        )
    except VerificationError as exc:
        failure = {"code": exc.code, "message": str(exc)}
    except Exception as exc:  # pragma: no cover - defensive fail-closed boundary
        failure = {"code": "E_UNEXPECTED", "message": f"{type(exc).__name__}: {exc}"}

    all_pass = failure is None and bool(checks.rows) and all(row["pass"] for row in checks.rows)
    receipt: dict[str, Any] = {
        "schema_name": RECEIPT_SCHEMA_NAME,
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "task_type": "exploratory_pilot",
        "claim_status": "PARTIAL",
        "validation_evidence_eligible": False,
        "sample": contract.sample,
        "scope": {
            "contigs": list(AUTOSOMES),
            "claim_limit": "HCC1395 input binding only; no extraction and no cohort claim",
        },
        "manifest": manifest_identity,
        "source_authority": source_identity,
        "binding_contract": {
            "bam": {
                "path": str(contract.bam_path),
                "size_bytes": contract.bam_size_bytes,
                "identity_policy": "fixed_metadata_plus_three_sampled_chunks_not_full_bam_hash",
                "is_full_content_hash": False,
                "chunk_size_bytes": contract.chunk_size_bytes,
                "chunks": expected_storage_identity(contract)["chunks"] if failure is None or contract.chunks else [],
            },
            "bai": {
                "path": str(contract.bai_path),
                "size_bytes": contract.bai_size_bytes,
                "identity_policy": "full_sha256",
                "sha256": contract.bai_sha256,
            },
        },
        "observed": {
            "bam": bam_observed,
            "bai": bai_observed,
            "samtools_quickcheck": quickcheck,
            "frozen_authority_artifacts": authority_artifacts,
        },
        "checks": checks.rows,
        "all_pass": all_pass,
        "failure": failure,
    }
    integrity_payload = copy.deepcopy(receipt)
    receipt["receipt_integrity"] = {
        "policy": "semantic_json_sha256_without_receipt_integrity",
        "sha256": canonical_sha256(integrity_payload),
    }
    return receipt


def write_receipt(path: Path, receipt: Mapping[str, Any]) -> None:
    if not path.is_absolute():
        raise VerificationError("E_OUTPUT_PATH", "receipt output path must be absolute")
    if path.exists() or path.is_symlink():
        raise VerificationError("E_OUTPUT_EXISTS", f"refusing to overwrite receipt: {path}")
    if not path.parent.is_dir():
        raise VerificationError("E_OUTPUT_PARENT", f"receipt parent does not exist: {path.parent}")
    with path.open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    manifest_path = args.manifest.absolute()
    output_path = args.output.absolute()
    receipt = verify_binding(manifest_path, samtools=args.samtools)
    try:
        write_receipt(output_path, receipt)
    except VerificationError as exc:
        print(f"{exc.code}: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "all_pass": receipt["all_pass"],
                "failure": receipt["failure"],
                "receipt": str(output_path),
                "bam": str(PRODUCTION_CONTRACT.bam_path),
                "bai": str(PRODUCTION_CONTRACT.bai_path),
                "middle_chunk_offset": PRODUCTION_CONTRACT.chunks[1].offset,
                "full_bam_sha256_computed": False,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
