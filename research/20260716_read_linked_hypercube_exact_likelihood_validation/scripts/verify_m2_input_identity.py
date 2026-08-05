#!/usr/bin/env python3
"""Verify the immutable direct-input identity for the M2 comprehensive run.

The seven BAM payloads are *not* hashed end-to-end.  They are verified using
the manifest's ``storage_identity_v1`` contract: logical/resolved metadata,
three manifest-defined sampled chunks, and a full SHA-256 of the BAI.  The
tree VCF/index and read-tag sidecar/index are full-hashed.  A POST invocation
can authenticate and compare against a PRE receipt.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import stat
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence


DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
SCHEMA_NAME = "intersubmod.m2_input_identity_verification"
SCHEMA_VERSION = "1.0.0"
EXPECTED_ARTIFACTS_PER_DATASET = 6
EXPECTED_ARTIFACTS = len(DATASETS) * EXPECTED_ARTIFACTS_PER_DATASET
AUTOSOMES = tuple(f"chr{value}" for value in range(1, 23))
STORAGE_CHUNK_SIZE = 1024 * 1024
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
CANONICAL_MANIFEST_SHA256 = "16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45"
CANONICAL_SCHEMA_SHA256 = "47f54d86159beae66b68719e034c7755eba59b2d1719d4546153fa0ae74ff19f"
CANONICAL_SCHEMA_RELATIVE = Path(
    "docs/methodology/_assets/20260627_subclone_4axis_teaching/"
    "schemas/layered_input_manifest_v3.schema.json"
)


def _discover_repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / CANONICAL_SCHEMA_RELATIVE
        if candidate.is_file():
            return parent
    return Path("/big7_disk/liaoyoyo2001/InterSubMod")


CANONICAL_SCHEMA_PATH = _discover_repo_root() / CANONICAL_SCHEMA_RELATIVE


class IdentityError(RuntimeError):
    """A fail-closed manifest, filesystem, receipt, or comparison error."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise IdentityError(message)


def _strict_int(value: Any, label: str, *, minimum: int = 0) -> int:
    require(
        isinstance(value, int) and not isinstance(value, bool) and value >= minimum,
        f"{label} must be an integer >= {minimum}",
    )
    return int(value)


def _validate_full_identity(expected: Mapping[str, Any], label: str) -> None:
    require(set(expected) == {"policy", "size_bytes", "sha256"}, f"{label} full identity schema mismatch")
    require(expected["policy"] == "full_sha256", f"{label} policy is not full_sha256")
    _strict_int(expected["size_bytes"], f"{label}.size_bytes")
    require(isinstance(expected["sha256"], str) and SHA256_RE.fullmatch(expected["sha256"]), f"{label}.sha256 is invalid")


class _MiniSchemaError(Exception):
    def __init__(self, validator: str, path: tuple[Any, ...], message: str) -> None:
        super().__init__(message)
        self.validator = validator
        self.path = path
        self.message = message


def _json_type_matches(value: Any, expected: str) -> bool:
    return {
        "object": isinstance(value, dict),
        "array": isinstance(value, list),
        "string": isinstance(value, str),
        "integer": isinstance(value, int) and not isinstance(value, bool),
        "number": isinstance(value, (int, float)) and not isinstance(value, bool),
        "boolean": isinstance(value, bool),
        "null": value is None,
    }.get(expected, False)


def _resolve_local_ref(root: Mapping[str, Any], reference: str) -> Any:
    require(reference.startswith("#/"), f"only local JSON Schema refs are supported: {reference}")
    value: Any = root
    for token in reference[2:].split("/"):
        token = token.replace("~1", "/").replace("~0", "~")
        require(isinstance(value, dict) and token in value, f"unresolved JSON Schema ref: {reference}")
        value = value[token]
    return value


def _mini_schema_validate(
    instance: Any,
    schema: Any,
    root: Mapping[str, Any],
    path: tuple[Any, ...] = (),
) -> None:
    """Validate every keyword used by the frozen canonical v3 schema."""
    if schema is True:
        return
    if schema is False:
        raise _MiniSchemaError("falseSchema", path, "value is forbidden")
    require(isinstance(schema, dict), f"JSON Schema node at {path} is not an object/boolean")
    supported = {
        "$schema", "$id", "$defs", "$ref", "title", "type", "additionalProperties",
        "required", "properties", "const", "enum", "pattern", "minLength", "minimum",
        "maximum", "format", "minItems", "maxItems", "items", "prefixItems", "oneOf",
        "allOf",
    }
    unknown = set(schema) - supported
    require(not unknown, f"unsupported JSON Schema keywords at {path}: {sorted(unknown)}")
    if "$ref" in schema:
        _mini_schema_validate(instance, _resolve_local_ref(root, schema["$ref"]), root, path)
    for branch in schema.get("allOf", []):
        _mini_schema_validate(instance, branch, root, path)
    if "oneOf" in schema:
        matches = 0
        for branch in schema["oneOf"]:
            try:
                _mini_schema_validate(instance, branch, root, path)
                matches += 1
            except _MiniSchemaError:
                pass
        if matches != 1:
            raise _MiniSchemaError("oneOf", path, f"expected one matching branch, observed {matches}")
    if "const" in schema and instance != schema["const"]:
        raise _MiniSchemaError("const", path, f"{instance!r} is not const {schema['const']!r}")
    if "enum" in schema and instance not in schema["enum"]:
        raise _MiniSchemaError("enum", path, f"{instance!r} is not in {schema['enum']!r}")
    expected_type = schema.get("type")
    if expected_type is not None and not _json_type_matches(instance, expected_type):
        raise _MiniSchemaError("type", path, f"{instance!r} is not type {expected_type}")
    if isinstance(instance, dict):
        required = schema.get("required", [])
        missing = [key for key in required if key not in instance]
        if missing:
            raise _MiniSchemaError("required", path, f"missing required properties: {missing}")
        properties = schema.get("properties", {})
        if schema.get("additionalProperties") is False:
            extras = sorted(set(instance) - set(properties))
            if extras:
                raise _MiniSchemaError("additionalProperties", path, f"unknown properties: {extras}")
        for key, subschema in properties.items():
            if key in instance:
                _mini_schema_validate(instance[key], subschema, root, path + (key,))
    if isinstance(instance, list):
        if "minItems" in schema and len(instance) < schema["minItems"]:
            raise _MiniSchemaError("minItems", path, f"array has {len(instance)} items")
        if "maxItems" in schema and len(instance) > schema["maxItems"]:
            raise _MiniSchemaError("maxItems", path, f"array has {len(instance)} items")
        prefixes = schema.get("prefixItems", [])
        for index, subschema in enumerate(prefixes[: len(instance)]):
            _mini_schema_validate(instance[index], subschema, root, path + (index,))
        if "items" in schema:
            remaining_start = len(prefixes) if prefixes else 0
            for index in range(remaining_start, len(instance)):
                _mini_schema_validate(instance[index], schema["items"], root, path + (index,))
    if isinstance(instance, str):
        if "minLength" in schema and len(instance) < schema["minLength"]:
            raise _MiniSchemaError("minLength", path, "string is too short")
        if "pattern" in schema and re.search(schema["pattern"], instance) is None:
            raise _MiniSchemaError("pattern", path, f"string does not match {schema['pattern']}")
        if schema.get("format") == "date-time":
            try:
                parsed = dt.datetime.fromisoformat(instance.replace("Z", "+00:00"))
            except ValueError as exc:
                raise _MiniSchemaError("format", path, "invalid date-time") from exc
            if parsed.tzinfo is None:
                raise _MiniSchemaError("format", path, "date-time lacks timezone")
    if isinstance(instance, (int, float)) and not isinstance(instance, bool):
        if "minimum" in schema and instance < schema["minimum"]:
            raise _MiniSchemaError("minimum", path, f"{instance} is below minimum")
        if "maximum" in schema and instance > schema["maximum"]:
            raise _MiniSchemaError("maximum", path, f"{instance} is above maximum")


def _validate_against_frozen_schema(
    manifest: Mapping[str, Any], schema_path: Path, expected_schema_sha256: str
) -> dict[str, Any]:
    require(SHA256_RE.fullmatch(expected_schema_sha256) is not None, "expected schema SHA-256 is invalid")
    schema_payload = schema_path.resolve(strict=True).read_bytes()
    observed_sha256 = hashlib.sha256(schema_payload).hexdigest()
    require(observed_sha256 == expected_schema_sha256, "canonical v3 schema SHA-256 mismatch")
    schema = strict_json_load_bytes(schema_payload, schema_path)
    require(isinstance(schema, dict), "canonical v3 schema root is not an object")
    try:
        _mini_schema_validate(manifest, schema, schema)
    except _MiniSchemaError as exc:
        location = ".".join(str(item) for item in exc.path) or "$"
        raise IdentityError(f"canonical v3 schema violation at {location}: {exc.message}") from exc
    return {"path": str(schema_path.resolve()), "sha256": observed_sha256}


def _validate_manifest_contract(manifest: Mapping[str, Any]) -> None:
    require(manifest.get("schema_name") == "intersubmod.layered_input_manifest", "manifest schema_name is not canonical v3")
    require(manifest.get("schema_version") == "3.0.0", "manifest schema_version is not 3.0.0")
    require(
        isinstance(manifest.get("$schema"), str)
        and manifest["$schema"].endswith("/schemas/layered_input_manifest_v3.schema.json"),
        "manifest $schema does not name layered_input_manifest_v3",
    )
    require(_strict_int(manifest.get("dataset_count"), "manifest.dataset_count") == 7, "manifest dataset_count is not 7")
    require(
        _strict_int(manifest.get("biological_sample_count"), "manifest.biological_sample_count") == 6,
        "manifest biological_sample_count is not 6",
    )
    contract = manifest.get("analysis_contract")
    require(isinstance(contract, dict), "manifest analysis_contract is missing")
    expected_contract = {
        "task_type": "comprehensive_validation",
        "read_tag_mode": "external_sidecar",
        "embedded_tag_policy": "ignore",
        "require_exact_join": True,
        "sidecar_identity_schema": "coordinate_join_v1",
        "bam_identity_policy": "storage_identity_v1",
        "bam_identity_assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
        "longphase_input_contract": "normalized_ClairS_raw_all",
        "tree_input_contract": "longphase_s_recalibrated_FILTER_PASS",
    }
    for key, expected in expected_contract.items():
        require(type(contract.get(key)) is type(expected) and contract.get(key) == expected, f"analysis_contract.{key} mismatch")
    scope = contract.get("scope")
    require(
        isinstance(scope, dict)
        and scope.get("name") == "whole_autosomes_chr1_22"
        and scope.get("contigs") == list(AUTOSOMES),
        "analysis contract is not exact chr1-22 scope",
    )
    summary = manifest.get("production_summary")
    require(isinstance(summary, dict), "production_summary missing")
    require(
        _strict_int(summary.get("expected_dataset_count"), "production_summary.expected_dataset_count") == 7
        and _strict_int(summary.get("completed_dataset_count"), "production_summary.completed_dataset_count") == 7
        and _strict_int(summary.get("passed_dataset_count"), "production_summary.passed_dataset_count") == 7
        and summary.get("all_pass") is True,
        "production_summary is not 7/7 PASS",
    )
    summary_rows = summary.get("datasets")
    require(isinstance(summary_rows, list) and len(summary_rows) == 7, "production_summary datasets are invalid")
    require(
        {row.get("sample") for row in summary_rows if isinstance(row, dict) and row.get("pass") is True}
        == set(DATASETS),
        "production_summary dataset PASS set differs",
    )


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for payload in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(payload)
    return digest.hexdigest()


def strict_json_load_bytes(payload: bytes, path: Path) -> Any:
    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise IdentityError(f"duplicate JSON key in {path}: {key}")
            value[key] = item
        return value

    def reject_constant(value: str) -> None:
        raise IdentityError(f"non-finite JSON constant in {path}: {value}")

    try:
        return json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=pairs_hook,
            parse_constant=reject_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise IdentityError(f"cannot parse JSON {path}: {exc}") from exc


def strict_json_load(path: Path) -> Any:
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise IdentityError(f"cannot read JSON {path}: {exc}") from exc
    return strict_json_load_bytes(payload, path)


def _regular_file(path: Path, label: str) -> Path:
    try:
        resolved = path.resolve(strict=True)
        mode = resolved.stat().st_mode
    except OSError as exc:
        raise IdentityError(f"unavailable {label}: {path}: {exc}") from exc
    require(stat.S_ISREG(mode), f"{label} is not a regular file: {path}")
    return resolved


def _stat_identity(path: Path) -> tuple[int, int, int, int, int]:
    value = path.stat()
    return (
        int(value.st_dev), int(value.st_ino), int(value.st_size),
        int(value.st_mtime_ns), int(value.st_ctime_ns),
    )


def _fd_stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        int(value.st_dev), int(value.st_ino), int(value.st_size),
        int(value.st_mtime_ns), int(value.st_ctime_ns),
    )


def _open_verified_fd(
    path: Path,
    *,
    expected_dev_ino: tuple[int, int],
    forbidden_dev_inos: frozenset[tuple[int, int]],
    label: str,
    allow_forbidden_identity: bool = False,
) -> tuple[int, os.stat_result]:
    """Open once, verify the payload identity, then permit any byte reads."""
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    except OSError as exc:
        raise IdentityError(f"cannot open {label}: {path}: {exc}") from exc
    try:
        metadata = os.fstat(descriptor)
        require(stat.S_ISREG(metadata.st_mode), f"{label} is not a regular file: {path}")
        observed_dev_ino = (int(metadata.st_dev), int(metadata.st_ino))
        require(
            observed_dev_ino == expected_dev_ino,
            f"{label} path target changed after role inventory: {path}",
        )
        require(
            allow_forbidden_identity or observed_dev_ino not in forbidden_dev_inos,
            f"refusing to full-hash or sample a BAM through non-BAM role {label}: {path}",
        )
        return descriptor, metadata
    except Exception:
        os.close(descriptor)
        raise


def _sha256_fd(descriptor: int) -> str:
    digest = hashlib.sha256()
    os.lseek(descriptor, 0, os.SEEK_SET)
    while True:
        payload = os.read(descriptor, 1024 * 1024)
        if not payload:
            break
        digest.update(payload)
    return digest.hexdigest()


def full_identity(
    path: Path,
    *,
    expected_dev_ino: tuple[int, int],
    forbidden_dev_inos: frozenset[tuple[int, int]],
    cache: dict[tuple[int, int, int, int, int], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    try:
        logical_before = path.lstat()
    except OSError as exc:
        raise IdentityError(f"cannot lstat full-sha artifact before hashing: {path}: {exc}") from exc
    descriptor, opened = _open_verified_fd(
        path,
        expected_dev_ino=expected_dev_ino,
        forbidden_dev_inos=forbidden_dev_inos,
        label="full-sha artifact",
    )
    before = _fd_stat_identity(opened)
    try:
        # Deliberately re-hash even when the same inode was seen earlier.  A
        # cached digest cannot prove that this verified open still exposes the
        # same bytes, and an early return would skip the logical-path and
        # boundary rechecks below.
        digest = _sha256_fd(descriptor)
        after = _fd_stat_identity(os.fstat(descriptor))
        require(before == after, f"full-sha artifact changed while hashing: {path}")
    finally:
        os.close(descriptor)
    try:
        logical_after = path.lstat()
    except OSError as exc:
        raise IdentityError(f"cannot lstat full-sha artifact after hashing: {path}: {exc}") from exc
    logical_before_identity = (
        int(logical_before.st_dev), int(logical_before.st_ino), int(logical_before.st_size),
        int(logical_before.st_mtime_ns), int(logical_before.st_ctime_ns), stat.S_IFMT(logical_before.st_mode),
    )
    logical_after_identity = (
        int(logical_after.st_dev), int(logical_after.st_ino), int(logical_after.st_size),
        int(logical_after.st_mtime_ns), int(logical_after.st_ctime_ns), stat.S_IFMT(logical_after.st_mode),
    )
    require(
        logical_before_identity == logical_after_identity,
        f"full-sha logical path changed while hashing: {path}",
    )
    boundary_descriptor, _ = _open_verified_fd(
        path,
        expected_dev_ino=expected_dev_ino,
        forbidden_dev_inos=forbidden_dev_inos,
        label="full-sha boundary recheck",
    )
    os.close(boundary_descriptor)
    value = {
        "policy": "full_sha256",
        "size_bytes": before[2],
        "sha256": digest,
    }
    if cache is not None:
        cache[before] = dict(value)
    return value


def _sample_chunk_fd(descriptor: int, path: Path, raw: Mapping[str, Any]) -> dict[str, Any]:
    require(set(raw) == {"label", "offset", "length", "sha256"}, "invalid storage chunk schema")
    label = raw["label"]
    offset = raw["offset"]
    length = raw["length"]
    require(isinstance(label, str) and label, "storage chunk label is invalid")
    require(isinstance(offset, int) and not isinstance(offset, bool) and offset >= 0, "storage chunk offset is invalid")
    require(isinstance(length, int) and not isinstance(length, bool) and length >= 0, "storage chunk length is invalid")
    try:
        payload = os.pread(descriptor, length, offset)
    except OSError as exc:
        raise IdentityError(f"cannot sample BAM chunk {path}@{offset}+{length}: {exc}") from exc
    require(len(payload) == length, f"short BAM chunk {path}@{offset}: {len(payload)} != {length}")
    return {
        "label": label,
        "offset": offset,
        "length": length,
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def _validate_storage_contract(
    bam_path: Path, bai_path: Path, expected: Mapping[str, Any], label: str
) -> None:
    required = {
        "policy", "assurance", "is_full_content_hash", "requested_path", "realpath",
        "logical_is_symlink", "logical_size_bytes", "logical_mtime_ns", "st_dev", "st_ino",
        "size_bytes", "mtime_ns", "ctime_ns", "chunk_size_bytes", "chunks", "index",
        "identity_sha256",
    }
    require(set(expected) == required, f"{label} storage_identity_v1 schema mismatch")
    require(expected["policy"] == "storage_identity_v1", f"{label} identity policy mismatch")
    require(
        expected["assurance"] == "metadata_plus_sampled_chunks_not_full_content_hash"
        and expected["is_full_content_hash"] is False,
        f"{label} identity assurance mismatch",
    )
    require(expected["requested_path"] == str(bam_path), f"{label} requested_path mismatch")
    require(isinstance(expected["realpath"], str) and expected["realpath"], f"{label} realpath invalid")
    require(isinstance(expected["logical_is_symlink"], bool), f"{label} logical_is_symlink invalid")
    for key in (
        "logical_size_bytes", "logical_mtime_ns", "st_dev", "st_ino", "size_bytes",
        "mtime_ns", "ctime_ns",
    ):
        _strict_int(expected[key], f"{label}.{key}")
    require(expected["size_bytes"] > 0, f"{label} BAM size is zero")
    require(
        _strict_int(expected["chunk_size_bytes"], f"{label}.chunk_size_bytes", minimum=1)
        == STORAGE_CHUNK_SIZE,
        f"{label} chunk_size_bytes is not 1 MiB",
    )
    chunks = expected["chunks"]
    require(isinstance(chunks, list) and len(chunks) == 3, f"{label} must contain three chunks")
    length = min(STORAGE_CHUNK_SIZE, int(expected["size_bytes"]))
    offsets = [0, max(0, (int(expected["size_bytes"]) - length) // 2), max(0, int(expected["size_bytes"]) - length)]
    for index, (chunk, expected_label, expected_offset) in enumerate(
        zip(chunks, ("first", "middle", "last"), offsets)
    ):
        require(isinstance(chunk, dict), f"{label}.chunks[{index}] is not an object")
        require(set(chunk) == {"label", "offset", "length", "sha256"}, f"{label}.chunks[{index}] schema mismatch")
        require(chunk["label"] == expected_label, f"{label}.chunks[{index}] label mismatch")
        require(_strict_int(chunk["offset"], f"{label}.chunks[{index}].offset") == expected_offset, f"{label}.chunks[{index}] offset mismatch")
        require(_strict_int(chunk["length"], f"{label}.chunks[{index}].length", minimum=1) == length, f"{label}.chunks[{index}] length mismatch")
        require(isinstance(chunk["sha256"], str) and SHA256_RE.fullmatch(chunk["sha256"]), f"{label}.chunks[{index}] sha256 invalid")
    index = expected["index"]
    require(isinstance(index, dict) and set(index) == {"path", "identity"}, f"{label} index schema mismatch")
    require(index["path"] == str(bai_path), f"{label} BAI path mismatch")
    require(isinstance(index["identity"], dict), f"{label} BAI identity missing")
    _validate_full_identity(index["identity"], f"{label}.index")
    require(isinstance(expected["identity_sha256"], str) and SHA256_RE.fullmatch(expected["identity_sha256"]), f"{label}.identity_sha256 invalid")
    without_digest = dict(expected)
    declared_digest = without_digest.pop("identity_sha256")
    require(canonical_sha256(without_digest) == declared_digest, f"{label} canonical identity digest mismatch")


def storage_identity(
    bam_path: Path,
    bai_path: Path,
    expected: Mapping[str, Any],
    *,
    expected_bam_dev_ino: tuple[int, int],
    expected_bai_dev_ino: tuple[int, int],
    forbidden_bam_dev_inos: frozenset[tuple[int, int]],
    full_hash_cache: dict[tuple[int, int, int, int, int], dict[str, Any]],
) -> dict[str, Any]:
    _validate_storage_contract(bam_path, bai_path, expected, "BAM")
    try:
        logical = bam_path.lstat()
        resolved = _regular_file(bam_path, "BAM")
    except OSError as exc:
        raise IdentityError(f"cannot observe BAM identity {bam_path}: {exc}") from exc
    descriptor, target = _open_verified_fd(
        bam_path,
        expected_dev_ino=expected_bam_dev_ino,
        forbidden_dev_inos=forbidden_bam_dev_inos,
        label="BAM",
        allow_forbidden_identity=True,
    )
    chunks_raw = expected["chunks"]
    index_raw = expected["index"]
    require(isinstance(index_raw, dict) and set(index_raw) == {"path", "identity"}, "invalid BAM index identity")
    require(Path(str(index_raw["path"])) == bai_path, "BAM index path differs from alignment index_path")
    before = _fd_stat_identity(target)
    try:
        observed: dict[str, Any] = {
            "policy": "storage_identity_v1",
            "assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
            "is_full_content_hash": False,
            "requested_path": str(bam_path),
            "realpath": str(resolved),
            "logical_is_symlink": stat.S_ISLNK(logical.st_mode),
            "logical_size_bytes": int(logical.st_size),
            "logical_mtime_ns": int(logical.st_mtime_ns),
            "st_dev": int(target.st_dev),
            "st_ino": int(target.st_ino),
            "size_bytes": int(target.st_size),
            "mtime_ns": int(target.st_mtime_ns),
            "ctime_ns": int(target.st_ctime_ns),
            "chunk_size_bytes": expected["chunk_size_bytes"],
            "chunks": [_sample_chunk_fd(descriptor, resolved, raw) for raw in chunks_raw],
            "index": {
                "path": str(bai_path),
                "identity": full_identity(
                    bai_path,
                    expected_dev_ino=expected_bai_dev_ino,
                    forbidden_dev_inos=forbidden_bam_dev_inos,
                    cache=full_hash_cache,
                ),
            },
        }
        after = _fd_stat_identity(os.fstat(descriptor))
        require(before == after, f"BAM changed while sampling: {resolved}")
    finally:
        os.close(descriptor)
    try:
        logical_after = bam_path.lstat()
        resolved_after = bam_path.resolve(strict=True)
    except OSError as exc:
        raise IdentityError(f"BAM logical path changed while sampling: {bam_path}: {exc}") from exc
    boundary_descriptor, _ = _open_verified_fd(
        bam_path,
        expected_dev_ino=expected_bam_dev_ino,
        forbidden_dev_inos=forbidden_bam_dev_inos,
        label="BAM boundary recheck",
        allow_forbidden_identity=True,
    )
    os.close(boundary_descriptor)
    require(
        resolved_after == resolved
        and stat.S_ISLNK(logical_after.st_mode) == stat.S_ISLNK(logical.st_mode)
        and int(logical_after.st_dev) == int(logical.st_dev)
        and int(logical_after.st_ino) == int(logical.st_ino)
        and int(logical_after.st_size) == int(logical.st_size)
        and int(logical_after.st_mtime_ns) == int(logical.st_mtime_ns)
        and int(logical_after.st_ctime_ns) == int(logical.st_ctime_ns),
        f"BAM logical path changed while sampling: {bam_path}",
    )
    observed["identity_sha256"] = canonical_sha256(observed)
    return observed


def _identity_artifact(
    dataset: str,
    role: str,
    path: Path,
    expected: Mapping[str, Any],
    observed: Mapping[str, Any],
    policy: str,
) -> dict[str, Any]:
    return {
        "dataset": dataset,
        "role": role,
        "policy": policy,
        "path": str(path),
        "expected": expected,
        "observed": observed,
        "match": observed == expected,
    }


def _build_role_inventory(by_dataset: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Validate direct-role schemas and aliases before any payload hashing."""
    inventory: list[dict[str, Any]] = []
    for dataset in DATASETS:
        sample = by_dataset[dataset]
        expected_biological = "HCC1395" if dataset == "HCC1395_DORADO" else dataset
        expected_platform = "ONT_DORADO" if dataset == "HCC1395_DORADO" else "ONT"
        expected_role = "platform_replica" if dataset == "HCC1395_DORADO" else "canonical"
        require(sample.get("biological_id") == expected_biological, f"{dataset} biological_id mismatch")
        require(sample.get("platform") == expected_platform, f"{dataset} platform mismatch")
        require(sample.get("replicate_role") == expected_role, f"{dataset} replicate_role mismatch")
        alignment = sample.get("alignment_payload")
        somatic = sample.get("somatic")
        read_tags = sample.get("read_tags")
        require(isinstance(alignment, dict) and isinstance(somatic, dict) and isinstance(read_tags, dict), f"{dataset} input blocks missing")
        require(
            alignment.get("identity_policy") == "storage_identity_v1"
            and alignment.get("embedded_hp_ps_policy") == "ignore",
            f"{dataset} alignment identity/tag policy mismatch",
        )
        require(somatic.get("backbone_role") == "longphase_s_recalibrated_filter_pass", f"{dataset} somatic backbone mismatch")
        require(
            read_tags.get("authority") == "external_coordinate_sidecar"
            and read_tags.get("identity_schema") == "coordinate_join_v1"
            and read_tags.get("fallback_policy") == "forbidden",
            f"{dataset} read-tag authority contract mismatch",
        )
        bam = Path(str(alignment.get("path", "")))
        bai = Path(str(alignment.get("index_path", "")))
        expected_storage = alignment.get("storage_identity_v1")
        require(isinstance(expected_storage, dict), f"{dataset} storage identity missing")
        _validate_storage_contract(bam, bai, expected_storage, f"{dataset}.alignment_payload")
        tree = somatic.get("tree_vcf")
        sidecar = read_tags.get("sidecar")
        sidecar_index = read_tags.get("index")
        require(isinstance(tree, dict) and isinstance(sidecar, dict) and isinstance(sidecar_index, dict), f"{dataset} direct input roles missing")
        tree_index = tree.get("index")
        require(isinstance(tree_index, dict), f"{dataset} tree index missing")
        role_specs = (
            ("alignment_bam", bam, expected_storage),
            ("alignment_bai", bai, expected_storage["index"]["identity"]),
            ("tree_vcf", Path(str(tree.get("path", ""))), tree.get("identity")),
            ("tree_vcf_index", Path(str(tree_index.get("path", ""))), tree_index.get("identity")),
            ("read_tags_sidecar", Path(str(sidecar.get("path", ""))), sidecar.get("identity")),
            ("read_tags_index", Path(str(sidecar_index.get("path", ""))), sidecar_index.get("identity")),
        )
        for role, path, expected_identity in role_specs:
            require(isinstance(expected_identity, dict), f"{dataset}.{role} expected identity missing")
            if role != "alignment_bam":
                _validate_full_identity(expected_identity, f"{dataset}.{role}")
            resolved = _regular_file(path, f"{dataset}.{role}")
            metadata = resolved.stat()
            inventory.append({
                "dataset": dataset,
                "role": role,
                "path": str(path),
                "resolved_path": str(resolved),
                "st_dev": int(metadata.st_dev),
                "st_ino": int(metadata.st_ino),
            })
    require(len(inventory) == EXPECTED_ARTIFACTS, "role inventory does not contain 42 entries")
    identities = {(row["st_dev"], row["st_ino"]) for row in inventory}
    require(len(identities) == EXPECTED_ARTIFACTS, "direct input roles alias the same resolved file")
    return inventory


def build_snapshot(
    manifest_path: Path,
    *,
    schema_path: Path,
    expected_manifest_sha256: str,
    expected_schema_sha256: str,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    manifest_path = manifest_path.resolve(strict=True)
    try:
        manifest_payload = manifest_path.read_bytes()
    except OSError as exc:
        raise IdentityError(f"cannot read manifest {manifest_path}: {exc}") from exc
    manifest = strict_json_load_bytes(manifest_payload, manifest_path)
    require(isinstance(manifest, dict), "manifest root is not an object")
    manifest_sha256 = hashlib.sha256(manifest_payload).hexdigest()
    require(SHA256_RE.fullmatch(expected_manifest_sha256) is not None, "expected manifest SHA-256 is invalid")
    require(manifest_sha256 == expected_manifest_sha256, "canonical v5 manifest SHA-256 mismatch")
    schema_identity = _validate_against_frozen_schema(
        manifest, schema_path, expected_schema_sha256
    )
    _validate_manifest_contract(manifest)
    samples = manifest.get("samples")
    require(isinstance(samples, list) and len(samples) == 7, "manifest samples are not seven rows")
    by_dataset: dict[str, Mapping[str, Any]] = {}
    for raw in samples:
        require(isinstance(raw, dict), "manifest sample row is not an object")
        dataset = raw.get("sample")
        require(isinstance(dataset, str) and dataset not in by_dataset, "invalid/duplicate manifest sample")
        by_dataset[dataset] = raw
    require(set(by_dataset) == set(DATASETS), "manifest dataset scope differs from frozen seven datasets")

    role_inventory = _build_role_inventory(by_dataset)
    role_identity = {
        (row["dataset"], row["role"]): (row["st_dev"], row["st_ino"])
        for row in role_inventory
    }
    forbidden_bam_dev_inos = frozenset(
        (row["st_dev"], row["st_ino"])
        for row in role_inventory if row["role"] == "alignment_bam"
    )
    require(len(forbidden_bam_dev_inos) == 7, "BAM forbidden-inode set does not contain seven identities")

    artifacts: list[dict[str, Any]] = []
    full_hash_cache: dict[tuple[int, int, int, int, int], dict[str, Any]] = {}
    for dataset in DATASETS:
        sample = by_dataset[dataset]
        alignment = sample.get("alignment_payload")
        somatic = sample.get("somatic")
        read_tags = sample.get("read_tags")
        require(isinstance(alignment, dict) and isinstance(somatic, dict) and isinstance(read_tags, dict), f"{dataset} input blocks missing")

        bam_path = Path(str(alignment.get("path", "")))
        bai_path = Path(str(alignment.get("index_path", "")))
        expected_storage = alignment.get("storage_identity_v1")
        require(isinstance(expected_storage, dict), f"{dataset} storage_identity_v1 missing")
        observed_storage = storage_identity(
            bam_path,
            bai_path,
            expected_storage,
            expected_bam_dev_ino=role_identity[(dataset, "alignment_bam")],
            expected_bai_dev_ino=role_identity[(dataset, "alignment_bai")],
            forbidden_bam_dev_inos=forbidden_bam_dev_inos,
            full_hash_cache=full_hash_cache,
        )
        artifacts.append(_identity_artifact(
            dataset, "alignment_bam", bam_path, expected_storage, observed_storage, "storage_identity_v1"
        ))
        expected_bai = expected_storage["index"]["identity"]
        observed_bai = full_identity(
            bai_path,
            expected_dev_ino=role_identity[(dataset, "alignment_bai")],
            forbidden_dev_inos=forbidden_bam_dev_inos,
            cache=full_hash_cache,
        )
        artifacts.append(_identity_artifact(
            dataset, "alignment_bai", bai_path, expected_bai, observed_bai, "full_sha256"
        ))

        tree_vcf = somatic.get("tree_vcf")
        require(isinstance(tree_vcf, dict), f"{dataset} tree_vcf missing")
        tree_path = Path(str(tree_vcf.get("path", "")))
        expected_tree = tree_vcf.get("identity")
        tree_index = tree_vcf.get("index")
        require(isinstance(expected_tree, dict) and isinstance(tree_index, dict), f"{dataset} tree VCF identity missing")
        observed_tree = full_identity(
            tree_path,
            expected_dev_ino=role_identity[(dataset, "tree_vcf")],
            forbidden_dev_inos=forbidden_bam_dev_inos,
            cache=full_hash_cache,
        )
        artifacts.append(_identity_artifact(
            dataset, "tree_vcf", tree_path, expected_tree, observed_tree, "full_sha256"
        ))
        tree_index_path = Path(str(tree_index.get("path", "")))
        expected_tree_index = tree_index.get("identity")
        require(isinstance(expected_tree_index, dict), f"{dataset} tree VCF index identity missing")
        observed_tree_index = full_identity(
            tree_index_path,
            expected_dev_ino=role_identity[(dataset, "tree_vcf_index")],
            forbidden_dev_inos=forbidden_bam_dev_inos,
            cache=full_hash_cache,
        )
        artifacts.append(_identity_artifact(
            dataset, "tree_vcf_index", tree_index_path, expected_tree_index, observed_tree_index, "full_sha256"
        ))

        sidecar = read_tags.get("sidecar")
        sidecar_index = read_tags.get("index")
        require(isinstance(sidecar, dict) and isinstance(sidecar_index, dict), f"{dataset} read-tag inputs missing")
        sidecar_path = Path(str(sidecar.get("path", "")))
        expected_sidecar = sidecar.get("identity")
        require(isinstance(expected_sidecar, dict), f"{dataset} sidecar identity missing")
        observed_sidecar = full_identity(
            sidecar_path,
            expected_dev_ino=role_identity[(dataset, "read_tags_sidecar")],
            forbidden_dev_inos=forbidden_bam_dev_inos,
            cache=full_hash_cache,
        )
        artifacts.append(_identity_artifact(
            dataset, "read_tags_sidecar", sidecar_path, expected_sidecar, observed_sidecar, "full_sha256"
        ))
        sidecar_index_path = Path(str(sidecar_index.get("path", "")))
        expected_sidecar_index = sidecar_index.get("identity")
        require(isinstance(expected_sidecar_index, dict), f"{dataset} sidecar index identity missing")
        observed_sidecar_index = full_identity(
            sidecar_index_path,
            expected_dev_ino=role_identity[(dataset, "read_tags_index")],
            forbidden_dev_inos=forbidden_bam_dev_inos,
            cache=full_hash_cache,
        )
        artifacts.append(_identity_artifact(
            dataset, "read_tags_index", sidecar_index_path, expected_sidecar_index, observed_sidecar_index, "full_sha256"
        ))

    require(len(artifacts) == EXPECTED_ARTIFACTS, "direct input artifact count is not 42")
    snapshot = {
        "manifest_sha256": manifest_sha256,
        "schema_sha256": schema_identity["sha256"],
        "datasets": list(DATASETS),
        "artifacts": [
            {
                "dataset": row["dataset"],
                "role": row["role"],
                "policy": row["policy"],
                "path": row["path"],
                "observed": row["observed"],
            }
            for row in artifacts
        ],
    }
    audit = {
        "artifacts": artifacts,
        "role_inventory": role_inventory,
        "n_artifacts": len(artifacts),
        "n_unique_resolved_files": len({(row["st_dev"], row["st_ino"]) for row in role_inventory}),
        "n_storage_identity_v1": sum(row["policy"] == "storage_identity_v1" for row in artifacts),
        "n_full_sha256": sum(row["policy"] == "full_sha256" for row in artifacts),
        "n_sampled_bam_chunks": sum(
            len(row["observed"]["chunks"])
            for row in artifacts if row["policy"] == "storage_identity_v1"
        ),
        "n_mismatches": sum(row["match"] is not True for row in artifacts),
    }
    return snapshot, audit, schema_identity


def authenticate_pre_receipt(path: Path) -> tuple[dict[str, Any], str]:
    path = path.resolve(strict=True)
    sidecar = path.with_name(f"{path.name}.sha256")
    require(sidecar.is_file(), f"PRE receipt sidecar missing: {sidecar}")
    fields = sidecar.read_text(encoding="ascii").strip().split()
    require(len(fields) == 2 and fields[1] == path.name, "PRE receipt sidecar is malformed")
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise IdentityError(f"cannot read PRE receipt {path}: {exc}") from exc
    digest = hashlib.sha256(payload).hexdigest()
    require(fields[0] == digest, "PRE receipt sidecar SHA mismatch")
    receipt = strict_json_load_bytes(payload, path)
    require(isinstance(receipt, dict), "PRE receipt is not an object")
    require(receipt.get("schema_name") == SCHEMA_NAME and receipt.get("schema_version") == SCHEMA_VERSION, "PRE receipt schema mismatch")
    require(receipt.get("all_pass") is True, "PRE receipt did not pass")
    require(receipt.get("mode") == "PRE", "compare-to receipt is not PRE mode")
    snapshot = receipt.get("input_identity_snapshot")
    require(isinstance(snapshot, dict), "PRE receipt lacks input snapshot")
    require(receipt.get("input_identity_snapshot_sha256") == canonical_sha256(snapshot), "PRE snapshot digest mismatch")
    return receipt, digest


def verify(
    manifest_path: Path,
    compare_to: Path | None = None,
    *,
    schema_path: Path = CANONICAL_SCHEMA_PATH,
    expected_manifest_sha256: str = CANONICAL_MANIFEST_SHA256,
    expected_schema_sha256: str = CANONICAL_SCHEMA_SHA256,
    _test_authority_override: bool = False,
) -> dict[str, Any]:
    canonical_authority = (
        expected_manifest_sha256 == CANONICAL_MANIFEST_SHA256
        and expected_schema_sha256 == CANONICAL_SCHEMA_SHA256
        and schema_path.resolve(strict=True) == CANONICAL_SCHEMA_PATH.resolve(strict=True)
    )
    require(
        canonical_authority or _test_authority_override,
        "production verification authority is frozen; caller-supplied manifest/schema authority is forbidden",
    )
    snapshot, audit, schema_identity = build_snapshot(
        manifest_path,
        schema_path=schema_path,
        expected_manifest_sha256=expected_manifest_sha256,
        expected_schema_sha256=expected_schema_sha256,
    )
    authenticated = authenticate_pre_receipt(compare_to) if compare_to is not None else None
    pre_receipt = None if authenticated is None else authenticated[0]
    pre_receipt_sha = None if authenticated is None else authenticated[1]
    snapshot_sha = canonical_sha256(snapshot)
    pre_snapshot_sha = None if pre_receipt is None else pre_receipt["input_identity_snapshot_sha256"]
    pre_snapshot_equal = pre_receipt is None or pre_receipt["input_identity_snapshot"] == snapshot
    if pre_receipt is not None:
        pre_eligible = pre_receipt.get("validation_evidence_eligible")
        require(
            isinstance(pre_eligible, bool),
            "PRE receipt validation_evidence_eligible is not a boolean",
        )
        require(
            pre_eligible is canonical_authority,
            "PRE receipt authority mode differs from POST verification authority",
        )
    verifier_path = Path(__file__).resolve()
    checks = {
        "manifest_matches_selected_authority_sha256": snapshot["manifest_sha256"] == expected_manifest_sha256,
        "schema_matches_selected_authority_sha256": schema_identity["sha256"] == expected_schema_sha256,
        "manifest_passes_complete_canonical_v3_json_schema": True,
        "manifest_scope_is_exactly_seven_datasets": snapshot["datasets"] == list(DATASETS),
        "direct_input_artifact_count_is_42": audit["n_artifacts"] == EXPECTED_ARTIFACTS,
        "all_42_direct_input_roles_resolve_to_distinct_files": (
            audit["n_unique_resolved_files"] == EXPECTED_ARTIFACTS
        ),
        "seven_bams_use_storage_identity_v1": audit["n_storage_identity_v1"] == 7,
        "thirty_five_non_bam_artifact_roles_use_full_sha256": audit["n_full_sha256"] == 35,
        "twenty_one_manifest_defined_bam_chunks_recomputed": audit["n_sampled_bam_chunks"] == 21,
        "all_manifest_identities_match_observed_files": audit["n_mismatches"] == 0,
        "post_snapshot_exactly_matches_authenticated_pre": pre_snapshot_equal,
    }
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "mode": "POST_COMPARE" if compare_to is not None else "PRE",
        "authority_mode": "CANONICAL_V5_FROZEN" if canonical_authority else "TEST_ONLY_UNFROZEN",
        "validation_evidence_eligible": canonical_authority,
        "authority": {
            "canonical_manifest_sha256": CANONICAL_MANIFEST_SHA256,
            "canonical_schema_sha256": CANONICAL_SCHEMA_SHA256,
            "canonical_schema_path": str(CANONICAL_SCHEMA_PATH.resolve(strict=True)),
            "selected_authority_is_canonical": canonical_authority,
            "test_only_override": _test_authority_override,
        },
        "scope": {
            "technical_datasets": 7,
            "biological_samples": 6,
            "chromosomes": "chr1-chr22",
            "tasks": 154,
            "datasets": list(DATASETS),
            "direct_input_artifacts": EXPECTED_ARTIFACTS,
        },
        "manifest": {
            "path": str(manifest_path.resolve()),
            "sha256": snapshot["manifest_sha256"],
            "expected_sha256": expected_manifest_sha256,
        },
        "canonical_schema": {
            **schema_identity,
            "expected_sha256": expected_schema_sha256,
        },
        "verifier": {"path": str(verifier_path), "sha256": sha256_path(verifier_path)},
        "assurance": {
            "bam": "storage_identity_v1 metadata plus three manifest-defined sampled chunks; no full BAM hash",
            "other_direct_inputs": "full SHA-256",
            "pre_post": (
                "canonical boundary snapshot exact equality; this detects persistent drift but does not prove "
                "that no transient mutation was made and restored between PRE and POST"
            ),
            "temporal_immutability_proven": False,
        },
        "input_identity_snapshot": snapshot,
        "input_identity_snapshot_sha256": snapshot_sha,
        "artifact_audit": audit,
        "compare_to": None if compare_to is None else {
            "path": str(compare_to.resolve()),
            "sha256": pre_receipt_sha,
            "pre_snapshot_sha256": pre_snapshot_sha,
            "exact_snapshot_equal": pre_snapshot_equal,
        },
        "checks": checks,
    }
    receipt["all_pass"] = bool(checks) and all(value is True for value in checks.values())
    return receipt


def _atomic_write_new(path: Path, payload: bytes) -> None:
    """Publish a fully flushed file atomically and refuse a concurrent overwrite."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    temp_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temp_path, path)
        except FileExistsError as exc:
            raise IdentityError(f"refusing to overwrite input identity evidence: {path}") from exc
    finally:
        if temp_path.exists():
            temp_path.unlink()


def write_receipt(path: Path, receipt: dict[str, Any]) -> str:
    path = path.resolve()
    sidecar = path.with_name(f"{path.name}.sha256")
    require(not path.exists() and not sidecar.exists(), f"refusing to overwrite input identity evidence: {path}")
    receipt["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": sidecar.name,
        "covers": path.name,
    }
    payload = json.dumps(receipt, ensure_ascii=False, allow_nan=False, indent=2).encode("utf-8") + b"\n"
    _atomic_write_new(path, payload)
    digest = sha256_path(path)
    _atomic_write_new(sidecar, f"{digest}  {path.name}\n".encode("ascii"))
    return digest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--compare-to", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    failure: str | None = None
    try:
        receipt = verify(
            args.manifest,
            args.compare_to,
        )
    except (IdentityError, OSError, ValueError, TypeError, KeyError, OverflowError) as exc:
        failure = f"{type(exc).__name__}: {exc}"
        receipt = {
            "schema_name": SCHEMA_NAME,
            "schema_version": SCHEMA_VERSION,
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "mode": "POST_COMPARE" if args.compare_to is not None else "PRE",
            "manifest": {"path": str(args.manifest.resolve())},
            "compare_to": None if args.compare_to is None else {"path": str(args.compare_to.resolve())},
            "failure": failure,
            "checks": {"verification_completed_without_contract_violation": False},
            "all_pass": False,
        }
    digest = write_receipt(args.output, receipt)
    print(json.dumps({
        "output": str(args.output.resolve()),
        "sha256": digest,
        "mode": receipt.get("mode"),
        "all_pass": receipt["all_pass"],
        "n_artifacts": (receipt.get("artifact_audit") or {}).get("n_artifacts"),
        "failure": failure,
    }, ensure_ascii=False))
    return 0 if receipt["all_pass"] and receipt.get("validation_evidence_eligible") is True else 1


if __name__ == "__main__":
    raise SystemExit(main())
