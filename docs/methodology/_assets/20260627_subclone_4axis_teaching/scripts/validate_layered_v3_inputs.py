#!/usr/bin/env python3
"""Fail-closed validator and frozen-lock writer for layered input manifest v3."""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import stat
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


SCHEMA_NAME = "intersubmod.layered_input_manifest"
SCHEMA_VERSION = "3.0.0"
LOCK_SCHEMA_NAME = "intersubmod.layered_input_lock"
LOCK_SCHEMA_VERSION = "1.0.0"
SIDECAR_VALIDATION_SCHEMA = "intersubmod.longphase_coordinate_sidecar_validation"
SIDECAR_VALIDATION_VERSION = "1.0.0"
SIDECAR_SUBJECT_SCHEMA = "intersubmod.coordinate_join_subject"
SIDECAR_SUBJECT_VERSION = "1.0.0"
PRODUCER_RECEIPT_SCHEMA = "intersubmod.longphase_production_capture_receipt"
PRODUCER_RECEIPT_VERSION = "1.0.0"
RAW_ALL_PRODUCER_RECEIPT_SCHEMA = "intersubmod.longphase_raw_all_capture_receipt"
RAW_ALL_PRODUCER_RECEIPT_VERSION = "2.0.0"
RAW_ALL_INPUT_CONTRACT = "normalized_ClairS_raw_all"
TREE_INPUT_CONTRACT = "longphase_s_recalibrated_FILTER_PASS"
SENSITIVITY_TREE_INPUT_CONTRACT = "clairs_FILTER_PASS_sensitivity"
RAW_ALL_TAGGING_SEMANTICS = "longphase_s_raw_all_bidirectional_recalibration_tags"
REDUNDANT_IDENTITY_POLICY = "collapse_redundant_rows_with_identical_HP_PS"
METHOD_RECEIPT_REPO_ROOT: Path | None = None
COORDINATE_IDENTITY_COLUMNS = ["QNAME", "CHROM", "START0", "END0", "FLAG", "CIGAR_B2"]
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
EXPECTED_DATASETS = {
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
STORAGE_CHUNK_SIZE = 1024 * 1024
PRODUCTION_POLICY = {
    "name": "production_default_filter",
    "longphase_internal_filter_enabled": True,
    "disable_filter_flag_present": False,
    "truth_flags_present": False,
    "forbidden_truth_options_absent": {"--truth-vcf": True, "--truth-bed": True},
}
PRODUCTION_EFFECTIVE_OPTIONS = {
    "threads": {"value": 12, "origin": "explicit"},
    "mapq_min": {"value": 20, "origin": "explicit"},
    "tag_supplementary": {"value": True, "origin": "explicit"},
    "output_somatic_vcf": {"value": True, "origin": "explicit"},
    "truth_vcf": {"value": False, "origin": "forbidden_absent"},
    "truth_bed": {"value": False, "origin": "forbidden_absent"},
    "disable_filter": {"value": False, "origin": "absent_default"},
    "purity_read_assignment": {"value": 0.6, "origin": "default"},
}


class ContractError(Exception):
    """An expected contract rejection with a stable machine code."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        exit_code: int = 2,
        stage: str = "schema",
        sample: str | None = None,
        artifact: str | None = None,
        expected: Any = None,
        observed: Any = None,
    ) -> None:
        super().__init__(message)
        self.code = code
        self.exit_code = exit_code
        self.stage = stage
        self.sample = sample
        self.artifact = artifact
        self.expected = expected
        self.observed = observed

    def record(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "stage": self.stage,
            "sample": self.sample,
            "artifact": self.artifact,
            "expected": self.expected,
            "observed": self.observed,
            "message": str(self),
        }


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def full_identity(path: Path) -> dict[str, Any]:
    ensure_regular_file(path, "artifact")
    return {
        "policy": "full_sha256",
        "size_bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def _sampled_chunk(path: Path, label: str, offset: int, length: int) -> dict[str, Any]:
    with path.open("rb") as handle:
        handle.seek(offset)
        payload = handle.read(length)
    if len(payload) != length:
        raise ContractError(
            "E_STORAGE_IDENTITY",
            f"short sampled chunk at {offset}: expected {length}, observed {len(payload)}",
            exit_code=3,
            stage="artifact_identity",
            artifact=str(path),
        )
    return {
        "label": label,
        "offset": offset,
        "length": length,
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def storage_identity_v1(path: Path, index_path: Path) -> dict[str, Any]:
    """Return a metadata + sampled-chunk identity. This is explicitly not a full BAM hash."""
    ensure_regular_file(path, "alignment_payload")
    ensure_regular_file(index_path, "alignment_payload_index")
    logical = path.lstat()
    resolved = path.resolve(strict=True)
    target = resolved.stat()
    size = target.st_size
    length = min(STORAGE_CHUNK_SIZE, size)
    middle = max(0, (size - length) // 2)
    last = max(0, size - length)
    value: dict[str, Any] = {
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
        "size_bytes": size,
        "mtime_ns": target.st_mtime_ns,
        "ctime_ns": target.st_ctime_ns,
        "chunk_size_bytes": STORAGE_CHUNK_SIZE,
        "chunks": [
            _sampled_chunk(resolved, "first", 0, length),
            _sampled_chunk(resolved, "middle", middle, length),
            _sampled_chunk(resolved, "last", last, length),
        ],
        "index": {"path": str(index_path), "identity": full_identity(index_path)},
    }
    value["identity_sha256"] = canonical_sha256(value)
    return value


def strict_json_load(path: Path) -> Any:
    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ContractError("E_JSON_DUPLICATE_KEY", f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise ContractError("E_JSON_NONFINITE", f"non-finite JSON number: {value}")

    try:
        return json.loads(
            path.read_text(encoding="utf-8", errors="strict"),
            object_pairs_hook=pairs_hook,
            parse_constant=reject_constant,
        )
    except ContractError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ContractError("E_JSON_PARSE", f"cannot parse {path}: {error}") from error


def ensure_regular_file(path: Path, label: str) -> None:
    try:
        mode = path.stat().st_mode
    except OSError as error:
        raise ContractError(
            "E_REQUIRED_ARTIFACT",
            f"{label} is unavailable: {path}: {error}",
            exit_code=3,
            stage="artifact_identity",
            artifact=str(path),
        ) from error
    if not stat.S_ISREG(mode):
        raise ContractError(
            "E_ARTIFACT_TYPE",
            f"{label} is not a regular file: {path}",
            exit_code=3,
            stage="artifact_identity",
            artifact=str(path),
        )


def atomic_write_json(path: Path, value: Any, *, mode: int = 0o644) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(json.dumps(value, ensure_ascii=False, allow_nan=False, indent=2).encode("utf-8"))
            handle.write(b"\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary_path, mode)
        os.replace(temporary_path, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except Exception:
        # Deliberately retain a .partial file for audit instead of deleting it.
        raise


def atomic_publish_frozen_lock(path: Path, lock: dict[str, Any], lock_schema_path: Path) -> None:
    candidate = path.parent / f".{path.name}.{os.getpid()}.{time.time_ns()}.candidate"
    atomic_write_json(candidate, lock, mode=0o444)
    readback = strict_json_load(candidate)
    apply_json_schema(readback, lock_schema_path)
    digest_payload = dict(readback)
    expected_digest = digest_payload.pop("lock_sha256")
    if canonical_sha256(digest_payload) != expected_digest:
        raise ContractError("E_LOCK_READBACK", "frozen lock candidate digest readback failed", exit_code=7, stage="lifecycle")
    if path.exists():
        raise ContractError("E_OUTPUT_EXISTS", f"frozen lock appeared during validation: {path}", exit_code=7, stage="lifecycle")
    os.replace(candidate, path)
    try:
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except OSError as error:
        failed_path = path.parent / f".{path.name}.{os.getpid()}.{time.time_ns()}.publish_failed"
        os.replace(path, failed_path)
        raise ContractError(
            "E_LOCK_PUBLISH", f"directory fsync failed; lock quarantined as {failed_path}: {error}",
            exit_code=7, stage="lifecycle",
        ) from error


def artifact(path: Path) -> dict[str, Any]:
    return {"path": str(path), "identity": full_identity(path)}


def find_index(path: Path) -> Path:
    candidates = [Path(f"{path}.tbi"), Path(f"{path}.csi")]
    existing = [candidate for candidate in candidates if candidate.is_file()]
    if len(existing) != 1:
        raise ContractError(
            "E_INDEX_MISMATCH",
            f"expected exactly one index for {path}; found {[str(x) for x in existing]}",
            exit_code=3,
            stage="artifact_identity",
            artifact=str(path),
        )
    return existing[0]


def indexed_artifact(path: Path) -> dict[str, Any]:
    index_path = find_index(path)
    result = artifact(path)
    result["index"] = artifact(index_path)
    return result


def inspect_coordinate_sidecar(path: Path, sample: str) -> dict[str, int]:
    expected_header = ["#CHROM", "START0", "END0", "QNAME", "FLAG", "MAPQ", "CIGAR_B2", "HP", "PS"]
    seen_at_locus: dict[tuple[str, ...], tuple[str, str]] = {}
    closed_chromosomes: set[str] = set()
    current_chrom: str | None = None
    current_start = -1
    rows = unique = duplicates = conflicts = malformed = 0
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            if header != expected_header:
                raise ContractError(
                    "E_SIDECAR_PARSE", f"unexpected 9-column sidecar header: {header}",
                    exit_code=5, stage="sidecar_join", sample=sample, artifact=str(path),
                )
            for line_number, line in enumerate(handle, start=2):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    malformed += 1
                    continue
                try:
                    start = int(fields[1]); end = int(fields[2]); int(fields[4]); int(fields[5])
                except ValueError:
                    malformed += 1
                    continue
                chrom = fields[0]
                if start < 0 or end < start:
                    malformed += 1
                    continue
                if chrom != current_chrom:
                    if chrom in closed_chromosomes:
                        raise ContractError(
                            "E_SIDECAR_PARSE",
                            f"coordinate sidecar returns to closed chromosome at line {line_number}: {chrom}",
                            exit_code=5, stage="sidecar_join", sample=sample, artifact=str(path),
                        )
                    if current_chrom is not None:
                        closed_chromosomes.add(current_chrom)
                    current_chrom = chrom
                    current_start = start
                    seen_at_locus = {}
                elif start < current_start:
                    raise ContractError(
                        "E_SIDECAR_PARSE",
                        f"coordinate sidecar start decreases at line {line_number}: {start} < {current_start}",
                        exit_code=5, stage="sidecar_join", sample=sample, artifact=str(path),
                    )
                elif start > current_start:
                    current_start = start
                    seen_at_locus = {}
                key = (fields[3], fields[2], fields[4], fields[6])
                tag = (fields[7], fields[8])
                if key in seen_at_locus:
                    duplicates += 1
                    conflicts += seen_at_locus[key] != tag
                else:
                    seen_at_locus[key] = tag
                    unique += 1
                rows += 1
    except ContractError:
        raise
    except (OSError, UnicodeError) as error:
        raise ContractError(
            "E_SIDECAR_PARSE", f"cannot read coordinate sidecar: {error}",
            exit_code=5, stage="sidecar_join", sample=sample, artifact=str(path),
        ) from error
    if malformed:
        raise ContractError(
            "E_SIDECAR_PARSE", f"sidecar contains {malformed} malformed rows",
            exit_code=5, stage="sidecar_join", sample=sample, artifact=str(path),
        )
    return {
        "mapped_alignment_count": rows,
        "identity_unique_count": unique,
        "duplicate_count": duplicates,
        "conflict_count": conflicts,
    }


def validate_redundant_coordinate_counts(counts: dict[str, Any], sample: str) -> None:
    mapped = counts.get("mapped_alignment_count")
    unique = counts.get("identity_unique_count")
    redundant = counts.get("duplicate_count")
    conflicts = counts.get("conflict_count")
    if (
        isinstance(mapped, bool)
        or not isinstance(mapped, int)
        or mapped <= 0
        or isinstance(unique, bool)
        or not isinstance(unique, int)
        or unique <= 0
        or isinstance(redundant, bool)
        or not isinstance(redundant, int)
        or redundant < 0
        or conflicts != 0
        or mapped != unique + redundant
    ):
        raise ContractError(
            "E_JOIN_MULTIPLICITY" if conflicts == 0 else "E_JOIN_CONFLICT",
            "coordinate rows must conserve as unique + redundant with zero HP/PS conflicts",
            exit_code=5,
            stage="sidecar_join",
            sample=sample,
            expected={"mapped_equals_unique_plus_redundant": True, "conflict_count": 0},
            observed=counts,
        )


def _pre_schema_checks(manifest: Any) -> None:
    if not isinstance(manifest, dict):
        raise ContractError("E_SCHEMA_TYPE", "manifest root must be an object")
    if manifest.get("schema_name") != SCHEMA_NAME or manifest.get("schema_version") != SCHEMA_VERSION:
        raise ContractError(
            "E_SCHEMA_VERSION",
            f"requires {SCHEMA_NAME} {SCHEMA_VERSION}",
            expected={"schema_name": SCHEMA_NAME, "schema_version": SCHEMA_VERSION},
            observed={"schema_name": manifest.get("schema_name"), "schema_version": manifest.get("schema_version")},
        )
    samples = manifest.get("samples")
    if manifest.get("dataset_count") == 0 or samples == []:
        raise ContractError("E_SCHEMA_EMPTY", "zero-dataset manifests are forbidden")
    if not isinstance(samples, list):
        raise ContractError("E_SCHEMA_TYPE", "samples must be an array")
    ids = [item.get("sample") for item in samples if isinstance(item, dict)]
    if len(ids) != len(set(ids)):
        raise ContractError("E_SAMPLE_DUPLICATE", "sample IDs must be unique")
    unsafe = [value for value in ids if not isinstance(value, str) or SAFE_ID.fullmatch(value) is None]
    if unsafe:
        raise ContractError("E_SAMPLE_ID_UNSAFE", f"unsafe sample IDs: {unsafe}")


def _schema_error_code(error: Any) -> str:
    path = ".".join(str(item) for item in error.absolute_path)
    if error.validator == "additionalProperties":
        return "E_UNKNOWN_FIELD"
    if "schema_version" in path or "schema_name" in path:
        return "E_SCHEMA_VERSION"
    if "scope" in path:
        return "E_SCOPE_NOT_COMPREHENSIVE"
    if "truth_flags_present" in path or "forbidden_truth_options_absent" in path:
        return "E_TRUTH_FLAG_PRESENT"
    if "production_filter_policy" in path or "producer_policy" in path:
        return "E_PRODUCTION_FILTER_POLICY"
    if "production_summary" in path:
        return "E_PRODUCTION_SUMMARY"
    if error.validator == "required":
        return "E_REQUIRED_FIELD"
    return "E_SCHEMA"


class _MiniSchemaError(Exception):
    def __init__(self, validator: str, path: tuple[Any, ...], message: str) -> None:
        super().__init__(message)
        self.validator = validator
        self.absolute_path = path
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


def _resolve_local_ref(root: dict[str, Any], reference: str) -> Any:
    if not reference.startswith("#/"):
        raise ContractError("E_SCHEMA_DEFINITION", f"only local schema refs are supported: {reference}")
    value: Any = root
    for token in reference[2:].split("/"):
        token = token.replace("~1", "/").replace("~0", "~")
        if not isinstance(value, dict) or token not in value:
            raise ContractError("E_SCHEMA_DEFINITION", f"unresolved schema ref: {reference}")
        value = value[token]
    return value


def _mini_schema_validate(instance: Any, schema: Any, root: dict[str, Any], path: tuple[Any, ...] = ()) -> None:
    if schema is True:
        return
    if schema is False:
        raise _MiniSchemaError("falseSchema", path, "value is forbidden")
    if not isinstance(schema, dict):
        raise ContractError("E_SCHEMA_DEFINITION", f"schema at {path} is not an object/boolean")
    supported = {
        "$schema", "$id", "$defs", "$ref", "title", "type", "additionalProperties", "required",
        "properties", "const", "enum", "pattern", "minLength", "minimum", "maximum", "format",
        "minItems", "maxItems", "items", "prefixItems", "oneOf", "allOf",
    }
    unknown = set(schema) - supported
    if unknown:
        raise ContractError("E_SCHEMA_DEFINITION", f"unsupported schema keywords at {path}: {sorted(unknown)}")
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
            raise _MiniSchemaError("oneOf", path, f"expected exactly one schema branch, observed {matches}")
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
                dt.datetime.fromisoformat(instance.replace("Z", "+00:00"))
            except ValueError as error:
                raise _MiniSchemaError("format", path, "invalid date-time") from error
    if isinstance(instance, (int, float)) and not isinstance(instance, bool):
        if "minimum" in schema and instance < schema["minimum"]:
            raise _MiniSchemaError("minimum", path, f"{instance} is below minimum")
        if "maximum" in schema and instance > schema["maximum"]:
            raise _MiniSchemaError("maximum", path, f"{instance} is above maximum")


def apply_json_schema(manifest: Any, schema_path: Path) -> None:
    schema = strict_json_load(schema_path)
    try:
        _mini_schema_validate(manifest, schema, schema)
    except _MiniSchemaError as error:
        path = ".".join(str(item) for item in error.absolute_path) or "$"
        code = _schema_error_code(error)
        exit_code = 4 if code in {"E_TRUTH_FLAG_PRESENT", "E_PRODUCTION_FILTER_POLICY"} else 2
        stage = "producer_scope" if exit_code == 4 else "schema"
        raise ContractError(code, f"schema violation at {path}: {error.message}", exit_code=exit_code, stage=stage)


def validate_semantics(manifest: dict[str, Any]) -> None:
    samples = manifest["samples"]
    observed = {item["sample"]: item["biological_id"] for item in samples}
    if manifest["dataset_count"] != 7 or len(samples) != 7:
        raise ContractError(
            "E_COUNT_MISMATCH", "comprehensive contract requires exactly 7 datasets", expected=7, observed=len(samples)
        )
    if observed != EXPECTED_DATASETS:
        raise ContractError(
            "E_DATASET_SET_MISMATCH",
            "dataset/biological-ID mapping is not the canonical 7-dataset/6-biological set",
            expected=EXPECTED_DATASETS,
            observed=observed,
        )
    if manifest["biological_sample_count"] != 6 or len(set(observed.values())) != 6:
        raise ContractError("E_COUNT_MISMATCH", "comprehensive contract requires exactly 6 biological IDs")
    scope = manifest["analysis_contract"]["scope"]
    if scope["name"] != "whole_autosomes_chr1_22" or scope["contigs"] != AUTOSOMES:
        raise ContractError("E_SCOPE_NOT_COMPREHENSIVE", "scope must be ordered chr1 through chr22")
    if manifest["analysis_contract"]["production_filter_policy"] != PRODUCTION_POLICY:
        raise ContractError(
            "E_PRODUCTION_FILTER_POLICY",
            "analysis production policy is not production_default_filter",
            exit_code=4,
            stage="producer_scope",
        )
    analysis = manifest["analysis_contract"]
    tree_contract = analysis.get("tree_input_contract")
    task_contract_ok = (
        (analysis.get("task_type") == "comprehensive_validation" and tree_contract == TREE_INPUT_CONTRACT)
        or (
            analysis.get("task_type") == "backbone_sensitivity"
            and tree_contract == SENSITIVITY_TREE_INPUT_CONTRACT
        )
    )
    if (
        analysis.get("longphase_input_contract") != RAW_ALL_INPUT_CONTRACT
        or not task_contract_ok
        or analysis.get("tagging_semantics") != RAW_ALL_TAGGING_SEMANTICS
    ):
        raise ContractError(
            "E_LONGPHASE_INPUT_CONTRACT",
            "layered v3 requires raw-all LongPhase input and an explicit canonical/sensitivity tree role",
            exit_code=4,
            stage="producer_scope",
        )
    for meta in samples:
        somatic = meta["somatic"]
        selected = somatic["tree_vcf"]
        if tree_contract == TREE_INPUT_CONTRACT:
            expected_tree = somatic["longphase_recalibrated_pass_vcf"]
            expected_role = "longphase_s_recalibrated_filter_pass"
        else:
            expected_tree = somatic["caller_pass_baseline_vcf"]
            expected_role = "clairs_filter_pass_sensitivity"
        if somatic.get("backbone_role") != expected_role or selected != expected_tree:
            raise ContractError(
                "E_TREE_VCF_ROLE",
                f"{meta['sample']} selected tree artifact does not match {tree_contract}",
                exit_code=4,
                stage="producer_scope",
                sample=meta["sample"],
            )
        if somatic["longphase_recalibrated_pass_vcf"] == somatic["longphase_recalibrated_all_vcf"]:
            raise ContractError(
                "E_TREE_VCF_ROLE",
                f"{meta['sample']} LongPhase-S PASS artifact must be distinct from recalibrated-all",
                exit_code=4,
                stage="producer_scope",
                sample=meta["sample"],
            )
    if (
        manifest["analysis_contract"]["sidecar_identity_schema"] != "coordinate_join_v1"
        or manifest["analysis_contract"]["sidecar_assurance"]
        != "bounded_coordinate_equivalence_not_global_payload_identity"
        or manifest["analysis_contract"].get("duplicate_identity_policy")
        != REDUNDANT_IDENTITY_POLICY
    ):
        raise ContractError(
            "E_SIDECAR_CONTRACT_UNSUPPORTED",
            "coordinate_join_v1 requires bounded assurance and explicit identical-tag collapse",
            exit_code=5,
            stage="sidecar_join",
        )
    summary = manifest["production_summary"]
    summary_ids = [item["sample"] for item in summary["datasets"]]
    if (
        summary["expected_dataset_count"] != 7
        or summary["completed_dataset_count"] != 7
        or summary["passed_dataset_count"] != 7
        or summary["all_pass"] is not True
        or len(summary_ids) != 7
        or set(summary_ids) != set(EXPECTED_DATASETS)
        or any(item.get("longphase_input_contract") != RAW_ALL_INPUT_CONTRACT for item in summary["datasets"])
    ):
        raise ContractError("E_PRODUCTION_SUMMARY", "production summary is not an exact 7/7 PASS")


def _compare(expected: Any, observed: Any, *, code: str, message: str, exit_code: int, stage: str, sample: str, artifact_name: str) -> None:
    if expected != observed:
        raise ContractError(
            code,
            message,
            exit_code=exit_code,
            stage=stage,
            sample=sample,
            artifact=artifact_name,
            expected=expected,
            observed=observed,
        )


def verify_artifact(spec: dict[str, Any], sample: str, name: str) -> dict[str, Any]:
    path = Path(spec["path"])
    observed = full_identity(path)
    _compare(
        spec["identity"], observed, code="E_HASH_MISMATCH", message=f"{name} identity changed",
        exit_code=3, stage="artifact_identity", sample=sample, artifact_name=str(path),
    )
    return {"path": str(path), "identity": observed}


def verify_indexed_artifact(spec: dict[str, Any], sample: str, name: str) -> dict[str, Any]:
    result = verify_artifact(spec, sample, name)
    result["index"] = verify_artifact(spec["index"], sample, f"{name}.index")
    return result


def _argv_contains(argv: list[str], option: str) -> bool:
    return any(token == option or token.startswith(f"{option}=") for token in argv)


def _argv_value(argv: list[str], options: tuple[str, ...]) -> str | None:
    matches: list[str] = []
    for index, token in enumerate(argv):
        for option in options:
            if token == option:
                if index + 1 >= len(argv):
                    raise ContractError(
                        "E_PRODUCER_OPTION_MISMATCH", f"missing value after {option}",
                        exit_code=4, stage="producer_scope",
                    )
                matches.append(argv[index + 1])
            elif token.startswith(f"{option}="):
                matches.append(token.split("=", 1)[1])
    if len(matches) > 1:
        raise ContractError(
            "E_PRODUCER_OPTION_MISMATCH", f"duplicate option aliases: {options}",
            exit_code=4, stage="producer_scope",
        )
    return matches[0] if matches else None


def _validate_sidecar_document(
    validation: dict[str, Any],
    sample: str,
    expected_paths: dict[str, str],
    expected_hashes: dict[str, str],
    expected_storage_sha256: str,
) -> dict[str, Any]:
    if (
        validation.get("schema_name") != SIDECAR_VALIDATION_SCHEMA
        or validation.get("schema_version") != SIDECAR_VALIDATION_VERSION
    ):
        raise ContractError(
            "E_SIDECAR_CONTRACT_UNSUPPORTED",
            "capture validation must implement coordinate-sidecar validation v1",
            exit_code=5,
            stage="sidecar_join",
            sample=sample,
        )
    if validation.get("sample") != sample or validation.get("region") != "all" or validation.get("pass") is not True:
        raise ContractError(
            "E_SIDECAR_VALIDATION",
            "sidecar validation is not an all-region PASS for this sample",
            exit_code=5,
            stage="sidecar_join",
            sample=sample,
        )
    if (
        validation.get("identity_schema") != "coordinate_join_v1"
        or validation.get("coordinate_identity_columns") != COORDINATE_IDENTITY_COLUMNS
    ):
        raise ContractError(
            "E_SIDECAR_CONTRACT_UNSUPPORTED",
            "9-column sidecar must declare the six coordinate/CIGAR identity columns",
            exit_code=5,
            stage="sidecar_join",
            sample=sample,
        )
    duplicate_count = validation.get("duplicate_exact_alignment_rows")
    conflict_count = validation.get("exact_alignment_identity_conflicts")
    mapped_count = validation.get("mapped_alignment_count")
    unique_count = validation.get("identity_unique_count")
    if duplicate_count != 0 or conflict_count != 0:
        raise ContractError(
            "E_JOIN_DUPLICATE" if duplicate_count != 0 else "E_JOIN_CONFLICT",
            "coordinate identities must be globally unique and conflict-free",
            exit_code=5,
            stage="sidecar_join",
            sample=sample,
        )
    if not isinstance(mapped_count, int) or mapped_count <= 0 or unique_count != mapped_count:
        raise ContractError(
            "E_JOIN_MULTIPLICITY",
            "capture validation must prove positive one-to-one global identity counts",
            exit_code=5,
            stage="sidecar_join",
            sample=sample,
        )
    producer = validation.get("producer")
    if not isinstance(producer, dict):
        raise ContractError(
            "E_PRODUCER_AMBIGUOUS", "capture validation lacks producer argv/input binding",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    policy = producer.get("production_policy")
    if not isinstance(policy, dict):
        raise ContractError(
            "E_PRODUCTION_POLICY_MISSING", "sidecar validation lacks machine-readable production policy",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if policy.get("truth_flags_present") is not False or policy.get("forbidden_truth_options_absent") != PRODUCTION_POLICY["forbidden_truth_options_absent"]:
        raise ContractError(
            "E_TRUTH_FLAG_PRESENT", "truth flags are present or not proven absent",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if policy != PRODUCTION_POLICY:
        raise ContractError(
            "E_PRODUCTION_FILTER_POLICY", "producer did not use production_default_filter",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    argv = producer.get("command_argv")
    if not isinstance(argv, list) or not argv or not all(isinstance(token, str) and token for token in argv):
        raise ContractError(
            "E_PRODUCER_AMBIGUOUS", "producer command_argv must be a non-empty string array",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if any(_argv_contains(argv, option) for option in ("--truth-vcf", "--truth-bed")):
        raise ContractError(
            "E_TRUTH_FLAG_PRESENT", "producer argv contains a forbidden truth option",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if _argv_contains(argv, "--disableFilter"):
        raise ContractError(
            "E_PRODUCTION_FILTER_POLICY", "producer argv disables the production default filter",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    effective_options = producer.get("effective_options")
    if effective_options != PRODUCTION_EFFECTIVE_OPTIONS:
        raise ContractError(
            "E_PRODUCER_OPTION_MISMATCH",
            "producer effective options/origins do not match production_default_filter",
            exit_code=4,
            stage="producer_scope",
            sample=sample,
            expected=PRODUCTION_EFFECTIVE_OPTIONS,
            observed=effective_options,
        )
    if (
        _argv_value(argv, ("-t", "--threads")) != "12"
        or _argv_value(argv, ("-q", "--mapq")) != "20"
        or not _argv_contains(argv, "--tagSupplementary")
        or not _argv_contains(argv, "--output-somatic-vcf")
    ):
        raise ContractError(
            "E_PRODUCER_OPTION_MISMATCH",
            "producer argv does not encode explicit threads=12, MAPQ=20, supplementary tags, and somatic VCF output",
            exit_code=4,
            stage="producer_scope",
            sample=sample,
        )
    argv_sha256 = canonical_sha256(argv)
    if producer.get("command_argv_sha256") != argv_sha256:
        raise ContractError(
            "E_PRODUCER_COMMAND_DIGEST", "producer argv digest mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    input_bindings = producer.get("input_bindings")
    if not isinstance(input_bindings, dict) or producer.get("input_binding_sha256") != canonical_sha256(input_bindings):
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "producer input-binding digest mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    expected_inputs = {
        "tumor_bam": {
            "path": expected_paths["tumor_bam"],
            "storage_identity_v1_sha256": expected_storage_sha256,
        },
        "caller_raw_vcf": {"path": expected_paths["caller_raw_vcf"], "sha256": expected_hashes["caller_raw_vcf"]},
        "caller_pass_vcf": {"path": expected_paths["caller_pass_vcf"], "sha256": expected_hashes["caller_pass_vcf"]},
    }
    if input_bindings != expected_inputs:
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "producer inputs do not bind the consumed BAM/VCFs",
            exit_code=4, stage="producer_scope", sample=sample,
            expected=expected_inputs, observed=input_bindings,
        )
    artifacts = validation.get("artifacts")
    expected_artifacts = {
        "sidecar_sha256": expected_hashes["sidecar"],
        "sidecar_index_sha256": expected_hashes["sidecar_index"],
    }
    if artifacts != expected_artifacts:
        raise ContractError(
            "E_SIDECAR_SUBJECT_MISMATCH", "capture validation does not bind sidecar/index bytes",
            exit_code=5, stage="sidecar_join", sample=sample,
            expected=expected_artifacts, observed=artifacts,
        )
    checks = validation.get("checks", {})
    if checks.get("truth_flags_absent") is not True or checks.get("sidecar_no_exact_identity_conflicts") is not True:
        raise ContractError(
            "E_SIDECAR_VALIDATION", "capture validation required checks are not true",
            exit_code=5, stage="sidecar_join", sample=sample,
        )
    return {
        "command_argv_sha256": argv_sha256,
        "input_binding_sha256": canonical_sha256(input_bindings),
        "effective_options_sha256": canonical_sha256(effective_options),
        "mapped_alignment_count": mapped_count,
        "identity_unique_count": unique_count,
        "duplicate_count": duplicate_count,
        "conflict_count": conflict_count,
    }


def _validate_native_sidecar_document(validation: dict[str, Any], sample: str) -> None:
    if validation.get("sample") not in (None, sample) or validation.get("region") != "all" or validation.get("pass") is not True:
        raise ContractError(
            "E_SIDECAR_VALIDATION", "native validation is not an all-region PASS for this sample",
            exit_code=5, stage="sidecar_join", sample=sample,
        )
    if "duplicate_exact_alignment_rows" not in validation:
        raise ContractError(
            "E_SIDECAR_VALIDATION", "native validation lacks duplicate identity count",
            exit_code=5, stage="sidecar_join", sample=sample,
        )
    redundant = validation["duplicate_exact_alignment_rows"]
    conflicts = validation.get("duplicate_exact_alignment_conflicts")
    capture_rows = validation.get("capture", {}).get("rows_mapped")
    if (
        isinstance(redundant, bool)
        or not isinstance(redundant, int)
        or redundant < 0
        or isinstance(capture_rows, bool)
        or not isinstance(capture_rows, int)
        or capture_rows <= redundant
    ):
        raise ContractError(
            "E_JOIN_MULTIPLICITY", "native validation reports invalid redundant-row conservation",
            exit_code=5, stage="sidecar_join", sample=sample,
        )
    if conflicts != 0:
        raise ContractError(
            "E_JOIN_CONFLICT", "native validation reports conflicting HP/PS for one identity",
            exit_code=5, stage="sidecar_join", sample=sample,
        )
    required_checks = (
        "truth_flags_absent", "parser_count_matches_input", "sidecar_row_count_matches_capture",
        "tagged_count_matches_execution", "sidecar_no_unknown_HP",
        "sidecar_no_exact_identity_conflicts", "recalibrated_preserves_all_input_keys",
    )
    failed = [name for name in required_checks if validation.get("checks", {}).get(name) is not True]
    if failed:
        raise ContractError(
            "E_SIDECAR_VALIDATION", f"native validation failed/missed checks: {failed}",
            exit_code=5, stage="sidecar_join", sample=sample,
        )


def _verify_storage_spec(spec: dict[str, Any], sample: str, name: str) -> dict[str, Any]:
    path = Path(spec["path"])
    index_path = Path(spec["index_path"])
    observed = storage_identity_v1(path, index_path)
    _compare(
        spec["storage_identity_v1"], observed, code="E_STORAGE_IDENTITY",
        message=f"{name} storage identity changed", exit_code=3, stage="artifact_identity",
        sample=sample, artifact_name=str(path),
    )
    return {"path": str(path), "index_path": str(index_path), "storage_identity_v1": observed}


def _validate_producer_receipt(
    receipt: dict[str, Any], sample: str, receipt_schema_path: Path, expected: dict[str, Any]
) -> dict[str, Any]:
    apply_json_schema(receipt, receipt_schema_path)
    if receipt["schema_name"] != PRODUCER_RECEIPT_SCHEMA or receipt["schema_version"] != PRODUCER_RECEIPT_VERSION:
        raise ContractError(
            "E_PRODUCER_RECEIPT_UNSUPPORTED", "producer receipt schema/version is unsupported",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if (
        receipt["sample"] != sample
        or receipt["evidence_origin"] != "post_run_normalization_from_frozen_execution_artifacts"
        or receipt["identity_schema"] != "coordinate_join_v1"
        or receipt["assurance"] != "bounded_coordinate_equivalence_not_global_payload_identity"
    ):
        raise ContractError(
            "E_PRODUCER_RECEIPT_SCOPE", "producer receipt sample/origin/assurance mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if receipt["production_policy"] != PRODUCTION_POLICY:
        raise ContractError(
            "E_PRODUCTION_FILTER_POLICY", "producer receipt is not production_default_filter",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    argv = receipt["command_argv"]
    if receipt["command_argv_sha256"] != canonical_sha256(argv):
        raise ContractError(
            "E_PRODUCER_COMMAND_DIGEST", "producer argv digest mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if any(_argv_contains(argv, option) for option in ("--truth-vcf", "--truth-bed")):
        raise ContractError(
            "E_TRUTH_FLAG_PRESENT", "producer argv contains forbidden truth input",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if _argv_contains(argv, "--disableFilter"):
        raise ContractError(
            "E_PRODUCTION_FILTER_POLICY", "producer argv disables LongPhase internal filtering",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    effective_options = receipt["effective_options"]
    if effective_options != PRODUCTION_EFFECTIVE_OPTIONS:
        raise ContractError(
            "E_PRODUCER_OPTION_MISMATCH", "producer effective option values/origins mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if (
        _argv_value(argv, ("-t", "--threads")) != "12"
        or _argv_value(argv, ("-q", "--mapq")) != "20"
        or not _argv_contains(argv, "--tagSupplementary")
        or not _argv_contains(argv, "--output-somatic-vcf")
    ):
        raise ContractError(
            "E_PRODUCER_OPTION_MISMATCH", "producer explicit options do not match the production profile",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    inputs = receipt["producer_inputs"]
    if receipt["producer_input_binding_sha256"] != canonical_sha256(inputs):
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "producer input-binding digest mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    observed_inputs = {
        "germline_phased_vcf": verify_indexed_artifact(
            inputs["germline_phased_vcf"], sample, "producer.germline_phased_vcf"
        ),
        "normal_bam": _verify_storage_spec(inputs["normal_bam"], sample, "producer.normal_bam"),
        "tumor_bam": _verify_storage_spec(inputs["tumor_bam"], sample, "producer.tumor_bam"),
        "caller_pass_vcf": verify_indexed_artifact(
            inputs["caller_pass_vcf"], sample, "producer.caller_pass_vcf"
        ),
        "reference": {
            "path": inputs["reference"]["path"],
            "fai_path": inputs["reference"]["fai_path"],
            "storage_identity_v1": _verify_storage_spec(
                {
                    "path": inputs["reference"]["path"],
                    "index_path": inputs["reference"]["fai_path"],
                    "storage_identity_v1": inputs["reference"]["storage_identity_v1"],
                },
                sample,
                "producer.reference",
            )["storage_identity_v1"],
            "dict": verify_artifact(inputs["reference"]["dict"], sample, "producer.reference.dict"),
        },
    }
    if observed_inputs != inputs:
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "producer input identities changed",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if inputs["tumor_bam"] != expected["tumor_bam"] or inputs["caller_pass_vcf"] != expected["caller_pass_vcf"]:
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "producer receipt does not bind manifest tumor BAM/caller PASS",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    required_paths = [
        inputs["germline_phased_vcf"]["path"], inputs["normal_bam"]["path"], inputs["tumor_bam"]["path"],
        inputs["caller_pass_vcf"]["path"], inputs["reference"]["path"],
    ]
    if any(path not in argv for path in required_paths):
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "not every hashed producer input appears as an exact argv token",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    raw = verify_indexed_artifact(
        receipt["layered_ledger_artifacts"]["caller_raw_vcf"], sample, "ledger.caller_raw_vcf"
    )
    if raw != expected["caller_raw_vcf"]:
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "caller raw VCF ledger binding mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    outputs = {
        "sidecar": verify_artifact(receipt["capture_outputs"]["sidecar"], sample, "capture.sidecar"),
        "sidecar_index": verify_artifact(
            receipt["capture_outputs"]["sidecar_index"], sample, "capture.sidecar_index"
        ),
        "native_validation": verify_artifact(
            receipt["capture_outputs"]["native_validation"], sample, "capture.native_validation"
        ),
        "longphase_recalibrated_all_vcf": verify_indexed_artifact(
            receipt["capture_outputs"]["longphase_recalibrated_all_vcf"], sample, "capture.recalibrated_vcf"
        ),
        "longphase_recalibrated_pass_vcf": verify_indexed_artifact(
            receipt["capture_outputs"]["longphase_recalibrated_pass_vcf"], sample, "capture.recalibrated_pass_vcf"
        ),
    }
    if outputs != expected["capture_outputs"]:
        raise ContractError(
            "E_SIDECAR_SUBJECT_MISMATCH", "producer receipt outputs do not bind consumed artifacts",
            exit_code=5, stage="sidecar_join", sample=sample, expected=expected["capture_outputs"], observed=outputs,
        )
    executable = {
        **verify_artifact(receipt["longphase_executable"], sample, "longphase_executable"),
        "version": receipt["longphase_executable"]["version"],
    }
    if not executable["version"]:
        raise ContractError(
            "E_PRODUCER_AMBIGUOUS", "LongPhase version is empty", exit_code=4, stage="producer_scope", sample=sample
        )
    code = {
        key: verify_artifact(value, sample, f"producer_code.{key}")
        for key, value in receipt["producer_code"].items()
    }
    environment_spec = receipt["environment_lock"]
    environment = {
        "production_manifest": verify_artifact(
            environment_spec["production_manifest"], sample, "environment.production_manifest"
        ),
        "run_params": verify_artifact(environment_spec["run_params"], sample, "environment.run_params"),
        "code_hash_manifest": verify_artifact(
            environment_spec["code_hash_manifest"], sample, "environment.code_hash_manifest"
        ),
        "sample_input_inventory": verify_artifact(
            environment_spec["sample_input_inventory"], sample, "environment.sample_input_inventory"
        ),
        "sample_input_hash_manifest": verify_artifact(
            environment_spec["sample_input_hash_manifest"], sample, "environment.sample_input_hash_manifest"
        ),
        "sample_output_hash_manifest": verify_artifact(
            environment_spec["sample_output_hash_manifest"], sample, "environment.sample_output_hash_manifest"
        ),
        "python_executable": verify_artifact(
            environment_spec["python_executable"], sample, "environment.python_executable"
        ),
        "python_version": environment_spec["python_version"],
        "platform": environment_spec["platform"],
        "normalizer_source": verify_artifact(
            environment_spec["normalizer_source"], sample, "environment.normalizer_source"
        ),
        "normalizer_argv": environment_spec["normalizer_argv"],
    }
    counts = receipt["global_coordinate_counts"]
    if (
        counts["mapped_alignment_count"] <= 0
        or counts["identity_unique_count"] != counts["mapped_alignment_count"]
        or counts["duplicate_count"] != 0
        or counts["conflict_count"] != 0
    ):
        raise ContractError(
            "E_JOIN_MULTIPLICITY", "producer receipt global counts are not one-to-one",
            exit_code=5, stage="sidecar_join", sample=sample,
        )
    return {
        "command_argv_sha256": canonical_sha256(argv),
        "input_binding_sha256": canonical_sha256(inputs),
        "effective_options_sha256": canonical_sha256(effective_options),
        "producer_receipt_sha256": expected["producer_receipt_sha256"],
        "mapped_alignment_count": counts["mapped_alignment_count"],
        "identity_unique_count": counts["identity_unique_count"],
        "duplicate_count": counts["duplicate_count"],
        "conflict_count": counts["conflict_count"],
        "producer_inputs": observed_inputs,
        "longphase_executable": executable,
        "producer_code": code,
        "environment_lock": environment,
    }


def _validate_raw_all_producer_receipt(
    receipt: dict[str, Any], sample: str, receipt_schema_path: Path, expected: dict[str, Any]
) -> dict[str, Any]:
    apply_json_schema(receipt, receipt_schema_path)
    if (
        receipt.get("schema_name") != RAW_ALL_PRODUCER_RECEIPT_SCHEMA
        or receipt.get("schema_version") != RAW_ALL_PRODUCER_RECEIPT_VERSION
        or receipt.get("sample") != sample
        or receipt.get("evidence_origin")
        != "post_run_normalization_from_frozen_raw_all_execution_artifacts"
        or receipt.get("identity_schema") != "coordinate_join_v1"
        or receipt.get("assurance") != "bounded_coordinate_equivalence_not_global_payload_identity"
        or receipt.get("longphase_input_contract") != RAW_ALL_INPUT_CONTRACT
        or receipt.get("tree_input_contract") != TREE_INPUT_CONTRACT
        or receipt.get("duplicate_identity_policy") != REDUNDANT_IDENTITY_POLICY
    ):
        raise ContractError(
            "E_PRODUCER_RECEIPT_SCOPE",
            "raw-all producer receipt sample/origin/contracts mismatch",
            exit_code=4,
            stage="producer_scope",
            sample=sample,
        )
    if receipt["production_policy"] != PRODUCTION_POLICY:
        raise ContractError(
            "E_PRODUCTION_FILTER_POLICY", "raw-all receipt production policy mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if receipt["effective_options"] != PRODUCTION_EFFECTIVE_OPTIONS:
        raise ContractError(
            "E_PRODUCER_OPTION_MISMATCH", "raw-all receipt effective options mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    argv = receipt["command_argv"]
    if receipt["command_argv_sha256"] != canonical_sha256(argv):
        raise ContractError(
            "E_PRODUCER_COMMAND_DIGEST", "raw-all producer argv digest mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if any(_argv_contains(argv, option) for option in ("--truth-vcf", "--truth-bed")):
        raise ContractError(
            "E_TRUTH_FLAG_PRESENT", "raw-all producer argv contains forbidden truth input",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if _argv_contains(argv, "--disableFilter"):
        raise ContractError(
            "E_PRODUCTION_FILTER_POLICY", "raw-all producer disables internal filtering",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if (
        _argv_value(argv, ("-t", "--threads")) != "12"
        or _argv_value(argv, ("-q", "--mapq")) != "20"
        or not _argv_contains(argv, "--tagSupplementary")
        or not _argv_contains(argv, "--output-somatic-vcf")
    ):
        raise ContractError(
            "E_PRODUCER_OPTION_MISMATCH", "raw-all producer explicit options mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )

    inputs = receipt["producer_inputs"]
    if receipt["producer_input_binding_sha256"] != canonical_sha256(inputs):
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "raw-all producer input digest mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    observed_inputs = {
        "germline_phased_vcf": verify_indexed_artifact(
            inputs["germline_phased_vcf"], sample, "producer.germline_phased_vcf"
        ),
        "normal_bam": _verify_storage_spec(inputs["normal_bam"], sample, "producer.normal_bam"),
        "tumor_bam": _verify_storage_spec(inputs["tumor_bam"], sample, "producer.tumor_bam"),
        "caller_raw_vcf": verify_indexed_artifact(
            inputs["caller_raw_vcf"], sample, "producer.caller_raw_vcf"
        ),
        "longphase_input_vcf": verify_indexed_artifact(
            inputs["longphase_input_vcf"], sample, "producer.longphase_input_vcf"
        ),
        "caller_pass_baseline_vcf": verify_indexed_artifact(
            inputs["caller_pass_baseline_vcf"], sample, "producer.caller_pass_baseline_vcf"
        ),
        "reference": {
            "path": inputs["reference"]["path"],
            "fai_path": inputs["reference"]["fai_path"],
            "storage_identity_v1": _verify_storage_spec(
                {
                    "path": inputs["reference"]["path"],
                    "index_path": inputs["reference"]["fai_path"],
                    "storage_identity_v1": inputs["reference"]["storage_identity_v1"],
                },
                sample,
                "producer.reference",
            )["storage_identity_v1"],
        },
    }
    if observed_inputs != inputs:
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "raw-all producer input identities changed",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    for role in ("tumor_bam", "caller_raw_vcf", "longphase_input_vcf", "caller_pass_baseline_vcf"):
        if inputs[role] != expected[role]:
            raise ContractError(
                "E_INPUT_BINDING_MISMATCH", f"raw-all receipt does not bind manifest role {role}",
                exit_code=4, stage="producer_scope", sample=sample,
            )
    if _argv_value(argv, ("--tumor-snv-file",)) != inputs["longphase_input_vcf"]["path"]:
        raise ContractError(
            "E_LONGPHASE_INPUT_CONTRACT", "LongPhase-S was not fed the normalized raw-all VCF",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    required_argv_paths = (
        inputs["germline_phased_vcf"]["path"], inputs["normal_bam"]["path"],
        inputs["tumor_bam"]["path"], inputs["longphase_input_vcf"]["path"], inputs["reference"]["path"],
    )
    if any(path not in argv for path in required_argv_paths):
        raise ContractError(
            "E_INPUT_BINDING_MISMATCH", "not every consumed raw-all input appears in argv",
            exit_code=4, stage="producer_scope", sample=sample,
        )

    output_specs = receipt["capture_outputs"]
    outputs = {
        "sidecar": verify_artifact(output_specs["sidecar"], sample, "capture.sidecar"),
        "sidecar_index": verify_artifact(output_specs["sidecar_index"], sample, "capture.sidecar_index"),
        "native_validation": verify_artifact(
            output_specs["native_validation"], sample, "capture.native_validation"
        ),
        "stream_capture_summary": verify_artifact(
            output_specs["stream_capture_summary"], sample, "capture.stream_capture_summary"
        ),
        "normalization_audit": verify_artifact(
            output_specs["normalization_audit"], sample, "capture.normalization_audit"
        ),
        "filter_transition_audit": verify_artifact(
            output_specs["filter_transition_audit"], sample, "capture.filter_transition_audit"
        ),
        "sample_verification": verify_artifact(
            output_specs["sample_verification"], sample, "capture.sample_verification"
        ),
        "longphase_recalibrated_all_vcf": verify_indexed_artifact(
            output_specs["longphase_recalibrated_all_vcf"], sample, "capture.recalibrated_all_vcf"
        ),
        "longphase_recalibrated_pass_vcf": verify_indexed_artifact(
            output_specs["longphase_recalibrated_pass_vcf"], sample, "capture.recalibrated_pass_vcf"
        ),
    }
    if outputs != expected["capture_outputs"]:
        raise ContractError(
            "E_SIDECAR_SUBJECT_MISMATCH", "raw-all receipt outputs do not bind consumed artifacts",
            exit_code=5, stage="sidecar_join", sample=sample,
        )

    executable = {
        **verify_artifact(receipt["longphase_executable"], sample, "longphase_executable"),
        "version": receipt["longphase_executable"]["version"],
    }
    patch = receipt["patch_evidence"]
    patch_evidence = {
        "approval_scope": patch["approval_scope"],
        "patch_receipt": verify_artifact(patch["patch_receipt"], sample, "patch_receipt"),
        "source_patch": verify_artifact(patch["source_patch"], sample, "source_patch"),
    }
    patch_document = strict_json_load(Path(patch_evidence["patch_receipt"]["path"]))
    if (
        patch_evidence["approval_scope"] != "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY"
        or patch_document.get("status") != "APPROVED_FOR_FAIL_CLOSED_7_DATASET_VALIDATION"
        or patch_document.get("approval", {}).get("scope") != patch_evidence["approval_scope"]
        or patch_document.get("patched_binary_sha256") != executable["identity"]["sha256"]
    ):
        raise ContractError(
            "E_PATCH_EVIDENCE", "raw-all fail-closed patch evidence mismatch",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    code = {
        role: verify_artifact(spec, sample, f"producer_code.{role}")
        for role, spec in receipt["producer_code"].items()
    }
    environment_spec = receipt["environment_lock"]
    environment = {
        role: verify_artifact(environment_spec[role], sample, f"environment.{role}")
        for role in (
            "production_manifest", "run_params", "code_hash_manifest", "sample_input_inventory",
            "sample_input_hash_manifest", "sample_output_hash_manifest", "python_executable", "normalizer_source",
        )
    }
    environment.update({
        "python_version": environment_spec["python_version"],
        "platform": environment_spec["platform"],
        "normalizer_argv": environment_spec["normalizer_argv"],
    })

    counts = receipt["global_coordinate_counts"]
    validate_redundant_coordinate_counts(counts, sample)
    transition = receipt["filter_transition_summary"]
    if (
        transition["pass"] is not True
        or transition["input_record_count"] != transition["output_record_count"]
        or sum(transition["transition_counts"].values()) != transition["input_record_count"]
    ):
        raise ContractError(
            "E_FILTER_TRANSITION_AUDIT", "raw-all FILTER transition counts do not conserve records",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    normalization_doc = strict_json_load(Path(outputs["normalization_audit"]["path"]))
    transition_doc = strict_json_load(Path(outputs["filter_transition_audit"]["path"]))
    sample_doc = strict_json_load(Path(outputs["sample_verification"]["path"]))
    expected_transition = {
        "input_record_count": transition_doc.get("input", {}).get("record_count"),
        "output_record_count": transition_doc.get("output", {}).get("record_count"),
        "rescued_nonpass_to_pass": transition_doc.get("rescued_nonpass_to_pass"),
        "removed_pass_to_nonpass": transition_doc.get("removed_pass_to_nonpass"),
        "transition_counts": transition_doc.get("filter_transition_counts"),
        "pass": transition_doc.get("pass"),
    }
    if transition != expected_transition:
        raise ContractError(
            "E_FILTER_TRANSITION_AUDIT", "raw-all receipt summary differs from audit artifact",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    if (
        normalization_doc.get("pass") is not True
        or normalization_doc.get("input", {}).get("record_count") != transition["input_record_count"]
        or normalization_doc.get("output", {}).get("record_count") != transition["input_record_count"]
        or normalization_doc.get("rescued_nonpass_to_pass") != 0
        or normalization_doc.get("removed_pass_to_nonpass") != 0
        or sample_doc.get("sample") != sample
        or sample_doc.get("pass") is not True
        or sample_doc.get("normalization") != normalization_doc
        or sample_doc.get("filter_transitions") != transition_doc
    ):
        raise ContractError(
            "E_SAMPLE_VERIFICATION", "raw-all normalization/sample evidence is inconsistent",
            exit_code=4, stage="producer_scope", sample=sample,
        )

    bam_policy = receipt["bam_output_policy"]
    fifo = Path(bam_policy["consumed_fifo_path"])
    try:
        fifo_mode = fifo.stat().st_mode
    except OSError as error:
        raise ContractError(
            "E_BAM_OUTPUT_POLICY", f"cannot stat consumed tagged-BAM FIFO: {fifo}",
            exit_code=4, stage="producer_scope", sample=sample,
        ) from error
    regular_bams = sorted(path for path in fifo.parent.glob("*.bam") if path.is_file())
    if (
        bam_policy != {
            "transport": "named_fifo",
            "persisted_bam": False,
            "regular_bam_count": 0,
            "consumed_fifo_path": str(fifo),
            "is_fifo_at_closeout": True,
        }
        or not stat.S_ISFIFO(fifo_mode)
        or regular_bams
    ):
        raise ContractError(
            "E_BAM_OUTPUT_POLICY", "raw-all run persisted a BAM or lost its FIFO evidence",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    return {
        "command_argv_sha256": canonical_sha256(argv),
        "input_binding_sha256": canonical_sha256(inputs),
        "effective_options_sha256": canonical_sha256(receipt["effective_options"]),
        "producer_receipt_sha256": expected["producer_receipt_sha256"],
        "mapped_alignment_count": counts["mapped_alignment_count"],
        "identity_unique_count": counts["identity_unique_count"],
        "duplicate_count": counts["duplicate_count"],
        "conflict_count": counts["conflict_count"],
        "producer_inputs": observed_inputs,
        "longphase_executable": executable,
        "patch_evidence": patch_evidence,
        "producer_code": code,
        "environment_lock": environment,
        "filter_transition_summary": transition,
        "bam_output_policy": bam_policy,
    }


def set_method_receipt_repo_root(root: Path | None) -> None:
    global METHOD_RECEIPT_REPO_ROOT
    if root is not None and (not root.is_absolute() or root.name != "InterSubMod"):
        raise ContractError(
            "E_JOIN_METHOD_RECEIPT", f"invalid explicit receipt repo root: {root}",
            exit_code=5, stage="sidecar_join",
        )
    METHOD_RECEIPT_REPO_ROOT = root


def _resolve_receipt_reference(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    if path.parts and path.parts[0] == "InterSubMod":
        repo_root = METHOD_RECEIPT_REPO_ROOT or Path(__file__).resolve().parents[5]
        return repo_root.parent / path
    repo_root = METHOD_RECEIPT_REPO_ROOT or Path(__file__).resolve().parents[5]
    return repo_root / path


def _validate_method_receipt(document: dict[str, Any], probe_kind: str) -> None:
    if document.get("schema_version") != "1.0" or document.get("pass") is not True:
        raise ContractError(
            "E_JOIN_METHOD_RECEIPT", f"invalid {probe_kind} join-method receipt",
            exit_code=5, stage="sidecar_join",
        )
    if probe_kind == "real_data_bounded":
        comparisons = document.get("comparisons", {})
        expected_comparisons = {
            "payload_equal",
            "payload_sha256_equal",
            "calls_sha256_equal",
            "hp_sha256_equal",
            "ps_sha256_equal",
            "sidecar_exact_matches_all_exposures",
            "sidecar_missing_zero",
            "sidecar_conflicts_zero",
        }
        sidecar = document.get("sidecar", {})
        direct_sha = document.get("direct", {}).get("payload_sha256")
        joined_sha = document.get("joined", {}).get("payload_sha256")
        if (
            document.get("scope") != "PARTIAL_BOUNDED_REAL_DATA_PROBE"
            or document.get("claim_limit") != "join equivalence only; historical tags are not clean production-tag evidence"
            or document.get("identity_schema") != "coordinate_join_v1"
            or document.get("region", {}).get("chrom") != "chr1"
            or not isinstance(sidecar.get("counts", {}).get("rows"), int)
            or sidecar["counts"]["rows"] <= 0
            or set(comparisons) != expected_comparisons
            or any(comparisons[key] is not True for key in expected_comparisons)
            or not isinstance(direct_sha, str)
            or direct_sha != joined_sha
        ):
            raise ContractError(
                "E_JOIN_METHOD_RECEIPT", "COLO829 bounded equivalence checks are incomplete",
                exit_code=5, stage="sidecar_join",
            )
        for path_key, hash_key in (("path", "sha256"), ("index_path", "index_sha256")):
            referenced = _resolve_receipt_reference(sidecar[path_key])
            ensure_regular_file(referenced, f"method_receipt.sidecar.{path_key}")
            if file_sha256(referenced) != sidecar[hash_key]:
                raise ContractError(
                    "E_JOIN_METHOD_RECEIPT_HASH", f"bounded receipt {path_key} hash changed",
                    exit_code=5, stage="sidecar_join", artifact=str(referenced),
                )
    else:
        expected_tests = {
            "test_exact_sidecar_matches_direct_tagged_bam",
            "test_missing_alignment_is_detected",
            "test_duplicate_identity_conflict_is_detected",
            "test_redundant_identical_sidecar_row_collapses_without_conflict",
            "test_duplicate_bam_identity_with_conflicting_allele_is_excluded",
        }
        tests = document.get("tests")
        if (
            document.get("scope") != "SYNTHETIC_COORDINATE_JOIN_V1_CONTRACT"
            or document.get("claim_limit") != "join edge-case behavior only; not production-data validation"
            or not isinstance(tests, list)
            or len(tests) != 5
            or {case.get("name") for case in tests if isinstance(case, dict)} != expected_tests
            or any(case.get("pass") is not True for case in tests if isinstance(case, dict))
            or document.get("tests_run") != 5
            or any(document.get(key) != 0 for key in ("failures", "errors", "exit_code"))
        ):
            raise ContractError(
                "E_JOIN_METHOD_RECEIPT", "synthetic receipt must contain the exact five PASS cases",
                exit_code=5, stage="sidecar_join",
            )
        for key in ("test_source", "consumer_source"):
            source = document.get(key, {})
            referenced = _resolve_receipt_reference(source.get("path", ""))
            ensure_regular_file(referenced, f"method_receipt.{key}")
            if file_sha256(referenced) != source.get("sha256"):
                raise ContractError(
                    "E_JOIN_METHOD_RECEIPT_HASH", f"synthetic receipt {key} hash changed",
                    exit_code=5, stage="sidecar_join", artifact=str(referenced),
                )


def verify_join_method(contract: dict[str, Any]) -> dict[str, Any]:
    method = contract["join_method_validation"]
    observed: dict[str, Any] = {
        "assurance": method["assurance"],
        "claim_limit": method["claim_limit"],
    }
    for key, kind in (
        ("real_data_bounded_receipt", "real_data_bounded"),
        ("synthetic_three_case_receipt", "synthetic_three_case"),
    ):
        verified = verify_artifact(method[key], "_GLOBAL_", key)
        _validate_method_receipt(strict_json_load(Path(verified["path"])), kind)
        observed[key] = verified
    return observed


def verify_sample(meta: dict[str, Any]) -> dict[str, Any]:
    sample = meta["sample"]
    alignment = meta["alignment_payload"]
    bam = Path(alignment["path"])
    bam_index = Path(alignment["index_path"])
    observed_storage = storage_identity_v1(bam, bam_index)
    _compare(
        alignment["storage_identity_v1"], observed_storage,
        code="E_STORAGE_IDENTITY", message="BAM storage_identity_v1 changed",
        exit_code=3, stage="artifact_identity", sample=sample, artifact_name=str(bam),
    )
    observed_somatic = {
        name: verify_indexed_artifact(spec, sample, name)
        for name, spec in meta["somatic"].items()
        if name != "backbone_role"
    }
    read_tags = meta["read_tags"]
    observed_read_tags = {
        "sidecar": verify_artifact(read_tags["sidecar"], sample, "read_tags.sidecar"),
        "index": verify_artifact(read_tags["index"], sample, "read_tags.index"),
        "validation": verify_artifact(read_tags["validation"], sample, "read_tags.validation"),
        "producer_capture_receipt": verify_artifact(
            read_tags["producer_capture_receipt"], sample, "read_tags.producer_capture_receipt"
        ),
    }
    observed_coordinate_counts = inspect_coordinate_sidecar(Path(read_tags["sidecar"]["path"]), sample)
    validation = strict_json_load(Path(read_tags["validation"]["path"]))
    _validate_native_sidecar_document(validation, sample)
    receipt = strict_json_load(Path(read_tags["producer_capture_receipt"]["path"]))
    expected_receipt = {
        "tumor_bam": {
            "path": str(bam), "index_path": str(bam_index), "storage_identity_v1": observed_storage,
        },
        "caller_raw_vcf": observed_somatic["caller_raw_vcf"],
        "longphase_input_vcf": observed_somatic["longphase_input_vcf"],
        "caller_pass_baseline_vcf": observed_somatic["caller_pass_baseline_vcf"],
        "capture_outputs": {
            "sidecar": observed_read_tags["sidecar"],
            "sidecar_index": observed_read_tags["index"],
            "native_validation": observed_read_tags["validation"],
            "stream_capture_summary": verify_artifact(
                receipt["capture_outputs"]["stream_capture_summary"], sample, "capture.stream_capture_summary"
            ),
            "normalization_audit": verify_artifact(
                receipt["capture_outputs"]["normalization_audit"], sample, "capture.normalization_audit"
            ),
            "filter_transition_audit": verify_artifact(
                receipt["capture_outputs"]["filter_transition_audit"], sample, "capture.filter_transition_audit"
            ),
            "sample_verification": verify_artifact(
                receipt["capture_outputs"]["sample_verification"], sample, "capture.sample_verification"
            ),
            "longphase_recalibrated_all_vcf": observed_somatic["longphase_recalibrated_all_vcf"],
            "longphase_recalibrated_pass_vcf": observed_somatic["longphase_recalibrated_pass_vcf"],
        },
        "producer_receipt_sha256": observed_read_tags["producer_capture_receipt"]["identity"]["sha256"],
    }
    validation_evidence = _validate_raw_all_producer_receipt(
        receipt,
        sample,
        Path(__file__).resolve().parent.parent / "schemas" / "longphase_raw_all_capture_receipt_v2.schema.json",
        expected_receipt,
    )
    for key in ("mapped_alignment_count", "identity_unique_count", "duplicate_count", "conflict_count"):
        if validation_evidence[key] != observed_coordinate_counts[key]:
            raise ContractError(
                "E_SIDECAR_SUBJECT_MISMATCH",
                f"capture validation {key} differs from sidecar scan",
                exit_code=5,
                stage="sidecar_join",
                sample=sample,
                expected=observed_coordinate_counts[key],
                observed=validation_evidence[key],
            )
    if read_tags["producer_policy"] != PRODUCTION_POLICY:
        raise ContractError(
            "E_PRODUCTION_FILTER_POLICY", "manifest producer policy differs from required policy",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    binding = read_tags["subject_binding"]
    expected_binding = {
        "schema_name": SIDECAR_SUBJECT_SCHEMA,
        "schema_version": SIDECAR_SUBJECT_VERSION,
        "sample": sample,
        "duplicate_identity_policy": REDUNDANT_IDENTITY_POLICY,
        "coordinate_identity_columns": COORDINATE_IDENTITY_COLUMNS,
        "sidecar_sha256": observed_read_tags["sidecar"]["identity"]["sha256"],
        "sidecar_index_sha256": observed_read_tags["index"]["identity"]["sha256"],
        "validation_sha256": observed_read_tags["validation"]["identity"]["sha256"],
        "producer_capture_receipt_sha256": observed_read_tags["producer_capture_receipt"]["identity"]["sha256"],
        "alignment_payload_storage_identity_sha256": observed_storage["identity_sha256"],
        "producer_command_argv_sha256": validation_evidence["command_argv_sha256"],
        "producer_input_binding_sha256": validation_evidence["input_binding_sha256"],
        "producer_effective_options_sha256": validation_evidence["effective_options_sha256"],
        "caller_raw_vcf_sha256": observed_somatic["caller_raw_vcf"]["identity"]["sha256"],
        "longphase_input_vcf_sha256": observed_somatic["longphase_input_vcf"]["identity"]["sha256"],
        "caller_pass_baseline_vcf_sha256": observed_somatic["caller_pass_baseline_vcf"]["identity"]["sha256"],
        "longphase_recalibrated_all_vcf_sha256": observed_somatic["longphase_recalibrated_all_vcf"]["identity"]["sha256"],
        "longphase_recalibrated_pass_vcf_sha256": observed_somatic[
            "longphase_recalibrated_pass_vcf"
        ]["identity"]["sha256"],
        "mapped_alignment_count": validation_evidence["mapped_alignment_count"],
        "identity_unique_count": validation_evidence["identity_unique_count"],
        "duplicate_count": validation_evidence["duplicate_count"],
        "conflict_count": validation_evidence["conflict_count"],
    }
    _compare(
        binding, expected_binding, code="E_SIDECAR_SUBJECT_MISMATCH",
        message="sidecar subject does not bind the consumed artifacts", exit_code=5,
        stage="sidecar_join", sample=sample, artifact_name=read_tags["sidecar"]["path"],
    )
    validate_redundant_coordinate_counts(binding, sample)
    if (
        read_tags.get("duplicate_identity_policy") != REDUNDANT_IDENTITY_POLICY
        or read_tags.get("require_unique_identity") is not False
    ):
        raise ContractError(
            "E_SIDECAR_CONTRACT_UNSUPPORTED",
            "manifest must explicitly collapse redundant identical HP/PS rows",
            exit_code=5,
            stage="sidecar_join",
            sample=sample,
        )
    copy_number: dict[str, Any] = {
        key: meta["copy_number"][key]
        for key in (
            "availability", "source", "semantics", "coordinate_system", "unlisted_position_semantics",
            "allowed_states", "overlap_policy", "reason",
        )
    }
    copy_number.update({
        key: verify_artifact(value, sample, f"copy_number.{key}") if value is not None else None
        for key, value in meta["copy_number"].items()
        if key in {"cn_bed", "cn_int_gain", "cn_int_loss", "integration_json"}
    })
    return {
        "sample": sample,
        "biological_id": meta["biological_id"],
        "pass": True,
        "alignment_payload": {"path": str(bam), "storage_identity_v1": observed_storage},
        "somatic": observed_somatic,
        "read_tags": {
            **observed_read_tags,
            "duplicate_identity_policy": REDUNDANT_IDENTITY_POLICY,
            "subject_binding": binding,
            "producer_policy": PRODUCTION_POLICY,
            "producer_evidence": validation_evidence,
        },
        "copy_number": copy_number,
    }


def validate_and_build_lock(manifest_path: Path, schema_path: Path, validator_path: Path | None = None) -> dict[str, Any]:
    manifest = strict_json_load(manifest_path)
    _pre_schema_checks(manifest)
    apply_json_schema(manifest, schema_path)
    validate_semantics(manifest)
    first_method_observation = verify_join_method(manifest["analysis_contract"])
    first_observation = [verify_sample(meta) for meta in manifest["samples"]]
    # Critical identities are observed twice before freeze to close validation-time TOCTOU.
    second_method_observation = verify_join_method(manifest["analysis_contract"])
    second_observation = [verify_sample(meta) for meta in manifest["samples"]]
    if (
        canonical_sha256(first_observation) != canonical_sha256(second_observation)
        or canonical_sha256(first_method_observation) != canonical_sha256(second_method_observation)
    ):
        raise ContractError(
            "E_PATH_RETARGETED", "critical input identities changed during validation",
            exit_code=3, stage="artifact_identity",
        )
    validator_path = validator_path or Path(__file__).resolve()
    lock: dict[str, Any] = {
        "$schema": "InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/schemas/layered_input_lock_v1.schema.json",
        "schema_name": LOCK_SCHEMA_NAME,
        "schema_version": LOCK_SCHEMA_VERSION,
        "lock_id": f"{manifest['manifest_id']}.frozen",
        "created_at_utc": utc_now(),
        "source_manifest": {
            "path": str(manifest_path.resolve()),
            "byte_sha256": file_sha256(manifest_path),
            "canonical_sha256": canonical_sha256(manifest),
        },
        "validator": {
            "path": str(validator_path),
            "sha256": file_sha256(validator_path),
            "schema_path": str(schema_path.resolve()),
            "schema_sha256": file_sha256(schema_path),
        },
        "dataset_count": 7,
        "biological_sample_count": 6,
        "analysis_contract": {**manifest["analysis_contract"], "join_method_validation": second_method_observation},
        "production_summary": manifest["production_summary"],
        "all_pass": True,
        "samples": second_observation,
    }
    lock["lock_sha256"] = canonical_sha256(lock)
    apply_json_schema(lock, schema_path.parent / "layered_input_lock_v1.schema.json")
    return lock


def failure_document(error: ContractError, manifest: Path, frozen_lock: Path) -> dict[str, Any]:
    return {
        "schema_name": "intersubmod.layered_input_failure",
        "schema_version": "1.0.0",
        "created_at_utc": utc_now(),
        "manifest": str(manifest),
        "frozen_lock_target": str(frozen_lock),
        "valid_lock_written": False,
        "error": error.record(),
    }


def main(argv: list[str] | None = None) -> int:
    script_dir = Path(__file__).resolve().parent
    default_schema = script_dir.parent / "schemas" / "layered_input_manifest_v3.schema.json"
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--schema", type=Path, default=default_schema)
    parser.add_argument("--frozen-lock", required=True, type=Path)
    parser.add_argument("--failure-report", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.frozen_lock.exists() or args.failure_report.exists():
        error = ContractError("E_OUTPUT_EXISTS", "lock/failure targets must both be new", exit_code=7, stage="lifecycle")
        print(f"{error.code}: {error}", file=sys.stderr)
        return error.exit_code
    try:
        lock = validate_and_build_lock(args.manifest, args.schema)
        atomic_publish_frozen_lock(
            args.frozen_lock, lock, args.schema.parent / "layered_input_lock_v1.schema.json"
        )
        print(f"LAYERED V3 PREFLIGHT: PASS 7/7 -> {args.frozen_lock}")
        return 0
    except ContractError as error:
        atomic_write_json(args.failure_report, failure_document(error, args.manifest, args.frozen_lock))
        print(f"{error.code}: {error} -> {args.failure_report}", file=sys.stderr)
        return error.exit_code
    except Exception as error:  # pragma: no cover - last-resort fail-closed boundary
        internal = ContractError("E_INTERNAL", repr(error), exit_code=70, stage="internal")
        atomic_write_json(args.failure_report, failure_document(internal, args.manifest, args.frozen_lock))
        print(f"{internal.code}: {internal} -> {args.failure_report}", file=sys.stderr)
        return internal.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
