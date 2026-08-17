#!/usr/bin/env python3
"""Normalize and verify the seven-dataset BAM intake manifest.

This collector is intentionally bounded.  It reads BAM headers and EOF markers,
three frozen 1 MiB BAM chunks, BAM indexes, and small producer receipts.  It does
not run flagstat, mosdepth, read-length scans, methylation scans, or variant
benchmarking.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import stat
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parents[1]
DEFAULT_SOURCE_MANIFEST = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/"
    "input_manifest.snapshot.json"
)
DEFAULT_TOPOLOGY = PROJECT_ROOT / "results" / "cohort_topology_metrics.csv"
DEFAULT_SCHEMA = PROJECT_ROOT / "multi_bam_input_manifest.schema.json"
DEFAULT_SOURCE_SCHEMA = (
    REPO_ROOT
    / "docs"
    / "methodology"
    / "_assets"
    / "20260627_subclone_4axis_teaching"
    / "schemas"
    / "layered_input_manifest_v3.schema.json"
)
CANONICAL_SOURCE_SCHEMA_ID = (
    "InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/"
    "schemas/layered_input_manifest_v3.schema.json"
)
CANONICAL_SOURCE_SCHEMA_SHA256 = (
    "47f54d86159beae66b68719e034c7755eba59b2d1719d4546153fa0ae74ff19f"
)
DEFAULT_OUTPUT = PROJECT_ROOT / "multi_bam_input_manifest.json"
DEFAULT_RECEIPT = PROJECT_ROOT / "results" / "multi_bam_input_manifest_validation.json"
FROZEN_CHUNK_SIZE_BYTES = 1024 * 1024
SUBPROCESS_TIMEOUT_SECONDS = 120
EXPECTED_DATASETS = {
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
}
EXPECTED_BIOLOGICAL_N = 6
PAYLOAD_MATCH_FIELDS = (
    "requested_path",
    "realpath",
    "logical_is_symlink",
    "logical_size_bytes",
    "logical_mtime_ns",
    "st_ino",
    "size_bytes",
    "mtime_ns",
    "ctime_ns",
    "chunk_size_bytes",
    "chunks",
    "index",
)


class ContractError(RuntimeError):
    """Raised when a source or output violates the intake contract."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-manifest", type=Path, default=DEFAULT_SOURCE_MANIFEST)
    parser.add_argument("--topology-csv", type=Path, default=DEFAULT_TOPOLOGY)
    parser.add_argument("--source-schema", type=Path, default=DEFAULT_SOURCE_SCHEMA)
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument(
        "--samtools",
        default="samtools",
        help="samtools executable used only for quickcheck and header readback",
    )
    return parser.parse_args()


def strict_json_load_with_sha256(path: Path) -> tuple[Any, str]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise ContractError(f"duplicate JSON key {key!r}: {path}")
            output[key] = value
        return output

    def reject_nonfinite(value: str) -> None:
        raise ContractError(f"non-finite JSON value {value!r}: {path}")

    try:
        raw = path.read_bytes()
        text = raw.decode("utf-8", errors="strict")
        value = json.loads(
            text,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_nonfinite,
        )
        return value, hashlib.sha256(raw).hexdigest()
    except ContractError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"cannot parse JSON {path}: {exc}") from exc


def strict_json_load(path: Path) -> Any:
    value, _ = strict_json_load_with_sha256(path)
    return value


SUPPORTED_SCHEMA_KEYWORDS = {
    "$schema",
    "$id",
    "$defs",
    "$ref",
    "title",
    "description",
    "type",
    "properties",
    "required",
    "additionalProperties",
    "const",
    "enum",
    "minItems",
    "maxItems",
    "prefixItems",
    "items",
    "uniqueItems",
    "minLength",
    "maxLength",
    "pattern",
    "minimum",
    "maximum",
    "oneOf",
    "allOf",
    "format",
}


def validate_supported_schema_keywords(schema: Any, label: str = "$") -> None:
    """Reject schema features this bounded validator cannot evaluate."""

    if isinstance(schema, bool):
        return
    if not isinstance(schema, dict):
        raise ContractError(f"JSON Schema node {label} must be an object or boolean")
    unsupported = sorted(set(schema) - SUPPORTED_SCHEMA_KEYWORDS)
    if unsupported:
        raise ContractError(f"unsupported JSON Schema keyword(s) at {label}: {unsupported}")
    for key in ("$defs", "properties"):
        mapping = schema.get(key, {})
        if not isinstance(mapping, dict):
            raise ContractError(f"JSON Schema {label}.{key} must be an object")
        for name, child in mapping.items():
            validate_supported_schema_keywords(child, f"{label}.{key}.{name}")
    for key in ("oneOf", "allOf", "prefixItems"):
        children = schema.get(key, [])
        if not isinstance(children, list):
            raise ContractError(f"JSON Schema {label}.{key} must be an array")
        for index, child in enumerate(children):
            validate_supported_schema_keywords(child, f"{label}.{key}[{index}]")
    for key in ("items", "additionalProperties"):
        if key in schema and isinstance(schema[key], dict):
            validate_supported_schema_keywords(schema[key], f"{label}.{key}")


def validate_canonical_source_schema(
    schema: Any,
    observed_sha256: str,
    *,
    path: Path | None = None,
) -> None:
    """Pin comprehensive intake to the reviewed layered-manifest v3 contract."""

    label = str(path) if path is not None else "source schema"
    if not isinstance(schema, dict):
        raise ContractError(f"canonical source JSON schema must be an object: {label}")
    if schema.get("$schema") != "https://json-schema.org/draft/2020-12/schema":
        raise ContractError(f"unexpected source JSON schema dialect: {label}")
    if schema.get("$id") != CANONICAL_SOURCE_SCHEMA_ID:
        raise ContractError(
            "canonical source schema $id mismatch: "
            f"expected {CANONICAL_SOURCE_SCHEMA_ID!r}, got {schema.get('$id')!r}"
        )
    if observed_sha256 != CANONICAL_SOURCE_SCHEMA_SHA256:
        raise ContractError(
            "canonical source schema SHA256 mismatch: "
            f"expected {CANONICAL_SOURCE_SCHEMA_SHA256}, got {observed_sha256}"
        )
    validate_supported_schema_keywords(schema)


def _json_schema_ref(root_schema: dict[str, Any], ref: str) -> Any:
    if not ref.startswith("#/"):
        raise ContractError(f"only local JSON Schema references are supported: {ref!r}")
    current: Any = root_schema
    for raw_part in ref[2:].split("/"):
        part = raw_part.replace("~1", "/").replace("~0", "~")
        if not isinstance(current, dict) or part not in current:
            raise ContractError(f"unresolvable JSON Schema reference: {ref!r}")
        current = current[part]
    return current


def _json_type_matches(value: Any, expected: str) -> bool:
    return {
        "object": isinstance(value, dict),
        "array": isinstance(value, list),
        "string": isinstance(value, str),
        "integer": isinstance(value, int) and not isinstance(value, bool),
        "number": isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value)),
        "boolean": isinstance(value, bool),
        "null": value is None,
    }.get(expected, False)


def validate_json_schema(
    value: Any,
    schema: Any,
    *,
    root_schema: dict[str, Any] | None = None,
    path: str = "$",
) -> None:
    """Evaluate the bounded Draft 2020-12 keyword subset used by both contracts."""

    if root_schema is None:
        if not isinstance(schema, dict):
            raise ContractError("root JSON Schema must be an object")
        root_schema = schema
        validate_supported_schema_keywords(schema)
    if schema is True:
        return
    if schema is False:
        raise ContractError(f"JSON Schema rejected {path}")
    if not isinstance(schema, dict):
        raise ContractError(f"invalid JSON Schema node at {path}")
    if "$ref" in schema:
        validate_json_schema(
            value,
            _json_schema_ref(root_schema, required_string(schema["$ref"], f"{path}.$ref")),
            root_schema=root_schema,
            path=path,
        )
    for child in schema.get("allOf", []):
        validate_json_schema(value, child, root_schema=root_schema, path=path)
    if "oneOf" in schema:
        matches = 0
        for child in schema["oneOf"]:
            try:
                validate_json_schema(value, child, root_schema=root_schema, path=path)
            except ContractError:
                continue
            matches += 1
        if matches != 1:
            raise ContractError(f"JSON Schema oneOf expected exactly one match at {path}; got {matches}")
    if "const" in schema and canonical_bytes(value) != canonical_bytes(schema["const"]):
        raise ContractError(f"JSON Schema const mismatch at {path}: {value!r}")
    if "enum" in schema and not any(
        canonical_bytes(value) == canonical_bytes(candidate) for candidate in schema["enum"]
    ):
        raise ContractError(f"JSON Schema enum mismatch at {path}: {value!r}")
    if "type" in schema:
        expected_types = schema["type"] if isinstance(schema["type"], list) else [schema["type"]]
        if not expected_types or not all(isinstance(item, str) for item in expected_types):
            raise ContractError(f"invalid JSON Schema type declaration at {path}")
        if not any(_json_type_matches(value, expected) for expected in expected_types):
            raise ContractError(
                f"JSON Schema type mismatch at {path}: expected {expected_types}, got {type(value).__name__}"
            )
    if isinstance(value, dict):
        required = schema.get("required", [])
        if not isinstance(required, list) or not all(isinstance(item, str) for item in required):
            raise ContractError(f"invalid JSON Schema required list at {path}")
        missing = [name for name in required if name not in value]
        if missing:
            raise ContractError(f"JSON Schema required field(s) missing at {path}: {missing}")
        properties = schema.get("properties", {})
        if not isinstance(properties, dict):
            raise ContractError(f"invalid JSON Schema properties at {path}")
        for name, child in properties.items():
            if name in value:
                validate_json_schema(
                    value[name], child, root_schema=root_schema, path=f"{path}.{name}"
                )
        extras = sorted(set(value) - set(properties))
        additional = schema.get("additionalProperties", True)
        if extras and additional is False:
            raise ContractError(f"JSON Schema additionalProperties at {path}: {extras}")
        if extras and isinstance(additional, dict):
            for name in extras:
                validate_json_schema(
                    value[name], additional, root_schema=root_schema, path=f"{path}.{name}"
                )
    if isinstance(value, list):
        if "minItems" in schema and len(value) < int(schema["minItems"]):
            raise ContractError(f"JSON Schema minItems mismatch at {path}")
        if "maxItems" in schema and len(value) > int(schema["maxItems"]):
            raise ContractError(f"JSON Schema maxItems mismatch at {path}")
        if schema.get("uniqueItems") is True:
            encoded = [canonical_bytes(item) for item in value]
            if len(encoded) != len(set(encoded)):
                raise ContractError(f"JSON Schema uniqueItems mismatch at {path}")
        prefix = schema.get("prefixItems", [])
        for index, child in enumerate(prefix[: len(value)]):
            validate_json_schema(
                value[index], child, root_schema=root_schema, path=f"{path}[{index}]"
            )
        items = schema.get("items")
        start = len(prefix) if prefix else 0
        if items is False and len(value) > start:
            raise ContractError(f"JSON Schema additional array items at {path}")
        if isinstance(items, dict):
            for index in range(start, len(value)):
                validate_json_schema(
                    value[index], items, root_schema=root_schema, path=f"{path}[{index}]"
                )
    if isinstance(value, str):
        if "minLength" in schema and len(value) < int(schema["minLength"]):
            raise ContractError(f"JSON Schema minLength mismatch at {path}")
        if "maxLength" in schema and len(value) > int(schema["maxLength"]):
            raise ContractError(f"JSON Schema maxLength mismatch at {path}")
        if "pattern" in schema and re.search(str(schema["pattern"]), value) is None:
            raise ContractError(f"JSON Schema pattern mismatch at {path}: {value!r}")
        if schema.get("format") == "date-time":
            try:
                parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
            except ValueError as exc:
                raise ContractError(f"JSON Schema date-time mismatch at {path}: {value!r}") from exc
            if parsed.tzinfo is None:
                raise ContractError(f"JSON Schema date-time lacks timezone at {path}: {value!r}")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < schema["minimum"]:
            raise ContractError(f"JSON Schema minimum mismatch at {path}: {value}")
        if "maximum" in schema and value > schema["maximum"]:
            raise ContractError(f"JSON Schema maximum mismatch at {path}: {value}")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return resolved.as_posix()


def require_regular(path: Path, label: str) -> Path:
    try:
        resolved = path.resolve(strict=True)
        mode = resolved.stat().st_mode
    except OSError as exc:
        raise ContractError(f"{label} is unavailable: {path}: {exc}") from exc
    if not stat.S_ISREG(mode):
        raise ContractError(f"{label} is not a regular file: {path} -> {resolved}")
    return resolved


def validate_frozen_storage_identity(
    frozen: dict[str, Any], *, label: str, payload_path: Path, index_path: Path
) -> None:
    if frozen.get("policy") != "storage_identity_v1":
        raise ContractError(f"{label}: frozen identity policy mismatch")
    if frozen.get("assurance") != "metadata_plus_sampled_chunks_not_full_content_hash":
        raise ContractError(f"{label}: frozen identity assurance mismatch")
    if frozen.get("is_full_content_hash") is not False:
        raise ContractError(f"{label}: bounded identity cannot claim a full-content hash")
    if frozen.get("requested_path") != str(payload_path):
        raise ContractError(f"{label}: frozen requested_path differs from declared payload path")
    size_bytes = exact_integer(frozen.get("size_bytes"), f"{label}.size_bytes", minimum=1)
    chunk_size = exact_integer(
        frozen.get("chunk_size_bytes"), f"{label}.chunk_size_bytes", minimum=1
    )
    if chunk_size != FROZEN_CHUNK_SIZE_BYTES:
        raise ContractError(
            f"{label}: chunk_size_bytes must be exactly {FROZEN_CHUNK_SIZE_BYTES}, got {chunk_size}"
        )
    if size_bytes < chunk_size:
        raise ContractError(f"{label}: payload is smaller than the fixed sampled chunk")
    chunks = frozen.get("chunks")
    if not isinstance(chunks, list) or len(chunks) != 3:
        raise ContractError(f"{label}: frozen storage identity must contain exactly three chunks")
    expected_layout = (
        ("first", 0),
        ("middle", (size_bytes - chunk_size) // 2),
        ("last", size_bytes - chunk_size),
    )
    for index, ((expected_label, expected_offset), chunk) in enumerate(
        zip(expected_layout, chunks)
    ):
        if not isinstance(chunk, dict):
            raise ContractError(f"{label}.chunks[{index}] must be an object")
        if chunk.get("label") != expected_label:
            raise ContractError(
                f"{label}.chunks[{index}].label must be {expected_label!r}, got {chunk.get('label')!r}"
            )
        if exact_integer(chunk.get("offset"), f"{label}.chunks[{index}].offset") != expected_offset:
            raise ContractError(
                f"{label}.chunks[{index}].offset does not match the fixed {expected_label} position"
            )
        length = exact_integer(chunk.get("length"), f"{label}.chunks[{index}].length", minimum=1)
        if length != chunk_size:
            raise ContractError(
                f"{label}.chunks[{index}].length must be exactly {chunk_size}, got {length}"
            )
        sha = required_string(chunk.get("sha256"), f"{label}.chunks[{index}].sha256")
        if re.fullmatch(r"[0-9a-f]{64}", sha) is None:
            raise ContractError(f"{label}.chunks[{index}].sha256 is malformed")
    frozen_index = frozen.get("index")
    if not isinstance(frozen_index, dict) or frozen_index.get("path") != str(index_path):
        raise ContractError(f"{label}: frozen index path differs from the declared index path")
    identity = frozen_index.get("identity")
    if not isinstance(identity, dict) or identity.get("policy") != "full_sha256":
        raise ContractError(f"{label}: frozen index must use full_sha256")
    required_string(identity.get("sha256"), f"{label}.index.sha256")
    exact_integer(identity.get("size_bytes"), f"{label}.index.size_bytes", minimum=1)
    expected_identity_sha = required_string(
        frozen.get("identity_sha256"), f"{label}.identity_sha256"
    )
    identity_payload = {key: value for key, value in frozen.items() if key != "identity_sha256"}
    if canonical_sha256(identity_payload) != expected_identity_sha:
        raise ContractError(f"{label}: frozen identity_sha256 does not match its canonical fields")


def sampled_chunk(path: Path, expected: dict[str, Any]) -> dict[str, Any]:
    offset = exact_integer(expected.get("offset"), f"{path}.chunk.offset", minimum=0)
    length = exact_integer(expected.get("length"), f"{path}.chunk.length", minimum=1)
    if length != FROZEN_CHUNK_SIZE_BYTES:
        raise ContractError(
            f"{path}: sampled chunk length must be exactly {FROZEN_CHUNK_SIZE_BYTES}, got {length}"
        )
    with path.open("rb") as handle:
        handle.seek(offset)
        payload = handle.read(length)
    if len(payload) != length:
        raise ContractError(
            f"short BAM chunk read: {path} offset={offset} expected={length} observed={len(payload)}"
        )
    return {
        "label": required_string(expected.get("label"), f"{path}.chunk.label"),
        "offset": offset,
        "length": length,
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def storage_identity(path: Path, index_path: Path, frozen: dict[str, Any]) -> dict[str, Any]:
    validate_frozen_storage_identity(
        frozen, label=str(path), payload_path=path, index_path=index_path
    )
    resolved = require_regular(path, "BAM")
    require_regular(index_path, "BAM index")
    logical = path.lstat()
    target = resolved.stat()
    chunks_raw = frozen.get("chunks")
    if not isinstance(chunks_raw, list) or len(chunks_raw) != 3:
        raise ContractError(f"{path}: frozen storage identity must contain three chunks")
    index_identity = {
        "policy": "full_sha256",
        "size_bytes": index_path.stat().st_size,
        "sha256": file_sha256(index_path),
    }
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
        "size_bytes": target.st_size,
        "mtime_ns": target.st_mtime_ns,
        "ctime_ns": target.st_ctime_ns,
        "chunk_size_bytes": FROZEN_CHUNK_SIZE_BYTES,
        "chunks": [sampled_chunk(resolved, item) for item in chunks_raw],
        "index": {"path": str(index_path), "identity": index_identity},
    }
    value["identity_sha256"] = canonical_sha256(value)
    return value


def exact_integer(value: Any, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{label} must be an integer, got {value!r}")
    if not math.isfinite(float(value)) or not float(value).is_integer():
        raise ContractError(f"{label} must be a finite integer, got {value!r}")
    output = int(value)
    if output < minimum:
        raise ContractError(f"{label} must be >= {minimum}, got {output}")
    return output


def required_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ContractError(f"{label} must be a non-empty string")
    return value


def ratio(numerator: int, denominator: int, label: str) -> float:
    if denominator <= 0 or numerator < 0 or numerator > denominator:
        raise ContractError(
            f"{label} requires 0 <= numerator <= positive denominator, got {numerator}/{denominator}"
        )
    return numerator / denominator


def run_command(argv: list[str], label: str) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            argv,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=SUBPROCESS_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as exc:
        raise ContractError(
            f"{label} exceeded the {SUBPROCESS_TIMEOUT_SECONDS}s timeout"
        ) from exc
    except OSError as exc:
        raise ContractError(f"cannot execute {label}: {exc}") from exc


def samtools_version(executable: str) -> str:
    result = run_command([executable, "--version"], "samtools --version")
    if result.returncode != 0 or not result.stdout.strip():
        raise ContractError(f"samtools --version failed: {result.stderr.strip()}")
    return result.stdout.splitlines()[0].strip()


def bam_quickcheck(executable: str, path: Path) -> tuple[str, str]:
    result = run_command([executable, "quickcheck", "-v", str(path)], "samtools quickcheck")
    diagnostic = "\n".join(
        part.strip() for part in (result.stdout, result.stderr) if part.strip()
    )
    return ("PASS" if result.returncode == 0 and not diagnostic else "FAIL", diagnostic)


def bam_header(executable: str, path: Path) -> dict[str, Any]:
    result = run_command([executable, "view", "-H", str(path)], "samtools view -H")
    if result.returncode != 0:
        raise ContractError(f"cannot read BAM header {path}: {result.stderr.strip()}")
    header = result.stdout
    if not header.startswith("@"):
        raise ContractError(f"BAM header is empty or malformed: {path}")
    sort_order: str | None = None
    read_group_samples: set[str] = set()
    dictionary: list[dict[str, Any]] = []
    for raw in header.splitlines():
        fields = raw.split("\t")
        tags = {
            field.split(":", 1)[0]: field.split(":", 1)[1]
            for field in fields[1:]
            if ":" in field
        }
        if fields[0] == "@HD":
            sort_order = tags.get("SO")
        elif fields[0] == "@RG" and tags.get("SM"):
            read_group_samples.add(tags["SM"])
        elif fields[0] == "@SQ":
            if not tags.get("SN") or not tags.get("LN", "").isdigit():
                raise ContractError(f"malformed @SQ line in {path}: {raw}")
            dictionary.append({"SN": tags["SN"], "LN": int(tags["LN"])})
    if not dictionary:
        raise ContractError(f"BAM header contains no @SQ dictionary: {path}")
    return {
        "sort_order": sort_order,
        "read_group_samples": sorted(read_group_samples),
        "contig_count": len(dictionary),
        "header_sha256": hashlib.sha256(header.encode("utf-8")).hexdigest(),
        "reference_dictionary_sha256": canonical_sha256(dictionary),
    }


def reference_dictionary(fai_path: Path) -> dict[str, Any]:
    require_regular(fai_path, "reference FAI")
    dictionary: list[dict[str, Any]] = []
    with fai_path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0] or not fields[1].isdigit():
                raise ContractError(f"malformed FAI row {fai_path}:{line_number}")
            dictionary.append({"SN": fields[0], "LN": int(fields[1])})
    if not dictionary:
        raise ContractError(f"reference FAI has no contigs: {fai_path}")
    return {
        "contig_count": len(dictionary),
        "reference_dictionary_sha256": canonical_sha256(dictionary),
    }


def storage_comparison(
    frozen: dict[str, Any], observed: dict[str, Any]
) -> tuple[list[str], bool, str]:
    differing_fields = sorted(
        key for key in set(frozen) | set(observed) if frozen.get(key) != observed.get(key)
    )
    payload_match = all(frozen.get(key) == observed.get(key) for key in PAYLOAD_MATCH_FIELDS)
    if not differing_fields:
        strict_status = "MATCH"
    elif set(differing_fields) == {"identity_sha256", "st_dev"} and payload_match:
        strict_status = "MOUNT_DEVICE_DRIFT_ONLY"
    else:
        strict_status = "MISMATCH"
    return differing_fields, payload_match, strict_status


def build_bam_record(
    spec: dict[str, Any], samtools: str, label: str
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    path = Path(required_string(spec.get("path"), f"{label}.path"))
    index_path = Path(required_string(spec.get("index_path"), f"{label}.index_path"))
    frozen = spec.get("storage_identity_v1")
    if not isinstance(frozen, dict):
        raise ContractError(f"{label}: missing frozen storage_identity_v1")
    observed_before = storage_identity(path, index_path, frozen)
    header = bam_header(samtools, path)
    quickcheck_status, quickcheck_diagnostic = bam_quickcheck(samtools, path)
    observed = storage_identity(path, index_path, frozen)
    if observed != observed_before:
        raise ContractError(f"{label}: BAM or index changed during bounded verification")
    differing_fields, payload_match, strict_status = storage_comparison(frozen, observed)
    index_expected_sha = required_string(
        frozen.get("index", {}).get("identity", {}).get("sha256"),
        f"{label}.frozen.index.sha256",
    )
    index_observed_sha = observed["index"]["identity"]["sha256"]
    record = {
        "path": str(path),
        "realpath": observed["realpath"],
        "index_path": str(index_path),
        "size_bytes": observed["size_bytes"],
        "logical_is_symlink": observed["logical_is_symlink"],
        **header,
        "quickcheck_status": quickcheck_status,
        "identity_policy": "storage_identity_v1",
        "identity_assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
        "frozen_identity_sha256": required_string(
            frozen.get("identity_sha256"), f"{label}.frozen.identity_sha256"
        ),
        "observed_identity_sha256": observed["identity_sha256"],
        "strict_storage_identity_status": strict_status,
        "strict_differing_fields": differing_fields,
        "bounded_payload_status": "MATCH" if payload_match else "MISMATCH",
        "index_sha256": index_observed_sha,
        "index_sha256_status": (
            "MATCH" if index_observed_sha == index_expected_sha else "MISMATCH"
        ),
        "full_content_sha256": None,
    }
    evidence = {
        "quickcheck_diagnostic": quickcheck_diagnostic,
        "frozen_st_dev": frozen.get("st_dev"),
        "observed_st_dev": observed.get("st_dev"),
        "payload_match_fields": list(PAYLOAD_MATCH_FIELDS),
    }
    return record, observed, evidence


def build_reference_record(
    spec: dict[str, Any], label: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    path = Path(required_string(spec.get("path"), f"{label}.path"))
    fai_path = Path(required_string(spec.get("fai_path"), f"{label}.fai_path"))
    frozen = spec.get("storage_identity_v1")
    if not isinstance(frozen, dict):
        raise ContractError(f"{label}: missing frozen storage_identity_v1")
    observed_before = storage_identity(path, fai_path, frozen)
    expected_fai_sha = required_string(
        frozen.get("index", {}).get("identity", {}).get("sha256"),
        f"{label}.frozen.fai.sha256",
    )
    dictionary = reference_dictionary(fai_path)
    observed = storage_identity(path, fai_path, frozen)
    if observed != observed_before:
        raise ContractError(f"{label}: reference or FAI changed during bounded verification")
    observed_fai_sha = observed["index"]["identity"]["sha256"]
    differing_fields, payload_match, strict_status = storage_comparison(frozen, observed)
    record = {
        "path": str(path),
        "realpath": observed["realpath"],
        "fai_path": str(fai_path),
        "size_bytes": observed["size_bytes"],
        "identity_policy": "storage_identity_v1",
        "identity_assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
        "strict_storage_identity_status": strict_status,
        "strict_differing_fields": differing_fields,
        "bounded_payload_status": "MATCH" if payload_match else "MISMATCH",
        "fai_sha256": observed_fai_sha,
        "fai_sha256_status": "MATCH" if observed_fai_sha == expected_fai_sha else "MISMATCH",
        **dictionary,
        "full_content_sha256": None,
    }
    evidence = {
        "frozen_st_dev": frozen.get("st_dev"),
        "observed_st_dev": observed.get("st_dev"),
        "payload_match_fields": list(PAYLOAD_MATCH_FIELDS),
    }
    return record, evidence


def producer_inputs(
    sample: str, read_tags: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    receipt_spec = read_tags.get("producer_capture_receipt")
    if not isinstance(receipt_spec, dict):
        raise ContractError(f"{sample}: missing producer_capture_receipt spec")
    receipt_path = Path(
        required_string(receipt_spec.get("path"), f"{sample}.producer_capture_receipt.path")
    )
    require_regular(receipt_path, f"{sample} producer capture receipt")
    expected_sha = required_string(
        receipt_spec.get("identity", {}).get("sha256"),
        f"{sample}.producer_capture_receipt.sha256",
    )
    receipt, observed_sha = strict_json_load_with_sha256(receipt_path)
    if observed_sha != expected_sha:
        raise ContractError(f"{sample}: producer capture receipt SHA changed")
    if not isinstance(receipt, dict) or receipt.get("sample") != sample:
        raise ContractError(f"{sample}: producer capture receipt sample mismatch")
    if (
        receipt.get("schema_name") != "intersubmod.longphase_raw_all_capture_receipt"
        or receipt.get("schema_version") != "2.0.0"
    ):
        raise ContractError(f"{sample}: producer capture receipt schema mismatch")
    if receipt.get("production_policy") != read_tags.get("producer_policy"):
        raise ContractError(f"{sample}: producer policy differs between source and receipt")
    inputs = receipt.get("producer_inputs")
    if not isinstance(inputs, dict):
        raise ContractError(f"{sample}: producer receipt lacks producer_inputs")
    for required in (
        "tumor_bam",
        "normal_bam",
        "reference",
        "caller_raw_vcf",
        "longphase_input_vcf",
        "caller_pass_baseline_vcf",
    ):
        if not isinstance(inputs.get(required), dict):
            raise ContractError(f"{sample}: producer receipt lacks {required}")
    if canonical_sha256(inputs) != receipt.get("producer_input_binding_sha256"):
        raise ContractError(f"{sample}: producer input binding SHA is invalid")
    if canonical_sha256(receipt.get("command_argv")) != receipt.get("command_argv_sha256"):
        raise ContractError(f"{sample}: producer command argv SHA is invalid")
    outputs = receipt.get("capture_outputs")
    if not isinstance(outputs, dict):
        raise ContractError(f"{sample}: producer receipt lacks capture_outputs")
    expected_output_keys = {
        "sidecar",
        "sidecar_index",
        "native_validation",
        "stream_capture_summary",
        "normalization_audit",
        "filter_transition_audit",
        "sample_verification",
        "longphase_recalibrated_all_vcf",
        "longphase_recalibrated_pass_vcf",
    }
    if set(outputs) != expected_output_keys:
        raise ContractError(
            f"{sample}: producer capture output set mismatch: {sorted(outputs)}"
        )
    output_bindings = {
        "sidecar": "sidecar",
        "sidecar_index": "index",
        "native_validation": "validation",
    }
    for output_name, read_tag_name in output_bindings.items():
        output_spec = outputs.get(output_name)
        source_spec = read_tags.get(read_tag_name)
        if not isinstance(output_spec, dict) or output_spec != source_spec:
            raise ContractError(
                f"{sample}: producer capture_outputs.{output_name} differs from read_tags.{read_tag_name}"
            )
        artifact_path = Path(
            required_string(output_spec.get("path"), f"{sample}.{output_name}.path")
        )
        require_regular(artifact_path, f"{sample} {output_name}")
        expected_size = exact_integer(
            output_spec.get("identity", {}).get("size_bytes"),
            f"{sample}.{output_name}.size_bytes",
            minimum=1,
        )
        if artifact_path.stat().st_size != expected_size:
            raise ContractError(f"{sample}: {output_name} size changed from producer receipt")
    subject = read_tags.get("subject_binding")
    if not isinstance(subject, dict) or subject.get("sample") != sample:
        raise ContractError(f"{sample}: read-tag subject binding sample mismatch")
    if (
        read_tags.get("identity_schema") != receipt.get("identity_schema")
        or read_tags.get("identity_schema") != "coordinate_join_v1"
    ):
        raise ContractError(f"{sample}: coordinate identity schema mismatch")
    if (
        read_tags.get("duplicate_identity_policy") != receipt.get("duplicate_identity_policy")
        or subject.get("duplicate_identity_policy") != read_tags.get("duplicate_identity_policy")
    ):
        raise ContractError(f"{sample}: duplicate identity policy mismatch")
    if subject.get("coordinate_identity_columns") != [
        "QNAME",
        "CHROM",
        "START0",
        "END0",
        "FLAG",
        "CIGAR_B2",
    ]:
        raise ContractError(f"{sample}: coordinate identity columns mismatch")
    subject_expectations = {
        "sidecar_sha256": read_tags["sidecar"]["identity"]["sha256"],
        "sidecar_index_sha256": read_tags["index"]["identity"]["sha256"],
        "validation_sha256": read_tags["validation"]["identity"]["sha256"],
        "producer_capture_receipt_sha256": observed_sha,
        "producer_command_argv_sha256": receipt["command_argv_sha256"],
        "producer_input_binding_sha256": receipt["producer_input_binding_sha256"],
        "producer_effective_options_sha256": canonical_sha256(receipt.get("effective_options")),
        "caller_raw_vcf_sha256": inputs["caller_raw_vcf"]["identity"]["sha256"],
        "longphase_input_vcf_sha256": inputs["longphase_input_vcf"]["identity"]["sha256"],
        "caller_pass_baseline_vcf_sha256": inputs["caller_pass_baseline_vcf"]["identity"]["sha256"],
        "longphase_recalibrated_all_vcf_sha256": outputs[
            "longphase_recalibrated_all_vcf"
        ]["identity"]["sha256"],
        "longphase_recalibrated_pass_vcf_sha256": outputs[
            "longphase_recalibrated_pass_vcf"
        ]["identity"]["sha256"],
    }
    for field, expected in subject_expectations.items():
        if subject.get(field) != expected:
            raise ContractError(f"{sample}: subject binding {field} mismatch")
    counts = receipt.get("global_coordinate_counts")
    count_fields = {
        "mapped_alignment_count",
        "identity_unique_count",
        "duplicate_count",
        "conflict_count",
    }
    if not isinstance(counts, dict) or set(counts) != count_fields:
        raise ContractError(f"{sample}: producer receipt lacks global_coordinate_counts")
    for field in ("mapped_alignment_count", "identity_unique_count", "duplicate_count", "conflict_count"):
        if subject.get(field) != counts.get(field):
            raise ContractError(f"{sample}: subject binding {field} differs from producer receipt")
    if (
        counts["mapped_alignment_count"]
        != counts["identity_unique_count"] + counts["duplicate_count"]
        or counts["conflict_count"] != 0
    ):
        raise ContractError(f"{sample}: producer coordinate counts do not conserve")
    context = {
        "capture_outputs": outputs,
        "global_coordinate_counts": counts,
        "subject_binding": subject,
        "receipt_sha256": observed_sha,
    }
    evidence = {
        "path": str(receipt_path),
        "sha256": observed_sha,
        "schema_name": receipt.get("schema_name"),
        "schema_version": receipt.get("schema_version"),
        "same_sample_output_binding": "PASS",
    }
    return inputs, context, evidence


def topology_contract(path: Path) -> dict[str, dict[str, Any]]:
    require_regular(path, "topology cohort CSV")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"sample", "biological_id", "technical_replicate"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ContractError(f"topology CSV is missing columns: {sorted(missing)}")
        rows = list(reader)
    output: dict[str, dict[str, Any]] = {}
    for row in rows:
        sample = required_string(row.get("sample"), "topology.sample")
        if sample in output:
            raise ContractError(f"duplicate topology dataset: {sample}")
        technical_raw = str(row.get("technical_replicate", "")).lower()
        if technical_raw not in {"true", "false"}:
            raise ContractError(f"{sample}: technical_replicate must be true/false")
        output[sample] = {
            "biological_id": required_string(
                row.get("biological_id"), f"{sample}.biological_id"
            ),
            "technical_replicate": technical_raw == "true",
        }
    if set(output) != EXPECTED_DATASETS:
        raise ContractError(
            f"topology sample set mismatch: expected={sorted(EXPECTED_DATASETS)} observed={sorted(output)}"
        )
    return output


def validate_sidecar_receipt(
    sample: str, read_tags: dict[str, Any], producer_context: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    validation_spec = read_tags.get("validation")
    if not isinstance(validation_spec, dict):
        raise ContractError(f"{sample}: missing read_tags.validation")
    validation_path = Path(required_string(validation_spec.get("path"), f"{sample}.validation.path"))
    require_regular(validation_path, f"{sample} sidecar validation")
    expected_sha = required_string(
        validation_spec.get("identity", {}).get("sha256"), f"{sample}.validation.sha256"
    )
    validation, observed_sha = strict_json_load_with_sha256(validation_path)
    if not isinstance(validation, dict):
        raise ContractError(f"{sample}: sidecar validation must be a JSON object")
    native_validation = producer_context.get("capture_outputs", {}).get("native_validation")
    if native_validation != validation_spec:
        raise ContractError(f"{sample}: sidecar validation is not bound to producer output")
    subject = producer_context.get("subject_binding")
    if not isinstance(subject, dict) or subject.get("sample") != sample:
        raise ContractError(f"{sample}: sidecar subject binding mismatch")
    capture = validation.get("capture")
    execution = validation.get("execution")
    checks = validation.get("checks")
    if not isinstance(capture, dict) or not isinstance(execution, dict) or not isinstance(checks, dict):
        raise ContractError(f"{sample}: sidecar validation lacks capture/execution/checks")
    required_checks = (
        "truth_flags_absent",
        "parser_count_matches_input",
        "capture_pass",
        "execution_alignment_count_matches_capture",
        "sidecar_row_count_matches_capture",
        "tagged_count_matches_execution",
        "sidecar_no_malformed_rows",
        "sidecar_no_unknown_HP",
        "sidecar_no_exact_identity_conflicts",
        "recalibrated_preserves_all_input_keys",
    )
    failed_checks = [name for name in required_checks if checks.get(name) is not True]
    classes = capture.get("alignment_classes")
    if not isinstance(classes, dict):
        raise ContractError(f"{sample}: missing alignment_classes")
    total = exact_integer(classes.get("total"), f"{sample}.alignment.total", minimum=1)
    primary = exact_integer(classes.get("primary"), f"{sample}.alignment.primary")
    secondary = exact_integer(classes.get("secondary"), f"{sample}.alignment.secondary")
    supplementary = exact_integer(
        classes.get("supplementary"), f"{sample}.alignment.supplementary"
    )
    if primary + secondary + supplementary != total:
        raise ContractError(f"{sample}: alignment classes do not sum to total")
    if exact_integer(execution.get("total_alignment"), f"{sample}.execution.total_alignment") != total:
        raise ContractError(f"{sample}: execution total differs from capture total")
    tagged = exact_integer(execution.get("total_tagged"), f"{sample}.total_tagged")
    hp_ps_counts = validation.get("HP_with_PS_counts")
    if not isinstance(hp_ps_counts, dict):
        raise ContractError(f"{sample}: missing HP_with_PS_counts")
    hp_ps = sum(
        exact_integer(value, f"{sample}.HP_with_PS_counts.{key}")
        for key, value in hp_ps_counts.items()
    )
    record_missing = exact_integer(
        validation.get("record_key_missing"), f"{sample}.record_key_missing"
    )
    record_extra = exact_integer(
        validation.get("record_key_extra"), f"{sample}.record_key_extra"
    )
    duplicate_rows = exact_integer(
        validation.get("duplicate_exact_alignment_rows"),
        f"{sample}.duplicate_exact_alignment_rows",
    )
    conflicts = exact_integer(
        validation.get("duplicate_exact_alignment_conflicts"),
        f"{sample}.duplicate_exact_alignment_conflicts",
    )
    unknown_hp = validation.get("unknown_HP_counts")
    if not isinstance(unknown_hp, dict):
        raise ContractError(f"{sample}: unknown_HP_counts must be an object")
    producer_counts = producer_context.get("global_coordinate_counts")
    if not isinstance(producer_counts, dict):
        raise ContractError(f"{sample}: producer coordinate counts are unavailable")
    expected_counts = {
        "mapped_alignment_count": total,
        "identity_unique_count": total - duplicate_rows,
        "duplicate_count": duplicate_rows,
        "conflict_count": conflicts,
    }
    if any(producer_counts.get(field) != value for field, value in expected_counts.items()):
        raise ContractError(f"{sample}: sidecar validation counts differ from producer receipt")
    if any(subject.get(field) != value for field, value in expected_counts.items()):
        raise ContractError(f"{sample}: sidecar validation counts differ from subject binding")
    validation_pass = (
        validation.get("pass") is True
        and observed_sha == expected_sha
        and not failed_checks
        and record_missing == 0
        and record_extra == 0
        and conflicts == 0
        and not unknown_hp
    )
    row = {
        "authority": required_string(read_tags.get("authority"), f"{sample}.read_tags.authority"),
        "assurance": required_string(read_tags.get("assurance"), f"{sample}.read_tags.assurance"),
        "denominator_population": (
            "all_captured_mapped_alignment_records_including_primary_secondary_supplementary"
        ),
        "validation_path": str(validation_path),
        "validation_sha256": observed_sha,
        "validation_sha256_status": "MATCH" if observed_sha == expected_sha else "MISMATCH",
        "validation_status": "PASS" if validation_pass else "FAIL",
        "total_alignment_records": total,
        "primary_alignment_records": primary,
        "secondary_alignment_records": secondary,
        "supplementary_alignment_records": supplementary,
        "hp_assigned_alignment_records": tagged,
        "hp_assigned_rate_all_alignment_records": ratio(
            tagged, total, f"{sample}.HP assigned rate"
        ),
        "hp_and_ps_alignment_records": hp_ps,
        "hp_and_ps_rate_all_alignment_records": ratio(
            hp_ps, total, f"{sample}.HP+PS/all rate"
        ),
        "hp_and_ps_rate_among_hp_assigned_alignment_records": ratio(
            hp_ps, tagged, f"{sample}.HP+PS/HP rate"
        ),
        "record_key_missing": record_missing,
        "record_key_extra": record_extra,
        "duplicate_identity_rows": duplicate_rows,
        "duplicate_identity_rate_all_alignment_records": ratio(
            duplicate_rows, total, f"{sample}.duplicate/all rate"
        ),
        "duplicate_identity_conflicts": conflicts,
        "unknown_hp_category_n": len(unknown_hp),
    }
    evidence = {
        "failed_checks": failed_checks,
        "source_schema_version": validation.get("schema_version"),
        "same_sample_evidence_chain": "PASS",
    }
    return row, evidence


def build_dataset(
    source: dict[str, Any], topology: dict[str, Any], samtools: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    sample = required_string(source.get("sample"), "source.sample")
    biological_id = required_string(source.get("biological_id"), f"{sample}.biological_id")
    if biological_id != topology["biological_id"]:
        raise ContractError(
            f"{sample}: biological_id mismatch source={biological_id} topology={topology['biological_id']}"
        )
    replicate_role = required_string(source.get("replicate_role"), f"{sample}.replicate_role")
    technical = replicate_role == "platform_replica"
    if replicate_role not in {"canonical", "platform_replica"}:
        raise ContractError(f"{sample}: unsupported replicate_role {replicate_role!r}")
    if technical != topology["technical_replicate"]:
        raise ContractError(f"{sample}: technical replicate mismatch between manifests")
    alignment = source.get("alignment_payload")
    if not isinstance(alignment, dict):
        raise ContractError(f"{sample}: missing alignment_payload")
    read_tags = source.get("read_tags")
    if not isinstance(read_tags, dict):
        raise ContractError(f"{sample}: missing read_tags")
    receipt_inputs, producer_context, receipt_evidence = producer_inputs(sample, read_tags)
    if receipt_inputs["tumor_bam"] != {
        "path": alignment.get("path"),
        "index_path": alignment.get("index_path"),
        "storage_identity_v1": alignment.get("storage_identity_v1"),
    }:
        raise ContractError(f"{sample}: source manifest tumor BAM differs from producer receipt")
    subject = producer_context["subject_binding"]
    if (
        subject.get("alignment_payload_storage_identity_sha256")
        != alignment.get("storage_identity_v1", {}).get("identity_sha256")
    ):
        raise ContractError(f"{sample}: read-tag subject is bound to a different tumor BAM")
    somatic = source.get("somatic")
    if not isinstance(somatic, dict):
        raise ContractError(f"{sample}: source manifest lacks somatic inputs")
    for artifact_name in (
        "caller_raw_vcf",
        "longphase_input_vcf",
        "caller_pass_baseline_vcf",
    ):
        if somatic.get(artifact_name) != receipt_inputs.get(artifact_name):
            raise ContractError(
                f"{sample}: source manifest {artifact_name} differs from producer receipt"
            )
    for artifact_name in (
        "longphase_recalibrated_all_vcf",
        "longphase_recalibrated_pass_vcf",
    ):
        if somatic.get(artifact_name) != producer_context["capture_outputs"].get(artifact_name):
            raise ContractError(
                f"{sample}: source manifest {artifact_name} differs from producer output"
            )
    tumor_record, _, tumor_evidence = build_bam_record(
        receipt_inputs["tumor_bam"], samtools, f"{sample}.tumor_bam"
    )
    normal_record, _, normal_evidence = build_bam_record(
        receipt_inputs["normal_bam"], samtools, f"{sample}.normal_bam"
    )
    reference_record, reference_evidence = build_reference_record(
        receipt_inputs["reference"], f"{sample}.reference"
    )
    reference_dictionary_sha = reference_record["reference_dictionary_sha256"]
    if tumor_record["reference_dictionary_sha256"] != reference_dictionary_sha:
        raise ContractError(f"{sample}: tumor BAM dictionary differs from reference FAI")
    if normal_record["reference_dictionary_sha256"] != reference_dictionary_sha:
        raise ContractError(f"{sample}: normal BAM dictionary differs from reference FAI")
    if tumor_record["reference_dictionary_sha256"] != normal_record["reference_dictionary_sha256"]:
        raise ContractError(f"{sample}: tumor and normal BAM dictionaries differ")
    caller_spec = receipt_inputs["caller_raw_vcf"]
    caller_path = Path(required_string(caller_spec.get("path"), f"{sample}.caller.path"))
    require_regular(caller_path, f"{sample} caller VCF")
    caller_expected_sha = required_string(
        caller_spec.get("identity", {}).get("sha256"), f"{sample}.caller.sha256"
    )
    caller_observed_sha = file_sha256(caller_path)
    if caller_observed_sha != caller_expected_sha:
        raise ContractError(f"{sample}: caller VCF SHA changed")
    source_caller = somatic.get("caller_raw_vcf")
    if source_caller != caller_spec:
        raise ContractError(f"{sample}: source manifest caller VCF differs from producer receipt")
    alignment_source_directory = str(Path(tumor_record["path"]).parent)
    caller_source_directory = str(caller_path.parent.parent)
    source_pairing_status = (
        "SAME_DIRECTORY_FAMILY"
        if Path(alignment_source_directory) == Path(caller_source_directory)
        else "CROSS_DIRECTORY_REVIEW_REQUIRED"
    )
    read_group_status = (
        "AVAILABLE_UNVERIFIED"
        if tumor_record["read_group_samples"] or normal_record["read_group_samples"]
        else "UNAVAILABLE_NO_RG"
    )
    read_tag_row, read_tag_evidence = validate_sidecar_receipt(
        sample, read_tags, producer_context
    )
    output = {
        "dataset_id": sample,
        "biological_id": biological_id,
        "replicate_type": "technical_platform_replica" if technical else "canonical",
        "technical_replicate": technical,
        "platform": required_string(source.get("platform"), f"{sample}.platform"),
        "bam": tumor_record,
        "paired_normal_bam": normal_record,
        "reference": reference_record,
        "input_compatibility": {
            "tumor_reference_dictionary_status": "PASS",
            "normal_reference_dictionary_status": "PASS",
            "tumor_normal_dictionary_status": "PASS",
            "read_group_identity_status": read_group_status,
            "caller_input_path": str(caller_path),
            "caller_input_sha256": caller_observed_sha,
            "alignment_source_directory": alignment_source_directory,
            "caller_source_directory": caller_source_directory,
            "source_directory_pairing_status": source_pairing_status,
            "interpretation_note": (
                "Source directories differ; this may be an intentional variant-transfer design, "
                "but the dataset comparison cannot be attributed to platform/basecaller alone."
                if source_pairing_status == "CROSS_DIRECTORY_REVIEW_REQUIRED"
                else "Alignment and caller inputs resolve under the same source directory family."
            ),
        },
        "producer_read_tags": read_tag_row,
    }
    evidence = {
        "dataset_id": sample,
        "producer_capture_receipt": receipt_evidence,
        "tumor_bam": tumor_evidence,
        "normal_bam": normal_evidence,
        "reference": reference_evidence,
        "caller_input": {
            "path": str(caller_path),
            "sha256": caller_observed_sha,
            "source_directory_pairing_status": source_pairing_status,
        },
        "read_tag_validation": read_tag_evidence,
    }
    return output, evidence


def validate_output(manifest: dict[str, Any], schema: dict[str, Any] | None = None) -> None:
    if schema is None:
        loaded = strict_json_load(DEFAULT_SCHEMA)
        if not isinstance(loaded, dict):
            raise ContractError(f"output JSON Schema must be an object: {DEFAULT_SCHEMA}")
        schema = loaded
    validate_json_schema(manifest, schema)
    if manifest.get("schema_name") != "intersubmod.multi_bam_dashboard_input_manifest":
        raise ContractError("output schema_name mismatch")
    if manifest.get("schema_version") != "1.1.0":
        raise ContractError("output schema_version mismatch")
    datasets = manifest.get("datasets")
    if not isinstance(datasets, list) or len(datasets) != 7:
        raise ContractError("output must contain exactly seven datasets")
    if any(not isinstance(row, dict) for row in datasets):
        raise ContractError("output dataset rows must be JSON objects")
    ids = [row.get("dataset_id") for row in datasets]
    if len(set(ids)) != 7 or set(ids) != EXPECTED_DATASETS:
        raise ContractError(f"output dataset IDs mismatch: {ids}")
    biological_ids = {row.get("biological_id") for row in datasets}
    if len(biological_ids) != EXPECTED_BIOLOGICAL_N:
        raise ContractError(f"expected six biological IDs, got {sorted(biological_ids)}")
    technical_n = sum(row.get("technical_replicate") is True for row in datasets)
    if technical_n != 1:
        raise ContractError(f"expected one technical replicate, got {technical_n}")
    for row in datasets:
        sample = required_string(row.get("dataset_id"), "output.dataset_id")
        biological_id = required_string(
            row.get("biological_id"), f"{sample}.biological_id"
        )
        technical = row.get("technical_replicate")
        if not isinstance(technical, bool):
            raise ContractError(f"{sample}: technical_replicate must be boolean")
        expected_technical = sample == "HCC1395_DORADO"
        expected_biological_id = "HCC1395" if expected_technical else sample
        if technical != expected_technical or biological_id != expected_biological_id:
            raise ContractError(f"{sample}: biological/technical identity contract mismatch")
        records = {
            key: row.get(key)
            for key in (
                "bam",
                "paired_normal_bam",
                "reference",
                "input_compatibility",
                "producer_read_tags",
            )
        }
        if any(not isinstance(value, dict) for value in records.values()):
            raise ContractError(f"{sample}: output row lacks required record objects")
        bam = records["bam"]
        normal = records["paired_normal_bam"]
        reference = records["reference"]
        compatibility = records["input_compatibility"]
        tags = records["producer_read_tags"]
        for role, record in (("tumor", bam), ("normal", normal)):
            if record["bounded_payload_status"] != "MATCH":
                raise ContractError(f"{sample}: sampled {role} BAM payload mismatch")
            if record["index_sha256_status"] != "MATCH":
                raise ContractError(f"{sample}: {role} BAM index SHA mismatch")
            if record["quickcheck_status"] != "PASS":
                raise ContractError(f"{sample}: {role} samtools quickcheck failed")
            strict_status = record.get("strict_storage_identity_status")
            strict_fields = set(record.get("strict_differing_fields", []))
            if not (
                (strict_status == "MATCH" and not strict_fields)
                or (
                    strict_status == "MOUNT_DEVICE_DRIFT_ONLY"
                    and strict_fields == {"identity_sha256", "st_dev"}
                )
            ):
                raise ContractError(f"{sample}: {role} strict identity is not bounded-pass")
            if record.get("full_content_sha256") is not None:
                raise ContractError(
                    f"{sample}: full BAM SHA must remain null until a full-content scan is run"
                )
        if reference["bounded_payload_status"] != "MATCH" or reference["fai_sha256_status"] != "MATCH":
            raise ContractError(f"{sample}: reference sampled payload or FAI mismatch")
        reference_strict_status = reference.get("strict_storage_identity_status")
        reference_strict_fields = set(reference.get("strict_differing_fields", []))
        if not (
            (
                reference_strict_status == "MATCH"
                and not reference_strict_fields
            )
            or (
                reference_strict_status == "MOUNT_DEVICE_DRIFT_ONLY"
                and reference_strict_fields == {"identity_sha256", "st_dev"}
            )
        ) or reference.get("full_content_sha256") is not None:
            raise ContractError(f"{sample}: reference identity is not bounded-pass")
        if any(
            compatibility[field] != "PASS"
            for field in (
                "tumor_reference_dictionary_status",
                "normal_reference_dictionary_status",
                "tumor_normal_dictionary_status",
            )
        ):
            raise ContractError(f"{sample}: BAM/reference dictionary mismatch")
        if len({
            bam.get("reference_dictionary_sha256"),
            normal.get("reference_dictionary_sha256"),
            reference.get("reference_dictionary_sha256"),
        }) != 1:
            raise ContractError(f"{sample}: dictionary SHA values do not reconcile")
        expected_read_group_status = (
            "AVAILABLE_UNVERIFIED"
            if bam.get("read_group_samples") or normal.get("read_group_samples")
            else "UNAVAILABLE_NO_RG"
        )
        if compatibility.get("read_group_identity_status") != expected_read_group_status:
            raise ContractError(f"{sample}: read-group identity status does not match BAM headers")
        expected_pairing = (
            "CROSS_DIRECTORY_REVIEW_REQUIRED"
            if sample == "HCC1395"
            else "SAME_DIRECTORY_FAMILY"
        )
        if compatibility.get("source_directory_pairing_status") != expected_pairing:
            raise ContractError(f"{sample}: caller/alignment source-directory pairing changed")
        if tags["validation_status"] != "PASS":
            raise ContractError(f"{sample}: producer tag validation failed")
        if tags.get("validation_sha256_status") != "MATCH":
            raise ContractError(f"{sample}: producer tag validation SHA mismatch")
        class_counts = {
            field: exact_integer(tags.get(field), f"{sample}.{field}")
            for field in (
                "total_alignment_records",
                "primary_alignment_records",
                "secondary_alignment_records",
                "supplementary_alignment_records",
                "hp_assigned_alignment_records",
                "hp_and_ps_alignment_records",
                "record_key_missing",
                "record_key_extra",
                "duplicate_identity_rows",
                "duplicate_identity_conflicts",
                "unknown_hp_category_n",
            )
        }
        if (
            class_counts["primary_alignment_records"]
            + class_counts["secondary_alignment_records"]
            + class_counts["supplementary_alignment_records"]
            != class_counts["total_alignment_records"]
        ):
            raise ContractError(f"{sample}: alignment-class counts do not sum to total")
        if class_counts["hp_assigned_alignment_records"] > class_counts["total_alignment_records"]:
            raise ContractError(f"{sample}: tagged alignments exceed total")
        if class_counts["hp_and_ps_alignment_records"] > class_counts["hp_assigned_alignment_records"]:
            raise ContractError(f"{sample}: HP+PS alignments exceed HP-assigned")
        expected_hp_rate = ratio(
            class_counts["hp_assigned_alignment_records"],
            class_counts["total_alignment_records"],
            f"{sample}.HP/all",
        )
        expected_hp_ps_rate = ratio(
            class_counts["hp_and_ps_alignment_records"],
            class_counts["total_alignment_records"],
            f"{sample}.HP+PS/all",
        )
        expected_hp_ps_among_hp = ratio(
            class_counts["hp_and_ps_alignment_records"],
            class_counts["hp_assigned_alignment_records"],
            f"{sample}.HP+PS/HP",
        )
        expected_duplicate_rate = ratio(
            class_counts["duplicate_identity_rows"],
            class_counts["total_alignment_records"],
            f"{sample}.duplicate/all",
        )
        for field, expected in (
            ("hp_assigned_rate_all_alignment_records", expected_hp_rate),
            ("hp_and_ps_rate_all_alignment_records", expected_hp_ps_rate),
            (
                "hp_and_ps_rate_among_hp_assigned_alignment_records",
                expected_hp_ps_among_hp,
            ),
            ("duplicate_identity_rate_all_alignment_records", expected_duplicate_rate),
        ):
            observed = tags.get(field)
            if (
                isinstance(observed, bool)
                or not isinstance(observed, (int, float))
                or not math.isclose(float(observed), expected, rel_tol=0, abs_tol=1e-12)
            ):
                raise ContractError(f"{sample}: {field} does not reconcile with exact counts")
    summary = manifest.get("verification_summary")
    if not isinstance(summary, dict):
        raise ContractError("output lacks verification_summary")
    identity_roles = [
        (row["bam"], row["paired_normal_bam"], row["reference"])
        for row in datasets
    ]
    expected_summary = {
        "dataset_count": 7,
        "bam_present_n": sum(Path(row["bam"]["path"]).is_file() for row in datasets),
        "bai_present_n": sum(Path(row["bam"]["index_path"]).is_file() for row in datasets),
        "quickcheck_pass_n": sum(row["bam"]["quickcheck_status"] == "PASS" for row in datasets),
        "tumor_bounded_payload_match_n": sum(
            row["bam"]["bounded_payload_status"] == "MATCH" for row in datasets
        ),
        "tumor_strict_storage_identity_match_n": sum(
            row["bam"]["strict_storage_identity_status"] == "MATCH" for row in datasets
        ),
        "tumor_mount_device_drift_only_n": sum(
            row["bam"]["strict_storage_identity_status"] == "MOUNT_DEVICE_DRIFT_ONLY"
            for row in datasets
        ),
        "all_input_roles_strict_storage_identity_match_n": sum(
            all(record["strict_storage_identity_status"] == "MATCH" for record in roles)
            for roles in identity_roles
        ),
        "all_input_roles_mount_device_drift_only_n": sum(
            all(
                record["strict_storage_identity_status"]
                in {"MATCH", "MOUNT_DEVICE_DRIFT_ONLY"}
                for record in roles
            )
            and any(
                record["strict_storage_identity_status"] == "MOUNT_DEVICE_DRIFT_ONLY"
                for record in roles
            )
            for roles in identity_roles
        ),
        "bai_sha256_match_n": sum(row["bam"]["index_sha256_status"] == "MATCH" for row in datasets),
        "full_bam_sha256_n": sum(row["bam"]["full_content_sha256"] is not None for row in datasets),
        "normal_bam_present_n": sum(Path(row["paired_normal_bam"]["path"]).is_file() for row in datasets),
        "normal_bai_present_n": sum(Path(row["paired_normal_bam"]["index_path"]).is_file() for row in datasets),
        "normal_quickcheck_pass_n": sum(
            row["paired_normal_bam"]["quickcheck_status"] == "PASS" for row in datasets
        ),
        "normal_bounded_payload_match_n": sum(
            row["paired_normal_bam"]["bounded_payload_status"] == "MATCH" for row in datasets
        ),
        "reference_bounded_payload_match_n": sum(
            row["reference"]["bounded_payload_status"] == "MATCH" for row in datasets
        ),
        "reference_dictionary_compatible_n": sum(
            row["input_compatibility"]["tumor_reference_dictionary_status"] == "PASS"
            and row["input_compatibility"]["normal_reference_dictionary_status"] == "PASS"
            and row["input_compatibility"]["tumor_normal_dictionary_status"] == "PASS"
            for row in datasets
        ),
        "no_read_group_n": sum(
            row["input_compatibility"]["read_group_identity_status"] == "UNAVAILABLE_NO_RG"
            for row in datasets
        ),
        "cross_directory_pairing_n": sum(
            row["input_compatibility"]["source_directory_pairing_status"]
            == "CROSS_DIRECTORY_REVIEW_REQUIRED"
            for row in datasets
        ),
        "header_read_n": sum(bool(row["bam"]["header_sha256"]) for row in datasets),
        "sidecar_validation_pass_n": sum(
            row["producer_read_tags"]["validation_status"] == "PASS" for row in datasets
        ),
        "read_tag_conservation_pass_n": sum(
            row["producer_read_tags"]["record_key_missing"] == 0
            and row["producer_read_tags"]["record_key_extra"] == 0
            and row["producer_read_tags"]["duplicate_identity_conflicts"] == 0
            and row["producer_read_tags"]["unknown_hp_category_n"] == 0
            for row in datasets
        ),
        "sample_set_match": True,
    }
    for field, expected in expected_summary.items():
        if summary.get(field) != expected:
            raise ContractError(
                f"verification_summary.{field} expected {expected!r}, got {summary.get(field)!r}"
            )
    expected_status = (
        "PASS_BOUNDED"
        if expected_summary["all_input_roles_strict_storage_identity_match_n"] == 7
        else "PASS_BOUNDED_WITH_METADATA_DRIFT"
    )
    if summary.get("status") != expected_status:
        raise ContractError(
            f"verification summary status expected {expected_status}, got {summary.get('status')}"
        )
    canonical_sha256(manifest)


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except BaseException:
        if temporary_path.exists():
            temporary_path.rename(temporary_path.with_suffix(".failed"))
        raise


def build_manifest(args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, Any]]:
    source_manifest_path = args.source_manifest.resolve(strict=True)
    topology_path = args.topology_csv.resolve(strict=True)
    source_schema_path = args.source_schema.resolve(strict=True)
    schema_path = args.schema.resolve(strict=True)
    schema, schema_sha = strict_json_load_with_sha256(schema_path)
    if not isinstance(schema, dict) or schema.get("$schema") != "https://json-schema.org/draft/2020-12/schema":
        raise ContractError(f"unexpected JSON schema document: {schema_path}")
    source_schema, source_schema_sha = strict_json_load_with_sha256(source_schema_path)
    validate_canonical_source_schema(
        source_schema,
        source_schema_sha,
        path=source_schema_path,
    )
    source, source_manifest_sha = strict_json_load_with_sha256(source_manifest_path)
    if not isinstance(source, dict):
        raise ContractError("source manifest must be a JSON object")
    validate_json_schema(source, source_schema)
    if (
        source.get("schema_name") != "intersubmod.layered_input_manifest"
        or source.get("schema_version") != "3.0.0"
        or source.get("dataset_count") != 7
        or source.get("biological_sample_count") != 6
    ):
        raise ContractError("source manifest schema or declared scope mismatch")
    topology = topology_contract(topology_path)
    samples = source.get("samples")
    if not isinstance(samples, list) or len(samples) != 7:
        raise ContractError("source manifest must contain seven sample rows")
    source_ids = [item.get("sample") for item in samples if isinstance(item, dict)]
    if len(source_ids) != 7 or len(set(source_ids)) != 7 or set(source_ids) != EXPECTED_DATASETS:
        raise ContractError(f"source sample set mismatch: {source_ids}")
    samtools_release = samtools_version(args.samtools)
    datasets: list[dict[str, Any]] = []
    evidence: list[dict[str, Any]] = []
    for row in samples:
        if not isinstance(row, dict):
            raise ContractError("source sample rows must be JSON objects")
        dataset, row_evidence = build_dataset(row, topology[row["sample"]], args.samtools)
        datasets.append(dataset)
        evidence.append(row_evidence)
    receipt_paths = [row["producer_capture_receipt"]["path"] for row in evidence]
    if len(set(receipt_paths)) != len(receipt_paths):
        raise ContractError("producer capture receipt paths must be unique per dataset")
    production_summary = source.get("production_summary")
    if not isinstance(production_summary, dict):
        raise ContractError("source manifest lacks production_summary")
    production_rows = production_summary.get("datasets")
    if not isinstance(production_rows, list):
        raise ContractError("source production_summary.datasets must be an array")
    production_by_sample = {
        row.get("sample"): row for row in production_rows if isinstance(row, dict)
    }
    if set(production_by_sample) != EXPECTED_DATASETS:
        raise ContractError("source production summary sample set mismatch")
    source_by_sample = {row["sample"]: row for row in samples}
    for sample in EXPECTED_DATASETS:
        production_row = production_by_sample[sample]
        source_row = source_by_sample[sample]
        if production_row.get("pass") is not True:
            raise ContractError(f"{sample}: source production summary is not PASS")
        if (
            production_row.get("validation_sha256")
            != source_row["read_tags"]["validation"]["identity"]["sha256"]
        ):
            raise ContractError(f"{sample}: production summary validation SHA mismatch")
    identity_roles = [
        (row["bam"], row["paired_normal_bam"], row["reference"])
        for row in datasets
    ]
    all_input_roles_strict_n = sum(
        all(record["strict_storage_identity_status"] == "MATCH" for record in roles)
        for roles in identity_roles
    )
    all_input_roles_drift_n = sum(
        all(
            record["strict_storage_identity_status"]
            in {"MATCH", "MOUNT_DEVICE_DRIFT_ONLY"}
            for record in roles
        )
        and any(
            record["strict_storage_identity_status"] == "MOUNT_DEVICE_DRIFT_ONLY"
            for record in roles
        )
        for roles in identity_roles
    )
    verification_summary = {
        "status": (
            "PASS_BOUNDED"
            if all_input_roles_strict_n == len(datasets)
            else "PASS_BOUNDED_WITH_METADATA_DRIFT"
        ),
        "dataset_count": len(datasets),
        "bam_present_n": sum(Path(row["bam"]["path"]).is_file() for row in datasets),
        "bai_present_n": sum(Path(row["bam"]["index_path"]).is_file() for row in datasets),
        "quickcheck_pass_n": sum(row["bam"]["quickcheck_status"] == "PASS" for row in datasets),
        "tumor_bounded_payload_match_n": sum(
            row["bam"]["bounded_payload_status"] == "MATCH" for row in datasets
        ),
        "tumor_strict_storage_identity_match_n": sum(
            row["bam"]["strict_storage_identity_status"] == "MATCH" for row in datasets
        ),
        "tumor_mount_device_drift_only_n": sum(
            row["bam"]["strict_storage_identity_status"] == "MOUNT_DEVICE_DRIFT_ONLY"
            for row in datasets
        ),
        "all_input_roles_strict_storage_identity_match_n": all_input_roles_strict_n,
        "all_input_roles_mount_device_drift_only_n": all_input_roles_drift_n,
        "bai_sha256_match_n": sum(
            row["bam"]["index_sha256_status"] == "MATCH" for row in datasets
        ),
        "full_bam_sha256_n": sum(row["bam"]["full_content_sha256"] is not None for row in datasets),
        "normal_bam_present_n": sum(
            Path(row["paired_normal_bam"]["path"]).is_file() for row in datasets
        ),
        "normal_bai_present_n": sum(
            Path(row["paired_normal_bam"]["index_path"]).is_file() for row in datasets
        ),
        "normal_quickcheck_pass_n": sum(
            row["paired_normal_bam"]["quickcheck_status"] == "PASS" for row in datasets
        ),
        "normal_bounded_payload_match_n": sum(
            row["paired_normal_bam"]["bounded_payload_status"] == "MATCH"
            for row in datasets
        ),
        "reference_bounded_payload_match_n": sum(
            row["reference"]["bounded_payload_status"] == "MATCH" for row in datasets
        ),
        "reference_dictionary_compatible_n": sum(
            row["input_compatibility"]["tumor_reference_dictionary_status"] == "PASS"
            and row["input_compatibility"]["normal_reference_dictionary_status"] == "PASS"
            and row["input_compatibility"]["tumor_normal_dictionary_status"] == "PASS"
            for row in datasets
        ),
        "no_read_group_n": sum(
            row["input_compatibility"]["read_group_identity_status"] == "UNAVAILABLE_NO_RG"
            for row in datasets
        ),
        "cross_directory_pairing_n": sum(
            row["input_compatibility"]["source_directory_pairing_status"]
            == "CROSS_DIRECTORY_REVIEW_REQUIRED"
            for row in datasets
        ),
        "header_read_n": sum(bool(row["bam"]["header_sha256"]) for row in datasets),
        "sidecar_validation_pass_n": sum(
            row["producer_read_tags"]["validation_status"] == "PASS" for row in datasets
        ),
        "read_tag_conservation_pass_n": sum(
            row["producer_read_tags"]["record_key_missing"] == 0
            and row["producer_read_tags"]["record_key_extra"] == 0
            and row["producer_read_tags"]["duplicate_identity_conflicts"] == 0
            and row["producer_read_tags"]["unknown_hp_category_n"] == 0
            for row in datasets
        ),
        "sample_set_match": set(source_ids) == set(topology) == EXPECTED_DATASETS,
    }
    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    manifest = {
        "schema_name": "intersubmod.multi_bam_dashboard_input_manifest",
        "schema_version": "1.1.0",
        "manifest_id": f"multi_bam_input_20260813_{source_manifest_sha[:12]}",
        "generated_at_utc": generated_at,
        "task_type": "comprehensive_validation",
        "claim_ceiling": "bounded_input_identity_and_existing_producer_receipts_only",
        "source_manifest": {
            "path": str(source_manifest_path),
            "sha256": source_manifest_sha,
            "schema_name": source["schema_name"],
            "schema_version": source["schema_version"],
            "manifest_id": required_string(source.get("manifest_id"), "source.manifest_id"),
            "created_at_utc": required_string(
                source.get("created_at_utc"), "source.created_at_utc"
            ),
        },
        "dataset_count": len(datasets),
        "biological_sample_count": len({row["biological_id"] for row in datasets}),
        "verification_summary": verification_summary,
        "datasets": datasets,
    }
    validate_output(manifest, schema)
    receipt = {
        "schema_name": "intersubmod.multi_bam_dashboard_input_manifest_validation",
        "schema_version": "1.0.0",
        "generated_at_utc": generated_at,
        "status": verification_summary["status"],
        "claim_ceiling": manifest["claim_ceiling"],
        "command_argv": sys.argv,
        "samtools_version": samtools_release,
        "builder": {
            "path": source_path(Path(__file__)),
            "sha256": file_sha256(Path(__file__)),
        },
        "schema": {
            "path": source_path(schema_path),
            "sha256": schema_sha,
            "validation_engine": (
                "intersubmod bounded Draft 2020-12 keyword validator plus semantic validation"
            ),
            "supported_keywords": sorted(SUPPORTED_SCHEMA_KEYWORDS),
        },
        "source_schema": {
            "path": source_path(source_schema_path),
            "sha256": source_schema_sha,
            "validation_engine": "intersubmod bounded Draft 2020-12 keyword validator",
        },
        "inputs": {
            "source_manifest": {
                "path": str(source_manifest_path),
                "sha256": manifest["source_manifest"]["sha256"],
            },
            "topology_csv": {
                "path": source_path(topology_path),
                "sha256": file_sha256(topology_path),
            },
        },
        "verification_summary": verification_summary,
        "dataset_evidence": evidence,
    }
    return manifest, receipt


def validate_io_path_separation(args: argparse.Namespace) -> None:
    inputs = {
        "source_manifest": args.source_manifest.resolve(),
        "topology_csv": args.topology_csv.resolve(),
        "source_schema": args.source_schema.resolve(),
        "output_schema": args.schema.resolve(),
    }
    outputs = {
        "output": args.output.resolve(),
        "receipt": args.receipt.resolve(),
    }
    if outputs["output"] == outputs["receipt"]:
        raise ContractError("--output and --receipt must resolve to different paths")
    for output_name, output_path in outputs.items():
        collisions = [name for name, input_path in inputs.items() if output_path == input_path]
        if collisions:
            raise ContractError(
                f"--{output_name.replace('_', '-')} collides with input path(s): {collisions}"
            )


def main() -> int:
    args = parse_args()
    try:
        validate_io_path_separation(args)
        manifest, receipt = build_manifest(args)
        output_path = args.output.resolve()
        receipt_path = args.receipt.resolve()
        atomic_write_json(output_path, manifest)
        receipt["output"] = {
            "path": source_path(output_path),
            "sha256": file_sha256(output_path),
        }
        atomic_write_json(receipt_path, receipt)
        summary = {
            "status": receipt["status"],
            "input": str(args.source_manifest.resolve()),
            "output": source_path(output_path),
            "receipt": source_path(receipt_path),
            "dataset_count": manifest["dataset_count"],
            "biological_sample_count": manifest["biological_sample_count"],
            **manifest["verification_summary"],
        }
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 0
    except ContractError as exc:
        print(f"MULTI-BAM INPUT MANIFEST ERROR: {exc}", file=sys.stderr)
        return 7


if __name__ == "__main__":
    raise SystemExit(main())
