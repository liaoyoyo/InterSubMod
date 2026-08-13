#!/usr/bin/env python3
"""Shared fail-closed contract helpers for the public claim registry."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any


REGISTRY_SCHEMA_NAME = "intersubmod.public_claim_remediation_registry"
REGISTRY_SCHEMA_VERSION = "2.0.0"
REGISTRY_ID = "intersubmod-public-claims-remediation-20260813"
GENERATED_AT = "2026-08-13T00:00:00+08:00"
TASK_TYPE = ["B_COMPREHENSIVE_VALIDATION", "D_EXTERNAL_HANDOFF"]
ALLOWED_CURRENT_VERDICTS = ["CONFIRMED", "CONFIRMED_WITH_LIMITS", "UNVERIFIED"]
EXPECTED_PUBLIC_SOURCE_COUNT = 34
EXPECTED_INVENTORY_SHA256 = "d5d8794c2467ac6ac466711f91890666f3dbc571d9b3e27b07ba5f10d642ddbc"
PUBLICATION_CANDIDATE_STATUS = "LOCAL_WORKING_TREE_CANDIDATE_NOT_PUBLICATION_EVIDENCE"

BUILD_RECEIPT_SCHEMA_NAME = "intersubmod.public_claim_remediation_build_receipt"
BUILD_RECEIPT_SCHEMA_VERSION = "2.0.0"
BUILD_RECEIPT_ID = "intersubmod-public-claims-remediation-build-20260813"

C066_CANONICAL_WORDING = (
    "MathUtils::cramers_v returns raw Cramer's V plus a separate reliability flag. "
    "RegionProcessor projects unreliable values to 0 only in the summary CramersV field. "
    "Legacy GlobalTest::passed_gate and per-region significance JSON use raw V. These are "
    "three distinct contracts and must not be described as one reliability gate."
)

EXPECTED_GATES = {
    "P0_SOURCE_READY": {
        "status": "PASS",
        "reason": (
            "33 local P0 claims pass fail-closed source guards; C108 About was updated and "
            "re-fetched with a bounded receipt."
        ),
    },
    "SOURCE_READY": {
        "status": "BLOCKED",
        "reason": "P1/P2 source edits still require dedicated guards/evidence closure.",
    },
    "PUBLICATION_READY": {
        "status": "BLOCKED",
        "reason": (
            "Default branch, Wiki and Pages have not been published and re-fetched; About "
            "C108 alone is confirmed."
        ),
    },
    "RELEASE_READY": {
        "status": "BLOCKED",
        "reason": "Publication gate and broader research-handoff release gates remain open.",
    },
}

EXPECTED_STATUS_SEMANTICS = {
    "inventory_proposition": (
        "The proposition observed in the immutable 2026-08-12 audit; it may be the false wording "
        "that triggered remediation."
    ),
    "remediated_claim": (
        "The bounded corrected claim to which current_verdict applies; null means source remediation "
        "is still open."
    ),
    "inventory_evidence": (
        "Evidence wording copied byte-for-byte from the immutable 2026-08-12 audit. It is historical "
        "audit input, may describe a superseded branch label, and is not the evidence basis for "
        "current_verdict."
    ),
    "inventory_minimum_rewrite": (
        "Minimum rewrite copied byte-for-byte from the immutable 2026-08-12 audit. It is preserved "
        "for provenance only; remediated_claim is the current bounded wording."
    ),
    "source_status": "State of the checked-in source; it does not describe GitHub live bytes.",
    "live_status": "State obtained by re-fetching the published surface; local edits cannot change it.",
    "SOURCE_READY": "All declared local P0 anchors pass, but publication remains unverified.",
    "VALIDATED_DERIVED_WITH_LIMITS": (
        "A bounded remediated claim is supported by checked-in derived evidence and a replayable "
        "receipt, but remains explicitly non-authority and non-final."
    ),
}


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode("utf-8")


def canonical_object_sha256(value: Any) -> str:
    return sha256_bytes(canonical_json_bytes(value))


def atomic_write_json(path: Path, value: Any) -> None:
    """Atomically replace one JSON file with fsync-backed canonical bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(canonical_json_bytes(value))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _git(root: Path, *arguments: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *arguments],
        cwd=root,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _index_entry(root: Path, relative: str) -> dict[str, str | None]:
    result = _git(root, "ls-files", "--stage", "--", relative)
    lines = [line for line in result.stdout.decode("utf-8").splitlines() if line]
    if not lines:
        return {"state": "ABSENT", "git_blob_oid": None}
    parsed = []
    for line in lines:
        metadata, _path = line.split("\t", 1)
        mode, oid, stage = metadata.split()
        parsed.append((mode, oid, stage))
    stage_zero = [entry for entry in parsed if entry[2] == "0"]
    if len(stage_zero) != 1 or len(parsed) != 1:
        return {"state": "UNMERGED", "git_blob_oid": None}
    return {"state": "PRESENT", "git_blob_oid": stage_zero[0][1]}


def _head_entry(root: Path, relative: str) -> dict[str, str | None]:
    result = _git(root, "ls-tree", "HEAD", "--", relative)
    output = result.stdout.decode("utf-8").strip()
    if not output:
        return {"state": "ABSENT", "git_blob_oid": None}
    metadata, _path = output.split("\t", 1)
    _mode, object_type, oid = metadata.split()
    if object_type != "blob":
        return {"state": f"NON_BLOB_{object_type.upper()}", "git_blob_oid": oid}
    return {"state": "PRESENT", "git_blob_oid": oid}


def _working_entry(root: Path, relative: str, index: dict[str, str | None]) -> dict[str, str]:
    result = _git(root, "hash-object", "--", relative)
    oid = result.stdout.decode("ascii").strip()
    if index["state"] == "ABSENT":
        state = "UNTRACKED_CANDIDATE"
    elif index["state"] != "PRESENT":
        state = "UNMERGED_CANDIDATE"
    elif oid == index["git_blob_oid"]:
        state = "MATCHES_INDEX"
    else:
        state = "MODIFIED_CANDIDATE"
    return {"state": state, "git_blob_oid": oid}


def expand_public_sources(root: Path, globs: object) -> list[Path]:
    if not isinstance(globs, list) or not globs or any(not isinstance(item, str) or not item for item in globs):
        raise ValueError("cross_document_guards.public_source_globs must be a non-empty string list")
    paths: set[Path] = set()
    for pattern in globs:
        matches = sorted(path for path in root.glob(pattern) if path.is_file())
        if not matches:
            raise ValueError(f"public source glob matched no files: {pattern}")
        for path in matches:
            resolved = path.resolve()
            try:
                resolved.relative_to(root.resolve())
            except ValueError as error:
                raise ValueError(f"public source escapes repository: {path}") from error
            paths.add(path)
    ordered = sorted(paths, key=lambda path: path.relative_to(root).as_posix())
    if len(ordered) != EXPECTED_PUBLIC_SOURCE_COUNT:
        raise ValueError(
            f"public source expansion must contain exactly {EXPECTED_PUBLIC_SOURCE_COUNT} files, got {len(ordered)}"
        )
    return ordered


def build_public_source_manifest(root: Path, p0_registry: dict[str, Any]) -> dict[str, Any]:
    cross_guards = p0_registry.get("cross_document_guards")
    if not isinstance(cross_guards, dict):
        raise ValueError("P0 registry cross_document_guards must be an object")
    globs = cross_guards.get("public_source_globs")
    paths = expand_public_sources(root, globs)
    files = []
    for path in paths:
        relative = path.relative_to(root).as_posix()
        raw = path.read_bytes()
        index = _index_entry(root, relative)
        head = _head_entry(root, relative)
        working = _working_entry(root, relative, index)
        if index["state"] == "PRESENT":
            index["state"] = "MATCHES_HEAD" if index["git_blob_oid"] == head["git_blob_oid"] else "STAGED_CANDIDATE"
        files.append(
            {
                "path": relative,
                "bytes": len(raw),
                "sha256": sha256_bytes(raw),
                "publication_status": PUBLICATION_CANDIDATE_STATUS,
                "git": {"working": working, "index": index, "HEAD": head},
            }
        )
    working_states = Counter(item["git"]["working"]["state"] for item in files)
    index_states = Counter(item["git"]["index"]["state"] for item in files)
    return {
        "source": "p0_claim_registry.cross_document_guards.public_source_globs",
        "source_globs": list(globs),
        "publication_semantics": (
            "Working-tree bytes are local publication candidates only; index and HEAD blob identities "
            "are provenance and do not prove any live GitHub surface."
        ),
        "counts": {
            "files": len(files),
            "working_states": dict(sorted(working_states.items())),
            "index_states": dict(sorted(index_states.items())),
        },
        "files": files,
    }


def build_about_snapshot_manifest(root: Path, about_receipt: dict[str, Any]) -> dict[str, Any]:
    verification = about_receipt.get("verification")
    commands = verification.get("commands") if isinstance(verification, dict) else None
    if not isinstance(commands, list):
        raise ValueError("validated GitHub About receipt has no command snapshots")
    files = []
    for command in sorted(commands, key=lambda item: str(item.get("command"))):
        relative = command.get("response_snapshot")
        if not isinstance(relative, str):
            raise ValueError("GitHub About snapshot path must be a string")
        path = (root / relative).resolve()
        try:
            path.relative_to(root.resolve())
        except ValueError as error:
            raise ValueError(f"GitHub About snapshot escapes repository: {relative}") from error
        raw = path.read_bytes()
        files.append(
            {
                "command": command.get("command"),
                "path": relative,
                "bytes": len(raw),
                "sha256": sha256_bytes(raw),
            }
        )
    if len(files) != 3:
        raise ValueError(f"GitHub About snapshot manifest must contain exactly 3 files, got {len(files)}")
    return {"counts": {"files": len(files)}, "files": files}


def build_counts(claims: list[dict[str, Any]]) -> dict[str, Any]:
    def count(field: str) -> dict[str, int]:
        return dict(sorted(Counter(str(item.get(field)) for item in claims).items()))

    return {
        "claims": len(claims),
        "by_priority": count("priority"),
        "by_current_verdict": count("current_verdict"),
        "by_source_status": count("source_status"),
    }


def build_receipt_payload(registry: dict[str, Any], output_sha256: str) -> dict[str, Any]:
    return {
        "schema_name": BUILD_RECEIPT_SCHEMA_NAME,
        "schema_version": BUILD_RECEIPT_SCHEMA_VERSION,
        "receipt_id": BUILD_RECEIPT_ID,
        "generated_at": GENERATED_AT,
        "output": {
            "path": registry["output_path"],
            "sha256": output_sha256,
            "schema_name": registry["schema_name"],
            "schema_version": registry["schema_version"],
            "registry_id": registry["registry_id"],
        },
        "source_inventory": {
            "path": registry["source_inventory"],
            "sha256": registry["source_inventory_sha256"],
        },
        "source_scope": {
            "path": registry["source_scope"],
            "sha256": registry["source_scope_sha256"],
        },
        "p0_guard": {
            "path": registry["p0_guard_registry"],
            "sha256": registry["p0_guard_registry_sha256"],
            "validation_summary": registry["p0_guard_validation_summary"],
        },
        "github_about": {
            "path": registry["github_about_receipt"],
            "sha256": registry["github_about_receipt_sha256"],
            "semantic_status": "PASS",
        },
        "manifests": {
            "public_sources": {
                "sha256": registry["public_source_manifest_sha256"],
                "files": registry["public_source_manifest"]["counts"]["files"],
            },
            "github_about_snapshots": {
                "sha256": registry["github_about_snapshot_manifest_sha256"],
                "files": registry["github_about_snapshot_manifest"]["counts"]["files"],
            },
        },
        "allowed_current_verdicts": registry["allowed_current_verdicts"],
        "counts": registry["counts"],
        "gates": registry["gates"],
    }
