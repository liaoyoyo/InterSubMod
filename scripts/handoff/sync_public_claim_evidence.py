#!/usr/bin/env python3
"""Atomically synchronize validated public-claim evidence into the handoff package.

The public-claim registry and validation receipt remain the source of truth.  This
script copies their exact bytes into the handoff package and derives the package
pointer and evidence-manifest metadata from those bytes.  It intentionally does
not build or validate the source registry itself.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any


PACKAGE_RELATIVE = Path("docs/handoff/20260813_完整研究資料與軟體交接_01")
SOURCE_REGISTRY_RELATIVE = Path("research/20260813_public_docs_p0_correction/claim_remediation_registry.json")
SOURCE_RECEIPT_RELATIVE = Path("research/20260813_public_docs_p0_correction/validation_receipt.json")
SOURCE_BROWSER_RECEIPT_RELATIVE = Path(
    "research/20260813_public_docs_p0_correction/html_browser_qa_receipt.json"
)
SOURCE_VALIDATION_RUNNER_RELATIVE = Path(
    "research/20260813_public_docs_p0_correction/scripts/run_claim_validation.py"
)
TARGET_REGISTRY_RELATIVE = Path("evidence/full_claim_registry.json")
TARGET_RECEIPT_RELATIVE = Path("evidence/public_claim_validation_receipt.json")
POINTER_RELATIVE = Path("registries/claim_registry.json")
MANIFEST_RELATIVE = Path("evidence/EVIDENCE_MANIFEST.json")

EXPECTED_RECORDS_URI = "../evidence/full_claim_registry.json"
EXPECTED_RECORD_COUNT = 158
EXPECTED_GATE_STATUSES = {
    "P0_SOURCE_READY": "PASS",
    "SOURCE_READY": "BLOCKED",
    "PUBLICATION_READY": "BLOCKED",
    "RELEASE_READY": "BLOCKED",
}
EVIDENCE_RECORDS = {
    "full_claim_registry_20260813": {
        "path": TARGET_REGISTRY_RELATIVE.as_posix(),
        "source_path": f"InterSubMod/{SOURCE_REGISTRY_RELATIVE.as_posix()}",
    },
    "public_claim_validation_20260813": {
        "path": TARGET_RECEIPT_RELATIVE.as_posix(),
        "source_path": f"InterSubMod/{SOURCE_RECEIPT_RELATIVE.as_posix()}",
    },
}


class SyncError(ValueError):
    """Raised when a source or package contract would be violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SyncError(message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def load_json_bytes(payload: bytes, label: str) -> Any:
    try:
        return json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise SyncError(f"{label} is not valid UTF-8 JSON: {error}") from error


def load_json(path: Path, label: str) -> Any:
    try:
        return load_json_bytes(path.read_bytes(), label)
    except OSError as error:
        raise SyncError(f"cannot read {label}: {path}: {error}") from error


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode("utf-8")


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write_bytes(path: Path, payload: bytes) -> None:
    """Publish bytes with same-directory temporary file, fsync and replace."""

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        fsync_directory(path.parent)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def validate_source_registry(registry: Any) -> tuple[dict[str, Any], dict[str, str]]:
    require(isinstance(registry, dict), "source claim registry must be an object")
    claims = registry.get("claims")
    require(isinstance(claims, list), "source claim registry claims must be a list")
    require(len(claims) == EXPECTED_RECORD_COUNT, f"source claim registry must contain {EXPECTED_RECORD_COUNT} claims")
    claim_ids = [row.get("claim_id") for row in claims if isinstance(row, dict)]
    require(len(claim_ids) == len(claims), "source claim registry contains a non-object claim")
    require(all(isinstance(value, str) and value for value in claim_ids), "claim_id must be non-empty")
    require(len(set(claim_ids)) == EXPECTED_RECORD_COUNT, "source claim registry claim_id values are not unique")

    counts = registry.get("counts")
    require(isinstance(counts, dict), "source claim registry counts must be an object")
    require(counts.get("claims") == EXPECTED_RECORD_COUNT, "source claim registry counts.claims must be 158")
    distributions = {key: value for key, value in counts.items() if key != "claims"}
    require(distributions, "source claim registry contains no count distributions")
    for category, distribution in distributions.items():
        require(isinstance(distribution, dict), f"claim count category must be an object: {category}")
        require(
            all(isinstance(value, int) and not isinstance(value, bool) and value >= 0 for value in distribution.values()),
            f"claim count category contains a non-negative-integer violation: {category}",
        )
        require(sum(distribution.values()) == EXPECTED_RECORD_COUNT, f"claim count category does not sum to 158: {category}")

    gates = registry.get("gates")
    require(isinstance(gates, dict), "source claim registry gates must be an object")
    require(set(gates) == set(EXPECTED_GATE_STATUSES), "source claim registry gate identities differ")
    statuses: dict[str, str] = {}
    for gate, expected in EXPECTED_GATE_STATUSES.items():
        row = gates.get(gate)
        require(isinstance(row, dict), f"source claim gate must be an object: {gate}")
        status = row.get("status")
        require(status == expected, f"source claim gate is not fail-closed: {gate}={status!r}")
        statuses[gate] = status
    return distributions, statuses


def validate_source_receipt(receipt: Any, repo_root: Path, registry_hash: str) -> str:
    require(isinstance(receipt, dict), "source validation receipt must be an object")
    require(
        receipt.get("schema_name") == "intersubmod.public_claim_validation_receipt"
        and receipt.get("schema_version") == "2.0.0",
        "source validation receipt must use the freshness-bound schema 2.0.0",
    )
    require(receipt.get("verdict") == "PASS", "source validation receipt verdict must be PASS")
    publication_status = receipt.get("publication_status")
    require(
        isinstance(publication_status, str) and publication_status.startswith("BLOCKED_"),
        "source validation receipt publication_status must remain BLOCKED",
    )
    require(receipt.get("release_status") == "BLOCKED", "source validation receipt release_status must remain BLOCKED")
    require(
        receipt.get("scope")
        in {
            "LOCAL_SOURCE_PLUS_C108_GITHUB_ABOUT_LIVE_REFETCH",
            "LOCAL_SOURCE_PLUS_C108_LIVE_RECEIPT_PLUS_CHROMIUM_QA",
        },
        "source validation receipt scope does not include the bounded C108 live re-fetch",
    )
    registry_contract = receipt.get("claim_registry_contract")
    browser_contract = receipt.get("browser_qa_contract")
    require(
        isinstance(registry_contract, dict)
        and registry_contract.get("verdict") == "PASS"
        and registry_contract.get("public_source_files") == 34
        and registry_contract.get("registry_sha256") == registry_hash,
        "source validation receipt does not bind the current 34-file claim registry",
    )
    browser_receipt = repo_root / SOURCE_BROWSER_RECEIPT_RELATIVE
    require(browser_receipt.is_file(), "source browser QA receipt is missing")
    require(
        isinstance(browser_contract, dict)
        and browser_contract.get("verdict") == "PASS"
        and browser_contract.get("html_files") == 21
        and browser_contract.get("standalone_svg_files") == 21
        and browser_contract.get("browser_cases") == 84
        and browser_contract.get("receipt_sha256") == sha256_bytes(browser_receipt.read_bytes()),
        "source validation receipt does not bind the exact current 84-case browser QA",
    )
    runner = receipt.get("runner")
    validation_runner = repo_root / SOURCE_VALIDATION_RUNNER_RELATIVE
    require(validation_runner.is_file(), "source validation runner is missing")
    require(
        isinstance(runner, dict)
        and runner.get("path") == SOURCE_VALIDATION_RUNNER_RELATIVE.as_posix()
        and runner.get("sha256") == sha256_bytes(validation_runner.read_bytes())
        and isinstance(runner.get("python"), str)
        and tuple(int(part) for part in runner["python"].split(".")[:2]) >= (3, 10),
        "source validation receipt does not bind the current runner executed with Python >=3.10",
    )
    created_at = receipt.get("created_at")
    require(isinstance(created_at, str) and created_at, "source validation receipt created_at must be non-empty")
    return created_at


def update_pointer(
    pointer: Any,
    *,
    registry_hash: str,
    distributions: dict[str, Any],
    gate_statuses: dict[str, str],
    verified_at: str,
) -> dict[str, Any]:
    require(isinstance(pointer, dict), "handoff claim pointer must be an object")
    require(pointer.get("records_uri") == EXPECTED_RECORDS_URI, "handoff claim pointer records_uri differs")
    updated = dict(pointer)
    updated["records_sha256"] = registry_hash
    updated["records_count"] = EXPECTED_RECORD_COUNT
    updated["counts"] = distributions
    updated["gates"] = gate_statuses
    updated["verified_at"] = verified_at
    materialization = updated.get("materialization")
    require(isinstance(materialization, dict), "handoff claim pointer materialization must be an object")
    materialization = dict(materialization)
    materialization["source_patch_delivery"] = "PUBLIC_CLAIM_CORRECTIONS_WORKING_STACK_VALIDATED"
    materialization["known_limit"] = (
        "The synchronized evidence describes the validated public-claim-corrections working stack. "
        "P0_SOURCE_READY is PASS; P1/P2 closure plus default-branch, Wiki and Pages publication "
        "checks keep SOURCE_READY, PUBLICATION_READY and RELEASE_READY BLOCKED."
    )
    updated["materialization"] = materialization
    return updated


def update_manifest(manifest: Any, payloads: dict[str, bytes]) -> dict[str, Any]:
    require(isinstance(manifest, dict), "handoff evidence manifest must be an object")
    records = manifest.get("records")
    require(isinstance(records, list), "handoff evidence manifest records must be a list")
    updated_records = []
    found: set[str] = set()
    for raw_row in records:
        require(isinstance(raw_row, dict), "handoff evidence manifest contains a non-object record")
        evidence_id = raw_row.get("evidence_id")
        if evidence_id not in EVIDENCE_RECORDS:
            updated_records.append(raw_row)
            continue
        require(evidence_id not in found, f"duplicate evidence record: {evidence_id}")
        expected = EVIDENCE_RECORDS[evidence_id]
        require(raw_row.get("copy_status") == "EXACT_COPY", f"evidence record is not EXACT_COPY: {evidence_id}")
        require(raw_row.get("path") == expected["path"], f"evidence target path differs: {evidence_id}")
        require(raw_row.get("source_path") == expected["source_path"], f"evidence source path differs: {evidence_id}")
        payload = payloads[evidence_id]
        row = dict(raw_row)
        row["sha256"] = sha256_bytes(payload)
        row["size_bytes"] = len(payload)
        updated_records.append(row)
        found.add(evidence_id)
    require(found == set(EVIDENCE_RECORDS), f"handoff evidence manifest lacks records: {sorted(set(EVIDENCE_RECORDS) - found)}")
    updated = dict(manifest)
    updated["records"] = updated_records
    return updated


def sync_evidence(repo_root: Path) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    require(repo_root.is_dir(), f"repository root is not a directory: {repo_root}")
    package_root = repo_root / PACKAGE_RELATIVE
    require(package_root.is_dir(), f"handoff package is missing: {package_root}")

    source_registry_path = repo_root / SOURCE_REGISTRY_RELATIVE
    source_receipt_path = repo_root / SOURCE_RECEIPT_RELATIVE
    try:
        source_registry_bytes = source_registry_path.read_bytes()
        source_receipt_bytes = source_receipt_path.read_bytes()
    except OSError as error:
        raise SyncError(f"cannot read public-claim source evidence: {error}") from error

    source_registry = load_json_bytes(source_registry_bytes, "source claim registry")
    source_receipt = load_json_bytes(source_receipt_bytes, "source validation receipt")
    distributions, gate_statuses = validate_source_registry(source_registry)
    registry_hash = sha256_bytes(source_registry_bytes)
    verified_at = validate_source_receipt(source_receipt, repo_root, registry_hash)

    pointer_path = package_root / POINTER_RELATIVE
    manifest_path = package_root / MANIFEST_RELATIVE
    pointer = load_json(pointer_path, "handoff claim pointer")
    manifest = load_json(manifest_path, "handoff evidence manifest")

    updated_pointer = update_pointer(
        pointer,
        registry_hash=registry_hash,
        distributions=distributions,
        gate_statuses=gate_statuses,
        verified_at=verified_at,
    )
    payloads = {
        "full_claim_registry_20260813": source_registry_bytes,
        "public_claim_validation_20260813": source_receipt_bytes,
    }
    updated_manifest = update_manifest(manifest, payloads)

    # Re-read before publication so a concurrent source rewrite fails closed.
    require(source_registry_path.read_bytes() == source_registry_bytes, "source claim registry changed during synchronization")
    require(source_receipt_path.read_bytes() == source_receipt_bytes, "source validation receipt changed during synchronization")

    target_registry_path = package_root / TARGET_REGISTRY_RELATIVE
    target_receipt_path = package_root / TARGET_RECEIPT_RELATIVE
    atomic_write_bytes(target_registry_path, source_registry_bytes)
    atomic_write_bytes(target_receipt_path, source_receipt_bytes)
    atomic_write_bytes(pointer_path, canonical_json_bytes(updated_pointer))
    atomic_write_bytes(manifest_path, canonical_json_bytes(updated_manifest))

    require(target_registry_path.read_bytes() == source_registry_bytes, "registry exact-copy verification failed")
    require(target_receipt_path.read_bytes() == source_receipt_bytes, "receipt exact-copy verification failed")

    return {
        "schema_version": "1.0.0",
        "operation": "sync_public_claim_evidence",
        "source_registry": SOURCE_REGISTRY_RELATIVE.as_posix(),
        "source_receipt": SOURCE_RECEIPT_RELATIVE.as_posix(),
        "package_root": PACKAGE_RELATIVE.as_posix(),
        "records_count": EXPECTED_RECORD_COUNT,
        "counts": distributions,
        "gates": gate_statuses,
        "outputs": {
            TARGET_REGISTRY_RELATIVE.as_posix(): {
                "size_bytes": len(source_registry_bytes),
                "sha256": registry_hash,
                "exact_copy": True,
            },
            TARGET_RECEIPT_RELATIVE.as_posix(): {
                "size_bytes": len(source_receipt_bytes),
                "sha256": sha256_bytes(source_receipt_bytes),
                "exact_copy": True,
            },
            POINTER_RELATIVE.as_posix(): {"sha256": sha256_bytes(canonical_json_bytes(updated_pointer))},
            MANIFEST_RELATIVE.as_posix(): {"sha256": sha256_bytes(canonical_json_bytes(updated_manifest))},
        },
    }


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=repository_root())
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        result = sync_evidence(arguments.repo_root)
    except (OSError, SyncError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    json.dump(result, sys.stdout, ensure_ascii=False, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
