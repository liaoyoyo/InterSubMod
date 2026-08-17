#!/usr/bin/env python3
"""Replay every hash recorded by the frozen authority manifest.

This verifier proves that frozen inputs are still byte-identical. It does not
claim that the current source tree can regenerate the scientific results.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def authority_records(manifest: dict[str, Any]) -> Iterable[dict[str, str]]:
    for artifact in manifest.get("artifacts", []):
        yield {
            "record_kind": "artifact",
            "artifact_id": str(artifact["artifact_id"]),
            "path": str(artifact["path"]),
            "expected_sha256": str(artifact["sha256"]),
        }

    implementation = manifest.get("implementation", {})
    binary = implementation.get("frozen_binary")
    if binary:
        yield {
            "record_kind": "binary",
            "artifact_id": "frozen_binary",
            "path": str(binary["path"]),
            "expected_sha256": str(binary["sha256"]),
        }

    for index, source in enumerate(implementation.get("source_snapshots", []), start=1):
        yield {
            "record_kind": "source_snapshot",
            "artifact_id": f"source_snapshot_{index}",
            "path": str(source["path"]),
            "expected_sha256": str(source["sha256_at_handoff"]),
        }


def replay(manifest_path: Path) -> dict[str, Any]:
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    results: list[dict[str, Any]] = []
    tally = {"MATCH": 0, "MISSING": 0, "HASH_MISMATCH": 0}

    for record in authority_records(manifest):
        path = Path(record["path"])
        result: dict[str, Any] = dict(record)
        if not path.is_file():
            result.update(status="MISSING", size_bytes=None, actual_sha256=None)
        else:
            actual = sha256_file(path)
            status = "MATCH" if actual == record["expected_sha256"] else "HASH_MISMATCH"
            result.update(status=status, size_bytes=path.stat().st_size, actual_sha256=actual)
        tally[result["status"]] += 1
        results.append(result)

    return {
        "schema_version": "1.0.0",
        "receipt_type": "frozen_authority_hash_replay",
        "claim_ceiling": "Frozen byte-integrity only; not a clean-source scientific rerun.",
        "verified_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "authority_as_of_date": manifest.get("as_of_date"),
        "tally": tally,
        "total": len(results),
        "pass": tally["MISSING"] == 0 and tally["HASH_MISMATCH"] == 0,
        "results": results,
    }


def write_json(path: Path, receipt: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_tsv(path: Path, receipt: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "status",
        "record_kind",
        "artifact_id",
        "size_bytes",
        "expected_sha256",
        "actual_sha256",
        "path",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(receipt["results"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-tsv", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.manifest.is_file():
        print(f"[ERROR] manifest not found: {args.manifest}", file=sys.stderr)
        return 2

    print(f"[INPUT] authority_manifest={args.manifest.resolve()}")
    receipt = replay(args.manifest)
    if args.output_json:
        write_json(args.output_json, receipt)
        print(f"[OUTPUT] receipt_json={args.output_json.resolve()}")
    if args.output_tsv:
        write_tsv(args.output_tsv, receipt)
        print(f"[OUTPUT] receipt_tsv={args.output_tsv.resolve()}")
    print(f"[RESULT] pass={str(receipt['pass']).lower()} tally={json.dumps(receipt['tally'], sort_keys=True)}")
    for result in receipt["results"][:3]:
        print(f"[SAMPLE] {result['status']} {result['artifact_id']} {result['path']}")
    return 0 if receipt["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
