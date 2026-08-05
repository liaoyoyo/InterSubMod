#!/usr/bin/env python3
"""Verify that frozen all-sSNV inputs were not modified by downstream analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


STAT_ONLY_ROLES = {"raw_alignment", "raw_alignment_index", "latest_read_tag_sidecar"}
AUDITED_ROLES = (
    "all_ssnv_vcf",
    "all_ssnv_vcf_index",
    "raw_alignment",
    "raw_alignment_index",
    "latest_read_tag_sidecar",
    "latest_read_tag_sidecar_index",
    "site_ledger",
    "site_ledger_index",
    "layered_reconstruction",
    "layered_region_view",
    "site_manifest",
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_artifact(role: str, frozen: dict[str, Any]) -> dict[str, Any]:
    path = Path(frozen["path"])
    observed = {
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else None,
        "mtime_ns": path.stat().st_mtime_ns if path.exists() else None,
    }
    checks = {
        "exists": observed["exists"],
        "size_bytes_equal": observed["size_bytes"] == frozen.get("size_bytes"),
        "mtime_ns_equal": observed["mtime_ns"] == frozen.get("mtime_ns"),
    }
    hash_policy = "stat_only_for_large_alignment_artifact"
    if path.exists() and "sha256" in frozen:
        observed["sha256"] = sha256(path)
        checks["sha256_equal"] = observed["sha256"] == frozen["sha256"]
        hash_policy = "sha256"
    elif role not in STAT_ONLY_ROLES:
        checks["frozen_sha256_present"] = False
    return {
        "role": role,
        "path": str(path),
        "hash_policy": hash_policy,
        "frozen": frozen,
        "observed": observed,
        "checks": checks,
        "pass": all(checks.values()),
    }


def main() -> None:
    topic_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=topic_root / "results" / "all_ssnv_input_manifest.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite immutability audit: {args.output}")

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    samples = []
    for entry in manifest["samples"]:
        artifacts = [audit_artifact(role, entry[role]) for role in AUDITED_ROLES]
        samples.append(
            {
                "sample": entry["sample"],
                "artifacts": artifacts,
                "pass": all(artifact["pass"] for artifact in artifacts),
            }
        )
    payload = {
        "schema_name": "intersubmod.all_ssnv_frozen_input_immutability",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "command": sys.argv,
        "source_code": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__)),
        },
        "manifest": {
            "path": str(args.manifest.resolve()),
            "sha256": sha256(args.manifest),
        },
        "samples": samples,
        "totals": {
            "n_samples": len(samples),
            "n_sample_pass": sum(sample["pass"] for sample in samples),
            "n_artifacts": sum(len(sample["artifacts"]) for sample in samples),
            "n_artifact_pass": sum(
                artifact["pass"] for sample in samples for artifact in sample["artifacts"]
            ),
        },
        "pass": len(samples) == 7 and all(sample["pass"] for sample in samples),
        "interpretation": (
            "Large BAM and read-tag sidecar identity uses frozen size plus nanosecond mtime; smaller "
            "manifest artifacts additionally use SHA-256. This verifies no observed input mutation, "
            "not biological correctness of the source files."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "totals": payload["totals"], "pass": payload["pass"]}, indent=2))
    raise SystemExit(0 if payload["pass"] else 1)


if __name__ == "__main__":
    main()
