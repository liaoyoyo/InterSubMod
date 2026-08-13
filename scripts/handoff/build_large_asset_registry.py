#!/usr/bin/env python3
"""Inventory tracked assets above a size threshold without inferring scientific finality."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime
from pathlib import Path


DEFAULT_THRESHOLD = 1024 * 1024
HARD_GIT_WARNING = 50 * 1024 * 1024


def run(repo: Path, *args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=repo)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact_type(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".html", ".htm"}:
        return "DERIVED_HTML"
    if suffix in {".png", ".jpg", ".jpeg", ".svg", ".pdf"}:
        return "FIGURE_OR_MEDIA"
    if suffix in {".json", ".jsonl", ".tsv", ".csv", ".parquet"}:
        return "DERIVED_TABLE_OR_DATA"
    return "OTHER_TRACKED_ASSET"


def tracked_files(repo: Path) -> list[Path]:
    payload = run(repo, "ls-files", "-z")
    return [repo / os.fsdecode(item) for item in payload.split(b"\0") if item]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold-bytes", type=int, default=DEFAULT_THRESHOLD)
    parser.add_argument("--verified-at", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo = args.repo_root.resolve()
    if args.threshold_bytes < 1:
        raise SystemExit("--threshold-bytes must be positive")
    try:
        datetime.fromisoformat(args.verified_at)
    except ValueError as exc:
        raise SystemExit("--verified-at must be ISO-8601") from exc

    head = run(repo, "rev-parse", "HEAD").decode().strip()
    source_state = (
        "COMMIT_PINNED"
        if subprocess.run(["git", "diff", "--quiet", "HEAD", "--"], cwd=repo, check=False).returncode == 0
        else "WORKTREE_HASH_BOUND_PENDING_COMMIT"
    )
    assets = []
    for path in tracked_files(repo):
        if path.is_symlink() or not path.is_file():
            continue
        size = path.stat().st_size
        if size <= args.threshold_bytes:
            continue
        rel = path.relative_to(repo).as_posix()
        assets.append(
            {
                "artifact_id": "tracked-large-" + hashlib.sha256(rel.encode()).hexdigest()[:16],
                "artifact_type": artifact_type(path),
                "semantic_description": f"Tracked large asset requiring producer/finality/license review: {rel}",
                "sample_id": None,
                "technical_dataset_id": None,
                "mode": "TRACKED_ASSET_AUDIT",
                "genome_build": None,
                "producer": "UNKNOWN_NEEDS_REVIEW",
                "producer_commit": head,
                "schema_version": "large-asset-audit.v1",
                "inputs": [f"git-tracked working-tree file:{rel}"],
                "derived_from": [],
                "validates": [],
                "supersedes": [],
                "used_by": [],
                "created_at": None,
                "generated_at": None,
                "verified_at": args.verified_at,
                "scope": "PARTIAL",
                "evidence_status": "IN_PROGRESS",
                "finality": "NON_FINAL",
                "availability": "GIT",
                "size_bytes": size,
                "sha256": sha256(path),
                "license": "NEEDS_REVIEW",
                "claim_ceiling": "Inventory only; no scientific finality, redistribution, or publication claim.",
                "known_limits": [
                    "Producer, finality, and redistribution license have not yet been adjudicated.",
                    "A tracked path and checksum do not make this artifact a current result."
                ],
                "regeneration_command": "NEEDS_PRODUCER_REVIEW",
                "public_uri": rel,
                "machine_locations": [str(path)],
                "policy_target": (
                    "REMOVE_FROM_REGULAR_GIT_AFTER_ARCHIVE"
                    if size > HARD_GIT_WARNING
                    else "REVIEW_FOR_RELEASE_ASSET_OR_EXPLICIT_ALLOWLIST"
                ),
                "migration_status": "BLOCKED_PENDING_ASSET_ADJUDICATION"
            }
        )

    assets.sort(key=lambda row: (-row["size_bytes"], row["public_uri"]))
    total = sum(row["size_bytes"] for row in assets)
    over_50 = sum(row["size_bytes"] > HARD_GIT_WARNING for row in assets)
    payload = {
        "schema_version": "1.0.0",
        "verified_at": args.verified_at,
        "source_git_head": head,
        "source_state": source_state,
        "threshold_bytes_exclusive": args.threshold_bytes,
        "summary": {
            "tracked_assets_over_threshold": len(assets),
            "total_size_bytes": total,
            "tracked_assets_over_50_mib": over_50,
            "verdict": "PASS" if not assets else "LARGE_ASSET_MIGRATION_BLOCKED"
        },
        "policy": {
            "regular_git_new_file_limit_bytes": 50 * 1024 * 1024,
            "derived_asset_review_threshold_bytes": args.threshold_bytes,
            "interpretation": "Every row remains non-final until producer, license, finality, and publication target are reviewed."
        },
        "assets": assets
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
