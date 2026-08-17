#!/usr/bin/env python3
"""Re-fetch public GitHub About evidence and write hash-bound snapshots atomically."""

from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from github_about_evidence import (
    EXPECTED_ACTOR,
    EXPECTED_COMMANDS,
    EXPECTED_DESCRIPTION,
    EXPECTED_REPOSITORY,
    EXPECTED_TOPICS,
    sha256_bytes,
    validate_about_receipt,
)


ROOT = Path(__file__).resolve().parents[3]
PROJECT = ROOT / "research/20260813_public_docs_p0_correction"
EVIDENCE_DIR = PROJECT / "github_about_evidence"
RECEIPT = PROJECT / "github_about_c108_receipt.json"
SNAPSHOTS = {
    "repository": EVIDENCE_DIR / "repository_response.json",
    "topics": EVIDENCE_DIR / "topics_response.json",
    "actor": EVIDENCE_DIR / "actor_response.json",
}


def atomic_write(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def split_include_response(raw: bytes) -> tuple[int, bytes]:
    normalized = raw.replace(b"\r\n", b"\n")
    header, separator, body = normalized.partition(b"\n\n")
    if not separator:
        raise RuntimeError("gh api -i response did not contain an HTTP header separator")
    status_match = re.search(rb"^HTTP/[^ ]+ ([0-9]{3})\b", header, flags=re.MULTILINE)
    if not status_match:
        raise RuntimeError("gh api -i response did not contain an HTTP status")
    return int(status_match.group(1)), body


def fetch(command_text: str, label: str) -> tuple[dict[str, Any], dict[str, Any]]:
    command = command_text.split()
    result = subprocess.run(command, cwd=ROOT, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"{command_text} failed: {result.stderr.decode(errors='replace')}")
    status, body = split_include_response(result.stdout)
    if status != 200:
        raise RuntimeError(f"{command_text} returned HTTP {status}")
    payload = json.loads(body.decode("utf-8"))
    atomic_write(SNAPSHOTS[label], body)
    record = {
        "command": command_text,
        "exit_code": result.returncode,
        "http_status": status,
        "response_snapshot": str(SNAPSHOTS[label].relative_to(ROOT)),
        "response_sha256": sha256_bytes(body),
    }
    return payload, record


def main() -> int:
    payloads: dict[str, Any] = {}
    records: list[dict[str, Any]] = []
    for command, label in EXPECTED_COMMANDS.items():
        payloads[label], record = fetch(command, label)
        records.append(record)
    repository = payloads["repository"]
    topics = payloads["topics"]
    actor = payloads["actor"]
    if repository.get("full_name") != EXPECTED_REPOSITORY:
        raise RuntimeError("repository identity mismatch")
    if repository.get("description") != EXPECTED_DESCRIPTION:
        raise RuntimeError("live repository description does not match bounded C108 wording")
    if set(topics.get("names", [])) != EXPECTED_TOPICS:
        raise RuntimeError("live repository topics do not match bounded C108 topic set")
    if {key: actor.get(key) for key in EXPECTED_ACTOR} != EXPECTED_ACTOR:
        raise RuntimeError("authenticated GitHub actor identity mismatch")

    previous = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {}
    previous_action = previous.get("action", {}) if isinstance(previous, dict) else {}
    now = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    receipt = {
        "schema_name": "intersubmod.github_about_claim_receipt",
        "schema_version": "2.0.0",
        "claim_id": "C108",
        "repository": {
            "id": repository.get("id"),
            "node_id": repository.get("node_id"),
            "name_with_owner": repository.get("full_name"),
            "html_url": repository.get("html_url"),
            "visibility": repository.get("visibility"),
            "default_branch": repository.get("default_branch"),
        },
        "action": {
            "api_route": previous_action.get(
                "api_route", "PATCH /repos/liaoyoyo/InterSubMod and PUT /repos/liaoyoyo/InterSubMod/topics"
            ),
            "executed_at": previous_action.get("executed_at"),
            "verified_at": now,
            "before_description": previous_action.get("before_description"),
            "after_description": EXPECTED_DESCRIPTION,
            "removed_topics": ["phylogeny", "subclone"],
            "verified_topics": sorted(EXPECTED_TOPICS),
        },
        "verification": {
            "actor": {key: actor.get(key) for key in ("login", "id", "node_id")},
            "refetched_at": now,
            "commands": records,
        },
        "evidence_status": "VALIDATED_DERIVED",
        "scope": "FULL",
        "finality": "FINAL_FOR_SCOPE",
        "live_status": "CONFIRMED_WITH_LIMITS_AFTER_REFETCH",
        "claim_ceiling": (
            "The About surface describes local mutation-state candidate analysis and methylation "
            "association only; it does not claim confirmed cellular subclones, biological phylogeny, "
            "production readiness, or release readiness."
        ),
        "publication_effect": (
            "C108 is closed on the GitHub About surface only. Default-branch, Wiki, Pages, and "
            "release gates remain BLOCKED."
        ),
    }
    atomic_write(RECEIPT, (json.dumps(receipt, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    _, errors = validate_about_receipt(ROOT, RECEIPT)
    if errors:
        raise RuntimeError("generated About evidence failed validation: " + "; ".join(errors))
    print(json.dumps({"receipt": str(RECEIPT.relative_to(ROOT)), "commands": records}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
