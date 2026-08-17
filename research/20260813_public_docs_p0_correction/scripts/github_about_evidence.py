#!/usr/bin/env python3
"""Validate hash-bound GitHub About API response snapshots for claim C108."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any


EXPECTED_DESCRIPTION = (
    "ONT read-level somatic-mutation integration for local mutation-state candidate analysis "
    "and methylation association; research software."
)
EXPECTED_TOPICS = {
    "bioinformatics",
    "clustering",
    "dna",
    "haplotype",
    "long-reads",
    "methylation",
    "phasing",
    "read-level",
    "somatic-variants",
    "tumor",
}
EXPECTED_REPOSITORY = "liaoyoyo/InterSubMod"
EXPECTED_ACTOR = {"login": "liaoyoyo", "id": 51171925}
EXPECTED_COMMANDS = {
    "gh api -i repos/liaoyoyo/InterSubMod": "repository",
    "gh api -i repos/liaoyoyo/InterSubMod/topics": "topics",
    "gh api -i user": "actor",
}


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def resolve_snapshot(root: Path, relative: object) -> Path:
    if not isinstance(relative, str) or not relative:
        raise ValueError("response_snapshot must be a non-empty repository-relative path")
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError(f"response snapshot escapes repository: {relative}") from error
    if not candidate.is_file():
        raise ValueError(f"response snapshot is missing: {relative}")
    return candidate


def validate_about_receipt(root: Path, receipt_path: Path) -> tuple[dict[str, Any], list[str]]:
    """Return the receipt and all offline-verifiable semantic/hash errors."""
    errors: list[str] = []
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return {}, [f"GitHub About receipt cannot be parsed: {error}"]

    if receipt.get("schema_name") != "intersubmod.github_about_claim_receipt":
        errors.append("GitHub About receipt schema_name mismatch")
    if receipt.get("schema_version") != "2.0.0":
        errors.append("GitHub About receipt schema_version must be 2.0.0")
    if receipt.get("claim_id") != "C108":
        errors.append("GitHub About receipt claim_id must be C108")
    if receipt.get("live_status") != "CONFIRMED_WITH_LIMITS_AFTER_REFETCH":
        errors.append("GitHub About receipt live_status is not bounded/refetched")

    repository = receipt.get("repository", {})
    action = receipt.get("action", {})
    verification = receipt.get("verification", {})
    if not isinstance(repository, dict) or repository.get("name_with_owner") != EXPECTED_REPOSITORY:
        errors.append("GitHub About repository identity mismatch")
    if not isinstance(action, dict) or action.get("after_description") != EXPECTED_DESCRIPTION:
        errors.append("GitHub About bounded description mismatch")
    if not isinstance(action, dict) or set(action.get("verified_topics", [])) != EXPECTED_TOPICS:
        errors.append("GitHub About bounded topic set mismatch")
    if not isinstance(action, dict) or set(action.get("removed_topics", [])) != {"phylogeny", "subclone"}:
        errors.append("GitHub About receipt must record removal of phylogeny/subclone")
    if not isinstance(verification, dict) or not verification.get("refetched_at"):
        errors.append("GitHub About receipt requires verification.refetched_at")

    commands = verification.get("commands", []) if isinstance(verification, dict) else []
    if not isinstance(commands, list) or {item.get("command") for item in commands if isinstance(item, dict)} != set(EXPECTED_COMMANDS):
        errors.append("GitHub About receipt does not contain the three exact API refetch commands")
        commands = []

    payloads: dict[str, Any] = {}
    for item in commands:
        command = item.get("command")
        label = EXPECTED_COMMANDS.get(command)
        if label is None:
            continue
        if item.get("exit_code") != 0 or item.get("http_status") != 200:
            errors.append(f"GitHub API command was not HTTP 200/exit 0: {command}")
        try:
            snapshot = resolve_snapshot(root, item.get("response_snapshot"))
            raw = snapshot.read_bytes()
            expected_hash = item.get("response_sha256")
            if not re.fullmatch(r"[0-9a-f]{64}", str(expected_hash or "")):
                errors.append(f"invalid API response SHA-256: {command}")
            elif sha256_bytes(raw) != expected_hash:
                errors.append(f"API response snapshot hash mismatch: {command}")
            payloads[label] = json.loads(raw.decode("utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
            errors.append(f"invalid API response snapshot for {command}: {error}")

    repo_payload = payloads.get("repository", {})
    topic_payload = payloads.get("topics", {})
    actor_payload = payloads.get("actor", {})
    if repo_payload:
        if repo_payload.get("full_name") != EXPECTED_REPOSITORY:
            errors.append("repository API snapshot full_name mismatch")
        if repo_payload.get("description") != EXPECTED_DESCRIPTION:
            errors.append("repository API snapshot description mismatch")
        if repo_payload.get("visibility") != "public" or repo_payload.get("private") is not False:
            errors.append("repository API snapshot does not record public visibility")
        if repo_payload.get("default_branch") != "main":
            errors.append("repository API snapshot default_branch mismatch")
    else:
        errors.append("repository API snapshot was not validated")
    if topic_payload:
        if set(topic_payload.get("names", [])) != EXPECTED_TOPICS:
            errors.append("topics API snapshot does not match the bounded topic set")
    else:
        errors.append("topics API snapshot was not validated")
    if actor_payload:
        if {key: actor_payload.get(key) for key in EXPECTED_ACTOR} != EXPECTED_ACTOR:
            errors.append("actor API snapshot identity mismatch")
    else:
        errors.append("actor API snapshot was not validated")

    recorded_actor = verification.get("actor", {}) if isinstance(verification, dict) else {}
    if {key: recorded_actor.get(key) for key in EXPECTED_ACTOR} != EXPECTED_ACTOR:
        errors.append("GitHub About receipt actor identity mismatch")
    return receipt, errors
