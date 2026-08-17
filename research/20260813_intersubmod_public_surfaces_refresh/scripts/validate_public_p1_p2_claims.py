#!/usr/bin/env python3
"""Validate the frozen P1/P2 public-document correction registry.

This guard checks claim-set closure and bounded wording in local source files.  It
does not assert that GitHub, Wiki, or Pages have already been published.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
from pathlib import Path
from typing import Any


PROBLEM_VERDICTS = {"NEEDS_CORRECTION", "CONTRADICTED", "UNVERIFIABLE"}
TARGET_PRIORITIES = {"P1", "P2"}
ALLOWED_DISPOSITIONS = {
    "local_source_corrected_publication_pending",
    "unsupported_public_claim_removed_publication_pending",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the 24 frozen InterSubMod P1/P2 public-document corrections."
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--inventory",
        type=Path,
        default=Path(
            "research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv"
        ),
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path(
            "research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json"
        ),
    )
    parser.add_argument(
        "--simulate-drop-claim",
        help="Negative-control probe: drop one claim from the loaded registry before validation.",
    )
    return parser.parse_args()


def normalize_document(path: Path) -> str:
    raw = path.read_text(encoding="utf-8")
    if path.suffix.lower() != ".html":
        return raw
    visible = re.sub(r"<script\b[^>]*>.*?</script>", " ", raw, flags=re.I | re.S)
    visible = re.sub(r"<style\b[^>]*>.*?</style>", " ", visible, flags=re.I | re.S)
    visible = re.sub(r"<[^>]+>", " ", visible)
    visible = html.unescape(visible)
    visible = re.sub(r"\s+", " ", visible)
    # Keep both raw attributes/data contracts and reader-visible text searchable.
    return raw + "\n" + visible


def load_inventory(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        row["claim_id"]: row
        for row in rows
        if row["verdict"] in PROBLEM_VERDICTS and row["priority"] in TARGET_PRIORITIES
    }


def compile_regex(pattern: Any, location: str, errors: list[str]) -> re.Pattern[str] | None:
    if not isinstance(pattern, str) or not pattern:
        errors.append(f"{location}: regex must be a non-empty string")
        return None
    try:
        return re.compile(pattern, flags=re.I | re.S)
    except re.error as exc:
        errors.append(f"{location}: invalid regex: {exc}")
        return None


def validate(
    registry: dict[str, Any], inventory: dict[str, dict[str, str]], repo_root: Path
) -> tuple[list[str], dict[str, Any]]:
    errors: list[str] = []
    claims = registry.get("claims")
    if not isinstance(claims, list):
        return ["registry.claims must be an array"], {"verdict": "FAIL", "errors": 1}

    registry_ids = [item.get("claim_id") for item in claims if isinstance(item, dict)]
    if len(registry_ids) != len(set(registry_ids)):
        errors.append("registry contains duplicate claim_id values")
    expected_ids = set(inventory)
    actual_ids = {item for item in registry_ids if isinstance(item, str)}
    if actual_ids != expected_ids:
        missing = sorted(expected_ids - actual_ids)
        extra = sorted(actual_ids - expected_ids)
        errors.append(f"claim-set mismatch: missing={missing}, extra={extra}")

    cache: dict[Path, str] = {}
    checked_rules = 0
    required_anchors = 0
    forbidden_anchors = 0
    dispositions: dict[str, int] = {}

    for claim_index, claim in enumerate(claims):
        location = f"claims[{claim_index}]"
        if not isinstance(claim, dict):
            errors.append(f"{location}: claim must be an object")
            continue
        claim_id = claim.get("claim_id")
        if claim_id not in inventory:
            errors.append(f"{location}: unknown claim_id {claim_id!r}")
            continue
        source = inventory[claim_id]
        if claim.get("source_verdict") != source["verdict"]:
            errors.append(f"{claim_id}: source_verdict does not match inventory")
        if claim.get("priority") != source["priority"]:
            errors.append(f"{claim_id}: priority does not match inventory")
        disposition = claim.get("disposition")
        if disposition not in ALLOWED_DISPOSITIONS:
            errors.append(f"{claim_id}: invalid disposition {disposition!r}")
        else:
            dispositions[disposition] = dispositions.get(disposition, 0) + 1
        if claim.get("live_status") != "publication_pending":
            errors.append(f"{claim_id}: live_status must remain publication_pending")
        if not isinstance(claim.get("evidence"), str) or not claim["evidence"].strip():
            errors.append(f"{claim_id}: evidence must be non-empty")
        if not isinstance(claim.get("external_actions"), list) or not claim["external_actions"]:
            errors.append(f"{claim_id}: external_actions must be a non-empty array")

        checks = claim.get("checks")
        if not isinstance(checks, list) or not checks:
            errors.append(f"{claim_id}: checks must be a non-empty array")
            continue
        for check_index, check in enumerate(checks):
            check_location = f"{claim_id}.checks[{check_index}]"
            if not isinstance(check, dict):
                errors.append(f"{check_location}: check must be an object")
                continue
            target = check.get("target")
            required = check.get("required_all")
            forbidden = check.get("forbidden_any")
            if not isinstance(target, str) or not target:
                errors.append(f"{check_location}: target must be a non-empty path")
                continue
            if not isinstance(required, list) or not required:
                errors.append(f"{check_location}: required_all must be non-empty")
                required = []
            if not isinstance(forbidden, list) or not forbidden:
                errors.append(f"{check_location}: forbidden_any must be non-empty")
                forbidden = []
            path = repo_root / target
            if not path.is_file():
                errors.append(f"{check_location}: target does not exist: {path}")
                continue
            checked_rules += 1
            if path not in cache:
                cache[path] = normalize_document(path)
            text = cache[path]
            for pattern_index, pattern in enumerate(required):
                required_anchors += 1
                regex = compile_regex(pattern, f"{check_location}.required_all[{pattern_index}]", errors)
                if regex is not None and regex.search(text) is None:
                    errors.append(f"{check_location}: missing required anchor {pattern!r}")
            for pattern_index, pattern in enumerate(forbidden):
                forbidden_anchors += 1
                regex = compile_regex(pattern, f"{check_location}.forbidden_any[{pattern_index}]", errors)
                if regex is not None and regex.search(text) is not None:
                    errors.append(f"{check_location}: residual forbidden claim matched {pattern!r}")

    summary = {
        "inventory_problem_p1_p2": len(inventory),
        "registry_problem_p1_p2": len(actual_ids),
        "checked_target_rules": checked_rules,
        "unique_documents": len(cache),
        "required_anchors": required_anchors,
        "forbidden_anchors": forbidden_anchors,
        "dispositions": dispositions,
        "live_status": "publication_pending",
        "errors": len(errors),
        "verdict": "PASS" if not errors else "FAIL",
    }
    return errors, summary


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    inventory_path = args.inventory if args.inventory.is_absolute() else repo_root / args.inventory
    registry_path = args.registry if args.registry.is_absolute() else repo_root / args.registry
    try:
        inventory = load_inventory(inventory_path)
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except (OSError, csv.Error, json.JSONDecodeError) as exc:
        print(json.dumps({"verdict": "FAIL", "errors": [str(exc)]}, ensure_ascii=False, indent=2))
        return 1
    if args.simulate_drop_claim:
        registry["claims"] = [
            claim for claim in registry.get("claims", []) if claim.get("claim_id") != args.simulate_drop_claim
        ]
    errors, summary = validate(registry, inventory, repo_root)
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    if errors:
        print("ERRORS:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
